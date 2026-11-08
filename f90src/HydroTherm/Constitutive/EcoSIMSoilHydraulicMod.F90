module EcoSIMSoilHydraulicMod
  ! Description:
  ! Implement the existing anchored EcoSIM retention curve and Green-Corey table,
  ! including historical auxiliary expressions needed for unchanged default runs.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use DebugToolMod, only: PrintInfo
  use SoilHydraulicModelMod, only: SoilHydraulicModel,SoilHydraulicParameters,NumSoilMoistureClasses
  implicit none
  private
  type, extends(SoilHydraulicModel), public :: EcoSIMSoilHydraulicModel
    ! Description:
    ! Stateless current scheme; all spatial parameters are passed in by the caller.
  contains
    procedure :: matric_potential => EcoSIMMatricPotential
    procedure :: water_content => EcoSIMWaterContent
    procedure :: conductivity_table => EcoSIMConductivityTable
    procedure :: dry_content => EcoSIMDryContent
    procedure :: ksat_reference_content => EcoSIMKsatReferenceContent
    procedure :: equilibrium_content => EcoSIMEquilibriumContent
  end type EcoSIMSoilHydraulicModel
contains
!------------------------------------------------------------------------------------------
  function EcoSIMMatricPotential(this,parameters,theta) result(potential)
    ! Description:
    ! Evaluate the existing piecewise log retention curve without regrouping its arithmetic.
    class(EcoSIMSoilHydraulicModel), intent(in) :: this ! current constitutive scheme [-]
    type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
    real(r8), intent(inout) :: theta ! liquid content; set to porosity near saturation [m3 m-3]
    real(r8) :: potential ! matric potential [MPa]
    character(len=*), parameter :: subname='EcoSIMMatricPotential'

    call PrintInfo('beg '//subname)
    associate(                                               &
      fc            => parameters%field_capacity           , &
      porosity      => parameters%porosity                 , &
      log_fc        => parameters%log_field_capacity       , &
      log_phi       => parameters%log_porosity             , &
      fcd           => parameters%log_dry_content_ratio    , &
      psd           => parameters%log_wet_content_ratio    , &
      log_psifc     => parameters%log_field_potential      , &
      log_psisat    => parameters%log_saturation_potential , &
      log_dry_ratio => parameters%log_dry_potential_ratio. , &
      log_wet_ratio => parameters%log_wet_potential_ratio. , &
      srp           => parameters%curvature                , &
      dry_limit     => parameters%dry_potential            , &
      sat_potential => parameters%saturation_potential     , &
      tolerance     => parameters%saturation_tolerance       &
    )
    IF(theta.LT.fc)THEN
      potential=AMAX1(dry_limit,-EXP(log_psifc+((log_fc-LOG(theta))/fcd*log_dry_ratio)))
    ELSEIF(theta.LT.porosity-tolerance)THEN
      potential=-EXP(log_psisat+(((log_phi-LOG(theta))/psd)**srp*log_wet_ratio))
    ELSE
      theta=porosity
      potential=sat_potential
    ENDIF
    end associate
    call PrintInfo('end '//subname)
  end function EcoSIMMatricPotential
!------------------------------------------------------------------------------------------
  function EcoSIMWaterContent(this,parameters,potential) result(theta)
    ! Description:
    ! Invert both branches of the current retention curve, bounded at saturation
    ! and the dry-pressure limit. This is distinct from its legacy diagnostic approximation.
    class(EcoSIMSoilHydraulicModel), intent(in) :: this ! current constitutive scheme [-]
    type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
    real(r8), intent(in) :: potential ! matric potential [MPa]
    real(r8) :: theta ! bounded liquid content [m3 m-3]
    real(r8) :: bounded_potential
    character(len=*), parameter :: subname='EcoSIMWaterContent'

    call PrintInfo('beg '//subname)
    associate(porosity => parameters%porosity,log_fc => parameters%log_field_capacity, &
      log_phi => parameters%log_porosity,fcd => parameters%log_dry_content_ratio, &
      psd => parameters%log_wet_content_ratio,log_psifc => parameters%log_field_potential, &
      log_psisat => parameters%log_saturation_potential,log_dry_ratio => parameters%log_dry_potential_ratio, &
      log_wet_ratio => parameters%log_wet_potential_ratio,srp => parameters%curvature, &
      dry_limit => parameters%dry_potential,sat_potential => parameters%saturation_potential)
      bounded_potential=AMAX1(potential,dry_limit)
      if(bounded_potential.GE.sat_potential)then
        theta=porosity
      elseif(bounded_potential.GE.-EXP(log_psifc))then
        theta=EXP(log_phi-psd*((LOG(-bounded_potential)-log_psisat)/log_wet_ratio)**(1._r8/srp))
      else
        theta=EXP(log_fc-(LOG(-bounded_potential)-log_psifc)*fcd/log_dry_ratio)
      endif
    end associate
    call PrintInfo('end '//subname)
  end function EcoSIMWaterContent
!------------------------------------------------------------------------------------------
  subroutine EcoSIMConductivityTable(this,parameters,ksat_vertical,ksat_horizontal,water_classes,pressure_classes,conductivity)
    ! Description:
    ! Preserve the original 100-class conductivity sums and multiplication order.
    ! Drainage thresholds are applied by the caller, separately from this scheme.
    class(EcoSIMSoilHydraulicModel), intent(in) :: this ! current constitutive scheme [-]
    type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
    real(r8), intent(in) :: ksat_vertical,ksat_horizontal ! saturated scales [m2 MPa-1 h-1]
    real(r8), intent(out) :: water_classes(NumSoilMoistureClasses) ! liquid samples [m3 m-3]
    real(r8), intent(out) :: pressure_classes(NumSoilMoistureClasses) ! corresponding potentials [MPa]
    real(r8), intent(out) :: conductivity(3,NumSoilMoistureClasses) ! horizontal 1:2, vertical 3 [m2 MPa-1 h-1]
    real(r8) :: xk,yk,sum1,sum2
    integer :: k,m,n,n100
    character(len=*), parameter :: subname='EcoSIMConductivityTable'

    call PrintInfo('beg '//subname)
    associate(porosity => parameters%porosity)
      sum2=0._r8;n100=NumSoilMoistureClasses
      DO k=1,n100
        xk=(k-1._r8)/n100
        water_classes(k)=porosity*(1._r8-xk)
        pressure_classes(k)=this%matric_potential(parameters,water_classes(k))
        sum2=sum2+(2._r8*k-1)/(pressure_classes(k)**2)
      ENDDO
      DO k=1,n100
        xk=(k-1._r8)/n100
        yk=(1._r8-xk)**1.33_r8
        sum1=0._r8
        DO m=k,n100
          sum1=sum1+(2*(m-k)+1._r8)/(pressure_classes(m)**2)
        ENDDO
        DO n=1,3
          IF(n.EQ.3)THEN
            conductivity(n,k)=ksat_vertical*yk*sum1/sum2
          ELSE
            conductivity(n,k)=ksat_horizontal*yk*sum1/sum2
          ENDIF
        ENDDO
      ENDDO
    end associate
    call PrintInfo('end '//subname)
  end subroutine EcoSIMConductivityTable
!------------------------------------------------------------------------------------------
  function EcoSIMDryContent(this,parameters) result(theta)
    ! Description:
    ! Preserve the original dry-bound expression and evaluation order.
    class(EcoSIMSoilHydraulicModel), intent(in) :: this ! current constitutive scheme [-]
    type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
    real(r8) :: theta ! dry-bound liquid content [m3 m-3]
    character(len=*), parameter :: subname='EcoSIMDryContent'

    call PrintInfo('beg '//subname)
    associate(log_psifc => parameters%log_field_potential,dry_limit => parameters%dry_potential, &
      fcd => parameters%log_dry_content_ratio,log_dry_ratio => parameters%log_dry_potential_ratio, &
      log_fc => parameters%log_field_capacity)
      theta=EXP((log_psifc-LOG(-dry_limit))*fcd/log_dry_ratio+log_fc)
    end associate
    call PrintInfo('end '//subname)
  end function EcoSIMDryContent
!------------------------------------------------------------------------------------------
  function EcoSIMKsatReferenceContent(this,parameters) result(theta)
    ! Description:
    ! Preserve the historical -0.033 MPa Ksat reference, which omits inverse SRP.
    class(EcoSIMSoilHydraulicModel), intent(in) :: this ! current constitutive scheme [-]
    type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
    real(r8) :: theta ! reference liquid content [m3 m-3]
    character(len=*), parameter :: subname='EcoSIMKsatReferenceContent'

    call PrintInfo('beg '//subname)
    associate(porosity => parameters%porosity,log_phi => parameters%log_porosity, &
      log_fc => parameters%log_field_capacity,log_psisat => parameters%log_saturation_potential, &
      log_wet_ratio => parameters%log_wet_potential_ratio)
      theta=AMIN1(porosity,EXP((log_psisat-LOG(0.033_r8))*(log_phi-log_fc)/log_wet_ratio+log_phi))
    end associate
    call PrintInfo('end '//subname)
  end function EcoSIMKsatReferenceContent
!------------------------------------------------------------------------------------------
  function EcoSIMEquilibriumContent(this,parameters,potential,max_content) result(theta)
    ! Description:
    ! Preserve the current water-table diagnostic's SRP-free wet-branch extrapolation.
    class(EcoSIMSoilHydraulicModel), intent(in) :: this ! current constitutive scheme [-]
    type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
    real(r8), intent(in) :: potential ! equilibrium matric potential [MPa]
    real(r8), intent(in) :: max_content ! liquid capacity after ice/air allowance [m3 m-3]
    real(r8) :: theta ! bounded diagnostic equilibrium content [m3 m-3]
    character(len=*), parameter :: subname='EcoSIMEquilibriumContent'

    call PrintInfo('beg '//subname)
    associate(log_psisat => parameters%log_saturation_potential,psd => parameters%log_wet_content_ratio, &
      log_wet_ratio => parameters%log_wet_potential_ratio,log_phi => parameters%log_porosity)
      IF(potential.GE.-EXP(log_psisat))THEN
        theta=max_content
      ELSE
        theta=AMIN1(max_content,EXP((log_psisat-LOG(-potential))*psd/log_wet_ratio+log_phi))
      ENDIF
    end associate
    call PrintInfo('end '//subname)
  end function EcoSIMEquilibriumContent
end module EcoSIMSoilHydraulicMod
