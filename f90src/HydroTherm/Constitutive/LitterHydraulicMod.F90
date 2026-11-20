module LitterHydraulicMod
  !
  ! Description:
  ! Own litter retention, conductivity and hydraulic limits independently of the
  ! selected soil model. Adapt existing layer-zero arrays without advancing mass
  ! or energy. Physical pore volume and hydraulic water-holding capacity differ.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use DebugToolMod, only: PrintInfo
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use abortutils, only: endrun,iulog
  use PorousSurfaceResistanceMod, only: PorousSurfaceResistance
  use BETResistanceMod, only: BETResistanceParameters,MakeBETResistanceParameters,BETResistanceState, &
    BETReferencePressure,LitterBETCutoffFCFraction
  use PhysPars, only: litter_bet_c
  use EcoSimConst, only: PSIHY,RGASC,TFice
  use EcoSiMParDataMod, only: micpar
  use SoilPhysDataType, only: FieldCapacity_vr,WiltPoint_vr,LOGPOROS_vr,LOGFldCapacity_vr, &
    LOGWiltPoint_vr,PSD_vr,FCD_vr,LOGPSIAtSat,LOGPSIFLD_col,LOGPSIMXD_col,LOGPSIMND_col, &
    SRP_vr,SatHydroCondVert_vr
  use SoilPropertyDataType, only: POROS_vr,SoilBulkDensity_vr,VLSoilMicPMass_vr,VGeomLayer_vr
  use SoilWaterDataType, only: PSI_at_Sautration_vr,HydroCond_3D,PSISoilAirEntry_vr,ThetaSat_vr,SoilWatAirDry_vr
  use EcoSIMCtrlDataType, only: ZEROS2
  use SurfLitterDataType, only: fLitRPore4WatHold_col,BulkDensLitR
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  integer, parameter, public :: NumLitterMoistureClasses=100
  real(r8), parameter, public :: LitterSaturatedConductivity=10.0_r8*0.098_r8 ! [m2 MPa-1 h-1]

  type(BETResistanceParameters), allocatable, save :: resistance_cache(:,:)
  real(r8), allocatable, save :: resistance_signature(:,:,:) ! native input signatures, (9,row,column)

  type, public :: LitterHydraulicParameters
    ! Description:
    ! Scalar litter material view. Contents use dry bulk litter volume; only
    ! holding_fraction bounds retention. physical_porosity is geometric pore space.
    real(r8) :: holding_fraction           ! retained capacity [m3 water m-3 bulk litter]
    real(r8) :: physical_porosity          ! physical pores [m3 pore m-3 bulk litter]
    real(r8) :: field_capacity             ! field-capacity liquid fraction [m3 m-3]
    real(r8) :: log_holding_fraction       ! ln(holding_fraction) [-]
    real(r8) :: log_field_capacity         ! ln(field_capacity) [-]
    real(r8) :: log_wet_content_ratio      ! ln(holding_fraction/field_capacity) [-]
    real(r8) :: log_dry_content_ratio      ! ln(field_capacity/wilting_point) [-]
    real(r8) :: log_saturation_potential   ! ln(abs(saturation anchor in MPa)) [-]
    real(r8) :: log_field_potential        ! ln(abs(field anchor in MPa)) [-]
    real(r8) :: log_wet_potential_ratio    ! ln(abs(field/saturation anchor)) [-]
    real(r8) :: log_dry_potential_ratio    ! ln(abs(wilting/field anchor)) [-]
    real(r8) :: curvature                  ! wet-branch exponent [-]
    real(r8) :: saturation_potential       ! saturated matric potential [MPa]
    real(r8) :: dry_potential              ! lower matric-potential bound [MPa]
    type(BETResistanceParameters) :: resistance=BETResistanceParameters() ! derived resistance-only material
  end type LitterHydraulicParameters

  public :: GetLitterHydraulicParameters,UpdateLitterHydraulicLimits,LitterHydroproperty
  public :: LitterMatricPotential,LitterWaterContent,LitterMatricPotentialDerivative,LitterDryContent
  public :: LitterConductivityTable,LitterMoistureClass,LitterResistance,LitterResistanceRetention
contains
!------------------------------------------------------------------------------------------
  function GetLitterHydraulicParameters(NY,NX) result(parameters)
  !
  ! Description:
  ! Read native litter parameters directly; never dispatch through a soil model.
  integer, intent(in) :: NY,NX                  ! column indices, row and column [-]
  type(LitterHydraulicParameters) :: parameters ! litter material view
  real(r8) :: signature(9)                      ! native anchors and adjustable adsorption C [mixed units as parameter view]
  character(len=*), parameter :: subname='GetLitterHydraulicParameters'

  call PrintInfo('beg '//subname)
  parameters=LitterHydraulicParameters( &
    fLitRPore4WatHold_col(NY,NX),POROS_vr(0,NY,NX),FieldCapacity_vr(0,NY,NX), &
    LOGPOROS_vr(0,NY,NX),LOGFldCapacity_vr(0,NY,NX),PSD_vr(0,NY,NX),FCD_vr(0,NY,NX), &
    LOGPSIAtSat(NY,NX),LOGPSIFLD_col(NY,NX),LOGPSIMXD_col(NY,NX),LOGPSIMND_col(NY,NX), &
    SRP_vr(0,NY,NX),PSI_at_Sautration_vr(0,NY,NX),PSIHY)
  if(parameters%holding_fraction>0._r8.and.parameters%log_dry_content_ratio>0._r8)then
    if(allocated(resistance_cache))then
      if(any(shape(resistance_cache)/=shape(fLitRPore4WatHold_col))) &
        deallocate(resistance_cache,resistance_signature)
    endif
    if(.not.allocated(resistance_cache))then
      allocate(resistance_cache(size(fLitRPore4WatHold_col,1),size(fLitRPore4WatHold_col,2)))
      allocate(resistance_signature(9,size(fLitRPore4WatHold_col,1),size(fLitRPore4WatHold_col,2)))
      resistance_signature=0._r8
    endif
    signature=[parameters%holding_fraction,parameters%log_field_capacity,parameters%log_dry_content_ratio, &
      parameters%log_field_potential,parameters%log_dry_potential_ratio,parameters%dry_potential, &
      parameters%curvature,litter_bet_c,parameters%field_capacity]
    if(.not.resistance_cache(NY,NX)%ready.or.any(signature/=resistance_signature(:,NY,NX)))then
      resistance_cache(NY,NX)=LitterResistanceRetention(parameters)
      resistance_signature(:,NY,NX)=signature
    endif
    parameters%resistance=resistance_cache(NY,NX)
  endif
  call PrintInfo('end '//subname)
  end function GetLitterHydraulicParameters
!------------------------------------------------------------------------------------------
  subroutine UpdateLitterHydraulicLimits(NY,NX,holding_fraction)
  !
  ! Description:
  ! Refresh the existing litter retention anchors from its water-holding fraction.
  ! Geometry, water/ice inventories and heat capacity are owned by the caller.
  integer, intent(in) :: NY,NX ! column indices, row and column [-]
  real(r8), intent(in) :: holding_fraction ! retained capacity per bulk volume [m3 m-3]
  character(len=*), parameter :: subname='UpdateLitterHydraulicLimits'

  call PrintInfo('beg '//subname)
  fLitRPore4WatHold_col(NY,NX) = holding_fraction
  FieldCapacity_vr(0,NY,NX)   = 0.500_r8*holding_fraction
  WiltPoint_vr(0,NY,NX)       = 0.125_r8*holding_fraction
  LOGPOROS_vr(0,NY,NX)        = LOG(fLitRPore4WatHold_col(NY,NX))
  LOGFldCapacity_vr(0,NY,NX)  = LOG(FieldCapacity_vr(0,NY,NX))
  LOGWiltPoint_vr(0,NY,NX)    = LOG(WiltPoint_vr(0,NY,NX))
  PSD_vr(0,NY,NX)            = LOGPOROS_vr(0,NY,NX)-LOGFldCapacity_vr(0,NY,NX)
  FCD_vr(0,NY,NX)            = LOGFldCapacity_vr(0,NY,NX)-LOGWiltPoint_vr(0,NY,NX)
  SRP_vr(0,NY,NX)            = 1.00_r8
  call PrintInfo('end '//subname)
  end subroutine UpdateLitterHydraulicLimits
!------------------------------------------------------------------------------------------
  function LitterMatricPotential(parameters,theta) result(potential)
  !
  ! Description:
  ! Evaluate the existing piecewise log litter retention curve. Bound saturation
  ! by retained capacity, independently of geometric pore space; dry zero is safe.
  type(LitterHydraulicParameters), intent(in) :: parameters ! litter material view
  real(r8), intent(in) :: theta ! liquid content [m3 water m-3 bulk litter]
  real(r8) :: potential ! matric potential [MPa]
  character(len=*), parameter :: subname='LitterMatricPotential'

  call PrintInfo('beg '//subname)
  if(theta<=0._r8)then
    potential=parameters%dry_potential
  elseif(theta<parameters%field_capacity)then
    potential=MAX(parameters%dry_potential,-EXP(parameters%log_field_potential &
      +((parameters%log_field_capacity-LOG(theta))/parameters%log_dry_content_ratio &
      *parameters%log_dry_potential_ratio)))
  elseif(theta<parameters%holding_fraction)then
    potential=-EXP(parameters%log_saturation_potential &
      +(((parameters%log_holding_fraction-LOG(theta))/parameters%log_wet_content_ratio) &
      **parameters%curvature*parameters%log_wet_potential_ratio))
  else
    potential=parameters%saturation_potential
  endif
  call PrintInfo('end '//subname)
  end function LitterMatricPotential
!------------------------------------------------------------------------------------------
  function LitterWaterContent(parameters,potential) result(theta)
  !
  ! Description:
  ! Invert the litter curve between its dry pressure limit and retained capacity.
  type(LitterHydraulicParameters), intent(in) :: parameters ! litter material view
  real(r8), intent(in) :: potential ! matric potential [MPa]
  real(r8) :: theta ! liquid content [m3 water m-3 bulk litter]
  real(r8) :: bounded_potential
  character(len=*), parameter :: subname='LitterWaterContent'

  call PrintInfo('beg '//subname)
  bounded_potential=MAX(potential,parameters%dry_potential)
  if(bounded_potential>=parameters%saturation_potential)then
    theta=parameters%holding_fraction
  elseif(bounded_potential>=-EXP(parameters%log_field_potential))then
    theta=EXP(parameters%log_holding_fraction-parameters%log_wet_content_ratio &
      *((LOG(-bounded_potential)-parameters%log_saturation_potential) &
      /parameters%log_wet_potential_ratio)**(1._r8/parameters%curvature))
  else
    theta=EXP(parameters%log_field_capacity-(LOG(-bounded_potential)-parameters%log_field_potential) &
      *parameters%log_dry_content_ratio/parameters%log_dry_potential_ratio)
  endif
  call PrintInfo('end '//subname)
  end function LitterWaterContent
!------------------------------------------------------------------------------------------
  function LitterDryContent(parameters) result(theta)
  !
  ! Description:
  ! Retain the legacy dry-branch litter content at the prescribed dry pressure.
  type(LitterHydraulicParameters), intent(in) :: parameters ! litter material view
  real(r8) :: theta     ! air-dry liquid fraction [m3 water m-3 bulk litter]
  character(len=*), parameter :: subname='LitterDryContent'

  call PrintInfo('beg '//subname)
  theta=EXP((parameters%log_field_potential-LOG(-parameters%dry_potential)) &
    *parameters%log_dry_content_ratio/parameters%log_dry_potential_ratio+parameters%log_field_capacity)
  call PrintInfo('end '//subname)
  end function LitterDryContent
!------------------------------------------------------------------------------------------
  function LitterMatricPotentialDerivative(parameters,theta) result(slope)
    ! Description:
    ! Differentiate only the litter curve with a bounded finite difference. Use
    ! one-sided endpoint secants and zero outside the unclipped retention interval.
    type(LitterHydraulicParameters), intent(in) :: parameters ! litter material view
    real(r8), intent(in) :: theta ! liquid fraction [m3 water m-3 bulk litter]
    real(r8) :: slope             ! dpsi/dtheta [MPa per (m3 m-3)]
    real(r8) :: dry,step,lo,hi
    character(len=*), parameter :: subname='LitterMatricPotentialDerivative'

    call PrintInfo('beg '//subname)
    slope = 0._r8
    dry   = MAX(tiny(1._r8),LitterDryContent(parameters))
    if(theta>=dry.and.theta<=parameters%holding_fraction.and.parameters%holding_fraction>dry)then
      step=epsilon(1._r8)**(1._r8/3._r8)*parameters%holding_fraction
      lo=MAX(dry,theta-step);hi=MIN(parameters%holding_fraction,theta+step)
      if(hi>lo)slope=(LitterMatricPotential(parameters,hi)-LitterMatricPotential(parameters,lo))/(hi-lo)
    endif
    call PrintInfo('end '//subname)
  end function LitterMatricPotentialDerivative
!------------------------------------------------------------------------------------------
  subroutine LitterConductivityTable(parameters,ksat,water_classes,pressure_classes,conductivity)
    ! Description:
    ! Build the existing 100-class litter conductivity sums; lateral flow is zero.
    type(LitterHydraulicParameters), intent(in) :: parameters ! litter material view
    real(r8), intent(in) :: ksat ! saturated vertical conductivity [m2 MPa-1 h-1]
    real(r8), intent(out) :: water_classes(NumLitterMoistureClasses) ! liquid fractions [m3 m-3]
    real(r8), intent(out) :: pressure_classes(NumLitterMoistureClasses) ! matric pressures [MPa]
    real(r8), intent(out) :: conductivity(3,NumLitterMoistureClasses) ! directions, classes [m2 MPa-1 h-1]
    integer :: k,m
    real(r8) :: xk,yk,sum1,sum2
    character(len=*), parameter :: subname='LitterConductivityTable'

    call PrintInfo('beg '//subname)
    sum2=0._r8
    do k=1,NumLitterMoistureClasses
      xk                  = (k-1._r8)/NumLitterMoistureClasses
      water_classes(k)    = parameters%holding_fraction*(1._r8-xk)
      pressure_classes(k) = LitterMatricPotential(parameters,water_classes(k))
      sum2                = sum2+(2*k-1)/(pressure_classes(k)**2)
    enddo
    conductivity=0._r8
    do k=1,NumLitterMoistureClasses
      sum1 = 0._r8
      xk   = (k-1._r8)/NumLitterMoistureClasses
      yk   = (1._r8-xk)**1.33_r8
      do m=k,NumLitterMoistureClasses
        sum1=sum1+(2*m+1-2*k)/(pressure_classes(m)**2)
      enddo
      conductivity(3,k)=ksat*yk*sum1/sum2
    enddo
    call PrintInfo('end '//subname)
  end subroutine LitterConductivityTable
!------------------------------------------------------------------------------------------
  function LitterMoistureClass(theta,holding_fraction) result(k)
    ! Description:
    ! Map retained liquid saturation to the existing wet-to-dry table indexing.
    real(r8), intent(in) :: theta ! liquid fraction [m3 water m-3 bulk litter]
    real(r8), intent(in) :: holding_fraction ! retained capacity per bulk volume [m3 m-3]
    integer :: k ! conductivity class, wet=1 and dry=100 [-]
    character(len=*), parameter :: subname='LitterMoistureClass'

    call PrintInfo('beg '//subname)
    k=NumLitterMoistureClasses
    if(holding_fraction>0._r8)k=MAX(1,MIN(NumLitterMoistureClasses, &
      NumLitterMoistureClasses+1-INT(100.0_r8*MAX(0._r8,MIN(holding_fraction,theta))/holding_fraction)))
    call PrintInfo('end '//subname)
  end function LitterMoistureClass
!------------------------------------------------------------------------------------------
  subroutine LitterHydroproperty(NY,NX,use_drainage_threshold)
    ! Description:
    ! Refresh litter density, dry content, conductivity cache and hydraulic limits.
    ! The caller passes its threshold policy explicitly; soil scheme selection is
    ! never consulted. No water, ice or energy inventory is modified.
    integer, intent(in) :: NY,NX ! column indices, row and column [-]
    logical, intent(in) :: use_drainage_threshold ! select the existing 0.1-Ksat threshold policy [-]
    real(r8) :: water_classes(NumLitterMoistureClasses),pressure_classes(NumLitterMoistureClasses)
    type(LitterHydraulicParameters) :: parameters
    integer :: k
    character(len=*), parameter :: subname='LitterHydroproperty'

    call PrintInfo('beg '//subname)
    if(VGeomLayer_vr(0,NY,NX)>ZEROS2(NY,NX))then
      SoilBulkDensity_vr(0,NY,NX)=VLSoilMicPMass_vr(0,NY,NX)/VGeomLayer_vr(0,NY,NX)
    else
      SoilBulkDensity_vr(0,NY,NX)=BulkDensLitR(micpar%k_fine_comp)
    endif
    parameters                = GetLitterHydraulicParameters(NY,NX)
    SoilWatAirDry_vr(0,NY,NX) = LitterDryContent(parameters)
    call LitterConductivityTable(parameters,SatHydroCondVert_vr(0,NY,NX),water_classes, &
      pressure_classes,HydroCond_3D(:,:,0,NY,NX))
    do k=2,NumLitterMoistureClasses
      if(use_drainage_threshold)then
        if(HydroCond_3D(3,k,0,NY,NX)<0.1_r8*SatHydroCondVert_vr(0,NY,NX).and. &
           HydroCond_3D(3,k-1,0,NY,NX)>=0.1_r8*SatHydroCondVert_vr(0,NY,NX))then
          PSISoilAirEntry_vr(0,NY,NX) = pressure_classes(k-1)
          ThetaSat_vr(0,NY,NX)        = water_classes(k-1)
        endif
      elseif(pressure_classes(k)<PSISoilAirEntry_vr(0,NY,NX))then
        if(pressure_classes(k-1)>=PSISoilAirEntry_vr(0,NY,NX))ThetaSat_vr(0,NY,NX)=water_classes(k)
      endif
    enddo
    call PrintInfo('end '//subname)
  end subroutine LitterHydroproperty
!------------------------------------------------------------------------------------------
  function LitterResistanceRetention(parameters) result(retention)
    ! Description:
    ! Derive a cubic/BET material at theta=0.5*field_capacity, independent of
    ! physical pore porosity. The getter caches it and refreshes on material change.
    type(LitterHydraulicParameters), intent(in) :: parameters ! current litter retention anchors
    type(BETResistanceParameters) :: retention ! reconstructible resistance-only material
    real(r8) :: theta,potential,slope,candidate
    character(len=*), parameter :: subname='LitterResistanceRetention'

    call PrintInfo('beg '//subname)
    associate(context => parameters) ! input: native litter anchors on dry bulk-volume basis [-]
    theta=LitterBETCutoffFCFraction*context%field_capacity
    potential=LitterMatricPotential(context,theta)
    slope=LitterMatricPotentialDerivative(context,theta)
    candidate=LitterWaterContent(context,BETReferencePressure(potential))
    retention=MakeBETResistanceParameters(theta,potential,slope,candidate,litter_bet_c)
    end associate
    call PrintInfo('end '//subname)
  end function LitterResistanceRetention
!------------------------------------------------------------------------------------------
  function LitterResistance(parameters,liquid_volume,ice_volume,bulk_volume,litter_area, &
    hydraulic_conductivity,free_vapor_diffusivity,tortuosity_factor,vapor_density,temperature,source_potential) result(resistance)
    ! Description:
    ! Apply Tang-Riley to the exposed litter-atmosphere branch. Retained liquid
    ! supplies capillary transport; actual liquid and physical ice displace pore
    ! air. All volume fractions use dry bulk litter volume. Bulk volume divided
    ! by total litter-covered area gives full patch thickness, without the 2.5-mm
    ! numerical floor or any excess ponded volume. Coverage/time enter flux once.
    ! Mobile liquid beyond holding capacity is treated as a wet surface film.
    ! Water/ice filling all physical pores is treated as exposed free water/ice.
    ! Those limits, and absent litter, have zero internal resistance; atmospheric
    ! film resistance and the accepted liquid/ice donor budgets remain external.
    ! Temperature enables cached cubic/BET resistance-only dry extrapolation.
    ! Global litter_bet_c defaults to 13; the wet cutoff is theta/field_capacity=0.5.
    ! Native litter pressure, hydraulic tables and the vapor flux gradient do not
    ! change. Without temperature, retain the native resistance API behavior.
    type(LitterHydraulicParameters), intent(in) :: parameters ! litter material view
    real(r8), intent(in) :: liquid_volume          ! current inner-step liquid [m3 per cell]
    real(r8), intent(in) :: ice_volume             ! current physical ice, not water equivalent [m3 per cell]
    real(r8), intent(in) :: bulk_volume            ! dry bulk litter volume [m3 per cell]
    real(r8), intent(in) :: litter_area            ! total litter-covered area, including snow-covered part [m2]
    real(r8), intent(in) :: hydraulic_conductivity ! K at retained liquid content [m2 MPa-1 h-1]
    real(r8), intent(in) :: free_vapor_diffusivity ! temperature-adjusted pore-air diffusivity [m2 h-1]
    real(r8), intent(in) :: tortuosity_factor      ! existing dimensionless pore-connectivity factor [-]
    real(r8), intent(in) :: vapor_density          ! common equilibrium source vapor density [Mg H2O m-3]
    real(r8), optional, intent(in) :: temperature ! current litter temperature [K]; enables BET resistance
    real(r8), optional, intent(in) :: source_potential ! matric pressure used for source vapor [MPa]
    real(r8) :: resistance ! internal litter resistance [h m-1]
    real(r8) :: theta,air_fraction,bulk_diffusivity,thickness,slope,bunsen
    real(r8) :: baseline_potential,log_slope,log_theta_rh,log_prefactor,potential
    logical :: extrapolate
    type(BETResistanceParameters) :: retention
    character(len=*), parameter :: subname='LitterResistance'

    call PrintInfo('beg '//subname)
    if(.not.all(ieee_is_finite([liquid_volume,ice_volume,bulk_volume,litter_area, &
      hydraulic_conductivity,free_vapor_diffusivity,tortuosity_factor,vapor_density, &
      parameters%holding_fraction,parameters%physical_porosity])))then
      write(iulog,*)subname,': nonfinite input'
      call endrun('Invalid litter resistance input')
    endif
    if(MIN(liquid_volume,ice_volume,bulk_volume,litter_area,hydraulic_conductivity, &
      free_vapor_diffusivity,tortuosity_factor,vapor_density)<0._r8)then
      write(iulog,*)subname,': negative input'
      call endrun('Invalid litter resistance input')
    endif
    resistance=0._r8
    if(bulk_volume>0._r8.and.litter_area>0._r8)then
      if(parameters%holding_fraction<=0._r8.or.parameters%physical_porosity<=0._r8)then
        write(iulog,*)subname,': invalid holding fraction or physical porosity'
        call endrun('Invalid litter resistance parameters')
      endif
      if(liquid_volume<=parameters%holding_fraction*bulk_volume.and. &
        liquid_volume+ice_volume<parameters%physical_porosity*bulk_volume)then
        theta=MIN(parameters%holding_fraction,liquid_volume/bulk_volume)
        air_fraction=MAX(0._r8,parameters%physical_porosity-(liquid_volume+ice_volume)/bulk_volume)
        bulk_diffusivity=free_vapor_diffusivity*air_fraction**2*tortuosity_factor/parameters%physical_porosity
        thickness=bulk_volume/litter_area
        slope=0._r8
        if(theta>0._r8.and.hydraulic_conductivity>0._r8)slope=LitterMatricPotentialDerivative(parameters,theta)
        ! Native liquid density is 1 Mg m-3; use the same source vapor as the flux.
        bunsen=1._r8/MAX(vapor_density,tiny(1._r8))
        extrapolate=.false.
        if(present(temperature).and.theta>0._r8.and.hydraulic_conductivity>0._r8)then
          if(.not.ieee_is_finite(temperature).or.temperature<=0._r8)call endrun('Invalid BET litter temperature')
          retention=parameters%resistance
          if(.not.retention%ready)retention=LitterResistanceRetention(parameters)
          extrapolate=theta<retention%theta_wet
          if(extrapolate)then
            call BETResistanceState(retention,theta,temperature,potential,log_slope,log_theta_rh)
            baseline_potential=LitterMatricPotential(parameters,theta)
            if(present(source_potential))baseline_potential=source_potential
            if(.not.ieee_is_finite(baseline_potential))call endrun('Invalid BET litter source potential')
            if(ice_volume>0._r8.and.temperature<TFice)then
              log_prefactor=LOG(bunsen)+LOG(theta)+log_slope
            else
              log_prefactor=LOG(bunsen)+baseline_potential*18._r8/(RGASC*temperature)+log_theta_rh+log_slope
            endif
          endif
        endif
        if(extrapolate)then
          resistance=PorousSurfaceResistance(theta,thickness,hydraulic_conductivity,0._r8, &
            bulk_diffusivity,bunsen,log_prefactor)
        else
          resistance=PorousSurfaceResistance(theta,thickness,hydraulic_conductivity,slope,bulk_diffusivity,bunsen)
        endif
      endif
    endif
    call PrintInfo('end '//subname)
  end function LitterResistance
end module LitterHydraulicMod
