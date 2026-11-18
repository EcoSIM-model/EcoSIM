module SoilHydraulicModelMod
  ! Description:
  ! Define the constitutive contract and scalar parameter view for soil hydraulic
  ! schemes, including numerical retention slopes and surface soil resistance.
  ! Spatial parameter storage and water inventories remain with EcoSIM.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use abortutils, only: endrun,iulog
  use DebugToolMod, only: PrintInfo
  implicit none
  private
  integer, parameter, public :: NumSoilMoistureClasses=100

  type, public :: SoilHydraulicParameters
    ! Description:
    ! View of one layer's existing parameters, including cached logarithms so the
    ! current scheme retains its arithmetic. Contents are per micropore soil volume.
    ! Other schemes may own additional parameter arrays indexed by layer/row/column,
    ! without adding their individual coefficients to this shared parameter view.
    real(r8) :: porosity                    ! pore capacity [m3 m-3]
    real(r8) :: field_capacity              ! field-capacity liquid content [m3 m-3]
    real(r8) :: log_porosity                ! ln(porosity) [-]
    real(r8) :: log_field_capacity          ! ln(field_capacity) [-]
    real(r8) :: log_wet_content_ratio       ! ln(porosity/field_capacity) [-]
    real(r8) :: log_dry_content_ratio       ! ln(field_capacity/wilting_point) [-]
    real(r8) :: log_saturation_potential    ! ln(abs(saturation anchor in MPa)) [-]
    real(r8) :: log_field_potential         ! ln(abs(field-capacity anchor in MPa)) [-]
    real(r8) :: log_wet_potential_ratio     ! ln(abs(field/saturation anchor)) [-]
    real(r8) :: log_dry_potential_ratio     ! ln(abs(wilting/field anchor)) [-]
    real(r8) :: curvature                   ! wet-branch exponent [-]
    real(r8) :: saturation_potential        ! returned potential at saturation [MPa]
    real(r8) :: dry_potential               ! lower matric-potential bound [MPa]
    real(r8) :: saturation_tolerance        ! near-saturation content tolerance [m3 m-3]
    integer :: layer=0,row=0,column=0       ! spatial indices for scheme-specific parameter arrays [-]
  end type SoilHydraulicParameters

  type, abstract, public :: SoilHydraulicModel
    ! Description:
    ! A paired retention/conductivity implementation. The inverse and auxiliary
    ! methods distinguish exact constitutive inversion from historical diagnostics.
  contains
    procedure(matric_interface), deferred :: matric_potential
    procedure :: matric_potential_derivative => MatricPotentialDerivative
    procedure :: soil_resistance             => SoilResistance
    procedure(inverse_interface), deferred :: water_content
    procedure(table_interface), deferred :: conductivity_table
    procedure :: relative_conductivity  => LookupRelativeConductivity
    procedure :: dry_content            => DefaultDryContent
    procedure :: ksat_reference_content => DefaultKsatReferenceContent
    procedure :: equilibrium_content    => DefaultEquilibriumContent
  end type SoilHydraulicModel

  abstract interface
    function matric_interface(this,parameters,theta) result(potential)
      ! Description:
      ! Evaluate matric potential; retain the scheme's saturation adjustment of theta.
      import :: SoilHydraulicModel,SoilHydraulicParameters,r8
      class(SoilHydraulicModel), intent(in) :: this ! constitutive implementation [-]
      type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
      real(r8), intent(inout) :: theta ! liquid content; may be adjusted at saturation [m3 m-3]
      real(r8) :: potential ! matric potential [MPa]
    end function matric_interface

    function inverse_interface(this,parameters,potential) result(theta)
      ! Description:
      ! Invert the retention relation within its saturated and dry limits.
      import :: SoilHydraulicModel,SoilHydraulicParameters,r8
      class(SoilHydraulicModel), intent(in) :: this ! constitutive implementation [-]
      type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
      real(r8), intent(in) :: potential ! matric potential [MPa]
      real(r8) :: theta ! bounded liquid content [m3 m-3]
    end function inverse_interface

    subroutine table_interface(this,parameters,ksat_vertical,ksat_horizontal,water_classes,pressure_classes,conductivity)
      ! Description:
      ! Build a scheme's paired pressure/conductivity samples at the fixed moisture classes.
      ! Class k must use theta=porosity*(1-(k-1)/100), for k=1..100. Existing consumers
      ! index total liquid saturation theta/porosity, not residual-adjusted saturation.
      import :: SoilHydraulicModel,SoilHydraulicParameters,r8,NumSoilMoistureClasses
      class(SoilHydraulicModel), intent(in) :: this                     ! constitutive implementation [-]
      type(SoilHydraulicParameters), intent(in) :: parameters           ! one-layer parameter view
      real(r8), intent(in) :: ksat_vertical,ksat_horizontal             ! saturated scales [m2 MPa-1 h-1]
      real(r8), intent(out) :: water_classes(NumSoilMoistureClasses)    ! wet-to-dry contents [m3 m-3]
      real(r8), intent(out) :: pressure_classes(NumSoilMoistureClasses) ! matric potentials [MPa]
      real(r8), intent(out) :: conductivity(3,NumSoilMoistureClasses)   ! directions 1:2 horizontal, 3 vertical [m2 MPa-1 h-1]
    end subroutine table_interface
  end interface
contains
!------------------------------------------------------------------------------------------
  function MatricPotentialDerivative(this,parameters,theta) result(dpsi_dtheta)
    ! Description:
    ! Differentiate the selected retention curve with a bounded finite difference.
    ! Use centered samples in the interior and one-sided samples at the dry and
    ! saturated endpoints. Return zero outside that interval, where pressure is
    ! clipped. A stencil crossing a branch join returns its interval-mean slope;
    ! an endpoint slope is a finite secant, not a potentially singular limit.
    ! Never expose matric_potential's inout adjustment to the caller or denominator.
    class(SoilHydraulicModel), intent(in) :: this           ! selected constitutive implementation [-]
    type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
    real(r8), intent(in) :: theta                           ! liquid content [m3 water m-3 micropore soil]
    real(r8) :: dpsi_dtheta                                 ! retention slope [MPa per (m3 water m-3 micropore soil)]
    real(r8) :: theta_min,theta_lo,theta_hi,theta_work,step
    real(r8) :: psi_lo,psi_hi
    character(len=*), parameter :: subname='MatricPotentialDerivative'

    call PrintInfo('beg '//subname)
    dpsi_dtheta=0._r8
    theta_min=MAX(tiny(1._r8),this%dry_content(parameters))
    if(theta>=theta_min.and.theta<=parameters%porosity.and.parameters%porosity>theta_min)then
      ! Balance central-difference truncation and roundoff on the content scale.
      step=epsilon(1._r8)**(1._r8/3._r8)*parameters%porosity
      theta_lo=MAX(theta_min,theta-step)
      theta_hi=MIN(parameters%porosity,theta+step)
      if(theta_hi>theta_lo)then
        theta_work=theta_lo
        psi_lo=this%matric_potential(parameters,theta_work)
        theta_work=theta_hi
        psi_hi=this%matric_potential(parameters,theta_work)
        dpsi_dtheta=(psi_hi-psi_lo)/(theta_hi-theta_lo)
      endif
    endif
    call PrintInfo('end '//subname)
  end function MatricPotentialDerivative
!------------------------------------------------------------------------------------------
  function SoilResistance(this,parameters,theta,thickness,hydraulic_conductivity, &
    bulk_vapor_diffusivity,bunsen_coefficient) result(resistance)
    !
    ! Description:
    ! Compute the soil-internal resistance of Tang and Riley (2013), Eqs. 13-14,
    ! doi:10.5194/hess-17-873-2013. Vapor and liquid transport act in parallel:
    ! r_s=(thickness/2)/(epsilon*D_g+B*theta*K*dpsi/dtheta).
    ! Supply epsilon*D_g as bulk_vapor_diffusivity, including air-filled porosity
    ! and tortuosity exactly once. B=rho_liquid/rho_soil_vapor uses equal density
    ! units. Native K [m2 MPa-1 h-1] times the pressure slope [MPa] is D_w [m2 h-1];
    ! no conversion to hydraulic head is needed. K must match the selected curve.
    ! No atmospheric/film resistance, empirical floor or water flux is included.
    ! Zero thickness has zero resistance; blocked transport returns HUGE, and
    ! resistances below the normal floating-point range return zero. Evaluate
    ! conductance sums in log space because B can be large in very dry soil.
    !
    class(SoilHydraulicModel), intent(in) :: this           ! selected retention implementation [-]
    type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
    real(r8), intent(in) :: theta                           ! liquid content [m3 water m-3 micropore soil]
    real(r8), intent(in) :: thickness                       ! topsoil control-volume thickness [m]
    real(r8), intent(in) :: hydraulic_conductivity          ! current vertical K [m2 MPa-1 h-1]
    real(r8), intent(in) :: bulk_vapor_diffusivity          ! epsilon*D_g on the same soil-volume basis [m2 h-1]
    real(r8), intent(in) :: bunsen_coefficient              ! liquid/soil-vapor density ratio B [-]
    real(r8) :: resistance                                  ! soil-internal resistance [h m-1]
    real(r8) :: slope,log_liquid,log_vapor,log_transport,log_resistance
    logical :: liquid_transport
    character(len=*), parameter :: subname='SoilResistance'

    call PrintInfo('beg '//subname)
    if(.not.all(ieee_is_finite([theta,thickness,hydraulic_conductivity, &
      bulk_vapor_diffusivity,bunsen_coefficient,parameters%porosity])))then
      write(iulog,*) subname,': nonfinite soil-resistance input'
      call endrun('Invalid soil resistance input')
    endif

    if(parameters%porosity<=0._r8 .or. theta<0._r8 .or. theta>parameters%porosity .or. &
      thickness<0._r8.or.hydraulic_conductivity<0._r8.or. &
      bulk_vapor_diffusivity<0._r8.or.bunsen_coefficient<0._r8)then
      write(iulog,*) subname,': invalid theta, thickness, K, vapor diffusivity, B or porosity', &
        theta,thickness,hydraulic_conductivity,bulk_vapor_diffusivity,bunsen_coefficient,parameters%porosity
      call endrun('Invalid soil resistance input')
    endif

    resistance=0._r8
    if(thickness>0._r8)then
      slope=0._r8
      if(theta>0._r8.and.hydraulic_conductivity>0._r8.and.bunsen_coefficient>0._r8)then
        slope=this%matric_potential_derivative(parameters,theta)
        if(.not.ieee_is_finite(slope).or.slope<0._r8)then
          write(iulog,*) subname,': invalid retention slope',theta,slope
          call endrun('Invalid soil resistance retention slope')
        endif
      endif
      liquid_transport=slope>0._r8
      if(liquid_transport)log_liquid=LOG(bunsen_coefficient)+LOG(theta) &
        +LOG(hydraulic_conductivity)+LOG(slope)
      if(bulk_vapor_diffusivity>0._r8)then
        log_vapor=LOG(bulk_vapor_diffusivity)
        log_transport=log_vapor
        if(liquid_transport)log_transport=MAX(log_liquid,log_vapor) &
          +LOG(1._r8+EXP(-ABS(log_liquid-log_vapor)))
      elseif(liquid_transport)then
        log_transport=log_liquid
      endif

      resistance=huge(1._r8)
      if(bulk_vapor_diffusivity>0._r8.or.liquid_transport)then
        log_resistance=LOG(thickness)-LOG(2._r8)-log_transport
        if(log_resistance<LOG(tiny(1._r8)))then
          resistance=0._r8
        elseif(log_resistance<LOG(huge(1._r8)))then
          resistance=EXP(log_resistance)
        endif
      endif
    endif
    call PrintInfo('end '//subname)
  end function SoilResistance
!------------------------------------------------------------------------------------------
  function LookupRelativeConductivity(this,parameters,theta) result(relative_k)
    ! Description:
    ! Select relative conductivity from a unit-scale table using EcoSIM's existing
    ! discrete lookup. Runtime flow continues to use its precomputed layer tables.
    class(SoilHydraulicModel), intent(in) :: this           ! constitutive implementation [-]
    type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
    real(r8), intent(in) :: theta                           ! liquid content [m3 m-3]
    real(r8) :: relative_k                                  ! vertical K/Ksat at the selected class [-]
    real(r8) :: water(NumSoilMoistureClasses),pressure(NumSoilMoistureClasses)
    real(r8) :: conductivity(3,NumSoilMoistureClasses)
    integer :: k
    character(len=*), parameter :: subname='LookupRelativeConductivity'

    call PrintInfo('beg '//subname)
    associate(porosity => parameters%porosity)
      call this%conductivity_table(parameters,1._r8,1._r8,water,pressure,conductivity)
      k=MAX(1,MIN(NumSoilMoistureClasses,NumSoilMoistureClasses+1-INT(NumSoilMoistureClasses*theta/porosity)))
      relative_k=conductivity(3,k)
    end associate
    call PrintInfo('end '//subname)
  end function LookupRelativeConductivity
!------------------------------------------------------------------------------------------
  function DefaultDryContent(this,parameters) result(theta)
    ! Description:
    ! Use a scheme's inverse at the dry-potential bound unless it overrides this rule.
    class(SoilHydraulicModel), intent(in) :: this           ! constitutive implementation [-]
    type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
    real(r8) :: theta                                       ! dry-limit liquid content [m3 m-3]
    character(len=*), parameter :: subname='DefaultDryContent'

    call PrintInfo('beg '//subname)
    associate(dry_potential => parameters%dry_potential)
      theta=this%water_content(parameters,dry_potential)
    end associate
    call PrintInfo('end '//subname)
  end function DefaultDryContent
!------------------------------------------------------------------------------------------
  function DefaultKsatReferenceContent(this,parameters) result(theta)
    ! Description:
    ! Obtain the -0.033 MPa reference content from a scheme's inverse by default.
    class(SoilHydraulicModel), intent(in) :: this           ! constitutive implementation [-]
    type(SoilHydraulicParameters), intent(in) :: parameters ! one-layer parameter view
    real(r8) :: theta                                       ! reference liquid content [m3 m-3]
    character(len=*), parameter :: subname='DefaultKsatReferenceContent'

    call PrintInfo('beg '//subname)
    theta=this%water_content(parameters,-.033_r8)
    call PrintInfo('end '//subname)
  end function DefaultKsatReferenceContent
!------------------------------------------------------------------------------------------
  function DefaultEquilibriumContent(this,parameters,potential,max_content) result(theta)
    ! Description:
    ! Bound inverse retention by the liquid capacity supplied by water-table diagnosis.
    class(SoilHydraulicModel), intent(in) :: this             ! constitutive implementation [-]
    type(SoilHydraulicParameters), intent(in) :: parameters   ! one-layer parameter view
    real(r8), intent(in) :: potential                         ! equilibrium matric potential [MPa]
    real(r8), intent(in) :: max_content                       ! liquid capacity after ice/air allowance [m3 m-3]
    real(r8) :: theta                                         ! bounded equilibrium content [m3 m-3]
    character(len=*), parameter :: subname='DefaultEquilibriumContent'

    call PrintInfo('beg '//subname)
    theta=AMIN1(max_content,this%water_content(parameters,potential))
    call PrintInfo('end '//subname)
  end function DefaultEquilibriumContent
end module SoilHydraulicModelMod
