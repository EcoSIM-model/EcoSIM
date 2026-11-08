module SoilHydraulicModelMod
  ! Description:
  ! Define the constitutive contract and scalar parameter view for soil hydraulic
  ! schemes. Spatial parameter storage and water inventories remain with EcoSIM.
  use data_kind_mod, only: r8 => DAT_KIND_R8
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
