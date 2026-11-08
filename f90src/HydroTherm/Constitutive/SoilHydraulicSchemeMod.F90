module SoilHydraulicSchemeMod
  ! Description:
  ! Allocate the simulation-wide hydraulic scheme and adapt existing spatial
  ! parameter arrays to its scalar interface. Selection never changes inventories.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use abortutils, only: endrun,iulog
  use DebugToolMod, only: PrintInfo
  use EcoSIMCtrlMod, only: soil_hydraulic_scheme,DefaultSoilHydraulicScheme
  use EcoSimConst, only: PSIHY
  use SoilPhysDataType, only: FieldCapacity_vr,LOGPOROS_vr,LOGFldCapacity_vr,PSD_vr,FCD_vr, &
    LOGPSIAtSat,LOGPSIFLD_col,LOGPSIMXD_col,LOGPSIMND_col,SRP_vr
  use SoilPropertyDataType, only: POROS_vr
  use SoilWaterDataType, only: PSI_at_Sautration_vr
  use SoilHydraulicModelMod, only: SoilHydraulicModel,SoilHydraulicParameters
  use EcoSIMSoilHydraulicMod, only: EcoSIMSoilHydraulicModel
  implicit none
  private
  class(SoilHydraulicModel), allocatable, target, save :: active_model
  character(len=32), save :: active_name=''
  integer, save :: active_revision=0
  public :: InitializeSoilHydraulicScheme,FinalizeSoilHydraulicScheme,GetSoilHydraulicModel
  public :: GetSoilHydraulicSchemeName,GetSoilHydraulicSchemeRevision,GetSoilHydraulicParameters
contains
!------------------------------------------------------------------------------------------
  subroutine InitializeSoilHydraulicScheme()
    ! Description:
    ! Validate/canonicalize the configured selector and allocate one implementation.
    ! Explicit initialization starts a new scheme lifetime; existing pointers expire.
    character(len=32) :: name
    integer :: i,code
    character(len=*), parameter :: subname='InitializeSoilHydraulicScheme'

    call PrintInfo('beg '//subname)
    name=adjustl(soil_hydraulic_scheme)
    do i=1,len_trim(name)
      code=iachar(name(i:i))
      if(code>=iachar('A').and.code<=iachar('Z'))name(i:i)=achar(code+iachar('a')-iachar('A'))
    enddo
    select case(trim(name))
    case(DefaultSoilHydraulicScheme)
      call FinalizeSoilHydraulicScheme()
      allocate(EcoSIMSoilHydraulicModel :: active_model)
      active_name           = DefaultSoilHydraulicScheme
      active_revision       = 1
      soil_hydraulic_scheme = active_name
    case default
      write(iulog,*) 'Unsupported soil_hydraulic_scheme: ',trim(name),'; available: ecosim'
      call endrun('Unsupported soil hydraulic scheme')
    end select
    call PrintInfo('end '//subname)
  end subroutine InitializeSoilHydraulicScheme
!------------------------------------------------------------------------------------------
  subroutine FinalizeSoilHydraulicScheme()
    ! Description:
    ! Release the selected implementation before explicit reinitialization or cleanup.
    character(len=*), parameter :: subname='FinalizeSoilHydraulicScheme'

    call PrintInfo('beg '//subname)
    if(allocated(active_model))deallocate(active_model)
    active_name     = ''
    active_revision = 0
    call PrintInfo('end '//subname)
  end subroutine FinalizeSoilHydraulicScheme
!------------------------------------------------------------------------------------------
  function GetSoilHydraulicModel() result(model)
    ! Description:
    ! Return the selected implementation. Lazy initialization also supports coupling
    ! adapters and small drivers that do not read the standalone control namelist.
    class(SoilHydraulicModel), pointer :: model ! borrowed simulation-wide implementation [-]
    character(len=*), parameter :: subname='GetSoilHydraulicModel'

    call PrintInfo('beg '//subname)
    if(.not.allocated(active_model))call InitializeSoilHydraulicScheme()
    model=>active_model
    call PrintInfo('end '//subname)
  end function GetSoilHydraulicModel
!------------------------------------------------------------------------------------------
  function GetSoilHydraulicSchemeName() result(name)
    ! Description:
    ! Report the canonical identity of the allocated implementation for output/restarts.
    character(len=32) :: name ! canonical scheme identifier [-]
    character(len=*), parameter :: subname='GetSoilHydraulicSchemeName'

    call PrintInfo('beg '//subname)
    if(.not.allocated(active_model))call InitializeSoilHydraulicScheme()
    name=active_name
    call PrintInfo('end '//subname)
  end function GetSoilHydraulicSchemeName
!------------------------------------------------------------------------------------------
  function GetSoilHydraulicSchemeRevision() result(revision)
    ! Description:
    ! Report the implementation revision used to guard restart cache compatibility.
    integer :: revision ! constitutive implementation revision [-]
    character(len=*), parameter :: subname='GetSoilHydraulicSchemeRevision'

    call PrintInfo('beg '//subname)
    if(.not.allocated(active_model))call InitializeSoilHydraulicScheme()
    revision=active_revision
    call PrintInfo('end '//subname)
  end function GetSoilHydraulicSchemeRevision
!------------------------------------------------------------------------------------------
  function GetSoilHydraulicParameters(L,NY,NX,saturation_tolerance) result(parameters)
    ! Description:
    ! Copy one layer's existing hydraulic parameters/caches into a scalar view.
    ! No per-layer polymorphic objects or duplicate persistent arrays are allocated.
    integer, intent(in) :: L       ! vertical layer index, including litter at zero [-]
    integer, intent(in) :: NY,NX   ! row and column indices [-]
    real(r8), intent(in) :: saturation_tolerance ! near-saturation tolerance [m3 m-3]; unused by inverse/auxiliary methods
    type(SoilHydraulicParameters) :: parameters ! one-layer parameter view with units defined by its type
    character(len=*), parameter :: subname='GetSoilHydraulicParameters'

    call PrintInfo('beg '//subname)
    parameters=SoilHydraulicParameters(                         &
      porosity                 = POROS_vr(L,NY,NX),             &
      field_capacity           = FieldCapacity_vr(L,NY,NX),     &
      log_porosity             = LOGPOROS_vr(L,NY,NX),          &
      log_field_capacity       = LOGFldCapacity_vr(L,NY,NX),    &
      log_wet_content_ratio    = PSD_vr(L,NY,NX),               &
      log_dry_content_ratio    = FCD_vr(L,NY,NX),               &
      log_saturation_potential = LOGPSIAtSat(NY,NX),            &
      log_field_potential      = LOGPSIFLD_col(NY,NX),          &
      log_wet_potential_ratio  = LOGPSIMXD_col(NY,NX),          &
      log_dry_potential_ratio  = LOGPSIMND_col(NY,NX),          &
      curvature                = SRP_vr(L,NY,NX),               &
      saturation_potential     = PSI_at_Sautration_vr(L,NY,NX), &
      dry_potential            = PSIHY,                         &
      saturation_tolerance     = saturation_tolerance,          &
      layer                    = L,                             &
      row                      = NY,                            &
      column                   = NX                             &
      )
    call PrintInfo('end '//subname)
  end function GetSoilHydraulicParameters
end module SoilHydraulicSchemeMod
