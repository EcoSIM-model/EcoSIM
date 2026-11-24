module ATSStateRegistryMod
  !
  ! Single list of the EcoSIM-private data that ATS stores on EcoSIM's behalf:
  !
  !   kRolePrivate: state carried between advances; ATS checkpoints it so that
  !                 restart is exact. Unpacked at the start of every advance.
  !   kRoleOutput:  EcoSIM-only outputs (not used by ATS physics); ATS
  !                 visualizes them. Packed after every advance, never unpacked.
  !
  ! Each entry maps the ATS field name (ATS key "surface-<ats_name>") to the
  ! EcoSIM variable it holds, with units and number of components per column.
  ! Units are copied from the EcoSIM declarations ("d-2" = per grid cell area).
  ! This list is the only place these variables are named: pack, unpack and
  ! the layout reported to ATS all loop over it. Adding, removing, reordering
  ! or resizing an entry changes the checkpoint layout, so bump kLayoutVersion.
  !
  ! Description:
  ! Register native carried state.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use, intrinsic :: iso_c_binding
  use EcoContainers_module, only : EcoSizes, EcoInternalState
  use abortutils, only : endrun
  implicit none
  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

  integer, parameter, public :: kRolePrivate = 0
  integer, parameter, public :: kRoleOutput = 1
  integer, parameter, public :: kLayoutVersion = 2

  integer, parameter, public :: kNameLen = 64
  integer, parameter, public :: kUnitsLen = 32
  integer, parameter, public :: kDescLen = 128

  type :: state_entry
    character(len=kNameLen)  :: ats_name
    character(len=kNameLen)  :: ecosim_name
    character(len=kUnitsLen) :: units
    character(len=kDescLen)  :: description
    integer :: ncomp
    integer :: role
    real(r8), pointer :: p(:,:) => null()   ! (ncomp, ncol) view of the EcoSIM data
  end type state_entry

  type(state_entry), allocatable :: entries(:)
  integer :: nentries = 0
  logical :: bind_arrays = .false.
  integer :: ncol_bind = 0

  public :: NumStateEntries
  public :: GetStateEntry
  public :: PackInternalState
  public :: UnpackInternalState

contains
!------------------------------------------------------------------------------------------

  subroutine DefineEntries(sizes, bind)
  ! The list. With bind = .false. only the layout is built, so it can be
  ! queried before EcoSIM allocates its arrays.
  use GridConsts,         only : JS, JP
  use ElmIDMod,           only : NumPlantChemElms
  use EcoSIMConfig,       only : jcplx => jcplxc, NumLiveMicrbCompts
  use SnowDataType,       only : VLDrySnoWE_snvr, VLWatSnow_snvr, VLIceSnow_snvr, &
    TKSnow_snvr, TCSnow_snvr, SnoDens_snvr, SnowThickL_snvr, VLSnoDWIprev_snvr, &
    VLHeatCapSnow_snvr
  use ChemTranspDataType, only : H2OVapDifsc_snvr
  use SoilWaterDataType,  only : VLWatMicP_vr, VLiceMicP_vr
  use SoilHeatDataType,   only : TKS_vr, VHeatCapacity_vr
  use CanopyDataType,     only : WatHeldOnCanopy_pft, &
    LWRadCanGPrev_col, WatHeldOnCanopy_col, PSICanopy_pft, TKC_pft, TKCanopy_pft, &
    ENGYX_pft, DeltaTKC_pft, CanopyBiomWater_pft, Transpiration_pft
  use ClimForcDataType,   only : TLEX_col, TSHX_col
  use SharedDataMod,      only : a_Transpiration, a_EvapCan, a_EvapGrnd, a_EvapLitr, &
    a_EvapSnow, a_Sublim, a_PlantInit
  use MicrobialDataType,  only : mBiomeHeter_vr, mBiomeAutor_vr
  implicit none
  type(EcoSizes), intent(in) :: sizes
  logical, intent(in) :: bind
  integer :: npft, nHeterBioms, nAutorBioms, nlay

  npft = min(sizes%num_pfts, JP)
  nHeterBioms = NumLiveMicrbCompts * sizes%num_heter_guilds
  nAutorBioms = NumLiveMicrbCompts * sizes%num_autor_guilds
  nlay = sizes%ncells_per_col_ + 1  ! +1 for the surface litter layer (index 0)
  bind_arrays = bind
  ncol_bind = sizes%num_columns
  if (allocated(entries)) deallocate(entries)
  allocate(entries(80))   ! capacity; new_entry stops the run if exceeded
  nentries = 0

  ! state carried between advances (checkpointed by ATS)
  !         ATS name                      EcoSIM variable       units       description, array, first index, components
  call add3('snow_layer_dry_swe', 'VLDrySnoWE_snvr', 'm3 d-2', kRolePrivate, &
    'dry snow water equivalent per snow layer', VLDrySnoWE_snvr, 1, JS)
  call add3('snow_layer_liquid', 'VLWatSnow_snvr', 'm3 d-2', kRolePrivate, &
    'liquid water per snow layer', VLWatSnow_snvr, 1, JS)
  call add3('snow_layer_ice', 'VLIceSnow_snvr', 'm3 d-2', kRolePrivate, &
    'ice per snow layer', VLIceSnow_snvr, 1, JS)
  call add3('snow_layer_temperature', 'TKSnow_snvr', 'K', kRolePrivate, &
    'temperature per snow layer', TKSnow_snvr, 1, JS)
  call add3('snow_layer_temperature_c', 'TCSnow_snvr', 'C', kRolePrivate, &
    'temperature per snow layer (kept separately, can differ from TKSnow after melt-out)', TCSnow_snvr, 1, JS)
  call add3('snow_layer_density', 'SnoDens_snvr', 'Mg m-3', kRolePrivate, &
    'density per snow layer', SnoDens_snvr, 1, JS)
  call add3('snow_layer_thickness', 'SnowThickL_snvr', 'm', kRolePrivate, &
    'thickness per snow layer', SnowThickL_snvr, 1, JS)
  call add3('snow_layer_volume', 'VLSnoDWIprev_snvr', 'm3 d-2', kRolePrivate, &
    'volume per snow layer', VLSnoDWIprev_snvr, 1, JS)
  call add3('snow_layer_heat_capacity', 'VLHeatCapSnow_snvr', 'MJ m-3 K-1', kRolePrivate, &
    'heat capacity per snow layer', VLHeatCapSnow_snvr, 1, JS)
  call add3('snow_layer_vapor_diffusivity', 'H2OVapDifsc_snvr', 'm2 h-1', kRolePrivate, &
    'water vapor diffusivity per snow layer', H2OVapDifsc_snvr, 1, JS)
  call add3('litter_water', 'VLWatMicP_vr(0)', 'm3 d-2', kRolePrivate, &
    'surface litter (soil layer 0) liquid water', VLWatMicP_vr, 0, 1)
  call add3('litter_ice', 'VLiceMicP_vr(0)', 'm3 d-2', kRolePrivate, &
    'surface litter (soil layer 0) ice', VLiceMicP_vr, 0, 1)
  call add3('litter_temperature', 'TKS_vr(0)', 'K', kRolePrivate, &
    'surface litter (soil layer 0) temperature', TKS_vr, 0, 1)
  call add3('litter_heat_capacity', 'VHeatCapacity_vr(0)', 'MJ m-3 K-1', kRolePrivate, &
    'surface litter (soil layer 0) heat capacity', VHeatCapacity_vr, 0, 1)
  call add3('canopy_water_pft', 'WatHeldOnCanopy_pft', 'm3 d-2', kRolePrivate, &
    'water held on the canopy per PFT', WatHeldOnCanopy_pft, 1, npft)
  ! canopy snow (SnowOnCanopy_pft) is exchanged with ATS through
  ! feedback%canopy_snow (a_CanSnow), so it is not registered here
  ! Prescribed phenology plant water and energy state, carried between hours
  ! (set at the first advance by InitPrescribedPlantState; the flag keeps a
  ! restarted run from initializing it again)
  call add2('plant_state_initialized', 'a_PlantInit', '-', kRolePrivate, &
    'prescribed phenology: 1 once the plant state below has been initialized', a_PlantInit)
  call add3('canopy_water_potential', 'PSICanopy_pft', 'MPa', kRolePrivate, &
    'canopy total water potential per PFT', PSICanopy_pft, 1, npft)
  call add3('plant_water', 'CanopyBiomWater_pft', 'm3 d-2', kRolePrivate, &
    'water in the plant tissue (canopy dry matter) per PFT', CanopyBiomWater_pft, 1, npft)
  call add3('canopy_temperature', 'TKC_pft', 'K', kRolePrivate, &
    'canopy temperature after the energy iteration per PFT', TKC_pft, 1, npft)
  call add3('canopy_temperature_start', 'TKCanopy_pft', 'K', kRolePrivate, &
    'canopy temperature that starts the next energy iteration per PFT', TKCanopy_pft, 1, npft)
  call add3('canopy_heat_storage', 'ENGYX_pft', 'MJ d-2', kRolePrivate, &
    'canopy heat storage from the previous step per PFT', ENGYX_pft, 1, npft)
  call add3('canopy_temperature_change', 'DeltaTKC_pft', 'K', kRolePrivate, &
    'canopy minus air temperature from the previous step per PFT', DeltaTKC_pft, 1, npft)
  ! The next three are EcoSIM carry-over values, not fluxes for ATS to use:
  ! they are what EcoSIM needs from the previous step to start the next one.
  call add2('canopy_longwave_radiation', 'LWRadCanGPrev_col', 'MJ h-1', kRolePrivate, &
    'carry-over, not an ATS flux: canopy longwave emission of the previous step, total per grid cell', &
    LWRadCanGPrev_col)
  call add2('canopy_latent_heat', 'TLEX_col', 'MJ m-1', kRolePrivate, &
    'carry-over, not a heat flux: latent heat flux x boundary-layer resistance, summed over previous step', &
    TLEX_col)
  call add2('canopy_sensible_heat', 'TSHX_col', 'MJ m-1', kRolePrivate, &
    'carry-over, not a heat flux: sensible heat flux x boundary-layer resistance, summed over previous step', &
    TSHX_col)

  ! Microbe biomass state (only when microbe model is active)
  if (sizes%microbe_model > 0) then
    call add6('heter_biomass', 'mBiomeHeter_vr', 'g d-2', kRolePrivate, &
      'heterotrophic microbial biomass: element x guild x complex x layer', &
      mBiomeHeter_vr, NumPlantChemElms, nHeterBioms, jcplx, nlay)
    call add5('autor_biomass', 'mBiomeAutor_vr', 'g d-2', kRolePrivate, &
      'autotrophic microbial biomass: element x guild x layer', &
      mBiomeAutor_vr, NumPlantChemElms, nAutorBioms, nlay)
  endif

  ! EcoSIM-only outputs (visualized by ATS, not used by ATS physics). Water vapor
  ! terms use EcoSIM's air -> surface sign: <0 = water leaving to the air.
  call add1('transpiration', 'a_Transpiration', 'm3 d-2 h-1', kRoleOutput, &
    'transpiration summed over PFTs; <0 = water to air', a_Transpiration)
  call add3('transpiration_pft', 'Transpiration_pft', 'm3 d-2 h-1', kRoleOutput, &
    'transpiration of individual PFTs; <0 = water to air', Transpiration_pft, 1, npft)
  call add1('evaporation_canopy', 'a_EvapCan', 'm2 d-2 h-1', kRoleOutput, &
    'canopy evaporation summed over PFTs, this hour; <0 = water to air', a_EvapCan)
  call add1('evaporation_ground', 'a_EvapGrnd', 'unannotated', kRoleOutput, &
    'top soil layer evaporation only, no snow evaporation/sublimation; <0 = water to air (TEvapXAir2Toplay_col)', a_EvapGrnd)
  call add1('evaporation_litter', 'a_EvapLitr', 'unannotated', kRoleOutput, &
    'litter evaporation; <0 = water to air (TEvapXAir2LitR_col)', a_EvapLitr)
  call add1('evaporation_snow', 'a_EvapSnow', 'm3 d-2 h-1', kRoleOutput, &
    'evaporation from snow, this hour; <0 = water to air (EVAPW_col added over substeps)', a_EvapSnow)
  call add1('sublimation_snow', 'a_Sublim', 'm3 d-2 h-1', kRoleOutput, &
    'sublimation from snow, this hour; <0 = water to air (EVAPS_col added over substeps)', a_Sublim)
  call add2('canopy_surface_water', 'WatHeldOnCanopy_col', 'm3 d-2', kRoleOutput, &
    'water held on the canopy, column total, this hour', WatHeldOnCanopy_col)

  end subroutine DefineEntries

!------------------------------------------------------------------------------------------

  subroutine new_entry(ats_name, ecosim_name, units, role, description, ncomp)
  ! Append one entry (metadata only)
  implicit none
  character(len=*), intent(in) :: ats_name, ecosim_name, units, description
  integer, intent(in) :: role, ncomp

  if (nentries == size(entries)) call endrun(trim(mod_filename)//' at line', __LINE__)
  nentries = nentries + 1
  entries(nentries)%ats_name = ats_name
  entries(nentries)%ecosim_name = ecosim_name
  entries(nentries)%units = units
  entries(nentries)%description = description
  entries(nentries)%role = role
  entries(nentries)%ncomp = ncomp
  nullify(entries(nentries)%p)
  end subroutine new_entry

!------------------------------------------------------------------------------------------

  subroutine add3(ats_name, ecosim_name, units, role, description, a, first, ncomp)
  ! (component, NY, NX) array; entry covers components first..first+ncomp-1
  implicit none
  character(len=*), intent(in) :: ats_name, ecosim_name, units, description
  integer, intent(in) :: role, first, ncomp
  real(r8), allocatable, target, intent(inout) :: a(:,:,:)

  call new_entry(ats_name, ecosim_name, units, role, description, ncomp)
  if (bind_arrays) then
    if (.not. allocated(a)) call endrun(trim(mod_filename)//' at line', __LINE__)
    entries(nentries)%p => a(first:first+ncomp-1, 1:ncol_bind, 1)
  endif
  end subroutine add3

!------------------------------------------------------------------------------------------

  subroutine add2(ats_name, ecosim_name, units, role, description, a)
  ! (NY, NX) column array, one component
  implicit none
  character(len=*), intent(in) :: ats_name, ecosim_name, units, description
  integer, intent(in) :: role
  real(r8), allocatable, target, intent(inout) :: a(:,:)

  call new_entry(ats_name, ecosim_name, units, role, description, 1)
  if (bind_arrays) then
    if (.not. allocated(a)) call endrun(trim(mod_filename)//' at line', __LINE__)
    entries(nentries)%p(1:1, 1:ncol_bind) => a(1:ncol_bind, 1)
  endif
  end subroutine add2

!------------------------------------------------------------------------------------------

  subroutine add1(ats_name, ecosim_name, units, role, description, a)
  ! (ncol) array, one component
  implicit none
  character(len=*), intent(in) :: ats_name, ecosim_name, units, description
  integer, intent(in) :: role
  real(r8), allocatable, target, intent(inout) :: a(:)

  call new_entry(ats_name, ecosim_name, units, role, description, 1)
  if (bind_arrays) then
    if (.not. allocated(a)) call endrun(trim(mod_filename)//' at line', __LINE__)
    entries(nentries)%p(1:1, 1:ncol_bind) => a(1:ncol_bind)
  endif
  end subroutine add1

!------------------------------------------------------------------------------------------

  subroutine add5(ats_name, ecosim_name, units, role, description, a, n1, n2, n3)
  ! (d1, d2, d3, NY, NX) array; first 3 dims flattened into ncomp = n1*n2*n3
  implicit none
  character(len=*), intent(in) :: ats_name, ecosim_name, units, description
  integer, intent(in) :: role, n1, n2, n3
  real(r8), allocatable, target, intent(inout) :: a(:,:,:,:,:)
  integer :: ncomp

  ncomp = n1 * n2 * n3
  call new_entry(ats_name, ecosim_name, units, role, description, ncomp)
  if (bind_arrays) then
    if (.not. allocated(a)) call endrun(trim(mod_filename)//' at line', __LINE__)
    if (size(a,1)*size(a,2)*size(a,3) /= ncomp) then
      write(*,*) 'add5 dimension mismatch for ', trim(ats_name), ': expected ', &
        ncomp, ' got ', size(a,1)*size(a,2)*size(a,3)
      call endrun(trim(mod_filename)//' at line', __LINE__)
    endif
    entries(nentries)%p(1:ncomp, 1:ncol_bind) => a(:,:,:, 1:ncol_bind, 1)
  endif
  end subroutine add5

!------------------------------------------------------------------------------------------

  subroutine add6(ats_name, ecosim_name, units, role, description, a, n1, n2, n3, n4)
  ! (d1, d2, d3, d4, NY, NX) array; first 4 dims flattened into ncomp = n1*n2*n3*n4
  implicit none
  character(len=*), intent(in) :: ats_name, ecosim_name, units, description
  integer, intent(in) :: role, n1, n2, n3, n4
  real(r8), allocatable, target, intent(inout) :: a(:,:,:,:,:,:)
  integer :: ncomp

  ncomp = n1 * n2 * n3 * n4
  call new_entry(ats_name, ecosim_name, units, role, description, ncomp)
  if (bind_arrays) then
    if (.not. allocated(a)) call endrun(trim(mod_filename)//' at line', __LINE__)
    if (size(a,1)*size(a,2)*size(a,3)*size(a,4) /= ncomp) then
      write(*,*) 'add6 dimension mismatch for ', trim(ats_name), ': expected ', &
        ncomp, ' got ', size(a,1)*size(a,2)*size(a,3)*size(a,4)
      call endrun(trim(mod_filename)//' at line', __LINE__)
    endif
    entries(nentries)%p(1:ncomp, 1:ncol_bind) => a(:,:,:,:, 1:ncol_bind, 1)
  endif
  end subroutine add6

!------------------------------------------------------------------------------------------

  integer function NumStateEntries(sizes)
  implicit none
  type(EcoSizes), intent(in) :: sizes

  call DefineEntries(sizes, .false.)
  NumStateEntries = nentries
  end function NumStateEntries

!------------------------------------------------------------------------------------------

  subroutine GetStateEntry(i, sizes, ats_name, ecosim_name, units, description, ncomp, role)
  ! Layout of entry i (1-based)
  implicit none
  integer, intent(in) :: i
  type(EcoSizes), intent(in) :: sizes
  character(len=*), intent(out) :: ats_name, ecosim_name, units, description
  integer, intent(out) :: ncomp, role

  call DefineEntries(sizes, .false.)
  if (i < 1 .or. i > nentries) call endrun(trim(mod_filename)//' at line', __LINE__)
  ats_name = entries(i)%ats_name
  ecosim_name = entries(i)%ecosim_name
  units = entries(i)%units
  description = entries(i)%description
  ncomp = entries(i)%ncomp
  role = entries(i)%role
  end subroutine GetStateEntry

!------------------------------------------------------------------------------------------

  subroutine CheckLayout(istate, sizes)
  ! Stop if the container ATS allocated does not match this list
  implicit none
  type(EcoInternalState), intent(in) :: istate
  type(EcoSizes), intent(in) :: sizes

  if (istate%layout_version /= kLayoutVersion .or. istate%num_entries /= nentries .or. &
      istate%values_per_column /= sum(entries(1:nentries)%ncomp) .or. &
      istate%num_columns /= sizes%num_columns) then
    write(*,*) "EcoSIM internal state layout mismatch: version ", istate%layout_version, &
      " (expected ", kLayoutVersion, "), entries ", istate%num_entries, " (expected ", nentries, &
      "), values per column ", istate%values_per_column, " (expected ", &
      sum(entries(1:nentries)%ncomp), ")"
    call endrun(trim(mod_filename)//' at line', __LINE__)
  endif
  end subroutine CheckLayout

!------------------------------------------------------------------------------------------

  subroutine PackInternalState(istate, sizes, include_outputs)
  ! EcoSIM -> container. Outputs are skipped at setup, before they are computed.
  implicit none
  type(EcoInternalState), intent(in) :: istate
  type(EcoSizes), intent(in) :: sizes
  logical, intent(in) :: include_outputs

  real(r8), pointer :: v(:,:)
  integer :: i, off

  call DefineEntries(sizes, .true.)
  call CheckLayout(istate, sizes)
  call c_f_pointer(istate%values%data, v, [istate%values_per_column, istate%num_columns])
  off = 0
  do i = 1, nentries
    if (include_outputs .or. entries(i)%role == kRolePrivate) then
      v(off+1:off+entries(i)%ncomp, :) = entries(i)%p
    endif
    off = off + entries(i)%ncomp
  enddo
  end subroutine PackInternalState

!------------------------------------------------------------------------------------------

  subroutine UnpackInternalState(istate, sizes)
  ! container -> EcoSIM, carried state only
  implicit none
  type(EcoInternalState), intent(in) :: istate
  type(EcoSizes), intent(in) :: sizes

  real(r8), pointer :: v(:,:)
  integer :: i, off

  call DefineEntries(sizes, .true.)
  call CheckLayout(istate, sizes)
  call c_f_pointer(istate%values%data, v, [istate%values_per_column, istate%num_columns])
  off = 0
  do i = 1, nentries
    if (entries(i)%role == kRolePrivate) then
      entries(i)%p = v(off+1:off+entries(i)%ncomp, :)
    endif
    off = off + entries(i)%ncomp
  enddo

  call RecomputeSnowTotals(sizes%num_columns)
  end subroutine UnpackInternalState

!------------------------------------------------------------------------------------------

  subroutine RecomputeSnowTotals(num_cols)
  ! Column snow totals derived from the layers, as computed in SnowMassUpdate
  use GridConsts,   only : JS
  use EcosimConst,  only : DENSICE
  use SnowDataType, only : VLDrySnoWE_snvr, VLWatSnow_snvr, VLIceSnow_snvr, &
    SnowThickL_snvr, VLSnoDWIprev_snvr, VcumDrySnoWE_col, VcumWatSnow_col, &
    VcumIceSnow_col, VcumSnoDWI_col, VcumSnowWE_col, cumSnowDepz_col
  implicit none
  integer, intent(in) :: num_cols
  integer :: NY, L

  do NY = 1, num_cols
    VcumDrySnoWE_col(NY,1) = sum(VLDrySnoWE_snvr(1:JS,NY,1))
    VcumWatSnow_col(NY,1)  = sum(VLWatSnow_snvr(1:JS,NY,1))
    VcumIceSnow_col(NY,1)  = sum(VLIceSnow_snvr(1:JS,NY,1))
    VcumSnoDWI_col(NY,1)   = sum(VLSnoDWIprev_snvr(1:JS,NY,1))
    VcumSnowWE_col(NY,1)   = VcumDrySnoWE_col(NY,1)+VcumIceSnow_col(NY,1)*DENSICE+VcumWatSnow_col(NY,1)
    cumSnowDepz_col(0,NY,1) = 0.0_r8
    do L = 1, JS
      cumSnowDepz_col(L,NY,1) = cumSnowDepz_col(L-1,NY,1)+SnowThickL_snvr(L,NY,1)
    enddo
  enddo
  end subroutine RecomputeSnowTotals

!------------------------------------------------------------------------------------------
end module ATSStateRegistryMod
