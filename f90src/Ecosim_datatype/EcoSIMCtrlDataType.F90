module EcoSIMCtrlDataType

! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  logical :: ats_shared_surface_layer=.false. ! ATS owns combined pond/litter inventories

  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  public
  real(r8) :: ZERO
  real(r8) :: ZERO2
  real(r8),target,allocatable :: ZEROS(:,:)
  real(r8),target,allocatable :: ZERO4Groth_pft(:,:,:)
  real(r8),target,allocatable :: ZERO4Uptk_pft(:,:,:)
  real(r8),target,allocatable :: ZERO4LeafVar_pft(:,:,:)
  real(r8),target,allocatable :: ZEROS2(:,:)

  integer :: NPX             !number of cycles per hour for water, heat, and solute flux calcns
  integer :: NPY             !number of cycles per NPX for gas flux calculations
  integer :: IGO             !flag for first scenario
  integer :: IDAYR           !day of recovery from earlier run
  integer :: IYRR            !year of recovery from earlier run
  integer :: ITERM           !end date for reading weather data
  integer :: IFIN            !end date for reading weather data
  integer :: IOLD            !last day of previous scenario
  integer :: ILAST           !last day of previous scenario
  integer :: IRUN            !start date of current scenario
  integer :: IBEGIN          !start date of model run
  integer :: ISTART          !start date of model run
  integer :: IEND            !end date of  current scenario
  integer :: LYRX            !last year
  integer :: DazCurrYear            !number of days in current year
  integer :: LYRG            !num_of_simdays, defined for regression tests

  integer :: iYearCurrent     !current year
  integer :: iyear_pre        !previous year
  integer :: iyear_rest       !restart year
  integer :: current_month, day_of_month, total_days_in_month
  contains

  subroutine InitEcoSIMCtrlData

  implicit none

  allocate(ZEROS(JY,JX)); ZEROS(:,:) = 0._r8
  allocate(ZERO4Groth_pft(JP,JY,JX)); ZERO4Groth_pft(:,:,:)=0._r8
  allocate(ZERO4Uptk_pft(JP,JY,JX)); ZERO4Uptk_pft(:,:,:)=0._r8
  allocate(ZERO4LeafVar_pft(JP,JY,JX)); ZERO4LeafVar_pft(:,:,:)=0._r8
  allocate(ZEROS2(JY,JX));  ZEROS2(:,:)=0._r8

  end subroutine InitEcoSIMCtrlData

  subroutine DestructEcoSIMCtrlData
  use abortutils, only : destroy

  implicit none

  call destroy(ZEROS)
  call destroy(ZERO4Groth_pft)
  call destroy(ZERO4Uptk_pft)
  call destroy(ZERO4LeafVar_pft)
  call destroy(ZEROS2)

  end subroutine DestructEcoSIMCtrlData

end module EcoSIMCtrlDataType
