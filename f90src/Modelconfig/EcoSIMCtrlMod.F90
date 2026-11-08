module EcoSIMCtrlMod
  ! Description:
  ! Store model controls, input selections and runtime configuration shared by drivers.
  use ncdio_pio,       only: file_desc_t
  use data_kind_mod,   only: r8 => DAT_KIND_R8
  use ecosim_Time_Mod, only: ecosim_time_type
  use fileUtil,        only: datestrlen
  use abortutils,      only: endrun, iulog
implicit none
  save
  character(len=*),private, parameter :: mod_filename =&
   __FILE__
  logical :: lmicrobeMLdiag = .false.
  logical :: LMacPoreOff          =.false.     !turnoff macropore, by default on
! --- Jing Tao ---
! numerical guard for a DYING stand. DEFAULT .false. so a switch-off build is
! bit-for-bit identical to the R3 baseline binary (exp/sposc-input-param @ 12a47976).
!
! WHY. Two independent guards in the plant BGC are scaled by a quantity that VANISHES
! as a stand dies, so both stop protecting exactly when they are needed:
!   (a) PlantBranchMod::WithdrawNutBranchLeaves divides by CPOOLT under the guard
!       CPOOLT > ZERO4Groth_pft, and ZERO4Groth_pft = ZERO*PlantPopuLive_pft
!       (InitPlantMod.F90:125) is proportional to the LIVE PLANT POPULATION. As the
!       population collapses the threshold collapses with it and the test degenerates
!       to "CPOOLT > 0", admitting an absolutely-tiny denominator -> SIGFPE.
!   (b) PlantBalMod normalises its carbon-balance relative error by total vegetation
!       carbon, which also collapses, so rounding-level residuals read as a large
!       relative error and abort the run.
! Measured on the EcoSIM_BioCON R3 Sobol prefix, 2026-08-20: 436 of 7071 attempted
! cases (7.5%) died this way -- 336 SIGFPE, 100 balance aborts -- and 12 of 12 balance
! aborts inspected had total vegetation C < 1 gC/m2 with GPP = 0.
!
! SCOPE, deliberately MINIMAL. This does NOT change ZERO4Groth_pft itself: that symbol
! has 260 references across 12+ files, 226 of them threshold comparisons, and it is
! used correctly elsewhere as an "is there any plant here" test. Only the DIVISION and
! the BALANCE TEST are protected, at their own sites.
!
! WHAT IT DOES NOT DO: it does not make these configurations viable. The stands still
! die; they die quietly and are counted as dead rather than as crashes. That matters
! for sensitivity analysis (a dead-but-complete case is a valid Y value, a crashed one
! is a hole in the design) but it is NOT a calibration improvement.
! Full record: A2MC memory/model_logs/20260820a_Two_EcoSIM_Guards_Scaled_By_The_Dying_Quantity.md
  logical :: guard_floor_dying_stand = .false.

! --- Jing Tao ---
! absolute floor [gC] on the carbon-pool DENOMINATOR in the leaf nutrient-withdrawal
! division. Far below any physically meaningful pool (a healthy stand carries
! O(1e2-1e3) gC/m2) and far above the denormal range, so it bounds the quotient
! without perturbing any live-plant arithmetic.
  real(r8), parameter :: CPOOL_DIV_FLOOR_dying = 1.0e-12_r8

! --- Jing Tao ---
! total vegetation carbon [gC/m2] below which the carbon-balance RELATIVE test is
! meaningless and is skipped (the absolute test still applies). Chosen as 1.0 because
! the 12 inspected aborts span 7e-4 to 0.54 gC/m2 while a healthy stand is 3e2-2e3,
! so 1.0 is ~0.05% of a live stand yet above every observed degenerate case.
! PI REVIEW WANTED on this value -- it decides which aborts are suppressed.
  real(r8), parameter :: VEGC_BAL_FLOOR_dying = 1.0_r8
  integer :: idebug_day           = -1         !specify the day being debugged
  integer :: idebug_year          = -1         !specify the year being debugged
  integer :: iselect_plantZ       = -1
  logical :: salt_model           = .false.    !toggle for salt model
  logical :: erosion_model        = .false.
  logical :: first_topou          = .false.  !only simulate first topo unit
  logical :: first_pft            = .false.  !only simulate first pft
  integer :: iErosionMode         = -1         !erosion option
  logical :: use_ecosys_drainage = .false.
  character(len=*), parameter :: DefaultSoilHydraulicScheme='ecosim'
  character(len=32) :: soil_hydraulic_scheme=DefaultSoilHydraulicScheme
  logical :: plant_model          = .true.
  logical :: microbial_model      = .true.
  logical :: soichem_model        = .true.
  logical :: snowRedist_model     = .true.
  logical :: llignification       = .false.
  logical :: lcoarseroot          =.false.
  logical :: ATS_cpl_mode         = .false.
  logical :: plantOM4Heat         = .false.
  logical :: fixWaterLevel        = .false.
  logical :: lsoilCompaction      = .false.
  logical :: ldo_sp_mode          = .false.   !do satellite phenology mode
  logical :: mod_snow_albedo      = .false.
  logical :: ldo_radiation_test   = .false.   !flag to turn on radiation test for prescribed phenolgoy
  logical :: ldo_transpt_bubbling = .true.   !flag to turn on/off ebullition during transport
  integer :: yearf1       !first year of daily climate forcing
  integer :: yearf2       !first year of hourly climate forcing
  integer :: nyeardal1    !number of daily climate forcing
  real(r8) :: aco2_ppm        = 280._r8
  real(r8) :: ach4_ppm        = 1.144_r8
  real(r8) :: an2o_ppm        = 0.270_r8
  real(r8) :: ao2_ppm         = 0.209e6_r8
  real(r8) :: arg_ppm         = 0.00934e6_r8
  real(r8) :: an2_ppm         = 0.78e6_r8
  real(r8) :: anh3_ppm        = 5.e-3_r8
  real(r8) :: ah2_ppm         = 0.55_r8
  real(r8) :: atm_co2_fix     = -100._r8
  real(r8) :: atm_ch4_fix     = -100._r8
  real(r8) :: atm_n2o_fix     = -100._r8
  character(len=256) :: warming_exp = ''
  character(len=300) :: pft_file_in
  character(len=300) :: pft_mgmt_in
  character(len=300) :: grid_file_in
  character(len=300) :: clm_hour_file_in =''    !file for hourly climate forcing
  character(len=300) :: clm_day_file_in  =''    !file for daily climate forcing
  character(len=300) :: soil_mgmt_in            !file for soil management information
  ! Opt in to dated soil schedules independent of recycled weather years.
  ! False preserves legacy forcing-year routing; no management is synthesized.
  logical :: soil_mgmt_use_model_year = .false.
  character(len=300) :: clm_factor_in           !file for climate change factors
  character(len=300) :: atm_ghg_in              !file for atmospheric GHG concentrations
  character(len=300) :: sp_phenol_in  =''    !file for monthly LAI, SAI and canopy height in prescribed phenology mode
  character(len=300) :: micpar_file_in=''       !file of microbial parameters 
  integer :: iVerbLevel = 0        
  logical :: do_budgets = .false.
  type(file_desc_t)  :: pft_nfid
  type(ecosim_time_type) :: etimer

  logical :: fixClime =.false.    !If true, do fixed climate forcing
  character(len=8) :: diag_opt ='nsteps'
  integer :: diag_frq=-999999999
  logical :: continue_run
  logical :: restart_out
  logical :: visual_out
  logical :: hist_yrclose
  character(len=16) :: hist_config(10)
  character(len=8)  :: sim_yyyymmdd
  integer :: forc_periods(15)
  integer :: NPXS(5),NPYS(5)
  integer :: NCYC_LITR  !number of subcycles for litr
  integer :: NCYC_SNOW  !number of subcycles for snow
  logical :: lverb           !logical switch for verbose output
  logical :: disp_planttrait =.true.
  logical :: disp_modelconfig=.true.
  logical :: do_rgres        !logical switch for regression tests
  integer :: grid_mode = 3  !vertical only
  
  type, public :: forc_data_rec_type
  integer :: pft_rec
  integer :: yearclm
  integer :: yearacc      !year passed
  integer :: yearcur      !current year
  integer :: yearpre      !previous year
  integer :: yearrst      !restart year
  logical :: lskip_loop   !logical switch to skip a loop
  character(len=datestrlen) :: ymdhs0  !the beginning yyyymmddhhmmss
  contains
  procedure, public :: Init  => init_frectyp
  end type forc_data_rec_type

  type(forc_data_rec_type), public :: frectyp

  contains

  subroutine Init_frectyp(this)
  implicit none
  class(forc_data_rec_type)  :: this

  this%pft_rec=0    !rec index for pft data
  this%yearclm=0    !year index for climate data
  this%yearacc=0    !accumulated years for the whole simulation
  this%yearcur=0    !current year
  this%yearpre=0    !previous year
  this%yearrst=0    !restarting year
  this%lskip_loop=.true.
  this%ymdhs0='00000000000000'
  end subroutine Init_frectyp

  !-----------------------------------------------------------------------
  integer function get_sim_len(forc_periods,nperiods)
  implicit none
  integer, dimension(:), intent(in) :: forc_periods
  integer :: nperiods
  integer :: nn1,id,nelms

  nelms = size(forc_periods)
  nperiods   = nelms/3

  get_sim_len=0
  DO nn1=0,nperiods-1
    id=nn1*3+1
    get_sim_len=get_sim_len+(abs(forc_periods(id+1)-forc_periods(id))+1)*forc_periods(id+2)
  enddo
  if(get_sim_len<0)then
  call endrun('Negative simulation length, check forc_periods set up in '//mod_filename,__LINE__)
  endif
  end function get_sim_len


end module EcoSIMCtrlMod
