module EcoSIMAPI
  use timings,          only: start_timer, end_timer, timer_stamp_type
  use data_kind_mod,    only: yearIJ_type
  use MicBGCAPI,        only: MicrobeModel,        MicAPI_Init,              MicAPI_cleanup
  use TracerIDMod,      only: ids_NO2B,            ids_NO2,                  idg_O2
  use EcosysWarmingMod, only: check_warming_year, apply_soil_cable_warming, config_soil_warming
  use ErosionMod,       only: erosion
  use Hour1Mod,         only: hour1
  use abortutils,       only: iulog
  use RedistMod,         only: redist
  use GeochemAPI,        only: soluteModel
  use PlantMod,          only: PlantModel, UpdateDailyGPP
  use TranspNoSaltMod,   only: TranspNoSalt
  use TranspSaltMod,     only: TranspSalt
  use WatsubMod,         only: watsub
  use PlantMgmtDataType, only: NP_col
  use FireMod,           only: config_fire
  use EcoSIMSolverPar,   only: oscal_test
  use DebugToolMod ,     only: PrintInfo,DebugPrint
  USE EcoSIMCtrlDataType
  use SoilWaterDataType
  use EcoSIMCtrlMod  
  use BalancesMod  
  use SoilDiagsMod 
  use MLDataDiagMod 
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: AdvanceModelOneYear
  public :: readnamelist
  public :: regressiontest,write_modelconfig
  public :: do_timing, do_timing_detail
  logical :: do_timing, do_timing_detail
contains

  subroutine Run_EcoSIM_one_step(yearIJ,NHW,NHE,NVN,NVS)

  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  integer, intent(in) :: NHW,NHE,NVN,NVS
  type(timer_stamp_type) :: t1
  character(len=*), parameter :: subname='Run_EcoSIM_one_step'

  call PrintInfo('beg '//subname)
  if(do_timing)call start_timer(t1)
  
  !do surface energy and water processes
  CALL HOUR1(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)

  if(do_timing)call end_timer('HOUR1',t1)
  !
  !   CALCULATE SOIL ENERGY BALANCE, WATER AND HEAT FLUXES IN 'WATSUB'
  !
  if(do_timing)call start_timer(t1)
  CALL WATSUB(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)
  if(do_timing)call end_timer('WAT',t1)
  !
  !   CALCULATE SOIL BIOLOGICAL TRANSFORMATIONS IN 'NITRO'
  !     
  if(microbial_model)then
    if(do_timing)call start_timer(t1)
    CALL MicrobeModel(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)
    if(do_timing)call end_timer('NIT',t1)
    
    if(lmicrobeMLdiag)call MLDataDiag(yearIJ,NHW,NHE,NVN,NVS)
  endif
  !
  !   UPDATE PLANT biogeochemistry
  !
  if(plant_model .and. (.not.ldo_radiation_test))then
    if(do_timing)call start_timer(t1)  
    call PlantModel(yearIJ,NHW,NHE,NVN,NVS)
    if(do_timing)call end_timer('PlantModel',t1)    
  endif

  call UpdateDailyGPP(yearIJ%J,NHW,NHE,NVN,NVS, &
    plant_model .and. .not.ldo_radiation_test .and. .not.ldo_sp_mode)

  !
  !   CALCULATE SOLUTE EQUILIBRIA IN 'SOLUTE'
  !
  if(soichem_model)then
    if(do_timing)call start_timer(t1)
    CALL soluteModel(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)
    if(do_timing)call end_timer('soluteModel',t1)
  endif

  !
  !   CALCULATE GAS AND SOLUTE FLUXES IN 'TranspNoSalt'
  !

  if(do_timing)call start_timer(t1)
  CALL TranspNoSalt(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)
  if(do_timing)call end_timer('TranspNoSalt',t1)

  !
  !   CALCULATE ADDITIONAL SOLUTE FLUXES IN 'TranspSalt' IF SALT OPTION SELECTED
  !

  if(salt_model)then
    if(do_timing)call start_timer(t1)
    CALL TranspSalt(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)
    if(do_timing)call end_timer('TranspSalt',t1)
  endif

  !
  !   CALCULATE SOIL SEDIMENT TRANSPORT IN 'EROSION'
  !
  if(lverb)WRITE(*,334)'EROSION'

  if(do_timing)call start_timer(t1)
  CALL EROSION(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)
  if(do_timing)call end_timer('EROSION',t1)
  !
  !   UPDATE ALL SOIL STATE VARIABLES FOR WATER, HEAT, GAS, SOLUTE
  !   AND SEDIMENT FLUXES IN 'REDIST'
  !
  if(do_timing)call start_timer(t1)    
  CALL REDIST(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)
  if(do_timing)call end_timer('REDIST',t1)
334   FORMAT(A8)

  if(do_timing)call start_timer(t1)
  call DiagSoilGasPressure(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)    

  call EndCheckBalances(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)

  if(do_timing)call end_timer('Diagnostics',t1)
  call PrintInfo('end '//subname)
  end subroutine Run_EcoSIM_one_step
! ----------------------------------------------------------------------
 subroutine readnamelist(nml_buffer, case_name, prefix,LYRG,nmicbguilds)
!!
! Description:
! read control namelist
  use abortutils     , only : endrun
  use EcoSIMConfig   , only : transport_on,column_mode, do_instequil,ref_date
  use EcoSIMConfig   , only : finidat,restartFileFullPath,brnch_retain_casename,start_date
  use ForcWriterMod  , only : bgc_forc_conf,do_bgcforc_write
  use fileUtil       , only : ecosim_namelist_buffer_size
  use EcoSIMHistMod  , only : DATAC
  use readsmod       , only : clim_var  
  use MicrobeConfigMod,only : ReadMicrobeNamelist
  use HistFileMod
  implicit none
  character(len=*), parameter :: mod_filename = &
  __FILE__
  character(len=ecosim_namelist_buffer_size), intent(in) :: nml_buffer
  character(len=36)    , intent(out) :: case_name
  character(len=80)    , intent(out) :: prefix
  integer              , intent(out) :: LYRG
  integer              , intent(out) :: nmicbguilds

  logical :: do_regression_test
  integer :: num_of_simdays
  logical :: lverbose
  integer :: num_microbial_guilds
  integer :: do_doy,do_year,do_layer
  character(len=64) :: bgc_fname
  character(len=256):: FireEvents
  real(r8) :: airT_C,Wind_ms,vap_Kpa,Rain_mmhr,SRAD_Wm2,Atm_kPa
  namelist /ecosim/case_name, prefix, do_regression_test, &
    num_of_simdays,lverbose,num_microbial_guilds,transport_on,column_mode,&
    do_instequil,salt_model, pft_file_in,grid_file_in,pft_mgmt_in, clm_factor_in,&
    clm_hour_file_in,clm_day_file_in,soil_mgmt_in,forc_periods,NCYC_LITR,NCYC_SNOW,&
    NPXS,NPYS,continue_run,restart_out,lsoilCompaction,LMacPoreOff,use_ecosys_drainage,&
    clm_hour_file_in,clm_day_file_in,soil_mgmt_in,soil_mgmt_use_model_year,forc_periods,NCYC_LITR,NCYC_SNOW,&
    NPXS,NPYS,continue_run,restart_out,lsoilCompaction,LMacPoreOff,&
    finidat,restartFileFullPath,brnch_retain_casename,plant_model,microbial_model,&
    soichem_model,atm_ghg_in,aco2_ppm,ao2_ppm,an2_ppm,ach4_ppm,anh3_ppm,&
    snowRedist_model,disp_planttrait,iErosionMode,grid_mode,atm_ch4_fix,atm_n2o_fix,&
    atm_co2_fix,first_topou,first_pft,fixWaterLevel,arg_ppm,idebug_day,idebug_year,ldo_sp_mode,sp_phenol_in,iverblevel,&
    ldo_radiation_test,ldo_transpt_bubbling,plantOM4Heat,micpar_file_in,iselect_plantZ,llignification
! --- Jing Tao ---
! guard_floor_dying_stand is appended to this group so ONE binary serves both arms of the
! verification, instead of needing a second build with the default flipped. It mirrors
! lmicrobeMLdiag below, which is likewise an EcoSIMCtrlMod logical surfaced through /ecosim/.
! DEFAULT-OFF IS PRESERVED BY THE READ ITSELF: a Fortran namelist read does not require every
! group member to be present, so the 14,849 R3 runfiles written before this variable existed
! simply leave it at its .false. initialiser. No existing case file needs rewriting, and a run
! against an old runfile is bit-for-bit unchanged. To enable it, add one line to the runfile:
!     guard_floor_dying_stand = .true.
  namelist /ecosim/hist_nhtfrq,hist_mfilt,hist_fincl1,hist_fincl2,hist_fincl3,hist_yrclose, lmicrobeMLdiag,&
    do_budgets,ref_date,start_date,do_timing,do_timing_detail,warming_exp,fixClime,FireEvents,oscal_test,&
    lcoarseroot,guard_floor_dying_stand

  logical :: laddband
  namelist /bbgcforc/do_bgcforc_write,do_year,do_doy,laddband,do_layer,&
    bgc_fname

  namelist /FixClimForc/airT_C,Wind_ms,vap_Kpa,Rain_mmhr,SRAD_Wm2,Atm_kPa
  !local variables
  character(len=256) :: ioerror_msg
  integer :: rc, fu
  integer :: nml_error
  integer :: year0

  do_timing    = .false.
  do_timing_detail = .false.
  continue_run = .false.
  NPXS         = 30   !number of cycles per hour for water, heat, solute flux calcns
  NPYS         = 10   !number of cycles per NPX for gas flux calculations

  lsoilCompaction=.false.
  idebug_day  =-1
  idebug_year =-1
  NCYC_LITR             = 20
  NCYC_SNOW             = 20
  grid_mode             = 3
  iErosionMode          = -1
  lmicrobeMLdiag        = .false.
  restart_out           = .false.
  do_budgets            = .false.
  plant_model           = .true.
  soichem_model         = .true.
  microbial_model       = .true.
  llignification        = .false.
  lcoarseroot           = .false.
  disp_planttrait       = .false.
  ref_date              = '18000101000000'   !place holder for future
  start_date            = '18000101000000'   !start date of the simulation, differ from the forcing date
  finidat               = ''
  restartFileFullPath   = ''                 !contains the full pathname of the restart file
  warming_exp           = ''
  brnch_retain_casename = .false.
  hist_yrclose          = .false.
  fixWaterLevel         = .false.
  FireEvents            =''
  forc_periods      = 0
  forc_periods(1:9) = (/1980,1980,1,1981,1988,2,1989,2008,1/)

  num_of_simdays       = -1
  do_year              = -1
  do_doy               = 0
  do_layer             = 1
  LMacPoreOff          = .false.
  use_ecosys_drainage  = .false.
  salt_model           = .false.
  laddband             = .false.
  do_regression_test   = .false.
  lverbose             = .false.
  num_microbial_guilds = 1
  do_bgcforc_write     = .false.
  bgc_fname            = 'bbforc.nc'
  do_instequil         = .false.
  plantOM4Heat         = .false.
  clm_factor_in      = 'NO'
  pft_file_in        = ''
  grid_file_in       = ''
  pft_mgmt_in        = 'NO'
  clm_hour_file_in   = ''
  clm_day_file_in    = ''
  soil_mgmt_in       = 'NO'
  soil_mgmt_use_model_year = .false. ! legacy default; opt in for dated schedules
  atm_ghg_in         = ''
  sp_phenol_in       = ''
  aco2_ppm           = 280._r8
  ach4_ppm           = 1.144_r8
  an2o_ppm           = 0.270_r8
  ao2_ppm            = 0.209e6_r8
  arg_ppm            = 0.00934e6_r8
  an2_ppm            = 0.78e6_r8
  anh3_ppm           = 5.e-3_r8
  atm_co2_fix        = -100._r8
  atm_n2o_fix        = -100._r8
  atm_ch4_fix        = -100._r8
  micpar_file_in     =''  
  first_topou        = .false.
  iselect_plantZ     = -1
  ldo_radiation_test = .false.
  
  read(nml_buffer, nml=ecosim, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     write(iulog,'(a)')"ERROR reading ecosim namelist ",nml_error,ioerror_msg
     call endrun('stopped in '//trim(mod_filename), __LINE__)
  end if

  read(nml_buffer, nml=bbgcforc, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     write(iulog,'(a)')"ERROR reading bbgcforc namelist "
     call endrun('stopped in '//trim(mod_filename), __LINE__)
  end if

  if (.true.) then
    write(iulog, *)
    write(iulog, *) '--------------------'
    write(iulog,ecosim)
    write(iulog, *)
    write(iulog, *) '--------------------'
    write(iulog,bbgcforc)
    write(iulog, *)
    write(iulog, *) '--------------------'
  endif
  erosion_model=iErosionMode>=0
  read(start_date,'(I4)')year0
  call etimer%Init(nml_buffer,year0=year0)
  if(do_bgcforc_write)then
    bgc_forc_conf%doy =do_doy
    bgc_forc_conf%year=do_year
    bgc_forc_conf%laddband=laddband
    bgc_forc_conf%Layer=do_layer
    bgc_forc_conf%bgc_fname=bgc_fname
  endif
  do_rgres=do_regression_test
  LYRG=num_of_simdays
  lverb=lverbose  
  nmicbguilds=num_microbial_guilds
  if(.not.soichem_model)then
    salt_model=.false.
  endif
  !radiation test can only be on for sp mode
  if(ldo_radiation_test)ldo_sp_mode=.true.
  if(ldo_sp_mode)then
    ! when using prescribed phenolgoy, turn on plant model for memory allocation, turn off soil chemistry model and microbial bgc
    plant_model           = .true.
    soichem_model         = .false.
    microbial_model       = .false.
  endif
  if(.not.plant_model)plantOM4Heat=.false.
  
  if (len_trim(warming_exp) > 10)call config_soil_warming(warming_exp)
  call config_fire(FireEvents)
  if(fixClime)then

    read(nml_buffer, nml=FixClimForc, iostat=nml_error, iomsg=ioerror_msg)
    if (nml_error /= 0) then
      write(iulog,'(a)')"ERROR reading FixClimForc namelist ",nml_error,ioerror_msg
      call endrun('stopped in '//trim(mod_filename), __LINE__)
    end if
    if(.true.)then
      write(iulog, *) '--------------------'
      write(iulog,FixClimForc)
      write(iulog, *)
      write(iulog, *) '--------------------'
    endif
    clim_var%airT_C    = airT_C
    clim_var%Wind_ms   = Wind_ms
    clim_var%vap_kpa   = vap_Kpa
    clim_var%Rain_mmhr = Rain_mmhr
    clim_var%SRAD_Wm2  = SRAD_Wm2
    clim_var%Atm_kPa   = Atm_kPa
  endif

  call ReadMicrobeNamelist(nml_buffer)

end subroutine readnamelist
! ----------------------------------------------------------------------

subroutine AdvanceModelOneYear(NHW,NHE,NVN,NVS,nlend)
!!
! Description:
! THIS IS THE MAIN SUBROUTINE FROM WHICH ALL OTHERS ARE CALLED
!
  use data_kind_mod,   only: r8 => DAT_KIND_R8
  use DayMod,          only: day, CompleteDailySetup
  use StarteMod,       only: starte
  use StartqMod,       only: startq
  use StartsMod,       only: starts,InitControlParms
  use WthrMod,         only: PrepHourlyWeather
  use RestartMod,      only: restFile
  use PlantInfoMod,    only: ReadPlantInfo
  use readsmod,        only: ReadClimSoilForcing
  use YearMod,         only: SetAnnualAccumlators
  use timings,         only: end_timer_loop, flush_timer_detail
  use ForcWriterMod,   only: do_bgcforc_write,   WriteBBGCForc
  use HistDataType,    only: hist_ecosim
  use HistFileMod,     only: hist_update_hbuf,   hist_htapes_wrapup
  use ClimReadMod,     only: read_soil_warming_Tref, ReadAtmosComposition
  use EcoSIMCtrlMod
  use GridConsts
  use EcoSIMCtrlDataType
  use EcoSIMHistMod
  use EcoSIMConfig

  implicit none
  integer :: yearc, yeari
  integer, intent(in) :: NHW,NHE,NVN,NVS
  logical, intent(out) :: nlend
  character(len=*), parameter :: mod_filename = &
  __FILE__
  type(timer_stamp_type) :: t1
  type(timer_stamp_type) :: step_stamp
  type(timer_stamp_type) :: setup_stamp
  type(yearIJ_type) :: yearIJ  

  character(len=*), parameter :: subname='AdvanceModelOneYear'
  integer :: I,J
  integer :: NX,NY
  integer :: NP_previous(JY,JX)
  integer :: idaz
  character(len=14) :: ymdhs
  logical :: rstwr, lnyr
  logical :: lverb0
  logical :: DeferDailySetup
  logical :: RefreshRainChem

! begin_execution
!
! READ INPUT DATA FOR SITE, SOILS AND MANAGEMENT IN 'READS'
! AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTS'
!
333   FORMAT(A8)
  call PrintInfo('beg '//subname)
  is_first_year=frectyp%yearacc.EQ.0

  if(do_timing)call start_timer(setup_stamp)

  !temporary set up for setting mass balance check
  IBEGIN=1;ISTART=1;ILAST=0

  call etimer%get_ymdhs(ymdhs)

  iYearCurrent=frectyp%yearcur

  !  Initialize controlling parameters
  if(do_timing)call start_timer(t1)
  call InitControlParms
  if(do_timing)call end_timer('InitControlParms',t1)

  ! Read annual atmospheric mixing ratios before STARTS converts them to mass concentrations.
  if(do_timing)call start_timer(t1)
  call ReadAtmosComposition(iYearCurrent,NHW,NHE,NVN,NVS)
  if(do_timing)call end_timer('ReadAtmosComposition',t1)

  IF(ymdhs(1:4)==frectyp%ymdhs0(1:4))THEN
    
    if(do_timing)call start_timer(t1)
    CALL STARTS(NHW,NHE,NVN,NVS)
    if(do_timing)call end_timer('STARTS',t1)
!
!   RECOVER VALUES OF ALL SOIL STATE VARIABLES FROM EARLIER RUN
!   IN 'ROUTS' IF NEEDED      
  ENDIF  
!
  if(plant_model)then
    !plant information is read in every year, but the active flags
    !are set using the checkpoint file.
    NP_previous=NP_col
    if(do_timing)call start_timer(t1)
    call ReadPlantInfo(frectyp%yearcur,frectyp%yearclm,NHW,NHE,NVN,NVS)
    if(do_timing)call end_timer('ReadPlantInfo',t1)
  endif

  if(do_timing)call start_timer(t1)
  CALL ReadClimSoilForcing(frectyp%yearcur,frectyp%yearclm,NHW,NHE,NVN,NVS,RefreshRainChem)
  if(do_timing)call end_timer('ReadClimSoilForcing',t1)

! INITIALIZE ALL PLANT VARIABLES IN 'STARTQ'
!
  if(ldo_sp_mode)then
    ! Prescribed PFT slots can appear in any forcing year. Initialize only
    ! newly named slots; existing canopies retain their water and root state.
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        call STARTQ(NX,NX,NY,NY,NP_previous(NY,NX)+1,NP_col(NY,NX))
      ENDDO
    ENDDO
  ELSEIF(ymdhs(1:4)==frectyp%ymdhs0(1:4) .and. plant_model)THEN
    !initialize by year
    if(do_timing)call start_timer(t1)
    CALL STARTQ(NHW,NHE,NVN,NVS,1,JP)
    if(do_timing)call end_timer('STARTQ',t1)
  ENDIF

  if(soichem_model)then
    ! Refresh irrigation chemistry after annual management reading; process
    ! rainfall only when its raw inputs were reloaded. Soil initialization
    ! remains restricted to the first cold-start year within STARTE.
    if(do_timing)call start_timer(t1)
    CALL STARTE(NHW,NHE,NVN,NVS,RefreshRainChem)
    if(do_timing)call end_timer('STARTE',t1)
  endif

  !load the reference whenever any day of the year is warmed, not only when Jan 1 is
  if(check_warming_year(iYearCurrent))then
    if(do_timing)call start_timer(t1)
    call read_soil_warming_Tref(iYearCurrent,NHW,NHE,NVN,NVS)    
    if(do_timing)call end_timer('ReadSoilWarmingTref',t1)
  endif
  lverb0      = lverb
  DazCurrYear = etimer%get_days_cur_year()
  yearIJ%year = frectyp%yearcur
  
  if(do_timing)call end_timer('AnnualSetup',setup_stamp)
  ! Keep setup out of the per-timestep detail records.
  if(do_timing)call flush_timer_detail()

  DO I  = 1, DazCurrYear
    call DebugPrint("beg step",I*1000)
    yearIJ%I=I
    if(idebug_day==I)then
      lverb=(idebug_year==iYearCurrent .and. idebug_year.GT.0) .or. (idebug_year.LT.0)      
    else
      lverb=lverb0
    endif  
    IF(do_rgres .and. I.eq.LYRG)RETURN
    !   UPDATE DAILY VARIABLES SUCH AS MANAGEMENT INPUTS
    !
    
    if(do_timing)call start_timer(t1)
    associate(RestartPending => frectyp%lskip_loop)
      DeferDailySetup=RestartPending .and. (is_restart() .or. is_branch())
    end associate
    CALL DAY(I,NHW,NHE,NVN,NVS,DeferStateChanges=DeferDailySetup)
    if(do_timing)call end_timer('DAY',t1)
    
    if(do_timing)call start_timer(t1)
    call SetAnnualAccumlators(I, NHW, NHE, NVN, NVS)
    if(do_timing)call end_timer('SetAnnualAccumlators',t1)

    DO J=1,24
      yearIJ%J=J
      call DebugPrint("beg step",I*1000+J)
      call etimer%get_ymdhs(ymdhs)
      
      if(ymdhs==frectyp%ymdhs0)then
        frectyp%lskip_loop=.false.        
        if(is_restart() .or. is_branch())then
          if(do_timing)call start_timer(t1)
          call restFile(flag='read')
          if(do_timing)call end_timer('RestartRead',t1)
          ! A January checkpoint precedes the new year's plant-management read.
          ! Restore last year's state first, then apply this year's slot counts,
          ! planting records and traits, as in an uninterrupted annual transition.
          if(I==1 .and. J==1)then
            call ReadAtmosComposition(iYearCurrent,NHW,NHE,NVN,NVS)
            if(plant_model)then
              NP_previous=NP_col
              if(do_timing)call start_timer(t1)
              call ReadPlantInfo(frectyp%yearcur,frectyp%yearclm,NHW,NHE,NVN,NVS)
              if(do_timing)call end_timer('ReadPlantInfo',t1)
              if(ldo_sp_mode)then
                DO NX=NHW,NHE
                  DO NY=NVN,NVS
                    call STARTQ(NX,NX,NY,NY,NP_previous(NY,NX)+1,NP_col(NY,NX))
                  ENDDO
                ENDDO
              endif
            endif
            ! Checkpoints also contain the previous year's soil-management
            ! state. Re-read the annual schedule against the restored geometry.
            call ReadClimSoilForcing(frectyp%yearcur,frectyp%yearclm,NHW,NHE,NVN,NVS,RefreshRainChem)
            if(soichem_model)call STARTE(NHW,NHE,NVN,NVS,RefreshRainChem)
          endif
          call CompleteDailySetup(I,J,NHW,NHE,NVN,NVS)
          if (j==1)then
            if(do_timing)call start_timer(t1)
            call SetAnnualAccumlators(I, NHW, NHE, NVN, NVS)
            if(do_timing)call end_timer('SetAnnualAccumlators',t1)
          endif

          if(do_timing)call start_timer(t1)
          call SummarizeTracerMass(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)          
          if(do_timing)call end_timer('SummarizeTracerMass',t1)
        endif
      endif
      if(frectyp%lskip_loop)then
        call etimer%update_time_stamp()
        cycle
      endif
    !
    !   UPDATE HOURLY VARIABLES IN 'HOUR1'
    !   set up climate forcing for the new hour

      if(do_timing)call start_timer(step_stamp)
      if(do_timing)call start_timer(t1)
      call PrepHourlyWeather(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)
      if(do_timing)call end_timer('WTHR',t1)

      if(do_timing)call start_timer(t1)
      call Run_EcoSIM_one_step(yearIJ,NHW,NHE,NVN,NVS)
      if(do_timing)call end_timer('EcoSIMOneStep',t1)

      if(do_timing)call start_timer(t1)
      call hist_ecosim%hist_update(yearIJ%I,yearIJ%J,bounds)
      if(do_timing)call end_timer('HistUpdate',t1)

      if(do_timing)call start_timer(t1)
      call hist_update_hbuf(bounds)
      if(do_timing)call end_timer('HistUpdateHbuf',t1)

      if(do_timing)call start_timer(t1)
      call etimer%update_time_stamp()

      nlend=etimer%its_time_to_exit()
      rstwr=etimer%its_time_to_write_restart(nlend)
      lnyr=etimer%its_a_new_year()
      if(do_timing)call end_timer('ClockUpdate',t1)

      if(do_timing)call start_timer(t1)
      call hist_htapes_wrapup( rstwr, nlend, bounds, lnyr )
      if(do_timing)call end_timer('HistoryIO',t1)
      
      if(rstwr)then
        
        if(do_timing)call start_timer(t1)
        call restFile(flag='write')
        if(do_timing)call end_timer('RestartWrite',t1)
      endif      
      if(do_timing)call end_timer('Timestep',step_stamp)
      if(do_timing)call end_timer_loop()
      call DebugPrint("end step",I*1000+J)
      if(nlend)exit
    END DO

!
! PERFORM MASS AND ENERGY BALANCE CHECKS IN 'EXEC'
!
    if(frectyp%lskip_loop)cycle

    if(do_bgcforc_write)then
      if(do_timing)call start_timer(t1)
      call WriteBBGCFORC(I,IYRR)
      if(do_timing)call end_timer('WriteBBGCForc',t1)
    endif
    if(nlend)exit
  END DO
  call PrintInfo('end '//subname)
  RETURN
END subroutine AdvanceModelOneYear
! ----------------------------------------------------------------------

subroutine regressiontest(nmfile,case_name, NX, NY)
!!
! Description:
! subroutine to conduct regression tests

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use TestMod      , only : regression
  use minimathmod  , only : safe_adb
  use GridConsts
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use PlantDataRateType
  use SoilBGCDataType
  use GridDataType
  implicit none

  character(len=*), parameter :: mod_filename = &
  __FILE__
  character(len=*), intent(in) :: nmfile
  character(len=*), intent(in) :: case_name
  integer, intent(in) :: NX
  integer, intent(in) :: NY
  !local variables
  character(len=128) :: category
  character(len=128) :: name
  integer :: nz, ll
  real(r8) :: datv(12)

  call regression%Init(trim(nmfile),case_name)

  if(regression%write_regression_output)then
    write(*,*)'write regression file'
    call regression%OpenOutput()

    do NZ=1,NP_col(NY,NX)
      IF(IsPlantActive_pft(NZ,NY,NX).EQ.iTrue)THEN

        category = 'flux'
        name = 'NH4_UPTK (g m^-3 h^-1)'
        datv=0._r8
        do ll=1,12
          if(AREA_3D(3,ll,NY,NX)>0._r8)then
            datv(ll)=safe_adb(RootNutUptake_pvr(ids_NH4,1,ll,NZ,NY,NX)+RootNutUptake_pvr(ids_NH4,2,ll,NZ,NY,NX) &
              +RootNutUptake_pvr(ids_NH4B,1,ll,NZ,NY,NX)+RootNutUptake_pvr(ids_NH4B,2,ll,NZ,NY,NX),AREA_3D(3,ll,NY,NX))
          endif
        enddo
        call regression%writedata(category,name, datv)
        exit
      endif

    enddo

    category = 'state'
    name = 'aqueous soil O2 (g m^3)'
    datv=trc_solcl_vr(idg_O2,1:12,NY,NX)
    call regression%writedata(category,name,datv)

    category = 'state'
    name = 'liquid soil water (m^3 m^-3)'
    datv=ThetaH2OZ_vr(1:12,NY,NX)
    call regression%writedata(category,name,datv)

    category = 'state'
    name = 'soil temperature (oC)'
    datv=TCS_vr(1:12,NY,NX)
    call regression%writedata(category,name,datv)
    write(*,*)'Finish regression file writing'
    call regression%CloseOutput()
  endif
end subroutine regressiontest

!------------------------------------------------------------------------------------------
  subroutine write_modelconfig()
  use fileUtil, only : getavu, relavu, opnfil
  use readiMod
  implicit none
  character(len=64) :: fnm_loc
  integer :: nu_cfg

  if(disp_modelconfig)then
    nu_cfg=getavu()
    write(fnm_loc,'(A,I4,A)')'ecosim.setup'
    call opnfil(fnm_loc,nu_cfg,'f')    

    write(nu_cfg,*)'microbial_model=',microbial_model
    write(nu_cfg,*)'plant_model=',plant_model
    write(nu_cfg,*)'salt_model=',salt_model
    write(nu_cfg,*)'erosion_model=',erosion_model_status(iErosionMode)
    write(nu_cfg,*)'grid_mode=',GridConectionMode(grid_mode)
    write(nu_cfg,*)'soil_mgmt_use_model_year=',soil_mgmt_use_model_year
    call relavu(nu_cfg)
  endif
  end subroutine write_modelconfig

end module EcoSIMAPI
