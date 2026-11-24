PROGRAM main
!!
! Description:
! THIS SUBROUTINE READS THE RUNSCRIPT AND ENTERS FILENAMES INTO DATA ARRAYS
! FOR USE IN 'READS' AND 'READQ'. WHEN FINISHED THIS SUBROUTINE CALLS
! 'SOIL' WHICH IS THE MAIN SUBROUTINE FROM WHICH ALL OTHERS ARE CALLED
!
  use data_kind_mod     , only : r8 => DAT_KIND_R8
  use TestMod           , only : regression
  use InitEcoSIM        , only : InitModules
  use EcoSIMDesctruct   , only : DestructEcoSIM
  use GridMod           , only : SetMesh
  use readiMod          , only : readi
  use abortutils        , only : iulog
  USE fileUtil          , ONLY : ecosim_namelist_buffer_size,namelist_to_buffer
  use HistFileMod       , only : hist_htapes_build
  use EcoSIMConfig      , only : case_name,set_sim_type,start_date,is_restart,datestrlen
  use StartsMod         , only : set_ecosim_solver
  use RestartMod        , only : get_restart_date
  use ClimReadMod       , only : get_clm_years
  use EcosysWarmingMod   , only : config_soil_warming
  use EcoSIMAPI         , only : AdvanceModelOneYear,readnamelist,regressiontest,write_modelconfig, &
                              do_timing,do_timing_detail
  use timings, only : init_timer,start_timer,end_timer,finalize_timer,timer_stamp_type,flush_timer_detail, &
                              timer_launch_label
  use EcoSIMCtrlMod
  use EcoSIMCtrlDataType
  use EcoSIMHistMod
  use EcosimConst
  implicit none

  character(len=*), parameter :: mod_filename = &
  __FILE__

  type(timer_stamp_type) :: initialization_stamp
  character(len=15) :: run_start_label
  integer :: NA(250),ND(250)
  integer :: NAX,NDX,NEX,NAY,NDY,NE,N,NTX,NT
  integer :: NHW,NVN,NHE,NVS
  integer :: nn1,nn2,nn3
  integer :: year_ini,nyr1,yeari,nstopyr
  CHARACTER(len=640):: BUF
  character(len=256):: nmlfile
  character(len=14) :: ymdhs
  character(len=14) :: ymdhs0
  logical :: is_dos,nlend
  character(len=datestrlen) :: curr_date
  integer :: nmicbguilds
  integer :: idy, nperiods
  character(len=ecosim_namelist_buffer_size) :: nml_buffer

!!
! begin_execution

  run_start_label=timer_launch_label()
  is_dos=.false.
!  open(111,file='fort11',status='unknown')
  write(iulog,*)'obtain working directory'
  CALL GETCWD(BUF)
!
! IDENTIFY OPERATING SYSTEM: DOS OR UNIX

  IF((.NOT.(BUF(1:1).EQ.'/'.OR.BUF(1:1).EQ.'~')).AND.BUF(2:2).EQ.':')THEN
    write(iulog,*)'Running ecosim on dos system'
    is_dos=.true.
    CALL GETARG(2,BUF)
    CALL CHDIR(BUF)
    PREFIX='.\\'  ! location where input files are put
!   make output directory
  ELSE
    write(iulog,*)'Running ecosim on unix system'
!   make output directory
  ENDIF
!
  CALL GETARG(1,nmlfile)

  call namelist_to_buffer(nmlfile,nml_buffer)

  write(iulog,*)'read namelist'
  call readnamelist(nml_buffer, case_name, prefix, LYRG, nmicbguilds)

  if(is_dos)then
    outdir=trim(buf)//'\\'//trim(case_name)//'_outputs\\'
  else
    outdir=trim(buf)//'/'//trim(case_name)//'_outputs/'
  endif
  call system('mkdir -p '//trim(outdir))

  if(do_timing)then
    call init_timer(outdir,do_timing_detail,run_start_label)
    call start_timer(initialization_stamp)
  endif

  read(start_date,'(I4)')year_ini

! NUMBER OF COLUMNS AND ROWS for the whole land scape
!
  call SetMesh(NHW,NVN,NHE,NVS)

  call  InitModules()

  if(lverb)WRITE(*,*)'read initialization information READI'
  CALL readi(NHW,NHE,NVN,NVS)

  call write_modelconfig()

  call set_sim_type()

  nstopyr=get_sim_len(forc_periods,nperiods)

  call etimer%update_sim_len(nstopyr)
  
  call frectyp%Init()

  yeari=year_ini

  if(is_restart())then
    call get_restart_date(curr_date)
    frectyp%ymdhs0=curr_date    !the actual simulation beginning year
    print*,'read restart/checkpoint info file: ecosim_rst ',curr_date    
    read(curr_date,'(I4)')frectyp%yearrst
  else
    frectyp%ymdhs0=start_date   !the actual simulation beginning year
  endif
  
  call hist_htapes_build()
  
  !prepare climate forcing
  call get_clm_years()
  
  if(do_timing)then
    call end_timer('Initialization',initialization_stamp)
    call flush_timer_detail()
  endif

  IGO=0

!  print*,frectyp%ymdhs0,yeari

  DO nn1=1,nperiods
    call set_ecosim_solver(NPXS(NN1),NPYS(NN1),NCYC_LITR,NCYC_SNOW)

    do nn2=1,forc_periods(nn1*3)
      nn3=(nn1-1)*3
      !determine the step size
      idy=1      
      if(forc_periods(nn3+1)>forc_periods(nn3+2))idy=-1

      !do the simulation loop for the period
      do nyr1=forc_periods(nn3+1),forc_periods(nn3+2),idy
        frectyp%yearclm=nyr1
        frectyp%yearcur=etimer%get_curr_yearAD()
        nlend=.false.
        IGO=yeari-year_ini
        if(frectyp%yearcur==yeari)then
          call AdvanceModelOneYear(NHW,NHE,NVN,NVS,nlend)
        endif
        if(nlend)exit
        frectyp%yearacc=frectyp%yearacc+1
        call etimer%get_ymdhs(ymdhs)
        if(.not.frectyp%lskip_loop)print*,frectyp%yearcur,nyr1,ymdhs
        yeari=yeari+1
      end do
      if(nlend)exit
    end do
    if(nlend)exit
  end do

  if(do_rgres)then
    call regressiontest(trim(nmlfile),trim(case_name),NHW,NVN)
  endif
  call DestructEcoSIM
  if(do_timing)call finalize_timer()
!close(111)  
END program main
