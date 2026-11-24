module ecosim_Time_Mod
!
! DESCRIPTION
! the module contains subroutine to march the time

  use data_kind_mod  , only : r8 => DAT_kind_r8
  use data_kind_mod  , only : i8 => DAT_kind_i8
  use fileUtil       , only : stdout, ecosim_string_length_long
  use EcosimConst    , only : secspday,secspyear
  implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
  __FILE__


  real(r8), parameter :: uninit_r8  = -999999999.0
  integer , parameter :: daz(12)=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
  integer , parameter :: cdaz(12)=(/31,59,90,120,151,181,212,243,273,304,334,365/)
  integer , parameter :: cal_str_len=18   !YYYY-MM-DD-HHMMSS
  integer , parameter :: optstrlen=8
  type, public :: ecosim_time_dat_type
    integer :: year0   !year0 of the simulation
    integer :: nstep   !number of steps has passed
    integer :: tstep   !number of steps in the current year
  end type ecosim_time_dat_type

  type, public:: ecosim_time_type
     ! NOTE(bja, 201603) all real variables have units of seconds!
     real(r8) :: delta_time
     real(r8) :: curr_time       !time in seconds
     real(r8) :: time0
     real(r8) :: cur_timef
     real(r8) :: toy                 !time of year
     integer  :: curr_stop_count
     integer  :: stop_count
     integer  :: tstep               !number of time steps
     integer  :: nelapstep
     integer  :: dow, dom, doy       !day of week, day month, day of year
     integer  :: tod_prev,dom_prev,moy_prev,year_prev
     integer  :: moy                 !mon of year
     integer  :: cyears              !cumulative years
     integer  :: cdays               !cumulative days
     integer  :: year0   !beginning year AD
     real(r8) :: tod     !time of day in seconds
     integer  :: hist_freq   !negative number, steps, positive number, 1: day, 30:mon, 365:year
     integer  :: stop_opt
     integer  :: leap_yr
     integer  :: rest_n
     integer  :: rest_frq
     integer  :: diag_n     
     character(len=optstrlen)  :: diag_opt
     integer  :: diag_frq
     character(len=optstrlen) :: rest_opt
   contains
     procedure, public :: Init
     procedure, public :: its_time_to_write_restart
     procedure, public :: its_time_to_diag
     procedure, public :: its_time_to_exit
     procedure, public :: update_time_stamp
     procedure, public :: set_nstep
     procedure, public :: get_nstep
     procedure, public :: set_time_offset
     procedure, public :: get_prev_date
     procedure, public :: get_curr_date
     procedure, public :: get_curr_doy
     procedure, public :: get_curr_dom
     procedure, public :: get_days_per_year
     procedure, public :: get_step_size
     procedure, public :: get_prev_time
     procedure, public :: get_curr_time
     procedure, public :: get_curr_mon_days
     procedure, public :: getdatetime
     procedure, public :: get_curr_timeful  !current time + offset
     procedure, public :: get_days_cur_year
     procedure, public :: proc_initstep
     procedure, public :: print_curr_time
     procedure, public :: its_time_to_histflush
     procedure, public :: get_calendar
     procedure, public :: setClock
     procedure, public :: its_a_new_hour
     procedure, public :: its_a_new_day
     procedure, public :: its_a_new_week
     procedure, public :: its_a_new_month
     procedure, public :: its_a_new_year
     procedure, public :: get_ymdhs
     procedure, public :: get_curr_year
     procedure, public :: get_curr_yearAD
     procedure, public :: get_curr_day
     procedure, public :: get_curr_mon
     procedure, public :: is_first_step
     procedure, public :: config_restart
     procedure, public :: config_diag
     procedure, public :: print_model_time_stamp
     procedure, public :: update_sim_len
     procedure, private:: proc_nextstep
     procedure, private:: ReadNamelist
  end type ecosim_time_type
  public :: getdow
  public :: get_steps_from_ymdhs
contains

  subroutine setClock(this, dtime, nelapstep)

  implicit none
  class(ecosim_time_type), intent(inout) :: this
  real(r8), intent(in) :: dtime
  integer, intent(in) :: nelapstep

  this%delta_time = dtime
  this%nelapstep = nelapstep

  if(this%year0>0)then
    if(isLeap(this%year0))then
      this%leap_yr=1
    else
      this%leap_yr=0
    endif
  endif
  end subroutine setClock
  !-------------------------------------------------------------------------------
  function get_curr_time(this)result(ans)

  implicit none
  class(ecosim_time_type), intent(in) :: this
  real(r8) :: ans

  ans = this%curr_time
  end function get_curr_time
  !-------------------------------------------------------------------------------
  function get_curr_timeful(this)result(ans)

  implicit none
  class(ecosim_time_type), intent(in) :: this
  real(r8) :: ans

  ans = this%curr_time+this%time0
  end function get_curr_timeful
  !-------------------------------------------------------------------------------
  subroutine Init(this, namelist_buffer, masterproc, year0, nyears)

    implicit none

    class(ecosim_time_type), intent(inout) :: this
    character(len=*), optional, intent(in) :: namelist_buffer
    logical, optional, intent(in) :: masterproc
    integer, optional, intent(in) :: year0    !beginning year, used to record in the AD format
    integer, optional, intent(in) :: nyears   !run the timer for nyears
    integer :: N
    this%rest_n=0
    this%diag_n=0
    this%tstep = 0
    this%time0 = 0._r8
    this%curr_time  = 0._r8
    this%tod   = 0._r8
    this%toy   = 0._r8
    this%cyears = 0
    this%year0 = 0
    this%cdays  = 0
    this%dow    = 1
    this%dom    = 1
    this%doy    = 1
    this%moy    = 1
    this%hist_freq=-1
    this%nelapstep=0
    this%stop_opt=3    !stop by year
    this%leap_yr=0
    this%curr_stop_count=0
    this%stop_count=0
    if(present(namelist_buffer))then
      if(present(masterproc))then
        call this%ReadNamelist(namelist_buffer, masterproc)
      else
        call this%ReadNamelist(namelist_buffer)
      endif
    else
      if(present(nyears))then
        this%stop_count=nyears
      endif
    endif
    if(present(year0))this%year0=year0
    this%leap_yr=isLeapi(this%year0)
  end subroutine Init

  ! ----------------------------------------------------------------------
  subroutine config_diag(this,diag_frq,diag_opt)
  implicit none
  ! !ARGUMENTS:
  class(ecosim_time_type), intent(inout) :: this
  character(len=*), intent(in) :: diag_opt
  integer, intent(in) :: diag_frq

  this%diag_frq=diag_frq
  this%diag_opt=diag_opt
  end subroutine config_diag
  ! ----------------------------------------------------------------------
  subroutine config_restart(this,rest_frq,rest_opt)
  implicit none
  ! !ARGUMENTS:
  class(ecosim_time_type), intent(inout) :: this
  character(len=*), intent(in) :: rest_opt
  integer, intent(in) :: rest_frq

  this%rest_frq=rest_frq
  this%rest_opt=rest_opt

  end subroutine config_restart
  ! ----------------------------------------------------------------------

  subroutine ReadNamelist(this, namelist_buffer, masterproc)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use abortutils    , only : endrun
    use ecosim_log_mod   , only : errMsg => shr_log_errMsg
    use fileUtil , only : stdout, ecosim_string_length_long

    implicit none
    ! !ARGUMENTS:
    class(ecosim_time_type), intent(inout) :: this
    character(len=*), intent(in) :: namelist_buffer
    logical, optional, intent(in) :: masterproc
    !
    ! !LOCAL VARIABLES:
    integer :: nml_error
    character(len=*), parameter :: subname = trim(mod_filename)//'::ReadNamelist'
    character(len=ecosim_string_length_long) :: ioerror_msg
    real(r8) :: delta_time         !model time step
    integer  :: stop_n
    integer  :: year0
    integer  :: rest_frq
    integer :: diag_frq
    character(len=optstrlen)  :: rest_opt
    character(len=optstrlen)  :: diag_opt
    character(len=optstrlen)  :: stop_option        !1: step, 2:day, 3:year
    logical :: masterproc_loc
    integer :: jj
    !-----------------------------------------------------------------------

    namelist / ecosim_time / delta_time, stop_n, stop_option,  &
      rest_frq, rest_opt, diag_frq, diag_opt

    ! FIXME(bja, 201603) Only reading time variables in seconds!
    ! Should enable other values with unit coversions.

    ! FIXME(bja, 201603) assign some defaults, should eventually remove
    ! when all input files are updated.
    masterproc_loc=.true.
    if(present(masterproc))masterproc_loc=masterproc
    delta_time = 1800._r8                !half hourly time step
    stop_n=2                             !by default 2 cycle
    stop_option='nyears'                 !by default years
    this%stop_opt=3
    rest_opt='never'
    this%hist_freq=-1                    !not used
    year0=0
    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( index(trim(namelist_buffer),'ecosim_time')/=0 )then
       ioerror_msg=''
       read(namelist_buffer, nml=ecosim_time, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          call endrun(msg="ERROR reading secosim_driver namelist "//errmsg(mod_filename, __LINE__))
       end if
    end if

    if (masterproc_loc) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' ecosim time :'
       write(stdout, *)
       write(stdout, *) ' ecosim_time namelist settings :'
       write(stdout, *)
       write(stdout, ecosim_time)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif
    this%diag_frq=diag_frq
    this%diag_opt=diag_opt
    this%delta_time = delta_time
    this%year0=year0
    this%rest_opt=rest_opt
    this%rest_frq=rest_frq
    this%stop_count=stop_n
    this%curr_stop_count=0
    select case (trim(stop_option))
    case ('nsteps')
      this%stop_opt=1
    case ('ndays')
      !day
      this%stop_opt=2
    case ('nmonths')
      this%stop_opt=3
    case ('nyears')
      !year
      this%stop_opt=4
    case default
      call endrun(msg="ERROR setting up stop_option "//errmsg(mod_filename, __LINE__))
    end select

  end subroutine ReadNamelist

  !-------------------------------------------------------------------------------
    function its_time_to_write_restart(this,nlend)result(ans)
     !
     ! DESCRIPTION
     ! decide if to write restart file
     !

     implicit none
     class(ecosim_time_type), intent(inout) :: this
     logical, optional, intent(in) :: nlend
     logical :: ans

     character(len=80) :: subname = trim(mod_filename)//'::its_time_to_write_restart'

     ans=.false.
     select case (this%rest_opt)
     case ('nsteps')
       this%rest_n=mod(this%rest_n+1,this%rest_frq)
       ans=(this%rest_n==0)
     case ('nhours')
       if(this%its_a_new_hour())then
         this%rest_n=mod(this%rest_n+1,this%rest_frq)
         ans=(this%rest_n==0)
       endif
     case ('ndays')
       if(this%its_a_new_day())then
         this%rest_n=mod(this%rest_n+1,this%rest_frq)
         ans=(this%rest_n==0)
       endif
     case ('nweeks')
       if(this%its_a_new_week())then
         this%rest_n=mod(this%rest_n+1,this%rest_frq)
         ans=(this%rest_n==0)
       endif
     case ('nmonths')
       if(this%its_a_new_month())then
         this%rest_n=mod(this%rest_n+1,this%rest_frq)
         ans=(this%rest_n==0)
       endif
     case ('nyears')
       if(this%its_a_new_year())then
         this%rest_n=mod(this%rest_n+1,this%rest_frq)
         ans=(this%rest_n==0)
       endif
     case default
       ans=.false.
     end select
     if(present(nlend))ans=ans .or. nlend
     end function its_time_to_write_restart
  !-------------------------------------------------------------------------------
    function its_time_to_diag(this)result(ans)
     !
     ! DESCRIPTION
     ! decide if to write restart file
     !

     implicit none
     class(ecosim_time_type), intent(inout) :: this
     logical :: ans

     character(len=80) :: subname = trim(mod_filename)//'::its_time_to_diag'

     ans=.false.
     select case (this%diag_opt)
     case ('nsteps')
       this%diag_n=this%diag_n+1
       this%diag_n=mod(this%diag_n,this%diag_frq)
       ans=(this%diag_n==0)
     case ('nhours')
       if(this%its_a_new_hour())then
         this%diag_n=mod(this%diag_n+1,this%diag_frq)
         ans=(this%diag_n==0)
       endif
     case ('ndays')
       if(this%its_a_new_day())then
         this%diag_n=mod(this%diag_n+1,this%diag_frq)
         ans=(this%diag_n==0)
       endif
     case ('nweeks')
       if(this%its_a_new_week())then
         this%diag_n=mod(this%diag_n+1,this%diag_frq)
         ans=(this%diag_n==0)
       endif
     case ('nmonths')
       if(this%its_a_new_month())then
         this%diag_n=mod(this%diag_n+1,this%diag_frq)
         ans=(this%diag_n==0)
       endif
     case ('nyears')
       if(this%its_a_new_year())then
         this%diag_n=mod(this%diag_n+1,this%diag_frq)
         ans=(this%diag_n==0)
       endif
     end select

     end function its_time_to_diag
  !-------------------------------------------------------------------------------
  function its_time_to_exit(this) result(ans)
    !
    ! DESCRIPTION
    ! decide if to exit the loop
    !

    implicit none
    class(ecosim_time_type), intent(in) :: this
    logical :: ans

    character(len=80) :: subname = 'its_time_to_exit'

    ans = (this%curr_stop_count >= this%stop_count)

  end function its_time_to_exit

  !-------------------------------------------------------------------------------
  subroutine update_time_stamp(this)
    !
    ! DESCRIPTION
    !

    implicit none

    class(ecosim_time_type), intent(inout) :: this

    character(len=80) :: subname='update_time_stamp'

    this%curr_time = this%curr_time + this%delta_time
    this%toy  = this%toy + this%delta_time

    this%tstep = this%tstep + 1
    !
    ! reset the clock every year, and assuming the time step
    ! size is always
    if(mod(this%toy, secspyear+86400._r8*this%leap_yr) == 0) then
       this%tstep = 0
    end if

    this%tod_prev=this%tod
    this%dom_prev=this%dom
    this%moy_prev=this%moy
    this%year_prev=this%cyears

    !update time of the day
    this%tod=this%tod+this%delta_time

    !update varaibles when it is to start a new day
    if(this%its_a_new_day())then
      !it reaches 24:00, reset to 00:00
      this%tod=0._r8
      this%dom=this%dom+1
      this%dow = mod(this%dow + 1, 7)
      if(this%dow==0)this%dow=7
      this%doy = this%doy + 1
      this%cdays= this%cdays + 1
      if(this%stop_opt==2)then
        this%curr_stop_count=this%curr_stop_count+1
      endif
    endif

    if(this%its_a_new_month())then
      this%moy=this%moy+1
      if(this%moy==13)then
        this%moy=1
        this%doy=mod(this%doy,365+this%leap_yr)
        if(this%doy==0)this%doy=1
      endif
      this%dom=1
      if(this%stop_opt==3)then
        this%curr_stop_count=this%curr_stop_count+1
      endif
    endif

    if(this%its_a_new_year())then
      this%cyears=this%cyears+1
      if(this%year0>0 .and. isLeap(this%cyears+this%year0))then
        this%leap_yr=1
      else
        this%leap_yr=0
      endif
      if(this%stop_opt==4)then
        this%curr_stop_count=this%curr_stop_count+1
      endif
    endif
    if(this%stop_opt==1)then
      this%curr_stop_count=this%curr_stop_count+1
    endif

    call this%proc_nextstep()
  end subroutine update_time_stamp

  !-------------------------------------------------------------------------------
  real(r8) function get_step_size(this) result(rc)

    ! Return the step size in seconds.

    implicit none

    class(ecosim_time_type), intent(in) :: this

    rc = this%delta_time
    return

  end function get_step_size

  !-------------------------------------------------------------------------------
  integer function get_nstep(this)

    ! Return the timestep number.
    implicit none

    class(ecosim_time_type), intent(in) :: this

    character(len=*), parameter :: sub = 'betr::get_nstep'

    get_nstep = this%nelapstep

  end function get_nstep


  !-------------------------------------------------------------------------------
  subroutine set_nstep(this, nstep)

    ! Return the timestep number.
    implicit none
    class(ecosim_time_type), intent(inout) :: this

    character(len=*), parameter :: sub = 'betr::get_nstep'
    integer, intent(in) :: nstep

    this%nelapstep = nstep

    if(this%its_a_new_year())then
      this%tstep = 0
    endif

  end subroutine set_nstep


  !-------------------------------------------------------------------------------
  subroutine set_time_offset(this, nstep, tstep, continue_run)

    ! Return the timestep number.
    implicit none
    class(ecosim_time_type), intent(inout) :: this

    character(len=*), parameter :: sub = 'betr::get_nstep'
    integer, intent(in) :: nstep
    integer, intent(in) :: tstep
    logical, intent(in) :: continue_run
    real(r8) :: cdtime
    integer :: nn,cdays

    this%time0  = nstep*this%delta_time
    if(continue_run)then
      this%nelapstep = nstep
    else
      this%nelapstep = 0
    endif
    this%tstep=tstep
    if(this%its_a_new_year())then
      this%tstep = 0
    endif
    cdtime=this%nelapstep*this%delta_time
    this%cdays=int(cdtime/secspday)

    if(cdtime>0._r8)then
      if(this%year0>0)then
        if(this%cdays>365)then
          nn=0
          do while(.true.)
            if(isLeap(this%year0+nn))then
              if(cdays>366)then
                nn=nn+1
                cdays=cdays-366
              else
                exit
              endif
            else
              nn=nn+1
              cdays=cdays-365
              if(cdays<365)exit
            endif
          enddo
          this%cyears=nn
          this%doy=cdays
          this%dow=getdow(this%year0+this%cyears,this%doy)
          do nn=1,12
            if(this%doy<=cdaz(nn))then
              this%moy=nn
              exit
            endif
          enddo
        endif
      else
        this%cyears=int(cdtime/secspyear)

        this%doy=mod(this%cdays,365)
        this%dow=mod(this%cdays,7)
        do nn=1,12
          if(this%doy<=cdaz(nn))then
            this%moy=nn
            exit
          endif
        enddo
        this%tod=mod(cdtime,secspday)
        this%toy=mod(cdtime,secspyear)
      endif
    endif
  end subroutine set_time_offset

  !-------------------------------------------------------------------------------
  subroutine proc_nextstep(this)

    implicit none

    class(ecosim_time_type), intent(inout) :: this

    this%nelapstep = this%nelapstep + 1
  end subroutine proc_nextstep

  !-------------------------------------------------------------------------------
  subroutine proc_initstep(this)

    implicit none
    class(ecosim_time_type), intent(inout) :: this

    this%nelapstep = 0
  end subroutine proc_initstep

  !-------------------------------------------------------------------------------
  integer function get_days_per_year(this)

    implicit none

    class(ecosim_time_type), intent(in) :: this

    ! Positive for future times, negative
    ! for previous times.

    ! remove unused dummy arg compiler warning
    if (this%tstep > 0) continue

    !hardwire at the moment
    get_days_per_year = 365
  end function get_days_per_year

  !-------------------------------------------------------------------------------
  integer function get_days_cur_year(this)
    implicit none
    class(ecosim_time_type), intent(in) :: this
    integer :: year
    year=this%get_curr_yearAD()    
    get_days_cur_year=365
    if(isLeap(year))get_days_cur_year=366
  end function  get_days_cur_year
  !-------------------------------------------------------------------------------
  integer function get_curr_mon_days(this)
  implicit none
  class(ecosim_time_type), intent(in) :: this

  get_curr_mon_days=daz(this%moy)
  if(this%moy==2)get_curr_mon_days=get_curr_mon_days+this%leap_yr
      
  end function get_curr_mon_days
  !-------------------------------------------------------------------------------

  subroutine print_curr_time(this)
  implicit none
  class(ecosim_time_type), intent(in) :: this

  print*,'time=',this%curr_time
  print*,'tod =',this%tod
  print*,'dow =',this%dow
  print*,'dom =',this%dom
  print*,'doy =',this%doy
  print*,'moy =',this%moy
  print*,'cyears=',this%cyears

  end subroutine print_curr_time


  !-------------------------------------------------------------------------------
  subroutine get_ymdhs(this, ymdhs)

  implicit none
  class(ecosim_time_type), intent(in) :: this
  character(len=*), intent(out) :: ymdhs

  write(ymdhs,'(I4.4,I2.2,I2.2,I2.2,I4.4)')this%cyears+this%year0,this%moy,this%dom,&
     int(this%tod/3600.0),int(mod(this%tod,3600._r8))
  end subroutine get_ymdhs
  !-------------------------------------------------------------------------------
  function its_a_new_hour(this)result(yesno)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: yesno


  yesno= abs(mod(this%tod, 3600._r8))<1.e-3_r8

  end function its_a_new_hour

  !-------------------------------------------------------------------------------

  function its_a_new_day(this)result(yesno)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: yesno


  yesno= abs(mod(this%tod, 86400._r8))<1.e-3_r8

  end function its_a_new_day

  !-------------------------------------------------------------------------------

  function its_a_new_week(this)result(yesno)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: yesno


  yesno= (this%dow == 0 .and. this%tod<1.e-3_r8)

  end function its_a_new_week

  !-------------------------------------------------------------------------------

  function its_a_new_month(this)result(yesno)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: yesno

  if(this%moy/=2)then
    !not Feb
    yesno = ((this%dom == daz(this%moy)+1 .or. this%dom==1) .and. this%tod < 1.e-3_r8)
  else
    !Feb
    yesno = ((this%dom == daz(this%moy)+this%leap_yr+1 .or. this%dom==1) .and. this%tod < 1.e-3_r8)
  endif
  end function its_a_new_month

  !-------------------------------------------------------------------------------
  function its_a_new_year(this)result(yesno)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: yesno

  yesno = (this%moy == 1 .and. this%dom == 1 .and. this%tod < 1.e-3_r8)

  end function its_a_new_year
  !-------------------------------------------------------------------------------

  function its_time_to_histflush(this)result(yesno)

  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: yesno

  if(this%hist_freq<0)then
    !by time steps
    yesno = (mod(this%nelapstep,this%hist_freq)==0)
  elseif(this%hist_freq==0)then
    !no history file output until the last time step
    yesno=this%its_time_to_exit()
  elseif(this%hist_freq==1)then
    !by day
    yesno=this%its_a_new_day()
  elseif(this%hist_freq==30)then
    !by month
    yesno=this%its_a_new_month()
  elseif(this%hist_freq==365)then
    !by year
    yesno=this%its_a_new_year()
  elseif(this%hist_freq==9999)then
    yesno=.false.
  endif
  end function its_time_to_histflush
  !-------------------------------------------------------------------------------
  function get_curr_year(this)result(ans)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  integer :: ans
  ans = this%cyears
  end function get_curr_year
  !-------------------------------------------------------------------------------
  function get_curr_yearAD(this)result(ans)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  integer :: ans
  ans = this%cyears+this%year0
  end function get_curr_yearAD
  !-------------------------------------------------------------------------------
  function get_curr_mon(this)result(ans)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  integer :: ans

  ans = this%moy
  return
  end function get_curr_mon

  !-------------------------------------------------------------------------------
  function get_curr_day(this)result(ans)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  integer :: ans
  ans = this%cdays
  end function get_curr_day
  !-------------------------------------------------------------------------------
  function is_first_step(this)result(ans)
  !
  !is it the beginning of the first time step?
  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: ans
  ans = (this%nelapstep==0)
  end function is_first_step

  !-------------------------------------------------------------------------------
  subroutine print_model_time_stamp(this, iulog)

  implicit none
  class(ecosim_time_type), intent(in) :: this
  integer, intent(in) :: iulog

  select case(this%stop_opt)
  case (1)
    write(iulog,*)'step', this%get_nstep()
  case (2)
    if (this%its_a_new_day()) then
      write(iulog,*)'day', this%get_curr_day()
    end if
  case (3)
    if (this%its_a_new_day()) then
      write(iulog,*)'month', this%get_curr_mon()
    end if
  case (4)
    if (this%its_a_new_year()) then
      write(iulog,*)'year', this%get_curr_year()
    end if
  end select

  end subroutine print_model_time_stamp

!------------------------------------------------------------------------------------------

  function isLeap(year)result(ans)
!
! Description
! Determine if it is a leap year

  implicit none
  integer, intent(in) :: year
  logical :: ans

  ans =(mod(year,400)== 0) .or. (mod(year,4)==0 .and. mod(year,100)/=0)
  end function isLeap
!------------------------------------------------------------------------------------------
  function isLeapi(year)result(ans)
  implicit none
  integer, intent(in) :: year

  integer :: ans

  ans=0
  if(year<=0)then
    ans=0
  else if(isLeap(year))then
    ans=1
  endif
  return
  end function isLeapi
!------------------------------------------------------------------------------------------

  function getdow(year,doy)result(ans)
  !
  !get day of week
  implicit none
  integer, intent(in) :: year    !inquired year
  integer, intent(in) :: doy     !inquired day of year

  integer :: ans,nd,nn

  nd=0
  if (year<1900)then
    do nn=1900,year,-1
      if(isLeap(nn))then
        nd=mod(-366+nd,7)
      else
        nd=mod(-365+nd,7)
      endif
    enddo
    ans=mod(nd+doy+1,7)
    if(ans==0)ans=7
    if(ans<0)ans=ans+7
  elseif (year==1900)then
    ans=mod(doy,7)
    if(ans==0)ans=7
  else
    do nn=1900,year-1
      if(isLeap(nn))then
        nd=mod(366+nd,7)
      else
        nd=mod(365+nd,7)
      endif
    enddo
    ans=mod(nd+doy,7)
    if(ans==0)ans=7
  endif
  end function getdow

!------------------------------------------------------------------------------------------
  subroutine get_steps_from_ymdhs(ymdhs,dtime_sec,etime_dat,year0)

  implicit none
  character(len=*), intent(in) :: ymdhs
  integer, intent(in) :: dtime_sec
  type(ecosim_time_dat_type), intent(out) :: etime_dat
  integer, optional, intent(in) :: year0
  integer :: year,mon,day,hh,ss
  integer :: cdays1,cdays2,nn

  read(ymdhs,'(I4,I2,I2,I2,I4)')year,mon,day,hh,ss

  cdays1=sum(daz(1:mon-1))
  if(mon>2)cdays1=cdays1+isLeapi(year)
  cdays1=cdays1+day-1

  cdays2=cdays1
  if (present(year0))then
    etime_dat%year0=year0
    do nn=year0,year-1
      cdays1=cdays1+365+isLeapi(nn)
    enddo
  else
    etime_dat%year0=year
  endif

  etime_dat%nstep=int((cdays1*secspday+hh*3600._r8+ss)/dtime_sec)
  etime_dat%tstep=int((cdays2*secspday+hh*3600._r8+ss)/dtime_sec)

  end subroutine get_steps_from_ymdhs
  !-------------------------------------------------------------------------------
  subroutine get_prev_time(this,days, seconds)
  implicit none
  class(ecosim_time_type), intent(inout) :: this
  integer, intent(out) ::&
        days,   &! number of whole days in time interval
        seconds  ! remaining seconds in time interval
  real(r8) :: diff_time

  diff_time=this%curr_time-this%delta_time
  if(diff_time<0._r8)then
    days=0
    seconds=0
  else
    days=int(diff_time/86400._r8)
    seconds=diff_time-days*86400
  endif
  end subroutine get_prev_time
  !-------------------------------------------------------------------------------
  subroutine get_prev_date(this,yr, mon, day, tod)
  implicit none
  class(ecosim_time_type), intent(inout) :: this
  integer, intent(out) ::&
        yr,    &! year
        mon,   &! month
        day,   &! day of month
        tod     ! time of day (seconds past 0Z)

  tod = this%tod_prev
  day = this%dom_prev
  mon = this%moy_prev
  yr  = this%year_prev+this%year0

  end subroutine get_prev_date
  !-------------------------------------------------------------------------------
  subroutine get_curr_date(this,yr,mon,day,tod)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  integer, intent(out) ::&
        yr,    &! year
        mon,   &! month
        day,   &! day of month
        tod     ! time of day (seconds past 0Z)

  tod=this%tod
  day=this%dom
  mon=this%moy
  yr=this%cyears+this%year0

  end subroutine get_curr_date
  !-------------------------------------------------------------------------------

  integer function get_curr_doy(this)
  !
  !get current day of year
  implicit none
  class(ecosim_time_type), intent(in) :: this

  get_curr_doy=this%doy
  end function get_curr_doy

  !-----------------------------------------------------------------------
  integer function get_curr_dom(this)
  implicit none
  class(ecosim_time_type), intent(in) :: this

  get_curr_dom=this%dom
  end function get_curr_dom
  !-----------------------------------------------------------------------
  function get_calendar(this,prev)

  implicit none
  class(ecosim_time_type)  :: this
  logical, optional, intent(in) :: prev
  character(len=cal_str_len) :: get_calendar

  integer :: year,mon,day,hh,mm,ss,tod
  logical :: prev_loc

  prev_loc=.false.
  if(present(prev))prev_loc=prev
  if(prev_loc)then
    call this%get_prev_date(year, mon, day, tod)
  else
    call this%get_curr_date(year, mon, day, tod)
  endif
  hh=tod/3600
  ss=mod(tod,60)
  mm=(tod-hh*3600-ss)/60
  write(get_calendar,'(I4.4,I2.2,I2.2,I2.2,I2.2,I2.2)')year,mon,day,hh,mm,ss

  end function get_calendar

!-----------------------------------------------------------------------
  subroutine getdatetime (this, cdate, ctime)
!
! !DESCRIPTION:
! A generic Date and Time routine
!
! !USES:
!  use spmdMod      , only : mpicom, masterproc, MPI_CHARACTER
! !ARGUMENTS:
  implicit none
  class(ecosim_time_type)  :: this
  character(len=8), intent(out) :: cdate  !current date
  character(len=8), intent(out) :: ctime  !current time
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  character(len=8)      :: date       !current date
  character(len=10)     :: time       !current time
  character(len=5)      :: zone       !zone
  integer, dimension(8) :: values !temporary
  integer               :: ier    !MPI error code
!-----------------------------------------------------------------------
!  if (masterproc) then

    call date_and_time (date, time, zone, values)

    cdate(1:2) = date(5:6)
    cdate(3:3) = '/'
    cdate(4:5) = date(7:8)
    cdate(6:6) = '/'
    cdate(7:8) = date(3:4)

    ctime(1:2) = time(1:2)
    ctime(3:3) = ':'
    ctime(4:5) = time(3:4)
    ctime(6:6) = ':'
    ctime(7:8) = time(5:6)

!  endif
  end subroutine getdatetime
!-----------------------------------------------------------------------
  subroutine update_sim_len(this, nstopyr)
  implicit none
  class(ecosim_time_type)  :: this
  integer, intent(in) :: nstopyr

  if(this%stop_opt==4)this%stop_count=min0(nstopyr,this%stop_count)

  end subroutine update_sim_len
end module ecosim_Time_Mod
