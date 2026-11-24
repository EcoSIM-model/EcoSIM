module EcoSIMConfig
  use fileUtil, only :   datestrlen
implicit none
  character(len=*),private, parameter :: mod_filename =&
   __FILE__
   
  logical :: is_first_year=.false.
  logical :: transport_on=.true.
  logical :: column_mode=.false.
  logical :: do_instequil=.false.
  logical :: brnch_retain_casename = .false.
  integer :: iFlagRaiseZ0GbyVeg                 !flag for raising Z0G with vegn

  character(len=256), public :: rpntdir = '.'
  character(len=256), public :: rpntfil = 'rpointer.esim'
  character(len=16) , public :: inst_suffix=''

  integer, parameter :: NumDeadMicrbCompts = 2 !# of microbial residue components
  integer, parameter :: NumLiveMicrbCompts = 3 !# of living biomass components
  integer, parameter :: jskenc    = 4 !# of kinetic components of the substrates
  integer, parameter :: jskenp1   = 5
  integer, parameter :: jcplxc    = 5 !# of microbe-substrate complexes
  integer, parameter :: jcplxcm1   = jcplxc-1
  integer, parameter :: NumMicbHFunGrupsPerCmplx  = 8 !# of heterotrophic microbial functional groups in each complex
  integer, parameter :: NumMicbAFunGrupsPerCmplx  = 6 !# of autotrohpic microbial functional groups 
  
  character(len=datestrlen)  :: ref_date  = '18000101000000'
  character(len=datestrlen)  :: start_date= '18000101000000'
  character(len=36)  :: case_name
  character(len=256) :: finidat =' '
  character(len=256) :: ctitle  = ' '
  character(len=256) :: restartFileFullPath 
  character(len=256), public :: hostname = ' '
  character(len=256), public :: username = ' '          ! Current user
  ! description of this source
  character(len=256), public :: source   = "ECOSIM"
  ! version of program
  character(len=256), public :: version  = " "
 
  integer, public :: nsrest  
  integer, public, parameter :: nsrBranch   = 2        ! Branch from restart files 
  integer, public, parameter :: nsrStartup  = 0        ! Startup from initial conditions
  integer, public, parameter :: nsrContinue = 1        ! Continue from restart files
  public :: is_restart,   &
            is_branch,    &
            set_sim_type, &
            is_cold_run
contains

  logical function is_restart( )
  !
  ! Determine if restart run
  implicit none
  
  if (nsrest == nsrContinue) then
    is_restart = .true.
  else
    is_restart = .false.
  end if
  end function is_restart
!-----------------------------------------------------------------------

  logical function is_branch( )
  !
  ! Determine if it is branch run
  implicit none
  
  if (nsrest == nsrBranch) then
    is_branch = .true.
  else
    is_branch = .false.
  end if
  end function is_branch

!-----------------------------------------------------------------------

  subroutine set_sim_type()
  use EcoSIMCtrlMod, only : continue_run
  implicit none
  !determine simulation type

  if(continue_run)then
    nsrest=nsrContinue
  else
    if(finidat ==' ')then
      nsrest=nsrStartup
    else
      nsrest=nsrBranch
    endif
  endif
  end subroutine set_sim_type
!-----------------------------------------------------------------------

logical function is_cold_run()
implicit none

is_cold_run=(nsrest==nsrStartup)

end function is_cold_run
end module EcoSIMConfig
