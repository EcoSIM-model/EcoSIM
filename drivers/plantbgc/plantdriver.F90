program main
!
  use abortutils, only : endrun,iulog
  use fileUtil
  use PlantTestMod
implicit none
  character(len=*), parameter :: mod_filename = &
  __FILE__

  integer :: arg_count
  character(len=ecosim_filename_length)      :: namelist_filename
  character(len=ecosim_namelist_buffer_size) :: namelist_buffer
  character(len=ecosim_filename_length)      :: base_filename

  arg_count = command_argument_count()
  if (arg_count /= 1) then
     write(stdout, '(a, i3)') 'ERROR: must pass exactly one arguement one the command line, received ', arg_count
     call usage()
     call endrun()
  end if

  call get_command_argument(1, namelist_filename)

  write(stdout, '(a, a)') 'Reading namelist filename : ', trim(namelist_filename)
  namelist_buffer = ''

  call namelist_to_buffer(namelist_filename, namelist_buffer)

  call RunModel(namelist_buffer)

end program main



! ----------------------------------------------------------------------
subroutine usage()
  !DESCRIPTION
  !display something
  use fileUtil, only : stdout
 implicit none
  write(stdout, *) 'plantmodel - standalone driver for plantmodel.'
  write(stdout, *) 'usage: plantmodel namelist_filename'
end subroutine usage


! ----------------------------------------------------------------------

subroutine RunModel(namelist_buffer)
  use ecosim_Time_Mod, only: ecosim_time_type
  use ModelStatusType, only: model_status_type
  use data_kind_mod,   only: r8 => DAT_KIND_R8
  use abortutils,      only: iulog
  use PlantTestMod
  use bhistMod
  implicit none
  character(len=*), intent(in) :: namelist_buffer

  character(len=hist_var_str_len) , allocatable :: varl(:)
  character(len=hist_var_lon_str_len) , allocatable :: varlnml(:)
  character(len=hist_unit_str_len), allocatable :: unitl(:)
  character(len=hist_freq_str_len), allocatable :: freql(:)
  integer                         , allocatable :: vartypes(:)
  real(r8), allocatable :: ystatesf(:,:)
  real(r8), allocatable :: ystates0l(:)
  real(r8), allocatable :: ystatesfl(:)
  integer :: nvars,ncols,jj
  real(r8) :: dtime
  type(ecosim_time_type) :: timer
  type(model_status_type) :: err_status
  type(histf_type) :: hist
  character(len=256) :: ioerror_msg
  character(len=256):: gname
  character(len=14) :: yymmddhhss
  integer :: rc, fu
  integer :: nml_error
  character(len=hist_freq_str_len) :: hist_freq
  character(len=32) :: model_name
  character(len=32) :: case_id
  namelist /driver_nml/model_name,case_id,hist_freq
  character(len=*), parameter :: mod_filename=&
  __FILE__
  hist_freq='day'
  model_name='plant'
  case_id='exp0'
  if ( .true. )then
     ioerror_msg=''
     read(namelist_buffer, nml=driver_nml, iostat=nml_error, iomsg=ioerror_msg)
     if (nml_error /= 0) then
        write(*,*)ioerror_msg
        call endrun(msg="ERROR reading driver_nml namelist "//errmsg(mod_filename, __LINE__))
     end if
  end if
  if (.true.) then
    write(stdout, *)
    write(stdout, *) '--------------------'
    write(stdout,driver_nml)
    write(stdout, *)
    write(stdout, *) '--------------------'
  endif

  ncols=1
!!============================================================
! customized model varlist creation
  nvars=getvarllen()
  allocate(varl(nvars)); allocate(varlnml(nvars));
  allocate(unitl(nvars)); allocate(freql(nvars)); allocate(vartypes(nvars))
  call getvarlist(nvars, varl, varlnml, unitl, vartypes)
!============================================================


  freql(:) = hist_freq
  allocate(ystatesf(1:ncols,1:nvars));ystatesf(:,:)=0._r8
  allocate(ystates0l(1:nvars));ystates0l(:) = 0._r8
  allocate(ystatesfl(1:nvars));ystatesfl(:) = 0._r8
!!============================================================
! customized model initialization
  call initmodel(nvars, ystates0l, err_status)

  if(err_status%check_status())then
    call endrun(msg=err_status%print_msg())
  endif

!============================================================
  !initialize timer
  call timer%Init(namelist_buffer=namelist_buffer)
  dtime=timer%get_step_size()

! initialize history file
  if(len(trim(case_id))==0)then
    write(gname,'(A)')'plantmodel'//'.'//trim(model_name)
  else
    write(gname,'(A)')'plantmodel'//'.'//trim(case_id)//'.'//trim(model_name)
  endif
  call hist%init(ncols, varl, varlnml, unitl, vartypes, freql, gname, dtime)
  !print*,'run the model'
  do

    call timer%update_time_stamp()

!!============================================================
! customized model run
    call runplant(nvars,ystates0l, ystatesfl, err_status)

    if(err_status%check_status())then
      call endrun(msg=err_status%print_msg())
    endif
!============================================================

    !print*,'hist_wrap'
    do jj = 1, nvars
      ystatesf(1,jj)=ystatesfl(jj)
    enddo

    call hist%hist_wrap(ystatesf, timer)

    if(timer%its_a_new_year())then
      write(iulog,*)'year ',timer%get_curr_year()
    endif
    if(timer%its_time_to_exit())exit
  enddo
  call timer%get_ymdhs(yymmddhhss)

  call hist%histrst('plantmodel', 'write', yymmddhhss)

end subroutine RunModel
