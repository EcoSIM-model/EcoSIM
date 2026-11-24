module TestMod
! DESCRIPTION
! codes to do regression tests
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : endrun
  use fileUtil, only : error_errmsg_len
implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  save

  integer, parameter :: ecosys_filename_length=256

  integer, parameter :: stdout = 6
  integer, parameter :: ecosim_var_name_length = 36
  integer, parameter :: ecosim_string_length = 128
  integer, parameter :: ecosim_string_length_long = 256
  integer, parameter :: ecosim_namelist_buffer_size = 4096
  integer, parameter :: ecosim_namelist_buffer_size_ext = 12288
  real(r8), parameter :: tiny_val=1.e-50_r8
  type, public :: error_status_type
    integer :: error
    character(len=error_errmsg_len) :: msg
  contains
    procedure, public :: reset
    procedure, public :: set_msg
    procedure, public :: check_status
    procedure, public :: print_err
    procedure, public :: print_msg
  end type error_status_type


  type, public :: ecosys_regression_type
     character(len=ecosys_filename_length), private :: filename
     integer :: output

     logical, public  :: write_regression_output
     integer, private :: num_cells
     ! FIXME(bja, 201603) specifying cell ids requires more careful
     !thought. Maybe just hard code a max....

     !X!integer, allocatable, private :: cell_ids(:)
   contains
     procedure, public  :: Init
     procedure, public  :: OpenOutput
     procedure, public  :: CloseOutput
     procedure, public  :: WriteData
     procedure, private :: ReadNamelist
     procedure, private :: CheckInput
  end type ecosys_regression_type

  type(ecosys_regression_type), public :: regression
  public :: create_error_status_type
  public :: errMsg
contains

    subroutine Init(this, namelist_file, case_name)

    implicit none

    class(ecosys_regression_type)              , intent(inout) :: this
    character(len=*), intent(in)     :: namelist_file
    character(len=*), intent(in)     :: case_name

    !local variables
    type(error_status_type) :: err_status
    this%output = 16
    call this%ReadNamelist(trim(namelist_file),err_status)
    if(err_status%check_status())then
      call endrun('stopped in '//mod_filename, __LINE__)
    endif
    call this%CheckInput(err_status)
    this%filename = trim(case_name) // '.regression'
  end subroutine Init



  subroutine ReadNamelist(this, namelist_file, err_status)
    implicit none

    class(ecosys_regression_type), intent(inout) :: this
    character(len=*), intent(in) :: namelist_file
    class(error_status_type)        , intent(out)   :: err_status

    character(len=*), parameter :: subname = 'ecosys_regression:ReadNamelist'
    ! !LOCAL VARIABLES:
    integer :: nml_error
    integer :: cells
    character(len=ecosim_string_length_long) :: ioerror_msg
    logical :: write_regression_output

    integer :: fu, rc
    namelist / regression_test / cells, write_regression_output


    call err_status%reset()
    cells = 0
    write_regression_output = .false.
    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------
    inquire (file=namelist_file, iostat=rc)
    if (rc /= 0) then
       write (stdout, '(3a)') 'Error: input file "', trim(namelist_file), '" does not exist.'
       call endrun('stopped in '//mod_filename, __LINE__)
     end if

     open (action='read', file=trim(namelist_file), iostat=rc, newunit=fu)
     if (rc /= 0) then
        write (*, '(3a)') 'Error openning input file "', trim(namelist_file)
        call endrun('stopped in '//mod_filename, __LINE__)
     end if

    if ( rc==0 )then
       ioerror_msg=''
       read(unit=fu, nml=regression_test, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
         call err_status%set_msg(msg="ERROR reading ecosys_regression_test namelist "//errmsg(mod_filename, __LINE__),err=-1)
         call endrun('stopped in '//mod_filename, __LINE__)
       end if
    end if

    close(fu)
    if (.true.) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' ecosim regression test type :'
       write(stdout, *)
       write(stdout, *) ' regression_test namelist settings :'
       write(stdout, *)
       write(stdout, regression_test)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif

    this%write_regression_output = write_regression_output
    this%num_cells = cells
  end subroutine ReadNamelist

  !---------------------------------------------------------------------------------

  subroutine CheckInput(this, err_status)

    implicit none

    class(ecosys_regression_type), intent(in) :: this
    class(error_status_type)        , intent(out)   :: err_status
    character(len=ecosim_string_length_long) :: msg

    if (this%write_regression_output) then
       ! sanity check on input
       if (this%num_cells < 0) then
          msg = 'ERROR num cells must be >= 0. '
          call err_status%set_msg(msg=msg//errmsg(mod_filename, __LINE__),err=-1)
          return
       end if
    end if
  end subroutine CheckInput

  !---------------------------------------------------------------------------------

  subroutine OpenOutput(this)


    implicit none

    class(ecosys_regression_type), intent(inout) :: this

    integer :: output = 16

    write(stdout, '(a, a)') 'Writing regression output to ', this%filename

    open(this%output, file=this%filename, status='REPLACE')

  end subroutine OpenOutput


  !---------------------------------------------------------------------------------

  subroutine CloseOutput(this)


    implicit none

    class(ecosys_regression_type), intent(inout) :: this

    close(this%output)

  end subroutine CloseOutput


  !---------------------------------------------------------------------------------

  subroutine WriteData(this, category, name, data)
    !
    ! Write regression data in a cfg/ini format that can easily be
    ! parsed by external processors.
    !
    ! sections names are the species name (and possible
    ! phase/location?)
    !
    ! section data is just keyword = value for things like min, max,
    ! mean, vector norms, cell point data.


    implicit none

    class(ecosys_regression_type)         , intent(inout) :: this
    character(len=ecosim_string_length)   , intent(in)    :: category
    character(len=ecosim_var_name_length) , intent(in)    :: name
    real(r8)                            , intent(in)    :: data(:)

    integer :: cell_increment, num_cells, cell
    real(r8) :: val, local_val

    ! FIXME(bja, 201603) cfg/ini format limits the characters that can
    ! be use in section names. We need to sanitize the names!
    write(this%output, '("[",a,"]")') trim(name)

    write(this%output, '("category = ",a)') trim(category)

    val = minval(data(:))
    if(abs(val)<tiny_val)val=0.
    write(this%output, '("min = ",e21.13)') val

    val = maxval(data(:))
    if(abs(val)<tiny_val)val=0.
    write(this%output, '("max = ",e21.13)') val

    val = sum(data(:)) / size(data)
    if(abs(val)<tiny_val)val=0.
    write(this%output, '("mean = ",e21.13)') val

    if (this%num_cells > 0) then
       num_cells = this%num_cells
       if (num_cells > size(data)) then
          ! need to truncate
          num_cells = size(data)
       end if

       cell_increment = int(size(data) / num_cells)

       do cell = 1, size(data), cell_increment
          local_val = data(cell)
          if(abs(local_val)<tiny_val)local_val=0.
          write(this%output, '("cell ", i4, " = ", e21.13)') cell, local_val
       end do
    write(this%output, '(a)')
    end if
  end subroutine WriteData




  function create_error_status_type()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(error_status_type), pointer :: create_error_status_type
    class(error_status_type), pointer :: error_status

    allocate(error_status)
    create_error_status_type => error_status

  end function create_error_status_type
!-------------------------------------------------------------------------------
  subroutine reset(this)
  implicit none
  class(error_status_type)  :: this

  this%msg = ''
  this%error = 0
  end subroutine reset
!-------------------------------------------------------------------------------
  subroutine set_msg(this, msg, err, c)
  implicit none
  class(error_status_type), intent(inout)  :: this
  character(len=*), intent(in) :: msg
  integer, intent(in) :: err
  integer, optional, intent(in) :: c

  if(present(c))then
    write(this%msg,'(A,A,I6.6)')trim(msg),' col=',c
  else
    write(this%msg,'(A)')trim(msg)
  endif
  this%error = err
  end subroutine set_msg
!-------------------------------------------------------------------------------
  function check_status(this)result(ans)
  implicit none
  class(error_status_type), intent(in)  :: this
  logical :: ans
  if(this%error < 0)then
    ans = .true.
  else
    ans = .false.
  endif
  end function check_status
!-------------------------------------------------------------------------------
  function print_msg(this)
  implicit none
  class(error_status_type), intent(in)  :: this
  character(len=error_errmsg_len) :: print_msg

  print_msg = this%msg
  end function print_msg

!-------------------------------------------------------------------------------
  function print_err(this)
  implicit none
  class(error_status_type), intent(in)  :: this
  integer :: print_err

  print_err = this%error
  end function print_err

!-------------------------------------------------------------------------------

  function errMsg(file, line)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(len=256)   :: errMsg
  character(len=*), intent(in) :: file
  integer         , intent(in) :: line

!EOP

  write(errMsg, '(a, a, a, i0)') 'ERROR in ', trim(file), ' at line ', line

  end function errMsg

end module TestMod
