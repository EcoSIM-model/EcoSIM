module abortutils

  !-----------------------------------------------------------------------
  ! !MODULE: abortutils
  !
  ! !DESCRIPTION:
  ! Abort the model for abnormal termination
  !-----------------------------------------------------------------------
  use data_kind_mod, only : r8 => DAT_KIND_R8

  implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename =&
   __FILE__
  integer, public, parameter :: iulog =6   
  public :: endrun
  public :: check_bool
  interface endrun
     module procedure endrun_vanilla
     module procedure endrun_line
     module procedure endrun_globalindex
  end interface
  public :: destroy
  interface destroy
    module procedure destroy_1d_r
    module procedure destroy_2d_r
    module procedure destroy_3d_r
    module procedure destroy_4d_r
    module procedure destroy_5d_r
    module procedure destroy_6d_r
    module procedure destroy_7d_r
    module procedure destroy_1d_i
    module procedure destroy_2d_i
    module procedure destroy_3d_i
    module procedure destroy_4d_i
    module procedure destroy_5d_i
    module procedure destroy_6d_i
    module procedure destroy_7d_i
    module procedure destroy_1d_l
    module procedure destroy_2d_l
    module procedure destroy_3d_l
    module procedure destroy_4d_l
    module procedure destroy_1d_char
    module procedure destroy_2d_char
    module procedure destroy_3d_char
    module procedure destroy_4d_char

    module procedure destroy_1d_rpt
    module procedure destroy_2d_rpt
    module procedure destroy_3d_rpt
    module procedure destroy_4d_rpt
    module procedure destroy_5d_rpt
    module procedure destroy_6d_rpt
    module procedure destroy_7d_rpt
    module procedure destroy_1d_ipt
    module procedure destroy_2d_ipt
    module procedure destroy_3d_ipt
    module procedure destroy_4d_ipt
    module procedure destroy_5d_ipt
    module procedure destroy_6d_ipt
    module procedure destroy_7d_ipt
    module procedure destroy_1d_lpt
    module procedure destroy_1d_charpt
    module procedure destroy_2d_charpt
    module procedure destroy_3d_charpt
    module procedure destroy_4d_charpt


  end interface destroy
  public :: padl, padr
  public :: print_info
  interface print_info
    module procedure print_info_arr
    module procedure print_info_msg
  end interface

CONTAINS
  !-----------------------------------------------------------------------
  subroutine endrun_vanilla(msg)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    !
    !
    ! !ARGUMENTS:

    implicit none
    character(len=*), intent(in), optional :: msg    ! string to be printed
    !-----------------------------------------------------------------------

    if (present (msg)) then
       write(iulog,*)'ENDRUN:', msg
    else
       write(iulog,*)'ENDRUN: called without a message string'
    end if

    stop

  end subroutine endrun_vanilla

  !-----------------------------------------------------------------------
  subroutine endrun_line(msg,line)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    !
    !
    ! !ARGUMENTS:

    implicit none
    character(len=*), intent(in) :: msg    ! string to be printed
    integer, intent(in) :: line
    !-----------------------------------------------------------------------

    write(iulog,*)'ENDRUN:', msg, 'at line ',line
    stop

  end subroutine endrun_line

  !-----------------------------------------------------------------------
  subroutine endrun_globalindex(decomp_index, elmlevel, msg)

    !-----------------------------------------------------------------------
    ! Description:
    ! Abort the model for abnormal termination
    !

    !
    ! Arguments:
    implicit none
    integer          , intent(in)           :: decomp_index
    character(len=*) , intent(in)           :: elmlevel
    character(len=*) , intent(in), optional :: msg    ! string to be printed
    !
    ! Local Variables:
    integer :: igrc, ilun, icol
    !-----------------------------------------------------------------------

    write(6,*)'calling getglobalwrite with decomp_index= ',decomp_index,' and elmlevel= ',trim(elmlevel)


    if (present (msg)) then
       write(iulog,*)'ENDRUN:', msg
    else
       write(iulog,*)'ENDRUN: called without a message string'
    end if

    stop

  end subroutine endrun_globalindex

  !-----------------------------------------------------------------------

  function padl(str,width,symbol) result(padded_str)
    !-----------------------------------------------------------------------
    ! Description:
    ! pad a string from left with symbol to width
    !
    !
    ! !ARGUMENTS:

    implicit none
    character(len=*), intent(in) :: str
    integer, intent(in) :: width
    character(len=1), intent(in), optional :: symbol
    character(len=max(len(str),width)) :: padded_str

    integer :: ns, jj, kk, kk1
    character(len=1) :: symboll

    if(present(symbol))then
      symboll=symbol
    else
      symboll=' '
    endif
    ns = len(str)
    if (ns >= width)then
      padded_str=str
    else
      do jj = width, 1, -1
        kk = width-jj
        if (kk <= ns-1) then
          padded_str(jj:jj) = str(ns-kk:ns-kk)
        else
          padded_str(jj:jj) = symboll
        endif
      enddo
    endif
    return
  end function padl
  !-----------------------------------------------------------------------

  ! right - padding
  function padr(str,width,symbol) result(padded_str)
    !-----------------------------------------------------------------------
    ! Description:
    ! pad a string from right with symbol to width
    !
    !
    ! !ARGUMENTS:

    implicit none
    character(len=*), intent(in) :: str
    integer, intent(in) :: width
    character(len=1), intent(in), optional :: symbol
    character(len=max(len(str),width)) :: padded_str

    integer :: ns, jj
    character(len=1) :: symboll

    if(present(symbol))then
      symboll=symbol
    else
      symboll=' '
    endif
    ns = len(str)
    if (ns >= width)then
      padded_str=str
    else
      do jj = 1, width
        if (jj <= ns) then
          padded_str(jj:jj) = str(jj:jj)
        else
          padded_str(jj:jj) = symboll
        endif
      enddo
    endif
    return

  end function padr

  !-----------------------------------------------------------------------

  subroutine print_info_arr(msg,strarr,valarr)

    !-----------------------------------------------------------------------
    ! Description:
    ! print information
    !
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: msg
    character(len=*), dimension(:), intent(in) :: strarr(:)
    real(r8), dimension(:), intent(in) :: valarr

    integer :: ns1
    integer :: jj

    ns1=size(valarr)
    write(*,*)trim(msg)
    do jj = 1 , ns1
      write(*,*)trim(strarr(jj))//'=',valarr(jj)
    enddo
  end subroutine print_info_arr

  !-----------------------------------------------------------------------
  subroutine print_info_msg(msg,lineno)

  implicit none
  character(len=*), intent(in) :: msg
  integer, intent(in) :: lineno

  write(*,*)''
  write(*,*)trim(msg)
  write(*,*)'at line',lineno
  end subroutine print_info_msg

  !-----------------------------------------------------------------------
  subroutine check_bool(bool_expr, msg, lineno, modfile)
  implicit none
  logical, intent(in) :: bool_expr
  character(len=*), intent(in) :: msg
  integer, intent(in) :: lineno
  character(len=*), intent(in) :: modfile

  if(bool_expr)then
    write(iulog,*)msg, 'at line',lineno, ' in file ',modfile
    call endrun()
  endif
  end subroutine check_bool

  !-----------------------------------------------------------------------
  subroutine destroy_1d_r(arr)
  implicit none
  real(r8), allocatable :: arr(:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_1d_r
  !-----------------------------------------------------------------------
  subroutine destroy_2d_r(arr)
  implicit none
  real(r8), allocatable :: arr(:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_2d_r
  !-----------------------------------------------------------------------
  subroutine destroy_3d_r(arr)
  implicit none
  real(r8), allocatable :: arr(:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_3d_r
  !-----------------------------------------------------------------------
  subroutine destroy_4d_r(arr)
  implicit none
  real(r8), allocatable :: arr(:,:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_4d_r
  !-----------------------------------------------------------------------
  subroutine destroy_5d_r(arr)
  implicit none
  real(r8), allocatable :: arr(:,:,:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_5d_r
  !-----------------------------------------------------------------------
  subroutine destroy_6d_r(arr)
  implicit none
  real(r8), allocatable :: arr(:,:,:,:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_6d_r
  !-----------------------------------------------------------------------
  subroutine destroy_7d_r(arr)
  implicit none
  real(r8), allocatable :: arr(:,:,:,:,:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_7d_r
  !-----------------------------------------------------------------------
  subroutine destroy_1d_i(arr)
  implicit none
  integer, allocatable :: arr(:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_1d_i
  !-----------------------------------------------------------------------
  subroutine destroy_2d_i(arr)
  implicit none
  integer, allocatable :: arr(:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_2d_i
  !-----------------------------------------------------------------------
  subroutine destroy_3d_i(arr)
  implicit none
  integer, allocatable :: arr(:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_3d_i
  !-----------------------------------------------------------------------
  subroutine destroy_4d_i(arr)
  implicit none
  integer, allocatable :: arr(:,:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_4d_i
  !-----------------------------------------------------------------------
  subroutine destroy_5d_i(arr)
  implicit none
  integer, allocatable :: arr(:,:,:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_5d_i
  !-----------------------------------------------------------------------
  subroutine destroy_6d_i(arr)
  implicit none
  integer, allocatable :: arr(:,:,:,:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_6d_i
  !-----------------------------------------------------------------------
  subroutine destroy_7d_i(arr)
  implicit none
  integer, allocatable :: arr(:,:,:,:,:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_7d_i
  !-----------------------------------------------------------------------
  subroutine destroy_1d_l(arr)
  implicit none
  logical, allocatable :: arr(:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_1d_L
  !-----------------------------------------------------------------------
  subroutine destroy_2d_l(arr)
  implicit none
  logical, allocatable :: arr(:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_2d_L
  !-----------------------------------------------------------------------
  subroutine destroy_3d_l(arr)
  implicit none
  logical, allocatable :: arr(:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_3d_L
  !-----------------------------------------------------------------------
  subroutine destroy_4d_l(arr)
  implicit none
  logical, allocatable :: arr(:,:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_4d_L


  !-----------------------------------------------------------------------
  subroutine destroy_1d_char(arr)
  implicit none
  character(len=*), allocatable :: arr(:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_1d_char

  !-----------------------------------------------------------------------
  subroutine destroy_2d_char(arr)
  implicit none
  character(len=*), allocatable :: arr(:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_2d_char
  !-----------------------------------------------------------------------
  subroutine destroy_3d_char(arr)
  implicit none
  character(len=*), allocatable :: arr(:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_3d_char

  !-----------------------------------------------------------------------
  subroutine destroy_4d_char(arr)
  implicit none
  character(len=*), allocatable :: arr(:,:,:,:)

  if(allocated(arr))deallocate(arr)

  end subroutine destroy_4d_char


  !-----------------------------------------------------------------------
  subroutine destroy_1d_rpt(arr)
  implicit none
  real(r8), pointer :: arr(:)

  deallocate(arr)

  end subroutine destroy_1d_rpt
  !-----------------------------------------------------------------------
  subroutine destroy_2d_rpt(arr)
  implicit none
  real(r8), pointer :: arr(:,:)

  deallocate(arr)

  end subroutine destroy_2d_rpt
  !-----------------------------------------------------------------------
  subroutine destroy_3d_rpt(arr)
  implicit none
  real(r8), pointer :: arr(:,:,:)

  deallocate(arr)

  end subroutine destroy_3d_rpt
  !-----------------------------------------------------------------------
  subroutine destroy_4d_rpt(arr)
  implicit none
  real(r8), pointer :: arr(:,:,:,:)

  deallocate(arr)

  end subroutine destroy_4d_rpt
  !-----------------------------------------------------------------------
  subroutine destroy_5d_rpt(arr)
  implicit none
  real(r8), pointer :: arr(:,:,:,:,:)

  deallocate(arr)

  end subroutine destroy_5d_rpt
  !-----------------------------------------------------------------------
  subroutine destroy_6d_rpt(arr)
  implicit none
  real(r8), pointer :: arr(:,:,:,:,:,:)

  deallocate(arr)

  end subroutine destroy_6d_rpt
  !-----------------------------------------------------------------------
  subroutine destroy_7d_rpt(arr)
  implicit none
  real(r8), pointer :: arr(:,:,:,:,:,:,:)

  deallocate(arr)

  end subroutine destroy_7d_rpt
  !-----------------------------------------------------------------------
  subroutine destroy_1d_ipt(arr)
  implicit none
  integer, pointer :: arr(:)

  deallocate(arr)

  end subroutine destroy_1d_ipt
  !-----------------------------------------------------------------------
  subroutine destroy_2d_ipt(arr)
  implicit none
  integer, pointer :: arr(:,:)

  deallocate(arr)

  end subroutine destroy_2d_ipt
  !-----------------------------------------------------------------------
  subroutine destroy_3d_ipt(arr)
  implicit none
  integer, pointer :: arr(:,:,:)

  deallocate(arr)

  end subroutine destroy_3d_ipt
  !-----------------------------------------------------------------------
  subroutine destroy_4d_ipt(arr)
  implicit none
  integer, pointer :: arr(:,:,:,:)

  deallocate(arr)

  end subroutine destroy_4d_ipt
  !-----------------------------------------------------------------------
  subroutine destroy_5d_ipt(arr)
  implicit none
  integer, pointer :: arr(:,:,:,:,:)

  deallocate(arr)

  end subroutine destroy_5d_ipt
  !-----------------------------------------------------------------------
  subroutine destroy_6d_ipt(arr)
  implicit none
  integer, pointer :: arr(:,:,:,:,:,:)

  deallocate(arr)

  end subroutine destroy_6d_ipt
  !-----------------------------------------------------------------------
  subroutine destroy_7d_ipt(arr)
  implicit none
  integer, pointer :: arr(:,:,:,:,:,:,:)

  deallocate(arr)

  end subroutine destroy_7d_ipt
  !-----------------------------------------------------------------------
  subroutine destroy_1d_lpt(arr)
  implicit none
  logical, pointer :: arr(:)

  deallocate(arr)

  end subroutine destroy_1d_Lpt


  !-----------------------------------------------------------------------
  subroutine destroy_1d_charpt(arr)
  implicit none
  character(len=*), pointer :: arr(:)

  deallocate(arr)

  end subroutine destroy_1d_charpt

  !-----------------------------------------------------------------------
  subroutine destroy_2d_charpt(arr)
  implicit none
  character(len=*), pointer :: arr(:,:)

  deallocate(arr)

  end subroutine destroy_2d_charpt
  !-----------------------------------------------------------------------
  subroutine destroy_3d_charpt(arr)
  implicit none
  character(len=*), pointer :: arr(:,:,:)

  deallocate(arr)

  end subroutine destroy_3d_charpt

  !-----------------------------------------------------------------------
  subroutine destroy_4d_charpt(arr)
  implicit none
  character(len=*), pointer :: arr(:,:,:,:)

  deallocate(arr)

  end subroutine destroy_4d_charpt



end module abortutils
