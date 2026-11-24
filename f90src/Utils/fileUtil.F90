module fileUtil
!!
! subroutines for file open with error check
  use abortutils, only : endrun,iulog
  use data_kind_mod
  implicit none
  private
  character(len=*),  parameter :: mod_filename = &
  __FILE__
  public :: open_safe
  public :: check_read
  public :: remove_filename_extension
  public :: file_exists
  public :: getfil,getavu
  public :: int2str
  public :: strip_null,print_ichar,strip_space
  public :: namelist_to_buffer
  public :: opnfil,relavu
  integer, public, parameter :: ecosim_filename_length=128
  integer, public, parameter :: stdout=6
  integer, public, parameter :: error_errmsg_len=256
  integer, public, parameter :: ecosim_string_length_long=256
  integer, public, parameter :: var_flux_type =1
  integer, public, parameter :: var_state_type=2
  integer, public, parameter :: ecosim_namelist_buffer_size = 4096
  logical, public, save :: continue_run = .false.
  integer, public, parameter :: datestrlen=14
  integer(DAT_KIND_IN) :: s_loglev = 0
  integer(DAT_KIND_IN),parameter :: file_maxUnit = 99 
  integer(DAT_KIND_IN),parameter :: file_minUnit = 10  
  logical, save :: UnitTag(0:file_maxUnit) = .false. ! Logical units in use
  contains
!------------------------------------------------------------------------------------------

  function file_exists(filename) result(res)
!
!! DESCRIPTION
! check existence of a file
  implicit none
  character(len=*),intent(in) :: filename
  logical                     :: res

  ! Check if the file exists
  inquire( file=trim(filename), exist=res )
  end function
!------------------------------------------------------------------------------------------

  subroutine open_safe(lun,prefix,fname,status,location,lineno,lverb)
!
! Description
! safely open a file for IO
  implicit none
  integer, intent(in) :: lun
  character(len=*), intent(in) :: prefix
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: status
  character(len=*), intent(in) :: location
  integer, intent(in) :: lineno
  logical, optional, intent(in) :: lverb
  character(len=2560) :: pathfile
  integer :: ierr

  call getfilenamef(prefix,fname,pathfile)

  if(trim(status)/='UNKNOWN' .and. .not.file_exists(pathfile))then
    call endrun(msg='Fail to find file '//trim(pathfile)//' in ' &
      //mod_filename,line=__LINE__)
  endif

  if(present(lverb))then
    if(lverb)write(*,*)'open file ',trim(pathfile)
  endif
  OPEN(lun,FILE=trim(pathfile),STATUS=status,iostat=ierr)
  if(ierr/=0)then
    call endrun(msg="error in "//location//" while reading file " &
      //TRIM(pathfile), line=lineno)
  endif
  end subroutine open_safe

!------------------------------------------------------------------------------------------

  subroutine check_read(ios,nelm,ifile,lineno,modfile)
  implicit none
  integer, intent(in) :: ios
  integer, intent(in) :: nelm
  character(len=*),intent(in) :: ifile
  integer, intent(in) :: lineno
  character(len=*), intent(in) :: modfile

  if(ios>0)then
    write(*,*)'input error'
    call endrun('read error when reading '//trim(ifile)//' in'//trim(modfile),lineno)
  elseif(ios<0)then
    write(*,*)'number of inputs is less than ',nelm
    call endrun('read error when reading '//trim(ifile)//' in'//trim(modfile),lineno)
  endif
  end subroutine check_read

!------------------------------------------------------------------------------------------

  function remove_filename_extension(filename) result(basename)
    !
    ! Remove the extension from a file name to get a base filename.
    !
    ! We start at the end of the filename and assume that the extension
    ! is marked by a period.
    !
    implicit none

    character(len=ecosim_filename_length), intent(in) :: filename
    character(len=ecosim_filename_length) :: basename
    integer :: ext_index

    ext_index = scan(filename, '.', .true.)
    if (ext_index == 0) then
       ! no period marking an extension...
       ext_index = len(trim(filename)) + 1

    end if
    basename = filename(1:ext_index-1)
  end function remove_filename_extension

!------------------------------------------------------------------------------------------

  subroutine namelist_to_buffer(namelist_filename, namelist_buffer)
  !DESCRIPTION
  !read in namelist
  !USES
  implicit none
  character(len=*)                         , intent(in)  :: namelist_filename
  character(len=ecosim_namelist_buffer_size) , intent(out) :: namelist_buffer

  character(len=*), parameter                            :: subname = 'namelist_to_buffer'
  character(len=ecosim_string_length_long)                 :: ioerror_msg
  integer :: nml_unit, nml_error,jj

  nml_unit = 16

  ! read the namelist file into a buffer.
  open(unit=nml_unit, file=trim(namelist_filename), action='read', access='stream', &
       form='unformatted', iostat=nml_error)
  if (nml_error == 0) then
     read(unit=nml_unit, iostat=nml_error, iomsg=ioerror_msg) namelist_buffer

     ! we should always reach the EOF to capture the entire file...
     if (nml_error/=-1) then
        write(stdout, '(a, a, i8)') subname, &
             ": IO ERROR reading namelist file into buffer: ", nml_error
        write(stdout, '(a)') ioerror_msg
        call abort()
     else
        write(stdout, '(a, a, a)') "Read '", trim(namelist_filename), "' until EOF."
     end if
     !remove junk
     do jj=ecosim_namelist_buffer_size,1,-1
       if(namelist_buffer(jj:jj)=='/')exit
       namelist_buffer(jj:jj)=''
     enddo
     write(stdout, '(a, a, i7, a)') subname, ": Read buffer of ", &
          len_trim(namelist_buffer), " characters."
     write(stdout, '(a)') "------------------------------"
     write(stdout, '(a)') "  If it looks like part of the namelist is missing, "
     write(stdout, '(a)') "  compare the number of characters read to the actual "
     write(stdout, '(a,a,a)') "  size of your file ($ wc -c ", trim(namelist_filename), ") and increase "
     write(stdout, '(a)') "  the buffer size if necessary."
     write(stdout, '(a)') "------------------------------"
     write(stdout, '(a)') trim(namelist_buffer)
     write(stdout, '(a)') "------------------------------"
     
  else
     write(stdout, '(a, a, i8, a, a)') subname, ": IO ERROR ", nml_error, &
          " opening namelist file : ", trim(namelist_filename)
     call abort()
  end if
  close(nml_unit)
  end subroutine namelist_to_buffer

!------------------------------------------------------------------------------------------
  subroutine getfilenamef(s1,s2,s3)
  implicit none
  character(len=*), intent(in) :: s1
  character(len=*), intent(in) :: s2
  character(len=*), intent(out) :: s3

  integer :: k,k1,l3,l2
  
  l3=len(s3)
  k=1
  do while(.true.)
    if(s1(k:k)/=' '.and. ichar(s1(k:k))/=0 .and. k<=l3)then
       s3(k:k)=s1(k:k)
    else
       exit
    endif
    k=k+1
  enddo
  if(k>l3)then
    call endrun("Not sufficient memory allocated for string s3 in " &
      //trim(mod_filename),__LINE__)
  endif
  k1=1;l2=len(s2)
  do while(k1<=l2)
    if(s2(k1:k1)/=' '.and. ichar(s2(k1:k1))/=0 .and. k<=l3)then
       s3(k:k)=s2(k1:k1)
    else
       exit
    endif
    k=k+1
    k1=k1+1
  enddo
  s3(k:)=' '
  end subroutine getfilenamef
  
!------------------------------------------------------------------------------------------

  logical function isprtablesymb(c)
  implicit none
  character(len=1), intent(in) :: c
  
  isprtablesymb=ICHAR(c) >= ICHAR(' ') .AND. ICHAR(c) <= ICHAR('~')
  end function isprtablesymb

  !------------------------------------------------------------------------
  subroutine getfil (fulpath, locfn, iflag)
    !
    ! !DESCRIPTION:
    ! Obtain local copy of file
    ! First check current working directory
    ! Next check full pathname[fulpath] on disk
    ! 
    ! !USES:
    implicit none
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)  :: fulpath !Archival or permanent disk full pathname
    character(len=*), intent(out) :: locfn   !output local file name
    integer,          intent(in)  :: iflag   !0=>abort if file not found 1=>do not abort
    !
    ! !LOCAL VARIABLES:
    integer :: klen            !length of fulpath character string
    logical :: lexist          !true if local file exists
    !------------------------------------------------------------------------

    ! get local file name from full name

    locfn = get_filename( fulpath )
    if (len_trim(locfn) == 0) then
        write(iulog,*)'(GETFIL): local filename has zero length'
       call endrun()
    else
        write(iulog,*)'(GETFIL): attempting to find local file ',  &
            trim(locfn)
    endif

    ! first check if file is in current working directory.

    inquire (file=locfn,exist=lexist)
    if (lexist) then
       write(iulog,*) '(GETFIL): using ',trim(locfn), &
            ' in current working directory'
       RETURN
    endif

    ! second check for full pathname on disk
    locfn = fulpath

    inquire (file=fulpath,exist=lexist)
    if (lexist) then
        write(iulog,*) '(GETFIL): using ',trim(fulpath)
       RETURN
    else
       write(iulog,*)'(GETFIL): failed getting file from full path: ', fulpath
       if (iflag==0) then
          call endrun('GETFIL: FAILED to get '//trim(fulpath)//' in ' &
      //mod_filename,line=__LINE__)
       else
          RETURN
       endif
    endif

  end subroutine getfil  

  !-----------------------------------------------------------------------
  character(len=512) function get_filename (fulpath)
    !
    ! !DESCRIPTION:
    ! Returns filename given full pathname
    !
    implicit none
    ! !ARGUMENTS:
    character(len=*), intent(in)  :: fulpath !full pathname
    !
    ! !LOCAL VARIABLES:
    integer :: i               !loop index
    integer :: klen            !length of fulpath character string
    !------------------------------------------------------------------------

    klen = len_trim(fulpath)
    do i = klen, 1, -1
       if (fulpath(i:i) == '/') go to 10
    end do
    i = 0
10  get_filename = fulpath(i+1:klen)

    return
  end function get_filename  

  !-----------------------------------------------------------------------
  subroutine strip_null(str)
  implicit none
  character(len=*), intent(inout) :: str
  integer :: i	
  do i=1,len(str)
      if(ichar(str(i:i))==0) str(i:i)=' '
  end do
  end subroutine strip_null  

  !------------------------------------------------------------------------
  integer function getavu()
    !
    ! !DESCRIPTION:
    ! Get next available Fortran unit number.
    !
    ! !USES:
    !------------------------------------------------------------------------
  implicit none

  getavu = file_getunit()

  end function getavu
  !------------------------------------------------------------------------
  subroutine relavu (iunit)
  implicit none
    !
    ! !DESCRIPTION:
    ! Close and release Fortran unit no longer in use!
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: iunit    !Fortran unit number
    !------------------------------------------------------------------------

    close(iunit)
    call file_freeUnit(iunit)

  end subroutine relavu  
  !------------------------------------------------------------------------
INTEGER FUNCTION file_getUnit ( unit )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(DAT_KIND_IN),intent(in),optional :: unit ! desired unit number

!EOP

   !----- local -----
   integer(DAT_KIND_IN)   :: n      ! loop index
   logical                :: opened ! If unit opened or not

   !----- formats -----
   character(*),parameter :: subName = '(file_getUnit) '
   character(*),parameter :: F00   = "('(file_getUnit) ',A,I4,A)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if (present (unit)) then
      inquire( unit, opened=opened )
      if (unit < 0 .or. unit > file_maxUnit) then
         write(iulog,F00) 'invalid unit number request:', unit
         call endrun( 'ERROR: bad input unit number '//trim(mod_filename),__LINE__)
      else if (opened .or. UnitTag(unit) .or. unit == 0 .or. unit == 5 &
               .or. unit == 6) then
         write(iulog,F00) 'unit number ', unit, ' is already in use'
         call endrun( 'ERROR: Input unit number already in use '//trim(mod_filename),__LINE__)
      else
         file_getUnit = unit
         UnitTag (unit)   = .true.
         return
      end if

   else
      ! --- Choose first available unit other than 0, 5, or 6  ------
      do n=file_maxUnit, file_minUnit, -1
         inquire( n, opened=opened )
         if (n == 5 .or. n == 6 .or. opened) then
            cycle
         end if
         if ( .not. UnitTag(n) ) then
            file_getUnit = n
            UnitTag(n)       = .true.
            return
         end if
      end do
   end if

   call endrun( subName//': Error: no available units found '//trim(mod_filename),__LINE__ )

END FUNCTION file_getUnit


!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: file_freeUnit -- Free up a FORTRAN unit number
!
! !DESCRIPTION: Free up the given unit number
!
! !REVISION HISTORY:
!     2005-Dec-14 - E. Kluzek - creation
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE file_freeUnit ( unit)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(DAT_KIND_IN),intent(in) :: unit  ! unit number to be freed

!EOP

   !----- local -----

   !----- formats -----
   character(*), parameter :: subName = '(file_freeUnit) '
   character(*), parameter :: F00 =   "('(file_freeUnit) ',A,I4,A)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if (unit < 0 .or. unit > file_maxUnit) then
      write(iulog,F00) 'invalid unit number request:', unit
   else if (unit == 0 .or. unit == 5 .or. unit == 6) then
      call endrun( subName//': Error: units 0, 5, and 6 must not be freed '//trim(mod_filename),__LINE__ )
   else if (UnitTag(unit)) then
      UnitTag (unit) = .false.
   else
      write(iulog,F00) 'unit ', unit, ' was not in use'
   end if

   return

END SUBROUTINE file_freeUnit

  !------------------------------------------------------------------------
  subroutine opnfil (locfn, iun, form)
    !
    ! !DESCRIPTION:
    ! Open file locfn in unformatted or formatted form on unit iun
    !
    ! !ARGUMENTS:
    character(len=*), intent(in):: locfn  !file name
    integer, intent(in):: iun             !fortran unit number
    character(len=1), intent(in):: form   !file format: u = unformatted, f = formatted
    !
    ! !LOCAL VARIABLES:
    integer ioe             !error return from fortran open
    character(len=11) ft    !format type: formatted. unformatted
    !------------------------------------------------------------------------

    if (len_trim(locfn) == 0) then
       call endrun('(OPNFIL): local filename has zero length in '//trim(mod_filename),__LINE__)
    endif
    if (form=='u' .or. form=='U') then
       ft = 'unformatted'
    else
       ft = 'formatted  '
    end if
    open (unit=iun,file=locfn,status='unknown',form=ft,iostat=ioe)
    if (ioe /= 0) then
       write(iulog,*)'(OPNFIL): failed to open file ',trim(locfn),        &
            &     ' on unit ',iun,' ierr=',ioe
       call endrun('Stopped in '//trim(mod_filename), __LINE__)
    else !if ( masterproc )then
       write(iulog,*)'(OPNFIL): Successfully opened file ',trim(locfn),   &
            &     ' on unit= ',iun
    end if

  end subroutine opnfil

  !------------------------------------------------------------------------
  subroutine print_ichar(instr)
  implicit none
  character(len=*), intent(in) ::instr 
  integer :: j
  do j=1,len(instr)
    print*,ichar(instr(j:j))
  enddo
  end subroutine print_ichar
  !------------------------------------------------------------------------
  subroutine strip_space(instr)
  implicit none
  character(len=*), intent(inout) ::instr 
  character(len=len(instr)) :: tstr
  integer :: j,k
  tstr=' '
  k=1
  do j=1,len(instr)
    if (ichar(instr(j:j))/=32)then
      tstr(k:k)=instr(j:j)
      k=k+1
    endif
  enddo
  instr=tstr
  end subroutine strip_space

  !------------------------------------------------------------------------
  function int2str(num)
  implicit none
  integer, intent(in) :: num
  character(len=256) :: int2str

  write(int2str,'(I0)')num
  end function int2str
end module fileUtil
