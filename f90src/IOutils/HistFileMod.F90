module HistFileMod
  ! Description:
  ! Register, accumulate and persist model histories and their restart buffers.

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use ncdio_pio
  use netcdf
  use fileUtil          , only : strip_null
  use abortutils        , only : endrun,iulog
  use TestMod           , only : errMsg
  use GridConsts        , only : JZ,JS,MaxNumBranches,bounds,bounds_type,NumOfPlantMorphUnits
  use GridConsts        , only : NumCanopyLayers,JP,NumGrowthStages,MaxNodesPerBranch,MaxNumRootAxes  
  use ElmIDMod          , only : NumPlantChemElms  
  use data_const_mod    , only : spval => DAT_CONST_SPVAL
  use EcosimConst       , only : secspday
  use EcoSIMCtrlMod     , only : etimer,lverb  
  use EcoSIMConfig      , only : case_name,hostname,version,source,username
  use DebugToolMod
implicit none
  save
  private
  character(len=*),  parameter :: mod_filename =&
  __FILE__
  integer, parameter :: max_chars = 256
  integer :: nfmaster = 0                        ! number of fields in master field list
  integer , public, parameter :: max_tapes = 6          ! max number of history tapes
  integer , public, parameter :: max_flds = 2500        ! max number of history fields
  integer , public, parameter :: max_namlen = 64        ! maximum number of characters for field name
  integer , parameter :: hist_dim_name_length = 16 ! lenngth of character strings in dimension names

  character(len=16), parameter :: nameg  = 'gridcell'     ! name of gridcells
  character(len=16), parameter :: namet  = 'topounit'     ! name of topographic units
  character(len=16), parameter :: namec  = 'column'       ! name of columns
  character(len=16), parameter :: namep  = 'pft'          ! name of patches

  integer :: ntapes = 0         ! index of max history file requested

  ! Pointers into datatype  arrays
  integer, parameter :: max_mapflds = 2500     ! Maximum number of fields to track

  integer :: ni                          ! implicit index below

  logical, public :: &
       hist_empty_htapes  = .false.      ! namelist: flag indicates no default history fields

  integer, public :: &
       hist_ndens(max_tapes) = 2         ! namelist: output density of netcdf history files

  integer, public :: &
       hist_mfilt(max_tapes) = (/ 1, (30, ni=2, max_tapes)/)        ! namelist: number of time samples per tape

  integer, public :: &
       hist_nhtfrq(max_tapes) = (/0, (-24, ni=2,max_tapes)/)        ! namelist: history write freq(0=monthly)

  character(len=1), public :: &
       hist_avgflag_pertape(max_tapes) = (/(' ',ni=1,max_tapes)/)   ! namelist: per tape averaging flag

  character(len=max_namlen), public :: &
       hist_type1d_pertape(max_tapes)  = (/(' ',ni=1,max_tapes)/)   ! namelist: per tape type1d

  character(len=max_namlen+2), public :: &
       fexcl(max_flds,max_tapes)         ! namelist-equivalence list of fields to remove

  character(len=max_namlen+2), public :: &
       fincl(max_flds,max_tapes)         ! namelist-equivalence list of fields to add

  character(len=max_namlen+2), public :: &
       hist_fincl1(max_flds) = ' '       ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl2(max_flds) = ' '       ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl3(max_flds) = ' '       ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl4(max_flds) = ' '       ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl5(max_flds) = ' '       ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl6(max_flds) = ' '       ! namelist: list of fields to add

  character(len=max_namlen+2), public :: &
       hist_fexcl1(max_flds) = ' ' ! namelist: list of fields to remove
  character(len=max_namlen+2), public :: &
       hist_fexcl2(max_flds) = ' ' ! namelist: list of fields to remove
  character(len=max_namlen+2), public :: &
       hist_fexcl3(max_flds) = ' ' ! namelist: list of fields to remove
  character(len=max_namlen+2), public :: &
       hist_fexcl4(max_flds) = ' ' ! namelist: list of fields to remove
  character(len=max_namlen+2), public :: &
       hist_fexcl5(max_flds) = ' ' ! namelist: list of fields to remove
  character(len=max_namlen+2), public :: &
       hist_fexcl6(max_flds) = ' ' ! namelist: list of fields to remove

  type ecosimpoint_rs                             ! Pointer to real scalar data (1D)
     real(r8), pointer :: ptr(:) => null()
  end type ecosimpoint_rs

  type ecosimpoint_ra1                             ! Pointer to real array data (2D)
     real(r8), pointer :: ptr(:,:) => null()
  end type ecosimpoint_ra1

  type ecosimpoint_ra2                             ! Pointer to real array data (2D)
     real(r8), pointer :: ptr(:,:,:) => null()
  end type ecosimpoint_ra2

  logical, private :: if_disphist(max_tapes)   ! restart, true => save history file

  type (ecosimpoint_rs), public :: esmptr_rs(max_mapflds) ! Real scalar data (1D)
  type (ecosimpoint_ra1), public :: esmptr_ra1(max_mapflds) ! Real array data (2D)

  type field_info
     character(len=max_namlen) :: name         ! field name
     character(len=max_chars)  :: long_name    ! long name
     character(len=max_chars)  :: standard_name  ! CF standard name
     character(len=max_chars)  :: units        ! units
     character(len=hist_dim_name_length) :: type1d                ! pointer to first dimension type from data type (nameg, etc)
     character(len=hist_dim_name_length) :: type1d_out            ! hbuf first dimension type from data type (nameg, etc)
     character(len=hist_dim_name_length) :: type2d                ! hbuf second dimension type ["levgrnd","levlak","numrad","ltype","natpft","cft","glc_nec","elevclas","subname(n)","month"]
     integer :: beg1d                          ! on-node 1d clm pointer start index
     integer :: end1d                          ! on-node 1d clm pointer end index
     integer :: num1d                          ! size of clm pointer first dimension (all nodes)
     integer :: numdims                        ! the actual number of dimensions, this allows
                                               ! for 2D arrays, where the second dimension is allowed
                                               ! to be 1
     integer :: beg1d_out                      ! on-node 1d hbuf pointer start index
     integer :: end1d_out                      ! on-node 1d hbuf pointer end index
     integer :: num1d_out                      ! size of hbuf first dimension (all nodes)
     integer :: num2d                          ! size of hbuf second dimension (e.g. number of vertical levels)
     integer :: hpindex                        ! history pointer index
     integer :: no_snow_behavior               ! for multi-layer snow fields, flag saying how to treat times when a given snow layer is absent
  end type field_info

  type history_entry
     type (field_info) :: field                ! field information
     character(len=1)  :: avgflag              ! time averaging flag
     real(r8), pointer :: hbuf(:,:)            ! history buffer (dimensions: dim1d x num2d)
     integer , pointer :: nacs(:,:)            ! accumulation counter (dimensions: dim1d x num2d)
  end type history_entry

  type history_tape
     integer  :: nflds                         ! number of active fields on tape
     integer  :: ntimes                        ! current number of time samples on tape
     integer  :: mfilt                         ! maximum number of time samples per tape
     integer  :: nhtfrq                        ! number of time samples per tape
     integer  :: ncprec                        ! netcdf output precision
     logical  :: is_endhist                    ! true => current time step is end of history interval
     real(r8) :: begtime                       ! time at beginning of history averaging interval
     type (history_entry) :: hlist(max_flds)   ! array of active history tape entries
  end type history_tape

  character(len=max_chars) :: locfnh(max_tapes)  ! local history file names
  ! Master list: an array of master_entry entities
  !
  type master_entry
     type (field_info)  :: field               ! field information
     logical            :: actflag(max_tapes)  ! active/inactive flag
     character(len=1)   :: avgflag(max_tapes)  ! time averaging flag ("X","A","M" or "I",)
  end type master_entry

  type (master_entry) :: masterlist(max_flds)  ! master field list

  type (history_tape), public :: tape(max_tapes)       ! array history tapes

  character(len=max_chars) :: locfnhr(max_tapes) ! local history restart file names
  logical :: htapes_defined = .false.            ! flag indicates history contents have been defined

! ids for netcdf files
  integer :: time_dimid                      ! time dimension id
  integer :: hist_interval_dimid             ! time bounds dimension id
  integer :: strlen_dimid                    ! string dimension id
  type(file_desc_t) :: nfid(max_tapes)       ! file ids
  type(file_desc_t) :: ncid_hist(max_tapes)  ! file ids for history restart files

  public :: hist_addfld1d        ! Add a 1d single-level field to the master field list
  public :: hist_addfld2d        ! Add a 2d multi-level field to the master field list
  public :: hist_printflds       ! Print summary of master field list
  public :: hist_htapes_build    ! Initialize history file handler for initial or continue run
  public :: hist_update_hbuf     ! Updates history buffer for all fields and tapes
  public :: hist_htapes_wrapup   ! Write history tape(s)
  public :: hist_restart_ncd     ! Read/write history file restart data
  public :: htapes_fieldlist     ! Define the contents of each history file based on namelist

  contains

  subroutine htape_create(t,histrest)

  !
  !DESCRIPTION
  !create history file for output

  use EcoSIMConfig      , only : case_name
  use EcoSIMConfig      , only : jcplx=>jcplxc,jsken=>jskenc
  implicit none
  integer, intent(in) :: t        ! tape index
  logical, optional, intent(in) :: histrest  !if creating the history restart file

  integer :: ncprec              ! output netCDF write precision
  logical :: lhistrest           ! local history restart flag
  type(file_desc_t) :: lnfid     ! local file id
  character(len=  8) :: curdate  ! current date
  character(len=  8) :: curtime  ! current time
  character(len=256) :: name     ! name of attribute
  character(len=256) :: units    ! units of attribute
  character(len=256) :: str      ! global attribute string
  character(len=  1) :: avgflag  ! time averaging flag
  integer :: numc
  integer :: nump
  integer :: numg,numl,numt
  integer :: dimid               ! dimension id temporary
  character(len=*),parameter :: subname = 'htape_create'
  integer :: ier

  if ( present(histrest) )then
    lhistrest = histrest
  else
    lhistrest = .false.
  end if

  call get_grid_info(ng=numg,nt=numt,nc=numc, np=nump)

  ncprec = tape(t)%ncprec

  if ( .not. lhistrest )then
    call ncd_pio_createfile(lnfid, trim(locfnh(t)))
    call check_ret(ncd_putatt(lnfid, ncd_global, 'title', 'EcoSIM History file information' ),&
      trim(subname)//' title1')
  else
    call ncd_pio_createfile(lnfid, trim(locfnhr(t)))
    call check_ret(ncd_putatt(lnfid, ncd_global, 'title', &
       'EcoSIM Restart History information, required to continue a simulation'),trim(subname)//' tilte2' )
    call check_ret(ncd_putatt(lnfid, ncd_global, 'comment', &
                   "This entire file NOT needed for startup or branch simulations"),trim(subname)//' comment1')
  endif

  call check_ret(ncd_putatt(lnfid, ncd_global, 'source'     , trim(source)),trim(subname)//' source')
  call check_ret(ncd_putatt(lnfid, ncd_global, 'source_id'  , trim(version)),trim(subname)//' source_id')
  call check_ret(ncd_putatt(lnfid, ncd_global, 'product'  , 'model-output'),trim(subname)//' product')
  call check_ret(ncd_putatt(lnfid, ncd_global, 'case', trim(case_name)),trim(subname)//'case')
  call check_ret(ncd_putatt(lnfid, ncd_global, 'username', trim(username)),trim(subname)//' username')
  call check_ret(ncd_putatt(lnfid, ncd_global, 'hostname', trim(hostname)),trim(subname)//' hostname')
  call check_ret(ncd_putatt(lnfid, ncd_global, 'git_version' , trim(version)),trim(subname)//' git_version')
  call getdatetime(curdate, curtime)

  ! Global compressed dimensions (not including non-land points)
  call ncd_defdim(lnfid, trim(nameg), numg, dimid)
  call ncd_defdim(lnfid, trim(namet), numt, dimid)
  call ncd_defdim(lnfid, trim(namec), numc, dimid)
  call ncd_defdim(lnfid, trim(namep), nump, dimid)

  ! "level" dimensions
  call ncd_defdim(lnfid, 'node', MaxNodesPerBranch, dimid)
  call ncd_defdim(lnfid, 'levsoi', JZ, dimid)
  call ncd_defdim(lnfid, 'levsno',  JS,dimid)
  call ncd_defdim(lnfid, 'levcan',NumCanopyLayers,dimid)
  call ncd_defdim(lnfid, 'npfts',  JP,dimid)
  call ncd_defdim(lnfid, 'nbranches',MaxNumBranches,dimid)
  call ncd_defdim(lnfid, 'ngrstages',NumGrowthStages,dimid)
  call ncd_defdim(lnfid, 'elements',NumPlantChemElms,dimid)
  call ncd_defdim(lnfid,'rootaxs',MaxNumRootAxes,dimid)
  call ncd_defdim(lnfid, 'nkinecomp',jsken,dimid)
  call ncd_defdim(lnfid, 'nomcomplx',jcplx,dimid)
  call ncd_defdim(lnfid, 'pmorphunits',NumOfPlantMorphUnits,dimid)
  if ( .not. lhistrest )then
    call ncd_defdim(lnfid, 'hist_interval', 2, hist_interval_dimid)
    call ncd_defdim(lnfid, 'time', ncd_unlimited, time_dimid)
    nfid(t) = lnfid
    write(iulog,*) trim(subname), &
                    ' : Successfully defined netcdf history file ',t
    flush(iulog)

  else
    ncid_hist(t) = lnfid
    write(iulog,*) trim(subname), &
                      ' : Successfully defined netcdf restart history file ',t
    flush(iulog)
  end if

  end subroutine htape_create

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: getdatetime
!
! !INTERFACE:
  subroutine getdatetime (cdate, ctime)
!
! !DESCRIPTION:
! A generic Date and Time routine
!
! !USES:
! !ARGUMENTS:
  implicit none
  character(len=8), intent(out) :: cdate  !current date
  character(len=8), intent(out) :: ctime  !current time
! !LOCAL VARIABLES:
!EOP
  character(len=8)      :: date       !current date
  character(len=10)     :: time       !current time
  character(len=5)      :: zone       !zone
  integer, dimension(8) :: values !temporary
  integer               :: ier    !MPI error code

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

  end subroutine getdatetime
!-----------------------------------------------------------------------
  subroutine get_grid_info(ng,nt,nc,np)
  implicit none
  integer, intent(out) :: ng,nt,nc,np

  nc=bounds%ncols
  np=bounds%npfts
  ng=bounds%ngrid
  nt=bounds%ntopou
  end subroutine get_grid_info

  !-----------------------------------------------------------------------
  subroutine hist_printflds()
    !
    ! !DESCRIPTION:
    ! Print summary of master field list.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
  implicit none
  integer :: nf
  character(len=*),parameter :: subname = 'hist_printflds'
    !-----------------------------------------------------------------------

!    if (masterproc) then
     write(iulog,*) trim(subname),' : number of master fields = ',nfmaster
     write(iulog,*)' ******* MASTER FIELD LIST *******'
     do nf = 1,nfmaster
       write(iulog,9000)nf, masterlist(nf)%field%name, masterlist(nf)%field%units
9000      format (i5,1x,a32,1x,a16)
     end do
     flush(iulog)
!    end if
  end subroutine hist_printflds
!-----------------------------------------------------------------------

  subroutine hist_addfld1d(fname,units,avgflag,long_name,type1d_out, standard_name,&
    ptr_gcell,ptr_topo,ptr_col,ptr_patch, default)
  ! Add a 1d single-level field to the master field list
  implicit none
  character(len=*), intent(in)           :: fname          ! field name
  character(len=*), intent(in)           :: units          ! units of field
  character(len=*), intent(in)           :: long_name      ! long name of field
  character(len=1), intent(in)           :: avgflag        ! time averaging flag
  character(len=*), optional, intent(in) :: type1d_out     ! output type (from data type)
  character(len=*), optional, intent(in) :: standard_name  ! CF standard name
  real(r8)        , optional, pointer    :: ptr_gcell(:)   ! pointer to gridcell array
  real(r8)        , optional, pointer    :: ptr_topo(:)    ! pointer to topounit array
  real(r8)        , optional, pointer    :: ptr_col(:)     ! pointer to column array
  real(r8)        , optional, pointer    :: ptr_patch(:)   ! pointer to pft array
  character(len=*), optional, intent(in) :: default        ! if set to 'inactive, field will not appear on primary tape

  character(len=*),parameter :: subname = 'hist_addfld1d'
  character(len=hist_dim_name_length) :: l_type1d       ! 1d data type
  character(len=hist_dim_name_length) :: l_type1d_out   ! 1d output type
  integer :: hpindex                 ! history buffer pointer index
  character(len=16):: l_default      ! local version of 'default'
  character(len=max_namlen) :: lstandard_name  ! local standard name

  hpindex = pointer_index()

  if (present(ptr_gcell)) then
    l_type1d = nameg
    l_type1d_out = nameg
    esmptr_rs(hpindex)%ptr => ptr_gcell
  else if (present(ptr_topo)) then
    l_type1d = namet
    l_type1d_out = namet
    esmptr_rs(hpindex)%ptr => ptr_topo
  else if (present(ptr_col)) then
    l_type1d = namec
    l_type1d_out = namec
    esmptr_rs(hpindex)%ptr => ptr_col
  else if (present(ptr_patch)) then
    l_type1d = namep
    l_type1d_out = namep
    esmptr_rs(hpindex)%ptr => ptr_patch
  else
    write(iulog,*) trim(subname),' ERROR: must specify a valid pointer index,', &
          ' choices are [ptr_atm, ptr_lnd, ptr_gcell, ptr_topo, ptr_lunit, ptr_col, ptr_patch] '
    call endrun(msg=errMsg(__FILE__, __LINE__))
  end if

  if (present(standard_name)) then
    lstandard_name = standard_name
  else
    lstandard_name = ' '
  end if

  call masterlist_addfld (fname=trim(fname), type1d=l_type1d, type1d_out=l_type1d_out, &
    type2d='unset', numdims=1, num2d=1, units=units, avgflag=avgflag, long_name=long_name,&
    standard_name=lstandard_name, hpindex= hpindex)

  l_default = 'active'
  if (present(default)) then
    l_default = default
  end if

  if (trim(l_default) == 'inactive') then
    return
  else
    call masterlist_make_active (name=trim(fname), tape_index=1)
  end if

  end subroutine hist_addfld1d

!-----------------------------------------------------------------------

  subroutine masterlist_addfld (fname, type1d, type1d_out,&
        type2d, numdims, num2d, units, avgflag, long_name, standard_name, hpindex)


  use GridConsts, only : bounds
  implicit none
  character(len=*), intent(in)  :: fname          ! field name
  character(len=*), intent(in)  :: type1d           ! 1d data type
  character(len=*), intent(in)  :: type1d_out       ! 1d output type
  character(len=*), intent(in)  :: type2d           ! 2d output type
  integer         , intent(in)  :: numdims          ! number of dimensions
  integer         , intent(in)  :: num2d            ! size of second dimension (e.g. number of vertical levels)
  character(len=*), intent(in)  :: units            ! units of field
  character(len=1), intent(in)  :: avgflag          ! time averaging flag
  character(len=*), intent(in)  :: long_name        ! long name of field
  character(len=*), intent(in)  :: standard_name        ! long name of field
  integer         , intent(in)  :: hpindex          ! data type index for history buffer output

  character(len=*), parameter :: subname='masterlist_addfld'
  integer :: n,f
  integer :: numg         ! total number of gridcells across all processors
  integer :: numt         ! total number of topounits across all processors
  integer :: numc         ! total number of columns across all processors
  integer :: nump         ! total number of pfts across all processors

  numg=bounds%ngrid;numt=bounds%ntopou;numc=bounds%ncols;nump=bounds%npfts
  if (fname == ' ') then
    write(iulog,*) trim(subname),' ERROR: blank field name not allowed'
    call endrun(msg=errMsg(__FILE__, __LINE__))
  end if

  do n = 1,nfmaster
    if (masterlist(n)%field%name == fname) then
      write(iulog,*) trim(subname),' ERROR:', fname, ' already on list'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
  end do

  ! Increase number of fields on master field list
  nfmaster = nfmaster + 1
  f = nfmaster

    ! Check number of fields in master list against maximum number for master list

  if (nfmaster > max_flds) then
    write(iulog,*) trim(subname),' ERROR: too many fields for primary history file ', &
          '-- max_flds,nfmaster=', max_flds, nfmaster
    call endrun(msg=errMsg(__FILE__, __LINE__))
  end if

  masterlist(f)%field%name           = fname
  masterlist(f)%field%long_name      = long_name
  masterlist(f)%field%standard_name  = standard_name
  masterlist(f)%field%units          = units
  masterlist(f)%field%type1d         = type1d
  masterlist(f)%field%type1d_out     = type1d_out
  masterlist(f)%field%type2d         = type2d
  masterlist(f)%field%numdims        = numdims
  masterlist(f)%field%num2d          = num2d
  masterlist(f)%field%hpindex        = hpindex

  select case (type1d)
  case (nameg)
      masterlist(f)%field%beg1d = bounds%begg
      masterlist(f)%field%end1d = bounds%endg
      masterlist(f)%field%num1d = numg
  case (namet)
      masterlist(f)%field%beg1d = bounds%begt
      masterlist(f)%field%end1d = bounds%endt
      masterlist(f)%field%num1d = numt
  case (namec)
      masterlist(f)%field%beg1d = bounds%begc
      masterlist(f)%field%end1d = bounds%endc
      masterlist(f)%field%num1d = numc
  case (namep)
      masterlist(f)%field%beg1d = bounds%begp
      masterlist(f)%field%end1d = bounds%endp
      masterlist(f)%field%num1d = nump
  case default
      write(iulog,*) trim(subname),' ERROR: unknown 1d output type= ',type1d
      call endrun(msg=errMsg(__FILE__, __LINE__))
  end select
  masterlist(f)%avgflag(:) = avgflag
  masterlist(f)%actflag(:) = .false.

  end subroutine masterlist_addfld
  !-----------------------------------------------------------------------
  integer function pointer_index ()
    !
    ! !DESCRIPTION:
    ! Set the current pointer index and increment the value of the index.
    !
    ! !ARGUMENTS:
    !
    integer, save :: lastindex = 1
    character(len=*),parameter :: subname = 'pointer_index'
    !-----------------------------------------------------------------------

    pointer_index = lastindex
    lastindex = lastindex + 1
    if (lastindex > max_mapflds) then
       write(iulog,*) trim(subname),' ERROR: ',&
            ' lastindex = ',lastindex,' greater than max_mapflds= ',max_mapflds
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end function pointer_index
!-----------------------------------------------------------------------

  subroutine hist_addfld2d(fname, units, type2d, avgflag, long_name, type1d_out, standard_name,&
    ptr_gcell,ptr_topo,ptr_col,ptr_patch, default)
  ! Add a 1d single-level field to the master field list
  implicit none
  character(len=*), intent(in)           :: fname            ! field name
  character(len=*), intent(in)           :: units            ! units of field
  character(len=*), intent(in)           :: type2d           ! 2d output type
  character(len=1), intent(in)           :: avgflag          ! time averaging flag
  character(len=*), intent(in)           :: long_name        ! long name of field
  character(len=*), optional, intent(in) :: type1d_out       ! output type (from data type)
  character(len=*), optional, intent(in) :: standard_name    ! output type (from data type)
  real(r8)        , optional, pointer    :: ptr_gcell(:,:)   ! pointer to gridcell array
  real(r8)        , optional, pointer    :: ptr_topo(:,:)    ! pointer to topounit array
  real(r8)        , optional, pointer    :: ptr_col(:,:)     ! pointer to column array
  real(r8)        , optional, pointer    :: ptr_patch(:,:)   ! pointer to pft array
  character(len=*), optional, intent(in) :: default        ! if set to 'inactive, field will not appear on primary tape
  character(len=*),parameter :: subname = 'hist_addfld2d'

  integer :: num2d                   ! size of second dimension (e.g. number of vertical levels)
  character(len=hist_dim_name_length) :: l_type1d       ! 1d data type
  character(len=hist_dim_name_length) :: l_type1d_out   ! 1d output type
  character(len=max_namlen) :: lstandard_name  ! local standard name
  integer :: hpindex                 ! history buffer pointer index
  character(len=16):: l_default      ! local version of 'default'

  hpindex = pointer_index()
  
  select case (type2d)
  case ('levsoi')
      num2d = JZ
  case ('nbranches')
      num2d = MaxNumBranches
  case ('rootaxs')    
      num2d=MaxNumRootAxes
  case ('pmorphunits')
      num2d=NumOfPlantMorphUnits    
  case ('elements')    
      num2d=NumPlantChemElms
  case ('levcan')
      num2d=NumCanopyLayers
  case ('node')
      num2d=MaxNodesPerBranch    
  case default
      write(iulog,*) trim(subname),' ERROR: unsupported 2d type ',type2d, &
        ' currently supported types for multi level fields are: ', &
        '[levgrnd,levlak,numrad,nmonthlevdcmp,levtrc,ltype,natpft,cft,'&
          //'glc_nec,elevclas,levsno,levsoi,nbranches,pmorphunits]'
      call endrun(msg=errMsg(__FILE__, __LINE__))
  end select
  if (present(ptr_gcell)) then
    l_type1d = nameg
    l_type1d_out = nameg
    esmptr_ra1(hpindex)%ptr => ptr_gcell
  else if (present(ptr_topo)) then
    l_type1d = namet
    l_type1d_out = namet
    esmptr_ra1(hpindex)%ptr => ptr_topo
  else if (present(ptr_col)) then
    l_type1d = namec
    l_type1d_out = namec
    esmptr_ra1(hpindex)%ptr => ptr_col
  else if (present(ptr_patch)) then
    l_type1d = namep
    l_type1d_out = namep
    print*,'fname',fname,hpindex,size(ptr_patch,2),hpindex
    esmptr_ra1(hpindex)%ptr => ptr_patch
  else
    write(iulog,*) trim(subname),' ERROR: must specify a valid pointer index,', &
          ' choices are [ptr_atm, ptr_lnd, ptr_gcell, ptr_topo, ptr_lunit, ptr_col, ptr_patch] '
    call endrun(msg=errMsg(__FILE__, __LINE__))
  end if

  if (present(standard_name)) then
    lstandard_name = standard_name
  else
    lstandard_name = ' '
  end if

  call masterlist_addfld (fname=trim(fname), type1d=l_type1d, type1d_out=l_type1d_out, &
    type2d=type2d, numdims=2, num2d=num2d, units=units, avgflag=avgflag, long_name=long_name, &
    standard_name=lstandard_name, hpindex=hpindex)

  l_default = 'active'
  if (present(default)) then
    l_default = default
  end if

  if (trim(l_default) == 'inactive') then
    return
  else
    call masterlist_make_active (name=trim(fname), tape_index=1)
  end if

  end subroutine hist_addfld2d

  !-----------------------------------------------------------------------
  subroutine masterlist_make_active (name, tape_index, avgflag)
    !
    ! !DESCRIPTION:
    ! Add a field to the default ``on'' list for a given history file.
    ! Also change the default time averaging flag if requested.
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: name          ! field name
    integer, intent(in) :: tape_index             ! history tape index
    character(len=1), intent(in), optional :: avgflag  ! time averaging flag
    !
    ! !LOCAL VARIABLES:
    integer :: f            ! field index
    logical :: found        ! flag indicates field found in masterlist
    character(len=*),parameter :: subname = 'masterlist_make_active'

    ! Check validity of input arguments

    if (tape_index > max_tapes) then
       write(iulog,*) trim(subname),' ERROR: tape index=', tape_index, ' is too big'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (present(avgflag)) then
       if ( avgflag /= ' ' .and. &   !is not defined
            avgflag /= 'A' .and. &   !is not temporal average
            avgflag /= 'P' .and. &   !is not temporal average over positive values           
            avgflag /= 'I' .and. &   !is not instantaneous
            avgflag /= 'X' .and. &   !is not the maximum over the time period
            avgflag /= 'M') then     !is not the minimum over the time period
          write(iulog,*) trim(subname),' ERROR: unknown averaging flag=', avgflag
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
    end if

    ! Look through master list for input field name.
    ! When found, set active flag for that tape to true.
    ! Also reset averaging flag if told to use other than default.

    found = .false.
    do f = 1,nfmaster
       if (trim(name) == trim(masterlist(f)%field%name)) then
          masterlist(f)%actflag(tape_index) = .true.
          if (present(avgflag)) then
             if (avgflag/= ' ') masterlist(f)%avgflag(tape_index) = avgflag
          end if
          found = .true.
          exit
       end if
    end do
    if (.not. found) then
       write(iulog,*) trim(subname),' ERROR: field=', name, ' not found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine masterlist_make_active


  !-----------------------------------------------------------------------
  subroutine hist_htapes_build ()
    !
    ! !DESCRIPTION:
    ! Initialize history file for initial or continuation run.  For example,
    ! on an initial run, this routine initializes ``ntapes'' history files.
    ! On a restart run, this routine only initializes history files declared
    ! beyond what existed on the previous run.  Files which already existed on
    ! the previous run have already been initialized (i.e. named and opened)
    ! in routine restart\_history.  Loop over tapes and fields per tape setting
    ! appropriate variables and calling appropriate routines
    !
    ! !USES:
    !
    implicit none
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: i                   ! index
    integer :: ier                 ! error code
    integer :: t, f                ! tape, field indices
    integer :: day, sec            ! day and seconds from base date
    character(len=*),parameter :: subname = trim(mod_filename)//'::hist_htapes_build'

    !-----------------------------------------------------------------------

!    if (masterproc) then
       write(iulog,*)  trim(subname),' Initializing ecosim history files'
       write(iulog,'(72a1)') ("-",i=1,60)
       call flush(iulog)
!    endif

    ! Define field list information for all history files.
    ! Update ntapes to reflect number of active history files
    ! Note - branch runs can have additional auxiliary history files
    ! declared).

    call htapes_fieldlist()

    ! Set number of time samples in each history file and
    ! Note - the following entries will be overwritten by history restart
    ! Note - with netcdf, only 1 (ncd_double) and 2 (ncd_float) are allowed

    do t=1,ntapes
       tape(t)%ntimes = 0
       tape(t)%nhtfrq = hist_nhtfrq(t)
       tape(t)%mfilt = hist_mfilt(t)
       if (hist_ndens(t) == 1) then
          tape(t)%ncprec = ncd_double
       else
          tape(t)%ncprec = ncd_float
       endif
    end do

    ! Set time of beginning of current averaging interval
    ! First etermine elapased time since reference date

    call etimer%get_prev_time(day, sec)
    do t=1,ntapes
       tape(t)%begtime = day + sec/secspday
    end do

!    if (masterproc) then
       write(iulog,*)  trim(subname),' Successfully initialized elm history files'
       write(iulog,'(72a1)') ("-",i=1,60)
       call flush(iulog)
 !   endif

  end subroutine hist_htapes_build

  !-----------------------------------------------------------------------
  subroutine htapes_fieldlist()
    !
    ! !DESCRIPTION:
    ! Define the contents of each history file based on namelist
    ! input for initial or branch run, and restart data if a restart run.
    ! Use arrays fincl and fexcl to modify default history tape contents.
    ! Then sort the result alphanumerically.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: t, f                         ! tape, field indices
    integer :: ff                           ! index into include, exclude and fprec list
    character(len=max_namlen) :: name       ! field name portion of fincl (i.e. no avgflag separator)
    character(len=max_namlen) :: mastername ! name from masterlist field
    character(len=1)  :: avgflag            ! averaging flag
    character(len=1)  :: prec_acc           ! history buffer precision flag
    character(len=1)  :: prec_wrt           ! history buffer write precision flag
    type (history_entry) :: tmp             ! temporary used for swapping
    character(len=*),parameter :: subname = 'htapes_fieldlist'
    !-----------------------------------------------------------------------

    ! Override averaging flag for all fields on a particular tape
    ! if namelist input so specifies

    do t=1,max_tapes
       if (hist_avgflag_pertape(t) /= ' ') then
          call masterlist_change_timeavg (t)
       end if
    end do

    fincl(:,1) = hist_fincl1(:)
    fincl(:,2) = hist_fincl2(:)
    fincl(:,3) = hist_fincl3(:)
    fincl(:,4) = hist_fincl4(:)
    fincl(:,5) = hist_fincl5(:)
    fincl(:,6) = hist_fincl6(:)

    fexcl(:,1) = hist_fexcl1(:)
    fexcl(:,2) = hist_fexcl2(:)
    fexcl(:,3) = hist_fexcl3(:)
    fexcl(:,4) = hist_fexcl4(:)
    fexcl(:,5) = hist_fexcl5(:)
    fexcl(:,6) = hist_fexcl6(:)


    ! First ensure contents of fincl and fexcl are valid names

    do t = 1,max_tapes
       f = 1
       do while (f < max_flds .and. fincl(f,t) /= ' ')
          name = getname (fincl(f,t))
          do ff = 1,nfmaster
             mastername = masterlist(ff)%field%name
             if (name == mastername) exit             
          end do

          if (name /= mastername) then
             write(iulog,*) trim(subname),' ERROR: ', trim(name), ' in fincl(', f, ') ',&
                  'for history tape ',t,' not found'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if
          f = f + 1
       end do
       
       f = 1
       do while (f < max_flds .and. fexcl(f,t) /= ' ')
          do ff = 1,nfmaster
             mastername = masterlist(ff)%field%name
             if (fexcl(f,t) == mastername) exit
          end do
          if (fexcl(f,t) /= mastername) then
             write(iulog,*) trim(subname),' ERROR: ', fexcl(f,t), ' in fexcl(', f, ') ', &
                  'for history tape ',t,' not found'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if
          f = f + 1
       end do
    end do

    tape(:)%nflds = 0
    do t = 1,max_tapes

       ! Loop through the masterlist set of field names and determine if any of those
       ! are in the FINCL or FEXCL arrays
       ! The call to list_index determines the index in the FINCL or FEXCL arrays
       ! that the masterlist field corresponds to
       ! Add the field to the tape if specified via namelist (FINCL[1-max_tapes]),
       ! or if it is on by default and was not excluded via namelist (FEXCL[1-max_tapes]).

       do f = 1,nfmaster
          mastername = masterlist(f)%field%name
          call list_index (fincl(1,t), mastername, ff)

          if (ff > 0) then

             ! if field is in include list, ff > 0 and htape_addfld
             ! will not be called for field

             avgflag = getflag (fincl(ff,t))
             call htape_addfld (t, f, avgflag)

          else if (.not. hist_empty_htapes) then

             ! find index of field in exclude list

             call list_index (fexcl(1,t), mastername, ff)

             ! if field is in exclude list, ff > 0 and htape_addfld
             ! will not be called for field
             ! if field is not in exclude list, ff =0 and htape_addfld
             ! will be called for field (note that htape_addfld will be
             ! called below only if field is not in exclude list OR in
             ! include list

             if (ff == 0 .and. masterlist(f)%actflag(t)) then
                call htape_addfld (t, f, ' ')
             end if

          end if
       end do

       ! Specification of tape contents now CO2CompenPoint_nodeete.
       ! Sort each list of active entries

       do f = tape(t)%nflds-1,1,-1
          do ff = 1,f
             if (tape(t)%hlist(ff)%field%name > tape(t)%hlist(ff+1)%field%name) then

                tmp = tape(t)%hlist(ff)
                tape(t)%hlist(ff  ) = tape(t)%hlist(ff+1)
                tape(t)%hlist(ff+1) = tmp

             else if (tape(t)%hlist(ff)%field%name == tape(t)%hlist(ff+1)%field%name) then

                write(iulog,*) trim(subname),' ERROR: Duplicate field ', &
                   tape(t)%hlist(ff)%field%name, &
                   't,ff,name=',t,ff,tape(t)%hlist(ff+1)%field%name
                call endrun(msg=errMsg(__FILE__, __LINE__))

             end if
          end do
       end do

!       if (masterproc) then
          if (tape(t)%nflds > 0) then
             write(iulog,*) trim(subname),' : Included fields tape ',t,'=',tape(t)%nflds
          end if
          do f = 1,tape(t)%nflds
             write(iulog,*) f,' ',tape(t)%hlist(f)%field%name, &
                  tape(t)%hlist(f)%field%num2d,' ',tape(t)%hlist(f)%avgflag
          end do
          call flush(iulog)
!       end if
    end do

    ! Determine total number of active history tapes

    ntapes = 0
    do t = max_tapes,1,-1
       if (tape(t)%nflds > 0) then
          ntapes = t
          exit
       end if
    end do

    ! Ensure there are no "holes" in tape specification, i.e. empty tapes.
    ! Enabling holes should not be difficult if necessary.

    do t = 1,ntapes
       if (tape(t)%nflds  ==  0) then
          write(iulog,*) trim(subname),' ERROR: Tape ',t,' is empty'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    end do

    ! Check that the number of history files declared does not exceed
    ! the maximum allowed.

    if (ntapes > max_tapes) then
       write(iulog,*) trim(subname),' ERROR: Too many history files declared, max_tapes=',max_tapes
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Change 1d output per tape output flag if requested - only for history
    ! tapes where 2d xy averaging is not enabled

    do t = 1,ntapes
       if (hist_type1d_pertape(t) /= ' ') then
          select case (trim(hist_type1d_pertape(t)))
          case ('PFTS','COLS', 'TOPO', 'GRID')
!             if ( masterproc ) &
             write(iulog,*)'history tape ',t,' will have 1d output type of ',hist_type1d_pertape(t)
          case default
             write(iulog,*) trim(subname),' ERROR: unknown namelist type1d per tape=',hist_type1d_pertape(t)
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
       end if
    end do

!    if (masterproc) then
       write(iulog,*) 'There will be a total of ',ntapes,' history tapes'
       do t=1,ntapes
          write(iulog,*)
          if (hist_nhtfrq(t) == 0) then
             write(iulog,*)'History tape ',t,' write frequency is MONTHLY'
          else
             write(iulog,*)'History tape ',t,' write frequency = ',hist_nhtfrq(t)
          endif
          write(iulog,*)'Number of time samples on history tape ',t,' is ',hist_mfilt(t)
          write(iulog,*)'Output precision on history tape ',t,'=',hist_ndens(t)
          write(iulog,*)
       end do
       call flush(iulog)
 !   end if

    ! Set flag indicating h-tape contents are now defined (needed by masterlist_addfld)

    htapes_defined = .true.

  end subroutine htapes_fieldlist


  !-----------------------------------------------------------------------
  subroutine list_index (list, name, index)
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: list(max_flds)  ! input list of names, possibly ":" delimited
    character(len=max_namlen), intent(in) :: name   ! name to be searched for
    integer, intent(out) :: index                   ! index of "name" in "list"
    !
    ! !LOCAL VARIABLES:
    !EOP
    character(len=max_namlen) :: listname           ! input name with ":" stripped off.
    integer :: f                                    ! field index
    character(len=*),parameter :: subname = 'list_index'
    !-----------------------------------------------------------------------

    ! Only list items

    index = 0
    do f=1,max_flds
      listname = getname (list(f))
      if (listname == ' ') exit
      if (listname == name) then
          index = f
          exit
      end if
    end do

  end subroutine list_index

  !-----------------------------------------------------------------------
  character(len=max_namlen) function getname (inname)
    !
  ! !DESCRIPTION:
  ! Retrieve name portion of inname. If an averaging flag separater character
  ! is present (:) in inname, lop it off.
  !
  ! !ARGUMENTS:
  character(len=*), intent(in) :: inname
  !
  ! !LOCAL VARIABLES:
  integer :: length
  integer :: i
  character(len=*),parameter :: subname = 'getname'
  !-----------------------------------------------------------------------

    length = len (inname)

    if (length < max_namlen .or. length > max_namlen+2) then
      write(iulog,*) trim(subname),' ERROR: bad length=',length
      call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    getname = ' '
    do i = 1,max_namlen
      if (inname(i:i) == ':') exit
      getname(i:i) = inname(i:i)
    end do
  end function getname

  !-----------------------------------------------------------------------
  character(len=1) function getflag (inname)
    !
    ! !DESCRIPTION:
    ! Retrieve flag portion of inname. If an averaging flag separater character
    ! is present (:) in inname, return the character after it as the flag
    !
    ! !ARGUMENTS:
    character(len=*) :: inname   ! character string
    !
    ! !LOCAL VARIABLES:
    integer :: length         ! length of inname
    integer :: i              ! loop index
    character(len=*),parameter :: subname = 'getflag'
    !-----------------------------------------------------------------------

    length = len (inname)

    if (length < max_namlen .or. length > max_namlen+2) then
      write(iulog,*) trim(subname),' ERROR: bad length=',length
      call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    getflag = ' '
    do i = 1,length
      if (inname(i:i) == ':') then
          getflag = inname(i+1:i+1)
          exit
      end if
    end do

  end function getflag

  !-----------------------------------------------------------------------
  subroutine htape_addfld (t, f, avgflag)
    !
    ! !DESCRIPTION:
    ! Add a field to the active list for a history tape. Copy the data from
    ! the master field list to the active list for the tape.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t                 ! history tape index
    integer, intent(in) :: f                 ! field index from master field list
    character(len=1), intent(in) :: avgflag  ! time averaging flag
    !
    ! !LOCAL VARIABLES:
    integer :: n                    ! field index on defined tape
    character(len=hist_dim_name_length) :: type1d      ! clm pointer 1d type
    character(len=hist_dim_name_length) :: type1d_out  ! history buffer 1d type
    integer :: numg                 ! total number of gridcells across all processors
    integer :: numt                 ! total number of topounits across all processors
    integer :: numc                 ! total number of columns across all processors
    integer :: nump                 ! total number of pfts across all processors
    integer :: num2d                ! size of second dimension (e.g. .number of vertical levels)
    integer :: beg1d_out,end1d_out  ! history output per-proc 1d beginning and ending indices
    integer :: num1d_out            ! history output 1d size

    character(len=*),parameter :: subname = 'htape_addfld'
    !-----------------------------------------------------------------------

    ! Ensure that it is not to late to add a field to the history tape

    if (htapes_defined) then
       write(iulog,*) trim(subname),' ERROR: attempt to add field ', &
            masterlist(f)%field%name, ' after history files are set'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    tape(t)%nflds = tape(t)%nflds + 1
    n = tape(t)%nflds

    ! Copy field information

    tape(t)%hlist(n)%field = masterlist(f)%field

    ! Determine bounds

    numg=bounds%ngrid;numt=bounds%ntopou;numc=bounds%ncols;nump=bounds%npfts

    ! Modify type1d_out if necessary

    if (hist_type1d_pertape(t) /= ' ') then

       ! Set output 1d type  based on namelist setting of  hist_type1d_pertape
       ! Only applies to tapes when xy output is not required

       type1d = tape(t)%hlist(n)%field%type1d

       select case (trim(hist_type1d_pertape(t)))
       case('GRID')
          tape(t)%hlist(n)%field%type1d_out = nameg
       case('TOPO')
          tape(t)%hlist(n)%field%type1d_out = namet
       case('COLS')
          tape(t)%hlist(n)%field%type1d_out = namec
       case ('PFTS')
          tape(t)%hlist(n)%field%type1d_out = namep
       case default
          write(iulog,*) trim(subname),' ERROR: unknown input hist_type1d_pertape= ', hist_type1d_pertape(t)
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select

    endif

    ! Determine output 1d dimensions

    type1d_out = tape(t)%hlist(n)%field%type1d_out
    if (type1d_out == nameg) then
       beg1d_out = bounds%begg
       end1d_out = bounds%endg
       num1d_out = numg
    else if (type1d_out == namet) then
       beg1d_out = bounds%begt
       end1d_out = bounds%endt
       num1d_out = numt
    else if (type1d_out == namec) then
       beg1d_out = bounds%begc
       end1d_out = bounds%endc
       num1d_out = numc
    else if (type1d_out == namep) then
       beg1d_out = bounds%begp
       end1d_out = bounds%endp
       num1d_out = nump
    else
       write(iulog,*) trim(subname),' ERROR: incorrect value of type1d_out= ',type1d_out
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    tape(t)%hlist(n)%field%beg1d_out = beg1d_out
    tape(t)%hlist(n)%field%end1d_out = end1d_out
    tape(t)%hlist(n)%field%num1d_out = num1d_out

    ! Allocate and initialize history buffer and related info

    num2d = tape(t)%hlist(n)%field%num2d
    allocate (tape(t)%hlist(n)%hbuf(beg1d_out:end1d_out,num2d))
    allocate (tape(t)%hlist(n)%nacs(beg1d_out:end1d_out,num2d))
    tape(t)%hlist(n)%hbuf(:,:) = 0._r8
    tape(t)%hlist(n)%nacs(:,:) = 0

    ! Set time averaging flag based on masterlist setting or
    ! override the default averaging flag with namelist setting

    select case (avgflag)
    case (' ')
       tape(t)%hlist(n)%avgflag = masterlist(f)%avgflag(t)
    case ('A','I','X','M','P')
       tape(t)%hlist(n)%avgflag = avgflag
    case default
       write(iulog,*) trim(subname),' ERROR: unknown avgflag=', avgflag
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine htape_addfld

  !-----------------------------------------------------------------------
  subroutine masterlist_change_timeavg (t)
    !
    ! !DESCRIPTION:
    ! Override default history tape contents for a specific tape.
    ! Copy the flag into the master field list.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t         ! history tape index
    !
    ! !LOCAL VARIABLES:
    integer :: f                     ! field index
    character(len=1) :: avgflag      ! lcl equiv of hist_avgflag_pertape(t)
    character(len=*),parameter :: subname = 'masterlist_change_timeavg'
    !-----------------------------------------------------------------------

    avgflag = hist_avgflag_pertape(t)

    do f = 1,nfmaster
       select case (avgflag)
       case ('A')  !average
          masterlist(f)%avgflag(t) = avgflag
       case ('P')  !average over positive values
          masterlist(f)%avgflag(t) = avgflag
       case ('I')  !instantaneous
          masterlist(f)%avgflag(t) = avgflag
       case ('X')  !maximum
          masterlist(f)%avgflag(t) = avgflag
       case ('M')  !minimum
          masterlist(f)%avgflag(t) = avgflag
       case default
          write(iulog,*) trim(subname),' ERROR: unknown avgflag=',avgflag
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    end do

  end subroutine masterlist_change_timeavg

  !-----------------------------------------------------------------------
  subroutine hist_update_hbuf(bounds)
    !
    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.
    !
    ! !ARGUMENTS:
    use GridConsts, only : bounds_type
    implicit none
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: t                   ! tape index
    integer :: f                   ! field index
    integer :: numdims             ! number of dimensions
    integer :: num2d               ! size of second dimension (e.g. number of vertical levels)
    character(len=*),parameter :: subname = 'hist_update_hbuf'
    !-----------------------------------------------------------------------

    do t = 1,ntapes
!!$OMP PARALLEL DO PRIVATE (f, numdims, num2d)
       do f = 1,tape(t)%nflds
          numdims = tape(t)%hlist(f)%field%numdims
          if ( numdims == 1) then
             call hist_update_hbuf_field_1d (t, f, bounds)
          else
             num2d = tape(t)%hlist(f)%field%num2d
             call hist_update_hbuf_field_2d (t, f, bounds, num2d)
          end if
       end do
!!$OMP END PARALLEL DO
    end do

  end subroutine hist_update_hbuf

  !-----------------------------------------------------------------------
  subroutine hist_update_hbuf_field_1d (t, f, bounds)
    !
    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.
    !
    ! This canNOT be called from within a threaded region (see comment below regarding the
    ! call to p2g, and the lack of explicit bounds on its arguments; see also bug 1786)
    !
    ! !USES:
    use GridConsts, only : bounds_type
    use EcoSIMCtrlMod, only: lverb    
    implicit none
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t            ! tape index
    integer, intent(in) :: f            ! field index
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: hpindex                 ! history pointer index
    integer  :: k                       ! gridcell, landunit, column or pft index
    integer  :: beg1d,end1d             ! beginning and ending indices
    logical  :: check_active            ! true => check 'active' flag of each point (this refers to a point being active, NOT a history field being active)
    logical  :: valid                   ! true => history operation is valid
    logical  :: map2gcell               ! true => map clm pointer field to gridcell
    character(len=hist_dim_name_length)  :: type1d         ! 1d clm pointerr type   ["gridcell","landunit","column","pft"]
    character(len=hist_dim_name_length)  :: type1d_out     ! 1d history buffer type ["gridcell","landunit","column","pft"]
    character(len=1)  :: avgflag        ! time averaging flag
    real(r8), pointer :: hbuf(:,:)      ! history buffer
    integer , pointer :: nacs(:,:)      ! accumulation counter
    real(r8), pointer :: field(:)       ! clm 1d pointer field
    logical , pointer :: active(:)      ! flag saying whether each point is active (used for type1d = landunit/column/pft) (this refers to a point being active, NOT a history field being active)
    real(r8) :: field_gcell(bounds%begg:bounds%endg)  ! gricell level field (used if mapping to gridcell is done)
    integer :: j
    character(len=*),parameter :: subname = 'hist_update_hbuf_field_1d'
    integer :: k_offset                    ! offset for mapping sliced subarray pointers when outputting variables in PFT/col vector form
    !-----------------------------------------------------------------------

    avgflag        =  tape(t)%hlist(f)%avgflag
    nacs           => tape(t)%hlist(f)%nacs
    hbuf           => tape(t)%hlist(f)%hbuf
    beg1d          =  tape(t)%hlist(f)%field%beg1d
    end1d          =  tape(t)%hlist(f)%field%end1d
    type1d         =  tape(t)%hlist(f)%field%type1d
    type1d_out     =  tape(t)%hlist(f)%field%type1d_out
    hpindex        =  tape(t)%hlist(f)%field%hpindex
    field          => esmptr_rs(hpindex)%ptr
    call PrintInfo('beg '//subname)    
    ! set variables to check weights when allocate all pfts

       ! For data defined on the pft, col, and landunit we need to check if a point is active
       ! to determine whether that point should be assigned spval
!       if (type1d == namep) then
!          check_active = .true.
!          active => veg_pp%active
!       else if (type1d == namec) then
!          check_active = .true.
!          active => col_pp%active
!       else if (type1d == namel) then
!          check_active = .true.
!          active =>lun_pp%active
!       else if (type1d == namet) then
!          check_active = .true.
!          active =>top_pp%active
!       else
          check_active = .false.
!       end if

       select case (avgflag)
       case ('I') ! Instantaneous
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   hbuf(k,1) = field(k)
                else
                   hbuf(k,1) = spval
                end if
             else
                hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case ('A') ! Time average
          ! create mappings for array slice pointers (which go from 1 to size(field) rather than beg1d to end1d)
          if ( end1d .eq. ubound(field,1) ) then
             k_offset = 0
          else
             k_offset = 1 - beg1d
          endif
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k+k_offset) /= spval) then   ! add k_offset
                   if (nacs(k,1) == 0) hbuf(k,1) = 0._r8
                   hbuf(k,1) = hbuf(k,1) + field(k+k_offset)   ! add k_offset
                   nacs(k,1) = nacs(k,1) + 1
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
          end do
       case ('P') !average over positive numbers
          if ( end1d .eq. ubound(field,1) ) then
             k_offset = 0
          else
             k_offset = 1 - beg1d
          endif
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k+k_offset) /= spval .and. field(k+k_offset)>1.e-8_r8) then   ! add k_offset
                   if (nacs(k,1) == 0) hbuf(k,1) = 0._r8
                   hbuf(k,1) = hbuf(k,1) + field(k+k_offset)   ! add k_offset
                   nacs(k,1) = nacs(k,1) + 1
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
          end do             
       case ('X') ! Maximum over time
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   if (nacs(k,1) == 0) hbuf(k,1) = -1.e50_r8
                   hbuf(k,1) = max( hbuf(k,1), field(k) )
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case ('M') ! Minimum over time
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   if (nacs(k,1) == 0) hbuf(k,1) = +1.e50_r8
                   hbuf(k,1) = min( hbuf(k,1), field(k) )
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case default
          write(iulog,*) trim(subname),' ERROR: invalid time averaging flag ', avgflag
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
  call PrintInfo('end '//subname)
  end subroutine hist_update_hbuf_field_1d
  !-----------------------------------------------------------------------
  subroutine hist_update_hbuf_field_2d (t, f, bounds, num2d)
    !
    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.
    !
    ! This canNOT be called from within a threaded region (see comment below regarding the
    ! call to p2g, and the lack of explicit bounds on its arguments; see also bug 1786)
    !
    !
    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: t            ! tape index
    integer, intent(in) :: f            ! field index
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num2d        ! size of second dimension
    !
    ! !LOCAL VARIABLES:
    integer  :: hpindex                 ! history pointer index
    integer  :: k                       ! gridcell, landunit, column or pft index
    integer  :: j                       ! level index
    integer  :: beg1d,end1d             ! beginning and ending indices
    logical  :: check_active            ! true => check 'active' flag of each point (this refers to a point being active, NOT a history field being active)
    logical  :: valid                   ! true => history operation is valid
    logical  :: map2gcell               ! true => map clm pointer field to gridcell
    character(len=hist_dim_name_length)  :: type1d         ! 1d clm pointerr type   ["gridcell","landunit","column","pft"]
    character(len=hist_dim_name_length)  :: type1d_out     ! 1d history buffer type ["gridcell","landunit","column","pft"]
    character(len=1)  :: avgflag        ! time averaging flag
    integer  :: no_snow_behavior        ! for multi-layer snow fields, behavior to use when a given layer is absent
    real(r8), pointer :: hbuf(:,:)      ! history buffer
    integer , pointer :: nacs(:,:)      ! accumulation counter
    real(r8), pointer :: field(:,:)     ! clm 2d pointer field
    logical           :: field_allocated! whether 'field' was allocated here
    logical , pointer :: active(:)      ! flag saying whether each point is active (used for type1d = landunit/column/pft)
                                        !(this refers to a point being active, NOT a history field being active)
    real(r8) :: field_gcell(bounds%begg:bounds%endg,num2d) ! gricell level field (used if mapping to gridcell is done)
    character(len=*),parameter :: subname = 'hist_update_hbuf_field_2d'
    !-----------------------------------------------------------------------

    avgflag             =  tape(t)%hlist(f)%avgflag
    nacs                => tape(t)%hlist(f)%nacs
    hbuf                => tape(t)%hlist(f)%hbuf
    beg1d               =  tape(t)%hlist(f)%field%beg1d
    end1d               =  tape(t)%hlist(f)%field%end1d
    type1d              =  tape(t)%hlist(f)%field%type1d
    type1d_out          =  tape(t)%hlist(f)%field%type1d_out
    no_snow_behavior    =  tape(t)%hlist(f)%field%no_snow_behavior
    hpindex             =  tape(t)%hlist(f)%field%hpindex

    call PrintInfo('beg '//subname)
!    if (no_snow_behavior /= no_snow_unset) then
       ! For multi-layer snow fields, build a special output variable that handles
       ! missing snow layers appropriately

       ! Note, regarding bug 1786: The following allocation is not what we would want if
       ! this routine were operating in a threaded region (or, more generally, within a
       ! loop over nclumps) - in that case we would want to use the bounds information for
       ! this clump. But currently that's not possible because the bounds of some fields
       ! have been reset to 1 - see also bug 1786. Similarly, if we wanted to allow
       ! operation within a loop over clumps, we would need to pass 'bounds' to
       ! hist_set_snow_field_2d rather than relying on beg1d & end1d (which give the proc,
       ! bounds not the clump bounds)

!       allocate(field(lbound(elmptr_ra(hpindex)%ptr, 1) : ubound(elmptr_ra(hpindex)%ptr, 1), 1:num2d))
!       field_allocated = .true.

!       call hist_set_snow_field_2d(field, esmptr_ra(hpindex)%ptr, no_snow_behavior, type1d, &
!            beg1d, end1d)
!    else
       field => esmptr_ra1(hpindex)%ptr(:,1:num2d)
       field_allocated = .false.
!    end if

    ! set variables to check weights when allocate all pfts

    ! Do not map to gridcell

       ! For data defined on the pft, col or landunit, we need to check if a point is active
       ! to determine whether that point should be assigned spval
!       if (type1d == namep) then
!          check_active = .true.
!          active => veg_pp%active
!       else if (type1d == namec) then
!          check_active = .true.
!          active => col_pp%active
!       else if (type1d == namel) then
!          check_active = .true.
!          active =>lun_pp%active
!       else if (type1d == namet) then
!          check_active = .true.
!          active =>top_pp%active
!       else
          check_active = .false.
!       end if

       ! Note that since field points to an array section the
       ! bounds are field(1:end1d-beg1d+1, num2d) - therefore
       ! need to do the shifting below

       select case (avgflag)
       case ('I') ! Instantaneous
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      hbuf(k,j) = field(k-beg1d+1,j)
                   else
                      hbuf(k,j) = spval
                   end if
                else
                   hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case ('A') ! Time average
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                      hbuf(k,j) = hbuf(k,j) + field(k-beg1d+1,j)
                      nacs(k,j) = nacs(k,j) + 1
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
             end do
          end do
       case ('P') ! Time average over values > 1.e-8
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval .and. field(k-beg1d+1,j)>1.e-8_r8) then
                      if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                      hbuf(k,j) = hbuf(k,j) + field(k-beg1d+1,j)
                      nacs(k,j) = nacs(k,j) + 1
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
             end do
          end do          
       case ('X') ! Maximum over time
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = -1.e50_r8
                      hbuf(k,j) = max( hbuf(k,j), field(k-beg1d+1,j) )
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case ('M') ! Minimum over time
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = +1.e50_r8
                      hbuf(k,j) = min( hbuf(k,j), field(k-beg1d+1,j))
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case default
          write(iulog,*) trim(subname),' ERROR: invalid time averaging flag ', avgflag
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select


    if (field_allocated) then
       deallocate(field)
    end if
  call PrintInfo('end '//subname)
  end subroutine hist_update_hbuf_field_2d

  !-----------------------------------------------------------------------
  subroutine hist_htapes_wrapup( rstwr, nlend, bounds, lnyr )

    ! Description:
    ! Write history tape(s)
    ! Determine if next time step is beginning of history interval and if so:
    !   increment the current time sample counter, open a new history file
    !   and if needed (i.e., when ntim = 1), write history data to current
    !   history file, reset field accumulation counters to zero.
    ! If primary history file is full or at the last time step of the simulation,
    !   write restart dataset and close all history fiels.
    ! If history file is full or at the last time step of the simulation:
    !   close history file
    !   and reset time sample counter to zero if file is full.
    ! Daily-averaged data for the first day in September are written on
    !   date = 00/09/02 with mscur = 0.
    ! Daily-averaged data for the first day in month mm are written on
    !   date = yyyy/mm/02 with mscur = 0.
    ! Daily-averaged data for the 30th day (last day in September) are written
    !   on date = 0000/10/01 mscur = 0.
    ! Daily-averaged data for the last day in month mm are written on
    !   date = yyyy/mm+1/01 with mscur = 0.
    !
    ! !USES:

!    use perf_mod        , only : t_startf, t_stopf
    !
    implicit none
    ! !ARGUMENTS:
    logical, intent(in) :: rstwr    ! true => write restart file this step
    logical, intent(in) :: nlend    ! true => end of run on this step
    logical, intent(in) :: lnyr     ! true => close current hist file and open a new one
    type(bounds_type) , intent(in) :: bounds
!    real(r8)          , intent(in) :: watsat_col( bounds%begc:,1: )
!    real(r8)          , intent(in) :: sucsat_col( bounds%begc:,1: )
!    real(r8)          , intent(in) :: bsw_col( bounds%begc:,1: )
!    real(r8)          , intent(in) :: hksat_col( bounds%begc:,1: )
    !
    ! !LOCAL VARIABLES:
    integer :: t                          ! tape index
    integer :: f                          ! field index
    integer :: ier                        ! error code
    integer :: nstep                      ! current step
    integer :: day                        ! current day (1 -> 31)
    integer :: mon                        ! current month (1 -> 12)
    integer :: yr                         ! current year (0 -> ...)
    integer :: mdcur                      ! current day
    integer :: mscur                      ! seconds of current day
    integer :: mcsec                      ! current time of day [seconds]
    integer :: daym1                      ! nstep-1 day (1 -> 31)
    integer :: monm1                      ! nstep-1 month (1 -> 12)
    integer :: yrm1                       ! nstep-1 year (0 -> ...)
    integer :: mcsecm1                    ! nstep-1 time of day [seconds]
    real(r8):: time                       ! current time
    character(len=256) :: str             ! global attribute string
    logical :: if_stop                    ! true => last time step of run
    logical, save :: do_3Dtconst = .true. ! true => write out 3D time-constant data
    logical :: rollover_tape(ntapes)       !close this tape and start its successor [-]
    character(len=14) :: ymdhs    
    character(len=*),parameter :: subname = trim(mod_filename)//'::hist_htapes_wrapup'
    !-----------------------------------------------------------------------

    call PrintInfo('beg '//subname)
    ! Monthly tapes skip step zero, so initialize their rollover decision too.
    rollover_tape(:)=.false.
    ! get current step
    
    nstep = etimer%get_nstep()

    ! Set calendar for current time step

    call etimer%get_curr_date (yr, mon, day, mcsec)
    time=etimer%get_curr_time ()/secspday

    ! Set calendar for current for previous time step

    call etimer%get_prev_date (yrm1, monm1, daym1, mcsecm1)

    ! Loop over active history tapes, create new history files if necessary
    ! and write data to history files if end of history interval.
    
    do t = 1, ntapes

       ! Skip nstep=0 if monthly average

       if (nstep==0 .and. tape(t)%nhtfrq==0) cycle

       ! Determine if end of history interval
       tape(t)%is_endhist = .false.
       if (tape(t)%nhtfrq==0) then   !monthly average
          !by year for monthly output
          if (mon /= monm1) tape(t)%is_endhist = .true.
       else
          if (mod(nstep,tape(t)%nhtfrq) == 0) tape(t)%is_endhist = .true.
       end if

       ! If end of history interval
       if (tape(t)%is_endhist) then

          ! Normalize history buffer if time averaged

          call hfields_normalize(t)

          ! Increment current time sample counter.

          tape(t)%ntimes = tape(t)%ntimes + 1

          ! Create history file if appropriate and build time comment

          ! If first time sample, generate unique history file name, open file,
          ! define dims, vars, etc.

          if (tape(t)%ntimes == 1) then
!             call t_startf('hist_htapes_wrapup_define')
             locfnh(t) = set_hist_filename (hist_freq=tape(t)%nhtfrq,hist_mfilt=tape(t)%mfilt, hist_file=t)
             print*,'locfnh(t), t=',t,locfnh(t)
!             if (masterproc) then
                write(iulog,*) trim(subname),' : Creating history file ', trim(locfnh(t)), &
                     ' at nstep = ',etimer%get_nstep()
!             endif
             !'call htape_create'
             call htape_create (t)

             ! Define time-constant field variables
             call htape_timeconst(t, mode='define')

             ! Define 3D time-constant field variables only to first primary tape
!             if ( do_3Dtconst .and. t == 1 ) then
!                call htape_timeconst3D(t, &
!                     bounds, watsat_col, sucsat_col, bsw_col, hksat_col, mode='define')
!                TimeConst3DVars_Filename = trim(locfnh(t))
!             end if

             !' Define model field variables'
             call hfields_write(t, mode='define')

             !' Exit define model'
             call ncd_enddef(nfid(t))
!             call t_stopf('hist_htapes_wrapup_define')
          endif

!          call t_startf('hist_htapes_wrapup_tconst')
          ! Write time constant history variables
          call htape_timeconst(t, mode='write')

          ! Write 3D time constant history variables only to first primary tape
!          if ( do_3Dtconst .and. t == 1 .and. tape(t)%ntimes == 1 )then
!             call htape_timeconst3D(t, &
!                  bounds, watsat_col, sucsat_col, bsw_col, hksat_col, mode='write')
!             do_3Dtconst = .false.
!          end if
!         if (masterproc) then
             call etimer%get_ymdhs(ymdhs)
             write(iulog,*)
             write(iulog,*) trim(subname),' : Writing current time sample to local history file ', &
                  trim(locfnh(t)),' at = ',ymdhs, ' doy=',etimer%get_curr_doy(),&
                  ' for history time interval beginning at ', tape(t)%begtime, &
                  ' and ending at ',time
             write(iulog,*)
             call flush(iulog)
!          endif

          ! Update beginning time of next interval
          tape(t)%begtime = time
!          call t_stopf('hist_htapes_wrapup_tconst')

          ! Write history time samples
!          call t_startf('hist_htapes_wrapup_write')
          call hfields_write(t, mode='write')
!          call t_stopf('hist_htapes_wrapup_write')

          ! Zero necessary history buffers
          call hfields_zero(t)

       end if
      ! Use one rollover decision for closing, reopening and resetting counts.
      ! Non-negative frequencies start a fresh tape each year even when the
      ! preceding file is partial (e.g. the last 24 hours of a leap year).
      rollover_tape(t)=tape(t)%ntimes>0 .and. &
        ((tape(t)%nhtfrq>=0 .and. (tape(t)%ntimes>=tape(t)%mfilt .or. lnyr)) .or. &
         (tape(t)%nhtfrq<0 .and. tape(t)%ntimes>=tape(t)%mfilt .and. lnyr))
    end do  ! end loop over history tapes
    ! Determine if file needs to be closed
    
    if(ntapes>0)then
      call hist_do_disp (ntapes, rollover_tape, if_stop, if_disphist, rstwr, nlend)
    endif

    ! Close open history file
    ! Auxilary files may have been closed and saved off without being full,
    ! must reopen the files

    do t = 1, ntapes
       if (if_disphist(t)) then
          if (tape(t)%ntimes /= 0) then
!             if (masterproc) then
                write(iulog,*)
                write(iulog,*)  trim(subname),' : Closing local history file ',&
                     trim(locfnh(t)),' at nstep = ', etimer%get_nstep()
                write(iulog,*)
!             endif
             call ncd_pio_closefile(nfid(t))
             ! Restart flushes continue the same tape; a rolled-over tape is
             ! already complete and must stay closed when its count is reset.
             if (.not.if_stop .and. .not.rollover_tape(t)) then
                call ncd_pio_openfile (nfid(t), trim(locfnh(t)), ncd_write)
             end if
          else
 !            if (masterproc) then
!                write(iulog,*) trim(subname),' : history tape ',t,': no open file to close'
!                print*,tape(t)%ntimes,tape(t)%mfilt,rstwr, nlend
 !            end if
          endif
       endif
    end do
    ! Reset number of time samples to zero if file is full
    ! make sure for high frequency output, daily or hourly, are recorded by whole year 
    do t = 1, ntapes
       if (rollover_tape(t)) then
          tape(t)%ntimes = 0
       end if
    end do
  call PrintInfo('end '//subname)
  end subroutine hist_htapes_wrapup

  !-----------------------------------------------------------------------
  subroutine htape_timeconst(t, mode)  
  implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: t              ! tape index
    character(len=*), intent(in) :: mode  ! 'define' or 'write'


  real(r8),pointer:: timedata(:)                ! time interval boundaries
  integer :: dim2id(2)           ! netCDF dimension id
  integer :: varid                             ! variable id
  integer :: nt

  if (tape(t)%ntimes == 1) then
    if (mode == 'define') then
    elseif (mode == 'write') then
    endif
  endif
   
  ! For define mode -- only do this for first time-sample
  if (mode == 'define' .and. tape(t)%ntimes == 1) then
    
    call ncd_defvar(nfid(t), 'time_bounds', ncd_double, dim1name='hist_interval', &
      dim2name= 'time', long_name = 'history time interval endpoints')   

  elseif (mode == 'write') then
    allocate(timedata(2))
    timedata(1) = tape(t)%begtime
    timedata(2) = etimer%get_curr_time ()/secspday
    nt         = tape(t)%ntimes

    call ncd_io(flag='write', varname='time_bounds', &
      dim1name='hist_interval', data=timedata, ncid=nfid(t), nt=nt)

    deallocate(timedata)
  endif       
  end subroutine htape_timeconst
  !-----------------------------------------------------------------------
  subroutine hfields_write(t, mode)
    !
    ! !DESCRIPTION:
    ! Write history tape.  Issue the call to write the variable.
    !
    ! !USES:
    !
    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: t                ! tape index
    character(len=*), intent(in) :: mode    ! 'define' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: f                         ! field index
    integer :: k                         ! 1d index
    integer :: topo,c,l,p                ! indices
    integer :: beg1d_out                 ! on-node 1d hbuf pointer start index
    integer :: end1d_out                 ! on-node 1d hbuf pointer end index
    integer :: num1d_out                 ! size of hbuf first dimension (overall all nodes)
    integer :: num2d                     ! hbuf second dimension size
    integer :: nt                        ! time index
    integer :: ier                       ! error status
    integer :: numdims                   ! number of dimensions
    character(len=1)         :: avgflag  ! time averaging flag
    character(len=max_chars) :: long_name! long name
    character(len=max_chars) :: standard_name! standard name
    character(len=max_chars) :: units    ! units
    character(len=max_namlen):: varname  ! variable name
    character(len=32) :: avgstr          ! time averaging type
    character(len=hist_dim_name_length)  :: type1d_out      ! history output 1d type
    character(len=hist_dim_name_length)  :: type2d          ! history output 2d type
    character(len=32) :: dim1name        ! temporary
    character(len=32) :: dim2name        ! temporary
    real(r8), pointer :: histo(:,:)      ! temporary
    real(r8), pointer :: hist1do(:)      ! temporary
    character(len=*),parameter :: subname = trim(mod_filename)//'::hfields_write'
!-----------------------------------------------------------------------
    ! Write/define 1d topological info

    call PrintInfo('beg '//subname)
!    if (mode == 'define') then
!      call hfields_1dinfo(t, mode='define')
!    else if (mode == 'write') then
!      call hfields_1dinfo(t, mode='write')
!    end if

    ! Define time-dependent variables create variables and attributes for field list

    do f = 1,tape(t)%nflds

       ! Set history field variables

       varname    = tape(t)%hlist(f)%field%name
       long_name  = tape(t)%hlist(f)%field%long_name
       standard_name  = tape(t)%hlist(f)%field%standard_name
       units      = tape(t)%hlist(f)%field%units
       avgflag    = tape(t)%hlist(f)%avgflag
       type1d_out = tape(t)%hlist(f)%field%type1d_out
       beg1d_out  = tape(t)%hlist(f)%field%beg1d_out
       end1d_out  = tape(t)%hlist(f)%field%end1d_out
       num1d_out  = tape(t)%hlist(f)%field%num1d_out
       type2d     = tape(t)%hlist(f)%field%type2d
       numdims    = tape(t)%hlist(f)%field%numdims
       num2d      = tape(t)%hlist(f)%field%num2d
       nt         = tape(t)%ntimes
       
       if (mode == 'define') then

          select case (avgflag)
          case ('A')
             avgstr = 'mean'
          case ('P')
             avgstr = 'positive_mean'   
          case ('I')
             avgstr = 'point'
          case ('X')
             avgstr = 'maximum'
          case ('M')
             avgstr = 'minimum'
          case default
             write(iulog,*) trim(subname),' ERROR: unknown time averaging flag (avgflag)=',avgflag
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select

          dim1name = type1d_out ; dim2name = 'undefined'

          if (dim2name == 'undefined') then

             if (numdims == 1) then
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name='time', &
                     long_name=long_name, standard_name=standard_name,units=units,&
                     cell_method=avgstr,  missing_value=spval, fill_value=spval)
             else
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name=type2d, dim3name='time', &
                     long_name=long_name,standard_name=standard_name, units=units, &
                     cell_method=avgstr,  missing_value=spval, fill_value=spval)
             end if
          else
             if (numdims == 1) then
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name=dim2name, dim3name='time', &
                     long_name=long_name, standard_name=standard_name, &
                     units=units, cell_method=avgstr, missing_value=spval, fill_value=spval)
             else
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name=dim2name, dim3name=type2d, dim4name='time', &
                     long_name=long_name, standard_name=standard_name, units=units,&
                     cell_method=avgstr, missing_value=spval, fill_value=spval)
             end if
          endif

       else if (mode == 'write') then

          ! Determine output buffer

          histo => tape(t)%hlist(f)%hbuf
          ! Allocate dynamic memory

          if (numdims == 1) then
             allocate(hist1do(beg1d_out:end1d_out), stat=ier)
             if (ier /= 0) then
                write(iulog,*) trim(subname),' ERROR: allocation'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end if
             hist1do(beg1d_out:end1d_out) = histo(beg1d_out:end1d_out,1)
          end if

          ! Write history output.  Always output land and ocean runoff on xy grid.
          call PrintInfo(varname)
          
          if (numdims == 1) then
!             print*,varname,hist1do(beg1d_out:end1d_out)
             call ncd_io(flag='write', varname=varname, &
                  dim1name=type1d_out, data=hist1do, ncid=nfid(t), nt=nt)
          else
!             if(lverb)then
!               write(*,*)varname,histo(1,:)
!             endif
             call ncd_io(flag='write', varname=varname, &
                  dim1name=type1d_out, data=histo, ncid=nfid(t), nt=nt)
          end if

          ! Deallocate dynamic memory

          if (numdims == 1) then
             deallocate(hist1do)
          end if

       end if

   end do
   call PrintInfo('end '//subname//' '//mode)
  end subroutine hfields_write

  !------------------------------------------------------------------------
  subroutine hist_do_disp (ntapes, rollover_tape, if_stop, if_disphist, rstwr, nlend)
    ! Description:
    ! Close rolled-over tapes and flush all open tapes for restart or run end.
    implicit none
    integer, intent(in) :: ntapes               !number of active history tapes [-]
    logical, intent(in) :: rollover_tape(ntapes) !true => begin a successor file [-]
    logical, intent(out) :: if_stop              !true => final simulation step [-]
    logical, intent(out) :: if_disphist(ntapes)   !true => close/flush this tape [-]
    logical, intent(in) :: rstwr                 !true => write a restart this step [-]
    logical, intent(in) :: nlend                 !true => end the run this step [-]
    character(len=*), parameter :: subname='hist_do_disp'
    call PrintInfo('beg '//subname)
    if_stop=nlend
    if_disphist(:)=rollover_tape(:).or.rstwr.or.nlend
    call PrintInfo('end '//subname)
  end subroutine hist_do_disp

   !-----------------------------------------------------------------------
   character(len=256) function set_hist_filename (hist_freq, hist_mfilt, hist_file)
     !
     ! !DESCRIPTION:
     ! Determine history dataset filenames.
     !
     ! !USES:
     !
     implicit none
     ! !ARGUMENTS:
     integer, intent(in)  :: hist_freq   !history file frequency
     integer, intent(in)  :: hist_mfilt  !history file number of time-samples
     integer, intent(in)  :: hist_file   !history file index

     !
     ! !LOCAL VARIABLES:
     !EOP
     character(len=256) :: cdate       !date char string
     character(len=  1) :: hist_index  !p,1 or 2 (currently)
     integer :: day                    !day (1 -> 31)
     integer :: mon                    !month (1 -> 12)
     integer :: yr                     !year (0 -> ...)
     integer :: sec                    !seconds into current day
     character(len=1) :: inst_suffix=' '
     character(len=*),parameter :: subname = 'set_hist_filename'

   if (hist_freq == 0) then   !monthly
     call etimer%get_prev_date (yr, mon, day, sec)

     if(hist_mfilt /= 1)mon=mon-1       
     write(cdate,'(i4.4,"-",i2.2)') yr,mon
   else                        !other
      call etimer%get_curr_date (yr, mon, day, sec)
      if(sec/=0)sec=sec-etimer%get_step_size()
      write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day-1,sec
   endif
   write(hist_index,'(i1.1)') hist_file - 1
   set_hist_filename = "./"//trim(case_name)//".ecosim"//trim(inst_suffix)//&
                       ".h"//hist_index//"."//trim(cdate)//".nc"

  end function set_hist_filename

!-----------------------------------------------------------------------
  subroutine hist_restart_ncd (bounds, ncid, flag, rdate)
   !
   ! !DESCRIPTION:
   ! Read/write history file restart data.
   ! If the current history file(s) are not full, file(s) are opened
   ! so that subsequent time samples are added until the file is full.
   ! A new history file is used on a branch run.
   !
   ! !USES:
   use EcoSIMConfig     , only : nsrest, nsrStartup, nsrBranch, is_restart
   use fileutil         , only : getfil
   !
   ! !ARGUMENTS:
   type(bounds_type), intent(in)    :: bounds
   type(file_desc_t), intent(inout) :: ncid     ! netcdf file
   character(len=*) , intent(in)    :: flag     !'read' or 'write'
   character(len=*) , intent(in), optional :: rdate    ! restart file time stamp for name
   !
   ! !LOCAL VARIABLES:
   integer :: max_nflds                     ! Max number of fields
   integer :: num1d,beg1d,end1d             ! 1d size, beginning and ending indices
   integer :: num1d_out,beg1d_out,end1d_out ! 1d size, beginning and ending indices
   integer :: num2d                         ! 2d size (e.g. number of vertical levels)
   integer :: numg                 ! total number of gridcells across all processors
   integer :: numt                 ! total number of topounits across all processors
   integer :: numc                 ! total number of columns across all processors
   integer :: nump                 ! total number of pfts across all processors
   character(len=max_namlen) :: name            ! variable name
   character(len=max_namlen) :: name_acc        ! accumulator variable name
   character(len=max_namlen) :: long_name       ! long name of variable
   character(len=max_namlen) :: standard_name       ! standard_name of var
   character(len=max_chars)  :: long_name_acc   ! long name for accumulator
   character(len=max_chars)  :: units           ! units of variable
   character(len=max_chars)  :: units_acc       ! accumulator units
   character(len=max_chars)  :: fname           ! full name of history file
   character(len=max_chars)  :: locrest(max_tapes) ! local history restart file names

   character(len=max_namlen),allocatable :: tname(:)
   character(len=max_chars), allocatable :: tunits(:),tlongname(:)
   character(len=hist_dim_name_length), allocatable :: tmpstr(:,:)
   character(len=1), allocatable :: tavgflag(:)
   integer :: start(2)

    character(len=1)   :: hnum                   ! history file index
    character(len=hist_dim_name_length)   :: type1d                 ! clm pointer 1d type
    character(len=hist_dim_name_length)   :: type1d_out             ! history buffer 1d type
    character(len=hist_dim_name_length)   :: type2d                 ! history buffer 2d type
    character(len=32)  :: dim1name               ! temporary
    character(len=32)  :: dim2name               ! temporary
    type(var_desc_t)   :: name_desc              ! variable descriptor for name
    type(var_desc_t)   :: longname_desc          ! variable descriptor for long_name
    type(var_desc_t)   :: units_desc             ! variable descriptor for units
    type(var_desc_t)   :: type1d_desc            ! variable descriptor for type1d
    type(var_desc_t)   :: type1d_out_desc        ! variable descriptor for type1d_out
    type(var_desc_t)   :: type2d_desc            ! variable descriptor for type2d
    type(var_desc_t)   :: avgflag_desc           ! variable descriptor for avgflag
    integer :: status                            ! error status
    integer :: dimid                             ! dimension ID
    integer :: k                                 ! 1d index
    integer :: ntapes_onfile                     ! number of history tapes on the restart file
    integer :: nflds_onfile                      ! number of history fields on the restart file
    integer :: t                                 ! tape index
    integer :: f                                 ! field index
    integer :: varid                             ! variable id
    integer, allocatable :: itemp2d(:)           ! 2D temporary
    real(r8), pointer :: hbuf(:,:)               ! history buffer
    real(r8), pointer :: hbuf1d(:)               ! 1d history buffer
    integer , pointer :: nacs(:,:)               ! accumulation counter
    integer , pointer :: nacs1d(:)               ! 1d accumulation counter
    integer           :: ier                     ! error code
    logical           :: readvar
    type(Var_desc_t)  :: vardesc                 ! netCDF variable description
    character(len=*),parameter :: subname = 'hist_restart_ncd'
    character(len=1) :: inst_suffix=' '
!------------------------------------------------------------------------

  numg=bounds%ngrid;numt=bounds%ntopou;numc=bounds%ncols;nump=bounds%npfts

  ! If branch run, initialize file times and return

  if (flag == 'read') then
    if (nsrest == nsrBranch) then
      do t = 1,ntapes
        tape(t)%ntimes = 0
      end do
      return
     end if
   ! If startup run just return
     if (nsrest == nsrStartup) then
       RETURN
     end if
   endif

    ! Read history file data only for restart run (not for branch run)

    !
    ! First when writing out and in define mode, create files and define all variables
    !
    !================================================
    if (flag == 'define') then
    !================================================

       if (.not. present(rdate)) then
          call endrun(msg=' variable rdate must be present for writing restart files'//&
               errMsg(__FILE__, __LINE__))
       end if

       !
       ! On master restart file add ntapes/max_chars dimension
       ! and then add the history and history restart filenames
       !
       call ncd_defdim( ncid, 'ntapes'       , ntapes      , dimid)
       call ncd_defdim( ncid, 'max_chars'    , max_chars   , dimid)

       call ncd_defvar(ncid=ncid, varname='locfnh', xtype=ncd_char, &
            long_name="History filename",     &
            comment="This variable NOT needed for startup or branch simulations", &
            dim1name='max_chars', dim2name="ntapes" )
       ier = ncd_inq_varid(ncid, 'locfnh', vardesc)
 !      ier = ncd_putatt(ncid, vardesc%varid, 'interpinic_flag', iflag_skip)

       call ncd_defvar(ncid=ncid, varname='locfnhr', xtype=ncd_char, &
            long_name="Restart history filename",     &
            comment="This variable NOT needed for startup or branch simulations", &
            dim1name='max_chars', dim2name="ntapes" )
       ier = ncd_inq_varid(ncid, 'locfnhr', vardesc)
!       ier = ncd_putatt(ncid, vardesc%varid, 'interpinic_flag', iflag_skip)

       ! max_nflds is the maximum number of fields on any tape
       ! max_flds is the maximum number possible number of fields

       max_nflds = max_nFields()

       ! Loop over tapes - write out namelist information to each restart-history tape
       ! only read/write accumulators and counters if needed

       do t = 1,ntapes

          ! Create the restart history filename and open it
          write(hnum,'(i1.1)') t-1
          locfnhr(t) = "./" // trim(case_name) //".ecosim"// trim(inst_suffix) &
                        // ".rh" // hnum //"."// trim(rdate) //".nc"

          call htape_create( t, histrest=.true. )

          ! Add read/write accumultators and counters if needed
          if (.not. tape(t)%is_endhist) then
             do f = 1,tape(t)%nflds
                name           =  tape(t)%hlist(f)%field%name
                long_name      =  tape(t)%hlist(f)%field%long_name
                standard_name      =  tape(t)%hlist(f)%field%standard_name
                units          =  tape(t)%hlist(f)%field%units
                name_acc       =  trim(name) // "_acc"
                units_acc      =  "unitless positive integer"
                long_name_acc  =  trim(long_name) // " accumulator number of samples"
                type1d_out     =  tape(t)%hlist(f)%field%type1d_out
                type2d         =  tape(t)%hlist(f)%field%type2d
                num2d          =  tape(t)%hlist(f)%field%num2d
                nacs           => tape(t)%hlist(f)%nacs
                hbuf           => tape(t)%hlist(f)%hbuf

                dim1name = type1d_out ; dim2name = 'undefined'

                if (dim2name == 'undefined') then
                   if (num2d == 1) then
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, &
                           dim1name=dim1name, &
                           long_name=trim(long_name), &
                           standard_name=standard_name, units=trim(units))
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                           dim1name=dim1name, &
                           long_name=trim(long_name_acc), units=trim(units_acc))
                   else
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, &
                           dim1name=dim1name, dim2name=type2d, &
                           long_name=trim(long_name), &
                           standard_name=standard_name, units=trim(units))
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                           dim1name=dim1name, dim2name=type2d, &
                           long_name=trim(long_name_acc), units=trim(units_acc))
                   end if
                else
                   if (num2d == 1) then
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, &
                           dim1name=dim1name, dim2name=dim2name, &
                           long_name=trim(long_name), &
                           standard_name=standard_name, units=trim(units))
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                           dim1name=dim1name, dim2name=dim2name, &
                           long_name=trim(long_name_acc), units=trim(units_acc))
                   else
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, &
                           dim1name=dim1name, dim2name=dim2name, dim3name=type2d, &
                           long_name=trim(long_name), &
                           standard_name=standard_name, units=trim(units))
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                           dim1name=dim1name, dim2name=dim2name, dim3name=type2d, &
                           long_name=trim(long_name_acc), units=trim(units_acc))
                   end if
                endif
             end do
          endif

          !
          ! Add namelist information to each restart history tape
          !
          call ncd_defdim( ncid_hist(t), 'fname_lenp2'  , max_namlen+2, dimid)
          call ncd_defdim( ncid_hist(t), 'fname_len'    , max_namlen  , dimid)
          call ncd_defdim( ncid_hist(t), 'len1'         , 1           , dimid)
          call ncd_defdim( ncid_hist(t), 'scalar'       , 1           , dimid)
          call ncd_defdim( ncid_hist(t), 'max_chars'    , max_chars   , dimid)
          call ncd_defdim( ncid_hist(t), 'max_nflds'    , max_nflds   ,  dimid)
          call ncd_defdim( ncid_hist(t), 'max_flds'     , max_flds    , dimid)
          call ncd_defdim( ncid_hist(t), 'string_length', 64          , dimid)
          call ncd_defvar(ncid=ncid_hist(t), varname='nhtfrq', xtype=ncd_int, &
               long_name="Frequency of history writes",               &
               comment="Namelist item", &
               units="absolute value of negative is in hours, 0=monthly, positive is time-steps",     &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='mfilt', xtype=ncd_int, &
               long_name="Number of history time samples on a file", units="1",     &
               comment="Namelist item", &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='ncprec', xtype=ncd_int, &
               long_name="Flag for data precision", flag_values=(/1,2/), &
               comment="Namelist item", &
               nvalid_range=(/1,2/), &
               flag_meanings=(/"single-precision", "double-precision"/), &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='fincl', xtype=ncd_char, &
               comment="Namelist item", &
               long_name="Fieldnames to include", &
               dim1name='fname_lenp2', dim2name='max_flds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='fexcl', xtype=ncd_char, &
               comment="Namelist item", &
               long_name="Fieldnames to exclude",  &
               dim1name='fname_lenp2', dim2name='max_flds' )

          call ncd_defvar(ncid=ncid_hist(t), varname='nflds', xtype=ncd_int, &
               long_name="Number of fields on file", units="1",        &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='ntimes', xtype=ncd_int, &
               long_name="Number of time steps on file", units="time-step",     &
               dim1name='scalar')

          call ncd_defvar(ncid=ncid_hist(t), varname='is_endhist', xtype=ncd_log, &
               long_name="End of history file", dim1name='scalar')

          call ncd_defvar(ncid=ncid_hist(t), varname='begtime', xtype=ncd_double, &
               long_name="Beginning time", units="time units",     &
               dim1name='scalar')

          call ncd_defvar(ncid=ncid_hist(t), varname='num2d', xtype=ncd_int, &
               long_name="Size of second dimension", units="1",     &
               dim1name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='hpindex', xtype=ncd_int, &
               long_name="History pointer index", units="1",     &
               dim1name='max_nflds' )

          call ncd_defvar(ncid=ncid_hist(t), varname='avgflag', xtype=ncd_char, &
               long_name="Averaging flag", &
               units="A=Average, X=Maximum, M=Minimum, I=Instantaneous", &
               dim1name='len1', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='name', xtype=ncd_char, &
               long_name="Fieldnames",  &
               dim1name='fname_len', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='long_name', xtype=ncd_char, &
               long_name="Long descriptive names for fields", &
               dim1name='max_chars', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='units', xtype=ncd_char, &
               long_name="Units for each history field output", &
               dim1name='max_chars', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='type1d', xtype=ncd_char, &
               long_name="1st dimension type", &
               dim1name='string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='type1d_out', xtype=ncd_char, &
               long_name="1st output dimension type", &
               dim1name='string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='type2d', xtype=ncd_char, &
               long_name="2nd dimension type", &
               dim1name='string_length', dim2name='max_nflds' )

          call ncd_enddef(ncid_hist(t))

       end do   ! end of ntapes loop

       RETURN

    !
    ! First write out namelist information to each restart history file
    !
    !================================================
    else if (flag == 'write') then
    !================================================

       ! Add history filenames to master restart file
       do t = 1,ntapes
          call ncd_putvar(ncid,'locfnh', t, locfnh(t))
          call ncd_putvar(ncid,'locfnhr', t, locfnhr(t))          
       end do

       fincl(:,1) = hist_fincl1(:)
       fincl(:,2) = hist_fincl2(:)
       fincl(:,3) = hist_fincl3(:)
       fincl(:,4) = hist_fincl4(:)
       fincl(:,5) = hist_fincl5(:)
       fincl(:,6) = hist_fincl6(:)

       fexcl(:,1) = hist_fexcl1(:)
       fexcl(:,2) = hist_fexcl2(:)
       fexcl(:,3) = hist_fexcl3(:)
       fexcl(:,4) = hist_fexcl4(:)
       fexcl(:,5) = hist_fexcl5(:)
       fexcl(:,6) = hist_fexcl6(:)

       max_nflds = max_nFields()

       start(1)=1

       allocate(itemp2d(max_nflds))

       !
       ! Add history namelist data to each history restart tape
       !
       
       do t = 1,ntapes
          
          call ncd_io(varname='fincl', data=fincl(:,t), ncid=ncid_hist(t), flag='write')

          call ncd_io(varname='fexcl', data=fexcl(:,t), ncid=ncid_hist(t), flag='write')

          call ncd_io(varname='is_endhist', data=tape(t)%is_endhist, ncid=ncid_hist(t), flag='write')

!          itemp2d(:,:) = 0
!          do f=1,tape(t)%nflds
!             itemp2d(f,t) = tape(t)%hlist(f)%field%num2d
!          end do
!          call ncd_io(varname='num2d', data=itemp2d(:,t), ncid=ncid_hist(t), flag='write')

!          itemp2d(:,:) = 0
!          do f=1,tape(t)%nflds
!             itemp2d(f,t) = tape(t)%hlist(f)%field%hpindex
!          end do
!          call ncd_io(varname='hpindex', data=itemp2d(:,t), ncid=ncid_hist(t), flag='write')

          itemp2d(:) = 0
          do f=1,tape(t)%nflds
             itemp2d(f) = tape(t)%hlist(f)%field%num2d
          end do
          call ncd_io(varname='num2d', data=itemp2d(:), ncid=ncid_hist(t), flag='write')

          itemp2d(:) = 0
          do f=1,tape(t)%nflds
             itemp2d(f) = tape(t)%hlist(f)%field%hpindex
          end do
          call ncd_io(varname='hpindex', data=itemp2d(:), ncid=ncid_hist(t), flag='write')

          call ncd_io('nflds',        tape(t)%nflds,   'write', ncid_hist(t) )
          call ncd_io('ntimes',       tape(t)%ntimes,  'write', ncid_hist(t) )
          call ncd_io('nhtfrq',  tape(t)%nhtfrq,  'write', ncid_hist(t) )
          call ncd_io('mfilt',   tape(t)%mfilt,   'write', ncid_hist(t) )
          call ncd_io('ncprec',  tape(t)%ncprec,  'write', ncid_hist(t) )
          call ncd_io('begtime',      tape(t)%begtime, 'write', ncid_hist(t) )
          allocate(tmpstr(tape(t)%nflds,7 ),tname(tape(t)%nflds), &
                   tavgflag(tape(t)%nflds),tunits(tape(t)%nflds),tlongname(tape(t)%nflds))
          do f=1,tape(t)%nflds
             tname(f)  = tape(t)%hlist(f)%field%name
             tunits(f) = tape(t)%hlist(f)%field%units
             tlongname(f) = tape(t)%hlist(f)%field%long_name
             tmpstr(f,1) = tape(t)%hlist(f)%field%type1d
             tmpstr(f,2) = tape(t)%hlist(f)%field%type1d_out
             tmpstr(f,3) = tape(t)%hlist(f)%field%type2d
             tavgflag(f) = tape(t)%hlist(f)%avgflag
          end do
          call ncd_io( 'name', tname, 'write',ncid_hist(t))
          call ncd_io('long_name', tlongname, 'write', ncid_hist(t))
          call ncd_io('units', tunits, 'write',ncid_hist(t))
          call ncd_io('type1d', tmpstr(:,1), 'write', ncid_hist(t))
          call ncd_io('type1d_out', tmpstr(:,2), 'write', ncid_hist(t))
          call ncd_io('type2d', tmpstr(:,3), 'write', ncid_hist(t))
          call ncd_io('avgflag',tavgflag , 'write', ncid_hist(t))
          deallocate(tname,tlongname,tunits,tmpstr,tavgflag)
       enddo
       deallocate(itemp2d)

    !
    ! Read in namelist information
    !
    !================================================
    else if (flag == 'read') then
    !================================================

       call ncd_inqdlen(ncid,dimid,ntapes_onfile, name='ntapes')
       
       if ( is_restart() .and. ntapes_onfile /= ntapes )then
          write(iulog,*) 'ntapes = ', ntapes, ' ntapes_onfile = ', ntapes_onfile
          call endrun(msg=' ERROR: number of ntapes different than on restart file!,'// &
               ' you can NOT change history options on restart!' //&
               errMsg(__FILE__, __LINE__))
       end if
       if ( is_restart() .and. ntapes > 0 )then
          do t = 1,ntapes
            call ncd_getvar(ncid,'locfnh', t, locfnh(t))
            call ncd_getvar(ncid,'locfnhr', t, locrest(t))
            call strip_null(locrest(t))
            call strip_null(locfnh(t))
          end do
       end if

       ! Determine necessary indices - the following is needed if model decomposition is different on restart

       start(1)=1

       if ( is_restart() )then
          D100: do t = 1,ntapes

             call getfil( locrest(t), locfnhr(t), 0 )

             call ncd_pio_openfile (ncid_hist(t), trim(locfnhr(t)), ncd_nowrite)
             
             if ( t == 1 )then

                call ncd_inqdlen(ncid_hist(1),dimid,max_nflds,name='max_nflds')

                allocate(itemp2d(max_nflds))
             end if

             call ncd_inqvid(ncid_hist(t), 'name',           varid, name_desc)
             call ncd_inqvid(ncid_hist(t), 'long_name',      varid, longname_desc)
             call ncd_inqvid(ncid_hist(t), 'units',          varid, units_desc)
             call ncd_inqvid(ncid_hist(t), 'type1d',         varid, type1d_desc)
             call ncd_inqvid(ncid_hist(t), 'type1d_out',     varid, type1d_out_desc)
             call ncd_inqvid(ncid_hist(t), 'type2d',         varid, type2d_desc)
             call ncd_inqvid(ncid_hist(t), 'avgflag',        varid, avgflag_desc)

             call ncd_io(varname='fincl', data=fincl(:,t), ncid=ncid_hist(t), flag='read')

             call ncd_io(varname='fexcl', data=fexcl(:,t), ncid=ncid_hist(t), flag='read')

             call ncd_io('nflds',   nflds_onfile, 'read', ncid_hist(t) )

             if ( nflds_onfile /= tape(t)%nflds )then
                write(iulog,*) 'nflds = ', tape(t)%nflds, ' nflds_onfile = ', nflds_onfile
                call endrun(msg=' ERROR: number of fields different than on restart file!,'// &
                     ' you can NOT change history options on restart!' //&
                     errMsg(__FILE__, __LINE__))
             end if
             call ncd_io('ntimes',  tape(t)%ntimes, 'read', ncid_hist(t) )
             call ncd_io('nhtfrq',  tape(t)%nhtfrq, 'read', ncid_hist(t) )
             call ncd_io('mfilt',   tape(t)%mfilt, 'read', ncid_hist(t) )
             call ncd_io('ncprec',  tape(t)%ncprec, 'read', ncid_hist(t) )
             call ncd_io('begtime', tape(t)%begtime, 'read', ncid_hist(t) )

             call ncd_io(varname='is_endhist', data=tape(t)%is_endhist, ncid=ncid_hist(t), flag='read')
             
             call ncd_io(varname='num2d', data=itemp2d(:), ncid=ncid_hist(t), flag='read')
             do f=1,tape(t)%nflds
                tape(t)%hlist(f)%field%num2d = itemp2d(f)
             end do

             call ncd_io(varname='hpindex', data=itemp2d(:), ncid=ncid_hist(t), flag='read')
             do f=1,tape(t)%nflds
                tape(t)%hlist(f)%field%hpindex = itemp2d(f)
             end do

             D101: do f=1,tape(t)%nflds
                start(2) = f
                
                call ncd_io( name_desc,           tape(t)%hlist(f)%field%name,       &
                             'read', ncid_hist(t), start )
                call ncd_io( longname_desc,       tape(t)%hlist(f)%field%long_name,  &
                             'read', ncid_hist(t), start )
                call ncd_io( units_desc,          tape(t)%hlist(f)%field%units,      &
                             'read', ncid_hist(t), start )
                call ncd_io( type1d_desc,         tape(t)%hlist(f)%field%type1d,     &
                             'read', ncid_hist(t), start )
                call ncd_io( type1d_out_desc,     tape(t)%hlist(f)%field%type1d_out, &
                             'read', ncid_hist(t), start )
                call ncd_io( type2d_desc,         tape(t)%hlist(f)%field%type2d,     &
                             'read', ncid_hist(t), start )
                call ncd_io( avgflag_desc,        tape(t)%hlist(f)%avgflag,          &
                             'read', ncid_hist(t), start )
                call strip_null(tape(t)%hlist(f)%field%name)
                call strip_null(tape(t)%hlist(f)%field%long_name)
                call strip_null(tape(t)%hlist(f)%field%units)
                call strip_null(tape(t)%hlist(f)%field%type1d)
                call strip_null(tape(t)%hlist(f)%field%type1d_out)
                call strip_null(tape(t)%hlist(f)%field%type2d)
                call strip_null(tape(t)%hlist(f)%avgflag)

                type1d_out = trim(tape(t)%hlist(f)%field%type1d_out)
                select case (trim(type1d_out))
                case (nameg)
                   num1d_out = numg
                   beg1d_out = bounds%begg
                   end1d_out = bounds%endg
                case (namet)
                   num1d_out = numt
                   beg1d_out = bounds%begt
                   end1d_out = bounds%endt
                case (namec)
                   num1d_out = numc
                   beg1d_out = bounds%begc
                   end1d_out = bounds%endc
                case (namep)
                   num1d_out = nump
                   beg1d_out = bounds%begp
                   end1d_out = bounds%endp
                case default
                   write(iulog,*) trim(subname),' ERROR: read unknown 1d output type=',trim(type1d_out)
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                end select

                tape(t)%hlist(f)%field%num1d_out = num1d_out
                tape(t)%hlist(f)%field%beg1d_out = beg1d_out
                tape(t)%hlist(f)%field%end1d_out = end1d_out

                num2d  = tape(t)%hlist(f)%field%num2d
                allocate (tape(t)%hlist(f)%hbuf(beg1d_out:end1d_out,num2d), &
                          tape(t)%hlist(f)%nacs(beg1d_out:end1d_out,num2d), &
                          stat=status)
                if (status /= 0) then
                   write(iulog,*) trim(subname),' ERROR: allocation error for hbuf,nacs at t,f=',t,f
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                endif
                tape(t)%hlist(f)%hbuf(:,:) = 0._r8
                tape(t)%hlist(f)%nacs(:,:) = 0

                type1d = tape(t)%hlist(f)%field%type1d
                select case (type1d)
                case (nameg)
                   num1d = numg
                   beg1d = bounds%begg
                   end1d = bounds%endg
                case (namet)
                   num1d = numt
                   beg1d = bounds%begt
                   end1d = bounds%endt
                case (namec)
                   num1d = numc
                   beg1d = bounds%begc
                   end1d = bounds%endc
                case (namep)
                   num1d = nump
                   beg1d = bounds%begp
                   end1d = bounds%endp
                case default
                   write(iulog,*) trim(subname),' ERROR: read unknown 1d type=',type1d
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                end select

                tape(t)%hlist(f)%field%num1d = num1d
                tape(t)%hlist(f)%field%beg1d = beg1d
                tape(t)%hlist(f)%field%end1d = end1d

             end do D101   ! end of flds loop

             ! If history file is not full, open it
             
             if (tape(t)%ntimes /= 0) then
                call ncd_pio_openfile (nfid(t), trim(locfnh(t)), ncd_write)
             end if

          end do  D100 ! end of tapes loop
          
          hist_fincl1(:) = fincl(:,1)
          hist_fincl2(:) = fincl(:,2)
          hist_fincl3(:) = fincl(:,3)
          hist_fincl4(:) = fincl(:,4)
          hist_fincl5(:) = fincl(:,5)
          hist_fincl6(:) = fincl(:,6)

          hist_fexcl1(:) = fexcl(:,1)
          hist_fexcl2(:) = fexcl(:,2)
          hist_fexcl3(:) = fexcl(:,3)
          hist_fexcl4(:) = fexcl(:,4)
          hist_fexcl5(:) = fexcl(:,5)
          hist_fexcl6(:) = fexcl(:,6)

       end if

       if ( allocated(itemp2d) ) deallocate(itemp2d)

    end if

    !======================================================================
    !  Read/write history file restart data.
    ! If the current history file(s) are not full, file(s) are opened
    ! so that subsequent time samples are added until the file is full.
    ! A new history file is used on a branch run.
    !======================================================================

    if (flag == 'write') then

       do t = 1,ntapes
          if (.not. tape(t)%is_endhist) then

             do f = 1,tape(t)%nflds
                name       =  tape(t)%hlist(f)%field%name
                name_acc   =  trim(name) // "_acc"
                type1d_out =  tape(t)%hlist(f)%field%type1d_out
                type2d     =  tape(t)%hlist(f)%field%type2d
                num2d      =  tape(t)%hlist(f)%field%num2d
                beg1d_out  =  tape(t)%hlist(f)%field%beg1d_out
                end1d_out  =  tape(t)%hlist(f)%field%end1d_out
                nacs       => tape(t)%hlist(f)%nacs
                hbuf       => tape(t)%hlist(f)%hbuf

                if (num2d == 1) then
                   allocate(hbuf1d(beg1d_out:end1d_out), &
                            nacs1d(beg1d_out:end1d_out), stat=status)
                   if (status /= 0) then
                      write(iulog,*) trim(subname),' ERROR: allocation'
                      call endrun(msg=errMsg(__FILE__, __LINE__))
                   end if

                   hbuf1d(beg1d_out:end1d_out) = hbuf(beg1d_out:end1d_out,1)
                   nacs1d(beg1d_out:end1d_out) = nacs(beg1d_out:end1d_out,1)
                   
                   call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name), &
                        dim1name=type1d_out, data=hbuf1d)
                   call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name_acc), &
                        dim1name=type1d_out, data=nacs1d)

                   deallocate(hbuf1d)
                   deallocate(nacs1d)
                else
                   call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name), &
                        dim1name=type1d_out, data=hbuf)
                   call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name_acc), &
                        dim1name=type1d_out, data=nacs)
                end if

             end do

          end if  ! end of is_endhist block

          call ncd_pio_closefile(ncid_hist(t))

       end do   ! end of ntapes loop

    else if (flag == 'read') then

       do t = 1,ntapes

          if (.not. tape(t)%is_endhist) then

             do f = 1,tape(t)%nflds
                name       =  tape(t)%hlist(f)%field%name
                name_acc   =  trim(name) // "_acc"
                type1d_out =  tape(t)%hlist(f)%field%type1d_out
                type2d     =  tape(t)%hlist(f)%field%type2d
                num2d      =  tape(t)%hlist(f)%field%num2d
                beg1d_out  =  tape(t)%hlist(f)%field%beg1d_out
                end1d_out  =  tape(t)%hlist(f)%field%end1d_out
                nacs       => tape(t)%hlist(f)%nacs
                hbuf       => tape(t)%hlist(f)%hbuf
                
                if (num2d == 1) then
                   allocate(hbuf1d(beg1d_out:end1d_out), &
                        nacs1d(beg1d_out:end1d_out), stat=status)
                   if (status /= 0) then
                      write(iulog,*) trim(subname),' ERROR: allocation'
                      call endrun(msg=errMsg(__FILE__, __LINE__))
                   end if
                
                   call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name), &
                        dim1name=type1d_out, data=hbuf1d)
                   call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name_acc), &
                        dim1name=type1d_out, data=nacs1d)

                   hbuf(beg1d_out:end1d_out,1) = hbuf1d(beg1d_out:end1d_out)
                   nacs(beg1d_out:end1d_out,1) = nacs1d(beg1d_out:end1d_out)

                   deallocate(hbuf1d)
                   deallocate(nacs1d)
                else
                   call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name), &
                        dim1name=type1d_out, data=hbuf)
                   call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name_acc), &
                        dim1name=type1d_out, data=nacs)
                end if
             end do

          end if

          call ncd_pio_closefile(ncid_hist(t))

       end do

    end if

  end subroutine hist_restart_ncd

  !-----------------------------------------------------------------------
  integer function max_nFields()
    !
    ! !DESCRIPTION:
    ! Get the maximum number of fields on all tapes.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: t  ! index
    character(len=*),parameter :: subname = 'max_nFields'
    !-----------------------------------------------------------------------

    max_nFields = 0
    do t = 1,ntapes
       max_nFields = max(max_nFields, tape(t)%nflds)
    end do
    return
  end function max_nFields

  !-----------------------------------------------------------------------
  subroutine hfields_normalize (t)
    !
    ! !DESCRIPTION:
    ! Normalize fields on a history file by the number of accumulations.
    ! Loop over fields on the tape.  Need averaging flag and number of
    ! accumulations to perform normalization.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t       ! tape index
    !
    ! !LOCAL VARIABLES:
    integer :: f                   ! field index
    integer :: k                   ! 1d index
    integer :: j                   ! 2d index
    logical :: aflag               ! averaging flag
    integer :: beg1d_out,end1d_out ! hbuf 1d beginning and ending indices
    integer :: num2d               ! hbuf size of second dimension (e.g. number of vertical levels)
    character(len=1)  :: avgflag   ! averaging flag
    real(r8), pointer :: hbuf(:,:) ! history buffer
    integer , pointer :: nacs(:,:) ! accumulation counter
    character(len=*),parameter :: subname = 'hfields_normalize'
    !-----------------------------------------------------------------------

    ! Normalize by number of accumulations for time averaged case

    do f = 1,tape(t)%nflds
       avgflag   =  tape(t)%hlist(f)%avgflag
       beg1d_out =  tape(t)%hlist(f)%field%beg1d_out
       end1d_out =  tape(t)%hlist(f)%field%end1d_out
       num2d     =  tape(t)%hlist(f)%field%num2d
       nacs      => tape(t)%hlist(f)%nacs
       hbuf      => tape(t)%hlist(f)%hbuf

       if (avgflag == 'A' .or. avgflag == 'P') then
          aflag = .true.
       else
          aflag = .false.
       end if

       do j = 1, num2d
          do k = beg1d_out, end1d_out
             if (aflag .and. nacs(k,j) /= 0) then
                hbuf(k,j) = hbuf(k,j) / float(nacs(k,j))
             end if
          end do
       end do
    end do

  end subroutine hfields_normalize

  !-----------------------------------------------------------------------
  subroutine hfields_zero (t)
    !
    ! !DESCRIPTION:
    ! Zero out accumulation and history buffers for a given history tape.
    ! Loop through fields on the tape.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t     ! tape index
    !
    ! !LOCAL VARIABLES:
    integer :: f                 ! field index
    character(len=*),parameter :: subname = 'hfields_zero'
    !-----------------------------------------------------------------------

    do f = 1,tape(t)%nflds
       tape(t)%hlist(f)%hbuf(:,:) = 0._r8
       tape(t)%hlist(f)%nacs(:,:) = 0
    end do

  end subroutine hfields_zero

end module HistFileMod
