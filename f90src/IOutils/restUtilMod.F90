module restUtilMod
   use EcoSIMConfig
   use ElmIDMod
   use ncdio_pio
   use netcdf
   use data_kind_mod, only : r8 => DAT_Kind_r8
   use data_const_mod,only : spval => DAT_CONST_SPVAL, ispval => DAT_CONST_ISPVAL
   use abortutils   , only : endrun
   use GridMod
  implicit none
  private
  character(len=*), parameter :: mod_filename=&
  __FILE__
  public :: restartvar
  public :: cppft
  public :: cpcol

  interface cpcol
    module procedure cpcol_i_1d,cpcol_r_1d,cpcol_r_2d
    module procedure cpcol_r_3d,cpcol_r_4d,cpcol_r_5d
  end interface cpcol

  interface cppft
    module procedure cppft_i_1d,cppft_i_2d,cppft_i_3d,cppft_r_1d
    module procedure cppft_r_2d,cppft_r_3d,cppft_r_4d,cppft_r_5d
  end interface cppft

  interface restartvar
    module procedure restartvar_int
    module procedure restartvar_int_1d
    module procedure restartvar_int_2d
    module procedure restartvar_int_3d
    module procedure restartvar_real_sp
    module procedure restartvar_real_sp_1d
    module procedure restartvar_real_sp_2d
    module procedure restartvar_real_sp_3d
    module procedure restartvar_real_sp_4d
    module procedure restartvar_real_sp_5d
    module procedure restartvar_logical_1d    
  end interface restartvar

  integer,parameter, public :: iflag_interp = 1
  integer,parameter, public :: iflag_copy   = 2
  integer,parameter, public :: iflag_skip   = 3

  contains

  subroutine restartvar_int(ncid, flag, varname, &
       long_name, units, interpinic_flag, data, &
       comment, flag_meanings, missing_value, fill_value, &
       flag_values, nvalid_range )

    !----------------------------------------------------
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    integer           , intent(inout)        :: data
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    integer           , intent(in), optional :: missing_value    ! attribute for real
    integer           , intent(in), optional :: fill_value       ! attribute for real
    integer           , intent(in), optional :: flag_values(:)   ! attribute for int
    integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int
    !
    ! Local variables
    character(len=*), parameter :: sub=trim(mod_filename)//'::'//'restartvar_int'
    logical           :: readvar          ! was var read?
    integer          :: ivalue
    type(var_desc_t) :: vardesc  ! local vardesc
    integer          :: status   ! return error code
    integer          :: varid
    !----------------------------------------------------

    readvar = .false.
    if (flag == 'define') then

       call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=ncd_int, &
            long_name=trim(long_name), units=units)

       status = nf90_inq_varid(ncid%fh,trim(varname),varid)

       if (trim(interpinic_flag) == 'interp') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_interp)
       else if (trim(interpinic_flag) == 'copy') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_copy)
       else if (trim(interpinic_flag) == 'skip') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_skip)
       end if

       status = ncd_putatt(ncid, varid, 'interpinic_flag_meanings', &
            "1=nearest neighbor, 2=copy directly, 3=skip")

       if (present(comment)) then
          call check_ret(ncd_putatt(ncid, varid, 'comment', trim(comment)),sub)
       end if
       if (present(units)) then
          call check_ret(ncd_putatt(ncid, varid, 'units', trim(units)),sub)
       end if
       if (present(fill_value)) then
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', fill_value),sub)
       end if
       if (present(missing_value)) then
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', missing_value),sub)
       end if

    else if (flag == 'read' .or. flag == 'write') then

       call ncd_io(varname=trim(varname), data=data, ncid=ncid, flag=flag, readvar=readvar)

    end if

    if (flag == 'read') then
       if (.not. readvar .and. is_restart()) call endrun('Reading restart file failed in '//trim(sub),__LINE__)
    end if

  end subroutine restartvar_int

  !------------------------------------------------------------------------
  subroutine restartvar_real_sp(ncid, flag, varname, &
       long_name, units, interpinic_flag, data, &
       comment, flag_meanings, missing_value, fill_value)

    !----------------------------------------------------
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    real(r8)          , intent(inout)        :: data

    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)          , intent(in), optional :: missing_value    ! attribute for real
    real(r8)          , intent(in), optional :: fill_value       ! attribute for real
    !
    ! Local variables
    character(len=*), parameter :: sub=trim(mod_filename)//'::'//'restartvar_real_sp'
    logical           :: readvar          ! was var read?
    integer          :: ivalue
    type(var_desc_t) :: vardesc  ! local vardesc
    integer          :: status   ! return error code
    integer          :: varid
    !----------------------------------------------------

    readvar = .false.
    if (flag == 'define') then

       call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=ncd_double, &
            long_name=trim(long_name), units=units)

       status = nf90_inq_varid(ncid%fh,trim(varname),varid)

       if (trim(interpinic_flag) == 'interp') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_interp)
       else if (trim(interpinic_flag) == 'copy') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_copy)
       else if (trim(interpinic_flag) == 'skip') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_skip)
       end if
       status = nf90_put_att(ncid%fh, varid, 'interpinic_flag_meanings', &
            "1=nearest neighbor, 2=copy directly, 3=skip")

       if (present(comment)) then
          call check_ret(ncd_putatt(ncid, varid, 'comment', trim(comment)),sub)
       end if
       if (present(units)) then
          call check_ret(ncd_putatt(ncid, varid, 'units', trim(units)),sub)
       end if
       if (present(fill_value)) then
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', fill_value,ncd_double),sub)
       end if
       if (present(missing_value)) then
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', missing_value,ncd_double),sub)
       end if

    else if (flag == 'read' .or. flag == 'write') then

       call ncd_io(varname=trim(varname), data=data, ncid=ncid, flag=flag, readvar=readvar)

    end if

    if (flag == 'read') then

       if (.not. readvar .and. is_restart()) &
         call endrun('Reading '//trim(varname)//' from restart file failed in '//trim(sub),__LINE__)
    end if

  end subroutine restartvar_real_sp
!-----------------------------------------------------------------------
  subroutine restartvar_int_1d(ncid, flag, varname, dim1name,  &
       long_name, units, interpinic_flag, data, &
       comment, flag_meanings, missing_value, fill_value, &
       flag_values, nvalid_range )

   implicit none
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    integer           , pointer              :: data(:)
    character(len=*)  , intent(in)           :: dim1name         ! dimension name
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    integer           , intent(in), optional :: missing_value   ! attribute for int
    integer           , intent(in), optional :: fill_value      ! attribute for int
    integer           , intent(in), optional :: flag_values(:)   ! attribute for int
    integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int


    ! Local variables
    ! Local variables
    character(len=*), parameter :: sub=trim(mod_filename)//'::'//'restartvar_int_arr'
    logical           :: readvar          ! was var read?
    integer          :: ivalue
    type(var_desc_t) :: vardesc  ! local vardesc
    integer          :: status   ! return error code
    integer          :: varid
    integer          :: lxtype   ! local external type (in case logical variable)
    !----------------------------------------------------

    readvar = .false.
    if (flag == 'define') then

      lxtype = ncd_int

      call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=lxtype, &
         dim1name=trim(dim1name), &
         long_name=trim(long_name), units=units)

      status = nf90_inq_varid(ncid%fh,trim(varname),varid)

       if (trim(interpinic_flag) == 'interp') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_interp)
       else if (trim(interpinic_flag) == 'copy') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_copy)
       else if (trim(interpinic_flag) == 'skip') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_skip)
       end if

       status = nf90_put_att(ncid%fh, varid, 'interpinic_flag_meanings', &
            "1=nearest neighbor, 2=copy directly, 3=skip")

       if (present(comment)) then
          call check_ret(ncd_putatt(ncid, varid, 'comment', trim(comment)),sub)
       end if
       if (present(units)) then
          call check_ret(ncd_putatt(ncid, varid, 'units', trim(units)),sub)
       end if

       if (present(fill_value)) then
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', fill_value),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', ispval),sub)
       end if
       if (present(missing_value)) then
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', missing_value),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', ispval),sub)
       end if

       if (present(nvalid_range)) then
          status = ncd_putatt(ncid,varid,'valid_range', nvalid_range )
       end if

    else if (flag == 'read' .or. flag == 'write') then

      call ncd_io(varname=trim(varname), data=data, &
        dim1name=trim(dim1name), ncid=ncid, flag=flag, readvar=readvar)
    end if

    if (flag == 'read') then
       if (.not. readvar .and. is_restart()) &
         call endrun('Reading '//trim(varname)//' from restart file failed in '//trim(sub),__LINE__)
    end if

  end subroutine restartvar_int_1d
!-----------------------------------------------------------------------
  subroutine restartvar_logical_1d(ncid, flag, varname, dim1name,  &
       long_name, units, interpinic_flag, data, &
       comment, flag_meanings, missing_value, fill_value, &
       flag_values, nvalid_range )

   implicit none
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    logical           , pointer              :: data(:)
    character(len=*)  , intent(in)           :: dim1name         ! dimension name
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    integer           , intent(in), optional :: missing_value   ! attribute for int
    integer           , intent(in), optional :: fill_value      ! attribute for int
    integer           , intent(in), optional :: flag_values(:)   ! attribute for int
    integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int


    ! Local variables
    ! Local variables
    character(len=*), parameter :: sub=trim(mod_filename)//'::'//'restartvar_int_arr'
    logical           :: readvar          ! was var read?
    integer          :: ivalue
    type(var_desc_t) :: vardesc  ! local vardesc
    integer          :: status   ! return error code
    integer          :: varid
    integer          :: lxtype   ! local external type (in case logical variable)
    !----------------------------------------------------

    readvar = .false.
    if (flag == 'define') then

      lxtype = ncd_int

      call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=lxtype, &
         dim1name=trim(dim1name), &
         long_name=trim(long_name), units=units)

      status = nf90_inq_varid(ncid%fh,trim(varname),varid)

       if (trim(interpinic_flag) == 'interp') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_interp)
       else if (trim(interpinic_flag) == 'copy') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_copy)
       else if (trim(interpinic_flag) == 'skip') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_skip)
       end if

       status = nf90_put_att(ncid%fh, varid, 'interpinic_flag_meanings', &
            "1=nearest neighbor, 2=copy directly, 3=skip")

       if (present(comment)) then
          call check_ret(ncd_putatt(ncid, varid, 'comment', trim(comment)),sub)
       end if
       if (present(units)) then
          call check_ret(ncd_putatt(ncid, varid, 'units', trim(units)),sub)
       end if

       if (present(fill_value)) then
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', fill_value),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', ispval),sub)
       end if
       if (present(missing_value)) then
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', missing_value),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', ispval),sub)
       end if

       if (present(nvalid_range)) then
          status = ncd_putatt(ncid,varid,'valid_range', nvalid_range )
       end if

    else if (flag == 'read' .or. flag == 'write') then

      call ncd_io(varname=trim(varname), data=data, &
        dim1name=trim(dim1name), ncid=ncid, flag=flag, readvar=readvar)
    end if

    if (flag == 'read') then
       if (.not. readvar .and. is_restart()) &
         call endrun('Reading '//trim(varname)//' from restart file failed in '//trim(sub),__LINE__)
    end if

  end subroutine restartvar_logical_1d  
!-----------------------------------------------------------------------
  subroutine restartvar_int_2d(ncid, flag, varname, dim1name, dim2name, &
       long_name, units, interpinic_flag, data, &
       comment, flag_meanings, missing_value, fill_value, &
       flag_values, nvalid_range )

   implicit none
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    integer           , pointer              :: data(:,:)
    character(len=*)  , intent(in)           :: dim1name         ! dimension name
    character(len=*)  , intent(in)           :: dim2name         ! dimension name
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    integer           , intent(in), optional :: missing_value   ! attribute for int
    integer           , intent(in), optional :: fill_value      ! attribute for int
    integer           , intent(in), optional :: flag_values(:)   ! attribute for int
    integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int


    ! Local variables
    ! Local variables
    character(len=*), parameter :: sub=trim(mod_filename)//'::'//'restartvar_int_2d'
    logical           :: readvar          ! was var read?
    integer          :: ivalue
    type(var_desc_t) :: vardesc  ! local vardesc
    integer          :: status   ! return error code
    integer          :: varid
    integer          :: lxtype   ! local external type (in case logical variable)
    !----------------------------------------------------

    readvar = .false.
    if (flag == 'define') then

      lxtype = ncd_int

      call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=lxtype, &
         dim1name=trim(dim1name), dim2name=trim(dim2name),&
         long_name=trim(long_name), units=units)

      status = nf90_inq_varid(ncid%fh,trim(varname),varid)

       if (trim(interpinic_flag) == 'interp') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_interp)
       else if (trim(interpinic_flag) == 'copy') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_copy)
       else if (trim(interpinic_flag) == 'skip') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_skip)
       end if

       status = nf90_put_att(ncid%fh, varid, 'interpinic_flag_meanings', &
            "1=nearest neighbor, 2=copy directly, 3=skip")

       if (present(comment)) then
          call check_ret(ncd_putatt(ncid, varid, 'comment', trim(comment)),sub)
       end if
       if (present(units)) then
          call check_ret(ncd_putatt(ncid, varid, 'units', trim(units)),sub)
       end if

       if (present(fill_value)) then
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', fill_value),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', ispval),sub)
       end if
       if (present(missing_value)) then
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', missing_value),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', ispval),sub)
       end if

       if (present(nvalid_range)) then
          status = ncd_putatt(ncid,varid,'valid_range', nvalid_range )
       end if

    else if (flag == 'read' .or. flag == 'write') then

      call ncd_io(varname=trim(varname), data=data, &
        dim1name=trim(dim1name), ncid=ncid, flag=flag, readvar=readvar)
    end if

    if (flag == 'read') then
       if (.not. readvar .and. is_restart()) &
         call endrun('Reading '//trim(varname)//' from restart file failed in '//trim(sub),__LINE__)
    end if

  end subroutine restartvar_int_2d
!-----------------------------------------------------------------------
  subroutine restartvar_int_3d(ncid, flag, varname, dim1name, dim2name, dim3name, &
       long_name, units, interpinic_flag, data, &
       comment, flag_meanings, missing_value, fill_value, &
       flag_values, nvalid_range )

   implicit none
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    integer           , pointer              :: data(:,:,:)
    character(len=*)  , intent(in)           :: dim1name         ! dimension name
    character(len=*)  , intent(in)           :: dim2name         ! dimension name
    character(len=*)  , intent(in)           :: dim3name         ! dimension name
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    integer           , intent(in), optional :: missing_value   ! attribute for int
    integer           , intent(in), optional :: fill_value      ! attribute for int
    integer           , intent(in), optional :: flag_values(:)   ! attribute for int
    integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int


    ! Local variables
    ! Local variables
    character(len=*), parameter :: sub=trim(mod_filename)//'::'//'restartvar_int_3d'
    logical           :: readvar          ! was var read?
    integer          :: ivalue
    type(var_desc_t) :: vardesc  ! local vardesc
    integer          :: status   ! return error code
    integer          :: varid
    integer          :: lxtype   ! local external type (in case logical variable)
    !----------------------------------------------------

    readvar = .false.
    if (flag == 'define') then

      lxtype = ncd_int

      call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=lxtype, &
         dim1name=trim(dim1name), dim2name=trim(dim2name), dim3name=trim(dim3name), &
         long_name=trim(long_name), units=units)

      status = nf90_inq_varid(ncid%fh,trim(varname),varid)

       if (trim(interpinic_flag) == 'interp') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_interp)
       else if (trim(interpinic_flag) == 'copy') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_copy)
       else if (trim(interpinic_flag) == 'skip') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_skip)
       end if

       status = nf90_put_att(ncid%fh, varid, 'interpinic_flag_meanings', &
            "1=nearest neighbor, 2=copy directly, 3=skip")

       if (present(comment)) then
          call check_ret(ncd_putatt(ncid, varid, 'comment', trim(comment)),sub)
       end if
       if (present(units)) then
          call check_ret(ncd_putatt(ncid, varid, 'units', trim(units)),sub)
       end if

       if (present(fill_value)) then
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', fill_value),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', ispval),sub)
       end if
       if (present(missing_value)) then
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', missing_value),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', ispval),sub)
       end if

       if (present(nvalid_range)) then
          status = ncd_putatt(ncid,varid,'valid_range', nvalid_range )
       end if

    else if (flag == 'read' .or. flag == 'write') then

      call ncd_io(varname=trim(varname), data=data, &
        dim1name=trim(dim1name), ncid=ncid, flag=flag, readvar=readvar)
    end if

    if (flag == 'read') then
       if (.not. readvar .and. is_restart()) &
         call endrun('Reading '//trim(varname)//' from restart file failed in '//trim(sub),__LINE__)
    end if

  end subroutine restartvar_int_3d
!-----------------------------------------------------------------------
  subroutine restartvar_real_sp_1d(ncid, flag, varname,  dim1name,  &
       long_name, units, interpinic_flag, data, &
       comment, flag_meanings, missing_value, fill_value)

   implicit none
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    real(r8)          , pointer              :: data(:)
    character(len=*)  , intent(in)           :: dim1name         ! dimension name
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)           , intent(in), optional :: missing_value   ! attribute for real
    real(r8)           , intent(in), optional :: fill_value      ! attribute for real


    ! Local variables
    ! Local variables
    character(len=*), parameter :: sub=trim(mod_filename)//'::'//'restartvar_real_sp_1d'
    logical           :: readvar          ! was var read?
    integer          :: ivalue
    type(var_desc_t) :: vardesc  ! local vardesc
    integer          :: status   ! return error code
    integer          :: varid
    integer          :: lxtype   ! local external type (in case logical variable)
    !----------------------------------------------------

    readvar = .false.
    if (flag == 'define') then

      lxtype = ncd_double

      call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=lxtype, &
         dim1name=trim(dim1name), &
         long_name=trim(long_name), units=units)

      status = nf90_inq_varid(ncid%fh,trim(varname),varid)

       if (trim(interpinic_flag) == 'interp') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_interp)
       else if (trim(interpinic_flag) == 'copy') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_copy)
       else if (trim(interpinic_flag) == 'skip') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_skip)
       end if

       status = nf90_put_att(ncid%fh, varid, 'interpinic_flag_meanings', &
            "1=nearest neighbor, 2=copy directly, 3=skip")

       if (present(comment)) then
          call check_ret(ncd_putatt(ncid, varid, 'comment', trim(comment)),sub)
       end if
       if (present(units)) then
          call check_ret(ncd_putatt(ncid, varid, 'units', trim(units)),sub)
       end if

       if (present(fill_value)) then
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', fill_value,ncd_double),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', spval,ncd_double),sub)
       end if
       if (present(missing_value)) then
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', missing_value,ncd_double),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', spval,ncd_double),sub)
       end if

    else if (flag == 'read' .or. flag == 'write') then

      call ncd_io(varname=trim(varname), data=data, &
        dim1name=trim(dim1name), ncid=ncid, flag=flag, readvar=readvar)
    end if

    if (flag == 'read') then
       if (.not. readvar .and. is_restart()) &
         call endrun('Reading restart file failed in '//trim(sub),__LINE__)
    end if

  end subroutine restartvar_real_sp_1d

!-----------------------------------------------------------------------
  subroutine restartvar_real_sp_2d(ncid, flag, varname,  dim1name,  dim2name, &
       long_name, units, interpinic_flag, data, &
       comment, flag_meanings, missing_value, fill_value)

   implicit none
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    real(r8)          , pointer              :: data(:,:)
    character(len=*)  , intent(in)           :: dim1name         ! dimension name
    character(len=*)  , intent(in)           :: dim2name         ! dimension name
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)           , intent(in), optional :: missing_value   ! attribute for real
    real(r8)           , intent(in), optional :: fill_value      ! attribute for real


    ! Local variables
    ! Local variables
    character(len=*), parameter :: sub=trim(mod_filename)//'::'//'restartvar_real_sp_2d'
    logical           :: readvar          ! was var read?
    integer          :: ivalue
    type(var_desc_t) :: vardesc  ! local vardesc
    integer          :: status   ! return error code
    integer          :: varid
    integer          :: lxtype   ! local external type (in case logical variable)
    !----------------------------------------------------

    readvar = .false.
    if (flag == 'define') then

      lxtype = ncd_double

      call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=lxtype, &
         dim1name=trim(dim1name), dim2name=trim(dim2name),&
         long_name=trim(long_name), units=units)

      status = nf90_inq_varid(ncid%fh,trim(varname),varid)

       if (trim(interpinic_flag) == 'interp') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_interp)
       else if (trim(interpinic_flag) == 'copy') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_copy)
       else if (trim(interpinic_flag) == 'skip') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_skip)
       end if

       status = nf90_put_att(ncid%fh, varid, 'interpinic_flag_meanings', &
            "1=nearest neighbor, 2=copy directly, 3=skip")

       if (present(comment)) then
          call check_ret(ncd_putatt(ncid, varid, 'comment', trim(comment)),sub)
       end if
       if (present(units)) then
          call check_ret(ncd_putatt(ncid, varid, 'units', trim(units)),sub)
       end if

       if (present(fill_value)) then
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', fill_value,ncd_double),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', spval,ncd_double),sub)
       end if
       if (present(missing_value)) then
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', missing_value,ncd_double),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', spval,ncd_double),sub)
       end if

    else if (flag == 'read' .or. flag == 'write') then

      call ncd_io(varname=trim(varname), data=data, &
        dim1name=trim(dim1name), ncid=ncid, flag=flag, readvar=readvar)
    end if

    if (flag == 'read') then
       if (.not. readvar .and. is_restart()) &
         call endrun('Reading restart file failed in '//trim(sub),__LINE__)
    end if

  end subroutine restartvar_real_sp_2d

!-----------------------------------------------------------------------
  subroutine restartvar_real_sp_3d(ncid, flag, varname,  dim1name,  dim2name, &
       dim3name, long_name, units, interpinic_flag, data,  &
       comment, flag_meanings, missing_value, fill_value)

   implicit none
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    real(r8)          , pointer              :: data(:,:,:)
    character(len=*)  , intent(in)           :: dim1name         ! dimension name
    character(len=*)  , intent(in)           :: dim2name         ! dimension name
    character(len=*)  , intent(in)           :: dim3name         ! dimension name
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)           , intent(in), optional :: missing_value   ! attribute for int
    real(r8)           , intent(in), optional :: fill_value      ! attribute for int


    ! Local variables
    ! Local variables
    character(len=*), parameter :: sub=trim(mod_filename)//'::'//'restartvar_real_sp_3d'
    logical           :: readvar          ! was var read?
    integer          :: ivalue
    type(var_desc_t) :: vardesc  ! local vardesc
    integer          :: status   ! return error code
    integer          :: varid
    integer          :: lxtype   ! local external type (in case logical variable)
    !----------------------------------------------------

    readvar = .false.
    if (flag == 'define') then

      lxtype = ncd_double

      call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=lxtype, &
         dim1name=trim(dim1name), dim2name=trim(dim2name), dim3name=trim(dim3name), &
         long_name=trim(long_name), units=units)

      status = nf90_inq_varid(ncid%fh,trim(varname),varid)

       if (trim(interpinic_flag) == 'interp') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_interp)
       else if (trim(interpinic_flag) == 'copy') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_copy)
       else if (trim(interpinic_flag) == 'skip') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_skip)
       end if

       status = nf90_put_att(ncid%fh, varid, 'interpinic_flag_meanings', &
            "1=nearest neighbor, 2=copy directly, 3=skip")

       if (present(comment)) then
          call check_ret(ncd_putatt(ncid, varid, 'comment', trim(comment)),sub)
       end if
       if (present(units)) then
          call check_ret(ncd_putatt(ncid, varid, 'units', trim(units)),sub)
       end if

       if (present(fill_value)) then
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', fill_value,ncd_double),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', spval,ncd_double),sub)
       end if
       if (present(missing_value)) then
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', missing_value,ncd_double),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', spval,ncd_double),sub)
       end if


    else if (flag == 'read' .or. flag == 'write') then

      call ncd_io(varname=trim(varname), data=data, &
        dim1name=trim(dim1name),  flag=flag, ncid=ncid, readvar=readvar)
    end if

    if (flag == 'read') then
       if (.not. readvar .and. is_restart()) &
         call endrun('Reading '//trim(varname)//' from restart file failed in '//trim(sub),__LINE__)
    end if

  end subroutine restartvar_real_sp_3d
!-----------------------------------------------------------------------
  subroutine restartvar_real_sp_4d(ncid, flag, varname,  dim1name,  dim2name, &
       dim3name, dim4name, long_name, units, interpinic_flag, data,  &
       comment, flag_meanings, missing_value, fill_value)

   implicit none
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    real(r8)          , pointer              :: data(:,:,:,:)
    character(len=*)  , intent(in)           :: dim1name         ! dimension name
    character(len=*)  , intent(in)           :: dim2name         ! dimension name
    character(len=*)  , intent(in)           :: dim3name         ! dimension name
    character(len=*)  , intent(in)           :: dim4name         ! dimension name
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)           , intent(in), optional :: missing_value   ! attribute for int
    real(r8)           , intent(in), optional :: fill_value      ! attribute for int


    ! Local variables
    ! Local variables
    character(len=*), parameter :: sub=trim(mod_filename)//'::'//'restartvar_real_sp_4d'
    logical           :: readvar          ! was var read?
    integer          :: ivalue
    type(var_desc_t) :: vardesc  ! local vardesc
    integer          :: status   ! return error code
    integer          :: varid
    integer          :: lxtype   ! local external type (in case logical variable)
    !----------------------------------------------------

    readvar = .false.
    if (flag == 'define') then

      lxtype = ncd_double

      call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=lxtype, &
         dim1name=trim(dim1name), dim2name=trim(dim2name), dim3name=trim(dim3name), &
         dim4name=trim(dim4name), long_name=trim(long_name), units=units)

      status = nf90_inq_varid(ncid%fh,trim(varname),varid)

       if (trim(interpinic_flag) == 'interp') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_interp)
       else if (trim(interpinic_flag) == 'copy') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_copy)
       else if (trim(interpinic_flag) == 'skip') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_skip)
       end if

       status = nf90_put_att(ncid%fh, varid, 'interpinic_flag_meanings', &
            "1=nearest neighbor, 2=copy directly, 3=skip")

       if (present(comment)) then
          call check_ret(ncd_putatt(ncid, varid, 'comment', trim(comment)),sub)
       end if
       if (present(units)) then
          call check_ret(ncd_putatt(ncid, varid, 'units', trim(units)),sub)
       end if

       if (present(fill_value)) then
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', fill_value,ncd_double),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', spval,ncd_double),sub)
       end if
       if (present(missing_value)) then
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', missing_value,ncd_double),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', spval,ncd_double),sub)
       end if


    else if (flag == 'read' .or. flag == 'write') then

      call ncd_io(varname=trim(varname), data=data, &
        dim1name=trim(dim1name), ncid=ncid, flag=flag, readvar=readvar)
    end if

    if (flag == 'read') then
       if (.not. readvar .and. is_restart()) &
         call endrun('Reading restart file failed in '//trim(sub),__LINE__)
    end if

  end subroutine restartvar_real_sp_4d
!-----------------------------------------------------------------------
  subroutine restartvar_real_sp_5d(ncid, flag, varname,  dim1name,  dim2name, &
       dim3name, dim4name, dim5name, long_name, units, interpinic_flag, data,  &
       comment, flag_meanings, missing_value, fill_value)

   implicit none
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    real(r8)          , pointer              :: data(:,:,:,:,:)
    character(len=*)  , intent(in)           :: dim1name         ! dimension name
    character(len=*)  , intent(in)           :: dim2name         ! dimension name
    character(len=*)  , intent(in)           :: dim3name         ! dimension name
    character(len=*)  , intent(in)           :: dim4name         ! dimension name
    character(len=*)  , intent(in)           :: dim5name         ! dimension name
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)           , intent(in), optional :: missing_value   ! attribute for int
    real(r8)           , intent(in), optional :: fill_value      ! attribute for int


    ! Local variables
    ! Local variables
    character(len=*), parameter :: sub=trim(mod_filename)//'::'//'restartvar_real_sp_5d'
    logical           :: readvar          ! was var read?
    integer          :: ivalue
    type(var_desc_t) :: vardesc  ! local vardesc
    integer          :: status   ! return error code
    integer          :: varid
    integer          :: lxtype   ! local external type (in case logical variable)
    !----------------------------------------------------

    readvar = .false.
    if (flag == 'define') then

      lxtype = ncd_double

      call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=lxtype, &
         dim1name=trim(dim1name), dim2name=trim(dim2name), dim3name=trim(dim3name), &
         dim4name=trim(dim4name), dim5name=trim(dim5name),long_name=trim(long_name), units=units)

      status = nf90_inq_varid(ncid%fh,trim(varname),varid)

       if (trim(interpinic_flag) == 'interp') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_interp)
       else if (trim(interpinic_flag) == 'copy') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_copy)
       else if (trim(interpinic_flag) == 'skip') then
          status = nf90_put_att(ncid%fh, varid, 'interpinic_flag', iflag_skip)
       end if

       status = nf90_put_att(ncid%fh, varid, 'interpinic_flag_meanings', &
            "1=nearest neighbor, 2=copy directly, 3=skip")

       if (present(comment)) then
          call check_ret(ncd_putatt(ncid, varid, 'comment', trim(comment)),sub)
       end if
       if (present(units)) then
          call check_ret(ncd_putatt(ncid, varid, 'units', trim(units)),sub)
       end if

       if (present(fill_value)) then
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', fill_value,ncd_double),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, '_FillValue', spval,ncd_double),sub)
       end if
       if (present(missing_value)) then
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', missing_value,ncd_double),sub)
       else
          call check_ret(ncd_putatt(ncid, varid, 'missing_value', spval,ncd_double),sub)
       end if


    else if (flag == 'read' .or. flag == 'write') then

      call ncd_io(varname=trim(varname), data=data, &
        dim1name=trim(dim1name), ncid=ncid, flag=flag, readvar=readvar)
    end if

    if (flag == 'read') then
       if (.not. readvar .and. is_restart()) &
         call endrun('Reading restart file failed in '//trim(sub),__LINE__)
    end if

  end subroutine restartvar_real_sp_5d
!------------------------------------------------------------------------------------------
  subroutine cpcol_i_1d(flag,NHW,NHE,NVN,NVS,dat_arc,datic_1d)
  implicit none
  character(len=*), intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(inout) :: dat_arc(:,:)
  integer, intent(inout):: datic_1d(:)
  integer :: NY,NX,icol

  if (flag=='read')then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        dat_arc(NY,NX)=datic_1d(icol)
      ENDDO
    ENDDO
  else if(flag=='write')then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        datic_1d(icol)=dat_arc(NY,NX)
      ENDDO
    ENDDO
  endif
  end subroutine cpcol_i_1d
!------------------------------------------------------------------------------------------
  subroutine cpcol_r_1d(flag,NHW,NHE,NVN,NVS,dat_arc,datrc_1d)
  implicit none
  character(len=*), intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), intent(inout) :: dat_arc(:,:)
  real(r8), intent(inout):: datrc_1d(:)
  integer :: NY,NX,icol

  if (flag=='read')then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        dat_arc(NY,NX)=datrc_1d(icol)
      ENDDO
    ENDDO
  else if(flag=='write')then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        datrc_1d(icol)=dat_arc(NY,NX)
      ENDDO
    ENDDO
  endif
  end subroutine cpcol_r_1d
!------------------------------------------------------------------------------------------
  subroutine cpcol_r_2d(flag,NHW,NHE,NVN,NVS,dat_arc,datrc_2d)
  implicit none
  character(len=*), intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), intent(inout) :: dat_arc(:,:,:)
  real(r8), intent(inout):: datrc_2d(:,:)
  integer :: NY,NX,icol,N

  if (flag=='read')then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        DO N=1,size(dat_arc,1)
          dat_arc(N,NY,NX)=datrc_2d(icol,N)
        enddo
      ENDDO
    ENDDO
  else if(flag=='write')then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        DO n=1,size(dat_arc,1)
          datrc_2d(icol,N)=dat_arc(N,NY,NX)
        enddo
      ENDDO
    ENDDO
  endif
  end subroutine cpcol_r_2d
!------------------------------------------------------------------------------------------
  subroutine cpcol_r_3d(flag,NHW,NHE,NVN,NVS,dat_arc,datrc_3d)
  implicit none
  character(len=*), intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), intent(inout) :: dat_arc(:,:,:,:)
  real(r8), intent(inout):: datrc_3d(:,:,:)
  integer :: NY,NX,icol,N1,N2

  if (flag=='read')then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        DO N2=1,size(dat_arc,2)
          DO N1=1,size(dat_arc,1)
            dat_arc(N1,N2,NY,NX)=datrc_3d(icol,N1,N2)
          enddo
        ENDDO
      ENDDO
    ENDDO
  else if(flag=='write')then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        DO n2=1,size(dat_arc,2)
          DO n1=1,size(dat_arc,1)
            datrc_3d(icol,N1,N2)=dat_arc(N1,N2,NY,NX)
          enddo
        enddo
      ENDDO
    ENDDO
  endif
  end subroutine cpcol_r_3d
!------------------------------------------------------------------------------------------
  subroutine cpcol_r_4d(flag,NHW,NHE,NVN,NVS,dat_arc,datrc_4d)
  implicit none
  character(len=*), intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), intent(inout) :: dat_arc(:,:,:,:,:)
  real(r8), intent(inout):: datrc_4d(:,:,:,:)
  integer :: NY,NX,icol,N1,N2,N3

  if (flag=='read')then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        DO N3=1,size(dat_arc,3)
          DO N2=1,size(dat_arc,2)
            DO N1=1,size(dat_arc,1)
              dat_arc(N1,N2,N3,NY,NX)=datrc_4d(icol,N1,N2,N3)
            enddo
          enddo
        ENDDO
      ENDDO
    ENDDO
  else if(flag=='write')then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        DO n3=1,size(dat_arc,3)
          DO n2=1,size(dat_arc,2)
            DO n1=1,size(dat_arc,1)
              datrc_4d(icol,N1,N2,N3)=dat_arc(N1,N2,N3,NY,NX)
            enddo
          enddo
        enddo
      ENDDO
    ENDDO
  endif
  end subroutine cpcol_r_4d
!------------------------------------------------------------------------------------------
  subroutine cpcol_r_5d(flag,NHW,NHE,NVN,NVS,dat_arc,datrc_5d)
  implicit none
  character(len=*), intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), intent(inout) :: dat_arc(:,:,:,:,:,:)
  real(r8), intent(inout):: datrc_5d(:,:,:,:,:)
  integer :: NY,NX,icol,N1,N2,N3,N4

  if (flag=='read')then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        DO n4=1,size(dat_arc,4)
          DO N3=1,size(dat_arc,3)
            DO N2=1,size(dat_arc,2)
              DO N1=1,size(dat_arc,1)
                dat_arc(N1,N2,N3,N4,NY,NX)=datrc_5d(icol,N1,N2,N3,N4)
              enddo
            enddo
          enddo
        ENDDO
      ENDDO
    ENDDO
  else if(flag=='write')then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        DO n4=1,size(dat_arc,4)
          DO n3=1,size(dat_arc,3)
            DO n2=1,size(dat_arc,2)
              DO n1=1,size(dat_arc,1)
                datrc_5d(icol,N1,N2,N3,N4)=dat_arc(N1,N2,N3,N4,NY,NX)
              enddo
            enddo
          enddo
        enddo
      ENDDO
    ENDDO
  endif
  end subroutine cpcol_r_5d

!------------------------------------------------------------------------------------------
  subroutine cppft_i_1d(flag,NHW,NHE,NVN,NVS,NP,dat_arp,datip_1d,NumActivePlants,IsPlantActive_pft)
  implicit none
  character(len=*),intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(in) :: NP(:,:)
  integer, intent(inout) :: dat_arp(:,:,:)
  integer, intent(inout):: datip_1d(:)
  integer, intent(in), optional :: NumActivePlants(:,:)
  integer, intent(in), optional :: IsPlantActive_pft(:,:,:)
  integer :: NZ,NY,NX,ip

  if (flag=='read')then
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  dat_arp(NZ,NY,NX)=datip_1d(ip)
                endif
              else
                ip=get_pft(NZ,NY,NX)
                dat_arp(NZ,NY,NX)=datip_1d(ip)
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          DO NZ=1,NP(NY,NX)
            ip=get_pft(NZ,NY,NX)
            dat_arp(NZ,NY,NX)=datip_1d(ip)
          ENDDO
        ENDDO
      ENDDO
    endif
  else if(flag=='write') then
    datip_1d=ispval
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  datip_1d(ip)=dat_arp(NZ,NY,NX)
                endif
              else
                ip=get_pft(NZ,NY,NX)
                datip_1d(ip)=dat_arp(NZ,NY,NX)
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
         DO NY=NVN,NVS
           DO NZ=1,NP(NY,NX)
             ip=get_pft(NZ,NY,NX)
             datip_1d(ip)=dat_arp(NZ,NY,NX)
           ENDDO
         ENDDO
      ENDDO
    endif
  endif
  end subroutine cppft_i_1d

!------------------------------------------------------------------------------------------
  subroutine cppft_i_2d(flag,NHW,NHE,NVN,NVS,NP,dat_arp,datip_2d,NumActivePlants,IsPlantActive_pft)
  implicit none
  character(len=*),intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(in) :: NP(:,:)
  integer, intent(inout) :: dat_arp(:,:,:,:)
  integer, intent(inout):: datip_2d(:,:)
  integer, intent(in), optional :: NumActivePlants(:,:)
  integer, intent(in), optional :: IsPlantActive_pft(:,:,:)
  integer :: NZ,NY,NX,ip,NN

  if (flag=='read')then
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  DO NN=1,SIZE(dat_arp,1)
                    dat_arp(NN,NZ,NY,NX)=datip_2d(ip,NN)
                  enddo
                endif
              else
                ip=get_pft(NZ,NY,NX)
                DO NN=1,SIZE(dat_arp,1)
                  dat_arp(NN,NZ,NY,NX)=datip_2d(ip,NN)
                ENDDO
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          DO NZ=1,NP(NY,NX)
            ip=get_pft(NZ,NY,NX)
            DO NN=1,SIZE(dat_arp,1)
              dat_arp(NN,NZ,NY,NX)=datip_2d(ip,NN)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    endif
  else if(flag=='write') then
    datip_2d=ispval
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS

          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  DO NN=1,SIZE(dat_arp,1)
                    datip_2d(ip,NN)=dat_arp(NN,NZ,NY,NX)
                  enddo
                endif
              else
                ip=get_pft(NZ,NY,NX)
                DO NN=1,SIZE(dat_arp,1)
                  datip_2d(ip,NN)=dat_arp(NN,NZ,NY,NX)
                enddo
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
         DO NY=NVN,NVS
           DO NZ=1,NP(NY,NX)
             ip=get_pft(NZ,NY,NX)
             do nn=1,size(dat_arp,1)
               datip_2d(ip,NN)=dat_arp(nn,nz,NY,NX)
             enddo
           ENDDO
         ENDDO
      ENDDO
    endif
  endif
  end subroutine cppft_i_2d

!------------------------------------------------------------------------------------------
  subroutine cppft_i_3d(flag,NHW,NHE,NVN,NVS,NP,dat_arp,datip_3d,NumActivePlants,IsPlantActive_pft)
  implicit none
  character(len=*),intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(in) :: NP(:,:)
  integer, intent(inout) :: dat_arp(:,:,:,:,:)
  integer, intent(inout):: datip_3d(:,:,:)
  integer, intent(in), optional :: NumActivePlants(:,:)
  integer, intent(in), optional :: IsPlantActive_pft(:,:,:)
  integer :: NZ,NY,NX,ip,N1,N2

  if (flag=='read')then
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  DO N2=1,SIZE(dat_arp,2)
                    DO N1=1,SIZE(dat_arp,1)
                      dat_arp(N1,N2,NZ,NY,NX)=datip_3d(ip,N1,N2)
                    ENDDO
                  enddo
                endif
              else
                ip=get_pft(NZ,NY,NX)
                DO N2=1,SIZE(dat_arp,2)
                  DO N1=1,SIZE(dat_arp,1)
                    dat_arp(N1,N2,NZ,NY,NX)=datip_3d(ip,N1,N2)
                  ENDDO
                ENDDO
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          DO NZ=1,NP(NY,NX)
            ip=get_pft(NZ,NY,NX)
            DO N2=1,SIZE(dat_arp,2)
              DO N1=1,SIZE(dat_arp,1)
                dat_arp(N1,N2,NZ,NY,NX)=datip_3d(ip,N1,N2)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    endif
  else if(flag=='write') then
    datip_3d=ispval
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  DO N2=1,SIZE(dat_arp,2)
                    DO N1=1,SIZE(dat_arp,1)
                      datip_3d(ip,N1,N2)=dat_arp(N1,N2,NZ,NY,NX)
                    ENDDO
                  enddo
                endif
              else
                ip=get_pft(NZ,NY,NX)
                DO N2=1,SIZE(dat_arp,2)
                  DO N1=1,SIZE(dat_arp,1)
                    datip_3d(ip,N1,N2)=dat_arp(N1,N2,NZ,NY,NX)
                  ENDDO
                enddo
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
         DO NY=NVN,NVS
           DO NZ=1,NP(NY,NX)
             ip=get_pft(NZ,NY,NX)
             DO N2=1,SIZE(dat_arp,2)
               DO N1=1,SIZE(dat_arp,1)
                 datip_3d(ip,N1,N2)=dat_arp(n1,n2,nz,NY,NX)
               ENDDO
             enddo
           ENDDO
         ENDDO
      ENDDO
    endif
  endif
  end subroutine cppft_i_3d

!------------------------------------------------------------------------------------------
  subroutine cppft_r_1d(flag,NHW,NHE,NVN,NVS,NP,dat_arp,datrp_1d,NumActivePlants,IsPlantActive_pft)
  implicit none
  character(len=*),intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(in) :: NP(:,:)
  real(r8), intent(inout) :: dat_arp(:,:,:)
  real(r8), intent(inout):: datrp_1d(:)
  integer, intent(in), optional :: NumActivePlants(:,:)
  integer, intent(in), optional :: IsPlantActive_pft(:,:,:)
  integer :: NZ,NY,NX,ip

  if (flag=='read')then
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  dat_arp(NZ,NY,NX)=datrp_1d(ip)
                endif
              else
                ip=get_pft(NZ,NY,NX)
                dat_arp(NZ,NY,NX)=datrp_1d(ip)
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          DO NZ=1,NP(NY,NX)
            ip=get_pft(NZ,NY,NX)
            dat_arp(NZ,NY,NX)=datrp_1d(ip)
          ENDDO
        ENDDO
      ENDDO
    endif
  else if(flag=='write') then
    datrp_1d=spval
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  datrp_1d(ip)=dat_arp(NZ,NY,NX)
                endif
              else
                ip=get_pft(NZ,NY,NX)
                datrp_1d(ip)=dat_arp(NZ,NY,NX)
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
         DO NY=NVN,NVS
         DO NZ=1,NP(NY,NX)
            ip=get_pft(NZ,NY,NX)
            datrp_1d(ip)=dat_arp(NZ,NY,NX)
         ENDDO
         ENDDO
      ENDDO
    endif
  endif
  end subroutine cppft_r_1d

!------------------------------------------------------------------------------------------
  subroutine cppft_r_2d(flag,NHW,NHE,NVN,NVS,NP,dat_arp,datrp_2d,NumActivePlants,IsPlantActive_pft)
  implicit none
  character(len=*),intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(in) :: NP(:,:)
  real(r8), intent(inout) :: dat_arp(:,:,:,:)    !data to read or write
  real(r8), intent(inout):: datrp_2d(:,:)
  integer, intent(in), optional :: NumActivePlants(:,:)
  integer, intent(in), optional :: IsPlantActive_pft(:,:,:)
  integer :: NZ,NY,NX,ip,NN

  if (flag=='read')then
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  DO NN=1,SIZE(dat_arp,1)
                    dat_arp(NN,NZ,NY,NX)=datrp_2d(ip,NN)
                  enddo
                endif
              else
                ip=get_pft(NZ,NY,NX)
                DO NN=1,SIZE(dat_arp,1)
                  dat_arp(NN,NZ,NY,NX)=datrp_2d(ip,NN)
                enddo
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          DO NZ=1,NP(NY,NX)
            ip=get_pft(NZ,NY,NX)
            DO NN=1,SIZE(dat_arp,1)
              dat_arp(NN,NZ,NY,NX)=datrp_2d(ip,NN)
            enddo
          ENDDO
        ENDDO
      ENDDO
    endif
  else if(flag=='write') then
    datrp_2d=spval
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  DO NN=1,SIZE(dat_arp,1)
                    datrp_2d(ip,NN)=dat_arp(NN,NZ,NY,NX)
                  enddo
                endif
              else
                ip=get_pft(NZ,NY,NX)
                DO NN=1,SIZE(dat_arp,1)
                  datrp_2d(ip,NN)=dat_arp(NN,NZ,NY,NX)
                enddo
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
         DO NY=NVN,NVS
           DO NZ=1,NP(NY,NX)
             ip=get_pft(NZ,NY,NX)
             DO NN=1,SIZE(dat_arp,1)
               datrp_2d(ip,NN)=dat_arp(NN,NZ,NY,NX)
             enddo
           ENDDO
         ENDDO
      ENDDO
    endif
  endif
  end subroutine cppft_r_2d

!------------------------------------------------------------------------------------------
  subroutine cppft_r_3d(flag,NHW,NHE,NVN,NVS,NP,dat_arp,datrp_3d,NumActivePlants,IsPlantActive_pft)
  implicit none
  character(len=*),intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(in) :: NP(:,:)
  real(r8), intent(inout) :: dat_arp(:,:,:,:,:)
  real(r8), intent(inout):: datrp_3d(:,:,:)
  integer, intent(in), optional :: NumActivePlants(:,:)
  integer, intent(in), optional :: IsPlantActive_pft(:,:,:)
  integer :: NZ,NY,NX,ip,N1,N2

  if (flag=='read')then
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                DO N2=1,SIZE(dat_arp,2)
                  DO N1=1,SIZE(dat_arp,1)
                    dat_arp(N1,N2,NZ,NY,NX)=datrp_3d(ip,N1,N2)
                  enddo
                enddo
                endif
              else
                ip=get_pft(NZ,NY,NX)
                DO N2=1,SIZE(dat_arp,2)
                  DO N1=1,SIZE(dat_arp,1)
                    dat_arp(N1,N2,NZ,NY,NX)=datrp_3d(ip,N1,N2)
                  enddo
                enddo
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          DO NZ=1,NP(NY,NX)
            ip=get_pft(NZ,NY,NX)
             DO N2=1,SIZE(dat_arp,2)
               DO N1=1,SIZE(dat_arp,1)
                 dat_arp(N1,N2,NZ,NY,NX)=datrp_3d(ip,N1,N2)
               enddo
             enddo
          ENDDO
        ENDDO
      ENDDO
    endif
  else if(flag=='write') then
    datrp_3d=spval
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  DO N2=1,SIZE(dat_arp,2)
                    DO N1=1,SIZE(dat_arp,1)
                      datrp_3d(ip,N1,N2)=dat_arp(N1,N2,NZ,NY,NX)
                    enddo
                  enddo
                endif
              else
                ip=get_pft(NZ,NY,NX)
                DO N2=1,SIZE(dat_arp,2)
                  DO N1=1,SIZE(dat_arp,1)
                    datrp_3d(ip,N1,N2)=dat_arp(N1,N2,NZ,NY,NX)
                  enddo
                enddo
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
         DO NY=NVN,NVS
           DO NZ=1,NP(NY,NX)
             ip=get_pft(NZ,NY,NX)
             DO N2=1,SIZE(dat_arp,2)
               DO N1=1,SIZE(dat_arp,1)
                 datrp_3d(ip,N1,N2)=dat_arp(N1,N2,NZ,NY,NX)
               enddo
             enddo
           ENDDO
         ENDDO
      ENDDO
    endif
  endif
  end subroutine cppft_r_3d
!------------------------------------------------------------------------------------------
  subroutine cppft_r_4d(flag,NHW,NHE,NVN,NVS,NP,dat_arp,datrp_4d,NumActivePlants,IsPlantActive_pft)
  implicit none
  character(len=*),intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(in) :: NP(:,:)
  real(r8), intent(inout) :: dat_arp(:,:,:,:,:,:)
  real(r8), intent(inout):: datrp_4d(:,:,:,:)
  integer, intent(in), optional :: NumActivePlants(:,:)
  integer, intent(in), optional :: IsPlantActive_pft(:,:,:)
  integer :: NZ,NY,NX,ip,N1,N2,N3

  if (flag=='read')then
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  DO N3=1,SIZE(dat_arp,3)
                    DO N2=1,SIZE(dat_arp,2)
                      DO N1=1,SIZE(dat_arp,1)
                        dat_arp(N1,N2,N3,NZ,NY,NX)=datrp_4d(ip,N1,N2,N3)
                      enddo
                  enddo
                enddo
                endif
              else
                ip=get_pft(NZ,NY,NX)
                DO N3=1,SIZE(dat_arp,3)
                  DO N2=1,SIZE(dat_arp,2)
                    DO N1=1,SIZE(dat_arp,1)
                      dat_arp(N1,N2,N3,NZ,NY,NX)=datrp_4d(ip,N1,N2,N3)
                    enddo
                  enddo
                enddo
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          DO NZ=1,NP(NY,NX)
            ip=get_pft(NZ,NY,NX)
             DO N3=1,SIZE(dat_arp,3)
               DO N2=1,SIZE(dat_arp,2)
                 DO N1=1,SIZE(dat_arp,1)
                   dat_arp(N1,N2,N3,NZ,NY,NX)=datrp_4d(ip,N1,N2,N3)
                 enddo
               enddo
             enddo
          ENDDO
        ENDDO
      ENDDO
    endif
  else if(flag=='write') then
    datrp_4d=spval
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  DO N3=1,SIZE(dat_arp,3)
                    DO N2=1,SIZE(dat_arp,2)
                      DO N1=1,SIZE(dat_arp,1)
                        datrp_4d(ip,N1,N2,N3)=dat_arp(N1,N2,N3,NZ,NY,NX)
                      enddo
                    enddo
                  enddo
                endif
              else
                ip=get_pft(NZ,NY,NX)
                DO N3=1,SIZE(dat_arp,3)
                  DO N2=1,SIZE(dat_arp,2)
                    DO N1=1,SIZE(dat_arp,1)
                      datrp_4d(ip,N1,N2,N3)=dat_arp(N1,N2,N3,NZ,NY,NX)
                    enddo
                  enddo
                enddo
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
         DO NY=NVN,NVS
           DO NZ=1,NP(NY,NX)
             ip=get_pft(NZ,NY,NX)
             DO N3=1,SIZE(dat_arp,3)
               DO N2=1,SIZE(dat_arp,2)
                  DO N1=1,SIZE(dat_arp,1)
                    datrp_4d(ip,N1,N2,N3)=dat_arp(N1,N2,N3,NZ,NY,NX)
                  enddo
               enddo
             enddo
           ENDDO
         ENDDO
      ENDDO
    endif
  endif
  end subroutine cppft_r_4d
!------------------------------------------------------------------------------------------
  subroutine cppft_r_5d(flag,NHW,NHE,NVN,NVS,NP,dat_arp,datrp_5d,NumActivePlants,IsPlantActive_pft)
  implicit none
  character(len=*),intent(in) :: flag
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(in) :: NP(:,:)
  real(r8), intent(inout) :: dat_arp(:,:,:,:,:,:,:)
  real(r8), intent(inout):: datrp_5d(:,:,:,:,:)
  integer, intent(in), optional :: NumActivePlants(:,:)
  integer, intent(in), optional :: IsPlantActive_pft(:,:,:)
  integer :: NZ,NY,NX,ip,N1,N2,N3,N4

  if (flag=='read')then
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  DO N4=1,SIZE(dat_arp,4)
                    DO N3=1,SIZE(dat_arp,3)
                      DO N2=1,SIZE(dat_arp,2)
                        DO N1=1,SIZE(dat_arp,1)
                          dat_arp(N1,N2,N3,N4,NZ,NY,NX)=datrp_5d(ip,N1,N2,N3,N4)
                        enddo
                      enddo
                    enddo
                  enddo
                endif
              else
                ip=get_pft(NZ,NY,NX)
                DO N4=1,SIZE(dat_arp,4)
                  DO N3=1,SIZE(dat_arp,3)
                      DO N2=1,SIZE(dat_arp,2)
                        DO N1=1,SIZE(dat_arp,1)
                          dat_arp(N1,N2,N3,N4,NZ,NY,NX)=datrp_5d(ip,N1,N2,N3,N4)
                        enddo
                      enddo
                  enddo
                enddo
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          DO NZ=1,NP(NY,NX)
            ip=get_pft(NZ,NY,NX)
             DO N4=1,SIZE(dat_arp,4)
               DO N3=1,SIZE(dat_arp,3)
                  DO N2=1,SIZE(dat_arp,2)
                    DO N1=1,SIZE(dat_arp,1)
                      dat_arp(N1,N2,N3,N4,NZ,NY,NX)=datrp_5d(ip,N1,N2,N3,N4)
                    enddo
                  enddo
               enddo
             enddo
          ENDDO
        ENDDO
      ENDDO
    endif
  else if(flag=='write') then
    datrp_5d=spval
    if(present(NumActivePlants))then
      DO NX=NHW,NHE
        DO NY=NVN,NVS
          IF(NumActivePlants(NY,NX)>0)THEN
            DO NZ=1,NP(NY,NX)
              if(present(IsPlantActive_pft))then
                IF(IsPlantActive_pft(NZ,NY,NX)==iTrue)THEN
                  ip=get_pft(NZ,NY,NX)
                  DO N4=1,SIZE(dat_arp,4)
                    DO N3=1,SIZE(dat_arp,3)
                      DO N2=1,SIZE(dat_arp,2)
                        DO N1=1,SIZE(dat_arp,1)
                          datrp_5d(ip,N1,N2,N3,N4)=dat_arp(N1,N2,N3,N4,NZ,NY,NX)
                        enddo
                      enddo
                    enddo
                  enddo
                endif
              else
                ip=get_pft(NZ,NY,NX)
                DO N4=1,SIZE(dat_arp,4)
                  DO N3=1,SIZE(dat_arp,3)
                      DO N2=1,SIZE(dat_arp,2)
                        DO N1=1,SIZE(dat_arp,1)
                          datrp_5d(ip,N1,N2,N3,N4)=dat_arp(N1,N2,N3,N4,NZ,NY,NX)
                        enddo
                      enddo
                  enddo
                enddo
              endif
            ENDDO
          endif
        enddo
      enddo
    else
      DO NX=NHW,NHE
         DO NY=NVN,NVS
           DO NZ=1,NP(NY,NX)
             ip=get_pft(NZ,NY,NX)
             DO N4=1,SIZE(dat_arp,4)
               DO N3=1,SIZE(dat_arp,3)
                  DO N2=1,SIZE(dat_arp,2)
                    DO N1=1,SIZE(dat_arp,1)
                      datrp_5d(ip,N1,N2,N3,N4)=dat_arp(N1,N2,N3,N4,NZ,NY,NX)
                    enddo
                  enddo
               enddo
             enddo
           ENDDO
         ENDDO
      ENDDO
    endif
  endif
  end subroutine cppft_r_5d

end module restUtilMod
