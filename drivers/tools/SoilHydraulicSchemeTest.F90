program SoilHydraulicSchemeTest
  ! Description:
  ! Exercise scheme selection and checkpoint metadata independently of full state.
  ! Negative cases must stop with their designated diagnostic before success is printed.
  use abortutils, only: iulog
  use EcoSIMCtrlMod, only: lverb,soil_hydraulic_scheme,use_ecosys_drainage
  use SoilHydraulicSchemeMod, only: InitializeSoilHydraulicScheme,FinalizeSoilHydraulicScheme, &
    GetSoilHydraulicSchemeName,GetSoilHydraulicSchemeRevision
  use RestartMod, only: restartnc_soil_hydraulic_scheme
  use ncdio_pio, only: file_desc_t,ncd_pio_createfile,ncd_pio_openfile,ncd_pio_closefile,ncd_enddef,ncd_putatt,check_ret,check_att
  use netcdf, only: nf90_global,nf90_nowrite
  implicit none
  character(len=32) :: mode
  type(file_desc_t) :: checkpoint
  character(len=128) :: filename
  logical :: found

  lverb=.false.
  use_ecosys_drainage=.false.
  call get_command_argument(1,mode)
  if(mode=='unknown_option')then
    soil_hydraulic_scheme='unknown'
    call InitializeSoilHydraulicScheme()
    error stop 'unknown scheme was accepted'
  endif

  ! Explicit selection is case insensitive; reinitialization resets the lifetime.
  soil_hydraulic_scheme='ECOSIM'
  call InitializeSoilHydraulicScheme()
  if(GetSoilHydraulicSchemeName()/='ecosim')error stop 'scheme was not canonicalized'
  if(GetSoilHydraulicSchemeRevision()/=1)error stop 'unexpected scheme revision'
  call FinalizeSoilHydraulicScheme()
  if(GetSoilHydraulicSchemeName()/='ecosim')error stop 'lazy initialization failed'

  filename='soil-hydraulic-scheme-'//trim(mode)//'.nc'
  call ncd_pio_createfile(checkpoint,trim(filename))
  if(mode=='incomplete')then
    call check_ret(ncd_putatt(checkpoint,nf90_global,'soil_hydraulic_scheme','ecosim'),trim(mode))
  elseif(mode/='legacy')then
    call restartnc_soil_hydraulic_scheme(checkpoint,'define')
    select case(trim(mode))
    case('scheme_mismatch')
      call check_ret(ncd_putatt(checkpoint,nf90_global,'soil_hydraulic_scheme','other'),trim(mode))
    case('revision_mismatch')
      call check_ret(ncd_putatt(checkpoint,nf90_global,'soil_hydraulic_scheme_revision',2),trim(mode))
    case('policy_mismatch')
      call check_ret(ncd_putatt(checkpoint,nf90_global,'soil_hydraulic_ecosys_drainage',1),trim(mode))
    case('matching')
      ! No mutation: the recorded scheme and policy must be accepted.
    case default
      write(iulog,*) 'Unknown test mode: ',trim(mode)
      error stop 1
    end select
  endif
  call ncd_enddef(checkpoint)
  call ncd_pio_closefile(checkpoint)
  call ncd_pio_openfile(checkpoint,trim(filename),nf90_nowrite)
  if(mode=='matching')then
    call check_att(checkpoint,nf90_global,'soil_hydraulic_scheme',found)
    if(.not.found)error stop 'scheme metadata was not written'
  endif
  call restartnc_soil_hydraulic_scheme(checkpoint,'read')
  call ncd_pio_closefile(checkpoint)
  if(mode/='legacy'.and.mode/='matching')error stop 'incompatible metadata was accepted'
  call FinalizeSoilHydraulicScheme()
  write(*,*) 'SoilHydraulicSchemeTest passed'
end program SoilHydraulicSchemeTest
