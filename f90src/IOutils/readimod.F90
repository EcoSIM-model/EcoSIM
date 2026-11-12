module readiMod
! Description:
!!
! code to read site, topographic data
!
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use abortutils,       only: endrun,iulog
  use fileUtil,         only: open_safe, check_read
  use minimathmod,      only: isclose,   AZMAX1, safe_adb
  use MiniFuncMod,      only: GetDayLength
  use EcoSIMConfig,     only: column_mode
  use EcoSiMParDataMod, only: micpar
  use SoilHydroParaMod, only: ComputeSoilHydroPars
  use SoilPhysParaMod,  only: SetDeepSoil
  use DebugToolMod,     only: PrintInfo
  use ncdio_pio
  use EcoSIMCtrlMod
  use SOMDataType
  use CanopyRadDataType
  use EcosimConst
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use ChemTranspDataType
  use LandSurfDataType
  use ClimForcDataType
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SurfLitterDataType
  use SnowDataType
  use EcoSIMHistMod
  use SoilPropertyDataType
  use SoilBGCDataType
  use AqueChemDatatype
  use GridDataType
  implicit none
  private

  character(len=*), parameter :: mod_filename = &
  __FILE__
  real(r8) :: datav(40)
  CHARACTER(len=16) :: OUTW,OUTI,OUTT,OUTN,OUTF
  CHARACTER(len=4) :: CHARY
  CHARACTER(len=1) :: TTYPE,CTYPE,IVAR(20),VAR(50),TYP(50)
  character(len=3), parameter :: model_status(0:1)=(/'off','on '/)

  integer :: ll
  real(r8) :: DAT(50),DATK(50)
  real(r8) :: ALATG,ATCAG,AZI,ASPX,WTBLDepz_nat,WTBLDepz_tile
  real(r8) :: DTBLGG,DEC,initSnowDepth,RCHQNG,RCHQEG
  real(r8) :: RCHQSG,RCHQWG,RCHGNUG,RCHGEUG,RCHGSUG,RCHGWUG
  real(r8) :: RCHGNTG,RCHGETG,RCHGSTG,RCHGWTG,RCHGDG
  real(r8) :: SL0,SLX,SL1,SL2

  integer :: L,NH1,NH2,NV1,NV2,NL1
  integer :: NL2
  type(file_desc_t) :: grid_nfid

  public :: readi,erosion_model_status
  public :: GridConectionMode
  contains

  SUBROUTINE readi(NHW,NHE,NVN,NVS)
!!
! Description:
! THIS SUBROUTINE READS ALL SOIL AND TOPOGRAPHIC INPUT FILES
!
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: jj,NX,NY
  integer :: ierr
  character(len=200) :: tline
!
!
! OPEN OUTPUT LOGFILES,AND SITE,TOPOGRAPHY FILES FROM
! FILE NAMES IN DATA ARRAYS LOADED IN 'MAIN'
!

  call ncd_pio_openfile(grid_nfid, grid_file_in, ncd_nowrite)

  call readsiteNC(NHW,NHE,NVN,NVS)

  call readTopoNC

  if(.not.column_mode)then
    DO  NX=NHW,NHE
      NL_col(NVS+1,NX)=NL_col(NVS,NX)
    ENDDO
    DO  NY=NVN,NVS
      NL_col(NY,NHE+1)=NL_col(NY,NHE)
    ENDDO
    NL_col(NVS+1,NHE+1)=NL_col(NVS,NHE)
  endif
  IOLD=0
  call ncd_pio_closefile(grid_nfid)


  END SUBROUTINE readi

!------------------------------------------------------------------------------------------

  function erosion_model_status(flag)result(status)

  implicit none
  integer, intent(in) :: flag
  character(len=40) :: status
  !0 means allowing  freeze- thaw to change elevation
  !1 means allowing freeze-thaw plus erosion to change elevation
  !2 means allowing freeze-thaw plus SOC accumulation to change elevation
  !3 means allowing freeze-thaw plus SOC accumulation, plus erosion to change elevation
  !-1 means no change in elevation.
  select case (flag)
  case (-1)
    status='no elv change'
  case (0)
    status='freeze-thaw elv change'
  case (1)
    status='freeze-thaw+erosion elv change'
  case (2)
    status='freeze-thaw+SOC accum elv change'
  case (3)
    status='freeze-thaw+SOC+erosion elv change'
  case default
    status=''
    call endrun('wrong erosion model option in '//trim(mod_filename)//' at line',__LINE__)
  end select
  end function erosion_model_status

!------------------------------------------------------------------------------------------

  function GridConectionMode(NCNG)result(status)
  implicit none
  integer, intent(in) :: NCNG
  character(len=40) :: status

  select case(NCNG)
  case (1)
    status='3D grid'   !3D
  case (2)
    status='2D Along the north-south direction' !2d, no x direction
  case (3)
    status='1D vertical column'           !1d, vertical only
  case default
    status=''
    call endrun('wrong option for NCNG in '//trim(mod_filename)//' at line',__LINE__)
  end select

  end function GridConectionMode

!------------------------------------------------------------------------------------------
  function WaterTableStatus(iWaterTabelMode)result(status)

  implicit none
  integer, intent(in) :: iWaterTabelMode

  character(len=64) :: status

  select case(iWaterTabelMode)
  case (0)
    status='No external water table'
  case (1)
    status='Natural stationary external water table'
  case (2)
    status='Natural mobile external water table'
  case (3)
    status='Artificial stationary extneral water table'
  case (4)
    status='Artificial mobile external water table'
  case default
    call endrun('wrong option for iWaterTabelMode in '//trim(mod_filename)//' at line',__LINE__)
  end select
  end function WaterTableStatus

!------------------------------------------------------------------------------------------

  subroutine read_aerosol_vars(loc,varname,var_val)
  implicit none
  integer, intent(in) :: loc
  character(len=*), intent(in) :: varname
  real(r8), intent(out) :: var_val
  type(Var_desc_t) :: vardesc
  logical :: readvar

  call check_var(grid_nfid, varname, vardesc, readvar)
  if(readvar)then
    !Read the optional bioaerosol carbon fraction.
    call ncd_getvar(grid_nfid,varname,loc,var_val)
  else
    var_val=0._r8
  endif

  end subroutine read_aerosol_vars

!------------------------------------------------------------------------------------------
  subroutine ReadVanGenuchtenProfile(ntp,NY,NX)
  !
  ! Description:
  ! Read layer parameters on the FC/WP volume basis. Alpha must declare pressure
  ! or head units; convert head-based coefficients using the model gravity factor.
  use netcdf, only: nf90_global,nf90_inq_varid,nf90_inquire_attribute,nf90_get_att,nf90_noerr,nf90_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, intent(in) :: ntp     ! topographic-unit record in the input dataset [-]
  integer, intent(in) :: NY,NX   ! source grid row and column [-]
  type(Var_desc_t) :: vardesc
  logical :: found
  integer :: variable,status,attribute_type,attribute_length,L
  real(r8) :: alpha_factor
  character(len=64) :: units
  character(len=*), parameter :: subname='ReadVanGenuchtenProfile'

  call PrintInfo('beg '//subname)
  associate(handle => grid_nfid%fh)
  if(allocated(VGThetaR_vr))then
    call check_var(grid_nfid,'VG_THETA_R',vardesc,found)
    if(.not.found)call endrun('van_genuchten requires grid variable VG_THETA_R')
    call ncd_getvar(grid_nfid,'VG_THETA_R',ntp,VGThetaR_vr(1:JZ,NY,NX))
    call check_var(grid_nfid,'VG_ALPHA',vardesc,found)
    if(.not.found)call endrun('van_genuchten requires grid variable VG_ALPHA')
    call ncd_getvar(grid_nfid,'VG_ALPHA',ntp,VGAlpha_vr(1:JZ,NY,NX))
    call check_var(grid_nfid,'VG_N',vardesc,found)
    if(.not.found)call endrun('van_genuchten requires grid variable VG_N')
    call ncd_getvar(grid_nfid,'VG_N',ntp,VGN_vr(1:JZ,NY,NX))
    call check_var(grid_nfid,'VG_L',vardesc,found)
    VGL_vr(1:JZ,NY,NX)=.5_r8
    if(found)call ncd_getvar(grid_nfid,'VG_L',ntp,VGL_vr(1:JZ,NY,NX))


    call check_ret(nf90_inq_varid(handle,'VG_ALPHA',variable),subname)
    status=nf90_inquire_attribute(handle,variable,'units',xtype=attribute_type,len=attribute_length)
    if(status/=nf90_noerr)call endrun('VG_ALPHA requires explicit units: MPa-1, m-1, or cm-1')
    if(attribute_type/=nf90_char.or.attribute_length>len(units))call endrun('Invalid VG_ALPHA units attribute')
    units=''
    call check_ret(nf90_get_att(handle,variable,'units',units),subname)

    select case(trim(adjustl(units)))
    case('MPa-1','MPa^-1','1/MPa')
      alpha_factor=1._r8
    case('m-1','m^-1','1/m')
      alpha_factor=1._r8/mGravAccelerat
    case('cm-1','cm^-1','1/cm')
      alpha_factor=100._r8/mGravAccelerat
    case default
      write(iulog,*) 'Unsupported VG_ALPHA units: ',trim(units)
      call endrun('VG_ALPHA requires explicit units: MPa-1, m-1, or cm-1')
    end select
    do L=1,JZ
      if(ieee_is_finite(VGAlpha_vr(L,NY,NX)))then
        if(ABS(VGAlpha_vr(L,NY,NX))>huge(1._r8)/alpha_factor)call endrun('VG_ALPHA unit conversion overflows')
        VGAlpha_vr(L,NY,NX)=VGAlpha_vr(L,NY,NX)*alpha_factor
      endif
    enddo
  endif
  end associate
  call PrintInfo('end '//subname)
  end subroutine ReadVanGenuchtenProfile
!------------------------------------------------------------------------------------------

  subroutine readsiteNC(NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8) :: DHI(JX),DVI(JY)
  character(len=*), parameter :: subname='readsiteNC'
  character(len=200) :: tline
  integer :: XI
  integer :: NY,NX
  integer :: IETYPG
  integer :: ierr,jj,loc
  integer :: iWaterTabelMode
  real(r8) :: fMoss,fLich,fLiveMB,fDeadMB,fDeadNMB,fDOM
!
! READ SITE DATA
!
  call PrintInfo('beg '//subname)
! ALATG,ALTIG,ATCAG=latitude,altitude,MAT(oC)
! iWaterTabelMode=water table flag
! :0=none
! :1,2=natural stationary,mobile
! :3,4=artificial stationary,mobile
! IETYPG,iErosionMode=Koppen climate zone,erosion options
! NCNG=1:lateral connections between grid cells,3:no connections
! WTBLDepz_nat,WTBLDepz_tile=depth of natural,artificial (tile) water table (iWaterTabelMode)
! DTBLGG=slope of natural water table relative to landscape surface
! RCHQNG,RCHQEG,RCHQSG,RCHQWG=boundary condns for N,E,S,W surface runoff
! RCHGNUG,RCHGEUG,RCHGSUG,RCHGWUG=bound condns for N,E,S,W subsurf flow
! RCHGNTG,RCHGETG,RCHGSTG,RCHGWTG=N,E,S,W distance to water table (m)
! RCHGDG=lower boundary conditions for water flow
! DHI=width of each W-E landscape column
! DVI=width of each N-S landscape row

! assuming only one grid, at the moment
  loc=1
  call ncd_getvar(grid_nfid,'ALATG',loc,ALATG)
  call ncd_getvar(grid_nfid,'ALTIG',loc,ALTIG)
  call ncd_getvar(grid_nfid,'ATCAG',loc,ATCAG)
  call ncd_getvar(grid_nfid,'IDTBLG',loc,iWaterTabelMode)

  call ncd_getvar(grid_nfid,'IETYPG',loc,IETYPG)
  call ncd_getvar(grid_nfid,'DTBLIG',loc,WTBLDepz_nat)
  call ncd_getvar(grid_nfid,'DTBLDIG',loc,WTBLDepz_tile)
  call ncd_getvar(grid_nfid,'DTBLGG',loc,DTBLGG)

  call ncd_getvar(grid_nfid,'RCHQNG',loc,RCHQNG)
  call ncd_getvar(grid_nfid,'RCHQEG',loc,RCHQEG)
  call ncd_getvar(grid_nfid,'RCHQSG',loc,RCHQSG)
  call ncd_getvar(grid_nfid,'RCHQWG',loc,RCHQWG)
  call ncd_getvar(grid_nfid,'RCHGNUG',loc,RCHGNUG)
  call ncd_getvar(grid_nfid,'RCHGEUG',loc,RCHGEUG)
  call ncd_getvar(grid_nfid,'RCHGSUG',loc,RCHGSUG)
  call ncd_getvar(grid_nfid,'RCHGWUG',loc,RCHGWUG)
  call ncd_getvar(grid_nfid,'RCHGNTG',loc,RCHGNTG)
  call ncd_getvar(grid_nfid,'RCHGETG',loc,RCHGETG)
  call ncd_getvar(grid_nfid,'RCHGSTG',loc,RCHGSTG)
  call ncd_getvar(grid_nfid,'RCHGWTG',loc,RCHGWTG)
  call ncd_getvar(grid_nfid,'RCHGDG',loc,RCHGDG)

  call ncd_getvar(grid_nfid,'DHI',loc,DHI(1:NHE))
  call ncd_getvar(grid_nfid,'DVI',loc,DVI(1:NVS))

  call read_aerosol_vars(loc,'fAeroDOM',fDOM)
  call read_aerosol_vars(loc,'fAeroDeadMB',fDeadMB)
  call read_aerosol_vars(loc,'fAeroLiveMB',fLiveMB)
  call read_aerosol_vars(loc,'fAeroMoss',fMoss)
  call read_aerosol_vars(loc,'fAeroLich',fLich)
  call read_aerosol_vars(loc,'fAeroDeadNMB',fDeadNMB)

  if(lverb)then
    write(*,*)'read site data file: ',DATA1(1)
    write(*,'(40A)')('-',ll=1,40)
    write(*,*)'Latitude (o): ALATG',ALATG
    write(*,*)'Altitude (m): ALTIG',ALTIG
    write(*,*)'Mean annual temperaure (oC): ATCAG',ATCAG
    write(*,'(40A)')('-',ll=1,40)
    write(*,*)'Koppen climate zone: IETYPG',IETYPG
    write(*,'(40A)')('-',ll=1,40)
    write(*,*)'width of each W-E landscape column: DHI'
    write(*,*)(DHI(NX),NX=1,NHE)
    write(*,*)'width of each N-S landscape row: DVI'
    write(*,*)(DVI(NY),NY=1,NVS)
    write(*,'(100A)')('=',ll=1,100)
  endif
  write(iulog,*)'water table config:',WaterTableStatus(iWaterTabelMode)
  write(iulog,*)'depth of natural water table: DTBLIG',WTBLDepz_nat
  write(iulog,*)'depth of artificial water table: DTBLDIG',WTBLDepz_tile
  write(iulog,*)'slope of natural water table relative to landscape surface: DTBLGG',DTBLGG

  write(iulog,*)'boundary condns for N surface runoff: RCHQNG',RCHQNG
  write(iulog,*)'boundary condns for E surface runoff: RCHQEG',RCHQEG
  write(iulog,*)'boundary condns for S surface runoff: RCHQSG',RCHQSG
  write(iulog,*)'boundary condns for W surface runoff: RCHQWG',RCHQWG

  write(iulog,*)'bound condns for N subsurf flow: RCHGNUG',RCHGNUG
  write(iulog,*)'bound condns for E subsurf flow: RCHGEUG',RCHGEUG
  write(iulog,*)'bound condns for S subsurf flow: RCHGSUG',RCHGSUG
  write(iulog,*)'bound condns for W subsurf flow: RCHGWUG',RCHGWUG

  write(iulog,*)'N distance to water table (m): RCHGNTG',RCHGNTG
  write(iulog,*)'E distance to water table (m): RCHGETG',RCHGETG
  write(iulog,*)'S distance to water table (m): RCHGSTG',RCHGSTG
  write(iulog,*)'W distance to water table (m): RCHGWTG',RCHGWTG
  write(iulog,*)'Lower boundary scaling rate for water flow:RCHGDG', RCHGDG
  
  D9895: DO NX=NHW,NHE
    D9890: DO NY=NVN,NVS
      !Partition bioaerosol carbon among organic matter and live biomass.
      f_aerosol_DOM_col(NY,NX)     = fDOM     !dissolved organic matter
      f_aerosol_DeadMB_col(NY,NX)  = fDeadMB  !dead microbial biomass
      f_aerosol_DeadNMB_col(NY,NX) = fDeadNMB !dead nonmicrobial biomass
      f_aerosol_LiveMB_col(NY,NX)  = fLiveMB  !live microbial biomass
      f_aerosol_MossB_col(NY,NX)   = fMoss    !live moss
      f_aerosol_LichB_col(NY,NX)   = fLich    !live lichen

      ALAT_col(NY,NX)                  = ALATG
      AltitudeLandscape_col(NY,NX)     = ALTIG
      PBOT_col(NY,NX)                  = PBOT_col(NY,NX)*exp(-AltitudeLandscape_col(NY,NX)/hpresc)
      ATCAI_col(NY,NX)                 = ATCAG
      TairKClimMean_col(NY,NX)         = ATCAG+273.15_r8
      IDWaterTable_col(NY,NX)          = iWaterTabelMode
      OXYE_col(NY,NX)                  = ao2_ppm
      Z2GE_col(NY,NX)                  = an2_ppm
      CO2EI_col(NY,NX)                 = aco2_ppm
      CH4E_col(NY,NX)                  = ach4_ppm
      Z2OE_col(NY,NX)                  = an2o_ppm
      ARGE_col(NY,NX)                  = arg_ppm
      ZNH3E_col(NY,NX)                 = anh3_ppm
      KoppenClimZone_col(NY,NX)        = IETYPG
      FlowDirIndicator_col(NY,NX)      = grid_mode
      NaturalWaterTableDepz_col(NY,NX) = WTBLDepz_nat
      TileWaterTableDepz_col(NY,NX)    = WTBLDepz_tile
      WaterTable2GroundSlope_col(NY,NX)       = DTBLGG

      RechargNorthSurf_col(NY,NX) = RCHQNG
      RechargEastSurf_col(NY,NX)  = RCHQEG
      RechargSouthSurf_col(NY,NX) = RCHQSG
      RechargWestSurf_col(NY,NX)  = RCHQWG 

      RechrgDistNorthSubSurf_col(NY,NX) = RCHGNTG
      RechrgDistEastSubSurf_col(NY,NX)  = RCHGETG
      RechrgDistSouthSubSurf_col(NY,NX) = RCHGSTG
      RechrgDistWestSubSurf_col(NY,NX)  = RCHGWTG

      RechargRateNorthWTBL_col(NY,NX) = RCHGNUG
      RechargRateEastWTBL_col(NY,NX)  = RCHGEUG
      RechargRateSouthWTBL_col(NY,NX) = RCHGSUG
      RechargRateWestWTBL_col(NY,NX)  = RCHGWUG      
      RechargBottom_col(NY,NX) = RCHGDG

      DH_col(NY,NX)   = DHI(NX)
      DV_col(NY,NX)   = DVI(NY)
      CO2E_col(NY,NX) = CO2EI_col(NY,NX)
      H2GE_col(NY,NX) = ah2_ppm
!
!     CALCULATE MAXIMUM DAYLENTH FOR PLANT PHENOLOGY
!
!     DayLenthMax_col=maximum daylength (h)
!
      IF(ALAT_col(NY,NX).GT.0.0_r8)THEN
        XI=173
      ELSE
        XI=356
      ENDIF
      DayLenthMax_col(NY,NX)=GetDayLength(ALAT_col(NY,NX),XI)

    ENDDO D9890
  ENDDO D9895
  call PrintInfo('end '//subname)
  end subroutine readsiteNC

!------------------------------------------------------------------------------------------

  subroutine readTopoNC()
!
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer :: jj,NY,NX,ll
  character(len=200) :: tline
  character(len=*), parameter :: subname='readTopoNC'
  integer :: NM(JY,JX),ntp,ntopus
  real(r8) :: dat1(1:JZ)
  real(r8) :: corrector,OrgVolFrac
  real(r8) :: fAeroScalar
  logical :: readvar_scalar
  type(Var_desc_t) :: vardesc
! begin_execution
  associate(                            &
  k_woody_comp => micpar%k_woody_comp , &
  k_fine_comp  => micpar%k_fine_comp  , &
  k_manure     => micpar%k_manure       &
  )
  call PrintInfo('beg '//subname)
!
! READ TOPOGRAPHY DATA AND SOIL FILE NAME FOR EACH GRID CELL
!
! for each unit within the landscape:
! NH1,NV1,NH2,NV2=NW,SE column,row
! ASPX=N,E,S,W aspect (o)
! SL0= slope (o)
! DPTHSX=initial snowpack depth
! SLX: not used

  ntopus=get_dim_len(grid_nfid, 'ntopou')
  if(first_topou)ntopus=1
  call check_var(grid_nfid, 'fAeroScalar', vardesc, readvar_scalar)
  fAeroScalar=0._r8
  DO ntp=1,ntopus
    call ncd_getvar(grid_nfid, 'NH1', ntp, NH1)
    call ncd_getvar(grid_nfid, 'NV1', ntp, NV1)
    call ncd_getvar(grid_nfid, 'NH2', ntp, NH2)
    call ncd_getvar(grid_nfid, 'NV2', ntp, NV2)
    call ncd_getvar(grid_nfid, 'ASPX', ntp,ASPX) !clockwise degrees from north
    call ncd_getvar(grid_nfid, 'SL0', ntp,SL0)
    call ncd_getvar(grid_nfid, 'DPTHSX', ntp,initSnowDepth)
    
    if(readvar_scalar)then
      !Read the optional bioaerosol carbon fraction.
      call ncd_getvar(grid_nfid,'fAeroScalar',ntp,fAeroScalar)      
    endif

!
! OPEN AND READ SOIL FILE
!

    DO  NX=NH1,NH2
      DO  NY=NV1,NV2
!
!     SURFACE SLOPES AND ASPECTS
!
        fAeroScalar_col(NY,NX)=fAeroScalar
        ASP_col(NY,NX)       = ASPX
        SL_col(NY,NX)        = SL0
        SnowDepth_col(NY,NX) = initSnowDepth
        !
        !     CONVERT ASPECT from geographic format TO GEOMETRIC FORMAT
        !
        !     what is geometric format mean? geographic format 0 is north, 90 east, 180 south, (clockwise)
        !     geometric or mathematical format: 0/360 east, 90 north, 180 west, 270 south (counter-clockwise)
        ASP_col(NY,NX)=450.0_r8-ASP_col(NY,NX)
        IF(ASP_col(NY,NX).GE.360.0_r8)ASP_col(NY,NX)=ASP_col(NY,NX)-360.0_r8
      ENDDO
    ENDDO

    call ncd_getvar(grid_nfid, 'PSIFC', ntp,PSIAtFldCapacity_col(NV1,NH1))
    call ncd_getvar(grid_nfid, 'PSIWP', ntp,PSIAtWiltPoint_col(NV1,NH1))
    call ncd_getvar(grid_nfid, 'ALBS',  ntp,SoilAlbedo_col(NV1,NH1))
    call ncd_getvar(grid_nfid, 'PH0',   ntp,PH_vr(0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSCf',  ntp,RSC_vr(k_fine_comp,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSNf',  ntp,RSN_vr(k_fine_comp,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSPf',  ntp,RSP_vr(k_fine_comp,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSCw',  ntp,RSC_vr(k_woody_comp,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSNw',  ntp,RSN_vr(k_woody_comp,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSPw',  ntp,RSP_vr(k_woody_comp,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSCm',  ntp,RSC_vr(k_manure,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSNm',  ntp,RSN_vr(k_manure,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'RSPm',  ntp,RSP_vr(k_manure,0,NV1,NH1))
    call ncd_getvar(grid_nfid, 'IXTYP1',ntp,iLitrType_col(1,NV1,NH1))
    call ncd_getvar(grid_nfid, 'IXTYP2',ntp,iLitrType_col(2,NV1,NH1))
    call ncd_getvar(grid_nfid, 'NUI'   ,ntp,NUI_col(NV1,NH1))
    call ncd_getvar(grid_nfid, 'NJ'    ,ntp,MaxNumRootLays_col(NV1,NH1))
    call ncd_getvar(grid_nfid, 'NL1'   ,ntp,NL1)
    call ncd_getvar(grid_nfid, 'NL2'   ,ntp,NL2)
    call ncd_getvar(grid_nfid, 'ISOILR',ntp,ISOILR_col(NV1,NH1))

    NU_col(NV1,NH1) = NUI_col(NV1,NH1)    
    NK_col(NV1,NH1) = MaxNumRootLays_col(NV1,NH1)+1
    NM(NV1,NH1)     = MaxNumRootLays_col(NV1,NH1)+NL1
    
    

!  the extra soil layer below root zone cannot be greater than what is allowed
    NL2=min0(JZ-NM(NV1,NH1),NL2)
    NLI_col(NV1,NH1) = NM(NV1,NH1)+NL2
    NL_col(NV1,NH1)  = NLI_col(NV1,NH1)

    call ncd_getvar(grid_nfid, 'CDPTH',ntp,CumDepz2LayBottom_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'BKDSI',ntp,SoiBulkDensityt0_vr(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'FC', ntp,FieldCapacity_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'WP', ntp,WiltPoint_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'SCNV', ntp,SatHydroCondVert_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'SCNH', ntp,SatHydroCondHrzn_vr(1:JZ,NV1,NH1))
    call ReadVanGenuchtenProfile(ntp,NV1,NH1)

    call ncd_getvar(grid_nfid, 'CSAND',ntp,CSAND_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CSILT',ntp,CSILT_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'FHOL',ntp,SoilFracAsMacP_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'ROCK',ntp,ROCK_vr(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'PH',ntp,PH_vr(1:JZ,NV1,NH1))
    !meq/100g means 1.e-3 mol/100 g =1 cmol/kg , cmol ~ 1.e-2 mol
    call ncd_getvar(grid_nfid, 'CEC',ntp,CEC_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'AEC',ntp,AEC_vr(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'CORGC',ntp,CSoilOrgM_vr(ielmc,1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CORGR',ntp,COMLitrC_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CORGN',ntp,CSoilOrgM_vr(ielmn,1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CORGP',ntp,CSoilOrgM_vr(ielmp,1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'CNH4',ntp,CNH4_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CNO3',ntp,CNO3_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CPO4',ntp,CPO4_vr(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'CAL' ,ntp,CAL_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CFE' ,ntp,CFE_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CCA' ,ntp,CCA_vr(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'CMG' ,ntp,CMG_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CNA' ,ntp,CNA_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CKA' ,ntp,CKA_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CSO4',ntp,CSO4_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CCL',ntp,CCL_vr(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'CALPO',ntp,CALPO_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CFEPO',ntp,CFEPO_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CCAPD',ntp,CCAPD_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CCAPH',ntp,CCAPH_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CALOH',ntp,CALOH_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CFEOH',ntp,CFEOH_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CCACO',ntp,CCACO_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'CCASO',ntp,CCASO_vr(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'GKC4',ntp,GKC4_vr(1:JZ,NV1,NH1))  !Gapon coeff for NH4(+)
    call ncd_getvar(grid_nfid, 'GKCH',ntp,GKCH_vr(1:JZ,NV1,NH1))  !Gapon coeff for H(+)
    call ncd_getvar(grid_nfid, 'GKCA',ntp,GKCA_vr(1:JZ,NV1,NH1))  !Gapon coeff for Al(3+)
    call ncd_getvar(grid_nfid, 'GKCM',ntp,GKCM_vr(1:JZ,NV1,NH1))  !Gapon coeff for Mg(2+)
    call ncd_getvar(grid_nfid, 'GKCN',ntp,GKCN_vr(1:JZ,NV1,NH1))  !Gapon coeff for Na(+)
    call ncd_getvar(grid_nfid, 'GKCK',ntp,GKCK_vr(1:JZ,NV1,NH1))  !Gapon coeff for K(+)

    call ncd_getvar(grid_nfid, 'THW',ntp,THW_vr(1:JZ,NV1,NH1))
    call ncd_getvar(grid_nfid, 'THI',ntp,THI_vr(1:JZ,NV1,NH1))

    call ncd_getvar(grid_nfid, 'RSCfL',ntp,dat1(1:JZ));RSC_vr(k_fine_comp,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSNfL',ntp,dat1(1:JZ));RSN_vr(k_fine_comp,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSPfL',ntp,dat1(1:JZ));RSP_vr(k_fine_comp,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSCwL',ntp,dat1(1:JZ));RSC_vr(k_woody_comp,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSNwL',ntp,dat1(1:JZ));RSN_vr(k_woody_comp,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSPwL',ntp,dat1(1:JZ));RSP_vr(k_woody_comp,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSCmL',ntp,dat1(1:JZ));RSC_vr(k_manure,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSNmL',ntp,dat1(1:JZ));RSN_vr(k_manure,1:JZ,NV1,NH1)=dat1(1:JZ)
    call ncd_getvar(grid_nfid, 'RSPmL',ntp,dat1(1:JZ));RSP_vr(k_manure,1:JZ,NV1,NH1)=dat1(1:JZ)

    DO  NX=NH1,NH2
      DO  NY=NV1,NV2
        IF (NX/=NH1 .OR. NY/=NV1) THEN
          !
          !     SURFACE PROPERTIES
          !
          !     PSIAtFldCapacity,PSIAtWiltPoint=water potentials at field capacity,wilting point (MPa)
          !     SoilAlbedo=wet soil albedo
          !     PH=litter pH
          !     RSC,RSC,RSP=C,N,P in fine(1,0),woody(0,0),manure(2,0) surface litter (g m-2)
          !     IXTYP=surface litter type:1=plant,2=manure
          !     NUI,MaxNumRootLays=number of soil surface layer,maximum rooting layer
          !     NL1,NL2=number of additional layers below NJ with,without data in file
          !     ISOILR_col=natural(0),reconstructed(1) soil profile
          !
          PSIAtFldCapacity_col(NY,NX)   = PSIAtFldCapacity_col(NV1,NH1)
          PSIAtWiltPoint_col(NY,NX)     = PSIAtWiltPoint_col(NV1,NH1)
          SoilAlbedo_col(NY,NX)     = SoilAlbedo_col(NV1,NH1)
          PH_vr(0,NY,NX)               = PH_vr(0,NV1,NH1)
          RSC_vr(k_fine_comp,0,NY,NX)  = RSC_vr(k_fine_comp,0,NV1,NH1)
          RSN_vr(k_fine_comp,0,NY,NX)  = RSN_vr(k_fine_comp,0,NV1,NH1)
          RSP_vr(k_fine_comp,0,NY,NX)  = RSP_vr(k_fine_comp,0,NV1,NH1)
          RSC_vr(k_woody_comp,0,NY,NX) = RSC_vr(k_woody_comp,0,NV1,NH1)
          RSN_vr(k_woody_comp,0,NY,NX) = RSN_vr(k_woody_comp,0,NV1,NH1)
          RSP_vr(k_woody_comp,0,NY,NX) = RSP_vr(k_woody_comp,0,NV1,NH1)
          RSC_vr(k_manure,0,NY,NX)     = RSC_vr(k_manure,0,NV1,NH1)
          RSN_vr(k_manure,0,NY,NX)     = RSN_vr(k_manure,0,NV1,NH1)
          RSP_vr(k_manure,0,NY,NX)     = RSP_vr(k_manure,0,NV1,NH1)
          iLitrType_col(1,NY,NX)            = iLitrType_col(1,NV1,NH1)
          iLitrType_col(2,NY,NX)            = iLitrType_col(2,NV1,NH1)
          NUI_col(NY,NX)                = NUI_col(NV1,NH1)
          MaxNumRootLays_col(NY,NX)     = MaxNumRootLays_col(NV1,NH1)
          ISOILR_col(NY,NX)             = ISOILR_col(NV1,NH1)
          NU_col(NY,NX)                 =NU_col(NV1,NH1)
          NK_col(NY,NX)                 = NK_col(NV1,NH1)
          NM(NY,NX)                 = NM(NV1,NH1)
          !  the extra soil layer below root zone cannot be greater than what is allowed
          NLI_col(NY,NX) = NLI_col(NV1,NH1)
          NL_col(NY,NX)  = NLI_col(NV1,NH1)
        ENDIF

        !
        !     PHYSICAL PROPERTIES
        !
        !     CDPTH=depth to bottom (m) > 0
        !     SoiBulkDensityt0_vr=initial bulk density (Mg m-3,0=water), it refers to solid matter
        !
        !
        iPondFlag_col(NY,NX)=.false.
        IF (NX/=NH1 .OR. NY/=NV1) THEN
          DO L=NU_col(NY,NX),NM(NY,NX)
            CumDepz2LayBottom_vr(L,NY,NX) = CumDepz2LayBottom_vr(L,NV1,NH1)
            SoiBulkDensityt0_vr(L,NY,NX)  = SoiBulkDensityt0_vr(L,NV1,NH1)
            if(.not.iPondFlag_col(NY,NX))iPondFlag_col(NY,NX)=SoiBulkDensityt0_vr(L,NY,NX).LE.ZERO
            if(SoiBulkDensityt0_vr(L,NY,NX).LE.ZERO)iPondBotLev_col(NY,NX)=L
            FieldCapacity_vr(L,NY,NX)     = FieldCapacity_vr(L,NV1,NH1)
            WiltPoint_vr(L,NY,NX)         = WiltPoint_vr(L,NV1,NH1)
            if(allocated(VGThetaR_vr))then
              VGThetaR_vr(L,NY,NX)=VGThetaR_vr(L,NV1,NH1)
              VGAlpha_vr(L,NY,NX)=VGAlpha_vr(L,NV1,NH1)
              VGN_vr(L,NY,NX)=VGN_vr(L,NV1,NH1)
              VGL_vr(L,NY,NX)=VGL_vr(L,NV1,NH1)
            endif
            SatHydroCondVert_vr(L,NY,NX)  = SatHydroCondVert_vr(L,NV1,NH1)
            SatHydroCondHrzn_vr(L,NY,NX)  = SatHydroCondHrzn_vr(L,NV1,NH1)
            CSAND_vr(L,NY,NX)             = CSAND_vr(L,NV1,NH1)
            CSILT_vr(L,NY,NX)              = CSILT_vr(L,NV1,NH1)
            SoilFracAsMacP_vr(L,NY,NX)    = SoilFracAsMacP_vr(L,NV1,NH1)
            ROCK_vr(L,NY,NX)              = ROCK_vr(L,NV1,NH1)
            PH_vr(L,NY,NX)                = PH_vr(L,NV1,NH1)
            CEC_vr(L,NY,NX)               = CEC_vr(L,NV1,NH1)
            AEC_vr(L,NY,NX)               = AEC_vr(L,NV1,NH1)
            CSoilOrgM_vr(ielmc,L,NY,NX)   = CSoilOrgM_vr(ielmc,L,NV1,NH1)
            COMLitrC_vr(L,NY,NX)          = COMLitrC_vr(L,NV1,NH1)
            CSoilOrgM_vr(ielmn,L,NY,NX)   = CSoilOrgM_vr(ielmn,L,NV1,NH1)
            CSoilOrgM_vr(ielmp,L,NY,NX)   = CSoilOrgM_vr(ielmp,L,NV1,NH1)

            CNH4_vr(L,NY,NX) = CNH4_vr(L,NV1,NH1)
            CNO3_vr(L,NY,NX) = CNO3_vr(L,NV1,NH1)
            CPO4_vr(L,NY,NX) = CPO4_vr(L,NV1,NH1)

            CAL_vr(L,NY,NX)  = CAL_vr(L,NV1,NH1)
            CFE_vr(L,NY,NX)  = CFE_vr(L,NV1,NH1)
            CCA_vr(L,NY,NX)  = CCA_vr(L,NV1,NH1)
            CMG_vr(L,NY,NX)  = CMG_vr(L,NV1,NH1)
            CNA_vr(L,NY,NX)  = CNA_vr(L,NV1,NH1)
            CKA_vr(L,NY,NX)  = CKA_vr(L,NV1,NH1)
            CSO4_vr(L,NY,NX) = CSO4_vr(L,NV1,NH1)
            CCL_vr(L,NY,NX)  = CCL_vr(L,NV1,NH1)

            CALPO_vr(L,NY,NX) = CALPO_vr(L,NV1,NH1)
            CFEPO_vr(L,NY,NX) = CFEPO_vr(L,NV1,NH1)
            CCAPD_vr(L,NY,NX) = CCAPD_vr(L,NV1,NH1)
            CCAPH_vr(L,NY,NX) = CCAPH_vr(L,NV1,NH1)
            CALOH_vr(L,NY,NX) = CALOH_vr(L,NV1,NH1)
            CFEOH_vr(L,NY,NX) = CFEOH_vr(L,NV1,NH1)
            CCACO_vr(L,NY,NX) = CCACO_vr(L,NV1,NH1)
            CCASO_vr(L,NY,NX) = CCASO_vr(L,NV1,NH1)

            GKC4_vr(L,NY,NX) = GKC4_vr(L,NV1,NH1)
            GKCH_vr(L,NY,NX) = GKCH_vr(L,NV1,NH1)
            GKCA_vr(L,NY,NX) = GKCA_vr(L,NV1,NH1)
            GKCM_vr(L,NY,NX) = GKCM_vr(L,NV1,NH1)
            GKCN_vr(L,NY,NX) = GKCN_vr(L,NV1,NH1)
            GKCK_vr(L,NY,NX) = GKCK_vr(L,NV1,NH1)

            THW_vr(L,NY,NX) = THW_vr(L,NV1,NH1)
            THI_vr(L,NY,NX) = THI_vr(L,NV1,NH1)

            RSC_vr(k_fine_comp,L,NY,NX)  = RSC_vr(k_fine_comp,L,NV1,NH1)
            RSN_vr(k_fine_comp,L,NY,NX)  = RSN_vr(k_fine_comp,L,NV1,NH1)
            RSP_vr(k_fine_comp,L,NY,NX)  = RSP_vr(k_fine_comp,L,NV1,NH1)
            RSC_vr(k_woody_comp,L,NY,NX) = RSC_vr(k_woody_comp,L,NV1,NH1)
            RSN_vr(k_woody_comp,L,NY,NX) = RSN_vr(k_woody_comp,L,NV1,NH1)
            RSP_vr(k_woody_comp,L,NY,NX) = RSP_vr(k_woody_comp,L,NV1,NH1)
            RSC_vr(k_manure,L,NY,NX)     = RSC_vr(k_manure,L,NV1,NH1)
            RSN_vr(k_manure,L,NY,NX)     = RSN_vr(k_manure,L,NV1,NH1)
            RSP_vr(k_manure,L,NY,NX)     = RSP_vr(k_manure,L,NV1,NH1)
          ENDDO
        ENDIF
                
        DO L=1,NL_col(NY,NX)

          if(CSoilOrgM_vr(ielmc,L,NY,NX) > 0._r8)then 
            if(CSoilOrgM_vr(ielmn,L,NY,NX)*1.e-3_r8>CSoilOrgM_vr(ielmc,L,NY,NX))then
              write(iulog,*)'Likely too larger N/C ratio',1.e-3_r8*safe_adb(CSoilOrgM_vr(ielmn,L,NY,NX),CSoilOrgM_vr(ielmc,L,NY,NX)), 'in L,NY,NX',L,NY,NX
              call endrun(trim(mod_filename)//' at line',__LINE__)  
            endif
            if(CSoilOrgM_vr(ielmp,L,NY,NX)>CSoilOrgM_vr(ielmn,L,NY,NX))then
              write(iulog,*)'Likely too larger P/N ratio',safe_adb(CSoilOrgM_vr(ielmp,L,NY,NX),CSoilOrgM_vr(ielmn,L,NY,NX)), 'in L,NY,NX',L,NY,NX
              call endrun(trim(mod_filename)//' at line',__LINE__)  
            endif
          endif  
        ENDDO

        if(lverb)then
          CALL Disp_topo_charc(NY,NX,NU_col(NY,NX),NM(NY,NX))
        endif
        RSC_vr(k_fine_comp,0,NY,NX)     = AMAX1(1.00E-06_r8,RSC_vr(k_fine_comp,0,NY,NX))
        RSN_vr(k_fine_comp,0,NY,NX)     = AMAX1(0.04E-06_r8,RSN_vr(k_fine_comp,0,NY,NX))
        RSP_vr(k_fine_comp,0,NY,NX)     = AMAX1(0.004E-06_r8,RSP_vr(k_fine_comp,0,NY,NX))

        SatHydroCondVert_vr(0,NY,NX) = 10.0_r8*0.098_r8
        !
        !     SET FLAGS FOR ESTIMATING FC,WP,SCNV,SCNH IF UNKNOWN
        !
        !     ISOIL=flag for calculating FC(1),WiltPoint_vr(2),SatHydroCondVert_vr(3),SatHydroCondHrzn_vr(4)
        !
        call ComputeSoilHydroPars(NY,NX,NU_col(NY,NX),NM(NY,NX))

        !
        !     FILL OUT SOIL BOUNDARY LAYERS ABOVE ROOTING ZONE (NOT USED)
        !     below is for soil repacking, whence NU>1
        !     root zone from NU to NM, why use 0.025?
        !
        IF(NU_col(NY,NX).GT.1)THEN
          DO  L=NU_col(NY,NX)-1,0,-1
            IF(SoiBulkDensityt0_vr(L+1,NY,NX).GT.0.025_r8)THEN 
              !next layer is not water
              CumDepz2LayBottom_vr(L,NY,NX)=CumDepz2LayBottom_vr(L+1,NY,NX)-0.01_r8
              !next layer is Likely water
            ELSE
              CumDepz2LayBottom_vr(L,NY,NX)=CumDepz2LayBottom_vr(L+1,NY,NX)-0.02_r8     !place holder for water
            ENDIF

            IF(L.GT.0)THEN
              SoiBulkDensityt0_vr(L,NY,NX) = SoiBulkDensityt0_vr(L+1,NY,NX)
              FieldCapacity_vr(L,NY,NX)    = FieldCapacity_vr(L+1,NY,NX)
              WiltPoint_vr(L,NY,NX)        = WiltPoint_vr(L+1,NY,NX)
              if(allocated(VGThetaR_vr))then
                VGThetaR_vr(L,NY,NX)=VGThetaR_vr(L+1,NY,NX)
                VGAlpha_vr(L,NY,NX)=VGAlpha_vr(L+1,NY,NX)
                VGN_vr(L,NY,NX)=VGN_vr(L+1,NY,NX)
                VGL_vr(L,NY,NX)=VGL_vr(L+1,NY,NX)
              endif
              SatHydroCondVert_vr(L,NY,NX) = SatHydroCondVert_vr(L+1,NY,NX)
              SatHydroCondHrzn_vr(L,NY,NX) = SatHydroCondHrzn_vr(L+1,NY,NX)
              CSAND_vr(L,NY,NX)            = CSAND_vr(L+1,NY,NX)
              CSILT_vr(L,NY,NX)            = CSILT_vr(L+1,NY,NX)
              CCLAY_vr(L,NY,NX)            = CCLAY_vr(L+1,NY,NX)
              SoilFracAsMacP_vr(L,NY,NX)   = SoilFracAsMacP_vr(L+1,NY,NX)
              ROCK_vr(L,NY,NX)             = ROCK_vr(L+1,NY,NX)
              PH_vr(L,NY,NX)               = PH_vr(L+1,NY,NX)
              CEC_vr(L,NY,NX)              = CEC_vr(L+1,NY,NX)
              AEC_vr(L,NY,NX)              = AEC_vr(L+1,NY,NX)
              CSoilOrgM_vr(ielmc,L,NY,NX)  = 1.0_r8*CSoilOrgM_vr(ielmc,L+1,NY,NX)
              COMLitrC_vr(L,NY,NX)         = 1.0_r8*COMLitrC_vr(L+1,NY,NX)
              CSoilOrgM_vr(ielmn,L,NY,NX)  = 1.0_r8*CSoilOrgM_vr(ielmn,L+1,NY,NX)
              CSoilOrgM_vr(ielmp,L,NY,NX)  = 1.0_r8*CSoilOrgM_vr(ielmp,L+1,NY,NX)
              CNH4_vr(L,NY,NX)             = CNH4_vr(L+1,NY,NX)
              CNO3_vr(L,NY,NX)             = CNO3_vr(L+1,NY,NX)
              CPO4_vr(L,NY,NX)             = CPO4_vr(L+1,NY,NX)
              CAL_vr(L,NY,NX)              = CAL_vr(L+1,NY,NX)
              CFE_vr(L,NY,NX)              = CFE_vr(L+1,NY,NX)
              CCA_vr(L,NY,NX)              = CCA_vr(L+1,NY,NX)
              CMG_vr(L,NY,NX)              = CMG_vr(L+1,NY,NX)
              CNA_vr(L,NY,NX)              = CNA_vr(L+1,NY,NX)
              CKA_vr(L,NY,NX)              = CKA_vr(L+1,NY,NX)
              CSO4_vr(L,NY,NX)             = CSO4_vr(L+1,NY,NX)
              CCL_vr(L,NY,NX)              = CCL_vr(L+1,NY,NX)
              CALOH_vr(L,NY,NX)            = CALOH_vr(L+1,NY,NX)
              CFEOH_vr(L,NY,NX)            = CFEOH_vr(L+1,NY,NX)
              CCACO_vr(L,NY,NX)            = CCACO_vr(L+1,NY,NX)
              CCASO_vr(L,NY,NX)            = CCASO_vr(L+1,NY,NX)
              CALPO_vr(L,NY,NX)            = CALPO_vr(L+1,NY,NX)
              CFEPO_vr(L,NY,NX)            = CFEPO_vr(L+1,NY,NX)
              CCAPD_vr(L,NY,NX)            = CCAPD_vr(L+1,NY,NX)
              CCAPH_vr(L,NY,NX)            = CCAPH_vr(L+1,NY,NX)
              GKC4_vr(L,NY,NX)             = GKC4_vr(L+1,NY,NX)
              GKCH_vr(L,NY,NX)             = GKCH_vr(L+1,NY,NX)
              GKCA_vr(L,NY,NX)             = GKCA_vr(L+1,NY,NX)
              GKCM_vr(L,NY,NX)             = GKCM_vr(L+1,NY,NX)
              GKCN_vr(L,NY,NX)             = GKCN_vr(L+1,NY,NX)
              GKCK_vr(L,NY,NX)             = GKCK_vr(L+1,NY,NX)
              THW_vr(L,NY,NX)              = THW_vr(L+1,NY,NX)
              THI_vr(L,NY,NX)              = THI_vr(L+1,NY,NX)
              ISOIL_vr(1:4,L,NY,NX)           = ISOIL_vr(1:4,L+1,NY,NX)
              RSC_vr(k_fine_comp,L,NY,NX)     = 0.0_r8
              RSN_vr(k_fine_comp,L,NY,NX)     = 0.0_r8
              RSP_vr(k_fine_comp,L,NY,NX)     = 0.0_r8
              RSC_vr(k_woody_comp,L,NY,NX)    = 0.0_r8
              RSN_vr(k_woody_comp,L,NY,NX)    = 0.0_r8
              RSP_vr(k_woody_comp,L,NY,NX)    = 0.0_r8
              RSC_vr(k_manure,L,NY,NX)        = 0.0_r8
              RSN_vr(k_manure,L,NY,NX)        = 0.0_r8
              RSP_vr(k_manure,L,NY,NX)        = 0.0_r8
            ENDIF
          ENDDO
        ENDIF
        !
        !     ADD SOIL BOUNDARY LAYERS BELOW SOIL ZONE
        !     depth of layer (L-1) is at the middle between that of layer L-2 and L
        if(lverb)write(*,*)'SetDeepSoil'
        call SetDeepSoil(NY,NX,NM(NY,NX),JZ)

        !
        !   CALCULATE DERIVED SOIL PROPERTIES FROM INPUT SOIL PROPERTIES
        !
        !   FracSoiAsMicP_vr=micropore fraction excluding macropore,rock
        !   SCNV,SCNH=vertical,lateral Ksat converted to m2 MPa-1 h-1
        !   CSAND_vr,CSILT,CCLAY_vr=sand,silt,clay content converted to g Mg-1
        !   CORGC,CORGR=SOC,POC converted to g Mg-1
        !   CEC,AEC=cation,anion exchange capacity converted to mol Mg-1
        !   CNH4...=solute concentrations converted to mol Mg-1
        !   SoiBulkDensityt0_vr: initial bulk density
        
        DO  L=1,NL_col(NY,NX)
        !   SoilFracAsMacP: macropore fraction
  !     SoiBulkDensityt0_vr(L,NY,NX)=SoiBulkDensityt0_vr(L,NY,NX)/(1.0_r8-SoilFracAsMacP_vr(L,NY,NX))
          SoilBulkDensity_vr(L,NY,NX)=SoiBulkDensityt0_vr(L,NY,NX)
          
          IF(isclose(SoilBulkDensity_vr(L,NY,NX),0.0_r8))SoilFracAsMacP_vr(L,NY,NX)=0.0_r8
          ! fraction of soil as micropore
          FracSoiAsMicP_vr(L,NY,NX)=(1.0_r8-ROCK_vr(L,NY,NX))*(1.0_r8-SoilFracAsMacP_vr(L,NY,NX))
        !  Macropore correction is off, when reporting from measurements, FieldCapacity includes contribution from
        !  both macropores and micropores    
  !     FieldCapacity_vr(L,NY,NX)=FieldCapacity_vr(L,NY,NX)/(1.0-SoilFracAsMacP_vr(L,NY,NX))
  !     WiltPoint_vr(L,NY,NX)=WiltPoint_vr(L,NY,NX)/(1.0-SoilFracAsMacP_vr(L,NY,NX))
  !
          SatHydroCondVert_vr(L,NY,NX) = 0.098_r8*SatHydroCondVert_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
          SatHydroCondHrzn_vr(L,NY,NX) = 0.098_r8*SatHydroCondHrzn_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
          CCLAY_vr(L,NY,NX)            = AZMAX1(1.0E+03_r8-(CSAND_vr(L,NY,NX)+CSILT_vr(L,NY,NX)))
          !convert from Kg to g C (C is input as kgC/Mg soil, N and P are input as g/Mg soil)
          CSoilOrgM_vr(ielmc,L,NY,NX)  = CSoilOrgM_vr(ielmc,L,NY,NX)*1.0E+03_r8  !from kg C to g C
          COMLitrC_vr(L,NY,NX)         = COMLitrC_vr(L,NY,NX)*1.0E+03_r8         !from kg C to g C
          CORGCI_vr(L,NY,NX)           = CSoilOrgM_vr(ielmc,L,NY,NX)
          SoilFracAsMacPt0_vr(L,NY,NX) = SoilFracAsMacP_vr(L,NY,NX)
          ! soil texture is reported based on mass basis soley for mineral component of the soil
          !volume of organic matter
          OrgVolFrac = CSoilOrgM_vr(ielmc,L,NY,NX)/orgcden
          corrector  = 1.0E-03_r8*AZMAX1((1.0_r8-OrgVolFrac))
          
          !convert soil texture into mass scale [0,1]
          CSAND_vr(L,NY,NX) = CSAND_vr(L,NY,NX)*corrector 
          CSILT_vr(L,NY,NX) = CSILT_vr(L,NY,NX)*corrector
          CCLAY_vr(L,NY,NX) = CCLAY_vr(L,NY,NX)*corrector
          CEC_vr(L,NY,NX)   = CEC_vr(L,NY,NX)*10.0_r8   !convert from meq/100g to cmol/kg
          AEC_vr(L,NY,NX)   = AEC_vr(L,NY,NX)*10.0_r8   !convert from meq/100g to cmol/kg
          CNH4_vr(L,NY,NX)  = CNH4_vr(L,NY,NX)/natomw
          CNO3_vr(L,NY,NX)  = CNO3_vr(L,NY,NX)/natomw
          CPO4_vr(L,NY,NX)  = CPO4_vr(L,NY,NX)/patomw
          CAL_vr(L,NY,NX)   = CAL_vr(L,NY,NX)/27.0_r8
          CFE_vr(L,NY,NX)   = CFE_vr(L,NY,NX)/56.0_r8
          CCA_vr(L,NY,NX)   = CCA_vr(L,NY,NX)/40.0_r8
          CMG_vr(L,NY,NX)   = CMG_vr(L,NY,NX)/24.3_r8
          CNA_vr(L,NY,NX)   = CNA_vr(L,NY,NX)/23.0_r8
          CKA_vr(L,NY,NX)   = CKA_vr(L,NY,NX)/39.1_r8
          CSO4_vr(L,NY,NX)  = CSO4_vr(L,NY,NX)/32.0_r8
          CCL_vr(L,NY,NX)   = CCL_vr(L,NY,NX)/35.5_r8
          CALPO_vr(L,NY,NX) = CALPO_vr(L,NY,NX)/patomw
          CFEPO_vr(L,NY,NX)=CFEPO_vr(L,NY,NX)/patomw
          CCAPD_vr(L,NY,NX)=CCAPD_vr(L,NY,NX)/patomw
          CCAPH_vr(L,NY,NX)=CCAPH_vr(L,NY,NX)/(patomw*3.0_r8)
          CALOH_vr(L,NY,NX)=CALOH_vr(L,NY,NX)/27.0_r8
          CFEOH_vr(L,NY,NX)=CFEOH_vr(L,NY,NX)/56.0_r8
          CCACO_vr(L,NY,NX)=CCACO_vr(L,NY,NX)/40.0_r8
          CCASO_vr(L,NY,NX)=CCASO_vr(L,NY,NX)/40.0_r8
          !
          !     ESTIMATE SON,SOP,CEC IF UNKNOWN
          ! Tipping et al (2008)   The C:N:P:S stoichiometry of soil organic matter, BIOCHEMISTRY 130:117-131
          !
          IF(CSoilOrgM_vr(ielmn,L,NY,NX).LT.0.0_r8)THEN
            !  default ORGN parameterization
            CSoilOrgM_vr(ielmn,L,NY,NX)=AMIN1(0.125_r8*CSoilOrgM_vr(ielmc,L,NY,NX),&
              8.9E+02_r8*(CSoilOrgM_vr(ielmc,L,NY,NX)/1.0E+04_r8)**0.80_r8)
          ENDIF
          IF(CSoilOrgM_vr(ielmp,L,NY,NX).LT.0.0_r8)THEN
            CSoilOrgM_vr(ielmp,L,NY,NX)=AMIN1(0.0125_r8*CSoilOrgM_vr(ielmc,L,NY,NX),&
              1.2E+02_r8*(CSoilOrgM_vr(ielmc,L,NY,NX)/1.0E+04_r8)**0.52_r8)
          ENDIF
          IF(CEC_vr(L,NY,NX).LT.0.0_r8)THEN
            !estimate from input data
            CEC_vr(L,NY,NX)=10.0_r8*(200.0_r8*2.0_r8*CSoilOrgM_vr(ielmc,L,NY,NX)/1.0E+06_r8 &
              +80.0_r8*CCLAY_vr(L,NY,NX)+20.0_r8*CSILT_vr(L,NY,NX) &
              +5.0_r8*CSAND_vr(L,NY,NX))
          ENDIF
        ENDDO
        CSoilOrgM_vr(ielmc,0,NY,NX)=orgcden
        FracSoiAsMicP_vr(0,NY,NX)=1.0_r8
      ENDDO
    ENDDO
  ENDDO
  Call PrintInfo('end '//subname)
  end associate
  end subroutine readTopoNC

!------------------------------------------------------------------------------------------

  subroutine Disp_topo_charc(NY,NX,NU,NM)

  implicit none
  integer, intent(in) :: NY,NX,NU,NM
  integer :: LL,L

  associate(                            &
  k_woody_comp => micpar%k_woody_comp , &
  k_fine_comp  => micpar%k_fine_comp  , &
  k_manure     => micpar%k_manure       &
  )

  write(*,'(40A)')('-',ll=1,40)
  write(*,*)'Topographic characterization'
  write(*,'(40A)')('-',ll=1,40)
  write(*,*)'NY, NX =',NY,NX
  write(*,*)'Aspect (o): ASPX',ASP_col(NY,NX)
  write(*,*)'Slope (o): SL0',SL_col(NY,NX)
  write(*,*)'Initial snowpack depth: initSnowDepth',SnowDepth_col(NY,NX)
  write(*,'(100A)')('=',ll=1,100)

  write(*,*)''
  write(*,*)'NY,NX=',NY,NX
  write(*,*)'Water potential at field capacity (MPa)',PSIAtFldCapacity_col(NY,NX)
  write(*,*)'Water potential at wilting point (MPa)',PSIAtWiltPoint_col(NY,NX)
  write(*,*)'Wet soil albedo',SoilAlbedo_col(NY,NX)

  write(*,*)'Litter pH',PH_vr(0,NY,NX)
  write(*,*)'C in surface fine litter (g m-2)',RSC_vr(k_fine_comp,0,NY,NX)
  write(*,*)'N in surface fine litter (g m-2)',RSN_vr(k_fine_comp,0,NY,NX)
  write(*,*)'P in surface fine litter (g m-2)',RSP_vr(k_fine_comp,0,NY,NX)
  write(*,*)'C in surface woody litter (g m-2)',RSC_vr(k_woody_comp,0,NY,NX)
  write(*,*)'N in surface woody litter (g m-2)',RSN_vr(k_woody_comp,0,NY,NX)
  write(*,*)'P in surface woody litter (g m-2)',RSP_vr(k_woody_comp,0,NY,NX)
  write(*,*)'C in surface manure litter (g m-2)',RSC_vr(k_manure,0,NY,NX)
  write(*,*)'N in surface manure litter (g m-2)',RSN_vr(k_manure,0,NY,NX)
  write(*,*)'P in surface manure litter (g m-2)',RSP_vr(k_manure,0,NY,NX)
  write(*,*)'surface litter type:1=plant,2=manure',iLitrType_col(1,NY,NX),iLitrType_col(2,NY,NX)
  write(*,*)'layer number of soil surface layer NUI',NUI_col(NY,NX)
  write(*,*)'layer number of maximum rooting layer NJ',MaxNumRootLays_col(NY,NX)
  write(*,*)'number of layers involved in model calculation',NL_col(NY,NX)
  write(*,*)'Flag for natural(0),reconstructed(1) soil profile', ISOILR_col(NY,NX)
  write(*,*)
  write(*,'(A,I2,A,I2)')'read data for layers from layer NU ',NU,' to layer NM ',NM

  write(*,*)'Depth to bottom of soil layer (m): CDPTH'
  write(*,*)(CumDepz2LayBottom_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Initial bulk density (Mg m-3, 0=water): SoiBulkDensityt0_vr'
  write(*,*)(SoiBulkDensityt0_vr(L,NY,NX),L=NU,NM)
!
!     HYDROLOGIC PROPERTIES
!
!     FC,WP=field capacity,wilting point:<0=unknown (m3 m-3)
!     SCNV,SCNH=vertical,lateral Ksat:<0=unknown (mm h-1)
!

  write(*,*)''
  write(*,*)'NY,NX=',NY,NX
  write(*,*)'L=',NU,NM
  write(*,*)'Field capacity (m3 m-3): FC'
  write(*,*)(FieldCapacity_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Wilting point (m3 m-3): WP'
  write(*,*)(WiltPoint_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Vertical Ksat (mm h-1): SCNV'
  write(*,*)(SatHydroCondVert_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Lateral Ksat (mm h-1): SCNH'
  write(*,*)(SatHydroCondHrzn_vr(L,NY,NX),L=NU,NM)
!
!     PHYSICAL PROPERTIES
!
!     CSAND_vr,CSILT=sand,silt contents (kg Mg-1)
!     SoilFracAsMacP_vr,ROCK=macropore,rock fraction
!

  write(*,*)''
  write(*,*)'Sand (kg Mg-1): CSAND_vr'
  write(*,*)(CSAND_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Silt (kg Mg-1): CSILT'
  write(*,*)(CSILT_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Macropore fraction (0-1): FOHL'
  write(*,*)(SoilFracAsMacP_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Rock fraction (0-1): ROCK'
  write(*,*)(ROCK_vr(L,NY,NX),L=NU,NM)
!
!     CHEMICAL PROPERTIES
!
!     PH=pH
!     CEC,AEC=cation,anion exchange capacity:CEC<0=unknown (cmol Kg-1)
!

  write(*,*)''
  write(*,*)'pH'
  write(*,*)(PH_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Cation exchange capacity (cmol Kg-1): CEC'
  write(*,*)(CEC_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Anion exchange capacity (cmol Kg-1): AEC'
  write(*,*)(AEC_vr(L,NY,NX),L=NU,NM)
!
!     ORGANIC C, N AND P CONCENTRATIONS
!
!     CORGC,COMLitrC_vr=total SOC,POC(part of SOC) (kg Mg-1)
!     CORGN,CORGP=SON,SOP:<0=unknown (g Mg-1)
!
  write(*,*)''
  write(*,*)'Total SOC (kg C/Mg soil): CORGC'
  write(*,*)(CSoilOrgM_vr(ielmc,L,NY,NX),L=NU,NM)
  write(*,*)'POC (part of SOC) (kg C/Mg soil): CORGR '
  write(*,*)(COMLitrC_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Total SON (g N/Mg soil): CORGN '
  write(*,*)(CSoilOrgM_vr(ielmn,L,NY,NX),L=NU,NM)
  write(*,*)'Total SOP (g P/Mg soil): CORGP'
  write(*,*)(CSoilOrgM_vr(ielmp,L,NY,NX),L=NU,NM)
!
!     INORGANIC N AND P CONCENTRATIONS
!
!     CNH4,CNO3,CPO4=soluble+exchangeable NH4,NO3,H2PO4 (g Mg-1)
!
  write(*,*)''
  write(*,*)'Total soil NH4 concentration (gN Mg-1): CNH4 '
  write(*,*)(CNH4_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Total soil NO3 concentration (gN Mg-1): CNO3'
  write(*,*)(CNO3_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Total soil H2PO4 concentration (gP Mg-1): CPO4 '
  write(*,*)(CPO4_vr(L,NY,NX),L=NU,NM)
!
!     CATION AND ANION CONCENTRATIONS
!
!     C*=soluble concentration from sat. paste extract (g Mg-1)
!     AL,FE,CA,MG,NA,KA,SO4,CL=Al,Fe,Ca,Mg,Na,K,SO4-S,Cl
!
  write(*,*)''
  write(*,*)'Soluble soil Al content (g Mg-1): CAL'
  write(*,*)(CAL_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soluble soil Fe content (g Mg-1): CFE'
  write(*,*)(CFE_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soluble soil Ca content (g Mg-1): CCA'
  write(*,*)(CCA_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soluble soil Na content (g Mg-1): CNA'
  write(*,*)(CNA_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soluble soil K content (g Mg-1): CKA'
  write(*,*)(CKA_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soluble soil SO4 content (gS Mg-1): CSO4'
  write(*,*)(CSO4_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soluble soil Cl content (g Mg-1): CCL'
  write(*,*)(CCL_vr(L,NY,NX),L=NU,NM)
!
!     PRECIPITATED MINERAL CONCENTRATIONS
!
!     CALPO,CFEPO,CCAPD,CCAPH=AlPO4,FePO4,CaHPO4,apatite (g Mg-1)
!     CALOH,CFEOH,CCACO,CCASO=AlOH3,FeOH3,CaSO4,CaCO3 (g Mg-1)
!
  write(*,*)''
  write(*,*)'Soil AlPO4 content (g Mg-1): CALPO'
  write(*,*)(CALPO_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soil FePO4 content (g Mg-1): CFEPO'
  write(*,*)(CFEPO_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soil CaHPO4 content (g Mg-1): CCAPD'
  write(*,*)(CCAPD_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soil apatite content (g Mg-1): CCAPH '
  write(*,*)(CCAPH_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soil Al(OH)3 content (g Mg-1): CALOH '
  write(*,*)(CALOH_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soil Fe(OH)3 content (g Mg-1): CFEOH '
  write(*,*)(CFEOH_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soil CaCO3 content (g Mg-1): CCACO'
  write(*,*)(CCACO_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Soil CaSO4 content (g Mg-1): CCASO'
  write(*,*)(CCASO_vr(L,NY,NX),L=NU,NM)
!
!     GAPON SELECTIVITY CO-EFFICIENTS
!
!     GKC4,GKCH,GKCA,GKCM,GKCN,GKCK=Gapon selectivity coefficients for
!     Ca-NH4,Ca-H,Ca-Al,Ca-Mg,Ca-Na,Ca-K
!
  write(*,*)''
  write(*,*)'Ca-NH4 Gapon selectivity coefficient: GKC4'
  write(*,*)(GKC4_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Ca-H Gapon selectivity coefficient: GKCH'
  write(*,*)(GKCH_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Ca-Al Gapon selectivity coefficient: GKCA'
  write(*,*)(GKCA_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Ca-Mg Gapon selectivity coefficient: GKCM'
  write(*,*)(GKCM_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Ca-Na Gapon selectivity coefficient: GKCN'
  write(*,*)(GKCN_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Ca-K Gapon selectivity coefficient: GKCK'
  write(*,*)(GKCK_vr(L,NY,NX),L=NU,NM)

!
!     INITIAL WATER, ICE CONTENTS
!
  write(*,*)'Initial soil water content (m3/m3): THW'
  write(*,*)(THW_vr(L,NY,NX),L=NU,NM)
  write(*,*)'Initial soil ice content (m3/m3): THI'
  write(*,*)(THI_vr(L,NY,NX),L=NU,NM)
!
!     THW,THI=initial water,ice:>1=satd,1=FC,0=WP,<0=0,0-1=m3 m-3
!
!     INITIAL PLANT AND ANIMAL RESIDUE C, N AND P
!
!     RSC,RSC,RSP=C,N,P in fine(1),woody(0),manure(2) litter (g m-2)
!
  write(*,*)''
  write(*,*)'Initial fine litter C (gC m-2): RSC_vr(k_fine_comp)'
  write(*,*)(RSC_vr(k_fine_comp,L,NY,NX),L=NU,NM)
  write(*,*)'Initial fine litter N (gN m-2): RSN_vr(k_fine_comp)'
  write(*,*)(RSN_vr(k_fine_comp,L,NY,NX),L=NU,NM)
  write(*,*)'Initial fine litter P (gP m-2): RSP_vr(k_fine_comp)'
  write(*,*)(RSP_vr(k_fine_comp,L,NY,NX),L=NU,NM)
  write(*,*)'Initial woody liter C (gC m-2): RSC_vr(k_woody_comp)'
  write(*,*)(RSC_vr(k_woody_comp,L,NY,NX),L=NU,NM)
  write(*,*)'Initial woody litter N (gN m-2): RSN_vr(k_woody_comp)'
  write(*,*)(RSN_vr(k_woody_comp,L,NY,NX),L=NU,NM)
  write(*,*)'Initial woody litter P (gP m-2): RSP_vr(k_woody_comp)'
  write(*,*)(RSP_vr(k_woody_comp,L,NY,NX),L=NU,NM)
  write(*,*)'Initial manure liter C (gC m-2): RSC_vr(k_manure)'
  write(*,*)(RSC_vr(k_manure,L,NY,NX),L=NU,NM)
  write(*,*)'Initial manure litter N (gN m-2): RSN_vr(k_manure)'
  write(*,*)(RSN_vr(k_manure,L,NY,NX),L=NU,NM)
  write(*,*)'Initial manure litter P (gP m-2): RSP_vr(k_manure)'
  write(*,*)(RSP_vr(k_manure,L,NY,NX),L=NU,NM)
  write(*,'(100A)')('=',ll=1,100)

  end associate
  end subroutine Disp_topo_charc

end module readiMod
