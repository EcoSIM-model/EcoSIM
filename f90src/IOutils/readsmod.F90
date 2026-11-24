module readsmod

  use data_kind_mod, only: r8 => DAT_KIND_R8
  use abortutils,    only: endrun
  use fileUtil,      only: open_safe
  use minimathmod,   only: isLeap, AZMAX1
  use EcoSIMCtrlMod, only: lverb,  salt_model, fixClime
  use FireMod
  use DebugToolMod
  use GridConsts
  use FlagDataType
  use FertilizerDataType
  use ClimForcDataType
  use SoilWaterDataType
  use LandSurfDataType
  use EcoSIMCtrlDataType
  use EcosimConst
  use EcoSIMHistMod
  use IrrigationDataType
  use GridDataType
  use EcoSIMConfig
  use ClimReadMod
  use AqueChemDatatype
  implicit none
  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: ReadClimSoilForcing
  integer, save :: yearhi  =0

  type, public :: clime_type
    real(r8) :: AirT_C      !air temperature, C
    real(r8) :: vap_Kpa     !vapor pressure, kPa
    real(r8) :: Wind_ms     !wind speed, m/s
    real(r8) :: Atm_kPa     !atmospheric vapor pressure kPa
    real(r8) :: Rain_mmhr   !rainfall mm/m2/hr
    real(r8) :: SRAD_Wm2    !solar radiation, W m-2
  end type clime_type

  type(clime_type), public :: clim_var
  contains

!------------------------------------------------------------------------------------------

  subroutine SetConstClimeForcing()

  implicit none
  
  IWTHR=2          !set as hourly forcing
  TMP_hrly(:,:)   = clim_var%airT_C
  WINDH(:,:)      = clim_var%Wind_ms * 3600._r8    ![m/s]->[m/hour]
  DWPTH(:,:)      = clim_var%vap_kpa               !
  RAINH(:,:)      = clim_var%Rain_mmhr * 1.e-3_r8  ![mm/hr]->[m/hr]
  SWRad_hrly(:,:) = clim_var%SRAD_Wm2 *3600.e-6_r8 ![W/m2]->[MJ/m2]
  PBOT_hrly(:,:)  = clim_var%Atm_kPa

  WindMesureHeight_col(:,:) = 0.5_r8
  SolarNoonHourYM_col(:,:)    = 12._r8
  pH_rain_col(:,:)          = 7._r8
  CN4RI_col(:,:)                = 0._r8
  CNORI_col(:,:)                = 0._r8
  NH4_rain_mole_conc(:,:)   = 0._r8
  NO3_rain_mole_conc(:,:)   = 0._r8
  H2PO4_rain_mole_conc(:,:) = 0._r8

  if(salt_model)then
    trcsalt_rain_mole_conc_col(idsalt_Al,:,:)  = 0._r8
    trcsalt_rain_mole_conc_col(idsalt_Fe,:,:)  = 0._r8
    trcsalt_rain_mole_conc_col(idsalt_Ca,:,:)  = 0._r8
    trcsalt_rain_mole_conc_col(idsalt_Mg,:,:)  = 0._r8
    trcsalt_rain_mole_conc_col(idsalt_Na,:,:)  = 0._r8
    trcsalt_rain_mole_conc_col(idsalt_K,:,:)   = 0._r8
    trcsalt_rain_mole_conc_col(idsalt_SO4,:,:) = 0._r8
    trcsalt_rain_mole_conc_col(idsalt_Cl,:,:)  = 0._r8
  endif

  end subroutine SetConstClimeForcing

!------------------------------------------------------------------------------------------

  SUBROUTINE ReadClimSoilForcing(yearc,yeari,NHW,NHE,NVN,NVS,RefreshRainChem)
!
! THIS SUBROUTINE read ALL SOIL AND PLANT MANAGEMENT INPUT FILES
!
  use ReadManagementMod, only : ReadSoilManagementForYear,ReadFire
  use EcoSIMCtrlMod, only : soil_mgmt_in
  implicit none
  integer, intent(in) :: yearc   !current model year
  integer, intent(in) :: yeari   !current data year
  integer, intent(in) :: NHW,NHE,NVN,NVS
  logical, intent(out) :: RefreshRainChem !raw rainfall chemistry was reloaded

  character(len=*), parameter :: subname='ReadClimSoilForcing'
  CHARACTER(len=1) :: TTYPE
  integer :: kk,N,L,NY,NX,NZ,K
  integer :: LL,J
  integer :: LPY,IX
  CHARACTER(len=16) :: OUTW,OUTI,OUTT,OUTN,OUTF
  CHARACTER(len=4) :: CHARY
  character(len=21) :: fire_event_entry
  integer :: IDY,IFLG3,I,ICHECK

! begin_execution
  call PrintInfo('beg '//subname)
  RefreshRainChem=.false.

  ICLM=0
  call ReadClimateCorrections(yeari)
!
! PREFIX=path for files in current or higher level directory
!
! ARTIFICIAL SOIL WARMING
!
! soiltemp=file with hourly soil temperatures from baseline run
! OUT=hourly soil temperatures from baseline run (oC)
! TKSZ=temperature used to calculate additional heat flux
! for warming in watsub.f
!
! generating checkpoint files,resuming from earlier checkpt files
! DRAD,DTMPX,DTMPN,DHUM,DPREC,DIRRI,DWIND,DCNR4,DCNOR
! =annual changes in radiation,max+min temperature,humidity,
! precip,irrign,windspeed,atm CO2 concn,NH4,NO3 concn in precip
! NPX=number of cycles per hour for water,heat,solute flux calcns
! NPY=number of cycles per NPX for gas flux calcns
! JOUT,IOUT,KOUT=output frequency for hourly,daily,checkpoint data
! ICLM=changes to weather data (0=none,1=step,2=transient)

! this reads the option file that sets up frequency of model output and # of
! iterations used by the ecosim solvers

! determine whether to read checkpoing file (i.e. an actual restart run)
! read changing factor for climate variables
! read from netcdf file

  if(lverb)then
    write(*,*)'annual changes in radiation: DRAD(1:4)',DRAD(1:4)
    write(*,*)'annual changes in max temperature: DTMPX(1:4)',DTMPX(1:4)
    write(*,*)'annual changes in min temperature: DTMPN(1:4)',DTMPN(1:4)
    write(*,*)'annual changes in humidity: DHUM(1:4)',DHUM(1:4)
    write(*,*)'annual changes in precip: DPREC(1:4)',DPREC(1:4)
    write(*,*)'annual changes in irrigation: DIRRI(1:4)',DIRRI(1:4)
    write(*,*)'annual changes in wind speed: DWIND(1:4)',DWIND(1:4)
    write(*,*)'annual changes in atm NH4 conc: DCN4R(1:4)',DCN4R(1:4)
    write(*,*)'annual changes in atm NO3 conc: DCNOR(1:4)',DCNOR(1:4)
    write(*,*)'number of cycles per hour for water, heat, and '// &
      'solute flux calcns: NPX ',NPX
    write(*,*)'number of cycles per NPX for gas flux calcns: NPY',NPY
    write(*,*)'changes to weather data (0=none,1=step,'// &
      '2=transient): ICLM ',ICLM
  endif

!
! INCREMENTS IN START AND END DATES FOR SUCCESSIVE SCENARIOS
! FROM LOOPS FOR SCENES, SCENARIOS IN RUNSCRIPT SET IN MAIN.F
!
! IDATA(3),IDATA(6)=start,end year of current scene
!
  iYearCurrent=yearc

!
! OPEN CHECKPOINT FILES FOR SOIL VARIABLES
!
! IDATE=year label for checkpoint files
! DATA(1)=site file name
! W,N=water+heat,nutrient checkpoint files
!

! READ WEATHER DATA
  if(fixClime)then
    !set constatnt climate forcing
    call SetConstClimeForcing()
    RefreshRainChem=.true.
  else
    IF(is_first_year)THEN
      L=1
    ELSE
      L=2
    ENDIF
    if(yearhi/=yeari)then
      call ReadClimateForc(yearc,yeari,L,NHW,NHE,NVN,NVS)
      RefreshRainChem=.true.
    endif
  endif
  yearhi=yeari

!
!     READ LAND MANAGEMENT FILE NAMES FOR EACH GRID CELL
!
  D9980: DO NX=NHW,NHE
    D9985: DO NY=NVN,NVS
      Lirri_auto_col(NY,NX) = .false.
      ROWSpaceNH4_col(NY,NX) = 0.0_r8
      ROWSpaceNO3Band_col(NY,NX) = 0.0_r8
      ROWSpacePO4_col(NY,NX) = 0.0_r8
      D325: DO I=1,366
        iSoilDisturbType_col(I,NY,NX) = 0
        DepzCorp_col(I,NY,NX)         = 0.0_r8
      ENDDO D325
      D40: DO I=1,366
        D45: DO N=1,20
          FERT(N,I,NY,NX)=0.0_r8
        ENDDO D45
        D35: DO N=1,3
          IYTYP(N,I,NY,NX)=0
        ENDDO D35
        ExternalManureC_col(I,NY,NX)=0._r8
        FDPTH(I,NY,NX)=0.0_r8
      ENDDO D40
      D125: DO I=1,366
        RRIG(1:24,I,NY,NX)             = 0.0_r8
        PHQ(I,NY,NX)                   = 7.0_r8
        NH4_irrig_mole_conc(I,NY,NX)   = 0.0_r8
        NO3_irrig_mole_conc(I,NY,NX)   = 0.0_r8
        H2PO4_irrig_mole_conc(I,NY,NX) = 0.0_r8
        WDPTH(I,NY,NX)                 = 0.0_r8
        ROWI(I,NY,NX)                  = 0.0_r8
      ENDDO D125
    ENDDO D9985
  ENDDO D9980
  if(salt_model)then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        DO I=1,366
          trcsalt_irrig_mole_conc_col(idsalt_Al,I,NY,NX)=0.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Fe,I,NY,NX)=0.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Ca,I,NY,NX)=0.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Mg,I,NY,NX)=0.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Na,I,NY,NX)=0.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_K,I,NY,NX)=0.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_SO4,I,NY,NX)=0.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Cl,I,NY,NX)=0.0_r8  
        ENDDO
      ENDDO  
    ENDDO   
  endif 
!
! READ LAND MANAGEMENT FILE DATAC(9 LOADED IN 'MAIN'.
! THIS FILE CONTAINS NAMES OF TILLAGE, IRRIGATION
! AND FERTILIZER FILES
!
  if(use_fire)then
    if(check_fire(yearc,fire_event_entry))then
      call ReadFire(yearc,fire_event_entry,NHW,NHE,NVN,NVS)
    endif
  endif

  IF(trim(soil_mgmt_in).NE.'NO')THEN
    ! Soil events may follow the simulation calendar while climate still
    ! cycles yeari. Fire, PFT scheduling, weather and GHG routing are unchanged.
    call ReadSoilManagementForYear(yearc,yeari)
  ENDIF
  IMNG=1

  call PrintInfo('end '//subname)
  
  RETURN
  END subroutine ReadClimSoilForcing
!------------------------------------------------------------------------------------------
  subroutine ReadClimateCorrections(yeari)
  !
  !DESCRIPTION
  !read climate change factors
  use EcoSIMCtrlMod, only : clm_factor_in
  use netcdf
  use ncdio_pio
  use fileUtil, only : file_exists
  implicit none
  integer, intent(in) :: yeari
  integer :: iyear,year
  type(file_desc_t) :: clm_factor_nfid
!  type(Var_desc_t) :: vardesc
!  logical :: readvar
  INTEGER :: N

  if(clm_factor_in=='NO')return
  if(.not.file_exists(clm_factor_in))then
    call endrun("Do not find clm_factor_in file "//trim(clm_factor_in)//" in "//mod_filename, __LINE__)
  endif
  call ncd_pio_openfile(clm_factor_nfid, clm_factor_in, ncd_nowrite)

  iyear=1
  DO while(.true.)
    call ncd_getvar(clm_factor_nfid,'year',iyear,year)
    if(year==yeari)exit
    iyear=iyear+1
  ENDDO
  call ncd_getvar(clm_factor_nfid,'DRAD',iyear,DRAD(1))
  call ncd_getvar(clm_factor_nfid,'DTMPX',iyear,DTMPX(1))
  call ncd_getvar(clm_factor_nfid,'DTMPN',iyear,DTMPN(1))
  call ncd_getvar(clm_factor_nfid,'DHUM',iyear,DHUM(1))
  call ncd_getvar(clm_factor_nfid,'DPREC',iyear,DPREC(1))
  call ncd_getvar(clm_factor_nfid,'DIRRI',iyear,DIRRI(1))
  call ncd_getvar(clm_factor_nfid,'DWIND',iyear,DWIND(1))
  call ncd_getvar(clm_factor_nfid,'DCN4R',iyear,DCN4R(1))
  call ncd_getvar(clm_factor_nfid,'DCNOR',iyear,DCNOR(1))
  call ncd_getvar(clm_factor_nfid,'ICLM',iyear,ICLM)

  call ncd_pio_closefile(clm_factor_nfid)

  DO N=2,12
    DRAD(N)  = DRAD(1)
    DTMPX(N) = DTMPX(1)
    DTMPN(N) = DTMPN(1)
    DHUM(N)  = DHUM(1)
    DPREC(N) = DPREC(1)
    DIRRI(N) = DIRRI(1)
    DWIND(N) = DWIND(1)
    DCN4R(N) = DCN4R(1)
    DCNOR(N) = DCNOR(1)
  ENDDO

  end subroutine ReadClimateCorrections
!------------------------------------------------------------------------------------------
  subroutine ReadClimateForc(yearc,yeari,L,NHW,NHE,NVN,NVS) 

  implicit none
  integer, intent(in) :: yearc  ! current year of simulation
  integer, intent(in) :: yeari  ! current year of forcing data
  integer, intent(in) :: L
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8) :: Z0G,ZNOONG
  real(r8) :: PHRG             !pH of precipitation
  real(r8) :: CN4RIG,CNORIG,CN4RG,CNORG,CPORG,CALRG
  real(r8) :: CFERG,CCARG,CMGRG,CNARG,CKARG,CSORG,CCLRG
  type(atm_forc_type) :: atmf
  integer :: NY,NX

  call ReadClimNC(yearc,yeari,L,atmf)    
!
! Z0G,IFLGW=windspeed mean height,flag for raising Z0G with vegn
! ZNOONG=time of solar noon
! PHRG,CN4RIG,CNORIG,CPORG,CALRG,CFERG,CCARG,CMGRG,CNARG,CKARG,
! CSORG,CCLRG=pH,NH4,NO3,H2PO4,Al,Fe,Ca,Mg,Na,K,SO4,Cl
! concentration in precipitation
!
! CALCULATE PRECIPITATION CONCENTRATIONS IN MOLE UNITS
!
  CN4RIG = atmf%CN4RIG/natomw
  CNORIG = atmf%CNORIG/natomw
  CN4RG  = CN4RIG
  CNORG  = CNORIG
  CPORG  = atmf%CPORG/31.0_r8
  CALRG  = atmf%CALRG/27.0_r8
  CFERG  = atmf%CFERG/55.8_r8
  CCARG  = atmf%CCARG/40.0_r8
  CMGRG  = atmf%CMGRG/24.3_r8
  CNARG  = atmf%CNARG/23.0_r8
  CKARG  = atmf%CKARG/39.1_r8
  CSORG  = atmf%CSORG/32.0_r8
  CCLRG  = atmf%CCLRG/35.5_r8
  PHRG   = atmf%PHRG
  Z0G    = atmf%Z0G
  ZNOONG = atmf%ZNOONG

  D8970: DO NX=NHW,NHE
    D8975: DO NY=NVN,NVS
      WindMesureHeight_col(NY,NX) = Z0G         !windspeed meast height
      SolarNoonHourYM_col(NY,NX)    = ZNOONG
      pH_rain_col(NY,NX)          = PHRG
      CN4RI_col(NY,NX)            = CN4RIG
      CNORI_col(NY,NX)            = CNORIG
      NH4_rain_mole_conc(NY,NX)   = CN4RIG
      NO3_rain_mole_conc(NY,NX)   = CNORIG
      H2PO4_rain_mole_conc(NY,NX) = CPORG
    ENDDO D8975
  ENDDO D8970

  if(salt_model)then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        trcsalt_rain_mole_conc_col(idsalt_Al,NY,NX)=CALRG
        trcsalt_rain_mole_conc_col(idsalt_Fe,NY,NX)=CFERG
        trcsalt_rain_mole_conc_col(idsalt_Ca,NY,NX)=CCARG
        trcsalt_rain_mole_conc_col(idsalt_Mg,NY,NX)=CMGRG
        trcsalt_rain_mole_conc_col(idsalt_Na,NY,NX)=CNARG
        trcsalt_rain_mole_conc_col(idsalt_K,NY,NX)=CKARG
        trcsalt_rain_mole_conc_col(idsalt_SO4,NY,NX)=CSORG
        trcsalt_rain_mole_conc_col(idsalt_Cl,NY,NX)=CCLRG  
      ENDDO
    ENDDO    
  endif
  end subroutine ReadClimateForc    

end module readsmod
