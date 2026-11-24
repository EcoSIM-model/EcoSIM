module StartsMod
  use AtmosCompositionMod, only: UpdateAtmosGasMassConcentrations
!!
! Description:
! code to initalize soil variables

  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use abortutils,       only: padr,   print_info, check_bool
  use minimathMod,      only: isclose,  AZMAX1,   AZMIN1, real_truncate, AZERO
  use EcoSiMParDataMod, only: micpar
  use SoilHydroParaMod, only: SetColdRunSoilStates
  use SnowPhysMod,      only: InitSnowLayers
  use InitSOMBGCMod,    only: InitSOMConsts, InitSOMProfile, InitSOMVars
  use InitVegBGC,       only: InitIrradianceGeometry
  use TracerPropMod,    only: gas_solubility
  use MicrobialDiagMod, only: sumSurfOMCK
  use UnitMod,          only: units
  use DebugToolMod
  use EcosimConst
  use TracerIDMod
  use SnowDataType
  use MicrobialDataType
  use EcoSIMSolverPar
  use SOMDataType
  use EcoSIMCtrlMod
  use ChemTranspDataType
  use FertilizerDataType
  use CanopyRadDataType
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use LandSurfDataType
  use PlantTraitDataType
  use PlantDataRateType
  use SurfLitterDataType
  use CanopyDataType
  use SurfSoilDataType
  use SoilBGCDataType
  use PlantMgmtDataType
  use EcoSimSumDataType
  use RootDataType
  use EcosimBGCFluxType
  use SoilPropertyDataType
  use IrrigationDataType
  use SedimentDataType
  use GridDataType
  use MiniFuncMod

  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  !
  !
  !     BulkDensLitR=dry bulk density of woody(0),fine(1),manure(2) litter
  !     FORGC=minimum SOC for organic soil (g Mg-1)
  !      VolMaxSoilMoist4Fire,FrcAsCH4byFire=maximum SWC,CH4 emission fraction for combustion
  !     PSIHY=hygroscopic water potential (MPa)
  !     FCI,WPI=FC,WP for water retention by ice (MPa)
  !     CDPTHSI=depth to bottom of snowpack layers
  !     POROQ=Penman Water Linear Reduction tortuosity used in gas flux calculations
  !

  public :: starts
  public :: set_ecosim_solver
  public :: startsim
  public :: InitControlParms
  contains

  SUBROUTINE starts(NHW,NHE,NVN,NVS)
  !
  !     THIS SUBROUTINE INITIALIZES ALL SOIL VARIABLES

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  character(len=*), parameter :: subname='starts'
  integer :: NY,NX,L,idg
  REAL(R8) :: SurfElevLandscapeMax
  real(r8) :: SurfElevLandscapeMin
  real(r8) :: XWS
  real(r8) :: LandScape1stSoiLayDepth
  real(r8) :: YSIN(NumOfSkyAzimuthSects)
  real(r8) :: YCOS(NumOfSkyAzimuthSects)
  real(r8) :: SkyAzimuthAngle(NumOfSkyAzimuthSects)

! begin_execution

  call PrintInfo('beg '//subname)
  
  call InitControlParms  
  !
  !  IRRADIANCE INTERCEPTION GEOMETRY, plant model
  call InitIrradianceGeometry(YSIN,YCOS,SkyAzimuthAngle)
  !  CALCULATE ELEVATION OF EACH GRID CELL
  !
  call InitGridElevation(NHW,NHE,NVN,NVS,YSIN,YCOS,SkyAzimuthAngle,SurfElevLandscapeMax)
  !
  !     INITIALIZE ACCUMULATORS AND MASS BALANCE CHECKS
  !     OF EACH GRID CELL
  !
  call InitAccumulators()

  call UpdateAtmosGasMassConcentrations(NHW,NHE,NVN,NVS)

! this assumes the whole landscape is a grid
  SurfElevLandscapeMin                   = 0.0_r8
  LandScape1stSoiLayDepth = 0.0_r8

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      !
      !     MINIMUM SURFACE ELEVATION IN LANDSCAPE
      !
      !     ALT=surface elevation relative to maximum
      !     SurfElevLandscapeMin=minimum surface elevation in landscape
      !CDPTH: depth to bottom of soil layer
      SurfElevation_col(NY,NX)=SurfElevation_col(NY,NX)-SurfElevLandscapeMax
      IF(NX.EQ.NHW.AND.NY.EQ.NVN)THEN
        SurfElevLandscapeMin=SurfElevation_col(NY,NX)
      ELSE
        SurfElevLandscapeMin=MIN(SurfElevLandscapeMin,SurfElevation_col(NY,NX))
      ENDIF
      !
      LandScape1stSoiLayDepth=AMAX1(LandScape1stSoiLayDepth,CumDepz2LayBottom_vr(NU_col(NY,NX),NY,NX)) !topsoil layer depth
      DO idg=idg_beg,idg_NH3
        trcs_solcoef_col(idg,NY,NX)= AMAX1(AtmGasCgperm3_col(idg,NY,NX)*gas_solubility(idg,ATCA_col(NY,NX))*1.e-2_r8,1.e-10_r8)
      ENDDO  

!
!     MICROBIAL THERMAL ADAPTATION
!
!     OFFSET=shift in Arrhenius curve used in nitro.f (oC)
!     ATCS=mean annual soil temperature (OC)
!
      TempOffset_col(NY,NX)=f_microbe_TOffset(ATCS_col(NY,NX))
!
!     INITIALIZE WATER POTENTIAL VARIABLES FOR SOIL LAYERS
!
!     LOGPSIMX,LOGPSIMN,LOGPSIAtSat=log water potential at FC,WP,POROS
!     LOGPSIMXD,LOGPSIMND=LOGPSIMX-LOGPSIAtSat,LOGPSIMN-LOGPSIMX
!
      LOGPSIAtSat(NY,NX) = LOG(-PSIPS)
      LOGPSIFLD_col(NY,NX)   = LOG(-PSIAtFldCapacity_col(NY,NX))
      LOGPSIMN_col(NY,NX)    = LOG(-PSIAtWiltPoint_col(NY,NX))
      LOGPSIMXD_col(NY,NX)   = LOGPSIFLD_col(NY,NX)-LOGPSIAtSat(NY,NX)
      LOGPSIMND_col(NY,NX)   = LOGPSIMN_col(NY,NX)-LOGPSIFLD_col(NY,NX)
!
!     DISTRIBUTION OF OM AMONG FRACTIONS OF DIFFERING
!     BIOLOGICAL ACTIVITY
!
      
      call InitHGrid(NY,NX)
    
      call InitLayerDepths(NY,NX)

    ! ActiveLayDepZ_col=active layer depth (m)
      ActiveLayDepZ_col(NY,NX)=9999.0_r8
!
!     VHeatCapLitR,VHCPNX=minimum heat capacities for solving
!      surface litter,soil layer water and heat fluxes
      VHeatCapLitRMin_col(NY,NX) = VLHeatCapLitRMin*AREA_3D(3,NU_col(NY,NX),NY,NX)
      VHCPNX_col(NY,NX)              = VLHeatCapSoiMin*AREA_3D(3,NU_col(NY,NX),NY,NX)

!
!     SURFACE WATER STORAGE AND LOWER HEAT SINK

!      SoilHeatSrcDepth_col=depth at which soil heat sink-source calculated
!     TCNDG=assumed thermal conductivity below lower soil boundary
!     (MJ m-1 K-1 h-1)
!     TKSD=deep source/sink temperature from geothermal flux(K)

      SoilHeatSrcDepth_col(NY,NX) = AMAX1(10.0_r8,CumDepz2LayBottom_vr(NL_col(NY,NX),NY,NX)+1.0_r8)
      TKSD_col(NY,NX)             = ATKS_col(NY,NX)+2.052E-04_r8*SoilHeatSrcDepth_col(NY,NX)/TCNDG
!
    ENDDO
  ENDDO
  !
  !     INITIALIZE COMMUNITY CANOPY
  !
  CanopyHeight_col(:,:)                        = 0.0_r8
  CanopyHeightZ_col(0:NumCanopyLayers,:,:)   = 0.0_r8
  CanopyLeafAareZ_col(1:NumCanopyLayers,:,:) = 0.0_r8
  CanopyStemAareZ_col(1:NumCanopyLayers,:,:) = 0.0_r8
  tCanLeafC_clyr(1:NumCanopyLayers,:,:)        = 0.0_r8
!
  
  call InitSoilVars(NHW,NHE,NVN,NVS,SurfElevLandscapeMin,LandScape1stSoiLayDepth)

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
!     INITIALIZE SNOWPACK LAYERS
      call InitSnowLayers(NY,NX,XWS)

      WatMass_col(NY,NX) = WatMass_col(NY,NX)+XWS

    ENDDO  
  ENDDO

  call PrintInfo('end '//subname)
  RETURN
  END subroutine starts
!------------------------------------------------------------------------------------------
  subroutine InitSoilVars(NHW,NHE,NVN,NVS,SurfElevLandscapeMin,LandScape1stSoiLayDepth)
  !     N3,N2,N1=L,NY,NX of source grid cell
  !     N6,N5,N4=L,NY,NX of destination grid cell
  !      SurfElevLandscapeMin=minimum surface elevation in landscape
  !     DTBLI,TileWaterTableDepz_col=depth of natural,artificial water table
  !     WaterTable2GroundSlope_col=slope of natural water table relative to landscape surface
  !     in geography, slope =rise/run
  !     ExtWaterTableLocDepz_col,DTBLD=depth of natural,artificial water table adjusted for elevn
  !     WaterTableDepZi_col=depth to internal water table
  !
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(in) :: SurfElevLandscapeMin,LandScape1stSoiLayDepth
  integer :: NY,NX
  character(len=*), parameter :: subname='InitSoilVars'

  ! begin_execution
  call PrintInfo('beg '//subname)
  !     INITIALIZE SEDIMENT LOAD IN EROSION MODEL
  !
  IF(iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
    SED_col(:,:)=0.0_r8
  ENDIF
  !
  !     INITIALIZE GRID CELL DIMENSIONS
  ! altz: topographic altitude
  ! ALT: grid altitude
  ! DTBLI: external water table depth, before applying the altitude correction
  ! WaterTable2GroundSlope_col: slope of water table relative to surface slope
  ! ExtWaterTableLocDepz_col: external water table depth
  ! TileWaterTableDepz_col: depth of artificial water table
  ! WaterTableDepZi_col: internal water table depth
  ! DTBLD: artifical water table depth, before applying the altitude correction
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      RefSurfElevation_col(NY,NX)=SurfElevLandscapeMin
      IF(SoilBulkDensity_vr(NU_col(NY,NX),NY,NX).GT.0.0_r8)THEN
        ExtWaterTableLocDepz_col(NY,NX)  = NaturalWaterTableDepz_col(NY,NX)+ &
          (SurfElevation_col(NY,NX)-RefSurfElevation_col(NY,NX))*(1.0_r8-WaterTable2GroundSlope_col(NY,NX))
        TileWaterTableLocDepz_col(NY,NX) = AZMAX1(TileWaterTableDepz_col(NY,NX)+ &
          (SurfElevation_col(NY,NX)-RefSurfElevation_col(NY,NX))*(1.0_r8-WaterTable2GroundSlope_col(NY,NX)))
      ELSE
        ExtWaterTableLocDepz_col(NY,NX)  = 0.0_r8
        TileWaterTableLocDepz_col(NY,NX) = 0.0_r8
      ENDIF
      WaterTableDepZi_col(NY,NX)=ExtWaterTableLocDepz_col(NY,NX)
    ENDDO
  ENDDO

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      !
      !     ALLOCATE LITTER,SOC TO WOODY,NON-WOODY,MANURE,POC AND HUMUS
      !
      call InitSoilProfile(NY,NX,LandScape1stSoiLayDepth)
      !
      !     SURFACE LITTER HEAT CAPACITY
      !
      SoilMicPMassLayerMn(NY,NX)=AZMAX1(SAND_vr(NU_col(NY,NX),NY,NX)+SILT_vr(NU_col(NY,NX),NY,NX)+CLAY_vr(NU_col(NY,NX),NY,NX))
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine InitSoilVars
!------------------------------------------------------------------------------------------
  subroutine InitSoilProfile(NY,NX,LandScape1stSoiLayDepth)

  implicit none
  integer, intent(in) :: NY,NX
  REAL(R8),INTENT(IN) :: LandScape1stSoiLayDepth
  character(len=*), parameter :: subname='InitSoilProfile'
  integer  :: L,M,K,N,KK,NN
  real(r8) :: CORGCM          !organic-matter mass fraction [Mg OM/Mg soil]
  real(r8) :: HCX,TORGC
  real(r8) :: CORGL,TORGLL,FCX
  REAL(R8) :: PTDS  !particle density [g cm-3]=[Mg m-3]
  real(r8) :: TORGM
  real(r8) :: VSAND
  real(r8) :: TORGL(JZ)
  real(r8) :: VMINL
  real(r8) :: VORGC,XS
  real(r8) :: VSolidSoil !volume of solid soil in one m3 volume, including mineral and organic particles
! begin_execution
! RSC,RSC,RSP=C,N,P in fine(1),woody(0),manure(2) litter (g m-2)
! CORGC,CORGR,CORGN,CORGP=SOC,POC,SON,SOP (g Mg-1)
! BKVL=soil mass (Mg)
!

!
!     INITIALIZE SOM FROM ORGANIC INPUTS IN SOIL FILE FROM 'READS'
!
!     CORGC,CORGR,CORGN,CORGP=SOC,POC,SON,SOP (g Mg-1)
!
  call PrintInfo('beg '//subname)
  TORGC=0.0_r8;TORGL=0._r8
  D1190: DO L=NU_col(NY,NX),NL_col(NY,NX)
    !     CORGCZ=CSoilOrgM_vr(ielmc,L,NY,NX)
    !     CORGRZ=CORGR(L,NY,NX)
    !     CORGNZ=CSoilOrgM_vr(ielmn,L,NY,NX)
    !     CORGPZ=CSoilOrgM_vr(ielmp,L,NY,NX)
    !
    !     ALLOCATE SOC TO POC(3) AND HUMUS(4)
    !
    !     CORGCX(3)=CORGRZ
    !     CORGCX(4)=AZMAX1(CORGCZ-CORGCX(3))
    !     CORGNX(3)=AMIN1(CNRH(3)*CORGCX(3),CORGNZ)
    !     CORGNX(4)=AZMAX1(CORGNZ-CORGNX(3))
    !     CORGPX(3)=AMIN1(CPRH(3)*CORGCX(3),CORGPZ)
    !     CORGPX(4)=AZMAX1(CORGPZ-CORGPX(3))
    CORGL    = AZMAX1(CSoilOrgM_vr(ielmc,L,NY,NX)-COMLitrC_vr(L,NY,NX))
    TORGL(L) = TORGC+CORGL*VLSoilMicPMass_vr(L,NY,NX)/AREA_3D(3,L,NY,NX)*0.5_r8
    TORGC    = TORGC+CORGL*VLSoilMicPMass_vr(L,NY,NX)/AREA_3D(3,L,NY,NX)
  ENDDO D1190
!
!     PARAMETERS TO ALLOCATE HUMUS TO LESS OR MORE RECALCITRANT FRACTIONS
!
!     TORGL=accumulated humus down to soil layer (g m-2)
!     TORGM=TORGL used to calculate allocation (g m-2)
!     HCX=shape parameter for depth effect on allocation
!
  TORGM=AMAX1(2.0E+03_r8,AMIN1(5.0E+03_r8,0.25_r8*TORGL(MaxNumRootLays_col(NY,NX))))
  IF(TORGM.GT.ZERO)THEN
    HCX=LOG(0.5_r8)/TORGM
  ELSE
    HCX=0.0_r8
  ENDIF
  WatMass_col(NY,NX) = 0._r8

  XS  = 0._r8
  D1200: DO L=0,NL_col(NY,NX)
    !
    if(L==0)then
      TORGLL=0.0_r8
    else
      TORGLL=TORGL(L)
    endif

    call InitSOMProfile(L,NY,NX,HCX,TORGLL,LandScape1stSoiLayDepth,FCX)
    !
    !     LAYER WATER, ICE, AIR CONTENTS
    !
    !     PSISE,PSISA=water potential at saturation,air entry (MPa)
    !     PTDS=particle density (Mg m-3)
    !     POROS=total porosity
    !     VOLA,VOLAH=micropore,macropore volume
    !     THW,THI=initial soil water,ice content
    !     VOLW,VOLWH=micropore,macropore water volume(m3)
    !     VOLI,VOLIH=micropore,macropore ice volume(m3)
    !     VOLP=total air volume (m3)
    !     PSISoilAirEntry_vr: soil matric potential at air entry point [MPa]
    PSI_at_Sautration_vr(L,NY,NX)         = PSIPS
    PSISoilAirEntry_vr(L,NY,NX)           = -1.5E-03_r8
    RGasTranspFlxPrev_vr(idg_O2,L,NY,NX)  = 0.0_r8
    RGasTranspFlxPrev_vr(idg_CO2,L,NY,NX) = 0.0_r8
    RO2AquaSourcePrev_vr(L,NY,NX)         = 0.0_r8
    RGasTranspFlxPrev_vr(idg_CH4,L,NY,NX) = 0.0_r8
    RCH4PhysexchPrev_vr(L,NY,NX)          = 0.0_r8

    IF(L.GT.0)THEN
      IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
        ! Use the organic-matter mass fraction and density mixture from Hour1Mod.
        CORGCM            = AZMAX1(AMIN1(1.0_r8,gC2MgOM*CSoilOrgM_vr(ielmc,L,NY,NX)))
        PTDS              = DensitySolidOM*CORGCM+DensitySolidMineral*(1.0_r8-CORGCM)
        VSolidSoil        = SoilBulkDensity_vr(L,NY,NX)/PTDS      ! volume of soil as solid [Mg m-3]/[Mg m-3]
        POROS_vr(L,NY,NX) = AZMAX1(1.0_r8-VSolidSoil)             ! volume of pores in one unit soil
      ELSE
        !for ponding water
        VSolidSoil        = 0._r8
        PTDS              = 0.0_r8
        POROS_vr(L,NY,NX) = 1.0_r8
      ENDIF
      POROSI_vr(L,NY,NX)    = POROS_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
      VLMicP_vr(L,NY,NX)    = POROS_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
      VLMacroPores_vr(L,NY,NX)    = SoilFracAsMacP_vr(L,NY,NX)*VGeomLayert0_vr(L,NY,NX)
      !
      !     LAYER HEAT CONTENTS
      !
      !     SAND,SILT,CLAY=sand,silt,clay mass (Mg)
      !     VORGC,VMINL,VSAND=volumetric fractions of SOC,non-sand,sand
      !     VHCM,VHCP=volumetric dry,wet soil heat capacity (MJ m-3 K-1)
      !     TKS_vr,TCS=soil temperature (oC,K)
      !     THETW,THETI,THETP=micropore water,ice,air concentration (m3 m-3)
      !
      SAND_vr(L,NY,NX)=CSAND_vr(L,NY,NX)*VLSoilMicPMass_vr(L,NY,NX) ![g sand/g soil ]*[Mg soil d-2]=[Mg sand d-2]
      SILT_vr(L,NY,NX)=CSILT_vr(L,NY,NX)*VLSoilMicPMass_vr(L,NY,NX) ![Mg silt d-2]
      CLAY_vr(L,NY,NX)=CCLAY_vr(L,NY,NX)*VLSoilMicPMass_vr(L,NY,NX) ![Mg clay d-2]
      
      IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
        ! PTDS=particle density [Mg m-3]
        ! soil volumetric heat capacity, [SoilBulkDensity_vr]=[Mg soil m-3]
        VORGC = CORGCM*VSolidSoil                                 ![m3 OM m-3 soil], same convention as Hour1Mod
        VMINL = (CSILT_vr(L,NY,NX)+CCLAY_vr(L,NY,NX))*VSolidSoil   ![m3 mineral m-3 soil]
        VSAND = CSAND_vr(L,NY,NX)*VSolidSoil                       ![m3 sand m-3 soil]
        VHeatCapMineralSoil_vr(L,NY,NX)=((cpmins*VMINL+cpsand*VSAND) &
          *FracSoiAsMicP_vr(L,NY,NX)+cpsand*ROCK_vr(L,NY,NX))*VGeomLayer_vr(L,NY,NX)     ![MJ m-3 K-1]*[m3]=[MJ K-1]          
        VHeatCapSolidSoil_vr(L,NY,NX)=VHeatCapMineralSoil_vr(L,NY,NX)+cpo*VORGC*FracSoiAsMicP_vr(L,NY,NX)*VGeomLayer_vr(L,NY,NX)     ![MJ m-3 K-1]*[m3]=[MJ K-1]          
      ELSE
        VHeatCapSolidSoil_vr(L,NY,NX)=0.0_r8
        VHeatCapMineralSoil_vr(L,NY,NX)=0._r8
      ENDIF
      !
      !     INITIAL SOIL WATER AND ICE CONTENTS
      ! the use of POROS in defining relative saturation is not very appropriate.
      ! Resolve supplied or estimated FC/WP and initialize all pore volumes
      ! before STARTE uses them to construct chemical inventories.
      call SetColdRunSoilStates(L,NY,NX)
      XS=XS+VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)+(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))*DENSICE
      HBAmin_vr(L,NY,NX) = 7._r8 * CSAND_vr(L,NY,NX)+3._r8*CCLAY_vr(L,NY,NX)+4.5_r8*CSILT_vr(L,NY,NX)
      !
      HBAconst_vr(L,NY,NX)=AMAX1(HBAmin_vr(L,NY,NX)*exp(-0.035_r8*SolidOMPercent_vr(L,NY,NX)),0.5_r8)
      
    ELSEIF(L.eq.0)THEN 
    
      IF(VLitR_col(NY,NX).GT.ZEROS(NY,NX))THEN
        VHeatCapacity_vr(0,NY,NX)=cpo*gOC_to_m3_OM(SoilOrgM_vr(ielmc,0,NY,NX))+cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)
      else
        VHeatCapacity_vr(0,NY,NX)=0._r8
      endif
      POROSI_vr(0,NY,NX)             = 1._r8  !this is added for numerical fixing
      VHeatCapSolidSoil_vr(0,NY,NX) = 0.0_r8
!      VLMicPt0_col(0,NY,NX)          = 0.0_r8
      XS=XS+VLWatMicP_vr(0,NY,NX)+VLiceMicP_vr(0,NY,NX)*DENSICE
    ENDIF

    TKS_vr(L,NY,NX) = ATKS_col(NY,NX)
    TCS_vr(L,NY,NX)    = ATCS_col(NY,NX)
    !
    !     INITIALIZE SOM VARIABLES
    call InitSOMVars(L,NY,NX,FCX)
    
  ENDDO D1200

  WatMass_col(NY,NX) = WatMass_col(NY,NX)+XS

  call sumSurfOMCK(NY,NX,SurfLitterC_col(:,NY,NX),RC0ff_col(NY,NX))

  !
  !  INITIALIZE FERTILIZER ARRAYS
  call initFertArrays(NY,NX)
  
  call PrintInfo('end '//subname)
  end subroutine InitSoilProfile
!------------------------------------------------------------------------------------------

  subroutine initFertArrays(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L2

! begin_execution
  L2=NL_col(NY,NX)
  FertN_mole_soil_vr(ifertn_beg:ifertn_end,0:L2,NY,NX)   = 0._r8
  FertN_mole_Band_vr(ifertnb_beg:ifertnb_end,1:L2,NY,NX) = 0._r8
  trcs_VLN_vr(ids_NH4,0:L2,NY,NX)                   = 1.0_r8
  trcs_VLN_vr(idg_NH3,0:L2,NY,NX)                   = trcs_VLN_vr(ids_NH4,0:L2,NY,NX)
  trcs_VLN_vr(ids_NO3,0:L2,NY,NX)                   = 1.0_r8
  trcs_VLN_vr(ids_NO2,0:L2,NY,NX)                   = trcs_VLN_vr(ids_NO3,0:L2,NY,NX)
  trcs_VLN_vr(ids_H1PO4,0:L2,NY,NX)                 = 1.0_r8
  trcs_VLN_vr(ids_H2PO4,0:L2,NY,NX)                 = trcs_VLN_vr(ids_H1PO4,0:L2,NY,NX)
  trcs_VLN_vr(ids_NH4B,0:L2,NY,NX)                  = 0.0_r8
  trcs_VLN_vr(idg_NH3B,0:L2,NY,NX)                  = trcs_VLN_vr(ids_NH4B,0:L2,NY,NX)
  trcs_VLN_vr(ids_NO3B,0:L2,NY,NX)                  = 0.0_r8
  trcs_VLN_vr(ids_NO2B,0:L2,NY,NX)                  = trcs_VLN_vr(ids_NO3B,0:L2,NY,NX)
  trcs_VLN_vr(ids_H1PO4B,0:L2,NY,NX)                = 0.0_r8
  trcs_VLN_vr(ids_H2PO4B,0:L2,NY,NX)                = trcs_VLN_vr(ids_H1PO4B,0:L2,NY,NX)

  REcoO2DmndResp_vr(0:L2,NY,NX)    = 0.0_r8
  REcoNH4DmndSoil_vr(0:L2,NY,NX)   = 0.0_r8
  REcoNO3DmndSoil_vr(0:L2,NY,NX)   = 0.0_r8
  RNO2EcoUptkSoil_vr(0:L2,NY,NX)   = 0.0_r8
  RN2OEcoUptkSoil_vr(0:L2,NY,NX)   = 0.0_r8
  REcoH2PO4DmndSoil_vr(0:L2,NY,NX) = 0.0_r8
  REcoH1PO4DmndSoil_vr(0:L2,NY,NX) = 0.0_r8
  RNO2DmndSoilChemo_vr(0:L2,NY,NX) = 0.0_r8
  REcoNH4DmndBand_vr(0:L2,NY,NX)   = 0.0_r8
  REcoNO3DmndBand_vr(0:L2,NY,NX)   = 0.0_r8
  RNO2EcoUptkBand_vr(0:L2,NY,NX)   = 0.0_r8
  REcoH2PO4DmndBand_vr(0:L2,NY,NX) = 0.0_r8
  REcoH1PO4DmndBand_vr(0:L2,NY,NX) = 0.0_r8
  RNO2DmndBandChemo_vr(0:L2,NY,NX) = 0.0_r8
  ZNHUI_vr(0:L2,NY,NX)                = 0.0_r8
  ZNHU0_vr(0:L2,NY,NX)                = 0.0_r8
  ZNFNI_vr(0:L2,NY,NX)                = 0.0_r8
  ZNFN0_vr(0:L2,NY,NX)                = 0.0_r8

  BandWidthNH4_vr(1:L2,NY,NX)     = 0.0_r8
  BandThicknessNH4_vr(1:L2,NY,NX) = 0.0_r8
  BandWidthNO3_vr(1:L2,NY,NX)     = 0.0_r8
  BandThicknessNO3_vr(1:L2,NY,NX) = 0.0_r8
  BandWidthPO4_vr(1:L2,NY,NX)     = 0.0_r8
  BandThicknessPO4_vr(1:L2,NY,NX) = 0.0_r8

  end subroutine initFertArrays

!------------------------------------------------------------------------------------------
  subroutine InitGridElevation(NHW,NHE,NVN,NVS,YSIN,YCOS,SkyAzimuthAngle,SurfElevLandscapeMax)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(in) :: YSIN(NumOfSkyAzimuthSects)
  real(r8),intent(in) :: YCOS(NumOfSkyAzimuthSects)
  real(r8),intent(in) :: SkyAzimuthAngle(NumOfSkyAzimuthSects)
  REAL(R8),INTENT(OUT):: SurfElevLandscapeMax  !peak elevation of the landscape [m]
  integer :: NY,NX,N,NN
  REAL(R8) :: DGAZI,GridArea
  real(r8) :: SineGrndSurfAzimuth_col(JY,JX),CosineGrndSurfAzimuth_col(JY,JX)  !diagnostic

! begin_execution
! GroundSurfaceAzimuth_col=ground surface azimuth, aspect in radians
! SineGrndSlope_col,CosineGrndSlope_col=sine,cosine of ground surface
! OMEGA2Ground=incident sky angle at ground surface
! SLOPE=tangent of ground surface slope in (0)aspect, (1)EW,(2)NS directions
! ALT=ground surface elevation
! SurfElevLandscapeMax=maximum surface elevation in landscape
! XGridRunoffFlag_2DH=runoff boundary flags:0=not possible,1=possible
! ASP_col=geometric aspect angle in degree
  SurfElevLandscapeMax=0.0_r8
  write(*,1112)'NY','NX','east','west','south','north','altitude','Dist(m):E-W','Dist(m):N-S',&
   'Geom-ASP(o)','SLOPE_col(o)','slope0','slope-east','slope-north','SineGrndSlope_col','CosineGrndSlope_col','SineGrndSurfAzimuth_col'

1112    FORMAT(2A4,4A6,25A12)
  D9985: DO NX=NHW,NHE
    D9980: DO NY=NVN,NVS
      GridArea=DH_col(NY,NX)*DV_col(NY,NX)
      ZEROS(NY,NX)  = ZERO*GridArea
      ZEROS2(NY,NX) = ZERO2*GridArea
!     compute slopes
      GroundSurfaceAzimuth_col(NY,NX)       = ASP_col(NY,NX)*RadianPerDegree   !radian
      SineGrndSurfAzimuth_col(NY,NX)        = ABS(SIN(GroundSurfaceAzimuth_col(NY,NX)))
      CosineGrndSurfAzimuth_col(NY,NX)      = ABS(COS(GroundSurfaceAzimuth_col(NY,NX)))
      SLOPE_col(0,NY,NX)                    = AMAX1(1.745E-04_r8,TAN(SL_col(NY,NX)*RadianPerDegree))  !minimum slope is 1.745E-4~0.01o
      SLOPE_col(iWestEastDirection,NY,NX)   = -AZERO(SLOPE_col(0,NY,NX)*COS(GroundSurfaceAzimuth_col(NY,NX)))   !elevation increases from west to east, if west facing
      SLOPE_col(iNorthSouthDirection,NY,NX) = AZERO(SLOPE_col(0,NY,NX)*SIN(GroundSurfaceAzimuth_col(NY,NX)))    !elevation increases from north to south, if north facing

      !aspect angle 
      IF(ASP_col(NY,NX).GE.0.0_r8 .AND. ASP_col(NY,NX).LT.90.0_r8)THEN
      ! face the northeast
        XGridRunoffFlag_2DH(1,1,NY,NX) = .true.    !west->east
        XGridRunoffFlag_2DH(2,1,NY,NX) = .false.
        XGridRunoffFlag_2DH(1,2,NY,NX) = .false.
        XGridRunoffFlag_2DH(2,2,NY,NX) = .true.    !south->north
      ELSEIF(ASP_col(NY,NX).GE.90.0_r8.AND.ASP_col(NY,NX).LT.180.0_r8)THEN
      ! face the northwest
        XGridRunoffFlag_2DH(1,1,NY,NX) = .false.
        XGridRunoffFlag_2DH(2,1,NY,NX) = .true.   !east -> west
        XGridRunoffFlag_2DH(1,2,NY,NX) = .false.
        XGridRunoffFlag_2DH(2,2,NY,NX) = .true.   !south -> north
      ELSEIF(ASP_col(NY,NX).GE.180.0_r8.AND.ASP_col(NY,NX).LT.270.0_r8)THEN
      !face the southwest
        XGridRunoffFlag_2DH(1,1,NY,NX) = .false.
        XGridRunoffFlag_2DH(2,1,NY,NX) = .true.  !east -> west
        XGridRunoffFlag_2DH(1,2,NY,NX) = .true.  !north -> south
        XGridRunoffFlag_2DH(2,2,NY,NX) = .false.
      ELSEIF(ASP_col(NY,NX).GE.270.0_r8.AND.ASP_col(NY,NX).LE.360.0_r8)THEN
      ! face the southeast
        XGridRunoffFlag_2DH(1,1,NY,NX) = .true.  !west->east
        XGridRunoffFlag_2DH(2,1,NY,NX) = .false.
        XGridRunoffFlag_2DH(1,2,NY,NX) = .true.  !north->south
        XGridRunoffFlag_2DH(2,2,NY,NX) = .false.
      ENDIF
      !In the vertical direction, for general treatment of flow against slope. 
      SLOPE_col(3,NY,NX)=-1.0_r8

      IF(.not.isclose(SLOPE_col(iWestEastDirection,NY,NX),0.0_r8) .OR. (.not.isclose(SLOPE_col(iNorthSouthDirection,NY,NX),0.0_r8)))THEN
        FSLOPE_2DH(iWestEastDirection,NY,NX)   = ABS(SLOPE_col(iWestEastDirection,NY,NX))/(ABS(SLOPE_col(iWestEastDirection,NY,NX))+ABS(SLOPE_col(iNorthSouthDirection,NY,NX)))  !
        FSLOPE_2DH(iNorthSouthDirection,NY,NX) = 1._r8-FSLOPE_2DH(iWestEastDirection,NY,NX)
      ELSE
        FSLOPE_2DH(iWestEastDirection,NY,NX)   = 0.5_r8
        FSLOPE_2DH(iNorthSouthDirection,NY,NX) = 0.5_r8
      ENDIF

!    Compute radiation geometry with sine, independently of the tangent used for SLOPE_col.
      SineGrndSlope_col(NY,NX)   = AMAX1(1.745E-04_r8,SIN(SL_col(NY,NX)*RadianPerDegree))    !preserve the radiation minimum slope
      CosineGrndSlope_col(NY,NX) = SQRT(1.0_r8-SineGrndSlope_col(NY,NX)**2)
      D240: DO N=1,NumOfSkyAzimuthSects
        DGAZI                 = COS(GroundSurfaceAzimuth_col(NY,NX)-SkyAzimuthAngle(N))
        OMEGA2Ground(N,NY,NX) = AZMAX1(AMIN1(1.0_r8,CosineGrndSlope_col(NY,NX)*YSIN(N)+ &
          SineGrndSlope_col(NY,NX)*YCOS(N)*DGAZI))
      ENDDO D240

!     compute ground surface elevation
!     DH, length in e-w direction
      IF(NX.EQ.NHW)THEN
        IF(NY.EQ.NVN)THEN
          !(west, north) corner as the anchor point
          SurfElevation_col(NY,NX)=0.5_r8*DH_col(NY,NX)*SLOPE_col(iWestEastDirection,NY,NX)+0.5_r8*DV_col(NY,NX)*SLOPE_col(iNorthSouthDirection,NY,NX)
        ELSE
          !west boundary: integrate the north-south slope between cell centers
          SurfElevation_col(NY,NX)=SurfElevation_col(NY-1,NX) &
            +0.5_r8*DV_col(NY-1,NX)*(SLOPE_col(iNorthSouthDirection,NY-1,NX)) &
            +0.5_r8*DV_col(NY,NX)*SLOPE_col(iNorthSouthDirection,NY,NX)
        ENDIF
      ELSE
        IF(NY.EQ.NVN)THEN
          !north boundary: integrate the east-west slope between cell centers
          SurfElevation_col(NY,NX)=SurfElevation_col(NY,NX-1) &
            +0.5_r8*DH_col(NY,NX-1)*SLOPE_col(iWestEastDirection,NY,NX-1) &
            +0.5_r8*DH_col(NY,NX)*SLOPE_col(iWestEastDirection,NY,NX)
        ELSE
          SurfElevation_col(NY,NX)=(SurfElevation_col(NY,NX-1) &
            +0.5_r8*DH_col(NY,NX-1)*SLOPE_col(iWestEastDirection,NY,NX-1) &
            +0.5_r8*DH_col(NY,NX)*SLOPE_col(iWestEastDirection,NY,NX) &
            +SurfElevation_col(NY-1,NX) &
            +0.5_r8*DV_col(NY-1,NX)*SLOPE_col(iNorthSouthDirection,NY-1,NX) &
            +0.5_r8*DV_col(NY,NX)*SLOPE_col(iNorthSouthDirection,NY,NX))/2.0_r8
        ENDIF
      ENDIF

      IF(NX.EQ.NHW .AND. NY.EQ.NVN)THEN
        !(west,north) corner, anchor point
        SurfElevLandscapeMax=SurfElevation_col(NY,NX)
      ELSE
        SurfElevLandscapeMax=MAX(SurfElevLandscapeMax,SurfElevation_col(NY,NX))
      ENDIF
      WRITE(iulog,1111)NY,NX,((XGridRunoffFlag_2DH(NN,N,NY,NX),NN=1,2),N=1,2) &
        ,SurfElevation_col(NY,NX),DH_col(NY,NX),DV_col(NY,NX),ASP_col(NY,NX),SL_col(NY,NX) &
        ,SLOPE_col(0,NY,NX),SLOPE_col(iWestEastDirection,NY,NX),SLOPE_col(iNorthSouthDirection,NY,NX) &
        ,SineGrndSlope_col(NY,NX),CosineGrndSurfAzimuth_col(NY,NX),SineGrndSurfAzimuth_col(NY,NX)
1111  FORMAT(2I4,4L6,20E12.4)
    ENDDO D9980
  ENDDO D9985
  end subroutine InitGridElevation
!------------------------------------------------------------------------------------------
  subroutine InitControlParms
  implicit none
  !     begin_execution
  !
  !     NPH=no. of cycles per hour for water, heat and solute flux calculns
  !     NPT=number of cycles per water iteration for gas flux calculations
  !     NPG=number of cycles per hour for gas flux calculations
  !     NPR,NPS=number of cycles NPH-1 for litter, snowpack flux calculns
  !     AirFillPore_Min=minimum air-filled porosity for gas flux calculations
  !     THETPI,DENSICE=ice porosity,density
  !     BulkDensLitR=surface litter bulk density, Mg m-3

  BulkDensLitR=(/0.0333_r8,0.0167_r8,0.0167_r8/)

  call InitSOMConsts

  ZERO  = 1.0E-16_r8
  ZERO2 = 1.0E-08_r8
  end subroutine InitControlParms
!------------------------------------------------------------------------------------------
  subroutine InitAccumulators()
  implicit none

!     begin_execution
  TDTPX(:,:,:) = 0.0_r8
  TDTPN(:,:,:) = 0.0_r8
  TDRAD(:,:,:) = 1.0_r8
  TDWND(:,:,:) = 1.0_r8
  TDHUM(:,:,:) = 1.0_r8
  TDPRC(:,:,:) = 1.0_r8
  TDIRI(:,:,:) = 1.0_r8
  TDCN4(:,:,:) = 1.0_r8
  TDCNO(:,:,:) = 1.0_r8

  iUreaHydInhibitorType_col(:,:) = 0
  iFertNH4Band_col(:,:)          = ifert_off
  iFertNO3Band_col(:,:)          = ifert_off
  iFertPO4Band_col(:,:)          = ifert_off
  iResetSoilProf_col(:,:)        = itrue
  NumActivePlants_col(:,:)       = 0
  ATCA_col(:,:)                  = ATCAI_col(:,:)
  ATCS_col(:,:)                  = ATCAI_col(:,:)
  TairKClimMean_col(:,:)         = units%Celcius2Kelvin(ATCA_col)
  ATKS_col(:,:)                  = units%Celcius2Kelvin(ATCS_col)
  QRain_CumYr_col(:,:)           = 0.0_r8

  CO2byFire_CumYr_col(:,:)       = 0.0_r8
  CH4byFire_CumYr_col(:,:)       = 0.0_r8
  O2byFire_CumYr_col(:,:)        = 0.0_r8
  N2ObyFire_CumYr_col(:,:)       = 0.0_r8
  NH3byFire_CumYr_col(:,:)       = 0.0_r8
  PO4byFire_CumYr_col(:,:)       = 0.0_r8
  AmendC_CumYr_flx_col(:,:)       = 0.0_r8
  FertN_Flx_CumYr_col(:,:)        = 0.0_r8
  FerP_Flx_CumYr_col(:,:)         = 0.0_r8
  H2OLoss_CumYr_col(:,:)         = 0.0_r8
  QEvap_CumYr_col(:,:)           = 0.0_r8
  Qrunoff_CumYr_col(:,:)         = 0.0_r8
  SedmErossLoss_CumYr_col(:,:)   = 0.0_r8
  RootResp_CumYr_col(:,:)        = 0.0_r8
  HydroSufDOCFlx_col(:,:)        = 0.0_r8
  HydroSubsDOCFlx_col(:,:)       = 0.0_r8
  HydroSufDONFlx_col(:,:)  = 0.0_r8
  HydroSubsDONFlx_col(:,:)       = 0.0_r8
  HydroSufDOPFlx_col(:,:)  = 0.0_r8
  HydroSubsDOPFlx_col(:,:)       = 0.0_r8
  HydroSufDICFlx_col(:,:)        = 0.0_r8
  HydroSubsDICFlx_col(:,:)       = 0.0_r8
  HydroSufDINFlx_CumYr_col(:,:)  = 0.0_r8
  HydroSubsDINFlx_col(:,:)       = 0.0_r8
  HydroSufDIPFlx_CumYr_col(:,:)  = 0.0_r8
  HydroSubsDIPFlx_col(:,:)       = 0.0_r8
  HydroIonFlx_CumYr_col(:,:)     = 0.0_r8
  LiterfalOrgM_col(ielmc,:,:)    = 0.0_r8
  LiterfalOrgM_col(ielmn,:,:)    = 0.0_r8
  LiterfalOrgM_col(ielmp,:,:)    = 0.0_r8
  QDrain_col(:,:)                = 0.0_r8
  ZDRAIN_col(:,:)                    = 0.0_r8
  PDRAIN_col(:,:)                    = 0.0_r8
  BandDepthNH4_col(:,:)          = 0.0_r8
  BandDepthNO3_col(:,:)          = 0.0_r8
  BandDepthPO4_col(:,:)          = 0.0_r8
  trcs_solml_vr(idg_O2,0,:,:)     = 0.0_r8
  FracSWRad2Grnd_col(:,:)        = 1.0_r8
  LWRadBySurf_col(:,:)           = 0.0_r8
  LWRadCanG_col(:,:)                 = 0.0_r8
  Eco_NetRad_col(:,:)            = 0.0_r8
  Eco_Heat_Latent_col(:,:)       = 0.0_r8
  Eco_Heat_Sens_col(:,:)         = 0.0_r8
  Eco_Heat_GrndSurf_col(:,:)     = 0.0_r8
  Air_Heat_Latent_store_col(:,:)    = 0.0_r8
  Air_Heat_Sens_store_col(:,:)      = 0.0_r8
  dSensHeatdTKQ_col(:,:)             = 0._r8
  TKQHeatRef_col(:,:)                = 0._r8
  dLatHeatdVPQ_col(:,:)              = 0._r8
  VPQLatHeatRef_col(:,:)             = -1._r8
  TLEX_col(:,:)                      = 0.0_r8
  TSHX_col(:,:)                      = 0.0_r8
  Eco_NEE_col(:,:)               = 0.0_r8
  WatHeldOnCanopy_col(:,:)          = 0.0_r8
  CanopyLeafArea_col(:,:)        = 0.0_r8
  StemArea_col(:,:)              = 0.0_r8
  RainIntceptByCanopy_col(:,:)   = 0.0_r8
  PlantPopu_col(:,:)             = 0.0_r8
  DayLensCurr_col(:,:)           = 12.0_r8
  SurfAlbedo_col(:,:)            = SoilAlbedo_col(:,:)
  EcoHavstElmnt_CumYr_col(:,:,:) = 0.0_r8
  EnergyImpact4Erosion_col(:,:)      = 0.0_r8
  end subroutine InitAccumulators

!------------------------------------------------------------------------------------------
  subroutine InitHGrid(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L

  call PrintInfo('beg InitHGrid')
  DO  L=0,NL_col(NY,NX)
    DLYRI_3D(1,L,NY,NX) = DH_col(NY,NX)        !east-west direction
    DLYRI_3D(2,L,NY,NX) = DV_col(NY,NX)        !north-south direction
    DLYR_3D(1,L,NY,NX)  = DLYRI_3D(1,L,NY,NX)
    DLYR_3D(2,L,NY,NX)  = DLYRI_3D(2,L,NY,NX)
    AREA_3D(3,L,NY,NX)     = DLYR_3D(1,L,NY,NX)*DLYR_3D(2,L,NY,NX)  !grid horizontal area
  ENDDO
  call PrintInfo('end InitHGrid')
  end subroutine InitHGrid
!------------------------------------------------------------------------------------------
  subroutine InitLayerDepths(NY,NX)

  implicit none
  integer, intent(in) :: NY, NX
  integer :: L,K,j
  real(r8) :: VLitR0
  character(len=*), parameter :: subname='InitLayerDepths'

  associate(                                   &
    NumOfLitrCmplxs => micpar%NumOfLitrCmplxs, &
    k_woody_comp    => micpar%k_woody_comp,    &
    k_fine_comp     => micpar%k_fine_comp,     &
    k_manure        => micpar%k_manure         &
  )
!     begin_execution
  call PrintInfo('beg '//subname)

  D111: DO   L=0,NL_col(NY,NX)
    !
    ! LAYER DEPTHS AND THEIR PHYSICAL PROPERTIES

    ! surface litter:L=0,soil layer:L>0
    ! DLYR,AREA=layer thickness,x-sectional area:1=EW,2=NS,3=vertical
    ! ORGC=organic C content (gC d-2)
    ! VOLT,VOLX=volume including,excluding macropores+rock
    ! BKVL=mass
    ! CDPTH,DPTH=depth to bottom,midpoint
    !
    IF(L.EQ.0)THEN
      ! surface litter residue layer
      TAREA                        = TAREA+AREA_3D(3,L,NY,NX)
      CumSoilThickness_vr(L,NY,NX) = 0.0_r8
      SoilSurfDepZ_col(NY,NX)      = 0._r8
      SoilOrgM_vr(ielmc,L,NY,NX)   = SUM(RSC_vr(1:NumOfLitrCmplxs,L,NY,NX))*AREA_3D(3,L,NY,NX)
      ORGCX_vr(L,NY,NX)            = SoilOrgM_vr(ielmc,L,NY,NX)
      VLitR0                       = 0._r8
      DO K=1,NumOfLitrCmplxs
        VLitR0=VLitR0+RSC_vr(K,L,NY,NX)/BulkDensLitR(K)
      ENDDO
      !volume of litter layer
      VLitR_col(NY,NX)           = VLitR0*ppmc*AREA_3D(3,L,NY,NX)
      
      VGeomLayer_vr(L,NY,NX)     = VLitR_col(NY,NX)
      VLSoilPoreMicP_vr(L,NY,NX) = VGeomLayer_vr(L,NY,NX)
      VLSoilMicP_vr(L,NY,NX)     = VLSoilPoreMicP_vr(L,NY,NX)
      VGeomLayert0_vr(L,NY,NX)   = VGeomLayer_vr(L,NY,NX)
      VLSoilMicPMass_vr(L,NY,NX) = gC2MgOM*SoilOrgM_vr(ielmc,L,NY,NX)  !mass of soil layer, Mg/d2
      !thickness of litter layer 

      DLYRI_3D(3,L,NY,NX) = VLSoilPoreMicP_vr(L,NY,NX)/AREA_3D(3,L,NY,NX)
      DLYR_3D(3,L,NY,NX)  = DLYRI_3D(3,L,NY,NX)
    ELSE

      ! if it is a standing water, no macropore fraction
      ! DPTH=depth of layer middle
      ! CumSoilThickness_vr=soil thickness from the bottom of surface litter layer to to bottom of layer L, [m]
      ! FracSoiAsMicP_vr=micropore fraction
      ! DPTHZ=depth to middle of soil layer from  surface of grid cell [m]
      ! VOLT=total volume
      ! VOLX=total micropore volume

      IF(SoiBulkDensityt0_vr(L,NY,NX).LE.ZERO)SoilFracAsMacP_vr(L,NY,NX)=0.0_r8
      ! layer thickness:=bottom depth-upper depth
      DLYRI_3D(3,L,NY,NX)=(CumDepz2LayBottom_vr(L,NY,NX)-CumDepz2LayBottom_vr(L-1,NY,NX))
      call check_bool(DLYRI_3D(3,L,NY,NX)<0._r8,'negative soil layer thickness',__LINE__,mod_filename)

      DLYR_3D(3,L,NY,NX)           = DLYRI_3D(3,L,NY,NX)
      SoilDepthMidLay_vr(L,NY,NX)  = 0.5_r8*(CumDepz2LayBottom_vr(L,NY,NX)+CumDepz2LayBottom_vr(L-1,NY,NX))
      ! Use the input upper boundary; the active top-layer thickness may not yet
      ! be calculated while visiting reserved layers above NU_col.
      CumSoilThickness_vr(L,NY,NX) = CumDepz2LayBottom_vr(L,NY,NX) &
        -CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX)
      CumSoilThickness_vr(L,NY,NX) = real_truncate(CumSoilThickness_vr(L,NY,NX),1.e-6_r8)
      !identify the first soil layer
      IF(SoiBulkDensityt0_vr(L,NY,NX).LE.ZERO)SoilSurfDepZ_col(NY,NX)=CumSoilThickness_vr(L,NY,NX)
      
      IF(L.EQ.NU_col(NY,NX))THEN
        CumSoilThickMidL_vr(L,NY,NX) = 0.5_r8*CumSoilThickness_vr(L,NY,NX)
      ELSE
        CumSoilThickMidL_vr(L,NY,NX) = 0.5_r8*(CumSoilThickness_vr(L,NY,NX)+CumSoilThickness_vr(L-1,NY,NX))
      ENDIF
      VGeomLayer_vr(L,NY,NX)       = AREA_3D(3,L,NY,NX)*DLYR_3D(3,L,NY,NX)     ![m3]
      VLSoilPoreMicP_vr(L,NY,NX)   = VGeomLayer_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)  !volume of soil occupied by micropores
      VLSoilMicP_vr(L,NY,NX)       = VLSoilPoreMicP_vr(L,NY,NX)
      VGeomLayert0_vr(L,NY,NX)     = VGeomLayer_vr(L,NY,NX)
!     bulk density is defined only for soil with micropores
!     bulk soil mass evaluated as micropore volume
      VLSoilMicPMass_vr(L,NY,NX) = SoilBulkDensity_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)

      totRootLenDens_vr(L,NY,NX) = 0.0_r8
    ENDIF
    AREA_3D(1,L,NY,NX) = DLYR_3D(3,L,NY,NX)*DLYR_3D(2,L,NY,NX)
    AREA_3D(2,L,NY,NX) = DLYR_3D(3,L,NY,NX)*DLYR_3D(1,L,NY,NX)
  ENDDO D111
  CumDepz2LayBottom_vr(0,NY,NX) = CumDepz2LayBottom_vr(NU_col(NY,NX),NY,NX)-DLYR_3D(3,NU_col(NY,NX),NY,NX)
  CumLitRDepzInit_col(NY,NX)    = CumDepz2LayBottom_vr(0,NY,NX)
  AREA_3D(3,NL_col(NY,NX)+1:JZ,NY,NX)  = DLYR_3D(1,NL_col(NY,NX),NY,NX)*DLYR_3D(2,NL_col(NY,NX),NY,NX)

  call PrintInfo('end '//subname)
  end associate
  end subroutine InitLayerDepths

!------------------------------------------------------------------------------------------

  subroutine set_ecosim_solver(NPXS1,NPYS1,NCYC_LITR,NCYC_SNOW)

  implicit none
  integer, intent(in) :: NPXS1,NPYS1
  integer, intent(in) :: NCYC_LITR,NCYC_SNOW
  !     begin_execution
  real(r8) :: XNPV

  NPX = NPXS1   !number of cycles per hour for water, heat, solute flux calcns
  NPY = NPYS1   !number of cycles per NPX for gas flux calcns

  NPH = NPX
  NPT = NPY
  NPG = NPH*NPT

  NPR           = NCYC_LITR     !sub-cycles of litter
  NPS           = NCYC_SNOW     !sub-cycles of snow iteration
  dts_HeatWatTP = 0.99999_r8/NPH
  dt_GasCyc     = 0.99999_r8/NPT
  dts_gas       = 0.99999_r8/NPG
  XNPR          = 0.99999_r8/NPR
  XNPS          = 0.99999_r8/NPS

  XNPV      = XNPR*XNPS
  XNPD      = 600.0_r8*dts_gas                     !600. is adjustable
  dts_wat   = AMIN1(1.0_r8,10.0_r8*dts_HeatWatTP)*0.99999_r8  !adjust/recompute the time step for water/heat update, no greater than 1 hour
  dts_sno   = dts_wat*XNPS
  XNPB      = dts_wat*XNPR                         !vapor flux in litter iteration
  dt_watvap = dts_wat*XNPV

  end subroutine set_ecosim_solver

!------------------------------------------------------------------------------------------
  subroutine InitLandScapeVars
  implicit none
  TAREA = 0.0_r8  !land scape area
  !
  !     INITIALIZE MASS BALANCE CHECKS
  !
  CRAIN_lnd           = 0.0_r8
  HEATIN_lnd          = 0.0_r8
  SurfGas_lnd(:)      = 0.0_r8
  TZIN                = 0.0_r8
  TPIN                = 0.0_r8
  tAmendOrgC_lnd      = 0.0_r8
  TORGN               = 0.0_r8
  TORGP               = 0.0_r8
  QH2OLoss_lnds       = 0.0_r8
  CEVAP               = 0.0_r8
  CRUN                = 0.0_r8
  HeatOut_lnds        = 0.0_r8
  OXYGOU              = 0.0_r8
  H2GOU               = 0.0_r8
  TSedmErossLoss_lnds = 0.0_r8
  TOMOU_lnds(:)       = 0.0_r8
  Litrfall_lnds(:)    = 0.0_r8
  TIONIN              = 0.0_r8
  TIONOU              = 0.0_r8

  end subroutine InitLandScapeVars  

!------------------------------------------------------------------------------------------
  subroutine startsim(NHW,NHE,NVN,NVS)
  use SoilHydroParaMod, only : ComputeSoilHydroPars
  use SoilPhysParaMod, only : SetDeepSoil
  implicit none
  integer,intent(in) :: NHW,NHE,NVN,NVS
  REAL(R8) :: SurfElevLandscapeMax ![m], peak altitude of the landscape
  real(r8) :: SurfElevLandscapeMin
  integer :: NY,NX,NM
  real(r8) :: LandScape1stSoiLayDepth
  real(r8) :: YSIN(NumOfSkyAzimuthSects)
  real(r8) :: YCOS(NumOfSkyAzimuthSects)
  real(r8) :: SkyAzimuthAngle(NumOfSkyAzimuthSects)

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      NM=MaxNumRootLays_col(NY,NX)+1
      !NM=2
      call ComputeSoilHydroPars(NY,NX,NU_col(NY,NX),NM)
      call SetDeepSoil(NY,NX,NM,JZ)
    enddo
  enddo

  call InitLandScapeVars

  !  IRRADIANCE INTERCEPTION GEOMETRY, plant model
  call InitIrradianceGeometry(YSIN,YCOS,SkyAzimuthAngle)
  !  CALCULATE ELEVATION OF EACH GRID CELL
  !
  call InitGridElevation(NHW,NHE,NVN,NVS,YSIN,YCOS,SkyAzimuthAngle,SurfElevLandscapeMax)

  call InitAccumulators()
  call UpdateAtmosGasMassConcentrations(NHW,NHE,NVN,NVS,PressureRatio=1._r8)

  !For whatever reason this sets a temerature, so I need to reset it:
  ATKS_col = 242.0
  SurfElevLandscapeMin=0.0_r8
  LandScape1stSoiLayDepth=0.0_r8   !pay attention to how it is set for many-grid simulations

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      !
      !     MINIMUM SURFACE ELEVATION IN LANDSCAPE
      !
      !     ALT=surface elevation relative to maximum
      !     SurfElevLandscapeMin=minimum surface elevation in landscape
      !CDPTH: depth to bottom of soil layer
      SurfElevation_col(NY,NX)=SurfElevation_col(NY,NX)-SurfElevLandscapeMax
      IF(NX.EQ.NHW .AND. NY.EQ.NVN)THEN
        SurfElevLandscapeMin=SurfElevation_col(NY,NX)
      ELSE
        SurfElevLandscapeMin=MIN(SurfElevLandscapeMin,SurfElevation_col(NY,NX))
      ENDIF
      !
      LandScape1stSoiLayDepth=AMAX1(LandScape1stSoiLayDepth,CumDepz2LayBottom_vr(NU_col(NY,NX),NY,NX)) !topsoil layer depth
!
!     MICROBIAL THERMAL ADAPTATION
!
!     OFFSET=shift in Arrhenius curve used in nitro.f (oC)
!     ATCS=mean annual soil temperature (OC)
!
      TempOffset_col(NY,NX)=f_microbe_TOffset(ATCS_col(NY,NX))

!
!     INITIALIZE WATER POTENTIAL VARIABLES FOR SOIL LAYERS
!
!     LOGPSIMX,PSIMN,LOGPSIAtSat=log water potential at FC,WP,POROS
!     LOGPSIMXD,LOGPSIMND=LOGPSIMX-LOGPSIAtSat,PSIMN-LOGPSIMX
!
      LOGPSIAtSat(NY,NX)   = LOG(-PSIPS)
      LOGPSIFLD_col(NY,NX) = LOG(-PSIAtFldCapacity_col(NY,NX))
      LOGPSIMN_col(NY,NX)  = LOG(-PSIAtWiltPoint_col(NY,NX))

      !for moisture interpolation between field capacity and saturation
      LOGPSIMXD_col(NY,NX)   = LOGPSIFLD_col(NY,NX)-LOGPSIAtSat(NY,NX)
      !for moisture interpolation between wilting point and field capacity
      LOGPSIMND_col(NY,NX)   = LOGPSIMN_col(NY,NX)-LOGPSIFLD_col(NY,NX)

!     VLSoilMicPMass_vr(0,NY,NX)=0.0_r8
!
!     DISTRIBUTION OF OM AMONG FRACTIONS OF DIFFERING
!     BIOLOGICAL ACTIVITY
!      
      call InitHGrid(NY,NX)

      call InitLayerDepths(NY,NX)

    ! ActiveLayDepZ_col=active layer depth (m)
      ActiveLayDepZ_col(NY,NX)=9999.0_r8
!
!     VHeatCapLitR,VHCPNX=minimum heat capacities for solving
!      surface litter,soil layer water and heat fluxes
      VHeatCapLitRMin_col(NY,NX) = VLHeatCapLitRMin*AREA_3D(3,NU_col(NY,NX),NY,NX)
      VHCPNX_col(NY,NX)          = VLHeatCapSoiMin*AREA_3D(3,NU_col(NY,NX),NY,NX)

!
!     SURFACE WATER STORAGE AND LOWER HEAT SINK

!      SoilHeatSrcDepth_col=depth at which soil heat sink-source calculated
!     TCNDG=assumed thermal conductivity below lower soil boundary
!     (MJ m-1 K-1 h-1)
!     TKSD=deep source/sink temperature from geothermal flux(K)

      SoilHeatSrcDepth_col(NY,NX) = AMAX1(10.0_r8,CumDepz2LayBottom_vr(NL_col(NY,NX),NY,NX)+1.0_r8)
      TCS_vr(0,NY,NX) = ATCS_col(NY,NX)
      TKS_vr(0,NY,NX) = ATKS_col(NY,NX)
      TKSD_col(NY,NX) = ATKS_col(NY,NX)+2.052E-04_r8*SoilHeatSrcDepth_col(NY,NX)/TCNDG
      TKS_vr(0,NY,NX) = ATCS_col(NY,NX)+273.15_r8
!
    ENDDO
  ENDDO

!
!     INITIALIZE COMMUNITY CANOPY
!
  CanopyHeight_col(:,:)                        = 0.0_r8
  CanopyHeightZ_col(0:NumCanopyLayers,:,:)   = 0.0_r8
  CanopyLeafAareZ_col(1:NumCanopyLayers,:,:) = 0.0_r8
  CanopyStemAareZ_col(1:NumCanopyLayers,:,:) = 0.0_r8
  tCanLeafC_clyr(1:NumCanopyLayers,:,:)        = 0.0_r8
!
  call InitSoilVars(NHW,NHE,NVN,NVS,SurfElevLandscapeMin,LandScape1stSoiLayDepth)

 !do NY=1,NYS
 !   DO L=NU_col(NY,NX),NL_col(NY,NX)
 !     AREA_3D(3,L,NY,NX)=a_AREA3(L,NY)
 !     write(*,*) "Before InitSnowLayers AREA_3D(3,L,NY,NX) = ", AREA_3D(3,L,NY,NX), ", a_AREA3(L,NY) = ", a_AREA3(L,NY)
 !   ENDDO
 ! ENDDO

!     INITIALIZE SNOWPACK LAYERS
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      call InitSnowLayers(NY,NX)
    ENDDO  
  ENDDO
  end subroutine startsim

end module StartsMod
