module SurfPhysMod
!
!Description
! code for doing surface physics
!
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use abortutils,    only: endrun
  use MiniFuncMod, only: gOC_to_m3_OM
  use DebugToolMod
  use GridDataType
  use HydroThermData
  use CanopyDataType
  use SoilPropertyDataType
  use ClimForcDataType
  use SurfSoilDataType
  use LandSurfDataType
  use SoilWaterDataType
  use SoilHeatDataType
  use SurfLitterDataType
  use MiniFuncMod
  use minimathmod
  use SurfPhysData
  use EcoSIMCtrlDataType
  use SnowDataType
  use SOMDataType
  USE ChemTranspDataType
  use HydroThermData
  use SnowPhysMod
  use PlantTraitDataType
  use SoilPhysDataType
  use SoilBGCDataType
  use IrrigationDataType
  use SnowPhysData
  use PhysPars
  use SurfLitterPhysMod
  use EcoSIMSolverPar
  use EcoSIMCtrlMod
  use SoilPhysParaMod
  use SnowBalanceMod
implicit none
  private
  character(len=*), parameter :: mod_filename=&
  __FILE__
  real(r8), parameter :: tiny_wat=1.e-13_r8

  !surface model
  public :: StageSurfacePhysModel
  public :: CorrectSoilSurfaceHeatM
  public :: RunSurfacePhysModelM
  public :: AggregateSurfRunoffFluxM
  public :: SetHourlyAccumulatorsATS
  public :: writeSurfDiagnosis
  !
  public :: XBoundSurfaceRunoffM
  public :: UpdateSurfaceAtM

! FEnergyImpact4Erosion_col=rate constant for restoring surface Ksat
! SoilEmisivity,SnowEmisivity,SurfLitREmisivity=emissivities of surface soil, snow and litter
! RACX,ResistanceLitrSurf=minimum boundary layer resistances of canopy,litter (h m-1)

  real(r8), parameter :: FEnergyImpact4Erosion_col = 1.0E-03_r8
  real(r8), parameter :: SatHydroCondLitR      = 25.0_r8      !saturated hydraulic conductivity of surface litter [x h-1], unit needs double check
  real(r8), parameter :: ResistanceLitrSurf    = 0.0139_r8    !minimum boundary layer resistances of litter [h m-1]
  real(r8), parameter :: SoilEmisivity         = 0.97_r8      !soil emissivity
  real(r8), parameter :: SnowEmisivity         = 0.97_r8      !snowpack emissivity
  real(r8), parameter :: SurfLitREmisivity     = 0.97_r8      !surfce litter emissivity
  real(r8), parameter :: RACX                  = 0.0139_r8    !total canopy boundary later resistance [h/m]

  real(r8) :: HeatSensVapAir2Grnd
  real(r8) :: HeatSensAir2Grnd
  real(r8) :: Radnet2Grnd   !net radiation onto ground [MJ]
  real(r8) :: LatentHeatEvapAir2Grnd,NetWatFlxAir2SoiMacP
  real(r8) :: CumHeatSensAir2LitR
  real(r8) :: cumHeatSensAir2Soil,NetWatFlxAir2LitR
  real(r8) :: PrecHeat2SoiNet,RainPrecAir2LitR

contains

  subroutine SetHourlyAccumulatorsATS(NY,NX)
!     implicit none
  integer, intent(in) :: NX,NY
  character(len=*), parameter :: subname='SetHourlyAccumulatorsATS'
  integer :: L
!     begin_execution
  call PrintInfo('beg '//subname)
  WatFLo2Litr_col(NY,NX)         = 0.0_r8
  HeatFLoByWat2LitR_col(NY,NX)   = 0.0_r8
  TLitrIceFlxThaw_col(NY,NX)     = 0.0_r8
  TLitrIceHeatFlxFrez_col(NY,NX) = 0.0_r8
  HeatByRad2Surf_col(NY,NX)      = 0.0_r8
  HeatSensAir2Surf_col(NY,NX)    = 0.0_r8
  HeatEvapAir2Surf_col(NY,NX)    = 0.0_r8
  HeatSensVapAir2Surf_col(NY,NX) = 0.0_r8
  HeatNet2Surf_col(NY,NX)        = 0.0_r8
  VapXAir2GSurf_col(NY,NX)       = 0.0_r8

  !TFLWCI(NY,NX)           = 0.0_r8
  RainIntceptByCanopy_col(NY,NX) = 0.0_r8

! zero arrays in the snow layers
  WatConvSno2MicP_snvr(1:JS,NY,NX)   = 0.0_r8
  WatConvSno2MacP_snvr(1:JS,NY,NX)   = 0.0_r8
  HeatConvSno2Soi_snvr(1:JS,NY,NX)   = 0.0_r8
  WatConvSno2LitR_snvr(1:JS,NY,NX)   = 0.0_r8
  HeatConvSno2LitR_snvr(1:JS,NY,NX)  = 0.0_r8
  SnoXfer2SnoLay_snvr(1:JS,NY,NX)    = 0.0_r8
  WatXfer2SnoLay_snvr(1:JS,NY,NX)    = 0.0_r8
  IceXfer2SnoLay_snvr(1:JS,NY,NX)    = 0.0_r8
  HeatXfer2SnoLay_snvr(1:JS,NY,NX)   = 0.0_r8
  XPhaseChangeHeatL_snvr(1:JS,NY,NX) = 0.0_r8
  
! zero plant precip accumulation arrays
  SnowIntceptByCanopy_col(NY,NX)      = 0.0_r8
  RainIntceptByCanopy_col(NY,NX)     = 0.0_r8 
  call PrintInfo('end '//subname)
  end subroutine SetHourlyAccumulatorsATS
!------------------------------------------------------------------------------------------

  subroutine StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistanceLitRLay)

  use SnowPhysMod, only : StageSnowModel
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8),dimension(:,:),intent(OUT) :: ResistanceLitRLay(JY,JX)
  character(len=*), parameter :: subname=trim(mod_filename)//'::StageSurfacePhysModel'
  integer :: NY,NX

  call PrintInfo('beg '//subname)
  watflw =0._r8;waticefl=0._r8

  !be careful about the following, consider move to another location.
  if(ATS_cpl_mode) then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
         NUM_col(NY,NX)=1
      enddo
    enddo
  endif

  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
      !reset accumulators to zero
      TEvapXAir2Toplay_col(NY,NX)       = 0._r8
      TEvapXAir2LitR_col(NY,NX)         = 0._r8
      TEvapXAir2Snow_col(NY,NX)         = 0._r8
      TXGridSurfRunoff_2DH(NY,NX)       = 0.0_r8
      THeatXGridBySurfRunoff_2DH(NY,NX) = 0.0_r8

      !make a local copy of the upper boundary index
      !
      !     ADJUST SURFACE ELEVATION USED IN RUNOFF FOR FREEZE-THAW, EROSION
      !     AND SOC
      !
      !     Altitude_col,ALT=current,initial elevation of ground surface
      !     CumDepz2LayBottom_vr(NUM_col(NY,NX)-1),=depth of ground surface
      !     EnergyImpact4Erosion_col=cumulative rainfall energy impact on soil surface
      !
      Altitude_col(NY,NX)             = SurfElevation_col(NY,NX)-CumDepz2LayBottom_vr(NUM_col(NY,NX)-1,NY,NX)
      EnergyImpact4Erosion_col(NY,NX) = EnergyImpact4Erosion_col(NY,NX)*(1.0_r8-FEnergyImpact4Erosion_col)

      call StageSnowModel(I,J,NY,NX)

      call CopySurfaceVars(I,J,NY,NX)

      call PartionSurfaceFraction(I,J,NY,NX)

      call PartitionPrecip(I,J,NY,NX)

      call SurfaceResistances(I,J,NY,NX,ResistanceLitRLay(NY,NX))

      call Radiation2Surface(I,J,NY,NX)

      call SetCanopyProperty(I,J,NY,NX)
      ! The prior-hour weighted derivative was consumed above. Build this
      ! hour's unweighted conductance sum from the actual exchange substeps.
      SurfSensHeatConduct_col(NY,NX)=0._r8
      SurfLatHeatVapDeriv_col(NY,NX)=0._r8

    ENDDO D9990
  ENDDO D9995
  call PrintInfo('end '//subname)
  end subroutine StageSurfacePhysModel

!------------------------------------------------------------------------------------------

  subroutine CopySurfaceVars(I,J,NY,NX)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX

  character(len=*), parameter :: subname='CopySurfaceVars'
  real(r8) :: VWatLitrZ,TVWatIceLitR,VOLIRZ
!
! SET INITIAL SOIL VALUES
!

  call PrintInfo('beg '//subname)
  LitrIceFlxThaw_col(NY,NX)     = 0.0_r8
  LitrIceHeatFlxFrez_col(NY,NX) = 0.0_r8
!
! ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
!     FOR USE AT INTERNAL TIME STEP IN SURFACE LITTER
!
  m3OM_col(NY,NX) = gOC_to_m3_OM(SoilOrgM_vr(ielmc,0,NY,NX))
  LWRadBySurf_col(NY,NX)        = 0.0_r8
  VHeatCapacity1_vr(0,NY,NX)    = cpo*m3OM_col(NY,NX)+cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)
  VLPoreLitR_col(NY,NX)         = VLMicP_vr(0,NY,NX)
  VLWatMicP1_vr(0,NY,NX)        = AZMAX1(VLWatMicP_vr(0,NY,NX))
  VLiceMicP1_vr(0,NY,NX)        = AZMAX1(VLiceMicP_vr(0,NY,NX))
  VLairMicP1_vr(0,NY,NX)        = AZMAX1(VLPoreLitR_col(NY,NX)-VLWatMicP1_vr(0,NY,NX)-VLiceMicP1_vr(0,NY,NX))
  VLWatMicPM_vr(1,0,NY,NX)      = VLWatMicP1_vr(0,NY,NX)
  VLsoiAirPM_vr(1,0,NY,NX)      = VLairMicP1_vr(0,NY,NX)
  TVWatIceLitR                  = VLWatMicP1_vr(0,NY,NX)+VLiceMicP1_vr(0,NY,NX)
  XVLMobileWaterLitR_col(NY,NX) = AZMAX1(TVWatIceLitR-VWatLitRHoldCapcity_col(NY,NX))
  IF(TVWatIceLitR.GT.ZEROS(NY,NX))THEN
    VWatLitrZ               = VLWatMicP1_vr(0,NY,NX)/TVWatIceLitR*VWatLitRHoldCapcity_col(NY,NX)
    VOLIRZ                  = VLiceMicP1_vr(0,NY,NX)/TVWatIceLitR*VWatLitRHoldCapcity_col(NY,NX)
    XVLMobileWatMicP(NY,NX) = AZMAX1(VLWatMicP1_vr(0,NY,NX)-VWatLitrZ)
    XVLiceMicP_col(NY,NX)   = AZMAX1(VLiceMicP1_vr(0,NY,NX)-VOLIRZ)
  ELSE
    XVLMobileWatMicP(NY,NX) = 0.0_r8
    XVLiceMicP_col(NY,NX)   = 0.0_r8
  ENDIF

  XVLMobileWaterLitRM(1,NY,NX) = XVLMobileWaterLitR_col(NY,NX)
  XVLMobileWatMicPM(1,NY,NX)   = XVLMobileWatMicP(NY,NX)
  XVLiceMicPM(1,NY,NX)         = XVLiceMicP_col(NY,NX)

  IF(VLitR_col(NY,NX).GT.ZEROS2(NY,NX))THEN
    FracSoiPAsWat_vr(0,NY,NX)     = AZMAX1t(VLWatMicP1_vr(0,NY,NX)/VLitR_col(NY,NX))
    FracSoiPAsIce_vr(0,NY,NX)     = AZMAX1t(VLiceMicP1_vr(0,NY,NX)/VLitR_col(NY,NX))
    FracAirFilledSoilPore_vr(0,NY,NX) = AZMAX1t(VLairMicP1_vr(0,NY,NX)/VLitR_col(NY,NX)) * &
      AZMAX1t((1.0_r8-XVLMobileWaterLitR_col(NY,NX)/VLWatHeldCapSurf_col(NY,NX)))
  ELSE
    FracSoiPAsWat_vr(0,NY,NX)     = 0.0_r8
    FracSoiPAsIce_vr(0,NY,NX)     = 0.0_r8
    FracAirFilledSoilPore_vr(0,NY,NX) = 1.0
  ENDIF

  FracAirFilledSoilPoreM_vr(1,0,NY,NX)   = FracAirFilledSoilPore_vr(0,NY,NX)
  PSISM1_vr(0,NY,NX)  = PSISoilMatricP_vr(0,NY,NX)
  TKSoil1_vr(0,NY,NX) = TKS_vr(0,NY,NX)
  
  call PrintInfo('end '//subname)
  end subroutine CopySurfaceVars

!------------------------------------------------------------------------------------------

  subroutine PartionSurfaceFraction(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  real(r8) :: depth
!     SNOW AND RESIDUE COVERAGE OF SOIL SURFACE
!     FSNW,FSNX=fractions of snow,snow-free cover
!     SnowDepth=snowpack depth
!     MinSnowDepth=minimum snowpack depth for full cover
!     BARE,CVRD=fractions of soil,litter cover

  FracSurfAsSnow_col(NY,NX)  = AMIN1(1.0_r8,SQRT((SnowDepth_col(NY,NX)/MinSnowDepth)))
  FracSurfSnoFree_col(NY,NX) = 1.0_r8-FracSurfAsSnow_col(NY,NX)

  call PartLitSoilFractionM(I,J,NY,NX)

  end subroutine PartionSurfaceFraction

!------------------------------------------------------------------------------------------
  subroutine SetCanopyProperty(I,J,NY,NX)
  !
  ! Update canopy vapor pressure and air temperature
  implicit none
  integer, intent(in) :: I,J,NY,NX
  character(len=*), parameter :: subname='SetCanopyProperty'
  real(r8), parameter :: SpecHeatCapAir=1.25E-03_r8    !heat capacity of 1m3 air, [MJ/(m3 K)]
  real(r8) :: RI,RAGX,AirHeatArea,CanopyAirHeatDenom
  real(r8) :: LatHeatArea,CanopyAirVapDenom
  real(r8) :: SensHeatDeriv,LatHeatDeriv
  ! Compile-time slope multiplier: 0 disables stabilization, 1 uses the
  ! full frozen-state estimate; values between 0 and 1 weaken the damping.
  real(r8), parameter :: cfeedback=1._r8

  !TLEX_col=total latent heat flux x boundary layer resistance, [MJ m-1]
  !TSHX_col=total sensible heat flux x boundary layer resistance, [MJ m-1]
  !VPQ_col=vapor pressure in canopy air,
  !TKQ=temperature in canopy air, Kelvin
  call PrintInfo('beg '//subname)
  !
  LatHeatArea              = EvapLHTC*AREA_3D(3,NUM_col(NY,NX),NY,NX)
  ! Keep restart-persisted derivatives unscaled; tune only this prediction.
  LatHeatDeriv            = cfeedback*dLatHeatdVPQ_col(NY,NX)
  SensHeatDeriv           = cfeedback*dSensHeatdTKQ_col(NY,NX)
  CanopyAirVapDenom        = 1._r8+LatHeatDeriv/LatHeatArea

  if(LatHeatDeriv/=0._r8 .and. VPQLatHeatRef_col(NY,NX)>=0._r8 .and. &
    abs(CanopyAirVapDenom)>1.e-2_r8)then
    ! S_L(q) ~= TLEX + K_V*(q-qref), paired from the same completed hour.
    ! This scalar predictor freezes thermal cross-responses; no saturation
    ! clipping or change to accepted water/energy transfers is introduced.
    VPQ_col(NY,NX)=(VPA_col(NY,NX) &
      -(TLEX_col(NY,NX)-LatHeatDeriv*VPQLatHeatRef_col(NY,NX))/LatHeatArea) &
      /CanopyAirVapDenom
  else
    VPQ_col(NY,NX)=VPA_col(NY,NX)-TLEX_col(NY,NX)/LatHeatArea
  endif

  AirHeatArea        = SpecHeatCapAir*AREA_3D(3,NUM_col(NY,NX),NY,NX)
  CanopyAirHeatDenom = 1._r8+SensHeatDeriv/AirHeatArea
  if(SensHeatDeriv/=0._r8 .and. TKQHeatRef_col(NY,NX)>0._r8 .and. &
    abs(CanopyAirHeatDenom)>1.e-2_r8)then
    ! S_H(Tq) ~= TSHX + K_H*(Tq-Tref), with TSHX/K_H/Tref from
    ! the same completed hour: x_new = (x_a-(F-D*x_old)/C)/(1+D/C).
    TKQ_col(NY,NX)=(TairK_col(NY,NX) &
      -(TSHX_col(NY,NX)-SensHeatDeriv*TKQHeatRef_col(NY,NX))/AirHeatArea) &
      /CanopyAirHeatDenom
  else
    ! No completed-hour estimate, zero conductance, or singular denominator.
    TKQ_col(NY,NX)=TairK_col(NY,NX)-TSHX_col(NY,NX)/AirHeatArea
  endif

  RI   = RichardsonNumber(RIB_col(NY,NX),TairK_col(NY,NX),TKQ_col(NY,NX))
  RAGX = RawIsoTAtm2CanopySinkZ_col(NY,NX)

  RawTAtm2CanopySinkZ_col(NY,NX) = AMAX1(RAM,0.8_r8*RAGX,AMIN1(1.2_r8*RAGX,RawIsoTAtm2CanopySinkZ_col(NY,NX)/(1.0_r8-10.0_r8*RI)))

  !
  call PrintInfo('end '//subname)
  end subroutine SetCanopyProperty
!------------------------------------------------------------------------------------------

  subroutine Radiation2Surface(I,J,NY,NX)

  implicit none
  integer, intent(in) :: I,J,NY,NX
  character(len=*), parameter :: subname='Radiation2Surface'

  real(r8) :: tScaledLWRad2Grnd !long-wave radiation incident on ground [MJ]
  real(r8) :: tScaledRAD2Grnd   !short-wave radiation incident on ground, [MJ]
  real(r8) :: Stefboltz_area    !area-scaled stefan-boltzman coefficient

  !
  !     INITIALIZE PARAMETERS, FLUXES FOR ENERGY EXCHANGE
  !     AT SNOW, RESIDUE AND SOIL SURFACES
  !
  !     tScaledRAD2Grnd=shortwave radiation at ground surface
  !     RadSWonSno,RadSW2Soil_col,RadSW2LitR_col= shortwave radn at snowpack,soil,litter, [MJ]
  !     FracSWRad2Grnd_col=fraction of shortwave radiation at ground surface
  !     FSNW,FSNX=fractions of snow,snow-free cover
  !     BARE,CVRD=fractions of soil,litter cover
  !     XNPS=internal time step for fluxes through snowpack
  !     tScaledLWRad2Grnd=longwave radiation at ground surface
  !     LWRad2Snow_col,LWRad2Soil_col,LWRad2LitR_col=longwave radn incident at snowpack,soil,litter
  !     LWEmscefSnow_col,LWEmscefSoil_col,LWEmscefLitR_col=longwave radn emitted by snowpack,soil,litter
  !     SnowEmisivity,SoilEmisivity,SurfLitREmisivity=emissivity of snowpack,soil,litter surfaces
  !     LWRadCanGPrev_col=longwave radiation emitted by canopy
  call PrintInfo('beg '//subname)
  tScaledRAD2Grnd      = RadSWGrnd_col(NY,NX)*dts_HeatWatTP
  RadSW2Sno_col(NY,NX) = tScaledRAD2Grnd*FracSurfAsSnow_col(NY,NX)*XNPS

  tScaledLWRad2Grnd     = (LWRadSky_col(NY,NX)*FracSWRad2Grnd_col(NY,NX)+LWRadCanGPrev_col(NY,NX))*dts_HeatWatTP
  LWRad2Snow_col(NY,NX) = tScaledLWRad2Grnd*FracSurfAsSnow_col(NY,NX)*XNPS

  ! SoilEmisivity,SnowEmisivity,SurfLitREmisivity=emissivities of surface soil, snow and litter
  ! The OTC emissivity modifier applies to canopy emission only.
  !stefboltz_const is stefan-boltzman constant converted into [MJ /(m^2 K^4 h)]
  Stefboltz_area          = stefboltz_const*AREA_3D(3,NUM_col(NY,NX),NY,NX)
  LWEmscefSnow_col(NY,NX) = SnowEmisivity*Stefboltz_area*FracSurfAsSnow_col(NY,NX)*dts_HeatWatTP

  LWEmscefSoil_col(NY,NX) = SoilEmisivity*Stefboltz_area*FracSurfSnoFree_col(NY,NX)*FracSurfBareSoil_col(NY,NX)*dts_HeatWatTP
  LWEmscefLitR_col(NY,NX) = SurfLitREmisivity*Stefboltz_area*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*dts_HeatWatTP

  call PrintInfo('end '//subname)
  end subroutine Radiation2Surface
!------------------------------------------------------------------------------------------

  subroutine LiterSoilRadiationM(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  character(len=*), parameter :: subname='LiterSoilRadiationM'
  real(r8) :: tScaledLWRad2Grnd,Stefboltz_area,tScaledRAD2Grnd

  call PrintInfo('beg '//subname)
  tScaledRAD2Grnd       = RadSWGrnd_col(NY,NX)*dts_HeatWatTP
  RadSW2Soil_col(NY,NX) = tScaledRAD2Grnd*FracSurfSnoFree_col(NY,NX)*FracSurfBareSoil_col(NY,NX)
  RadSW2LitR_col(NY,NX) = tScaledRAD2Grnd*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*XNPR

  tScaledLWRad2Grnd     = (LWRadSky_col(NY,NX)*FracSWRad2Grnd_col(NY,NX)+LWRadCanGPrev_col(NY,NX))*dts_HeatWatTP
  LWRad2Soil_col(NY,NX) = tScaledLWRad2Grnd*FracSurfSnoFree_col(NY,NX)*FracSurfBareSoil_col(NY,NX)
  LWRad2LitR_col(NY,NX) = tScaledLWRad2Grnd*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*XNPR

  Stefboltz_area          = stefboltz_const*AREA_3D(3,NUM_col(NY,NX),NY,NX)
  LWEmscefSoil_col(NY,NX) = SoilEmisivity*Stefboltz_area*FracSurfSnoFree_col(NY,NX)*FracSurfBareSoil_col(NY,NX)*dts_HeatWatTP
  LWEmscefLitR_col(NY,NX) = SurfLitREmisivity*Stefboltz_area*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*dts_HeatWatTP

  call PrintInfo('end '//subname)
  end subroutine LiterSoilRadiationM

!------------------------------------------------------------------------------------------
  subroutine SurfaceResistances(I,J,NY,NX,ResistanceLitRLay)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  real(r8), intent(out):: ResistanceLitRLay     !litter layer resistance [h/m]
  character(len=*), parameter :: subname='SurfaceResistances'

  real(r8) :: FracLitrPoreAsAir       !
  real(r8) :: VaporTortuosityLitr     !porosity limitation to diffusion through litter
  real(r8) :: TFACR
  real(r8) :: PAREX  !area scaled conductances for latent heat fluxes, m^2/h
  real(r8) :: PARSX  !area scaled conductances for sensible heat fluxes, m^2/h
  real(r8) :: RAS    !Snow layer resistance, m
  real(r8) :: ALFZ,WindSpeedGrnd

  call PrintInfo('beg '//subname)
  !
  !     BOUNDARY LAYER CONDUCTANCES FOR EXPORT TO TranspNoSalt.F
  !
  !     AERODYNAMIC RESISTANCE OF CANOPY TO SNOW/RESIDUE/SOIL
  !     SURFACE ENERGY EXCHANGE WITH ATMOSPHERE
  !
  ! ALFZ=attenuation coefficient
  !
  ALFZ=2.0_r8*(1.0_r8-FracSWRad2Grnd_col(NY,NX))
  IF(RawIsoTSurf2CanopyHScal_col(NY,NX).GT.ZERO .AND. CanopyHeight_col(NY,NX).GT.SoilSurfRoughness_col(NY,NX) &
    .AND. ALFZ.GT.ZERO)THEN
    !If there are plants
    WindSpeedGrnd               = WindSpeedAtm_col(NY,NX)*EXP(-ALFZ)
    !Choudhury and Monteith (1988)
    RawIsoTSurf2SinkZ_col(NY,NX) = AMIN1(RACX,AZMAX1(CanopyHeight_col(NY,NX)*EXP(ALFZ) &
      /(ALFZ/RawIsoTSurf2CanopyHScal_col(NY,NX))*AZMAX1(EXP(-ALFZ*SoilSurfRoughness_col(NY,NX)/CanopyHeight_col(NY,NX)) &
      -EXP(-ALFZ*(ZeroPlaneDisplacem_col(NY,NX)+RoughnessLength_col(NY,NX))/CanopyHeight_col(NY,NX)))))
  ELSE
    WindSpeedGrnd               = WindSpeedAtm_col(NY,NX)
    RawIsoTSurf2SinkZ_col(NY,NX) = 0.0_r8
  ENDIF
  !
  !     AERODYNAMIC RESISTANCE OF SNOWPACK, RESIDUE AND SOIL
  !     SURFACES TO ENERGY EXCHANGE WITH ATMOSPHERE
  !     Soil Sci. Soc. Am. J. 48:25-32
  !
  !     RAR=porosity-unlimited litter blr
  !     DLYRR=litter depth
  !     VaporDiffusivityLitR_col=vapor diffusivity in litter
  !     RAG,RAGW,RAGR=isothermal blrs at ground,snowpack,litter surfaces
  !     ResistanceLitrSurf=boundary layer resistance (blr) of litter surface
  !     FracSoiPAsAir*=air-filled porosity of litter
  !     POROQ=litter tortuosity
  !     ResistanceLitRLay=porosity-limited litter blr [h/m]
  !     PAREX,PARSX=
  !     PAREW,PARSW=conductances for snowpack latent,sensible heatfluxes
  !     PAREG,PARSG=conductances for soil latent,sensible heat fluxes
  !     PARER,PARSR=conductances for litter latent,sensible heat fluxes
  !     XNPR=internal time step for fluxes through litter
  !     dts_HeatWatTP=time step for heat calculation in watsub 1/h

  ! Note: VaporDiffusivityLitR_col is computed in Hour1 which is not used in the
  !       ATS coupler. However the calculation is simple so I'm just reproducing it
  !       here. This depends on the temperature of the liter (TKS) which is the same
  !       as TKSoi1 at the surface
  TFACR                           = TEFGASDIF(TKS_vr(0,NY,NX))
  VaporDiffusivityLitR_col(NY,NX) = TFACR*7.70E-02_r8
  VapDiffusResistanceLitR(NY,NX)  = DLYRR_COL(NY,NX)/VaporDiffusivityLitR_col(NY,NX)

  RawIsoTSoil2ATM_col(NY,NX)       = RawIsoTSurf2SinkZ_col(NY,NX)+RawIsoTAtm2CanopySinkZ_col(NY,NX)
  ResistAreodynOverSnow_col(NY,NX) = RawIsoTSoil2ATM_col(NY,NX)
  ResistAreodynOverLitr_col(NY,NX) = RawIsoTSoil2ATM_col(NY,NX)+ResistanceLitrSurf
  RawIsoTLitr2ATM_col(NY,NX)       = ResistAreodynOverLitr_col(NY,NX)
  FracLitrPoreAsAir                = AMAX1(ZERO2,FracAirFilledSoilPore_vr(0,NY,NX))
  VaporTortuosityLitr              = FracLitrPoreAsAir*POROQ*FracLitrPoreAsAir/POROS_vr(0,NY,NX)
  ResistanceLitRLay                = RawIsoTSoil2ATM_col(NY,NX)+VapDiffusResistanceLitR(NY,NX)/VaporTortuosityLitr
  PAREX                            = AREA_3D(3,NUM_col(NY,NX),NY,NX)*dts_HeatWatTP               !conductance for latent heat flux, m^2 h
  ! assuming volumetric heat capacity of air is 1.25E-3_r8 MJ/(m^3 K)
  PARSX                            = 1.25E-03_r8*AREA_3D(3,NUM_col(NY,NX),NY,NX)*dts_HeatWatTP   !conductance for sensible heat flux, MJ h /(m K)

  AScaledCdWOverSnow_col(NY,NX) = PAREX*FracSurfAsSnow_col(NY,NX)*XNPS
  AScaledCdHOverSnow_col(NY,NX) = PARSX*FracSurfAsSnow_col(NY,NX)*XNPS
  AScaledCdWOverSoil_col(NY,NX) = PAREX*FracSurfSnoFree_col(NY,NX)
  AScaledCdWOverLitr_col(NY,NX) = PAREX*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*XNPR
  AScaledCdHOverSoil_col(NY,NX) = PARSX*FracSurfSnoFree_col(NY,NX)
  AScaledCdHOverLitr_col(NY,NX) = PARSX*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)*XNPR
  !
  !     PARR=boundary layer conductance above litter,soil surfaces, (m^2/h)/(h/m)
  !
  RAS             = CalcSnowBNDResistance(NY,NX)
  PARR_col(NY,NX) = AREA_3D(3,NUM_col(NY,NX),NY,NX)*dts_HeatWatTP/(ResistAreodynOverLitr_col(NY,NX)+RAS)   !this includes snow layer resistance
  call PrintInfo('end '//subname)
  end subroutine SurfaceResistances

  !------------------------------------------------------------------------------------------

  subroutine SoilSRFEnerbyBalanceM(M,I,J,NY,NX,PSISV1,LWRadGrnd_col,ResistanceLitRLay,TopLayWatVol,&
    VapXAir2TopLay,HeatFluxAir2Soi)
  !
  !Description
  !Compute surface energy partition
  !Inputs:
  !Topsoil variables:
  !SoilWatAirDry_vr,POROS_vr,VLWatMicP1_vr,VLSoilMicP_vr,PSISoilOsmotic_vr
  !TKSoil1_vr,NUM
  !Atmospheric variables:
  !VPQ_col,TKQ_col
  implicit none
  integer, intent(in) :: M,I,J,NY,NX
  real(r8), intent(out):: PSISV1      !surface soil suction pressure [MPa]
  real(r8), intent(out) :: LWRadGrnd_col  !outgoing long wave radiation from soil surface [MJ]
  real(r8), intent(inout) :: ResistanceLitRLay
  real(r8), intent(inout) :: TopLayWatVol
  real(r8), intent(out) :: VapXAir2TopLay
  real(r8), intent(out) :: HeatFluxAir2Soi   !total heat flux into soil after surface energy budget [MJ]
  
  real(r8) :: RAa            ! [h/m]
  real(r8) :: VaporSoi1      ! soil vapor pressure [ton water/m3]
  real(r8) :: CdSoiEvap      ! scaled conductance for evaporation   [m^3]
  real(r8) :: RSoilEvap      ! soil surface resistance to evaporation [h/m]
  real(r8) :: tRadIncid      !total incoming radiation to soil or liter surface, short + long [MJ]
  real(r8) :: RI,WPLX,WPX,FracLitrPoreAsAir
  real(r8) :: VLWatGrnd,VLIceGrnd,VaporTortuosityLitr
  real(r8) :: TKX1,THETA1S
  real(r8) :: tHeatAir2Grnd   !residual heat flux into soil from incoming radiation minus sensible and latent heat [MJ]
  real(r8) :: AlbedoGrnd      !albedo at the ground
  real(r8) :: RadSWbySoil     !shortwave radiation absorbed by exposed soil [MJ]
  real(r8) :: ResistBndlSurf0
  character(len=*), parameter :: subname='SoilSRFEnerbyBalanceM'

  ! begin_execution
  !
  ! PHYSICAL AND HYDRAULIC PROPERTIES OF SOIL SURFACE INCLUDING
  ! AIR AND WATER-FILLED POROSITY, AND WATER POTENTIAL USED IN
  ! FLUX CALCULATIONS
  !
  ! THETY=current,hygroscopic water concentration
  ! POROS=porosity
  ! VLWatMicP1,VLSoilPoreMicP_vrI=volume of micropore water
  ! VLSoilPoreMicP_vrI=soil volume less macropore,rock
  ! FC,WP=water contents at field capacity,wilting point
  ! FCL,WPL=log FC,WP
  ! FCD,PSD=FCL-WPL,log(POROS)-FCL
  ! PSISM1,PSIHY,PSISE=soil matric,hygroscopic,air entry potential
  ! LOGPSIMX,LOGPSIMN,LOGPSIAtSat=log water potential at FC,WP,POROS
  ! PSISD,PSIMD=LOGPSIMX-LOGPSIAtSat,LOGPSIMN-LOGPSIMX
  ! SRP=parameter for deviation from linear log-log water retention
  ! function from hour1.f
  ! PSISO=osmotic potential
  ! BKVL=bulk soil mass for a given layer

  call PrintInfo('beg '//subname)

  !In ATS mode PSISM1_vr already holds the ATS matric potential; recomputing
  !it from EcoSIM's own water content and retention curve hides surface drying
  !from the soil vapor pressure and keeps evaporation at the potential rate
  if(.not.ats_cpl_mode)then
    call CalcSoilWatPotential(NY,NX,NX,NY,NUM_col(NY,NX),PSISM1_vr(NUM_col(NY,NX),NY,NX),THETA1S)
  endif

  PSISV1=PSISM1_vr(NUM_col(NY,NX),NY,NX)+PSISoilOsmotic_vr(NUM_col(NY,NX),NY,NX)
!
! SOIL SURFACE ALBEDO, NET RADIATION
!
! VLWatMicP1,VLiceMicP1=water,ice volume in micopores
! VLWatMacP1,VLiceMacP1=water,ice volume in macopores
! AlbedoGrnd,SoilAlbedo=albedo of ground surface,soil
! BKVL=soil mass
! RadSW2Soil_col,LWRad2Soil_col,Radnet2Grnd=incoming shortwave,longwave,net radiation
! LWRadGrnd_col,LWEmscefSoil_col=emitted longwave radiation, emissivity
! TK1=soil temperature
! albedo of water and ice are set to 0.06, and 0.30 respectively

  VLWatGrnd=VLWatMicP1_vr(NUM_col(NY,NX),NY,NX)+VLWatMacP1_vr(NUM_col(NY,NX),NY,NX)
  VLIceGrnd=VLiceMicP1_vr(NUM_col(NY,NX),NY,NX)+VLiceMacP1_vr(NUM_col(NY,NX),NY,NX)

  IF((VLWatGrnd+VLIceGrnd).GT.ZEROS2(NY,NX))THEN
    !top soil has water or ice
    !ice albedo seems too low.
    !write(*,*) "Albedo recomputed"
    AlbedoGrnd=(SoilAlbedo_col(NY,NX)*VLSoilMicPMass_vr(NUM_col(NY,NX),NY,NX)+0.06_r8*VLWatGrnd &
      +0.30_r8*VLIceGrnd)/(VLSoilMicPMass_vr(NUM_col(NY,NX),NY,NX)+VLWatGrnd+VLIceGrnd)
  ELSE
    !write(*,*) "Albedo from soil"
    AlbedoGrnd=SoilAlbedo_col(NY,NX)
  ENDIF
  !absorbed radiation
  !Radnet2Grnd=net radiation, after taking out outgoing surface layer radiation
  !LWRadGrnd_col=emitted longwave radiation
  RadSWbySoil          = (1.0_r8-AlbedoGrnd)*RadSW2Soil_col(NY,NX)
  tRadIncid            = RadSWbySoil+LWRad2Soil_col(NY,NX)
  LWRadGrnd_col        = LWEmscefSoil_col(NY,NX)*TKSoil1_vr(NUM_col(NY,NX),NY,NX)**4
  Radnet2Grnd          = tRadIncid-LWRadGrnd_col
  Eco_RadSW_col(NY,NX) = Eco_RadSW_col(NY,NX) + RadSWbySoil
!
! AERODYNAMIC RESISTANCE ABOVE SOIL SURFACE INCLUDING
! RESISTANCE IMPOSED BY PLANT CANOPY
!
! FracSoiPAsAir*=air-filled porosity of soil
! VaporTortuosityLitr=porosity limitation to diffusion through soil
! POROQ=soil tortuosity
! ResistanceLitRLay=porosity-limited litter blr
! RAGZ=combined soil+litter blr
! RI=Richardsons number
! RIB=isothermal RI
! TKQ,TK1=canopy air,soil temperature
! RAGZ,RAa=soil+litter blr
! ResistBndlSurf_col=isothermal blr at ground surface
!
  FracLitrPoreAsAir   = AMAX1(ZERO,FracAirFilledSoilPore_vr(0,NY,NX))
  VaporTortuosityLitr = FracLitrPoreAsAir*POROQ*FracLitrPoreAsAir/POROS_vr(0,NY,NX)
  ResistanceLitRLay   = RawIsoTSoil2ATM_col(NY,NX)+VapDiffusResistanceLitR(NY,NX)/VaporTortuosityLitr
  RI                  = RichardsonNumber(RIB_col(NY,NX),TKQ_col(NY,NX),TKSoil1_vr(NUM_col(NY,NX),NY,NX))

  ResistBndlSurf0  = 1.0_r8/((1._r8-FracSurfByLitR_col(NY,NX))/ResistAreodynOverLitr_col(NY,NX)+FracSurfByLitR_col(NY,NX)/ResistanceLitRLay)

!  ResistBndlSurf_col(NY,NX) = AMAX1(RAM,0.8_r8*ResistBndlSurf_col(NY,NX),AMIN1(1.2_r8*ResistBndlSurf_col(NY,NX), &
!    ResistanceLitRLay/(1.0_r8-10.0_r8*RI)))

  ResistBndlSurf_col(NY,NX) = AMAX1(RAM,0.8_r8*ResistBndlSurf0,AMIN1(1.2_r8*ResistBndlSurf0, ResistanceLitRLay/(1.0_r8-10.0_r8*RI)))
  RAa                       = ResistAreodynOverLitr_col(NY,NX)+ResistBndlSurf_col(NY,NX)
  !
  ! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
  !
  ! CdSoiEvap,CdSoiHSens=conductance for latent,sensible heat fluxes over soil
  ! PAREG,PARSG=conductances for latent,sensible heat fluxes
  ! RZ=minimum surface resistance
  ! VaporSoi1,VPQ_col=vapor pressure at soil surface, canopy air
  ! VapXAir2TopLay=negative of evaporation, or water vapor transfer from canopy air to top layer soil/lake
  ! LatentHeatEvapAir2Grnd=latent heat flux
  ! XH=rate constant
  ! VOLW2=soil water volume
  ! VAP=latent heat of evaporation
  ! HeatSensVapAir2Grnd=convective heat of evaporation flux
  !
  if(ats_cpl_mode)then
    !A coarse ATS top cell cannot resolve the dry surface crust that limits
    !bare-soil evaporation, so use a surface resistance from its liquid
    !saturation (THETW_vr = a_LSAT): r = exp(8.206-4.255*s) [s m-1] (Sellers et
    !al., 1992, SiB2). It is about RZ when saturated and ~1500 s m-1 at s=0.2.
    RSoilEvap = AMAX1(RZ, EXP(8.206_r8-4.255_r8*AMIN1(1._r8,AZMAX1(THETW_vr(NUM_col(NY,NX),NY,NX))))/3600._r8)
  else
    RSoilEvap = RZ
  endif
  CdSoiEvap             = AScaledCdWOverSoil_col(NY,NX)/(RAa+RSoilEvap)
  CdSoiHSens_col(NY,NX) = AScaledCdHOverSoil_col(NY,NX)/RAa
  TKX1                  = TKSoil1_vr(NUM_col(NY,NX),NY,NX)
  VaporSoi1             = vapsat(TKX1)*EXP(18.0_r8*PSISV1/(RGASC*TKX1))

  !evaporation, no more than what is available, g H2O
  if(TopLayWatVol<1.0e-30) TopLayWatVol=0.0
  VapXAir2TopLay=AMAX1(CdSoiEvap*(VPQ_col(NY,NX)-VaporSoi1),-AZMAX1(TopLayWatVol*dts_wat))
  ! Freeze water availability: a donor-limited evaporation branch has zero
  ! local derivative. CdSoiEvap already includes area and substep duration.
  if(CdSoiEvap*(VPQ_col(NY,NX)-VaporSoi1)>-AZMAX1(TopLayWatVol*dts_wat))then
    SurfLatHeatVapDeriv_col(NY,NX)=SurfLatHeatVapDeriv_col(NY,NX)+EvapLHTC*CdSoiEvap
  endif
  !VapXAir2TopLay=0.0

  !latent heat > 0 into soil/ground
  LatentHeatEvapAir2Grnd=VapXAir2TopLay*EvapLHTC
  IF(VapXAir2TopLay.LT.0.0_r8)THEN
    !evaporation (<0 into atmosphere)
    HeatSensVapAir2Grnd=VapXAir2TopLay*cpw*TKSoil1_vr(NUM_col(NY,NX),NY,NX)*HeatAdv_scal
  ELSE
    !condensation (>0 into ground)
    HeatSensVapAir2Grnd=VapXAir2TopLay*cpw*TKQ_col(NY,NX)*HeatAdv_scal
  ENDIF
  !
  TopLayWatVol=TopLayWatVol+VapXAir2TopLay
  !
  !
  ! SOLVE FOR SOIL SURFACE TEMPERATURE AT WHICH ENERGY
  ! BALANCE OCCURS, SOLVE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
  !
  ! HeatSensAir2Grnd,LatentHeatEvapAir2Grnd,Radnet2Grnd=sensible,latent heat fluxes, net radiation
  ! HeatSensVapAir2Grnd=convective heat flux from LatentHeatEvapAir2Grnd, > 0 into ground
  ! HeatFluxAir2Soi=storage heat flux
  !
  HeatSensAir2Grnd = CdSoiHSens_col(NY,NX)*(TKQ_col(NY,NX)-TKSoil1_vr(NUM_col(NY,NX),NY,NX))
  SurfSensHeatConduct_col(NY,NX)=SurfSensHeatConduct_col(NY,NX)+CdSoiHSens_col(NY,NX)
  tHeatAir2Grnd    = Radnet2Grnd+LatentHeatEvapAir2Grnd+HeatSensAir2Grnd !net energy into soil, subtracting latent heat and sensible heat
  HeatFluxAir2Soi  = tHeatAir2Grnd+HeatSensVapAir2Grnd !total heat plus convective heat > 0 to ground
  
  HeatSensAir2Grnd_col(NY,NX) = HeatSensAir2Grnd_col(NY,NX) + HeatSensAir2Grnd
  Radnet2Grnd_col(NY,NX) = Radnet2Grnd_col(NY,NX) + Radnet2Grnd
  LatentHeatEvapAir2Grnd_col(NY,NX) = LatentHeatEvapAir2Grnd_col(NY,NX) + LatentHeatEvapAir2Grnd
  HeatSensVapAir2Grnd_col(NY,NX) = HeatSensVapAir2Grnd_col(NY,NX) + HeatSensVapAir2Grnd
  
  call PrintInfo('end '//subname)

  end subroutine SoilSRFEnerbyBalanceM

!------------------------------------------------------------------------------------------

  subroutine ExposedSoilFluxM(I,J,M,NY,NX,CumSnowWatFlow2LitR,CumSnoHeatFlow2LitR,&
    CumWatFlx2SoiMicP,CumWatFlx2SoiMacP,cumSnoHeatFlow2Soil,Prec2LitR2,RainHeat2LitR2,&
    PrecNet2SoiMicP,PrecNet2SoiMacP, RainPrecHeatAir2LitR,ResistanceLitRLay,TopLayWatVol,&
    VapXAir2TopLay,HeatFluxAir2Soi,NetWatFlxAir2SoiMicP,NetWatXFlxAir2SoiMicP,NetWatFlxAir2SoiMacP)
  implicit none
  integer, intent(in) :: M,I,J   !soil heat-flow iteration id
  integer, intent(in) :: NY,NX
  real(r8), intent(in)  :: CumSnowWatFlow2LitR,CumSnoHeatFlow2LitR
  real(r8), intent(in) :: CumWatFlx2SoiMicP,CumWatFlx2SoiMacP,cumSnoHeatFlow2Soil
  real(r8), intent(in) :: Prec2LitR2,RainHeat2LitR2
  real(r8), intent(in) :: PrecNet2SoiMicP
  real(r8), intent(in) :: PrecNet2SoiMacP
  real(r8), intent(in) :: RainPrecHeatAir2LitR
  real(r8) ,intent(inout) :: ResistanceLitRLay
  real(r8), intent(inout) :: TopLayWatVol   !water content in topsoil layer
  real(r8), intent(out) :: VapXAir2TopLay
  real(r8), intent(out) :: HeatFluxAir2Soi    !heat flux from air to top soil [MJ]
  real(r8), intent(out) :: NetWatFlxAir2SoiMacP
  real(r8), intent(out) :: NetWatXFlxAir2SoiMicP  !
  real(r8), intent(out) :: NetWatFlxAir2SoiMicP

  character(len=*), parameter :: subname='ExposedSoilFluxM'
  real(r8) :: VapXAir2LitR
  real(r8) :: PSISV1
  real(r8) :: LWRadGrnd_col
  real(r8) :: HeatSensLitR2Soi1,HeatSensVapLitR2Soi1,EvapLitR2Soi1,HeatFluxAir2LitR

! begin_execution
  call PrintInfo('beg '//subname)

! Watch out for L, is its value defined?
  call SoilSRFEnerbyBalanceM(M,I,J,NY,NX,PSISV1,LWRadGrnd_col,ResistanceLitRLay,TopLayWatVol,&
    VapXAir2TopLay,HeatFluxAir2Soi)
  !
  call SurfLitREnergyBalanceM(I,J,M,NY,NX,PSISV1,TopLayWatVol,Prec2LitR2,RainHeat2LitR2,CumSnowWatFlow2LitR,&
    CumSnoHeatFlow2LitR,CumWatFlx2SoiMicP,CumWatFlx2SoiMacP,cumSnoHeatFlow2Soil,&
    HeatSensLitR2Soi1,HeatSensVapLitR2Soi1,EvapLitR2Soi1,VapXAir2LitR,HeatFluxAir2LitR)
  !
  call SumAftEnergyBalanceM(I,J,M,NY,NX,LWRadGrnd_col,VapXAir2TopLay,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,EvapLitR2Soi1,VapXAir2LitR,HeatFluxAir2LitR,HeatFluxAir2Soi,PrecNet2SoiMicP,&
    PrecNet2SoiMacP,RainPrecHeatAir2LitR,NetWatFlxAir2SoiMicP,NetWatXFlxAir2SoiMicP,NetWatFlxAir2SoiMacP)

  call PrintInfo('end '//subname)
  end subroutine ExposedSoilFluxM

!------------------------------------------------------------------------------------------

  subroutine AtmLandSurfExchangeM(I,J,M,NY,NX,Prec2LitR2,RainHeat2LitR2,PrecNet2SoiMicP,&
    PrecNet2SoiMacP,RainPrecHeatAir2LitR,ResistanceLitRLay,TopLayWatVol,LatentHeatAir2Sno,&
    HeatSensEvapAir2Snow,HeatSensAir2Snow,Radnet2Snow,VapXAir2TopLay,HeatFluxAir2Soi1)
  implicit none
  integer, intent(in) :: M           !soil heat-flow iteration id
  integer, intent(in) :: NY,NX,I,J
  real(r8), intent(in) :: Prec2LitR2,RainHeat2LitR2
  real(r8), intent(in):: PrecNet2SoiMicP
  real(r8), intent(in):: PrecNet2SoiMacP
  real(r8), intent(in):: RainPrecHeatAir2LitR
  real(r8), intent(inout) :: ResistanceLitRLay
  real(r8), intent(inout) :: TopLayWatVol
  real(r8), intent(out) :: LatentHeatAir2Sno,HeatSensAir2Snow,Radnet2Snow,HeatSensEvapAir2Snow
  real(r8), intent(out) :: VapXAir2TopLay
  real(r8), intent(out) :: HeatFluxAir2Soi1   !MJ/d2

  character(len=*), parameter :: subname='AtmLandSurfExchangeM'
  real(r8) :: cumSnoHeatFlow2Soil
  real(r8) :: CumWatFlx2SoiMicP,CumWatFlx2SoiMacP
  real(r8) :: CumSnoHeatFlow2LitR
  real(r8) :: CumWatXFlx2SoiMicP
  real(r8) :: CumSnowWatFlow2LitR
  real(r8) :: HeatNetFlx2Snow
  real(r8) :: NetWatXFlxAir2SoiMicP
  real(r8) :: NetWatFlxAir2SoiMicP
  real(r8) :: NetWatFlxAir2SoiMacP
  real(r8) :: TMX,deltaT
  integer  :: L

! begin_execution
  call PrintInfo('beg '//subname)
  VapXAir2TopLay       = 0._r8
  CumWatXFlx2SoiMicP   = 0._r8
  CumSnowWatFlow2LitR  = 0._r8
  CumSnoHeatFlow2LitR  = 0._r8
  CumWatFlx2SoiMicP    = 0._r8
  CumWatFlx2SoiMacP    = 0._r8
  HeatSensEvapAir2Snow = 0._r8
  HeatNetFlx2Snow      = 0._r8
  LatentHeatAir2Sno    = 0._r8
  HeatSensAir2Snow     = 0._r8
  Radnet2Snow          = 0._r8
  cumSnoHeatFlow2Soil  = 0._r8
  HeatFluxAir2Soi1     = 0._r8

  !solve for energy balance over significant snow layer
  IF(VLSnowHeatCapM_snvr(M,1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
    !   VHCPW,VLHeatCapSnowMin_col=current, minimum snowpack heat capacities
    call SolveSnowpackM(I,J,M,NY,NX,LatentHeatAir2Sno,Radnet2Snow,HeatSensEvapAir2Snow,HeatSensAir2Snow,&
      HeatNetFlx2Snow,CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,CumWatXFlx2SoiMicP,CumSnowWatFlow2LitR,&
      CumSnoHeatFlow2LitR,cumSnoHeatFlow2Soil)
  ENDIF

  !AGGREGATE fluxes from solving the snow model
  WaterFlow2Micpt_3D(3,NUM_col(NY,NX),NY,NX)  = CumWatFlx2SoiMicP
  WaterFlow2MicptX_3D(3,NUM_col(NY,NX),NY,NX) = CumWatXFlx2SoiMicP
  WaterFlow2Macpt_3D(3,NUM_col(NY,NX),NY,NX)  = CumWatFlx2SoiMacP
  HeatFlow2Soili_3D(3,NUM_col(NY,NX),NY,NX)   = cumSnoHeatFlow2Soil
  WatFLo2LitrM_col(NY,NX)                     = CumSnowWatFlow2LitR
  HeatFLoByWat2LitRM_col(NY,NX)               = CumSnoHeatFlow2LitR
  HeatSnow2LitrM_col(NY,NX)                   = CumSnoHeatFlow2LitR
  !!!!!
  !
  ! ENERGY EXCHANGE AT SOIL SURFACE IF EXPOSED UNDER SNOWPACK
  !
  IF(FracSurfSnoFree_col(NY,NX).GT.0.0_r8 .AND. (SoilBulkDensity_vr(NUM_col(NY,NX),NY,NX).GT.ZERO .OR. &
    VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX).GT.VHCPNX_col(NY,NX)))THEN

    !Ground partially covered by snow, focus on litter-soil interaction
    call ExposedSoilFluxM(I,J,M,NY,NX,CumSnowWatFlow2LitR,CumSnoHeatFlow2LitR,&
      CumWatFlx2SoiMicP,CumWatFlx2SoiMacP,cumSnoHeatFlow2Soil,Prec2LitR2,RainHeat2LitR2,&
      PrecNet2SoiMicP,PrecNet2SoiMacP,RainPrecHeatAir2LitR,ResistanceLitRLay,TopLayWatVol,&
      VapXAir2TopLay,HeatFluxAir2Soi1,NetWatFlxAir2SoiMicP,NetWatXFlxAir2SoiMicP,NetWatFlxAir2SoiMacP)
  ELSE
!   ground is fully snow covered, thus no flux from soil & litter
    call ZeroFlxOverSnowCoveredSoilM(I,J,M,NY,NX,NetWatFlxAir2SoiMicP,NetWatXFlxAir2SoiMicP,NetWatFlxAir2SoiMacP)
  ENDIF

!
  if(lverb)write(*,*)'AGGREGATE RESIDUE AND SOIL SURFACE FLUXES BENEATH SNOW AND ATMOSPHERE'
!
! WaterFlow2Micpt_3D,WaterFlow2MicptX_3D=total water flux into soil micropores
! FLWHL=total water flux into soil macropores
! HFLWL=total heat flux into soil
! FLWRL,WaterFlow2MicptX_3D=total water flux into litter
! HFLWRL=total heat flux into litter
! FLWV*=total internal vapor flux in soil
!
  WaterFlow2Micpt_3D(3,NUM_col(NY,NX),NY,NX)  = WaterFlow2Micpt_3D(3,NUM_col(NY,NX),NY,NX)+NetWatFlxAir2SoiMicP
  WaterFlow2MicptX_3D(3,NUM_col(NY,NX),NY,NX) = WaterFlow2MicptX_3D(3,NUM_col(NY,NX),NY,NX)+NetWatXFlxAir2SoiMicP
  WaterFlow2Macpt_3D(3,NUM_col(NY,NX),NY,NX)  = WaterFlow2Macpt_3D(3,NUM_col(NY,NX),NY,NX)+NetWatFlxAir2SoiMacP
  HeatFlow2Soili_3D(3,NUM_col(NY,NX),NY,NX)   = HeatFlow2Soili_3D(3,NUM_col(NY,NX),NY,NX)+cumHeatSensAir2Soil
  !
  !fluxes to litter
  WatFLo2LitrM_col(NY,NX)       = WatFLo2LitrM_col(NY,NX)+NetWatFlxAir2LitR
  HeatFLoByWat2LitRM_col(NY,NX) = HeatFLoByWat2LitRM_col(NY,NX)+CumHeatSensAir2LitR
  HeatAir2LitrM_col(NY,NX)      = CumHeatSensAir2LitR
  call PrintInfo('end '//subname)
  end subroutine AtmLandSurfExchangeM
!------------------------------------------------------------------------------------------

  subroutine ZeroFlxOverSnowCoveredSoilM(I,J,M,NY,NX,NetWatFlxAir2SoiMicP,NetWatXFlxAir2SoiMicP,NetWatFlxAir2SoiMacP)
  !
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NY,NX
  real(r8),intent(out):: NetWatFlxAir2SoiMicP
  real(r8),intent(out):: NetWatXFlxAir2SoiMicP
  real(r8),intent(out):: NetWatFlxAir2SoiMacP

! begin_execution
  Radnet2Grnd            = 0.0_r8        !net radiation into soil
  LatentHeatEvapAir2Grnd = 0.0_r8        !latent heat flux from air and topsoil
  HeatSensVapAir2Grnd    = 0.0_r8        !convective heat flux from air and topsoil
  HeatSensAir2Grnd       = 0.0_r8        !sensible heat flux from air to topsoil

  NetWatFlxAir2SoiMicP  = 0.0_r8       !total water flux from air into top soil
  NetWatXFlxAir2SoiMicP = 0.0_r8       !total water flux from air into soil micropores
  NetWatFlxAir2SoiMacP  = 0.0_r8       !total water flux from air into macropores
  cumHeatSensAir2Soil   = 0.0_r8       !total water associated heat flux from air into soil
  NetWatFlxAir2LitR     = 0.0_r8       !total water flux from air to litter
  CumHeatSensAir2LitR   = 0.0_r8       !total water associated heat flux from air to litter
  CdLitRHSens_col(NY,NX)= 0._r8       !no exposed-litter sensible exchange beneath snow
  end subroutine ZeroFlxOverSnowCoveredSoilM
!------------------------------------------------------------------------------------------

  subroutine SurfLitrSoilWaterExchangeM(I,J,M,NY,NX,RainEkReducedKsat)
  implicit none
  integer, intent(in) :: I,J,M,NY,NX
  real(r8),intent(in) :: RainEkReducedKsat
  character(len=*), parameter :: subname='SurfLitrSoilWaterExchangeM'
  real(r8) :: THETW1,ThetaWLitR,PSIST1
  real(r8) :: PSIST0,HeatFlxLitR2Soi,FLQZ,DarcyFlxLitR2Soil
  real(r8) :: WatDarcyFloLitR2SoiMicP,FLQ2,CNDR,DarcyCondLitR2Soil
  real(r8) :: HeatFlowLitR2MacP,WatFlowLitR2MacP
  real(r8) :: AvailableLitRWater,AvailableSoilWater,AvailableMicPSpace,AvailableMacPSpace
  real(r8) :: MobileLitRWater,TotalLitRWaterIce
  real(r8) :: CND1
  integer :: K0,K1
  real(r8) :: TKRt,VHCapLitR,ENGYR,VLWatLitR,deltaT,TMX

! begin_execution
! CNDR,HCNDLitr_col=current,saturated litter hydraulic conductivity
! PSISE,PSISM1_vr(0,=air entry,current litter water potential
! VLWatMicP1_vr(0,VWatLitRHoldCapcity=current,maximum litter water volume
! CND1,HydroCond_3D=soil hydraulic conductivity
! RainEkReducedKsat=reduction in soil surface Ksat from rainfall energy impact
! K1=soil relative water-filled porosity
! THETWX,POROS=soil water content,porosity
! DarcyCondLitR2Soil=litter-soil hydraulic conductance
! DLYRR,DLYR=litter,soil thicknesses
! WatDarcyFloLitR2SoiMicP,FLQL=litter-soil water flux unltd,ltd by water
! dts_HeatWatTP=time step of flux calculations
! CVRD=fraction of litter cover
! PSISM1_vr(NUM=soil matric water potential
! VLWatMicP1_vr(NUM=soil water volume
! HFLQL=convective heat from litter-soil water flux
! FLWL,HFLWL=micropore water,heat flux
! WatFLow2LitR,HFLWRL=total litter water,heat flux
! FLWRM=litter-soil water flux for solute transfer in TranspNoSalt.f
! CND1,CNDL=hydraulic conductivity of source,destination layer
! HydroCond_3D=lateral(1,2),vertical(3) micropore hydraulic conductivity
!
  call PrintInfo('beg '//subname)
  HeatSoil2LitrM_col(NY,NX) = 0._r8
  ! Budget accepted surface fluxes before either infiltration pathway. Keep
  ! micropore exchange first; macropores receive only remaining mobile liquid.
  WatFLoLitr2SoilM_col(M,NY,NX)=0._r8
  WatFLoLitr2MacPM_col(M,NY,NX)=0._r8
  AvailableLitRWater=AZMAX1(VLWatMicP1_vr(0,NY,NX)+WatFLo2LitrM_col(NY,NX))
  VLWatLitR=AvailableLitRWater
  AvailableSoilWater=AZMAX1(VLWatMicP1_vr(NUM_col(NY,NX),NY,NX) &
    +WaterFlow2Micpt_3D(3,NUM_col(NY,NX),NY,NX))
  AvailableMicPSpace=AZMAX1(VLairMicP1_vr(NUM_col(NY,NX),NY,NX) &
    -WaterFlow2Micpt_3D(3,NUM_col(NY,NX),NY,NX))

  IF(SoilBulkDensity_vr(NUM_col(NY,NX),NY,NX).GT.ZERO)THEN
    !top layer is soil
    IF(VWatLitRHoldCapcity_col(NY,NX).GT.ZEROS2(NY,NX))THEN
      !surface litter holds water
      ThetaWLitR=AMIN1(VWatLitRHoldCapcity_col(NY,NX),VLWatLitR)/VLitR_col(NY,NX)
    ELSE
      ThetaWLitR=POROS0_col(NY,NX)
    ENDIF
    THETW1=AMAX1(SoilWatAirDry_vr(NUM_col(NY,NX),NY,NX),AMIN1(POROS_vr(NUM_col(NY,NX),NY,NX) &
      ,safe_adb(VLWatMicP1_vr(NUM_col(NY,NX),NY,NX),VLSoilMicP_vr(NUM_col(NY,NX),NY,NX))))
    !K0 litter layer
    !K1 topsoil layer
    !DarcyFlxLitR2Soil = water flux from litter layer into the topsoil
    
    K0   = getMoistK(ThetaWLitR,POROS0_col(NY,NX))
    K1   = getMoistK(THETW1,POROS_vr(NUM_col(NY,NX),NY,NX))
    CNDR = HydroCond_3D(3,K0,0,NY,NX)
    CND1 = HydroCond_3D(3,K1,NUM_col(NY,NX),NY,NX)*RainEkReducedKsat
    if(ats_cpl_mode .or. CNDR.LE.0._r8 .or. CND1.LE.0._r8)then
      DarcyCondLitR2Soil = 0._r8
    else
      DarcyCondLitR2Soil = 2.0_r8*CNDR*CND1/(CNDR*DLYR_3D(3,NUM_col(NY,NX),NY,NX)+CND1*DLYRR_COL(NY,NX))
    endif
    PSIST0             = PSISM1_vr(0,NY,NX)+PSIGrav_vr(0,NY,NX)+PSISoilOsmotic_vr(0,NY,NX)
    PSIST1             = PSISM1_vr(NUM_col(NY,NX),NY,NX)+PSIGrav_vr(NUM_col(NY,NX),NY,NX)+PSISoilOsmotic_vr(NUM_col(NY,NX),NY,NX)
    DarcyFlxLitR2Soil  = DarcyCondLitR2Soil*(PSIST0-PSIST1)*AREA_3D(3,NUM_col(NY,NX),NY,NX)*FracEffAsLitR_col(NY,NX)*dts_HeatWatTP

    IF(DarcyFlxLitR2Soil.GE.0.0_r8)THEN
      !flow from liter to soil
      !more water than saturation
      !note dts_wat:=AMIN1(1.0_r8,20.0_r8*dts_HeatWatTP)
      IF(ThetaWLitR.GT.ThetaSat_vr(0,NY,NX))THEN
        FLQZ=DarcyFlxLitR2Soil+AMIN1((ThetaWLitR-ThetaSat_vr(0,NY,NX))*VLitR_col(NY,NX),&
          AZMAX1((ThetaSat_vr(NUM_col(NY,NX),NY,NX)-THETW1)*VLSoilMicP_vr(NUM_col(NY,NX),NY,NX)))*dts_wat
      ELSE
        FLQZ=DarcyFlxLitR2Soil
      ENDIF
      WatDarcyFloLitR2SoiMicP = AZMAX1(AMIN1(FLQZ,VLWatLitR*dts_wat,AvailableMicPSpace))
      FLQ2                    = AZMAX1(AMIN1(DarcyFlxLitR2Soil,VLWatLitR*dts_wat,AvailableMicPSpace))
    ELSE
      !from soil to liter
      IF(THETW1.GT.ThetaSat_vr(NUM_col(NY,NX),NY,NX))THEN
        FLQZ=DarcyFlxLitR2Soil+AMAX1((ThetaSat_vr(NUM_col(NY,NX),NY,NX)-THETW1)*VLSoilMicP_vr(NUM_col(NY,NX),NY,NX),&
          AZMIN1((ThetaWLitR-ThetaSat_vr(0,NY,NX))*VLitR_col(NY,NX)))*dts_wat
      ELSE
        FLQZ=DarcyFlxLitR2Soil
      ENDIF
      WatDarcyFloLitR2SoiMicP = AZMIN1(AMAX1(FLQZ,-AvailableSoilWater*dts_wat,-VLairMicP1_vr(0,NY,NX)))
      FLQ2                    = AZMIN1(AMAX1(DarcyFlxLitR2Soil,-AvailableSoilWater*dts_wat,-VLairMicP1_vr(0,NY,NX)))
    ENDIF

    !top soil layer has negative air volume
    IF(VLairMicP_vr(NUM_col(NY,NX),NY,NX).LT.0.0_r8)THEN
      WatDarcyFloLitR2SoiMicP = WatDarcyFloLitR2SoiMicP+AZMIN1(AMAX1(-AvailableSoilWater*dts_wat,VLairMicP_vr(NUM_col(NY,NX),NY,NX)))
      FLQ2                    = FLQ2+AZMIN1(AMAX1(-AvailableSoilWater*dts_wat,VLairMicP_vr(NUM_col(NY,NX),NY,NX)))
    ENDIF

    ! The oversaturation correction above shares the same donor budget.
    WatDarcyFloLitR2SoiMicP=AMAX1(-AvailableSoilWater, &
      AMIN1(AvailableLitRWater,WatDarcyFloLitR2SoiMicP))
    IF(WatDarcyFloLitR2SoiMicP.GT.0.0_r8)THEN
      !litter layer to soil
      HeatFlxLitR2Soi=cpw*TKSoil1_vr(0,NY,NX)*WatDarcyFloLitR2SoiMicP*HeatAdv_scal
    ELSE
      !soil to litter layer
      HeatFlxLitR2Soi=cpw*TKSoil1_vr(NUM_col(NY,NX),NY,NX)*WatDarcyFloLitR2SoiMicP*HeatAdv_scal
    ENDIF

    WaterFlow2Micpt_3D(3,NUM_col(NY,NX),NY,NX) = WaterFlow2Micpt_3D(3,NUM_col(NY,NX),NY,NX)+WatDarcyFloLitR2SoiMicP
    HeatFlow2Soili_3D(3,NUM_col(NY,NX),NY,NX)  = HeatFlow2Soili_3D(3,NUM_col(NY,NX),NY,NX)+HeatFlxLitR2Soi

    !soil to litter
  ELSE
    !top layer is water
    TotalLitRWaterIce=AvailableLitRWater+VLiceMicP1_vr(0,NY,NX)
    MobileLitRWater=0._r8
    IF(TotalLitRWaterIce.GT.ZEROS(NY,NX))MobileLitRWater=AvailableLitRWater &
      *AZMAX1(1._r8-VWatLitRHoldCapcity_col(NY,NX)/TotalLitRWaterIce)
    WatDarcyFloLitR2SoiMicP=AMIN1(AvailableLitRWater,MobileLitRWater*dts_wat)
    HeatFlxLitR2Soi                            = cpw*TKSoil1_vr(0,NY,NX)*WatDarcyFloLitR2SoiMicP*HeatAdv_scal
    WaterFlow2Micpt_3D(3,NUM_col(NY,NX),NY,NX) = WaterFlow2Micpt_3D(3,NUM_col(NY,NX),NY,NX)+WatDarcyFloLitR2SoiMicP
    HeatFlow2Soili_3D(3,NUM_col(NY,NX),NY,NX)  = HeatFlow2Soili_3D(3,NUM_col(NY,NX),NY,NX)+HeatFlxLitR2Soi
  ENDIF

  AvailableLitRWater=AZMAX1(AvailableLitRWater-WatDarcyFloLitR2SoiMicP)
  WatFLo2LitrM_col(NY,NX)       = WatFLo2LitrM_col(NY,NX)-WatDarcyFloLitR2SoiMicP
  HeatFLoByWat2LitRM_col(NY,NX) = HeatFLoByWat2LitRM_col(NY,NX)-HeatFlxLitR2Soi
  WatFLoLitr2SoilM_col(M,NY,NX) = WatDarcyFloLitR2SoiMicP
  HeatSoil2LitrM_col(NY,NX)     = HeatSoil2LitrM_col(NY,NX) -HeatFlxLitR2Soi
!     OVERLAND FLOW INTO SOIL MACROPORES WHEN WATER STORAGE CAPACITY
!     OF THE LITTER IS EXCEEDED
!
!     VOLPH1=air-filled macroporosity
!     FINHR,HFINHR=water,convective heat from litter to macropores
!     VLWatMicP1_vr(0,VWatLitRHoldCapcity=current,maximum litter water volume
!     WaterFlow2Micpt_3D,HFLWL=micropore water,heat flux
!     WatFLow2LitR,HFLWRL=total litter water,heat flux
!
  TotalLitRWaterIce=AvailableLitRWater+VLiceMicP1_vr(0,NY,NX)
  MobileLitRWater=0._r8
  IF(TotalLitRWaterIce.GT.ZEROS(NY,NX))MobileLitRWater=AvailableLitRWater &
    *AZMAX1(1._r8-VWatLitRHoldCapcity_col(NY,NX)/TotalLitRWaterIce)
  AvailableMacPSpace=AZMAX1(VLairMacP1_vr(NUM_col(NY,NX),NY,NX) &
    -WaterFlow2Macpt_3D(3,NUM_col(NY,NX),NY,NX))
  IF(AvailableMacPSpace.GT.0._r8 .AND. MobileLitRWater.GT.0._r8)THEN
    WatFlowLitR2MacP=AMIN1(MobileLitRWater*dts_wat,AvailableLitRWater,AvailableMacPSpace)
    WatFLoLitr2MacPM_col(M,NY,NX)=WatFlowLitR2MacP
    HeatFlowLitR2MacP                      = WatFlowLitR2MacP*cpw*TKSoil1_vr(0,NY,NX)*HeatAdv_scal
    WaterFlow2Macpt_3D(3,NUM_col(NY,NX),NY,NX) = WaterFlow2Macpt_3D(3,NUM_col(NY,NX),NY,NX)+WatFlowLitR2MacP
    HeatFlow2Soili_3D(3,NUM_col(NY,NX),NY,NX)  = HeatFlow2Soili_3D(3,NUM_col(NY,NX),NY,NX)+HeatFlowLitR2MacP

    WatFLo2LitrM_col(NY,NX)       = WatFLo2LitrM_col(NY,NX)-WatFlowLitR2MacP
    HeatFLoByWat2LitRM_col(NY,NX) = HeatFLoByWat2LitRM_col(NY,NX)-HeatFlowLitR2MacP
    HeatSoil2LitrM_col(NY,NX)     = HeatSoil2LitrM_col(NY,NX) -HeatFlowLitR2MacP
  ENDIF


  call PrintInfo('end '//subname)
  end subroutine SurfLitrSoilWaterExchangeM
!------------------------------------------------------------------------------------------

  subroutine UpdateLitterPhaseChangeM(I,J,M,NY,NX)
  !
  !Calculate litter phase-change fluxes and water-film thicknesses
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NY,NX
  character(len=*), parameter :: subname='UpdateLitterPhaseChangeM'
  real(r8) :: LitrIceHeatFlxFrezPt,TK1X,ENGYR,VLWatMicP1X,VLHeatCapacityX
  real(r8) :: TFREEZ,TFLX


  !     begin_execution
  !
  !     FREEZE-THAW IN RESIDUE SURFACE FROM NET CHANGE IN RESIDUE
  !     SURFACE HEAT STORAGE
  !
  !     TFREEZ=litter freezing temperature
  !     PSISM1=litter water potential
  !     VLWatMicP1*,VLiceMicP1=litter water,ice volume
  !     VLHeatCapacity*=litter volumetric heat capacity
  !     TK1*=litter temperature
  !     ORGC_vr=litter organic C
  !     HFLWRL=total litter conductive, convective heat flux
  !     LitrIceHeatFlxFrezPt,TFLX=unltd,ltd latent heat from freeze-thaw
  !     LitrIceHeatFlxFrez,LitrIceFlxThaw=litter water,latent heat flux from freeze-thaw
  ! using Clausius-Clapeyron equation for freezing temeprature calculation
  call PrintInfo('beg '//subname)
  !TFREEZ          = -9.0959E+04_r8/(PSISM1_vr(0,NY,NX)-LtHeatIceMelt)
  TFREEZ          = get_Tfrez(PSISM1_vr(0,NY,NX))
  VLWatMicP1X     = AZMAX1(VLWatMicP1_vr(0,NY,NX)+WatFLo2LitrM_col(NY,NX))
  ENGYR           = VHeatCapacity1_vr(0,NY,NX)*TKSoil1_vr(0,NY,NX)
  VLHeatCapacityX = cpo*m3OM_col(NY,NX)+cpw*VLWatMicP1X+cpi*VLiceMicP1_vr(0,NY,NX)

  IF(VLHeatCapacityX.GT.ZEROS(NY,NX))THEN
    TK1X=(ENGYR+HeatFLoByWat2LitRM_col(NY,NX))/VLHeatCapacityX
  ELSE
    TK1X=TKSoil1_vr(0,NY,NX)
  ENDIF

  !do freeze thaw of litter water
  IF((TK1X.LT.TFREEZ .AND. VLWatMicP1_vr(0,NY,NX).GT.ZERO*VGeomLayer_vr(0,NY,NX)) &
    .OR. (TK1X.GT.TFREEZ .AND. VLiceMicP1_vr(0,NY,NX).GT.ZERO*VGeomLayer_vr(0,NY,NX)))THEN
    LitrIceHeatFlxFrezPt = VHeatCapacity1_vr(0,NY,NX)*(TFREEZ-TK1X) &
      /((1.0_r8+(TFREEZ-TFICE)*6.2913E-03_r8)*(1.0_r8-0.10_r8*PSISM1_vr(0,NY,NX)))*dts_wat

    IF(LitrIceHeatFlxFrezPt.LT.0.0_r8)THEN
      !ice thaw, current temperature above freezing TFREEZ<TK1X
      TFLX=AMAX1(-LtHeatIceMelt*DENSICE*VLiceMicP1_vr(0,NY,NX)*dts_wat,LitrIceHeatFlxFrezPt)
    ELSE
      !water freeze (TFLX > 0._r8)
      TFLX=AMIN1(LtHeatIceMelt*VLWatMicP1X*dts_wat,LitrIceHeatFlxFrezPt)
    ENDIF

    LitrIceHeatFlxFrez_col(NY,NX) = TFLX
    LitrIceFlxThaw_col(NY,NX)     = -TFLX/LtHeatIceMelt
  ELSE
    LitrIceFlxThaw_col(NY,NX)     = 0.0_r8
    LitrIceHeatFlxFrez_col(NY,NX) = 0.0_r8
  ENDIF

!
!     THICKNESS OF WATER FILMS IN LITTER AND SOIL SURFACE
!     FROM WATER POTENTIALS FOR GAS EXCHANGE IN TranspNoSalt.F
!
  IF(VHeatCapacity1_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX))THEN
    FILMM_vr(M,0,NY,NX)=FilmThickness(PSISM1_vr(0,NY,NX), is_top_layer=.true.)
  ELSE
    FILMM_vr(M,0,NY,NX)=1.0E-03_r8
  ENDIF
  FILMM_vr(M,NUM_col(NY,NX),NY,NX)=FilmThickness(PSISM1_vr(NUM_col(NY,NX),NY,NX),is_top_layer=.true.)
  call PrintInfo('end '//subname)
  end subroutine UpdateLitterPhaseChangeM
!------------------------------------------------------------------------------------------

  subroutine CalcSurfaceRunoffPotentialM(I,J,M,NY,NX)
  !Calculate runoff potential from the locally updated, pre-runoff state.
  implicit none
  integer, intent(in) :: I,J,M,NY,NX
  character(len=*), parameter :: subname='CalcSurfaceRunoffPotentialM'
  integer :: N1,N2
  real(r8) :: VWatExcess,WatExcessDetph,Q,VLWatMicP1X
  real(r8) :: HydraulicRadius,CrossSectVelocity

  call PrintInfo('beg '//subname)
  !
  !     OVERLAND FLOW WHEN WATER STORAGE CAPACITY
  !     OF THE SOIL SURFACE PLUS MACROPORES IS EXCEEDED
  !
  N1=NX;N2=NY
  !
  !     SURFACE WATER FLUX
  !
  !     N2,N1=NY,NX of source grid cell
  !     XVOLT,XVOLW=excess water+ice,water in source grid cell
  !     VWatStoreCapSurf_col=ground surface water retention capacity
  !     VWatExcess=ponded water volume above surface retention capacity
  !     Q=runoff from Mannings equation

  ! there is mobile water
  IF(XVLMobileWaterLitR_col(N2,N1).GT.VWatStoreCapSurf_col(N2,N1))THEN
    VWatExcess                      = XVLMobileWaterLitR_col(N2,N1)-VWatStoreCapSurf_col(N2,N1)
    WatExcessDetph                  = VWatExcess/AREA_3D(3,0,N2,N1)
    HydraulicRadius                 = WatExcessDetph/2.828_r8
    CrossSectVelocity               = GaucklerManningVelocity(HydraulicRadius,SLOPE_col(0,N2,N1))/SoiSurfRoughness(N2,N1)     ![1/s]
    Q                               = CrossSectVelocity*WatExcessDetph*AREA_3D(3,NUM_col(N2,N1),N2,N1)*secsphour*dts_HeatWatTP  ![kg/h/d2]
    VLWatMicP1X                     = AZMAX1(VLWatMicP1_vr(0,N2,N1))
    RunoffVelocityM_col(M,N2,N1)    = CrossSectVelocity
    SurfRunoffPotentM_col(M,N2,N1) = AMIN1(Q,VWatExcess*dts_wat,VLWatMicP1X*dts_wat) &
      *XVLMobileWatMicP(N2,N1)/XVLMobileWaterLitR_col(N2,N1)                             !>0
  ELSE
    RunoffVelocityM_col(M,N2,N1)   = 0.0_r8
    SurfRunoffPotentM_col(M,N2,N1) = 0.0_r8
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine CalcSurfaceRunoffPotentialM
!------------------------------------------------------------------------------------------

  subroutine AccumWaterVaporHeatFluxesM(I,J,M,NY,NX,LatentHeatAir2Sno,HeatSensEvapAir2Snow,HeatSensAir2Snow,&
    Radnet2Snow,VapXAir2TopLay)
  implicit none
  integer , intent(in) :: I,J
  integer , intent(in) :: M,NY,NX
  real(r8), intent(in) :: LatentHeatAir2Sno,Radnet2Snow
  real(r8), intent(in) :: HeatSensAir2Snow,HeatSensEvapAir2Snow
  real(r8), intent(in) :: VapXAir2TopLay

  character(len=*), parameter :: subname='AccumWaterVaporHeatFluxesM'
  real(r8) :: VLWatLitR,WaterRoundoffTol

  ! begin_execution
  ! HOURLY-ACCUMULATED WATER, VAPOR AND HEAT FLUXES THROUGH
  ! SURFACE RESIDUE AND SOIL SURFACE
  !

  call PrintInfo('beg '//subname)

  VLWatLitR=VLWatMicP1_vr(0,NY,NX)+WatFLo2LitrM_col(NY,NX)+LitrIceFlxThaw_col(NY,NX)

  ! A material deficit is a flux-budget error; never create water to hide it.
  WaterRoundoffTol=MAX(tiny_wat,100._r8*EPSILON(1._r8)*MAX( &
    ABS(VLWatMicP1_vr(0,NY,NX)),ABS(WatFLo2LitrM_col(NY,NX)),ABS(LitrIceFlxThaw_col(NY,NX))))
  IF(VLWatLitR.LT.-WaterRoundoffTol)THEN
    write(*,*)subname,': negative litter liquid at I,J,M,NY,NX=',I,J,M,NY,NX,VLWatLitR
    call endrun(trim(mod_filename)//' at line',__LINE__)
  ENDIF

  TLitrIceFlxThaw_col(NY,NX)                   = TLitrIceFlxThaw_col(NY,NX)+LitrIceFlxThaw_col(NY,NX)
  TLitrIceHeatFlxFrez_col(NY,NX)               = TLitrIceHeatFlxFrez_col(NY,NX)+LitrIceHeatFlxFrez_col(NY,NX)
  WaterFlowSoiMicP_3D(3,NUM_col(NY,NX),NY,NX)  = WaterFlowSoiMicP_3D(3,NUM_col(NY,NX),NY,NX)+WaterFlow2Micpt_3D(3,NUM_col(NY,NX),NY,NX)
  WaterFlowSoiMicPX_3D(3,NUM_col(NY,NX),NY,NX) = WaterFlowSoiMicPX_3D(3,NUM_col(NY,NX),NY,NX)+WaterFlow2MicptX_3D(3,NUM_col(NY,NX),NY,NX)
  WaterFlowSoiMacP_3D(3,NUM_col(NY,NX),NY,NX)  = WaterFlowSoiMacP_3D(3,NUM_col(NY,NX),NY,NX)+WaterFlow2Macpt_3D(3,NUM_col(NY,NX),NY,NX)
  HeatFlow2Soil_3D(3,NUM_col(NY,NX),NY,NX)     = HeatFlow2Soil_3D(3,NUM_col(NY,NX),NY,NX)+HeatFlow2Soili_3D(3,NUM_col(NY,NX),NY,NX)

  WatFLo2LitR_col(NY,NX)         = WatFLo2LitR_col(NY,NX)+WatFLo2LitrM_col(NY,NX)
  HeatFLoByWat2LitR_col(NY,NX)   = HeatFLoByWat2LitR_col(NY,NX)+HeatFLoByWat2LitRM_col(NY,NX)
  HeatByRad2Surf_col(NY,NX)      = HeatByRad2Surf_col(NY,NX)+Radnet2Grnd+Radnet2Snow
  HeatSensAir2Surf_col(NY,NX)    = HeatSensAir2Surf_col(NY,NX)+HeatSensAir2Grnd+HeatSensAir2Snow
  HeatEvapAir2Surf_col(NY,NX)    = HeatEvapAir2Surf_col(NY,NX)+LatentHeatEvapAir2Grnd+LatentHeatAir2Sno
  HeatSensVapAir2Surf_col(NY,NX) = HeatSensVapAir2Surf_col(NY,NX)+HeatSensVapAir2Grnd+HeatSensEvapAir2Snow

  HeatNet2Surf_col(NY,NX)                 = HeatNet2Surf_col(NY,NX)          &
    +Radnet2Grnd+HeatSensAir2Grnd+LatentHeatEvapAir2Grnd+HeatSensVapAir2Grnd &
    +Radnet2Snow+HeatSensAir2Snow+LatentHeatAir2Sno+HeatSensEvapAir2Snow

  !EVAPG=negative evaporation from ground/top soil layer
  !EVAPR=evaporation from litter layer
  !EVAPSN=evaporation from snow, sublimation+evaporation

  VapXAir2GSurf_col(NY,NX)                 = VapXAir2GSurf_col(NY,NX)+VapXAir2TopLay+VapXAir2Sno_col(NY,NX)   !>0 into ground
  WaterFlow2MicPM_3D(M,3,NUM_col(NY,NX),NY,NX) = WaterFlow2Micpt_3D(3,NUM_col(NY,NX),NY,NX)
  WaterFlow2MacPM_3D(M,3,NUM_col(NY,NX),NY,NX) = WaterFlow2Macpt_3D(3,NUM_col(NY,NX),NY,NX)
  TEvapXAir2Toplay_col(NY,NX)              = TEvapXAir2Toplay_col(NY,NX)+VapXAir2TopLay
  TEvapXAir2Snow_col(NY,NX)                = TEvapXAir2Snow_col(NY,NX)+VapXAir2Sno_col(NY,NX)

  call PrintInfo('end '//subname)
  end subroutine AccumWaterVaporHeatFluxesM

!------------------------------------------------------------------------------------------

  subroutine InitSurfModelM(I,J,M,NY,NX,ResistanceLitRLay,RainEkReducedKsat,Prec2LitR2,RainHeat2LitR2,&
    PrecNet2SoiMicP,PrecNet2SoiMacP,RainPrecHeatAir2LitR)
  !
  !Initialize surface model for iteration M
  implicit none
  integer, intent(in) :: M   !soil heat-flow iteration id
  integer, intent(in) :: NY,NX,I,J
  real(r8),intent(in) :: ResistanceLitRLay
  real(r8),intent(out):: Prec2LitR2,RainHeat2LitR2
  real(r8),intent(out):: RainEkReducedKsat
  real(r8),intent(out):: PrecNet2SoiMicP
  real(r8),intent(out):: PrecNet2SoiMacP
  real(r8),intent(out):: RainPrecHeatAir2LitR
  integer  :: L
  real(r8) :: scalar,THETWT
  real(r8) :: VOLAT0,ENGYD
  real(r8) :: Rain2MicPflow2LitR    !rainfall overflow from micropore to litter
  real(r8) :: Rain2MacPflow2LitR    !rainfall overflow from macropore to litter
  real(r8) :: RainHeat2Soiflow2LitR !rainfall heat overflow from soil pores to litter
  real(r8) :: ENGYB,RAS,TScal4Aquadifsvity,THETWA
  real(r8) :: MobileH2OSurf   !mm H2O exceed surface water holding capacity

! begin_execution
! INITIALIZE NET SURFACE FLUX ACCUMULATORS
!
! TQS1,TQW1,TQI1=net water and snowpack snow,water,ice runoff
! cumHeatFlx2LitRByRunoff_col,THQS1=net convective heat from surface water and snow runoff
! FracAsExposedSoil_col,FracEffAsLitR_col=fractions of soil,litter cover including free water+ice
! ResistBndlSurf_col= boundary layer resistance at soil surface
! CondGasXSnowM_col=boundary layer conductance above soil surface
!
  call ZeroSnowFluxM(NY,NX)

  call PartLitSoilFractionM(I,J,NY,NX)

  call LiterSoilRadiationM(I,J,NY,NX)

  IF(VHeatCapacity1_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX))THEN
    FracAsExposedSoil_col(NY,NX)=AZMAX1(FracSurfBareSoil_col(NY,NX)-AMIN1(1.0_r8,&
      AZMAX1(XVLMobileWaterLitR_col(NY,NX)/VLWatHeldCapSurf_col(NY,NX))))
  ELSE
    FracAsExposedSoil_col(NY,NX)=1.0_r8
  ENDIF
  FracEffAsLitR_col(NY,NX)   = 1.0_r8-FracAsExposedSoil_col(NY,NX)

  ResistBndlSurf_col(NY,NX)  = 1.0_r8/(FracAsExposedSoil_col(NY,NX)/ResistAreodynOverLitr_col(NY,NX)+FracEffAsLitR_col(NY,NX)/ResistanceLitRLay)

  RAS                        = CalcSnowBNDResistance(NY,NX)
  CondGasXSnowM_col(M,NY,NX) = AREA_3D(3,NUM_col(NY,NX),NY,NX)*dts_HeatWatTP/(ResistBndlSurf_col(NY,NX)+RAS)  !m^2 h/(h/m) = m3
  CondGasXSurf_col(NY,NX)    = 1._r8/(ResistBndlSurf_col(NY,NX)+RAS)  !m/h
  !
  ! REDISTRIBUTE INCOMING PRECIPITATION
  ! BETWEEN RESIDUE AND SOIL SURFACE
  !
  ! The model always assumes there is a litter layer
  IF(SoilBulkDensity_vr(NUM_col(NY,NX),NY,NX).GT.ZERO)THEN
    !soil bulk density significant
    !get flow to litter
    if(FracSurfByLitR_col(NY,NX)>ZEROL)then
      Rain2MicPflow2LitR    = AZMAX1(Rain2SoiMicP1_col(NY,NX)-VLairMicP1_vr(NUM_col(NY,NX),NY,NX))
      Rain2MacPflow2LitR    = AZMAX1(Rain2SoiMacP1_col(NY,NX)-VLairMacP1_vr(NUM_col(NY,NX),NY,NX))
      RainHeat2Soiflow2LitR = cpw*TairK_col(NY,NX)*(Rain2MicPflow2LitR+Rain2MacPflow2LitR)
      RainPrecAir2LitR      = Rain2LitR1_col(NY,NX)+Rain2MicPflow2LitR+Rain2MacPflow2LitR
      RainPrecHeatAir2LitR  = RainHeat2LitR1_col(NY,NX)+RainHeat2Soiflow2LitR
    else
      Rain2MicPflow2LitR    = 0._r8
      Rain2MacPflow2LitR    = 0._r8
      RainHeat2Soiflow2LitR = 0._r8
      RainPrecAir2LitR      = 0._r8
      RainPrecHeatAir2LitR  = 0._r8
    endif
    !get flow to soil
    PrecNet2SoiMicP  = Rain2SoiMicP1_col(NY,NX)-Rain2MicPflow2LitR
    PrecNet2SoiMacP  = Rain2SoiMacP1_col(NY,NX)-Rain2MacPflow2LitR
    PrecHeat2SoiNet  = RainHeat2SoilP1_col(NY,NX)-RainHeat2Soiflow2LitR
  ELSE
    RainPrecAir2LitR     = Rain2LitR1_col(NY,NX)
    RainPrecHeatAir2LitR = RainHeat2LitR1_col(NY,NX)
    PrecNet2SoiMicP      = Rain2SoiMicP1_col(NY,NX)
    PrecNet2SoiMacP      = Rain2SoiMacP1_col(NY,NX)
    PrecHeat2SoiNet      = RainHeat2SoilP1_col(NY,NX)
  ENDIF

  Rain2LitR_col(NY,NX) = Rain2LitR_col(NY,NX)+RainPrecAir2LitR
  Rain2Soil_col(NY,NX) = Rain2Soil_col(NY,NX)+PrecNet2SoiMicP+PrecNet2SoiMacP

  Prec2LitR2             = RainPrecAir2LitR*XNPR
  RainHeat2LitR2         = RainPrecHeatAir2LitR*XNPR
  !
  ! WATER GAS EXCHANGE COEFFICIENTS IN SURFACE LITTER
  !
  ! VOLA1,VLiceMicP1,VLWatMicP1,VsoiPM=total,ice-,water-,air-filled porosity
  ! TScal4Aquadifsvity=temperature effect on gas diffusivity
  ! DiffusivitySolutEff=rate constant for air-water gas exchange
  ! Z1R,Z2RW,Z2RD,Z3RX=parameters for litter air-water gas transfers
  ! XNPD=time step for gas transfer calculations, it is tunable parameter
  ! TORT=tortuosity for aqueous diffusivity
  ! VOLAT0=ice-excluded porosity in litter
  !
  VOLAT0=VLPoreLitR_col(NY,NX)-VLiceMicP1_vr(0,NY,NX)
  IF(VOLAT0.GT.ZEROS2(NY,NX).AND.VLsoiAirPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    !litter layer is not saturated
    THETWA                             = AZMAX1(AMIN1(1.0_r8,VLWatMicP1_vr(0,NY,NX)/VOLAT0))
    TScal4Aquadifsvity                 = TEFAQUDIF(TKSoil1_vr(0,NY,NX))
    scalar                             = TScal4Aquadifsvity*XNPD
    DiffusivitySolutEffM_vr(M,0,NY,NX) = fDiffusivitySolutEff(scalar,THETWA,0.0_r8,is_litter=.true.)
  ELSE
    !litter layer saturated
    DiffusivitySolutEffM_vr(M,0,NY,NX)=0.0_r8
  ENDIF
! VWatLitRHoldCapcity=surface litter water holding capacity, [m3 d-2]
  IF(VWatLitRHoldCapcity_col(NY,NX).GT.ZEROS(NY,NX))THEN
    !litter is able to hold water
    THETWT=AMIN1(1.0_r8,VLWatMicP_vr(0,NY,NX)/VWatLitRHoldCapcity_col(NY,NX))
  ELSE
    THETWT=1.0
  ENDIF
  !TORT=tortuosity in litter (treated as micropore)
  TortMicPM_vr(M,0,NY,NX)=TortMicporeW(THETWT)
  !
  ! KINETIC ENERGY OF DIRECT RAINFALL AND THROUGHFALL
  !
  ! PRECD,PRECB=direct,indirect precipn+irrign at soil surface
  ! ENGYD,ENGYB=energy impact of direct,indirect precipn+irrign at soil surface
  ! VWatStoreCapSurf_col=ground surface water retention capacity
  ! XVOLW=free surface water
  ! ZT=canopy height
  ! EnergyImpact4Erosion_colM=total energy impact for use in erosion.f
  ! EnergyImpact4Erosion_col=cumulative rainfall energy impact on soil surface
  ! RainEkReducedKsat=reduction in soil surface Ksat from rainfall energy impact
  ! Note: A good reference for the following formula and alternatives
  ! is "Rainfall intensity-kinetic energy relationships for soil loss prediction",
  ! Kinnell, 1981

  IF(PrecDirect2Grnd_col(NY,NX).GT.ZERO)THEN
    ENGYD=AZMAX1(8.95_r8+8.44_r8*LOG(PRECM_col(NY,NX)))
  ELSE
    ENGYD=0.0_r8
  ENDIF
  IF(PrecIndirect2Grnd_col(NY,NX).GT.ZERO)THEN
    ENGYB=AZMAX1(15.8_r8*SQRT(AMIN1(2.5_r8,CanopyHeight_col(NY,NX)))-5.87_r8)
  ELSE
    ENGYB=0.0_r8
  ENDIF

  IF(ENGYD+ENGYB.GT.ZERO)THEN
    MobileH2OSurf                  = AMIN1(1.0E+03_r8*AZMAX1(XVLMobileWaterLitR_col(NY,NX)-VWatStoreCapSurf_col(NY,NX))/AREA_3D(3,NUM_col(NY,NX),NY,NX),20._r8)
    EnergyImpact4Erosion_colM(M,NY,NX) = (ENGYD*PrecDirect2Grnd_col(NY,NX)+ENGYB*PrecIndirect2Grnd_col(NY,NX))*EXP(-2.0_r8*MobileH2OSurf) &
      *FracSurfBareSoil_col(NY,NX)*dts_HeatWatTP
    EnergyImpact4Erosion_col(NY,NX)=EnergyImpact4Erosion_col(NY,NX)+AZERO(EnergyImpact4Erosion_colM(M,NY,NX))
  ELSE
    EnergyImpact4Erosion_colM(M,NY,NX)=0.0_r8
  ENDIF
  RainEkReducedKsat=EXP(-2.0E-03_r8*(CSILT_vr(NUM_col(NY,NX),NY,NX)+CCLAY_vr(NUM_col(NY,NX),NY,NX)) &
    *EnergyImpact4Erosion_col(NY,NX))
  !
  !  SNOWPACK FLUX ACCUMULATORS
  !
   call InitSnowAccumsM(I,J,M,NY,NX)

  call PrepIterSnowLayerM(I,J,M,NY,NX)

!
  end subroutine InitSurfModelM
!------------------------------------------------------------------------------------------

  subroutine XBoundSurfaceRunoffM(I,J,M,N,NN,N1,N2,M4,M5,RCHQF,XN)
  !
  !Description:
  !Surface runoff X the Northern, western, southern, and eastern boundaries,
  !It is driven by surface ponding water elevation difference
  !On the east or south, (M5, M4) is off the boundary.
  !On the west or north, (M5,M4) is next to the boundary facing the east or south.
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M      !iteration index
  integer, intent(in) :: N      !Flow direction, N=1: west-east, N=2:north-south
  integer, intent(in) :: NN     !front/behind indicator
  integer, intent(in) :: N1,N2  !source grid, always on the boundary
  integer, intent(in) :: M4,M5  !dest grid, differ from (N1,N2) at the southern and eastern boundary, same as (N1,N2) on the west and north
  real(r8),intent(in) :: RCHQF  !surface water flux scalar for X-boundary flow
  real(r8),intent(in) :: XN     !flow direction: XN=1 when NN=iBehind/2, XN=-1 when NN=iFront/1

  character(len=*), parameter :: subname='XBoundSurfaceRunoffM'
  real(r8) :: ALT1     !elevation at the source grid
  real(r8) :: ALT2     !elevation in the dest center
  real(r8) :: DPTHW1,DPTHW2
  real(r8) :: VX, PondDepz
  !

  call PrintInfo('beg '//subname)
  !
  ! DPTHW1,DPTHW2=surface water depth of source,destination
  ! ALT1,ALT2=elevation of source,destination
  ! XVOLT=excess surface water+ice
  ! VWatStoreCapSurf_col=ground surface water retention capacity
  ! ExtWaterTable_col=natural water table depth, > 0 into soil
  ! XN=direction
  !
  ! RUNOFF
  !

  DPTHW1 = XVLMobileWaterLitR_col(N2,N1)/AREA_3D(3,NUM_col(N2,N1),N2,N1)
  DPTHW2 = VWatStoreCapSurf_col(N2,N1)/AREA_3D(3,NUM_col(N2,N1),N2,N1)

  ALT1     = Altitude_col(N2,N1)+DPTHW1
  !linear extrapolation of surface elevation over the grid
  ALT2     = Altitude_col(N2,N1)+DPTHW2-XN*SLOPE_col(N,N2,N1)*DLYR_3D(N,NUM_col(N2,N1),N2,N1)
  PondDepz = CumDepz2LayBottom_vr(NU_col(N2,N1)-1,N2,N1)-DPTHW1    !can be < 0, if the upper edge of topsoil is at 0 m,

  !depth is counting downward
  !Grid elevation is higher than outside the grid, and in grid water layer higher than external water table
  !grid discharges water to external water table
  IF(ALT1.GT.ALT2 .AND. PondDepz.LT.ExtWaterTable_col(N2,N1))THEN
    !say y=WatFlx2LitRByRunoff_2DH: on the west/north, XN >0, so y <0, +y to update (N2,N1) ;
    !on the east/south, XN < 0, so y > 0, -y to update (N2,N1).

    WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)   = -XN*SurfRunoffPotentM_col(M,N2,N1)*FSLOPE_2DH(N,N2,N1)*RCHQF   !<0 if NN=iBehind, W/N; > 0 if NN=iFront, E/S
    HeatFlx2LitRByRunoff_2DH(N,NN,M5,M4)  = cpw*TKSoil1_vr(0,N2,N1)*WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)
    XGridSurfRunoff_2DH(N,NN,M5,M4)       = XGridSurfRunoff_2DH(N,NN,M5,M4)+WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)
    HeatXGridBySurfRunoff_2DH(N,NN,M5,M4) = HeatXGridBySurfRunoff_2DH(N,NN,M5,M4)+HeatFlx2LitRByRunoff_2DH(N,NN,M5,M4)

    !when NN=2, (M5,M4) = (N2,N1), so plus
    !when NN=1, (M5,M4) differs (N2,N1), so minus
   ! RUNON
   ! grid water table in lower than external water table, so water is recharged from external water table
  ELSEIF(PondDepz.GT.ExtWaterTable_col(N2,N1))THEN
    !elevation difference
    VX                                    = AZMIN1((ExtWaterTable_col(N2,N1)-PondDepz)*AREA_3D(3,NUM_col(N2,N1),N2,N1))  !<0
    SurfRunoffPotentM_col(M,N2,N1)        = VX*dts_wat                                                            !<0
    RunoffVelocityM_col(M,N2,N1)          = 0.0_r8

    WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)   = -XN*SurfRunoffPotentM_col(M,N2,N1)*FSLOPE_2DH(N,N2,N1)*RCHQF          !>0 if NN=iBehind, W/N; < 0 if NN=iFront, E/S
    !should air temperature be used instead of TKSoil1_vr(0,N2,N1)?
    HeatFlx2LitRByRunoff_2DH(N,NN,M5,M4)  = cpw*TKSoil1_vr(0,N2,N1)*WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)
    XGridSurfRunoff_2DH(N,NN,M5,M4)       = XGridSurfRunoff_2DH(N,NN,M5,M4)+WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)
    HeatXGridBySurfRunoff_2DH(N,NN,M5,M4) = HeatXGridBySurfRunoff_2DH(N,NN,M5,M4)+HeatFlx2LitRByRunoff_2DH(N,NN,M5,M4)

    !when NN=iBehind, (M5,M4)=(N2,N1), so plus (>0._r8)
    !when NN=iFront, (M5,M4) differs from (N2,N1), so minus (<0)
  ELSE
    WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)  = 0.0_r8
    HeatFlx2LitRByRunoff_2DH(N,NN,M5,M4) = 0.0_r8
  ENDIF
  QflxSurfRunoffM_2DH(M,N,NN,M5,M4) = WatFlx2LitRByRunoff_2DH(N,NN,M5,M4)
  IFLBM_2DH(M,N,NN,M5,M4)           = 0
  IFLB_2DH(N,NN,M5,M4)              = IFLBM_2DH(M,N,NN,M5,M4)
  call PrintInfo('end '//subname)
  end subroutine XBoundSurfaceRunoffM
!------------------------------------------------------------------------------------------

  subroutine PartitionPrecip(I,J,NY,NX)
  !
  !Do precipitation partitioning over exposed soil, exposed litter and snow.
  !
  implicit none
  integer, intent(in) :: I,J,NY,NX

  character(len=*), parameter :: subname ='PartitionPrecip'
  real(r8) :: SnoFall,PrecHeat2Sno
  real(r8) :: RainThrufall2SoiMicP,RainThrufall2SoiMacP,Rain2Snow
  real(r8) :: RainThrufall2LitR,RainThrufall2Soil,RainHeat2LitR,RainHeat2Soil
  real(r8) :: Rain2ExposedSurf,RainFrac,IrrigFrac
  real(r8), parameter :: m2mm=1.e3_r8

  call PrintInfo('beg '//subname)

  !convert water flux from m/hour to mm/hour
  PRECM_col(NY,NX)             = m2mm*PrecRainAndIrrig_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
  PrecDirect2Grnd_col(NY,NX)   = m2mm*(PrecRainAndIrrig_col(NY,NX)-Rain2Canopy_col(NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
  PrecIndirect2Grnd_col(NY,NX) = m2mm*(Rain2Canopy_col(NY,NX)-RainIntceptByCanopy_col(NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
!
!     RESIDUE WATER ABSORPTION CAPACITY
!
!     HCNDLitr_col=litter saturated hydraulic conductivity
!     DLYRR=litter depth
!
!  HCNDLitr_col(NY,NX) = SatHydroCondLitR
  DLYRR_COL(NY,NX)    = AMAX1(2.5E-03_r8,DLYR_3D(3,0,NY,NX))
!
!     DISTRIBUTION OF PRECIPITATION AND ITS HEAT AMONG SURFACE
!     RESIDUE, SOIL SURFACE, AND MACROPORES
!
!     PRECA,SnoFalPrec=rainfall+irrigation,snowfall (water equiv)
!     Rain2Snow=rainfall to snowpack
!     SnoFall=snowfall to snowpack
!     PrecHeat2Sno=convective heat flux to snowpack
!     Rain2ExposedSurf=precip to litter+soil surfaces
!     RainThrufall2Soil,RainThrufall2LitR=precip to soil,litter surfaces
!     RainHeat2Soil,RainHeat2LitR=convective heat flux to soil,litter surfaces
!     RainThrufall2SoiMicP,RainThrufall2SoiMacP=precip to soil micropores,macropores
!     TFLWC=canopy intercepted precipitation
!     FSNW=fraction of snow cover
!     RainIntceptByCanopy_col=precipitation intercepted by plant canopy

! partition throughfall
  IF(PrecRainAndIrrig_col(NY,NX).GT.0.0_r8 .OR. SnowPrecThrufall_col(NY,NX).GT.0.0_r8)THEN
  ! there is precipitation
    Rain2Snow            = RainPrecThrufall_col(NY,NX)*FracSurfAsSnow_col(NY,NX)
    SnoFall              = SnowPrecThrufall_col(NY,NX) !snowfall
    PrecHeat2Sno         = cps*TKSnowThrufall_col(NY,NX)*SnoFall+cpw*TairK_col(NY,NX)*Rain2Snow                                  !incoming heat flux from precipitations to snow-covered surface
    Rain2ExposedSurf     = RainPrecThrufall_col(NY,NX) *FracSurfSnoFree_col(NY,NX)       !incoming precipitation to snow-free surface
    RainThrufall2LitR    = Rain2ExposedSurf*FracSurfByLitR_col(NY,NX)                             !water flux to snow-free coverd by litter
    RainHeat2LitR        = cpw*TairK_col(NY,NX)*RainThrufall2LitR                                 !heat flux to snow-free surface covered by litter
    RainThrufall2Soil    = Rain2ExposedSurf*FracSurfBareSoil_col(NY,NX)                          !heat flux to snow-free surface not covered by litter
    RainHeat2Soil        = cpw*TairK_col(NY,NX)*RainThrufall2Soil*HeatAdv_scal
    RainThrufall2SoiMicP = RainThrufall2Soil*SoilFracAsMicP_vr(NUM_col(NY,NX),NY,NX)          !water flux to micropore
    RainThrufall2SoiMacP = RainThrufall2Soil*SoilFracAsMacP1_vr(NUM_col(NY,NX),NY,NX)         !water flux to macropore
  ELSE
  ! no rainfall and irrigation, interception should be zero too, thus no throughfall and evertying is zero
    Rain2Snow            = 0._r8      !should be zero
    SnoFall              = 0._r8
    PrecHeat2Sno         = 0._r8
    Rain2ExposedSurf     = 0._r8
    RainThrufall2LitR    = 0._r8
    RainHeat2LitR        = 0._r8
    RainThrufall2Soil    = 0._r8
    RainHeat2Soil        = 0._r8
    RainThrufall2SoiMicP = 0._r8
    RainThrufall2SoiMacP = 0._r8
  ENDIF
  Rain2ExposedSurf_col(NY,NX) = Rain2ExposedSurf
  PrecHeat_col(NY,NX)         = PrecHeat_col(NY,NX)+(PrecHeat2Sno+RainHeat2Soil)*dts_HeatWatTP
!
!     RECIPIENT WATER RATES EXPORTED TO SALT AND NON-SALT DEPOSITION
!
!     Rain2SnowSurf_col includes snowfall water equivalent and liquid rain.
!     Snow/litter/soil exports are rates; deposition applies dts_HeatWatTP once.
!
  ! Export the same recipient split used by the water calculation. Liquid
  ! throughfall carries the rain/irrigation mixture; snowfall goes to snow.
  RainFrac  = 0._r8
  IrrigFrac = 0._r8
  IF(RainFalPrec_col(NY,NX)+IrrigSurface_col(NY,NX)>0._r8)THEN
    RainFrac  = RainFalPrec_col(NY,NX)/(RainFalPrec_col(NY,NX)+IrrigSurface_col(NY,NX))
    IrrigFrac = IrrigSurface_col(NY,NX)/(RainFalPrec_col(NY,NX)+IrrigSurface_col(NY,NX))
  ENDIF
  Rain2SnowSurf_col(NY,NX)  = SnoFall+Rain2Snow*RainFrac
  Irrig2SnowSurf_col(NY,NX) = Rain2Snow*IrrigFrac
  Rain2LitRSurf_col(NY,NX)  = RainThrufall2LitR*RainFrac
  Irrig2LitRSurf_col(NY,NX) = RainThrufall2LitR*IrrigFrac
  Rain2SoilSurf_col(NY,NX)  = RainThrufall2Soil*RainFrac
  Irrig2SoilSurf_col(NY,NX) = RainThrufall2Soil*IrrigFrac
  !
  !     GATHER PRECIPITATION AND MELTWATER FLUXES AND THEIR HEATS
  !     AMONG ATMOSPHERE, SNOWPACK, RESIDUE AND SOIL SURFACES
  !     INTO LOCAL ARRAYS FOR USE IN MASS AND ENERGY EXCHANGE
  !     ALGORITHMS
  !
  !     dts_HeatWatTP=internal time step for fluxes through soil profile
  !
  !     FLW0S,Ice2Snowt_col,Rain2Snowt_col=snow,ice,water input to snowpack
  !     PrecHeat2Snowt_col=convective heat flux to snowpack
  !     Rain2SoiMicP1_col,Rain2SoiMacP1_col,Rain2LitR1_col=rain+irrigation to micropores,macropores,litter
  !     RainHeat2SoilP1_col,HWFLY1=convective heat flux to soil,litter surfaces
  !
  SnowFallt_col(NY,NX)      = SnoFall*dts_HeatWatTP
  Ice2Snowt_col(NY,NX)      = 0.0_r8
  Rain2Snowt_col(NY,NX)     = Rain2Snow*dts_HeatWatTP
  PrecHeat2Snowt_col(NY,NX) = PrecHeat2Sno*dts_HeatWatTP

  Rain2SoiMicP1_col(NY,NX)   = RainThrufall2SoiMicP*dts_HeatWatTP*HeatAdv_scal
  Rain2SoiMacP1_col(NY,NX)   = RainThrufall2SoiMacP*dts_HeatWatTP*HeatAdv_scal
  RainHeat2SoilP1_col(NY,NX) = RainHeat2Soil*dts_HeatWatTP

  Rain2LitR1_col(NY,NX)     = RainThrufall2LitR*dts_HeatWatTP
  RainHeat2LitR1_col(NY,NX) = RainHeat2LitR*dts_HeatWatTP

  call PrintInfo('end '//subname)
  end subroutine PartitionPrecip
!------------------------------------------------------------------------------------------

  subroutine UpdateSurfaceAtM(I,J,M,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: M,I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer  :: NY,NX

! begin_execution
  D9795: DO NX=NHW,NHE
    D9790: DO NY=NVN,NVS

      !update snow state variables
      call UpdateSnowAtM(I,J,M,NY,NX)

      ! Snow cover suppresses atmospheric exchange, not accepted melt/runoff transfers.
      call UpdateLitRAftRunoffM(I,J,M,NY,NX)

    ENDDO D9790
  ENDDO D9795

  end subroutine UpdateSurfaceAtM
!------------------------------------------------------------------------------------------

  subroutine SumAftEnergyBalanceM(I,J,M,NY,NX,LWRadGrnd_col,VapXAir2TopLay,HeatSensLitR2Soi1,&
    HeatSensVapLitR2Soi1,EvapLitR2Soi1,VapXAir2LitR,HeatFluxAir2LitR,HeatFluxAir2Soi,PrecNet2SoiMicP,&
    PrecNet2SoiMacP,RainPrecHeatAir2LitR,NetWatFlxAir2SoiMicP,NetWatXFlxAir2SoiMicP,NetWatFlxAir2SoiMacP)
  !
  !Sum up fluxes after doing surface energy balance calculation
  implicit none
  integer, intent(in)  :: I,J,M
  integer, intent(in)  :: NY,NX
  real(r8), intent(in) :: LWRadGrnd_col            !>0 into atmosphere, long wave radiation from ground surface
  real(r8), intent(in) :: VapXAir2TopLay
  real(r8), intent(in) :: HeatSensLitR2Soi1
  real(r8), intent(in) :: HeatSensVapLitR2Soi1
  real(r8), intent(in) :: EvapLitR2Soi1
  real(r8), intent(in) :: HeatFluxAir2LitR
  real(r8), intent(in) :: VapXAir2LitR
  real(r8), intent(in) :: HeatFluxAir2Soi
  real(r8), intent(in) :: PrecNet2SoiMicP
  real(r8), intent(in) :: PrecNet2SoiMacP
  real(r8), intent(in) :: RainPrecHeatAir2LitR
  real(r8), intent(out):: NetWatFlxAir2SoiMicP
  real(r8), intent(out):: NetWatXFlxAir2SoiMicP
  real(r8), intent(out):: NetWatFlxAir2SoiMacP

  character(len=*), parameter :: subname='SumAftEnergyBalanceM'
  real(r8) :: VLSnowHeatCap,VLSnowHeatCap0,TKSX,ENGY
  real(r8) :: TMX,deltaT

  ! begin_execution
  call PrintInfo('beg '//subname)
  !
  ! GATHER WATER, VAPOR AND HEAT FLUXES INTO FLUX ARRAYS
  ! FOR LATER UPDATES TO STATE VARIABLES
  !
  ! NetWatFlxAir2SoiMicP,NetWatFlxAir2SoiMacP=water flux from atm to soil micropores,macropores
  ! cumHeatSensAir2Soil=convective heat flux from atm to soil
  ! NetWatFlxAir2LitR=water flux from atm to litter
  ! CumHeatSensAir2LitR=convective heat flux from atm to litter
  !
  NetWatFlxAir2SoiMicP  = PrecNet2SoiMicP+VapXAir2TopLay+EvapLitR2Soi1
  NetWatXFlxAir2SoiMicP = PrecNet2SoiMicP+VapXAir2TopLay+EvapLitR2Soi1
  NetWatFlxAir2SoiMacP  = PrecNet2SoiMacP
  cumHeatSensAir2Soil   = PrecHeat2SoiNet+HeatFluxAir2Soi+HeatSensVapLitR2Soi1+HeatSensLitR2Soi1
  !write(*,*) "PrecHeat2SoiNet: ", PrecHeat2SoiNet
  !write(*,*) "HeatFluxAir2Soi: ", HeatFluxAir2Soi
  !write(*,*) "HeatSensVapLitR2Soi1: ", HeatSensVapLitR2Soi1
  !write(*,*) "HeatSensLitR2Soi1: ", HeatSensLitR2Soi1

  !total water and heat fluxes to litter
  NetWatFlxAir2LitR   = RainPrecAir2LitR+VapXAir2LitR-EvapLitR2Soi1
  CumHeatSensAir2LitR = RainPrecHeatAir2LitR+HeatFluxAir2LitR-HeatSensVapLitR2Soi1-HeatSensLitR2Soi1

  !
  ! GENERATE NEW SNOWPACK
  !
  ! XFLWS,WatXfer2SnoLay,IceXfer2SnoLay=hourly snow,water,ice transfer
  ! SnowFallt_col,Rain2Snowt_col,Ice2Snowt_col=snow,water,ice input to snowpack
  ! HeatXfer2SnoLay=hourly convective heat flux from snow,water,ice transfer
  ! PrecHeat2Snowt_col=convective heat flux from snow,water,ice to snowpack
  !

  IF(VLSnowHeatCapM_snvr(M,1,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX) .AND. SnowFallt_col(NY,NX)+Rain2Snowt_col(NY,NX).GT.ZEROS(NY,NX))THEN

    SnoXfer2SnoLay_snvr(1,NY,NX)  = SnoXfer2SnoLay_snvr(1,NY,NX)+SnowFallt_col(NY,NX)
    WatXfer2SnoLay_snvr(1,NY,NX)  = WatXfer2SnoLay_snvr(1,NY,NX)+Rain2Snowt_col(NY,NX)
    IceXfer2SnoLay_snvr(1,NY,NX)  = IceXfer2SnoLay_snvr(1,NY,NX)+Ice2Snowt_col(NY,NX)
    HeatXfer2SnoLay_snvr(1,NY,NX) = HeatXfer2SnoLay_snvr(1,NY,NX)+PrecHeat2Snowt_col(NY,NX)
    Prec2Snow_col(NY,NX)          = Prec2Snow_col(NY,NX) + SnowFallt_col(NY,NX)+Rain2Snowt_col(NY,NX)+Ice2Snowt_col(NY,NX)
    PrecHeat2Snow_col(NY,NX)      = PrecHeat2Snow_col(NY,NX)+PrecHeat2Snowt_col(NY,NX)
    RainPrec2Sno_col(NY,NX)       = RainPrec2Sno_col(NY,NX)+Rain2Snowt_col(NY,NX)

    !accumulate new snallfall
    VLSnowHeatCap0 = cps*VLDrySnoWE0_snvr(1,NY,NX)+cpw*VLWatSnow0_snvr(1,NY,NX)+cpi*VLIceSnow0_snvr(1,NY,NX)

    VLDrySnoWE0_snvr(1,NY,NX)      = VLDrySnoWE0_snvr(1,NY,NX)+SnowFallt_col(NY,NX)

    VLIceSnow0_snvr(1,NY,NX)       = VLIceSnow0_snvr(1,NY,NX)+Ice2Snowt_col(NY,NX)
    VLWatSnow0_snvr(1,NY,NX)       = VLWatSnow0_snvr(1,NY,NX)+Rain2Snowt_col(NY,NX)
    VLSnowHeatCap                  = cps*VLDrySnoWE0_snvr(1,NY,NX)+cpw*VLWatSnow0_snvr(1,NY,NX)+cpi*VLIceSnow0_snvr(1,NY,NX)
    VLSnowHeatCapM_snvr(M,1,NY,NX) = VLSnowHeatCap
    TKSX                           = TKSnow0_snvr(1,NY,NX)
    ENGY                           = VLSnowHeatCap0*TKSX
    TKSnow0_snvr(1,NY,NX)          = (ENGY+PrecHeat2Snowt_col(NY,NX))/VLSnowHeatCap

  ENDIF

  !LWRadBySurf_col=longwave emission from litter and surface soil into atmosphere
  LWRadBySurf_col(NY,NX) = LWRadBySurf_col(NY,NX)+LWRadGrnd_col

  VapXAir2GSurf_col(NY,NX)  = VapXAir2GSurf_col(NY,NX) + VapXAir2LitR
  TEvapXAir2LitR_col(NY,NX) = TevapXAir2LitR_col(NY,NX)+ VapXAir2LitR
  call PrintInfo('end '//subname)
  end subroutine SumAftEnergyBalanceM
!------------------------------------------------------------------------------------------
  subroutine CorrectSoilSurfaceHeatM(NY,NX,TKPrev,TKNew,dHeatdTK,NetHeatCorrection)
  ! Record the boundary exchange accepted by the implicit topsoil solve.
  implicit none
  integer, intent(in) :: NY,NX             !horizontal grid indices
  real(r8), intent(in) :: TKPrev,TKNew     !soil temperature before/after solve [K]
  real(r8), intent(in) :: dHeatdTK         !derivative of step-integrated net boundary heat [MJ/K]
  real(r8), intent(out) :: NetHeatCorrection !accepted change to soil net heat [MJ/grid cell]
  real(r8) :: DeltaT,RadHeatCorrection,SensHeatCorrection
  integer :: L

  NetHeatCorrection=0._r8
  if(dHeatdTK==0._r8)return
  L=NUM_col(NY,NX)
  DeltaT=TKNew-TKPrev
  NetHeatCorrection=dHeatdTK*DeltaT
  RadHeatCorrection=-4._r8*LWEmscefSoil_col(NY,NX)*TKPrev**3*DeltaT
  SensHeatCorrection=NetHeatCorrection-RadHeatCorrection

  ! Correct both the accepted surface face and its cell/hourly accumulators.
  HeatFlow2Soili_3D(3,L,NY,NX) = HeatFlow2Soili_3D(3,L,NY,NX)+NetHeatCorrection
  HeatFlow2Soil_3D(3,L,NY,NX)  = HeatFlow2Soil_3D(3,L,NY,NX)+NetHeatCorrection
  HeatFlx2Grnd_col(NY,NX)      = HeatFlx2Grnd_col(NY,NX)+NetHeatCorrection
  HeatByRad2Surf_col(NY,NX)    = HeatByRad2Surf_col(NY,NX)+RadHeatCorrection
  HeatSensAir2Surf_col(NY,NX)  = HeatSensAir2Surf_col(NY,NX)+SensHeatCorrection
  HeatNet2Surf_col(NY,NX)      = HeatNet2Surf_col(NY,NX)+NetHeatCorrection
  LWRadBySurf_col(NY,NX)       = LWRadBySurf_col(NY,NX)-RadHeatCorrection
  LWRadGrnd_col(NY,NX)         = LWRadGrnd_col(NY,NX)-RadHeatCorrection
  end subroutine CorrectSoilSurfaceHeatM

!------------------------------------------------------------------------------------------
  subroutine RunSurfacePhysModelM(I,J,M,NHE,NHW,NVS,NVN,ResistanceLitRLay,RainEkReducedKsat,&
    TopLayWatVol_col,HeatFluxAir2Soi,dHGrnddTKsoi1,Qinfl2MicP,HeatInfl2Soil,Qinfl2MacP)
  !
  !run surface energy/water model for iteration M
  implicit none
  integer, intent(in) :: I,J !day, hour
  integer, intent(in) :: M   !soil heat-water iteration id
  integer, intent(in) :: NHE,NHW,NVS,NVN
  real(r8), dimension(:,:),intent(inout) :: ResistanceLitRLay(JY,JX)
  REAL(R8), dimension(:,:),INTENT(OUT) :: RainEkReducedKsat
  real(r8), dimension(:,:),intent(inout) :: TopLayWatVol_col(JY,JX)
  
  ! Atmospheric heat predictor for soil conduction [MJ/grid cell/step];
  ! the soil temperature solve subsequently accepts its implicit correction.
  real(r8), dimension(:,:),intent(out) :: HeatFluxAir2Soi(JY,JX)

  ! Derivative of net soil radiation+sensible heat integrated over this step
  ! [MJ/K], evaluated at the current topsoil temperature; normally negative.
  ! UpdateSoilMoistTempM accepts its flux correction after solving temperature.
  real(r8), dimension(:,:),intent(out) :: dHGrnddTKsoi1(JY,JX)
  real(r8), dimension(:,:),optional,intent(out) :: Qinfl2MacP(JY,JX)       !water into soil macropores [m3/grid cell/step]
  real(r8), dimension(:,:),optional,intent(out) :: Qinfl2MicP(JY,JX)       !water into soil micropores [m3/grid cell/step]
  real(r8), dimension(:,:),optional,intent(out) :: HeatInfl2Soil(JY,JX)    !boundary heat before soil implicit correction [MJ/grid cell]

  character(len=*), parameter :: subname='RunSurfacePhysModelM'
  real(r8) :: LatentHeatAir2Sno,HeatSensAir2Snow,Radnet2Snow,HeatSensEvapAir2Snow,VapXAir2TopLay
  integer  :: N1,N2,NX,NY,L
  real(r8) :: QWatinfl2Mic, QHeatInfl2Soil

  call PrintInfo('beg '//subname)
  WatFlx2LitRByRunoff_2DH(:,:,:,:)  = 0._r8
  HeatFlx2LitRByRunoff_2DH(:,:,:,:) = 0._r8

  ! Finish every cell's local exchanges before reading neighboring water
  ! levels or temperatures. No runoff is applied during this first pass.
  D9895: DO  NX=NHW,NHE
    D9890: DO  NY=NVN,NVS

      ! No exposed-soil sensible derivative survives a snow-covered step.
      CdSoiHSens_col(NY,NX)=0._r8
      call SurfaceEnergyModelM(I,J,M,NX,NY,ResistanceLitRLay(NY,NX),RainEkReducedKsat(NY,NX),&
        HeatFluxAir2Soi(NY,NX),LatentHeatAir2Sno,HeatSensEvapAir2Snow,HeatSensAir2Snow,Radnet2Snow,&
        TopLayWatVol_col(NY,NX),VapXAir2TopLay)

      call SurfLitrSoilWaterExchangeM(I,J,M,NY,NX,RainEkReducedKsat(NY,NX))

      call UpdateLitterPhaseChangeM(I,J,M,NY,NX)

      call AccumWaterVaporHeatFluxesM(I,J,M,NY,NX,LatentHeatAir2Sno,HeatSensEvapAir2Snow,HeatSensAir2Snow,&
        Radnet2Snow,VapXAir2TopLay)

      ! Commit snowmelt and litter-soil exchange even under complete snow cover.
      call UpdateLitRBe4RunoffM(I,J,M,NY,NX)

      if(present(Qinfl2MicP))Qinfl2MicP(NY,NX)=WaterFlow2Micpt_3D(3,NUM_col(NY,NX),NY,NX)
      if(present(Qinfl2MacP))Qinfl2MacP(NY,NX)=WaterFlow2Macpt_3D(3,NUM_col(NY,NX),NY,NX)
      if(present(HeatInfl2Soil))HeatInfl2Soil(NY,NX)=HeatFlow2Soili_3D(3,NUM_col(NY,NX),NY,NX)
      LWRadGrnd_col(NY,NX)   = LWRadBySurf_col(NY,NX)
      
      dHGrnddTKsoi1(NY,NX)=0._r8
      ! Match the exposed-soil exchange gate in AtmLandSurfExchangeM.
      if(FracSurfSnoFree_col(NY,NX)>0._r8 .and. &
        (SoilBulkDensity_vr(NUM_col(NY,NX),NY,NX)>ZERO .or. &
        VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX)>VHCPNX_col(NY,NX)))then
        dHGrnddTKsoi1(NY,NX)=-4._r8*LWEmscefSoil_col(NY,NX) &
          *TKSoil1_vr(NUM_col(NY,NX),NY,NX)**3-CdSoiHSens_col(NY,NX)
      endif
    ENDDO D9890
  ENDDO D9895

  ! Clear this substep's snow faces, including unused boundary slots, before
  ! constructing accepted transfers from the common post-local donor state.
  IF(snowRedist_model)THEN
    DrySnoFlxByRedistM_2DH(M,:,:,:,:)=0._r8
    WatSnoFlxByRedistM_2DH(M,:,:,:,:)=0._r8
    IceSnoFlxByRedistM_2DH(M,:,:,:,:)=0._r8
    HeatSnoFlxByRedistM_2DH(M,:,:,:,:)=0._r8
    IFLBSM_2DH(M,:,:,:,:)=-1
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        call SnowRedistributionM(M,NY,NX,NHE,NHW,NVS,NVN)
      ENDDO
    ENDDO
  ENDIF

  ! All donors and receivers now share the same post-local-process state.
  ! Potential and face construction only read that state; the donor limiter
  ! and simultaneous runoff update are performed later in WatsubMod.
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      call CalcSurfaceRunoffPotentialM(I,J,M,NY,NX)
      if(.not.ATS_cpl_mode)call XGridsSurfRunoffM(I,J,M,NY,NX,NHE,NHW,NVS,NVN)
    ENDDO
  ENDDO

  call PrintInfo('end '//subname)
  end subroutine RunSurfacePhysModelM

!------------------------------------------------------------------------------------------
  subroutine SurfaceEnergyModelM(I,J,M,NX,NY,ResistanceLitRLay,RainEkReducedKsat,&
    HeatFluxAir2Soi1,LatentHeatAir2Sno,HeatSensEvapAir2Snow,HeatSensAir2Snow,Radnet2Snow,&
    TopLayWatVol,VapXAir2TopLay)
  !
  !Description
  !call surface energy balance model
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M     !soil heat-flow iteration
  integer, intent(in) :: NX,NY
  real(r8),intent(inout) :: ResistanceLitRLay
  REAL(R8), INTENT(OUT) :: RainEkReducedKsat,HeatFluxAir2Soi1
  real(r8), intent(out) :: Radnet2Snow,LatentHeatAir2Sno,HeatSensAir2Snow,HeatSensEvapAir2Snow
  real(r8), intent(inout) :: TopLayWatVol
  real(r8), intent(out) :: VapXAir2TopLay

  character(len=*), parameter :: subname = 'SurfaceEnergyModelM'
  integer :: N1,N2,L
  real(r8):: PrecNet2SoiMicP,PrecNet2SoiMacP,RainPrecHeatAir2LitR
  real(r8) :: Prec2LitR2,RainHeat2LitR2,SnowRadHeatCorrection

  call PrintInfo('beg '//subname)

  !ResistanceLitRLay is input
  call InitSurfModelM(I,J,M,NY,NX,ResistanceLitRLay,RainEkReducedKsat,Prec2LitR2,&
    RainHeat2LitR2,PrecNet2SoiMicP,PrecNet2SoiMacP,RainPrecHeatAir2LitR)

  ! updates ResistanceLitRLay
  call AtmLandSurfExchangeM(I,J,M,NY,NX,Prec2LitR2,RainHeat2LitR2,PrecNet2SoiMicP,&
    PrecNet2SoiMacP,RainPrecHeatAir2LitR,ResistanceLitRLay,TopLayWatVol,LatentHeatAir2Sno,&
    HeatSensEvapAir2Snow,HeatSensAir2Snow,Radnet2Snow,VapXAir2TopLay,HeatFluxAir2Soi1)

  !update snow pack before doing snow redistribution to avoid negative mass values
  call UpdateSnowPack1M(I,J,M,NY,NX,SnowRadHeatCorrection)
  ! Include the outer snow stabilizer before surface heat is accumulated.
  Radnet2Snow=Radnet2Snow+SnowRadHeatCorrection

  call PrintInfo('end '//subname)
  end subroutine SurfaceEnergyModelM

!------------------------------------------------------------------------------------------
  pure function GaucklerManningVelocity(HydraulicRadius, slope)result(CrossSectVelocity)
  !
  !Description
  !
  !ref: https://en.wikipedia.org/wiki/Manning_formula
  !V=k/n*Rh^(2/3)*S^(1/2)
  implicit none
  real(r8), intent(in) :: HydraulicRadius
  real(r8), intent(in) :: slope   !stream slope, or hydraulic gradient, or channel bed slope

  real(r8) :: CrossSectVelocity    !m/s

  CrossSectVelocity=HydraulicRadius**0.67_r8*SQRT(SLOPE)

  end function GaucklerManningVelocity

!------------------------------------------------------------------------------------------

  subroutine writeSurfDiagnosis(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,j
  integer, intent(in) :: NY,NX

  write(106,*)I+J/24.,TEvapXAir2Toplay_col(NY,NX)*1000._r8/AREA_3D(3,NU_col(NY,NX),NY,NX), &
    TEvapXAir2LitR_col(NY,NX)*1000._r8/AREA_3D(3,NU_col(NY,NX),NY,NX), &
    TEvapXAir2Snow_col(NY,NX)*1000._r8/AREA_3D(3,NU_col(NY,NX),NY,NX)

  end subroutine writeSurfDiagnosis

!------------------------------------------------------------------------------------------

  subroutine AggregateSurfRunoffFluxM(I,J,M,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J,M,NHW,NHE,NVN,NVS
  character(len=*), parameter :: subname='AggregateSurfRunoffFluxM'
  integer :: NX,NY,N,Side
  integer :: FaceNN(2,2),FaceNY(2,2),FaceNX(2,2)
  real(r8) :: OutwardSign(2,2)

  call PrintInfo('beg '//subname)

  !TODO: consider adaptive substeps for converging flows into depressions.
  !The donor limit preserves water/heat accounting and liquid positivity, but
  !does not prevent collective water-level overshoot from multiple donors.
  !Finalize every donor's outgoing faces before any receiving cell is credited.
  call LimitSurfaceRunoffM(M,NHW,NHE,NVN,NVS)

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      cumWatFlx2LitRByRunoff_col(NY,NX)=0._r8
      cumHeatFlx2LitRByRunoff_col(NY,NX)=0._r8
      call GetSurfaceRunoffFaces(NX,NY,NHW,NHE,NVN,NVS,FaceNN,FaceNY,FaceNX,OutwardSign)

      DO N=1,2
        !Internal inflows are stored in this cell's receiving slots.
        !External inflows use the boundary face and are included below instead.
        if((N.EQ.iWestEastDirection .and. NX.LT.NHE) .or. &
           (N.EQ.iNorthSouthDirection .and. NY.LT.NVS))then
          cumWatFlx2LitRByRunoff_col(NY,NX)=cumWatFlx2LitRByRunoff_col(NY,NX) &
            +WatFlx2LitRByRunoff_2DH(N,iFront,NY,NX)
          cumHeatFlx2LitRByRunoff_col(NY,NX)=cumHeatFlx2LitRByRunoff_col(NY,NX) &
            +HeatFlx2LitRByRunoff_2DH(N,iFront,NY,NX)
        endif
        if((N.EQ.iWestEastDirection .and. NX.GT.NHW) .or. &
           (N.EQ.iNorthSouthDirection .and. NY.GT.NVN))then
          cumWatFlx2LitRByRunoff_col(NY,NX)=cumWatFlx2LitRByRunoff_col(NY,NX) &
            +WatFlx2LitRByRunoff_2DH(N,iBehind,NY,NX)
          cumHeatFlx2LitRByRunoff_col(NY,NX)=cumHeatFlx2LitRByRunoff_col(NY,NX) &
            +HeatFlx2LitRByRunoff_2DH(N,iBehind,NY,NX)
        endif

        DO Side=iFront,iBehind
          cumWatFlx2LitRByRunoff_col(NY,NX)=cumWatFlx2LitRByRunoff_col(NY,NX) &
            -OutwardSign(N,Side)*WatFlx2LitRByRunoff_2DH(N,FaceNN(N,Side),FaceNY(N,Side),FaceNX(N,Side))
          cumHeatFlx2LitRByRunoff_col(NY,NX)=cumHeatFlx2LitRByRunoff_col(NY,NX) &
            -OutwardSign(N,Side)*HeatFlx2LitRByRunoff_2DH(N,FaceNN(N,Side),FaceNY(N,Side),FaceNX(N,Side))
        ENDDO
      ENDDO

      TXGridSurfRunoff_2DH(NY,NX)=TXGridSurfRunoff_2DH(NY,NX)+cumWatFlx2LitRByRunoff_col(NY,NX)
      THeatXGridBySurfRunoff_2DH(NY,NX)=THeatXGridBySurfRunoff_2DH(NY,NX)+cumHeatFlx2LitRByRunoff_col(NY,NX)
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine AggregateSurfRunoffFluxM

!------------------------------------------------------------------------------------------
  subroutine LimitSurfaceRunoffM(M,NHW,NHE,NVN,NVS)
  !Project all outgoing runoff faces onto each donor's beginning-of-pass liquid.
  !Do not credit concurrent inflow: limiting must be independent of cell order.
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  integer :: NX,NY,N,Side,NN,FY,FX
  integer :: FaceNN(2,2),FaceNY(2,2),FaceNX(2,2)
  real(r8) :: OutwardSign(2,2),ProposedOutflow,AvailableWater,RunoffScale
  real(r8) :: OldWater,OldHeat,AcceptedWater,AcceptedHeat

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      call GetSurfaceRunoffFaces(NX,NY,NHW,NHE,NVN,NVS,FaceNN,FaceNY,FaceNX,OutwardSign)
      ProposedOutflow=0._r8
      DO N=1,2
        DO Side=iFront,iBehind
          ProposedOutflow=ProposedOutflow+MAX(0._r8,OutwardSign(N,Side) &
            *WatFlx2LitRByRunoff_2DH(N,FaceNN(N,Side),FaceNY(N,Side),FaceNX(N,Side)))
        ENDDO
      ENDDO

      AvailableWater=MAX(0._r8,VLWatMicP1_vr(0,NY,NX)-tiny_wat)
      if(ProposedOutflow.LE.AvailableWater)cycle
      RunoffScale=AvailableWater/ProposedOutflow
      DO N=1,2
        DO Side=iFront,iBehind
          NN=FaceNN(N,Side);FY=FaceNY(N,Side);FX=FaceNX(N,Side)
          OldWater=WatFlx2LitRByRunoff_2DH(N,NN,FY,FX)
          if(OutwardSign(N,Side)*OldWater.LE.0._r8)cycle
          OldHeat=HeatFlx2LitRByRunoff_2DH(N,NN,FY,FX)
          AcceptedWater=OldWater*RunoffScale
          AcceptedHeat=OldHeat*RunoffScale
          WatFlx2LitRByRunoff_2DH(N,NN,FY,FX)=AcceptedWater
          HeatFlx2LitRByRunoff_2DH(N,NN,FY,FX)=AcceptedHeat
          QflxSurfRunoffM_2DH(M,N,NN,FY,FX)=AcceptedWater
          !Only replace this substep's contribution; retain earlier substeps.
          XGridSurfRunoff_2DH(N,NN,FY,FX)=XGridSurfRunoff_2DH(N,NN,FY,FX)+AcceptedWater-OldWater
          HeatXGridBySurfRunoff_2DH(N,NN,FY,FX)=HeatXGridBySurfRunoff_2DH(N,NN,FY,FX)+AcceptedHeat-OldHeat
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  end subroutine LimitSurfaceRunoffM

!------------------------------------------------------------------------------------------
  subroutine GetSurfaceRunoffFaces(NX,NY,NHW,NHE,NVN,NVS,FaceNN,FaceNY,FaceNX,OutwardSign)
  !Return this donor's E/S (iFront) and W/N (iBehind) face slots.
  !Internal transfers are positive in receiving slots. Boundary W/N outflows
  !are negative in the donor slot; boundary E/S outflows are positive outside.
  implicit none
  integer, intent(in) :: NX,NY,NHW,NHE,NVN,NVS
  integer, intent(out) :: FaceNN(2,2),FaceNY(2,2),FaceNX(2,2)
  real(r8), intent(out) :: OutwardSign(2,2)
  integer :: N,Position,FirstCell,LastCell

  DO N=1,2
    FaceNY(N,:)=NY;FaceNX(N,:)=NX
    if(N.EQ.iWestEastDirection)then
      Position=NX;FirstCell=NHW;LastCell=NHE
      FaceNX(N,iFront)=NX+1
      FaceNX(N,iBehind)=NX-1
    else
      Position=NY;FirstCell=NVN;LastCell=NVS
      FaceNY(N,iFront)=NY+1
      FaceNY(N,iBehind)=NY-1
    endif
    FaceNN(N,iFront)=iBehind
    FaceNN(N,iBehind)=iFront
    OutwardSign(N,:)=1._r8
    if(Position.EQ.LastCell)FaceNN(N,iFront)=iFront
    if(Position.EQ.FirstCell)then
      FaceNN(N,iBehind)=iBehind
      FaceNY(N,iBehind)=NY;FaceNX(N,iBehind)=NX
      OutwardSign(N,iBehind)=-1._r8
    endif
  ENDDO
  end subroutine GetSurfaceRunoffFaces
end module SurfPhysMod
