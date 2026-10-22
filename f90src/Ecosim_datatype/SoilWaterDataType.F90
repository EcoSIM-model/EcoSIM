module SoilWaterDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8),target,allocatable ::  TXGridSurfRunoff_2DH(:,:)                 !water flux into the grid due to runoff, [m3 d-2 h-1]
  real(r8),target,allocatable ::  THeatXGridBySurfRunoff_2DH(:,:)           !heat flux into the grid due to runoff,[MJ d-2 h-1]
  logical, target,allocatable ::  iPondFlag_col(:,:)                        !flag for the col if it is a pond,[-]
  integer, target,allocatable ::  iPondBotLev_col(:,:)                      !Bottom level ID, [-]
  real(r8),target,allocatable ::  DVLiceMicP_vr(:,:,:)                      !change of water+ice volume, [m3 d-2]
  real(r8),target,allocatable ::  ThetaAir_vr(:,:,:)                        !air concentration ,[m3 m-3]
  real(r8),target,allocatable ::  VLsoiAirP_vr(:,:,:)                       !soil air content ,[m3 d-2]
  real(r8),target,allocatable ::  THETW_vr(:,:,:)                           !volumetric water content ,[m3 m-3]
  real(r8),target,allocatable ::  THETI_vr(:,:,:)                           !volumetric ice content ,[m3 m-3]
  real(r8),target,allocatable ::  ThetaH2OZ_vr(:,:,:)                       !relative water saturation in soil [m3 water/m3 pore]; litter (L=0): surface depth [m]
  real(r8),target,allocatable ::  ThetaICEZ_vr(:,:,:)                       !relative ice saturation in soil [m3 ice/m3 pore]; litter (L=0): surface depth [m]
  real(r8),target,allocatable ::  VLWatMicP_vr(:,:,:)                       !soil micropore water content ,[m3 d-2]
  real(r8),target,allocatable ::  VLiceMicP_vr(:,:,:)                       !soil micropore ice content   ,[m3 d-2]
  real(r8),target,allocatable ::  VLWatMacP_vr(:,:,:)                       !soil macropore water content ,[m3 d-2]
  real(r8),target,allocatable ::  PSISoilMatricP_vr(:,:,:)                  !soil micropore matric water potential ,[MPa]
  real(r8),target,allocatable ::  ElvAdjstedSoilH2OPSIMPa_vr(:,:,:)         !elevation adjusted total soil micropore total water potential ,[MPa]
  real(r8),target,allocatable ::  VLWatMicPX_vr(:,:,:)                      !soil micropore water content before wetting front ,[m3 d-2]
  real(r8),target,allocatable ::  FWatExMacP2MicP_vr(:,:,:)                 !soil macropore - micropore water transfer ,[m3 d-2 h-1]
  real(r8),target,allocatable ::  VLiceMacP_vr(:,:,:)                       !soil macropore ice content ,[m3 d-2]
  real(r8),target,allocatable ::  VLWatMicPM_vr(:,:,:,:)                    !soil micropore water content, [m3 d-2]
  real(r8),target,allocatable ::  VLWatMacPM_vr(:,:,:,:)                    !soil macropore water content, [m3 d-2]
  real(r8),target,allocatable ::  VLsoiAirPM_vr(:,:,:,:)                    !soil air content, [m3 d-2]
  real(r8),target,allocatable ::  FILMM_vr(:,:,:,:)                         !soil water film thickness , [m]
  real(r8),target,allocatable ::  WaterTable2GroundSlope_col(:,:)                    !slope of water table relative to surface slope, [-]
  real(r8),target,allocatable ::  TileWaterTableDepz_col(:,:)                     !depth of artificial water table
  real(r8),target,allocatable ::  TileWaterTable_col(:,:)                   !artificial water table depth, [m]
  real(r8),target,allocatable ::  TileWaterTableLocDepz_col(:,:)                                !depth of artificial water table adjusted for elevation,[m]
  real(r8),target,allocatable ::  WaterTableDepZi_col(:,:)                      !internal water table depth, [m]
  real(r8),target,allocatable ::  NaturalWaterTableDepzTopoC_col(:,:)       !topography (slope and elevation) corrected local water table depth,[m]
  real(r8),target,allocatable ::  NaturalWaterTableDepzLitrC_col(:,:)       !litter-layer corrected (on top of topography) water table depth (>0 lower than soil surface), [m]
  real(r8),target,allocatable ::  NaturalWaterTableDepz_col(:,:)            !landscape external water table depth, [m]
  real(r8),target,allocatable ::  EnergyImpact4Erosion_colM(:,:,:)          !total energy impact for erosion
  real(r8),target,allocatable ::  XVLMobileWaterLitRM(:,:,:)                !excess water+ice, [m3 d-2]
  real(r8),target,allocatable ::  XVLMobileWatMicPM(:,:,:)                  !excess water,[m3 d-2]
  real(r8),target,allocatable ::  XVLiceMicPM(:,:,:)                        !excess ice, [m3 d-2]
  real(r8),target,allocatable ::  HydroCond_3D(:,:,:,:,:)                   !hydraulic conductivity at different moisture levels, [m MPa-1 h-1]
  real(r8),target,allocatable ::  HydroCondMacP_vr(:,:,:)                   !macropore hydraulic conductivity, [m MPa-1 h-1]
  real(r8),target,allocatable ::  HYCDMicP4RootUptake_vr(:,:,:)        !soil micropore hydraulic conductivity for root water uptake ,[m MPa-1 h-1]
  real(r8),target,allocatable ::  SurfRunoffPotentM_col(:,:,:)              !runoff water flux out of grid (>=0), [m3 d-2 t-1]
  real(r8),target,allocatable ::  RunoffVelocityM_col(:,:,:)                !runoff velocity, [m t-1], 1 W/N, 0 E/S, -1 ,[-]
  integer,target,allocatable ::  IFLBM_2DH(:,:,:,:,:)                       !flag for directional surface runoff,[-]
  logical,target,allocatable ::  XGridRunoffFlag_2DH(:,:,:,:)               !enables or disables boundary water flux depending on aspect, [-]
  integer,target,allocatable ::  IFLB_2DH(:,:,:,:)                           !flag for directional runoff, related to IFLBM_2DH
  real(r8),target,allocatable ::  RechrgDistNorthSubSurf_col(:,:)             !sclar for northern subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RechrgDistEastSubSurf_col(:,:)                 !sclar for eastern subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RechrgDistSouthSubSurf_col(:,:)                !sclar for southern subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RechrgDistWestSubSurf_col(:,:)                 !sclar for western subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargRateNorthWTBL_col(:,:)                 !northern subsurface boundary water flux rate constant, [h-1]
  real(r8),target,allocatable ::  RechargRateEastWTBL_col(:,:)                  !eastern subsurface boundary water flux  rate constant, [h-1]
  real(r8),target,allocatable ::  RechargRateSouthWTBL_col(:,:)                 !southern subsurface boundary water flux  rate constant, [h-1]
  real(r8),target,allocatable ::  RechargRateWestWTBL_col(:,:)                  !western subsurface boundary water flux  rate constant, [h-1]
  real(r8),target,allocatable ::  RechargNorthSurf_col(:,:)                     !sclar for northern surface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargEastSurf_col(:,:)                      !sclar for eastern surface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargSouthSurf_col(:,:)                     !sclar for southern surface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargWestSurf_col(:,:)                  !sclar for western surface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargBottom_col(:,:)                    !sclar for lower subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  WaterFlow2MicPM_3D(:,:,:,:,:)             !micropore water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  WaterFlow2MacPM_3D(:,:,:,:,:)             !macropore water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  ReductVLsoiAirPM_vr(:,:,:,:)              !change in soil air volume for layer from last to current iteration, >0, shrink,[g d-2 h-1]
  real(r8),target,allocatable ::  FWatExMacP2MicPM_vr(:,:,:,:)              !soil macropore - micropore water transfer, [g d-2 t-1]
  real(r8),target,allocatable ::  WatFlowSno2MicPM_col(:,:,:)               !meltwater flux into soil micropores, [m3 d-2 h-1]
  real(r8),target,allocatable ::  WatFlowSno2MacPM_col(:,:,:)               !meltwater flux into soil macropores,[m3 d-2 h-1]
  real(r8),target,allocatable ::  FracAirFilledSoilPoreM_vr(:,:,:,:)        !air-filled soil porosity, [m3 m-3]
  real(r8),target,allocatable ::  TortMicPM_vr(:,:,:,:)                     !soil micropore tortuosity, [m3 m-3]
  real(r8),target,allocatable ::  TortMacPM_vr(:,:,:,:)                     !soil macropore tortuosity, [m3 m-3]
  real(r8),target,allocatable ::  DiffusivitySolutEffM_vr(:,:,:,:)          !coefficient for dissolution - volatilization, []
  real(r8),target,allocatable ::  SoilBulkModulus4RootPent_vr(:,:,:)         !elastic modulus of the undisturbed soil, [MPa]
  real(r8),target,allocatable ::  SoilModulus4RootRadialexp_vr(:,:,:)       !soil modulus for root radial expansion, [MPa]    
  real(r8),target,allocatable ::  SoilWeightStress_vr(:,:,:)                  !soil weight stress on root thickening, [MPa]
  real(r8),target,allocatable ::  PSI_at_Sautration_vr(:,:,:)                           !soil water potential at saturation, [Mpa]
  real(r8),target,allocatable ::  PSISoilAirEntry_vr(:,:,:)                    !soil water potential at air entry, [Mpa]
  real(r8),target,allocatable ::  PSISoilOsmotic_vr(:,:,:)                  !osmotic soil water potential , [Mpa]
  real(r8),target,allocatable ::  PSIGrav_vr(:,:,:)                         !gravimetric soil water potential , [Mpa]
  real(r8),target,allocatable ::  SoilWatAirDry_vr(:,:,:)                   !air-dry water content, [m3 m-3]
  real(r8),target,allocatable ::  ThetaSat_vr(:,:,:)                        !micropore class water content
  real(r8),target,allocatable ::  WaterFlowSoiMicPX_3D(:,:,:,:)             !unsaturated water flux , [m3 d-2 h-1]
  real(r8),target,allocatable ::  EvapoTransp_col(:,:)                      !evapotranspiration,[m3 d-2 h-1]
  real(r8),target,allocatable ::  QEvap_CumYr_col(:,:)                      !total evaporation, [m3 d-2]
  real(r8),target,allocatable ::  QRain_CumYr_col(:,:)                      !total precipitation, [m3 d-2]
  real(r8),target,allocatable ::  Qrunoff_CumYr_col(:,:)                    !total surface runoff, [m3 d-2]
  real(r8),target,allocatable ::  WatMass_col(:,:)                          !total soil water content, [m3 d-2]
  real(r8),target,allocatable ::  H2OLoss_CumYr_col(:,:)                    !total subsurface water flux, [m3 d-2]
  real(r8),target,allocatable ::  QDrainloss_vr(:,:,:)                      !water drainage from different soil layers, [m3 d-2]
  real(r8),target,allocatable ::  QDrain_col(:,:)                           !total water drainage from soil, [m3 d-2]
  real(r8),target,allocatable ::  QDrain_cum_col(:,:)                       !cumulative total water drainage from soil, [m3 d-2]
  real(r8),target,allocatable ::  XGridSurfRunoff_2DH(:,:,:,:)              !soil surface runoff water, [m3 d-2 h-1]
  real(r8),target,allocatable ::  HeatXGridBySurfRunoff_2DH(:,:,:,:)        !soil surface runoff heat, [MJ d-2 h-1]
  real(r8),target,allocatable ::  QRunSurf_col(:,:)                         !runoff from surface water, [m3 d-2 h-1]
  real(r8),target,allocatable ::  QDischarg2WTBL_col(:,:)                   !water discharge to external water table, [m3 d-2 h-1]
  real(r8),target,allocatable ::  QflxSurfRunoffM_2DH(:,:,:,:,:)            !surface runoff in iteration M,
  real(r8),target,allocatable ::  Qinflx2Soil_col(:,:)                      !infiltration into soil [m3 d-2 h-1]
  real(r8),target,allocatable ::  SoilWatMassBeg_col(:,:)                   !soil water mass at the begnining of time step
  real(r8),target,allocatable ::  SoilWatMassEnd_col(:,:)                   !soil water mass at the end of time step
  real(r8),target,allocatable ::  Rain2Soil_col(:,:)                        !water flow into soil due to precipitation (+ surface irrigation), [m3 H2O/d2/h]
  real(r8),target,allocatable ::  QdewCanopy_CumYr_pft(:,:,:)               !cumulative dew deposition on canopy, [m3 d-2]
  real(r8),target,allocatable ::  QSnoWatXfer2Soil_col(:,:)                 !snow water flux to soil, [m3 d-2 h-1]
  real(r8),target,allocatable ::  QSnoIceXfer2Soil_col(:,:)                 !snow ice flux to soil, [m3 d-2 h-1]
  real(r8),target,allocatable ::  PrecipAtm2LandSurf_col(:,:)                !precipiation from atmosphere to land surface ,[m3 d-2 h-1]
  real(r8),target,allocatable ::  RainPrecThrufall_col(:,:)                  !rainfall precipitation through canopy, [m3 H2O d-2 h-1]
  real(r8),target,allocatable ::  SnowPrecThrufall_col(:,:)                  !snowfall through canopy, [m3 H2O d-2 h-1]
  real(r8),target,allocatable ::  TKSnowThrufall_col(:,:)                    !temperature of canopy throughfall snow, [K]
  real(r8),target,allocatable ::  RainPrec2Sno_col(:,:)                      !rainfall to snow, [m3 H2O d-2 h-1]
  real(r8),target,allocatable ::  Rain2ExposedSurf_col(:,:)                  !rainfall to exposed surface, [m3 H2O d-2 h-1]
  real(r8),target,allocatable ::  QLaterFlow2Cell_col(:,:)                  !Internal lateral flow between grids, [m3 H2O d-2 h-1]
  real(r8),target,allocatable ::  HydCondSoil_3D(:,:,:,:)                    !3D micropore hydraulic conductivity, [m MPa-1 h-1]
  real(r8),target,allocatable ::  TWatFlowCellMicP_vr(:,:,:)                 !water flow into cell through micropores, [m3 H2O d-2 h-1]
  real(r8),target,allocatable ::  TWatFlowCellMacP_vr(:,:,:)                 !water flow into cell through macropores, [m3 H2O d-2 h-1]  
  private :: InitAllocate
  contains

  subroutine InitSoilWater
  implicit none

  call InitAllocate
  end subroutine InitSoilWater

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none

  allocate(DVLiceMicP_vr(JZ,JY,JX)); DVLiceMicP_vr=0._r8
  allocate(TWatFlowCellMicP_vr(JZ,JY,JX));     TWatFlowCellMicP_vr=0._r8
  allocate(TWatFlowCellMacP_vr(JZ,JY,JX));    TWatFlowCellMacP_vr=0._r8
  allocate(iPondBotLev_col(JY,JX)); iPondBotLev_col=0
  allocate(iPondFlag_col(JY,JX)); iPondFlag_col =.false.
  allocate(QLaterFlow2Cell_col(JY,JX)); QLaterFlow2Cell_col=0._r8
  allocate(Rain2ExposedSurf_col(JY,JX)); Rain2ExposedSurf_col=0._r8
  allocate(SnowPrecThrufall_col(JY,JX)); SnowPrecThrufall_col=0._r8
  allocate(TKSnowThrufall_col(JY,JX)); TKSnowThrufall_col=0._r8
  allocate(RainPrec2Sno_col(JY,JX)); RainPrec2Sno_col = 0._r8
  allocate(RainPrecThrufall_col(JY,JX)); RainPrecThrufall_col=0._r8
  allocate(PrecipAtm2LandSurf_col(JY,JX)); PrecipAtm2LandSurf_col=0._r8
  allocate(QSnoIceXfer2Soil_col(JY,JX)); QSnoIceXfer2Soil_col=0._r8
  allocate(QSnoWatXfer2Soil_col(JY,JX)); QSnoWatXfer2Soil_col=0._r8
  allocate(Rain2Soil_col(JY,JX));  Rain2Soil_col=0._r8
  allocate(SoilWatMassBeg_col(JY,JX)); SoilWatMassBeg_col=0._r8
  allocate(SoilWatMassEnd_col(JY,JX)); SoilWatMassEnd_col=0._r8
  allocate(QdewCanopy_CumYr_pft(JZ,JY,JX)); QdewCanopy_CumYr_pft=0._r8
  allocate(EvapoTransp_col(JY,JX)); EvapoTransp_col=0._r8
  allocate(Qinflx2Soil_col(JY,JX)); Qinflx2Soil_col=0._r8
  allocate(ThetaAir_vr(0:JZ,JY,JX));  ThetaAir_vr=0._r8
  allocate(VLsoiAirP_vr(0:JZ,JY,JX));   VLsoiAirP_vr=0._r8
  allocate(THETW_vr(0:JZ,JY,JX));  THETW_vr=0._r8
  allocate(THETI_vr(0:JZ,JY,JX));  THETI_vr=0._r8
  allocate(ThetaH2OZ_vr(0:JZ,JY,JX)); ThetaH2OZ_vr=0._r8
  allocate(ThetaICEZ_vr(0:JZ,JY,JX)); ThetaICEZ_vr=0._r8
  allocate(VLWatMicP_vr(0:JZ,JY,JX));   VLWatMicP_vr=0._r8
  allocate(VLiceMicP_vr(0:JZ,JY,JX));   VLiceMicP_vr=0._r8
  allocate(VLWatMacP_vr(JZ,JY,JX));    VLWatMacP_vr=0._r8
  allocate(PSISoilMatricP_vr(0:JZ,JY,JX));  PSISoilMatricP_vr=0._r8
  allocate(ElvAdjstedSoilH2OPSIMPa_vr(0:JZ,JY,JX));  ElvAdjstedSoilH2OPSIMPa_vr=0._r8
  allocate(VLWatMicPX_vr(0:JZ,JY,JX));  VLWatMicPX_vr=0._r8
  allocate(FWatExMacP2MicP_vr(JZ,JY,JX));     FWatExMacP2MicP_vr=0._r8
  allocate(VLiceMacP_vr(JZ,JY,JX));    VLiceMacP_vr=0._r8
  allocate(VLWatMicPM_vr(60,0:JZ,JY,JX));VLWatMicPM_vr=0._r8
  allocate(VLWatMacPM_vr(60,JZ,JY,JX));VLWatMacPM_vr=0._r8
  allocate(VLsoiAirPM_vr(60,0:JZ,JY,JX));VLsoiAirPM_vr=0._r8
  allocate(FILMM_vr(60,0:JZ,JY,JX));FILMM_vr=0._r8
  allocate(WaterTable2GroundSlope_col(JY,JX));       WaterTable2GroundSlope_col=0._r8
  allocate(TileWaterTableDepz_col(JY,JX));      TileWaterTableDepz_col=0._r8
  allocate(TileWaterTable_col(JY,JX));       TileWaterTable_col=0._r8
  allocate(TileWaterTableLocDepz_col(JY,JX));       TileWaterTableLocDepz_col=0._r8
  allocate(WaterTableDepZi_col(JY,JX));       WaterTableDepZi_col=0._r8
  allocate(NaturalWaterTableDepzTopoC_col(JY,JX));       NaturalWaterTableDepzTopoC_col=0._r8
  allocate(NaturalWaterTableDepzLitrC_col(JY,JX));       NaturalWaterTableDepzLitrC_col=0._r8
  allocate(NaturalWaterTableDepz_col(JY,JX));       NaturalWaterTableDepz_col=0._r8
  allocate(EnergyImpact4Erosion_colM(60,JY,JX));   EnergyImpact4Erosion_colM=0._r8
  allocate(XVLMobileWaterLitRM(60,JY,JX));   XVLMobileWaterLitRM=0._r8
  allocate(XVLMobileWatMicPM(60,JY,JX));   XVLMobileWatMicPM=0._r8
  allocate(XVLiceMicPM(60,JY,JX));   XVLiceMicPM=0._r8
  allocate(HydroCond_3D(3,100,0:JZ,JY,JX));HydroCond_3D=0._r8
  allocate(HydroCondMacP_vr(JZ,JY,JX));     HydroCondMacP_vr=0._r8
  allocate(HYCDMicP4RootUptake_vr(JZ,JY,JX));     HYCDMicP4RootUptake_vr=0._r8
  allocate(SurfRunoffPotentM_col(60,JV,JH));      SurfRunoffPotentM_col=0._r8
  allocate(RunoffVelocityM_col(60,JY,JX));      RunoffVelocityM_col=0._r8
  allocate(IFLBM_2DH(60,2,2,JY,JX));IFLBM_2DH=-1
  allocate(XGridRunoffFlag_2DH(2,2,JY,JX));   XGridRunoffFlag_2DH=.false.
  allocate(IFLB_2DH(2,2,JY,JX));   IFLB_2DH=0
  allocate(RechrgDistNorthSubSurf_col(JY,JX));      RechrgDistNorthSubSurf_col=0._r8
  allocate(RechrgDistEastSubSurf_col(JY,JX));      RechrgDistEastSubSurf_col=0._r8
  allocate(RechrgDistSouthSubSurf_col(JY,JX));      RechrgDistSouthSubSurf_col=0._r8
  allocate(RechrgDistWestSubSurf_col(JY,JX));      RechrgDistWestSubSurf_col=0._r8
  allocate(RechargRateNorthWTBL_col(JY,JX));      RechargRateNorthWTBL_col=0._r8
  allocate(RechargRateEastWTBL_col(JY,JX));      RechargRateEastWTBL_col=0._r8
  allocate(RechargRateSouthWTBL_col(JY,JX));      RechargRateSouthWTBL_col=0._r8
  allocate(RechargRateWestWTBL_col(JY,JX));      RechargRateWestWTBL_col=0._r8
  allocate(RechargNorthSurf_col(JY,JX));       RechargNorthSurf_col=0._r8
  allocate(RechargEastSurf_col(JY,JX));       RechargEastSurf_col=0._r8
  allocate(RechargSouthSurf_col(JY,JX));       RechargSouthSurf_col=0._r8
  allocate(RechargWestSurf_col(JY,JX));       RechargWestSurf_col=0._r8
  allocate(RechargBottom_col(JY,JX));       RechargBottom_col=0._r8
  allocate(WaterFlow2MicPM_3D(60,3,JD,JV,JH));WaterFlow2MicPM_3D=0._r8
  allocate(WaterFlow2MacPM_3D(60,3,JD,JV,JH));WaterFlow2MacPM_3D=0._r8
  allocate(ReductVLsoiAirPM_vr(60,JZ,JY,JX));  ReductVLsoiAirPM_vr=0._r8
  allocate(FWatExMacP2MicPM_vr(60,JZ,JY,JX)); FWatExMacP2MicPM_vr=0._r8
  allocate(WatFlowSno2MicPM_col(60,JY,JX));    WatFlowSno2MicPM_col=0._r8
  allocate(WatFlowSno2MacPM_col(60,JY,JX));    WatFlowSno2MacPM_col=0._r8
  allocate(FracAirFilledSoilPoreM_vr(60,0:JZ,JY,JX));FracAirFilledSoilPoreM_vr=0._r8
  allocate(TortMicPM_vr(60,0:JZ,JY,JX));TortMicPM_vr=0._r8
  allocate(TortMacPM_vr(60,JZ,JY,JX)); TortMacPM_vr=0._r8
  allocate(DiffusivitySolutEffM_vr(60,0:JZ,JY,JX));DiffusivitySolutEffM_vr=0._r8
  allocate(SoilBulkModulus4RootPent_vr(JZ,JY,JX));     SoilBulkModulus4RootPent_vr=0._r8
  allocate(SoilModulus4RootRadialexp_vr(JZ,JY,JX)); SoilModulus4RootRadialexp_vr=0._r8
  allocate(SoilWeightStress_vr(JZ,JY,JX)); SoilWeightStress_vr=0._r8
  allocate(PSI_at_Sautration_vr(0:JZ,JY,JX));  PSI_at_Sautration_vr=0._r8
  allocate(PSISoilAirEntry_vr(0:JZ,JY,JX));  PSISoilAirEntry_vr=0._r8
  allocate(PSISoilOsmotic_vr(0:JZ,JY,JX));  PSISoilOsmotic_vr=0._r8
  allocate(PSIGrav_vr(0:JZ,JY,JX));  PSIGrav_vr=0._r8
  allocate(SoilWatAirDry_vr(0:JZ,JY,JX));  SoilWatAirDry_vr=0._r8
  allocate(ThetaSat_vr(0:JZ,JY,JX));  ThetaSat_vr=0._r8
  allocate(WaterFlowSoiMicPX_3D(3,JD,JV,JH));   WaterFlowSoiMicPX_3D=0._r8
  allocate(QEvap_CumYr_col(JY,JX));       QEvap_CumYr_col=0._r8
  allocate(QRain_CumYr_col(JY,JX));       QRain_CumYr_col=0._r8
  allocate(Qrunoff_CumYr_col(JY,JX));        Qrunoff_CumYr_col=0._r8
  allocate(WatMass_col(JY,JX));       WatMass_col=0._r8
  allocate(H2OLoss_CumYr_col(JY,JX));       H2OLoss_CumYr_col=0._r8
  allocate(QDrain_col(JY,JX));      QDrain_col=0._r8
  allocate(QDrainloss_vr(JZ,JY,JX));      QDrainloss_vr=0._r8
  allocate(XGridSurfRunoff_2DH(2,2,JV,JH));      XGridSurfRunoff_2DH=0._r8
  allocate(HeatXGridBySurfRunoff_2DH(2,2,JV,JH));     HeatXGridBySurfRunoff_2DH=0._r8
  allocate(QRunSurf_col(JY,JX));        QRunSurf_col=0._r8
  allocate(QDischarg2WTBL_col(JY,JX));       QDischarg2WTBL_col=0._r8
  allocate(QflxSurfRunoffM_2DH(60,2,2,JV,JH)); QflxSurfRunoffM_2DH=0._r8
  allocate(TXGridSurfRunoff_2DH(JY,JX));         TXGridSurfRunoff_2DH=0._r8
  allocate(THeatXGridBySurfRunoff_2DH(JY,JX));        THeatXGridBySurfRunoff_2DH=0._r8
  allocate(QDrain_cum_col(JY,JX));      QDrain_cum_col=0._r8  
  allocate(HydCondSoil_3D(3,JZ,JY,JX)); HydCondSoil_3D=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSoilWater
  use abortutils, only : destroy
  implicit none

  CALL destroy(DVLiceMicP_vr)
  call destroy(iPondBotLev_col)
  call destroy(iPondFlag_col)
  call destroy(QLaterFlow2Cell_col)
  call destroy(Rain2ExposedSurf_col)
  call destroy(SnowPrecThrufall_col)
  call destroy(TKSnowThrufall_col)
  call destroy(RainPrec2Sno_col)
  call destroy(RainPrecThrufall_col)
  call destroy(PrecipAtm2LandSurf_col)
  call destroy(QSnoIceXfer2Soil_col)
  call destroy(QSnoWatXfer2Soil_col)
  call destroy(Rain2Soil_col)
  call destroy(SoilWatMassBeg_col)
  call destroy(SoilWatMassEnd_col)
  call destroy(QdewCanopy_CumYr_pft)
  call destroy(ThetaAir_vr)
  call destroy(VLsoiAirP_vr)
  call destroy(THETW_vr)
  call destroy(THETI_vr)
  call destroy(ThetaH2OZ_vr)
  call destroy(ThetaICEZ_vr)
  call destroy(VLWatMicP_vr)
  call destroy(VLiceMicP_vr)
  call destroy(VLWatMacP_vr)
  call destroy(PSISoilMatricP_vr)
  call destroy(ElvAdjstedSoilH2OPSIMPa_vr)
  call destroy(VLWatMicPX_vr)
  call destroy(FWatExMacP2MicP_vr)
  call destroy(VLiceMacP_vr)
  call destroy(VLWatMicPM_vr)
  call destroy(VLWatMacPM_vr)
  call destroy(VLsoiAirPM_vr)
  call destroy(FILMM_vr)
  call destroy(WaterTable2GroundSlope_col)
  call destroy(TileWaterTableDepz_col)
  call destroy(TileWaterTable_col)
  call destroy(TileWaterTableLocDepz_col)
  call destroy(WaterTableDepZi_col)
  call destroy(NaturalWaterTableDepzTopoC_col)
  call destroy(NaturalWaterTableDepzLitrC_col)
  call destroy(NaturalWaterTableDepz_col)
  call destroy(EnergyImpact4Erosion_colM)
  call destroy(XVLMobileWaterLitRM)
  call destroy(XVLMobileWatMicPM)
  call destroy(XVLiceMicPM)
  call destroy(HydroCond_3D)
  call destroy(HydroCondMacP_vr)
  call destroy(HYCDMicP4RootUptake_vr)
  call destroy(SurfRunoffPotentM_col)
  call destroy(RunoffVelocityM_col)
  call destroy(IFLBM_2DH)
  call destroy(XGridRunoffFlag_2DH)
  call destroy(IFLB_2DH)
  call destroy(TWatFlowCellMicP_vr)
  call destroy(TWatFlowCellMacP_vr)  
  call destroy(RechrgDistNorthSubSurf_col)
  call destroy(RechrgDistEastSubSurf_col)
  call destroy(RechrgDistSouthSubSurf_col)
  call destroy(RechrgDistWestSubSurf_col)
  call destroy(RechargRateNorthWTBL_col)
  call destroy(RechargRateEastWTBL_col)
  call destroy(RechargRateSouthWTBL_col)
  call destroy(RechargRateWestWTBL_col)
  call destroy(RechargNorthSurf_col)
  call destroy(RechargEastSurf_col)
  call destroy(RechargSouthSurf_col)
  call destroy(RechargWestSurf_col)
  call destroy(RechargBottom_col)
  call destroy(WaterFlow2MicPM_3D)
  call destroy(WaterFlow2MacPM_3D)
  call destroy(ReductVLsoiAirPM_vr)
  call destroy(FWatExMacP2MicPM_vr)
  call destroy(WatFlowSno2MicPM_col)
  call destroy(WatFlowSno2MacPM_col)
  call destroy(FracAirFilledSoilPoreM_vr)
  call destroy(TortMicPM_vr)
  call destroy(TortMacPM_vr)
  call destroy(DiffusivitySolutEffM_vr)
  call destroy(SoilBulkModulus4RootPent_vr)
  call destroy(SoilModulus4RootRadialexp_vr)
  call destroy(SoilWeightStress_vr)
  call destroy(PSI_at_Sautration_vr)
  call destroy(PSISoilAirEntry_vr)
  call destroy(PSISoilOsmotic_vr)
  call destroy(PSIGrav_vr)
  call destroy(SoilWatAirDry_vr)
  call destroy(ThetaSat_vr)
  call destroy(WaterFlowSoiMicPX_3D)
  call destroy(QEvap_CumYr_col)
  call destroy(EvapoTransp_col)
  call destroy(QRain_CumYr_col)
  call destroy(Qrunoff_CumYr_col)
  call destroy(WatMass_col)
  call destroy(H2OLoss_CumYr_col)
  call destroy(QDrain_col)
  call destroy(QDrainloss_vr)  
  call destroy(QDrain_cum_col)  
  call destroy(XGridSurfRunoff_2DH)
  call destroy(HeatXGridBySurfRunoff_2DH)
  call destroy(QRunSurf_col)
  call destroy(QDischarg2WTBL_col)
  call destroy(QflxSurfRunoffM_2DH)
  call destroy(Qinflx2Soil_col)
  call destroy(TXGridSurfRunoff_2DH)
  call destroy(THeatXGridBySurfRunoff_2DH)
  call destroy(HydCondSoil_3D)
  end subroutine DestructSoilWater

end module SoilWaterDataType
