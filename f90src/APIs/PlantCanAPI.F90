module PlantCanAPI

! interface to integrate the plant model
  use data_kind_mod,              only: r8 => DAT_KIND_R8
  use EcoSiMParDataMod,           only: micpar, pltpar
  use SoilPhysDataType,           only: SurfAlbedo_col
  use MiniMathMod,                only: AZMAX1
  use EcoSIMCtrlMod,              only: ldo_sp_mode
  use PlantSiteAPIData,           only: plt_site
  use PlantPhotosynthesisAPIData, only: plt_photo
  use PlantRadiationAPIData,      only: plt_rad
  use PlantMorphologyAPIData,     only: plt_morph
  use PlantPhenologyAPIData,      only: plt_pheno
  use PlantSoilChemistryAPIData,  only: plt_soilchem
  use PlantBiomassAPIData,        only: plt_biom
  use PlantEnergyWaterAPIData,    only: plt_ew
  use EcoSIMSolverPar
  use EcoSIMHistMod
  use SnowDataType
  use TracerIDMod
  use SurfLitterDataType
  use LandSurfDataType
  use SoilPropertyDataType
  use ChemTranspDataType
  use EcoSimSumDataType
  use SoilHeatDataType
  use SOMDataType
  use ClimForcDataType
  use EcoSIMCtrlDataType
  use GridDataType
  use RootDataType
  use SoilWaterDataType
  use CanopyDataType
  use PlantDataRateType
  use PlantTraitDataType
  use CanopyRadDataType
  use FlagDataType
  use EcosimBGCFluxType
  use FertilizerDataType
  use SoilBGCDataType
  use PlantMgmtDataType
  use PlantAPICommonData

implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: PlantAPICanMSend
  public :: PlantAPICanMRecv

  contains
!------------------------------------------------------------------------------------------

  subroutine PlantAPICanMSend(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,N,M,NN,NZ,K,NB
  associate(CanopySnow => plt_ew%SnowOnCanopy_pft, &
    StandingDead => plt_biom%StandDeadStrutElms_pft, &
    GrowthThreshold => plt_biom%ZERO4Groth_pft)
!  Integers
  plt_site%KoppenClimZone        = KoppenClimZone_col(NY,NX)
  plt_site%NP                    = NP_col(NY,NX)
  plt_site%NU                     =NU_col(NY,NX)
  plt_morph%StemArea_col         = StemArea_col(NY,NX)
  plt_morph%CanopyLeafArea_col   = CanopyLeafArea_col(NY,NX)
  plt_site%ZEROS                 = ZEROS(NY,NX)
  plt_site%ZERO                  = ZERO
  plt_ew%SnowDepth               = SnowDepth_col(NY,NX)
  plt_ew%TairK                   = TairK_col(NY,NX)
  plt_morph%CanopyHeight_col     = CanopyHeight_col(NY,NX)
  plt_site%WindMesureHeight_col  = WindMesureHeight_col(NY,NX)
  plt_ew%ZeroPlaneDisplacem_col  = ZeroPlaneDisplacem_col(NY,NX)
  plt_site%WindSpeedAtm_col      = WindSpeedAtm_col(NY,NX)
  plt_ew%VLHeatCapSnowMin_col    = VLHeatCapSnowMin_col(NY,NX)
  plt_ew%VLHeatCapSurfSnow_col   = VLHeatCapSnow_snvr(1,NY,NX)
  plt_morph%CanopyHeightZ_col(0) = CanopyHeightZ_col(0,NY,NX)

  if(ldo_sp_mode)then
    DO NZ=1,NP_col(NY,NX)
      plt_morph%tlai_day_pft(NZ) = tlai_day_pft(NZ,NY,NX)
      plt_morph%tsai_day_pft(NZ) = tsai_day_pft(NZ,NY,NX)

      DO  N=1,NumLeafInclinationClasses
        plt_morph%LeafAngleClass_pft(N,NZ)=LeafAngleClass_pft(N,NZ,NY,NX)
      ENDDO
    ENDDO      
  endif
  DO L=1,NumCanopyLayers
    plt_morph%CanopyStemAareZ_col(L) = CanopyStemAareZ_col(L,NY,NX)
    plt_morph%CanopyLeafAareZ_col(L) = CanopyLeafAareZ_col(L,NY,NX)
    plt_morph%CanopyHeightZ_col(L)   = CanopyHeightZ_col(L,NY,NX)
    plt_rad%TAU_DirectSunLit(L)      = TAU_DirectSunLit(L,NY,NX)
  ENDDO
  DO L=1,NumCanopyLayers+1
    plt_rad%TAU_DirectSunLit(L) = TAU_DirectSunLit(L,NY,NX)
    plt_rad%TAU_DirectSunSha(L) = TAU_DirectSunSha(L,NY,NX)
  ENDDO  

  DO L=0,NL_col(NY,NX)
    plt_site%AREA3(L)                 = AREA_3D(3,L,NY,NX)
    plt_soilchem%VLSoilPoreMicP_vr(L) = VLSoilPoreMicP_vr(L,NY,NX)
    plt_soilchem%VLSoilMicP_vr(L)     = VLSoilMicP_vr(L,NY,NX)
    plt_soilchem%VLWatMicP_vr(L)      = VLWatMicP_vr(L,NY,NX)
  ENDDO

  plt_ew%SnowOnCanopy_col          = SnowOnCanopy_col(NY,NX)
  plt_ew%fSnowCanopy_col           = fSnowCanopy_col(NY,NX)
  plt_rad%SineSunInclinationAngle_col     = SineSunInclinationAngle_col(NY,NX)
  plt_site%SolarNoonHour_col       = SolarNoonHour_col(NY,NX)
  plt_morph%LeafStalkAreaAll_col      = LeafStalkAreaAll_col(NY,NX)
  plt_rad%GroundSurfaceAzimuth_col    = GroundSurfaceAzimuth_col(NY,NX)
  plt_rad%CosineGrndSlope_col      = CosineGrndSlope_col(NY,NX)
  plt_rad%SineGrndSlope_col        = SineGrndSlope_col(NY,NX)
  plt_rad%RadSWDiffus_col          = RadSWDiffus_col(NY,NX)
  plt_rad%RadPARDiffus_col         = RadPARDiffus_col(NY,NX)
  plt_rad%RadSWDirect_col          = RadSWDirect_col(NY,NX)
  plt_rad%RadDirectPAR_col         = RadDirectPAR_col(NY,NX)
  plt_rad%RadPARGrnd_col           = RadPARGrnd_col(NY,NX)
  plt_site%SoilSurfRoughness_col = SoilSurfRoughness_col(NY,NX)
  plt_ew%VcumWatSnow_col           = VcumWatSnow_col(NY,NX)
  plt_ew%VcumIceSnow_col           = VcumIceSnow_col(NY,NX)
  plt_ew%VcumDrySnoWE_col          = VcumDrySnoWE_col(NY,NX)
  plt_rad%TotSineSkyAngles_grd     = TotSineSkyAngles_grd
  plt_rad%SoilAlbedo               = SoilAlbedo_col(NY,NX)
  plt_rad%SurfAlbedo_col           = SurfAlbedo_col(NY,NX)
  plt_site%ZEROS2                  = ZEROS2(NY,NX)
  plt_ew%BulkFactor4Snow_col       = BulkFactor4Snow_col(NY,NX)
  plt_site%POROS1                  = POROS_vr(NU_col(NY,NX),NY,NX)
  DO NZ=1,NP_col(NY,NX)
    ! Radiation precedes PlantAPISend, including the first hour after restart.
    ! Supply these inputs here instead of relying on the preceding plant call.
    CanopySnow(NZ) = SnowOnCanopy_pft(NZ,NY,NX)
    StandingDead(:,NZ) = StandDeadStrutElms_pft(:,NZ,NY,NX)
    GrowthThreshold(NZ) = ZERO4Groth_pft(NZ,NY,NX)
    plt_site%PlantPopuDead_pft(NZ)          = PlantPopuDead_pft(NZ,NY,NX)
    ! Radiation runs before PlantModel refreshes the shared plant API state.
    ! Refill biomass, snow and leaf-angle inputs for this column and restart.
    plt_biom%StandDeadStrutElms_pft(1:NumPlantChemElms,NZ) = StandDeadStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_photo%LeafAreaSunlit_pft(NZ)       = LeafAreaSunlit_pft(NZ,NY,NX)
    plt_biom%StomatalStress_pft(NZ)        = StomatalStress_pft(NZ,NY,NX)
    plt_biom%ZERO4Groth_pft(NZ)            = ZERO4Groth_pft(NZ,NY,NX)
    plt_ew%SnowOnCanopy_pft(NZ)            = SnowOnCanopy_pft(NZ,NY,NX)
    plt_morph%LeafAngleClass_pft(:,NZ)     = LeafAngleClass_pft(:,NZ,NY,NX)
    plt_morph%StemSpecVolume_pft(NZ)        = StemSpecVolume_pft(NZ,NY,NX)
    plt_ew%SnowOnCanopy_pft(NZ)             = SnowOnCanopy_pft(NZ,NY,NX)
    plt_ew%fSnowCanopy_pft(NZ)              = fSnowCanopy_pft(NZ,NY,NX)
    plt_ew%BulkFactor4Snow_pft(NZ)          = BulkFactor4Snow_pft(NZ,NY,NX)
    plt_pheno%iPlant2ndGrothPattern_pft(NZ) = iPlant2ndGrothPattern_pft(NZ,NY,NX)
    plt_morph%CanopyLeafArea_pft(NZ)        = CanopyLeafArea_pft(NZ,NY,NX)
    plt_morph%CanopyHeightDead_pft(NZ)      = CanopyHeightDead_pft(NZ,NY,NX)
    plt_morph%CanopyHeightLive_pft(NZ)      = CanopyHeightLive_pft(NZ,NY,NX)
    plt_morph%ClumpFactorNow_pft(NZ)        = ClumpFactorNow_pft(NZ,NY,NX)
    plt_rad%LeafSWabsorptivity_pft(NZ)      = LeafSWabsorptivity_pft(NZ,NY,NX)
    plt_rad%LeafPARabsorptivity_pft(NZ)     = LeafPARabsorptivity_pft(NZ,NY,NX)
    plt_rad%RadSWLeafTransmitance_pft(NZ)   = RadSWLeafTransmitance_pft(NZ,NY,NX)
    plt_rad%RadSWLeafAlbedo_pft(NZ)         = RadSWLeafAlbedo_pft(NZ,NY,NX)
    plt_rad%RadPARLeafTransmitance_pft(NZ)  = RadPARLeafTransmitance_pft(NZ,NY,NX)
    plt_rad%RadPARLeafAlbedo_pft(NZ)        = RadPARLeafAlbedo_pft(NZ,NY,NX)
    plt_morph%NumOfBranches_pft(NZ)         = NumOfBranches_pft(NZ,NY,NX)
    plt_morph%ClumpFactor_pft(NZ)           = ClumpFactor_pft(NZ,NY,NX)

    DO  L=1,NumCanopyLayers
      plt_morph%CanopyStemSurfAreaZ_pft(L,NZ) = CanopyStemSurfAreaZ_pft(L,NZ,NY,NX)
      plt_morph%CanopyLeafAreaZ_pft(L,NZ) = CanopyLeafAreaZ_pft(L,NZ,NY,NX)
      plt_morph%CanopySurfAreaProfDead_pft(L,NZ) = CanopySurfAreaProfDead_pft(L,NZ,NY,NX)
    ENDDO

    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO K=0,MaxNodesPerBranch
        DO  L=1,NumCanopyLayers
          plt_morph%CanopyLeafArea_lnode(L,K,NB,NZ)=CanopyLeafArea_lnode(L,K,NB,NZ,NY,NX)
        ENDDO
      ENDDO
      DO  L=1,NumCanopyLayers
        plt_morph%CanopyStalkSurfArea_lbrch(L,NB,NZ)=CanopyStalkSurfArea_lbrch(L,NB,NZ,NY,NX)        
      ENDDO
      DO K=1,MaxNodesPerBranch
        DO  L=1,NumCanopyLayers
          DO N=1,NumLeafInclinationClasses
            plt_morph%LeafAreaZsec_brch(N,L,K,NB,NZ)=LeafAreaZsec_brch(N,L,K,NB,NZ,NY,NX)
          ENDDO
        ENDDO
      ENDDO
      DO  L=1,NumCanopyLayers
        DO N=1,NumLeafInclinationClasses
          plt_morph%StemAreaZsec_brch(N,L,NB,NZ)=StemAreaZsec_brch(N,L,NB,NZ,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO N=1,NumOfSkyAzimuthSects
    plt_rad%OMEGA2Ground(N)=OMEGA2Ground(N,NY,NX)
  ENDDO
  DO N=1,NumLeafInclinationClasses
    plt_rad%CosineLeafAngle(N) = CosineLeafAngle(N)
    plt_rad%SineLeafAngle(N)   = SineLeafAngle(N)
  ENDDO
  DO NN=1,NumOfLeafAzimuthSectors
    DO M=1,NumLeafInclinationClasses
      DO N=1,NumOfSkyAzimuthSects
        plt_rad%OMEGA2Leaf(N,M,NN)        = OMEGA2Leaf(N,M,NN)
        plt_rad%OMEGX(N,M,NN)             = OMEGX(N,M,NN)
        plt_rad%iScatteringDiffus(N,M,NN) = iScatteringDiffus(N,M,NN)
      ENDDO
    ENDDO
  ENDDO

  end associate
  end subroutine PlantAPICanMSend

!------------------------------------------------------------------------------------------

  subroutine PlantAPICanMRecv(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  integer :: N,M,NN,L,NZ,K,NB

  ZeroPlaneDisplacem_col(NY,NX)  = plt_ew%ZeroPlaneDisplacem_col
  RoughnessLength_col(NY,NX)         = plt_ew%RoughnessLength
  RawIsoTAtm2CanopySinkZ_col(NY,NX) = plt_ew%RawIsoTAtm2CanopySinkZ_col
  RawIsoTSurf2CanopyHScal_col(NY,NX)= plt_ew%RawIsoTSurf2CanopyHScal_col
  RawCanopyH2SinkZ_col(NY,NX)       = plt_ew%RawCanopyH2SinkZ_col
  RIB_col(NY,NX)                     = plt_ew%RIB
  
  CanopyHeight_col(NY,NX)    = plt_morph%CanopyHeight_col
  RadSWDirect_col(NY,NX)     = plt_rad%RadSWDirect_col
  RadSWDiffus_col(NY,NX)     = plt_rad%RadSWDiffus_col
  RadDirectPAR_col(NY,NX)    = plt_rad%RadDirectPAR_col
  RadPARDiffus_col(NY,NX)    = plt_rad%RadPARDiffus_col
  RadSWGrnd_col(NY,NX)       = plt_rad%RadSWGrnd_col
  RadPARGrnd_col(NY,NX)      = plt_rad%RadPARGrnd_col
  FracSWRad2Grnd_col(NY,NX)  = plt_rad%FracSWRad2Grnd_col
  RadSWSolarBeam_col(NY,NX)  = plt_rad%RadSWSolarBeam_col
  RadPARSolarBeam_col(NY,NX) = plt_rad%RadPARSolarBeam_col
  
  DO L=0,NumCanopyLayers
    CanopyHeightZ_col(L,NY,NX)=plt_morph%CanopyHeightZ_col(L)
  ENDDO
  DO L=1,NumCanopyLayers+1
    TAU_DirectSunLit(L,NY,NX) = plt_rad%TAU_DirectSunLit(L)
    TAU_DirectSunSha(L,NY,NX)         = plt_rad%TAU_DirectSunSha(L)
  ENDDO

  LeafStalkAreaAll_col(NY,NX)=plt_morph%LeafStalkAreaAll_col
  
  DO NZ=1,NP_col(NY,NX)
    DO L=1,NumCanopyLayers
      RadSWCanopyLAbsroption_pft(L,NZ,NY,NX)=plt_rad%RadSWCanopyLAbsroption_pft(L,NZ)      
    ENDDO
    StandDeadSurfArea_pft(NZ,NY,NX)      = plt_morph%StandDeadSurfArea_pft(NZ)
    LeafStalkAreaAct_pft(NZ,NY,NX)       = plt_morph%LeafStalkAreaAct_pft(NZ)
    RadSWCanopyAbsorption_pft(NZ,NY,NX)  = plt_rad%RadSWCanopyAbsorption_pft(NZ)
    RadPARCanopyAbsorption_pft(NZ,NY,NX) = plt_rad%RadPARCanopyAbsorption_pft(NZ)
    ClumpFactorNow_pft(NZ,NY,NX)         = plt_morph%ClumpFactorNow_pft(NZ)
    FracPARads2Canopy_pft(NZ,NY,NX)      = plt_rad%FracPARads2Canopy_pft(NZ)
    FracPARads2LiveCanopy_pft(NZ,NY,NX)  = plt_rad%FracPARads2LiveCanopy_pft(NZ)
    StomatalStress_pft(NZ,NY,NX)         = plt_biom%StomatalStress_pft(NZ)
    Eco_RadSW_col(NY,NX)                 = Eco_RadSW_col(NY,NX)+RadSWCanopyAbsorption_pft(NZ,NY,NX)
    RadSW_Canopy_col(NY,NX)         = RadSW_Canopy_col(NY,NX)+RadSWCanopyAbsorption_pft(NZ,NY,NX)
    LeafAreaSunlit_pft(NZ,NY,NX)    = plt_photo%LeafAreaSunlit_pft(NZ)
    PARSunlit_pft(NZ,NY,NX)         = plt_photo%PARSunlit_pft(NZ)
    PARSunsha_pft(NZ,NY,NX)         = plt_photo%PARSunsha_pft(NZ)
    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO K=1,MaxNodesPerBranch
        DO  L=1,NumCanopyLayers
          DO N=1,NumLeafInclinationClasses
            LeafAreaZsec_brch(N,L,K,NB,NZ,NY,NX)  = plt_morph%LeafAreaZsec_brch(N,L,K,NB,NZ)
            LeafEffArea_zsec(N,L,K,NB,NZ,NY,NX) = plt_photo%LeafEffArea_zsec(N,L,K,NB,NZ)                        
          ENDDO  
        ENDDO
      ENDDO  
    ENDDO
    DO L=1,NumCanopyLayers
      DO M=1,NumOfSkyAzimuthSects
        DO  N=1,NumLeafInclinationClasses
          RadDifPARAbsorption_zsec(N,M,L,NZ,NY,NX)=plt_rad%RadDifPARAbsorption_zsec(N,M,L,NZ)
          RadTotPARAbsorption_zsec(N,M,L,NZ,NY,NX)   =plt_rad%RadTotPARAbsorption_zsec(N,M,L,NZ)          
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  end subroutine PlantAPICanMRecv

end module PlantCanAPI
