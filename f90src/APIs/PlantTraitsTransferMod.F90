module PlantTraitsTransferMod
  ! Ordered traits transfers used by PlantAPISend/PlantAPIRecv.
  use data_kind_mod,              only: r8 => DAT_KIND_R8, yearIJ_type
  use EcoSiMParDataMod,           only: micpar,            pltpar
  use SoilPhysDataType,           only: SurfAlbedo_col,    SoilSurfDepZ_col
  use MiniMathMod,                only: AZMAX1,            safe_adb
  use DebugToolMod,               only: PrintInfo
  use PlantSiteAPIData,           only: plt_site
  use PlantPhotosynthesisAPIData, only: plt_photo
  use PlantRadiationAPIData,      only: plt_rad
  use PlantMorphologyAPIData,     only: plt_morph
  use PlantPhenologyAPIData,      only: plt_pheno
  use PlantAllometryAPIData,      only: plt_allom
  use PlantBiomassAPIData,        only: plt_biom
  use PlantEnergyWaterAPIData,    only: plt_ew
  use PlantDisturbanceAPIData,    only: plt_distb
  use PlantRootBGCAPIData,        only: plt_rbgc
  use NumericalAuxMod
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
  public :: SendPlantTraits
contains

  subroutine SendPlantTraits(NY,NX,NZ)
  use EcoSIMConfig, only : jsken=>jskenc,jcplx=>jcplxc
  integer, intent(in) :: NY,NX,NZ
  integer :: NB,L,M,N

    !plant properties begin
    plt_morph%StemSpecVolume_pft(NZ)             = StemSpecVolume_pft(NZ,NY,NX)
    plt_morph%Radius95pctMature_pft(NZ)          = Radius95pctMature_pft(NZ,NY,NX)
    plt_photo%iPlantPhotosynsType_pft(NZ)        = iPlantPhotosynsType_pft(NZ,NY,NX)
    plt_pheno%iPlantRootProfile_pft(NZ)           = iPlantRootProfile_pft(NZ,NY,NX)
    plt_morph%xylemPhi_min_pft(NZ)                = xylemPhi_min_pft(NZ,NY,NX)
    plt_morph%xylemPhi_mean_pft(NZ)                = xylemPhi_mean_pft(NZ,NY,NX)
    plt_morph%xylemPhi_max_pft(NZ)                = xylemPhi_max_pft(NZ,NY,NX)
    plt_pheno%iPlantPhenolPattern_pft(NZ)         = iPlantPhenolPattern_pft(NZ,NY,NX)
    plt_pheno%iPlantDevelopPattern_pft(NZ)        = iPlantDevelopPattern_pft(NZ,NY,NX)
    plt_pheno%iPlantPhenolType_pft(NZ)            = iPlantPhenolType_pft(NZ,NY,NX)
    plt_pheno%iEmbryophyteType_pft(NZ)             = iEmbryophyteType_pft(NZ,NY,NX)
    plt_pheno%iPlantPhotoperiodType_pft(NZ)       = iPlantPhotoperiodType_pft(NZ,NY,NX)
    plt_pheno%iPlantTurnoverPattern_pft(NZ)       = iPlantTurnoverPattern_pft(NZ,NY,NX)
    plt_pheno%iPlant2ndGrothPattern_pft(NZ)       = iPlant2ndGrothPattern_pft(NZ,NY,NX)
    plt_pheno%PlantInitThermoAdaptZone_pft(NZ)    = PlantInitThermoAdaptZone_pft(NZ,NY,NX)
    plt_morph%iPlantGrainType_pft(NZ)             = iPlantGrainType_pft(NZ,NY,NX)
    plt_morph%iPlantNfixType_pft(NZ)              = iPlantNfixType_pft(NZ,NY,NX)
    plt_morph%Myco_pft(NZ)                        = Myco_pft(NZ,NY,NX)
    plt_photo%VmaxSpecRubCarboxyRef_pft(NZ)       = VmaxSpecRubCarboxyRef_pft(NZ,NY,NX)
    plt_photo%VmaxRubOxyRef_pft(NZ)               = VmaxRubOxyRef_pft(NZ,NY,NX)
    plt_photo%VmaxPEPCarboxyRef_pft(NZ)           = VmaxPEPCarboxyRef_pft(NZ,NY,NX)
    plt_photo%XKCO2_pft(NZ)                       = XKCO2_pft(NZ,NY,NX)
    plt_photo%XKO2_pft(NZ)                        = XKO2_pft(NZ,NY,NX)
    plt_photo%Km4PEPCarboxy_pft(NZ)               = Km4PEPCarboxy_pft(NZ,NY,NX)
    plt_photo%LeafRubisco2Protein_pft(NZ)                = LeafRubisco2Protein_pft(NZ,NY,NX)
    plt_photo%LeafPEP2Protein_pft(NZ) = LeafPEP2Protein_pft(NZ,NY,NX)
    plt_photo%SpecLeafChlAct_pft(NZ)            = SpecLeafChlAct_pft(NZ,NY,NX)
    plt_photo%LeafProtein2Chl_pft(NZ)         = LeafProtein2Chl_pft(NZ,NY,NX)
    plt_photo%fMesophyllChlProtein_pft(NZ)         = fMesophyllChlProtein_pft(NZ,NY,NX)
    plt_photo%CanopyCi2CaRatio_pft(NZ)                  = CanopyCi2CaRatio_pft(NZ,NY,NX)

    plt_pheno%RefNodeInitRate_pft(NZ)        = RefNodeInitRate_pft(NZ,NY,NX)
    plt_pheno%RateRefLeafAppearance_pft(NZ)      = RateRefLeafAppearance_pft(NZ,NY,NX)
    plt_pheno%TCChill4Seed_pft(NZ)           = TCChill4Seed_pft(NZ,NY,NX)
    plt_morph%rLen2WidthLeaf_pft(NZ)         = rLen2WidthLeaf_pft(NZ,NY,NX)
    plt_pheno%NonstCMinConc2InitBranch_pft(NZ)   = NonstCMinConc2InitBranch_pft(NZ,NY,NX)
    plt_morph%ShootNodeNumAtPlanting_pft(NZ) = ShootNodeNumAtPlanting_pft(NZ,NY,NX)
    plt_pheno%CriticPhotoPeriod_pft(NZ)      = CriticPhotoPeriod_pft(NZ,NY,NX)
    plt_pheno%PhotoPeriodSens_pft(NZ)        = PhotoPeriodSens_pft(NZ,NY,NX)
    plt_morph%SLA1_pft(NZ)                   = SLA1_pft(NZ,NY,NX)
    plt_morph%PetolShethLen2Mass_pft(NZ)           = PetolShethLen2Mass_pft(NZ,NY,NX)
    plt_morph%NodeLenPergC_pft(NZ)               = NodeLenPergC_pft(NZ,NY,NX)
    DO  N=1,NumLeafInclinationClasses
      plt_morph%LeafAngleClass_pft(N,NZ)=LeafAngleClass_pft(N,NZ,NY,NX)
    ENDDO
    plt_morph%ClumpFactorInit_pft(NZ)     = ClumpFactorInit_pft(NZ,NY,NX)
    plt_morph%SineBranchAngle_pft(NZ)     = SineBranchAngle_pft(NZ,NY,NX)
    plt_morph%SinePetolShethAngle_pft(NZ)    = SinePetolShethAngle_pft(NZ,NY,NX)
    plt_morph%GrothStalkMaxSeedSites_pft(NZ) = GrothStalkMaxSeedSites_pft(NZ,NY,NX)
    plt_morph%MaxSeedNumPerSite_pft(NZ)   = MaxSeedNumPerSite_pft(NZ,NY,NX)
    plt_morph%SeedCMassMax_pft(NZ)        = SeedCMassMax_pft(NZ,NY,NX)
    plt_morph%SeedCMass_pft(NZ)           = SeedCMass_pft(NZ,NY,NX)
    plt_morph%SeedWidth2LenRatio_pft(NZ)  = SeedWidth2LenRatio_pft(NZ,NY,NX)
    plt_pheno%GrainFillRate25C_pft(NZ)    = GrainFillRate25C_pft(NZ,NY,NX)
    plt_biom%StandingDeadInitC_pft(NZ)    = StandingDeadInitC_pft(NZ,NY,NX)
    plt_morph%StalkAxialResist_pft(NZ) =StalkAxialResist_pft(NZ,NY,NX)
    plt_morph%RootSingleVesselRstaxial_pft(NZ)  = RootSingleVesselRstaxial_pft(NZ,NY,NX)
    plt_morph%RootSingleVesselArea_pft(NZ)    = RootSingleVesselArea_pft(NZ,NY,NX)
    plt_morph%RootVesselRadius_pft(NZ)        = RootVesselRadius_pft(NZ,NY,NX)
    plt_morph%AlphaVesselAxialResist_pft(NZ)  = AlphaVesselAxialResist_pft(NZ,NY,NX)
    !initial root values
    DO N=1,Myco_pft(NZ,NY,NX)
      plt_morph%Root1stMaxRadius_pft(N,NZ) = Root1stMaxRadius_pft(N,NZ,NY,NX)
      plt_morph%Root2ndMaxRadius_pft(N,NZ) = Root2ndMaxRadius_pft(N,NZ,NY,NX)
      plt_morph%RootPorosity_pft(N,NZ)     = RootPorosity_pft(N,NZ,NY,NX)
      plt_morph%RootRadialResist_pft(N,NZ) = RootRadialResist_pft(N,NZ,NY,NX)
      plt_morph%Root2ndAxialResist_pft(N,NZ)  = Root2ndAxialResist_pft(N,NZ,NY,NX)
      plt_rbgc%VmaxNH4Root_pft(N,NZ)       = VmaxNH4Root_pft(N,NZ,NY,NX)
      plt_rbgc%KmNH4Root_pft(N,NZ)         = KmNH4Root_pft(N,NZ,NY,NX)
      plt_rbgc%CMinNH4Root_pft(N,NZ)       = CMinNH4Root_pft(N,NZ,NY,NX)
      plt_rbgc%VmaxNO3Root_pft(N,NZ)       = VmaxNO3Root_pft(N,NZ,NY,NX)
      plt_rbgc%KmNO3Root_pft(N,NZ)         = KmNO3Root_pft(N,NZ,NY,NX)
      plt_rbgc%CminNO3Root_pft(N,NZ)       = CminNO3Root_pft(N,NZ,NY,NX)
      plt_rbgc%VmaxPO4Root_pft(N,NZ)       = VmaxPO4Root_pft(N,NZ,NY,NX)
      plt_rbgc%KmPO4Root_pft(N,NZ)         = KmPO4Root_pft(N,NZ,NY,NX)
      plt_rbgc%CMinPO4Root_pft(N,NZ)       = CMinPO4Root_pft(N,NZ,NY,NX)
    ENDDO
    plt_pheno%NonstCMinCon2InitRoot_pft(NZ)        = NonstCMinCon2InitRoot_pft(NZ,NY,NX)
    plt_pheno%ShootRootNonstElmConduts_pft(NZ) = ShootRootNonstElmConduts_pft(NZ,NY,NX)
    plt_morph%FineRootBranchFreq_pft(NZ)            = FineRootBranchFreq_pft(NZ,NY,NX)
    plt_morph%MediumRootBranchFreq_pft(NZ) = MediumRootBranchFreq_pft(NZ,NY,NX)
    plt_ew%OrganOsmoPsi0pt_pft(NZ)                = OrganOsmoPsi0pt_pft(NZ,NY,NX)
    plt_photo%RCS_pft(NZ)                           = RCS_pft(NZ,NY,NX)
    plt_photo%CuticleResist_pft(NZ)             = CuticleResist_pft(NZ,NY,NX)

    plt_allom%LeafBiomGrowthYld_pft(NZ)    = LeafBiomGrowthYld_pft(NZ,NY,NX)
    plt_allom%PetolShethBiomGrowthYld_pft(NZ) = PetolShethBiomGrowthYld_pft(NZ,NY,NX)
    plt_allom%StalkBiomGrowthYld_pft(NZ)   = StalkBiomGrowthYld_pft(NZ,NY,NX)
    plt_allom%ReserveBiomGrowthYld_pft(NZ) = ReserveBiomGrowthYld_pft(NZ,NY,NX)
    plt_allom%HuskBiomGrowthYld_pft(NZ)    = HuskBiomGrowthYld_pft(NZ,NY,NX)
    plt_allom%EarBiomGrowthYld_pft(NZ)     = EarBiomGrowthYld_pft(NZ,NY,NX)
    plt_allom%GrainBiomGrowthYld_pft(NZ)   = GrainBiomGrowthYld_pft(NZ,NY,NX)
    plt_allom%RootBiomGrosYld_pft(NZ)      = RootBiomGrosYld_pft(NZ,NY,NX)
    plt_allom%NoduGrowthYield_pft(NZ)      = NoduGrowthYield_pft(NZ,NY,NX)
    plt_allom%rNCLeaf_pft(NZ)              = rNCLeaf_pft(NZ,NY,NX)
    plt_allom%rNCSheath_pft(NZ)            = rNCSheath_pft(NZ,NY,NX)
    plt_allom%rNCStalk_pft(NZ)             = rNCStalk_pft(NZ,NY,NX)
    plt_allom%rECLiveCRoot_pft(:,NZ)         = rECLiveCRoot_pft(:,NZ,NY,NX)
    plt_allom%rECDeadCRoot_pft(:,NZ)       = rECDeadCRoot_pft(:,NZ,NY,NX)
    plt_allom%rNCReserve_pft(NZ)           = rNCReserve_pft(NZ,NY,NX)
    plt_allom%rNCHusk_pft(NZ)              = rNCHusk_pft(NZ,NY,NX)
    plt_allom%rNCEar_pft(NZ)               = rNCEar_pft(NZ,NY,NX)
    plt_allom%rNCGrain_pft(NZ)             = rNCGrain_pft(NZ,NY,NX)
    plt_allom%rNCRoot_pft(NZ)              = rNCRoot_pft(NZ,NY,NX)
    plt_allom%rNCNodule_pft(NZ)            = rNCNodule_pft(NZ,NY,NX)
    plt_allom%rNCLigRoot_pft(NZ)           = rNCLigRoot_pft(NZ,NY,NX)
    plt_allom%rPCLeaf_pft(NZ)              = rPCLeaf_pft(NZ,NY,NX)
    plt_allom%rPCSheath_pft(NZ)            = rPCSheath_pft(NZ,NY,NX)
    plt_allom%rPCStalk_pft(NZ)             = rPCStalk_pft(NZ,NY,NX)
    plt_biom%KLigMax_pft(NZ)               = KLigMax_pft(NZ,NY,NX)
    plt_biom%KLigMM_pft(NZ)                = KLigMM_pft(NZ,NY,NX)
    plt_allom%rPCReserve_pft(NZ)           = rPCReserve_pft(NZ,NY,NX)
    plt_allom%rPCHusk_pft(NZ)              = rPCHusk_pft(NZ,NY,NX)
    plt_allom%rPCEar_pft(NZ)               = rPCEar_pft(NZ,NY,NX)
    plt_allom%rPCGrain_pft(NZ)             = rPCGrain_pft(NZ,NY,NX)
    plt_allom%rPCRootr_pft(NZ)             = rPCRootr_pft(NZ,NY,NX)
    plt_allom%rPCNoduler_pft(NZ)           = rPCNoduler_pft(NZ,NY,NX)
    plt_allom%rPCLigRoot_pft(NZ)           = rPCLigRoot_pft(NZ,NY,NX)

    !plant properties end

    plt_morph%LeafStalkAreaAct_pft(NZ)                     = LeafStalkAreaAct_pft(NZ,NY,NX)
    plt_distb%iPlantingYear_pft(NZ)                     = iPlantingYear_pft(NZ,NY,NX)
    plt_distb%iPlantingDay_pft(NZ)                      = iPlantingDay_pft(NZ,NY,NX)
    plt_distb%iHarvestYear_pft(NZ)                      = iHarvestYear_pft(NZ,NY,NX)
    plt_rad%RadPARCanopyAbsorption_pft(NZ)              = RadPARCanopyAbsorption_pft(NZ,NY,NX)
    plt_rad%RadSWCanopyAbsorption_pft(NZ)               = RadSWCanopyAbsorption_pft(NZ,NY,NX)
    plt_ew%RainIntcptByCanopy_pft(NZ)                   = RainIntcptByCanopy_pft(NZ,NY,NX)
    plt_ew%SnowIntcptByCanopy_pft(NZ)                   = SnowIntcptByCanopy_pft(NZ,NY,NX)
    plt_site%PPatSeeding_pft(NZ)                        = PPatSeeding_pft(NZ,NY,NX)
    plt_distb%iHarvestDay_pft(NZ)                       = iHarvestDay_pft(NZ,NY,NX)
    plt_morph%ClumpFactorNow_pft(NZ)                    = ClumpFactorNow_pft(NZ,NY,NX)
    plt_site%DATAP(NZ)                                  = DATAP(NZ,NY,NX)
    plt_pheno%MatureGroup_pft(NZ)                       = MatureGroup_pft(NZ,NY,NX)
    plt_biom%AvgCanopyBiomC2Graze_pft(NZ)               = AvgCanopyBiomC2Graze_pft(NZ,NY,NX)

    DO NB=1,pltpar%MaxNumBranches
      plt_pheno%HourReq4LeafOut_brch(NB,NZ)=HourReq4LeafOut_brch(NB,NZ,NY,NX)
      plt_pheno%HourReq4LeafOff_brch(NB,NZ)=HourReq4LeafOff_brch(NB,NZ,NY,NX)
    ENDDO

    DO L=1,NumCanopyLayers
      plt_rad%RadSWCanopyLAbsroption_pft(L,NZ)  =RadSWCanopyLAbsroption_pft(L,NZ,NY,NX)
      DO  M=1,NumOfSkyAzimuthSects
        DO  N=1,NumLeafInclinationClasses
          plt_rad%RadTotPARAbsorption_zsec(N,M,L,NZ) = RadTotPARAbsorption_zsec(N,M,L,NZ,NY,NX)
          plt_rad%RadDifPARAbsorption_zsec(N,M,L,NZ) = RadDifPARAbsorption_zsec(N,M,L,NZ,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  end subroutine SendPlantTraits
end module PlantTraitsTransferMod
