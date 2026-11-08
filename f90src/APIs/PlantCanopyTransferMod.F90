module PlantCanopyTransferMod
  ! Ordered canopy transfers used by PlantAPISend/PlantAPIRecv.
  use data_kind_mod,              only: r8 => DAT_KIND_R8, yearIJ_type
  use EcoSiMParDataMod,           only: micpar,            pltpar
  use SoilPhysDataType,           only: SurfAlbedo_col,    SoilSurfDepZ_col
  use MiniMathMod,                only: AZMAX1,            safe_adb
  use DebugToolMod,               only: PrintInfo
  use PlantPhotosynthesisAPIData, only: plt_photo
  use PlantMorphologyAPIData,     only: plt_morph
  use PlantPhenologyAPIData,      only: plt_pheno
  use PlantAllometryAPIData,      only: plt_allom
  use PlantBiomassAPIData,        only: plt_biom
  use PlantEnergyWaterAPIData,    only: plt_ew
  use PlantBGCRatesAPIData,       only: plt_bgcr
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
  public :: ReceivePlantCanopy
  public :: SendPlantCanopy
contains

  subroutine ReceivePlantCanopy(NY,NX,NZ)
  use EcoSIMConfig, only : jsken=>jskenc,jcplx=>jcplxc
  integer, intent(in) :: NY,NX,NZ
  integer :: NB,NR,K,L,M,N,NE

    DO L=1,NK_col(NY,NX)
      RootFineFrac2Med_rpvr(:,L,NZ,NY,NX)                  = plt_morph%RootFineFrac2Med_rpvr(:,L,NZ)
      RootMediumLength_pvr(L,NZ,NY,NX)                  = plt_morph%RootMediumLength_pvr(L,NZ)
      RootSinkWeight_pvr(L,NZ,NY,NX)                    = plt_morph%RootSinkWeight_pvr(L,NZ)
      Root1stSinkWeight_pvr(L,NZ,NY,NX)                 = plt_morph%Root1stSinkWeight_pvr(L,NZ)
      RootMSinkWeight_pvr(L,NZ,NY,NX)                   = plt_morph%RootMSinkWeight_pvr(L,NZ)
      Root2ndSinkWeight_pvr(L,1:pltpar%jroots,NZ,NY,NX) = plt_morph%Root2ndSinkWeight_pvr(L,1:pltpar%jroots,NZ)
      RootNodulStrutElms_rpvr(1:NumPlantChemElms,L,NZ,NY,NX) = plt_biom%RootNodulStrutElms_rpvr(1:NumPlantChemElms,L,NZ)
      RootNodulNonstElms_rpvr(1:NumPlantChemElms,L,NZ,NY,NX) = plt_biom%RootNodulNonstElms_rpvr(1:NumPlantChemElms,L,NZ)
      RootN2Fix_pvr(L,NZ,NY,NX)                              = plt_bgcr%RootN2Fix_pvr(L,NZ)
      fTgrowRootP_vr(L,NZ,NY,NX)                             = plt_pheno%fTgrowRootP_vr(L,NZ)
      RootN2Fix_vr(L,NY,NX)                                  = RootN2Fix_vr(L,NY,NX)+RootN2Fix_pvr(L,NZ,NY,NX)
      DO NE=1,NumPlantChemElms
        RootMedStruct_pvr(NE,L,NZ,NY,NX)         = plt_biom%RootMedStruct_pvr(NE,L,NZ)
        Root1stActStruct_pvr(NE,L,NZ,NY,NX)      = plt_biom%Root1stActStruct_pvr(NE,L,NZ)
        Root1stLigStruct_pvr(NE,L,NZ,NY,NX)      = plt_biom%Root1stLigStruct_pvr(NE,L,NZ)

        DO N=1,Myco_pft(NZ,NY,NX)
          RootMycoMassElm_pvr(NE,N,L,NZ,NY,NX)                 = plt_biom%RootMycoMassElm_pvr(NE,N,L,NZ)
          RootMycoMassElm_vr(NE,N,L,NY,NX)= RootMycoMassElm_vr(NE,N,L,NY,NX)+ plt_biom%RootMycoMassElm_pvr(NE,N,L,NZ)
        ENDDO
      ENDDO
    ENDDO
    DO L=1,NumCanopyLayers
      CanopyLeafAreaZ_pft(L,NZ,NY,NX)        = plt_morph%CanopyLeafAreaZ_pft(L,NZ)
      CanopyLeafCLyr_pft(L,NZ,NY,NX)         = plt_biom%CanopyLeafCLyr_pft(L,NZ)
      CanopyStemSurfAreaZ_pft(L,NZ,NY,NX)    = plt_morph%CanopyStemSurfAreaZ_pft(L,NZ)
      CanopySurfAreaProfDead_pft(L,NZ,NY,NX) = plt_morph%CanopySurfAreaProfDead_pft(L,NZ)
    ENDDO

    DO L=0,NL_col(NY,NX)
      DO K=1,micpar%NumOfPlantLitrCmplxs
        DO M=1,jsken
          LitrfallElms_pvr(1:NumPlantChemElms,M,K,L,NZ,NY,NX)=plt_bgcr%LitrfallElms_pvr(1:NumPlantChemElms,M,K,L,NZ)
        ENDDO
      ENDDO
    ENDDO

    DO NB=1,pltpar%MaxNumBranches
      DO NE=1,NumPlantChemElms
        CanopyNonstElms_brch(NE,NB,NZ,NY,NX) = plt_biom%CanopyNonstElms_brch(NE,NB,NZ)
        EarStrutElms_brch(NE,NB,NZ,NY,NX)    = plt_biom%EarStrutElms_brch(NE,NB,NZ)
      ENDDO
      C4PhotoShootNonstC_brch(NB,NZ,NY,NX)=plt_biom%C4PhotoShootNonstC_brch(NB,NZ)
    ENDDO

    DO NR=1,pltpar%MaxNumRootAxes
      SapFlowVLinear_rpvr(:,NR,NZ,NY,NX) = plt_ew%SapFlowVLinear_rpvr(:,NR,NZ)
      RootSegAges_raxes(:,NR,NZ,NY,NX)      = plt_morph%RootSegAges_raxes(:,NR,NZ)
      RootSeglengths_raxes(:,NR,NZ,NY,NX)   = plt_morph%RootSeglengths_raxes(:,NR,NZ)
      NActiveRootSegs_raxes(NR,NZ,NY,NX)  = plt_morph%NActiveRootSegs_raxes(NR,NZ)
      RootSegBaseDepth_raxes(NR,NZ,NY,NX) = plt_morph%RootSegBaseDepth_raxes(NR,NZ)
      IndRootSegBase_raxes(NR,NZ,NY,NX)   = plt_morph%IndRootSegBase_raxes(NR,NZ)
      IndRootSegTip_raxes(NR,NZ,NY,NX)    = plt_morph%IndRootSegTip_raxes(NR,NZ)
    ENDDO
    DO NB=1,pltpar%MaxNumBranches
      DO NE=1,NumPlantChemElms
        CanopyNodulNonstElms_brch(NE,NB,NZ,NY,NX) = plt_biom%CanopyNodulNonstElms_brch(NE,NB,NZ)
        ShootElms_brch(NE,NB,NZ,NY,NX)            = plt_biom%ShootElms_brch(NE,NB,NZ)
        PetolShethStrutElms_brch(NE,NB,NZ,NY,NX)      = plt_biom%PetolShethStrutElms_brch(NE,NB,NZ)
        StalkStrutElms_brch(NE,NB,NZ,NY,NX)       = plt_biom%StalkStrutElms_brch(NE,NB,NZ)
        LeafPetoNonstElmConc_brch(NE,NB,NZ,NY,NX) = plt_biom%LeafPetoNonstElmConc_brch(NE,NB,NZ)
        LeafStrutElms_brch(NE,NB,NZ,NY,NX)        = plt_biom%LeafStrutElms_brch(NE,NB,NZ)
        StalkRsrvElms_brch(NE,NB,NZ,NY,NX)        = plt_biom%StalkRsrvElms_brch(NE,NB,NZ)
        HuskStrutElms_brch(NE,NB,NZ,NY,NX)        = plt_biom%HuskStrutElms_brch(NE,NB,NZ)
        GrainStrutElms_brch(NE,NB,NZ,NY,NX)       = plt_biom%GrainStrutElms_brch(NE,NB,NZ)
        CanopyNodulStrutElms_brch(NE,NB,NZ,NY,NX) = plt_biom%CanopyNodulStrutElms_brch(NE,NB,NZ)
      ENDDO
    ENDDO

    DO NB=1,pltpar%MaxNumBranches
      DO L=1,NumCanopyLayers
        CanopyStalkSurfArea_lbrch(L,NB,NZ,NY,NX)=plt_morph%CanopyStalkSurfArea_lbrch(L,NB,NZ)
      ENDDO
      Hours2LeafOut_brch(NB,NZ,NY,NX) = plt_pheno%Hours2LeafOut_brch(NB,NZ)
      LeafAreaLive_brch(NB,NZ,NY,NX)  = plt_morph%LeafAreaLive_brch(NB,NZ)
      LeafAreaDying_brch(NB,NZ,NY,NX) = plt_morph%LeafAreaDying_brch(NB,NZ)

      HourFailGrainFill_brch(NB,NZ,NY,NX)                         = plt_pheno%HourFailGrainFill_brch(NB,NZ)
      RubiscoActivity_brch(NB,NZ,NY,NX)                           = plt_photo%RubiscoActivity_brch(NB,NZ)
      GrainFillDowreg_brch(NB,NZ,NY,NX)                          = plt_photo%GrainFillDowreg_brch(NB,NZ)
      HoursDoingRemob_brch(NB,NZ,NY,NX)                           = plt_pheno%HoursDoingRemob_brch(NB,NZ)
      MatureGroup_brch(NB,NZ,NY,NX)                               = plt_pheno%MatureGroup_brch(NB,NZ)
      NodeNumNormByMatgrp_brch(NB,NZ,NY,NX)                       = plt_pheno%NodeNumNormByMatgrp_brch(NB,NZ)
      ReprodNodeNumNormByMatrgrp_brch(NB,NZ,NY,NX)                = plt_pheno%ReprodNodeNumNormByMatrgrp_brch(NB,NZ)
      PotentialSeedSites_brch(NB,NZ,NY,NX)                        = plt_morph%PotentialSeedSites_brch(NB,NZ)
      SetNumberSeeds_brch(NB,NZ,NY,NX)                                = plt_morph%SetNumberSeeds_brch(NB,NZ)
      SingleGrainMeanBiomC_brch(NB,NZ,NY,NX)                        = plt_allom%SingleGrainMeanBiomC_brch(NB,NZ)
      CanPBranchHeight(NB,NZ,NY,NX)                               = plt_morph%CanPBranchHeight(NB,NZ)
      doRemobilization_brch(NB,NZ,NY,NX)                          = plt_pheno%doRemobilization_brch(NB,NZ)
      isPlantBranchAlive_brch(NB,NZ,NY,NX)                         = plt_pheno%isPlantBranchAlive_brch(NB,NZ)
      doPlantLeaveOff_brch(NB,NZ,NY,NX)                           = plt_pheno%doPlantLeaveOff_brch(NB,NZ)
      EnablePlantLeafOut_brch(NB,NZ,NY,NX)                            = plt_pheno%EnablePlantLeafOut_brch(NB,NZ)
      doInitLeafOut_brch(NB,NZ,NY,NX)                             = plt_pheno%doInitLeafOut_brch(NB,NZ)
      doSenescence_brch(NB,NZ,NY,NX)                              = plt_pheno%doSenescence_brch(NB,NZ)
      Prep4Literfall_brch(NB,NZ,NY,NX)                            = plt_pheno%Prep4Literfall_brch(NB,NZ)
      Hours4LiterfalAftMature_brch(NB,NZ,NY,NX)                   = plt_pheno%Hours4LiterfalAftMature_brch(NB,NZ)
      KHiestGroLeafNode_brch(NB,NZ,NY,NX)                         = plt_pheno%KHiestGroLeafNode_brch(NB,NZ)
      KLeafNumber_brch(NB,NZ,NY,NX)                               = plt_morph%KLeafNumber_brch(NB,NZ)
      KMinNumLeaf4GroAlloc_brch(NB,NZ,NY,NX)                      = plt_morph%KMinNumLeaf4GroAlloc_brch(NB,NZ)
      KLowestGroLeafNode_brch(NB,NZ,NY,NX)                        = plt_pheno%KLowestGroLeafNode_brch(NB,NZ)
      BranchNumerID_brch(NB,NZ,NY,NX)                              = plt_morph%BranchNumerID_brch(NB,NZ)
      ShootNodeNum_brch(NB,NZ,NY,NX)                              = plt_morph%ShootNodeNum_brch(NB,NZ)
      ShootNodeNumAtInitFloral_brch(NB,NZ,NY,NX)                        = plt_morph%ShootNodeNumAtInitFloral_brch(NB,NZ)
      ShootNodeNumAtAnthesis_brch(NB,NZ,NY,NX)                      = plt_morph%ShootNodeNumAtAnthesis_brch(NB,NZ)
      LeafElmntRemobFlx_brch(1:NumPlantChemElms,NB,NZ,NY,NX)      = plt_pheno%LeafElmntRemobFlx_brch(1:NumPlantChemElms,NB,NZ)
      LeafSenescInitialElms_brch(1:NumPlantChemElms,NB,NZ,NY,NX)      = plt_pheno%LeafSenescInitialElms_brch(1:NumPlantChemElms,NB,NZ)
      PetolShethChemElmRemobFlx_brch(1:NumPlantChemElms,NB,NZ,NY,NX) = plt_pheno%PetolShethChemElmRemobFlx_brch(1:NumPlantChemElms,NB,NZ)
      PetolSenescInitialElms_brch(1:NumPlantChemElms,NB,NZ,NY,NX) = plt_pheno%PetolSenescInitialElms_brch(1:NumPlantChemElms,NB,NZ)

      TotalNodeNumNormByMatgrp_brch(NB,NZ,NY,NX)               = plt_pheno%TotalNodeNumNormByMatgrp_brch(NB,NZ)
      TotReproNodeNumNormByMatrgrp_brch(NB,NZ,NY,NX)           = plt_pheno%TotReproNodeNumNormByMatrgrp_brch(NB,NZ)
      Hours4LenthenPhotoPeriod_brch(NB,NZ,NY,NX)               = plt_pheno%Hours4LenthenPhotoPeriod_brch(NB,NZ)
      Hours4ShortenPhotoPeriod_brch(NB,NZ,NY,NX)               = plt_pheno%Hours4ShortenPhotoPeriod_brch(NB,NZ)
      Hours4Leafout_brch(NB,NZ,NY,NX)                          = plt_pheno%Hours4Leafout_brch(NB,NZ)
      Hours4LeafOff_brch(NB,NZ,NY,NX)                          = plt_pheno%Hours4LeafOff_brch(NB,NZ)
      NumOfLeaves_brch(NB,NZ,NY,NX)                            = plt_morph%NumOfLeaves_brch(NB,NZ)
      LeafNumberAtFloralInit_brch(NB,NZ,NY,NX)                 = plt_pheno%LeafNumberAtFloralInit_brch(NB,NZ)
      CanopyLeafSheathC_brch(NB,NZ,NY,NX)                      = plt_biom%CanopyLeafSheathC_brch(NB,NZ)
      dReproNodeNumNormByMatG_brch(NB,NZ,NY,NX)                = plt_pheno%dReproNodeNumNormByMatG_brch(NB,NZ)
!      LeafChemElmRemob_brch(1:NumPlantChemElms,NB,NZ,NY,NX)    = plt_biom%LeafChemElmRemob_brch(1:NumPlantChemElms,NB,NZ)
      SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ,NY,NX) = plt_biom%SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ)
      SapwoodBiomassC_brch(NB,NZ,NY,NX)                        = plt_biom%SapwoodBiomassC_brch(NB,NZ)
      CanopyNLimFactor_brch(NB,NZ,NY,NX)                       = plt_bgcr%CanopyNLimFactor_brch(NB,NZ)
      CanopyPLimFactor_brch(NB,NZ,NY,NX)                       = plt_bgcr%CanopyPLimFactor_brch(NB,NZ)

      DO K=0,MaxNodesPerBranch
        LeafArea_node(K,NB,NZ,NY,NX)                           = plt_morph%LeafArea_node(K,NB,NZ)
        StalkNodeVertLength_brch(K,NB,NZ,NY,NX)                = plt_morph%StalkNodeVertLength_brch(K,NB,NZ)
        StalkNodeHeight_brch(K,NB,NZ,NY,NX)                    = plt_morph%StalkNodeHeight_brch(K,NB,NZ)
        PetoleLength_node(K,NB,NZ,NY,NX)                     = plt_morph%PetoleLength_node(K,NB,NZ)
        StructInternodeElms_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX) = plt_biom%StructInternodeElms_brch(1:NumPlantChemElms,K,NB,NZ)
        LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)      = plt_biom%LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)
        LeafProteinC_node(K,NB,NZ,NY,NX)                      = plt_biom%LeafProteinC_node(K,NB,NZ)
        PetolShethElmntNode_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)   = plt_biom%PetolShethElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)
        PetoleProteinC_node(K,NB,NZ,NY,NX)                    = plt_biom%PetoleProteinC_node(K,NB,NZ)
      ENDDO
      DO  L=1,NumCanopyLayers
        DO N=1,NumLeafInclinationClasses
          StemAreaZsec_brch(N,L,NB,NZ,NY,NX)=plt_morph%StemAreaZsec_brch(N,L,NB,NZ)
        ENDDO
      ENDDO
      DO K=0,MaxNodesPerBranch
        DO  L=1,NumCanopyLayers
          CanopyLeafArea_lnode(L,K,NB,NZ,NY,NX)                         = plt_morph%CanopyLeafArea_lnode(L,K,NB,NZ)
          LeafLayerElms_node(1:NumPlantChemElms,L,K,NB,NZ,NY,NX) = plt_biom%LeafLayerElms_node(1:NumPlantChemElms,L,K,NB,NZ)
        ENDDO
      ENDDO
      DO M=1,pltpar%NumGrowthStages
        iPlantCalendar_brch(M,NB,NZ,NY,NX)=plt_pheno%iPlantCalendar_brch(M,NB,NZ)
      ENDDO

      DO K=1,MaxNodesPerBranch
        DO  L=1,NumCanopyLayers
          DO N=1,NumLeafInclinationClasses
            LeafAreaZsec_brch(N,L,K,NB,NZ,NY,NX)  = plt_morph%LeafAreaZsec_brch(N,L,K,NB,NZ)
          ENDDO
        ENDDO

        CPOOL3_node(K,NB,NZ,NY,NX)                   = plt_photo%CPOOL3_node(K,NB,NZ)
        CPOOL4_node(K,NB,NZ,NY,NX)                   = plt_photo%CPOOL4_node(K,NB,NZ)
        CMassCO2BundleSheath_node(K,NB,NZ,NY,NX)     = plt_photo%CMassCO2BundleSheath_node(K,NB,NZ)
        CO2CompenPoint_node(K,NB,NZ,NY,NX)           = plt_photo%CO2CompenPoint_node(K,NB,NZ)
        RubiscoCarboxyEff_node(K,NB,NZ,NY,NX)        = plt_photo%RubiscoCarboxyEff_node(K,NB,NZ)
        C4CarboxyEff_node(K,NB,NZ,NY,NX)             = plt_photo%C4CarboxyEff_node(K,NB,NZ)
        LigthSatCarboxyRate_node(K,NB,NZ,NY,NX)      = plt_photo%LigthSatCarboxyRate_node(K,NB,NZ)
        LigthSatC4CarboxyRate_node(K,NB,NZ,NY,NX)    = plt_photo%LigthSatC4CarboxyRate_node(K,NB,NZ)
        NutrientCtrlonC4Carboxy_node(K,NB,NZ,NY,NX)  = plt_photo%NutrientCtrlonC4Carboxy_node(K,NB,NZ)
        CMassHCO3BundleSheath_node(K,NB,NZ,NY,NX)    = plt_photo%CMassHCO3BundleSheath_node(K,NB,NZ)
        Vmax4RubiscoCarboxy_node(K,NB,NZ,NY,NX)       = plt_photo%Vmax4RubiscoCarboxy_node(K,NB,NZ)
        ProteinCperm2LeafArea_node(K,NB,NZ,NY,NX)     = plt_photo%ProteinCperm2LeafArea_node(K,NB,NZ)
        CO2lmtRubiscoCarboxyRate_node(K,NB,NZ,NY,NX) = plt_photo%CO2lmtRubiscoCarboxyRate_node(K,NB,NZ)
        Vmax4PEPCarboxy_node(K,NB,NZ,NY,NX)           = plt_photo%Vmax4PEPCarboxy_node(K,NB,NZ)
        CO2lmtPEPCarboxyRate_node(K,NB,NZ,NY,NX)     = plt_photo%CO2lmtPEPCarboxyRate_node(K,NB,NZ)
      ENDDO
    ENDDO
    DO M=1,jsken
      StandDeadCompKElms_pft(1:NumPlantChemElms,M,NZ,NY,NX)=plt_biom%StandDeadCompKElms_pft(1:NumPlantChemElms,M,NZ)
    ENDDO

  end subroutine ReceivePlantCanopy

  subroutine SendPlantCanopy(NY,NX,NZ)
  use EcoSIMConfig, only : jsken=>jskenc,jcplx=>jcplxc
  integer, intent(in) :: NY,NX,NZ
  integer :: NB,NR,K,L,M,N,NE

    plt_morph%RootSinkWeight_pvr(:,NZ)=RootSinkWeight_pvr(:,NZ,NY,NX)
    DO L=1,NK_col(NY,NX)
      plt_rbgc%GroSrcRootStress_pvr(L,NZ)  = GroSrcRootStress_pvr(L,NZ,NY,NX)
      plt_morph%RootMediumLength_pvr(L,NZ) = RootMediumLength_pvr(L,NZ,NY,NX)
      plt_morph%RootFineFrac2Med_rpvr(:,L,NZ) = RootFineFrac2Med_rpvr(:,L,NZ,NY,NX)
      DO K=1,jcplx
        DO N=1,Myco_pft(NZ,NY,NX)
          DO NE=1,NumPlantChemElms
            plt_rbgc%Soil2RootMycoExudE_pvr(NE,N,K,L,NZ)=Soil2RootMycoExudE_pvr(NE,N,K,L,NZ,NY,NX)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    DO NR=1,pltpar%MaxNumRootAxes
      plt_morph%RootSegAges_raxes(:,NR,NZ)   =RootSegAges_raxes(:,NR,NZ,NY,NX)
      plt_ew%SapFlowVLinear_rpvr(:,NR,NZ) =SapFlowVLinear_rpvr(:,NR,NZ,NY,NX)
      plt_morph%RootSeglengths_raxes(:,NR,NZ)  = RootSeglengths_raxes(:,NR,NZ,NY,NX)
      plt_morph%NActiveRootSegs_raxes(NR,NZ) = NActiveRootSegs_raxes(NR,NZ,NY,NX)
      plt_morph%RootSegBaseDepth_raxes(NR,NZ) = RootSegBaseDepth_raxes(NR,NZ,NY,NX)
      plt_morph%IndRootSegBase_raxes(NR,NZ)  = IndRootSegBase_raxes(NR,NZ,NY,NX)
      plt_morph%IndRootSegTip_raxes(NR,NZ)   = IndRootSegTip_raxes(NR,NZ,NY,NX)
    ENDDO
    DO NB=1,pltpar%MaxNumBranches
      ! Inactive branches retain their saved counter until phenology updates it.
      plt_morph%KLeafNumber_brch(NB,NZ) = KLeafNumber_brch(NB,NZ,NY,NX)
      plt_photo%RubiscoActivity_brch(NB,NZ) = RubiscoActivity_brch(NB,NZ,NY,NX)
      DO K=1,MaxNodesPerBranch
        DO  L=1,NumCanopyLayers
          DO N=1,NumLeafInclinationClasses
            plt_photo%LeafEffArea_zsec(N,L,K,NB,NZ)=LeafEffArea_zsec(N,L,K,NB,NZ,NY,NX)
          ENDDO
        ENDDO
      ENDDO

      DO NE=1,NumPlantChemElms
        plt_biom%CanopyNonstElms_brch(NE,NB,NZ)      = CanopyNonstElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%CanopyNodulNonstElms_brch(NE,NB,NZ) = CanopyNodulNonstElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%ShootElms_brch(NE,NB,NZ)            = ShootElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%LeafPetoNonstElmConc_brch(NE,NB,NZ) = LeafPetoNonstElmConc_brch(NE,NB,NZ,NY,NX)
        plt_biom%PetolShethStrutElms_brch(NE,NB,NZ)      = PetolShethStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%StalkStrutElms_brch(NE,NB,NZ)       = StalkStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%LeafStrutElms_brch(NE,NB,NZ)        = LeafStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%StalkRsrvElms_brch(NE,NB,NZ)        = StalkRsrvElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%HuskStrutElms_brch(NE,NB,NZ)        = HuskStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%GrainStrutElms_brch(NE,NB,NZ)       = GrainStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%EarStrutElms_brch(NE,NB,NZ)         = EarStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%CanopyNodulStrutElms_brch(NE,NB,NZ) = CanopyNodulStrutElms_brch(NE,NB,NZ,NY,NX)

      ENDDO
      plt_biom%C4PhotoShootNonstC_brch(NB,NZ)                     = C4PhotoShootNonstC_brch(NB,NZ,NY,NX)
      plt_biom%CanopyLeafSheathC_brch(NB,NZ)                      = CanopyLeafSheathC_brch(NB,NZ,NY,NX)
      plt_biom%SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ) = SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ,NY,NX)
      plt_biom%SapwoodBiomassC_brch(NB,NZ)                        = SapwoodBiomassC_brch(NB,NZ,NY,NX)

      plt_photo%GrainFillDowreg_brch(NB,NZ) = GrainFillDowreg_brch(NB,NZ,NY,NX)
      plt_pheno%Hours2LeafOut_brch(NB,NZ)    = Hours2LeafOut_brch(NB,NZ,NY,NX)
      plt_morph%LeafAreaDying_brch(NB,NZ)    = LeafAreaDying_brch(NB,NZ,NY,NX)
      plt_morph%LeafAreaLive_brch(NB,NZ)     = LeafAreaLive_brch(NB,NZ,NY,NX)

      plt_pheno%dReproNodeNumNormByMatG_brch(NB,NZ)    = dReproNodeNumNormByMatG_brch(NB,NZ,NY,NX)
      plt_pheno%HourFailGrainFill_brch(NB,NZ)          = HourFailGrainFill_brch(NB,NZ,NY,NX)
      plt_pheno%HoursDoingRemob_brch(NB,NZ)            = HoursDoingRemob_brch(NB,NZ,NY,NX)
      plt_pheno%MatureGroup_brch(NB,NZ)                = MatureGroup_brch(NB,NZ,NY,NX)
      plt_pheno%NodeNumNormByMatgrp_brch(NB,NZ)        = NodeNumNormByMatgrp_brch(NB,NZ,NY,NX)
      plt_pheno%ReprodNodeNumNormByMatrgrp_brch(NB,NZ) = ReprodNodeNumNormByMatrgrp_brch(NB,NZ,NY,NX)
      plt_morph%PotentialSeedSites_brch(NB,NZ)         = PotentialSeedSites_brch(NB,NZ,NY,NX)
      plt_morph%SetNumberSeeds_brch(NB,NZ)                 = SetNumberSeeds_brch(NB,NZ,NY,NX)
      plt_allom%SingleGrainMeanBiomC_brch(NB,NZ)         = SingleGrainMeanBiomC_brch(NB,NZ,NY,NX)
      plt_morph%CanPBranchHeight(NB,NZ)                = CanPBranchHeight(NB,NZ,NY,NX)
      plt_pheno%isPlantBranchAlive_brch(NB,NZ)          = isPlantBranchAlive_brch(NB,NZ,NY,NX)
      plt_pheno%doRemobilization_brch(NB,NZ)           = doRemobilization_brch(NB,NZ,NY,NX)
      plt_pheno%doPlantLeaveOff_brch(NB,NZ)            = doPlantLeaveOff_brch(NB,NZ,NY,NX)
      plt_pheno%EnablePlantLeafOut_brch(NB,NZ)             = EnablePlantLeafOut_brch(NB,NZ,NY,NX)
      plt_pheno%doInitLeafOut_brch(NB,NZ)              = doInitLeafOut_brch(NB,NZ,NY,NX)
      plt_pheno%doSenescence_brch(NB,NZ)               = doSenescence_brch(NB,NZ,NY,NX)
      plt_pheno%Prep4Literfall_brch(NB,NZ)             = Prep4Literfall_brch(NB,NZ,NY,NX)
      plt_pheno%Hours4LiterfalAftMature_brch(NB,NZ)    = Hours4LiterfalAftMature_brch(NB,NZ,NY,NX)
      plt_pheno%KHiestGroLeafNode_brch(NB,NZ)          = KHiestGroLeafNode_brch(NB,NZ,NY,NX)
      plt_pheno%KLowestGroLeafNode_brch(NB,NZ)         = KLowestGroLeafNode_brch(NB,NZ,NY,NX)
      plt_morph%BranchNumerID_brch(NB,NZ)               = BranchNumerID_brch(NB,NZ,NY,NX)
      plt_morph%ShootNodeNum_brch(NB,NZ)               = ShootNodeNum_brch(NB,NZ,NY,NX)
      plt_morph%ShootNodeNumAtInitFloral_brch(NB,NZ)                        = ShootNodeNumAtInitFloral_brch(NB,NZ,NY,NX)
      plt_morph%ShootNodeNumAtAnthesis_brch(NB,NZ)                      = ShootNodeNumAtAnthesis_brch(NB,NZ,NY,NX)
      plt_pheno%LeafElmntRemobFlx_brch(1:NumPlantChemElms,NB,NZ)      = LeafElmntRemobFlx_brch(1:NumPlantChemElms,NB,NZ,NY,NX)
      plt_pheno%LeafSenescInitialElms_brch(1:NumPlantChemElms,NB,NZ)      = LeafSenescInitialElms_brch(1:NumPlantChemElms,NB,NZ,NY,NX)
      plt_pheno%PetolShethChemElmRemobFlx_brch(1:NumPlantChemElms,NB,NZ) = PetolShethChemElmRemobFlx_brch(1:NumPlantChemElms,NB,NZ,NY,NX)
      plt_pheno%PetolSenescInitialElms_brch(1:NumPlantChemElms,NB,NZ) = PetolSenescInitialElms_brch(1:NumPlantChemElms,NB,NZ,NY,NX)
      plt_pheno%TotalNodeNumNormByMatgrp_brch(NB,NZ)                  = TotalNodeNumNormByMatgrp_brch(NB,NZ,NY,NX)
      plt_pheno%TotReproNodeNumNormByMatrgrp_brch(NB,NZ)              = TotReproNodeNumNormByMatrgrp_brch(NB,NZ,NY,NX)
      plt_pheno%LeafNumberAtFloralInit_brch(NB,NZ)                    = LeafNumberAtFloralInit_brch(NB,NZ,NY,NX)
      plt_morph%NumOfLeaves_brch(NB,NZ)                               = NumOfLeaves_brch(NB,NZ,NY,NX)
      plt_pheno%Hours4LenthenPhotoPeriod_brch(NB,NZ)                  = Hours4LenthenPhotoPeriod_brch(NB,NZ,NY,NX)
      plt_pheno%Hours4ShortenPhotoPeriod_brch(NB,NZ)                  = Hours4ShortenPhotoPeriod_brch(NB,NZ,NY,NX)
      plt_pheno%Hours4Leafout_brch(NB,NZ)                             = Hours4Leafout_brch(NB,NZ,NY,NX)
      plt_pheno%Hours4LeafOff_brch(NB,NZ)                             = Hours4LeafOff_brch(NB,NZ,NY,NX)
      DO M=1,pltpar%NumGrowthStages
        plt_pheno%iPlantCalendar_brch(M,NB,NZ)=iPlantCalendar_brch(M,NB,NZ,NY,NX)
      ENDDO
      DO K=1,MaxNodesPerBranch
        plt_photo%CPOOL3_node(K,NB,NZ)                = CPOOL3_node(K,NB,NZ,NY,NX)
        plt_photo%CPOOL4_node(K,NB,NZ)                = CPOOL4_node(K,NB,NZ,NY,NX)
        plt_photo%CMassCO2BundleSheath_node(K,NB,NZ)  = CMassCO2BundleSheath_node(K,NB,NZ,NY,NX)
        plt_photo%CO2CompenPoint_node(K,NB,NZ)        = CO2CompenPoint_node(K,NB,NZ,NY,NX)
        plt_photo%RubiscoCarboxyEff_node(K,NB,NZ)     = RubiscoCarboxyEff_node(K,NB,NZ,NY,NX)
        plt_photo%CMassHCO3BundleSheath_node(K,NB,NZ) = CMassHCO3BundleSheath_node(K,NB,NZ,NY,NX)

        DO  L=1,NumCanopyLayers
          DO N=1,NumLeafInclinationClasses
            plt_morph%LeafAreaZsec_brch(N,L,K,NB,NZ)=LeafAreaZsec_brch(N,L,K,NB,NZ,NY,NX)
          ENDDO
        ENDDO

      ENDDO
      DO K=0,MaxNodesPerBranch
        plt_morph%LeafArea_node(K,NB,NZ)                              = LeafArea_node(K,NB,NZ,NY,NX)
        plt_morph%StalkNodeVertLength_brch(K,NB,NZ)                   = StalkNodeVertLength_brch(K,NB,NZ,NY,NX)
        plt_morph%PetoleLength_node(K,NB,NZ)                        = PetoleLength_node(K,NB,NZ,NY,NX)
        plt_morph%StalkNodeHeight_brch(K,NB,NZ)                       = StalkNodeHeight_brch(K,NB,NZ,NY,NX)

        plt_biom%StructInternodeElms_brch(1:NumPlantChemElms,K,NB,NZ) = StructInternodeElms_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)
        plt_biom%LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)       = LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)
        plt_biom%LeafProteinC_node(K,NB,NZ)                           = LeafProteinC_node(K,NB,NZ,NY,NX)
        plt_biom%PetolShethElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)    = PetolShethElmntNode_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)
        plt_biom%PetoleProteinC_node(K,NB,NZ)                     = PetoleProteinC_node(K,NB,NZ,NY,NX)
      ENDDO

      DO K=0,MaxNodesPerBranch
        DO  L=1,NumCanopyLayers
          plt_morph%CanopyLeafArea_lnode(L,K,NB,NZ)                        = CanopyLeafArea_lnode(L,K,NB,NZ,NY,NX)
          plt_biom%LeafLayerElms_node(1:NumPlantChemElms,L,K,NB,NZ) = LeafLayerElms_node(1:NumPlantChemElms,L,K,NB,NZ,NY,NX)
        ENDDO
      ENDDO
      DO  L=1,NumCanopyLayers
        plt_morph%CanopyStalkSurfArea_lbrch(L,NB,NZ)=CanopyStalkSurfArea_lbrch(L,NB,NZ,NY,NX)
      ENDDO
    enddo

  end subroutine SendPlantCanopy
end module PlantCanopyTransferMod
