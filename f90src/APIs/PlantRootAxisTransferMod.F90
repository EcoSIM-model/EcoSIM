module PlantRootAxisTransferMod
  ! Ordered rootaxis transfers used by PlantAPISend/PlantAPIRecv.
  use data_kind_mod,             only: r8 => DAT_KIND_R8, yearIJ_type
  use EcoSiMParDataMod,          only: micpar,            pltpar
  use SoilPhysDataType,          only: SurfAlbedo_col,    SoilSurfDepZ_col
  use MiniMathMod,               only: AZMAX1,            safe_adb
  use DebugToolMod,              only: PrintInfo
  use PlantMorphologyAPIData,    only: plt_morph
  use PlantSoilChemistryAPIData, only: plt_soilchem
  use PlantBiomassAPIData,       only: plt_biom
  use PlantBGCRatesAPIData,      only: plt_bgcr
  use PlantRootBGCAPIData,       only: plt_rbgc
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
  public :: ReceivePlantRootAxes
  public :: SendPlantRootAxes
contains

  subroutine ReceivePlantRootAxes(NY,NX,NZ)
  use EcoSIMConfig, only : jsken=>jskenc,jcplx=>jcplxc
  integer, intent(in) :: NY,NX,NZ
  integer :: NR,K,L,M,N,NE

  associate(YoungNumber=>plt_morph%RootMediumYoungXNum_rpvr, &
    YoungLength=>plt_morph%RootMediumYoungLength_rpvr,YoungElements=>plt_morph%RootMediumYoungElms_rpvr)
    RootMediumYoungXNum_rpvr(1:NK_col(NY,NX),:,NZ,NY,NX) = YoungNumber(1:NK_col(NY,NX),:,NZ)
    RootMediumYoungLength_rpvr(1:NK_col(NY,NX),:,NZ,NY,NX) = YoungLength(1:NK_col(NY,NX),:,NZ)
    RootMediumYoungElms_rpvr(:,1:NK_col(NY,NX),:,NZ,NY,NX) = YoungElements(:,1:NK_col(NY,NX),:,NZ)
  end associate

    DO NR=1,pltpar%MaxNumRootAxes
      NRoot1stTipLay_raxes(NR,NZ,NY,NX) = plt_morph%NRoot1stTipLay_raxes(NR,NZ)
      Root1stDepz_raxes(NR,NZ,NY,NX)    = plt_morph%Root1stDepz_raxes(NR,NZ)
      RootMyco1stElm_raxs(1:NumPlantChemElms,NR,NZ,NY,NX) = plt_biom%RootMyco1stElm_raxs(1:NumPlantChemElms,NR,NZ)

      DO L=1,NK_col(NY,NX)
        RootMediumXNum_rpvr(L,NR,NZ,NY,NX)                          = plt_morph%RootMediumXNum_rpvr(L,NR,NZ)
        RootMediumRadius_rpvr(L,NR,NZ,NY,NX)                        = plt_morph%RootMediumRadius_rpvr(L,NR,NZ)
        RootMediumLength_rpvr(L,NR,NZ,NY,NX)                        = plt_morph%RootMediumLength_rpvr(L,NR,NZ)
        fctyok_scalar_rpvr(L,NR,NZ,NY,NX)                           = plt_morph%fctyok_scalar_rpvr(L,NR,NZ)
        RootCRRadius0_rpvr(L,NR,NZ,NY,NX)                           = plt_morph%RootCRRadius0_rpvr(L,NR,NZ)
        Root1stRadius_rpvr(L,NR,NZ,NY,NX)                           = plt_morph%Root1stRadius_rpvr(L,NR,NZ)
        RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,L,NR,NZ,NY,NX) = plt_biom%RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,L,NR,NZ)
        Root1stActStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ,NY,NX) = plt_biom%Root1stActStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ)
        Root1stLigStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ,NY,NX) = plt_biom%Root1stLigStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ)
        RootMediumStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ,NY,NX) = plt_biom%RootMediumStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ)

        Root1stLenPP_rpvr(L,NR,NZ,NY,NX)     = plt_morph%Root1stLenPP_rpvr(L,NR,NZ)
        RootAge_rpvr(L,NR,NZ,NY,NX)          = plt_morph%RootAge_rpvr(L,NR,NZ)
        RootMyco1stSinkC_rpvr(L,NR,NZ,NY,NX) = plt_rbgc%RootMyco1stSinkC_rpvr(L,NR,NZ)
        Cytokinin1stConc_rpvr(L,NR,NZ,NY,NX) = plt_rbgc%Cytokinin1stConc_rpvr(L,NR,NZ)
        CytokininMRConc_rpvr(L,NR,NZ,NY,NX)  = plt_rbgc%CytokininMRConc_rpvr(L,NR,NZ)
        CRootLumenArea_rpvr(L,NR,NZ,NY,NX)   = plt_morph%CRootLumenArea_rpvr(L,NR,NZ)
        DO N=1,Myco_pft(NZ,NY,NX)
          RootMyco2ndSinkC_rpvr(N,L,NR,NZ,NY,NX)  = plt_rbgc%RootMyco2ndSinkC_rpvr(N,L,NR,NZ)
          RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ,NY,NX) = plt_biom%RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ)
          Root2ndLen_rpvr(N,L,NR,NZ,NY,NX)                              = plt_morph%Root2ndLen_rpvr(N,L,NR,NZ)
          Root2ndXNum_rpvr(N,L,NR,NZ,NY,NX)                             = plt_morph%Root2ndXNum_rpvr(N,L,NR,NZ)
          Cytokinin2ndConc_rpvr(N,L,NR,NZ,NY,NX) = plt_rbgc%Cytokinin2ndConc_rpvr(N,L,NR,NZ)
        ENDDO
      ENDDO
    ENDDO

    DO L=NU_col(NY,NX),MaxSoilLays4Root_pft(NZ,NY,NX)
      DO K=1,jcplx
        DO N=1,Myco_pft(NZ,NY,NX)
          DO NE=1,NumPlantChemElms
            Soil2RootMycoExudE_pvr(NE,N,K,L,NZ,NY,NX)=plt_rbgc%Soil2RootMycoExudE_pvr(NE,N,K,L,NZ)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO M=1,jsken
      DO N=0,pltpar%NumLitterGroups
        DO NE=1,NumPlantChemElms
          PlantElmAllocMat4Litr(NE,N,M,NZ,NY,NX)=plt_soilchem%PlantElmAllocMat4Litr(NE,N,M,NZ)
        enddo
      enddo
    ENDDO

    !the following needs to be updated
    Root1stMaxRadius_pft(2,NZ,NY,NX) = plt_morph%Root1stMaxRadius_pft(2,NZ)
    Root2ndMaxRadius_pft(2,NZ,NY,NX) = plt_morph%Root2ndMaxRadius_pft(2,NZ)
    RootPorosity_pft(2,NZ,NY,NX)     = plt_morph%RootPorosity_pft(2,NZ)
    VmaxNH4Root_pft(2,NZ,NY,NX)      = plt_rbgc%VmaxNH4Root_pft(2,NZ)
    KmNH4Root_pft(2,NZ,NY,NX)        = plt_rbgc%KmNH4Root_pft(2,NZ)
    CMinNH4Root_pft(2,NZ,NY,NX)      = plt_rbgc%CMinNH4Root_pft(2,NZ)
    VmaxNO3Root_pft(2,NZ,NY,NX)      = plt_rbgc%VmaxNO3Root_pft(2,NZ)
    KmNO3Root_pft(2,NZ,NY,NX)        = plt_rbgc%KmNO3Root_pft(2,NZ)
    CminNO3Root_pft(2,NZ,NY,NX)      = plt_rbgc%CminNO3Root_pft(2,NZ)
    VmaxPO4Root_pft(2,NZ,NY,NX)      = plt_rbgc%VmaxPO4Root_pft(2,NZ)
    KmPO4Root_pft(2,NZ,NY,NX)        = plt_rbgc%KmPO4Root_pft(2,NZ)
    CMinPO4Root_pft(2,NZ,NY,NX)      = plt_rbgc%CMinPO4Root_pft(2,NZ)
    RootRadialResist_pft(2,NZ,NY,NX) = plt_morph%RootRadialResist_pft(2,NZ)
    Root2ndAxialResist_pft(2,NZ,NY,NX)  = plt_morph%Root2ndAxialResist_pft(2,NZ)
    !vessel axial resistances are recomputed at planting in StartPlants
    Root2ndAxialResist_pft(ipltroot,NZ,NY,NX) = plt_morph%Root2ndAxialResist_pft(ipltroot,NZ)
    RootSingleVesselRstaxial_pft(NZ,NY,NX)    = plt_morph%RootSingleVesselRstaxial_pft(NZ)
    RootSingleVesselArea_pft(NZ,NY,NX)        = plt_morph%RootSingleVesselArea_pft(NZ)
    CRootActVolPerMassC_pft(NZ,NY,NX)= plt_morph%CRootActVolPerMassC_pft(NZ)
    DO N=1,Myco_pft(NZ,NY,NX)
      RootMycoNonstElms_pft(1:NumPlantChemElms,N,NZ,NY,NX) = plt_biom%RootMycoNonstElms_pft(1:NumPlantChemElms,N,NZ)
      RootPoreTortu4Gas_pft(N,NZ,NY,NX)                     = plt_morph%RootPoreTortu4Gas_pft(N,NZ)
      RootRaidus_rpft(N,NZ,NY,NX)                          = plt_morph%RootRaidus_rpft(N,NZ)
      FineRootVolPerMassC_pft(N,NZ,NY,NX)                  = plt_morph%FineRootVolPerMassC_pft(N,NZ)
      Root1stSpecLen_pft(N,NZ,NY,NX)                       = plt_morph%Root1stSpecLen_pft(N,NZ)
      Root2ndSpecLen_pft(N,NZ,NY,NX)                       = plt_morph%Root2ndSpecLen_pft(N,NZ)
      Root1stMaxRadius1_pft(N,NZ,NY,NX)                    = plt_morph%Root1stMaxRadius1_pft(N,NZ)
      Root2ndMaxRadius1_pft(N,NZ,NY,NX)                    = plt_morph%Root2ndMaxRadius1_pft(N,NZ)
      Root1stXSecArea_pft(N,NZ,NY,NX)                      = plt_morph%Root1stXSecArea_pft(N,NZ)
      Root2ndXSecArea_pft(N,NZ,NY,NX)                      = plt_morph%Root2ndXSecArea_pft(N,NZ)
    enDDO
    ! Rebuild the native layer totals after receiving the authoritative axes.
    RootMyco2ndStrutElms_pvr(:,:,:,NZ,NY,NX)=sum( &
      RootMyco2ndStrutElms_rpvr(:,:,:,1:NumStructuralRootAxes_pft(NZ,NY,NX),NZ,NY,NX),dim=4)
  end subroutine ReceivePlantRootAxes

  subroutine SendPlantRootAxes(NY,NX,NZ)
  use EcoSIMConfig, only : jsken=>jskenc,jcplx=>jcplxc
  integer, intent(in) :: NY,NX,NZ
  integer :: NR,K,L,M,N,NE

    ! Native axis pools may have moved with the soil layers since the last step.
    RootMyco2ndStrutElms_pvr(:,:,:,NZ,NY,NX)=sum( &
      RootMyco2ndStrutElms_rpvr(:,:,:,1:NumStructuralRootAxes_pft(NZ,NY,NX),NZ,NY,NX),dim=4)
    plt_biom%RootMyco2ndStrutElms_pvr(:,:,:,NZ)=0._r8
    plt_biom%RootMyco2ndStrutElms_pvr(:,:,1:NL_col(NY,NX),NZ)= &
      RootMyco2ndStrutElms_pvr(:,:,1:NL_col(NY,NX),NZ,NY,NX)
  associate(YoungNumber=>plt_morph%RootMediumYoungXNum_rpvr, &
    YoungLength=>plt_morph%RootMediumYoungLength_rpvr,YoungElements=>plt_morph%RootMediumYoungElms_rpvr)
    YoungNumber(1:NK_col(NY,NX),:,NZ) = RootMediumYoungXNum_rpvr(1:NK_col(NY,NX),:,NZ,NY,NX)
    YoungLength(1:NK_col(NY,NX),:,NZ) = RootMediumYoungLength_rpvr(1:NK_col(NY,NX),:,NZ,NY,NX)
    YoungElements(:,1:NK_col(NY,NX),:,NZ) = RootMediumYoungElms_rpvr(:,1:NK_col(NY,NX),:,NZ,NY,NX)
  end associate

    DO NR=1,pltpar%MaxNumRootAxes
      plt_morph%NRoot1stTipLay_raxes(NR,NZ)=NRoot1stTipLay_raxes(NR,NZ,NY,NX)
      DO L=1,NK_col(NY,NX)
        plt_morph%RootCRRadius0_rpvr(L,NR,NZ) = RootCRRadius0_rpvr(L,NR,NZ,NY,NX)
        plt_morph%Root1stRadius_rpvr(L,NR,NZ) = Root1stRadius_rpvr(L,NR,NZ,NY,NX)
        plt_morph%RootMediumRadius_rpvr(L,NR,NZ) = RootMediumRadius_rpvr(L,NR,NZ,NY,NX)
        plt_morph%RootMediumLength_rpvr(L,NR,NZ) = RootMediumLength_rpvr(L,NR,NZ,NY,NX)
        plt_biom%RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,L,NR,NZ) = RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,L,NR,NZ,NY,NX)
        plt_biom%Root1stActStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ) = Root1stActStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ,NY,NX)
        plt_biom%Root1stLigStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ) = Root1stLigStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ,NY,NX)
        plt_biom%RootMediumStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ) = RootMediumStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ,NY,NX)

        plt_morph%Root1stLenPP_rpvr(L,NR,NZ)    = Root1stLenPP_rpvr(L,NR,NZ,NY,NX)
        plt_morph%RootAge_rpvr(L,NR,NZ)         = RootAge_rpvr(L,NR,NZ,NY,NX)
        plt_rbgc%RootMyco1stSinkC_rpvr(L,NR,NZ) = RootMyco1stSinkC_rpvr(L,NR,NZ,NY,NX)
        plt_rbgc%Cytokinin1stConc_rpvr(L,NR,NZ) = AZMAX1(Cytokinin1stConc_rpvr(L,NR,NZ,NY,NX))
        plt_rbgc%CytokininMRConc_rpvr(L,NR,NZ)  = CytokininMRConc_rpvr(L,NR,NZ,NY,NX)
        plt_morph%CRootLumenArea_rpvr(L,NR,NZ)   = CRootLumenArea_rpvr(L,NR,NZ,NY,NX)
        plt_morph%RootMediumXNum_rpvr(L,NR,NZ) = RootMediumXNum_rpvr(L,NR,NZ,NY,NX)
        DO N=1,Myco_pft(NZ,NY,NX)
          plt_morph%Root2ndLen_rpvr(N,L,NR,NZ)      = Root2ndLen_rpvr(N,L,NR,NZ,NY,NX)
          plt_morph%Root2ndXNum_rpvr(N,L,NR,NZ)     = Root2ndXNum_rpvr(N,L,NR,NZ,NY,NX)
          plt_rbgc%Cytokinin2ndConc_rpvr(N,L,NR,NZ)= Cytokinin2ndConc_rpvr(N,L,NR,NZ,NY,NX)
          plt_rbgc%RootMyco2ndSinkC_rpvr(N,L,NR,NZ) = RootMyco2ndSinkC_rpvr(N,L,NR,NZ,NY,NX)
          plt_biom%RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = &
            RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ,NY,NX)
        enddo
      enddo
      plt_morph%Root1stDepz_raxes(NR,NZ)    = Root1stDepz_raxes(NR,NZ,NY,NX)
      plt_biom%RootMyco1stElm_raxs(1:NumPlantChemElms,NR,NZ) = RootMyco1stElm_raxs(1:NumPlantChemElms,NR,NZ,NY,NX)
    enddo

    DO M=1,jsken
      DO NE=1,NumPlantChemElms
        plt_biom%StandDeadCompKElms_pft(NE,M,NZ)=StandDeadCompKElms_pft(NE,M,NZ,NY,NX)
      ENDDO
    ENDDO
!!!!  LitrfallElms_pvr in restart file?
    DO L=0,NK_col(NY,NX)
      DO K=1,micpar%NumOfPlantLitrCmplxs
        DO M=1,jsken
          plt_bgcr%LitrfallElms_pvr(1:NumPlantChemElms,M,K,L,NZ)=LitrfallElms_pvr(1:NumPlantChemElms,M,K,L,NZ,NY,NX)
        enddo
      enddo
    ENDDO

    DO M=1,jsken
      DO N=0,pltpar%NumLitterGroups
        DO NE=1,NumPlantChemElms
          plt_soilchem%PlantElmAllocMat4Litr(NE,N,M,NZ)=PlantElmAllocMat4Litr(NE,N,M,NZ,NY,NX)
        enddo
      enddo
    ENDDO

!!!
  end subroutine SendPlantRootAxes
end module PlantRootAxisTransferMod
