module PlantRootTransferMod
  ! Ordered root transfers used by PlantAPISend/PlantAPIRecv.
  use data_kind_mod,             only: r8 => DAT_KIND_R8, yearIJ_type
  use EcoSiMParDataMod,          only: micpar,            pltpar
  use SoilPhysDataType,          only: SurfAlbedo_col,    SoilSurfDepZ_col
  use MiniMathMod,               only: AZMAX1,            safe_adb
  use DebugToolMod,              only: PrintInfo
  use PlantMorphologyAPIData,    only: plt_morph
  use PlantPhenologyAPIData,     only: plt_pheno
  use PlantSoilChemistryAPIData, only: plt_soilchem
  use PlantBiomassAPIData,       only: plt_biom
  use PlantEnergyWaterAPIData,   only: plt_ew
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
  public :: ReceivePlantRoots
  public :: SendPlantRoots
contains

  subroutine ReceivePlantRoots(NY,NX,NZ)
  use EcoSIMConfig, only : jsken=>jskenc,jcplx=>jcplxc
  integer, intent(in) :: NY,NX,NZ
  integer :: K,L,N,NE,idg

    DO  L=NUI_col(NY,NX),NK_col(NY,NX)
      GroSrcRootStress_pvr(L,NZ,NY,NX) = plt_rbgc%GroSrcRootStress_pvr(L,NZ)
      DO K=1,jcplx
        DOM_MicP_vr(idom_doc:idom_dop,K,L,NY,NX)=plt_soilchem%DOM_MicP_vr(idom_doc:idom_dop,K,L)
        DOM_MicP_drib_vr(idom_doc:idom_dop,K,L,NY,NX)=plt_soilchem%DOM_MicP_drib_vr(idom_doc:idom_dop,K,L)
      ENDDO
      RootMediumXNum_pvr(L,NZ,NY,NX) = plt_morph%RootMediumXNum_pvr(L,NZ)
      Root1stXNumL_pvr(L,NZ,NY,NX)    = plt_morph%Root1stXNumL_pvr(L,NZ)
      DO NE=1,NumPlantChemElms
        RootShootExch_pvr(NE,L,NZ,NY,NX) = plt_bgcr%RootShootExch_pvr(NE,L,NZ)
      ENDDO
      DO N = 1, Myco_pft(NZ,NY,NX)
        ROOTNLim_rpvr(N,L,NZ,NY,NX)                                = plt_biom%ROOTNLim_rpvr(N,L,NZ)
        ROOTPLim_rpvr(N,L,NZ,NY,NX)                                = plt_biom%ROOTPLim_rpvr(N,L,NZ)
        RootMaintDef_CO2_pvr(N,L,NZ,NY,NX)                         = plt_bgcr%RootMaintDef_CO2_pvr(N,L,NZ)
        Nutruptk_fClim_rpvr(N,L,NZ,NY,NX)                          = plt_bgcr%Nutruptk_fClim_rpvr(N,L,NZ)
        Nutruptk_fNlim_rpvr(N,L,NZ,NY,NX)                          = plt_bgcr%Nutruptk_fNlim_rpvr(N,L,NZ)
        Nutruptk_fPlim_rpvr(N,L,NZ,NY,NX)                          = plt_bgcr%Nutruptk_fPlim_rpvr(N,L,NZ)
        Nutruptk_fProtC_rpvr(N,L,NZ,NY,NX)                         = plt_bgcr%Nutruptk_fProtC_rpvr(N,L,NZ)
        fRootGrowPSISense_pvr(N,L,NZ,NY,NX)                        = plt_pheno%fRootGrowPSISense_pvr(N,L,NZ)
        RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ,NY,NX)    = plt_biom%RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ)
        RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ,NY,NX) = plt_biom%RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ)
        RootProteinConc_rpvr(N,L,NZ,NY,NX)                         = plt_biom%RootProteinConc_rpvr(N,L,NZ)
        trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ,NY,NX)              = plt_rbgc%trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)
        trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ,NY,NX)              = plt_rbgc%trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)
        PSIRoot_pvr(N,L,NZ,NY,NX)                                  = plt_ew%PSIRoot_pvr(N,L,NZ)
        PSIRootOSMO_vr(N,L,NZ,NY,NX)                               = plt_ew%PSIRootOSMO_vr(N,L,NZ)
        PSIRootTurg_vr(N,L,NZ,NY,NX)                               = plt_ew%PSIRootTurg_vr(N,L,NZ)
        Root2ndXNumL_rpvr(N,L,NZ,NY,NX)                            = plt_morph%Root2ndXNumL_rpvr(N,L,NZ)
        RootTotLenPerPlant_pvr(N,L,NZ,NY,NX)                       = plt_morph%RootTotLenPerPlant_pvr(N,L,NZ)
        RootAbsorbLenPerPlant_pvr(N,L,NZ,NY,NX)                    = plt_morph%RootAbsorbLenPerPlant_pvr(N,L,NZ)
        RootLenPerPlant_pvr(N,L,NZ,NY,NX)                          = plt_morph%RootLenPerPlant_pvr(N,L,NZ)
        RootLenDensPerPlant_pvr(N,L,NZ,NY,NX)                      = plt_morph%RootLenDensPerPlant_pvr(N,L,NZ)
        RootPoreVol_pvr(N,L,NZ,NY,NX)                             = plt_morph%RootPoreVol_pvr(N,L,NZ)
        RootVH2O_pvr(N,L,NZ,NY,NX)                                 = plt_morph%RootVH2O_pvr(N,L,NZ)
        Root1stRadius_pvr(N,L,NZ,NY,NX)                            = plt_morph%Root1stRadius_pvr(N,L,NZ)
        Root2ndRadius_rpvr(N,L,NZ,NY,NX)                           = plt_morph%Root2ndRadius_rpvr(N,L,NZ)
        RootSAreaPerPlant_pvr(N,L,NZ,NY,NX)                         = plt_morph%RootSAreaPerPlant_pvr(N,L,NZ)
        RootArea1stPP_pvr(N,L,NZ,NY,NX)                            = plt_morph%RootArea1stPP_pvr(N,L,NZ)
        RootArea2ndPP_pvr(N,L,NZ,NY,NX)                            = plt_morph%RootArea2ndPP_pvr(N,L,NZ)
        Root2ndEffLen4uptk_rpvr(N,L,NZ,NY,NX)                      = plt_morph%Root2ndEffLen4uptk_rpvr(N,L,NZ)
        RootRespPotent_pvr(N,L,NZ,NY,NX)                           = plt_rbgc%RootRespPotent_pvr(N,L,NZ)
        RootCO2EmisPot_pvr(N,L,NZ,NY,NX)                           = plt_rbgc%RootCO2EmisPot_pvr(N,L,NZ)
        RootCO2Autor_pvr(N,L,NZ,NY,NX)                             = plt_rbgc%RootCO2Autor_pvr(N,L,NZ)
        RCO2Emis2Root_rpvr(N,L,NZ,NY,NX)                            = plt_rbgc%RCO2Emis2Root_rpvr(N,L,NZ)
        RootO2Uptk_pvr(N,L,NZ,NY,NX)                               = plt_rbgc%RootO2Uptk_pvr(N,L,NZ)
        RootAtmGasConductance_rpvr(idg_beg:idg_NH3,N,L,NZ,NY,NX)   = plt_rbgc%RootAtmGasConductance_rpvr(idg_beg:idg_NH3,N,L,NZ)
        RootUptkSoiSol_pvr(idg_CO2,N,L,NZ,NY,NX)                   = plt_rbgc%RootUptkSoiSol_pvr(idg_CO2,N,L,NZ)
        RootUptkSoiSol_pvr(idg_O2,N,L,NZ,NY,NX)                    = plt_rbgc%RootUptkSoiSol_pvr(idg_O2,N,L,NZ)
        RootUptkSoiSol_pvr(idg_CH4,N,L,NZ,NY,NX)                   = plt_rbgc%RootUptkSoiSol_pvr(idg_CH4,N,L,NZ)
        RootUptkSoiSol_pvr(idg_N2O,N,L,NZ,NY,NX)                   = plt_rbgc%RootUptkSoiSol_pvr(idg_N2O,N,L,NZ)
        RootUptkSoiSol_pvr(idg_NH3,N,L,NZ,NY,NX)                   = plt_rbgc%RootUptkSoiSol_pvr(idg_NH3,N,L,NZ)
        RootUptkSoiSol_pvr(idg_NH3B,N,L,NZ,NY,NX)                  = plt_rbgc%RootUptkSoiSol_pvr(idg_NH3B,N,L,NZ)
        RootUptkSoiSol_pvr(idg_H2,N,L,NZ,NY,NX)                    = plt_rbgc%RootUptkSoiSol_pvr(idg_H2,N,L,NZ)
        trcg_air2root_flx_pvr(idg_CO2,N,L,NZ,NY,NX)                = plt_rbgc%trcg_air2root_flx_pvr(idg_CO2,N,L,NZ)
        trcg_air2root_flx_pvr(idg_O2,N,L,NZ,NY,NX)                 = plt_rbgc%trcg_air2root_flx_pvr(idg_O2,N,L,NZ)
        trcg_air2root_flx_pvr(idg_CH4,N,L,NZ,NY,NX)                = plt_rbgc%trcg_air2root_flx_pvr(idg_CH4,N,L,NZ)
        trcg_air2root_flx_pvr(idg_N2O,N,L,NZ,NY,NX)                = plt_rbgc%trcg_air2root_flx_pvr(idg_N2O,N,L,NZ)
        trcg_air2root_flx_pvr(idg_NH3,N,L,NZ,NY,NX)                = plt_rbgc%trcg_air2root_flx_pvr(idg_NH3,N,L,NZ)
        trcg_air2root_flx_pvr(idg_H2,N,L,NZ,NY,NX)                 = plt_rbgc%trcg_air2root_flx_pvr(idg_H2,N,L,NZ)
        trcg_Root_gas2aqu_flx_vr(idg_CO2,N,L,NZ,NY,NX)             = plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_CO2,N,L,NZ)
        trcg_Root_gas2aqu_flx_vr(idg_O2,N,L,NZ,NY,NX)              = plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_O2,N,L,NZ)
        trcg_Root_gas2aqu_flx_vr(idg_CH4,N,L,NZ,NY,NX)             = plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_CH4,N,L,NZ)
        trcg_Root_gas2aqu_flx_vr(idg_N2O,N,L,NZ,NY,NX)             = plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_N2O,N,L,NZ)
        trcg_Root_gas2aqu_flx_vr(idg_NH3,N,L,NZ,NY,NX)             = plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_NH3,N,L,NZ)
        trcg_Root_gas2aqu_flx_vr(idg_H2,N,L,NZ,NY,NX)              = plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_H2,N,L,NZ)
        RootNH4DmndSoil_pvr(N,L,NZ,NY,NX)                          = plt_rbgc%RootNH4DmndSoil_pvr(N,L,NZ)
        VmaxNH4Root_pvr(N,L,NZ,NY,NX)                              = plt_rbgc%VmaxNH4Root_pvr(N,L,NZ)
        VmaxNO3Root_pvr(N,L,NZ,NY,NX)                              = plt_rbgc%VmaxNO3Root_pvr(N,L,NZ)
        RootRadialKond2H2O_pvr(N,L,NZ,NY,NX)                       = plt_ew%RootRadialKond2H2O_pvr(N,L,NZ)
        RootAXialKond2H2O_pvr(N,L,NZ,NY,NX)                        = plt_ew%RootAXialKond2H2O_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_NH4,N,L,NZ,NY,NX)                    = plt_rbgc%RootNutUptake_pvr(ids_NH4,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_NH4,N,L,NZ,NY,NX)                = plt_rbgc%RootOUlmNutUptake_pvr(ids_NH4,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_NH4,N,L,NZ,NY,NX)                = plt_rbgc%RootCUlmNutUptake_pvr(ids_NH4,N,L,NZ)
        RootNH4DmndBand_pvr(N,L,NZ,NY,NX)                          = plt_rbgc%RootNH4DmndBand_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_NH4B,N,L,NZ,NY,NX)                   = plt_rbgc%RootNutUptake_pvr(ids_NH4B,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_NH4B,N,L,NZ,NY,NX)               = plt_rbgc%RootOUlmNutUptake_pvr(ids_NH4B,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_NH4B,N,L,NZ,NY,NX)               = plt_rbgc%RootCUlmNutUptake_pvr(ids_NH4B,N,L,NZ)
        RootNO3DmndSoil_pvr(N,L,NZ,NY,NX)                          = plt_rbgc%RootNO3DmndSoil_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_NO3,N,L,NZ,NY,NX)                    = plt_rbgc%RootNutUptake_pvr(ids_NO3,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_NO3,N,L,NZ,NY,NX)                = plt_rbgc%RootOUlmNutUptake_pvr(ids_NO3,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_NO3,N,L,NZ,NY,NX)                = plt_rbgc%RootCUlmNutUptake_pvr(ids_NO3,N,L,NZ)
        RootNO3DmndBand_pvr(N,L,NZ,NY,NX)                          = plt_rbgc%RootNO3DmndBand_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_NO3B,N,L,NZ,NY,NX)                   = plt_rbgc%RootNutUptake_pvr(ids_NO3B,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_NO3B,N,L,NZ,NY,NX)               = plt_rbgc%RootOUlmNutUptake_pvr(ids_NO3B,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_NO3B,N,L,NZ,NY,NX)               = plt_rbgc%RootCUlmNutUptake_pvr(ids_NO3B,N,L,NZ)
        RootH2PO4DmndSoil_pvr(N,L,NZ,NY,NX)                        = plt_rbgc%RootH2PO4DmndSoil_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_H2PO4,N,L,NZ,NY,NX)                  = plt_rbgc%RootNutUptake_pvr(ids_H2PO4,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_H2PO4,N,L,NZ,NY,NX)              = plt_rbgc%RootOUlmNutUptake_pvr(ids_H2PO4,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_H2PO4,N,L,NZ,NY,NX)              = plt_rbgc%RootCUlmNutUptake_pvr(ids_H2PO4,N,L,NZ)
        RootH2PO4DmndBand_pvr(N,L,NZ,NY,NX)                        = plt_rbgc%RootH2PO4DmndBand_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_H2PO4B,N,L,NZ,NY,NX)                 = plt_rbgc%RootNutUptake_pvr(ids_H2PO4B,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ,NY,NX)             = plt_rbgc%RootOUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ,NY,NX)             = plt_rbgc%RootCUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ)
        RootH1PO4DmndSoil_pvr(N,L,NZ,NY,NX)                        = plt_rbgc%RootH1PO4DmndSoil_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_H1PO4,N,L,NZ,NY,NX)                  = plt_rbgc%RootNutUptake_pvr(ids_H1PO4,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_H1PO4,N,L,NZ,NY,NX)              = plt_rbgc%RootOUlmNutUptake_pvr(ids_H1PO4,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_H1PO4,N,L,NZ,NY,NX)              = plt_rbgc%RootCUlmNutUptake_pvr(ids_H1PO4,N,L,NZ)
        RootH1PO4DmndBand_pvr(N,L,NZ,NY,NX)                        = plt_rbgc%RootH1PO4DmndBand_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_H1PO4B,N,L,NZ,NY,NX)                 = plt_rbgc%RootNutUptake_pvr(ids_H1PO4B,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ,NY,NX)             = plt_rbgc%RootOUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ,NY,NX)             = plt_rbgc%RootCUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ)
        RootO2Dmnd4Resp_pvr(N,L,NZ,NY,NX)                          = plt_rbgc%RootO2Dmnd4Resp_pvr(N,L,NZ)
        RPlantRootH2OUptk_pvr(N,L,NZ,NY,NX)                        = plt_ew%RPlantRootH2OUptk_pvr(N,L,NZ)
        RootH2OUptkStress_pvr(N,L,NZ,NY,NX)                        = plt_ew%RootH2OUptkStress_pvr(N,L,NZ)
        RootMycoActiveBiomC_pvr(N,L,NZ,NY,NX)                      = plt_biom%RootMycoActiveBiomC_pvr(N,L,NZ)
        Root1stTransptArea_pvr(N,L,NZ,NY,NX)                       = plt_morph%Root1stTransptArea_pvr(N,L,NZ)
        RootMedTransptArea_pvr(N,L,NZ,NY,NX)                       = plt_morph%RootMedTransptArea_pvr(N,L,NZ)
        PopuRootMycoC_pvr(N,L,NZ,NY,NX)                            = AZMAX1(plt_biom%PopuRootMycoC_pvr(N,L,NZ))
        RootResist4H2O_pvr(N,L,NZ,NY,NX)                           = plt_ew%RootResist4H2O_pvr(N,L,NZ)
        RootProteinC_pvr(N,L,NZ,NY,NX)                             = plt_biom%RootProteinC_pvr(N,L,NZ)
        RAutoRootO2Limter_rpvr(N,L,NZ,NY,NX)                       = plt_rbgc%RAutoRootO2Limter_rpvr(N,L,NZ)
        RootCO2Autor_vr(L,NY,NX)                                   = RootCO2Autor_vr(L,NY,NX)+RootCO2Autor_pvr(N,L,NZ,NY,NX)
      ENDDO
      SapFlowVlinear_pvr(L,NZ,NY,NX) = plt_ew%SapFlowVlinear_pvr(L,NZ)
      CRootLumenArea_pvr(L,NZ,NY,NX) = plt_morph%CRootLumenArea_pvr(L,NZ)
      MRootLumenArea_pvr(L,NZ,NY,NX) = plt_morph%MRootLumenArea_pvr(L,NZ)
      RootCO2Ar2Soil_vr(L,NY,NX)     = RootCO2Ar2Soil_vr(L,NY,NX)+plt_rbgc%RootCO2Ar2Soil_pvr(L,NZ)
      RootCO2Ar2Root_vr(L,NY,NX)     = RootCO2Ar2Root_vr(L,NY,NX)+plt_rbgc%RootCO2Ar2RootX_pvr(L,NZ)
      do idg=idg_beg,idg_NH3
        trcs_deadroot2soil_vr(idg,L,NY,NX)    = trcs_deadroot2soil_vr(idg,L,NY,NX) + plt_rbgc%trcs_deadroot2soil_pvr(idg,L,NZ)
      ENDDO
    ENDDO

  end subroutine ReceivePlantRoots

  subroutine SendPlantRoots(NY,NX,NZ)
  use EcoSIMConfig, only : jsken=>jskenc,jcplx=>jcplxc
  integer, intent(in) :: NY,NX,NZ
  integer :: L,N,NE

    ! Unrooted layers retain these caches until roots regrow into them.
    plt_ew%RootRadialKond2H2O_pvr(:,:,NZ)=RootRadialKond2H2O_pvr(:,:,NZ,NY,NX)
    plt_ew%RootAxialKond2H2O_pvr(:,:,NZ)=RootAxialKond2H2O_pvr(:,:,NZ,NY,NX)
    plt_rbgc%RootAtmGasConductance_rpvr(:,:,:,NZ)=RootAtmGasConductance_rpvr(:,:,:,NZ,NY,NX)

    DO L=1,NL_col(NY,NX)
      DO NE=1,NumPlantChemElms
        plt_biom%RootNodulStrutElms_rpvr(NE,L,NZ) =RootNodulStrutElms_rpvr(NE,L,NZ,NY,NX)
      ENDDO
    ENDDO

    DO L=1,NK_col(NY,NX)
      plt_morph%RootMediumXNum_pvr(L,NZ)=RootMediumXNum_pvr(L,NZ,NY,NX)
      plt_morph%Root1stXNumL_pvr(L,NZ) = Root1stXNumL_pvr(L,NZ,NY,NX)
      plt_morph%CRootLumenArea_pvr(L,NZ)   = CRootLumenArea_pvr(L,NZ,NY,NX)
      plt_morph%MRootLumenArea_pvr(L,NZ)   = MRootLumenArea_pvr(L,NZ,NY,NX)
      DO N=1,Myco_pft(NZ,NY,NX)
        plt_biom%RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ) = RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ,NY,NX)
        plt_biom%RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ) = RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ,NY,NX)
        plt_biom%RootProteinConc_rpvr(N,L,NZ)                         = RootProteinConc_rpvr(N,L,NZ,NY,NX)

        plt_rbgc%trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)           = trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ,NY,NX)
        plt_rbgc%trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)           = trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ,NY,NX)

        plt_ew%PSIRoot_pvr(N,L,NZ)                = PSIRoot_pvr(N,L,NZ,NY,NX)
        plt_ew%PSIRootOSMO_vr(N,L,NZ)             = PSIRootOSMO_vr(N,L,NZ,NY,NX)
        plt_ew%PSIRootTurg_vr(N,L,NZ)             = PSIRootTurg_vr(N,L,NZ,NY,NX)
        plt_rbgc%RootRespPotent_pvr(N,L,NZ)       = RootRespPotent_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootCO2EmisPot_pvr(N,L,NZ)       = RootCO2EmisPot_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootCO2AutorX_pvr(N,L,NZ)        = RootCO2Autor_pvr(N,L,NZ,NY,NX)
        plt_morph%Root2ndXNumL_rpvr(N,L,NZ)       = Root2ndXNumL_rpvr(N,L,NZ,NY,NX)
        plt_morph%RootTotLenPerPlant_pvr(N,L,NZ)  = RootTotLenPerPlant_pvr(N,L,NZ,NY,NX)
        plt_morph%RootAbsorbLenPerPlant_pvr(N,L,NZ)=RootAbsorbLenPerPlant_pvr(N,L,NZ,NY,NX)
        plt_morph%RootLenDensPerPlant_pvr(N,L,NZ) = RootLenDensPerPlant_pvr(N,L,NZ,NY,NX)
        plt_morph%RootPoreVol_pvr(N,L,NZ)        = RootPoreVol_pvr(N,L,NZ,NY,NX)
        plt_morph%RootVH2O_pvr(N,L,NZ)            = RootVH2O_pvr(N,L,NZ,NY,NX)
        plt_morph%Root1stRadius_pvr(N,L,NZ)       = Root1stRadius_pvr(N,L,NZ,NY,NX)
        plt_morph%Root2ndRadius_rpvr(N,L,NZ)      = Root2ndRadius_rpvr(N,L,NZ,NY,NX)
        plt_morph%RootSAreaPerPlant_pvr(N,L,NZ)    = RootSAreaPerPlant_pvr(N,L,NZ,NY,NX)
        plt_morph%Root2ndEffLen4uptk_rpvr(N,L,NZ) = Root2ndEffLen4uptk_rpvr(N,L,NZ,NY,NX)

        plt_rbgc%RootO2Dmnd4Resp_pvr(N,L,NZ)       = RootO2Dmnd4Resp_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootNH4DmndSoilPrev_pvr(N,L,NZ)   = RootNH4DmndSoil_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootNH4DmndBandPrev_pvr(N,L,NZ)   = RootNH4DmndBand_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootNO3DmndSoilPrev_pvr(N,L,NZ)   = RootNO3DmndSoil_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootNO3DmndBandPrev_pvr(N,L,NZ)   = RootNO3DmndBand_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootH2PO4DmndSoilPrev_pvr(N,L,NZ) = RootH2PO4DmndSoil_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootH2PO4DmndBandPrev_pvr(N,L,NZ) = RootH2PO4DmndBand_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootH1PO4DmndSoilPrev_pvr(N,L,NZ) = RootH1PO4DmndSoil_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootH1PO4DmndBandPrev_pvr(N,L,NZ) = RootH1PO4DmndBand_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RAutoRootO2Limter_rpvr(N,L,NZ)    = RAutoRootO2Limter_rpvr(N,L,NZ,NY,NX)
        plt_biom%RootMycoActiveBiomC_pvr(N,L,NZ)   = RootMycoActiveBiomC_pvr(N,L,NZ,NY,NX)
        plt_morph%Root1stTransptArea_pvr(N,L,NZ)   = Root1stTransptArea_pvr(N,L,NZ,NY,NX)
        plt_morph%RootMedTransptArea_pvr(N,L,NZ)   = RootMedTransptArea_pvr(N,L,NZ,NY,NX)
        plt_biom%PopuRootMycoC_pvr(N,L,NZ)         = PopuRootMycoC_pvr(N,L,NZ,NY,NX)
        plt_biom%RootProteinC_pvr(N,L,NZ)          = RootProteinC_pvr(N,L,NZ,NY,NX)

      enddo
      plt_biom%RootNodulNonstElms_rpvr(1:NumPlantChemElms,L,NZ)=RootNodulNonstElms_rpvr(1:NumPlantChemElms,L,NZ,NY,NX)
    ENDDO

    DO L=1,NumCanopyLayers
      plt_morph%CanopyStemSurfAreaZ_pft(L,NZ) = CanopyStemSurfAreaZ_pft(L,NZ,NY,NX)
      plt_morph%CanopyLeafAreaZ_pft(L,NZ) = CanopyLeafAreaZ_pft(L,NZ,NY,NX)
      plt_biom%CanopyLeafCLyr_pft(L,NZ)   = CanopyLeafCLyr_pft(L,NZ,NY,NX)
      plt_morph%CanopySurfAreaProfDead_pft(L,NZ) = CanopySurfAreaProfDead_pft(L,NZ,NY,NX)
    ENDDO
    plt_morph%CRootActVolPerMassC_pft(NZ) = CRootActVolPerMassC_pft(NZ,NY,NX)
    DO N=1,Myco_pft(NZ,NY,NX)
      plt_morph%FineRootVolPerMassC_pft(N,NZ)   = FineRootVolPerMassC_pft(N,NZ,NY,NX)
      plt_morph%RootPoreTortu4Gas_pft(N,NZ)     = RootPoreTortu4Gas_pft(N,NZ,NY,NX)
      plt_morph%Root2ndXSecArea_pft(N,NZ)   = Root2ndXSecArea_pft(N,NZ,NY,NX)
      plt_morph%Root1stXSecArea_pft(N,NZ)   = Root1stXSecArea_pft(N,NZ,NY,NX)
      plt_morph%Root1stMaxRadius1_pft(N,NZ) = Root1stMaxRadius1_pft(N,NZ,NY,NX)
      plt_morph%Root2ndMaxRadius1_pft(N,NZ) = Root2ndMaxRadius1_pft(N,NZ,NY,NX)
      plt_morph%Root1stSpecLen_pft(N,NZ)    = Root1stSpecLen_pft(N,NZ,NY,NX)
      plt_morph%RootRaidus_rpft(N,NZ)       = RootRaidus_rpft(N,NZ,NY,NX)
      plt_morph%Root2ndSpecLen_pft(N,NZ)    = Root2ndSpecLen_pft(N,NZ,NY,NX)
    ENDDO
  end subroutine SendPlantRoots
end module PlantRootTransferMod
