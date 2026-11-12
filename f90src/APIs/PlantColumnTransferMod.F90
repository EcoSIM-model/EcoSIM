module PlantColumnTransferMod
  ! Ordered column transfers used by PlantAPISend/PlantAPIRecv.
  use data_kind_mod,             only: r8 => DAT_KIND_R8, yearIJ_type
  use EcoSiMParDataMod,          only: micpar,            pltpar
  use SoilPhysDataType,          only: SurfAlbedo_col,    SoilSurfDepZ_col
  use MiniMathMod,               only: AZMAX1,            safe_adb, AZERO
  use DebugToolMod,              only: PrintInfo
  use PlantSiteAPIData,          only: plt_site
  use PlantRadiationAPIData,     only: plt_rad
  use PlantMorphologyAPIData,    only: plt_morph
  use PlantPhenologyAPIData,     only: plt_pheno
  use PlantSoilChemistryAPIData, only: plt_soilchem
  use PlantAllometryAPIData,     only: plt_allom
  use PlantBiomassAPIData,       only: plt_biom
  use PlantEnergyWaterAPIData,   only: plt_ew
  use PlantDisturbanceAPIData,   only: plt_distb
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
  public :: ReceivePlantColumns
  public :: SendPlantColumnInputs
  public :: SendPlantColumnState
contains

  subroutine ReceivePlantColumns(I1,NY,NX)
  use EcoSIMConfig, only : jsken=>jskenc,jcplx=>jcplxc
  integer, intent(in) :: I1,NY,NX
  integer :: K,L,M,NE,NN

  NumActivePlants_col(NY,NX)                          = sum(plt_pheno%IsPlantActive_pft(1:NP0_col(NY,NX)))
  PlantPopu_col(NY,NX)                                = plt_site%PlantPopu_col
  ECO_ER_col(NY,NX)                                   = plt_bgcr%ECO_ER_col
  Eco_NBP_CumYr_col(NY,NX)                            = plt_bgcr%Eco_NBP_CumYr_col
  Air_Heat_Latent_store_col(NY,NX)                    = plt_ew%Air_Heat_Latent_store_col
  Air_Heat_Sens_store_col(NY,NX)                      = plt_ew%Air_Heat_Sens_store_col
  Eco_AutoR_CumYr_col(NY,NX)                          = AZERO(plt_bgcr%Eco_AutoR_CumYr_col)
  LitrFallStrutElms_col(1:NumPlantChemElms,NY,NX)     = plt_bgcr%LitrFallStrutElms_col(1:NumPlantChemElms)
  EcoHavstElmnt_CumYr_col(1:NumPlantChemElms,NY,NX)   = plt_distb%EcoHavstElmnt_CumYr_col(1:NumPlantChemElms)
  WatHeldOnCanopy_col(NY,NX)                          = plt_ew%WatHeldOnCanopy_col
  SnowOnCanopy_col(NY,NX)                             = plt_ew%SnowOnCanopy_col
  Eco_Heat_Sens_col(NY,NX)                            = plt_ew%Eco_Heat_Sens_col
  StandingDeadStrutElms_col(1:NumPlantChemElms,NY,NX) = plt_biom%StandingDeadStrutElms_col(1:NumPlantChemElms)
  H2OLoss_CumYr_col(NY,NX)                            = plt_ew%H2OLoss_CumYr_col
  StemArea_col(NY,NX)                                 = plt_morph%StemArea_col
  CanopyHeatLoss2Dist_col(NY,NX)                          = plt_ew%CanopyHeatLoss2Dist_col
  QCanopyWatLoss2Dist_col(NY,NX)                          = plt_ew%QCanopyWatLoss2Dist_col
  CanopyLeafArea_col(NY,NX)                           = plt_morph%CanopyLeafArea_col
  Eco_NetRad_col(NY,NX)                               = plt_rad%Eco_NetRad_col
  Eco_Heat_Latent_col(NY,NX)                          = plt_ew%Eco_Heat_Latent_col
  Eco_Heat_GrndSurf_col(NY,NX)                        = plt_ew%Eco_Heat_GrndSurf_col
  QVegET_col(NY,NX)                                   = plt_ew%QVegET_col
  LWRadCanG_col(NY,NX)                                = plt_ew%LWRadCanG
  VapXAir2Canopy_col(NY,NX)                           = plt_ew%VapXAir2Canopy_col
  HeatSensVapAir2Canopy_col(NY,NX)                           = plt_ew%HeatSensVapAir2Canopy_col
  HeatPrecAir2Canopy_col(NY,NX)                           = plt_ew%HeatPrecAir2Canopy_col
  HeatFlx2Canopy_col(NY,NX)                           = plt_ew%HeatFlx2Canopy_col
  CanopyBiomWater_col(NY,NX)                          = plt_ew%CanopyBiomWater_col
  CanopyHeatStor_col(NY,NX)                           = plt_ew%CanopyHeatStor_col
  TRootGasLossDisturb_col(idg_beg:idg_NH3,NY,NX)      = plt_rbgc%TRootGasLossDisturb_col(idg_beg:idg_NH3)
  Canopy_NEE_col(NY,NX)                               = plt_bgcr%Canopy_NEE_col
  TPlantRootH2OUptake_col(NY,NX)                      = plt_ew%TPlantRootH2OUptake_col
  FERT(ifert_plant_manuC:ifert_plant_manuP,I1,NY,NX)  = plt_distb%FERT(ifert_plant_manuC:ifert_plant_manuP)
  FERT(ifert_N_urea,I1,NY,NX)                         = plt_distb%FERT(ifert_N_urea)
  IYTYP(iAmendtyp_Manure,I1,NY,NX)                    = plt_distb%IYTYP
  FracWoodStalkElmAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs) = plt_allom%FracWoodStalkElmAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)
  FracRootElmAllocm(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)      = plt_allom%FracRootElmAllocm(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)
  FracLeafShethElmAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)    = plt_allom%FracLeafShethElmAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)
  FracPetolShethAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)   = plt_allom%FracPetolShethAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)
  QH2OLoss_lnds                                                         = plt_site%QH2OLoss_lnds

  DO L=1,NumCanopyLayers
    tCanLeafC_clyr(L,NY,NX)        = plt_biom%tCanLeafC_clyr(L)
    CanopyStemAareZ_col(L,NY,NX) = plt_morph%CanopyStemAareZ_col(L)
    CanopyLeafAareZ_col(L,NY,NX) = plt_morph%CanopyLeafAareZ_col(L)
  ENDDO

  DO L=NU_col(NY,NX),NL_col(NY,NX)

    DO K=1,jcplx
      DO NE=1,NumPlantChemElms
        REcoDOMProd_vr(NE,K,L,NY,NX)=plt_bgcr%REcoDOMProd_vr(NE,K,L)
      ENDDO
    ENDDO
    DO NN=1,NPH
      REcoUptkSoilO2M_vr(NN,L,NY,NX)=REcoUptkSoilO2M_vr(NN,L,NY,NX)+plt_rbgc%REcoUptkSoilO2M_vr(NN,L)
    ENDDO
  ENDDO

  DO L=0,NL_col(NY,NX)
    REcoH2PO4DmndBand_vr(L,NY,NX)  = plt_bgcr%REcoH2PO4DmndBand_vr(L)
    REcoH1PO4DmndBand_vr(L,NY,NX)  = plt_bgcr%REcoH1PO4DmndBand_vr(L)
    REcoNO3DmndBand_vr(L,NY,NX)    = plt_bgcr%REcoNO3DmndBand_vr(L)
    REcoNH4DmndBand_vr(L,NY,NX)    = plt_bgcr%REcoNH4DmndBand_vr(L)
    REcoH1PO4DmndSoil_vr(L,NY,NX)  = plt_bgcr%REcoH1PO4DmndSoil_vr(L)
    REcoH2PO4DmndSoil_vr(L,NY,NX)  = plt_bgcr%REcoH2PO4DmndSoil_vr(L)
    REcoNO3DmndSoil_vr(L,NY,NX)    = plt_bgcr%REcoNO3DmndSoil_vr(L)
    REcoNH4DmndSoil_vr(L,NY,NX)    = plt_bgcr%REcoNH4DmndSoil_vr(L)
    REcoO2DmndResp_vr(L,NY,NX)     = plt_bgcr%REcoO2DmndResp_vr(L)
    THeatLossRoot2Soil_vr(L,NY,NX) = plt_ew%THeatLossRoot2Soil_vr(L)

    DO  K=1,micpar%NumOfPlantLitrCmplxs
      DO  M=1,jsken
        DO NE=1,NumPlantChemElms
          LitrfalStrutElms_vr(NE,M,K,L,NY,NX)=plt_bgcr%LitrfalStrutElms_vr(NE,M,K,L)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO L=1,NK_col(NY,NX)
    TWaterPlantRoot2Soil_vr(L,NY,NX)              = plt_ew%TWaterPlantRoot2Soil_vr(L)
    totRootLenDens_vr(L,NY,NX)                    = plt_morph%totRootLenDens_vr(L)
    trcg_root_vr(idg_beg:idg_NH3,L,NY,NX)         = plt_rbgc%trcg_root_vr(idg_beg:idg_NH3,L)
    trcg_air2root_flx_vr(idg_beg:idg_NH3,L,NY,NX) = plt_rbgc%trcg_air2root_flx_vr(idg_beg:idg_NH3,L)
    RootCO2Emis2Root_vr(L,NY,NX)                  = plt_bgcr%RootCO2Emis2Root_vr(L)
    RUptkRootO2_vr(L,NY,NX)                       = plt_bgcr%RUptkRootO2_vr(L)
    RootO2_TotSink_vr(L,NY,NX)                    = plt_bgcr%RootO2_TotSink_vr(L)
    trcs_Soil2plant_uptake_vr(ids_beg:ids_end,L,NY,NX) =plt_rbgc%trcs_Soil2plant_uptake_vr(ids_beg:ids_end,L)

    associate(RootSoilCO2Debt_vr => plt_rbgc%RootSoilCO2Debt_vr)
      ! Materialize the predictor loan with equal aqueous mass and debt.
      ! Adding debt alone would charge fixed transport/root uptake twice.
      ! The pair leaves soil CO2 mass minus debt unchanged; normal transport
      ! stages this obligation once and SubstrateDribbling repays/carries it.
      ! The persistent dribble pool already participates in balances/restarts.
      trcs_solml_vr(idg_CO2,L,NY,NX)=trcs_solml_vr(idg_CO2,L,NY,NX)+RootSoilCO2Debt_vr(L)
      trcs_solml_drib_vr(idg_CO2,L,NY,NX)=trcs_solml_drib_vr(idg_CO2,L,NY,NX)+RootSoilCO2Debt_vr(L)
      RootSoilCO2Debt_vr(L)=0._r8
    end associate

    DO  K=1,jcplx
      tRootMycoExud2Soil_vr(1:NumPlantChemElms,K,L,NY,NX)=plt_bgcr%tRootMycoExud2Soil_vr(1:NumPlantChemElms,K,L)
    ENDDO
    RootMycoMassElm_vr(:,:,L,NY,NX)=0._r8
  ENDDO

  end subroutine ReceivePlantColumns

  subroutine SendPlantColumnInputs(I,NY,NX,I1)
  ! Description:
  ! Load current column forcing, surface exposure and geometry for the full
  ! plant model before its energy/water and biogeochemical processes.
  use EcoSIMConfig, only : jsken=>jskenc,jcplx=>jcplxc
  integer, intent(in) :: I,NY,NX
  integer, intent(out) :: I1
  integer :: NZ,K,L,N,ids

  IF((ALAT_col(NY,NX).GE.0.0_r8.AND.I.EQ.1) .OR. (ALAT_col(NY,NX).LT.0.0_r8.AND.I.EQ.1))THEN
    DO NZ=1,NP0_col(NY,NX)
      plt_morph%lreset_laimax_pft(NZ)=.true.
    ENDDO
  ELSE
    DO NZ=1,NP0_col(NY,NX)
      plt_morph%lreset_laimax_pft(NZ)=.false.
    ENDDO
  endif
  plt_site%NY=NY;plt_site%NX=NX
  plt_site%DazCurrYear=DazCurrYear
  I1=I+1;if(I1>DazCurrYear)I1=1
  plt_site%SoilSurfDepZ_col           = SoilSurfDepZ_col(NY,NX)
  plt_site%ZERO                       = ZERO
  plt_site%ZERO2                      = ZERO2
  plt_site%ALAT                       = ALAT_col(NY,NX)
  plt_site%ATCA                       = ATCA_col(NY,NX)
  plt_ew%BulkFactor4Snow_col          =BulkFactor4Snow_col(NY,NX)
  plt_morph%LeafStalkAreaAll_col         = LeafStalkAreaAll_col(NY,NX)
  plt_morph%CanopyLeafArea_col        = CanopyLeafArea_col(NY,NX)
  plt_site%ALT                        = SurfElevation_col(NY,NX)
  plt_site%CCO2EI_gperm3                     = CCO2EI_gperm3_col(NY,NX)
  plt_site%CO2EI                      = CO2EI_col(NY,NX)
  plt_bgcr%NetCO2Flx2Canopy_col       = NetCO2Flx2Canopy_col(NY,NX)
  plt_site%CO2E                       = CO2E_col(NY,NX)
  plt_site%AtmGasc(idg_beg:idg_NH3) = AtmGasCgperm3_col(idg_beg:idg_NH3,NY,NX)
  plt_site%DayLenthPrev               = DayLenthPrev_col(NY,NX)
  plt_site%DayLenthCurrent            = DayLensCurr_col(NY,NX)
  plt_ew%SnowDepth                    = SnowDepth_col(NY,NX)
  plt_ew%canopy_exposed_col           = canopy_exposed_col(NY,NX)
  plt_site%DayLenthMax_col                = DayLenthMax_col(NY,NX)
  plt_site%KoppenClimZone             = KoppenClimZone_col(NY,NX)
  plt_site%iYearCurrent               = iYearCurrent
  plt_site%NL                         = NL_col(NY,NX)
  plt_site%NP0                        = NP0_col(NY,NX)
  plt_site%MaxNumRootLays             = MaxNumRootLays_col(NY,NX)
  plt_site%NP                         = NP_col(NY,NX)
  plt_site%NU                         = NUM_col(NY,NX)
  plt_site%NK                         = NK_col(NY,NX)
  plt_site%OXYE                       = OXYE_col(NY,NX)
  plt_ew%RawCanopyH2SinkZ_col         = RawCanopyH2SinkZ_col(NY,NX)
  plt_ew%RawIsoTAtm2CanopySinkZ_col   = RawIsoTAtm2CanopySinkZ_col(NY,NX)
  plt_ew%RIB                          = RIB_col(NY,NX)
  plt_rad%SineSunInclAnglNxtHour_col  = SineSunInclAnglNxtHour_col(NY,NX)
  plt_rad%SineSunInclinationAngle_col        = SineSunInclinationAngle_col(NY,NX)
  plt_ew%TKSnow                       = TKSnow_snvr(1,NY,NX)  !surface layer snow temperature
  plt_ew%TairK                        = TairK_col(NY,NX)
  plt_rad%LWRadGrnd_col               = LWRadGrnd_col(NY,NX)
  plt_rad%LWRadSky_col                = LWRadSky_col(NY,NX)
  plt_ew%VPA                          = VPA_col(NY,NX)
  plt_ew%EMS_Modify_Scalar_col        = EMS_Modify_Scalar_col(NY,NX)
  plt_distb%XCORP                     = XTillCorp_col(NY,NX)
  plt_site%SolarNoonHour_col          = SolarNoonHour_col(NY,NX)
  plt_site%ZEROS2                     = ZEROS2(NY,NX)
  plt_site%ZEROS                      = ZEROS(NY,NX)
  plt_ew%RoughnessLength                  = RoughnessLength_col(NY,NX)
  plt_morph%CanopyHeight_col          = CanopyHeight_col(NY,NX)
  plt_ew%ZeroPlaneDisplacem_col       = ZeroPlaneDisplacem_col(NY,NX)
  plt_distb%DCORP                     = DepzCorp_col(I,NY,NX)
  plt_distb%iSoilDisturbType_col      = iSoilDisturbType_col(I,NY,NX)
  plt_morph%CanopyHeightZ_col(0)      = CanopyHeightZ_col(0,NY,NX)
  DO  L=1,NumCanopyLayers
    plt_morph%CanopyHeightZ_col(L) = CanopyHeightZ_col(L,NY,NX)
  ENDDO
  DO  L=1,NumCanopyLayers+1
    plt_rad%TAU_DirectSunLit(L) = TAU_DirectSunLit(L,NY,NX)
    plt_rad%TAU_DirectSunSha(L)         = TAU_DirectSunSha(L,NY,NX)
  ENDDO

  DO N=1,NumLeafInclinationClasses
    plt_rad%SineLeafAngle(N)=SineLeafAngle(N)
  ENDDO

  DO L=1,NL_col(NY,NX)
    plt_soilchem%HYCDMicP4RootUptake_vr(L) = HYCDMicP4RootUptake_vr(L,NY,NX)
    plt_soilchem%GasDifcT_vr(idg_beg:idg_end,L)  = GasDifcT_vr(idg_beg:idg_end,L,NY,NX)
    plt_soilchem%SoilBulkModulus4RootPent_vr(L)   = SoilBulkModulus4RootPent_vr(L,NY,NX)
    plt_soilchem%SoilModulus4RootRadialexp_vr(L) = SoilModulus4RootRadialexp_vr(L,NY,NX)
    plt_site%CumSoilThickMidL_vr(L)             = CumSoilThickMidL_vr(L,NY,NX)
  ENDDO

  DO L=1,NK_col(NY,NX)
    plt_soilchem%trcg_gasml_vr(idg_beg:idg_NH3,L) = trcg_gasml_vr(idg_beg:idg_NH3,L,NY,NX)
  ENDDO

  plt_site%CumSoilThickness_vr(0)                       = CumSoilThickness_vr(0,NY,NX)
  DO L=1,NK_col(NY,NX)
    plt_site%CumSoilThickness_vr(L)                       = CumSoilThickness_vr(L,NY,NX)
    plt_site%AREA3(L)                                     = AREA_3D(3,L,NY,NX)
    plt_soilchem%SoilBulkDensity_vr(L)                    = SoilBulkDensity_vr(L,NY,NX)
    plt_soilchem%trc_solcl_vr(ids_beg:ids_end,L)          = trc_solcl_vr(ids_beg:ids_end,L,NY,NX)
    plt_soilchem%SoluteDifusvtyT_vr(ids_beg:ids_end,L)     = SoluteDifusvtyT_vr(ids_beg:ids_end,L,NY,NX)
    plt_soilchem%trcg_gascl_vr(idg_beg:idg_NH3,L)         = trcg_gascl_vr(idg_beg:idg_NH3,L,NY,NX)
    plt_soilchem%CSoilOrgM_vr(ielmc,L)                    = CSoilOrgM_vr(ielmc,L,NY,NX)
    plt_site%FracSoiAsMicP_vr(L)                          = FracSoiAsMicP_vr(L,NY,NX)
    DO ids=ids_beg,ids_end
      plt_soilchem%trcs_solml_vr(ids,L)         =AZMAX1(trcs_solml_vr(ids,L,NY,NX)-trcs_solml_drib_vr(ids,L,NY,NX))
    ENDDO
    plt_soilchem%GasSolbility_vr(idg_beg:idg_NH3,L)       = GasSolbility_vr(idg_beg:idg_NH3,L,NY,NX)
!    plt_soilchem%trcs_RMicbUptake_vr(idg_beg:idg_NH3-1,L) = trcs_RMicbUptake_vr(idg_beg:idg_NH3-1,L,NY,NX)
    plt_ew%ElvAdjstedSoilH2OPSIMPa_vr(L)                  = ElvAdjstedSoilH2OPSIMPa_vr(L,NY,NX)
    plt_bgcr%RH2PO4EcoDmndSoilPrev_vr(L)                  = RH2PO4EcoDmndSoilPrev_vr(L,NY,NX)
    plt_bgcr%RH2PO4EcoDmndBandPrev_vr(L)                  = RH2PO4EcoDmndBandPrev_vr(L,NY,NX)
    plt_bgcr%RH1PO4EcoDmndSoilPrev_vr(L)                  = RH1PO4EcoDmndSoilPrev_vr(L,NY,NX)
    plt_bgcr%RH1PO4EcoDmndBandPrev_vr(L)                  = RH1PO4EcoDmndBandPrev_vr(L,NY,NX)
    plt_bgcr%RNO3EcoDmndSoilPrev_vr(L)                    = RNO3EcoDmndSoilPrev_vr(L,NY,NX)
    plt_bgcr%RNH4EcoDmndSoilPrev_vr(L)                    = RNH4EcoDmndSoilPrev_vr(L,NY,NX)
    plt_bgcr%RNH4EcoDmndBandPrev_vr(L)                    = RNH4EcoDmndBandPrev_vr(L,NY,NX)
    plt_bgcr%RNO3EcoDmndBandPrev_vr(L)                    = RNO3EcoDmndBandPrev_vr(L,NY,NX)
    plt_bgcr%RGasTranspFlxPrev_vr(idg_beg:idg_NH3,L)      = RGasTranspFlxPrev_vr(idg_beg:idg_NH3,L,NY,NX)
    plt_bgcr%RO2AquaSourcePrev_vr(L)                      = RO2AquaSourcePrev_vr(L,NY,NX)
    plt_bgcr%RO2EcoDmndPrev_vr(L)                         = RO2EcoDmndPrev_vr(L,NY,NX)
    plt_ew%TKS_vr(L)                                      = TKS_vr(L,NY,NX)
    plt_soilchem%THETW_vr(L)               = THETW_vr(L,NY,NX)
    plt_soilchem%SoilWatAirDry_vr(L)       = SoilWatAirDry_vr(L,NY,NX)
    plt_soilchem%TScal4Difsvity_vr(L)      = TScal4Difsvity_vr(L,NY,NX)
    plt_soilchem%VLSoilPoreMicP_vr(L)      = VLSoilPoreMicP_vr(L,NY,NX)
    plt_soilchem%trcs_VLN_vr(ids_H1PO4B,L) = trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    plt_soilchem%trcs_VLN_vr(ids_NO3,L)    = trcs_VLN_vr(ids_NO3,L,NY,NX)
    plt_soilchem%trcs_VLN_vr(ids_H1PO4,L)  = trcs_VLN_vr(ids_H1PO4,L,NY,NX)
    plt_soilchem%VLSoilMicP_vr(L)          = VLSoilMicP_vr(L,NY,NX)
    plt_soilchem%VLiceMicP_vr(L)           = VLiceMicP_vr(L,NY,NX)
    plt_soilchem%VLWatMicP_vr(L)           = VLWatMicP_vr(L,NY,NX)
    plt_soilchem%VLMicP_vr(L)              = VLMicP_vr(L,NY,NX)
    plt_soilchem%trcs_VLN_vr(ids_NO3B,L)   = trcs_VLN_vr(ids_NO3B,L,NY,NX)
    plt_soilchem%trcs_VLN_vr(ids_NH4,L)    = trcs_VLN_vr(ids_NH4,L,NY,NX)
    plt_soilchem%trcs_VLN_vr(ids_NH4B,L)   = trcs_VLN_vr(ids_NH4B,L,NY,NX)
    plt_site%DLYR3(L)                      = DLYR_3D(3,L,NY,NX)
    DO K=1,jcplx
      plt_soilchem%FracBulkSOMC_vr(K,L)          = FracBulkSOMC_vr(K,L,NY,NX)
      plt_soilchem%DOM_MicP_vr(idom_doc:idom_dop,K,L) = DOM_MicP_vr(idom_doc:idom_dop,K,L,NY,NX)
      plt_soilchem%DOM_MicP_drib_vr(idom_doc:idom_dop,K,L)=DOM_MicP_drib_vr(idom_doc:idom_dop,K,L,NY,NX)
    ENDDO
  ENDDO

  end subroutine SendPlantColumnInputs

  subroutine SendPlantColumnState(I1,NY,NX)
  use EcoSIMConfig, only : jsken=>jskenc,jcplx=>jcplxc
  integer, intent(in) :: I1,NY,NX
  integer :: K,L,M,NE

  DO L=1,NK_col(NY,NX)
    plt_site%SoilWeightStress_vr(L) = SoilWeightStress_vr(L,NY,NX)
    plt_site%SoilSuctStress_vr(L) = PSISoilMatricP_vr(L,NY,NX)+PSISoilOsmotic_vr(L,NY,NX)
    plt_site%rSat_vr(L)           = safe_adb(VLWatMicP_vr(L,NY,NX),VLSoilMicP_vr(L,NY,NX))
    DO M=1,NPH
      plt_site%VLWatMicPM_vr(M,L)               = VLWatMicPM_vr(M,L,NY,NX)
      plt_site%VLsoiAirPM_vr(M,L)               = VLsoiAirPM_vr(M,L,NY,NX)
      plt_site%TortMicPM_vr(M,L)                = TortMicPM_vr(M,L,NY,NX)
      plt_site%FILMM_vr(M,L)                    = FILMM_vr(M,L,NY,NX)
      plt_soilchem%DiffusivitySolutEffM_vr(M,L) = DiffusivitySolutEffM_vr(M,L,NY,NX)
    ENDDO
  ENDDO

  ! sent variables also modified
  plt_site%NumActivePlants                               = NumActivePlants_col(NY,NX)
  plt_site%QH2OLoss_lnds                                 = QH2OLoss_lnds
  plt_site%PlantPopu_col                                 = PlantPopu_col(NY,NX)
  plt_bgcr%ECO_ER_col                                    = ECO_ER_col(NY,NX)
  plt_biom%StandingDeadStrutElms_col(1:NumPlantChemElms) = StandingDeadStrutElms_col(1:NumPlantChemElms,NY,NX)
  plt_bgcr%LitrFallStrutElms_col(1:NumPlantChemElms)     = LitrFallStrutElms_col(1:NumPlantChemElms,NY,NX)
  plt_morph%StemArea_col                                 = StemArea_col(NY,NX)
  plt_ew%Eco_Heat_Sens_col                               = Eco_Heat_Sens_col(NY,NX)
  plt_ew%WatHeldOnCanopy_col                             = WatHeldOnCanopy_col(NY,NX)
  plt_ew%SnowOnCanopy_col                                = SnowOnCanopy_col(NY,NX)
  plt_bgcr%Eco_NBP_CumYr_col                             = Eco_NBP_CumYr_col(NY,NX)
  plt_ew%Air_Heat_Latent_store_col                       = Air_Heat_Latent_store_col(NY,NX)
  plt_ew%Air_Heat_Sens_store_col                         = Air_Heat_Sens_store_col(NY,NX)
  plt_bgcr%Eco_AutoR_CumYr_col                           = Eco_AutoR_CumYr_col(NY,NX)
  plt_ew%H2OLoss_CumYr_col                               = H2OLoss_CumYr_col(NY,NX)
  plt_distb%EcoHavstElmnt_CumYr_col(1:NumPlantChemElms)  = EcoHavstElmnt_CumYr_col(1:NumPlantChemElms,NY,NX)
  plt_rad%Eco_NetRad_col                                 = Eco_NetRad_col(NY,NX)
  plt_ew%VapXAir2Canopy_col                              = VapXAir2Canopy_col(NY,NX)
  plt_ew%Eco_Heat_Latent_col                             = Eco_Heat_Latent_col(NY,NX)
  plt_rbgc%TRootGasLossDisturb_col(idg_beg:idg_NH3)    = TRootGasLossDisturb_col(idg_beg:idg_NH3,NY,NX)
  plt_ew%Eco_Heat_GrndSurf_col                           = Eco_Heat_GrndSurf_col(NY,NX)
  plt_ew%QVegET_col                                        = QVegET_col(NY,NX)
  plt_ew%HeatSensVapAir2Canopy_col                              = HeatSensVapAir2Canopy_col(NY,NX)
  plt_ew%HeatPrecAir2Canopy_col                              = HeatPrecAir2Canopy_col(NY,NX)
  plt_ew%HeatFlx2Canopy_col                              = HeatFlx2Canopy_col(NY,NX)
  plt_ew%LWRadCanG                                       = LWRadCanG_col(NY,NX)
  plt_ew%CanopyBiomWater_col                                   = CanopyBiomWater_col(NY,NX)
  plt_ew%CanopyHeatStor_col                              = CanopyHeatStor_col(NY,NX)
  plt_bgcr%Canopy_NEE_col                                = Canopy_NEE_col(NY,NX)
  plt_distb%FERT(1:20)                                   = FERT(1:20,I1,NY,NX)
  plt_ew%CanopyHeatLoss2Dist_col                             = CanopyHeatLoss2Dist_col(NY,NX)
  plt_ew%QCanopyWatLoss2Dist_col                             = QCanopyWatLoss2Dist_col(NY,NX)
  DO  L=1,NumCanopyLayers
    plt_morph%CanopyStemAareZ_col(L) = CanopyStemAareZ_col(L,NY,NX)
    plt_biom%tCanLeafC_clyr(L)         = tCanLeafC_clyr(L,NY,NX)
    plt_morph%CanopyLeafAareZ_col(L) = CanopyLeafAareZ_col(L,NY,NX)
  ENDDO

  DO L=0,NL_col(NY,NX)
    DO K=1,jcplx
      DO NE=1,NumPlantChemElms
        plt_bgcr%REcoDOMProd_vr(NE,K,L)=REcoDOMProd_vr(NE,K,L,NY,NX)
      ENDDO
    ENDDO
  ENDDO


  DO L=0,NL_col(NY,NX)
    plt_bgcr%REcoH2PO4DmndBand_vr(L) = REcoH2PO4DmndBand_vr(L,NY,NX)
    plt_bgcr%REcoH1PO4DmndBand_vr(L) = REcoH1PO4DmndBand_vr(L,NY,NX)
    plt_bgcr%REcoNO3DmndBand_vr(L)   = REcoNO3DmndBand_vr(L,NY,NX)
    plt_bgcr%REcoNH4DmndBand_vr(L)   = REcoNH4DmndBand_vr(L,NY,NX)
    plt_bgcr%REcoH1PO4DmndSoil_vr(L) = REcoH1PO4DmndSoil_vr(L,NY,NX)
    plt_bgcr%REcoH2PO4DmndSoil_vr(L) = REcoH2PO4DmndSoil_vr(L,NY,NX)
    plt_bgcr%REcoNO3DmndSoil_vr(L)   = REcoNO3DmndSoil_vr(L,NY,NX)
    plt_bgcr%REcoNH4DmndSoil_vr(L)   = REcoNH4DmndSoil_vr(L,NY,NX)
    plt_bgcr%REcoO2DmndResp_vr(L)    = REcoO2DmndResp_vr(L,NY,NX)
    plt_ew%THeatLossRoot2Soil_vr(L)     = THeatLossRoot2Soil_vr(L,NY,NX)

    DO  K=1,micpar%NumOfPlantLitrCmplxs
      DO  M=1,jsken
        DO NE=1,NumPlantChemElms
          plt_bgcr%LitrfalStrutElms_vr(NE,M,K,L)=LitrfalStrutElms_vr(NE,M,K,L,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO L=1,NK_col(NY,NX)
    plt_ew%TWaterPlantRoot2Soil_vr(L)                = TWaterPlantRoot2Soil_vr(L,NY,NX)
    plt_morph%totRootLenDens_vr(L)                   = totRootLenDens_vr(L,NY,NX)
    plt_rbgc%trcg_root_vr(idg_beg:idg_NH3,L)         = trcg_root_vr(idg_beg:idg_NH3,L,NY,NX)
    plt_rbgc%trcg_air2root_flx_vr(idg_beg:idg_NH3,L) = trcg_air2root_flx_vr(idg_beg:idg_NH3,L,NY,NX)
    plt_bgcr%RootCO2Emis2Root_vr(L)                  = RootCO2Emis2Root_vr(L,NY,NX)
    plt_bgcr%RUptkRootO2_vr(L)                       = RUptkRootO2_vr(L,NY,NX)
    plt_bgcr%RootO2_TotSink_vr(L)                    = RootO2_TotSink_vr(L,NY,NX)
    DO  K=1,jcplx
      plt_bgcr%tRootMycoExud2Soil_vr(1:NumPlantChemElms,K,L)=tRootMycoExud2Soil_vr(1:NumPlantChemElms,K,L,NY,NX)
    ENDDO
  ENDDO

  end subroutine SendPlantColumnState
end module PlantColumnTransferMod
