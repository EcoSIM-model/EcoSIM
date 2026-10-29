submodule (HistDataType) HistUpdateColumn
  use GridConsts, only: JZ
  use ElmIDMod, only: ielmc, ielmn, ielmp, NumPlantChemElms
  use EcosimConst, only: natomw, patomw
  use data_const_mod, only: spval  => DAT_CONST_SPVAL
  use MicrobialDiagMod, only: SumMicbGroup, sumDOML, sumMicBiomLayL
  use MiniMathMod, only: safe_adb, AZMAX1, VapMass2KPa,AZERO
  use EcoSiMParDataMod, only: micpar
  use MLDataDiagType, only: AcetConc30cm_col, AcetConc60cm_col, AcetMGC30cm_col, AcetMGC60cm_col, &
    AeroHRBactC30cm_col, AeroHRBactC60cm_col, AeroHRFungC30cm_col, AeroHRFungC60cm_col, AeroMOC30cm_col, &
    AeroMOC60cm_col, CH4wConc30cm_col, CH4wConc60cm_col, DOCConc30cm_col, DOCConc60cm_col, FermC30cm_col, &
    FermC60cm_col, H2MGC30cm_col, H2MGC60cm_col, H2wConc30cm_col, H2wConc60cm_col, O2wConc30cm_col, &
    O2wConc60cm_col, RAeroCH4Oxi30cm_col, RAeroCH4Oxi60cm_col, RCH4ProdAcet30cm_col, RCH4ProdAcet60cm_col, &
    RCH4ProdHG30cm_col, RCH4ProdHG60cm_col, RCO2Ht30cm_col, RCO2Ht60cm_col, RFerment30cm_col, RFerment60cm_col, &
    TEMP30cm_col, TEMP60cm_col, THETW30cm_col, THETW60cm_col
  use LandSurfDataType, only: RoughnessLength_col, TKQ_col, VPQ_col, ZeroPlaneDisplacem_col
  use EcoSIMCtrlMod, only: lmicrobeMLdiag
  use MicrobialDataType, only: tRespGrossHeter_vr, tRespGrossHeterUlm_vr
  use CanopyRadDataType, only: RadSW_Canopy_col
  use GridDataType, only: AREA_3D, CumDepz2LayBottom_vr, DLYR_3D, NU_col
  use EcoSIMCtrlDataType, only: ZEROS
  use PlantTraitDataType, only: CanopyLeafArea_col, StemArea_col
  use PlantDataRateType, only: Eco_NEE_col, idg_Ar, idg_CH4, idg_CO2, idg_H2, idg_N2, idg_N2O, idg_NH3, idg_O2, &
    idom_acetate, idom_beg, idom_doc, idom_don, idom_dop, idom_end, ids_H2PO4, ids_NH4, ids_NO2, ids_NO3, &
    idx_H2PO4, idx_HPO4, idx_NH4, RootN2Fix_col, RUptkRootO2_col
  use ClimForcDataType, only: CH4E_col, CO2E_col, Eco_RadSW_col, LWRadSky_col, NWetDep_col, PBOT_col, &
    RadPARSolarBeam_col, RadSWSolarBeam_col, RainFalPrec_col, SnoFalPrec_col, TairK_col, TRAD_col, VPK_col, &
    WindSpeedAtm_col
  use CanopyDataType, only: QVegET_col, RadPAR2LitR_col, RadPAR2Soil_col, RadPARGrnd_col, RadSWGrnd_col, &
    SnowOnCanopy_col, VapXAir2Canopy_col
  use SOMDataType, only: DIC_mass_col, FracLitrMix_vr, SoilOrgM_vr, tHumOM_col, tHxPO4_col, tLitrOM_col, &
    tMicBiome_col, tNH4_col, tNO3_col, tOMActC_vr, tSoilOrgM_col, TSolidOMActC_vr, TSolidOMC_vr
  use SoilPhysDataType, only: ActiveLayDepZ_col, CondGasXSurf_col
  use SoilHeatDatatype, only: HeatDischar_col, HeatFlx2Grnd_col, HeatStore_col, TCS_vr, VHeatCapacity_vr
  use SoilWaterDataType, only: WaterTableDepZi_col, NaturalWaterTableDepzTopoC_col, PSISoilMatricP_vr, QDischarg2WTBL_col, QDrain_col, &
    QEvap_CumYr_col, Qinflx2Soil_col, QRain_CumYr_col, QRunSurf_col, ThetaH2OZ_vr, ThetaICEZ_vr, VLiceMicP_vr, &
    VLWatMicP_vr, WatMass_col
  use SnowDataType, only: SnowDepth_col, TCSnow_snvr, VcumSnowWE_col
  use PlantMgmtDataType, only: CH4byFire_CumYr_col, CO2byFire_CumYr_col, PO4byFire_CumYr_col
  use EcosimBGCFluxType, only: Canopy_NEE_col, CumDryDepoOM_col, Eco_AutoR_CumYr_col, Eco_GPP_CumYr_col, &
    Eco_Heat_GrndSurf_col, Eco_Heat_Latent_col, Eco_Heat_Sens_col, ECO_HR_CO2_col, ECO_HR_CO2_vr, &
    Eco_HR_CumYr_col, Eco_NBP_CumYr_col, Eco_NetRad_col, Eco_NPP_CumYr_col, EcoHavstElmnt_CumYr_col, &
    NetNH4Mineralize_CumYr_col, NetPO4Mineralize_CumYr_col
  use SoilPropertyDataType, only: VLSoilMicPMass_vr
  use SurfLitterDataType, only: FracSurfByLitR_col, VLitR_col
  use SoilBGCDataType, only: AmendC_CumYr_flx_col, FerP_Flx_CumYr_col, FertN_Flx_CumYr_col, &
    Gas_Prod_TP_cumRes_col, Gas_WetDeposit_flx_col, GasDiff2Surf_flx_col, GasHydroLoss_cumflx_col, &
    HydroIonFlx_CumYr_col, Hydroloss_NH4_cumflx_col, Hydroloss_NO3_cumflx_col, HydroSubsDICFlx_col, &
    HydroSubsDINFlx_col, HydroSubsDIPFlx_col, HydroSubsDOCFlx_col, HydroSubsDONFlx_col, HydroSubsDOPFlx_col, &
    HydroSufDICFlx_col, HydroSufDINFlx_CumYr_col, HydroSufDIPFlx_CumYr_col, HydroSufDOCFlx_col, &
    HydroSufDONFlx_col, HydroSufDOPFlx_col, LiterfalOrgM_col, Micb_N2Fixation_vr, MoistSensDecomp_vr, &
    OxyDecompLimiter_vr, RCH4Oxi_aero_vr, RCH4Oxi_anmo_vr, RCH4ProdAcetcl_vr, RCH4ProdHydrog_vr, &
    RDen_N2OtoN2_vr, RDen_NO2toN2O_vr, RDen_NO3toNO2_vr, RFerment_vr, RN2OChemoProd_vr, RN2ONitProd_vr, &
    RNit_NH3toNO2_vr, RNit_NO2toNO3_vr, RO2DecompUptk_vr, SedmErossLoss_CumYr_col, StandingDeadStrutElms_col, &
    SurfGasEmiss_all_flx_col, TempSensDecomp_vr, TMicHeterActivity_vr, trc_solcl_vr, trcg_air2root_flx_col, &
    trcg_ebu_flx_col, trcg_soilml_vr, trcg_soilMass_col, trcg_TotalMass_col, trcs_drainage_flx_col, trcs_solml_vr, &
    tRHydlyBioReSOM_vr, tRHydlySOM_vr, tRHydlySoprtOM_vr, tXPO4_col
  use AqueChemDatatype, only: TProd_CO2_geochem_soil_vr, trcg_mass_cumerr_col, trcx_solml_vr
  use SurfSoilDataType, only: FracSurfSnowCovered_col, HeatByRad2Surf_col, HeatEvapAir2Surf_col, HeatNet2Surf_col, &
    HeatSensAir2Surf_col, HeatSensVapAir2Surf_col, VapXAir2GSurf_col
  use TracerPropMod, only: GramPerHr2umolPerSec
  implicit none
contains

  module procedure update_hist_columns
    real(r8) :: micBE(1:NumPlantChemElms)
    real(r8) :: DOM(idom_beg:idom_end)
    real(r8), parameter :: secs1hour=3600._r8
    real(r8), parameter :: MJ2W=1.e6_r8/secs1hour
    real(r8), parameter :: m2mm=1000._r8
    real(r8), parameter :: million=1.e6_r8

      this%h1D_cumFIRE_CO2_col(ncol)        =  CO2byFire_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_cumFIRE_CH4_col(ncol)        =  CH4byFire_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_cNH4_LITR_col(ncol)        =  safe_adb(trcs_solml_vr(ids_NH4,0,NY,NX)+&
        natomw*trcx_solml_vr(idx_NH4,0,NY,NX),VLSoilMicPMass_vr(0,NY,NX)*million)
      this%h1D_cNO3_LITR_col(ncol)        =  safe_adb(trcs_solml_vr(ids_NO3,0,NY,NX)+&
        trcs_solml_vr(ids_NO2,0,NY,NX),VLSoilMicPMass_vr(0,NY,NX)*million)
      
      this%h1D_ECO_HVST_N_col(ncol)   = EcoHavstElmnt_CumYr_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_NET_N_MIN_col(ncol)    = -NetNH4Mineralize_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tLITR_P_col(ncol)      = tLitrOM_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_HUMUS_C_col(ncol)      = tHumOM_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_HUMUS_N_col(ncol)      = tHumOM_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_HUMUS_P_col(ncol)      = tHumOM_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_AMENDED_P_col(ncol)    = FerP_Flx_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tLITRf_C_FLX_col(ncol) = LiterfalOrgM_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tLITRf_N_FLX_col(ncol) = LiterfalOrgM_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tLITRf_P_FLX_col(ncol) = LiterfalOrgM_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tEXCH_PO4_col(ncol)        = tHxPO4_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUR_DOP_FLX_col(ncol)      = HydroSufDOPFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUB_DOP_FLX_col(ncol)      = HydroSubsDOPFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUR_DIP_FLX_col(ncol)      = HydroSufDIPFlx_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUB_DIP_FLX_col(ncol)      = HydroSubsDIPFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_HeatFlx2Grnd_col(ncol)     = HeatFlx2Grnd_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CumDryDepoOM_col(ncol) = CumDryDepoOM_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CanSWRad_col(ncol)         = MJ2W*RadSW_Canopy_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RadSW_Grnd_col(ncol)       = MJ2W*RadSWGrnd_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RadPAR_Grnd_col(ncol)      = RadPARGrnd_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RadPAR2Soil_col(ncol)      = RadPAR2Soil_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RadPAR2LitR_col(ncol)      = RadPAR2LitR_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Qinfl2soi_col(ncol)        = m2mm*Qinflx2Soil_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Qdrain_col(ncol)           = m2mm*QDrain_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

      this%h1D_SUR_DON_FLX_col(ncol)      = HydroSufDONFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUB_DON_FLX_col(ncol)      = HydroSubsDONFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSALT_DISCHG_FLX_col(ncol) = HydroIonFlx_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUR_DIN_FLX_col(ncol)      = HydroSufDINFlx_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUB_DIN_FLX_col(ncol)      = HydroSubsDINFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUR_DOC_FLX_col(ncol)      = HydroSufDOCFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUB_DOC_FLX_col(ncol)      = HydroSubsDOCFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUR_DIC_FLX_col(ncol)      = HydroSufDICFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUB_DIC_FLX_col(ncol)      = HydroSubsDICFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUR_DIP_FLX_col(ncol)      = HydroSufDIPFlx_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tPREC_P_col(ncol)          = tXPO4_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tMICRO_P_col(ncol)         = tMicBiome_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SnowCanopy_col(ncol)        =SnowOnCanopy_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_PO4_FIRE_col(ncol)         = PO4byFire_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_cPO4_LITR_col(ncol)        = safe_adb(trcs_solml_vr(ids_H2PO4,0,NY,NX),VLSoilMicPMass_vr(0,NY,NX)*million)
      this%h1D_cEXCH_P_LITR_col(ncol)     = patomw*safe_adb(trcx_solml_vr(idx_HPO4,0,NY,NX)+&
        trcx_solml_vr(idx_H2PO4,0,NY,NX),VLSoilMicPMass_vr(0,NY,NX)*million)
      this%h1D_ECO_HVST_P_col(ncol) = EcoHavstElmnt_CumYr_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_NET_P_MIN_col(ncol)  = -NetPO4Mineralize_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_PSI_SURF_col(ncol)   = PSISoilMatricP_vr(0,NY,NX)
      this%h1D_SURF_ELEV_col(ncol)  = -CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX)+DLYR_3D(3,0,NY,NX)
      this%h1D_tLITR_N_col(ncol)    = tLitrOM_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_AMENDED_N_col(ncol)  = FertN_Flx_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tNH4X_col(ncol)      = tNH4_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tNO3_col(ncol)       = tNO3_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tRAD_col(ncol)       = TRAD_col(NY,NX)
      if(this%h1D_tNH4X_col(ncol)<0._r8)then
        write(*,*)'negative tNH4X',this%h1D_tNH4X_col(ncol),this%h1D_tNO3_col(ncol)
        stop
      endif
      this%h1D_tMICRO_N_col(ncol)         = tMicBiome_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_TEMP_LITR_col(ncol)        = TCS_vr(0,NY,NX)
      this%h1D_TEMP_surf_col(ncol)        = FracSurfByLitR_col(NY,NX)*TCS_vr(0,NY,NX)+(1._r8-FracSurfByLitR_col(NY,NX))*TCS_vr(NU_col(NY,NX),NY,NX)
      if(VcumSnowWE_col(NY,NX)<=ZEROS(NY,NX))then
        this%h1D_TEMP_SNOW_col(ncol)   = spval
      else
        this%h1D_TEMP_SNOW_col(ncol)   = TCSnow_snvr(1,NY,NX)
      endif
      this%h1D_FracBySnow_col(ncol) = FracSurfSnowCovered_col(NY,NX)
      this%h1D_FracByLitr_col(ncol) = FracSurfByLitR_col(NY,NX)
      this%h1D_tLITR_C_col(ncol)    = tLitrOM_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

      if(lmicrobeMLdiag)then
        this%h1d_TEMP60cm_col(ncol)         = TEMP60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_THETW60cm_col(ncol)        = THETW60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_O2wConc60cm_col(ncol)      = O2wConc60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_AcetConc60cm_col(ncol)     = AcetConc60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_H2wConc60cm_col(ncol)      = H2wConc60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_CH4wConc60cm_col(ncol)     = CH4wConc60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_DOCConc60cm_col(ncol)      = DOCConc60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_AcetMGC60cm_col(ncol)      = AcetMGC60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_H2MGC60cm_col(ncol)        = H2MGC60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_FermC60cm_col(ncol)        = FermC60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_AeroMOC60cm_col(ncol)      = AeroMOC60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_AeroHRFungC60cm_col(ncol)  = AeroHRFungC60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_AeroHRBactC60cm_col(ncol)  = AeroHRBactC60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_RAeroCH4Oxi60cm_col(ncol)  = RAeroCH4Oxi60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_RCH4ProdAcet60cm_col(ncol) = RCH4ProdAcet60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_RCH4ProdHG60cm_col(ncol)   = RCH4ProdHG60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_RFerment60cm_col(ncol)     = RFerment60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_RCO2Ht60cm_col(ncol)       = RCO2Ht60cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_TEMP30cm_col(ncol)         = TEMP30cm_col(NY,NX) /AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_THETW30cm_col(ncol)        = THETW30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_O2wConc30cm_col(ncol)      = O2wConc30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_AcetConc30cm_col(ncol)     = AcetConc30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_H2wConc30cm_col(ncol)      = H2wConc30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_CH4wConc30cm_col(ncol)     = CH4wConc30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_DOCConc30cm_col(ncol)      = DOCConc30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_AcetMGC30cm_col(ncol)      = AcetMGC30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_H2MGC30cm_col(ncol)        = H2MGC30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_FermC30cm_col(ncol)        = FermC30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_AeroMOC30cm_col(ncol)      = AeroMOC30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_AeroHRFungC30cm_col(ncol)  = AeroHRFungC30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_AeroHRBactC30cm_col(ncol)  = AeroHRBactC30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_RAeroCH4Oxi30cm_col(ncol)  = RAeroCH4Oxi30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_RCH4ProdAcet30cm_col(ncol) = RCH4ProdAcet30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_RCH4ProdHG30cm_col(ncol)   = RCH4ProdHG30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_RFerment30cm_col(ncol)     = RFerment30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1d_RCO2Ht30cm_col(ncol)       = RCO2Ht30cm_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      endif


      this%h1D_AMENDED_C_col(ncol)        = AmendC_CumYr_flx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tMICRO_C_col(ncol)         = tMicBiome_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSoilOrgC_col(ncol)        = tSoilOrgM_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSoilOrgN_col(ncol)        = tSoilOrgM_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSoilOrgP_col(ncol)        = tSoilOrgM_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_OMC_LITR_col(ncol)         = SoilOrgM_vr(ielmc,0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_OMN_LITR_col(ncol)         = SoilOrgM_vr(ielmn,0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_OMP_LITR_col(ncol)         = SoilOrgM_vr(ielmp,0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ATM_CO2_col(ncol)          = CO2E_col(NY,NX)
      this%h1D_ATM_CH4_col(ncol)          = CH4E_col(NY,NX)
      this%h1D_NBP_col(ncol)              = Eco_NBP_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      
      this%h1D_ECO_HVST_C_col(ncol)       = EcoHavstElmnt_CumYr_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_LAI_col(ncol)          = CanopyLeafArea_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_SAI_col(ncol)          = StemArea_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Eco_GPP_CumYr_col(ncol)    = Eco_GPP_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_RA_col(ncol)           = Eco_AutoR_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Eco_NPP_CumYr_col(ncol)    = Eco_NPP_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Eco_HR_CumYr_col(ncol)     = Eco_HR_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Eco_HR_CO2_col(ncol)       = ECO_HR_CO2_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
!      this%h1D_Eco_HR_CH4_col(ncol)       = ECO_HR_CH4_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tDIC_col(ncol)             = DIC_mass_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSTANDING_DEAD_C_col(ncol) = StandingDeadStrutElms_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSTANDING_DEAD_N_col(ncol) = StandingDeadStrutElms_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSTANDING_DEAD_P_col(ncol) = StandingDeadStrutElms_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tPRECIP_col(ncol)           = m2mm*QRain_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_ET_col(ncol)           = m2mm*QEvap_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_Ar_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_CO2_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_CH4_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_CH4,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_O2_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_O2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_N2_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_N2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_NH3_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_NH3,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_H2_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_H2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

      this%h1D_ECO_RADSW_col(ncol)        = MJ2W*Eco_RadSW_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_N2O_LITR_col(ncol)         = trc_solcl_vr(idg_N2O,0,NY,NX)
      this%h1D_NH3_LITR_col(ncol)         = trc_solcl_vr(idg_NH3,0,NY,NX)
      this%h1D_SOL_RADN_col(ncol)         = RadSWSolarBeam_col(NY,NX)*MJ2W
      this%h1D_AIR_TEMP_col(ncol)         = TairK_col(NY,NX)-273.15_r8
      this%h1D_HUM_col(ncol)              = VPK_col(NY,NX)
      this%h1D_PATM_col(ncol)             = PBOT_col(NY,NX)
      this%h1D_WIND_col(ncol)             = WindSpeedAtm_col(NY,NX)/secs1hour
      this%h1D_PREC_col(ncol)             = (RainFalPrec_col(NY,NX)+SnoFalPrec_col(NY,NX))*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Snofall_col(ncol)          = SnoFalPrec_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SOIL_RN_col(ncol)          = HeatByRad2Surf_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_LWSky_col(ncol)            = LWRadSky_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SOIL_LE_col(ncol)          = HeatEvapAir2Surf_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SOIL_H_col(ncol)           = HeatSensAir2Surf_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SOIL_G_col(ncol)           = -(HeatNet2Surf_col(NY,NX)-HeatSensVapAir2Surf_col(NY,NX))*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_RN_col(ncol)           = Eco_NetRad_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_LE_col(ncol)           = Eco_Heat_Latent_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Eco_HeatSen_col(ncol)      = Eco_Heat_Sens_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_Heat2G_col(ncol)       = Eco_Heat_GrndSurf_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_O2_LITR_col(ncol)          = trc_solcl_vr(idg_O2,0,NY,NX)
      this%h1D_CO2_SEMIS_FLX_col(ncol)    = SurfGasEmiss_all_flx_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
      this%h1D_AR_SEMIS_FLX_col(ncol)     = SurfGasEmiss_all_flx_col(idg_AR,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_AR)
      this%h1D_ECO_CO2_FLX_col(ncol)      = Eco_NEE_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
      this%h1d_CAN_NEE_col(ncol)          = Canopy_NEE_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
      this%h1D_CH4_SEMIS_FLX_col(ncol)    = SurfGasEmiss_all_flx_col(idg_CH4,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
      this%h1D_O2_SEMIS_FLX_col(ncol)     = SurfGasEmiss_all_flx_col(idg_O2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_O2)
      this%h1D_CH4_EBU_flx_col(ncol)      = trcg_ebu_flx_col(idg_CH4,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CH4)
      this%h1D_Ar_EBU_flx_col(ncol)       = trcg_ebu_flx_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_Ar)      
      this%h1D_AR_PLTROOT_flx_col(ncol)   = trcg_air2root_flx_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_Ar)
      
      this%h1D_CH4_PLTROOT_flx_col(ncol)  = trcg_air2root_flx_col(idg_CH4,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CH4)
      this%h1D_CO2_PLTROOT_flx_col(ncol)  = trcg_air2root_flx_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
      this%h1D_O2_PLTROOT_flx_col(ncol)   = trcg_air2root_flx_col(idg_O2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_O2)
      this%h1D_VPQ_col(ncol)                = VapMass2KPa(VPQ_col(NY,NX),TKQ_col(NY,NX),fwd=.FALSE.)
      this%h1D_TKQ_col(ncol)                = TKQ_col(NY,NX)
      this%h1D_RoughnessLength_col(ncol)    = RoughnessLength_col(NY,NX)
      this%h1D_ZeroPlaneDisplacem_col(ncol) = ZeroPlaneDisplacem_col(NY,NX)
      this%h1D_CO2_DIF_flx_col(ncol)        = GasDiff2Surf_flx_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
      this%h1D_CH4_DIF_flx_col(ncol)      = GasDiff2Surf_flx_col(idg_CH4,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CH4)
      this%h1D_Ar_DIF_flx_col(ncol)       = GasDiff2Surf_flx_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_Ar)
      this%h1D_NH3_DIF_flx_col(ncol)      = GasDiff2Surf_flx_col(idg_NH3,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_NH3)
      this%h1D_O2_DIF_flx_col(ncol)        =GasDiff2Surf_flx_col(idg_O2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_O2)
      this%h1D_CO2_TPR_err_col(ncol)      = Gas_Prod_TP_cumRes_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Ar_TPR_err_col(ncol)       = Gas_Prod_TP_cumRes_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CO2_Drain_flx_col(ncol)    = trcs_drainage_flx_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CO2_hydloss_flx_col(ncol)  = GasHydroLoss_cumflx_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_NH3_hydloss_flx_col(ncol)  = (GasHydroLoss_cumflx_col(idg_NH3,NY,NX)+Hydroloss_NH4_cumflx_col(NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_NO3_hydloss_flx_col(ncol)  = (Hydroloss_NO3_cumflx_col(NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CO2_LITR_col(ncol)         = trc_solcl_vr(idg_CO2,0,NY,NX)
      this%h1D_CH4_LITR_col(ncol)         = trc_solcl_vr(idg_CH4,0,NY,NX)
      this%h1D_H2_LITR_col(ncol)          = trc_solcl_vr(idg_H2,0,NY,NX)
      this%h1D_EVAPG_col(ncol)            = VapXAir2GSurf_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CondGasXSurf_col(ncol)     = CondGasXSurf_col(NY,NX)
      this%h1D_CanopyEvap_col(ncol)       = VapXAir2Canopy_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CANET_col(ncol)            = QVegET_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RUNOFF_FLX_col(ncol)       = -QRunSurf_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SEDIMENT_FLX_col(ncol)     = SedmErossLoss_CumYr_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSWC_col(ncol)             = WatMass_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tHeat_col(ncol)            = HeatStore_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_QDISCHG_FLX_col(ncol)      = QDischarg2WTBL_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_HeatDISCHG_FLX_col(ncol)   = HeatDischar_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SNOWPACK_col(ncol)         = AZMAX1((VcumSnowWE_col(NY,NX))*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX))
      if(this%h1D_SNOWPACK_col(ncol)>0._r8)then
        this%h1D_SNOWDENS_col(ncol)         = this%h1D_SNOWPACK_col(ncol)/SnowDepth_col(NY,NX)
      else
        this%h1D_SNOWDENS_col(ncol)         =0._r8
      endif
      this%h1D_SURF_WTR_col(ncol)         = ThetaH2OZ_vr(0,NY,NX)
      this%h1D_SURF_ICE_col(ncol)         = ThetaICEZ_vr(0,NY,NX)
      this%h1D_ThetaW_litr_col(ncol)      = safe_adb(VLWatMicP_vr(0,NY,NX),VLitR_col(NY,NX))
      this%h1D_ThetaI_litr_col(ncol)      = safe_adb(VLiceMicP_vr(0,NY,NX),VLitR_col(NY,NX))
      this%h1D_ACTV_LYR_col(ncol)         = -(ActiveLayDepZ_col(NY,NX)-CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX))
      this%h1D_WTR_TBL_col(ncol)          = -(WaterTableDepZi_col(NY,NX)-CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX))
      this%h1D_WTBL_TopoC_col(ncol)       = -NaturalWaterTableDepzTopoC_col(NY,NX)
      this%h1D_N2O_SEMIS_FLX_col(ncol)         = SurfGasEmiss_all_flx_col(idg_N2O,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_N2_SEMIS_FLX_col(ncol)         = SurfGasEmiss_all_flx_col(idg_N2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_NH3_SEMIS_FLX_col(ncol)         = SurfGasEmiss_all_flx_col(idg_NH3,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_H2_SEMIS_FLX_col(ncol)          = SurfGasEmiss_all_flx_col(idg_H2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_PAR_col(ncol)            = RadPARSolarBeam_col(NY,NX)
      this%h1D_NWetDep_flx_col(ncol)  = NWetDep_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_VHeatCap_litr_col(ncol)  = VHeatCapacity_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_AR_WetDep_FLX_col(ncol)  = Gas_WetDeposit_flx_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CO2_WetDep_FLX_col(ncol) = Gas_WetDeposit_flx_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RootXO2_flx_col(ncol)    = RUptkRootO2_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RootN_Fix_col(ncol)      = RootN2Fix_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      call sumMicBiomLayL(0,NY,NX,micBE)
      this%h2D_MicroBiomeE_litr_col(ncol,1:NumPlantChemElms) =micBE/AREA_3D(3,NU_col(NY,NX),NY,NX)

      call SumMicbGroup(0,NY,NX,micpar%mid_HeterAerobBacter,MicbE)
      this%h2D_AeroHrBactE_litr_col(ncol,1:NumPlantChemElms) =micBE/AREA_3D(3,NU_col(NY,NX),NY,NX)     !aerobic heterotropic bacteria

      call SumMicbGroup(0,NY,NX,micpar%mid_Aerob_Fungi,MicbE)
      this%h2D_AeroHrFungE_litr_col(ncol,1:NumPlantChemElms) = micBE/AREA_3D(3,NU_col(NY,NX),NY,NX)   !aerobic heterotropic fungi

      call SumMicbGroup(0,NY,NX,micpar%mid_Facult_DenitBacter,MicbE)
      this%h2D_faculDenitE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)  !facultative denitrifier

      call SumMicbGroup(0,NY,NX,micpar%mid_fermentor,MicbE)
      this%h2D_fermentorE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)  !fermentor

      call SumMicbGroup(0,NY,NX,micpar%mid_HeterMixtCynoBacter,MicbE)
      this%h2D_cyanoBactC_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)  !cyanobacteria
      this%h1D_cyanoBactC_col(ncol) = MicbE(ielmc)

      call SumMicbGroup(0,NY,NX,micpar%mid_HeterAcetoCH4GenArchea,MicbE)
      this%h2D_acetometgE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)  !acetogenic methanogen

      call SumMicbGroup(0,NY,NX,micpar%mid_HeterAerobN2Fixer,MicbE)
      this%h2D_aeroN2fixE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)  !aerobic N2 fixer

      call SumMicbGroup(0,NY,NX,micpar%mid_HeterAnaerobN2Fixer,MicbE)
      this%h2D_anaeN2FixE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)  !anaerobic N2 fixer

      call SumMicbGroup(0,NY,NX,micpar%mid_AutoAmmoniaOxidBacter,MicbE,isauto=.true.)
      this%h2D_NH3OxiBactE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)

      call SumMicbGroup(0,NY,NX,micpar%mid_AutoNitriteOxidBacter,MicbE,isauto=.true.)
      this%h2D_NO2OxiBactE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)

      call SumMicbGroup(0,NY,NX,micpar%mid_AutoAeroCH4OxiBacter,MicbE,isauto=.true.)
      this%h2D_CH4AeroOxiE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)

      associate(mid_ANME2d => micpar%mid_AutoAMOANME2D, mid_NC10 => micpar%mid_AutoAMONC10)
        call SumMicbGroup(0,NY,NX,mid_ANME2d,MicbE,isauto=.true.)
        this%h2D_ANME2dE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)
        call SumMicbGroup(0,NY,NX,mid_NC10,MicbE,isauto=.true.)
        this%h2D_NC10E_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)
      end associate

      call SumMicbGroup(0,NY,NX,micpar%mid_AutoH2GenoCH4GenArchea,MicbE,isauto=.true.)
      this%h2D_H2MethogenE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)

      call sumDOML(0,NY,NX,DOM)

      this%h1D_DOC_LITR_col(ncol)     = safe_adb(DOM(idom_doc),VLitR_col(NY,NX))
      this%h1D_DON_LITR_col(ncol)     = safe_adb(DOM(idom_don),VLitR_col(NY,NX))
      this%h1D_DOP_LITR_col(ncol)     = safe_adb(DOM(idom_dop),VLitR_col(NY,NX))
      this%h1D_acetate_LITR_col(ncol) = safe_adb(DOM(idom_acetate),VLitR_col(NY,NX))

      this%h1D_RCH4ProdHydrog_litr_col(ncol) = RCH4ProdHydrog_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RCH4ProdAcetcl_litr_col(ncol) = RCH4ProdAcetcl_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RCH4Oxi_aero_litr_col(ncol)   = RCH4Oxi_aero_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RCH4Oxi_aero_col(ncol) =    RCH4Oxi_aero_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RCH4Oxi_ANMO_litr_col(ncol)   = RCH4Oxi_anmo_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RCH4Oxi_ANMO_col(ncol) =    RCH4Oxi_anmo_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RDen_NO3toNO2_col(ncol) = RDen_NO3toNO2_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      
      this%h1D_RFermen_litr_col(ncol)        = AZERO(RFerment_vr(0,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_NH3oxi_litr_col(ncol)         = RNit_NH3toNO2_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_NO2Oxi_litr_col(ncol)         = RNit_NO2toNO3_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_N2Oprod_litr_col(ncol)        = (RDen_NO2toN2O_vr(0,NY,NX)+RN2ONitProd_vr(0,NY,NX) &
                               +RN2OChemoProd_vr(0,NY,NX)-RDen_N2OtoN2_vr(0,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)

      this%h1D_decomp_OStress_litr_col(ncol)   = OxyDecompLimiter_vr(0,NY,NX)
      this%h1D_MicrobAct_litr_col(ncol)        = TMicHeterActivity_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RO2Decomp_litr_col(ncol)        = RO2DecompUptk_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RDECOMPC_SOM_litr_col(ncol)     = tRHydlySOM_vr(ielmc,0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RDECOMPC_BReSOM_litr_col(ncol)  = tRHydlyBioReSOM_vr(ielmc,0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RDECOMPC_SorpSOM_litr_col(ncol) = tRHydlySoprtOM_vr(ielmc,0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tRespGrossHete_litr_col(ncol)   = tRespGrossHeter_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tRespGrossHeteUlm_litr_col(ncol) = tRespGrossHeterUlm_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

      this%h1D_Decomp_temp_FN_litr_col(ncol)   = TempSensDecomp_vr(0,NY,NX)
      this%h1D_Decomp_moist_FN_litr_col(ncol)  = MoistSensDecomp_vr(0,NY,NX)
      this%h1D_FracLitMix_litr_col(ncol)       = FracLitrMix_vr(0,NY,NX)
      
      this%h1D_Eco_HR_CO2_litr_col(ncol)       = ECO_HR_CO2_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_TSolidOMActC_litr_col(ncol)     = TSolidOMActC_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_TSolidOMActCDens_litr_col(ncol) = safe_adb(TSolidOMActC_vr(0,NY,NX),TSolidOMC_vr(0,NY,NX))
      this%h1D_tOMActCDens_litr_col(ncol)      = safe_adb(tOMActC_vr(0,NY,NX),(TSolidOMC_vr(0,NY,NX)+tOMActC_vr(0,NY,NX)))
      this%h1D_Ar_mass_col(ncol)               = trcg_TotalMass_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CO2_MassLit_col(ncol) = trcg_soilml_vr(idg_CO2,0,NY,NX) &
        /AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_O2_MassLit_col(ncol) = trcg_soilml_vr(idg_O2,0,NY,NX) &
        /AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CH4_MassLit_col(ncol) = trcg_soilml_vr(idg_CH4,0,NY,NX) &
        /AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_H2_MassLit_col(ncol) = trcg_soilml_vr(idg_H2,0,NY,NX) &
        /AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Ar_soilMass_col(ncol)           = trcg_soilMass_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CO2_mass_col(ncol)               = trcg_TotalMass_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Gchem_CO2_prod_col(ncol)         = sum(TProd_CO2_geochem_soil_vr(1:JZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)

      this%h1D_tDOC_soil_col(ncol)=0._r8
      this%h1D_tDON_soil_col(ncol)=0._r8
      this%h1D_tDOP_soil_col(ncol)=0._r8
      this%h1D_tAcetate_soil_col(ncol)=0._r8
      
      this%h1D_FreeNFix_col(ncol)=Micb_N2Fixation_vr(0,NY,NX)

  end procedure update_hist_columns

end submodule HistUpdateColumn
