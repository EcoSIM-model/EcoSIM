submodule (HistDataType) HistUpdateSoil
  use GridConsts, only: JZ
  use ElmIDMod, only: ielmc, ielmn, ielmp, ipltroot, NumPlantChemElms
  use EcoSIMConfig, only: jcplx => jcplxc
  use EcosimConst, only: natomw, patomw
  use MicrobialDiagMod, only: SumMicbGroup, sumDOML, sumMicBiomLayL, SumSolidOML
  use MiniMathMod, only: safe_adb, AZMAX1
  use EcoSiMParDataMod, only: micpar
  use EcoSIMCtrlMod, only: plant_model
  use MicrobialDataType, only: FermOXYI_vr, tRespGrossHeter_vr, tRespGrossHeterUlm_vr
  use GridDataType, only: AREA_3D, CumDepz2LayBottom_vr, DLYR_3D, NU_col, NL_col
  use PlantDataRateType, only: idg_Ar, idg_CH4, idg_CO2, idg_H2, idg_N2, idg_N2O, idg_NH3, idg_O2, &
    idom_acetate, idom_beg, idom_doc, idom_don, idom_dop, idom_end, ids_H1PO4, ids_H1PO4B, ids_H2PO4, &
    ids_H2PO4B, ids_NH4, ids_NH4B, ids_NO2, ids_NO2B, ids_NO3, ids_NO3B, idx_H2PO4, idx_H2PO4B, idx_HPO4, &
    idx_HPO4B, idx_NH4, idx_NH4B, RootCO2Ar2Root_vr, RootCO2Ar2Soil_vr, RootCO2Autor_vr, THeatLossRoot2Soil_vr, &
    trcg_root_vr, TWaterPlantRoot2SoilPrev_vr
  use RootDataType, only: RootMycoMassElm_vr
  use SOMDataType, only: FracLitrMix_vr, litrOM_vr, RHydlySOCK_vr, SoilOrgM_vr, SolidOM_vr, SorbedOM_vr, &
    tOMActC_vr, TSolidOMActC_vr, TSolidOMC_vr
  use SoilPhysDataType, only: PAR_RAD_vr
  use SoilHeatDatatype, only: TCS_vr, THeatFlowCellSoil_vr, VHeatCapacity_vr
  use SoilWaterDataType, only: HydCondSoil_3D, PSISoilMatricP_vr, PSISoilOsmotic_vr, QDrainloss_vr, &
    SoilBulkModulus4RootPent_vr, ThetaH2OZ_vr, ThetaICEZ_vr, TWatFlowCellMacP_vr, TWatFlowCellMicP_vr, &
    VLWatMicP_vr
  use EcosimBGCFluxType, only: ECO_HR_CO2_vr
  use SoilPropertyDataType, only: POROS_vr, VLSoilMicPMass_vr
  use SoilBGCDataType, only: AeroBact_PrimeS_lim_vr, AeroFung_PrimeS_lim_vr, Micb_N2Fixation_vr, &
    MoistSensDecomp_vr, OxyDecompLimiter_vr, RCH4Oxi_aero_vr, RCH4Oxi_anmo_vr, RCH4ProdAcetcl_vr, &
    RCH4ProdHydrog_vr, RDen_N2OtoN2_vr, RDen_NO2toN2O_vr, RDen_NO3toNO2_vr, RFerment_vr, &
    RHydrolysisScalCmpK_vr, RN2OChemoProd_vr, RN2ONitProd_vr, RNit_NH3toNO2_vr, RNit_NO2toNO3_vr, &
    RO2DecompUptk_vr, ROQC4HeterMicActCmpK_vr, Soil_Gas_Frac_vr, Soil_Gas_pressure_vr, TempSensDecomp_vr, &
    TMicHeterActivity_vr, trc_solcl_vr, trcg_soilml_vr, trcs_solml_vr, tRHydlyBioReSOM_vr, tRHydlySOM_vr, tRHydlySoprtOM_vr
  use AqueChemDatatype, only: ElectricConductivity_vr, TProd_CO2_geochem_soil_vr, trcx_solml_vr
  implicit none
contains

  module procedure update_hist_soil
    integer :: L
    real(r8) :: micBE(1:NumPlantChemElms)
    real(r8) :: DOM(idom_beg:idom_end)
    real(r8), parameter :: m2mm=1000._r8
    real(r8) :: DVOLL
    real(r8) :: SOMC(jcplx)
    integer :: jj


    ! Layer 0 is litter. Sum only active soil layers, excluding unused slots,
    ! and convert cell-total C rates to rates per ground area exactly once.
    this%h1D_RCH4ProdAcet_col(ncol)=(RCH4ProdAcetcl_vr(0,NY,NX) &
      +SUM(RCH4ProdAcetcl_vr(NU_col(NY,NX):NL_col(NY,NX),NY,NX)))/AREA_3D(3,NU_col(NY,NX),NY,NX)
    this%h1D_RCH4ProdHG_col(ncol)=(RCH4ProdHydrog_vr(0,NY,NX) &
      +SUM(RCH4ProdHydrog_vr(NU_col(NY,NX):NL_col(NY,NX),NY,NX)))/AREA_3D(3,NU_col(NY,NX),NY,NX)
    this%h1D_RCH4Oxi_col(ncol)=(RCH4Oxi_aero_vr(0,NY,NX)+RCH4Oxi_anmo_vr(0,NY,NX) &
      +SUM(RCH4Oxi_aero_vr(NU_col(NY,NX):NL_col(NY,NX),NY,NX)) &
      +SUM(RCH4Oxi_anmo_vr(NU_col(NY,NX):NL_col(NY,NX),NY,NX)))/AREA_3D(3,NU_col(NY,NX),NY,NX)

      DO L=1,JZ
        call SumSolidOML(ielmc,L,NY,NX,SOMC)
        DO jj=1,jcplx
          this%h3D_SOC_Cps_vr(ncol,L,jj) = SOMC(jj)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        ENDDO
        this%h2D_CO2_MassL_vr(ncol,L) = trcg_soilml_vr(idg_CO2,L,NY,NX) &
          /AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_O2_MassL_vr(ncol,L) = trcg_soilml_vr(idg_O2,L,NY,NX) &
          /AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_CH4_MassL_vr(ncol,L) = trcg_soilml_vr(idg_CH4,L,NY,NX) &
          /AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_H2_MassL_vr(ncol,L) = trcg_soilml_vr(idg_H2,L,NY,NX) &
          /AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_Gas_Pressure_vr(ncol,L) = Soil_Gas_pressure_vr(L,NY,NX)
        this%h2D_CO2_Gas_ppmv_vr(ncol,L) = Soil_Gas_Frac_vr(idg_CO2,L,NY,NX)
        this%h2D_CH4_Gas_ppmv_vr(ncol,L) = Soil_Gas_Frac_vr(idg_CH4,L,NY,NX)
        this%h2D_Ar_Gas_ppmv_vr(ncol,L)  = Soil_Gas_Frac_vr(idg_Ar,L,NY,NX)
        this%h2D_O2_Gas_ppmv_vr(ncol,L)  = Soil_Gas_Frac_vr(idg_O2,L,NY,NX)
        this%h2D_H2_Gas_ppmv_vr(ncol,L)  = Soil_Gas_Frac_vr(idg_H2,L,NY,NX)
        this%h2D_N2O_Gas_ppmv_vr(ncol,L) = Soil_Gas_Frac_vr(idg_N2O,L,NY,NX)
        this%h2D_N2_Gas_ppmv_vr(ncol,L)  = Soil_Gas_Frac_vr(idg_N2,L,NY,NX)
        this%h2D_NH3_Gas_ppmv_vr(ncol,L) = Soil_Gas_Frac_vr(idg_NH3,L,NY,NX)

        DVOLL=DLYR_3D(3,L,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)

        this%h2D_Eco_HR_CO2_vr(ncol,L)    = ECO_HR_CO2_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_Gchem_CO2_prod_vr(ncol,L)= TProd_CO2_geochem_soil_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_SoilBulkStress_vr(ncol,L) = SoilBulkModulus4RootPent_vr(L,NY,NX)
        if(DVOLL<=1.e-8_r8)cycle

        call sumDOML(L,NY,NX,DOM)
        call sumMicBiomLayL(L,NY,NX,micBE,I,J)
        this%h1D_tDOC_soil_col(ncol)=this%h1D_tDOC_soil_col(ncol)+DOM(idom_doc)
        this%h1D_tDON_soil_col(ncol)=this%h1D_tDON_soil_col(ncol)+DOM(idom_don)
        this%h1D_tDOP_soil_col(ncol)=this%h1D_tDOP_soil_col(ncol)+DOM(idom_dop)
        this%h1D_tAcetate_soil_col(ncol)=this%h1D_tAcetate_soil_col(ncol)+DOM(idom_acetate)
        this%h2D_QDrainloss_vr(ncol,L)=QDrainloss_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_FermOXYI_vr(ncol,L)=FermOXYI_vr(L,NY,NX)
        this%h2D_DOC_vr(ncol,L)             = DOM(idom_doc)/DVOLL
        this%h2D_DON_vr(ncol,L)             = DOM(idom_don)/DVOLL
        this%h2D_DOP_vr(ncol,L)             = DOM(idom_dop)/DVOLL
        this%h2D_BotDEPZ_vr(ncol,L)         = CumDepz2LayBottom_vr(L,NY,NX)
        this%h2D_acetate_vr(ncol,L)         = DOM(idom_acetate)/DVOLL
        this%h2D_litrC_vr(ncol,L)           = litrOM_vr(ielmc,L,NY,NX)/DVOLL
        this%h2D_litrN_vr(ncol,L)           = litrOM_vr(ielmn,L,NY,NX)/DVOLL
        this%h2D_litrP_vr(ncol,L)           = litrOM_vr(ielmp,L,NY,NX)/DVOLL
        this%h2D_tSOC_vr(ncol,L)            = SoilOrgM_vr(ielmc,L,NY,NX)/DVOLL

        this%h2D_POM_C_vr(ncol,L)           = sum(SolidOM_vr(ielmc,:,micpar%k_POM,L,NY,NX))*safe_adb(1.e-3_r8,VLSoilMicPMass_vr(L,NY,NX))
        this%h2D_MAOM_C_vr(ncol,L)          = (sum(SorbedOM_vr(idom_doc,:,L,NY,NX))+sum(SorbedOM_vr(idom_acetate,:,L,NY,NX)))*safe_adb(1.e-3_r8,VLSoilMicPMass_vr(L,NY,NX))
        this%h2D_microbC_vr(ncol,L)         =  micBE(ielmc)/DVOLL
        this%h2D_microbN_vr(ncol,L)         =  micBE(ielmn)/DVOLL
        this%h2D_microbP_vr(ncol,L)         =  micBE(ielmp)/DVOLL
        this%h2D_AeroBact_PrimS_lim_vr(ncol,L)=AeroBact_PrimeS_lim_vr(L,NY,NX)
        this%h2D_AeroFung_PrimS_lim_vr(ncol,L)=AeroFung_PrimeS_lim_vr(L,NY,NX)
        this%h2D_tSOCL_vr(ncol,L)           = SoilOrgM_vr(ielmc,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_tSON_vr(ncol,L)            = SoilOrgM_vr(ielmn,L,NY,NX)/DVOLL
        this%h2D_tSOP_vr(ncol,L)            = SoilOrgM_vr(ielmp,L,NY,NX)/DVOLL
        this%h2D_VHeatCap_vr(ncol,L)        = VHeatCapacity_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_Root_CO2_vr(ncol,L)        = trcg_root_vr(idg_CO2,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_Root_O2_vr(ncol,L) = trcg_root_vr(idg_O2,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_Root_CH4_vr(ncol,L) = trcg_root_vr(idg_CH4,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_Root_H2_vr(ncol,L) = trcg_root_vr(idg_H2,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_NO3_vr(ncol,L)             = trc_solcl_vr(ids_NO3,L,NY,NX)
        this%h2D_NH4_vr(ncol,L)             = trc_solcl_vr(ids_NH4,L,NY,NX)
        this%h2D_Aqua_CO2_vr(ncol,L)        = trc_solcl_vr(idg_CO2,L,NY,NX)
        this%h2D_Aqua_CH4_vr(ncol,L)        = trc_solcl_vr(idg_CH4,L,NY,NX)
        this%h2D_Aqua_O2_vr(ncol,L)         = trc_solcl_vr(idg_O2,L,NY,NX)
        this%h2D_Aqua_N2O_vr(ncol,L)        = trc_solcl_vr(idg_N2O,L,NY,NX)
        this%h2D_Aqua_NH3_vr(ncol,L)        = trc_solcl_vr(idg_NH3,L,NY,NX)
        this%h2D_Aqua_N2_vr(ncol,L)         = trc_solcl_vr(idg_N2,L,NY,NX)
        this%h2D_Aqua_Ar_vr(ncol,L)         = trc_solcl_vr(idg_Ar,L,NY,NX)
        this%h2D_Aqua_H2_vr(ncol,L)         = trc_solcl_vr(idg_H2,L,NY,NX)

        this%h2D_TEMP_vr(ncol,L)            = TCS_vr(L,NY,NX)
        this%h2D_decomp_OStress_vr(ncol,L)  = OxyDecompLimiter_vr(L,NY,NX)
        this%h2D_RO2Decomp_vr(ncol,L)       = RO2DecompUptk_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_Decomp_temp_FN_vr(ncol,L)  = TempSensDecomp_vr(L,NY,NX)
        this%h2D_FracLitMix_vr(ncol,L)      = FracLitrMix_vr(L,NY,NX)
        this%h2D_Decomp_Moist_FN_vr(ncol,L) = MoistSensDecomp_vr(L,NY,NX)
        this%h2D_HeatFlow_vr(ncol,L)        = THeatFlowCellSoil_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_HeatUptk_vr(ncol,L)        = THeatLossRoot2Soil_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_VSPore_vr(ncol,L)          = POROS_vr(L,NY,NX)
        this%h2D_FLO_MICP_vr(ncol,L)        = m2mm*TWatFlowCellMicP_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_FLO_MACP_vr(ncol,L)        = m2mm*TWatFlowCellMacP_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_rVSM_vr    (ncol,L)        = ThetaH2OZ_vr(L,NY,NX)
        this%h2D_rVSICE_vr  (ncol,L)        = ThetaICEZ_vr(L,NY,NX)
        this%h2D_PSI_vr(ncol,L)             = PSISoilMatricP_vr(L,NY,NX)+PSISoilOsmotic_vr(L,NY,NX)
        this%h2D_PsiO_vr(ncol,L)            = PSISoilOsmotic_vr(L,NY,NX)
        this%h2D_RootH2OUP_vr(ncol,L)       = TWaterPlantRoot2SoilPrev_vr(L,NY,NX)
        this%h2D_cNH4t_vr(ncol,L)           = safe_adb(trcs_solml_vr(ids_NH4,L,NY,NX)+trcs_solml_vr(ids_NH4B,L,NY,NX) &
                                               +natomw*(trcx_solml_vr(idx_NH4,L,NY,NX)+trcx_solml_vr(idx_NH4B,L,NY,NX)),&
                                               VLSoilMicPMass_vr(L,NY,NX))

        this%h2D_cNO3t_vr(ncol,L)= safe_adb(trcs_solml_vr(ids_NO3,L,NY,NX)+trcs_solml_vr(ids_NO3B,L,NY,NX) &
                                               +trcs_solml_vr(ids_NO2,L,NY,NX)+trcs_solml_vr(ids_NO2B,L,NY,NX),&
                                               VLSoilMicPMass_vr(L,NY,NX))

        this%h2D_cPO4_vr(ncol,L) = safe_adb(trcs_solml_vr(ids_H1PO4,L,NY,NX)+trcs_solml_vr(ids_H1PO4B,L,NY,NX) &
                                               +trcs_solml_vr(ids_H2PO4,L,NY,NX)+trcs_solml_vr(ids_H2PO4B,L,NY,NX),&
                                               VLWatMicP_vr(L,NY,NX))
        this%h2D_cEXCH_P_vr(ncol,L)= patomw*safe_adb(trcx_solml_vr(idx_HPO4,L,NY,NX)+trcx_solml_vr(idx_H2PO4,L,NY,NX) &
                                               +trcx_solml_vr(idx_HPO4B,L,NY,NX)+trcx_solml_vr(idx_H2PO4B,L,NY,NX),&
                                               VLSoilMicPMass_vr(L,NY,NX))
        this%h2D_microb_N2fix_vr(ncol,L)         = 1.e6_r8*Micb_N2Fixation_vr(L,NY,NX)/DVOLL
        this%h1D_FreeNFix_col(ncol)              = this%h1D_FreeNFix_col(ncol)+ Micb_N2Fixation_vr(L,NY,NX)
        this%h2D_ElectricConductivity_vr(ncol,L) = ElectricConductivity_vr(L,NY,NX)

        this%h2D_HydCondSoil_vr(ncol,L) = HydCondSoil_3D(3,L,NY,NX)

        call SumMicbGroup(L,NY,NX,micpar%mid_HeterMixtCynoBacter,MicbE)
        this%h2D_cyanoBactC_vr(ncol,L) = MicbE(ielmc)/DVOLL
        this%h1D_cyanoBactC_col(ncol) = this%h1D_cyanoBactC_col(ncol)+MicbE(ielmc)
        !aerobic heterotropic bacteria
        call SumMicbGroup(L,NY,NX,micpar%mid_HeterAerobBacter,MicbE)
        this%h2D_AeroHrBactC_vr(ncol,L) = MicbE(ielmc)/DVOLL
        this%h2D_AeroHrBactN_vr(ncol,L) = MicbE(ielmn)/DVOLL
        this%h2D_AeroHrBactP_vr(ncol,L) = MicbE(ielmp)/DVOLL

        !facultative denitrifier
        call SumMicbGroup(L,NY,NX,micpar%mid_Facult_DenitBacter,MicbE)
        this%h2D_faculDenitC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_faculDenitN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_faculDenitP_vr(ncol,L) = micBE(ielmp)/DVOLL

        !aerobic heterotropic fungi
        call SumMicbGroup(L,NY,NX,micpar%mid_Aerob_Fungi,MicbE)
        this%h2D_AeroHrFungC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_AeroHrFungN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_AeroHrFungP_vr(ncol,L) = micBE(ielmp)/DVOLL

        !fermentor
        call SumMicbGroup(L,NY,NX,micpar%mid_fermentor,MicbE)
        this%h2D_fermentorC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_fermentorN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_fermentorP_vr(ncol,L) = micBE(ielmp)/DVOLL

          this%h2D_fermentor_frac_vr(ncol,L)=safe_adb(this%h2D_fermentorC_vr(ncol,L),this%h2D_microbC_vr(ncol,L))

        !acetogenic methanogen
        call SumMicbGroup(L,NY,NX,micpar%mid_HeterAcetoCH4GenArchea,MicbE)
        this%h2D_acetometgC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_acetometgN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_acetometgP_vr(ncol,L) = micBE(ielmp)/DVOLL

        this%h2D_acetometh_frac_vr(ncol,L)=safe_adb(this%h2D_acetometgC_vr(ncol,L),this%h2D_microbC_vr(ncol,L))

        !aerobic N2 fixer
        call SumMicbGroup(L,NY,NX,micpar%mid_HeterAerobN2Fixer,MicbE)
        this%h2D_aeroN2fixC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_aeroN2fixN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_aeroN2fixP_vr(ncol,L) = micBE(ielmp)/DVOLL

        this%h2D_RDECOMPC_SOM_vr(ncol,L)     = tRHydlySOM_vr(ielmc,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RDECOMPC_BReSOM_vr(ncol,L)  = tRHydlyBioReSOM_vr(ielmc,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RDECOMPC_SorpSOM_vr(ncol,L) = tRHydlySoprtOM_vr(ielmc,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_MicrobAct_vr(ncol,L)        = TMicHeterActivity_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

        DO jj=1,jcplx
          this%h3D_MicrobActCps_vr(ncol,L,jj) = ROQC4HeterMicActCmpK_vr(jj,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

          this%h3D_SOMHydrylScalCps_vr(ncol,L,jj) = AZMAX1(RHydrolysisScalCmpK_vr(jj,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX))

          this%h3D_HydrolCSOMCps_vr(ncol,L,jj) = RHydlySOCK_vr(jj,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        ENDDO
        !anaerobic N2 fixer
        call SumMicbGroup(L,NY,NX,micpar%mid_HeterAnaerobN2Fixer,MicbE)
        this%h2D_anaeN2FixC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_anaeN2FixN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_anaeN2FixP_vr(ncol,L) = micBE(ielmp)/DVOLL

        call SumMicbGroup(L,NY,NX,micpar%mid_AutoAmmoniaOxidBacter,MicbE,isauto=.true.)
        this%h2D_NH3OxiBactC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_NH3OxiBactN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_NH3OxiBactP_vr(ncol,L) = micBE(ielmp)/DVOLL

        call SumMicbGroup(L,NY,NX,micpar%mid_AutoNitriteOxidBacter,MicbE,isauto=.true.)
        this%h2D_NO2OxiBactC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_NO2OxiBactN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_NO2OxiBactP_vr(ncol,L) = micBE(ielmp)/DVOLL

        call SumMicbGroup(L,NY,NX,micpar%mid_AutoAeroCH4OxiBacter,MicbE,isauto=.true.)
        this%h2D_CH4AeroOxiC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_CH4AeroOxiN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_CH4AeroOxiP_vr(ncol,L) = micBE(ielmp)/DVOLL

        this%h2D_tRespGrossHeterUlm_vr(ncol,L)= tRespGrossHeterUlm_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_tRespGrossHeter_vr(ncol,L)= tRespGrossHeter_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

        call SumMicbGroup(L,NY,NX,micpar%mid_AutoH2GenoCH4GenArchea,MicbE,isauto=.true.)
        this%h2D_H2MethogenC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_H2MethogenN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_H2MethogenP_vr(ncol,L) = micBE(ielmp)/DVOLL

        this%h2D_hydrogMeth_frac_vr(ncol,L)=safe_adb(this%h2D_H2MethogenC_vr(ncol,L),this%h2D_microbC_vr(ncol,L))

        this%h2D_TSolidOMActC_vr(ncol,L)     = TSolidOMActC_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_TSolidOMActCDens_vr(ncol,L) = safe_adb(TSolidOMActC_vr(L,NY,NX),TSolidOMC_vr(L,NY,NX))
        this%h2D_tOMActCDens_vr(ncol,L)      = safe_adb(tOMActC_vr(L,NY,NX),(tOMActC_vr(L,NY,NX)+TSolidOMC_vr(L,NY,NX)))
        this%h2D_RCH4ProdAcetcl_vr(ncol,L)   = RCH4ProdAcetcl_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RCH4ProdHydrog_vr(ncol,L)   = RCH4ProdHydrog_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

        this%h2D_RCH4Oxi_aero_vr(ncol,L) = RCH4Oxi_aero_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RCH4Oxi_aero_col(ncol)  = this%h1D_RCH4Oxi_aero_col(ncol) + RCH4Oxi_aero_vr(L,NY,NX)
        this%h2D_RCH4Oxi_anmo_vr(ncol,L) = RCH4Oxi_anmo_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RCH4Oxi_anmo_col(ncol)  = this%h1D_RCH4Oxi_anmo_col(ncol) + RCH4Oxi_anmo_vr(L,NY,NX)

        this%h2D_RFerment_vr(ncol,L) = RFerment_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RDen_NO3toNO2_vr(ncol,L) = RDen_NO3toNO2_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RNit_NO2toNO3_vr(ncol,L)   = RNit_NO2toNO3_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_nh3oxi_vr(ncol,L)   = RNit_NH3toNO2_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_N2Oprod_vr(ncol,L)  = (RDen_NO2toN2O_vr(L,NY,NX)+RN2ONitProd_vr(L,NY,NX) &
                               +RN2OChemoProd_vr(L,NY,NX)-RDen_N2OtoN2_vr(L,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RootAR_vr(ncol,L) = -RootCO2Autor_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_PAR_RAD_vr(ncol,L)=PAR_RAD_vr(L,NY,NX)
        this%h2D_RootAR2soil_vr(ncol,L)=-RootCO2Ar2Soil_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RootAR2Root_vr(ncol,L)=-RootCO2Ar2Root_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        if(plant_model)then
          this%h2D_RootMassC_vr(ncol,L)     = RootMycoMassElm_vr(ielmc,ipltroot,L,NY,NX)/DVOLL
          this%h2D_RootMassN_vr(ncol,L)     = RootMycoMassElm_vr(ielmn,ipltroot,L,NY,NX)/DVOLL
          this%h2D_RootMassP_vr(ncol,L)     = RootMycoMassElm_vr(ielmp,ipltroot,L,NY,NX)/DVOLL
        endif
      ENDDO
      this%h1D_cyanoBactC_col(ncol) = this%h1D_cyanoBactC_col(ncol)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RCH4Oxi_aero_col(ncol) = this%h1D_RCH4Oxi_aero_col(ncol)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RCH4Oxi_anmo_col(ncol) = this%h1D_RCH4Oxi_anmo_col(ncol)/AREA_3D(3,NU_col(NY,NX),NY,NX)

  end procedure update_hist_soil

end submodule HistUpdateSoil
