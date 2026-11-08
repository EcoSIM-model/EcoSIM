submodule (HistDataType) HistUpdatePlant
  ! Populate plant and root history diagnostics from current model state.
  use GridConsts, only: JZ, MaxNodesPerBranch, MaxNumRootAxes, NumCanopyLayers, NumGrowthStages, &
    NumOfPlantMorphUnits
  use ElmIDMod, only: ielmc, ielmn, ielmp, imycorr_arbu, imycorrhz, ipltroot, iTrue, NumPlantChemElms
  use EcosimConst, only: natomw, BndlResistCO2ToH2ORatio
  use data_const_mod, only: spval  => DAT_CONST_SPVAL
  use MiniMathMod, only: safe_adb, AZERO
  use GridMod, only: get_pft
  use GridDataType, only: AREA_3D, DLYR_3D, NU_col
  use EcoSIMCtrlDataType, only: ZEROS
  use EcoSIMCtrlMod, only: ldo_sp_mode
  use FlagDataType, only: IsPlantActive_pft
  use PlantTraitDataType, only: CanopyHeightLive_pft, CanopyLeafArea_pft, CanopySeedNum_pft, &
    CanPhenoMoistStress_pft, CanPhenoTempStress_pft, ClumpFactorNow_pft, fTCanopyGroth_pft, &
    HoursTooLowPsiCan_pft, HypocotHeight_pft, iPlantCalendar_brch, LeafStalkAreaAct_pft, MainBranchNum_pft, &
    NumOfBranches_pft, NumOfLeaves_brch, PARTS_brch, PlantO2Stress_pft, PlantPopuLive_pft, PSICanPDailyMin_pft, &
    SeedDepth_pft, ShootNodeNum_brch, StalkAveRadius_pft, StalkHeight_pft, TCGroth_pft
  use PlantDataRateType, only: CanopyGrosRCO2_pft, CH4ByFire_CumYr_pft, CO2ByFire_CumYr_pft, &
    EcoHavstElmnt_CumYr_pft, ETCanopy_CumYr_pft, fRootGrowPSISense_pvr, GrossCO2Fix_CumYr_pft, GrossCO2Fix_pft, &
    GrossResp_pft, idg_CO2, idg_O2, ids_H2PO4, ids_H2PO4B, ids_NH4, ids_NH4B, ids_NO3, ids_NO3B, &
    LitrfalStrutElms_CumYr_pft, NetPrimProduct_pft, NH3byFire_CumYr_pft, NH3Dep2Can_pft, NH3Emis_CumYr_pft, &
    PlantElmBalCum_pft, PlantExudElm_CumYr_pft, PlantN2Fix_CumYr_pft, PO4byFire_CumYr_pft, PTSHTR_pft, &
    RAutoRootO2Limter_rpvr, RNFixCO2_pft, RootCO2Autor_col, RootCO2Autor_pvr, RootCO2Emis2Root_col, &
    RootH2PO4Uptake_pft, RootN2Fix_pft, RootNH4Uptake_pft, RootNO3Uptake_pft, RootNutUptake_pvr, &
    RootUptk_N_CumYr_pft, RootUptk_Nmin_cumYr_pft, RootUptk_P_CumYr_pft, RootUptk_Pmin_cumYr_pft, &
    ShootRootXferElm_pft, SurfLitrfalStrutElms_CumYr_pft, VmaxNH4Root_pvr, VmaxNO3Root_pvr
  use ClimForcDataType, only: RadPARSolarBeam_col
  use CanopyDataType, only: canopy_growth_pft, CanopyEvapTransLHeat_pft, CanopyGasCO2_pft, CanopyLeafAreaZ_pft, &
    CanopyMinStomaResistH2O_pft, CanopyNLimFactor_brch, CanopyNonstElmConc_pft, CanopyNonstOMConc_pft, CanopyNonstElms_pft, &
    CanopyPLimFactor_brch, CanopyVcMaxPEP25C_pft, CanopyVcMaxRubisco25C_pft, CanopyVoMaxRubisco25C_pft, &
    CanPStomaResistH2O_pft, CdH2ORootxSoil_pft, CH2OSunlit_pft, CH2OSunsha_pft, CO2FixCL_pft, CO2FixLL_pft, &
    CO2NetFix_pft, DynCi2CaRatio_pft, EarStrutElms_pft, ElectronTransptJmax25C_pft, fNCLFW_pft, fPCLFW_pft, &
    FracPARads2Canopy_pft, fSnowCanopy_pft, GrainStrutElms_pft, HeatStorCanopy_pft, HeatXAir2PCan_pft, &
    HuskStrutElms_pft, Km4RubiscoCarboxy_pft, LeafAreaSunlit_pft, LeafC3ChlCperm2LA_pft, LeafC4ChlCperm2LA_pft, &
    LeafPEPCperm2LA_pft, LeafProteinCperm2LA_pft, LeafRubiscoCperm2LA_pft, LeafStrutElms_pft, O2L_pft, &
    PARSunlit_pft, PARSunsha_pft, PetolShethStrutElms_pft, ProteinCperm2LeafArea_node, PSICanopy_pft, &
    PSICanopyTurg_pft, RadNet2Canopy_pft, RadPARCanopyAbsorption_pft, RawH2OCanopy2Atm_pft, RawHeatCanopy2Atm_pft, RCanMaintDef_CO2_pft, &
    RLeafAppear_pft, RNodeInitiate_pft, RubiscoActivity_brch, SeasonalNonstElms_pft, ShootElms_pft, &
    ShootNoduleElms_pft, SpecificLeafArea_pft, StalkRsrvElms_pft, StalkStrutElms_pft, StandDeadStrutElms_pft, &
    StomatalStress_pft, TdegCCanopy_pft, TFN_Carboxy_pft, TFN_eTranspt_pft, TFN_Oxygen_pft, Transpiration_pft, &
    TurgEff4CanopyResp_pft
  use RootDataType, only: CRootLumenArea_rpvr, Cytokinin1stConc_rpvr, fctyok_scalar_rpvr, fTgrowRootP_vr, &
    NumStructuralRootAxes_pft, Nutruptk_fClim_rpvr, Nutruptk_fNlim_rpvr, Nutruptk_fPlim_rpvr, &
    Nutruptk_fProtC_rpvr, PopuRootMycoC_pvr, PSIRoot_pvr, Root1stActStruct_pvr, Root1stDepz_raxes, &
    Root1stLenPP_rpvr, Root1stLigStruct_pvr, Root1stRadius_pvr, Root1stSinkWeight_pvr, &
    Root1stTipSinkWeight_pft, Root1stTransptArea_pvr, Root1stXNumL_pvr, Root2ndSinkWeight_pvr, &
    Root2ndXNumL_rpvr, Root2ndLen_rpvr, MeanSecondaryAxisLength_pft, &
    RootAtmGasConductance_rpvr, RootAxialKond2H2O_pvr, RootElms_pft, RootH2OUptkStress_pvr, &
    RootLenDensPerPlant_pvr, RootLenPerPlant_pvr, RootMaintDef_CO2_pvr, RootMediumLength_pvr, &
    RootMediumRadius_rpvr, RootMediumXNum_pvr, RootMedStruct_pvr, RootMSinkWeight_pvr, &
    RootMyco1stStrutElms_rpvr, RootMyco2ndStrutElms_pvr, RootMycoMassElm_pvr, RootMycoNonstElms_pft, &
    RootMycoNonstElms_rpvr, ROOTNLim_rpvr, RootNoduleElms_pft, RootNonstructElmConc_pvr, RootNonstOMConc_pvr, ROOTPLim_rpvr, &
    RootProteinC_pvr, RootRadialKond2H2O_pvr, RootResist4H2O_pvr, RootSAreaPerPlant_pvr, RootShootExch_pvr, &
    RootSinkWeight_pvr, RootStrutElms_pft, RPlantRootH2OUptk_pvr, SapFlowVlinear_pvr
  use SoilWaterDataType, only: QdewCanopy_CumYr_pft
  use PlantMgmtDataType, only: CanopyCutProxy_pft, NP0_col
  use TracerPropMod, only: GramPerHr2umolPerSec
  implicit none
contains

  module procedure update_hist_plants
    ! Fill instantaneous plant diagnostics for history averaging and output.
    integer :: nptc
    integer :: L
    integer :: NZ
    integer :: KN
    integer :: NB
    integer :: NR
    integer :: NB1
    integer :: K
    real(r8), parameter :: secs1hour=3600._r8
    real(r8), parameter :: MJ2W=1.e6_r8/secs1hour
    real(r8), parameter :: m2mm=1000._r8
    real(r8) :: DVOLL

      this%h1D_RootAR_col(ncol)       = -AZERO(RootCO2Autor_col(NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RootCO2Relez_col(ncol) = RootCO2Emis2Root_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1d_fPAR_col(ncol) = 0._r8
      this%h1D_QTRANSP_col(ncol)=0._r8

      DO NZ=1,NP0_col(NY,NX)
        nptc=get_pft(NZ,NY,NX)
        ! Reset live diagnostics before filling this step's branches and roots.
        ! Unused axes/layers must not retain a previous fallow year's zeros in
        ! one run and fresh-history fill values in a restarted run.
        if(PlantPopuLive_pft(NZ,NY,NX)>0._r8)call this%ZeroPlantHistVars(nptc)
        this%h1D_CanopyCutProxy_ptc(nptc)   = CanopyCutProxy_pft(NZ,I,NY,NX)
        this%h1D_POPN_ptc(nptc)             = PlantPopuLive_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_EXUD_CumYr_C_FLX_ptc(nptc) = PlantExudElm_CumYr_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_cumGPP_ptc(nptc)       = GrossCO2Fix_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_LITRf_C_FLX_ptc(nptc)      = LitrfalStrutElms_CumYr_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_SURF_LITRf_C_FLX_ptc(nptc) = SurfLitrfalStrutElms_CumYr_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_HVST_C_FLX_ptc(nptc)       = EcoHavstElmnt_CumYr_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_HVST_N_FLX_ptc(nptc)       = EcoHavstElmnt_CumYr_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_HVST_P_FLX_ptc(nptc)       = EcoHavstElmnt_CumYr_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_PLANT_BALANCE_C_ptc(nptc)  = PlantElmBalCum_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_PLANT_BALANCE_N_ptc(nptc)  = PlantElmBalCum_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_PLANT_BALANCE_P_ptc(nptc)  = PlantElmBalCum_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_FIREp_CO2_FLX_ptc(nptc)    = CO2ByFire_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_FIREp_CH4_FLX_ptc(nptc)    = CH4ByFire_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_NPP_ptc(nptc)              = NetPrimProduct_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_tTRANSPN_ptc(nptc)         = -ETCanopy_CumYr_pft(NZ,NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_EXUD_N_FLX_ptc(nptc)       = PlantExudElm_CumYr_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Uptk_NMIN_Flx_ptc(nptc)    = RootUptk_Nmin_cumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Uptk_PMIN_Flx_ptc(nptc)    = RootUptk_Pmin_cumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Uptk_N_Flx_ptc(nptc)       = RootUptk_N_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Uptk_P_Flx_ptc(nptc)       = RootUptk_P_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_LITRf_N_FLX_ptc(nptc)      = LitrfalStrutElms_CumYr_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_cum_N_FIXED_ptc(nptc)      = PlantN2Fix_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_NH3can_FLX_ptc(nptc)       = NH3Emis_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_FIREp_N_FLX_ptc(nptc)      = NH3byFire_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_SURF_LITRf_N_FLX_ptc(nptc) = SurfLitrfalStrutElms_CumYr_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CANDew_ptc(nptc)          = QdewCanopy_CumYr_pft(NZ,NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_EXUD_P_FLX_ptc(nptc)       = PlantExudElm_CumYr_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_LITRf_P_FLX_ptc(nptc)      = LitrfalStrutElms_CumYr_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_FIREp_P_FLX_ptc(nptc)      = PO4byFire_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_SURF_LITRf_P_FLX_ptc(nptc) = SurfLitrfalStrutElms_CumYr_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STANDING_DEAD_N_ptc(nptc)  = StandDeadStrutElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STANDING_DEAD_P_ptc(nptc)  = StandDeadStrutElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STANDING_DEAD_C_ptc(nptc)  = StandDeadStrutElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        if(MainBranchNum_pft(NZ,NY,NX)>0) then
          this%h1D_ShootNodeNumber_ptc(nptc)  = SUM(ShootNodeNum_brch(1:NumOfBranches_pft(NZ,NY,NX),NZ,NY,NX))
          this%h1D_MainBranchNodeNumber_ptc(nptc)=ShootNodeNum_brch(MainBranchNum_pft(NZ,NY,NX),NZ,NY,NX)
        else
          this%h1D_MainBranchNodeNumber_ptc(nptc)=0
          this%h1D_ShootNodeNumber_ptc(nptc)  = 0
        endif
        if (PlantPopuLive_pft(NZ,NY,NX) .LE. 0._r8)then
          call this%ZeroPlantHistVars(nptc)
          cycle
        endif
        this%h1D_PTSHTR_ptc(nptc) = PTSHTR_pft(NZ,NY,NX)
        this%h1D_fSnowCan_ptc(nptc)     = fSnowCanopy_pft(NZ,NY,NX)
        this%h1D_ROOT_NONSTC_ptc(nptc)  = RootMycoNonstElms_pft(ielmc,ipltroot,NZ,NY,NX)
        this%h1D_ROOT_NONSTN_ptc(nptc)  = RootMycoNonstElms_pft(ielmn,ipltroot,NZ,NY,NX)
        this%h1D_ROOT_NONSTP_ptc(nptc)  = RootMycoNonstElms_pft(ielmp,ipltroot,NZ,NY,NX)
        this%h1D_SHOOT_NONSTC_ptc(nptc) = CanopyNonstElms_pft(ielmc,NZ,NY,NX)
        this%h1D_SHOOT_NONSTN_ptc(nptc) = CanopyNonstElms_pft(ielmn,NZ,NY,NX)
        this%h1D_SHOOT_NONSTP_ptc(nptc) = CanopyNonstElms_pft(ielmp,NZ,NY,NX)

        if(CO2FixLL_pft(NZ,NY,NX)/=spval .and. CO2FixCL_pft(NZ,NY,NX)/=spval)then
          this%h1D_dCAN_GPP_CLIM_ptc(nptc) = AZERO(CO2FixCL_pft(NZ,NY,NX)-GrossCO2Fix_pft(NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
          this%h1D_dCAN_GPP_eLIM_ptc(nptc) = AZERO(CO2FixLL_pft(NZ,NY,NX)-GrossCO2Fix_pft(NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        endif

        this%h1D_LeafChlCperm2LA_ptc(nptc)   = (LeafC3ChlCperm2LA_pft(NZ,NY,NX)+LeafC4ChlCperm2LA_pft(NZ,NY,NX))*1.e3_r8
        this%h1D_LeafC4ChlCperm2LA_ptc(nptc)   = LeafC4ChlCperm2LA_pft(NZ,NY,NX)*1.e3_r8
        this%h1D_LeafRubiscoNperm2LA_ptc(nptc) = LeafRubiscoCperm2LA_pft(NZ,NY,NX)/3.3_r8
        this%h1D_LeafPEPCNperm2LA_ptc(nptc)    = LeafPEPCperm2LA_pft(NZ,NY,NX)/3.3_r8

        this%h1D_MIN_LWP_ptc(nptc)      = PSICanPDailyMin_pft(NZ,NY,NX)
        this%h1D_SLA_ptc(nptc)          = 1.e4_r8*SpecificLeafArea_pft(NZ,NY,NX)
        this%h1D_LEAF_PC_ptc(nptc)       = safe_adb(LeafStrutElms_pft(ielmp,NZ,NY,NX)+CanopyNonstElms_pft(ielmp,NZ,NY,NX), &
                                                 LeafStrutElms_pft(ielmc,NZ,NY,NX)+CanopyNonstElms_pft(ielmc,NZ,NY,NX))
        this%h1D_CAN_RN_ptc(nptc)        = MJ2W*RadNet2Canopy_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_LE_ptc(nptc)        = MJ2W*CanopyEvapTransLHeat_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_H_ptc(nptc)         = MJ2W*HeatXAir2PCan_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_G_ptc(nptc)         = MJ2W*HeatStorCanopy_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_TEMPC_ptc(nptc)     = TdegCCanopy_pft(NZ,NY,NX)
        this%h1D_CAN_TEMPFN_ptc(nptc)    = fTCanopyGroth_pft(NZ,NY,NX)
        this%h1D_CAN_CO2_FLX_ptc(nptc)   = CO2NetFix_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
        this%h1D_CAN_GPP_ptc(nptc)       = GrossCO2Fix_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_RA_ptc(nptc)        = CanopyGrosRCO2_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_TurgEff4CanopyResp_ptc(nptc) = TurgEff4CanopyResp_pft(NZ,NY,NX)
        this%h1D_CAN_GROWTH_ptc(nptc)    = canopy_growth_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_cTNC_ptc(nptc)          = CanopyNonstElmConc_pft(ielmc,NZ,NY,NX)
        this%h1D_cTNN_ptc(nptc)          = CanopyNonstElmConc_pft(ielmn,NZ,NY,NX)
        this%h1D_cTNP_ptc(nptc)          = CanopyNonstElmConc_pft(ielmp,NZ,NY,NX)
        this%h1D_CanopyNonstOMConc_ptc(nptc) = CanopyNonstOMConc_pft(NZ,NY,NX)
        this%h1D_STOML_RSC_CO2_ptc(nptc) = CanPStomaResistH2O_pft(NZ,NY,NX)*1.56_r8*secs1hour
        this%h1D_STOML_Min_RSC_CO2_ptc(nptc)=CanopyMinStomaResistH2O_pft(NZ,NY,NX)*1.56_r8*secs1hour
        this%h1D_Km_CO2_carboxy_ptc(nptc)= Km4RubiscoCarboxy_pft(NZ,NY,NX)
        this%h1D_DynCi2CaRatio_ptc(nptc) = DynCi2CaRatio_pft(NZ,NY,NX)

        this%h1D_BLYR_RSC_CO2_ptc(nptc)  = RawH2OCanopy2Atm_pft(NZ,NY,NX)*BndlResistCO2ToH2ORatio*secs1hour
        this%h1D_CAN_CO2_ptc(nptc)       = CanopyGasCO2_pft(NZ,NY,NX)
        this%h1D_O2L_ptc(nptc)           = O2L_pft(NZ,NY,NX)
        this%h1D_LAI_ptc(nptc)           = LeafStalkAreaAct_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CanPhenol_WSTRSS_ptc(nptc)=CanPhenoMoistStress_pft(NZ,NY,NX)
        this%h1D_CanPhenol_TSTRSS_ptc(nptc)=CanPhenoTempStress_pft(NZ,NY,NX)
        this%h1D_PSI_CAN_ptc(nptc)       = PSICanopy_pft(NZ,NY,NX)
        this%h1D_TURG_CAN_ptc(nptc)      = PSICanopyTurg_pft(NZ,NY,NX)
        this%h1D_STOML_RSC_H2O_ptc(nptc)  = CanPStomaResistH2O_pft(NZ,NY,NX)*secs1hour
        this%h1D_BLYR_RSC_H2O_ptc(nptc)  = RawH2OCanopy2Atm_pft(NZ,NY,NX)*secs1hour
        this%h1D_BLYR_RSC_HEAT_ptc(nptc)  = RawHeatCanopy2Atm_pft(NZ,NY,NX)*secs1hour
        this%h1D_CdH2ORootxSoil_ptc(nptc) = CdH2ORootxSoil_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_TRANSPN_ptc(nptc)       = Transpiration_pft(NZ,NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_QTRANSP_col(ncol)       = this%h1D_QTRANSP_col(ncol)+this%h1D_TRANSPN_ptc(nptc)
        this%h1D_NH4_UPTK_FLX_ptc(nptc)  = RootNH4Uptake_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_NO3_UPTK_FLX_ptc(nptc)  = RootNO3Uptake_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_N2_FIXN_FLX_ptc(nptc)   = RootN2Fix_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_cNH3_FLX_ptc(nptc)      = NH3Dep2Can_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RNFixCO2_ptc(nptc)      =RNFixCO2_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_TC_Groth_ptc(nptc)      = TCGroth_pft(NZ,NY,NX)
        this%h1D_PO4_UPTK_FLX_ptc(nptc)  = RootH2PO4Uptake_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_frcPARabs_ptc(nptc)     = FracPARads2Canopy_pft(NZ,NY,NX)
        this%h1D_PAR_CAN_ptc(nptc)       = RadPARCanopyAbsorption_pft(NZ,NY,NX)   !umol /m2/s
        this%h1d_fPAR_col(ncol)          = this%h1d_fPAR_col(ncol)+RadPARCanopyAbsorption_pft(NZ,NY,NX)
        this%h1D_SHOOT_C_ptc(nptc)       = ShootElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Plant_C_ptc(nptc)       = (ShootElms_pft(ielmc,NZ,NY,NX) &
          +RootElms_pft(ielmc,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RNodeInitiate_ptc(nptc) = RNodeInitiate_pft(NZ,NY,NX)*24._r8
        this%h1D_RLeafAppear_ptc(nptc) = RLeafAppear_pft(NZ,NY,NX)*24._r8
        this%h1D_LEAF_C_ptc(nptc)        = LeafStrutElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Petole_C_ptc(nptc)      = PetolShethStrutElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STALK_C_ptc(nptc)       = StalkStrutElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RESERVE_C_ptc(nptc)     = StalkRsrvElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_HUSK_C_ptc(nptc)        = (HuskStrutElms_pft(ielmc,NZ,NY,NX) &
          +EarStrutElms_pft(ielmc,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_GRAIN_C_ptc(nptc)       = GrainStrutElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ROOT_C_ptc(nptc)        = RootElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ROOTST_C_ptc(nptc)      = RootStrutElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ROOTST_N_ptc(nptc)      = RootStrutElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ROOTST_P_ptc(nptc)      = RootStrutElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RootNodule_C_ptc(nptc)  = RootNoduleElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ShootNodule_C_ptc(nptc)  = ShootNoduleElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ShootNodule_N_ptc(nptc)  = ShootNoduleElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ShootNodule_P_ptc(nptc)  = ShootNoduleElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

        this%h1D_STORED_C_ptc(nptc)      = SeasonalNonstElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_GRAIN_NO_ptc(nptc)      = CanopySeedNum_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_LAIb_ptc(nptc)          = CanopyLeafArea_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_AUTO_RESP_FLX_ptc(nptc)    = GrossResp_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

        this%h1D_CAN_HT_ptc(nptc)            = CanopyHeightLive_pft(NZ,NY,NX)
        this%h1D_Stalk_HT_ptc(nptc)          = StalkHeight_pft(NZ,NY,NX)
        this%h1D_EmergeHeight_ptc(nptc)     = HypocotHeight_pft(NZ,NY,NX)-SeedDepth_pft(NZ,NY,NX)
        this%h1D_WTR_STRESS_ptc(nptc)        = HoursTooLowPsiCan_pft(NZ,NY,NX)
        this%h1D_LeafProteinNperm2_ptc(nptc) = LeafProteinCperm2LA_pft(NZ,NY,NX)/3.3_r8
        this%h1D_VcMaxRubisco_ptc(nptc)      = CanopyVcMaxRubisco25C_pft(NZ,NY,NX)
        this%h1D_VoMaxRubisco_ptc(nptc)      = CanopyVoMaxRubisco25C_pft(NZ,NY,NX)
        this%h1D_VcMaxPEP_ptc(nptc)          = CanopyVcMaxPEP25C_pft(NZ,NY,NX)
        this%h1D_JMaxPhoto_ptc(nptc)         = ElectronTransptJmax25C_pft(NZ,NY,NX)
        this%h1D_TFN_Carboxy_ptc(nptc)       = TFN_Carboxy_pft(NZ,NY,NX)
        this%h1D_TFN_Oxygen_ptc(nptc)        = TFN_Oxygen_pft(NZ,NY,NX)
        this%h1D_TFN_eTranspt_ptc(nptc)      = TFN_eTranspt_pft(NZ,NY,NX)
        this%h1D_PARSunlit_ptc(nptc)         = PARSunlit_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_PARSunsha_ptc(nptc)         = PARSunsha_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CH2OSunlit_ptc(nptc)        = CH2OSunlit_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CH2OSunsha_ptc(nptc)        = CH2OSunsha_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

        this%h1D_fClump_ptc(nptc)  = ClumpFactorNow_pft(NZ,NY,NX)
        this%h1D_LeafAreaSunlit_ptc(nptc)=LeafAreaSunlit_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_OXY_STRESS_ptc(nptc)   = PlantO2Stress_pft(NZ,NY,NX)
        this%h1D_SHOOT_N_ptc(nptc)      = ShootElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Plant_N_ptc(nptc)      = (ShootElms_pft(ielmn,NZ,NY,NX)&
          +RootElms_pft(ielmn,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_fCNLFW_ptc(nptc) = safe_adb(1._r8,fNCLFW_pft(NZ,NY,NX))
        this%h1D_fCPLFW_ptc(nptc) = safe_adb(1._r8,fPCLFW_pft(NZ,NY,NX))
        this%h1D_LEAF_N_ptc(nptc)    = LeafStrutElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_LeafNperm2LAI_ptc(nptc) = safe_adb(LeafStrutElms_pft(ielmn,NZ,NY,NX),CanopyLeafArea_pft(NZ,NY,NX))
        this%h1D_Petole_N_ptc(nptc)  = PetolShethStrutElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STALK_N_ptc(nptc)   = StalkStrutElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RESERVE_N_ptc(nptc) = StalkRsrvElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_HUSK_N_ptc(nptc)    = (HuskStrutElms_pft(ielmn,NZ,NY,NX) &
          +EarStrutElms_pft(ielmn,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_GRAIN_N_ptc(nptc)          = GrainStrutElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ROOT_N_ptc(nptc)           = RootElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RootNodule_N_ptc(nptc)     = RootNoduleElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STORED_N_ptc(nptc)         = SeasonalNonstElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_TreeRingRadius_ptc(nptc)   = StalkAveRadius_pft(NZ,NY,NX)
        this%h1D_Root1stTipSinkWt_ptc(nptc) = Root1stTipSinkWeight_pft(NZ,NY,NX)
        this%h1D_SHOOT_P_ptc(nptc)          = ShootElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Plant_P_ptc(nptc)          = (ShootElms_pft(ielmp,NZ,NY,NX) &
          +RootElms_pft(ielmp,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_stomatal_stress_ptc(nptc) = StomatalStress_pft(NZ,NY,NX)
        this%h1D_LEAF_P_ptc(nptc)          = LeafStrutElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Petole_P_ptc(nptc)        = PetolShethStrutElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STALK_P_ptc(nptc)         = StalkStrutElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RESERVE_P_ptc(nptc)       = StalkRsrvElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_HUSK_P_ptc(nptc)          = (HuskStrutElms_pft(ielmp,NZ,NY,NX) &
          +EarStrutElms_pft(ielmp,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_GRAIN_P_ptc(nptc)          = GrainStrutElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ROOT_P_ptc(nptc)           = RootElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RootNodule_P_ptc(nptc)     = RootNoduleElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STORED_P_ptc(nptc)         = SeasonalNonstElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ShootRootXferC_ptc(nptc)   = ShootRootXferElm_pft(ielmc,NZ,NY,NX)/PlantPopuLive_pft(NZ,NY,NX)
        this%h1D_ShootRootXferN_ptc(nptc)   = ShootRootXferElm_pft(ielmn,NZ,NY,NX)/PlantPopuLive_pft(NZ,NY,NX)
        this%h1D_ShootRootXferP_ptc(nptc)   = ShootRootXferElm_pft(ielmp,NZ,NY,NX)/PlantPopuLive_pft(NZ,NY,NX)
        this%h1D_BRANCH_NO_ptc(nptc)        = NumOfBranches_pft(NZ,NY,NX)
        this%h1D_MainBranchNO_ptc(nptc)     = MainBranchNum_pft(NZ,NY,NX)
        this%h1D_RCanMaintDef_CO2_pft(nptc) = RCanMaintDef_CO2_pft(NZ,NY,NX)
        this%h1D_LEAF_NC_ptc(nptc)      = safe_adb(LeafStrutElms_pft(ielmn,NZ,NY,NX)+CanopyNonstElms_pft(ielmn,NZ,NY,NX),&
                                                 LeafStrutElms_pft(ielmc,NZ,NY,NX)+CanopyNonstElms_pft(ielmc,NZ,NY,NX))
        this%h1D_RootMaintDef_CO2_pft(nptc) = sum(RootMaintDef_CO2_pvr(ipltroot,1:JZ,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_NumPrimeRootAxes_ptc(nptc)      = NumStructuralRootAxes_pft(NZ,NY,NX)

        IF(NumOfBranches_pft(NZ,NY,NX)>0)then
          DO K=1,MaxNodesPerBranch
            this%h2D_ProteinNperm2LeafArea_pnd(nptc,K)=0._r8
            DO NB=1,NumOfBranches_pft(NZ,NY,NX)
              this%h2D_ProteinNperm2LeafArea_pnd(nptc,K)=this%h2D_ProteinNperm2LeafArea_pnd(nptc,K)+ProteinCperm2LeafArea_node(K,NB,NZ,NY,NX)
            ENDDO
            !assuming protein C to N mass ratio is 3.3 (median value)
            this%h2D_ProteinNperm2LeafArea_pnd(nptc,K)=this%h2D_ProteinNperm2LeafArea_pnd(nptc,K)/(NumOfBranches_pft(NZ,NY,NX)*3.3_r8)
          ENDDO
        ENDIF
        this%h1D_Growth_Stage_ptc(nptc) = 0
        this%h1D_RUB_ACTVN_ptc(nptc)    = 0._r8
        this%h1D_CanopyNLim_ptc(nptc)   = 0._r8
        this%h1D_CanopyPLim_ptc(nptc)   = 0._r8
        this%h1D_Num_Leaves_ptc(nptc)   = 0._r8
        if(MainBranchNum_pft(NZ,NY,NX)> 0)then
          DO KN=NumGrowthStages,0,-1
            IF(KN>0)THEN
              if(any(iPlantCalendar_brch(KN,:,NZ,NY,NX)>0))then
                this%h1D_Growth_Stage_ptc(nptc) =KN
                exit
              endif
            ENDIF
          ENDDO
          NB1=0
          DO NB=1,NumOfBranches_pft(NZ,NY,NX)
            this%h1D_Num_Leaves_ptc(nptc)  = this%h1D_Num_Leaves_ptc(nptc)+NumOfLeaves_brch(NB,NZ,NY,NX)
            if(RubiscoActivity_brch(NB,NZ,NY,NX)>0._r8)then
              NB1=NB1+1
              this%h1D_RUB_ACTVN_ptc(nptc)   = this%h1D_RUB_ACTVN_ptc(nptc)+ RubiscoActivity_brch(NB,NZ,NY,NX)
              this%h1D_CanopyNLim_ptc(nptc)  = this%h1D_CanopyNLim_ptc(nptc)+ CanopyNLimFactor_brch(NB,NZ,NY,NX)
              this%h1D_CanopyPLim_ptc(nptc)  = this%h1D_CanopyPLim_ptc(nptc)+ CanopyPLimFactor_brch(NB,NZ,NY,NX)
            endif
            this%h3D_PARTS_ptc(nptc,1:NumOfPlantMorphUnits,NB) = PARTS_brch(1:NumOfPlantMorphUnits,NB,NZ,NY,NX)
          ENDDO
          if(NB1>0)then
            this%h1D_RUB_ACTVN_ptc(nptc)=this%h1D_RUB_ACTVN_ptc(nptc)/real(NB1,kind=r8)
            this%h1D_CanopyNLim_ptc(nptc)=this%h1D_CanopyNLim_ptc(nptc)/real(NB1,kind=r8)
            this%h1D_CanopyPLim_ptc(nptc)=this%h1D_CanopyPLim_ptc(nptc)/real(NB1,kind=r8)
          endif
        endif

        DO L=1,NumCanopyLayers
          this%h2D_RootShootExchC_pvr(nptc,L)=RootShootExch_pvr(ielmc,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
          this%h2D_RootShootExchN_pvr(nptc,L)=RootShootExch_pvr(ielmn,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
          this%h2D_RootShootExchP_pvr(nptc,L)=RootShootExch_pvr(ielmp,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
          this%h2D_CanopyLAIZ_plyr(nptc,L)=CanopyLeafAreaZ_pft(L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        ENDDO
        this%h1D_RootAR_ptc(nptc)          = 0._r8
        this%h1D_RootLenPerPlant_ptc(nptc) = 0._r8
        if(IsPlantActive_pft(NZ,NY,NX).EQ.iTrue .and. PlantPopuLive_pft(NZ,NY,NX) .GT. ZEROS(NY,NX))then
          DO NR=1,NumStructuralRootAxes_pft(NZ,NY,NX)
            this%h2D_Root1stDepz_ptc(nptc,NR)      = Root1stDepz_raxes(NR,NZ,NY,NX)
          ENDDO
        else
          DO NR=1,MaxNumRootAxes
            this%h2D_Root1stDepz_ptc(nptc,NR)      = 0._r8
          ENDDO
        endif
        this%h1D_RootAbsorbAreaPP_pft(nptc)=0._r8
        this%h1D_MycorrizhalBiomC_ptc(nptc) = 0._r8
        this%h1D_Root1stStrutC_ptc(nptc)=0._r8
        this%h1D_RootMeDStrutC_ptc(nptc)=0._r8
        this%h1D_Root1stStrutN_ptc(nptc)=0._r8
        this%h1D_Root2ndStrutC_ptc(nptc)=0._r8
        this%h1D_RootAct1stC_ptc(nptc)=0._r8
        DO L=1,JZ
          this%h1D_RootAR_ptc(nptc)=this%h1D_RootAR_ptc(nptc)-RootCO2Autor_pvr(ipltroot,L,NZ,NY,NX)
          DVOLL                                  = DLYR_3D(3,L,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
          this%h2D_Cytokinin1stConc_pvr(nptc,L) = 0._r8
          this%h2D_Cytok_scalar_pvr(nptc,L)      = 0._r8
          this%h2D_Root1stStrutC_pvr(nptc,L)     = 0._r8
          this%h2D_MycoBiomC_pvr(nptc,L)         = 0._r8
          this%h2D_CRootLumenArea_pvr(nptc,L)    = 0._R8
          this%h2D_Root1stStrutN_pvr(nptc,L)     = 0._r8
          this%h2D_Root1stStrutP_pvr(nptc,L)     = 0._r8
          this%h2D_Root2ndStrutC_pvr(nptc,L)     = 0._r8
          this%h2D_Root2ndStrutN_pvr(nptc,L)     = 0._r8
          this%h2D_Root2ndStrutP_pvr(nptc,L)     = 0._r8
          this%h2D_Root2ndMeanLen_pvr(nptc,L)    = 0._r8
          this%h2D_RootAct1stC_pvr(nptc,L)       = 0._r8
          this%h2D_RootLig1stC_pvr(nptc,L)       = 0._r8
          this%h2D_NonstC_conc_pvr(nptc,L)       = RootNonstructElmConc_pvr(ielmc,ipltroot,L,NZ,NY,NX)
          if(DVOLL>1.e-8_r8)then
            this%h2d_RootPop_pvr(nptc,L)=PopuRootMycoC_pvr(ipltroot,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
            this%h2D_MycoPop_pvr(nptc,L)=PopuRootMycoC_pvr(imycorrhz,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
            this%h2D_RootRadialKond2H2O_pvr(nptc,L)=RootRadialKond2H2O_pvr(ipltroot,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
            this%h2D_RootAXialKond2H2O_pvr(nptc,L) =RootAxialKond2H2O_pvr(ipltroot,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
            this%h2D_VmaxNH4Root_pvr(nptc,L)       = VmaxNH4Root_pvr(ipltroot,L,NZ,NY,NX)*1.e6/natomw
            this%h2D_VmaxNO3Root_pvr(nptc,L)       = VmaxNO3Root_pvr(ipltroot,L,NZ,NY,NX)*1.e6/natomw
            this%h2D_RootMassC_pvr(nptc,L)         = RootMycoMassElm_pvr(ielmc,ipltroot,L,NZ,NY,NX)/DVOLL
            this%h2D_RootNutupk_fClim_pvr(nptc,L)  = Nutruptk_fClim_rpvr(ipltroot,L,NZ,NY,NX)
            this%h2D_RootNutupk_fNlim_pvr(nptc,L)  = Nutruptk_fNlim_rpvr(ipltroot,L,NZ,NY,NX)
            this%h2D_RootNutupk_fPlim_pvr(nptc,L)  = Nutruptk_fPlim_rpvr(ipltroot,L,NZ,NY,NX)
            this%h2D_RootNutupk_fProtC_pvr(nptc,L) = Nutruptk_fProtC_rpvr(ipltroot,L,NZ,NY,NX)
            this%h2D_Root1stSArea4GasTP_pvr(nptc,L)= Root1stTransptArea_pvr(ipltroot,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
            this%h2D_RootProteinC_pvr(nptc,L)      = RootProteinC_pvr(ipltroot,L,NZ,NY,NX)/DVOLL
            this%h2D_O2_rootconduct_pvr(nptc,L)    = RootAtmGasConductance_rpvr(idg_O2,ipltroot,L,NZ,NY,NX)
            this%h2D_CO2_rootconduct_pvr(nptc,L)   = RootAtmGasConductance_rpvr(idg_CO2,ipltroot,L,NZ,NY,NX)
            this%h2D_fTRootGro_pvr(nptc,L)         = fTgrowRootP_vr(L,NZ,NY,NX)
            this%h2D_fRootGrowPSISense_pvr(nptc,L) = fRootGrowPSISense_pvr(ipltroot,L,NZ,NY,NX)
            this%h1D_RootAct1stC_ptc(nptc) = this%h1D_RootAct1stC_ptc(nptc)+Root1stActStruct_pvr(ielmc,L,NZ,NY,NX)
            this%h2D_RootAbsorbAreaPP_pvr(nptc,L)  = RootSAreaPerPlant_pvr(ipltroot,L,NZ,NY,NX)
            this%h1D_RootAbsorbAreaPP_pft(nptc)    = this%h1D_RootAbsorbAreaPP_pft(nptc)+RootSAreaPerPlant_pvr(ipltroot,L,NZ,NY,NX)
            this%h2D_ROOT_OSTRESS_pvr(nptc,L)      = RAutoRootO2Limter_rpvr(ipltroot,L,NZ,NY,NX)
            this%h2D_PSI_RT_pvr(nptc,L)            = PSIRoot_pvr(ipltroot,L,NZ,NY,NX)
            this%h2D_RootH2OUptkStress_pvr(nptc,L) = RootH2OUptkStress_pvr(ipltroot,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
            this%h2D_SapFlowVlinear_pvr(nptc,L) = (SapFlowVlinear_pvr(L,NZ,NY,NX))
            this%h2D_RootH2OUptk_pvr(nptc,L) = AZERO(m2mm*RPlantRootH2OUptk_pvr(ipltroot,L,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
            this%h2D_RootMaintDef_CO2_pvr(nptc,L)=RootMaintDef_CO2_pvr(ipltroot,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
            this%h2D_prtUP_NH4_pvr(nptc,L)    = (sum(RootNutUptake_pvr(ids_NH4,:,L,NZ,NY,NX))+&
              sum(RootNutUptake_pvr(ids_NH4B,:,L,NZ,NY,NX)))/AREA_3D(3,L,NY,NX)
            this%h2D_prtUP_NO3_pvr(nptc,L)  = (sum(RootNutUptake_pvr(ids_NO3,:,L,NZ,NY,NX))+&
              sum(RootNutUptake_pvr(ids_NO3B,:,L,NZ,NY,NX)))/AREA_3D(3,L,NY,NX)
            this%h2D_prtUP_PO4_pvr(nptc,L)  = (sum(RootNutUptake_pvr(ids_H2PO4,:,L,NZ,NY,NX))+&
              sum(RootNutUptake_pvr(ids_H2PO4B,:,L,NZ,NY,NX)))/AREA_3D(3,L,NY,NX)
            this%h2D_DNS_RT_pvr(nptc,L) = RootLenDensPerPlant_pvr(ipltroot,L,NZ,NY,NX)* &
              PlantPopuLive_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*1.e-4_r8
            this%h1D_RootLenPerPlant_ptc(nptc)=this%h1D_RootLenPerPlant_ptc(nptc)+RootLenPerPlant_pvr(ipltroot,L,NZ,NY,NX)

            this%h2D_ROOTNLim_rpvr(nptc,L)         = ROOTNLim_rpvr(ipltroot,L,NZ,NY,NX)
            this%h2D_ROOTPLim_rpvr(nptc,L)         = ROOTPLim_rpvr(ipltroot,L,NZ,NY,NX)
            this%h2D_RootNonstC_rpvr(nptc,L)       = RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ,NY,NX)
            this%h2D_RootSinkWeight_pvr(nptc,L)    = RootSinkWeight_pvr(L,NZ,NY,NX)
            this%h2D_RootMSinkWeight_pvr(nptc,L)   = RootMSinkWeight_pvr(L,NZ,NY,NX)
            this%h2D_Root2ndSinkWeight_pvr(nptc,L) = Root2ndSinkWeight_pvr(L,ipltroot,NZ,NY,NX)
            this%h2D_Root1stSinkWeight_pvr(nptc,L) = Root1stSinkWeight_pvr(L,NZ,NY,NX)
            this%h2D_Root1stRadius_rpvr(nptc,L)    = Root1stRadius_pvr(ipltroot,L,NZ,NY,NX)*m2mm
            this%h2D_RootNonstBConc_pvr(nptc,L)    = RootNonstOMConc_pvr(ipltroot,L,NZ,NY,NX)

            this%h2D_Rootmedlength_pvr(nptc,L) = RootMediumLength_pvr(L,NZ,NY,NX)*RootMediumXNum_pvr(L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
            if(PlantPopuLive_pft(NZ,NY,NX).GT.0._r8)then
              this%h2D_Root1stAxesNumL_pvr(nptc,L)= Root1stXNumL_pvr(L,NZ,NY,NX)/PlantPopuLive_pft(NZ,NY,NX)
              this%h2D_RootMedAxesNumL_pvr(nptc,L)=RootMediumXNum_pvr(L,NZ,NY,NX)/PlantPopuLive_pft(NZ,NY,NX)
            else
              this%h2D_Root1stAxesNumL_pvr(nptc,L)= 0._r8
              this%h2D_RootMedAxesNumL_pvr(nptc,L)=0._r8
            endif
            this%h2D_Root2ndAxesNumL_pvr(nptc,L)= Root2ndXNumL_rpvr(ipltroot,L,NZ,NY,NX)
            if(Root2ndXNumL_rpvr(ipltroot,L,NZ,NY,NX)>0._r8)then
              if(ldo_sp_mode)then
                this%h2D_Root2ndMeanLen_pvr(nptc,L)=MeanSecondaryAxisLength_pft(NZ,NY,NX)
              else
                ! Count-weighted mean from actual secondary-root geometry;
                ! unlike hydraulic effective length, this has no tip-length floor.
                this%h2D_Root2ndMeanLen_pvr(nptc,L)= &
                  sum(Root2ndLen_rpvr(ipltroot,L,1:NumStructuralRootAxes_pft(NZ,NY,NX),NZ,NY,NX)) &
                  /Root2ndXNumL_rpvr(ipltroot,L,NZ,NY,NX)
              endif
            endif
            this%h2D_RootKond2H2O_pvr(nptc,L)= safe_adb(1._r8,RootResist4H2O_pvr(ipltroot,L,NZ,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX))*1.e7/3600._r8
            this%h2D_RootMedC_pvr(nptc,L) = RootMedStruct_pvr(ielmc,L,NZ,NY,NX)/DVOLL
            this%h2D_RootAct1stC_pvr(nptc,L) = Root1stActStruct_pvr(ielmc,L,NZ,NY,NX)/DVOLL
            this%h2D_RootLig1stC_pvr(nptc,L) = Root1stLigStruct_pvr(ielmc,L,NZ,NY,NX)/DVOLL
            this%h2D_Root1stLenPP_pvr(nptc,L)   = 0._r8
            this%h2D_RootmedRadius_pvr(nptc,L)  = 0._r8
            this%h2D_Root2ndStrutC_pvr(nptc,L)  = RootMyco2ndStrutElms_pvr(ielmc,ipltroot,L,NZ,NY,NX)
            this%h2D_Root2ndStrutN_pvr(nptc,L)  = RootMyco2ndStrutElms_pvr(ielmn,ipltroot,L,NZ,NY,NX)
            this%h2D_Root2ndStrutP_pvr(nptc,L)  = RootMyco2ndStrutElms_pvr(ielmp,ipltroot,L,NZ,NY,NX)
            this%h2D_MycoBiomC_pvr(nptc,L)      = RootMyco2ndStrutElms_pvr(ielmc,imycorr_arbu,L,NZ,NY,NX)
            this%h1D_MycorrizhalBiomC_ptc(nptc) = this%h1D_MycorrizhalBiomC_ptc(nptc)+this%h2D_MycoBiomC_pvr(nptc,L)
            DO NR=1,NumStructuralRootAxes_pft(NZ,NY,NX)
              this%h2D_RootmedRadius_pvr(nptc,L)=this%h2D_RootmedRadius_pvr(nptc,L)+RootMediumRadius_rpvr(L,NR,NZ,NY,NX)
              this%h2D_Root1stLenPP_pvr(nptc,L)=this%h2D_Root1stLenPP_pvr(nptc,L)+Root1stLenPP_rpvr(L,NR,NZ,NY,NX)
              this%h2D_CRootLumenArea_pvr(nptc,L)    = this%h2D_CRootLumenArea_pvr(nptc,L)+CRootLumenArea_rpvr(L,NR,NZ,NY,NX)
              this%h2D_Cytokinin1stConc_pvr(nptc,L)  = this%h2D_Cytokinin1stConc_pvr(nptc,L)+Cytokinin1stConc_rpvr(L,NR,NZ,NY,NX)
              this%h2D_Cytok_scalar_pvr(nptc,L)      = this%h2D_Cytok_scalar_pvr(nptc,L)+fctyok_scalar_rpvr(L,NR,NZ,NY,NX)
              this%h2D_Root1stStrutC_pvr(nptc,L)     = this%h2D_Root1stStrutC_pvr(nptc,L) + RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ,NY,NX)
              this%h2D_Root1stStrutN_pvr(nptc,L)     = this%h2D_Root1stStrutN_pvr(nptc,L) + RootMyco1stStrutElms_rpvr(ielmn,L,NR,NZ,NY,NX)
              this%h2D_Root1stStrutP_pvr(nptc,L)     = this%h2D_Root1stStrutP_pvr(nptc,L) + RootMyco1stStrutElms_rpvr(ielmp,L,NR,NZ,NY,NX)
            ENDDO

            if(NumStructuralRootAxes_pft(NZ,NY,NX).GT.0)THEN
              this%h2D_RootmedRadius_pvr(nptc,L)=this%h2D_RootmedRadius_pvr(nptc,L)*1.e3_r8/NumStructuralRootAxes_pft(NZ,NY,NX)
              this%h2D_Root1stLenPP_pvr(nptc,L)     = this%h2D_Root1stLenPP_pvr(nptc,L)/NumStructuralRootAxes_pft(NZ,NY,NX)
              this%h2D_Cytokinin1stConc_pvr(nptc,L) = AZERO(this%h2D_Cytokinin1stConc_pvr(nptc,L)/NumStructuralRootAxes_pft(NZ,NY,NX))
              this%h2D_Cytok_scalar_pvr(nptc,L)     = AZERO(this%h2D_Cytok_scalar_pvr(nptc,L)/NumStructuralRootAxes_pft(NZ,NY,NX))
              this%h2D_CRootLumenArea_pvr(nptc,L)   = AZERO(this%h2D_CRootLumenArea_pvr(nptc,L)/NumStructuralRootAxes_pft(NZ,NY,NX))
            ENDIF
            this%h1D_RootMeDStrutC_ptc(nptc) = this%h1D_RootMeDStrutC_ptc(nptc)+RootMedStruct_pvr(ielmc,L,NZ,NY,NX)
            this%h1D_Root1stStrutC_ptc(nptc) = this%h1D_Root1stStrutC_ptc(nptc)+this%h2D_Root1stStrutC_pvr(nptc,L)
            this%h1D_Root1stStrutN_ptc(nptc) = this%h1D_Root1stStrutN_ptc(nptc)+this%h2D_Root1stStrutN_pvr(nptc,L)
            this%h1D_Root2ndStrutC_ptc(nptc) = this%h1D_Root2ndStrutC_ptc(nptc)+this%h2D_Root2ndStrutC_pvr(nptc,L)
            this%h2D_MycoBiomC_pvr(nptc,L) = this%h2D_MycoBiomC_pvr(nptc,L)/DVOLL
            this%h2D_Root1stStrutC_pvr(nptc,L) = this%h2D_Root1stStrutC_pvr(nptc,L)/DVOLL
            this%h2D_Root1stStrutN_pvr(nptc,L) = this%h2D_Root1stStrutN_pvr(nptc,L)/DVOLL
            this%h2D_Root1stStrutP_pvr(nptc,L) = this%h2D_Root1stStrutP_pvr(nptc,L)/DVOLL
            this%h2D_Root2ndStrutC_pvr(nptc,L) = this%h2D_Root2ndStrutC_pvr(nptc,L)/DVOLL
            this%h2D_Root2ndStrutN_pvr(nptc,L) = this%h2D_Root2ndStrutN_pvr(nptc,L)/DVOLL
            this%h2D_Root2ndStrutP_pvr(nptc,L) = this%h2D_Root2ndStrutP_pvr(nptc,L)/DVOLL
          endif
        ENDDO
        this%h1D_RootAR_ptc(nptc)=this%h1D_RootAR_ptc(nptc)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_MycorrizhalBiomC_ptc(nptc)=this%h1D_MycorrizhalBiomC_ptc(nptc)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RootMeDStrutC_ptc(nptc) = this%h1D_RootMeDStrutC_ptc(nptc)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Root1stStrutC_ptc(nptc) = this%h1D_Root1stStrutC_ptc(nptc)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Root1stStrutN_ptc(nptc) = this%h1D_Root1stStrutN_ptc(nptc)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Root2ndStrutC_ptc(nptc) = this%h1D_Root2ndStrutC_ptc(nptc)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      ENDDO
      this%h1d_fPAR_col(ncol)=safe_adb(this%h1d_fPAR_col(ncol),RadPARSolarBeam_col(NY,NX))
  end procedure update_hist_plants

end submodule HistUpdatePlant
