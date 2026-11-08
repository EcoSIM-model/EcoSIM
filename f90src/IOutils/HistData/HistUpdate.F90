submodule (HistDataType) HistUpdate
  use GridConsts, only: JZ, NumOfPlantMorphUnits
  use GridMod, only: get_col
  use DebugToolMod, only: PrintInfo
  implicit none
contains

  module procedure hist_update
    integer :: NX,NY,ncol
    character(len=*), parameter :: subname='hist_update'

    call PrintInfo('beg '//subname)
    DO NX=bounds%NHW,bounds%NHE
      DO NY=bounds%NVN,bounds%NVS
        ncol=get_col(NY,NX)
        ! Keep column initialization, soil accumulation, and plant accumulation
        ! in their original order, including all resets and final normalization.
        call update_hist_columns(this,I,J,NY,NX,ncol)
        call update_hist_soil(this,I,J,NY,NX,ncol)
        call update_hist_plants(this,I,J,NY,NX,ncol)
      ENDDO
    ENDDO
    call PrintInfo('end '//subname)
  end procedure hist_update

  module procedure ZeroPlantHistVars
  implicit none
  this%h1D_PTSHTR_ptc(nptc) = 0._r8
  this%h1D_fSnowCan_ptc(nptc) = 0._r8
  this%h1D_RootAct1stC_ptc(nptc) = 0._r8
  this%h1D_RootMeDStrutC_ptc(nptc)=0._r8
  this%h1D_Root2ndStrutC_ptc(nptc)          = 0._r8
  this%h1D_Root1stStrutC_ptc(nptc)          = 0._r8
  this%h1D_Root1stStrutN_ptc(nptc)          = 0._r8
  this%h1D_MycorrizhalBiomC_ptc(nptc)        = 0._r8
  this%h1D_ShootRootXferC_ptc(nptc)          = 0._r8
  this%h1D_ShootRootXferN_ptc(nptc)          = 0._r8
  this%h1D_ShootRootXferP_ptc(nptc)          = 0._r8
  this%h1D_RootLenPerPlant_ptc(nptc)         = 0._r8
  this%h2d_RootPop_pvr(nptc,1:JZ)            = 0._r8
  this%h2D_MycoPop_pvr(nptc,1:JZ)            = 0._r8
  this%h2D_RootRadialKond2H2O_pvr(nptc,1:JZ) = 0._r8
  this%h2D_RootAXialKond2H2O_pvr(nptc,1:JZ)  = 0._r8
  this%h2D_VmaxNH4Root_pvr(nptc,1:JZ)        = 0._r8
  this%h2D_VmaxNO3Root_pvr(nptc,1:JZ)        = 0._r8
  this%h2D_RootMassC_pvr(nptc,1:JZ)          = 0._r8
  this%h2D_RootNutupk_fClim_pvr(nptc,1:JZ)   = 0._r8
  this%h2D_RootNutupk_fNlim_pvr(nptc,1:JZ)   = 0._r8
  this%h2D_RootNutupk_fPlim_pvr(nptc,1:JZ)   = 0._r8
  this%h2D_RootNutupk_fProtC_pvr(nptc,1:JZ)  = 0._r8
  this%h2D_Root1stSArea4GasTP_pvr(nptc,1:JZ)  = 0._r8
  this%h2D_RootProteinC_pvr(nptc,1:JZ)       = 0._r8
  this%h2D_O2_rootconduct_pvr(nptc,1:JZ)     = 0._r8
  this%h2D_CO2_rootconduct_pvr(nptc,1:JZ)    = 0._r8
  this%h2D_fTRootGro_pvr(nptc,1:JZ)          = 0._r8
  this%h2D_fRootGrowPSISense_pvr(nptc,1:JZ)  = 0._r8
  this%h2D_RootAbsorbAreaPP_pvr(nptc,1:JZ)     = 0._r8
  this%h1D_RootAbsorbAreaPP_pft(nptc)=0._r8
  this%h2D_ROOT_OSTRESS_pvr(nptc,1:JZ)       = 0._r8
  this%h2D_PSI_RT_pvr(nptc,1:JZ)             = 0._r8
  this%h2D_RootH2OUptkStress_pvr(nptc,1:JZ)  = 0._r8
  this%h2D_RootH2OUptk_pvr(nptc,1:JZ)        = 0._r8
  this%h2D_SapFlowVlinear_pvr(nptc,1:JZ)     = 0._R8
  this%h2D_RootMaintDef_CO2_pvr(nptc,1:JZ)   = 0._r8
  this%h2D_prtUP_NH4_pvr(nptc,1:JZ)          = 0._r8
  this%h2D_prtUP_NO3_pvr(nptc,1:JZ)          = 0._r8
  this%h2D_prtUP_PO4_pvr(nptc,1:JZ)          = 0._r8
  this%h2D_DNS_RT_pvr(nptc,1:JZ)             = 0._r8

  this%h2D_ROOTNLim_rpvr(nptc,1:JZ)       = 0._r8
  this%h2D_ROOTPLim_rpvr(nptc,1:JZ)       = 0._r8
  this%h2D_RootNonstC_rpvr(nptc,1:JZ)     = 0._r8
  this%h2D_RootSinkWeight_pvr(nptc,1:JZ)  = 0._r8
  this%h2D_RootMSinkWeight_pvr(nptc,1:JZ) = 0._r8
  this%h2D_Root1stSinkWeight_pvr(nptc,1:JZ) = 0._r8
  this%h2D_Root2ndSinkWeight_pvr(nptc,1:JZ) = 0._r8
  this%h2D_Root1stRadius_rpvr(nptc,1:JZ)  = 0._r8
  this%h2D_RootNonstBConc_pvr(nptc,1:JZ)  = 0._r8
  this%h2D_Root1stAxesNumL_pvr(nptc,1:JZ) = 0._r8
  this%h2D_RootMedAxesNumL_pvr(nptc,1:JZ) = 0._r8
  this%h2D_RootmedRadius_pvr(nptc,1:JZ) = 0._r8
  this%h2D_Rootmedlength_pvr(nptc,1:JZ) = 0._r8
  this%h2D_Root1stLenPP_pvr(nptc,1:JZ)    = 0._r8
  this%h2D_Root2ndAxesNumL_pvr(nptc,1:JZ) = 0._r8
  this%h2D_Root2ndMeanLen_pvr(nptc,1:JZ)  = 0._r8
  this%h2D_RootKond2H2O_pvr(nptc,1:JZ)    = 0._r8

  this%h2D_RootMedC_pvr(nptc,1:JZ) = 0._r8
  this%h2D_RootAct1stC_pvr(nptc,1:JZ)    = 0._r8
  this%h2D_RootLig1stC_pvr(nptc,1:JZ)    = 0._r8
  this%h2D_NonstC_conc_pvr(nptc,1:JZ)    = 0._r8
  this%h1D_ROOT_NONSTC_ptc(nptc)  = 0._r8
  this%h1D_ROOT_NONSTN_ptc(nptc)  = 0._r8
  this%h1D_ROOT_NONSTP_ptc(nptc)  = 0._r8
  this%h1D_SHOOT_NONSTC_ptc(nptc) = 0._r8
  this%h1D_SHOOT_NONSTN_ptc(nptc) = 0._r8
  this%h1D_SHOOT_NONSTP_ptc(nptc) = 0._r8

  this%h1D_dCAN_GPP_CLIM_ptc(nptc) = 0._r8
  this%h1D_dCAN_GPP_eLIM_ptc(nptc) = 0._r8

  this%h1D_LeafChlCperm2LA_ptc(nptc)   = 0._r8
  this%h1D_LeafC4ChlCperm2LA_ptc(nptc)   = 0._r8
  this%h1D_LeafRubiscoNperm2LA_ptc(nptc) = 0._r8
  this%h1D_LeafPEPCNperm2LA_ptc(nptc)     = 0._r8

  this%h1D_MIN_LWP_ptc(nptc)      = 0._r8
  this%h1D_SLA_ptc(nptc)          = 0._r8
  this%h1D_LEAF_PC_ptc(nptc)       = 0._r8
  this%h1D_CAN_RN_ptc(nptc)        = 0._r8
  this%h1D_CAN_LE_ptc(nptc)        = 0._r8
  this%h1D_CAN_H_ptc(nptc)         = 0._r8
  this%h1D_CAN_G_ptc(nptc)         = 0._r8
  this%h1D_CAN_TEMPC_ptc(nptc)     = 0._r8
  this%h1D_CAN_TEMPFN_ptc(nptc)    =0._r8
  this%h1D_CAN_CO2_FLX_ptc(nptc)   = 0._r8
  this%h1D_CAN_GPP_ptc(nptc)       = 0._r8
  this%h1D_CAN_RA_ptc(nptc)        = 0._r8
  this%h1D_TurgEff4CanopyResp_ptc(nptc) = 0._r8
  this%h1D_CAN_GROWTH_ptc(nptc)    = 0._r8
  this%h1D_cTNC_ptc(nptc)          = 0._r8
  this%h1D_cTNN_ptc(nptc)          = 0._r8
  this%h1D_cTNP_ptc(nptc)          = 0._r8
  this%h1D_CanopyNonstOMConc_ptc(nptc) = 0._r8
  this%h1D_STOML_RSC_CO2_ptc(nptc) = 0._r8
  this%h1D_STOML_Min_RSC_CO2_ptc(nptc)=0._r8
  this%h1D_Km_CO2_carboxy_ptc(nptc)= 0._r8
  this%h1D_DynCi2CaRatio_ptc(nptc) = 0._r8
  this%h1D_BLYR_RSC_CO2_ptc(nptc)  =0._r8
  this%h1D_CAN_CO2_ptc(nptc)       = 0._r8
  this%h1D_O2L_ptc(nptc)           = 0._r8
  this%h1D_LAI_ptc(nptc)           = 0._r8
  this%h1D_PSI_CAN_ptc(nptc)       = 0._r8
  this%h1D_TURG_CAN_ptc(nptc)      = 0._r8
  this%h1D_CanPhenol_WSTRSS_ptc(nptc)=0._r8
  this%h1D_CanPhenol_TSTRSS_ptc(nptc)=0._r8
  this%h1D_STOML_RSC_H2O_ptc(nptc)  = 0._r8
  this%h1D_BLYR_RSC_H2O_ptc(nptc)  = 0._r8
  this%h1D_BLYR_RSC_HEAT_ptc(nptc)  = 0._r8
  this%h1D_CdH2ORootxSoil_ptc(nptc) = 0._r8
  this%h1D_TRANSPN_ptc(nptc)       = 0._r8
  this%h1D_NH4_UPTK_FLX_ptc(nptc)  = 0._r8
  this%h1D_NO3_UPTK_FLX_ptc(nptc)  = 0._r8
  this%h1D_N2_FIXN_FLX_ptc(nptc)   = 0._r8
  this%h1D_cNH3_FLX_ptc(nptc)      = 0._r8
  this%h1D_RNFixCO2_ptc(nptc)      = 0._R8
  this%h1D_TC_Groth_ptc(nptc)      = 0._r8
  this%h1D_PO4_UPTK_FLX_ptc(nptc)  = 0._r8
  this%h1D_frcPARabs_ptc(nptc)     = 0._r8
  this%h1D_PAR_CAN_ptc(nptc)       = 0._r8
  this%h1D_SHOOT_C_ptc(nptc)       = 0._r8
  this%h1D_Plant_C_ptc(nptc)       = 0._r8
  this%h1D_RNodeInitiate_ptc(nptc) = 0._r8
  this%h1D_RLeafAppear_ptc(nptc) = 0._r8
  this%h1D_LEAF_C_ptc(nptc)        = 0._r8
  this%h1D_Petole_C_ptc(nptc)      = 0._r8
  this%h1D_STALK_C_ptc(nptc)       =  0._r8
  this%h1D_RESERVE_C_ptc(nptc)     = 0._r8
  this%h1D_HUSK_C_ptc(nptc)        = 0._r8
  this%h1D_GRAIN_C_ptc(nptc)       = 0._r8
  this%h1D_ROOT_C_ptc(nptc)        = 0._r8
  this%h1D_ROOTST_C_ptc(nptc)      = 0._r8
  this%h1D_ROOTST_N_ptc(nptc)      = 0._r8
  this%h1D_ROOTST_P_ptc(nptc)      =0._r8
  this%h1D_RootNodule_C_ptc(nptc)  =0._r8
  this%h1D_ShootNodule_C_ptc(nptc)  = 0._r8
  this%h1D_ShootNodule_N_ptc(nptc)  = 0._r8
  this%h1D_ShootNodule_P_ptc(nptc)  = 0._r8

  this%h1D_STORED_C_ptc(nptc)      = 0._r8
  this%h1D_GRAIN_NO_ptc(nptc)      = 0._r8
  this%h1D_LAIb_ptc(nptc)          = 0._r8
  this%h1D_LITRf_C_FLX_ptc(nptc)      =0._r8
  this%h1D_SURF_LITRf_C_FLX_ptc(nptc) = 0._r8
  this%h1D_AUTO_RESP_FLX_ptc(nptc)    = 0._r8
  this%h1D_HVST_C_FLX_ptc(nptc)       = 0._r8
  this%h1D_HVST_N_FLX_ptc(nptc)       = 0._r8
  this%h1D_HVST_P_FLX_ptc(nptc)       = 0._r8
  this%h1D_CAN_HT_ptc(nptc)           = 0._r8
  this%h1D_Stalk_HT_ptc(nptc)          = 0._r8
  this%h1D_EmergeHeight_ptc(nptc)      = 0._r8
  this%h1D_WTR_STRESS_ptc(nptc)       = 0._r8
  this%h1D_LeafProteinNperm2_ptc(nptc)= 0._r8
  this%h1D_VcMaxRubisco_ptc(nptc) =  0._r8
  this%h1D_VoMaxRubisco_ptc(nptc) =  0._r8
  this%h1D_VcMaxPEP_ptc(nptc)     =  0._r8
  this%h1D_JMaxPhoto_ptc(nptc)    =  0._r8
  this%h1D_TFN_Carboxy_ptc(nptc)  =  0._r8
  this%h1D_TFN_Oxygen_ptc(nptc)   =  0._r8
  this%h1D_TFN_eTranspt_ptc(nptc) =  0._r8
  this%h1D_PARSunlit_ptc(nptc)    =  0._r8
  this%h1D_PARSunsha_ptc(nptc)    =  0._r8
  this%h1D_CH2OSunlit_ptc(nptc)   =  0._r8
  this%h1D_CH2OSunsha_ptc(nptc)   =  0._r8

  this%h1D_fClump_ptc(nptc)  =  0._r8
  this%h1D_LeafAreaSunlit_ptc(nptc)= 0._r8
  this%h1D_OXY_STRESS_ptc(nptc)   =  0._r8
  this%h1D_SHOOT_N_ptc(nptc)      =  0._r8
  this%h1D_Plant_N_ptc(nptc)      =  0._r8
  this%h1D_fCNLFW_ptc(nptc) =  0._r8
  this%h1D_fCPLFW_ptc(nptc) =  0._r8
  this%h1D_LEAF_N_ptc(nptc)    =  0._r8
  this%h1D_LeafNperm2LAI_ptc(nptc) =  0._r8
  this%h1D_Petole_N_ptc(nptc)  =  0._r8
  this%h1D_STALK_N_ptc(nptc)   =  0._r8
  this%h1D_RESERVE_N_ptc(nptc) =  0._r8
  this%h1D_HUSK_N_ptc(nptc)    =  0._r8
  this%h1D_GRAIN_N_ptc(nptc)          =  0._r8
  this%h1D_ROOT_N_ptc(nptc)           =  0._r8
  this%h1D_RootNodule_N_ptc(nptc)     =  0._r8
  this%h1D_STORED_N_ptc(nptc)         =  0._r8
  this%h1D_TreeRingRadius_ptc(nptc)   =  0._r8
  this%h1D_SHOOT_P_ptc(nptc)          =  0._r8
  this%h1D_Root1stTipSinkWt_ptc(nptc) =  0._r8
  this%h1D_Plant_P_ptc(nptc)          =  0._r8
  this%h1D_stomatal_stress_ptc(nptc) =  0._r8
  this%h1D_LEAF_P_ptc(nptc)          =  0._r8
  this%h1D_Petole_P_ptc(nptc)        =  0._r8
  this%h1D_STALK_P_ptc(nptc)         = 0._r8
  this%h1D_RESERVE_P_ptc(nptc)       =  0._r8
  this%h1D_HUSK_P_ptc(nptc)          =  0._r8
  this%h1D_GRAIN_P_ptc(nptc)          =  0._r8
  this%h1D_ROOT_P_ptc(nptc)           =  0._r8
  this%h1D_RootNodule_P_ptc(nptc)     =  0._r8
  this%h1D_STORED_P_ptc(nptc)         =  0._r8
  this%h1D_BRANCH_NO_ptc(nptc)        = 0._r8
  this%h1D_MainBranchNO_ptc(nptc)     = 0._r8
  this%h1D_RCanMaintDef_CO2_pft(nptc) = 0._r8
  this%h1D_LEAF_NC_ptc(nptc)      = 0._r8
  this%h1D_RootMaintDef_CO2_pft(nptc) = 0._r8
  this%h1D_NumPrimeRootAxes_ptc(nptc) = 0._r8
  this%h2D_ProteinNperm2LeafArea_pnd(nptc,:)=0._r8
  this%h1D_Growth_Stage_ptc(nptc)   =0

  this%h1D_Num_Leaves_ptc(nptc)                     = 0._r8
  this%h1D_RUB_ACTVN_ptc(nptc)                      = 0._r8;
  this%h1D_CanopyNLim_ptc(nptc)                     = 0._r8
  this%h1D_CanopyPLim_ptc(nptc)                     = 0._r8
  this%h3D_PARTS_ptc(nptc,1:NumOfPlantMorphUnits,:) = 0._r8
  this%h2D_RootShootExchC_pvr(nptc,:)               = 0._r8
  this%h2D_RootShootExchN_pvr(nptc,:)               = 0._r8
  this%h2D_RootShootExchP_pvr(nptc,:)               = 0._r8
  this%h2D_CanopyLAIZ_plyr(nptc,:)                  = 0._r8
  this%h1D_RootAR_ptc(nptc)                         = 0._r8
  this%h1D_RootLenPerPlant_ptc(nptc)                = 0._r8
  this%h2D_Root1stDepz_ptc(nptc,:)                  = 0._r8
  this%h2D_Root1stStrutC_pvr(nptc,1:JZ)             = 0._r8
  this%h2D_MycoBiomC_pvr(nptc,1:JZ)                 = 0._R8
  this%h2D_Root1stStrutN_pvr(nptc,1:JZ)             = 0._r8
  this%h2D_Root1stStrutP_pvr(nptc,1:JZ)             = 0._r8
  this%h2D_Root2ndStrutC_pvr(nptc,1:JZ)             = 0._r8
  this%h2D_Root2ndStrutN_pvr(nptc,1:JZ)             = 0._r8
  this%h2D_Root2ndStrutP_pvr(nptc,1:JZ)             = 0._r8
  this%h2D_Cytokinin1stConc_pvr(nptc,1:JZ)         = 0._R8
  this%h2D_CRootLumenArea_pvr(nptc,1:JZ)           = 0._r8
  this%h2D_Cytok_scalar_pvr(nptc,1:JZ)              = 0._R8
  end procedure ZeroPlantHistVars
end submodule HistUpdate
