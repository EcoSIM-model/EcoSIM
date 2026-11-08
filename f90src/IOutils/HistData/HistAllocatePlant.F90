submodule (HistDataType) HistAllocatePlant
  ! Allocate and initialize plant history-output buffers.
  use data_const_mod, only: spval => DAT_CONST_SPVAL, ispval => DAT_CONST_ISPVAL
  implicit none
contains

  module procedure allocate_hist_plants
    ! Initialize plant history fields to the missing-value sentinel.
    integer :: beg_col, end_col, beg_ptc, end_ptc

    beg_col=1; end_col=bounds%ncols
    beg_ptc=1; end_ptc=bounds%npfts

  allocate(this%h1D_CAN_G_ptc(beg_ptc:end_ptc))           ;this%h1D_CAN_G_ptc(:)=spval
  allocate(this%h1D_CAN_TEMPC_ptc(beg_ptc:end_ptc))        ;this%h1D_CAN_TEMPC_ptc(:)=spval
  allocate(this%h1D_CAN_TEMPFN_ptc(beg_ptc:end_ptc))      ;this%h1D_CAN_TEMPFN_ptc(:)=spval
  allocate(this%h1D_CAN_CO2_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_CAN_CO2_FLX_ptc(:)=spval
  allocate(this%h1D_CAN_GPP_ptc(beg_ptc:end_ptc))         ;this%h1D_CAN_GPP_ptc(:)=spval
  allocate(this%h1D_CAN_cumGPP_ptc(beg_ptc:end_ptc))         ;this%h1D_CAN_cumGPP_ptc(:)=spval
  allocate(this%h1D_dCAN_GPP_CLIM_ptc(beg_ptc:end_ptc))    ;this%h1D_dCAN_GPP_CLIM_ptc(:)=spval
  allocate(this%h1D_dCAN_GPP_eLIM_ptc(beg_ptc:end_ptc))    ;this%h1D_dCAN_GPP_eLIM_ptc(:)=spval
  allocate(this%h1D_CAN_RA_ptc(beg_ptc:end_ptc))          ;this%h1D_CAN_RA_ptc(:)=spval
  allocate(this%h1D_TurgEff4CanopyResp_ptc(beg_ptc:end_ptc)); this%h1D_TurgEff4CanopyResp_ptc(:)=spval
  allocate(this%h1D_LEAF_PC_ptc(beg_ptc:end_ptc))         ;this%h1D_LEAF_PC_ptc(:)=spval
  allocate(this%h1D_CAN_RN_ptc(beg_ptc:end_ptc))          ;this%h1D_CAN_RN_ptc(:)=spval
  allocate(this%h1D_CAN_LE_ptc(beg_ptc:end_ptc))          ;this%h1D_CAN_LE_ptc(:)=spval
  allocate(this%h1D_CAN_H_ptc(beg_ptc:end_ptc))           ;this%h1D_CAN_H_ptc(:)=spval
  allocate(this%h1D_CAN_GROWTH_ptc(beg_ptc:end_ptc))      ;this%h1D_CAN_GROWTH_ptc(:)=spval
  allocate(this%h1D_cTNC_ptc(beg_ptc:end_ptc))            ;this%h1D_cTNC_ptc(:)=spval
  allocate(this%h1D_cTNN_ptc(beg_ptc:end_ptc))            ;this%h1D_cTNN_ptc(:)=spval
  allocate(this%h1D_cTNP_ptc(beg_ptc:end_ptc))            ;this%h1D_cTNP_ptc(:)=spval
  allocate(this%h1D_CanopyNonstOMConc_ptc(beg_ptc:end_ptc)) ;this%h1D_CanopyNonstOMConc_ptc(:)=spval
  allocate(this%h1D_STOML_RSC_CO2_ptc(beg_ptc:end_ptc))   ;this%h1D_STOML_RSC_CO2_ptc(:)=spval
  allocate(this%h1D_STOML_Min_RSC_CO2_ptc(beg_ptc:end_ptc));this%h1D_STOML_Min_RSC_CO2_ptc(:)=spval
  allocate(this%h1D_Km_CO2_carboxy_ptc(beg_ptc:end_ptc)); this%h1D_Km_CO2_carboxy_ptc(:)=spval
  allocate(this%h1D_DynCi2CaRatio_ptc(beg_ptc:end_ptc));this%h1D_DynCi2CaRatio_ptc(:)=spval
  allocate(this%h1D_BLYR_RSC_CO2_ptc(beg_ptc:end_ptc))    ;this%h1D_BLYR_RSC_CO2_ptc(:)=spval
  allocate(this%h1D_CAN_CO2_ptc(beg_ptc:end_ptc))         ;this%h1D_CAN_CO2_ptc(:)=spval
  allocate(this%h1D_O2L_ptc(beg_ptc:end_ptc)); this%h1D_O2L_ptc(:)=spval
  allocate(this%h1D_LAI_ptc(beg_ptc:end_ptc))             ;this%h1D_LAI_ptc(:)=spval
  allocate(this%h1D_PSI_CAN_ptc(beg_ptc:end_ptc))         ;this%h1D_PSI_CAN_ptc(:)=spval
  allocate(this%h1D_CanPhenol_WSTRSS_ptc(beg_ptc:end_ptc));this%h1D_CanPhenol_WSTRSS_ptc(:)=spval
  allocate(this%h1D_CanPhenol_TSTRSS_ptc(beg_ptc:end_ptc)); this%h1D_CanPhenol_TSTRSS_ptc(:)=spval
  allocate(this%h1D_RootAR_ptc(beg_ptc:end_ptc))          ;this%h1D_RootAR_ptc(:)=spval
  allocate(this%h1D_RootLenPerPlant_ptc(beg_ptc:end_ptc))         ;this%h1D_RootLenPerPlant_ptc(:)=spval
  allocate(this%h1D_TURG_CAN_ptc(beg_ptc:end_ptc))        ;this%h1D_TURG_CAN_ptc(:)=spval
  allocate(this%h1D_STOML_RSC_H2O_ptc(beg_ptc:end_ptc))    ;this%h1D_STOML_RSC_H2O_ptc(:)=spval
  allocate(this%h1D_BLYR_RSC_H2O_ptc(beg_ptc:end_ptc))    ;this%h1D_BLYR_RSC_H2O_ptc(:)=spval
  allocate(this%h1D_BLYR_RSC_HEAT_ptc(beg_ptc:end_ptc))    ;this%h1D_BLYR_RSC_HEAT_ptc(:)=spval
  allocate(this%h1D_CdH2ORootxSoil_ptc(beg_ptc:end_ptc)) ;this%h1D_CdH2ORootxSoil_ptc(:)=spval
  allocate(this%h1D_TRANSPN_ptc(beg_ptc:end_ptc))         ;this%h1D_TRANSPN_ptc(:)=spval
  allocate(this%h1D_QTRANSP_col(beg_col:end_col)); this%h1D_QTRANSP_col(:)=spval
  allocate(this%h1D_fSnowCan_ptc(beg_ptc:end_ptc));this%h1D_fSnowCan_ptc(:)=spval
  allocate(this%h1D_PTSHTR_ptc(beg_ptc:end_ptc)); this%h1D_PTSHTR_ptc(:)=spval
  allocate(this%h1D_NH4_UPTK_FLX_ptc(beg_ptc:end_ptc))    ;this%h1D_NH4_UPTK_FLX_ptc(:)=spval
  allocate(this%h1D_NO3_UPTK_FLX_ptc(beg_ptc:end_ptc))    ;this%h1D_NO3_UPTK_FLX_ptc(:)=spval
  allocate(this%h1D_N2_FIXN_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_N2_FIXN_FLX_ptc(:)=spval
  allocate(this%h1D_cNH3_FLX_ptc(beg_ptc:end_ptc))        ;this%h1D_cNH3_FLX_ptc(:)=spval
  allocate(this%h1D_RNFixCO2_ptc(beg_ptc:end_ptc)); this%h1D_RNFixCO2_ptc(:)=spval
  allocate(this%h1D_TC_Groth_ptc(beg_ptc:end_ptc))        ;this%h1D_TC_Groth_ptc(:)=spval
  allocate(this%h1D_PO4_UPTK_FLX_ptc(beg_ptc:end_ptc))    ;this%h1D_PO4_UPTK_FLX_ptc(:)=spval
  allocate(this%h1D_SHOOT_C_ptc(beg_ptc:end_ptc))         ;this%h1D_SHOOT_C_ptc(:)=spval
  allocate(this%h1D_RNodeInitiate_ptc(beg_ptc:end_ptc));this%h1D_RNodeInitiate_ptc(:)=spval
  allocate(this%h1D_RLeafAppear_ptc(beg_ptc:end_ptc)); this%h1D_RLeafAppear_ptc(:)=spval
  allocate(this%h1D_Plant_C_ptc(beg_ptc:end_ptc))         ;this%h1D_Plant_C_ptc(:)=spval
  allocate(this%h1D_frcPARabs_ptc(beg_ptc:end_ptc))       ;this%h1D_frcPARabs_ptc(:)=spval
  allocate(this%h1D_PAR_CAN_ptc(beg_ptc:end_ptc))         ;this%h1D_PAR_CAN_ptc(:)=spval
  allocate(this%h1D_LEAF_C_ptc(beg_ptc:end_ptc))          ;this%h1D_LEAF_C_ptc(:)=spval
  allocate(this%h1D_Petole_C_ptc(beg_ptc:end_ptc))        ;this%h1D_Petole_C_ptc(:)=spval
  allocate(this%h1D_STALK_C_ptc(beg_ptc:end_ptc))         ;this%h1D_STALK_C_ptc(:)=spval
  allocate(this%h1D_RESERVE_C_ptc(beg_ptc:end_ptc))       ;this%h1D_RESERVE_C_ptc(:)=spval
  allocate(this%h1D_HUSK_C_ptc(beg_ptc:end_ptc))          ;this%h1D_HUSK_C_ptc(:)=spval
  allocate(this%h1D_GRAIN_C_ptc(beg_ptc:end_ptc))         ;this%h1D_GRAIN_C_ptc(:)=spval
  allocate(this%h1D_ROOT_C_ptc(beg_ptc:end_ptc))          ;this%h1D_ROOT_C_ptc(:)=spval
  allocate(this%h1D_CytokininShootExport_ptc(beg_ptc:end_ptc));this%h1D_CytokininShootExport_ptc(:)=spval
  allocate(this%h1D_CytokininRootLoss_ptc(beg_ptc:end_ptc));this%h1D_CytokininRootLoss_ptc(:)=spval
  allocate(this%h1D_RootActiveC_ptc(beg_ptc:end_ptc));this%h1D_RootActiveC_ptc(:)=spval
  allocate(this%h1D_ROOTST_C_ptc(beg_ptc:end_ptc))        ;this%h1D_ROOTST_C_ptc(:)=spval
  allocate(this%h1D_ROOTST_N_ptc(beg_ptc:end_ptc))        ;this%h1D_ROOTST_N_ptc(:)=spval
  allocate(this%h1D_ROOTST_P_ptc(beg_ptc:end_ptc))        ;this%h1D_ROOTST_P_ptc(:)=spval
  allocate(this%h1D_ShootNodule_C_ptc(beg_ptc:end_ptc))    ;this%h1D_ShootNodule_C_ptc(:)=spval
  allocate(this%h1D_ShootNodule_N_ptc(beg_ptc:end_ptc))    ;this%h1D_ShootNodule_N_ptc(:)=spval
  allocate(this%h1D_ShootNodule_P_ptc(beg_ptc:end_ptc))    ;this%h1D_ShootNodule_P_ptc(:)=spval
  allocate(this%h1D_RootNodule_C_ptc(beg_ptc:end_ptc))    ;this%h1D_RootNodule_C_ptc(:)=spval
  allocate(this%h1D_STORED_C_ptc(beg_ptc:end_ptc))        ;this%h1D_STORED_C_ptc(:)=spval
  allocate(this%h1D_GRAIN_NO_ptc(beg_ptc:end_ptc))        ;this%h1D_GRAIN_NO_ptc(:)=spval
  allocate(this%h1D_LAIb_ptc(beg_ptc:end_ptc))            ;this%h1D_LAIb_ptc(:)=spval
  allocate(this%h1D_EXUD_CumYr_C_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_EXUD_CumYr_C_FLX_ptc(:)=spval
  allocate(this%h1D_LITRf_C_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_LITRf_C_FLX_ptc=spval
  allocate(this%h1D_LITRf_P_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_LITRf_P_FLX_ptc=spval
  allocate(this%h1D_SURF_LITRf_C_FLX_ptc(beg_ptc:end_ptc));this%h1D_SURF_LITRf_C_FLX_ptc(:)=spval
  allocate(this%h1D_AUTO_RESP_FLX_ptc(beg_ptc:end_ptc))   ;this%h1D_AUTO_RESP_FLX_ptc(:)=spval
  allocate(this%h1D_HVST_C_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_HVST_C_FLX_ptc(:)=spval
  allocate(this%h1D_RootAct1stC_ptc(beg_ptc:end_ptc)) ; this%h1D_RootAct1stC_ptc(:)=spval
  allocate(this%h1D_MainBranchNodeNumber_ptc(beg_ptc:end_ptc));this%h1D_MainBranchNodeNumber_ptc(:)=spval
  allocate(this%h1D_ShootNodeNumber_ptc(beg_ptc:end_ptc)) ;this%h1D_ShootNodeNumber_ptc(:)=spval
  allocate(this%h1D_STANDING_DEAD_C_ptc(beg_ptc:end_ptc)) ;this%h1D_STANDING_DEAD_C_ptc(:)=spval
  allocate(this%h1D_FIREp_CO2_FLX_ptc(beg_ptc:end_ptc))   ;this%h1D_FIREp_CO2_FLX_ptc(:)=spval
  allocate(this%h1D_FIREp_CH4_FLX_ptc(beg_ptc:end_ptc))   ;this%h1D_FIREp_CH4_FLX_ptc(:)=spval
  allocate(this%h1D_NPP_ptc(beg_ptc:end_ptc))             ;this%h1D_NPP_ptc(:)=spval
  allocate(this%h1D_CAN_HT_ptc(beg_ptc:end_ptc))          ;this%h1D_CAN_HT_ptc(:)=spval
  allocate(this%h1D_Stalk_HT_ptc(beg_ptc:end_ptc))        ;this%h1D_Stalk_HT_ptc(:)=spval
  allocate(this%h1D_EmergeHeight_ptc(beg_ptc:end_ptc)); this%h1D_EmergeHeight_ptc(:)=spval
  allocate(this%h1D_POPN_ptc(beg_ptc:end_ptc))            ;this%h1D_POPN_ptc(:)=spval
  allocate(this%h1D_CanopyCutProxy_ptc(beg_ptc:end_ptc))  ;this%h1D_CanopyCutProxy_ptc(:)=spval
  allocate(this%h1D_tTRANSPN_ptc(beg_ptc:end_ptc))        ;this%h1D_tTRANSPN_ptc(:)=spval
  allocate(this%h1D_WTR_STRESS_ptc(beg_ptc:end_ptc))      ;this%h1D_WTR_STRESS_ptc(:)=spval
  allocate(this%h1D_OXY_STRESS_ptc(beg_ptc:end_ptc))      ;this%h1D_OXY_STRESS_ptc(:)=spval
  allocate(this%h1D_VcMaxRubisco_ptc(beg_ptc:end_ptc))    ;this%h1D_VcMaxRubisco_ptc(:)=spval
  allocate(this%h1D_LeafProteinNperm2_ptc(beg_ptc:end_ptc));this%h1D_LeafProteinNperm2_ptc(:)=spval
  allocate(this%h1D_VoMaxRubisco_ptc(beg_ptc:end_ptc))    ;this%h1D_VoMaxRubisco_ptc(:)=spval
  allocate(this%h1D_VcMaxPEP_ptc(beg_ptc:end_ptc))        ;this%h1D_VcMaxPEP_ptc(:)=spval
  allocate(this%h1D_JMaxPhoto_ptc(beg_ptc:end_ptc))       ;this%h1D_JMaxPhoto_ptc(:)=spval
  allocate(this%h1D_TFN_Carboxy_ptc(beg_ptc:end_ptc))     ;this%h1D_TFN_Carboxy_ptc(:)=spval
  allocate(this%h1D_TFN_Oxygen_ptc(beg_ptc:end_ptc))      ;this%h1D_TFN_Oxygen_ptc(:)=spval
  allocate(this%h1D_TFN_eTranspt_ptc(beg_ptc:end_ptc))    ;this%h1D_TFN_eTranspt_ptc(:)=spval

  allocate(this%h1D_PARSunlit_ptc(beg_ptc:end_ptc))  ;this%h1D_PARSunlit_ptc(:)=spval
  allocate(this%h1D_PARSunsha_ptc(beg_ptc:end_ptc))  ;this%h1D_PARSunsha_ptc(:)=spval
  allocate(this%h1D_CH2OSunlit_ptc(beg_ptc:end_ptc)) ;this%h1D_CH2OSunlit_ptc(:)=spval
  allocate(this%h1D_CH2OSunsha_ptc(beg_ptc:end_ptc)); this%h1D_CH2OSunsha_ptc(:)=spval
  allocate(this%h1D_fClump_ptc(beg_ptc:end_ptc)); this%h1D_fClump_ptc(:)=spval
  allocate(this%h1D_LeafAreaSunlit_ptc(beg_ptc:end_ptc))  ;this%h1D_LeafAreaSunlit_ptc(:)=spval
  allocate(this%h1D_SHOOT_N_ptc(beg_ptc:end_ptc))         ;this%h1D_SHOOT_N_ptc(:)=spval
  allocate(this%h1D_Plant_N_ptc(beg_ptc:end_ptc))         ;this%h1D_Plant_N_ptc(:)=spval
  allocate(this%h1D_LEAF_N_ptc(beg_ptc:end_ptc))          ;this%h1D_LEAF_N_ptc(:)=spval
  allocate(this%h1D_fCNLFW_ptc(beg_ptc:end_ptc)) ; this%h1D_fCNLFW_ptc(:)=spval
  allocate(this%h1D_fCPLFW_ptc(beg_ptc:end_ptc)); this%h1D_fCPLFW_ptc(:)=spval
  allocate(this%h1D_LeafNperm2LAI_ptc(beg_ptc:end_ptc))      ;this%h1D_LeafNperm2LAI_ptc(:)=spval
  allocate(this%h1D_Petole_N_ptc(beg_ptc:end_ptc))       ;this%h1D_Petole_N_ptc(:)=spval
  allocate(this%h1D_STALK_N_ptc(beg_ptc:end_ptc))         ;this%h1D_STALK_N_ptc(:)=spval
  allocate(this%h1D_RESERVE_N_ptc(beg_ptc:end_ptc))       ;this%h1D_RESERVE_N_ptc(:)=spval
  allocate(this%h1D_HUSK_N_ptc(beg_ptc:end_ptc))          ;this%h1D_HUSK_N_ptc(:)=spval
  allocate(this%h1D_GRAIN_N_ptc(beg_ptc:end_ptc))         ;this%h1D_GRAIN_N_ptc(:)=spval
  allocate(this%h1D_ROOT_N_ptc(beg_ptc:end_ptc))          ;this%h1D_ROOT_N_ptc(:)=spval

  allocate(this%h1D_RootNodule_N_ptc(beg_ptc:end_ptc))        ;this%h1D_RootNodule_N_ptc(:)=spval
  allocate(this%h1D_STORED_N_ptc(beg_ptc:end_ptc))        ;this%h1D_STORED_N_ptc(:)=spval
  allocate(this%h1D_EXUD_N_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_EXUD_N_FLX_ptc(:)=spval
  allocate(this%h1D_Uptk_NMIN_Flx_ptc(beg_ptc:end_ptc))   ;this%h1D_Uptk_NMIN_Flx_ptc(:)=spval
  allocate(this%h1D_Uptk_PMIN_Flx_ptc(beg_ptc:end_ptc))   ;this%h1D_Uptk_PMIN_Flx_ptc(:)=spval
  allocate(this%h1D_Uptk_N_Flx_ptc(beg_ptc:end_ptc))      ;this%h1D_Uptk_N_Flx_ptc(:)=spval
  allocate(this%h1D_Uptk_P_Flx_ptc(beg_ptc:end_ptc))      ;this%h1D_Uptk_P_Flx_ptc(:)=spval
  allocate(this%h1D_LITRf_N_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_LITRf_N_FLX_ptc(:)=spval
  allocate(this%h1D_cum_N_FIXED_ptc(beg_ptc:end_ptc))  ;this%h1D_cum_N_FIXED_ptc(:)=spval
  allocate(this%h1D_TreeRingRadius_ptc(beg_ptc:end_ptc));this%h1D_TreeRingRadius_ptc(:)=spval
  allocate(this%h1D_HVST_N_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_HVST_N_FLX_ptc(:)=spval
  allocate(this%h1D_NH3can_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_NH3can_FLX_ptc(:)=spval
  allocate(this%h1D_PLANT_BALANCE_C_ptc(beg_ptc:end_ptc)) ;this%h1D_PLANT_BALANCE_C_ptc(:)=spval
  allocate(this%h1D_PLANT_BALANCE_N_ptc(beg_ptc:end_ptc)) ;this%h1D_PLANT_BALANCE_N_ptc(:)=spval
  allocate(this%h1D_PLANT_BALANCE_P_ptc(beg_ptc:end_ptc)) ;this%h1D_PLANT_BALANCE_P_ptc(:)=spval
  allocate(this%h1D_STANDING_DEAD_N_ptc(beg_ptc:end_ptc)) ;this%h1D_STANDING_DEAD_N_ptc(:)=spval
  allocate(this%h1D_FIREp_N_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_FIREp_N_FLX_ptc(:)=spval
  allocate(this%h1D_SURF_LITRf_N_FLX_ptc(beg_ptc:end_ptc));this%h1D_SURF_LITRf_N_FLX_ptc(:)=spval
  allocate(this%h1D_SHOOT_P_ptc(beg_ptc:end_ptc))         ;this%h1D_SHOOT_P_ptc(:)=spval
  allocate(this%h1D_Plant_P_ptc(beg_ptc:end_ptc))         ;this%h1D_Plant_P_ptc(:)=spval
  allocate(this%h1D_stomatal_stress_ptc(beg_ptc:end_ptc)) ;this%h1D_stomatal_stress_ptc(:)=spval
  allocate(this%h1D_CANDew_ptc(beg_ptc:end_ptc)); this%h1D_CANDew_ptc(:)=spval
  allocate(this%h1D_LEAF_P_ptc(beg_ptc:end_ptc))          ;this%h1D_LEAF_P_ptc(:)=spval
  allocate(this%h1D_Petole_P_ptc(beg_ptc:end_ptc))        ;this%h1D_Petole_P_ptc(:)=spval
  allocate(this%h1D_STALK_P_ptc(beg_ptc:end_ptc))         ;this%h1D_STALK_P_ptc(:)=spval
  allocate(this%h1D_RESERVE_P_ptc(beg_ptc:end_ptc))       ;this%h1D_RESERVE_P_ptc(:)=spval
  allocate(this%h1D_HUSK_P_ptc(beg_ptc:end_ptc))          ;this%h1D_HUSK_P_ptc(:)=spval
  allocate(this%h1D_GRAIN_P_ptc(beg_ptc:end_ptc))         ;this%h1D_GRAIN_P_ptc(:)=spval
  allocate(this%h1D_ROOT_P_ptc(beg_ptc:end_ptc))          ;this%h1D_ROOT_P_ptc(:)=spval
  allocate(this%h1D_RootNodule_P_ptc(beg_ptc:end_ptc))        ;this%h1D_RootNodule_P_ptc(:)=spval
  allocate(this%h1D_STORED_P_ptc(beg_ptc:end_ptc))        ;this%h1D_STORED_P_ptc(:)=spval
  allocate(this%h1D_EXUD_P_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_EXUD_P_FLX_ptc(:)=spval
  allocate(this%h1D_LITTERf_P_ptc(beg_ptc:end_ptc))       ;this%h1D_LITTERf_P_ptc(:)=spval
  allocate(this%h1D_HVST_P_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_HVST_P_FLX_ptc(:)=spval
  allocate(this%h1D_STANDING_DEAD_P_ptc(beg_ptc:end_ptc)) ;this%h1D_STANDING_DEAD_P_ptc(:)=spval
  allocate(this%h1D_FIREp_P_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_FIREp_P_FLX_ptc(:)=spval
  allocate(this%h1D_SURF_LITRf_P_FLX_ptc(beg_ptc:end_ptc));this%h1D_SURF_LITRf_P_FLX_ptc(:)=spval
  allocate(this%h1D_ShootRootXferC_ptc(beg_ptc:end_ptc)); this%h1D_ShootRootXferC_ptc(:)=spval
  allocate(this%h1D_ShootRootXferN_ptc(beg_ptc:end_ptc)); this%h1D_ShootRootXferN_ptc(:)=spval
  allocate(this%h1D_ShootRootXferP_ptc(beg_ptc:end_ptc)); this%h1D_ShootRootXferP_ptc(:)=spval
  allocate(this%h1D_BRANCH_NO_ptc(beg_ptc:end_ptc))       ;this%h1D_BRANCH_NO_ptc(:)=spval
  allocate(this%h1D_Growth_Stage_ptc(beg_ptc:end_ptc))    ;this%h1D_Growth_Stage_ptc(:)=spval
  allocate(this%h1D_LEAF_NC_ptc(beg_ptc:end_ptc));        ;this%h1D_LEAF_NC_ptc(:)=spval
  allocate(this%h1D_SHOOT_NONSTC_ptc(beg_ptc:end_ptc))     ;this%h1D_SHOOT_NONSTC_ptc(:)=spval
  allocate(this%h1D_SHOOT_NONSTN_ptc(beg_ptc:end_ptc))     ;this%h1D_SHOOT_NONSTN_ptc(:)=spval
  allocate(this%h1D_SHOOT_NONSTP_ptc(beg_ptc:end_ptc))     ;this%h1D_SHOOT_NONSTP_ptc(:)=spval

  allocate(this%h1D_LeafChlCperm2LA_ptc(beg_ptc:end_ptc));this%h1D_LeafChlCperm2LA_ptc=spval
  allocate(this%h1D_LeafC4ChlCperm2LA_ptc(beg_ptc:end_ptc));this%h1D_LeafC4ChlCperm2LA_ptc=spval
  allocate(this%h1D_LeafRubiscoNperm2LA_ptc(beg_ptc:end_ptc));this%h1D_LeafRubiscoNperm2LA_ptc=spval
  allocate(this%h1D_LeafPEPCNperm2LA_ptc(beg_ptc:end_ptc));this%h1D_LeafPEPCNperm2LA_ptc=spval
!  allocate(this%h1D_CFIX_lmtf_ptc(beg_ptc:end_ptc))        ;this%h1D_CFIX_lmtf_ptc(:)=spval
  allocate(this%h1D_MainBranchNO_ptc(beg_ptc:end_ptc))     ;this%h1D_MainBranchNO_ptc(:)=ispval
  allocate(this%h1D_RCanMaintDef_CO2_pft(beg_ptc:end_ptc)) ;this%h1D_RCanMaintDef_CO2_pft(:)=spval
  allocate(this%h1D_RootMaintDef_CO2_pft(beg_ptc:end_ptc)) ;this%h1D_RootMaintDef_CO2_pft(:)=spval
  allocate(this%h1D_NumPrimeRootAxes_ptc(beg_ptc:end_ptc)) ;this%h1D_NumPrimeRootAxes_ptc(:)=spval
  allocate(this%h1D_ROOT_NONSTC_ptc(beg_ptc:end_ptc))     ;this%h1D_ROOT_NONSTC_ptc(:)=spval
  allocate(this%h1D_MycorrizhalBiomC_ptc(beg_ptc:end_ptc));        this%h1D_MycorrizhalBiomC_ptc(:) = spval
  allocate(this%h1D_Root1stStrutC_ptc(beg_ptc:end_ptc)); this%h1D_Root1stStrutC_ptc(:)=spval
  allocate(this%h1D_RootMeDStrutC_ptc(beg_ptc:end_ptc)); this%h1D_RootMeDStrutC_ptc(:)=spval
  allocate(this%h1D_Root1stStrutN_ptc(beg_ptc:end_ptc)); this%h1D_Root1stStrutN_ptc(:)=spval
  allocate(this%h1D_Root2ndStrutC_ptc(beg_ptc:end_ptc)); this%h1D_Root2ndStrutC_ptc(:)=spval

  allocate(this%h1D_ROOT_NONSTN_ptc(beg_ptc:end_ptc))     ;this%h1D_ROOT_NONSTN_ptc(:)=spval
  allocate(this%h1D_ROOT_NONSTP_ptc(beg_ptc:end_ptc))     ;this%h1D_ROOT_NONSTP_ptc(:)=spval
  end procedure allocate_hist_plants

end submodule HistAllocatePlant
