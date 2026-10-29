submodule (HistDataType) HistAllocateProfile
  ! Allocate soil- and PFT-layer history storage and initialize missing values.
  use data_const_mod, only: spval => DAT_CONST_SPVAL
  use GridConsts, only: JZ, MaxNodesPerBranch, MaxNumBranches, MaxNumRootAxes, NumCanopyLayers, NumOfPlantMorphUnits
  use ElmIDMod, only: NumPlantChemElms
  use EcoSIMConfig, only: jcplx => jcplxc
  implicit none
contains

  module procedure allocate_hist_profiles
    ! Allocate the layer diagnostics registered with the history system.
    integer :: beg_col, end_col, beg_ptc, end_ptc

    beg_col=1; end_col=bounds%ncols
    beg_ptc=1; end_ptc=bounds%npfts

  allocate(this%h2D_CO2_MassL_vr(beg_col:end_col,1:JZ)); this%h2D_CO2_MassL_vr(:,:)=spval
  allocate(this%h2D_O2_MassL_vr(beg_col:end_col,1:JZ)); this%h2D_O2_MassL_vr(:,:)=spval
  allocate(this%h2D_CH4_MassL_vr(beg_col:end_col,1:JZ)); this%h2D_CH4_MassL_vr(:,:)=spval
  allocate(this%h2D_H2_MassL_vr(beg_col:end_col,1:JZ)); this%h2D_H2_MassL_vr(:,:)=spval
  allocate(this%h2D_tSOC_vr(beg_col:end_col,1:JZ))    ;this%h2D_tSOC_vr(:,:)=spval
  allocate(this%h2D_POM_C_vr(beg_col:end_col,1:JZ)) ; this%h2D_POM_C_vr(:,:)=spval
  allocate(this%h2D_MAOM_C_vr(beg_col:end_col,1:JZ)); this%h2D_MAOM_C_vr(:,:)=spval
  allocate(this%h2D_microbC_vr(beg_col:end_col,1:JZ))  ;this%h2D_microbC_vr(:,:)=spval
  allocate(this%h2D_microbN_vr(beg_col:end_col,1:JZ))  ;this%h2D_microbN_vr(:,:)=spval
  allocate(this%h2D_microbP_vr(beg_col:end_col,1:JZ))  ;this%h2D_microbP_vr(:,:)=spval
  allocate(this%h2D_AeroBact_PrimS_lim_vr(beg_col:end_col,1:JZ)); this%h2D_AeroBact_PrimS_lim_vr(:,:)=spval
  allocate(this%h2D_AeroFung_PrimS_lim_vr(beg_col:end_col,1:JZ)); this%h2D_AeroFung_PrimS_lim_vr(:,:)=spval
  allocate(this%h2D_tSOCL_vr(beg_col:end_col,1:JZ))    ;this%h2D_tSOCL_vr(:,:)=spval
  allocate(this%h2D_NO3_vr(beg_col:end_col,1:JZ)); this%h2D_NO3_vr(:,:)=spval
  allocate(this%h2D_NH4_vr(beg_col:end_col,1:JZ)); this%h2D_NH4_vr(:,:)=spval
  allocate(this%h2D_tSON_vr(beg_col:end_col,1:JZ))    ;this%h2D_tSON_vr(:,:)=spval
  allocate(this%h2D_tSOP_vr(beg_col:end_col,1:JZ))    ;this%h2D_tSOP_vr(:,:)=spval
  allocate(this%h2D_litrC_vr(beg_col:end_col,1:JZ))   ;this%h2D_litrC_vr(:,:)=spval
  allocate(this%h2D_litrN_vr(beg_col:end_col,1:JZ))   ;this%h2D_litrN_vr(:,:)=spval
  allocate(this%h2D_litrP_vr(beg_col:end_col,1:JZ))   ;this%h2D_litrP_vr(:,:)=spval
  allocate(this%h2D_VHeatCap_vr(beg_col:end_col,1:JZ));this%h2D_VHeatCap_vr(:,:)=spval
  allocate(this%h1D_Num_Leaves_ptc(beg_ptc:end_ptc));this%h1D_Num_Leaves_ptc(:)=spval
  allocate(this%h1D_RUB_ACTVN_ptc(beg_ptc:end_ptc));  this%h1D_RUB_ACTVN_ptc(:)=spval
  allocate(this%h1D_CanopyNLim_ptc(beg_ptc:end_ptc)); this%h1D_CanopyNLim_ptc(:)=spval
  allocate(this%h1D_CanopyPLim_ptc(beg_ptc:end_ptc)); this%h1D_CanopyPLim_ptc(:)=spval
  allocate(this%h2D_Aqua_N2_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_N2_vr(:,:)=spval
  allocate(this%h2D_Aqua_H2_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_H2_vr(:,:)=spval
  allocate(this%h2D_Aqua_Ar_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_Ar_vr(:,:)=spval
  allocate(this%h2D_Aqua_CO2_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_CO2_vr(:,:)=spval
  allocate(this%h2D_Root_CO2_vr(beg_col:end_col,1:JZ))  ; this%h2D_Root_CO2_vr(:,:)=spval
  allocate(this%h2D_Root_O2_vr(beg_col:end_col,1:JZ)); this%h2D_Root_O2_vr(:,:)=spval
  allocate(this%h2D_Root_CH4_vr(beg_col:end_col,1:JZ)); this%h2D_Root_CH4_vr(:,:)=spval
  allocate(this%h2D_Root_H2_vr(beg_col:end_col,1:JZ)); this%h2D_Root_H2_vr(:,:)=spval
  allocate(this%h2D_RootNutupk_fClim_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootNutupk_fClim_pvr=spval
  allocate(this%h2D_RootNutupk_fNlim_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootNutupk_fNlim_pvr=spval
  allocate(this%h2D_RootNutupk_fPlim_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootNutupk_fPlim_pvr=spval
  allocate(this%h2D_RootNutupk_fProtC_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootNutupk_fProtC_pvr=spval
  allocate(this%h2D_Root1stSArea4GasTP_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_Root1stSArea4GasTP_pvr=spval
  allocate(this%h2D_RootProteinC_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_RootProteinC_pvr=spval
  allocate(this%h2D_O2_rootconduct_pvr(beg_ptc:end_ptc,1:JZ))     ;this%h2D_O2_rootconduct_pvr(:,:)=spval
  allocate(this%h2D_CO2_rootconduct_pvr(beg_ptc:end_ptc,1:JZ))     ;this%h2D_CO2_rootconduct_pvr(:,:)=spval
  allocate(this%h2D_ProteinNperm2LeafArea_pnd(beg_ptc:end_ptc,1:MaxNodesPerBranch)); this%h2D_ProteinNperm2LeafArea_pnd(:,:)=spval
  allocate(this%h2D_Aqua_CH4_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_CH4_vr(:,:)=spval
  allocate(this%h2D_Aqua_O2_vr(beg_col:end_col,1:JZ))         ;this%h2D_Aqua_O2_vr(:,:)=spval
  allocate(this%h2D_Aqua_N2O_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_N2O_vr(:,:)=spval
  allocate(this%h2D_Aqua_NH3_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_NH3_vr(:,:)=spval
  allocate(this%h2D_TEMP_vr(beg_col:end_col,1:JZ))       ;this%h2D_TEMP_vr(:,:)=spval
  allocate(this%h1D_decomp_OStress_litr_col(beg_col:end_col)); this%h1D_decomp_OStress_litr_col(:)=spval
  allocate(this%h1D_Decomp_temp_FN_litr_col(beg_col:end_col)); this%h1D_Decomp_temp_FN_litr_col(:)=spval
  allocate(this%h1D_FracLitMix_litr_col(beg_col:end_col)); this%h1D_FracLitMix_litr_col(:)=spval
  allocate(this%h1D_Decomp_Moist_FN_litr_col(beg_col:end_col)); this%h1D_Decomp_Moist_FN_litr_col(:)=spval
  allocate(this%h1D_RO2Decomp_litr_col(beg_col:end_col)); this%h1D_RO2Decomp_litr_col(:)=spval
  allocate(this%h2D_decomp_OStress_vr(beg_col:end_col,1:JZ)); this%h2D_decomp_OStress_vr(:,:)=spval
  allocate(this%h2D_RO2Decomp_vr(beg_col:end_col,1:JZ)); this%h2D_RO2Decomp_vr(:,:)=spval
  allocate(this%h2D_Decomp_temp_FN_vr(beg_col:end_col,1:JZ)); this%h2D_Decomp_temp_FN_vr(:,:)=spval
  allocate(this%h2D_FracLitMix_vr(beg_col:end_col,1:JZ)); this%h2D_FracLitMix_vr(:,:)=spval
  allocate(this%h2D_Decomp_Moist_FN_vr(beg_col:end_col,1:JZ)); this%h2D_Decomp_Moist_FN_vr(:,:)=spval
  allocate(this%h2D_RootMassC_vr(beg_col:end_col,1:JZ)); this%h2D_RootMassC_vr(:,:)=spval
  allocate(this%h2D_RootMassN_vr(beg_col:end_col,1:JZ)); this%h2D_RootMassN_vr(:,:)=spval
  allocate(this%h2D_RootMassP_vr(beg_col:end_col,1:JZ)); this%h2D_RootMassP_vr(:,:)=spval
  allocate(this%h2D_RootMassC_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootMassC_pvr(:,:)=spval
  allocate(this%h2D_RootRadialKond2H2O_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootRadialKond2H2O_pvr(:,:)=spval
  allocate(this%h2d_RootPop_pvr(beg_ptc:end_ptc,1:JZ)); this%h2d_RootPop_pvr(:,:)=spval
  allocate(this%h2D_MycoPop_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_MycoPop_pvr(:,:)=spval
  allocate(this%h2D_Root1stDepz_ptc(beg_ptc:end_ptc,1:MaxNumRootAxes)); this%h2D_Root1stDepz_ptc(:,:)=spval
  allocate(this%h2D_RootAxialKond2H2O_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootAxialKond2H2O_pvr(:,:)=spval
  allocate(this%h2D_VmaxNH4Root_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_VmaxNH4Root_pvr(:,:)=spval
  allocate(this%h2D_VmaxNO3Root_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_VmaxNO3Root_pvr(:,:)=spval
  allocate(this%h2D_DOC_vr(beg_col:end_col,1:JZ)); this%h2D_DOC_vr(:,:)=spval
  allocate(this%h2D_DON_vr(beg_col:end_col,1:JZ)); this%h2D_DON_vr(:,:)=spval
  allocate(this%h2D_DOP_vr(beg_col:end_col,1:JZ)); this%h2D_DOP_vr(:,:)=spval
  allocate(this%h2D_SoilBulkStress_vr(beg_col:end_col,1:JZ)); this%h2D_SoilBulkStress_vr(:,:)=spval
  allocate(this%h2D_acetate_vr(beg_col:end_col,1:JZ)); this%h2D_acetate_vr(:,:)=spval
  allocate(this%h1D_tDOC_soil_col(beg_col:end_col)); this%h1D_tDOC_soil_col(:)=spval
  allocate(this%h1D_tDON_soil_col(beg_col:end_col)); this%h1D_tDON_soil_col(:)=spval
  allocate(this%h1D_tDOP_soil_col(beg_col:end_col)); this%h1D_tDOP_soil_col(:)=spval
  allocate(this%h1D_tAcetate_soil_col(beg_col:end_col)); this%h1D_tAcetate_soil_col(:)=spval
  allocate(this%h1D_DOC_litr_col(beg_col:end_col));  this%h1D_DOC_litr_col(:)=spval
  allocate(this%h1D_DON_litr_col(beg_col:end_col)); this%h1D_DON_litr_col(:)=spval
  allocate(this%h1D_DOP_litr_col(beg_col:end_col)); this%h1D_DOP_litr_col(:)=spval
  allocate(this%h1D_acetate_litr_col(beg_col:end_col)); this%h1D_acetate_litr_col(:)=spval
  allocate(this%h1D_VHeatCap_litr_col(beg_col:end_col)); this%h1D_VHeatCap_litr_col(:)=spval
  allocate(this%h2D_Gas_Pressure_vr(beg_col:end_col,1:JZ)); this%h2D_Gas_Pressure_vr(:,:)=spval
  allocate(this%h2D_CO2_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_CO2_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_CH4_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_CH4_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_H2_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_H2_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_Ar_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_Ar_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_N2_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_N2_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_N2O_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_N2O_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_NH3_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_NH3_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_O2_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_O2_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_cyanoBactC_vr(beg_col:end_col,1:JZ)); this%h2D_cyanoBactC_vr(:,:)=spval
  allocate(this%h2D_AeroHrBactC_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrBactC_vr(:,:)=spval
  allocate(this%h2D_AeroHrFungC_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrFungC_vr(:,:)=spval
  allocate(this%h2D_faculDenitC_vr(beg_col:end_col,1:JZ)); this%h2D_faculDenitC_vr(:,:)=spval
  allocate(this%h2D_fermentorC_vr(beg_col:end_col,1:JZ));  this%h2D_fermentorC_vr(:,:)=spval

  allocate(this%h2D_fermentor_frac_vr(beg_col:end_col,1:JZ)); this%h2D_fermentor_frac_vr(:,:)=spval
  allocate(this%h2D_acetometh_frac_vr(beg_col:end_col,1:JZ)); this%h2D_acetometh_frac_vr(:,:)=spval
  allocate(this%h2D_hydrogMeth_frac_vr(beg_col:end_col,1:JZ)); this%h2D_hydrogMeth_frac_vr(:,:)=spval
  allocate(this%h2D_HydCondSoil_vr(beg_col:end_col,1:JZ)); this%h2D_HydCondSoil_vr=spval
  allocate(this%h2D_acetometgC_vr(beg_col:end_col,1:JZ));  this%h2D_acetometgC_vr(:,:)=spval
  allocate(this%h2D_aeroN2fixC_vr(beg_col:end_col,1:JZ));  this%h2D_aeroN2fixC_vr(:,:)=spval
  allocate(this%h2D_anaeN2FixC_vr(beg_col:end_col,1:JZ));  this%h2D_anaeN2FixC_vr(:,:)=spval
  allocate(this%h2D_NH3OxiBactC_vr(beg_col:end_col,1:JZ)); this%h2D_NH3OxiBactC_vr(:,:)=spval
  allocate(this%h2D_NO2OxiBactC_vr(beg_col:end_col,1:JZ)); this%h2D_NO2OxiBactC_vr(:,:)=spval
  allocate(this%h2D_CH4AeroOxiC_vr(beg_col:end_col,1:JZ)); this%h2D_CH4AeroOxiC_vr(:,:)=spval
  allocate(this%h2D_H2MethogenC_vr(beg_col:end_col,1:JZ)); this%h2D_H2MethogenC_vr(:,:)=spval
  allocate(this%h2D_tOMActCDens_vr(beg_col:end_col,1:JZ)); this%h2D_tOMActCDens_vr(:,:)=spval
  allocate(this%h2D_TSolidOMActC_vr(beg_col:end_col,1:JZ));this%h2D_TSolidOMActC_vr(:,:)=spval
  allocate(this%h1D_tOMActCDens_litr_col(beg_col:end_col)); this%h1D_tOMActCDens_litr_col(:)=spval
  allocate(this%h1D_TSolidOMActC_litr_col(beg_col:end_col)); this%h1D_TSolidOMActC_litr_col(:)=spval
  allocate(this%h2D_TSolidOMActCDens_vr(beg_col:end_col,1:JZ));this%h2D_TSolidOMActCDens_vr(:,:)=spval
  allocate(this%h1D_TSolidOMActCDens_litr_col(beg_col:end_col)); this%h1D_TSolidOMActCDens_litr_col(:)=spval
  allocate(this%h2D_RCH4ProdHydrog_vr(beg_col:end_col,1:JZ));   this%h2D_RCH4ProdHydrog_vr(:,:)=spval
  allocate(this%h2D_RCH4ProdAcetcl_vr(beg_col:end_col,1:JZ));   this%h2D_RCH4ProdAcetcl_vr(:,:)=spval
  allocate(this%h2D_RCH4Oxi_aero_vr(beg_col:end_col,1:JZ)); this%h2D_RCH4Oxi_aero_vr(:,:)=spval
  allocate(this%h2D_RCH4Oxi_anmo_vr(beg_col:end_col,1:JZ)); this%h2D_RCH4Oxi_anmo_vr(:,:)=spval
  allocate(this%h2D_RFerment_vr(beg_col:end_col,1:JZ)); this%h2D_RFerment_vr(:,:)=spval
  allocate(this%h2D_nh3oxi_vr(beg_col:end_col,1:JZ));  this%h2D_nh3oxi_vr(:,:)=spval
  allocate(this%h2D_RNit_NO2toNO3_vr(beg_col:end_col,1:JZ));this%h2D_RNit_NO2toNO3_vr(:,:)=spval
  allocate(this%h2D_RDen_NO3toNO2_vr(beg_col:end_col,1:JZ)); this%h2D_RDen_NO3toNO2_vr(:,:)=spval
  allocate(this%h2D_n2oprod_vr(beg_col:end_col,1:JZ));  this%h2D_n2oprod_vr(:,:)=spval
  allocate(this%h2D_RootAR_vr(beg_col:end_col,1:JZ)); this%h2D_RootAR_vr(:,:)=spval
  allocate(this%h2D_PAR_RAD_vr(beg_col:end_col,1:JZ)); this%h2D_PAR_RAD_vr(:,:)=spval
  allocate(this%h2D_RootAR2soil_vr(beg_col:end_col,1:JZ)); this%h2D_RootAR2soil_vr(:,:)=spval
  allocate(this%h2D_RootAR2Root_vr(beg_col:end_col,1:JZ)); this%h2D_RootAR2Root_vr(:,:)=spval
  allocate(this%h1D_RCH4ProdHydrog_litr_col(beg_col:end_col));  this%h1D_RCH4ProdHydrog_litr_col(:)=spval
  allocate(this%h1D_RCH4ProdAcetcl_litr_col(beg_col:end_col));  this%h1D_RCH4ProdAcetcl_litr_col(:)=spval
  allocate(this%h1D_RCH4Oxi_aero_litr_col(beg_col:end_col)); this%h1D_RCH4Oxi_aero_litr_col(:)=spval
  allocate(this%h1D_RCH4Oxi_anmo_litr_col(beg_col:end_col)); this%h1D_RCH4Oxi_anmo_litr_col(:)=spval
  allocate(this%h1D_RCH4Oxi_aero_col(beg_col:end_col));this%h1D_RCH4Oxi_aero_col(:)=spval
  allocate(this%h1D_RCH4Oxi_anmo_col(beg_col:end_col));this%h1D_RCH4Oxi_anmo_col(:)=spval
  allocate(this%h1D_RFermen_litr_col(beg_col:end_col));  this%h1D_RFermen_litr_col(:)=spval
  allocate(this%h1D_nh3oxi_litr_col(beg_col:end_col)); this%h1D_nh3oxi_litr_col(:)=spval
  allocate(this%h1D_RDen_NO3toNO2_col(beg_col:end_col));this%h1D_RDen_NO3toNO2_col(:)=spval
  allocate(this%h1D_NO2Oxi_litr_col(beg_col:end_col)); this%h1D_NO2Oxi_litr_col(:)=spval
  allocate(this%h1D_n2oprod_litr_col(beg_col:end_col));  this%h1D_n2oprod_litr_col(:)=spval
  allocate(this%h2D_Gchem_CO2_prod_vr(beg_col:end_col,1:JZ)); this%h2D_Gchem_CO2_prod_vr(:,:)=spval
  allocate(this%h2D_Eco_HR_CO2_vr(beg_col:end_col,1:JZ)); this%h2D_Eco_HR_CO2_vr(:,:)=spval
  allocate(this%h2D_AeroHrBactN_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrBactN_vr(:,:)=spval
  allocate(this%h2D_AeroHrFungN_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrFungN_vr(:,:)=spval
  allocate(this%h2D_faculDenitN_vr(beg_col:end_col,1:JZ)); this%h2D_faculDenitN_vr(:,:)=spval
  allocate(this%h2D_fermentorN_vr(beg_col:end_col,1:JZ));  this%h2D_fermentorN_vr(:,:)=spval
  allocate(this%h2D_acetometgN_vr(beg_col:end_col,1:JZ));  this%h2D_acetometgN_vr(:,:)=spval
  allocate(this%h2D_aeroN2fixN_vr(beg_col:end_col,1:JZ));  this%h2D_aeroN2fixN_vr(:,:)=spval
  allocate(this%h2D_anaeN2FixN_vr(beg_col:end_col,1:JZ));  this%h2D_anaeN2FixN_vr(:,:)=spval
  allocate(this%h2D_NH3OxiBactN_vr(beg_col:end_col,1:JZ)); this%h2D_NH3OxiBactN_vr(:,:)=spval
  allocate(this%h2D_NO2OxiBactN_vr(beg_col:end_col,1:JZ)); this%h2D_NO2OxiBactN_vr(:,:)=spval
  allocate(this%h2D_CH4AeroOxiN_vr(beg_col:end_col,1:JZ)); this%h2D_CH4AeroOxiN_vr(:,:)=spval
  allocate(this%h2D_H2MethogenN_vr(beg_col:end_col,1:JZ)); this%h2D_H2MethogenN_vr(:,:)=spval
  allocate(this%h2D_tRespGrossHeterUlm_vr(beg_col:end_col,1:JZ)); this%h2D_tRespGrossHeterUlm_vr(:,:)=spval
  allocate(this%h2D_tRespGrossHeter_vr(beg_col:end_col,1:JZ)); this%h2D_tRespGrossHeter_vr(:,:)=spval

  allocate(this%h2D_RDECOMPC_SOM_vr(beg_col:end_col,1:JZ)); this%h2D_RDECOMPC_SOM_vr(:,:)=spval
  allocate(this%h2D_RDECOMPC_BReSOM_vr(beg_col:end_col,1:JZ));this%h2D_RDECOMPC_BReSOM_vr(:,:)=spval
  allocate(this%h2D_RDECOMPC_SorpSOM_vr(beg_col:end_col,1:JZ));this%h2D_RDECOMPC_SorpSOM_vr(:,:)=spval
  allocate(this%h2D_MicrobAct_vr(beg_col:end_col,1:JZ)); this%h2D_MicrobAct_vr(:,:)=spval

  allocate(this%h3D_MicrobActCps_vr(beg_col:end_col,1:JZ,jcplx)); this%h3D_MicrobActCps_vr(:,:,:)=spval

  allocate(this%h3D_SOMHydrylScalCps_vr(beg_col:end_col,1:JZ,jcplx));this%h3D_SOMHydrylScalCps_vr(:,:,:)=spval

  allocate(this%h3D_SOC_Cps_vr(beg_col:end_col,1:JZ,jcplx)); this%h3D_SOC_Cps_vr(:,:,:)=spval

  allocate(this%h3D_HydrolCSOMCps_vr(beg_col:end_col,1:JZ,jcplx)); this%h3D_HydrolCSOMCps_vr(:,:,:)=spval

  allocate(this%h1D_RDECOMPC_SOM_litr_col(beg_col:end_col)); this%h1D_RDECOMPC_SOM_litr_col(:)=spval
  allocate(this%h1D_RDECOMPC_BReSOM_litr_col(beg_col:end_col));this%h1D_RDECOMPC_BReSOM_litr_col(:)=spval
  allocate(this%h1D_RDECOMPC_SorpSOM_litr_col(beg_col:end_col));this%h1D_RDECOMPC_SorpSOM_litr_col(:)=spval
  allocate(this%h1D_MicrobAct_litr_col(beg_col:end_col)); this%h1D_MicrobAct_litr_col(:)=spval
  allocate(this%h1D_tRespGrossHeteUlm_litr_col(beg_col:end_col)); this%h1D_tRespGrossHeteUlm_litr_col=spval
  allocate(this%h1D_tRespGrossHete_litr_col(beg_col:end_col)); this%h1D_tRespGrossHete_litr_col=spval

  allocate(this%h2D_AeroHrBactP_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrBactP_vr(:,:)=spval
  allocate(this%h2D_AeroHrFungP_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrFungP_vr(:,:)=spval
  allocate(this%h2D_faculDenitP_vr(beg_col:end_col,1:JZ)); this%h2D_faculDenitP_vr(:,:)=spval
  allocate(this%h2D_fermentorP_vr(beg_col:end_col,1:JZ));  this%h2D_fermentorP_vr(:,:)=spval
  allocate(this%h2D_acetometgP_vr(beg_col:end_col,1:JZ));  this%h2D_acetometgP_vr(:,:)=spval
  allocate(this%h2D_aeroN2fixP_vr(beg_col:end_col,1:JZ));  this%h2D_aeroN2fixP_vr(:,:)=spval
  allocate(this%h2D_anaeN2FixP_vr(beg_col:end_col,1:JZ));  this%h2D_anaeN2FixP_vr(:,:)=spval
  allocate(this%h2D_NH3OxiBactP_vr(beg_col:end_col,1:JZ)); this%h2D_NH3OxiBactP_vr(:,:)=spval
  allocate(this%h2D_NO2OxiBactP_vr(beg_col:end_col,1:JZ)); this%h2D_NO2OxiBactP_vr(:,:)=spval
  allocate(this%h2D_CH4AeroOxiP_vr(beg_col:end_col,1:JZ)); this%h2D_CH4AeroOxiP_vr(:,:)=spval
  allocate(this%h2D_H2MethogenP_vr(beg_col:end_col,1:JZ)); this%h2D_H2MethogenP_vr(:,:)=spval

  allocate(this%h2D_MicroBiomeE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_MicroBiomeE_litr_col(:,:)=spval
  allocate(this%h2D_AeroHrBactE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_AeroHrBactE_litr_col(:,:)=spval
  allocate(this%h2D_AeroHrFungE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_AeroHrFungE_litr_col(:,:)=spval
  allocate(this%h2D_faculDenitE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_faculDenitE_litr_col(:,:)=spval
  allocate(this%h2D_fermentorE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_fermentorE_litr_col(:,:)=spval
  allocate(this%h2D_cyanoBactC_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_cyanoBactC_litr_col(:,:)=spval
  allocate(this%h2D_acetometgE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_acetometgE_litr_col(:,:)=spval
  allocate(this%h2D_aeroN2fixE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_aeroN2fixE_litr_col(:,:)=spval
  allocate(this%h2D_anaeN2FixE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_anaeN2FixE_litr_col(:,:)=spval
  allocate(this%h2D_NH3OxiBactE_litr_col(beg_col:end_col,1:NumPlantChemElms));this%h2D_NH3OxiBactE_litr_col(:,:)=spval
  allocate(this%h2D_NO2OxiBactE_litr_col(beg_col:end_col,1:NumPlantChemElms));this%h2D_NO2OxiBactE_litr_col(:,:)=spval
  allocate(this%h2D_CH4AeroOxiE_litr_col(beg_col:end_col,1:NumPlantChemElms));this%h2D_CH4AeroOxiE_litr_col(:,:)=spval
  allocate(this%h2D_ANME2dC_vr(beg_col:end_col,1:JZ)); this%h2D_ANME2dC_vr(:,:)=spval
  allocate(this%h2D_ANME2dN_vr(beg_col:end_col,1:JZ)); this%h2D_ANME2dN_vr(:,:)=spval
  allocate(this%h2D_ANME2dP_vr(beg_col:end_col,1:JZ)); this%h2D_ANME2dP_vr(:,:)=spval
  allocate(this%h2D_ANME2dE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_ANME2dE_litr_col(:,:)=spval
  allocate(this%h2D_NC10C_vr(beg_col:end_col,1:JZ)); this%h2D_NC10C_vr(:,:)=spval
  allocate(this%h2D_NC10N_vr(beg_col:end_col,1:JZ)); this%h2D_NC10N_vr(:,:)=spval
  allocate(this%h2D_NC10P_vr(beg_col:end_col,1:JZ)); this%h2D_NC10P_vr(:,:)=spval
  allocate(this%h2D_NC10E_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_NC10E_litr_col(:,:)=spval
  allocate(this%h2D_H2MethogenE_litr_col(beg_col:end_col,1:NumPlantChemElms));this%h2D_H2MethogenE_litr_col(:,:)=spval


  allocate(this%h2D_HeatFlow_vr(beg_col:end_col,1:JZ))   ;this%h2D_HeatFlow_vr(:,:)=spval
  allocate(this%h2D_HeatUptk_vr(beg_col:end_col,1:JZ))   ;this%h2D_HeatUptk_vr(:,:)=spval
  allocate(this%h2D_rVSM_vr    (beg_col:end_col,1:JZ))     ;this%h2D_rVSM_vr    (:,:)=spval
  allocate(this%h2D_VSPore_vr (beg_col:end_col,1:JZ))     ;this%h2D_VSPore_vr (:,:)=spval
  allocate(this%h2D_FLO_MICP_vr(beg_col:end_col,1:JZ)) ;this%h2D_FLO_MICP_vr(:,:)=spval
  allocate(this%h2D_FLO_MACP_vr(beg_col:end_col,1:JZ)) ;this%h2D_FLO_MACP_vr(:,:)=spval
  allocate(this%h2D_rVSICE_vr    (beg_col:end_col,1:JZ))       ;this%h2D_rVSICE_vr    (:,:)=spval
  allocate(this%h2D_PSI_vr(beg_col:end_col,1:JZ))        ;this%h2D_PSI_vr(:,:)=spval
  allocate(this%h2D_PsiO_vr(beg_col:end_col,1:JZ)) ; this%h2D_PsiO_vr(:,:)=spval
  allocate(this%h2D_RootH2OUP_vr(beg_col:end_col,1:JZ))  ;this%h2D_RootH2OUP_vr(:,:)=spval
  allocate(this%h2D_cNH4t_vr(beg_col:end_col,1:JZ))      ;this%h2D_cNH4t_vr(:,:)=spval
  allocate(this%h2D_cNO3t_vr(beg_col:end_col,1:JZ))      ;this%h2D_cNO3t_vr(:,:)=spval

  allocate(this%h2D_cPO4_vr(beg_col:end_col,1:JZ))       ;this%h2D_cPO4_vr(:,:)=spval
  allocate(this%h2D_cEXCH_P_vr(beg_col:end_col,1:JZ))    ;this%h2D_cEXCH_P_vr(:,:)=spval
  allocate(this%h2D_microb_N2fix_vr(beg_col:end_col,1:JZ)); this%h2D_microb_N2fix_vr(:,:)=spval
  allocate(this%h2D_ElectricConductivity_vr(beg_col:end_col,1:JZ))       ;this%h2D_ElectricConductivity_vr(:,:)=spval
  allocate(this%h2D_PSI_RT_pvr(beg_ptc:end_ptc,1:JZ))     ;this%h2D_PSI_RT_pvr(:,:)=spval
  allocate(this%h2D_RootH2OUptkStress_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootH2OUptkStress_pvr(:,:)=spval
  allocate(this%h2D_RootH2OUptk_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_RootH2OUptk_pvr(:,:)=spval
  allocate(this%h2D_RootAct1stC_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_RootAct1stC_pvr(:,:)=spval
  allocate(this%h2D_RootMedC_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_RootMedC_pvr(:,:)=spval
  allocate(this%h2D_NonstC_conc_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_NonstC_conc_pvr(:,:)=spval
  allocate(this%h2D_RootLig1stC_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_RootLig1stC_pvr(:,:)=spval
  allocate(this%h2D_RootShootExchC_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootShootExchC_pvr(:,:)=spval
  allocate(this%h2D_RootShootExchN_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootShootExchN_pvr(:,:)=spval
  allocate(this%h2D_RootShootExchP_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootShootExchP_pvr(:,:)=spval
  allocate(this%h2D_SapFlowVlinear_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_SapFlowVlinear_pvr(:,:)=spval
  allocate(this%h2D_RootMaintDef_CO2_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootMaintDef_CO2_pvr(:,:)=spval
  allocate(this%h2D_RootAbsorbAreaPP_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootAbsorbAreaPP_pvr(:,:)=spval
  allocate(this%h1D_RootAbsorbAreaPP_pft(beg_ptc:end_ptc)); this%h1D_RootAbsorbAreaPP_pft(:)=spval
  allocate(this%h2D_ROOTNLim_rpvr(beg_ptc:end_ptc,1:JZ)); this%h2D_ROOTNLim_rpvr(:,:)=spval
  allocate(this%h2D_ROOTPLim_rpvr(beg_ptc:end_ptc,1:JZ)); this%h2D_ROOTPLim_rpvr(:,:)=spval
  allocate(this%h2D_RootNonstC_rpvr(beg_ptc:end_ptc,1:JZ)); this%h2D_RootNonstC_rpvr(:,:)=spval
  allocate(this%h2D_RootMSinkWeight_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_RootMSinkWeight_pvr(:,:)=spval
  allocate(this%h2D_RootSinkWeight_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_RootSinkWeight_pvr(:,:)=spval
  allocate(this%h2D_Root2ndSinkWeight_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_Root2ndSinkWeight_pvr(:,:)=spval
  allocate(this%h2D_Root1stSinkWeight_pvr(beg_ptc:end_ptc,1:jz));this%h2D_Root1stSinkWeight_pvr(:,:)=spval
  allocate(this%h2D_Root1stRadius_rpvr(beg_ptc:end_ptc,1:JZ)); this%h2D_Root1stRadius_rpvr(:,:)=spval
  allocate(this%h2D_ROOT_OSTRESS_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_ROOT_OSTRESS_pvr(:,:)=spval
  allocate(this%h2D_prtUP_NH4_pvr(beg_ptc:end_ptc,1:JZ))  ;this%h2D_prtUP_NH4_pvr(:,:)=spval
  allocate(this%h2D_prtUP_NO3_pvr(beg_ptc:end_ptc,1:JZ))  ;this%h2D_prtUP_NO3_pvr(:,:)=spval
  allocate(this%h2D_prtUP_PO4_pvr(beg_ptc:end_ptc,1:JZ))  ;this%h2D_prtUP_PO4_pvr(:,:)=spval
  allocate(this%h2D_DNS_RT_pvr(beg_ptc:end_ptc,1:JZ))     ;this%h2D_DNS_RT_pvr(:,:)=spval
  allocate(this%h2D_RootNonstBConc_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootNonstBConc_pvr(:,:)=spval
  allocate(this%h2D_MycoBiomC_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_MycoBiomC_pvr(:,:)=spval
  allocate(this%h2D_Root1stStrutC_pvr(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root1stStrutC_pvr=spval
  allocate(this%h2D_Cytok_scalar_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_Cytok_scalar_pvr(:,:)=spval
  allocate(this%h2D_Cytokinin1stConc_pvr(beg_ptc:end_ptc,1:JZ));  this%h2D_Cytokinin1stConc_pvr(:,:)=spval
  allocate(this%h2D_CRootLumenArea_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_CRootLumenArea_pvr(:,:)=spval
  allocate(this%h2D_Root1stStrutN_pvr(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root1stStrutN_pvr=spval
  allocate(this%h2D_Root1stStrutP_pvr(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root1stStrutP_pvr=spval
  allocate(this%h2D_Root2ndStrutC_pvr(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root2ndStrutC_pvr=spval
  allocate(this%h2D_Root2ndStrutN_pvr(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root2ndStrutN_pvr=spval
  allocate(this%h2D_Root2ndStrutP_pvr(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root2ndStrutP_pvr=spval
  allocate(this%h2D_Root2ndAxesNumL_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_Root2ndAxesNumL_pvr=spval
  allocate(this%h2D_Root2ndMeanLen_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_Root2ndMeanLen_pvr=spval
  allocate(this%h2D_RootKond2H2O_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootKond2H2O_pvr=spval
  allocate(this%h2D_Root1stLenPP_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_Root1stLenPP_pvr=spval
  allocate(this%h2D_Rootmedlength_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_Rootmedlength_pvr=spval
  allocate(this%h2D_RootmedRadius_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootmedRadius_pvr=spval
  allocate(this%h2D_Root1stAxesNumL_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_Root1stAxesNumL_pvr=spval
  allocate(this%h2D_RootMedAxesNumL_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootMedAxesNumL_pvr=spval
  allocate(this%h2D_fTRootGro_pvr(beg_ptc:end_ptc,1:JZ)) ; this%h2D_fTRootGro_pvr=spval
  allocate(this%h2D_fRootGrowPSISense_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_fRootGrowPSISense_pvr=spval
  allocate(this%h3D_PARTS_ptc(beg_ptc:end_ptc,1:NumOfPlantMorphUnits,1:MaxNumBranches));this%h3D_PARTS_ptc(:,:,:)=spval
  allocate(this%h2D_CanopyLAIZ_plyr(beg_ptc:end_ptc,1:NumCanopyLayers)); this%h2D_CanopyLAIZ_plyr(:,:)=spval
  end procedure allocate_hist_profiles

end submodule HistAllocateProfile
