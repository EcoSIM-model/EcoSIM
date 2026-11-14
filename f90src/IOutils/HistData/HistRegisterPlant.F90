submodule (HistDataType) HistRegisterPlant
  ! Register PFT-level plant diagnostics with history metadata and pointers.
  use HistFileMod, only: hist_addfld1d
  implicit none
contains

  module procedure register_hist_plants
    ! Expose supported plant history fields to the output system.
    integer :: beg_col, end_col, beg_ptc, end_ptc
    real(r8), pointer :: data1d_ptr(:)

    beg_col=1; end_col=bounds%ncols
    beg_ptc=1; end_ptc=bounds%npfts

  data1d_ptr => this%h1D_CAN_G_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_G_pft',units='W/m2',avgflag='A',&
    long_name='Canopy storage heat flux',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CAN_TEMPC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_TEMPC_pft',units='oC',avgflag='A',&
    long_name='Canopy temperature',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_CAN_TEMPFN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CANGRO_TEMP_FN_pft',units='none',avgflag='A',&
    long_name='Canopy temperature growth function/stress',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CAN_CO2_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_CO2_FLX_pft',units='umol C/m2/s',avgflag='A',&
    long_name='Canopy net CO2 exchange',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CAN_GPP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_GPP_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant canopy gross CO2 fixation',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_CAN_cumGPP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_cumGPP_pft',units='gC/m2',avgflag='I',&
    long_name='Plant canopy cumulative gross CO2 fixation',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_dCAN_GPP_CLIM_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='dCAN_GPP_CLIM_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant canopy CO2-limited gross CO2 fixation minus acutal value (>0 light-limitation)',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_dCAN_GPP_eLIM_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='dCAN_GPP_eLIM_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant canopy light-limited gross CO2 fixation minus actual value (>0 C-limitation)',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CAN_RA_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_RA_pft',units='gC/m2/hr',avgflag='A',&
    long_name='total aboveground autotrophic respiration',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_TurgEff4CanopyResp_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='TurgEff4CanopyResp_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Turgor pressure effect on canopy autotrophic respiration',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CAN_GROWTH_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_GROWTH_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Canopy structural growth rate',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_cTNC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='cAbvNonstC_pft',units='gC/gC',avgflag='A',&
    long_name='Canopy nonstructural C concentration',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_cTNN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='cAbvNonstN_pft',units='gN/gC',avgflag='A',&
    long_name='Canopy nonstructural N concentration',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_cTNP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='cAbvNonstP_pft',units='gP/gC',avgflag='A',&
    long_name='Canopy nonstructural P concentration',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_fSnowCan_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='fSnowCanopy_pft',units='-',avgflag='A',&
    long_name='Canopy covered by snow',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_PTSHTR_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='PTSHTR_pft',units='h-1',avgflag='A',&
    long_name='Root-shoot C coupling rate',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_CanopyNonstOMConc_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CanopyNonstOMConc_pft',units='g/gC',avgflag='A',&
    long_name='Canopy nonstructural C+N+P concentration per unit leaf/sheath carbon',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_STOML_RSC_CO2_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='STOML_RSC_CO2_pft',units='s/m',avgflag='A',&
    long_name='Canopy stomatal resistance for CO2',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_STOML_Min_RSC_CO2_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='STOML_MinRSC_CO2_pft',units='s/m',avgflag='A',&
    long_name='Canopy minimal stomatal resistance for CO2',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_Km_CO2_carboxy_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Km_CO2_carboxy_pft',units='uM',avgflag='A',&
    long_name='MM parameter for CO2 carboxylation by Rubisco',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_DynCi2CaRatio_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='DynCi2CaRatio_pft',units='1',avgflag='A',&
    long_name='Dynamic intracellular-to-canopy CO2 ratio',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_BLYR_RSC_CO2_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='BLYR_RSC_CO2_pft',units='s/m',avgflag='A',&
    long_name='Canopy boundary layer resistance for CO2',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CAN_CO2_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_CO2_pft',units='umol/mol',avgflag='A',&
    long_name='Canopy gaesous CO2 concentration',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_O2L_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Leaf_O2_pft',units='umol/mol',avgflag='A',&
    long_name='Leaf aqueous O2 concentration',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_LAI_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LAIstk_pft',units='m2/m2',avgflag='A',&
    long_name='whole plant leaf area, including stalk',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_PSI_CAN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='PSI_CAN_pft',units='MPa',avgflag='A',&
    long_name='Canopy total water potential',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_CanPhenol_WSTRSS_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CanPhenol_WSTRESS_pft',units='-',avgflag='A',&
    long_name='Canopy water stress for phenology development',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CanPhenol_TSTRSS_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CanPhenol_TSTRESS_pft',units='-',avgflag='A',&
    long_name='Canopy temperature stress for phenology development',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_RootAR_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootAR_pft',units='gC/m2/h',avgflag='A',&
    long_name='Root autotrophic respiraiton',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_RootLenPerPlant_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootLen_pft',units='m plant-1',avgflag='A',&
    long_name='Root length per pft (excluding root hair)',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_TURG_CAN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='TURG_CAN_pft',units='MPa',avgflag='A',&
    long_name='Canopy turgor water potential',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_STOML_RSC_H2O_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='STOM_RSC_H2O_pft',units='s/m',avgflag='A',&
    long_name='Canopy stomatal resistance for H2O',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_BLYR_RS_H2O_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='BLYR_RS_H2O_pft',units='s/m',avgflag='A',&
    long_name='Canopy boundary layer resistance for H2O',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_BLYR_RSC_HEAT_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='BLYR_RSC_HEAT_pft',units='s/m',avgflag='A',&
    long_name='Canopy boundary layer resistance for sensible heat',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CdH2ORootxSoil_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Cd_Root2Soil_H2O_pft',units='kg H2O m-2 h-1 MPa-1',avgflag='A',&
    long_name='total root soil conductance for plant root H2O uptake',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_TRANSPN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='QVegTransp_pft',units='mmH2O/m2/h',avgflag='A',&
    long_name='Canopy transpiration (<0 into atmosphere)',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_NH4_UPTK_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='UPTK_NH4_FLX_pft',units='gN/m2/hr',&
    avgflag='A',long_name='total root uptake of NH4',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_NO3_UPTK_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='UPTK_NO3_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='total root uptake of NO3',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_N2_FIXN_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='N2_FIX_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='total root N2 fixation',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_cNH3_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_NH3_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='*canopy NH3 flux',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_RNFixCO2_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RNFixCO2_FLX_pft',units='gC/m2/hr',avgflag='A',&
    long_name='CO2 respired for N fixation',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_TC_Groth_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='TC_Groth_pft',units='oC',avgflag='A',&
    long_name='Plant growth temperature',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_PO4_UPTK_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='UPTK_PO4_FLX_pft',units='gP/m2/hr',avgflag='A',&
    long_name='total root uptake of PO4',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_SHOOT_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOT_C_pft',units='gC/m2',avgflag='A',&
    long_name='Live plant shoot C',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_frcPARabs_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='frcPARabs_pft',units='none',avgflag='A',&
    long_name='fraction of PAR absorbed by plant',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_PAR_CAN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Canopy_PAR_pft',units='umol m-2 s-1',avgflag='A',&
    long_name='PAR absorbed by plant canopy',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_Plant_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Plant_C_pft',units='gC/m2',avgflag='A',&
    long_name='Plant C',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_RNodeInitiate_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RNodeAppear_pft',units='1/day',avgflag='A',&
    long_name='Plant node initation rate',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_RLeafAppear_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RLeafAppear_pft',units='1/day',avgflag='A',&
    long_name='Plant leaf appearance rate',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_LEAF_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LEAF_C_pft',units='gC/m2',avgflag='A',&
    long_name='Canopy leaf C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_Petole_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='PetolSheth_C_pft',units='gC/m2',avgflag='A',&
    long_name='Canopy sheath C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_STALK_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='STALK_C_pft',units='gC/m2',avgflag='A',&
    long_name='Canopy stalk C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_RESERVE_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RESERVE_C_pft',units='gC/m2',avgflag='A',&
    long_name='Canopy reserve C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_HUSK_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='HUSK_C_pft',units='gC/m2',avgflag='A',&
    long_name='Canopy husk C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_GRAIN_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='GRAIN_C_pft',units='gC/m2',avgflag='A',&
    long_name='Canopy grain C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ROOT_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root_C_pft',units='gC/m2',avgflag='A',&
    long_name='Plant root C, exluding nodule',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_ROOTST_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootST_C_pft',units='gC/m2',avgflag='A',&
    long_name='Plant root structural C',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_ROOTST_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootST_N_pft',units='gN/m2',avgflag='A',&
    long_name='Plant root structural N',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ROOTST_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootST_P_pft',units='gP/m2',avgflag='A',&
    long_name='Plant root structural P',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_RootNodule_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootNodule_C_pft',units='gC/m2',avgflag='A',&
    long_name='Root total nodule C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ShootNodule_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ShootNodule_C_pft',units='gC/m2',avgflag='A',&
    long_name='Shoot total nodule C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ShootNodule_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ShootNodule_N_pft',units='gN/m2',avgflag='A',&
    long_name='Shoot total nodule N',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ShootNodule_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ShootNodule_P_pft',units='gP/m2',avgflag='A',&
    long_name='Shoot total nodule P',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_STORED_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SSTORED_C_pft',units='gC/m2',avgflag='A',&
    long_name='Plant seasonal storage of nonstructural C',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_ROOT_NONSTC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root_NONSTC_pft',units='gC/m2',avgflag='A',&
    long_name='Plant root storage of nonstructural C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_MycorrizhalBiomC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='MycoArbuBiomC_pft',units='gC/m2',avgflag='A',&
    long_name='arbuscular mycorrhizal structural C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_Root1stStrutC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root1stBiomC_pft',units='gC/m2',avgflag='A',&
    long_name='Primary root structural C',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_RootMeDStrutC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootMedBiomC_pft',units='gC/m2',avgflag='A',&
    long_name='Medium size root structural C',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_Root1stStrutN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root1stBiomN_pft',units='gN/m2',avgflag='A',&
    long_name='Primary root structural N',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_Root2ndStrutC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root2ndBiomC_pft',units='gC/m2',avgflag='A',&
    long_name='Secondary root structural C',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_ROOT_NONSTN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root_NONSTN_pft',units='gN/m2',avgflag='A',&
    long_name='Plant root storage of nonstructural N',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ROOT_NONSTP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root_NONSTP_pft',units='gP/m2',avgflag='A',&
    long_name='Plant root storage of nonstructural P',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SHOOT_NONSTC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOT_NONSTC_pft',units='gC/m2',avgflag='A',&
    long_name='Plant leaf storage of nonstructural C',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_SHOOT_NONSTN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOT_NONSTN_pft',units='gN/m2',avgflag='A',&
    long_name='Plant leaf storage of nonstructural N',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SHOOT_NONSTP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOT_NONSTP_pft',units='gP/m2',avgflag='A',&
    long_name='Plant leaf storage of nonstructural P',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_LeafChlCperm2LA_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafChlC_pft',units='mgC chl m-2 leaf area',avgflag='A',&
    long_name='Total cholorophyll carbon in leaves',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_LeafC4ChlCperm2LA_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafC4ChlC_pft',units='mgC chl m-2 leaf area',avgflag='A',&
    long_name='chlorophyll carbon in mesophyll for C4 plants',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_LeafRubiscoNperm2LA_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafRubiscoN_pft',units='gN rubisco m-2 leaf area',avgflag='A',&
    long_name='Rubisco nitrogen in mesophyll for C3/bundle sheath for C4 plants',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_LeafPEPCNperm2LA_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafPEPCN_pft',units='gN PEP m-2 leaf area',avgflag='A',&
    long_name='PEP carboxylase N in mesophyll for C4 plants',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_GRAIN_NO_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='GRAIN_Number_pft',units='1/m2',avgflag='A',&
    long_name='Number of grains in canopy',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_LAIb_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LAI_xstk_pft',units='m2/m2',avgflag='A',&
    long_name='whole plant leaf area, exclude stalk',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_EXUD_CumYr_C_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='EXUD_CumYr_C_FLX_pft',units='gC/m2',avgflag='I',&
    long_name='Cumulative root organic C uptake (<0 exudation into soil)',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_LITRf_C_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LITRf_C_pft',units='gC/m2/hr',avgflag='A',&
    long_name='total plant LitrFall C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SURF_LITRf_C_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SURF_LITRf_C_FLX_pft',units='gC/m2',avgflag='I',&
    long_name='Cumulative plant LitrFall C to the soil surface',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_AUTO_RESP_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='AUTO_RESP_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Whole plant autotrophic respiration',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_HVST_C_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='HVST_C_FLX_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant C harvest',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_RootAct1stC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root1ActC_pft',units='gC/m2',avgflag='A',&
    long_name='Primary root active C biomass',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_PLANT_BALANCE_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Plant_BALANCE_C_pft',units='gC/m2',avgflag='A',&
    long_name='Cumulative plant C conservation error',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_STANDING_DEAD_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='STANDING_DEAD_C_pft',units='gC/m2',avgflag='A',&
    long_name='pft Standing dead C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_MainBranchNodeNumber_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='MainBranchNodeNumber_pft',units='-',avgflag='I',&
    long_name='Plant main branch node number as a measure of phenoloigcal development',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ShootNodeNumber_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ShootNodeNumber_pft',units='-',avgflag='I',&
    long_name='Plant shoot total node number as a measure of phenoloigcal development',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_FIREp_CO2_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='FIREp_CO2_FLX_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant CO2 from fire',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_FIREp_CH4_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='FIREp_CH4_FLX_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant CH4 emission from fire',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_NPP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='NPP_pft',units='gC/m2',avgflag='A',&
    long_name='Plant net primary productivity',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_CAN_HT_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_HT_pft',units='m',avgflag='A',&
    long_name='Canopy height',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_Stalk_HT_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Stalk_HT_pft',units='m',avgflag='A',&
    long_name='Canopy stalk height',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_EmergeHeight_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='EmergHeight_pft',units='m',avgflag='A',&
    long_name='Hypocotyl height exceeds seeding depth for emergence (>0 yes)',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_POPN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='POPN_pft',units='1/m2',avgflag='A',&
    long_name='Plant population',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_CanopyCutProxy_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CutProxy_pft',units='-',avgflag='A',&
    long_name='Plant cut proxy',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_tTRANSPN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='tTRANSPN_pft',units='mmH2O/m2',avgflag='I',&
    long_name='Cumulative canopy evapotranspiration (>0 into atmosphere)',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_WTR_STRESS_ptc(beg_ptc:end_ptc)    !HoursTooLowPsiCan_pft(NZ,NY,NX)
  call hist_addfld1d(fname='WTR_STRESS_pft',units='hr',avgflag='A',&
    long_name='Canopy plant water stress indicator: number of ' &
    //'hours PSICanopy_pft(< PSILY)',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_OXY_STRESS_ptc(beg_ptc:end_ptc)    !OSTR(NZ,NY,NX)
  call hist_addfld1d(fname='OXY_STRESS_pft',units='none',avgflag='A',&
    long_name='Plant root O2 stress indicator [0->1 weaker stress]',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_VcMaxRubisco_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='VcMax25C_RUBISCO_pft',units='umol CO2 s-1 m-2 leaf area',avgflag='A',&
    long_name='Maximum carboxylation rate by Rubisco at 25oC',&
    ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_LeafProteinNperm2_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafProteinN_pft',units='gN protein m-2 leaf area',avgflag='A',&
    long_name='Protein nitrogen mass per unit of leaf area',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_VoMaxRubisco_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='VoMax25C_RUBISCO_pft',units='umol O2 s-1 m-2 leaf area',avgflag='A',&
    long_name='Maximum oxygenation rate by Rubisco at 25oC',&
    ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_VcMaxPEP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='VcMax25C_PEP_pft',units='umol CO2 s-1 m-2 leaf area',avgflag='A',&
    long_name='Maximum carboxylation rate by PEP at 25oC',&
    ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_JMaxPhoto_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='JMax25C_photo_pft',units='umol e- s-1 m-2 leaf area',avgflag='A',&
    long_name='Maximum electron transport rate at 25oC',&
    ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_TFN_Carboxy_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='TFN_Carboxy_pft',units='none',avgflag='A',&
    long_name='Temperature response of carboyxlation in photosynthesis',&
    ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_TFN_Oxygen_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='TFN_Oxygen_pft',units='none',avgflag='A',&
    long_name='Temperature response of oxygenation in photosynthesis',&
    ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_TFN_eTranspt_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='TFN_eTranspt_pft',units='none',avgflag='A',&
    long_name='Temperature response of electron transport in photosynthesis',&
    ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_PARSunlit_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='PARSunlit_pft',units='umol m-2 s-1',avgflag='A',&
    long_name='PAR absorbed by sunlit leaves',&
    ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_PARSunsha_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='PARSunsha_pft',units='umol m-2 s-1',avgflag='A',&
    long_name='PAR absorbed by sun-shaded leaves',&
    ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_CH2OSunlit_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CH2OSunlit_pft',units='gC m-2 h-1',avgflag='A',&
    long_name='Photosynthesis by sunlit leaves',&
    ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_CH2OSunsha_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CH2OSunsha_pft',units='gC m-2 h-1',avgflag='A',&
    long_name='Photosynthesis by sun-shaded leaves',&
    ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_LeafAreaSunlit_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafArea_sunlit_pft',units='m2 m-2',avgflag='A',&
    long_name='Irridiance-lit leaf area',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_fClump_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='fClump_pft',units='none',avgflag='A',&
    long_name='Clumping factor of leaf area',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SHOOT_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOT_N_pft',units='gN/m2',avgflag='A',&
    long_name='Live plant shoot N',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_Plant_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Plant_N_pft',units='gN/m2',avgflag='A',&
    long_name='Plant N',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_LEAF_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LEAF_N_pft',units='gN/m2',avgflag='A',&
    long_name='Canopy leaf N',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_fCNLFW_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='fLEAF_CN_pft',units='gC/gN',avgflag='A',&
    long_name='Canopy new leaf C:N mass ratio',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_fCPLFW_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='fLEAF_CP_pft',units='gC/gP',avgflag='A',&
    long_name='Canopy new leaf C:P mass ratio',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_LeafNperm2LAI_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafN_pft',units='gN m-2 LA',avgflag='A',&
    long_name='Canopy leaf structural N per m2 leaf area, excluding storage',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_Petole_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='PetolSheth_N_pft',units='gN/m2',avgflag='A',&
    long_name='Canopy sheath N',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_STALK_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='STALK_N_pft',units='gN/m2',avgflag='A',&
    long_name='Canopy stalk N',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_RESERVE_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RESERVE_N_pft',units='gN/m2',avgflag='A',&
    long_name='Canopy reserve N',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_HUSK_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='HUSK_N_pft',units='gN/m2',avgflag='A',&
    long_name='Canopy husk N',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_GRAIN_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='GRAIN_N_pft',units='gN/m2',avgflag='A',&
    long_name='Canopy grain C',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ROOT_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root_N_pft',units='gN/m2',avgflag='A',&
    long_name='Root nitrogen',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_RootNodule_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootNodule_N_pft',units='gN/m2',avgflag='A',&
    long_name='Root total nodule N',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_STORED_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SSTORED_N_pft',units='gN/m2',avgflag='A',&
    long_name='Plant seasonal storage of nonstructural N',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_EXUD_N_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='EXUD_CumYr_N_FLX_pft',units='gN/m2',avgflag='I',&
    long_name='Cumulative Root organic N uptake (<0 exudation into soil)',&
    ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_Uptk_N_Flx_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Uptk_N_CumYr_FLX_pft',units='gN/m2',avgflag='I',&
    long_name='Cumulative Root N uptake (including <0 exudation to soil)',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_Uptk_NMIN_Flx_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Uptk_NMin_CumYr_FLX_pft',units='gN/m2',avgflag='I',&
    long_name='Cumulative Root mineral N uptake',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_Uptk_PMIN_Flx_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Uptk_PMin_CumYr_FLX_pft',units='gP/m2',avgflag='I',&
    long_name='Cumulative Root mineral P uptake',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_Uptk_P_Flx_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Uptk_P_CumYr_FLX_pft',units='gP/m2',avgflag='I',&
    long_name='Cumulative Root P uptake (including <0 exudation to soil)',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_LITRf_N_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LITRf_N_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='total plant LitrFall N',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_cum_N_FIXED_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='cumN_FIXED_pft',units='gN/m2',avgflag='I',&
    long_name='cumulative plant N2 fixation',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_TreeRingRadius_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='TreeRingRadius_pft',units='m',avgflag='I',&
    long_name='Mean main stalk radius',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_HVST_N_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='HVST_N_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='Plant N harvest',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_NH3can_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='NH3can_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='total canopy NH3 flux',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_PLANT_BALANCE_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Plant_BALANCE_N_pft',units='gC/m2',avgflag='A',&
    long_name='Cumulative plant N conservation error',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_STANDING_DEAD_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='STANDING_DEAD_N_pft',units='gN/m2',avgflag='A',&
    long_name='pft standing dead N',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_FIREp_N_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='FIREp_N_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='Plant N emission from fire',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_SURF_LITRf_N_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SURF_LITRf_N_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='total surface LitrFall N',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_SHOOT_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOT_P_pft',units='gP/m2',avgflag='A',&
    long_name='Live plant shoot P',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_Root1stTipSinkWt_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root1stTipSinkwt_pft',units='-',avgflag='A',&
    long_name='Primary root Tip Sink weight',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_Plant_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Plant_P_pft',units='gP/m2',avgflag='A',&
    long_name='Plant P',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_stomatal_stress_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='STOMATAL_STRESS_pft',units='none',avgflag='A',&
    long_name='stomatal stress from root turogr [0->1 increasing stress]',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CANDew_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Canopy_DEW_pft',units='mm H2O/m2',avgflag='I',&
    long_name='Cumulative canopy dew deposition',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_LEAF_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LEAF_P_pft',units='gP/m2',avgflag='A',&
    long_name='Canopy leaf P',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_Petole_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='PetolSheth_P_pft',units='gP/m2',avgflag='A',&
    long_name='Canopy sheath P',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_STALK_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='STALK_P_pft',units='gP/m2',avgflag='A',&
    long_name='Plant stalk P',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_RESERVE_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RESERVE_P_pft',units='gP/m2',avgflag='A',&
    long_name='Plant reserve P',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_HUSK_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='HUSK_P_pft',units='gP/m2',avgflag='A',&
    long_name='Husk P',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_GRAIN_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='GRAIN_P_pft',units='gP/m2',avgflag='A',&
    long_name='Canopy grain P',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_ROOT_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root_P_pft',units='gP/m2',avgflag='A',&
    long_name='Plant root P',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_RootNodule_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootNodule_P_pft',units='gP/m2',avgflag='A',&
    long_name='Root total nodule P',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_STORED_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SSTORED_P_pft',units='gP/m2',avgflag='A',&
    long_name='Plant seasonal storage of nonstructural P',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_EXUD_P_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='EXUD_CumYr_P_FLX_pft',units='gP/m2',avgflag='I',&
    long_name='Cumulative root organic P uptake (<0 exudation into soil)',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_RDECOMPC_SOM_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RDecompC_SOM_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Hydrolysis of SOM C in litter layer',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_MicrobAct_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='MicrobAct_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Respiration-based micoribal activity for hydrolysis in litter layer',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_RDECOMPC_BReSOM_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RDecompC_BReSOM_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Hydrolysis of microbial residual OM C in litter layer',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_RDECOMPC_SorpSOM_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RDecompC_SorpSOM_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Hydrolysis of adsorbed OM C in litter layer',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tRespGrossHeteUlm_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='tRespGrossHeteUlm_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Oxygen unlimited gross heterotrophic respiraiton in litter layer',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tRespGrossHete_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='tRespGrossHete_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Oxygen-limited gross heterotrophic respiraiton in litter layer',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_LITRf_P_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LITRf_P_FLX_pft',units='gP/m2/hr',avgflag='A',&
    long_name='total plant LitrFall P',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_HVST_P_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='HVST_P_FLX_pft',units='gP/m2/hr',avgflag='A',&
    long_name='Plant P harvest',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_PLANT_BALANCE_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Plant_BALANCE_P_pft',units='gP/m2',avgflag='A',&
    long_name='Cumulative plant P conservation error',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_STANDING_DEAD_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='STANDING_DEAD_P_pft',units='gP/m2',avgflag='A',&
    long_name='pft Standing dead P',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_FIREp_P_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='FIREp_P_FLX_pft',units='gP/m2/hr',avgflag='A',&
    long_name='Plant PO4 emission from fire',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_SURF_LITRf_P_FLX_ptc(beg_ptc:end_ptc)         !SurfLitrfallElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
  call hist_addfld1d(fname='SURF_LITRf_P_FLX_pft',units='gP/m2/hr',avgflag='A',&
    long_name='Plant LitrFall P to the soil surface',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_ShootRootXferC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ShootRoot_XFER_C_pft',units='gC/hr/plant',avgflag='A',&
    long_name='Shoot C transfered to root via phloem (>0 to roots)',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_ShootRootXferN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ShootRoot_XFER_N_pft',units='gN/hr/plant',avgflag='A',&
    long_name='Shoot N transfered to root via phloem and xylem (>0 to roots)',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_ShootRootXferP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ShootRoot_XFER_P_pft',units='gP/hr/plant',avgflag='A',&
    long_name='Shoot P transfered to root via phloem and xylem (>0 to roots)',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_BRANCH_NO_ptc(beg_ptc:end_ptc)            !NumOfBranches_pft(NZ,NY,NX)
  call hist_addfld1d(fname='BRANCH_NO_pft',units='none',avgflag='I',&
    long_name='Plant branch number',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_Growth_Stage_ptc(beg_ptc:end_ptc)  !plant development stage, integer, 0-10, planting, emergence, floral_init, jointing,
                                                               !elongation, heading, anthesis, seed_fill, see_no_set, seed_mass_set, end_seed_fill
  call hist_addfld1d(fname='Growth_Stage_pft',units='none',avgflag='I',&
    long_name='Plant development stage, integer, 0-planting, 1-emergence, 2-floral_init, 3-jointing,'// &
    '4-elongation, 5-heading, 6-anthesis, 7-seed_fill, 8-see_no_set, 9-seed_mass_set, 10-end_seed_fill',&
    ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_LEAF_NC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LEAF_rNC_pft',units='gN/gC',avgflag='A',&
    long_name='Mass based plant leaf NC ratio',ptr_patch=data1d_ptr,default='inactive')

    data1d_ptr => this%h1D_MainBranchNO_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='MainBranchNO_pft',units='-',avgflag='A',&
    long_name='Main branch number',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_RCanMaintDef_CO2_pft(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RCanMaintDef_CO2_pft',units='gC m-2 h-1',avgflag='A',&
    long_name='Canopy maintenance respiraiton deficit as CO2 (<0 deficit)',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_RootMaintDef_CO2_pft(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootMaintDef_CO2_pft',units='gC m-2 h-1',avgflag='A',&
    long_name='Root maintenance respiraiton deficit as CO2 (<0 deficit)',ptr_patch=data1d_ptr,default='inactive')

  end procedure register_hist_plants

end submodule HistRegisterPlant
