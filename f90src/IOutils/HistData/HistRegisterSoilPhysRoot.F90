submodule (HistDataType) HistRegisterSoilPhysRoot
  ! Register soil physical and root diagnostics with history metadata and pointers.
  use GridConsts, only: JZ, MaxNumBranches, NumCanopyLayers, NumOfPlantMorphUnits
  use ElmIDMod, only: NumPlantChemElms
  use HistFileMod, only: hist_addfld1d, hist_addfld2d
  implicit none
contains

  module procedure register_hist_soil_phys_root
    ! Expose supported soil and root history fields to the output system.
    integer :: beg_col, end_col, beg_ptc, end_ptc
    real(r8), pointer :: data1d_ptr(:)
    real(r8), pointer :: data2d_ptr(:,:)
    integer :: nbr
    character(len=32) :: fieldname

    beg_col=1; end_col=bounds%ncols
    beg_ptc=1; end_ptc=bounds%npfts

  data2d_ptr =>  this%h2D_RDECOMPC_SOM_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RDecompC_SOM_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved Hydrolysis of solid OM C',ptr_col=data2d_ptr)

  data2d_ptr =>  this%h2D_RDECOMPC_BReSOM_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RDecompC_BReSOM_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved Hydrolysis of microbial residual OM C',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_RDECOMPC_SorpSOM_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RDecompC_SorpSOM_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved Hydrolysis of adsorbed OM C',ptr_col=data2d_ptr,default='inactive')

!------
  data2d_ptr =>  this%h2D_AeroHrBactP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrBacterP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic P biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_AeroHrFungP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrFungiP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic fungi P profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_faculDenitP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Facult_denitrifierP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Facultative denitrifier P biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_fermentorP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FermentorP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Fermentor P biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_acetometgP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Acetic_methanogenP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aceticlastic methanogen P biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_aeroN2fixP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_N2fixerP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic N2 fixer P biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_anaeN2FixP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Anaerobic_N2fixerP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Anaerobic N2 fixer P biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_NH3OxiBactP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Ammonia_OxidizerBactP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Ammonia oxidize bacteria P biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_NO2OxiBactP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Nitrie_OxidizerBactP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Nitrite oxidize bacteria P profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_CH4AeroOxiP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_methanotrophP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic methanotroph P biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_H2MethogenP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Hygrogen_methanogenP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Hydrogenotrophic methanogen P biomass profile',ptr_col=data2d_ptr,default='inactive')
!---

  data2d_ptr =>  this%h2D_MicroBiomeE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='MicroBiomE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Total micorobial elemental biomass in litter',ptr_col=data2d_ptr)

  data2d_ptr =>  this%h2D_AeroHrBactE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Aerobic_HetrBacterE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Aerobic heterotrophic bacterial elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_AeroHrFungE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Aerobic_HetrFungiE_litr',units='gP/m2',type2d='elements',avgflag='A',&
    long_name='Aerobic heterotrophic fungi elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_faculDenitE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Facult_denitrifierE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Facultative denitrifier elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_fermentorE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='FermentorE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Fermentor elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_cyanoBactC_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='CynoBacterE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Cynobacterial elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_acetometgE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Acetic_methanogenE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Aceticlastic methanogen elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_aeroN2fixE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Aerobic_N2fixerE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Aerobic N2 fixer elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_anaeN2FixE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Anaerobic_N2fixerE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Anaerobic N2 fixer elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_NH3OxiBactE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Ammonia_OxidizerBactE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Ammonia oxidize bacteria elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_NO2OxiBactE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Nitrie_OxidizerBactE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Nitrite oxidize bacteria elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_CH4AeroOxiE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Aerobic_methanotrophE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Aerobic methanotroph elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_H2MethogenE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Hygrogen_methanogenE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Hydrogenotrophic methanogen elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_HeatFlow_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='HeatFlow_vr',units='MJ m-3 hr-1',type2d='levsoi',avgflag='A',&
    long_name='soil heat flow profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_HeatUptk_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='HeatUptk_vr',units='MJ m-3 hr-1',type2d='levsoi',avgflag='A',&
    long_name='soil heat flow by plant water uptake (<0 into roots)',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_VSPore_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='VSPore_vr',units='m3 pore/m3 soil',type2d='levsoi',avgflag='A',&
    long_name='Volumetric soil porosity',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_rVSM_vr(beg_col:end_col,1:JZ)        !ThetaH2OZ_vr(1:JZ,NY,NX)
  call hist_addfld2d(fname='rWatFLP_vr',units='m3 H2O/m3 soil pore',type2d='levsoi',avgflag='A',&
    long_name='Fraction of soil porosity filled by water (relative saturation)',ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_FLO_MICP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='MicPFlo_vr',units='mm H2O h-1',type2d='levsoi',avgflag='A',&
    long_name='Micropore water flow (>0) into soil layer',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_FLO_MACP_vr(beg_col:end_col,1:JZ)        !ThetaH2OZ_vr(1:JZ,NY,NX)
  call hist_addfld2d(fname='MacPFlo_vr',units='mm H2O h-1',type2d='levsoi',avgflag='A',&
    long_name='Macropore water flow (>0) into soil',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_rVSICE_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='rIceFLP_vr',units='m3 ice/m3 soil pore',type2d='levsoi',avgflag='A',&
    long_name='fraction of soil porosity filled by ice',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_PSI_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='PSI_vr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='soil matric pressure+osmotic pressure',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_PsiO_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='PsiO_vr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='soil osmotic pressure',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootH2OUP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootH2OUptake_vr',units='mmH2O/hr',type2d='levsoi',avgflag='A',&
    long_name='soil water taken up by root (<0 into roots)',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_cNH4t_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='cNH4t_vr',units='gN/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='soil NH4x concentration',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_cNO3t_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='cNO3t_vr',units='gN/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='Soil NO3+NO2 concentration',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_cPO4_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='cPO4_vr',units='gP/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='soil dissolved PO4 concentration',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_cEXCH_P_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='cEXCH_P_vr',units='gP/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='total exchangeable soil PO4 concentration',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_microb_N2fix_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Free_N2Fix_vr',units='ugN/m3 h-1',type2d='levsoi',avgflag='A',&
    long_name='Free dizotrophic N2 fixation in soil',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_TEMP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='TMAX_SOIL_vr',units='oC',type2d='levsoi',avgflag='X',&
    long_name='Soil maximum temperature profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_TEMP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='TMIN_SOIL_vr',units='oC',type2d='levsoi',avgflag='M',&
    long_name='Soil minimum temperature profile',ptr_col=data2d_ptr,default='inactive')

  data1d_ptr => this%h1D_decomp_ostress_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='Decomp_OStress_LITR',units='none',avgflag='A',&
    long_name='Decomposition O2 stress in litter layer [0->1: weaker]',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_Decomp_temp_FN_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='Decomp_TEMP_FN_LITR',units='none',avgflag='A',&
    long_name='Decomposition temperature sensitivity in litter layer',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_FracLitMix_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='FracLitMix_LITR',units='none',avgflag='A',&
    long_name='Fraction of surface litter layer to be mixed downward',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_Decomp_Moist_FN_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='Decomp_Moist_FN_LITR',units='none',avgflag='A',&
    long_name='Decomposition moisture sensitivity in litter layer',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_RO2Decomp_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RO2Decomp_flx_LITR',units='gO2/m2/h',avgflag='A',&
    long_name='Decomposition O2 uptake in litter layer ',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_TSolidOMActC_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='SoildOM_Act_litr',units='gC/m2',avgflag='M',&
    long_name='Active solid OM in litter',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tOMActCDens_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='tOMActC_Dens_litr',units='gC/gC',avgflag='M',&
    long_name='Active heterotrophic microbial C density in litter',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_TSolidOMActCDens_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='SoildOM_Act_Dens_litr',units='gC/gC',avgflag='M',&
    long_name='Active solid OM density in litter',ptr_col=data1d_ptr,default='inactive')

  data2d_ptr => this%h2D_ElectricConductivity_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='ElectricConductivity_vr',units='dS m-1',type2d='levsoi',avgflag='A',&
    long_name='electrical conductivity',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_HydCondSoil_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='HydCondSoil_vr',units='m MPa-1 h-1',type2d='levsoi',avgflag='A',&
    long_name='Vertical hydraulic conductivity',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootShootExchC_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootShootCX_pvr',units='gC m-2 h-1',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved pft root shoot C exchange (>0 to root)',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootShootExchN_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootShootNX_pvr',units='gN m-2 h-1',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved pft root shoot N exchange (>0 to root)',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootShootExchP_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootShootPX_pvr',units='gP m-2 h-1',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved pft root shoot P exchange (>0 to root)',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_PSI_RT_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='PSI_RT_pvr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='Root total water potential of each pft',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootH2OUptkStress_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootH2OUptkStress_pvr',units='m3 m-2 h-1',type2d='levsoi',avgflag='A',&
    long_name='Rate indicated root water uptake stress of each pft (>0 hydraulic stress)',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootH2OUptk_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootH2OUptk_pvr',units='mm H2O h-1',type2d='levsoi',avgflag='A',&
    long_name='Plant root water uptake from soil (>0 release water to soil)',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_RootAct1stC_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootAct1stC_pvr',units='gC m-3',type2d='levsoi',avgflag='A',&
    long_name='Active zone C in primary roots',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_RootMedC_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootMedC_pvr',units='gC m-3',type2d='levsoi',avgflag='A',&
    long_name='Medium size root biomass C',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_RootLig1stC_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootLig1stC_pvr',units='gC m-3',type2d='levsoi',avgflag='A',&
    long_name='Lignified zone C in primary roots',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_NonstC_conc_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='NonstC_conc_pvr',units='gC nonst gC struct-1',type2d='levsoi',avgflag='A',&
    long_name='Nonstructural C concentration',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_SapFlowVlinear_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='SapFlowVlinear_pvr',units='m h-1',type2d='levsoi',avgflag='A',&
    long_name='Lumen area normalized mean linear sap flow velocity along the vessels of coarse roots',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_RootMaintDef_CO2_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootMaintDef_CO2_pvr',units='g CO2 m-2 h-1',type2d='levsoi',avgflag='A',&
    long_name='Plant root maintenance deficit each pft (<0 deficit)',ptr_patch=data2d_ptr,default='inactive')

  data1d_ptr => this%h1D_RootAbsorbAreaPP_pft(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootAbsorbAreaPP_pft',units='m2/plant',avgflag='M',&
    long_name='Root surface area per plant for nutrient and water absorption',ptr_col=data1d_ptr)

  data2d_ptr => this%h2D_RootAbsorbAreaPP_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootSurfAreaPP_pvr',units='m2 surface per plant',type2d='levsoi',avgflag='A',&
    long_name='Root surface area per plant (for nutrient uptake)',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_ROOTNLim_rpvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootNlim_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Plant root nitrogen limitation for each pft',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_ROOTPLim_rpvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootPlim_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Plant root phosphorus limitation for each pft',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootNonstC_rpvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootNonstC_pvr',units='gC',type2d='levsoi',avgflag='A',&
    long_name='Plant root nonstrucal C for each pft',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_RootSinkWeight_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootSinkWeight_pvr',units='d-2',type2d='levsoi',avgflag='A',&
    long_name='Root nonstructural allocation weight profile for each pft',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_RootMSinkWeight_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootMedSinkWeight_pvr',units='d-2',type2d='levsoi',avgflag='A',&
    long_name='Medium size root nonstructural allocation weight profile for each pft',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_Root2ndSinkWeight_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='Root2ndSinkWeight_pvr',units='d-2',type2d='levsoi',avgflag='A',&
    long_name='Root nonstructural allocation weight for fine roots of each pft',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Root1stSinkWeight_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='Root1stSinkWeight_pvr',units='d-2',type2d='levsoi',avgflag='A',&
    long_name='Root nonstructural allocation weight for primary roots of each pft (excluding tip)',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Root1stRadius_rpvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='Root1stRadius_pvr',units='mm',type2d='levsoi',avgflag='A',&
    long_name='Plant corase root radius for each pft',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_ROOT_OSTRESS_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='Root_OXYSTRESS_pvr',units='None',type2d='levsoi',avgflag='A',&
    long_name='Root Oxygen stress profile [0->1 weaker stress]',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_Root1stStrutC_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootC_1st_pvr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Primary root structural biomass C density',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_MycoBiomC_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='MycoBiomC_pvr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Mycorrhizal biomass C density',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Cytokinin1stConc_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='Cytokinin1stConc_pvr',units='1.e-3gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Primary root cytokinin mean concentration',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_CRootLumenArea_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='CRootLumenArea_pvr',units='m2',type2d='levsoi',avgflag='A',&
    long_name='Mean lumen area for primary root axis',ptr_patch=data2d_ptr,default='inactive')


  data2d_ptr => this%h2D_Cytok_scalar_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='Cytok_scalar_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Cytokinin scalar for corase root thickening',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootNonstBConc_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootNonstBConc_pvr',units='g/gC',type2d='levsoi',avgflag='A',&
    long_name='Plant-root nonstructural C+N+P concentration per unit active root carbon',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_Root1stStrutN_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootN_1st_pvr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Primary root structural biomass N density',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_Root1stStrutP_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootP_1st_pvr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Primary root structural biomass P density',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Root2ndStrutC_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootC_2nd_pvr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Secondary root structural biomass C density',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_Root2ndStrutN_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootN_2nd_pvr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Secondary root structural biomass N density',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Root2ndStrutP_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootP_2nd_pvr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Secondary root structural biomass P density',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Root1stAxesNumL_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='Root1st_AxesNumL_pvr',units='# plant-1',type2d='levsoi',avgflag='A',&
    long_name='Primary root axes number in soil layer',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_Rootmedlength_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootMed_length_pvr',units='m root m-2',type2d='levsoi',avgflag='A',&
    long_name='Total length of medium size root axes in layer',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootmedRadius_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootMed_radius_pvr',units='mm',type2d='levsoi',avgflag='A',&
    long_name='Mean medium size root radius',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootMedAxesNumL_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootMed_AxesNumL_pvr',units='# plant-1',type2d='levsoi',avgflag='A',&
    long_name='Medium size root axes number in soil layer',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Root1stLenPP_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='Root1stLenPP_pvr',units='m (plant)-1',type2d='levsoi',avgflag='A',&
    long_name='Mean primary root axes length in soil layer',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Root2ndAxesNumL_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='Root2nd_AxesNumL_pvr',units='1/d2',type2d='levsoi',avgflag='A',&
    long_name='Secondary root axes number in soil layer',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_Root2ndMeanLen_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='Root2ndMeanLen_pvr',units='m',type2d='levsoi',avgflag='A',&
    long_name='Mean secondary plant-root axis length in soil layer',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_RootKond2H2O_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootKond2H2O_pvr',units='x1.e7 m s-1 MPa-1',type2d='levsoi',avgflag='A',&
    long_name='Total root conductance to water uptake in soil layer',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_prtUP_NH4_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='prtUP_NH4_pvr',units='gN/m3/hr',type2d='levsoi',avgflag='A',&
    long_name='Root uptake of NH4',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_prtUP_NO3_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='prtUP_NO3_pvr',units='gN/m3/hr',type2d='levsoi',&
    avgflag='A',long_name='Root uptake of NO3',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_prtUP_PO4_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='prtUP_PO4_pvr',units='gP/m3/hr',type2d='levsoi',avgflag='A',&
    long_name='Root uptake of PO4',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_DNS_RT_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootLDS_pvr',units='cm/cm3',type2d='levsoi',avgflag='A',&
    long_name='Root length density (including root hair)',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootNutupk_fClim_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootNutUptk_fClim_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='C-availability for root nutrient uptake limitation, 0->1 stronger limitation',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootNutupk_fNlim_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootNutUptk_fNlim_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='N-limitation for root nutrient uptake limitation, 0->1 stronger limitation',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootNutupk_fPlim_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootNutUptk_fPlim_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='P-limitation for root nutrient uptake limitation, 0->1 stronger limitation',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_RootNutupk_fProtC_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootNutUptk_fProtC_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Transporter-limitation for root nutrient uptake capacity, 0->1 weaker limitation',&
    ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_Root1stSArea4GasTP_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='Root1stSA4GasTP_pvr',units='m2',type2d='levsoi',avgflag='A',&
    long_name='Primary root surface area for gas transport',ptr_patch=data2d_ptr,default='inactive')


  data2d_ptr => this%h2D_RootProteinC_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootProteinC_pvr',units='gC m-3',type2d='levsoi',avgflag='A',&
    long_name='Root proteinC concentration',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_fTRootGro_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootGRO_TEMP_FN_pvr',units='none',type2d='levsoi',avgflag='A',&
    long_name='Root growth temperature dependence function',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_CanopyLAIZ_plyr(beg_ptc:end_ptc,1:NumCanopyLayers)
  call hist_addfld2d(fname='CanopyLAIZ_plyr',units='m2/m2',type2d='levcan',avgflag='A',&
    long_name='Vertically distributed leaf area',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_fRootGrowPSISense_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootGRO_PSI_FN_pvr',units='none',type2d='levsoi',avgflag='A',&
    long_name='Root growth moisture dependence function',ptr_patch=data2d_ptr,default='inactive')
  ![terminate]
  do nbr=1,MaxNumBranches
    data2d_ptr => this%h3D_PARTS_ptc(beg_ptc:end_ptc,1:NumOfPlantMorphUnits,nbr)
    write(fieldname,'(I2.2)')nbr
    call hist_addfld2d(fname='C_PARTS_brch_'//trim(fieldname),units='none',&
      type2d='pmorphunits',avgflag='A',&
      long_name='C allocation to different morph unit in branch '//trim(fieldname),ptr_patch=data2d_ptr,default='inactive')
  enddo
  end procedure register_hist_soil_phys_root

end submodule HistRegisterSoilPhysRoot
