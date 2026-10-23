submodule (HistDataType) HistRegisterSoilBGC
  use GridConsts, only: JZ, MaxNumRootAxes
  use EcoSIMConfig, only: jcplx => jcplxc
  use EcoSiMParDataMod, only: micpar
  use HistFileMod, only: hist_addfld1d, hist_addfld2d
  implicit none
contains

  module procedure register_hist_soil_bgc
    integer :: beg_col, end_col, beg_ptc, end_ptc
    real(r8), pointer :: data1d_ptr(:)
    real(r8), pointer :: data2d_ptr(:,:)
    integer :: jj

    beg_col=1; end_col=bounds%ncols
    beg_ptc=1; end_ptc=bounds%npfts

  data2d_ptr => this%h2D_CO2_MassL_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CO2_MassL_vr',units='gC/m2',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved CO2 mass in soil, dissolved plus gas, excluding roots',&
    ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_O2_MassL_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='O2_MassL_vr',units='g/m2',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved O2 mass in soil, dissolved plus gas, excluding roots',&
    ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_CH4_MassL_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4_MassL_vr',units='gC/m2',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved CH4 mass in soil, dissolved plus gas, excluding roots',&
    ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_H2_MassL_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='H2_MassL_vr',units='g/m2',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved H2 mass in soil, dissolved plus gas, excluding roots',&
    ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_QDrainloss_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='QDrainloss_vr',units='mm H2O/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved water drainage',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_FermOXYI_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FermOXYI_vr',units='none',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved oxygen inhibitor of fermentation [0->1, weaker inhibition]',&
    ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_litrC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='litrC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Column-level Vertically resolved litter C',ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_litrN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='litrN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved litter N',ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_litrP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='litrP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved litter P',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_tSOC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tSOC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total soil organic C (everything organic)',&
    ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_POM_C_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='POM_C_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved particulate organic C (resulting from the humification process)',&
    ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_MAOM_C_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='MAOM_C_vr',units='gC (kg soil)-1',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved mineral associated organic C (resulting from sorption process)',&
    ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_microbC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tMicrobeC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total live microbial C',&
    ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_microbN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tMicrobeN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total live microbial N',&
    ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_microbP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tMicrobeP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total live microbial P',&
    ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_AeroBact_PrimS_lim_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='AeroBact_PrimS_lim_vr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved primary substrate limitation for aerobic heterotrophic bacteria',&
    ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_AeroFung_PrimS_lim_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='AeroFung_PrimS_lim_vr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved primary substrate limitation for aerobic heterotrophic fungi',&
    ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_tSOCL_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tSOCL_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved total soil organic C (everything organic)',ptr_col=data2d_ptr,&
    default='inactive')

  data2d_ptr => this%h2D_tSON_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tSON_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total soil organic N (everything organic)',ptr_col=data2d_ptr,&
    default='inactive')

  data2d_ptr => this%h2D_BotDEPZ_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='BOTDEPZ_vr',units='m',type2d='levsoi',avgflag='A',&
    long_name='Bottom depth of soil layer',ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_tSOP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tSOP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total soil organic P (everything organic)',&
    ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_NO3_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='NO3_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved dissolved NO3 concentration',&
    ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_NH4_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='NH4_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved dissolved NH4 concentration',&
    ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_VHeatCap_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='VHeatCap_vr',units='MJ/m3/K',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved Volumetric heat capacity',ptr_col=data2d_ptr,default='inactive')

  data1d_ptr => this%h1D_Num_Leaves_ptc(beg_ptc:end_ptc)        !NumOfLeaves_brch(MainBranchNum_pft(NZ,NY,NX),NZ,NY,NX), leaf NO
  call hist_addfld1d(fname='Num_Leaves_pft',units='none',avgflag='A',&
    long_name='Number of leaves',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_RUB_ACTVN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RUB_ACTVN_pft',units='none',avgflag='A',&
    long_name='mean rubisco activity for CO2 fixation across branches, 0-1',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_NumPrimeRootAxes_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='N1stRootAXes_pft',units='none',avgflag='A',&
    long_name='Mean number of seminal root axes of the plant population',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_CanopyNLim_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CanopyNLim_pft',units='none',avgflag='A',&
    long_name='mean canopy nitrogen limitation across branches, 0->1 weaker limitation',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_CanopyPLim_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CanopyPLim_pft',units='none',avgflag='A',&
    long_name='mean canopy phosphorus limitation across branches, 0->1 weaker limitation',ptr_patch=data1d_ptr,default='inactive')

  data2d_ptr => this%h2D_Root_CO2_vr(beg_col:end_col,1:JZ)          !trc_solcl_vr(idg_CO2,1:JZ,NY,NX)
  call hist_addfld2d(fname='Root_CO2_mass_vr',units='gC/m2',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved CO2 mass in roots',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Root_O2_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Root_O2_mass_vr',units='g/m2',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved O2 mass in roots',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Root_CH4_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Root_CH4_mass_vr',units='gC/m2',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved CH4 mass in roots',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Root_H2_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Root_H2_mass_vr',units='g/m2',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved H2 mass in roots',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_O2_rootconduct_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='O2_root_conductance_pvr',units='1/h',type2d='levsoi',avgflag='A',&
    long_name='Root conductance for O2 gaseous in soil layer',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_CO2_rootconduct_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='CO2_root_conductance_pvr',units='1/h',type2d='levsoi',avgflag='A',&
    long_name='Root conductance for CO2 gaseous in soil layer',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_ProteinNperm2LeafArea_pnd(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='ProteinNperm2LeafArea_pnd',units='gN/(m2 leaf area)',type2d='node',avgflag='A',&
    long_name='Areal leaf protein N concentration by node',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Aqua_CO2_vr(beg_col:end_col,1:JZ)          !trc_solcl_vr(idg_CO2,1:JZ,NY,NX)
  call hist_addfld2d(fname='CO2w_conc_vr',units='gC/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous CO2 concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Aqua_CH4_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4w_conc_vr',units='gC/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous CH4 concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Aqua_O2_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='O2w_conc_vr',units='gO/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous O2 concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Aqua_N2O_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='N2Ow_conc_vr',units='gN/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous N2O concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Aqua_NH3_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='NH3w_conc_vr',units='gN/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous NH3 concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Aqua_H2_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='H2w_conc_vr',units='gH/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous H2 concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Aqua_Ar_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Arw_conc_vr',units='gN/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous Ar concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Aqua_N2_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='N2w_conc_vr',units='gN/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous N2 concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_TEMP_vr(beg_col:end_col,1:JZ)         !TCS_vr(1:JZ,NY,NX)
  call hist_addfld2d(fname='Temp_vr',units='oC',type2d='levsoi',avgflag='A',&
    long_name='soil temperature profile',ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_decomp_OStress_vr(beg_col:end_col,1:JZ)         !
  call hist_addfld2d(fname='Decomp_OStress_vr',units='none',type2d='levsoi',avgflag='A',&
    long_name='decomposition oxygen stress [0->1 weaker]',ptr_col=data2d_ptr)

  data2d_ptr => this%h2D_RO2Decomp_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RO2Decomp_flx_vr',units='gO2/m2/h',type2d='levsoi',avgflag='A',&
    long_name='Decomposition O2 uptake in soil layers',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Decomp_temp_FN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Decomp_TEMP_FN_vr',units='none',type2d='levsoi',avgflag='A',&
    long_name='Temeprature dependence of microbial decomposition in soil layers',&
    ptr_col=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_FracLitMix_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FracLitMix_vr',units='none',type2d='levsoi',avgflag='A',&
    long_name='Fraction of litter to mixed with the next layer (>0 downward mixing)',ptr_col=data2d_ptr,&
    default='inactive')

  data2d_ptr => this%h2D_Decomp_Moist_FN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Decomp_Moist_FN_vr',units='none',type2d='levsoi',avgflag='A',&
    long_name='Moisture dependence of microbial decomposition in soil layers',ptr_col=data2d_ptr,&
    default='inactive')
!-----

  data2d_ptr =>  this%h2D_RootMassC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Root C density profile',ptr_col=data2d_ptr)

  data2d_ptr =>  this%h2D_RootMassC_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootC_pvr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Root C density profile of different pft',ptr_patch=data2d_ptr)

  data2d_ptr =>  this%h2D_RootRadialKond2H2O_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='KH2ORadial_pvr',units='m H2O h-1 MPa-1',type2d='levsoi',avgflag='A',&
    long_name='Radial root conductance for water uptake',ptr_patch=data2d_ptr)

  data2d_ptr =>  this%h2d_RootPop_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootPop_pvr',units='#/m2',type2d='levsoi',avgflag='A',&
    long_name='Root population',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_MycoPop_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='MycoPop_pvr',units='#/m2',type2d='levsoi',avgflag='A',&
    long_name='Mycorrhizal population',ptr_patch=data2d_ptr,default='inactive')

  data2d_ptr => this%h2D_Root1stDepz_ptc(beg_ptc:end_ptc,1:MaxNumRootAxes)
  call hist_addfld2d(fname='Root1stDepz_pft',units='m',type2d='rootaxs',avgflag='A',&
    long_name='Primary root tips depth',ptr_patch=data2d_ptr)

  data2d_ptr =>  this%h2D_RootAxialKond2H2O_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='KH2OAxial_pvr',units='m3 H2O h-1 MPa-1',type2d='levsoi',avgflag='A',&
    long_name='Axial root conductance for water uptake',ptr_patch=data2d_ptr)

  data2d_ptr =>  this%h2D_VmaxNH4Root_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='VmaxNH4Root_pvr',units='umolN h-1 (gC root)-1',type2d='levsoi',avgflag='A',&
    long_name='Maximum NH4 uptake rate for given pft',ptr_patch=data2d_ptr)

  data2d_ptr =>  this%h2D_VmaxNO3Root_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='VmaxNO3Root_pvr',units='umolN h-1 (gC root)-1',type2d='levsoi',avgflag='A',&
    long_name='Maximum NO3 uptake rate for given pft',ptr_patch=data2d_ptr)

  data2d_ptr =>  this%h2D_RootMassN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Root N density profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_RootMassP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Root P density profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_DOC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='DOC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='DOC profile',ptr_col=data2d_ptr)

  data2d_ptr =>  this%h2D_DON_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='DON_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='DON profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_DOP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='DOP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='DOP profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_SoilBulkStress_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootPenitStress_vr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='Soil resistance for root penetration',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_acetate_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Acetate_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Acetate profile',ptr_col=data2d_ptr,default='inactive')

  data1d_ptr => this%h1D_tDOC_soil_col(beg_col:end_col)
  call hist_addfld1d(fname='DOC_soil_col',units='gC/m2',avgflag='A',&
    long_name='DOC in soil',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_tDON_soil_col(beg_col:end_col)
  call hist_addfld1d(fname='DON_soil_col',units='gN/m2',avgflag='A',&
    long_name='DON in soil',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RCH4Oxi_aero_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RCH4_AMOX_litr_col',units='gC/m2/h',avgflag='A',&
    long_name='Aerobic CH4 oxidation in litter layer',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RCH4Oxi_anmo_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RCH4_ANMO_litr_col',units='gC/m2/h',avgflag='A',&
    long_name='Anaerobic CH4 oxidation in litter layer',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RFermen_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RFerment_litr_col',units='gC/m2/h',avgflag='A',&
    long_name='Anaerobic C fermentation in litter layer',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_NH3oxi_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='NH3Oxi_litr_col',units='gN/m2/h',avgflag='A',&
    long_name='Nitrifier NH3 to NO2(-) oxidation in litter layer',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_NO2Oxi_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='NO2Oxi_litr_col',units='gN/m2/h',avgflag='A',&
    long_name='Nitrifier NO2(-) to NO3(-) oxidation rate in litter layer',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RDen_NO3toNO2_col(beg_col:end_col)
  call hist_addfld1d(fname='NO3redux_litr_col',units='gN/m2/h',avgflag='A',&
    long_name='Denitrifier NO3(-) to NO2(-) reduction rate in litter layer',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RCH4ProdAcet_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4Prod_Acet_col',units='gC/m2/h',avgflag='A',&
    long_name='Acetoclastic CH4 production integrated over soil and surface litter',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RCH4ProdHG_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4Prod_HG_col',units='gC/m2/h',avgflag='A',&
    long_name='Hydrogenotrophic CH4 production integrated over soil and surface litter',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RCH4Oxi_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4Oxi_col',units='gC/m2/h',avgflag='A',&
    long_name='Aerobic plus anaerobic CH4 oxidation integrated over soil and surface litter',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RCH4Oxi_aero_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4Oxi_aero_col',units='gC/m2/h',avgflag='A',&
    long_name='Aerobic CH4 oxidation integrated over all layers',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RCH4Oxi_anmo_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4Oxi_anmo_col',units='gC/m2/h',avgflag='A',&
    long_name='Anaerobic CH4 oxidation integrated over all layers',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RCH4ProdHydrog_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4Prod_HG_litr_col',units='gC/m2/h',avgflag='A',&
    long_name='Hydrogenotrophic CH4 produciton in litter layer',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_RCH4ProdAcetcl_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4Prod_Acet_litr_col',units='gC/m2/h',avgflag='A',&
    long_name='Acetoclastic CH4 produciton in litter layer',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tDOP_soil_col(beg_col:end_col)
  call hist_addfld1d(fname='DOP_soil_col',units='gP/m2',avgflag='A',&
    long_name='DOP in soil',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_tAcetate_soil_col(beg_col:end_col)
  call hist_addfld1d(fname='Acetate_soil_col',units='gC/m2',avgflag='A',&
    long_name='Acetate in soil',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_DOC_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='DOC_litr_col',units='gC/m3',avgflag='A',&
    long_name='DOC in litter',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_DON_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='DON_litr_col',units='gN/m3',avgflag='A',&
    long_name='DON in litter',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_DOP_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='DOP_litr_col',units='gP/m3',avgflag='A',&
    long_name='DOP in litter',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_acetate_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='Acetate_litr_col',units='gC/m3',avgflag='A',&
    long_name='Acetate in litter',ptr_col=data1d_ptr,default='inactive')

!------
  data2d_ptr =>  this%h2D_AeroHrBactC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrBacterC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic bacteria C profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_cyanoBactC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CynoBacterC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Mixtrophic cyanobacteria C profile',ptr_col=data2d_ptr)

  data2d_ptr =>  this%h2D_AeroHrFungC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrFungiC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic fungi C profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_faculDenitC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Facult_denitrifierC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Facultative denitrifier C biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_fermentorC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FermentorC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Fermentor C biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_fermentor_frac_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Fermentor_frac_vr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Fraction of microbial C in fermentor',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_acetometh_frac_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='AcetoMethgen_frac_vr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Fraction of microbial C in acetoclastic methanogen',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_hydrogMeth_frac_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='HydrogMethgen_frac_vr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Fraction of microbial C in hydrogenotrohpic methanogen',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_acetometgC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Acetic_methanogenC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aceticlastic methanogen C biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_aeroN2fixC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_N2fixerC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic N2 fixer C biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_Gas_Pressure_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='GAS_PRESSURE_vr',units='Pa',type2d='levsoi',avgflag='A',&
    long_name='Soil gas pressure profile',ptr_col=data2d_ptr)

  data2d_ptr =>  this%h2D_CO2_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CO2_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous CO2 profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_CH4_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous CH4 profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_H2_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='H2_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous H2 profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_Ar_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Ar_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous Ar profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_O2_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='O2_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous O2 profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_N2_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='N2_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous N2 profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_N2O_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='N2O_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous N2O profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_NH3_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='NH3_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous NH3 profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_anaeN2FixC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Anaerobic_N2fixerC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Anaerobic N2 fixer C biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_NH3OxiBactC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Ammonia_OxidizerBactC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Ammonia oxidize bacteria C biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_NO2OxiBactC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Nitrie_OxidizerBactC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Nitrite oxidize bacteria C profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_CH4AeroOxiC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_methanotrophC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic methanotroph C biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_H2MethogenC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Hygrogen_methanogenC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Hydrogenotrophic methanogen C biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_TSolidOMActC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='SoildOM_Act_vr',units='gC/m2',type2d='levsoi',avgflag='A',&
    long_name='Active solid organic C in soil layer',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_TSolidOMActCDens_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='SoildOM_Act_Dens_vr',units='gC/gC',type2d='levsoi',avgflag='A',&
    long_name='Active solid organic C Density in soil layer',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_tOMActCDens_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tOMActC_Dens_vr',units='gC/gC',type2d='levsoi',avgflag='A',&
    long_name='Active heterotrophic microbial C Density in soil layer',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_RCH4ProdHydrog_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4Prod_hydro_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved hydrogenotrophic CH4 production rate',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_RCH4ProdAcetcl_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4Prod_aceto_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved acetoclastic CH4 production rate',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_RCH4Oxi_aero_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4Oxi_Aero_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved aerobic CH4 oxidation rate',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_RCH4Oxi_anmo_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4Oxi_ANMO_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved anaerobic CH4 oxidation rate',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_RFerment_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Fermentation_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved fermentation rate',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_RootAR_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootAR_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved root respiration rate',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_PAR_RAD_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='PAR_vr',units='umol photon m-2 s-1',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved PAR in soil',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_RootAR2soil_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootAR2Soil_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved root respiratory CO2 rate to soil',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_RootAR2Root_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootAR2Root_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved root respiratory CO2 rate to roots',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_nh3oxi_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Nit_NH3Oxid_vr',units='gN/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved NH3 to NO2(-) oxidation rate',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_RNit_NO2toNO3_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Nit_NO2Oxid_vr',units='gN/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved NO2(-) to NO3(-) oxidation rate',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_RDen_NO3toNO2_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Den_NO3toNO2_vr',units='gN/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved NO3(-) to NO2(-) reduction rate by denitrifcation',&
    ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_n2oprod_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='N2OProd_vr',units='gN/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total N2O production rate by nitrification and (chemo)denitrifcation',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_Eco_HR_CO2_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='HR_CO2_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved heterotrophic respiration rate',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_Gchem_CO2_prod_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Gchem_CO2_Prod_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved geochemical CO2 production rate',ptr_col=data2d_ptr,default='inactive')
!------
  data2d_ptr =>  this%h2D_AeroHrBactN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrBacterN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic bacteria N profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_AeroHrFungN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrFungiN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic fungi N profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_faculDenitN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Facult_denitrifierN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Facultative denitrifier N biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_fermentorN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FermentorN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Fermentor N biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_acetometgN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Acetic_methanogenN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aceticlastic methanogen N biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_aeroN2fixN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_N2fixerN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic N2 fixer N biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_anaeN2FixN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Anaerobic_N2fixerN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Anaerobic N2 fixer N biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_NH3OxiBactN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Ammonia_OxidizerBactN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Ammonia oxidize bacteria N biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_NO2OxiBactN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Nitrie_OxidizerBactN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Nitrite oxidize bacteria N profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_CH4AeroOxiN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_methanotrophN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic methanotroph N biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_tRespGrossHeterUlm_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tRespGrossHeterUlm_vr',units='gC/m2/h',type2d='levsoi',avgflag='A',&
    long_name='Total oxygen-unlimited gross respiraiton by heterotrophs',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_tRespGrossHeter_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tRespGrossHeter_vr',units='gC/m2/h',type2d='levsoi',avgflag='A',&
    long_name='Total oxygen gross respiraiton by heterotrophs',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_H2MethogenN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Hygrogen_methanogenN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Hydrogenotrophic methanogen N biomass profile',ptr_col=data2d_ptr,default='inactive')

  data2d_ptr =>  this%h2D_MicrobAct_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='MicrobAct_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved respiration-based microbial acitivity for hydrolysis',&
    ptr_col=data2d_ptr,default='inactive')

  DO jj=1,jcplx
    data2d_ptr =>  this%h3D_HydrolCSOMCps_vr(beg_col:end_col,1:JZ,jj)
    call hist_addfld2d(fname='HydrolCSOM_'//trim(micpar%cplxname(jj))//'_cplx_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
      long_name='Layer resolved carbon hydrolysis in '//trim(micpar%cplxname(jj))//' complex',&
      ptr_col=data2d_ptr,default='inactive')

    data2d_ptr =>  this%h3D_SOC_Cps_vr(beg_col:end_col,1:JZ,jj)
    call hist_addfld2d(fname='SOC_'//trim(micpar%cplxname(jj))//'_cplx_vr',units='gC/m2',type2d='levsoi',avgflag='A',&
      long_name='Layer resolved carbon in '//trim(micpar%cplxname(jj))//' complex',&
      ptr_col=data2d_ptr,default='inactive')

    data2d_ptr =>  this%h3D_SOMHydrylScalCps_vr(beg_col:end_col,1:JZ,jj)
    call hist_addfld2d(fname='SOMHydrlScal_'//trim(micpar%cplxname(jj))//'_cplx_vr',units='none',type2d='levsoi',avgflag='A',&
      long_name='Layer resolved SOM hydrolysis scalar in '//trim(micpar%cplxname(jj))//' complex',&
      ptr_col=data2d_ptr,default='inactive')

    data2d_ptr =>  this%h3D_MicrobActCps_vr(beg_col:end_col,1:JZ,jj)
    call hist_addfld2d(fname='MicrobAct_'//trim(micpar%cplxname(jj))//'_cplx_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
      long_name='Layer resolved respiration-based microbial acitivity for hydrolysis in '//trim(micpar%cplxname(jj))//' complex',&
      ptr_col=data2d_ptr,default='inactive')
  ENDDO

  DO jj=1,micpar%FG_guilds_heter(micpar%mid_HeterAerobBacter)

  enddo

  DO JJ=1,micpar%FG_guilds_heter(micpar%mid_Aerob_Fungi)
  ENDDO

  DO JJ=1,micpar%FG_guilds_heter(micpar%mid_Facult_DenitBacter)
  ENDDO

  DO JJ=1,micpar%FG_guilds_heter(micpar%mid_HeterAerobN2Fixer)

  ENDDO

  DO JJ=1,micpar%FG_guilds_heter(micpar%mid_HeterAnaerobN2Fixer)
  ENDDO

  DO JJ=1,micpar%FG_guilds_heter(micpar%mid_fermentor)
  ENDDO

  DO JJ=1,micpar%FG_guilds_heter(micpar%mid_HeterAcetoCH4GenArchea)

  ENDDO

  DO JJ=1,micpar%FG_guilds_autor(micpar%mid_AutoH2GenoCH4GenArchea)

  ENDDO

  DO JJ=1,micpar%FG_guilds_autor(micpar%mid_AutoAmmoniaOxidBacter)

  ENDDO

  DO JJ=1,micpar%FG_guilds_autor(micpar%mid_AutoNitriteOxidBacter)

  ENDDO

  DO JJ=1,micpar%FG_guilds_autor(micpar%mid_AutoAeroCH4OxiBacter)

  ENDDO

  end procedure register_hist_soil_bgc

end submodule HistRegisterSoilBGC
