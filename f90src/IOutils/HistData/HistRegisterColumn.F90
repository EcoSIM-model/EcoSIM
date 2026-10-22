submodule (HistDataType) HistRegisterColumn
  use EcoSIMCtrlMod, only: lmicrobeMLdiag, salt_model
  use HistFileMod, only: hist_addfld1d
  implicit none
contains

  module procedure register_hist_columns
    integer :: beg_col, end_col, beg_ptc, end_ptc
    real(r8), pointer :: data1d_ptr(:)

    beg_col=1; end_col=bounds%ncols
    beg_ptc=1; end_ptc=bounds%npfts

  !-----------------------------------------------------------------------
  ! initialize history fields
  !--------------------------------------------------------------------

  if(lmicrobeMLdiag)then
    data1d_ptr => this%h1d_TEMP30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='TEMP30cm_col',units='K',avgflag='A',&
      long_name='0-30 cm mean soil temperature',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_THETW30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='THETW30cm_col',units='m3 H2O m-3 soil',avgflag='A',&
      long_name='0-30 cm mean volumetric soil moisture',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_O2wConc30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='O2wConc30cm_col',units='gO m-3 water',avgflag='A',&
      long_name='0-30 cm mean dissolved O2 concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_AcetConc30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='AcetConc30cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-30 cm mean dissolved acetate concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_H2wConc30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='H2wConc30cm_col',units='gH m-3 water',avgflag='A',&
      long_name='0-30 cm mean dissolved H2 concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_CH4wConc30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='CH4wConc30cm_col',units='gC m-3 water',avgflag='A',&
      long_name='0-30 cm mean dissolved CH4 concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_DOCConc30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='DOCConc30cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-30 cm mean dissolved organic C concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_AcetMGC30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='AcetMGC30cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-30 cm mean acetoclastic methanogen C biomass concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_H2MGC30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='H2MGC30cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-30 cm mean hydrogenotrohpic methanogen C biomass concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_FermC30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='FermC30cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-30 cm mean fermentor C biomass concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_AeroMOC30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='AeroMOC30cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-30 cm mean aerobic methanotroph C biomass concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_AeroHRFungC30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='AeroHRFungC30cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-30 cm mean aerobic fungi heterotroph C biomass concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_AeroHRBactC30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='AeroHRBactC30cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-30 cm mean aerobic bacterial heterotroph C biomass concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_RAeroCH4Oxi30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='RAeroCH4Oxi30cm_col',units='gC h-1 m-3 soil',avgflag='A',&
      long_name='0-30 cm mean aerobic CH4 oxidation rate',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_RCH4ProdAcet30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='CH4Prod_Acet30cm_col',units='gC h-1 m-3 soil',avgflag='A',&
      long_name='0-30 cm mean acetoclastic CH4 produciton rate',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_RCH4ProdHG30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='CH4Prod_HG30cm_col',units='gC h-1 m-3 soil',avgflag='A',&
      long_name='0-30 cm mean hydrogenotrohpic CH4 produciton rate',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_RFerment30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='Ferment_30cm_col',units='gC h-1 m-3 soil',avgflag='A',&
      long_name='0-30 cm mean fermentation rate',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_RCO2Ht30cm_col(beg_col:end_col)
    call hist_addfld1d(fname='CO2Ht_30cm_col',units='gC h-1 m-3 soil',avgflag='A',&
      long_name='0-30 cm mean CO2 respiration rate',ptr_col=data1d_ptr,default='inactive')
!----------------------------------------------------------------------------------------------------
    data1d_ptr => this%h1d_TEMP60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='TEMP60cm_col',units='K',avgflag='A',&
      long_name='0-60 cm mean soil temperature',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_THETW60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='THETW60cm_col',units='m3 H2O m-3 soil',avgflag='A',&
      long_name='0-60 cm mean volumetric soil moisture',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_O2wConc60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='O2wConc60cm_col',units='gO m-3 water',avgflag='A',&
      long_name='0-60 cm mean dissolved O2 concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_AcetConc60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='AcetConc60cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-60 cm mean dissolved acetate concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_H2wConc60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='H2wConc60cm_col',units='gH m-3 water',avgflag='A',&
      long_name='0-60 cm mean dissolved H2 concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_CH4wConc60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='CH4wConc60cm_col',units='gC m-3 water',avgflag='A',&
      long_name='0-60 cm mean dissolved CH4 concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_DOCConc60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='DOCConc60cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-60 cm mean dissolved organic C concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_AcetMGC60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='AcetMGC60cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-60 cm mean acetoclastic methanogen C biomass concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_H2MGC60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='H2MGC60cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-60 cm mean hydrogenotrohpic methanogen C biomass concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_FermC60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='FermC60cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-60 cm mean fermentor C biomass concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_AeroMOC60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='AeroMOC60cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-60 cm mean aerobic methanotroph C biomass concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_AeroHRFungC60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='AeroHRFungC60cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-60 cm mean aerobic fungi heterotroph C biomass concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_AeroHRBactC60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='AeroHRBactC60cm_col',units='gC m-3 soil',avgflag='A',&
      long_name='0-60 cm mean aerobic bacterial heterotroph C biomass concentration',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_RAeroCH4Oxi60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='RAeroCH4Oxi60cm_col',units='gC h-1 m-3 soil',avgflag='A',&
      long_name='0-60 cm mean aerobic CH4 oxidation rate',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_RCH4ProdAcet60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='CH4Prod_Acet60cm_col',units='gC h-1 m-3 soil',avgflag='A',&
      long_name='0-60 cm mean acetoclastic CH4 produciton rate',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_RCH4ProdHG60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='CH4Prod_HG60cm_col',units='gC h-1 m-3 soil',avgflag='A',&
      long_name='0-60 cm mean hydrogenotrohpic CH4 produciton rate',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_RFerment60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='Ferment_60cm_col',units='gC h-1 m-3 soil',avgflag='A',&
      long_name='0-60 cm mean fermentation rate',ptr_col=data1d_ptr,default='inactive')

    data1d_ptr => this%h1d_RCO2Ht60cm_col(beg_col:end_col)
    call hist_addfld1d(fname='CO2Ht_60cm_col',units='gC h-1 m-3 soil',avgflag='A',&
      long_name='0-60 cm mean CO2 respiration rate',ptr_col=data1d_ptr,default='inactive')
  endif

  data1d_ptr => this%h1D_cumFIRE_CO2_col(beg_col:end_col)
  call hist_addfld1d(fname='cumFIRE_CO2_col',units='gC m-2',avgflag='I',&
    long_name='cumulative CO2 flux from fire (<0 into atmosphere)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_cumFIRE_CH4_col(beg_col:end_col)
  call hist_addfld1d(fname='cumFIRE_CH4_col',units='gC d-2',avgflag='I', &
    long_name='cumulative CH4 flux from fire (<0 into atmosphere)', ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_cNH4_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='cNH4_LITR_col',units='gN NH4/g litter',avgflag='A', &
    long_name='NH4 concentration in litter',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_cNO3_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='cNO3_LITR_col',units='gN NO3/g litter',avgflag='A',&
    long_name='NO3 concentration in litter',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ECO_HVST_C_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_HVST_C_col',units='gC/m2',avgflag='A',&
    long_name='Harvested C',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_ECO_HVST_N_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_HVST_N_col',units='gN/m2',avgflag='A',&
    long_name='Harvested N',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_ECO_HVST_P_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_HVST_P_col',units='gP/m2',avgflag='A',&
    long_name='Harvested P',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_NET_N_MIN_col(beg_col:end_col)
  call hist_addfld1d(fname='NET_N_MIN_col',units='gN/m2',avgflag='I',&
    long_name='Cumulative net microbial NH4 mineralization (<0 immobilization)',&
    ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tLITR_C_col(beg_col:end_col)
  call hist_addfld1d(fname='tLITR_C_col',units='gC/m2',avgflag='A',&
    long_name='Column integrated total (above+belowground) litter C',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_tRAD_col(beg_col:end_col)
  call hist_addfld1d(fname='RADN_col',units='MJ/m2/hr',avgflag='A',&
    long_name='Total incoming solar radiation',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_tLITR_N_col(beg_col:end_col)
  call hist_addfld1d(fname='tLITR_N_col',units='gN/m2',avgflag='A',&
    long_name='Column integrated total (above+belowground) litter N',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RootAR_col(beg_col:end_col)
  call hist_addfld1d(fname='Root_AR_col',units='gC/m2/h',avgflag='A',&
    long_name='Column integrated root autotrophic respiration',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RootCO2Relez_col(beg_col:end_col)
  call hist_addfld1d(fname='Root_CO2Relez_col',units='gC/m2/h',avgflag='A',&
    long_name='Column integrated root CO2 flux',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_tLITR_P_col(beg_col:end_col)
  call hist_addfld1d(fname='tLITR_P_col',units='gP/m2',avgflag='A',&
    long_name='Column integrated total (above+belowground) litter P',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_HUMUS_C_col(beg_col:end_col)
  call hist_addfld1d(fname='HUMUS_C_col',units='gC/m2',avgflag='A',&
    long_name='colum integrated humus C',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_HUMUS_N_col(beg_col:end_col)
  call hist_addfld1d(fname='HUMUS_N_col',units='gN/m2',avgflag='A',&
    long_name='colum integrated humus N',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_HUMUS_P_col(beg_col:end_col)
  call hist_addfld1d(fname='HUMUS_P_col',units='gP/m2',avgflag='A',&
    long_name='colum integrated humus P',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_AMENDED_C_col(beg_col:end_col)
  call hist_addfld1d(fname='AMENDED_C_col',units='gC/m2',avgflag='A',&
    long_name='Column-integrated total organic fertilizer C amendment',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_AMENDED_N_col(beg_col:end_col)
  call hist_addfld1d(fname='AMENDED_N_col',units='gN/m2',avgflag='A',&
    long_name='Column-integrated total organic fertilizer N amendment',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_AMENDED_P_col(beg_col:end_col)
  call hist_addfld1d(fname='AMENDED_P_col',units='gP/m2',avgflag='A',&
    long_name='Column-integrated total organic fertilizer P amendment',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tLITRf_C_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='tLITRf_C_col',units='gC/m2/hr',avgflag='A',&
    long_name='Column-integrated total (above+belowground) litrFall C flux',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tLITRf_N_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='tLITRf_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='Column-integrated total (above+belowground) litrFall N',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tLITRf_P_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='tLITRf_P',units='gP/m2/hr',avgflag='A',&
    long_name='Column-integrated total (above+belowground) litrFall P',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tEXCH_PO4_col(beg_col:end_col)
  call hist_addfld1d(fname='tEXCH_PO4_col',units='gP/m2',avgflag='A',&
    long_name='Column-integrated total bioavailable (exchangeable) mineral P: H2PO4+HPO4',&
    ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SUR_DOC_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SUR_DOC_FLX_col',units='gC/m2/hr',avgflag='A',&
    long_name='Column-integrated surface DOC flux',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SUR_DON_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SUR_DON_FLX_col',units='gN/m2',avgflag='A',&
    long_name='Column-integrated surface DON flux',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SUR_DOP_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SUR_DOP_FLX',units='gP/m2',avgflag='A',&
    long_name='Column-integrated surface DOP flux',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SUB_DOC_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SUB_DOC_FLX_col',units='gC/m2/hr',avgflag='A',&
    long_name='total subsurface DOC flux',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SUB_DON_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SUB_DON_FLX_col',units='gN/m2/hr',avgflag='A',&
    long_name='total subsurface DON flux',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SUB_DOP_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SUB_DOP_FLX_col',units='gP/m2/hr',avgflag='A',&
    long_name='total subsurface DOP flux',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SUR_DIC_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SUR_DIC_FLX_col',units='gC/m2/hr',avgflag='A',&
    long_name='surface carbon export (CO2, CH4 and carbonate species)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SUR_DIN_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SUR_DIN_FLX_col',units='gN/m2',avgflag='I',&
    long_name='annual cumulative surface NH3/NH4/NO2/NO3 export',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SUR_DIP_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SUR_DIP_FLX_col',units='gP/m2',avgflag='I',&
    long_name='annual cumulative surface phosphate export',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SUB_DIC_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SUB_DIC_FLX_col',units='gC/m2/hr',avgflag='A',&
    long_name='subsurface carbon export (CO2, CH4, carbonate and boundary gas transport)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SUB_DIN_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SUB_DIN_FLX_col',units='gN/m2/hr',avgflag='A',&
    long_name='landscape total subsurface DIN flux',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SUB_DIP_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SUB_DIP_FLX_col',units='gP/m2/hr',avgflag='A',&
    long_name='total subsurface DIP flux',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_HeatFlx2Grnd_col(beg_col:end_col)
  call hist_addfld1d(fname='HeatFlx2Grnd_col',units='MJ/m2/hr',avgflag='A',&
    long_name='Heat flux into the ground',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_CumDryDepoOM_col(beg_col:end_col)
  call hist_addfld1d(fname='CumDryDepoOM_col',units='gC/m2',avgflag='I',&
    long_name='Dry deposition C to ground',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_cyanoBactC_col(beg_col:end_col)
  call hist_addfld1d(fname='CynoBacterC_col',units='gC/m2',avgflag='A',&
    long_name='Mixotrophic cyanobacteria C in surface litter and soil column',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RadSW_Grnd_col(beg_col:end_col)
  call hist_addfld1d(fname='RadSW_Grnd_col',units='W/m2',avgflag='A',&
    long_name='Shortwave Radiation onto the ground',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RadPAR_Grnd_col(beg_col:end_col)
  call hist_addfld1d(fname='RadPAR_Grnd_col',units='umol m-2 s-1',avgflag='A',&
    long_name='PAR Radiation onto the ground',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RadPAR2Soil_col(beg_col:end_col)
  call hist_addfld1d(fname='RadPAR2Soil_col',units='umol m-2 s-1',avgflag='A',&
    long_name='PAR Radiation onto exposed soil',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RadPAR2LitR_col(beg_col:end_col)
  call hist_addfld1d(fname='RadPAR2LitR_col',units='umol m-2 s-1',avgflag='A',&
    long_name='PAR Radiation onto litter',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_CanSWRad_col(beg_col:end_col)
  call hist_addfld1d(fname='RadSW_Canopy_col',units='W/m2',avgflag='A',&
    long_name='Shortwave Radiation onto the grid canopy',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_Qinfl2soi_col(beg_col:end_col)
  call hist_addfld1d(fname='Qinfl2soi_col',units='mm H2O/hr',avgflag='A',&
    long_name='Water flux into the ground',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_QTRANSP_col(beg_col:end_col)
  call hist_addfld1d(fname='QTransp_col',units='mm H2O/hr',avgflag='A',&
    long_name='Soil water loss through transpiration (<0 into atmosphere)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_Qdrain_col(beg_col:end_col)
  call hist_addfld1d(fname='Qdrain_col',units='mm H2O/hr',avgflag='A',&
    long_name='Drainage water flux out (>0) of the soil column',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_Ar_mass_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_mass_col',units='g/m2',avgflag='A',&
    long_name='total Ar mass of the soil column, include that in snow and roots',&
    ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_CO2_mass_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_mass_col',units='g/m2',avgflag='A',&
    long_name='total CO2 mass of the soil column, include that in snow and roots',&
    ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_Gchem_CO2_prod_col(beg_col:end_col)
  call hist_addfld1d(fname='Gchem_CO2_prod_col',units='gC/m2',avgflag='A',&
    long_name='Column integrated CO2 production rate from geochemistry',&
    ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_Ar_soilMass_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_soil_mass_col',units='g/m2',avgflag='A',&
    long_name='total Ar mass of the soil column, excluding that in snow',&
    ptr_col=data1d_ptr,default='inactive')

  IF(salt_model)THEN
    data1d_ptr => this%h1D_tSALT_DISCHG_FLX_col(beg_col:end_col)
    call hist_addfld1d(fname='tSALT_DISCHG_FLX_col',units='mol/m2/hr',avgflag='A',&
      long_name='total subsurface ion flux',ptr_col=data1d_ptr)
  endif

  data1d_ptr => this%h1D_tPREC_P_col(beg_col:end_col)
  call hist_addfld1d(fname='tPREC_P_col',units='gP/m2',avgflag='A',&
    long_name='column integrated total soil precipited P',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tSoilOrgC_col(beg_col:end_col)
  call hist_addfld1d(fname='tSoilOrgC_col',units='gC/m2',avgflag='A', &
    long_name='Column-integrated total soil organic C',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tSoilOrgN_col(beg_col:end_col)
  call hist_addfld1d(fname='tSoilOrgN_col',units='gN/m2',avgflag='A', &
    long_name='Column-integrated total soil organic N',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tSoilOrgP_col(beg_col:end_col)
  call hist_addfld1d(fname='tSoilOrgP_col',units='gP/m2',avgflag='A', &
    long_name='Column-integrated total soil organic P',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tMICRO_C_col(beg_col:end_col)
  call hist_addfld1d(fname='tMICROB_C_col',units='gC/m2',avgflag='A', &
    long_name='Column-integrated micriobial C (include surface litter)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_tMICRO_N_col(beg_col:end_col)
  call hist_addfld1d(fname='tMICROB_N_col',units='gN/m2',avgflag='A', &
    long_name='Column-integrated micriobial N (include surface litter)',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tMICRO_P_col(beg_col:end_col)
  call hist_addfld1d(fname='tMICROB_P_col',units='gP/m2',avgflag='A', &
    long_name='Column-integrated micriobial P (include surface litter)',ptr_col=data1d_ptr, &
    default='inactive')

  data1d_ptr => this%h1D_SnowCanopy_col(beg_col:end_col)
  call hist_addfld1d(fname='SWECanopy_col',units='mmH2O/m2',avgflag='A', &
    long_name='Column-integrated canopy held  (water equivalent) snow',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_PO4_FIRE_col(beg_col:end_col)
  call hist_addfld1d(fname='PO4_FIRE_col',units='gP/m2',avgflag='I',&
    long_name='Cumulative PO4 flux from fire',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_cPO4_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='cPO4_LITR_col',units='gP/g litr',avgflag='A',&
    long_name='PO4 concentration in litter',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_cEXCH_P_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='cEXCH_P_LITR_col',units='gP/g litr',avgflag='A',&
    long_name='concentration of exchangeable inorganic P in litterr',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_NET_P_MIN_col(beg_col:end_col)
  call hist_addfld1d(fname='NET_P_MIN_col',units='gP/m2',avgflag='I',&
    long_name='Cumulative net microbial P mineralization (<0 immobilization)',&
    ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_HUM_col(beg_col:end_col)
  call hist_addfld1d(fname='HMAX_AIR_col',units='kPa',avgflag='X',&
    long_name='daily maximum vapor pressure',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_HUM_col(beg_col:end_col)
  call hist_addfld1d(fname='HMIN_AIR_col',units='kPa',avgflag='M',&
    long_name='daily maximum vapor pressure',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_PSI_SURF_col(beg_col:end_col)
  call hist_addfld1d(fname='PSI_LITR_col',units='MPa',avgflag='A',&
    long_name='Litter layer micropore matric water potential',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SURF_ELEV_col(beg_col:end_col)
  call hist_addfld1d(fname='SURF_ELEV_col',units='m',avgflag='A',&
    long_name='Surface elevation, including litter layer',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tNH4X_col(beg_col:end_col)
  call hist_addfld1d(fname='tNH4_col',units='gN/m2',avgflag='A', &
    long_name='Column-integrated NH4+NH3',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_tNO3_col(beg_col:end_col)
  call hist_addfld1d(fname='tNO3_col',units='gN/m2',avgflag='A',&
    long_name='Column integrated NO3+NO2 content',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_TEMP_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='TEMP_LITR_col',units='oC',avgflag='A',&
    long_name='Litter layer temperature',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_TEMP_surf_col(beg_col:end_col)
  call hist_addfld1d(fname='TEMP_SURF_col',units='oC',avgflag='A',&
    long_name='Ground surface temperature',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_TEMP_SNOW_col(beg_col:end_col)
  call hist_addfld1d(fname='TEMP_SNOW_col',units='oC',avgflag='A',&
    long_name='First snow layer temperature',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_FracBySnow_col(beg_col:end_col)
  call hist_addfld1d(fname='Frac_Snow_Ground_col',units='none',avgflag='A',&
    long_name='Fraction of ground covered by snow',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_FracByLitr_col(beg_col:end_col)
  call hist_addfld1d(fname='Frac_Litr_Ground_col',units='none',avgflag='A',&
    long_name='Fraction of ground covered by litter',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_OMC_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='Surf_LitrC_col',units='gC/m2',avgflag='A',&
    long_name='Total surface litter C, including microbes',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_OMN_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='Surf_LitrN_col',units='gN/m2',avgflag='A',&
    long_name='Total surface litter N, including microbes',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_OMP_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='Surf_LitrP_col',units='gP/m2',avgflag='A',&
    long_name='Total surface litter P, including microbes',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ATM_CO2_col(beg_col:end_col)
  call hist_addfld1d(fname='ATM_CO2_col',units='umol/mol',avgflag='A',&
    long_name='Atmospheric CO2 concentration',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_ATM_CH4_col(beg_col:end_col)
  call hist_addfld1d(fname='ATM_CH4_col',units='umol/mol',avgflag='A',&
    long_name='Atmospheric CH4 concentration',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_NBP_col(beg_col:end_col)
  call hist_addfld1d(fname='NBP_col',units='gC/m2',avgflag='I',&
    long_name='Cumulative net biosphere productivity (<0 into atmosphere)',&
    ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ECO_LAI_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_LAI_col',units='m2/m2',avgflag='A',&
    long_name='Ecosystem LAI',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_ECO_SAI_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_SAI_col',units='m2/m2',avgflag='A',&
    long_name='Ecosystem stem Area index for all live branches',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_Eco_GPP_CumYr_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_GPP_col',units='gC/m2',avgflag='I',&
    long_name='cumulative ecosystem GPP',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_ECO_RA_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_RA_col',units='gC/m2',avgflag='I',&
    long_name='cumulative ecosystem autotrophic respiration',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_Eco_NPP_CumYr_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_NPP_col',units='gC/m2',avgflag='I',&
    long_name='cumulative ecosystem NPP',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_Eco_HR_CumYr_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_RH_col',units='gC/m2',avgflag='I',&
    long_name='Cumulative ecosystem heterotrophic respiration (<0 into atmosphere)',&
    ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_Eco_HR_CO2_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_HR_CO2_col',units='gC/m2/hr',avgflag='A',&
    long_name='Ecosystem heterotrophic respiration as CO2 (<0 into atmosphere)',&
    ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_Eco_HR_CO2_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='HR_CO2_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Heterotrophic respiration as CO2 in litter (<0 into atmosphere)',&
    ptr_col=data1d_ptr)

!  data1d_ptr => this%h1D_Eco_HR_CH4_col(beg_col:end_col)
!  call hist_addfld1d(fname='ECO_RH_CH4',units='gC/m2/hr',avgflag='A',&
!    long_name='Ecosystem heterotrophic respiration as CH4',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_tDIC_col(beg_col:end_col)
  call hist_addfld1d(fname='tDIC_col',units='gC/m2',avgflag='A',&
    long_name='column integrated total soil DIC: CO2+CH4',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_tSTANDING_DEAD_C_col(beg_col:end_col)
  call hist_addfld1d(fname='tSTANDING_DEAD_C_col',units='gC/m2',avgflag='A',&
    long_name='total standing dead C',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_tSTANDING_DEAD_N_col(beg_col:end_col)
  call hist_addfld1d(fname='tSTANDING_DEAD_N_col',units='gN/m2',avgflag='A',&
    long_name='total standing dead N',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tSTANDING_DEAD_P_col(beg_col:end_col)
  call hist_addfld1d(fname='tSTANDING_DEAD_P_col',units='gP/m2',avgflag='A',&
    long_name='total standing dead P',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tPRECIP_col(beg_col:end_col)
  call hist_addfld1d(fname='tPRECIP_col',units='mm/m2',avgflag='I',&
    long_name='cumulative precipitation, including irrigation',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_ECO_ET_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_ET_col',units='mm H2O/m2',avgflag='I',&
    long_name='cumulative total evapotranspiration',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_trcg_Ar_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_cumerr_col',units='gAr/m2',avgflag='I',&
    long_name='cumulative mass error for Ar',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_trcg_O2_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='O2_cumerr_col',units='gO/m2',avgflag='I',&
    long_name='cumulative mass error for O2',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_trcg_N2_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='N2_cumerr_col',units='gN/m2',avgflag='I',&
    long_name='cumulative mass error for N2',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_trcg_NH3_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='NH3_cumerr_col',units='gN/m2',avgflag='I',&
    long_name='cumulative mass error for NH3',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_trcg_H2_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='H2_cumerr_col',units='gH/m2',avgflag='I',&
    long_name='cumulative mass error for H2',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_trcg_CO2_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_cumerr_col',units='gC/m2',avgflag='I',&
    long_name='cumulative mass error for CO2',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_trcg_CH4_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4_cumerr_col',units='gC/m2',avgflag='I',&
    long_name='cumulative mass error for CH4',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1d_CAN_NEE_col(beg_col:end_col)
  call hist_addfld1d(fname='CAN_NEE_col',units='umol C/m2/s',avgflag='A',&
    long_name='Canopy net CO2 exchange (<0 into atmosphere)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_ECO_RADSW_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_RADSW_col',units='W/m2',avgflag='A',&
    long_name='Shortwave radiation absorbed by the ecosystem',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_N2O_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='N2O_LITR_col',units='g/m3',avgflag='A',&
    long_name='N2O solute concentration in soil micropores',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_NH3_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='NH3_LITR_col',units='g/m3',avgflag='A',&
    long_name='NH3 solute concentration in soil micropores',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SOL_RADN_col(beg_col:end_col)
  call hist_addfld1d(fname='SOL_RADN_col',units='W/m2',avgflag='A',&
    long_name='Incoming shortwave radiation on the ecosystem',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_AIR_TEMP_col(beg_col:end_col)
  call hist_addfld1d(fname='AIR_TEMP_col',units='oC',avgflag='A',&
    long_name='air temperature',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_FreeNFix_col(beg_col:end_col)
  call hist_addfld1d(fname='FreeNFix_col',units='gN/hr/m2',avgflag='A',&
    long_name='N fixation by free-living microbes',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_PATM_col(beg_col:end_col)
  call hist_addfld1d(fname='PATM_col',units='kPa',avgflag='A',&
    long_name='atmospheric pressure',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_HUM_col(beg_col:end_col)
  call hist_addfld1d(fname='HUM_col',units='kPa',avgflag='A',&
    long_name='vapor pressure',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_WIND_col(beg_col:end_col)
  call hist_addfld1d(fname='WIND_col',units='m/s',avgflag='A',&
    long_name='wind speed',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_PREC_col(beg_col:end_col)
  call hist_addfld1d(fname='PREC_col',units='mm H2O/m2/hr',avgflag='A',&
    long_name='Total precipitation, excluding irrigation',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_Snofall_col(beg_col:end_col)
  call hist_addfld1d(fname='SNOFAL_col',units='mm H2O/m2/hr',avgflag='A',&
    long_name='Precipitation as snowfall',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SOIL_RN_col(beg_col:end_col)
  call hist_addfld1d(fname='SOIL_RN_col',units='W/m2',avgflag='A',&
    long_name='total net radiation at ground surface (incoming short/long wave - outgoing short/long wave at soil/snow/litter)',&
    ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_LWSky_col(beg_col:end_col)
  call hist_addfld1d(fname='LW_Sky_col',units='W/m2',avgflag='A',&
    long_name='Incoming sky long wave radiation',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SOIL_LE_col(beg_col:end_col)
  call hist_addfld1d(fname='SOIL_LE_col',units='W/m2',avgflag='A',&
    long_name='Latent heat flux into ground surface (exclude plant canopy)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SOIL_H_col(beg_col:end_col)
  call hist_addfld1d(fname='SOIL_H_col',units='W/m2',avgflag='A',&
    long_name='Sensible heat flux into ground surface (exclude plant canopy)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SOIL_G_col(beg_col:end_col)
  call hist_addfld1d(fname='SOIL_G_col',units='W/m2',avgflag='A',&
    long_name='total heat flux out of ground surface (>0 into atmosphere)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_ECO_RN_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_Radnet_col',units='W/m2',avgflag='A',&
    long_name='Ecosystem net radiation (>0 into ecosystem, short+sky_long - plant_long-surf_long)',&
    ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ECO_LE_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_LE_col',units='W/m2',avgflag='A',&
    long_name='Ecosystem latent heat flux (>0 into surface)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_Eco_HeatSen_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_HeatS_col',units='W/m2',avgflag='A',&
    long_name='Ecosystem sensible heat flux (>0 into surface)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_ECO_Heat2G_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_Heat2G_col',units='W/m2',avgflag='A',&
    long_name='Heat flux to warm the ecosystem (<0 into atmosphere),' &
    //' including canopy, snow, litter and exposed soil',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_O2_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='O2w_conc_LITR_col',units='g/m3',avgflag='A',&
    long_name='O2 solute concentration in litter layer',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_MIN_LWP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='MIN_LWP_pft',units='MPa',avgflag='A',&
    long_name='minimum daily canopy water potential',ptr_patch=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SLA_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SLA_pft',units='cm2 leaf (gC leaf)-1',avgflag='A',&
    long_name='Specific leaf area',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_CO2_SEMIS_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_SEMIS_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='Surface CO2 flux (< 0 into atmosphere), '// &
    'excluding wet deposition from rainfall and irrigation',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_AR_SEMIS_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_SEMIS_FLX_col',units='umol Ar/m2/s',avgflag='A',&
    long_name='soil Ar flux (< 0 into atmosphere), '// &
    'excluding wet deposition from rainfall and irrigation',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_ECO_CO2_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_NEE_CO2_col',units='umol C/m2/s',avgflag='A',&
    long_name='ecosystem net CO2 exchange (<0 into atmosphere)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_CH4_SEMIS_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4_SEMIS_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='Surface CH4 flux (<0 into atmosphere), '// &
    'excluding wet deposition from rainfall and surface irrigation',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_CH4_EBU_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4_EBU_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='soil CH4 ebullition flux (<0 into atmosphere)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_Ar_EBU_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_EBU_FLX_col',units='umol Ar/m2/s',avgflag='A',&
    long_name='soil Ar ebullition flux (<0 into atmosphere)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_CO2_TPR_err_col(beg_col:end_col)
  call hist_addfld1d(fname='CumCO2_Transpt_Residual_col',units='gC/m2',avgflag='I',&
    long_name='Cumulative difference between soil CO2 production and surface CO2 flux',&
    ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_CO2_Drain_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_DRAINLOSS_col',units='gC/m2/hr',avgflag='A',&
    long_name='CO2 loss flux through subsurface drainage',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_CO2_hydloss_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_Cum_Hyd_Loss_col',units='gC/m2',avgflag='I',&
    long_name='Cumulative hydrological CO2 loss flux, including subsurface drainage',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_NH3_hydloss_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='NH3_Cum_Hyd_Loss_col',units='gN/m2',avgflag='I',&
    long_name='Cumulative hydrological NH3/N4 loss flux, including subsurface drainage',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_NO3_hydloss_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='NO3_Cum_Hyd_Loss_col',units='gN/m2',avgflag='I',&
    long_name='Cumulative hydrological NO3 loss flux, including subsurface drainage',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_Ar_TPR_err_col(beg_col:end_col)
  call hist_addfld1d(fname='CumAr_Transpt_Residual_col',units='g/m2',avgflag='I',&
    long_name='Cumulative difference between soil Ar production and surface Ar flux',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CH4_PLTROOT_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4_PLTROOT_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='soil CH4 flux through plants(<0 into atmosphere)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_AR_PLTROOT_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_PLTROOT_FLX_col',units='umol Ar/m2/s',avgflag='A',&
    long_name='soil AR flux through plants(<0 into atmosphere)',ptr_col=data1d_ptr, &
    default='inactive')

  data1d_ptr => this%h1D_CO2_PLTROOT_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_PLTROOT_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='soil CO2 flux through plants(<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_TKQ_col(beg_col:end_col)
  call hist_addfld1d(fname='TKQ_col',units='K',avgflag='A',&
    long_name='Sink level atmospheric temperature',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_VPQ_col(beg_col:end_col)
  call hist_addfld1d(fname='VPQ_col',units='kPa',avgflag='A',&
    long_name='Sink level atmospheric vapor pressure',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_O2_PLTROOT_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='O2_PLTROOT_FLX_col',units='umol O2/m2/s',avgflag='A',&
    long_name='soil O2 flux through plants(<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CO2_DIF_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_DIF_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='soil CO2 flux through advection+diffusion (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_O2_DIF_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='O2_DIF_FLX_col',units='umol O2/m2/s',avgflag='A',&
    long_name='soil O2 flux through advection+diffusion (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CH4_DIF_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4_DIF_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='soil CH4 flux through advection+diffusion (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_NH3_DIF_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='NH3_DIF_FLX_col',units='umol N/m2/s',avgflag='A',&
    long_name='soil NH3 flux through advection+diffusion (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_Ar_DIF_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_DIF_FLX_col',units='umol Ar/m2/s',avgflag='A',&
    long_name='soil Ar flux through advection+diffusion (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_RoughnessLength_col(beg_col:end_col)
  call hist_addfld1d(fname='Z0_col',units='m',avgflag='A',&
    long_name='Roughness length of the grid (considering vegetation cover)',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_ZeroPlaneDisplacem_col(beg_col:end_col)
  call hist_addfld1d(fname='d_col',units='m',avgflag='A',&
    long_name='Zero plane displacement of the grid (considering vegetation cover)',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_O2_SEMIS_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='O2_SEMIS_FLX_col',units='umol O2/m2/s',avgflag='A',&
    long_name='Surface O2 flux (<0 into atmosphere)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_CO2_LITR_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_LITR_col',units='gC/m3',avgflag='A',&
    long_name='CO2 solute concentration in litter',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_EVAPG_col(beg_col:end_col)
  call hist_addfld1d(fname='EVAPGrnd_col',units='mm H2O/m2/hr',avgflag='A',&
    long_name='Column-integrated ground surface evaporation(>0 into soil)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_CondGasXSurf_col(beg_col:end_col)
  call hist_addfld1d(fname='GasXSurfConduct_col',units='m/hr',avgflag='A',&
    long_name='Conductance for soil-air gas exchange',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_CANET_col(beg_col:end_col)
  call hist_addfld1d(fname='QVegET_col',units='mm H2O/m2/hr',avgflag='A',&
    long_name='Column-integrated canopy evapotranspiration(<0 int atmosphere)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_CanopyEvap_col(beg_col:end_col)
  call hist_addfld1d(fname='QVegEvap_col',units='mm H2O/m2/hr',avgflag='A',&
    long_name='Column-integrated canopy evaporation(<0 int atmosphere)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_tSWC_col(beg_col:end_col)
  call hist_addfld1d(fname='tSWC_col',units='mmH2O/m2',avgflag='A', &
    long_name='column integrated water content (include snow)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_tHeat_col(beg_col:end_col)
  call hist_addfld1d(fname='tSHeat_col',units='MJ/m2',avgflag='A', &
    long_name='column integrated heat content (include snow)',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_SNOWPACK_col(beg_col:end_col)
  call hist_addfld1d(fname='SNOWPACK_col',units='mmH2O/m2',&
    avgflag='A',long_name='total water equivalent snow',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SNOWDENS_col(beg_col:end_col)
  call hist_addfld1d(fname='SNOWDENS_col',units='kg/m3',&
    avgflag='A',long_name='snow density',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SURF_WTR_col(beg_col:end_col)
  call hist_addfld1d(fname='SURF_WTR_col',units='m3/m3',avgflag='A',&
    long_name='Volumetric water content in surface litter layer',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_ThetaW_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='ThetaW_litr_col',units='none',avgflag='A',&
    long_name='Relative saturation of water content in surface litter layer [0-1]',&
    ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_ThetaI_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='ThetaI_litr_col',units='none',avgflag='A',&
    long_name='Relative volume of ice content in surface litter layer [0-1]',&
    ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_SURF_ICE_col(beg_col:end_col)
  call hist_addfld1d(fname='SURF_ICE_col',units='m3/m3',avgflag='A',&
    long_name='Volumetric ice content in surface litter layer',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_ACTV_LYR_col(beg_col:end_col)
  call hist_addfld1d(fname='ACTV_LYR_col',units='m',avgflag='A',&
    long_name='active layer depth',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_WTR_TBL_col(beg_col:end_col)
  call hist_addfld1d(fname='WTBLi_col',units='m',avgflag='A',&
    long_name='internal water table depth (<0 below soil surface)',&
    ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_WTBL_TopoC_col(beg_col:end_col)
  call hist_addfld1d(fname='WTBL_TopoC_col',units='m',avgflag='A',&
    long_name='Topography (slope and elevation) corrected natural water table depth (<0 below soil surface)',&
    ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_N2O_SEMIS_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='N2O_SEMIS_FLX_col',units='gN/m2/hr',&
    avgflag='A',long_name='Surface N2O flux (<0 into atmosphere), '// &
    'including wet deposition from rainfall and irrigation',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_PAR_col(beg_col:end_col)
  call hist_addfld1d(fname='PAR_col',units='umol m-2 s-1',avgflag='A',&
    long_name='Direct plus diffusive incoming photosynthetic photon flux density',&
    ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1d_fPAR_col(beg_col:end_col)
  call hist_addfld1d(fname='fPAR_col',units='-',avgflag='P',&
    long_name='Fraction of absorbed PAR by canopy',&
    ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_CO2_WetDep_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_WetDep_FLX_col',units='gC/m2/hr',&
    avgflag='A',long_name='Wet deposition CO2 flux to soil, '// &
    'from rainfall and irrigation (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_RootN_Fix_col(beg_col:end_col)
  call hist_addfld1d(fname='Root_N_FIX_col',units='gN/m2/hr',&
    avgflag='A',long_name='Root N2 fixation',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_AR_WetDep_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_WetDep_FLX_col',units='gAr/m2/hr',&
    avgflag='A',long_name='Wet deposition Ar flux to soil, '// &
    'from rainfall and irrigation (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_NWetDep_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='Ni_cumWetDep_col',units='gN/m2',&
    avgflag='I',long_name='Cumulative atmospheric wet inorganic N deposition',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_RootXO2_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='RootO2_X_Flx_col',units='gO2/m2/hr',&
    avgflag='A',long_name='O2 consumption rates in roots',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_N2_SEMIS_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='N2_SEMIS_FLX_col',units='gN/m2/hr',&
    avgflag='A',long_name='Surface N2 flux (<0 into atmosphere), '// &
    'including wet deposition from rainfall and irrigation',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_NH3_SEMIS_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='NH3_SEMIS_FLX_col',units='gN/m2/hr',avgflag='A',&
    long_name='Surface NH3 flux (<0 into atmosphere), '// &
    'including wet deposition from rainfall and irrigation',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_H2_SEMIS_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='H2_SEMIS_FLX_col',units='gH/m2/hr',avgflag='A',&
    long_name='Surface H2 flux (<0 into atmosphere), '// &
    'including wet deposition from rainfall and irrigation',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_VHeatCap_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='vHeatCap_litr_col',units='MJ/m3/K',avgflag='A',&
    long_name='surface litter heat capacity',ptr_col=data1d_ptr, &
    default='inactive')

  data1d_ptr => this%h1D_RUNOFF_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='RUNOFF_FLX_col',units='mmH2O/m2/hr',avgflag='A',&
    long_name='Surface runoff from surface water',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_SEDIMENT_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SEDIMENT_FLX_col',units='kg/m2/hr',avgflag='A',&
    long_name='total sediment subsurface flux',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_QDISCHG_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='QDischarge_FLX_col',units='mmH2O/m2/hr',avgflag='A',&
    long_name='grid water lateral discharge with respect external water table (>0 out of grid)',ptr_col=data1d_ptr)

  data1d_ptr => this%h1D_HeatDISCHG_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='HeatDischarge_FLX_col',units='MJ/m2/hr',avgflag='A',&
    long_name='Column-integrated heat flux through discharge',ptr_col=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_LEAF_PC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LEAF_rPC_pft',units='gP/gC',avgflag='I',&
    long_name='Mass based leaf PC ratio',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CAN_RN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_RN_pft',units='W/m2',avgflag='A',&
    long_name='Canopy net radiation',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CAN_LE_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_LE_pft',units='W/m2',avgflag='A',&
    long_name='Canopy latent heat flux (<0 to ATM)',ptr_patch=data1d_ptr,&
    default='inactive')

  data1d_ptr => this%h1D_CAN_H_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_H_pft',units='W/m2',avgflag='A',&
    long_name='Canopy sensible heat flux',ptr_patch=data1d_ptr,&
    default='inactive')

  end procedure register_hist_columns

end submodule HistRegisterColumn
