module MLDataDiagType

  use data_kind_mod, only: r8 => DAT_KIND_R8
  use EcoSIMConfig,  only: jcplx => jcplxc, nlbiomcp=>NumLiveMicrbCompts
  use ElmIDMod,      only: NumPlantChemElms
  use abortutils,    only: destroy
  use GridConsts
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__


  real(r8), allocatable :: TEMP30cm_col(:,:)         !0-30 cm mean temperature, [K]
  real(r8), allocatable :: THETW30cm_col(:,:)        !0-30 cm mean volumetric water content,  [m3 H2O m-3 soil]
  real(r8), allocatable :: O2wConc30cm_col(:,:)      !0-30 cm mean dissolved O2 concentration, [gO/m3 water] 
  real(r8), allocatable :: CO2wConc30cm_col(:,:)     !0-30 cm mean dissolved CO2 concentration, [gC/m3 water]   
  real(r8), allocatable :: AcetConc30cm_col(:,:)     !0-30 cm mean dissolved acetate concentration, [gC/m3]
  real(r8), allocatable :: H2wConc30cm_col(:,:)      !0-30 cm mean dissolved H2 concentration, [gH/m3]
  real(r8), allocatable :: CH4wConc30cm_col(:,:)     !0-30 cm mean dissolved CH4 concentration,[gC/m3]
  real(r8), allocatable :: DOCConc30cm_col(:,:)      !0-30 cm mean DOC concentration, [gC/m3]
  real(r8), allocatable :: AcetMGC30cm_col(:,:)      !0-30 cm mean acetoclastic methanogen C biomass concentration, [gC/m3]
  real(r8), allocatable :: H2MGC30cm_col(:,:)        !0-30 cm mean hydrogenotrohpic methanogen C biomass concentration, [gC/m3]
  real(r8), allocatable :: FermC30cm_col(:,:)        !0-30 cm mean fermentor C biomass concentration, [gC/m3]
  real(r8), allocatable :: AeroMOC30cm_col(:,:)      !0-30 cm mean aerobic CH4 methanotroph C biomass concentration, [gC/m3]
  real(r8), allocatable :: AeroHRFungC30cm_col(:,:)  !0-30 cm mean aerobic fungi heterotroph C biomass concentration, [gC/m3]
  real(r8), allocatable :: AeroHRBactC30cm_col(:,:)  !0-30 cm mean aerobic bacterial heterotroph C biomass concentration, [gC/m3]
  real(r8), allocatable :: RAeroCH4Oxi30cm_col(:,:)  !0-30 cm mean aerobic CH4 oxidation rate, [gC/h/m3]
  real(r8), allocatable :: RCH4ProdHG30cm_col(:,:)   !0-30 cm mean hydrogenotrohpic CH4 produciton rate, [gC/h/m3]
  real(r8), allocatable :: RCH4ProdAcet30cm_col(:,:) !0-30 cm mean acetoclastic CH4 produciton rate, [gC/h/m3]
  real(r8), allocatable :: RFerment30cm_col(:,:)     !0-30 cm mean fermentation rate, [gC/h/m3]
  real(r8), allocatable :: RCO2Ht30cm_col(:,:)       !0-30 cm mean CO2 respiration rate, [gC/h/m3]

  real(r8), allocatable :: TEMP60cm_col(:,:)         !0-60 cm mean temperature, [K]
  real(r8), allocatable :: THETW60cm_col(:,:)        !0-60 cm mean volumetric water content,  [m3 H2O m-3 soil]
  real(r8), allocatable :: O2wConc60cm_col(:,:)      !0-60 cm mean dissolved O2 concentration, [gO/m3 water] 
  real(r8), allocatable :: CO2wConc60cm_col(:,:)     !0-30 cm mean dissolved CO2 concentration, [gC/m3 water]     
  real(r8), allocatable :: AcetConc60cm_col(:,:)     !0-60 cm mean dissolved acetate concentration, [gC/m3]
  real(r8), allocatable :: H2wConc60cm_col(:,:)      !0-60 cm mean dissolved H2 concentration, [gH/m3]
  real(r8), allocatable :: CH4wConc60cm_col(:,:)     !0-60 cm mean dissolved CH4 concentration,[gC/m3]
  real(r8), allocatable :: DOCConc60cm_col(:,:)      !0-60 cm mean DOC concentration, [gC/m3]
  real(r8), allocatable :: AcetMGC60cm_col(:,:)      !0-60 cm mean acetoclastic methanogen C biomass concentration, [gC/m3]
  real(r8), allocatable :: H2MGC60cm_col(:,:)        !0-60 cm mean hydrogenotrohpic methanogen C biomass concentration, [gC/m3]
  real(r8), allocatable :: FermC60cm_col(:,:)        !0-60 cm mean fermentor C biomass concentration, [gC/m3]
  real(r8), allocatable :: AeroMOC60cm_col(:,:)      !0-60 cm mean aerobic CH4 methanotroph C biomass concentration, [gC/m3]
  real(r8), allocatable :: AeroHRFungC60cm_col(:,:)  !0-60 cm mean aerobic fungi heterotroph C biomass concentration, [gC/m3]
  real(r8), allocatable :: AeroHRBactC60cm_col(:,:)  !0-60 cm mean aerobic bacterial heterotroph C biomass concentration, [gC/m3]
  real(r8), allocatable :: RAeroCH4Oxi60cm_col(:,:)  !0-60 cm mean aerobic CH4 oxidation rate, [gC/h/m3]
  real(r8), allocatable :: RCH4ProdHG60cm_col(:,:)   !0-60 cm mean hydrogenotrohpic CH4 produciton rate, [gC/h/m3]
  real(r8), allocatable :: RCH4ProdAcet60cm_col(:,:) !0-60 cm mean acetoclastic CH4 produciton rate, [gC/h/m3]
  real(r8), allocatable :: RFerment60cm_col(:,:)     !0-60 cm mean fermentation rate, [gC/h/m3]
  real(r8), allocatable :: RCO2Ht60cm_col(:,:)       !0-60 cm mean CO2 respiration rate, [gC/h/m3]

  public :: InitMLDataType
  public :: DestructMLDataType
  contains

!----------------------------------------------------------------------
  subroutine InitMLDataType()
  implicit none

  allocate(TEMP30cm_col(JY,JX)); TEMP30cm_col=0._r8
  allocate(THETW30cm_col(JY,JX)); THETW30cm_col=0._r8
  allocate(O2wConc30cm_col(JY,JX)); O2wConc30cm_col=0._r8
  allocate(CO2wConc30cm_col(JY,JX)); CO2wConc30cm_col=0._r8  
  allocate(AcetConc30cm_col(JY,JX)); AcetConc30cm_col=0._r8
  allocate(H2wConc30cm_col(JY,JX)); H2wConc30cm_col=0._r8
  allocate(CH4wConc30cm_col(JY,JX)); CH4wConc30cm_col=0._r8
  allocate(DOCConc30cm_col(JY,JX)); DOCConc30cm_col=0._r8
  allocate(AcetMGC30cm_col(JY,JX)); AcetMGC30cm_col=0._r8
  allocate(H2MGC30cm_col(JY,JX)); H2MGC30cm_col=0._r8
  allocate(FermC30cm_col(JY,JX)); FermC30cm_col=0._r8
  allocate(AeroMOC30cm_col(JY,JX)); AeroMOC30cm_col=0._r8
  allocate(AeroHRFungC30cm_col(JY,JX)); AeroHRFungC30cm_col=0._r8
  allocate(AeroHRBactC30cm_col(JY,JX)); AeroHRBactC30cm_col=0._r8
  allocate(RAeroCH4Oxi30cm_col(JY,JX)); RAeroCH4Oxi30cm_col=0._r8
  allocate(RCH4ProdAcet30cm_col(JY,JX)); RCH4ProdAcet30cm_col=0._r8
  allocate(RCH4ProdHG30cm_col(JY,JX)); RCH4ProdHG30cm_col=0._r8
  allocate(RFerment30cm_col(JY,JX)); RFerment30cm_col=0._r8
  allocate(RCO2Ht30cm_col(JY,JX)); RCO2Ht30cm_col=0._r8

  allocate(TEMP60cm_col(JY,JX)); TEMP60cm_col=0._r8
  allocate(THETW60cm_col(JY,JX)); THETW60cm_col=0._r8
  allocate(O2wConc60cm_col(JY,JX)); O2wConc60cm_col=0._r8
  allocate(CO2wConc60cm_col(JY,JX)); CO2wConc60cm_col=0._r8  
  allocate(AcetConc60cm_col(JY,JX)); AcetConc60cm_col=0._r8
  allocate(H2wConc60cm_col(JY,JX)); H2wConc60cm_col=0._r8
  allocate(CH4wConc60cm_col(JY,JX)); CH4wConc60cm_col=0._r8
  allocate(DOCConc60cm_col(JY,JX)); DOCConc60cm_col=0._r8
  allocate(AcetMGC60cm_col(JY,JX)); AcetMGC60cm_col=0._r8
  allocate(H2MGC60cm_col(JY,JX)); H2MGC60cm_col=0._r8
  allocate(FermC60cm_col(JY,JX)); FermC60cm_col=0._r8
  allocate(AeroMOC60cm_col(JY,JX)); AeroMOC60cm_col=0._r8
  allocate(AeroHRFungC60cm_col(JY,JX)); AeroHRFungC60cm_col=0._r8
  allocate(AeroHRBactC60cm_col(JY,JX)); AeroHRBactC60cm_col=0._r8
  allocate(RAeroCH4Oxi60cm_col(JY,JX)); RAeroCH4Oxi60cm_col=0._r8
  allocate(RCH4ProdAcet60cm_col(JY,JX)); RCH4ProdAcet60cm_col=0._r8
  allocate(RCH4ProdHG60cm_col(JY,JX)); RCH4ProdHG60cm_col=0._r8
  allocate(RFerment60cm_col(JY,JX)); RFerment60cm_col=0._r8
  allocate(RCO2Ht60cm_col(JY,JX)); RCO2Ht60cm_col=0._r8

  end subroutine InitMLDataType

!----------------------------------------------------------------------

  subroutine DestructMLDataType
  use abortutils, only : destroy
  implicit none

  call destroy(TEMP60cm_col) 
  call destroy(THETW60cm_col)
  call destroy(O2wConc60cm_col) 
  call destroy(CO2wConc60cm_col)   
  call destroy(AcetConc60cm_col)
  call destroy(H2wConc60cm_col)
  call destroy(CH4wConc60cm_col)
  call destroy(DOCConc60cm_col)
  call destroy(AcetMGC60cm_col)
  call destroy(H2MGC60cm_col)
  call destroy(FermC60cm_col)
  call destroy(AeroMOC60cm_col)
  call destroy(AeroHRFungC60cm_col)
  call destroy(AeroHRBactC60cm_col)
  call destroy(RAeroCH4Oxi60cm_col)
  call destroy(RCH4ProdAcet60cm_col)
  call destroy(RCH4ProdHG60cm_col)
  call destroy(RFerment60cm_col)
  call destroy(RCO2Ht60cm_col)

  call destroy(TEMP30cm_col) 
  call destroy(THETW30cm_col)
  call destroy(O2wConc30cm_col) 
  call destroy(CO2wConc30cm_col)   
  call destroy(AcetConc30cm_col)
  call destroy(H2wConc30cm_col)
  call destroy(CH4wConc30cm_col)
  call destroy(DOCConc30cm_col)
  call destroy(AcetMGC30cm_col)
  call destroy(H2MGC30cm_col)
  call destroy(FermC30cm_col)
  call destroy(AeroMOC30cm_col)
  call destroy(AeroHRFungC30cm_col)
  call destroy(AeroHRBactC30cm_col)
  call destroy(RAeroCH4Oxi30cm_col)
  call destroy(RCH4ProdAcet30cm_col)
  call destroy(RCH4ProdHG30cm_col)
  call destroy(RFerment30cm_col)
  call destroy(RCO2Ht30cm_col)

  end subroutine DestructMLDataType

end module MLDataDiagType
