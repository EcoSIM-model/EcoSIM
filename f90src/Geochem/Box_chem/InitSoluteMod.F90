module InitSoluteMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : AZMAX1
  use SoluteChemDataType, only : solutedtype
  use AqueChemDatatype, only : trcSaltIonNumber
  use ChemTracerParsMod
  use SoluteParMod
  use TracerPropMod
  use TracerIDMod
  use EcoSiMParDataMod, only : micpar
  use EcoSIMCtrlMod, only : salt_model
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  real(r8), parameter :: COOH1=2.5E-02_r8
  real(r8), parameter :: TAD=5.0E-02_r8
  real(r8), parameter :: CALMX=10.0_r8
  real(r8), parameter :: CFEMX=10.0_r8
  real(r8), parameter :: ZERO=1.0E-15_r8

  real(r8) :: A1,A2,CEC_conc,XCOOH_conc,SPOH2,SPOH1,NO3_1e_aque_mole_conc
  real(r8), pointer :: H2CO3_aque_mole_conc
  real(r8), pointer :: CH4_aque_mole_conc
  real(r8), pointer :: O2_aque_mole_conc
  real(r8), pointer :: N2_aque_mole_conc
  real(r8), pointer :: N2O_aque_mole_conc
  real(r8), pointer :: NH4_1p_aque_mole_conc
  real(r8), pointer :: NH3_aque_mole_conc
  real(r8), pointer :: Al_3p_aque_mole_conc
  real(r8), pointer :: Fe_3p_aque_mole_conc
  real(r8), pointer :: H_1p_aque_mole_conc
  real(r8), pointer :: Ca_2p_aque_mole_conc
  real(r8), pointer :: Mg_2p_aque_mole_conc
  real(r8), pointer :: Na_1p_aque_mole_conc
  real(r8), pointer :: K_1p_aque_mole_conc
  real(r8), pointer :: OH_1e_aque_mole_conc
  real(r8), pointer :: SO4_2e_aque_mole_conc
  real(r8), pointer :: Cl_e_conc
  real(r8), pointer :: CO3_2e_aque_mole_conc
  real(r8), pointer :: HCO3_e_conc
  real(r8), pointer :: AlOH_2p_aque_mole_conc
  real(r8), pointer :: AlO2H2_1p_aque_mole_conc
  real(r8), pointer :: AlO3H3_conc
  real(r8), pointer :: AlO4H4_1e_aque_mole_conc
  real(r8), pointer :: AlSO4_1p_aque_mole_conc
  real(r8), pointer :: FeOH_2p_aque_mole_conc
  real(r8), pointer :: FeO2H2_p_conc
  real(r8), pointer :: FeO3H3_conc
  real(r8), pointer :: FeO4H4_1e_aque_mole_conc
  real(r8), pointer :: FeSO4_1p_aque_mole_conc
  real(r8), pointer :: CaO2H2_conc
  real(r8), pointer :: CaCO3_conc
  real(r8), pointer :: CaHCO3_1p_aque_mole_conc
  real(r8), pointer :: CaSO4_conc
  real(r8), pointer :: MgOH_1p_aque_mole_conc
  real(r8), pointer :: MgCO3_conc
  real(r8), pointer :: MgHCO3_1p_aque_mole_conc
  real(r8), pointer :: MgSO4_conc
  real(r8), pointer :: NaCO3_1e_aque_mole_conc
  real(r8), pointer :: NaSO4_1e_aque_mole_conc
  real(r8), pointer :: KSO4_1e_aque_mole_conc
  real(r8), pointer :: H0PO4_3e_conc
  real(r8), pointer :: H1PO4_2e_aque_mole_conc
  real(r8), pointer :: H2PO4_1e_aque_mole_conc
  real(r8), pointer :: H3PO4_conc
  real(r8), pointer :: FeHPO4_p_conc
  real(r8), pointer :: FeH2PO4_2p_aque_mole_conc
  real(r8), pointer :: CaPO4_1e_con
  real(r8), pointer :: CaHPO4_conc
  real(r8), pointer :: CaH4P2O8_1p_aque_mole_conc
  real(r8), pointer :: MgHPO4_conc
  real(r8), pointer :: CSTR1
  real(r8), pointer :: CCO2M
  real(r8), pointer :: CCH4M
  real(r8), pointer :: COXYM
  real(r8), pointer :: CZ2GM
  real(r8), pointer :: CZ2OM
  real(r8), pointer :: CN4Z
  real(r8), pointer :: CNOZ
  real(r8), pointer :: CNAZ
  real(r8), pointer :: CKAZ
  real(r8), pointer :: CSOZ
  real(r8), pointer :: CCLZ
  real(r8), pointer :: CNOX
  real(r8), pointer :: CCASOX
  real(r8), pointer :: CN4X
  real(r8), pointer :: CPOZ
  real(r8), pointer :: Al_mole_conc
  real(r8), pointer :: CFEZ
  real(r8), pointer :: CCAZ
  real(r8), pointer :: CMGZ
  real(r8), pointer :: CALX
  real(r8), pointer :: CFEX
  real(r8), pointer :: CaX_conc
  real(r8), pointer :: MgX_conc
  real(r8), pointer :: CNAX
  real(r8), pointer :: CKAX
  real(r8), pointer :: CSOX
  real(r8), pointer :: CCLX
  real(r8), pointer :: CALPOX
  real(r8), pointer :: CFEPOX
  real(r8), pointer :: CCAPDX
  real(r8), pointer :: CCAPHX
  real(r8), pointer :: CALOHX
  real(r8), pointer :: CFEOHX
  real(r8), pointer :: CCACOX
  real(r8), pointer :: XNH4_mole_conc
  real(r8), pointer :: XHY1
  real(r8), pointer :: XAl_conc
  real(r8), pointer :: XFe_conc
  real(r8), pointer :: XCa_conc
  real(r8), pointer :: Precp_Ca5P3O12O3H3_mole_conc
  real(r8), pointer :: XMg_conc
  real(r8), pointer :: XNa_conc
  real(r8), pointer :: XK_conc
  real(r8), pointer :: XHC1
  real(r8), pointer :: XAlO2H2_conc
  real(r8), pointer :: XFeO2H2_conc
  real(r8), pointer :: XOH_mole_conc
  real(r8), pointer :: XROH1_mole_conc
  real(r8), pointer :: XROH2_mole_conc
  real(r8), pointer :: XHPO4_mole_conc
  real(r8), pointer :: XH2PO4_mole_conc
  real(r8), pointer :: Precp_AlO3H3_conc
  real(r8), pointer :: Precp_FeO3H3_conc
  real(r8), pointer :: Precp_CaCO3_conc
  real(r8), pointer :: Precp_CaSO4_conc
  real(r8), pointer :: Precp_AlPO4_mole_conc
  real(r8), pointer :: Precp_FePO4_mole_conc
  real(r8), pointer :: Precp_CaHPO4_mole_conc
  real(r8), pointer :: FH2O
  real(r8), pointer :: ATCA
  real(r8), pointer :: XAEC
  real(r8), pointer :: CEC
  real(r8), pointer :: ORGC
  real(r8), pointer :: VLPO4
  real(r8), pointer :: XCEC
  real(r8), pointer :: GKC4
  real(r8), pointer :: GKCA
  real(r8), pointer :: GKCH
  real(r8), pointer :: GKCK
  real(r8), pointer :: GKCN
  real(r8), pointer :: GKCM
  real(r8), pointer :: ZEROS
  real(r8), pointer :: VLWatMicP

  public :: InitSoluteModel
  public :: InitSoluteProperty
  contains

  subroutine InitSoluteProperty
  implicit none

  trcSaltIonNumber(idsalt_Al)        = 1._r8
  trcSaltIonNumber(idsalt_Fe)        = 1._r8
  trcSaltIonNumber(idsalt_Hp)        = 1._r8
  trcSaltIonNumber(idsalt_Ca)        = 1._r8
  trcSaltIonNumber(idsalt_Mg)        = 1._r8
  trcSaltIonNumber(idsalt_Na)        = 1._r8
  trcSaltIonNumber(idsalt_K)         = 1._r8
  trcSaltIonNumber(idsalt_OH)        = 1._r8
  trcSaltIonNumber(idsalt_SO4)       = 1._r8
  trcSaltIonNumber(idsalt_Cl)        = 1._r8
  trcSaltIonNumber(idsalt_CO3)       = 1._r8
  trcSaltIonNumber(idsalt_H0PO4)     = 1._r8
  trcSaltIonNumber(idsalt_H0PO4B)    = 1._r8
  trcSaltIonNumber(idsalt_HCO3)      = 2._r8
  trcSaltIonNumber(idsalt_AlOH)      = 2._r8
  trcSaltIonNumber(idsalt_AlSO4)     = 2._r8
  trcSaltIonNumber(idsalt_FeOH)      = 2._r8
  trcSaltIonNumber(idsalt_FeSO4)     = 2._r8
  trcSaltIonNumber(idsalt_CaOH)      = 2._r8
  trcSaltIonNumber(idsalt_CaCO3)     = 2._r8
  trcSaltIonNumber(idsalt_CaSO4)     = 2._r8
  trcSaltIonNumber(idsalt_MgOH2)     = 3._r8
  trcSaltIonNumber(idsalt_MgCO3)     = 2._r8
  trcSaltIonNumber(idsalt_MgSO4)     = 2._r8
  trcSaltIonNumber(idsalt_NaCO3)     = 2._r8
  trcSaltIonNumber(idsalt_NaSO4)     = 2._r8
  trcSaltIonNumber(idsalt_KSO4)      = 2._r8
  trcSaltIonNumber(idsalt_CaPO4)     = 2._r8
  trcSaltIonNumber(idsalt_CaPO4B)    = 2._r8
  trcSaltIonNumber(idsalt_AlOH2)     = 3._r8
  trcSaltIonNumber(idsalt_FeOH2)     = 3._r8
  trcSaltIonNumber(idsalt_CaHCO3)    = 3._r8
  trcSaltIonNumber(idsalt_MgHCO3)    = 3._r8
  trcSaltIonNumber(idsalt_FeHPO4)    = 3._r8
  trcSaltIonNumber(idsalt_CaHPO4)    = 3._r8
  trcSaltIonNumber(idsalt_MgHPO4)    = 3._r8
  trcSaltIonNumber(idsalt_FeHPO4B)   = 3._r8
  trcSaltIonNumber(idsalt_CaHPO4B)   = 3._r8
  trcSaltIonNumber(idsalt_MgHPO4B)   = 3._r8
  trcSaltIonNumber(idsalt_AlOH3)     = 4._r8
  trcSaltIonNumber(idsalt_FeOH3)     = 4._r8
  trcSaltIonNumber(idsalt_H3PO4)     = 4._r8
  trcSaltIonNumber(idsalt_FeH2PO4)   = 4._r8
  trcSaltIonNumber(idsalt_CaH4P2O8)  = 4._r8
  trcSaltIonNumber(idsalt_H3PO4B)    = 4._r8
  trcSaltIonNumber(idsalt_FeH2PO4B)  = 4._r8
  trcSaltIonNumber(idsalt_CaH4P2O8B) = 4._r8
  trcSaltIonNumber(idsalt_AlOH4)     = 5._r8
  trcSaltIonNumber(idsalt_FeOH4)     = 5._r8

  end subroutine InitSoluteProperty
!------------------------------------------------------------------------------------------

  SUBROUTINE InitSoluteModel(K,BulkSoilMass,solutevar)
!
!     THIS SUBROUTINE INITIALIZES ALL SOIL CHEMISTRY VARIABLES
!

  implicit none
  integer, intent(in) :: K
  real(r8), intent(in) :: BulkSoilMass
  type(solutedtype), target, intent(inout) :: solutevar

  integer, parameter :: MRXN1=1000
  integer :: M
!     begin_execution
  H2CO3_aque_mole_conc       => solutevar%H2CO3_aque_mole_conc
  CH4_aque_mole_conc         => solutevar%CH4_aque_mole_conc
  O2_aque_mole_conc          => solutevar%O2_aque_mole_conc
  N2_aque_mole_conc          => solutevar%N2_aque_mole_conc
  N2O_aque_mole_conc         => solutevar%N2O_aque_mole_conc
  NH4_1p_aque_mole_conc      => solutevar%NH4_1p_aque_mole_conc
  NH3_aque_mole_conc         => solutevar%NH3_aque_mole_conc
  Al_3p_aque_mole_conc       => solutevar%Al_3p_aque_mole_conc
  Fe_3p_aque_mole_conc       => solutevar%Fe_3p_aque_mole_conc
  H_1p_aque_mole_conc        => solutevar%H_1p_aque_mole_conc
  Ca_2p_aque_mole_conc       => solutevar%Ca_2p_aque_mole_conc
  Mg_2p_aque_mole_conc       => solutevar%Mg_2p_aque_mole_conc
  Na_1p_aque_mole_conc       => solutevar%Na_1p_aque_mole_conc
  K_1p_aque_mole_conc        => solutevar%K_1p_aque_mole_conc
  OH_1e_aque_mole_conc       => solutevar%OH_1e_aque_mole_conc
  SO4_2e_aque_mole_conc      => solutevar%SO4_2e_aque_mole_conc
  Cl_e_conc                  => solutevar%Cl_e_conc
  CO3_2e_aque_mole_conc      => solutevar%CO3_2e_aque_mole_conc
  HCO3_e_conc                => solutevar%HCO3_e_conc
  AlOH_2p_aque_mole_conc     => solutevar%AlOH_2p_aque_mole_conc
  AlO2H2_1p_aque_mole_conc   => solutevar%AlO2H2_1p_aque_mole_conc
  AlO3H3_conc                => solutevar%AlO3H3_conc
  AlO4H4_1e_aque_mole_conc   => solutevar%AlO4H4_1e_aque_mole_conc
  AlSO4_1p_aque_mole_conc    => solutevar%AlSO4_1p_aque_mole_conc
  FeOH_2p_aque_mole_conc     => solutevar%FeOH_2p_aque_mole_conc
  FeO2H2_p_conc              => solutevar%FeO2H2_p_conc
  FeO3H3_conc                => solutevar%FeO3H3_conc
  FeO4H4_1e_aque_mole_conc   => solutevar%FeO4H4_1e_aque_mole_conc
  FeSO4_1p_aque_mole_conc    => solutevar%FeSO4_1p_aque_mole_conc
  CaO2H2_conc                => solutevar%CaO2H2_conc
  CaCO3_conc                 => solutevar%CaCO3_conc
  CaHCO3_1p_aque_mole_conc   => solutevar%CaHCO3_1p_aque_mole_conc
  CaSO4_conc                 => solutevar%CaSO4_conc
  MgOH_1p_aque_mole_conc     => solutevar%MgOH_1p_aque_mole_conc
  MgCO3_conc                 => solutevar%MgCO3_conc
  MgHCO3_1p_aque_mole_conc   => solutevar%MgHCO3_1p_aque_mole_conc
  MgSO4_conc                 => solutevar%MgSO4_conc
  NaCO3_1e_aque_mole_conc    => solutevar%NaCO3_1e_aque_mole_conc
  NaSO4_1e_aque_mole_conc    => solutevar%NaSO4_1e_aque_mole_conc
  KSO4_1e_aque_mole_conc     => solutevar%KSO4_1e_aque_mole_conc
  H0PO4_3e_conc              => solutevar%H0PO4_3e_conc
  H1PO4_2e_aque_mole_conc    => solutevar%H1PO4_2e_aque_mole_conc
  H2PO4_1e_aque_mole_conc    => solutevar%H2PO4_1e_aque_mole_conc
  H3PO4_conc                 => solutevar%H3PO4_conc
  FeHPO4_p_conc              => solutevar%FeHPO4_p_conc
  FeH2PO4_2p_aque_mole_conc  => solutevar%FeH2PO4_2p_aque_mole_conc
  CaPO4_1e_con               => solutevar%CaPO4_1e_con
  CaHPO4_conc                => solutevar%CaHPO4_conc
  CaH4P2O8_1p_aque_mole_conc => solutevar%CaH4P2O8_1p_aque_mole_conc
  MgHPO4_conc                => solutevar%MgHPO4_conc
  CSTR1                      => solutevar%CSTR1
  CCO2M                      => solutevar%CCO2M
  CCH4M                      => solutevar%CCH4M
  COXYM                      => solutevar%COXYM
  CZ2GM                      => solutevar%CZ2GM
  CZ2OM                      => solutevar%CZ2OM
  CN4Z                       => solutevar%CN4Z
  CNOZ                       => solutevar%CNOZ
  CNAZ                       => solutevar%CNAZ
  CKAZ                       => solutevar%CKAZ
  CSOZ                       => solutevar%CSOZ
  CCLZ                       => solutevar%CCLZ
  CNOX                       => solutevar%CNOX
  CCASOX                     => solutevar%CCASOX
  CN4X                       => solutevar%CN4X
  CPOZ                       => solutevar%CPOZ
  Al_mole_conc                       => solutevar%Al_mole_conc
  CFEZ                       => solutevar%CFEZ
  CCAZ                       => solutevar%CCAZ
  CMGZ                       => solutevar%CMGZ
  CALX                       => solutevar%CALX
  CFEX                       => solutevar%CFEX
  CaX_conc                   => solutevar%CaX_conc
  MgX_conc                   => solutevar%MgX_conc
  CNAX                       => solutevar%CNAX
  CKAX                       => solutevar%CKAX
  CSOX                       => solutevar%CSOX
  CCLX                       => solutevar%CCLX
  CALPOX                     => solutevar%CALPOX
  CFEPOX                     => solutevar%CFEPOX
  CCAPDX                     => solutevar%CCAPDX
  CCAPHX                     => solutevar%CCAPHX
  CALOHX                     => solutevar%CALOHX
  CFEOHX                     => solutevar%CFEOHX
  CCACOX                     => solutevar%CCACOX
  XNH4_mole_conc                  => solutevar%XNH4_mole_conc
  XHY1                       => solutevar%XHY1
  XAl_conc                   => solutevar%XAl_conc
  XFe_conc                   => solutevar%XFe_conc
  XCa_conc                   => solutevar%XCa_conc
  Precp_Ca5P3O12O3H3_mole_conc    => solutevar%Precp_Ca5P3O12O3H3_mole_conc
  XMg_conc                   => solutevar%XMg_conc
  XNa_conc                   => solutevar%XNa_conc
  XK_conc                    => solutevar%XK_conc
  XHC1                       => solutevar%XHC1
  XAlO2H2_conc               => solutevar%XAlO2H2_conc
  XFeO2H2_conc               => solutevar%XFeO2H2_conc
  XOH_mole_conc                   => solutevar%XOH_mole_conc
  XROH1_mole_conc                 => solutevar%XROH1_mole_conc
  XROH2_mole_conc                 => solutevar%XROH2_mole_conc
  XHPO4_mole_conc                 => solutevar%XHPO4_mole_conc
  XH2PO4_mole_conc                => solutevar%XH2PO4_mole_conc
  Precp_AlO3H3_conc          => solutevar%Precp_AlO3H3_conc
  Precp_FeO3H3_conc          => solutevar%Precp_FeO3H3_conc
  Precp_CaCO3_conc           => solutevar%Precp_CaCO3_conc
  Precp_CaSO4_conc           => solutevar%Precp_CaSO4_conc
  Precp_AlPO4_mole_conc           => solutevar%Precp_AlPO4_mole_conc
  Precp_FePO4_mole_conc           => solutevar%Precp_FePO4_mole_conc
  Precp_CaHPO4_mole_conc          => solutevar%Precp_CaHPO4_mole_conc
  FH2O                       => solutevar%FH2O
  ATCA                       => solutevar%ATCA
  XAEC                       => solutevar%XAEC
  CEC                        => solutevar%CEC
  ORGC                       => solutevar%ORGC
  VLPO4                      => solutevar%VLPO4
  XCEC                       => solutevar%XCEC
  GKC4                       => solutevar%GKC4
  GKCA                       => solutevar%GKCA
  GKCH                       => solutevar%GKCH
  GKCK                       => solutevar%GKCK
  GKCN                       => solutevar%GKCN
  GKCM                       => solutevar%GKCM
  ZEROS                      => solutevar%ZEROS
  VLWatMicP                  => solutevar%VLWatMicP

  call InitEquilibria(K,BulkSoilMass)
!
!     CONVERGE TOWARDS ALL SOLUBILITY EQUILIBRIA
!     IF SALT OPTION IS SELECTED
!
  IF(salt_model)THEN
    DO   M=1,MRXN1
      call SolubilityEquilibiriaSalt(K,M,BulkSoilMass)
    ENDDO
!
!     CONVERGE TOWARDS ALL SOLUBILITY EQUILIBRIA
!     IF SALT OPTION IS NOT SELECTED
!
  ELSE
    DO  M=1,MRXN1
      call SolubilityEquilibriaNoSalt(K,M,BulkSoilMass)
    ENDDO
  ENDIF
  end SUBROUTINE InitSoluteModel
!------------------------------------------------------------------------------------------

  subroutine InitEquilibria(K,BulkSoilMass)

  implicit none
  integer, intent(in) :: K
  real(r8), intent(in) :: BulkSoilMass
  integer :: MM
  real(r8) :: XNAQ,XPT,XN4Q,XMGQ,XKAQ,XHYQ,XFEQ
  real(r8) :: XHP,XCECQ,XTLQ,XALQ,XCAX
  real(r8) :: XCAQ,XOH,SPH1P,FX,FXH2,FXP1,FSTR2
  real(r8) :: FHPA,FXH1,FXP2,FHP1,FHP3,FHCO
  real(r8) :: FCO3,CX1,CX2,CSTRZ,FHP2,CCO2Y
  real(r8) :: CION2,A3,anion_1e_aque_mole_conc,anion_2e_aque_mole_conc,anion_3e_conc,cation_1p_aque_mole_conc,cation_2p_aque_mole_conc,cation_3p_aque_mole_conc
  real(r8) :: CCO2X,CCO2Z,CN,CSTR2,FHP0,FXH0,R,Z
  real(r8) :: OH_activity
  real(r8) :: AdsorpEquilTerm,XPTDiscriminant,XPTDenom
!
!     INITIALIZE SOLUTE EQUILIBRIA
!
  cation_3p_aque_mole_conc=AZMAX1(Al_mole_conc)+AZMAX1(CFEZ)
  anion_3e_conc=0._r8
  cation_2p_aque_mole_conc=AZMAX1(CCAZ)+CMGZ
  anion_2e_aque_mole_conc=CSOZ
  cation_1p_aque_mole_conc=CN4Z+CNAZ+CKAZ
  anion_1e_aque_mole_conc=CNOZ+CCLZ+CPOZ
  CX2=0._r8
  CX1=0._r8
  CN=0._r8
!
!     INITIALIZE ION STRENGTH AND ACTIVITIES
!
  CION2=AZMAX1(cation_3p_aque_mole_conc+anion_3e_conc+cation_2p_aque_mole_conc+anion_2e_aque_mole_conc+cation_1p_aque_mole_conc+anion_1e_aque_mole_conc+CN)
  CSTR1=0.5E-03*(9.0*(cation_3p_aque_mole_conc+anion_3e_conc)+4.0*(cation_2p_aque_mole_conc+anion_2e_aque_mole_conc)+cation_1p_aque_mole_conc+anion_1e_aque_mole_conc)
  CSTRZ=0.5E-03*(9.0*(cation_3p_aque_mole_conc+anion_3e_conc)+4.0*(cation_2p_aque_mole_conc+CX2)+cation_1p_aque_mole_conc+CX1)
  CSTR2=SQRT(CSTR1)
  FSTR2=CSTR2/(1.0_r8+CSTR2)
  FH2O=5.56E+04/(5.56E+04+CION2)

  IF(salt_model)THEN
    A1=AMIN1(1.0,10.0**(-0.509*1.0*FSTR2+0.20*CSTR2))
    A2=AMIN1(1.0,10.0**(-0.509*4.0*FSTR2+0.20*CSTR2))
    A3=AMIN1(1.0,10.0**(-0.509*9.0*FSTR2+0.20*CSTR2))
  ELSE
    A1=1.0_r8
    A2=1.0_r8
    A3=1.0_r8
  ENDIF
  ! Keep prescribed [H(+)], but enforce a_H*a_OH=DPH2O in salt chemistry.
  ! Both monovalent activity coefficients are A1; no-salt inputs are unchanged.
  IF(salt_model) OH_1e_aque_mole_conc=DPH2O/(A1*A1*H_1p_aque_mole_conc)
  ! Hydroxide reaction constants use a_OH = gamma_OH*[OH], including each power.
  OH_activity=OH_1e_aque_mole_conc*A1
!
!     INITIALIZE GASES
!

  ! CO2(g) <=> CO2(aq), represented by H2CO3 here.
  ! Henry partitioning includes temperature, salting-out and the coupled carbonate equilibria.
  CCO2X          = CCO2M*gas_solubility(idg_CO2,ATCA)/(EXP(GasSechenovConst(idg_CO2)*CSTRZ))*FH2O
  CCO2Y          = LOG(CCO2X)
  CCO2Z          = ABS(CCO2Y)
  H2CO3_aque_mole_conc = CCO2X
  FCO3           = DPCO3*A0/(H_1p_aque_mole_conc**2._r8*A2)
  FHCO           = DPCO2*A0/(H_1p_aque_mole_conc*A1)
  Z              = GasSechenovConst(idg_CO2)*(2.0E-03_r8*FCO3+0.5E-03_r8*FHCO)
  DO MM   = 1, 25
    R=(LOG(H2CO3_aque_mole_conc)+Z*H2CO3_aque_mole_conc-CCO2Y)/CCO2Z
    IF(R.LT.1.0E-03_r8)exit
    H2CO3_aque_mole_conc=H2CO3_aque_mole_conc/SQRT(1.0_r8+R)
  ENDDO

  ! CH4(g) <=> CH4(aq): initialize dissolved gas by Henry partitioning.
  CH4_aque_mole_conc = CCH4M*gas_solubility(idg_CH4,ATCA)/(EXP(GasSechenovConst(idg_CH4)*CSTR1))*FH2O
  ! O2(g) <=> O2(aq): initialize dissolved gas by Henry partitioning.
  O2_aque_mole_conc  = COXYM*gas_solubility(idg_O2,ATCA)/(EXP(GasSechenovConst(idg_O2)*CSTR1))*FH2O
  ! N2(g) <=> N2(aq): initialize dissolved gas by Henry partitioning.
  N2_aque_mole_conc  = CZ2GM*gas_solubility(idg_N2,ATCA)/(EXP(GasSechenovConst(idg_N2)*CSTR1))*FH2O
  ! N2O(g) <=> N2O(aq): initialize dissolved gas by Henry partitioning.
  N2O_aque_mole_conc = CZ2OM*gas_solubility(idg_N2O,ATCA)/(EXP(GasSechenovConst(idg_N2O)*CSTR1))*FH2O

  ! H2CO3(aq) <=> CO3(2-) + 2H(+): combine the two carbonate acid dissociations.
  CO3_2e_aque_mole_conc=H2CO3_aque_mole_conc*DPCO3*A0/(H_1p_aque_mole_conc**2._r8*A2)
  ! H2CO3(aq) <=> HCO3(-) + H(+): evaluate bicarbonate at prescribed pH.
  HCO3_e_conc=H2CO3_aque_mole_conc*DPCO2*A0/(H_1p_aque_mole_conc*A1)
  NO3_1e_aque_mole_conc=CNOZ
!
!     INITIALIZE ION PAIR EQUILIBRIA
!
  ! NH4(+) <=> NH3(aq) + H(+): partition total dissolved ammonia at prescribed pH.
  IF(K.NE.micpar%k_POM)THEN
    NH4_1p_aque_mole_conc=CN4Z/(1.0_r8+DPN4*A1/(H_1p_aque_mole_conc*A0))
    NH3_aque_mole_conc=NH4_1p_aque_mole_conc*DPN4*A1/(H_1p_aque_mole_conc*A0)
  ELSE
    NH4_1p_aque_mole_conc=ZERO
    NH3_aque_mole_conc=ZERO
  ENDIF
  ! Al(OH)3(s) <=> Al(3+) + 3OH(-): initialize the free ion subject to its mineral solubility cap.
  IF(Al_mole_conc.LE.0.0_r8)THEN
    Al_3p_aque_mole_conc=AMIN1(CALMX,SPALO/(OH_activity**3*A3))
  ELSE
    Al_3p_aque_mole_conc=AMIN1(Al_mole_conc,SPALO/(OH_activity**3*A3))
  ENDIF
  ! Fe(OH)3(s) <=> Fe(3+) + 3OH(-): initialize the free ion subject to its mineral solubility cap.
  IF(CFEZ.LE.0.0_r8)THEN
    Fe_3p_aque_mole_conc=AMIN1(CFEMX,SPFEO/(OH_activity**3*A3))
  ELSE
    Fe_3p_aque_mole_conc=AMIN1(CFEZ,SPFEO/(OH_activity**3*A3))
  ENDIF
  ! CaCO3(s) <=> Ca(2+) + CO3(2-): initialize the free ion subject to its mineral solubility cap.
  IF(CCAZ.LE.0.0_r8)THEN
    Ca_2p_aque_mole_conc=AMIN1(CCAMX,SPCAC/(CO3_2e_aque_mole_conc*A2**2))
  ELSE
    Ca_2p_aque_mole_conc=AMIN1(CCAZ,SPCAC/(CO3_2e_aque_mole_conc*A2**2))
  ENDIF
  Mg_2p_aque_mole_conc   = CMGZ
  Na_1p_aque_mole_conc             = CNAZ
  K_1p_aque_mole_conc              = CKAZ
  SO4_2e_aque_mole_conc            = CSOZ
  Cl_e_conc              = CCLZ

  ! Al(3+) + OH(-) <=> AlOH(2+): initialize the hydrolysis product.
  AlOH_2p_aque_mole_conc = Al_3p_aque_mole_conc*OH_activity*A3/(DPAL1*A2)

  ! AlOH(2+) + OH(-) <=> Al(OH)2(+): initialize the hydrolysis product.
  AlO2H2_1p_aque_mole_conc         = Al_3p_aque_mole_conc*OH_activity**2*A3/(DPAL1*DPAL2*A1)
  
  ! Al(OH)2(+) + OH(-) <=> Al(OH)3(aq): initialize the hydrolysis product.
  AlO3H3_conc            = Al_3p_aque_mole_conc*OH_activity**3*A3/(DPAL1*DPAL2*DPAL3*A0)
  
  ! Al(OH)3(aq) + OH(-) <=> Al(OH)4(-): initialize the hydrolysis product.
  AlO4H4_1e_aque_mole_conc         = Al_3p_aque_mole_conc*OH_activity**4*A3/(DPAL1*DPAL2*DPAL3*DPAL4*A1)
  AlSO4_1p_aque_mole_conc          = 0._r8
  
  ! Fe(3+) + OH(-) <=> FeOH(2+): initialize the hydrolysis product.
  FeOH_2p_aque_mole_conc = Fe_3p_aque_mole_conc*OH_activity*A3/(DPFE1*A2)
  
  ! FeOH(2+) + OH(-) <=> Fe(OH)2(+): initialize the hydrolysis product.
  FeO2H2_p_conc          = Fe_3p_aque_mole_conc*OH_activity**2*A3/(DPFE1*DPFE2*A1)
  
  ! Fe(OH)2(+) + OH(-) <=> Fe(OH)3(aq): initialize the hydrolysis product.
  FeO3H3_conc            = Fe_3p_aque_mole_conc*OH_activity**3*A3/(DPFE1*DPFE2*DPFE3*A0)
  
  ! Fe(OH)3(aq) + OH(-) <=> Fe(OH)4(-): initialize the hydrolysis product.
  FeO4H4_1e_aque_mole_conc         = Fe_3p_aque_mole_conc*OH_activity**4*A3/(DPFE1*DPFE2*DPFE3*DPFE4*A1)
  FeSO4_1p_aque_mole_conc          = 0._r8
  
  ! Ca(2+) + OH(-) <=> CaOH(+): initialize the aqueous complex.
  CaO2H2_conc            = Ca_2p_aque_mole_conc*OH_activity*A2/(DPCAO*A1)
  
  ! Ca(2+) + CO3(2-) <=> CaCO3(aq): initialize the aqueous complex.
  CaCO3_conc=Ca_2p_aque_mole_conc*CO3_2e_aque_mole_conc*A2**2/(DPCAC*A0)
  
  ! Ca(2+) + HCO3(-) <=> CaHCO3(+): initialize the aqueous complex.
  CaHCO3_1p_aque_mole_conc=Ca_2p_aque_mole_conc*HCO3_e_conc*A2/DPCAH
  CaSO4_conc=0._r8
  
  ! Mg(2+) + OH(-) <=> MgOH(+): initialize the aqueous complex.
  MgOH_1p_aque_mole_conc=Mg_2p_aque_mole_conc*OH_activity*A2/(DPMGO*A1)
  
  ! Mg(2+) + CO3(2-) <=> MgCO3(aq): initialize the aqueous complex.
  MgCO3_conc=Mg_2p_aque_mole_conc*CO3_2e_aque_mole_conc*A2**2/(DPMGC*A0)
  
  ! Mg(2+) + HCO3(-) <=> MgHCO3(+): initialize the aqueous complex.
  MgHCO3_1p_aque_mole_conc=Mg_2p_aque_mole_conc*HCO3_e_conc*A2/DPMGH
  MgSO4_conc=0._r8
  
  ! Na(+) + CO3(2-) <=> NaCO3(-): initialize the aqueous complex.
  NaCO3_1e_aque_mole_conc=Na_1p_aque_mole_conc*CO3_2e_aque_mole_conc*A2/DPNAC
  NaSO4_1e_aque_mole_conc=0._r8
  KSO4_1e_aque_mole_conc=0._r8
  FeHPO4_p_conc=0._r8
  FeH2PO4_2p_aque_mole_conc=0._r8
  CaPO4_1e_con=0._r8
  CaHPO4_conc=0._r8
  CaH4P2O8_1p_aque_mole_conc=0._r8
  MgHPO4_conc=0._r8
!
!     INITIALIZE PHOSPHORUS EQUILIBRIA AMONG SOLUBLE, ADSORBED
!     AND PRECIPITATED FORMS
  ! Phosphate acid equilibria: H3PO4 <=> H(+) + H2PO4(-),
  ! H2PO4(-) <=> H(+) + HPO4(2-), and HPO4(2-) <=> H(+) + PO4(3-).
  ! Salt mode represents all four forms. No-salt mode represents only H2PO4
  ! and HPO4, so normalize their fractions to retain all dissolved P.
  IF(salt_model)THEN
    FHP3=1.0_r8/(1.0_r8+DPH3P*A0/(H_1p_aque_mole_conc*A1)+DPH3P*DPH2P*A0 &
      /(H_1p_aque_mole_conc**2*A2)+DPH3P*DPH2P*DPH1P*A0/(H_1p_aque_mole_conc**3*A3))
    FHP2=FHP3*DPH3P*A0/(H_1p_aque_mole_conc*A1)
    FHP1=FHP3*DPH3P*DPH2P*A0/(H_1p_aque_mole_conc**2*A2)
    FHP0=FHP3*DPH3P*DPH2P*DPH1P*A0/(H_1p_aque_mole_conc**3*A3)
  ELSE
    FHP3=0._r8
    FHP2=H_1p_aque_mole_conc/(H_1p_aque_mole_conc+DPH2P)
    FHP1=DPH2P/(H_1p_aque_mole_conc+DPH2P)
    FHP0=0._r8
  ENDIF
  IF(K.NE.micpar%k_POM)THEN
    ! Boundary water: distribute the supplied P among represented dissolved forms.
    H3PO4_conc=CPOZ*FHP3
    H2PO4_1e_aque_mole_conc=CPOZ*FHP2
    H1PO4_2e_aque_mole_conc=CPOZ*FHP1
    H0PO4_3e_conc=CPOZ*FHP0
  ELSE
    ! Solve the same phosphate acid equilibria together with surface-site protonation:
    ! R-OH2(+) <=> R-OH + H(+); R-OH <=> R-O(-) + H(+).
    ! Ligand exchange: R-OH + H2PO4(-) <=> R-H2PO4 + OH(-),
    ! and R-OH + HPO4(2-) <=> R-HPO4(-) + OH(-).
    XHP=CPOZ
    XOH=XAEC/BulkSoilMass
    SPOH2=SXOH2/A1
    SPOH1=SXOH1/A1
    SPH1P=SXH1P*DPH2O*A1/A2
    FXH2=1.0_r8/(1.0_r8+SPOH2/H_1p_aque_mole_conc+SPOH2*SPOH1/H_1p_aque_mole_conc**2)
    FXH1=FXH2*SPOH2/H_1p_aque_mole_conc
    FXH0=FXH1*SPOH1/H_1p_aque_mole_conc
    FXP2=1.0_r8/(1.0_r8+SXH2P*DPH2P/(SXH1P*H_1p_aque_mole_conc))
    FXP1=FXP2*SXH2P*DPH2P/(SXH1P*H_1p_aque_mole_conc)
    FHPA=FHP2*A1
    ! Coupled phosphate adsorption above: solve total occupied sites XPT,
    ! then distribute free sites and sorbed phosphate among their protonation states.
    ! Rationalize the smaller root of (XOH-XPT)*(XHP-XPT)=AdsorpEquilTerm*XPT.
    ! This retains tiny P inputs and gives exactly zero adsorption when XHP=0.
    AdsorpEquilTerm=SXH2P*FXP2*OH_1e_aque_mole_conc/(FXH1*FHPA)
    XPTDiscriminant=(XOH-XHP)**2+2._r8*AdsorpEquilTerm*(XOH+XHP)+AdsorpEquilTerm**2
    XPTDenom=XOH+XHP+AdsorpEquilTerm+SQRT(MAX(0._r8,XPTDiscriminant))
    IF(XPTDenom.GT.0._r8)THEN
      XPT=2._r8*(XOH/XPTDenom)*XHP
    ELSE
      XPT=0._r8
    ENDIF
    XROH2_mole_conc=(XOH-XPT)*FXH2
    XROH1_mole_conc=(XOH-XPT)*FXH1
    XOH_mole_conc=(XOH-XPT)*FXH0
    XHPO4_mole_conc=XPT*FXP1
    XH2PO4_mole_conc=XPT*FXP2
    H3PO4_conc=(XHP-XPT)*FHP3
    H2PO4_1e_aque_mole_conc=(XHP-XPT)*FHP2
    H1PO4_2e_aque_mole_conc=(XHP-XPT)*FHP1
    H0PO4_3e_conc=(XHP-XPT)*FHP0
!
!     INITIALIZE CATION EQILIBRIA BETWEEN SOLUBLE
!     AND EXCHANGEABLE FORMS
!
    ! Cation exchange uses Ca as the Gapon reference: 2X-NH4 + Ca(2+)
    ! <=> X2-Ca + 2NH4(+), with analogous H, Al, Fe, Mg, Na and K exchange.
    ! Initialize the exchangeable pools; X denotes one negative exchange-site charge.
    XCECQ=AMAX1(CN4X,CEC)
    XN4Q=CN4X
    XHYQ=0._r8
    XALQ=0._r8
    XFEQ=0._r8
    XCAQ=0._r8
    XMGQ=0._r8
    XNAQ=0._r8
    XKAQ=0._r8
    XHC1=0._r8
    XAlO2H2_conc=0._r8
    XFeO2H2_conc=0._r8
    XCOOH_conc=AZMAX1(COOH1*ORGC)
  ENDIF
  cation_3p_aque_mole_conc=Al_3p_aque_mole_conc+Fe_3p_aque_mole_conc
  anion_3e_conc=H0PO4_3e_conc
  cation_2p_aque_mole_conc=Ca_2p_aque_mole_conc+Mg_2p_aque_mole_conc+AlOH_2p_aque_mole_conc+FeOH_2p_aque_mole_conc+FeH2PO4_2p_aque_mole_conc
  anion_2e_aque_mole_conc=SO4_2e_aque_mole_conc+CO3_2e_aque_mole_conc+H1PO4_2e_aque_mole_conc
  cation_1p_aque_mole_conc=NH4_1p_aque_mole_conc+H_1p_aque_mole_conc+Na_1p_aque_mole_conc+K_1p_aque_mole_conc+AlO2H2_1p_aque_mole_conc+FeO2H2_p_conc+AlSO4_1p_aque_mole_conc &
    +FeSO4_1p_aque_mole_conc+CaO2H2_conc+CaHCO3_1p_aque_mole_conc+MgOH_1p_aque_mole_conc+MgHCO3_1p_aque_mole_conc+FeHPO4_p_conc+CaH4P2O8_1p_aque_mole_conc
  anion_1e_aque_mole_conc=NO3_1e_aque_mole_conc+OH_1e_aque_mole_conc+HCO3_e_conc+Cl_e_conc+AlO4H4_1e_aque_mole_conc+FeO4H4_1e_aque_mole_conc+NaCO3_1e_aque_mole_conc &
    +NaSO4_1e_aque_mole_conc+KSO4_1e_aque_mole_conc+H2PO4_1e_aque_mole_conc+CaPO4_1e_con
  CN=H2CO3_aque_mole_conc+CH4_aque_mole_conc+O2_aque_mole_conc+N2_aque_mole_conc+N2O_aque_mole_conc+NH3_aque_mole_conc &
    +AlO3H3_conc+FeO3H3_conc+CaCO3_conc+CaSO4_conc &
    +MgCO3_conc+MgSO4_conc+H3PO4_conc+CaHPO4_conc+MgHPO4_conc
  CX2=anion_2e_aque_mole_conc-CO3_2e_aque_mole_conc
  CX1=anion_1e_aque_mole_conc-HCO3_e_conc
!
!     INITIALIZE EQUILIBRIA BETWEEN SOLUBLE AND PRECIPITATED FORMS
! POM complex
  ! Seed the solid reservoirs for Al(OH)3, Fe(OH)3, CaCO3, CaSO4,
  ! AlPO4, FePO4, CaHPO4 and Ca5(PO4)3OH; their dissolution reactions
  ! are solved in the salt/no-salt iteration sections below.
  IF(K.EQ.micpar%k_POM)THEN
    Precp_AlO3H3_conc=CALOHX
    Precp_FeO3H3_conc=CFEOHX
    Precp_CaCO3_conc=CCACOX
    Precp_CaSO4_conc=CCASOX
    Precp_AlPO4_mole_conc=CALPOX*VLPO4
    Precp_FePO4_mole_conc=CFEPOX*VLPO4
    Precp_CaHPO4_mole_conc=CCAPDX*VLPO4
    Precp_Ca5P3O12O3H3_mole_conc=CCAPHX*VLPO4
    CEC_conc=AMAX1(ZERO,XCEC/BulkSoilMass)
    CALX=Al_3p_aque_mole_conc**0.333_r8
    CFEX=Fe_3p_aque_mole_conc**0.333_r8
    CaX_conc=Ca_2p_aque_mole_conc**0.500
    MgX_conc=Mg_2p_aque_mole_conc**0.500
    XCAX=CEC_conc/(1.0_r8+GKC4*NH4_1p_aque_mole_conc/CaX_conc &
      +GKCH*H_1p_aque_mole_conc/CaX_conc+GKCA*CALX/CaX_conc &
      +GKCA*CFEX/CaX_conc+GKCA*Mg_2p_aque_mole_conc/CaX_conc &
      +GKCN*Na_1p_aque_mole_conc/CaX_conc+GKCK*K_1p_aque_mole_conc/CaX_conc)
    XN4Q=CN4X
    XHYQ=XCAX*GKCH
    XALQ=XCAX*GKCA
    XFEQ=XCAX*GKCA
    XCAQ=XCAX*CaX_conc
    XMGQ=XCAX*GKCM
    XNAQ=XCAX*GKCN
    XKAQ=XCAX*GKCK
    XTLQ=XN4Q+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ
    IF(XTLQ.GT.ZERO)THEN
      FX=CEC_conc/XTLQ
    ELSE
      FX=0._r8
    ENDIF
    XNH4_mole_conc=CN4X
    XHY1=FX*XHYQ
    XAl_conc=FX*XALQ/3.0_r8
    XFe_conc=FX*XFEQ/3.0_r8
    XCa_conc=FX*XCAQ/2.0_r8
    XMg_conc=FX*XMGQ/2.0_r8
    XNa_conc=FX*XNAQ
    XK_conc=FX*XKAQ
  ENDIF
  end subroutine InitEquilibria
!------------------------------------------------------------------------------------------

  subroutine SolubilityEquilibiriaSalt(K,M,BulkSoilMass)

  implicit none
  integer, intent(in) :: K,M
  real(r8), intent(in) :: BulkSoilMass  !Mg/d-2, soil mass whole layer
  real(r8) :: XN4Q,XMGQ,XKAQ,XHYQ,XFEQ,XTLQ
  real(r8) :: XALQ,XCAX,SPH2P,XCAQ,SPH1P,FX
  real(r8) :: RN4S,RN3S,RNA,RNAC,RNH4
  real(r8) :: RPALOX,H2PO4_1e_AlPO4_dissol_flx,RPCACX,H2PO4_1e_CaHPO4_dissol_flx
  real(r8) :: H2PO4_e_to_HPO4_2e_flx,RHP1,RHP2,H2PO4_1e_apatite_dissol_flx
  real(r8) :: H2PO4_1e_FePO4_dissol_flx,H1PO4_to_XHPO4_ROH_flx
  real(r8) :: FSTR2,RCAO,FHCO,FHP2
  real(r8) :: H2PO4_1e_to_XH2PO4_ROH2_flx,RXN4,RF1P,FCO3,SP
  real(r8) :: XCOO,XNAQ,RHP0,RHP3,RKAS,RM1P
  real(r8) :: CX1,CX2,CSTRZ,CCO2Y,CION2
  real(r8) :: RFE1,RFE2,RFE3,RFE4
  real(r8) :: RFEO1,RFEO2,RFEO3,RFEO4
  real(r8) :: RHA4P1,RHA4P2,RHAL1
  real(r8) :: RHALO1,RHALO2,RHALO3,RHALO4
  real(r8) :: H2PO4_to_XH2PO4_ROH_flx,S0,S1,VLWatMicPBK,RCA,RCAC,RCAH
  real(r8) :: RCAS,RF2P,RFE,RFES,RH1P,RH3P
  real(r8) :: RHA0P1,RHA0P2,RHA1P1,RHA1P2,RHA2P1
  real(r8) :: RHA2P2,RHA3P1,RHCAC3,RHCACH,RHA3P2
  real(r8) :: RHCACO,RHCAH1,RHCAH2,RHF0P1,RHCAD2
  real(r8) :: RHF1P2,RHF2P2,RHF3P1,RHF0P2,RHF3P2
  real(r8) :: RHF1P1,RHF2P1,RHF4P1,RHF4P2,RHFE1
  real(r8) :: RHFEO1,RHFEO2,RHFEO3,RHFEO4,RKA,RMG
  real(r8) :: RMGC,RMGH,RMGO,RMGS,RNAS,RPCAD1,RPCASO
  real(r8) :: RPFEOX,RSO4,RXCA,RXFE,RXHY,RXKA,RXMG
  real(r8) :: RXAL,RXNA,RXOH1,RXOH2,SPX,R,Z
  real(r8) :: A3,Al_3p_activity,AlOH_2p_activity,AlO2H2_1p_activity,AlO3H3_activity,AlO4H4_1e_activity
  real(r8) :: AlSO4_1p_activity,AALX,CaPO4_1e_activity,CaHPO4_activity,CaH4P2O8_1p_activity,Ca_2p_activity,CaCO3_activity
  real(r8) :: CaHCO3_1p_activity,CaSO4_activity,CaO2H2_activity,ACAX,H2CO3_activity,CO3_2e_activity
  real(r8) :: FeHPO4_p_activity,FeH2PO4_2p_activity,Fe_3p_activity,FeOH_2p_activity
  real(r8) :: FeO2H2_p_activity,FeO3H3_activity,FeO4H4_1e_activity
  real(r8) :: FeSO4_1p_activity,AFEX,H0PO4_3e_activity,H1PO4_2e_activity,H2PO4_1e_activity,H3PO4_activity,HCO3_e_activity
  real(r8) :: H_1p_activity,K_1p_activity,KSO4_1e_activity,MgHPO4_activity,Mg_2p_activity,MgCO3_activity,MgOH_1p_activity
  real(r8) :: AMGX,MgHCO3_1p_activity,MgSO4_activity,NH3_activity,NH4_1p_activity,Na_1p_activity,NaCO3_1e_activity
  real(r8) :: OH_1e_activity,SO4_2e_activity,NaSO4_1e_activity,anion_1e_aque_mole_conc,anion_2e_aque_mole_conc
  real(r8) :: anion_3e_conc,cation_1p_aque_mole_conc,cation_2p_aque_mole_conc
  real(r8) :: cation_3p_aque_mole_conc,CCO2X,CCO2Z,CN,CSTR2,DP,P1,P2,P3,PX
  real(r8) :: PY,R1,RAL,RAL1,RAL2,RAL3,RAL4,RALO1
  real(r8) :: RALO2,RALO3,RALO4,RALS,RC0P,RC1P,RC2P
  integer :: NR1,NP2,NP3
  integer :: MM

! begin_execution
  H2CO3_aque_mole_conc       = AMAX1(ZERO,H2CO3_aque_mole_conc)
  CO3_2e_aque_mole_conc      = H2CO3_aque_mole_conc*DPCO3*A0/(H_1p_aque_mole_conc**2*A2)
  HCO3_e_conc                = H2CO3_aque_mole_conc*DPCO2*A0/(H_1p_aque_mole_conc*A1)
  NH4_1p_aque_mole_conc      = AMAX1(ZERO,NH4_1p_aque_mole_conc)
  NH3_aque_mole_conc         = AMAX1(ZERO,NH3_aque_mole_conc)
  Al_3p_aque_mole_conc       = AMAX1(ZERO,Al_3p_aque_mole_conc)
  Fe_3p_aque_mole_conc       = AMAX1(ZERO,Fe_3p_aque_mole_conc)
  Ca_2p_aque_mole_conc       = AMAX1(ZERO,Ca_2p_aque_mole_conc)
  Ca_2p_aque_mole_conc       = AMIN1(Ca_2p_aque_mole_conc,SPCAC/(CO3_2e_aque_mole_conc*A2**2))
  Mg_2p_aque_mole_conc       = AMAX1(ZERO,Mg_2p_aque_mole_conc)
  Na_1p_aque_mole_conc       = AMAX1(ZERO,Na_1p_aque_mole_conc)
  K_1p_aque_mole_conc        = AMAX1(ZERO,K_1p_aque_mole_conc)
  SO4_2e_aque_mole_conc      = AMAX1(ZERO,SO4_2e_aque_mole_conc)
  AlOH_2p_aque_mole_conc     = AMAX1(ZERO,AlOH_2p_aque_mole_conc)
  AlO2H2_1p_aque_mole_conc   = AMAX1(ZERO,AlO2H2_1p_aque_mole_conc)
  AlO3H3_conc                = AMAX1(ZERO,AlO3H3_conc)
  AlO4H4_1e_aque_mole_conc   = AMAX1(ZERO,AlO4H4_1e_aque_mole_conc)
  AlSO4_1p_aque_mole_conc    = AMAX1(ZERO,AlSO4_1p_aque_mole_conc)
  FeOH_2p_aque_mole_conc     = AMAX1(ZERO,FeOH_2p_aque_mole_conc)
  FeO2H2_p_conc              = AMAX1(ZERO,FeO2H2_p_conc)
  FeO3H3_conc                = AMAX1(ZERO,FeO3H3_conc)
  FeO4H4_1e_aque_mole_conc   = AMAX1(ZERO,FeO4H4_1e_aque_mole_conc)
  FeSO4_1p_aque_mole_conc    = AMAX1(ZERO,FeSO4_1p_aque_mole_conc)
  CaO2H2_conc                = AMAX1(ZERO,CaO2H2_conc)
  CaCO3_conc                 = AMAX1(ZERO,CaCO3_conc)
  CaHCO3_1p_aque_mole_conc   = AMAX1(ZERO,CaHCO3_1p_aque_mole_conc)
  CaSO4_conc                 = AMAX1(ZERO,CaSO4_conc)
  MgOH_1p_aque_mole_conc     = AMAX1(ZERO,MgOH_1p_aque_mole_conc)
  MgCO3_conc                 = AMAX1(ZERO,MgCO3_conc)
  MgHCO3_1p_aque_mole_conc   = AMAX1(ZERO,MgHCO3_1p_aque_mole_conc)
  MgSO4_conc                 = AMAX1(ZERO,MgSO4_conc)
  NaCO3_1e_aque_mole_conc    = AMAX1(ZERO,NaCO3_1e_aque_mole_conc)
  NaSO4_1e_aque_mole_conc    = AMAX1(ZERO,NaSO4_1e_aque_mole_conc)
  KSO4_1e_aque_mole_conc     = AMAX1(ZERO,KSO4_1e_aque_mole_conc)
  H0PO4_3e_conc              = AMAX1(ZERO,H0PO4_3e_conc)
  H1PO4_2e_aque_mole_conc    = AMAX1(ZERO,H1PO4_2e_aque_mole_conc)
  H2PO4_1e_aque_mole_conc    = AMAX1(ZERO,H2PO4_1e_aque_mole_conc)
  H3PO4_conc                 = AMAX1(ZERO,H3PO4_conc)
  FeHPO4_p_conc              = AMAX1(ZERO,FeHPO4_p_conc)
  FeH2PO4_2p_aque_mole_conc  = AMAX1(ZERO,FeH2PO4_2p_aque_mole_conc)
  CaPO4_1e_con               = AMAX1(ZERO,CaPO4_1e_con)
  CaHPO4_conc                = AMAX1(ZERO,CaHPO4_conc)
  CaH4P2O8_1p_aque_mole_conc = AMAX1(ZERO,CaH4P2O8_1p_aque_mole_conc)
  MgHPO4_conc                = AMAX1(ZERO,MgHPO4_conc)
!
!     ION ACTIVITY COEFFICIENTS
!
  cation_3p_aque_mole_conc=Al_3p_aque_mole_conc+Fe_3p_aque_mole_conc
  anion_3e_conc=H0PO4_3e_conc
  cation_2p_aque_mole_conc=Ca_2p_aque_mole_conc+Mg_2p_aque_mole_conc+AlOH_2p_aque_mole_conc+FeOH_2p_aque_mole_conc+FeH2PO4_2p_aque_mole_conc
  anion_2e_aque_mole_conc=SO4_2e_aque_mole_conc+CO3_2e_aque_mole_conc+H1PO4_2e_aque_mole_conc
  cation_1p_aque_mole_conc=NH4_1p_aque_mole_conc+H_1p_aque_mole_conc+Na_1p_aque_mole_conc+K_1p_aque_mole_conc+AlO2H2_1p_aque_mole_conc &
    +FeO2H2_p_conc+AlSO4_1p_aque_mole_conc+FeSO4_1p_aque_mole_conc+CaO2H2_conc &
    +CaHCO3_1p_aque_mole_conc+MgOH_1p_aque_mole_conc+MgHCO3_1p_aque_mole_conc+FeHPO4_p_conc+CaH4P2O8_1p_aque_mole_conc
  anion_1e_aque_mole_conc=NO3_1e_aque_mole_conc+OH_1e_aque_mole_conc+HCO3_e_conc+Cl_e_conc+AlO4H4_1e_aque_mole_conc &
    +FeO4H4_1e_aque_mole_conc+NaCO3_1e_aque_mole_conc+NaSO4_1e_aque_mole_conc+KSO4_1e_aque_mole_conc &
    +H2PO4_1e_aque_mole_conc+CaPO4_1e_con
  CN=H2CO3_aque_mole_conc+CH4_aque_mole_conc+O2_aque_mole_conc+N2_aque_mole_conc+N2O_aque_mole_conc+NH3_aque_mole_conc &
    +AlO3H3_conc+FeO3H3_conc+CaCO3_conc+CaSO4_conc &
    +MgCO3_conc+MgSO4_conc+H3PO4_conc+CaHPO4_conc+MgHPO4_conc
  CX2=anion_2e_aque_mole_conc-CO3_2e_aque_mole_conc
  CX1=anion_1e_aque_mole_conc-HCO3_e_conc
  CION2=AZMAX1(cation_3p_aque_mole_conc+anion_3e_conc+cation_2p_aque_mole_conc+anion_2e_aque_mole_conc+cation_1p_aque_mole_conc+anion_1e_aque_mole_conc+CN)
  CSTR1=0.5E-03_r8*(9.0_r8*(cation_3p_aque_mole_conc+anion_3e_conc)+4.0*(cation_2p_aque_mole_conc+anion_2e_aque_mole_conc)+cation_1p_aque_mole_conc+anion_1e_aque_mole_conc)
  CSTRZ=0.5E-03_r8*(9.0_r8*(cation_3p_aque_mole_conc+anion_3e_conc)+4.0*(cation_2p_aque_mole_conc+CX2)+cation_1p_aque_mole_conc+CX1)
  CSTR2=SQRT(CSTR1)
  FSTR2=CSTR2/(1.0_r8+CSTR2)
  FH2O=5.56E+04_r8/(5.56E+04_r8+CION2)
  A1=AMIN1(1.0_r8,10.0_r8**(-0.509_r8*1.0_r8*FSTR2+0.20_r8*CSTR2))
  A2=AMIN1(1.0_r8,10.0_r8**(-0.509_r8*4.0_r8*FSTR2+0.20_r8*CSTR2))
  A3=AMIN1(1.0_r8,10.0_r8**(-0.509_r8*9.0_r8*FSTR2+0.20_r8*CSTR2))
!
!     PRECIPITATION-DISSOLUTION EQUILIBRIA
!
  ! Refresh water dissociation after each ionic-strength/activity update.
  OH_1e_aque_mole_conc=DPH2O/(A1*A1*H_1p_aque_mole_conc)
  H_1p_activity        = H_1p_aque_mole_conc*A1
  OH_1e_activity       = OH_1e_aque_mole_conc*A1
  Al_3p_activity       = Al_3p_aque_mole_conc*A3
  AlOH_2p_activity     = AlOH_2p_aque_mole_conc*A2
  AlO2H2_1p_activity   = AlO2H2_1p_aque_mole_conc*A1
  AlO3H3_activity      = AlO3H3_conc*A0
  AlO4H4_1e_activity   = AlO4H4_1e_aque_mole_conc*A1
  Fe_3p_activity       = Fe_3p_aque_mole_conc*A3
  FeOH_2p_activity     = FeOH_2p_aque_mole_conc*A2
  FeO2H2_p_activity    = FeO2H2_p_conc*A1
  FeO3H3_activity      = FeO3H3_conc*A0
  FeO4H4_1e_activity   = FeO4H4_1e_aque_mole_conc*A1
  Ca_2p_activity       = Ca_2p_aque_mole_conc*A2
  CO3_2e_activity      = CO3_2e_aque_mole_conc*A2
  HCO3_e_activity      = HCO3_e_conc*A1
  H2CO3_activity       = H2CO3_aque_mole_conc*A0
  SO4_2e_activity      = SO4_2e_aque_mole_conc*A2
  H0PO4_3e_activity    = H0PO4_3e_conc*A3
  H1PO4_2e_activity    = H1PO4_2e_aque_mole_conc*A2
  H2PO4_1e_activity    = H2PO4_1e_aque_mole_conc*A1
  H3PO4_activity       = H3PO4_conc*A0
  FeHPO4_p_activity    = FeHPO4_p_conc*A1
  FeH2PO4_2p_activity  = FeH2PO4_2p_aque_mole_conc*A2
  CaPO4_1e_activity    = CaPO4_1e_con*A1
  CaHPO4_activity      = CaHPO4_conc*A0
  CaH4P2O8_1p_activity = CaH4P2O8_1p_aque_mole_conc*A1
  MgHPO4_activity      = MgHPO4_conc*A0
  NH4_1p_activity      = NH4_1p_aque_mole_conc*A1
  NH3_activity         = NH3_aque_mole_conc*A0
  Mg_2p_activity       = Mg_2p_aque_mole_conc*A2
  Na_1p_activity       = Na_1p_aque_mole_conc*A1
  K_1p_activity        = K_1p_aque_mole_conc*A1
  AALX                 = Al_3p_activity**0.333_r8
  AFEX                 = Fe_3p_activity**0.333_r8
  ACAX                 = Ca_2p_activity**0.500_r8
  AMGX                 = Mg_2p_activity**0.500_r8
  AlSO4_1p_activity    = AlSO4_1p_aque_mole_conc*A1
  FeSO4_1p_activity    = FeSO4_1p_aque_mole_conc*A1
  CaO2H2_activity      = CaO2H2_conc*A1
  CaCO3_activity       = CaCO3_conc*A0
  CaSO4_activity       = CaSO4_conc*A0
  CaHCO3_1p_activity   = CaHCO3_1p_aque_mole_conc*A1
  MgOH_1p_activity     = MgOH_1p_aque_mole_conc*A1
  MgCO3_activity       = MgCO3_conc*A0
  MgHCO3_1p_activity   = MgHCO3_1p_aque_mole_conc*A1
  MgSO4_activity       = MgSO4_conc*A0
  NaCO3_1e_activity    = NaCO3_1e_aque_mole_conc*A1
  NaSO4_1e_activity    = NaSO4_1e_aque_mole_conc*A1
  KSO4_1e_activity     = KSO4_1e_aque_mole_conc*A1
!
!     ALUMINUM HYDROXIDE (GIBBSITE)
!
!     Al(OH)3(s) + 3H(+) <=> Al(3+) + 3H2O.
!     The dominant dissolved metal hydrolysis form selects an equivalent pH-dependent reaction.
!     Positive mineral flux precipitates the solid; negative flux dissolves it.
  IF(K.EQ.micpar%k_POM)THEN
    ! PX/PY copy an AMAX1 operand, so exact equality identifies the selected
    ! species even at trace activity; the first exact tie keeps precedence.
    PX=AMAX1(Al_3p_activity,AlOH_2p_activity,AlO2H2_1p_activity,AlO3H3_activity,AlO4H4_1e_activity)
    IF(PX.EQ.Al_3p_activity)THEN
      R1=H_1p_activity
      P1=Al_3p_activity
      P2=OH_1e_activity
      NR1=3
      NP2=0
      SP=SHALO
    ELSEIF(PX.EQ.AlOH_2p_activity)THEN
      R1=H_1p_activity
      P1=AlOH_2p_activity
      P2=OH_1e_activity
      NR1=2
      NP2=0
      SP=SHAL1
    ELSEIF(PX.EQ.AlO2H2_1p_activity)THEN
      R1=H_1p_activity
      P1=AlO2H2_1p_activity
      P2=OH_1e_activity
      NR1=1
      NP2=0
      SP=SHAL2
    ELSEIF(PX.EQ.AlO3H3_activity)THEN
      R1=H_1p_activity
      P1=AlO3H3_activity
      P2=OH_1e_activity
      NR1=0
      NP2=0
      SP=SPAL3
    ELSEIF(PX.EQ.AlO4H4_1e_activity)THEN
      R1=OH_1e_activity
      P1=AlO4H4_1e_activity
      P2=H_1p_activity
      NR1=0
      NP2=1
      SP=SHAL4
    ENDIF
    RHAL1=0._r8
    RHALO1=0._r8
    RHALO2=0._r8
    RHALO3=0._r8
    RHALO4=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1/P2**NP2
    RPALOX=AMAX1(-Precp_AlO3H3_conc,TPD*(P1-SPX))
    IF(PX.EQ.Al_3p_activity)THEN
      RHAL1=RPALOX
    ELSEIF(PX.EQ.AlOH_2p_activity)THEN
      RHALO1=RPALOX
    ELSEIF(PX.EQ.AlO2H2_1p_activity)THEN
      RHALO2=RPALOX
    ELSEIF(PX.EQ.AlO3H3_activity)THEN
      RHALO3=RPALOX
    ELSEIF(PX.EQ.AlO4H4_1e_activity)THEN
      RHALO4=RPALOX
    ENDIF
!
!     IRON HYDROXIDE
!
!     Fe(OH)3(s) + 3H(+) <=> Fe(3+) + 3H2O.
!     The dominant dissolved metal hydrolysis form selects an equivalent pH-dependent reaction.
!     Positive mineral flux precipitates the solid; negative flux dissolves it.
    PX=AMAX1(Fe_3p_activity,FeOH_2p_activity,FeO2H2_p_activity,FeO3H3_activity,FeO4H4_1e_activity)
    IF(PX.EQ.Fe_3p_activity)THEN
      R1=H_1p_activity
      P1=Fe_3p_activity
      P2=OH_1e_activity
      NR1=3
      NP2=0
      SP=SHFEO
    ELSEIF(PX.EQ.FeOH_2p_activity)THEN
      R1=H_1p_activity
      P1=FeOH_2p_activity
      P2=OH_1e_activity
      NR1=2
      NP2=0
      SP=SHFE1
    ELSEIF(PX.EQ.FeO2H2_p_activity)THEN
      R1=H_1p_activity
      P1=FeO2H2_p_activity
      P2=OH_1e_activity
      NR1=1
      NP2=0
      SP=SHFE2
    ELSEIF(PX.EQ.FeO3H3_activity)THEN
      R1=H_1p_activity
      P1=FeO3H3_activity
      P2=OH_1e_activity
      NR1=0
      NP2=0
      SP=SPFE3
    ELSEIF(PX.EQ.FeO4H4_1e_activity)THEN
      R1=OH_1e_activity
      P1=FeO4H4_1e_activity
      P2=H_1p_activity
      NR1=0
      NP2=1
      SP=SHFE4
    ENDIF
    RHFE1=0._r8
    RHFEO1=0._r8
    RHFEO2=0._r8
    RHFEO3=0._r8
    RHFEO4=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1/P2**NP2
    RPFEOX=AMAX1(-Precp_FeO3H3_conc,TPD*(P1-SPX))
    IF(PX.EQ.Fe_3p_activity)THEN
      RHFE1=RPFEOX
    ELSEIF(PX.EQ.FeOH_2p_activity)THEN
      RHFEO1=RPFEOX
    ELSEIF(PX.EQ.FeO2H2_p_activity)THEN
      RHFEO2=RPFEOX
    ELSEIF(PX.EQ.FeO3H3_activity)THEN
      RHFEO3=RPFEOX
    ELSEIF(PX.EQ.FeO4H4_1e_activity)THEN
      RHFEO4=RPFEOX
    ENDIF
!
!     CALCITE
!
!     CaCO3(s) <=> Ca(2+) + CO3(2-); with HCO3(-), consume one H(+),
!     or with H2CO3(aq), consume two H(+).
!     Positive mineral flux precipitates the solid; negative flux dissolves it.
    PX=AMAX1(CO3_2e_activity,HCO3_e_activity,H2CO3_activity)
    R1=H_1p_activity
    P1=Ca_2p_activity
    IF(PX.EQ.CO3_2e_activity)THEN
      P2=CO3_2e_activity
      NR1=0
      SP=SPCAC
    ELSEIF(PX.EQ.HCO3_e_activity)THEN
      P2=HCO3_e_activity
      NR1=1
      SP=SHCAC1
    ELSEIF(PX.EQ.H2CO3_activity)THEN
      P2=H2CO3_activity
      NR1=2
      SP=SHCAC2
    ENDIF
    RHCAC3=0._r8
    RHCACH=0._r8
    RHCACO=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1
    S0=P1+P2
    ! Keep small positive discriminants: concentration cutoffs do not apply to squares.
    ! Mineral D=(P1-P2)**2+4*SPX; rationalize S0-sqrt(D) to retain small fluxes.
    S1=MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    RPCACX=AMAX1(-Precp_CaCO3_conc,TPD*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))
    IF(PX.EQ.CO3_2e_activity)THEN
      RHCAC3=RPCACX
    ELSEIF(PX.EQ.HCO3_e_activity)THEN
      RHCACH=RPCACX
    ELSEIF(PX.EQ.H2CO3_activity)THEN
      RHCACO=RPCACX
    ENDIF
!
!     GYPSUM
!
!     CaSO4(s) <=> Ca(2+) + SO4(2-).
!     Positive mineral flux precipitates the solid; negative flux dissolves it.
    P1=Ca_2p_activity
    P2=SO4_2e_activity
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SPCAS
    S0=P1+P2
    S1=MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    RPCASO=AMAX1(-Precp_CaSO4_conc,TPD*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))
!
!     PHOSPHORUS PRECIPITATION-DISSOLUTION IN NON-BAND SOIL ZONE
!
!
!     ALUMINUM PHOSPHATE (VARISCITE)
!
!     AlPO4(s) + 2H(+) <=> Al(3+) + H2PO4(-),
!     or AlPO4(s) + H(+) <=> Al(3+) + HPO4(2-).
!     The dominant dissolved metal hydrolysis form selects an equivalent pH-dependent reaction.
!     Positive mineral flux precipitates the solid; negative flux dissolves it.
    PX=AMAX1(Al_3p_activity,AlOH_2p_activity,AlO2H2_1p_activity,AlO3H3_activity,AlO4H4_1e_activity)
    PY=AMAX1(H1PO4_2e_activity,H2PO4_1e_activity)
    R1=H_1p_activity
    P3=H_1p_activity
    IF(PY.EQ.H1PO4_2e_activity)THEN
      P2=H1PO4_2e_activity
      IF(PX.EQ.Al_3p_activity)THEN
        P1=Al_3p_activity
        NR1=1
        NP3=0
        SP=SHA0P1
      ELSEIF(PX.EQ.AlOH_2p_activity)THEN
        P1=AlOH_2p_activity
        NR1=0
        NP3=0
        SP=SPA1P1
      ELSEIF(PX.EQ.AlO2H2_1p_activity)THEN
        P1=AlO2H2_1p_activity
        NR1=0
        NP3=1
        SP=SHA2P1
      ELSEIF(PX.EQ.AlO3H3_activity)THEN
        P1=AlO3H3_activity
        NR1=0
        NP3=2
        SP=SHA3P1
      ELSEIF(PX.EQ.AlO4H4_1e_activity)THEN
        P1=AlO4H4_1e_activity
        NR1=0
        NP3=3
        SP=SHA4P1
      ENDIF
    ELSE
      P2=H2PO4_1e_activity
      IF(PX.EQ.Al_3p_activity)THEN
        P1=Al_3p_activity
        NR1=2
        NP3=0
        SP=SHA0P2
      ELSEIF(PX.EQ.AlOH_2p_activity)THEN
        P1=AlOH_2p_activity
        NR1=1
        NP3=0
        SP=SHA1P2
      ELSEIF(PX.EQ.AlO2H2_1p_activity)THEN
        P1=AlO2H2_1p_activity
        NR1=0
        NP3=0
        SP=SPA2P2
      ELSEIF(PX.EQ.AlO3H3_activity)THEN
        P1=AlO3H3_activity
        NR1=0
        NP3=1
        SP=SHA3P2
      ELSEIF(PX.EQ.AlO4H4_1e_activity)THEN
        P1=AlO4H4_1e_activity
        NR1=0
        NP3=2
        SP=SHA4P2
      ENDIF
    ENDIF
    RHA0P1=0._r8
    RHA1P1=0._r8
    RHA2P1=0._r8
    RHA3P1=0._r8
    RHA4P1=0._r8
    RHA0P2=0._r8
    RHA1P2=0._r8
    RHA2P2=0._r8
    RHA3P2=0._r8
    RHA4P2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    P3=AMAX1(ZERO,P3)
    SPX=SP*R1**NR1/P3**NP3
    S0=P1+P2
    S1=MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    H2PO4_1e_AlPO4_dissol_flx=AMAX1(-Precp_AlPO4_mole_conc,TPD*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))
    IF(PY.EQ.H1PO4_2e_activity)THEN
      IF(PX.EQ.Al_3p_activity)THEN
        RHA0P1=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(PX.EQ.AlOH_2p_activity)THEN
        RHA1P1=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(PX.EQ.AlO2H2_1p_activity)THEN
        RHA2P1=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(PX.EQ.AlO3H3_activity)THEN
        RHA3P1=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(PX.EQ.AlO4H4_1e_activity)THEN
        RHA4P1=H2PO4_1e_AlPO4_dissol_flx
      ENDIF
    ELSE
      IF(PX.EQ.Al_3p_activity)THEN
        RHA0P2=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(PX.EQ.AlOH_2p_activity)THEN
        RHA1P2=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(PX.EQ.AlO2H2_1p_activity)THEN
        RHA2P2=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(PX.EQ.AlO3H3_activity)THEN
        RHA3P2=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(PX.EQ.AlO4H4_1e_activity)THEN
        RHA4P2=H2PO4_1e_AlPO4_dissol_flx
      ENDIF
    ENDIF
!
!     IRON PHOSPHATE (STRENGITE)
!
!     FePO4(s) + 2H(+) <=> Fe(3+) + H2PO4(-),
!     or FePO4(s) + H(+) <=> Fe(3+) + HPO4(2-).
!     The dominant dissolved metal hydrolysis form selects an equivalent pH-dependent reaction.
!     Positive mineral flux precipitates the solid; negative flux dissolves it.
    PX=AMAX1(Fe_3p_activity,FeOH_2p_activity,FeO2H2_p_activity,FeO3H3_activity,FeO4H4_1e_activity)
    PY=AMAX1(H1PO4_2e_activity,H2PO4_1e_activity)
    R1=H_1p_activity
    P3=H_1p_activity
    IF(PY.EQ.H1PO4_2e_activity)THEN
      P2=H1PO4_2e_activity
      IF(PX.EQ.Fe_3p_activity)THEN
        P1=Fe_3p_activity
        NR1=1
        NP3=0
        SP=SHF0P1
      ELSEIF(PX.EQ.FeOH_2p_activity)THEN
        P1=FeOH_2p_activity
        NR1=0
        NP3=0
        SP=SPF1P1
      ELSEIF(PX.EQ.FeO2H2_p_activity)THEN
        P1=FeO2H2_p_activity
        NR1=0
        NP3=1
        SP=SHF2P1
      ELSEIF(PX.EQ.FeO3H3_activity)THEN
        P1=FeO3H3_activity
        NR1=0
        NP3=2
        SP=SHF3P1
      ELSEIF(PX.EQ.FeO4H4_1e_activity)THEN
        P1=FeO4H4_1e_activity
        NR1=0
        NP3=3
        SP=SHF4P1
      ENDIF
    ELSE
      P2=H2PO4_1e_activity
      IF(PX.EQ.Fe_3p_activity)THEN
        P1=Fe_3p_activity
        NR1=2
        NP3=0
        SP=SHF0P2
      ELSEIF(PX.EQ.FeOH_2p_activity)THEN
        P1=FeOH_2p_activity
        NR1=1
        NP3=0
        SP=SHF1P2
      ELSEIF(PX.EQ.FeO2H2_p_activity)THEN
        P1=FeO2H2_p_activity
        NR1=0
        NP3=0
        SP=SPF2P2
      ELSEIF(PX.EQ.FeO3H3_activity)THEN
        P1=FeO3H3_activity
        NR1=0
        NP3=1
        SP=SHF3P2
      ELSEIF(PX.EQ.FeO4H4_1e_activity)THEN
        P1=FeO4H4_1e_activity
        NR1=0
        NP3=2
        SP=SHF4P2
      ENDIF
    ENDIF
    RHF0P1=0._r8
    RHF1P1=0._r8
    RHF2P1=0._r8
    RHF3P1=0._r8
    RHF4P1=0._r8
    RHF0P2=0._r8
    RHF1P2=0._r8
    RHF2P2=0._r8
    RHF3P2=0._r8
    RHF4P2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    P3=AMAX1(ZERO,P3)
    SPX=SP*R1**NR1/P3**NP3
    S0=P1+P2
    S1=MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    H2PO4_1e_FePO4_dissol_flx=AMAX1(-Precp_FePO4_mole_conc,TPD*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))
    IF(PY.EQ.H1PO4_2e_activity)THEN
      IF(PX.EQ.Fe_3p_activity)THEN
        RHF0P1=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(PX.EQ.FeOH_2p_activity)THEN
        RHF1P1=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(PX.EQ.FeO2H2_p_activity)THEN
        RHF2P1=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(PX.EQ.FeO3H3_activity)THEN
        RHF3P1=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(PX.EQ.FeO4H4_1e_activity)THEN
        RHF4P1=H2PO4_1e_FePO4_dissol_flx
      ENDIF
    ELSE
      IF(PX.EQ.Fe_3p_activity)THEN
        RHF0P2=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(PX.EQ.FeOH_2p_activity)THEN
        RHF1P2=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(PX.EQ.FeO2H2_p_activity)THEN
        RHF2P2=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(PX.EQ.FeO3H3_activity)THEN
        RHF3P2=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(PX.EQ.FeO4H4_1e_activity)THEN
        RHF4P2=H2PO4_1e_FePO4_dissol_flx
      ENDIF
    ENDIF
!
!     DICALCIUM PHOSPHATE
!
!     CaHPO4(s) <=> Ca(2+) + HPO4(2-),
!     or CaHPO4(s) + H(+) <=> Ca(2+) + H2PO4(-).
!     Positive mineral flux precipitates the solid; negative flux dissolves it.
    PX=AMAX1(H1PO4_2e_activity,H2PO4_1e_activity)
    R1=H_1p_activity
    P1=Ca_2p_activity
    IF(PX.EQ.H1PO4_2e_activity)THEN
      P2=H1PO4_2e_activity
      NR1=0
      SP=SPCAD
    ELSEIF(PX.EQ.H2PO4_1e_activity)THEN
      P2=H2PO4_1e_activity
      NR1=1
      SP=SHCAD2
    ENDIF
    RPCAD1=0._r8
    RHCAD2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1
    S0=P1+P2
    S1=MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    H2PO4_1e_CaHPO4_dissol_flx=AMAX1(-Precp_CaHPO4_mole_conc,TPD*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))
    IF(PX.EQ.H1PO4_2e_activity)THEN
      RPCAD1=H2PO4_1e_CaHPO4_dissol_flx
    ELSEIF(PX.EQ.H2PO4_1e_activity)THEN
      RHCAD2=H2PO4_1e_CaHPO4_dissol_flx
    ENDIF
!
!     HYDROXYAPATITE
!
!     Ca5(PO4)3OH(s) + 4H(+) <=> 5Ca(2+) + 3HPO4(2-) + H2O,
!     or with 7H(+), release 3H2PO4(-) instead of 3HPO4(2-).
!     Positive mineral flux precipitates the solid; negative flux dissolves it.
    PX=AMAX1(H1PO4_2e_activity,H2PO4_1e_activity)
    R1=H_1p_activity
    P1=Ca_2p_activity
    IF(PX.EQ.H1PO4_2e_activity)THEN
      P2=H1PO4_2e_activity
      NR1=4
      SP=SHCAH1
    ELSEIF(PX.EQ.H2PO4_1e_activity)THEN
      P2=H2PO4_1e_activity
      NR1=7
      SP=SHCAH2
    ENDIF
    RHCAH1=0._r8
    RHCAH2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=(SP*R1**NR1/P1**5_r8)**0.333_r8
    H2PO4_1e_apatite_dissol_flx=AMAX1(-Precp_Ca5P3O12O3H3_mole_conc,TPD*(P2-SPX))
    IF(PX.EQ.H1PO4_2e_activity)THEN
      RHCAH1=H2PO4_1e_apatite_dissol_flx
    ELSEIF(PX.EQ.H2PO4_1e_activity)THEN
      RHCAH2=H2PO4_1e_apatite_dissol_flx
    ENDIF
    Precp_AlO3H3_conc            = Precp_AlO3H3_conc+RPALOX
    Precp_FeO3H3_conc            = Precp_FeO3H3_conc+RPFEOX
    Precp_CaCO3_conc             = Precp_CaCO3_conc+RPCACX
    Precp_CaSO4_conc             = Precp_CaSO4_conc+RPCASO
    Precp_AlPO4_mole_conc        = Precp_AlPO4_mole_conc+H2PO4_1e_AlPO4_dissol_flx
    Precp_FePO4_mole_conc        = Precp_FePO4_mole_conc+H2PO4_1e_FePO4_dissol_flx
    Precp_CaHPO4_mole_conc       = Precp_CaHPO4_mole_conc+H2PO4_1e_CaHPO4_dissol_flx
    Precp_Ca5P3O12O3H3_mole_conc = Precp_Ca5P3O12O3H3_mole_conc+H2PO4_1e_apatite_dissol_flx
!
!     ANION EXCHANGE EQILIBRIA
!
    IF(VLWatMicP.GT.ZEROS)THEN
      VLWatMicPBK=AMIN1(1.0_r8,BulkSoilMass/VLWatMicP)
    ELSE
      VLWatMicPBK=1.0_r8
    ENDIF
    IF(XAEC.GT.ZEROS)THEN
      ! R-OH + H(+) <=> R-OH2(+); positive RXOH2 protonates the site.
      RXOH2=TAD*(XROH1_mole_conc*H_1p_activity-SXOH2*XROH2_mole_conc)/(XROH1_mole_conc+SPOH2)*VLWatMicPBK
      ! R-O(-) + H(+) <=> R-OH; positive RXOH1 protonates the site.
      RXOH1=TAD*(XOH_mole_conc*H_1p_activity-SXOH1*XROH1_mole_conc)/(XOH_mole_conc+SPOH1)*VLWatMicPBK
      ! Intended ligand exchange: R-OH2(+) + H2PO4(-) <=> R-H2PO4 + H2O.
      ! Positive flux binds dissolved H2PO4 to a protonated surface site.
      ! R-H2PO4 + H2O <=> R-OH2(+) + H2PO4(-).
      ! Combine neutral-site exchange (SXH2P), water dissociation (DPH2O),
      ! and R-OH2(+) <=> R-OH + H(+) (SXOH2) to close the equilibrium cycle.
      SPH2P = SXH2P*DPH2O/SXOH2
      H2PO4_1e_to_XH2PO4_ROH2_flx=TAD*(XROH2_mole_conc*H2PO4_1e_activity-SPH2P*XH2PO4_mole_conc)/(XROH2_mole_conc+SPH2P)*VLWatMicPBK
      ! R-OH + H2PO4(-) <=> R-H2PO4 + OH(-); positive flux binds phosphate.
      H2PO4_to_XH2PO4_ROH_flx=TAD*(XROH1_mole_conc*H2PO4_1e_activity-SXH2P*XH2PO4_mole_conc*OH_1e_activity) &
        /(XROH1_mole_conc+SXH2P*OH_1e_activity)*VLWatMicPBK
!
!     HPO4 EXCHANGE
!     Intended ligand exchange: R-OH + HPO4(2-) <=> R-HPO4(-) + OH(-).
!     Positive flux removes dissolved HPO4 and adds adsorbed R-HPO4.
!
      ! R-OH + HPO4(2-) <=> R-HPO4(-) + OH(-); positive flux binds phosphate.
      SPH1P=SXH1P*DPH2O/DPH2P
      H1PO4_to_XHPO4_ROH_flx=TAD*(XROH1_mole_conc*H1PO4_2e_activity-SPH1P*XHPO4_mole_conc)/(XROH1_mole_conc+SPH1P)*VLWatMicPBK
      XOH_mole_conc=XOH_mole_conc-RXOH1
      XROH1_mole_conc=XROH1_mole_conc+RXOH1-RXOH2-H2PO4_to_XH2PO4_ROH_flx-H1PO4_to_XHPO4_ROH_flx
      XROH2_mole_conc=XROH2_mole_conc+RXOH2-H2PO4_1e_to_XH2PO4_ROH2_flx
      XHPO4_mole_conc=XHPO4_mole_conc+H1PO4_to_XHPO4_ROH_flx
      XH2PO4_mole_conc=XH2PO4_mole_conc+H2PO4_1e_to_XH2PO4_ROH2_flx+H2PO4_to_XH2PO4_ROH_flx
    ELSE
      RXOH2=0._r8
      RXOH1=0._r8
      H2PO4_1e_to_XH2PO4_ROH2_flx=0._r8
      H2PO4_to_XH2PO4_ROH_flx=0._r8
      H1PO4_to_XHPO4_ROH_flx=0._r8
    ENDIF
!
!     CATION EXCHANGE
!
    ! Gapon exchange: 2X-NH4 + Ca(2+) <=> X2-Ca + 2NH4(+),
    ! with analogous H, Al, Fe, Mg, Na and K equilibria relative to Ca.
    ! X denotes one negative site charge; calculate targets before transfers.
    IF(XCEC.GT.ZEROS)THEN
      AALX=Al_3p_activity**0.333_r8
      AFEX=Fe_3p_activity**0.333_r8
      ACAX=Ca_2p_activity**0.500_r8
      AMGX=Mg_2p_activity**0.500_r8
      XCAX=CEC_conc/(1.0_r8+GKC4*NH4_1p_activity/ACAX &
       +GKCH*H_1p_activity/ACAX+GKCA*AALX/ACAX &
       +GKCA*AFEX/ACAX+GKCM*AMGX/ACAX &
       +GKCN*Na_1p_activity/ACAX+GKCK*K_1p_activity/ACAX)
      XN4Q=XCAX*NH4_1p_activity*GKC4
      XHYQ=XCAX*H_1p_activity*GKCH
      XALQ=XCAX*AALX*GKCA
      XFEQ=XCAX*AFEX*GKCA
      XCAQ=XCAX*ACAX
      XMGQ=XCAX*AMGX*GKCM
      XNAQ=XCAX*Na_1p_activity*GKCN
      XKAQ=XCAX*K_1p_activity*GKCK
      XTLQ=XN4Q+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ
      IF(XTLQ.GT.ZERO)THEN
        FX=CEC_conc/XTLQ
      ELSE
        FX=0._r8
      ENDIF
      XN4Q=FX*XN4Q
      XHYQ=FX*XHYQ
      XALQ=FX*XALQ/3.0_r8
      XFEQ=FX*XFEQ/3.0_r8
      XCAQ=FX*XCAQ/2.0_r8
      XMGQ=FX*XMGQ/2.0_r8
      XNAQ=FX*XNAQ
      XKAQ=FX*XKAQ
      ! Match hourly NH4 exchange: bound desorption before rate scaling.
      ! The same RXN4 enters the dissolved balance with the opposite sign.
      ! X-NH4 <=> X(-) + NH4(+); positive RXN4 adsorbs NH4 onto exchange sites.
      RXN4=TAD*AMAX1(AMIN1( &
        (XN4Q-XNH4_mole_conc)*NH4_1p_activity/XN4Q, &
        NH4_1p_aque_mole_conc),-XNH4_mole_conc)
      ! TODO: Other cation desorption can overdraw exchangeable donors when equilibrium
      ! targets are small. Other reactions also compete for aqueous donors, so
      ! independent reaction limits do not guarantee nonnegative shared pools.
      ! Retain these equations pending replacement with a geochemistry solver
      ! that handles all reactions concurrently and conserves elemental totals.

      ! X-H <=> X(-) + H(+); positive flux adsorbs the cation onto exchange sites.
      RXHY=TAD*AMIN1((XHYQ-XHY1)*H_1p_activity/XHYQ,H_1p_aque_mole_conc)

      ! X3-Al <=> 3X(-) + Al(3+); positive flux adsorbs the cation onto exchange sites.
      RXAL=TAD*AMIN1((XALQ-XAl_conc)*AALX/XALQ,Al_3p_aque_mole_conc)

      ! X3-Fe <=> 3X(-) + Fe(3+); positive flux adsorbs the cation onto exchange sites.
      RXFE=TAD*AMIN1((XFEQ-XFe_conc)*AFEX/XFEQ,Fe_3p_aque_mole_conc)

      ! X2-Ca <=> 2X(-) + Ca(2+); positive flux adsorbs the cation onto exchange sites.
      RXCA=TAD*AMIN1((XCAQ-XCa_conc)*ACAX/XCAQ,Ca_2p_aque_mole_conc)

      ! X2-Mg <=> 2X(-) + Mg(2+); positive flux adsorbs the cation onto exchange sites.
      RXMG=TAD*AMIN1((XMGQ-XMg_conc)*AMGX/XMGQ,Mg_2p_aque_mole_conc)

      ! X-Na <=> X(-) + Na(+); positive flux adsorbs the cation onto exchange sites.
      RXNA=TAD*AMIN1((XNAQ-XNa_conc)*Na_1p_activity/XNAQ,Na_1p_aque_mole_conc)

      ! X-K <=> X(-) + K(+); positive flux adsorbs the cation onto exchange sites.
      RXKA=TAD*AMIN1((XKAQ-XK_conc)*K_1p_activity/XKAQ,K_1p_aque_mole_conc)
      XNH4_mole_conc=XNH4_mole_conc+RXN4
      XHY1=XHY1+RXHY
      XAl_conc=XAl_conc+RXAL
      XFe_conc=XFe_conc+RXFE
      XCa_conc=XCa_conc+RXCA
      XMg_conc=XMg_conc+RXMG
      XNa_conc=XNa_conc+RXNA
      XK_conc=XK_conc+RXKA
    ELSE
      RXN4=0._r8
      RXHY=0._r8
      RXAL=0._r8
      RXFE=0._r8
      RXCA=0._r8
      RXMG=0._r8
      RXNA=0._r8
      RXKA=0._r8
    ENDIF
!
!     ORGANIC MATTER
!
    DP=DPCOH*DPALO
    ! R-COOH <=> R-COO(-) + H(+): solve protonated carboxyl sites.
    XHC1=H_1p_activity*(XCOOH_conc-XAlO2H2_conc-XFeO2H2_conc)/(H_1p_activity+DPCOH)

    ! R-COO(-) + Al(OH)2(+) <=> R-COO-Al(OH)2: solve Al binding.
    XAlO2H2_conc=AlO2H2_1p_activity*(XCOOH_conc-XHC1)/(AlO2H2_1p_activity+DPALO)

    ! R-COO(-) + Fe(OH)2(+) <=> R-COO-Fe(OH)2: solve Fe binding.
    XFeO2H2_conc=FeO2H2_p_activity*(XCOOH_conc-XHC1)/(FeO2H2_p_activity+DPFEO)
    XCOO=AZMAX1(XCOOH_conc-XHC1-XAlO2H2_conc-XFeO2H2_conc)
  ELSE
    RHAL1=0._r8
    RHALO1=0._r8
    RHALO2=0._r8
    RHALO3=0._r8
    RHALO4=0._r8
    RHFE1=0._r8
    RHFEO1=0._r8
    RHFEO2=0._r8
    RHFEO3=0._r8
    RHFEO4=0._r8
    RHCAC3=0._r8
    RHCACH=0._r8
    RHCACO=0._r8
    RPCACX=0._r8
    RPCASO=0._r8
    H2PO4_1e_CaHPO4_dissol_flx=0._r8
    H2PO4_1e_apatite_dissol_flx=0._r8
    RHA0P1=0._r8
    RHA1P1=0._r8
    RHA2P1=0._r8
    RHA3P1=0._r8
    RHA4P1=0._r8
    RHA0P2=0._r8
    RHA1P2=0._r8
    RHA2P2=0._r8
    RHA3P2=0._r8
    RHA4P2=0._r8
    RHF0P1=0._r8
    RHF1P1=0._r8
    RHF2P1=0._r8
    RHF3P1=0._r8
    RHF4P1=0._r8
    RHF0P2=0._r8
    RHF1P2=0._r8
    RHF2P2=0._r8
    RHF3P2=0._r8
    RHF4P2=0._r8
    RPCAD1=0._r8
    RHCAD2=0._r8
    RHCAH1=0._r8
    RHCAH2=0._r8
    RXOH2=0._r8
    RXOH1=0._r8
    H2PO4_1e_to_XH2PO4_ROH2_flx=0._r8
    H2PO4_to_XH2PO4_ROH_flx=0._r8
    H1PO4_to_XHPO4_ROH_flx=0._r8
    RXN4=0._r8
    RXHY=0._r8
    RXAL=0._r8
    RXFE=0._r8
    RXCA=0._r8
    RXMG=0._r8
    RXNA=0._r8
    RXKA=0._r8
  ENDIF
!
!     ION SPECIATION
!
  ! Rationalize S0-sqrt(S1) as 4*(A*B-K*C)/(S0+sqrt(S1)) to retain
  ! small reaction transfers. Preserve TSL, which includes the root factor 1/2.
  ! NH4(+) <=> NH3(aq) + H(+); positive RNH4 forms NH4(+).
  S0=H_1p_activity+NH3_activity+DPN4
  S1=MAX(0._r8,(H_1p_activity-NH3_activity)**2+2._r8*DPN4*(H_1p_activity+NH3_activity)+DPN4**2+4._r8*DPN4*NH4_1p_activity)
  RNH4=TSL*(4._r8*((H_1p_activity*NH3_activity-DPN4*NH4_1p_activity)/(S0+SQRT(S1))))

  ! AlOH(2+) <=> Al(3+) + OH(-); positive RALO1 forms AlOH(2+).
  S0=Al_3p_activity+OH_1e_activity+DPAL1
  S1=MAX(0._r8,(Al_3p_activity-OH_1e_activity)**2+2._r8*DPAL1*(Al_3p_activity+OH_1e_activity)+DPAL1**2+4._r8*DPAL1*AlOH_2p_activity)
  RALO1=TSL*(4._r8*((Al_3p_activity*OH_1e_activity-DPAL1*AlOH_2p_activity)/(S0+SQRT(S1))))

  ! Al(OH)2(+) <=> AlOH(2+) + OH(-); positive RALO2 forms Al(OH)2(+).
  S0=AlOH_2p_activity+OH_1e_activity+DPAL2
  S1=MAX(0._r8,(AlOH_2p_activity-OH_1e_activity)**2+2._r8*DPAL2*(AlOH_2p_activity+OH_1e_activity)+DPAL2**2+4._r8*DPAL2*AlO2H2_1p_activity)
  RALO2=TSL*(4._r8*((AlOH_2p_activity*OH_1e_activity-DPAL2*AlO2H2_1p_activity)/(S0+SQRT(S1))))

  ! Al(OH)3(aq) <=> Al(OH)2(+) + OH(-); positive RALO3 forms Al(OH)3(aq).
  ! Retain the existing concentration-based S0 for this reaction; its
  ! activity-based product is unchanged by the rationalization below.
  S0=AlO2H2_1p_aque_mole_conc+OH_1e_aque_mole_conc+DPAL3
  S1=MAX(0._r8,S0**2-4._r8*(AlO2H2_1p_activity*OH_1e_activity-DPAL3*AlO3H3_activity))
  RALO3=TSL*(4._r8*((AlO2H2_1p_activity*OH_1e_activity-DPAL3*AlO3H3_activity)/(S0+SQRT(S1))))

  ! Al(OH)4(-) <=> Al(OH)3(aq) + OH(-); positive RALO4 forms Al(OH)4(-).
  S0=AlO3H3_activity+OH_1e_activity+DPAL4
  S1=MAX(0._r8,(AlO3H3_activity-OH_1e_activity)**2+2._r8*DPAL4*(AlO3H3_activity+OH_1e_activity)+DPAL4**2+4._r8*DPAL4*AlO4H4_1e_activity)
  RALO4=TSL*(4._r8*((AlO3H3_activity*OH_1e_activity-DPAL4*AlO4H4_1e_activity)/(S0+SQRT(S1))))

  ! AlSO4(+) <=> Al(3+) + SO4(2-); positive RALS forms AlSO4(+).
  S0=Al_3p_activity+SO4_2e_activity+DPALS
  S1=MAX(0._r8,(Al_3p_activity-SO4_2e_activity)**2+2._r8*DPALS*(Al_3p_activity+SO4_2e_activity)+DPALS**2+4._r8*DPALS*AlSO4_1p_activity)
  RALS=TSL*(4._r8*((Al_3p_activity*SO4_2e_activity-DPALS*AlSO4_1p_activity)/(S0+SQRT(S1))))

  ! FeOH(2+) <=> Fe(3+) + OH(-); positive RFEO1 forms FeOH(2+).
  S0=Fe_3p_activity+OH_1e_activity+DPFE1
  S1=MAX(0._r8,(Fe_3p_activity-OH_1e_activity)**2+2._r8*DPFE1*(Fe_3p_activity+OH_1e_activity)+DPFE1**2+4._r8*DPFE1*FeOH_2p_activity)
  RFEO1=TSL*(4._r8*((Fe_3p_activity*OH_1e_activity-DPFE1*FeOH_2p_activity)/(S0+SQRT(S1))))

  ! Fe(OH)2(+) <=> FeOH(2+) + OH(-); positive RFEO2 forms Fe(OH)2(+).
  S0=FeOH_2p_activity+OH_1e_activity+DPFE2
  S1=MAX(0._r8,(FeOH_2p_activity-OH_1e_activity)**2+2._r8*DPFE2*(FeOH_2p_activity+OH_1e_activity)+DPFE2**2+4._r8*DPFE2*FeO2H2_p_activity)
  RFEO2=TSL*(4._r8*((FeOH_2p_activity*OH_1e_activity-DPFE2*FeO2H2_p_activity)/(S0+SQRT(S1))))

  ! Fe(OH)3(aq) <=> Fe(OH)2(+) + OH(-); positive RFEO3 forms Fe(OH)3(aq).
  S0=FeO2H2_p_activity+OH_1e_activity+DPFE3
  S1=MAX(0._r8,(FeO2H2_p_activity-OH_1e_activity)**2+2._r8*DPFE3*(FeO2H2_p_activity+OH_1e_activity)+DPFE3**2+4._r8*DPFE3*FeO3H3_activity)
  RFEO3=TSL*(4._r8*((FeO2H2_p_activity*OH_1e_activity-DPFE3*FeO3H3_activity)/(S0+SQRT(S1))))

  ! Fe(OH)4(-) <=> Fe(OH)3(aq) + OH(-); positive RFEO4 forms Fe(OH)4(-).
  S0=FeO3H3_activity+OH_1e_activity+DPFE4
  S1=MAX(0._r8,(FeO3H3_activity-OH_1e_activity)**2+2._r8*DPFE4*(FeO3H3_activity+OH_1e_activity)+DPFE4**2+4._r8*DPFE4*FeO4H4_1e_activity)
  RFEO4=TSL*(4._r8*((FeO3H3_activity*OH_1e_activity-DPFE4*FeO4H4_1e_activity)/(S0+SQRT(S1))))

  ! FeSO4(+) <=> Fe(3+) + SO4(2-); positive RFES forms FeSO4(+).
  S0=Fe_3p_activity+SO4_2e_activity+DPFES
  S1=MAX(0._r8,(Fe_3p_activity-SO4_2e_activity)**2+2._r8*DPFES*(Fe_3p_activity+SO4_2e_activity)+DPFES**2+4._r8*DPFES*FeSO4_1p_activity)
  RFES=TSL*(4._r8*((Fe_3p_activity*SO4_2e_activity-DPFES*FeSO4_1p_activity)/(S0+SQRT(S1))))

  ! CaOH(+) <=> Ca(2+) + OH(-); positive RCAO forms CaOH(+).
  S0=Ca_2p_activity+OH_1e_activity+DPCAO
  S1=MAX(0._r8,(Ca_2p_activity-OH_1e_activity)**2+2._r8*DPCAO*(Ca_2p_activity+OH_1e_activity)+DPCAO**2+4._r8*DPCAO*CaO2H2_activity)
  RCAO=TSL*(4._r8*((Ca_2p_activity*OH_1e_activity-DPCAO*CaO2H2_activity)/(S0+SQRT(S1))))

  ! CaCO3(aq) <=> Ca(2+) + CO3(2-); positive RCAC forms CaCO3(aq).
  S0=Ca_2p_activity+CO3_2e_activity+DPCAC
  S1=MAX(0._r8,(Ca_2p_activity-CO3_2e_activity)**2+2._r8*DPCAC*(Ca_2p_activity+CO3_2e_activity)+DPCAC**2+4._r8*DPCAC*CaCO3_activity)
  RCAC=TSL*(4._r8*((Ca_2p_activity*CO3_2e_activity-DPCAC*CaCO3_activity)/(S0+SQRT(S1))))

  ! CaHCO3(+) <=> Ca(2+) + HCO3(-); positive RCAH forms CaHCO3(+).
  S0=Ca_2p_activity+HCO3_e_activity+DPCAH
  S1=MAX(0._r8,(Ca_2p_activity-HCO3_e_activity)**2+2._r8*DPCAH*(Ca_2p_activity+HCO3_e_activity)+DPCAH**2+4._r8*DPCAH*CaHCO3_1p_activity)
  RCAH=TSL*(4._r8*((Ca_2p_activity*HCO3_e_activity-DPCAH*CaHCO3_1p_activity)/(S0+SQRT(S1))))

  ! CaSO4(aq) <=> Ca(2+) + SO4(2-); positive RCAS forms CaSO4(aq).
  S0=Ca_2p_activity+SO4_2e_activity+DPCAS
  S1=MAX(0._r8,(Ca_2p_activity-SO4_2e_activity)**2+2._r8*DPCAS*(Ca_2p_activity+SO4_2e_activity)+DPCAS**2+4._r8*DPCAS*CaSO4_activity)
  RCAS=TSL*(4._r8*((Ca_2p_activity*SO4_2e_activity-DPCAS*CaSO4_activity)/(S0+SQRT(S1))))

  ! MgOH(+) <=> Mg(2+) + OH(-); positive RMGO forms MgOH(+).
  S0=Mg_2p_activity+OH_1e_activity+DPMGO
  S1=MAX(0._r8,(Mg_2p_activity-OH_1e_activity)**2+2._r8*DPMGO*(Mg_2p_activity+OH_1e_activity)+DPMGO**2+4._r8*DPMGO*MgOH_1p_activity)
  RMGO=TSL*(4._r8*((Mg_2p_activity*OH_1e_activity-DPMGO*MgOH_1p_activity)/(S0+SQRT(S1))))

  ! MgCO3(aq) <=> Mg(2+) + CO3(2-); positive RMGC forms MgCO3(aq).
  S0=Mg_2p_activity+CO3_2e_activity+DPMGC
  S1=MAX(0._r8,(Mg_2p_activity-CO3_2e_activity)**2+2._r8*DPMGC*(Mg_2p_activity+CO3_2e_activity)+DPMGC**2+4._r8*DPMGC*MgCO3_activity)
  RMGC=TSL*(4._r8*((Mg_2p_activity*CO3_2e_activity-DPMGC*MgCO3_activity)/(S0+SQRT(S1))))

  ! MgHCO3(+) <=> Mg(2+) + HCO3(-); positive RMGH forms MgHCO3(+).
  S0=Mg_2p_activity+HCO3_e_activity+DPMGH
  S1=MAX(0._r8,(Mg_2p_activity-HCO3_e_activity)**2+2._r8*DPMGH*(Mg_2p_activity+HCO3_e_activity)+DPMGH**2+4._r8*DPMGH*MgHCO3_1p_activity)
  RMGH=TSL*(4._r8*((Mg_2p_activity*HCO3_e_activity-DPMGH*MgHCO3_1p_activity)/(S0+SQRT(S1))))

  ! MgSO4(aq) <=> Mg(2+) + SO4(2-); positive RMGS forms MgSO4(aq).
  S0=Mg_2p_activity+SO4_2e_activity+DPMGS
  S1=MAX(0._r8,(Mg_2p_activity-SO4_2e_activity)**2+2._r8*DPMGS*(Mg_2p_activity+SO4_2e_activity)+DPMGS**2+4._r8*DPMGS*MgSO4_activity)
  RMGS=TSL*(4._r8*((Mg_2p_activity*SO4_2e_activity-DPMGS*MgSO4_activity)/(S0+SQRT(S1))))

  ! NaCO3(-) <=> Na(+) + CO3(2-); positive RNAC forms NaCO3(-).
  S0=Na_1p_activity+CO3_2e_activity+DPNAC
  S1=MAX(0._r8,(Na_1p_activity-CO3_2e_activity)**2+2._r8*DPNAC*(Na_1p_activity+CO3_2e_activity)+DPNAC**2+4._r8*DPNAC*NaCO3_1e_activity)
  RNAC=TSL*(4._r8*((Na_1p_activity*CO3_2e_activity-DPNAC*NaCO3_1e_activity)/(S0+SQRT(S1))))

  ! NaSO4(-) <=> Na(+) + SO4(2-); positive RNAS forms NaSO4(-).
  S0=Na_1p_activity+SO4_2e_activity+DPNAS
  S1=MAX(0._r8,(Na_1p_activity-SO4_2e_activity)**2+2._r8*DPNAS*(Na_1p_activity+SO4_2e_activity)+DPNAS**2+4._r8*DPNAS*NaSO4_1e_activity)
  RNAS=TSL*(4._r8*((Na_1p_activity*SO4_2e_activity-DPNAS*NaSO4_1e_activity)/(S0+SQRT(S1))))

  ! KSO4(-) <=> K(+) + SO4(2-); positive RKAS forms KSO4(-).
  S0=K_1p_activity+SO4_2e_activity+DPKAS
  S1=MAX(0._r8,(K_1p_activity-SO4_2e_activity)**2+2._r8*DPKAS*(K_1p_activity+SO4_2e_activity)+DPKAS**2+4._r8*DPKAS*KSO4_1e_activity)
  RKAS=TSL*(4._r8*((K_1p_activity*SO4_2e_activity-DPKAS*KSO4_1e_activity)/(S0+SQRT(S1))))

  ! HPO4(2-) <=> PO4(3-) + H(+); positive RH1P forms HPO4(2-).
  S0=H0PO4_3e_activity+H_1p_activity+DPH1P
  S1=MAX(0._r8,(H0PO4_3e_activity-H_1p_activity)**2+2._r8*DPH1P*(H0PO4_3e_activity+H_1p_activity)+DPH1P**2+4._r8*DPH1P*H1PO4_2e_activity)
  RH1P=TSL*(4._r8*((H0PO4_3e_activity*H_1p_activity-DPH1P*H1PO4_2e_activity)/(S0+SQRT(S1))))

  ! H2PO4(-) <=> HPO4(2-) + H(+); positive H2PO4_e_to_HPO4_2e_flx forms H2PO4(-).
  S0=H1PO4_2e_activity+H_1p_activity+DPH2P
  S1=MAX(0._r8,(H1PO4_2e_activity-H_1p_activity)**2+2._r8*DPH2P*(H1PO4_2e_activity+H_1p_activity)+DPH2P**2+4._r8*DPH2P*H2PO4_1e_activity)
  H2PO4_e_to_HPO4_2e_flx=TSL*(4._r8*((H1PO4_2e_activity*H_1p_activity-DPH2P*H2PO4_1e_activity)/(S0+SQRT(S1))))

  ! H3PO4 <=> H2PO4(-) + H(+); positive RH3P forms H3PO4.
  S0=H2PO4_1e_activity+H_1p_activity+DPH3P
  S1=MAX(0._r8,(H2PO4_1e_activity-H_1p_activity)**2+2._r8*DPH3P*(H2PO4_1e_activity+H_1p_activity)+DPH3P**2+4._r8*DPH3P*H3PO4_activity)
  RH3P=TSL*(4._r8*((H2PO4_1e_activity*H_1p_activity-DPH3P*H3PO4_activity)/(S0+SQRT(S1))))

  ! FeHPO4(+) <=> Fe(3+) + HPO4(2-); positive RF1P forms FeHPO4(+).
  S0=Fe_3p_activity+H1PO4_2e_activity+DPF1P
  S1=MAX(0._r8,(Fe_3p_activity-H1PO4_2e_activity)**2+2._r8*DPF1P*(Fe_3p_activity+H1PO4_2e_activity)+DPF1P**2+4._r8*DPF1P*FeHPO4_p_activity)
  RF1P=TSL*(4._r8*((Fe_3p_activity*H1PO4_2e_activity-DPF1P*FeHPO4_p_activity)/(S0+SQRT(S1))))

  ! FeH2PO4(2+) <=> Fe(3+) + H2PO4(-); positive RF2P forms FeH2PO4(2+).
  S0=Fe_3p_activity+H2PO4_1e_activity+DPF2P
  S1=MAX(0._r8,(Fe_3p_activity-H2PO4_1e_activity)**2+2._r8*DPF2P*(Fe_3p_activity+H2PO4_1e_activity)+DPF2P**2+4._r8*DPF2P*FeH2PO4_2p_activity)
  RF2P=TSL*(4._r8*((Fe_3p_activity*H2PO4_1e_activity-DPF2P*FeH2PO4_2p_activity)/(S0+SQRT(S1))))

  ! CaPO4(-) <=> Ca(2+) + PO4(3-); positive RC0P forms CaPO4(-).
  S0=Ca_2p_activity+H0PO4_3e_activity+DPC0P
  S1=MAX(0._r8,(Ca_2p_activity-H0PO4_3e_activity)**2+2._r8*DPC0P*(Ca_2p_activity+H0PO4_3e_activity)+DPC0P**2+4._r8*DPC0P*CaPO4_1e_activity)
  RC0P=TSL*(4._r8*((Ca_2p_activity*H0PO4_3e_activity-DPC0P*CaPO4_1e_activity)/(S0+SQRT(S1))))

  ! CaHPO4(aq) <=> Ca(2+) + HPO4(2-); positive RC1P forms CaHPO4(aq).
  S0=Ca_2p_activity+H1PO4_2e_activity+DPC1P
  S1=MAX(0._r8,(Ca_2p_activity-H1PO4_2e_activity)**2+2._r8*DPC1P*(Ca_2p_activity+H1PO4_2e_activity)+DPC1P**2+4._r8*DPC1P*CaHPO4_activity)
  RC1P=TSL*(4._r8*((Ca_2p_activity*H1PO4_2e_activity-DPC1P*CaHPO4_activity)/(S0+SQRT(S1))))

  ! Legacy dissolved CaH4P2O8 variables represent the one-P CaH2PO4(+) complex:
  ! CaH2PO4(+) <-> Ca(2+) + H2PO4(-), distinct from solid Ca(H2PO4)2.
  ! Positive RC2P forms the dissolved one-P complex CaH2PO4(+).
  S0=Ca_2p_activity+H2PO4_1e_activity+DPC2P
  S1=MAX(0._r8,(Ca_2p_activity-H2PO4_1e_activity)**2+2._r8*DPC2P*(Ca_2p_activity+H2PO4_1e_activity)+DPC2P**2+4._r8*DPC2P*CaH4P2O8_1p_activity)
  RC2P=TSL*(4._r8*((Ca_2p_activity*H2PO4_1e_activity-DPC2P*CaH4P2O8_1p_activity)/(S0+SQRT(S1))))

  ! MgHPO4(aq) <=> Mg(2+) + HPO4(2-); positive RM1P forms MgHPO4(aq).
  S0=Mg_2p_activity+H1PO4_2e_activity+DPM1P
  S1=MAX(0._r8,(Mg_2p_activity-H1PO4_2e_activity)**2+2._r8*DPM1P*(Mg_2p_activity+H1PO4_2e_activity)+DPM1P**2+4._r8*DPM1P*MgHPO4_activity)
  RM1P=TSL*(4._r8*((Mg_2p_activity*H1PO4_2e_activity-DPM1P*MgHPO4_activity)/(S0+SQRT(S1))))
!
!     TOTAL ION FLUX FOR EACH ION SPECIES
!
  RN4S=RNH4-RXN4
  RN3S=-RNH4
  RAL=-RHAL1-RHA0P1-RHA0P2-RALO1-RALS-RXAL
  RFE=-RHFE1-RHF0P1-RHF0P2-RFEO1-RFES-RXFE-RF1P-RF2P
  RCA=-RPCACX-RPCASO-H2PO4_1e_CaHPO4_dissol_flx-5.0*H2PO4_1e_apatite_dissol_flx-RXCA-RCAO-RCAC-RCAH-RCAS-RC0P-RC1P-RC2P
  RMG=-RMGO-RMGC-RMGH-RMGS-RM1P-RXMG
  RNA=-RNAC-RNAS-RXNA
  RKA=-RKAS-RXKA
  RSO4=-RPCASO-RALS-RFES-RCAS-RMGS-RNAS-RKAS
  RAL1=-RHALO1-RHA1P1-RHA1P2+RALO1-RALO2
  RAL2=-RHALO2-RHA2P1-RHA2P2+RALO2-RALO3
  RAL3=-RHALO3-RHA3P1-RHA3P2+RALO3-RALO4
  RAL4=-RHALO4-RHA4P1-RHA4P2+RALO4
  RFE1=-RHFEO1-RHF1P1-RHF1P2+RFEO1-RFEO2
  RFE2=-RHFEO2-RHF2P1-RHF2P2+RFEO2-RFEO3
  RFE3=-RHFEO3-RHF3P1-RHF3P2+RFEO3-RFEO4
  RFE4=-RHFEO4-RHF4P1-RHF4P2+RFEO4
  RHP0=-RH1P-RC0P
  RHP1=-RHA0P1-RHA1P1-RHA2P1-RHA3P1 &
    -RHA4P1-RHF0P1-RHF1P1-RHF2P1 &
    -RHF3P1-RHF4P1-RPCAD1-3.0*RHCAH1 &
    -H1PO4_to_XHPO4_ROH_flx+RH1P-H2PO4_e_to_HPO4_2e_flx-RF1P-RC1P-RM1P
  RHP2=-RHA0P2-RHA1P2-RHA2P2-RHA3P2 &
    -RHA4P2-RHF0P2-RHF1P2-RHF2P2 &
    -RHF3P2-RHF4P2-RHCAD2-3.0*RHCAH2 &
    -H2PO4_1e_to_XH2PO4_ROH2_flx-H2PO4_to_XH2PO4_ROH_flx+H2PO4_e_to_HPO4_2e_flx-RH3P-RF2P-RC2P
  RHP3=RH3P
!
!     ION CONCENTRATIONS
!
  ! CO2(g) <=> H2CO3(aq), coupled to HCO3(-) and CO3(2-) formation.
  ! Iterate dissolved CO2 with the salting-out correction.
  CCO2X=CCO2M*SCO2X/(EXP(GasSechenovConst(idg_CO2)*CSTRZ))*EXP(0.843-0.0281*ATCA)*FH2O
  CCO2Y=LOG(CCO2X)
  CCO2Z=ABS(CCO2Y)
  H2CO3_aque_mole_conc=CCO2X
  FCO3=DPCO3*A0/(H_1p_activity**2*A2)
  FHCO=DPCO2*A0/(H_1p_activity*A1)
  Z=GasSechenovConst(idg_CO2)*(2.0E-03*FCO3+0.5E-03*FHCO)
  DO  MM=1,25
    R=(LOG(H2CO3_aque_mole_conc)+Z*H2CO3_aque_mole_conc-CCO2Y)/CCO2Z
    IF(R.LT.1.0E-03)exit
    H2CO3_aque_mole_conc=H2CO3_aque_mole_conc/SQRT(1.0_r8+R)
  ENDDO

  ! CH4(g) <=> CH4(aq): refresh Henry partitioning at the current ionic strength.
  CH4_aque_mole_conc    = CCH4M*SCH4X/(EXP(GasSechenovConst(idg_CH4)*CSTR1))*EXP(0.597-0.0199*ATCA)*FH2O

  ! O2(g) <=> O2(aq): refresh Henry partitioning at the current ionic strength.
  O2_aque_mole_conc     = COXYM*SOXYX/(EXP(GasSechenovConst(idg_O2)*CSTR1))*EXP(0.516-0.0172*ATCA)*FH2O

  ! N2(g) <=> N2(aq): refresh Henry partitioning at the current ionic strength.
  N2_aque_mole_conc     = CZ2GM*SN2GX/(EXP(GasSechenovConst(idg_N2)*CSTR1))*EXP(0.456-0.0152*ATCA)*FH2O

  ! N2O(g) <=> N2O(aq): refresh Henry partitioning at the current ionic strength.
  N2O_aque_mole_conc    = CZ2OM*SN2OX/(EXP(GasSechenovConst(idg_N2O)*CSTR1))*EXP(0.897-0.0299*ATCA)*FH2O
  NH4_1p_aque_mole_conc = NH4_1p_aque_mole_conc+RN4S
  NH3_aque_mole_conc    = NH3_aque_mole_conc+RN3S
  Al_3p_aque_mole_conc  = Al_3p_aque_mole_conc+RAL
  Fe_3p_aque_mole_conc  = Fe_3p_aque_mole_conc+RFE
  Ca_2p_aque_mole_conc  = Ca_2p_aque_mole_conc+RCA
  Mg_2p_aque_mole_conc  = Mg_2p_aque_mole_conc+RMG
  Na_1p_aque_mole_conc  = Na_1p_aque_mole_conc+RNA
  K_1p_aque_mole_conc   = K_1p_aque_mole_conc+RKA
  SO4_2e_aque_mole_conc    = SO4_2e_aque_mole_conc+RSO4
  CO3_2e_aque_mole_conc    = H2CO3_aque_mole_conc*DPCO3*A0/(H_1p_activity**2*A2)
  HCO3_e_conc              = H2CO3_aque_mole_conc*DPCO2*A0/(H_1p_activity*A1)
  AlOH_2p_aque_mole_conc   = AlOH_2p_aque_mole_conc+RAL1
  AlO2H2_1p_aque_mole_conc = AlO2H2_1p_aque_mole_conc+RAL2
  AlO3H3_conc              = AlO3H3_conc+RAL3
  AlO4H4_1e_aque_mole_conc = AlO4H4_1e_aque_mole_conc+RAL4
  AlSO4_1p_aque_mole_conc  = AlSO4_1p_aque_mole_conc+RALS
  FeOH_2p_aque_mole_conc   = FeOH_2p_aque_mole_conc+RFE1
  FeO2H2_p_conc            = FeO2H2_p_conc+RFE2
  FeO3H3_conc              = FeO3H3_conc+RFE3
  FeO4H4_1e_aque_mole_conc = FeO4H4_1e_aque_mole_conc+RFE4
  FeSO4_1p_aque_mole_conc  = FeSO4_1p_aque_mole_conc+RFES
  CaO2H2_conc              = CaO2H2_conc+RCAO
  CaCO3_conc               = CaCO3_conc+RCAC
  CaHCO3_1p_aque_mole_conc = CaHCO3_1p_aque_mole_conc+RCAH
  CaSO4_conc               = CaSO4_conc+RCAS
  MgOH_1p_aque_mole_conc   = MgOH_1p_aque_mole_conc+RMGO
  MgCO3_conc               = MgCO3_conc+RMGC
  MgHCO3_1p_aque_mole_conc   = MgHCO3_1p_aque_mole_conc+RMGH
  MgSO4_conc                 = MgSO4_conc+RMGS
  NaCO3_1e_aque_mole_conc    = NaCO3_1e_aque_mole_conc+RNAC
  NaSO4_1e_aque_mole_conc    = NaSO4_1e_aque_mole_conc+RNAS
  KSO4_1e_aque_mole_conc     = KSO4_1e_aque_mole_conc+RKAS
  H0PO4_3e_conc              = H0PO4_3e_conc+RHP0
  H1PO4_2e_aque_mole_conc    = H1PO4_2e_aque_mole_conc+RHP1
  H2PO4_1e_aque_mole_conc    = H2PO4_1e_aque_mole_conc+RHP2
  H3PO4_conc                 = H3PO4_conc+RHP3
  FeHPO4_p_conc              = FeHPO4_p_conc+RF1P
  FeH2PO4_2p_aque_mole_conc  = FeH2PO4_2p_aque_mole_conc+RF2P
  CaPO4_1e_con               = CaPO4_1e_con+RC0P
  CaHPO4_conc                = CaHPO4_conc+RC1P
  CaH4P2O8_1p_aque_mole_conc = CaH4P2O8_1p_aque_mole_conc+RC2P
  MgHPO4_conc                = MgHPO4_conc+RM1P

  end subroutine SolubilityEquilibiriaSalt
!------------------------------------------------------------------------------------------

  subroutine SolubilityEquilibriaNoSalt(K,M,BulkSoilMass)

  implicit none
  integer, intent(in) :: K, M
  real(r8), intent(in) :: BulkSoilMass
  real(r8) :: H2PO4_1e_AlPO4_eqv,H2PO4_1e_CaHPO4_eqv,XNAQ,XCAQ,SPH1P
  real(r8) :: H2PO4_1e_FePO4_eqv,H2PO4_1e_apatite_eqv,FX,RN4S,RN3S,RNH4
  real(r8) :: H2PO4_1e_AlPO4_dissol_flx,H2PO4_1e_CaHPO4_dissol_flx,H2PO4_e_to_HPO4_2e_flx,RHP1,RHP2
  real(r8) :: H2PO4_1e_apatite_dissol_flx,H2PO4_1e_FePO4_dissol_flx,H1PO4_to_XHPO4_ROH_flx,H2PO4_1e_to_XH2PO4_ROH2_flx
  real(r8) :: RXN4,H2PO4_to_XH2PO4_ROH_flx,SPH2P,VLWatMicPBK,S0,S1
  real(r8) :: XALQ,XCAX,XFEQ,XHYQ,XKAQ,XMGQ,XN4Q
  real(r8) :: XTLQ
! begin_execution
  H2CO3_aque_mole_conc=AMAX1(ZERO,H2CO3_aque_mole_conc)

  ! H2CO3(aq) <=> CO3(2-) + 2H(+): combine the two carbonate acid dissociations.
  CO3_2e_aque_mole_conc=H2CO3_aque_mole_conc*DPCO3*A0/(H_1p_aque_mole_conc**2*A2)

  ! H2CO3(aq) <=> HCO3(-) + H(+): evaluate bicarbonate at prescribed pH.
  HCO3_e_conc           = H2CO3_aque_mole_conc*DPCO2*A0/(H_1p_aque_mole_conc*A1)
  NH4_1p_aque_mole_conc = AMAX1(ZERO,NH4_1p_aque_mole_conc)
  NH3_aque_mole_conc    = AMAX1(ZERO,NH3_aque_mole_conc)
  Al_3p_aque_mole_conc  = AMAX1(ZERO,Al_3p_aque_mole_conc)
  Fe_3p_aque_mole_conc  = AMAX1(ZERO,Fe_3p_aque_mole_conc)
  Ca_2p_aque_mole_conc  = AMAX1(ZERO,Ca_2p_aque_mole_conc)

  ! CaCO3(s) <=> Ca(2+) + CO3(2-): cap free Ca at calcite saturation.
  Ca_2p_aque_mole_conc    = AMIN1(Ca_2p_aque_mole_conc,SPCAC/(CO3_2e_aque_mole_conc*A2**2))
  Mg_2p_aque_mole_conc    = AMAX1(ZERO,Mg_2p_aque_mole_conc)
  Na_1p_aque_mole_conc    = AMAX1(ZERO,Na_1p_aque_mole_conc)
  K_1p_aque_mole_conc     = AMAX1(ZERO,K_1p_aque_mole_conc)
  H1PO4_2e_aque_mole_conc = AMAX1(ZERO,H1PO4_2e_aque_mole_conc)
  H2PO4_1e_aque_mole_conc = AMAX1(ZERO,H2PO4_1e_aque_mole_conc)
!
!     PRECIPITATION-DISSOLUTION FLUXES
!
  IF(K.EQ.micpar%k_POM)THEN
    ! AlPO4(s) + 2H2O <=> Al(3+) + H2PO4(-) + 2OH(-).
    ! Positive mineral flux precipitates the solid; negative flux dissolves it.
    H2PO4_1e_AlPO4_eqv=SYA0P2/(Al_3p_aque_mole_conc*OH_1e_aque_mole_conc**2)
    H2PO4_1e_AlPO4_dissol_flx=AMAX1(-Precp_AlPO4_mole_conc,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_AlPO4_eqv))
 
    ! FePO4(s) + 2H2O <=> Fe(3+) + H2PO4(-) + 2OH(-).
    ! Positive mineral flux precipitates the solid; negative flux dissolves it.
    H2PO4_1e_FePO4_eqv=SYF0P2/(Fe_3p_aque_mole_conc*OH_1e_aque_mole_conc**2._r8)
    H2PO4_1e_FePO4_dissol_flx=AMAX1(-Precp_FePO4_mole_conc,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_FePO4_eqv))
 
    ! CaHPO4(s) + H2O <=> Ca(2+) + H2PO4(-) + OH(-).
    ! Positive mineral flux precipitates the solid; negative flux dissolves it.
    H2PO4_1e_CaHPO4_eqv=SYCAD2/(Ca_2p_aque_mole_conc*OH_1e_aque_mole_conc)
    H2PO4_1e_CaHPO4_dissol_flx=AMAX1(-Precp_CaHPO4_mole_conc,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_CaHPO4_eqv))
 
    ! Ca5(PO4)3OH(s) + 6H2O <=> 5Ca(2+) + 3H2PO4(-) + 7OH(-).
    ! Positive mineral flux precipitates the solid; negative flux dissolves it.
    H2PO4_1e_apatite_eqv=(SYCAH2/(Ca_2p_aque_mole_conc**5*OH_1e_aque_mole_conc**7))**0.333
    H2PO4_1e_apatite_dissol_flx=AMAX1(-Precp_Ca5P3O12O3H3_mole_conc,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_apatite_eqv))
    Precp_AlPO4_mole_conc=Precp_AlPO4_mole_conc+H2PO4_1e_AlPO4_dissol_flx
    Precp_FePO4_mole_conc=Precp_FePO4_mole_conc+H2PO4_1e_FePO4_dissol_flx
    Precp_CaHPO4_mole_conc=Precp_CaHPO4_mole_conc+H2PO4_1e_CaHPO4_dissol_flx
    Precp_Ca5P3O12O3H3_mole_conc=Precp_Ca5P3O12O3H3_mole_conc+H2PO4_1e_apatite_dissol_flx
!
!     ANION EXCHANGE FLUXES
!
    IF(VLWatMicP.GT.ZEROS)THEN
      VLWatMicPBK=AMIN1(1.0,BulkSoilMass/VLWatMicP)
    ELSE
      VLWatMicPBK=1.0
    ENDIF
    IF(XAEC.GT.ZEROS)THEN
      ! Intended ligand exchange: R-OH2(+) + H2PO4(-) <=> R-H2PO4 + H2O.
      ! Positive flux binds dissolved H2PO4 to a protonated surface site.
      ! R-H2PO4 + H2O <=> R-OH2(+) + H2PO4(-).
      ! Combine neutral-site exchange (SXH2P), water dissociation (DPH2O),
      ! and R-OH2(+) <=> R-OH + H(+) (SXOH2) to close the equilibrium cycle.
      SPH2P = SXH2P*DPH2O/SXOH2
      H2PO4_1e_to_XH2PO4_ROH2_flx=TAD*(XROH2_mole_conc*H2PO4_1e_aque_mole_conc-SPH2P*XH2PO4_mole_conc)/(XROH2_mole_conc+SPH2P)*VLWatMicPBK

      ! R-OH + H2PO4(-) <=> R-H2PO4 + OH(-); positive flux binds phosphate.
      H2PO4_to_XH2PO4_ROH_flx=TAD*(XROH1_mole_conc*H2PO4_1e_aque_mole_conc-SXH2P*OH_1e_aque_mole_conc*XH2PO4_mole_conc) &
        /(XROH1_mole_conc+SXH2P*OH_1e_aque_mole_conc)*VLWatMicPBK

      ! R-OH + HPO4(2-) <=> R-HPO4(-) + OH(-); positive flux binds phosphate.
      SPH1P=SXH1P*DPH2O/DPH2P
      H1PO4_to_XHPO4_ROH_flx=TAD*(XROH1_mole_conc*H1PO4_2e_aque_mole_conc-SPH1P*XHPO4_mole_conc)/(XROH1_mole_conc+SPH1P)*VLWatMicPBK
      XROH1_mole_conc=XROH1_mole_conc-H2PO4_to_XH2PO4_ROH_flx-H1PO4_to_XHPO4_ROH_flx
      XROH2_mole_conc=XROH2_mole_conc-H2PO4_1e_to_XH2PO4_ROH2_flx
      XHPO4_mole_conc=XHPO4_mole_conc+H1PO4_to_XHPO4_ROH_flx
      XH2PO4_mole_conc=XH2PO4_mole_conc+H2PO4_1e_to_XH2PO4_ROH2_flx+H2PO4_to_XH2PO4_ROH_flx
    ELSE
      H2PO4_1e_to_XH2PO4_ROH2_flx=0._r8
      H2PO4_to_XH2PO4_ROH_flx=0._r8
      H1PO4_to_XHPO4_ROH_flx=0._r8
    ENDIF
!
!     CATION EXCHANGE
!
    ! Gapon exchange: 2X-NH4 + Ca(2+) <=> X2-Ca + 2NH4(+),
    ! with analogous H, Al, Fe, Mg, Na and K equilibria relative to Ca.
    ! X denotes one negative site charge; calculate targets before transfers.
    IF(XCEC.GT.ZEROS)THEN
      CALX=Al_3p_aque_mole_conc**0.333
      CFEX=Fe_3p_aque_mole_conc**0.333
      CaX_conc=Ca_2p_aque_mole_conc**0.500_r8
      MgX_conc=Mg_2p_aque_mole_conc**0.500_r8
      XCAX=CEC_conc/(1.0_r8+GKC4*NH4_1p_aque_mole_conc/CaX_conc &
        +GKCH*H_1p_aque_mole_conc/CaX_conc+GKCA*CALX/CaX_conc &
        +GKCA*CFEX/CaX_conc+GKCM*MgX_conc/CaX_conc &
        +GKCN*Na_1p_aque_mole_conc/CaX_conc+GKCK*K_1p_aque_mole_conc/CaX_conc)
      XN4Q = XCAX*NH4_1p_aque_mole_conc*GKC4
      XHYQ = XCAX*H_1p_aque_mole_conc*GKCH
      XALQ = XCAX*CALX*GKCA
      XFEQ = XCAX*CFEX*GKCA
      XCAQ = XCAX*CaX_conc
      XMGQ = XCAX*MgX_conc*GKCM
      XNAQ = XCAX*Na_1p_aque_mole_conc*GKCN
      XKAQ = XCAX*K_1p_aque_mole_conc*GKCK
      XTLQ = XN4Q+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ
      IF(XTLQ.GT.ZERO)THEN
        FX=CEC_conc/XTLQ
      ELSE
        FX=0._r8
      ENDIF
      XN4Q=FX*XN4Q
      ! Match hourly NH4 exchange: bound desorption before rate scaling.
      ! Preserve the paired dissolved transfer; shared aqueous donors still need
      ! the concurrent reaction treatment noted in the salt initialization block.
      ! X-NH4 <=> X(-) + NH4(+); positive RXN4 adsorbs NH4 onto exchange sites.
      RXN4=TSL*AMAX1(AMIN1( &
        (XN4Q-XNH4_mole_conc)*NH4_1p_aque_mole_conc/XN4Q, &
        NH4_1p_aque_mole_conc),-XNH4_mole_conc)
      XNH4_mole_conc=XNH4_mole_conc+RXN4
    ELSE
      RXN4=0._r8
    ENDIF
  ELSE
    H2PO4_1e_AlPO4_dissol_flx   = 0._r8
    H2PO4_1e_FePO4_dissol_flx   = 0._r8
    H2PO4_1e_CaHPO4_dissol_flx  = 0._r8
    H2PO4_1e_apatite_dissol_flx = 0._r8
    H2PO4_1e_to_XH2PO4_ROH2_flx = 0._r8
    H2PO4_to_XH2PO4_ROH_flx     = 0._r8
    H1PO4_to_XHPO4_ROH_flx      = 0._r8
    RXN4=0._r8
  ENDIF
!
!     NH4 <-> NH3
!
  ! NH4(+) <=> NH3(aq) + H(+); positive RNH4 forms NH4(+).
  S0=H_1p_aque_mole_conc+NH3_aque_mole_conc+DPN4
  S1=MAX(0._r8,(H_1p_aque_mole_conc-NH3_aque_mole_conc)**2+2._r8*DPN4*(H_1p_aque_mole_conc+NH3_aque_mole_conc)+DPN4**2+4._r8*DPN4*NH4_1p_aque_mole_conc)
  RNH4=TSL*(4._r8*((H_1p_aque_mole_conc*NH3_aque_mole_conc-DPN4*NH4_1p_aque_mole_conc)/(S0+SQRT(S1))))
!
!     H2PO4 <-> HPO4
!
  ! H2PO4(-) <=> HPO4(2-) + H(+); positive flux forms H2PO4(-).
  S0=H1PO4_2e_aque_mole_conc+H_1p_aque_mole_conc+DPH2P
  S1=MAX(0._r8,(H1PO4_2e_aque_mole_conc-H_1p_aque_mole_conc)**2+2._r8*DPH2P*(H1PO4_2e_aque_mole_conc+H_1p_aque_mole_conc)+DPH2P**2+4._r8*DPH2P*H2PO4_1e_aque_mole_conc)
  H2PO4_e_to_HPO4_2e_flx=TSL*(4._r8*((H1PO4_2e_aque_mole_conc*H_1p_aque_mole_conc-DPH2P*H2PO4_1e_aque_mole_conc)/(S0+SQRT(S1))))
!
!     ION FLUXES
!
  RN4S=RNH4-RXN4
  RN3S=-RNH4
  RHP1=-H1PO4_to_XHPO4_ROH_flx-H2PO4_e_to_HPO4_2e_flx
  RHP2=-H2PO4_1e_to_XH2PO4_ROH2_flx-H2PO4_to_XH2PO4_ROH_flx+H2PO4_e_to_HPO4_2e_flx &
    -H2PO4_1e_AlPO4_dissol_flx-H2PO4_1e_FePO4_dissol_flx-H2PO4_1e_CaHPO4_dissol_flx-3.0_r8*H2PO4_1e_apatite_dissol_flx

  NH4_1p_aque_mole_conc=NH4_1p_aque_mole_conc+RN4S
  NH3_aque_mole_conc=NH3_aque_mole_conc+RN3S
  H1PO4_2e_aque_mole_conc=H1PO4_2e_aque_mole_conc+RHP1
  H2PO4_1e_aque_mole_conc=H2PO4_1e_aque_mole_conc+RHP2
  end subroutine SolubilityEquilibriaNoSalt

end module InitSoluteMod
