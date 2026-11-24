module ChemEquilibriaMod
!
! Description
! code to solve equilibrium chemistry, without the salt equilbrium
! the solver assumes prescirbed pH, and therefore [OH(-)]
! Note:
! the model does not include H3PO4, which seems problematic, July 28, 2022
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoluteParMod
  use EcosimConst
  use SoluteChemDataType
  use EcoSIMCtrlDataType
  use EcoSIMSolverPar
  implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__
  public ::  NoSaltChemEquilibria
  contains

!------------------------------------------------------------------
  subroutine NoSaltChemEquilibria(I,J,L,NY,NX,chemvar,solflx)
!
!  Solve for chemical concentrations using input pH
!     IF NO SALTS IS SELECTED IN SITE FILE THEN A SUBSET
!     OF THE EQUILIBRIA REACTIONS ARE SOLVED: MOSTLY THOSE
!     FOR PHOSPHORUS AND CO-REACTANTS
! Reference:
! Grant R.F. and Heaney D.J., 1997, Inorganic phosphorus transformation
! and transport in soils: mathematical modeling in ecosys, Soil Sci. Soc. Am. J.
! 61: 752-764.
  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  type(solute_flx_type),target, intent(inout) :: solflx
  type(chem_var_type), target, intent(inout) :: chemvar

  real(r8) :: Al_3p_aque_mole_conc,CALX,Ca_2p_aque_mole_conc
  real(r8) :: CaX_conc,CEC_conc
  real(r8) :: H2CO3_aque_mole_conc,CO3_2e_aque_mole_conc
  real(r8) :: CFEX,Fe_3p_aque_mole_conc
  real(r8) :: H2PO4_1e_AlPO4_eqv,H2PO4_1e_CaHPO4_eqv,H2PO4_1e_FePO4_eqv,H2PO4_1e_apatite_eqv,H2PO4_1e_CaH4P2O8_eqv
  real(r8) :: H_1p_aque_mole_conc,Na_1p_aque_mole_conc,K_1p_aque_mole_conc,Mg_2p_aque_mole_conc,MgX_conc
  real(r8) :: OH_1e_aque_mole_conc,DP,FX
  real(r8) :: RB1P,RB2P,RBH1,RBH2,RH2B
  real(r8) :: RHB1,RHB2,RHP1,RHP2
  real(r8) :: RN3B,RN3S,RNH4Net2Band,RNH4Net2Soil
  real(r8) :: RNHB,H2PO4_1e_AlPO4_dissolB_flx,H2PO4_1e_CaHPO4_dissolB_flx,H2PO4_1e_apatite_dissolB_flx
  real(r8) :: H2PO4_1e_CaH4P2O8_dissolB_flx,H2PO4_1e_FePO4_dissolB_flx,RX1P,RX2P
  real(r8) :: RXH1,RXH2,RXH1B,H2PO4_1e_to_XH2PO4_ROH2_Bflx,H1PO4_to_XHPO4_ROH_flx,H2PO4_1e_to_XH2PO4_ROH2_flx,H2PO4_to_XH2PO4_ROH_flx
  real(r8) :: RXNB,H2PO4_to_XH2PO4_ROH_Bflx
  real(r8) :: S0,S1,SPH1P,SPH2P
  real(r8) :: XALQ,XCAQ,XCAX,XFEQ,XHYQ,XKAQ,XMGQ
  real(r8) :: XN4Q,XTLQ,XNAQ,XNBQ
  real(r8) :: H2PO4_e_to_HPO4_2e_flx,RNH4,H2PO4_1e_AlPO4_dissol_flx,H2PO4_1e_CaHPO4_dissol_flx
  real(r8) :: H2PO4_1e_apatite_dissol_flx,H2PO4_1e_CaH4P2O8_dissol_flx,H2PO4_1e_FePO4_dissol_flx,RXN4
  real(r8) :: VLWatMicPBK
  real(r8), parameter :: ZEROS = 1.0E-015_r8
  real(r8), parameter :: ZEROS2= 1.0E-008_r8

  real(r8), pointer :: SoilMicPMassLayerX
  real(r8), pointer :: VLWatMicPM
  real(r8), pointer :: SoilMicPMassLayer
  real(r8), pointer :: GKCA
  real(r8), pointer :: GKCH
  real(r8), pointer :: GKCM
  real(r8), pointer :: GKCK
  real(r8), pointer :: GKCN
  real(r8), pointer :: GKC4
  real(r8), pointer :: VLWatMicPNB
  real(r8), pointer :: VLWatMicPNH
  real(r8), pointer :: VLWatMicPPB
  real(r8), pointer :: VLWatMicPPO
  real(r8), pointer :: VLNH4
  real(r8), pointer :: VLNHB

  real(r8), pointer :: XCEC
  real(r8), pointer :: PH
  real(r8), pointer :: XAEC

  real(r8), pointer :: CAL
  real(r8), pointer :: CFE
  real(r8), pointer :: ZMG
  real(r8), pointer :: ZNA
  real(r8), pointer :: ZKA
  real(r8), pointer :: CO2S
  real(r8), pointer :: CCA
  real(r8), pointer :: H1PO4_2e_aque_mole_conc
  real(r8), pointer :: H1PO4_2e_band_mole_conc
  real(r8), pointer :: H2PO4_1e_aque_mole_conc
  real(r8), pointer :: H2PO4_1e_band_mole_conc
  real(r8), pointer :: NH3_aque_mole_conc
  real(r8), pointer :: NH3_aqu_band_mole_conc
  real(r8), pointer :: NH4_1p_aque_mole_conc
  real(r8), pointer :: NH4_1p_band_mole_conc
  real(r8), pointer :: XNH4_mole_conc
  real(r8), pointer :: XNH4_band_mole_conc
  real(r8), pointer :: XHPO4_band_mole_conc
  real(r8), pointer :: XH2PO4_band_mole_conc
  real(r8), pointer :: XROH_band_mole_conc
  real(r8), pointer :: XHPO4_mole_conc
  real(r8), pointer :: XROH2_band_mole_conc
  real(r8), pointer :: XH2PO4_mole_conc
  real(r8), pointer :: XROH1_mole_conc
  real(r8), pointer :: XROH2_mole_conc
  real(r8), pointer :: Precp_AlPO4_mole_conc
  real(r8), pointer :: PrecpB_AlPO4_mole_conc
  real(r8), pointer :: Precp_CaHPO4_mole_conc
  real(r8), pointer :: PrecpB_CaHPO4_mole_conc
  real(r8), pointer :: Precp_Ca5P3O12O3H3_mole_conc
  real(r8), pointer :: PrecpB_Ca5P3O12O3H3_mole_conc
  real(r8), pointer :: Precp_CaH4P2O8_mole_conc
  real(r8), pointer :: PrecpB_CaH4P2O8_mole_conc
  real(r8), pointer :: Precp_FePO4_mole_conc
  real(r8), pointer :: PrecpB_FePO4_mole_conc


  real(r8), pointer :: TRChem_NH4_soil_mole
  real(r8), pointer :: TRChem_NH4_band_soil_mole
  real(r8), pointer :: TRChem_NH3_soil_mole
  real(r8), pointer :: TRChem_NH3_band_soil
  real(r8), pointer :: TRChem_H1PO4_soil
  real(r8), pointer :: TRChem_H2PO4_soil
  real(r8), pointer :: TRChem_H1PO4_band_soil
  real(r8), pointer :: TRChem_H2PO4_band_soil
  real(r8), pointer :: TRChem_NH4_sorbed_soil
  real(r8), pointer :: TRChem_NH4_sorbed_band_soil
  real(r8), pointer :: TRChem_ROH_sorbed_soil
  real(r8), pointer :: TRChem_ROH2_sorbed_soil
  real(r8), pointer :: TRChem_RHPO4_sorbed_soil
  real(r8), pointer :: TRChem_RH2PO4_sorbed_soil
  real(r8), pointer :: TRChem_ROH_sorbed_band_soil
  real(r8), pointer :: TRChem_ROH2_sorbed_band_soil
  real(r8), pointer :: TRChem_RHPO4_sorbed_band_soil
  real(r8), pointer :: TRChem_RH2PO4_sorbed_band_soil
  real(r8), pointer :: TRChem_AlPO4_precip_soil
  real(r8), pointer :: TRChem_FePO4_precip_soil
  real(r8), pointer :: TRChem_CaHPO4_precip_soil
  real(r8), pointer :: TRChem_apatite_precip_soil
  real(r8), pointer :: TRChem_CaH4P2O8_precip_soil
  real(r8), pointer :: TRChem_AlPO4_precip_band_soil
  real(r8), pointer :: TRChem_FePO4_precip_band_soil
  real(r8), pointer :: TRChem_CaHPO4_precip_band_soil
  real(r8), pointer :: TRChem_apatite_precip_band_soil
  real(r8), pointer :: TRChem_CaH4P2O8_precip_band_soil

!     begin_execution
  call solflx%SetZero()

  VLWatMicPNB              => chemvar%VLWatMicPNB
  VLWatMicPNH              => chemvar%VLWatMicPNH
  VLWatMicPPB              => chemvar%VLWatMicPPB
  VLWatMicPPO              => chemvar%VLWatMicPPO
  XROH1_mole_conc               => chemvar%XROH1_mole_conc
  XNH4_mole_conc           => chemvar%XNH4_mole_conc
  XNH4_band_mole_conc           => chemvar%XNH4_band_mole_conc
  H1PO4_2e_band_mole_conc       => chemvar%H1PO4_2e_band_mole_conc
  H1PO4_2e_aque_mole_conc  => chemvar%H1PO4_2e_aque_mole_conc
  H2PO4_1e_aque_mole_conc  => chemvar%H2PO4_1e_aque_mole_conc
  H2PO4_1e_band_mole_conc       => chemvar%H2PO4_1e_band_mole_conc
  XHPO4_band_mole_conc          => chemvar%XHPO4_band_mole_conc
  XH2PO4_band_mole_conc         => chemvar%XH2PO4_band_mole_conc
  XROH_band_mole_conc           => chemvar%XROH_band_mole_conc
  XHPO4_mole_conc               => chemvar%XHPO4_mole_conc
  XROH2_band_mole_conc          => chemvar%XROH2_band_mole_conc
  XH2PO4_mole_conc              => chemvar%XH2PO4_mole_conc
  XROH2_mole_conc               => chemvar%XROH2_mole_conc
  NH3_aque_mole_conc       => chemvar%NH3_aque_mole_conc
  NH3_aqu_band_mole_conc        => chemvar%NH3_aqu_band_mole_conc
  NH4_1p_aque_mole_conc    => chemvar%NH4_1p_aque_mole_conc
  NH4_1p_band_mole_conc         => chemvar%NH4_1p_band_mole_conc
  Precp_AlPO4_mole_conc         => chemvar%Precp_AlPO4_mole_conc
  PrecpB_AlPO4_mole_conc        => chemvar%PrecpB_AlPO4_mole_conc
  Precp_CaHPO4_mole_conc        => chemvar%Precp_CaHPO4_mole_conc
  PrecpB_CaHPO4_mole_conc       => chemvar%PrecpB_CaHPO4_mole_conc
  Precp_Ca5P3O12O3H3_mole_conc  => chemvar%Precp_Ca5P3O12O3H3_mole_conc
  PrecpB_Ca5P3O12O3H3_mole_conc => chemvar%PrecpB_Ca5P3O12O3H3_mole_conc
  Precp_CaH4P2O8_mole_conc      => chemvar%Precp_CaH4P2O8_mole_conc
  PrecpB_CaH4P2O8_mole_conc     => chemvar%PrecpB_CaH4P2O8_mole_conc
  Precp_FePO4_mole_conc         => chemvar%Precp_FePO4_mole_conc
  PrecpB_FePO4_mole_conc         => chemvar%PrecpB_FePO4_mole_conc
  SoilMicPMassLayerX       => chemvar%SoilMicPMassLayerX
  XCEC                     => chemvar%XCEC
  PH                       => chemvar%PH
  CAL                      => chemvar%CAL
  CFE                      => chemvar%CFE
  VLWatMicPM               => chemvar%VLWatMicPM
  ZMG                      => chemvar%ZMG
  ZNA                      => chemvar%ZNA
  ZKA                      => chemvar%ZKA
  CO2S                     => chemvar%CO2S
  CCA                      => chemvar%CCA
  SoilMicPMassLayer        => chemvar%SoilMicPMassLayer
  XAEC                     => chemvar%XAEC
  VLNH4                    => chemvar%VLNH4
  GKC4                     => chemvar%GKC4
  VLNHB                    => chemvar%VLNHB
  GKCA                     => chemvar%GKCA
  GKCH                     => chemvar%GKCH
  GKCM                     => chemvar%GKCM
  GKCK                     => chemvar%GKCK
  GKCN                     => chemvar%GKCN

  TRChem_NH4_soil_mole             => solflx%TRChem_NH4_soil_mole
  TRChem_NH4_band_soil_mole        => solflx%TRChem_NH4_band_soil_mole
  TRChem_NH3_soil_mole             => solflx%TRChem_NH3_soil_mole
  TRChem_NH3_band_soil             => solflx%TRChem_NH3_band_soil
  TRChem_H1PO4_soil                => solflx%TRChem_H1PO4_soil
  TRChem_H2PO4_soil                => solflx%TRChem_H2PO4_soil
  TRChem_H1PO4_band_soil           => solflx%TRChem_H1PO4_band_soil
  TRChem_H2PO4_band_soil           => solflx%TRChem_H2PO4_band_soil
  TRChem_NH4_sorbed_soil           => solflx%TRChem_NH4_sorbed_soil
  TRChem_NH4_sorbed_band_soil      => solflx%TRChem_NH4_sorbed_band_soil
  TRChem_ROH_sorbed_soil           => solflx%TRChem_ROH_sorbed_soil
  TRChem_ROH2_sorbed_soil          => solflx%TRChem_ROH2_sorbed_soil
  TRChem_RHPO4_sorbed_soil         => solflx%TRChem_RHPO4_sorbed_soil
  TRChem_RH2PO4_sorbed_soil        => solflx%TRChem_RH2PO4_sorbed_soil
  TRChem_ROH_sorbed_band_soil      => solflx%TRChem_ROH_sorbed_band_soil
  TRChem_ROH2_sorbed_band_soil     => solflx%TRChem_ROH2_sorbed_band_soil
  TRChem_RHPO4_sorbed_band_soil    => solflx%TRChem_RHPO4_sorbed_band_soil
  TRChem_RH2PO4_sorbed_band_soil   => solflx%TRChem_RH2PO4_sorbed_band_soil
  TRChem_AlPO4_precip_soil         => solflx%TRChem_AlPO4_precip_soil
  TRChem_FePO4_precip_soil         => solflx%TRChem_FePO4_precip_soil
  TRChem_CaHPO4_precip_soil        => solflx%TRChem_CaHPO4_precip_soil
  TRChem_apatite_precip_soil       => solflx%TRChem_apatite_precip_soil
  TRChem_CaH4P2O8_precip_soil      => solflx%TRChem_CaH4P2O8_precip_soil
  TRChem_AlPO4_precip_band_soil    => solflx%TRChem_AlPO4_precip_band_soil
  TRChem_FePO4_precip_band_soil    => solflx%TRChem_FePO4_precip_band_soil
  TRChem_CaHPO4_precip_band_soil   => solflx%TRChem_CaHPO4_precip_band_soil
  TRChem_apatite_precip_band_soil  => solflx%TRChem_apatite_precip_band_soil
  TRChem_CaH4P2O8_precip_band_soil => solflx%TRChem_CaH4P2O8_precip_band_soil

!
!     PRECIPITATION-DISSOLUTION CALCULATED FROM ACTIVITIES
!     OF REACTANTS AND PRODUCTS THROUGH SOLUTIONS
!     FOR THEIR EQUILIBRIUM CONSTANTS USING CURRENT
!     ION CONCENTRATION
!
!     CEC_conc,XCEC=cation exchange concentration [mol/m3],capacity
!     BLVLX=soil mass
!     C*,Z*=solute concentration, mass
!     DP*=dissociation constant from PARAMETER above
!     SP*=solubility product from PARAMETER above
!     C*<0.0=solve for C* from equilibrium with other solutes
!  call PrintChemVar(chemvar)
! SoilMicPMassLayerX=soil mass of the layer, Mg/d2
! XCEC=mol/d-2
! CEC_conc=mol/Mg, assuming Mg is for 1 m3 water

  IF(SoilMicPMassLayerX.GT.ZEROS)THEN
    CEC_conc=AMAX1(ZERO,XCEC/SoilMicPMassLayerX)
  ELSE
    CEC_conc=ZERO
  ENDIF
  !VLWatMicPM=m3/d2
  !ZMG=mol/d2

  IF(VLWatMicPM.GT.ZEROS2)THEN
    ! aqueous Mg(2+)
    Mg_2p_aque_mole_conc=AMAX1(ZERO,ZMG/VLWatMicPM)
    ! aqueous Na(+)
    Na_1p_aque_mole_conc=AMAX1(ZERO,ZNA/VLWatMicPM)
    ! aqueous K(+)
    K_1p_aque_mole_conc=AMAX1(ZERO,ZKA/VLWatMicPM)
    ! aqueous CO2 (H2CO3), mol/m3
    H2CO3_aque_mole_conc=AMAX1(ZERO,CO2S/(Catomw*VLWatMicPM))
    !volmetric rescaling factor
    VLWatMicPBK=AMIN1(1.0_r8,SoilMicPMassLayer/VLWatMicPM)
  ELSE
    Mg_2p_aque_mole_conc = 0._r8
    Na_1p_aque_mole_conc = 0._r8
    K_1p_aque_mole_conc  = 0._r8
    H2CO3_aque_mole_conc = 0._r8
    VLWatMicPBK          = 1._r8
  ENDIF


  ![H(+)] is prescribed, mol/m3, because pH is defined as log([H(+)]), with [H(+)] of unit moles per liter
  H_1p_aque_mole_conc=AMAX1(ZERO,10.0_r8**(-(PH-3.0_r8)))
  ! [OH(-)]
  OH_1e_aque_mole_conc=AMAX1(ZERO,DPH2O/H_1p_aque_mole_conc)
  ! Aluminum,[Al(3+)], Al(3+)+3OH(-)<->Al(OH)3
  IF(CAL.LT.0.0_r8)THEN
    Al_3p_aque_mole_conc=AMAX1(ZERO,SPALO/OH_1e_aque_mole_conc**3._r8)
  ELSE
    Al_3p_aque_mole_conc=AMAX1(ZERO,AMIN1(CAL,SPALO/OH_1e_aque_mole_conc**3._r8))
  ENDIF
! Fe(3+) + 3OH(-)<->Fe(OH)3
  IF(CFE.LT.0.0_r8)THEN
    Fe_3p_aque_mole_conc=AMAX1(ZERO,SPFEO/OH_1e_aque_mole_conc**3._r8)
  ELSE
    Fe_3p_aque_mole_conc=AMAX1(ZERO,AMIN1(CFE,SPFEO/OH_1e_aque_mole_conc**3._r8))
  ENDIF
!
!     CA CONCENTRATION FROM CURRENT CO2 CONCENTRATION
! H2CO3<-> CO3(2-)+2H(+) 
! ([CO3(2-)][H+]^2)/[H2CO3]=K
  CO3_2e_aque_mole_conc=AMAX1(ZERO,H2CO3_aque_mole_conc*DPCO3/H_1p_aque_mole_conc**2._r8)

! [Ca(2+)], CaCO3<->Ca(2+)+CO3(2-)
! [Ca(2+)][CO3(2-)]/[CaCO3(solid)]=K (SPCAC)
! 
  IF(CCA.LT.0.0_r8)THEN
    Ca_2p_aque_mole_conc=AMAX1(ZERO,AMIN1(CCAMX,SPCAC/CO3_2e_aque_mole_conc))
  ELSE
    Ca_2p_aque_mole_conc=AMAX1(ZERO,AMIN1(CCA,SPCAC/CO3_2e_aque_mole_conc))
  ENDIF

!
!     PHOSPHORUS TRANSFORMATIONS IN NON-BAND SOIL ZONE
!
  IF(VLWatMicPPO.GT.ZEROS2)THEN
!
!     ALUMINUM PHOSPHATE (VARISCITE, AlPO4)
!
!     H2PO4_1e_AlPO4_eqv,H2PO4_1e_aque_mole_conc=equilibrium, current H2PO4 concentration in non-band
!     SYA0P2=solubility product derived from SPALO
!     H2PO4_1e_AlPO4_dissol_flx=H2PO4 dissolution from AlPO4 in non-band
!
!   PALP01: precipitated Al(PO4)
!   AlPO4.(H2O)2<-> 2OH(-)+H2PO4(-)+Al(3+) 
!   ([OH(-)]^2[Al(3+)][H2PO4(-)])/[AlPO4]=K
    H2PO4_1e_AlPO4_eqv        = SYA0P2/(Al_3p_aque_mole_conc*OH_1e_aque_mole_conc**2._r8)
    H2PO4_1e_AlPO4_dissol_flx = AMAX1(-Precp_AlPO4_mole_conc,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_AlPO4_eqv))
!------------------------------------------
!     IRON PHOSPHATE (STRENGITE)
!
!     H2PO4_1e_FePO4_eqv,H2PO4_1e_aque_mole_conc=equilibrium,current H2PO4 concentration in non-band
!     SYF0P2=solubility product derived from SPALO
!     H2PO4_1e_FePO4_dissol_flx=H2PO4 dissolution from FePO4 in non-band
!
!   PFEP01: precipitated FePO4
!   FePO4.(H2O)2 <-> Fe(3+)+2OH(-) + H2PO4(-)
    H2PO4_1e_FePO4_eqv=SYF0P2/(Fe_3p_aque_mole_conc*OH_1e_aque_mole_conc**2._r8)
    H2PO4_1e_FePO4_dissol_flx=AMAX1(-Precp_FePO4_mole_conc,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_FePO4_eqv))
!------------------------------------------
!     DICALCIUM PHOSPHATE (CaHPO4)
!
!     H2PO4_1e_CaHPO4_eqv,H2PO4_1e_aque_mole_conc=equilibrium,current H2PO4 concentration in non-band
!     SYCAD2=solubility product derived from SPALO
!     H2PO4_1e_CaHPO4_dissol_flx=H2PO4 dissolution from CaHPO4 in non-band
!  Precp_CaHPO4_mole_conc: precipitated CaHPO4
!  CaHPO4 + H2O <-> Ca(2+) + OH(-)+H2PO4(-)
    H2PO4_1e_CaHPO4_eqv=SYCAD2/(Ca_2p_aque_mole_conc*OH_1e_aque_mole_conc)
    H2PO4_1e_CaHPO4_dissol_flx=AMAX1(-Precp_CaHPO4_mole_conc,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_CaHPO4_eqv))
!------------------------------------------
!     'HYDROXYAPATITE'
!
!     H2PO4_1e_apatite_eqv,H2PO4_1e_aque_mole_conc=equilibrium,current H2PO4 concentration in non-band
!     SYCAH2=solubility product derived from SPALO
!     H2PO4_1e_apatite_dissol_flx=H2PO4 dissolution from apatite in non-band
!
!  Precp_Ca5P3O12O3H3_mole_conc: precipitated Ca5(PO4)3OH
!  Ca5(PO4)3OH+6H2O<->5Ca(2+) +7OH(-)+3H2PO4(-)
!
    H2PO4_1e_apatite_eqv=(SYCAH2/(Ca_2p_aque_mole_conc**5._r8*OH_1e_aque_mole_conc**7))**0.333_r8
    H2PO4_1e_apatite_dissol_flx=AMAX1(-Precp_Ca5P3O12O3H3_mole_conc,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_apatite_eqv))
!
!------------------------------------------
!    'MONOCALCIUM PHOSPHATE'
!
!     H2PO4_1e_CaH4P2O8_eqv,H2PO4_1e_aque_mole_conc=equilibrium,current H2PO4 concentration in non-band
!     SPCAM=solubility product for Ca(H2PO4)2
!     H2PO4_1e_CaH4P2O8_dissol_flx=H2PO4 dissolution from Ca(H2PO4)2 in non-band
!   Precp_CaH4P2O8_mole_conc: precipitated Ca(H2PO4)2
!   Ca(H2PO4)2 <-> Ca(2+)+2H2PO4(-)
    H2PO4_1e_CaH4P2O8_eqv=SQRT(SPCAM/Ca_2p_aque_mole_conc)
    H2PO4_1e_CaH4P2O8_dissol_flx=AMAX1(-Precp_CaH4P2O8_mole_conc*SPPO4,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_CaH4P2O8_eqv))

!------------------------------------------
!
!     PHOSPHORUS ANION EXCHANGE IN NON-BAND SOIL ZONE
!     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
!     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
!     EXCHANGE SITES
!
!   'XAEC: anaion exchange capacity'
    IF(XAEC.GT.ZEROS)THEN
!     Anion adsorption equilibria
!     H2PO4 EXCHANGE IN NON-BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG H2PO4-, H+, OH-, X-OH
!     AND X-H2PO4
!
!     SPH2P,SXH2P=equilibrium constant for H2PO4 exchange with R-OH2,R-OH
!     H2PO4_1e_to_XH2PO4_ROH2_flx,H2PO4_to_XH2PO4_ROH_flx=H2PO4 exchange with R-OH2,R-OH in non-band
!     TADA: rate constant for adsorption equilibria
!     XROH2_mole_conc: X-OH2(+)
!     XROH1_mole_conc: X-OH
!     X-H2PO4+H2O   <-> X-OH2(+)+H2PO4(-)
!     [X-OH2(+)][H2PO4(-)]/[X-H2PO4]=K, SPH2P
!     X-H2PO4+OH(-) <-> X-OH+H2PO4(-)
!     [X-OH][H2PO4(-)]/([X-H2PO4][OH(-)])=K,SXH2P

      ! R-H2PO4 + H2O <=> R-OH2(+) + H2PO4(-).
      ! Combine neutral-site exchange (SXH2P), water dissociation (DPH2O),
      ! and R-OH2(+) <=> R-OH + H(+) (SXOH2) to close the equilibrium cycle.
      SPH2P = SXH2P*DPH2O/SXOH2
      H2PO4_1e_to_XH2PO4_ROH2_flx = TADA*(XROH2_mole_conc*H2PO4_1e_aque_mole_conc-SPH2P*XH2PO4_mole_conc)/(XROH2_mole_conc+SPH2P)*VLWatMicPBK
      H2PO4_to_XH2PO4_ROH_flx     = TADA*(XROH1_mole_conc*H2PO4_1e_aque_mole_conc-SXH2P*OH_1e_aque_mole_conc*XH2PO4_mole_conc)/(XROH1_mole_conc+SXH2P)*VLWatMicPBK
!------------------------------------------
!     HPO4 EXCHANGE IN NON-BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG HPO4--, H+, OH-, X-OH
!     AND X-HPO4
!
!     SPH1P=equilibrium constant for HPO4 exchange with R-OH
!     H1PO4_to_XHPO4_ROH_flx=HPO4 exchange with R-OH in non-band
!     X-HPO4(-)+ OH(-)<-> X-OH+HPO4(2-)
!     [X-OH][HPO4(2-)]/([X-HPO4(-)][OH(-)])=K
      SPH1P                  = SXH1P*DPH2O/DPH2P
      H1PO4_to_XHPO4_ROH_flx = TADA*(XROH1_mole_conc*H1PO4_2e_aque_mole_conc-SPH1P*XHPO4_mole_conc)/(XROH1_mole_conc+SPH1P)*VLWatMicPBK
    ELSE
      H2PO4_1e_to_XH2PO4_ROH2_flx = 0._r8
      H2PO4_to_XH2PO4_ROH_flx     = 0._r8
      H1PO4_to_XHPO4_ROH_flx      = 0._r8
    ENDIF
!------------------------------------------    
!
!    'H2PO4(-) <-> H(+)+HPO4(2-)'
!
!     DPH2P=dissociation constant
!     S1=equilibrium concentration in non-band
!     H2PO4_e_to_HPO4_2e_flx=H2PO4- <-> H(+) + HPO4(--) dissociation in non-band
!   The problem is equivalent to solving a quadratic equation, with S1 being the \Delta
!   of the following equation
!    x^2+(A+B+K)x+AB-KC=0
!    corresponding to
!  (AB) <-> A+B,  (A+x)(B+x)/((AB)-x)=K
!  X represents flux from H2PO4(-) dissolution
    ! Rationalize the smaller root to retain small phosphate transfers;
    ! TSL retains the existing relaxation and root factor 1/2.
    DP                     = DPH2P
    S0                     = H1PO4_2e_aque_mole_conc+H_1p_aque_mole_conc+DP
    S1=MAX(0._r8,(H1PO4_2e_aque_mole_conc-H_1p_aque_mole_conc)**2+2._r8*DP*(H1PO4_2e_aque_mole_conc+H_1p_aque_mole_conc)+DP**2+4._r8*DP*H2PO4_1e_aque_mole_conc)
    H2PO4_e_to_HPO4_2e_flx=TSL*(4._r8*((H1PO4_2e_aque_mole_conc*H_1p_aque_mole_conc-DP*H2PO4_1e_aque_mole_conc)/(S0+SQRT(S1))))
  ELSE
    H2PO4_1e_AlPO4_dissol_flx    = 0._r8
    H2PO4_1e_FePO4_dissol_flx    = 0._r8
    H2PO4_1e_CaHPO4_dissol_flx   = 0._r8
    H2PO4_1e_apatite_dissol_flx  = 0._r8
    H2PO4_1e_CaH4P2O8_dissol_flx = 0._r8
    H2PO4_1e_to_XH2PO4_ROH2_flx  = 0._r8
    H2PO4_to_XH2PO4_ROH_flx      = 0._r8
    H1PO4_to_XHPO4_ROH_flx       = 0._r8
    H2PO4_e_to_HPO4_2e_flx       = 0._r8
  ENDIF

!------------------------------------------
!     PHOSPHORUS PRECIPITATION-DISSOLUTION IN BAND SOIL ZONE
!
  IF(VLWatMicPPB.GT.ZEROS2)THEN
!------------------------------------------
!     ALUMINUM PHOSPHATE (VARISCITE)
!
!     H2PO4_1e_AlPO4_eqv,H2PO4_1e_band_mole_conc=equilibrium,current H2PO4 concentration in band
!     SYA0P2=solubility product derived from SPALO
!     H2PO4_1e_AlPO4_dissolB_flx=H2PO4 dissolution from AlPO4 in band
! Al(3+) and OH(-) are shared between band and non-band soils
!
    H2PO4_1e_AlPO4_eqv=SYA0P2/(Al_3p_aque_mole_conc*OH_1e_aque_mole_conc**2)
    H2PO4_1e_AlPO4_dissolB_flx=AMAX1(-PrecpB_AlPO4_mole_conc,TPD*(H2PO4_1e_band_mole_conc-H2PO4_1e_AlPO4_eqv))
!------------------------------------------
!     IRON PHOSPHATE (STRENGITE)
!
!     H2PO4_1e_FePO4_eqv,H2PO4_1e_band_mole_conc=equilibrium,current H2PO4 concentration in band
!     SYF0P2=solubility product derived from SPALO
!     H2PO4_1e_FePO4_dissolB_flx=H2PO4 dissolution from FePO4 in band
!   Fe(3+) is shared between band and non-band soils
    H2PO4_1e_FePO4_eqv=SYF0P2/(Fe_3p_aque_mole_conc*OH_1e_aque_mole_conc**2._r8)
    H2PO4_1e_FePO4_dissolB_flx=AMAX1(-PrecpB_FePO4_mole_conc,TPD*(H2PO4_1e_band_mole_conc-H2PO4_1e_FePO4_eqv))
!------------------------------------------
!     'DICALCIUM PHOSPHATE'
!
!     H2PO4_1e_CaHPO4_eqv,H2PO4_1e_band_mole_conc=equilibrium,current H2PO4 concentration in band
!     SYCAD2=solubility product derived from SPALO
!     H2PO4_1e_CaHPO4_dissolB_flx=H2PO4 dissolution from CaHPO4 in band
!
    H2PO4_1e_CaHPO4_eqv         = SYCAD2/(Ca_2p_aque_mole_conc*OH_1e_aque_mole_conc)
    H2PO4_1e_CaHPO4_dissolB_flx = AMAX1(-PrecpB_CaHPO4_mole_conc,TPD*(H2PO4_1e_band_mole_conc-H2PO4_1e_CaHPO4_eqv))
!------------------------------------------
!     HYDROXYAPATITE
!
!     H2PO4_1e_apatite_eqv,H2PO4_1e_band_mole_conc=equilibrium,current H2PO4 concentration in band
!     SYCAH2=solubility product derived from SPALO
!     H2PO4_1e_apatite_dissolB_flx=H2PO4 dissolution from apatite in band
!
    H2PO4_1e_apatite_eqv=(SYCAH2/(Ca_2p_aque_mole_conc**5*OH_1e_aque_mole_conc**7))**0.333_r8
    H2PO4_1e_apatite_dissolB_flx=AMAX1(-PrecpB_Ca5P3O12O3H3_mole_conc,TPD*(H2PO4_1e_band_mole_conc-H2PO4_1e_apatite_eqv))
!------------------------------------------
!     'MONOCALCIUM PHOSPHATE'
!
!     H2PO4_1e_CaH4P2O8_eqv,H2PO4_1e_band_mole_conc=equilibrium,current H2PO4 concentration in band
!     SPCAM=solubility product for Ca(H2PO4)2
!     H2PO4_1e_CaH4P2O8_dissolB_flx=H2PO4 dissolution from Ca(H2PO4)2 in band
!
    H2PO4_1e_CaH4P2O8_eqv         = SQRT(SPCAM/Ca_2p_aque_mole_conc)
    H2PO4_1e_CaH4P2O8_dissolB_flx = AMAX1(-PrecpB_CaH4P2O8_mole_conc*SPPO4,TPD*(H2PO4_1e_band_mole_conc-H2PO4_1e_CaH4P2O8_eqv))
!     IF(I.GT.315)THEN
!     WRITE(*,1117)'RPPOB',I,J,L,H2PO4_1e_CaH4P2O8_dissolB_flx,H2PO4_1e_CaH4P2O8_eqv,H2PO4_1e_band_mole_conc,SPCAM,Ca_2p_aque_mole_conc
!    2,PrecpB_CaH4P2O8_mole_conc,SPPO4,TPD,PCPMB,BKVLPB
!     ENDIF
!------------------------------------------
!     PHOSPHORUS ANION EXCHANGE IN BAND SOIL ZONE
!     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
!     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
!     EXCHANGE SITES
!  Anion exchange concentration is also shared among band and non-band soils

    IF(XAEC.GT.ZEROS)THEN
!
!     H2PO4 EXCHANGE IN BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG H2PO4-, H+, OH-, X-OH
!     AND X-H2PO4
!
!     SPH2P,SXH2P=equilibrium constant for H2PO4 exchange with R-OH2,R-OH
!     H2PO4_1e_to_XH2PO4_ROH2_Bflx,H2PO4_to_XH2PO4_ROH_Bflx=H2PO4 exchange with R-OH2,R-OH in band
!     R-H2PO4 <-> R-OH2+H2PO4(-)
!   H2PO4_1e_band_mole_conc: H2PO4(-)
!   XH2PO4_band_mole_conc: R-H2PO4
!   XROH2_band_mole_conc: R-OH2
      ! R-H2PO4 + H2O <=> R-OH2(+) + H2PO4(-).
      ! Combine neutral-site exchange (SXH2P), water dissociation (DPH2O),
      ! and R-OH2(+) <=> R-OH + H(+) (SXOH2) to close the equilibrium cycle.
      SPH2P = SXH2P*DPH2O/SXOH2
      H2PO4_1e_to_XH2PO4_ROH2_Bflx = TADA*(XROH2_band_mole_conc*H2PO4_1e_band_mole_conc-SPH2P*XH2PO4_band_mole_conc)/(XROH2_band_mole_conc+SPH2P)*VLWatMicPBK
      H2PO4_to_XH2PO4_ROH_Bflx     = TADA*(XROH_band_mole_conc*H2PO4_1e_band_mole_conc-SXH2P*XH2PO4_band_mole_conc*OH_1e_aque_mole_conc)/(XROH_band_mole_conc+SXH2P)*VLWatMicPBK
!
!     HPO4 EXCHANGE IN BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG HPO4--, H+, OH-, X-OH
!     AND X-HPO4
!
!     SPH1P=equilibrium constant for HPO4 exchange with R-OH
!     RXH1B=HPO4 exchange with R-OH in band
!
      SPH1P = SXH1P*DPH2O/DPH2P
      RXH1B = TADA*(XROH_band_mole_conc*H1PO4_2e_band_mole_conc-SPH1P*XHPO4_band_mole_conc)/(XROH_band_mole_conc+SPH1P)*VLWatMicPBK
    ELSE
      H2PO4_1e_to_XH2PO4_ROH2_Bflx = 0._r8
      H2PO4_to_XH2PO4_ROH_Bflx     = 0._r8
      RXH1B                        = 0._r8
    ENDIF
!
!     H2PO4-H+HPO4
!
!     DPH2P=dissociation constant
!     S1=eqilibriunm concentration in band
!     RH2B=H2PO4-H+HPO4 dissociation in band
!
    ! Use the same cancellation-resistant root in the fertilizer band.
    DP   = DPH2P
    S0   = H1PO4_2e_band_mole_conc+H_1p_aque_mole_conc+DP
    S1=MAX(0._r8,(H1PO4_2e_band_mole_conc-H_1p_aque_mole_conc)**2+2._r8*DP*(H1PO4_2e_band_mole_conc+H_1p_aque_mole_conc)+DP**2+4._r8*DP*H2PO4_1e_band_mole_conc)
    RH2B=TSLX*(4._r8*((H1PO4_2e_band_mole_conc*H_1p_aque_mole_conc-DP*H2PO4_1e_band_mole_conc)/(S0+SQRT(S1))))
  ELSE
    H2PO4_1e_AlPO4_dissolB_flx    = 0._r8
    H2PO4_1e_FePO4_dissolB_flx    = 0._r8
    H2PO4_1e_CaHPO4_dissolB_flx   = 0._r8
    H2PO4_1e_apatite_dissolB_flx  = 0._r8
    H2PO4_1e_CaH4P2O8_dissolB_flx = 0._r8
    H2PO4_1e_to_XH2PO4_ROH2_Bflx  = 0._r8
    H2PO4_to_XH2PO4_ROH_Bflx      = 0._r8
    RXH1B                         = 0._r8
    RH2B                          = 0._r8
  ENDIF
!
!------------------------------------------------------
!     CATION EXCHANGE FROM GAPON SELECTIVITY COEFFICIENTS
!     FOR CA-NH4, CA-H, CA-AL
! cation-exchange capacity
  IF(XCEC.GT.ZEROS)THEN
!
!     CATION CONCENTRATIONS
!
!     EQUILIBRIUM X-CA CONCENTRATION FROM CEC, GAPON COEFFICIENTS
!     AND CATION CONCENTRATIONS
!
!     CEC_conc,XCEC=cation exchange concentration,capacity
!     XCAX=equilibrium R-Ca concentration
!     GKC4,GKCH,GKCA,GKCM,GKCN,GKCK=Gapon selectivity coefficients for
!     CA-NH4,CA-H,CA-AL,CA-MG,CA-NA,CA-K
!     X*Q=equilibrium exchangeable concentrations
!     XTLQ=total equilibrium exchangeable concentration
!
!   non-band X-NH4 + H(+) <-> X-H + NH4(+)
    NH4_1p_aque_mole_conc=AMAX1(ZERO,NH4_1p_aque_mole_conc)

!   band X-NH4 + H(+) <-> X-H + NH4(+)
    NH4_1p_band_mole_conc=AMAX1(ZERO,NH4_1p_band_mole_conc)

!   X-Al(1/3) + H(+) <-> X-H + (1/3)Al(3+)
    CALX=AMAX1(ZERO,Al_3p_aque_mole_conc)**0.333_r8

!   X-Fe(1/3) + H(+) <-> X-H + (1/3)Fe(3+)
    CFEX=AMAX1(ZERO,Fe_3p_aque_mole_conc)**0.333_r8

!   X-Ca(1/2) + H(+) <-> X-H + (1/2)Ca(2+)
    CaX_conc=AMAX1(ZERO,Ca_2p_aque_mole_conc)**0.50_r8

!   X-Mg(1/2) + H(+) <-> X-H + (1/2)Mg(2+)
    MgX_conc=AMAX1(ZERO,Mg_2p_aque_mole_conc)**0.50_r8

!   X-Na + H(+) <-> X-H + Na(+)
    Na_1p_aque_mole_conc=AMAX1(ZERO,Na_1p_aque_mole_conc)

!   X-K + H(+) <-> X-H + K(+)
    K_1p_aque_mole_conc=AMAX1(ZERO,K_1p_aque_mole_conc)
!
!   EQUILIBRIUM X-CA CONCENTRATION FROM CEC AND CATION
!   CONCENTRATIONS
!
!   The following equation can be found in Robbins et al. (1980)
!   XCAX designates X-Ca
!   one X corresponds to one e.
!   Given a general reaction

!   X-H + 1/n A(n+) <-> X-A[(n+)/n] + H(+)
!   X-A[(n+)/n]*H(+)/[(A)^(1/n)*(X-H)]=GK
!   X-A[(n+)/n]=GK*(A)^(1/n)*(X-H)

    XCAX=CEC_conc/(1.0+(GKC4*NH4_1p_aque_mole_conc*VLNH4 &
      +GKC4*NH4_1p_band_mole_conc*VLNHB &
      +GKCH*H_1p_aque_mole_conc &
      +GKCA*CALX+GKCA*CFEX+GKCM*MgX_conc &
      +GKCN*Na_1p_aque_mole_conc &
      +GKCK*K_1p_aque_mole_conc)/CaX_conc)

!   X-H+NH4(+) <-> X-NH4 + H(+)
!   XN4Q: equilibrium exchangeable NH4(+) non-band soil
    XN4Q=XCAX*NH4_1p_aque_mole_conc*GKC4

!   XNBq: equilibrium exchangeable NH4(+) band soil
    XNBQ=XCAX*NH4_1p_band_mole_conc*GKC4

!   XHYQ: equilibrium exchangeable H(+) in soil
    XHYQ=XCAX*H_1p_aque_mole_conc*GKCH

!   3X-H+Al(3+) <-> X-Al+ 3H(+)
!   XALQ: equilibrium exchangeable Al(3+) in soil
    XALQ=XCAX*CALX*GKCA

!   3X-H+Fe(3+) <-> X-Fe + 3(H+)
!   XFEQ: equilibrium exchangeable Fe(3+) soil

    XFEQ = XCAX*CFEX*GKCA
    XCAQ = XCAX*CaX_conc
    XMGQ = XCAX*MgX_conc*GKCM
    XNAQ = XCAX*Na_1p_aque_mole_conc*GKCN
    XKAQ = XCAX*K_1p_aque_mole_conc*GKCK
    XTLQ = XN4Q*VLNH4+XNBQ*VLNHB+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ

    IF(XTLQ.GT.ZERO)THEN
      FX=CEC_conc/XTLQ
    ELSE
      FX=0._r8
    ENDIF
    XN4Q = FX*XN4Q
    XNBQ = FX*XNBQ
!
!  'NH4 EXCHANGE IN NON-BAND AND BAND SOIL ZONES'
!
!     RXN4,RXNB=NH4 adsorption in non-band,band
!     TADC=adsorption rate constant
!  if RNX4>0, XNH4_mole_conc will be increased
    RXN4 = TADC*AMAX1(AMIN1((XN4Q-XNH4_mole_conc)*NH4_1p_aque_mole_conc/XN4Q,NH4_1p_aque_mole_conc),-XNH4_mole_conc)
    RXNB = TADC*AMAX1(AMIN1((XNBQ-XNH4_band_mole_conc)*NH4_1p_band_mole_conc/XNBQ,NH4_1p_band_mole_conc),-XNH4_band_mole_conc)
  ELSE
    RXN4=0._r8
    RXNB=0._r8
  ENDIF
!----------------------------------------------------------------------------------
!     NH4-NH3+H IN NON-BAND AND BAND SOIL ZONES
!
! RNH4,RNHB=NH4-NH3+H dissociation in non-band,band
! DPN4=NH4 dissociation constant
! NH4(+) <-> NH3 + H(+)
! K=[NH3]*[H(+)]/[NH4(+)]
! H(+): H_1p_aque_mole_conc
! NH3_aque_mole_conc: NH3
! NH4_1p_aque_mole_conc: NH4(+)
! DPN4: K
! Rate function:=[NH3]*[H(+)]-K*[NH4(+)], >0, form more NH4(+), =0 equilbrium, <0, more formation of NH3
! RNH4: if RNH4 > 0., then the equilibrium should move to association, and NH4(+) will be increased
!
  IF(VLWatMicPNH.GT.ZEROS2)THEN
    RNH4 = (H_1p_aque_mole_conc*NH3_aque_mole_conc-DPN4*NH4_1p_aque_mole_conc)/(DPN4+H_1p_aque_mole_conc)
    ! Prescribed pH supplies H(+); only NH3 and NH4 are finite N donors.
    ! Match the band limiter and retain equal, opposite NH3/NH4 transfers.
    RNH4 = AMAX1(AMIN1(RNH4,NH3_aque_mole_conc),-NH4_1p_aque_mole_conc)

  ELSE
    RNH4=0._r8
  ENDIF

  IF(VLWatMicPNB.GT.ZEROS2)THEN
    RNHB = (H_1p_aque_mole_conc*NH3_aqu_band_mole_conc-DPN4*NH4_1p_band_mole_conc)/(DPN4+H_1p_aque_mole_conc)
    RNHB = AMAX1(AMIN1(RNHB,NH3_aqu_band_mole_conc),-NH4_1p_band_mole_conc)
  ELSE
    RNHB=0._r8
  ENDIF
!
!     TOTAL ION FLUXES FOR ALL REACTIONS ABOVE
!
!     RNH4Net2Soil,RNH4Net2Band=net NH4 flux in non-band,band
!     RN3S,RN3B=net NH3 flux in non-band,band
!     RAL,RFE,RHY,RCA,RMG,RNA,RKA,ROH=net Al,Fe,H,Ca,Mg,Na,K,OH flux
!     RHP1,RHP2=net HPO4,H2PO4 flux in non-band
!     RXH1,RXH2,RX1P,RX2P=net R-OH,R-OH2,R-HPO4,R-H2PO4 in non-band
!     RHB1,RHB2=net HPO4,H2PO4 flux in band
!     RBH1,RBH2,RB1P,RB2P=net R-OH,R-OH2,R-HPO4,R-H2PO4 in band
!

  RNH4Net2Soil = RNH4-RXN4
  RNH4Net2Band = RNHB-RXNB
  RN3S         = -RNH4
  RN3B         = -RNHB
!----------------------------------------------------------------------------------  
  RHP1 = -H2PO4_e_to_HPO4_2e_flx-H1PO4_to_XHPO4_ROH_flx
  RHP2 = H2PO4_e_to_HPO4_2e_flx-H2PO4_1e_to_XH2PO4_ROH2_flx-H2PO4_to_XH2PO4_ROH_flx-H2PO4_1e_AlPO4_dissol_flx &
    -H2PO4_1e_FePO4_dissol_flx-H2PO4_1e_CaHPO4_dissol_flx &
    -2.0_r8*H2PO4_1e_CaH4P2O8_dissol_flx-3.0_r8*H2PO4_1e_apatite_dissol_flx
  RHB1 = -RH2B-RXH1B
  RHB2 = RH2B-H2PO4_1e_to_XH2PO4_ROH2_Bflx-H2PO4_to_XH2PO4_ROH_Bflx-H2PO4_1e_AlPO4_dissolB_flx-H2PO4_1e_FePO4_dissolB_flx &
    -H2PO4_1e_CaHPO4_dissolB_flx-2.0_r8*H2PO4_1e_CaH4P2O8_dissolB_flx-3.0_r8*H2PO4_1e_apatite_dissolB_flx
  RXH1 = -H2PO4_to_XH2PO4_ROH_flx-H1PO4_to_XHPO4_ROH_flx
  RXH2 = -H2PO4_1e_to_XH2PO4_ROH2_flx
  RX1P = H1PO4_to_XHPO4_ROH_flx
  RX2P = H2PO4_1e_to_XH2PO4_ROH2_flx+H2PO4_to_XH2PO4_ROH_flx
  RBH1 = -H2PO4_to_XH2PO4_ROH_Bflx-RXH1B
  RBH2 = -H2PO4_1e_to_XH2PO4_ROH2_Bflx
  RB1P = RXH1B
  RB2P = H2PO4_1e_to_XH2PO4_ROH2_Bflx+H2PO4_to_XH2PO4_ROH_Bflx
!
!     CONVERT TOTAL ION FLUXES FROM CHANGES IN CONCENTRATION
!     TO CHANGES IN MASS PER UNIT AREA FOR USE IN 'REDIST'
!
!     TRChem_NH4_soil_mole,TRChem_NH4_band_soil_mole=total NH4 flux in non-band,band
!     TRChem_NH3_soil_mole,TRChem_NH3_band_soil=total NH3 flux in non-band,band
!     TRChem_H1PO4_soil,TRChem_H2PO4_soil=net HPO4,H2PO4 flux in non-band
!     TRChem_H1PO4_band_soil,TRChem_H2PO4_band_soil=net HPO4,H2PO4 flux in band
!     TRNX4,TRNXB=total NH4 adsorption in non-band,band
!     TRChem_ROH_sorbed_soil,TRChem_ROH2_sorbed_soil,TRChem_RHPO4_sorbed_soil,TRChem_RH2PO4_sorbed_soil
!     =total R-OH,R-OH2,R-HPO4,R-H2PO4 adsorption in non-band
!     TRChem_ROH_sorbed_band_soil,TRChem_ROH2_sorbed_band_soil,TRChem_RHPO4_sorbed_band_soil,TRChem_RH2PO4_sorbed_band_soil
!     =total R-OH,R-OH2,R-HPO4,R-H2PO4 adsorption in band
!     TRChem_AlPO4_precip_soil,TRChem_FePO4_precip_soil,TRChem_CaHPO4_precip_soil,TRChem_apatite_precip_soil,TRChem_CaH4P2O8_precip_soil
!     =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation in non-band
!     TRChem_AlPO4_precip_band_soil,TRChem_FePO4_precip_band_soil,TRChem_CaHPO4_precip_band_soil,TRChem_apatite_precip_band_soil,TRChem_CaH4P2O8_precip_band_soil
!     =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation in band
!

  TRChem_NH4_soil_mole        = TRChem_NH4_soil_mole+RNH4Net2Soil*VLWatMicPNH
  TRChem_NH3_soil_mole        = TRChem_NH3_soil_mole+RN3S*VLWatMicPNH
  TRChem_NH4_sorbed_soil      = TRChem_NH4_sorbed_soil+RXN4*VLWatMicPNH
  TRChem_NH4_band_soil_mole   = TRChem_NH4_band_soil_mole+RNH4Net2Band*VLWatMicPNB
  TRChem_NH3_band_soil        = TRChem_NH3_band_soil+RN3B*VLWatMicPNB
  TRChem_NH4_sorbed_band_soil = TRChem_NH4_sorbed_band_soil+RXNB*VLWatMicPNB

  
  TRChem_H1PO4_soil         = TRChem_H1PO4_soil+RHP1*VLWatMicPPO
  TRChem_H2PO4_soil         = TRChem_H2PO4_soil+RHP2*VLWatMicPPO
  TRChem_ROH_sorbed_soil    = TRChem_ROH_sorbed_soil+RXH1*VLWatMicPPO
  TRChem_ROH2_sorbed_soil   = TRChem_ROH2_sorbed_soil+RXH2*VLWatMicPPO
  TRChem_RHPO4_sorbed_soil  = TRChem_RHPO4_sorbed_soil+RX1P*VLWatMicPPO
  TRChem_RH2PO4_sorbed_soil = TRChem_RH2PO4_sorbed_soil+RX2P*VLWatMicPPO

  TRChem_AlPO4_precip_soil    = TRChem_AlPO4_precip_soil+H2PO4_1e_AlPO4_dissol_flx*VLWatMicPPO
  TRChem_FePO4_precip_soil    = TRChem_FePO4_precip_soil+H2PO4_1e_FePO4_dissol_flx*VLWatMicPPO
  TRChem_CaHPO4_precip_soil   = TRChem_CaHPO4_precip_soil+H2PO4_1e_CaHPO4_dissol_flx*VLWatMicPPO
  TRChem_apatite_precip_soil  = TRChem_apatite_precip_soil+H2PO4_1e_apatite_dissol_flx*VLWatMicPPO
  TRChem_CaH4P2O8_precip_soil = TRChem_CaH4P2O8_precip_soil+H2PO4_1e_CaH4P2O8_dissol_flx*VLWatMicPPO

  TRChem_H1PO4_band_soil         = TRChem_H1PO4_band_soil+RHB1*VLWatMicPPB
  TRChem_H2PO4_band_soil         = TRChem_H2PO4_band_soil+RHB2*VLWatMicPPB
  TRChem_ROH_sorbed_band_soil    = TRChem_ROH_sorbed_band_soil+RBH1*VLWatMicPPB
  TRChem_ROH2_sorbed_band_soil   = TRChem_ROH2_sorbed_band_soil+RBH2*VLWatMicPPB
  TRChem_RHPO4_sorbed_band_soil  = TRChem_RHPO4_sorbed_band_soil+RB1P*VLWatMicPPB
  TRChem_RH2PO4_sorbed_band_soil = TRChem_RH2PO4_sorbed_band_soil+RB2P*VLWatMicPPB
  
  TRChem_AlPO4_precip_band_soil    = TRChem_AlPO4_precip_band_soil+H2PO4_1e_AlPO4_dissolB_flx*VLWatMicPPB
  TRChem_FePO4_precip_band_soil    = TRChem_FePO4_precip_band_soil+H2PO4_1e_FePO4_dissolB_flx*VLWatMicPPB
  TRChem_CaHPO4_precip_band_soil   = TRChem_CaHPO4_precip_band_soil+H2PO4_1e_CaHPO4_dissolB_flx*VLWatMicPPB
  TRChem_apatite_precip_band_soil  = TRChem_apatite_precip_band_soil+H2PO4_1e_apatite_dissolB_flx*VLWatMicPPB
  TRChem_CaH4P2O8_precip_band_soil = TRChem_CaH4P2O8_precip_band_soil+H2PO4_1e_CaH4P2O8_dissolB_flx*VLWatMicPPB

  end subroutine NoSaltChemEquilibria

!--------------------------------------------------------------------------
  subroutine PrintChemVar(chemvar)
  implicit none
  type(chem_var_type), target, intent(in) :: chemvar

  write(*,*)'XCEC'  ,chemvar%XCEC
  write(*,*)'SoilMicPMassLayerX' ,chemvar%SoilMicPMassLayerX
  write(*,*)'SoilMicPMassLayer'  ,chemvar%SoilMicPMassLayer
  write(*,*)'VLWatMicPM' ,chemvar%VLWatMicPM
  write(*,*)'VLWatMicPPO',chemvar%VLWatMicPPO
  write(*,*)'VLWatMicPPB',chemvar%VLWatMicPPB
  write(*,*)'VLWatMicPNH',chemvar%VLWatMicPNH
  write(*,*)'VLWatMicPNB',chemvar%VLWatMicPNB
  write(*,*)'VLNH4' ,chemvar%VLNH4
  write(*,*)'VLNHB' ,chemvar%VLNHB
  write(*,*)'ZMG'   ,chemvar%ZMG
  write(*,*)'ZNA'   ,chemvar%ZNA
  write(*,*)'ZKA'   ,chemvar%ZKA
  write(*,*)'XAEC'  ,chemvar%XAEC
  write(*,*)'CO2S'  ,chemvar%CO2S
  write(*,*)'PH'    ,chemvar%PH
  write(*,*)'H_1p_aque_mole_conc'  ,10.0**(-(chemvar%PH-3.0))
  write(*,*)'CAL'   ,chemvar%CAL
  write(*,*)'CFE'   ,chemvar%CFE
  write(*,*)'CCA'   ,chemvar%CCA
  write(*,*)'H2PO4_1e_aque_mole_conc' ,chemvar%H2PO4_1e_aque_mole_conc
  write(*,*)'Precp_AlPO4_mole_conc',chemvar%Precp_AlPO4_mole_conc
  write(*,*)'Precp_FePO4_mole_conc',chemvar%Precp_FePO4_mole_conc
  write(*,*)'Precp_CaHPO4_mole_conc',chemvar%Precp_CaHPO4_mole_conc
  write(*,*)'Precp_Ca5P3O12O3H3_mole_conc',chemvar%Precp_Ca5P3O12O3H3_mole_conc
  write(*,*)'Precp_CaH4P2O8_mole_conc',chemvar%Precp_CaH4P2O8_mole_conc
  write(*,*)'XROH2_mole_conc' ,chemvar%XROH2_mole_conc
  write(*,*)'XH2PO4_mole_conc' ,chemvar%XH2PO4_mole_conc
  write(*,*)'XROH1_mole_conc' ,chemvar%XROH1_mole_conc
  write(*,*)'XHPO4_mole_conc' ,chemvar%XHPO4_mole_conc
  write(*,*)'H1PO4_2e_band_mole_conc' ,chemvar%H1PO4_2e_band_mole_conc
  write(*,*)'H1PO4_2e_aque_mole_conc' ,chemvar%H1PO4_2e_aque_mole_conc
  write(*,*)'H2PO4_1e_band_mole_conc' ,chemvar%H2PO4_1e_band_mole_conc
  write(*,*)'PrecpB_AlPO4_mole_conc',chemvar%PrecpB_AlPO4_mole_conc
  write(*,*)'PrecpB_FePO4_mole_conc',chemvar%PrecpB_FePO4_mole_conc
  write(*,*)'PrecpB_CaHPO4_mole_conc',chemvar%PrecpB_CaHPO4_mole_conc
  write(*,*)'PrecpB_Ca5P3O12O3H3_mole_conc',chemvar%PrecpB_Ca5P3O12O3H3_mole_conc
  write(*,*)'PrecpB_CaH4P2O8_mole_conc',chemvar%PrecpB_CaH4P2O8_mole_conc
  write(*,*)'XROH2_band_mole_conc' ,chemvar%XROH2_band_mole_conc
  write(*,*)'XH2PO4_band_mole_conc' ,chemvar%XH2PO4_band_mole_conc
  write(*,*)'XROH_band_mole_conc' ,chemvar%XROH_band_mole_conc
  write(*,*)'XHPO4_band_mole_conc' ,chemvar%XHPO4_band_mole_conc
  write(*,*)'H1PO4_2e_band_mole_conc' ,chemvar%H1PO4_2e_band_mole_conc
  write(*,*)'NH3_aque_mole_conc'  ,chemvar%NH3_aque_mole_conc
  write(*,*)'NH3_aqu_band_mole_conc'  ,chemvar%NH3_aqu_band_mole_conc
  write(*,*)'NH4_1p_aque_mole_conc'  ,chemvar%NH4_1p_aque_mole_conc
  write(*,*)'NH4_1p_band_mole_conc'  ,chemvar%NH4_1p_band_mole_conc
  write(*,*)'XNH4_mole_conc'  ,chemvar%XNH4_mole_conc
  write(*,*)'XNH4_band_mole_conc'  ,chemvar%XNH4_band_mole_conc
!  IF(chemvar%VLWatMicPNH>0. .AND. chemvar%VLWatMicPNB>0. &
!    .AND. chemvar%VLWatMicPPO>0. .AND. chemvar%VLWatMicPPB>0.)PAUSE
  end subroutine PrintChemVar

end module ChemEquilibriaMod
