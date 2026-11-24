module SaltChemEquilibriaMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoluteChemDataType, only: chem_var_type, solute_flx_type
  use minimathmod,        only: AZMAX1
  use DebugToolMod
  use SoluteParMod
  use EcosimConst
  use EcoSIMSolverPar
  implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME = &
  __FILE__

  real(r8), parameter :: ZERO=1.0E-15_r8
  real(r8), parameter :: ZEROS = 1.0E-015_r8
  real(r8), parameter :: ZEROS2= 1.0E-008_r8

  real(r8) :: RAL,RAL1,RAL2,RAL3,RAL4,RALS
  real(r8) :: RB1P,RB2P,RBH0,RBH1,RBH2,RC0B,RC0P
  real(r8) :: RC1B,RC1P,RC2P,RCA,RCAC,RCAH,RCAS,RMGO
  real(r8) :: RC2B,RCAO,RCO2,RCO3,RF1B,RF2B,RFE
  real(r8) :: RFEO1,RFEO2,RFEO3,RFEO4,RFES
  real(r8) :: RFE1,RFE2,RFE3,RFE4,RF1P,RF2P
  real(r8) :: RHB0,RHB3,RHB2,RHCO,RHCO3,RHB1,RHHY,RXH0
  real(r8) :: RHOH,RHP0,RHP1,RHP2,RHP3,RHY,RKA,RKAS
  ! Scratch targets/correction refreshed by IonStrengthActivity each iteration.
  real(r8) :: OH_1e_equil_conc,RHOH_equil
  real(r8) :: RM1B,RH1B,RM1P,RMG
  real(r8) :: RMGC,RMGH,RMGS,RH3B
  real(r8) :: XHY1,XAl_conc,XFe_conc
  real(r8) :: XCa_conc,XMg_conc,XNa_conc,XK_conc,XHC1,XAlO2H2_conc,XFeO2H2_conc,XCOOH_mole_conc,XCOO
  real(r8) :: RN3B,RN3S,RN4B,RN4S,ROH
  real(r8) :: RNA,RNAC,RNAS,H2PO4_1e_AlPO4_dissolB_flx,RPALOX
  real(r8) :: RHAL1,RHALO1,RHALO2,RHALO3,RHALO4,RHFE1
  real(r8) :: R1,Precp_CaSO4_conc,Precp_CaCO3_conc,Precp_FeO3H3_conc,Precp_AlO3H3_conc
  real(r8) :: H0PO4_3e_band_conc,H3PO4_band_conc,FeHPO4_1p_band_conc,FeH2PO4_2p_band_conc
  real(r8) :: CaPO4_1e_band_conc,CaHPO4_band_conc,CaH4P2O8_1p_band_conc,MgHPO4_band_conc
  real(r8) :: H0PO4_3e_conc,H3PO4_conc,FeHPO4_p_conc,FeH2PO4_2p_aque_mole_conc,CaPO4_1e_con,CaHPO4_conc,CaH4P2O8_1p_aque_mole_conc,MgHPO4_conc
  real(r8) :: CaSO4_mole_conc,MgOH_1p_aque_mole_conc,MgCO3_conc,MgHCO3_1p_aque_mole_conc,MgSO4_conc,NaCO3_1e_aque_mole_conc,NaSO4_1e_aque_mole_conc,KSO4_1e_aque_mole_conc
  real(r8) :: FeOH_2p_aque_mole_conc,FeO2H2_p_mole_conc,FeO3H3_mole_conc,FeO4H4_1e_aque_mole_conc,FeSO4_1p_aque_mole_conc,CaO2H2_mole_conc,CaCO3_mole_conc,CaHCO3_1p_aque_mole_conc
  real(r8) :: SO4_2e_aque_mole_conc,Cl_e_conc,HCO3_e_mole_conc,AlOH_2p_aque_mole_conc,AlO2H2_1p_aque_mole_conc,AlO3H3_mole_conc,AlO4H4_1e_aque_mole_conc,AlSO4_1p_aque_mole_conc
  real(r8) :: NO3_1e_aque_mole_conc,NO3_1e_band_conc,Al_3p_aque_mole_conc,CEC_conc,H2CO3_aque_mole_conc,Mg_2p_aque_mole_conc,Na_1p_aque_mole_conc,Fe_3p_aque_mole_conc,H_1p_aque_mole_conc
  real(r8) :: CO3_2e_aque_mole_conc,K_1p_aque_mole_conc,RPCACX,H2PO4_e_to_HPO4_2e_flx,RNH4
  real(r8) :: RHCACO,RPCASO,RHA0P1,RHA1P1,RHA2P1,RHA3P1
  real(r8) :: RHFEO1,RHFEO2,RHFEO3,RHFEO4,RPFEOX,RHCAC3,RHCACH
  real(r8) :: RHA4P1,RHA0P2,RHA1P2,RHA2P2,RHA3P2,RHA4P2,RHF0P1
  real(r8) :: RHF1P1,RHF2P1,RHF3P1,RHF4P1,RHF0P2,RHF1P2,RHF2P2
  real(r8) :: RHF3P2,RHF4P2,RPCAD1,RHCAD2,RHCAH1,RHCAH2,RHA0B1
  real(r8) :: RHA1B1,RHA2B1,RHA3B1,RHA4B1,RHA0B2,RHA1B2,RHA2B2
  real(r8) :: RHA3B2,RHA4B2,RHF0B1,RHF1B1,RHF2B1,RHF3B1,RHF4B1
  real(r8) :: RHF0B2,RHF1B2,RHF2B2,RHF3B2,RHF4B2,RPCDB1,RHCDB2
  real(r8) :: RHCHB1,RHCHB2,RXOH2,RXOH1,RXO2B,RXO1B,RXHY,RXAL
  real(r8) :: SO4_2e_activity,OH_1e_activity,RH1P,RH3P,H2PO4_1e_AlPO4_dissol_flx,H2PO4_1e_CaHPO4_dissol_flx
  real(r8) :: H2PO4_1e_apatite_dissol_flx,H2PO4_1e_CaH4P2O8_dissol_flx
  real(r8) :: H2CO3_activity,H0PO4_3e_activity,H1PO4_2e_activity,H2PO4_1e_activity,H3PO4_activity,FeHPO4_p_activity,FeH2PO4_2p_activity
  real(r8) :: CaPO4_1e_activity,CaHPO4_activity,CaH4P2O8_1p_activity,MgHPO4_activity,H0PO4_3e_band_activity
  real(r8) :: H1PO4_2e_band_activity,H2PO4_1e_band_activity,H3PO4_band_activity
  real(r8) :: A2,H_1p_activity,Al_3p_activity,AlOH_2p_activity,AlO2H2_1p_activity,AlO3H3_activity,AlO4H4_1e_activity
  real(r8) :: Fe_3p_activity,FeOH_2p_activity,FeO2H2_p_activity,FeO3H3_activity,FeO4H4_1e_activity,Ca_2p_activity,CO3_2e_activity,HCO3_e_activity
  real(r8) :: H2PO4_1e_CaHPO4_dissolB_flx,H2PO4_1e_apatite_dissolB_flx
  real(r8) :: H2PO4_1e_CaH4P2O8_dissolB_flx,H2PO4_1e_FePO4_dissolB_flx,RSO4,RX1P,RX2P
  real(r8) :: RXFE,RXCA,RXMG,RXNA,RXKA,RXHC,RXALO2,RXFEO2,RCO2Q
  real(r8) :: RXH1,RXH2,RXNB,Ca_2p_aque_mole_conc,OH_1e_aque_mole_conc,H1PO4_to_XHPO4_ROH_flx
  real(r8) :: H2PO4_1e_to_XH2PO4_ROH2_Bflx,H2PO4_1e_FePO4_dissol_flx
  real(r8) :: RALO1,RALO2,RALO3,RALO4,RNHB,RXH1B,RXN4
  real(r8) :: FeHPO4_1p_band_activity,FeH2PO4_2p_band_activity,CaPO4_1e_band_activity
  real(r8) :: CaHPO4_band_activity,CaH4P2O8_1p_band_activity,MgHPO4_band_activity,NH4_1p_activity,NH4_1p_band_activity
  real(r8) :: NH3_activity,NH3_band_activity,Mg_2p_activity,Na_1p_activity
  real(r8) :: K_1p_activity,AlSO4_1p_activity,FeSO4_1p_activity,CaO2H2_activity,CaCO3_activity
  real(r8) :: RH2B,H2PO4_1e_to_XH2PO4_ROH2_flx,H2PO4_to_XH2PO4_ROH_Bflx,H2PO4_to_XH2PO4_ROH_flx,AKAS1,AALX,AFEX,ACAX,AMGX
  real(r8) :: CaSO4_activity,CaHCO3_1p_activity,MgOH_1p_activity,AMGC1,MgHCO3_1p_activity,MgSO4_activity,NaCO3_1e_activity,NaSO4_1e_activity
  real(r8) :: VLWatMicPBK

  real(r8), pointer :: VLWatMicPM
  real(r8), pointer :: VLWatMicPNH
  real(r8), pointer :: VLWatMicPNB
  real(r8), pointer :: VLWatMicPPB
  real(r8), pointer :: VLWatMicPPO
  real(r8), pointer :: VLPOB
  real(r8), pointer :: VLPO4
  real(r8), pointer :: VLNHB
  real(r8), pointer :: VLNH4
  real(r8), pointer :: VLNOB
  real(r8), pointer :: VLNO3
  real(r8), pointer :: SoilMicPMassLayer
  real(r8), pointer :: SoilMicPMassLayerX
  real(r8), pointer :: VLWatMicPNO
  real(r8), pointer :: VLWatMicPNZ
  real(r8), pointer :: GKC4     !Ca-NH4 Gapon selectivity coefficient, [-]
  real(r8), pointer :: GKCA     !Ca-Al Gapon selectivity coefficient, [-]
  real(r8), pointer :: GKCH     !Ca-H Gapon selectivity coefficient, [-]
  real(r8), pointer :: GKCM     !Ca-Mg Gapon selectivity coefficient, [-]
  real(r8), pointer :: GKCK     !Ca-K Gapon selectivity coefficient, [-]
  real(r8), pointer :: GKCN     !Ca-Na Gapon selectivity coefficient, [-]
  real(r8), pointer :: XCEC     !cation exchange capacity, [mol d-2]
  real(r8), pointer :: XAEC     !anion exchange capacity, [mol d-2]
  real(r8), pointer :: ORGC     !total soil organic C [g d-2]
  real(r8), pointer :: PH       !pH value of the system
  real(r8), pointer :: RProd_Hp    !total H+ production, [flux]

  real(r8), pointer :: CO2S     !aqueous CO2  micropore	[g d-2]
  real(r8), pointer :: H1PO4_2e_aque_mole_conc    !soil aqueous HPO4 content micropore non-band, [mol m-3]
  real(r8), pointer :: H1PO4_2e_band_mole_conc    !soil aqueous HPO4 content micropore band, [mol m-3]
  real(r8), pointer :: H2PO4_1e_aque_mole_conc    !soil aqueous H2PO4 content micropore non-band, [mol m-3]
  real(r8), pointer :: H2PO4_1e_band_mole_conc    !soil aqueous H2PO4 content micropore  band, [mol m-3]
  real(r8), pointer :: NH3_aque_mole_conc     !soil NH3 concentration in non-band soil, [mol m-3]
  real(r8), pointer :: NH3_aqu_band_mole_conc     !soil NH3 concentration in band soil, [mol m-3]
  real(r8), pointer :: NH4_1p_aque_mole_conc     !soil NH4 concentration in non-band soil, [mol m-3]
  real(r8), pointer :: NH4_1p_band_mole_conc     !soil NH4 concentration in band soil, [mol m-3]
  real(r8), pointer :: ZNO3S    !NO3 mass non-band micropore, [g d-2]
  real(r8), pointer :: ZNO3B    !NO3 mass band micropore, [g d-2]
  real(r8), pointer :: ZHY      !soil aqueous H content micropore, [mol d-2]
  real(r8), pointer :: ZOH      !soil aqueous OH content micropore, [mol d-2]
  real(r8), pointer :: ZAL      !soil aqueous Al content micropore, [mol d-2]
  real(r8), pointer :: ZFE      !soil aqueous Fe content micropore, [mol d-2]
  real(r8), pointer :: ZCA      !soil aqueous Ca content micropore, [mol d-2]
  real(r8), pointer :: ZMG      !soil aqueous Mg content micropore, [mol d-2]
  real(r8), pointer :: ZNA      !soil aqueous Na content micropore, [mol d-2]
  real(r8), pointer :: ZKA      !soil aqueous K content micropore, [mol d-2]
  real(r8), pointer :: ZSO4     !soil aqueous SO4 content micropore, [mol d-2]
  real(r8), pointer :: ZCL      !soil aqueous Cl content micropore, [mol d-2]
  real(r8), pointer :: ZCO3     !soil aqueous CO3 content micropore, [mol d-2]
  real(r8), pointer :: ZHCO3    !soil aqueous HCO3 content micropore, [mol d-2]
  real(r8), pointer :: ZALOH1   !soil aqueous AlOH content micropore, [mol d-2]
  real(r8), pointer :: ZALOH2   !soil aqueous AlOH2 content micropore, [mol d-2]
  real(r8), pointer :: ZALOH3   !soil aqueous AlOH3 content micropore, [mol d-2]
  real(r8), pointer :: ZALOH4   !soil aqueous AlOH4 content micropore, [mol d-2]
  real(r8), pointer :: ZALS     !soil aqueous AlSO4 content micropore, [mol d-2]
  real(r8), pointer :: ZFEOH1   !soil aqueous FeOH content micropore, [mol d-2]
  real(r8), pointer :: ZFEOH2   !soil aqueous FeOH2 content micropore, [mol d-2]
  real(r8), pointer :: ZFEOH3   !soil aqueous FeOH3 content micropore, [mol d-2]
  real(r8), pointer :: ZFEOH4   !soil aqueous FeOH4 content micropore, [mol d-2]
  real(r8), pointer :: ZFES     !soil aqueous FeSO4 content micropore, [mol d-2]
  real(r8), pointer :: ZCAO     !soil aqueous CaOH2 content micropore, [mol d-2]
  real(r8), pointer :: ZCAC     !soil aqueous CACO3 content micropore, [mol d-2]
  real(r8), pointer :: ZCAH     !soil aqueous CaHCO3 content micropore, [mol d-2]
  real(r8), pointer :: ZCAS     !soil aqueous CaSO4 content micropore, [mol d-2]
  real(r8), pointer :: ZMGO     !soil aqueous MgOH content micropore, [mol d-2]
  real(r8), pointer :: ZMGC     !soil aqueous MgCO3 content micropore, [mol d-2]
  real(r8), pointer :: ZMGH     !soil aqueous MgHCO3 content micropore, [mol d-2]
  real(r8), pointer :: ZMGS     !soil aqueous MgSO4 content micropore, [mol d-2]
  real(r8), pointer :: ZNAC     !soil aqueous NaCO3 content micropore, [mol d-2]
  real(r8), pointer :: ZNAS     !soil aqueous NaSO4 content micropore, [mol d-2]
  real(r8), pointer :: ZKAS     !soil aqueous KSO4 content micropore, [mol d-2]
  real(r8), pointer :: H0PO4    !soil aqueous PO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: H3PO4    !soil aqueous H3PO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: ZFE1P    !soil aqueous FeHPO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: ZFE2P    !soil aqueous FeH2PO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: ZCA0P    !soil aqueous CaPO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: ZCA1P    !soil aqueous CaHPO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: ZCA2P    !soil aqueous CaH4P2O8 content micropore non-band, [mol d-2]
  real(r8), pointer :: ZMG1P    !soil aqueous MgHPO4 content micropore non-band, [mol d-2]
  real(r8), pointer :: H0POB    !soil aqueous PO4 content micropore band, [mol d-2]
  real(r8), pointer :: H3POB    !soil aqueous H3PO4 content micropore band, [mol d-2]
  real(r8), pointer :: ZFE1PB   !soil aqueous FeHPO4 content micropore band, [mol d-2]
  real(r8), pointer :: ZFE2PB   !soil aqueous FeH2PO4 content micropore band, [mol d-2
  real(r8), pointer :: ZCA0PB   !soil aqueous CaPO4 content micropore band, [mol d-2]
  real(r8), pointer :: ZCA1PB   !soil aqueous CaHPO4 content micropore band, [mol d-2]
  real(r8), pointer :: ZCA2PB   !soil aqueous CaH4P2O8 content micropore band, [mol d-2]
  real(r8), pointer :: ZMG1PB   !soil aqueous MgHPO4 content micropore band, [mol d-2]
  real(r8), pointer :: Precp_AlPO4_mole_conc   !precipitated AlPO4 non-band, [mol m-3]
  real(r8), pointer :: PrecpB_AlPO4_mole_conc   !precipitated AlPO4 band soil, [mol m-3]
  real(r8), pointer :: Precp_CaHPO4_mole_conc   !precipitated CaHPO4 non-band soil, [mol m-3]
  real(r8), pointer :: PrecpB_CaHPO4_mole_conc   !precipitated CaHPO4 band soil, [mol m-3]
  real(r8), pointer :: Precp_Ca5P3O12O3H3_mole_conc   !precipitated Ca5(PO4)3OH hydroxyapatite non-band soil, [mol m-3]
  real(r8), pointer :: PrecpB_Ca5P3O12O3H3_mole_conc   !precipitated Ca5(PO4)3OH hydroxyapatite band soil, [mol m-3]
  real(r8), pointer :: Precp_CaH4P2O8_mole_conc   !precipitated Ca(H2PO4)2 non-band soil, [mol m-3]
  real(r8), pointer :: PrecpB_CaH4P2O8_mole_conc   !precipitated CaH4P2O8 band soil, [mol m-3]
  real(r8), pointer :: Precp_FePO4_mole_conc   !precipitated FePO4 non-band soil, [mol m-3]
  real(r8), pointer :: PrecpB_FePO4_mole_conc   !precipitated FePO4 band soil, [mol m-3]
  real(r8), pointer :: PALOH    !precipitated Al(OH)3, [mol d-2]
  real(r8), pointer :: PFEOH    !precipitated Fe(OH)3, [mol d-2]
  real(r8), pointer :: PCACO    !precipitated CaCO3, [mol d-2]
  real(r8), pointer :: PCASO    !precipitated CaSO4, [mol d-2]
  real(r8), pointer :: XHY      !exchangeable H(+) , [mol d-2]
  real(r8), pointer :: XAL      !exchangeable Al, [mol d-2]
  real(r8), pointer :: XFE      !exchangeable Fe, [mol d-2]
  real(r8), pointer :: XCA      !exchangeable Ca, [mol d-2]
  real(r8), pointer :: XMG      !exchangeable Mg, [mol d-2]
  real(r8), pointer :: XNA      !exchangeable Na, [mol d-2]
  real(r8), pointer :: XKA      !exchangeable K, [mol d-2]
  real(r8), pointer :: XHC      !exchangeable COOH , [mol d-2]
  real(r8), pointer :: XROH1_mole_conc    !exchangeable OH  non-band, [mol d-2]
  real(r8), pointer :: XALO2    !exchangeable AlOH2 , [mol d-2]
  real(r8), pointer :: XFEO2    !exchangeable Fe(OH)2, [mol d-2]
  real(r8), pointer :: XNH4_mole_conc     !exchangeable NH4 non-band soil, [mol d-2]
  real(r8), pointer :: XNH4_band_mole_conc     !exchangeable NH4 band soil, [mol d-2]
  real(r8), pointer :: XROH1_band_mole_conc    !exchangeable OH- band, [mol d-2]
  real(r8), pointer :: XOH_mole_conc    !exchangeable OH- non-band, [mol d-2]
  real(r8), pointer :: XHPO4_band_mole_conc    !exchangeable HPO4 concentration band-soil, [mol m-3]
  real(r8), pointer :: XH2PO4_band_mole_conc    !exchangeable H2PO4 concentration band-soil, [mol m-3]
  real(r8), pointer :: XROH_band_mole_conc    !exchangeable OH band-soil, [mol m-3]
  real(r8), pointer :: XHPO4_mole_conc    !exchangeable HPO4  non-band, [mol m-3]
  real(r8), pointer :: XROH2_band_mole_conc    !exchangeable OH2 band-soil, [mol m-3]
  real(r8), pointer :: XROH2_mole_conc    !exchangeable OH2  non-band soil, [mol m-3]
  real(r8), pointer :: XH2PO4_mole_conc    !exchangeable H2PO4  non-band soil, [mol m-3]

! fluxes
  real(r8), pointer :: TRChem_CaCO3_precip_soil   !total precipitated CaCO3 transformation, [mol d-2 h-1]
  real(r8), pointer :: TRChem_NaCO3_soil
  real(r8), pointer :: TRChem_MgCO3_soil
  real(r8), pointer :: TRChem_CaCO3_soil
  real(r8), pointer :: TRChem_MgHCO3_soil
  real(r8), pointer :: TRChem_CaHCO3_soil
  real(r8), pointer :: TRChem_HCO3_soil
  real(r8), pointer :: TRChem_HCO3_sorbed_soil
  real(r8), pointer :: TRChem_CaH4P2O8_soil
  real(r8), pointer :: TRChem_FeH2PO4_soil
  real(r8), pointer :: TRChem_CaH4P2O8_band_soil
  real(r8), pointer :: TRChem_FeH2PO4_band_soil
  real(r8), pointer :: TRChem_MgHPO4_soil
  real(r8), pointer :: TRChem_CaHPO4_soil
  real(r8), pointer :: TRChem_FeHPO4_soil
  real(r8), pointer :: TRChem_MgHPO4_band_soil
  real(r8), pointer :: TRChem_CaHPO4_band_soil
  real(r8), pointer :: TRChem_FeHPO4_band_soil
  real(r8), pointer :: TRChem_CaPO4_band_soil
  real(r8), pointer :: TRChem_PO4_band_soil
  real(r8), pointer :: TRChem_CaPO4_soil
  real(r8), pointer :: TRChem_PO4_soil
  real(r8), pointer :: TRChem_FePO4_precip_soil
  real(r8), pointer :: TRChem_AlPO4_precip_soil
  real(r8), pointer :: TRChem_FePO4_precip_band_soil
  real(r8), pointer :: TRChem_AlPO4_precip_band_soil
  real(r8), pointer :: TRChem_CaH4P2O8_precip_soil
  real(r8), pointer :: TRChem_CaHPO4_precip_soil
  real(r8), pointer :: TRChem_CaH4P2O8_precip_band_soil
  real(r8), pointer :: TRChem_CaHPO4_precip_band_soil
  real(r8), pointer :: TRChem_apatite_precip_band_soil
  real(r8), pointer :: TRChem_apatite_precip_soil
  real(r8), pointer :: TRChem_FeOH_soil
  real(r8), pointer :: TRChem_AlOH_soil
  real(r8), pointer :: TRChem_FeO2H2_soil
  real(r8), pointer :: TRChem_AlO2H2_soil
  real(r8), pointer :: TRChem_FeO3H3_soil
  real(r8), pointer :: TRChem_AlO3H3_soil
  real(r8), pointer :: TRChem_FeO4H4_soil
  real(r8), pointer :: TRChem_AlO4H4_soil
  real(r8), pointer :: TRChem_MgOH_soil
  real(r8), pointer :: TRChem_CaOH_soil
  real(r8), pointer :: TRChem_FeOH3_precip_soil
  real(r8), pointer :: TRChem_AlOH3_precip_soil
  real(r8), pointer :: TRChem_FeO2H2_sorbed_soil
  real(r8), pointer :: TRChem_Al_3p_soil
  real(r8), pointer :: TRChem_AlSO4_soil
  real(r8), pointer :: TRChem_RHPO4_sorbed_band_soil
  real(r8), pointer :: TRChem_RH2PO4_sorbed_band_soil
  real(r8), pointer :: TRChem_RO_sorbed_band_soil
  real(r8), pointer :: TRChem_ROH_sorbed_band_soil
  real(r8), pointer :: TRChem_ROH2_sorbed_band_soil
  real(r8), pointer :: TRChem_Ca_2p_soil
  real(r8), pointer :: TRChem_CaSO4_soil
  real(r8), pointer :: TRChem_CaSO4_precip_soil
  real(r8), pointer :: TRChem_CO2_gchem_soil
  real(r8), pointer :: TRChem_CO3_2e_soil
  real(r8), pointer :: TRChem_Fe_3p_soil
  real(r8), pointer :: TRChem_FeSO4_soil
  real(r8), pointer :: TRChem_H1PO4_band_soil
  real(r8), pointer :: TRChem_H2PO4_band_soil
  real(r8), pointer :: TRChem_H3PO4_band_soil
  real(r8), pointer :: TRChem_H1PO4_soil
  real(r8), pointer :: TRChem_H2PO4_soil
  real(r8), pointer :: TRChem_H3PO4_sorbed_soil
  real(r8), pointer :: TRChem_H_p_soil
  real(r8), pointer :: TRChem_K_1p_soil
  real(r8), pointer :: TRChem_KSO4_soil
  real(r8), pointer :: TRChem_Mg_2p_soil
  real(r8), pointer :: TRChem_MgSO4_soil
  real(r8), pointer :: TRChem_NH3_band_soil
  real(r8), pointer :: TRChem_NH3_soil_mole
  real(r8), pointer :: TRChem_NH4_band_soil_mole
  real(r8), pointer :: TRChem_NH4_soil_mole
  real(r8), pointer :: TRChem_Na_p_soil
  real(r8), pointer :: TRChem_NaSO4_soil
  real(r8), pointer :: TRChem_OH_1e_soil
  real(r8), pointer :: TRChem_SO4_2e_soil
  real(r8), pointer :: TRChem_RHPO4_sorbed_soil
  real(r8), pointer :: TRChem_RH2PO4_sorbed_soil
  real(r8), pointer :: TRChem_Al_sorbed_soil
  real(r8), pointer :: TRChem_AlO2H2_sorbed_soil   !total adsorbed AlOH2 transformation, [mol d-2 h-1]
  real(r8), pointer :: TRChem_Ca_sorbed_soil
  real(r8), pointer :: TRChem_Fe_sorbed_soil
  real(r8), pointer :: TRChem_RO_sorbed_soil
  real(r8), pointer :: TRChem_ROH_sorbed_soil
  real(r8), pointer :: TRChem_ROH2_sorbed_soil
  real(r8), pointer :: TRChem_H_p_sorbed_soil
  real(r8), pointer :: TRChem_K_sorbed_soil
  real(r8), pointer :: TRChem_Mg_sorbed_soil
  real(r8), pointer :: TRChem_NH4_sorbed_soil
  real(r8), pointer :: TRChem_Na_sorbed_soil
  real(r8), pointer :: TRChem_NH4_sorbed_band_soil
  real(r8), pointer :: TRAquaChem_CO2_soil
  real(r8), pointer :: TBION_soil
  real(r8), pointer :: TRH2O_soil
  public :: SaltChemEquilibria
  contains

!--------------------------------------------------------------------------
  subroutine SaltChemEquilibria(chemvar,solflx)

  implicit none
  type(solute_flx_type),target, intent(inout) :: solflx
  type(chem_var_type),target, intent(inout) :: chemvar

!     begin_execution

  SoilMicPMassLayerX            => chemvar%SoilMicPMassLayerX
  VLWatMicPNZ                   => chemvar%VLWatMicPNZ
  VLWatMicPNO                   => chemvar%VLWatMicPNO
  VLWatMicPNB                   => chemvar%VLWatMicPNB
  VLWatMicPNH                   => chemvar%VLWatMicPNH
  VLWatMicPPB                   => chemvar%VLWatMicPPB
  VLWatMicPPO                   => chemvar%VLWatMicPPO
  XROH1_mole_conc               => chemvar%XROH1_mole_conc
  XROH2_mole_conc               => chemvar%XROH2_mole_conc
  XNH4_mole_conc                => chemvar%XNH4_mole_conc
  XNH4_band_mole_conc           => chemvar%XNH4_band_mole_conc
  H1PO4_2e_aque_mole_conc       => chemvar%H1PO4_2e_aque_mole_conc
  H1PO4_2e_band_mole_conc       => chemvar%H1PO4_2e_band_mole_conc
  H2PO4_1e_aque_mole_conc       => chemvar%H2PO4_1e_aque_mole_conc
  H2PO4_1e_band_mole_conc       => chemvar%H2PO4_1e_band_mole_conc
  XHPO4_band_mole_conc          => chemvar%XHPO4_band_mole_conc
  XH2PO4_band_mole_conc         => chemvar%XH2PO4_band_mole_conc
  XROH_band_mole_conc           => chemvar%XROH_band_mole_conc
  XHPO4_mole_conc               => chemvar%XHPO4_mole_conc
  XROH2_band_mole_conc          => chemvar%XROH2_band_mole_conc
  XH2PO4_mole_conc              => chemvar%XH2PO4_mole_conc
  NH3_aque_mole_conc            => chemvar%NH3_aque_mole_conc
  NH3_aqu_band_mole_conc        => chemvar%NH3_aqu_band_mole_conc
  NH4_1p_aque_mole_conc         => chemvar%NH4_1p_aque_mole_conc
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
  PrecpB_FePO4_mole_conc        => chemvar%PrecpB_FePO4_mole_conc
  ZNO3S                         => chemvar%ZNO3S
  ZNO3B                         => chemvar%ZNO3B
  VLWatMicPM                    => chemvar%VLWatMicPM
  RProd_Hp                      => chemvar%RProd_Hp
  ZHY                           => chemvar%ZHY
  XCEC                          => chemvar%XCEC
  ZOH                           => chemvar%ZOH
  ZAL                           => chemvar%ZAL
  ZFE                           => chemvar%ZFE
  ZCA                           => chemvar%ZCA
  ZMG                           => chemvar%ZMG
  ZNA                           => chemvar%ZNA
  ZKA                           => chemvar%ZKA
  ZSO4                          => chemvar%ZSO4
  ZCL                           => chemvar%ZCL
  ZCO3                          => chemvar%ZCO3
  ZHCO3                         => chemvar%ZHCO3
  CO2S                          => chemvar%CO2S
  ZALOH1                        => chemvar%ZALOH1
  ZALOH2                        => chemvar%ZALOH2
  ZALOH3                        => chemvar%ZALOH3
  ZALOH4                        => chemvar%ZALOH4
  ZALS                          => chemvar%ZALS
  ZFEOH1                        => chemvar%ZFEOH1
  ZFEOH2                        => chemvar%ZFEOH2
  ZFEOH3                        => chemvar%ZFEOH3
  ZFEOH4                        => chemvar%ZFEOH4
  ZFES                          => chemvar%ZFES
  ZCAO                          => chemvar%ZCAO
  ZCAC                          => chemvar%ZCAC
  ZCAH                          => chemvar%ZCAH
  ZCAS                          => chemvar%ZCAS
  ZMGO                          => chemvar%ZMGO
  ZMGC                          => chemvar%ZMGC
  ZMGH                          => chemvar%ZMGH
  ZMGS                          => chemvar%ZMGS
  ZNAC                          => chemvar%ZNAC
  ZNAS                          => chemvar%ZNAS
  ZKAS                          => chemvar%ZKAS
  H0PO4                         => chemvar%H0PO4
  H3PO4                         => chemvar%H3PO4
  ZFE1P                         => chemvar%ZFE1P
  ZFE2P                         => chemvar%ZFE2P
  ZCA0P                         => chemvar%ZCA0P
  ZCA1P                         => chemvar%ZCA1P
  ZCA2P                         => chemvar%ZCA2P
  ZMG1P                         => chemvar%ZMG1P
  H0POB                         => chemvar%H0POB
  H3POB                         => chemvar%H3POB
  ZFE1PB                        => chemvar%ZFE1PB
  ZFE2PB                        => chemvar%ZFE2PB
  ZCA0PB                        => chemvar%ZCA0PB
  ZCA1PB                        => chemvar%ZCA1PB
  ZCA2PB                        => chemvar%ZCA2PB
  ZMG1PB                        => chemvar%ZMG1PB
  XHY                           => chemvar%XHY
  XAL                           => chemvar%XAL
  XFE                           => chemvar%XFE
  XCA                           => chemvar%XCA
  XMG                           => chemvar%XMG
  XNA                           => chemvar%XNA
  XKA                           => chemvar%XKA
  XHC                           => chemvar%XHC
  XALO2                         => chemvar%XALO2
  XFEO2                         => chemvar%XFEO2
  ORGC                          => chemvar%ORGC
  PALOH                         => chemvar%PALOH
  PFEOH                         => chemvar%PFEOH
  PCACO                         => chemvar%PCACO
  PCASO                         => chemvar%PCASO
  PH                            => chemvar%PH
  VLPOB                         => chemvar%VLPOB
  VLPO4                         => chemvar%VLPO4
  VLNHB                         => chemvar%VLNHB
  VLNH4                         => chemvar%VLNH4
  VLNOB                         => chemvar%VLNOB
  VLNO3                         => chemvar%VLNO3
  SoilMicPMassLayer             => chemvar%SoilMicPMassLayer
  XAEC                          => chemvar%XAEC
  GKC4                          => chemvar%GKC4
  GKCA                          => chemvar%GKCA
  GKCH                          => chemvar%GKCH
  GKCM                          => chemvar%GKCM
  GKCK                          => chemvar%GKCK
  GKCN                          => chemvar%GKCN
  XROH1_band_mole_conc          => chemvar%XROH1_band_mole_conc
  XOH_mole_conc                 => chemvar%XOH_mole_conc

  TRChem_CaCO3_precip_soil         => solflx%TRChem_CaCO3_precip_soil
  TRChem_NaCO3_soil                => solflx%TRChem_NaCO3_soil
  TRChem_MgCO3_soil                => solflx%TRChem_MgCO3_soil
  TRChem_CaCO3_soil                => solflx%TRChem_CaCO3_soil
  TRChem_MgHCO3_soil               => solflx%TRChem_MgHCO3_soil
  TRChem_CaHCO3_soil               => solflx%TRChem_CaHCO3_soil
  TRChem_HCO3_soil                 => solflx%TRChem_HCO3_soil
  TRChem_HCO3_sorbed_soil          => solflx%TRChem_HCO3_sorbed_soil
  TRChem_CaH4P2O8_soil             => solflx%TRChem_CaH4P2O8_soil
  TRChem_FeH2PO4_soil              => solflx%TRChem_FeH2PO4_soil
  TRChem_CaH4P2O8_band_soil        => solflx%TRChem_CaH4P2O8_band_soil
  TRChem_FeH2PO4_band_soil         => solflx%TRChem_FeH2PO4_band_soil
  TRChem_MgHPO4_soil               => solflx%TRChem_MgHPO4_soil
  TRChem_CaHPO4_soil               => solflx%TRChem_CaHPO4_soil
  TRChem_FeHPO4_soil               => solflx%TRChem_FeHPO4_soil
  TRChem_MgHPO4_band_soil          => solflx%TRChem_MgHPO4_band_soil
  TRChem_CaHPO4_band_soil          => solflx%TRChem_CaHPO4_band_soil
  TRChem_FeHPO4_band_soil          => solflx%TRChem_FeHPO4_band_soil
  TRChem_CaPO4_band_soil           => solflx%TRChem_CaPO4_band_soil
  TRChem_PO4_band_soil             => solflx%TRChem_PO4_band_soil
  TRChem_CaPO4_soil                => solflx%TRChem_CaPO4_soil
  TRChem_PO4_soil                  => solflx%TRChem_PO4_soil
  TRChem_FePO4_precip_soil         => solflx%TRChem_FePO4_precip_soil
  TRChem_AlPO4_precip_soil         => solflx%TRChem_AlPO4_precip_soil
  TRChem_FePO4_precip_band_soil    => solflx%TRChem_FePO4_precip_band_soil
  TRChem_AlPO4_precip_band_soil    => solflx%TRChem_AlPO4_precip_band_soil
  TRChem_CaH4P2O8_precip_soil      => solflx%TRChem_CaH4P2O8_precip_soil
  TRChem_CaHPO4_precip_soil        => solflx%TRChem_CaHPO4_precip_soil
  TRChem_CaH4P2O8_precip_band_soil => solflx%TRChem_CaH4P2O8_precip_band_soil
  TRChem_CaHPO4_precip_band_soil   => solflx%TRChem_CaHPO4_precip_band_soil
  TRChem_apatite_precip_band_soil  => solflx%TRChem_apatite_precip_band_soil
  TRChem_apatite_precip_soil       => solflx%TRChem_apatite_precip_soil
  TRChem_FeOH_soil                 => solflx%TRChem_FeOH_soil
  TRChem_AlOH_soil                 => solflx%TRChem_AlOH_soil
  TRChem_FeO2H2_soil               => solflx%TRChem_FeO2H2_soil
  TRChem_AlO2H2_soil               => solflx%TRChem_AlO2H2_soil
  TRChem_FeO3H3_soil               => solflx%TRChem_FeO3H3_soil
  TRChem_AlO3H3_soil               => solflx%TRChem_AlO3H3_soil
  TRChem_FeO4H4_soil               => solflx%TRChem_FeO4H4_soil
  TRChem_AlO4H4_soil               => solflx%TRChem_AlO4H4_soil
  TRChem_MgOH_soil                 => solflx%TRChem_MgOH_soil
  TRChem_CaOH_soil                 => solflx%TRChem_CaOH_soil
  TRChem_FeOH3_precip_soil         => solflx%TRChem_FeOH3_precip_soil
  TRChem_AlOH3_precip_soil         => solflx%TRChem_AlOH3_precip_soil
  TRChem_FeO2H2_sorbed_soil        => solflx%TRChem_FeO2H2_sorbed_soil
  TRChem_Al_3p_soil                => solflx%TRChem_Al_3p_soil
  TRChem_AlSO4_soil                => solflx%TRChem_AlSO4_soil
  TRChem_RHPO4_sorbed_band_soil    => solflx%TRChem_RHPO4_sorbed_band_soil
  TRChem_RH2PO4_sorbed_band_soil   => solflx%TRChem_RH2PO4_sorbed_band_soil
  TRChem_RO_sorbed_band_soil       => solflx%TRChem_RO_sorbed_band_soil
  TRChem_ROH_sorbed_band_soil      => solflx%TRChem_ROH_sorbed_band_soil
  TRChem_ROH2_sorbed_band_soil     => solflx%TRChem_ROH2_sorbed_band_soil
  TRChem_Ca_2p_soil                => solflx%TRChem_Ca_2p_soil
  TRChem_CaSO4_soil                => solflx%TRChem_CaSO4_soil
  TRChem_CaSO4_precip_soil         => solflx%TRChem_CaSO4_precip_soil
  TRChem_CO2_gchem_soil            => solflx%TRChem_CO2_gchem_soil
  TRChem_CO3_2e_soil               => solflx%TRChem_CO3_2e_soil
  TRChem_Fe_3p_soil                => solflx%TRChem_Fe_3p_soil
  TRChem_FeSO4_soil                => solflx%TRChem_FeSO4_soil
  TRChem_H1PO4_band_soil           => solflx%TRChem_H1PO4_band_soil
  TRChem_H2PO4_band_soil           => solflx%TRChem_H2PO4_band_soil
  TRChem_H3PO4_band_soil           => solflx%TRChem_H3PO4_band_soil
  TRChem_H1PO4_soil                => solflx%TRChem_H1PO4_soil
  TRChem_H2PO4_soil                => solflx%TRChem_H2PO4_soil
  TRChem_H3PO4_sorbed_soil         => solflx%TRChem_H3PO4_sorbed_soil
  TRChem_H_p_soil                  => solflx%TRChem_H_p_soil
  TRChem_K_1p_soil                 => solflx%TRChem_K_1p_soil
  TRChem_KSO4_soil                 => solflx%TRChem_KSO4_soil
  TRChem_Mg_2p_soil                => solflx%TRChem_Mg_2p_soil
  TRChem_MgSO4_soil                => solflx%TRChem_MgSO4_soil
  TRChem_NH3_band_soil             => solflx%TRChem_NH3_band_soil
  TRChem_NH3_soil_mole             => solflx%TRChem_NH3_soil_mole
  TRChem_NH4_band_soil_mole        => solflx%TRChem_NH4_band_soil_mole
  TRChem_NH4_soil_mole             => solflx%TRChem_NH4_soil_mole
  TRChem_Na_p_soil                 => solflx%TRChem_Na_p_soil
  TRChem_NaSO4_soil                => solflx%TRChem_NaSO4_soil
  TRChem_OH_1e_soil                => solflx%TRChem_OH_1e_soil
  TRChem_SO4_2e_soil               => solflx%TRChem_SO4_2e_soil
  TRChem_RHPO4_sorbed_soil         => solflx%TRChem_RHPO4_sorbed_soil
  TRChem_RH2PO4_sorbed_soil        => solflx%TRChem_RH2PO4_sorbed_soil
  TRChem_Al_sorbed_soil            => solflx%TRChem_Al_sorbed_soil
  TRChem_AlO2H2_sorbed_soil        => solflx%TRChem_AlO2H2_sorbed_soil
  TRChem_Ca_sorbed_soil            => solflx%TRChem_Ca_sorbed_soil
  TRChem_Fe_sorbed_soil            => solflx%TRChem_Fe_sorbed_soil
  TRChem_RO_sorbed_soil            => solflx%TRChem_RO_sorbed_soil
  TRChem_ROH_sorbed_soil           => solflx%TRChem_ROH_sorbed_soil
  TRChem_ROH2_sorbed_soil          => solflx%TRChem_ROH2_sorbed_soil
  TRChem_H_p_sorbed_soil           => solflx%TRChem_H_p_sorbed_soil
  TRChem_K_sorbed_soil             => solflx%TRChem_K_sorbed_soil
  TRChem_Mg_sorbed_soil            => solflx%TRChem_Mg_sorbed_soil
  TRChem_NH4_sorbed_soil           => solflx%TRChem_NH4_sorbed_soil
  TRChem_Na_sorbed_soil            => solflx%TRChem_Na_sorbed_soil
  TRChem_NH4_sorbed_band_soil      => solflx%TRChem_NH4_sorbed_band_soil
  TRAquaChem_CO2_soil              => solflx%TRAquaChem_CO2_soil
  TBION_soil                       => solflx%TBION_soil
  TRH2O_soil                       => solflx%TRH2O_soil
!
  call PrepIonConcentrations
!
!     CONVERGENCE TOWARDS SOLUTE EQILIBRIA
!
  call SolveChemEquilibria
!
  call SummarizeIonFluxes

  end subroutine SaltChemEquilibria

!------------------------------------------------------------------------
  subroutine PrepIonConcentrations
  implicit none
  real(r8) :: VLWatMicPPX
!     begin_execution
!     SOLUBLE NO3 CONCENTRATIONS
!     IN NON-BAND AND BAND SOIL ZONES
!
!     VLWatMicPNO,VLWatMicPNZ=soil water volume in NO3 non-band,band
!     ZNO3S,ZNO3B=NO3 mass in non-band,band
!     NO3_1e_aque_mole_conc,NO3_1e_band_conc=NO3 concentrations in non-band,band
!
  IF(VLWatMicPNO.GT.ZEROS2)THEN
    NO3_1e_aque_mole_conc=AZMAX1(ZNO3S/(natomw*VLWatMicPNO))
  ELSE
    NO3_1e_aque_mole_conc=0._r8
  ENDIF

  IF(VLWatMicPNZ.GT.ZEROS2)THEN
    NO3_1e_band_conc=AZMAX1(ZNO3B/(natomw*VLWatMicPNZ))
  ELSE
    NO3_1e_band_conc=0._r8
  ENDIF
!
!     H CONCENTRATION
!
!     RProd_Hp=total H+ production from nitro.f
!
  H_1p_aque_mole_conc=AZMAX1(ZHY+RProd_Hp)/VLWatMicPM
!
!     SOLUTE ION AND ION PAIR CONCENTRATIONS
!
!     CEC_conc,XCEC=cation exchange concentration,capacity
!     SoilMicPMassLayerX=soil mass
!     VLWatMicPM=soil water volume
!
  IF(SoilMicPMassLayerX.GT.ZEROS)THEN
    CEC_conc=AMAX1(ZERO,XCEC/SoilMicPMassLayerX)
  ELSE
    CEC_conc=ZERO
  ENDIF

  IF(VLWatMicPM.GT.ZEROS2)THEN
    OH_1e_aque_mole_conc     = AZMAX1(ZOH/VLWatMicPM)
    Al_3p_aque_mole_conc     = AZMAX1(ZAL/VLWatMicPM)
    Fe_3p_aque_mole_conc     = AZMAX1(ZFE/VLWatMicPM)
    Ca_2p_aque_mole_conc     = AZMAX1(ZCA/VLWatMicPM)
    Mg_2p_aque_mole_conc     = AZMAX1(ZMG/VLWatMicPM)
    Na_1p_aque_mole_conc     = AZMAX1(ZNA/VLWatMicPM)
    K_1p_aque_mole_conc      = AZMAX1(ZKA/VLWatMicPM)
    SO4_2e_aque_mole_conc    = AZMAX1(ZSO4/VLWatMicPM)
    Cl_e_conc                = AZMAX1(ZCL/VLWatMicPM)
    CO3_2e_aque_mole_conc    = AZMAX1(ZCO3/VLWatMicPM)
    HCO3_e_mole_conc         = AZMAX1(ZHCO3/VLWatMicPM)
    H2CO3_aque_mole_conc     = AZMAX1(CO2S/(catomw*VLWatMicPM))
    AlOH_2p_aque_mole_conc   = AZMAX1(ZALOH1/VLWatMicPM)
    AlO2H2_1p_aque_mole_conc = AZMAX1(ZALOH2/VLWatMicPM)
    AlO3H3_mole_conc         = AZMAX1(ZALOH3/VLWatMicPM)
    AlO4H4_1e_aque_mole_conc = AZMAX1(ZALOH4/VLWatMicPM)
    AlSO4_1p_aque_mole_conc  = AZMAX1(ZALS/VLWatMicPM)
    FeOH_2p_aque_mole_conc   = AZMAX1(ZFEOH1/VLWatMicPM)
    FeO2H2_p_mole_conc       = AZMAX1(ZFEOH2/VLWatMicPM)
    FeO3H3_mole_conc         = AZMAX1(ZFEOH3/VLWatMicPM)
    FeO4H4_1e_aque_mole_conc = AZMAX1(ZFEOH4/VLWatMicPM)
    FeSO4_1p_aque_mole_conc  = AZMAX1(ZFES/VLWatMicPM)
    CaO2H2_mole_conc         = AZMAX1(ZCAO/VLWatMicPM)
    CaCO3_mole_conc          = AZMAX1(ZCAC/VLWatMicPM)
    CaHCO3_1p_aque_mole_conc = AZMAX1(ZCAH/VLWatMicPM)
    CaSO4_mole_conc               = AZMAX1(ZCAS/VLWatMicPM)
    MgOH_1p_aque_mole_conc   = AZMAX1(ZMGO/VLWatMicPM)
    MgCO3_conc               = AZMAX1(ZMGC/VLWatMicPM)
    MgHCO3_1p_aque_mole_conc = AZMAX1(ZMGH/VLWatMicPM)
    MgSO4_conc               = AZMAX1(ZMGS/VLWatMicPM)
    NaCO3_1e_aque_mole_conc  = AZMAX1(ZNAC/VLWatMicPM)
    NaSO4_1e_aque_mole_conc  = AZMAX1(ZNAS/VLWatMicPM)
    KSO4_1e_aque_mole_conc   = AZMAX1(ZKAS/VLWatMicPM)
  ELSE
    OH_1e_aque_mole_conc     = 0._r8
    Al_3p_aque_mole_conc     = 0._r8
    Fe_3p_aque_mole_conc     = 0._r8
    Ca_2p_aque_mole_conc     = 0._r8
    Mg_2p_aque_mole_conc     = 0._r8
    Na_1p_aque_mole_conc     = 0._r8
    K_1p_aque_mole_conc      = 0._r8
    SO4_2e_aque_mole_conc    = 0._r8
    Cl_e_conc      = 0._r8
    CO3_2e_aque_mole_conc    = 0._r8
    HCO3_e_mole_conc    = 0._r8
    H2CO3_aque_mole_conc = 0._r8
    AlOH_2p_aque_mole_conc   = 0._r8
    AlO2H2_1p_aque_mole_conc = 0._r8
    AlO3H3_mole_conc    = 0._r8
    AlO4H4_1e_aque_mole_conc = 0._r8
    AlSO4_1p_aque_mole_conc  = 0._r8
    FeOH_2p_aque_mole_conc   = 0._r8
    FeO2H2_p_mole_conc  = 0._r8
    FeO3H3_mole_conc    = 0._r8
    FeO4H4_1e_aque_mole_conc = 0._r8
    FeSO4_1p_aque_mole_conc  = 0._r8
    CaO2H2_mole_conc    = 0._r8
    CaCO3_mole_conc     = 0._r8
    CaHCO3_1p_aque_mole_conc = 0._r8
    CaSO4_mole_conc     = 0._r8
    MgOH_1p_aque_mole_conc   = 0._r8
    MgCO3_conc     = 0._r8
    MgHCO3_1p_aque_mole_conc = 0._r8
    MgSO4_conc     = 0._r8
    NaCO3_1e_aque_mole_conc  = 0._r8
    NaSO4_1e_aque_mole_conc  = 0._r8
    KSO4_1e_aque_mole_conc   = 0._r8
  ENDIF
!
!     PO4 CONCENTRATIONS IN NON-BAND AND BAND SOIL ZONES
!
!     VLWatMicPPO,VLWatMicPPB=water volume in PO4 non-band,band
!
  IF(VLWatMicPPO.GT.ZEROS2)THEN
    VLWatMicPPX      = patomw*VLWatMicPPO
    H0PO4_3e_conc    = AZMAX1(H0PO4/VLWatMicPPO)
    H3PO4_conc       = AZMAX1(H3PO4/VLWatMicPPO)
    FeHPO4_p_conc    = AZMAX1(ZFE1P/VLWatMicPPO)
    FeH2PO4_2p_aque_mole_conc  = AZMAX1(ZFE2P/VLWatMicPPO)
    CaPO4_1e_con     = AZMAX1(ZCA0P/VLWatMicPPO)
    CaHPO4_conc      = AZMAX1(ZCA1P/VLWatMicPPO)
    CaH4P2O8_1p_aque_mole_conc = AZMAX1(ZCA2P/VLWatMicPPO)
    MgHPO4_conc      = AZMAX1(ZMG1P/VLWatMicPPO)
  ELSE
    H0PO4_3e_conc    = 0._r8
    H3PO4_conc       = 0._r8
    FeHPO4_p_conc    = 0._r8
    FeH2PO4_2p_aque_mole_conc  = 0._r8
    CaPO4_1e_con     = 0._r8
    CaHPO4_conc      = 0._r8
    CaH4P2O8_1p_aque_mole_conc = 0._r8
    MgHPO4_conc      = 0._r8
  ENDIF
  IF(VLWatMicPPB.GT.ZEROS2)THEN
    H0PO4_3e_band_conc    = AZMAX1(H0POB/VLWatMicPPB)
    H3PO4_band_conc       = AZMAX1(H3POB/VLWatMicPPB)
    FeHPO4_1p_band_conc   = AZMAX1(ZFE1PB/VLWatMicPPB)
    FeH2PO4_2p_band_conc  = AZMAX1(ZFE2PB/VLWatMicPPB)
    CaPO4_1e_band_conc    = AZMAX1(ZCA0PB/VLWatMicPPB)
    CaHPO4_band_conc      = AZMAX1(ZCA1PB/VLWatMicPPB)
    CaH4P2O8_1p_band_conc = AZMAX1(ZCA2PB/VLWatMicPPB)
    MgHPO4_band_conc      = AZMAX1(ZMG1PB/VLWatMicPPB)
  ELSE
    H0PO4_3e_band_conc    = 0._r8
    H3PO4_band_conc       = 0._r8
    FeHPO4_1p_band_conc   = 0._r8
    FeH2PO4_2p_band_conc  = 0._r8
    CaPO4_1e_band_conc    = 0._r8
    CaHPO4_band_conc      = 0._r8
    CaH4P2O8_1p_band_conc = 0._r8
    MgHPO4_band_conc      = 0._r8
  ENDIF
!
!     EXCHANGEABLE ION CONCENTRATIONS
!
  IF(SoilMicPMassLayerX.GT.ZEROS)THEN
    XHY1         = AZMAX1(XHY/SoilMicPMassLayerX)
    XAl_conc     = AZMAX1(XAL/SoilMicPMassLayerX)
    XFe_conc     = AZMAX1(XFE/SoilMicPMassLayerX)
    XCa_conc     = AZMAX1(XCA/SoilMicPMassLayerX)
    XMg_conc     = AZMAX1(XMG/SoilMicPMassLayerX)
    XNa_conc     = AZMAX1(XNA/SoilMicPMassLayerX)
    XK_conc      = AZMAX1(XKA/SoilMicPMassLayerX)
    XHC1         = AZMAX1(XHC/SoilMicPMassLayerX)
    XAlO2H2_conc = AZMAX1(XALO2/SoilMicPMassLayerX)
    XFeO2H2_conc = AZMAX1(XFEO2/SoilMicPMassLayerX)
    XCOOH_mole_conc   = AZMAX1(COOH*ORGC/SoilMicPMassLayerX)
!
!     PRECIPITATE CONCENTRATIONS
!
    Precp_AlO3H3_conc = AZMAX1(PALOH/SoilMicPMassLayerX)
    Precp_FeO3H3_conc = AZMAX1(PFEOH/SoilMicPMassLayerX)
    Precp_CaCO3_conc  = AZMAX1(PCACO/SoilMicPMassLayerX)
    Precp_CaSO4_conc  = AZMAX1(PCASO/SoilMicPMassLayerX)
  ELSE
    XHY1              = 0._r8
    XAl_conc          = 0._r8
    XFe_conc          = 0._r8
    XCa_conc          = 0._r8
    XMg_conc          = 0._r8
    XNa_conc          = 0._r8
    XK_conc           = 0._r8
    XHC1              = 0._r8
    XAlO2H2_conc      = 0._r8
    XFeO2H2_conc      = 0._r8
    XCOOH_mole_conc        = 0._r8
    Precp_AlO3H3_conc = 0._r8
    Precp_FeO3H3_conc = 0._r8
    Precp_CaCO3_conc  = 0._r8
    Precp_CaSO4_conc  = 0._r8
  ENDIF
  end subroutine PrepIonConcentrations
!--------------------------------------------------------------------------
  subroutine SolveChemEquilibria

  implicit none
  real(r8) :: SP,SPX,S0,S1,P1,P2,PX
  integer :: M,NR1,NP2
!     begin_execution

  DO 1000 M=1,MRXN
!
!     SOLUTE CONCENTRATIONS
!
    call GetSoluteConcentrations
!
    call IonStrengthActivity
!
!     ALUMINUM HYDROXIDE (GIBBSITE)
!
    ! PX/PY copy an AMAX1 operand, so exact equality identifies the selected
    ! species even at trace activity; the first exact tie keeps precedence.
    PX=AMAX1(Al_3p_activity,AlOH_2p_activity,AlO2H2_1p_activity,AlO3H3_activity,AlO4H4_1e_activity)
    IF(PX.EQ.Al_3p_activity)THEN
      R1  = H_1p_activity
      P1  = Al_3p_activity
      P2  = OH_1e_activity
      NR1 = 3
      NP2 = 0
      SP  = SHALO
    ELSEIF(PX.EQ.AlOH_2p_activity)THEN
      R1  = H_1p_activity
      P1  = AlOH_2p_activity
      P2  = OH_1e_activity
      NR1 = 2
      NP2 = 0
      SP  = SHAL1
    ELSEIF(PX.EQ.AlO2H2_1p_activity)THEN
      R1  = H_1p_activity
      P1  = AlO2H2_1p_activity
      P2  = OH_1e_activity
      NR1 = 1
      NP2 = 0
      SP  = SHAL2
    ELSEIF(PX.EQ.AlO3H3_activity)THEN
      R1  = H_1p_activity
      P1  = AlO3H3_activity
      P2  = OH_1e_activity
      NR1 = 0
      NP2 = 0
      SP  = SPAL3
    ELSEIF(PX.EQ.AlO4H4_1e_activity)THEN
      R1  = OH_1e_activity
      P1  = AlO4H4_1e_activity
      P2  = H_1p_activity
      NR1 = 0
      NP2 = 1
      SP  = SHAL4
    ENDIF
    RHAL1  = 0._r8
    RHALO1 = 0._r8
    RHALO2 = 0._r8
    RHALO3 = 0._r8
    RHALO4 = 0._r8
! R1     = AMAX1(ZERO,R1)
!     P1=AMAX1(ZERO,P1)
!     P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1/P2**NP2
    RPALOX=AMAX1(-Precp_AlO3H3_conc,TPDX*(P1-SPX))
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
    PX=AMAX1(Fe_3p_activity,FeOH_2p_activity,FeO2H2_p_activity,FeO3H3_activity,FeO4H4_1e_activity)
    IF(PX.EQ.Fe_3p_activity)THEN
      R1  = H_1p_activity
      P1  = Fe_3p_activity
      P2  = OH_1e_activity
      NR1 = 3
      NP2 = 0
      SP  = SHFEO
    ELSEIF(PX.EQ.FeOH_2p_activity)THEN
      R1  = H_1p_activity
      P1  = FeOH_2p_activity
      P2  = OH_1e_activity
      NR1 = 2
      NP2 = 0
      SP  = SHFE1
    ELSEIF(PX.EQ.FeO2H2_p_activity)THEN
      R1  = H_1p_activity
      P1  = FeO2H2_p_activity
      P2  = OH_1e_activity
      NR1 = 1
      NP2 = 0
      SP  = SHFE2
    ELSEIF(PX.EQ.FeO3H3_activity)THEN
      R1  = H_1p_activity
      P1  = FeO3H3_activity
      P2  = OH_1e_activity
      NR1 = 0
      NP2 = 0
      SP  = SPFE3
    ELSEIF(PX.EQ.FeO4H4_1e_activity)THEN
      R1  = OH_1e_activity
      P1  = FeO4H4_1e_activity
      P2  = H_1p_activity
      NR1 = 0
      NP2 = 1
      SP  = SHFE4
    ENDIF
    RHFE1  = 0._r8
    RHFEO1 = 0._r8
    RHFEO2 = 0._r8
    RHFEO3 = 0._r8
    RHFEO4 = 0._r8
    R1     = AMAX1(ZERO,R1)
    P1     = AMAX1(ZERO,P1)
    P2     = AMAX1(ZERO,P2)
    SPX    = SP*R1**NR1/P2**NP2
    RPFEOX = AMAX1(-Precp_FeO3H3_conc,TPDX*(P1-SPX))
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
!     IF(I.EQ.180.AND.J.EQ.12)THEN
!     WRITE(*,1112)'FEOH',I,J,L,M,Precp_FeO3H3_conc,Fe_3p_activity,FeOH_2p_activity,FeO2H2_p_activity,FeO3H3_activity,FeO4H4_1e_activity
!    2,OH_1e_activity,R1,P1,P2,SP,SPX,RPFEOX,RHFE1,RHFEO1,RHFEO2,RHFEO3,RHFEO4
!    3,Fe_3p_activity*OH_1e_activity**3,SPFEO
!     ENDIF
!
!     CALCITE
!
    PX = AMAX1(CO3_2e_activity,HCO3_e_activity,H2CO3_activity)
    R1 = H_1p_activity
    P1 = Ca_2p_activity
    IF(PX.EQ.CO3_2e_activity)THEN
      P2  = CO3_2e_activity
      NR1 = 0
      SP  = SPCAC
    ELSEIF(PX.EQ.HCO3_e_activity)THEN
      P2  = HCO3_e_activity
      NR1 = 1
      SP  = SHCAC1
    ELSEIF(PX.EQ.H2CO3_activity)THEN
      P2  = H2CO3_activity
      NR1 = 2
      SP  = SHCAC2
    ENDIF
    RHCAC3 = 0._r8
    RHCACH = 0._r8
    RHCACO = 0._r8
    R1     = AMAX1(ZERO,R1)
    P1     = AMAX1(ZERO,P1)
    P2     = AMAX1(ZERO,P2)
    SPX    = SP*R1**NR1
    S0     = P1+P2
    ! Keep small positive discriminants: concentration cutoffs do not apply to squares.
    ! Mineral D=(P1-P2)**2+4*SPX; rationalize S0-sqrt(D) to retain small fluxes.
    S1     = MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    RPCACX = AMAX1(-Precp_CaCO3_conc,TPDX*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))  !negative flux, suction of CO2

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
    P1     = Ca_2p_activity
    P2     = SO4_2e_activity
    P1     = AMAX1(ZERO,P1)
    P2     = AMAX1(ZERO,P2)
    SPX    = SPCAS
    S0     = P1+P2
    S1     = MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    RPCASO = AMAX1(-Precp_CaSO4_conc,TPDX*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))
!     IF((M/10)*10.EQ.M)THEN
!     WRITE(*,1112)'CALC',I,J,L,M,Precp_CaSO4_conc,CO3_2e_activity,HCO3_e_activity,H2CO3_activity,H_1p_aque_mole_conc
!    2,OH_1e_aque_mole_conc,R1,P1,P2,P3,SP,Z,TX,RPCACX,RHCAC3,RHCACH,RHCACO
!    3,Ca_2p_aque_mole_conc*A2*CCO3*A2,SPCAC
!     ENDIF
!
!     PHOSPHORUS PRECIPITATION-DISSOLUTION IN NON-BAND SOIL ZONE
!
    call PhospPrecipDissolNonBand
!
!     PHOSPHORUS PRECIPITATION-DISSOLUTION IN BAND SOIL ZONE
!
    call PhospPrecipDissolBand
!
    call PhospAnionExchNoBand
!
    call PhospAnionExchBand
!
!     CATION EXCHANGE FROM GAPON SELECTIVITY COEFFICIENTS
!     FOR CA-NH4, CA-H, CA-AL, CA-MG, CA-NA, CA-K
!
    call CationExchange
!
!     SOLUTE DISSOCIATION REACTIONS
!
    call SoluteDissociation
!
!     TOTAL ION FLUXES FOR CURRENT ITERATION
!     FROM ALL REACTIONS ABOVE
!
    call UpdateIonFluxCurentIter

    call UpdateIonConcCurrentIter
!
    call AccumulateIonFlux
!
!     GO TO NEXT ITERATION
!
1000  CONTINUE
!     CONVERGENCE ITERATIONS CO2CompenPoint_nodeETED
!
!     IF(J.EQ.24)THEN
!     WRITE(*,1119)'GAPON',I,J,L,M,H0PO4_3e_conc,Al_3p_aque_mole_conc,Fe_3p_aque_mole_conc,H0PO4_3e_conc*A3*Al_3p_aque_mole_conc*A3
!    2,SPALP,H0PO4_3e_conc*A3*Fe_3p_aque_mole_conc*A3,SPFEP
!    6,SPOH2,XROH1_mole_conc*H_1p_aque_mole_conc*A1/XROH2_mole_conc,SPOH1,XOH_mole_conc*H_1p_aque_mole_conc*A1/XROH1_mole_conc
!    7,SPH2P,XROH2_mole_conc*H2PO4_1e_aque_mole_conc*A1/XH2PO4_mole_conc,SXH2P,XROH1_mole_conc*H2PO4_1e_aque_mole_conc/(XH2PO4_mole_conc*OH_1e_aque_mole_conc)
!    8,SPH1P,XROH1_mole_conc*H1PO4_2e_aque_mole_conc*A2/(XHPO4_mole_conc*OH_1e_aque_mole_conc*A1)
!    9,OH_1e_aque_mole_conc*A1,H_1p_aque_mole_conc*A1
!1119  FORMAT(A8,4I4,24E11.3)
!     WRITE(*,1119)'CATION',I,J,L,M,CEC_conc,XNH4_mole_conc+XHY1+3*XAl_conc+2*(XCa_conc+XMg_conc)
!    2+XNa_conc+XK_conc,XNH4_mole_conc,XHY1,XAl_conc,XCa_conc,XMg_conc,XNa_conc,XK_conc,NH4_1p_aque_mole_conc,H_1p_aque_mole_conc,Al_3p_aque_mole_conc,Ca_2p_aque_mole_conc
!    2,Mg_2p_aque_mole_conc,Na_1p_aque_mole_conc,K_1p_aque_mole_conc,(Ca_2p_aque_mole_conc*A2)**0.5*XNH4_mole_conc/(NH4_1p_aque_mole_conc*A1*XCa_conc*2)
!    3,(Ca_2p_aque_mole_conc*A2)**0.5*XHY1/(H_1p_aque_mole_conc*A1*XCa_conc*2)
!    2,(Ca_2p_aque_mole_conc*A2)**0.5*XAl_conc*3/((Al_3p_aque_mole_conc*A3)**0.333*XCa_conc*2)
!    3,(Ca_2p_aque_mole_conc*A2)**0.5*XMg_conc*2/((Mg_2p_aque_mole_conc*A2)**0.5*XCa_conc*2)
!    3,(Ca_2p_aque_mole_conc*A2)**0.5*XNa_conc/(Na_1p_aque_mole_conc*A1*XCa_conc*2)
!    5,(Ca_2p_aque_mole_conc*A2)**0.5*XK_conc/(K_1p_aque_mole_conc*A1*XCa_conc*2)
!    6,H_1p_aque_mole_conc*A1*XCOO/XHC1,AlO2H2_1p_aque_mole_conc*A1*XCOO/XAlO2H2_conc
!     ENDIF
  end subroutine SolveChemEquilibria
!------------------------------------------------------------------------------------------

  subroutine PhospPrecipDissolNonBand

  implicit none
  real(r8) :: SP,SPX,S0,S1,P1,P2,P3,PX,PY
  integer :: NR1,NP3
  IF(VLWatMicPPO.GT.ZEROS2)THEN
!
!     ALUMINUM PHOSPHATE (VARISCITE)
!
    ! PX/PY copy an AMAX1 operand, so exact equality identifies the selected
    ! species even at trace activity; the first exact tie keeps precedence.
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
    RHA0P1 = 0._r8
    RHA1P1 = 0._r8
    RHA2P1 = 0._r8
    RHA3P1 = 0._r8
    RHA4P1 = 0._r8
    RHA0P2 = 0._r8
    RHA1P2 = 0._r8
    RHA2P2 = 0._r8
    RHA3P2 = 0._r8
    RHA4P2 = 0._r8
    R1     = AMAX1(ZERO,R1)
    P1     = AMAX1(ZERO,P1)
    P2     = AMAX1(ZERO,P2)
    P3     = AMAX1(ZERO,P3)
    SPX    = SP*R1**NR1/P3**NP3
    S0     = P1+P2
    S1     = MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    H2PO4_1e_AlPO4_dissol_flx=AMAX1(-Precp_AlPO4_mole_conc,TPDX*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))
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
!     IF(I.EQ.180.AND.J.EQ.12)THEN
!     WRITE(*,1112)'ALPO4',I,J,L,M,Precp_AlPO4_mole_conc,Al_3p_activity,AlOH_2p_activity,AlO2H2_1p_activity,AlO3H3_activity,AlO4H4_1e_activity
!    2,H0PO4_3e_activity,H1PO4_2e_activity,H2PO4_1e_activity,H_1p_activity,OH_1e_activity,H2PO4_1e_AlPO4_dissol_flx,RHA0P1,RHA1P1,RHA2P1,RHA3P1
!    3,RHA4P1,RHA0P2,RHA1P2,RHA2P2,RHA3P2,RHA4P2,SP,SPX,Al_3p_activity*H0PO4_3e_activity
!    4,SPALP,H0PO4_3e_conc,H1PO4_2e_aque_mole_conc,H2PO4_1e_aque_mole_conc
!     ENDIF
!1112  FORMAT(A8,4I5,80E12.4)
!     ENDIF
!
!     IRON PHOSPHATE (STRENGITE)
!
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
    RHF0P1 = 0._r8
    RHF1P1 = 0._r8
    RHF2P1 = 0._r8
    RHF3P1 = 0._r8
    RHF4P1 = 0._r8
    RHF0P2 = 0._r8
    RHF1P2 = 0._r8
    RHF2P2 = 0._r8
    RHF3P2 = 0._r8
    RHF4P2 = 0._r8
    R1     = AMAX1(ZERO,R1)
    P1     = AMAX1(ZERO,P1)
    P2     = AMAX1(ZERO,P2)
    P3     = AMAX1(ZERO,P3)
    SPX    = SP*R1**NR1/P3**NP3
    S0     = P1+P2
    S1     = MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    H2PO4_1e_FePO4_dissol_flx=AMAX1(-Precp_FePO4_mole_conc,TPDX*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))
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
!     IF(I.EQ.180.AND.J.EQ.12)THEN
!     WRITE(*,1112)'FEPO4',I,J,L,M,Precp_FePO4_mole_conc,Fe_3p_activity,FeOH_2p_activity, &
!       FeO2H2_p_activity,FeO3H3_activity,FeO4H4_1e_activity,&
!       H0PO4_3e_activity,H1PO4_2e_activity,H2PO4_1e_activity,H_1p_activity,OH_1e_activity,&
!       H2PO4_1e_FePO4_dissol_flx,RHF0P1,RHF1P1,RHF2P1,RHF3P1,&
!       RHF4P1,RHF0P2,RHF1P2,RHF2P2,RHF3P2,RHF4P2,SP,SPX,Fe_3p_activity*H0PO4_3e_activity,&
!       SPFEP
!     ENDIF
!
!     DICALCIUM PHOSPHATE
!
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
    RPCAD1 = 0._r8
    RHCAD2 = 0._r8
    R1     = AMAX1(ZERO,R1)
    P1     = AMAX1(ZERO,P1)
    P2     = AMAX1(ZERO,P2)
    SPX    = SP*R1**NR1
    S0     = P1+P2
    S1     = MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    H2PO4_1e_CaHPO4_dissol_flx=AMAX1(-Precp_CaHPO4_mole_conc,TPDX*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))
    IF(PX.EQ.H1PO4_2e_activity)THEN
      RPCAD1=H2PO4_1e_CaHPO4_dissol_flx
    ELSEIF(PX.EQ.H2PO4_1e_activity)THEN
      RHCAD2=H2PO4_1e_CaHPO4_dissol_flx
    ENDIF
!     IF((M/10)*10.EQ.M)THEN
!     WRITE(*,1112)'CAPO4',I,J,L,M,Precp_CaH4P2O8_mole_conc,Precp_CaHPO4_mole_conc,Ca_2p_aque_mole_conc
!    2,H1PO4_2e_aque_mole_conc,H2PO4_1e_aque_mole_conc,H_1p_aque_mole_conc,OH_1e_aque_mole_conc,H2PO4_1e_CaHPO4_dissol_flx,RPCAD1,RHCAD2,R1,P1,P2,P3
!    3,SP,Z,FX,Y,X,TX,A2,Ca_2p_aque_mole_conc*A2*H1PO4_2e_aque_mole_conc*A2,SPCAD
!     ENDIF
!
!     HYDROXYAPATITE
!
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
    RHCAH1 = 0._r8
    RHCAH2 = 0._r8
    R1     = AMAX1(ZERO,R1)
    P1     = AMAX1(ZERO,P1)
    P2     = AMAX1(ZERO,P2)
    SPX    = (SP*R1**NR1/P1**5)**0.333_r8
    H2PO4_1e_apatite_dissol_flx=AMAX1(-Precp_Ca5P3O12O3H3_mole_conc,TPDX*(P2-SPX))
    IF(PX.EQ.H1PO4_2e_activity)THEN
      RHCAH1=H2PO4_1e_apatite_dissol_flx
    ELSEIF(PX.EQ.H2PO4_1e_activity)THEN
      RHCAH2=H2PO4_1e_apatite_dissol_flx
    ENDIF
!     IF((I/10)*10.EQ.I.AND.J.EQ.12)THEN
!     WRITE(*,1112)'A1',I,L,K,M,A1,A2,A3,FSTR2,CSTR1
!    2,CSTR2,cation_3p_aque_mole_conc,anion_3e_conc,cation_2p_aque_mole_conc,CA2,cation_1p_aque_mole_conc,CA1,VLWatMicPM
!     WRITE(*,1112)'APATITE',I,J,L,M,Precp_Ca5P3O12O3H3_mole_conc,Ca_2p_activity,XCa_conc
!    2,H0PO4_3e_activity,H1PO4_2e_activity,H2PO4_1e_activity,H_1p_activity,OH_1e_activity,H2PO4_1e_apatite_dissol_flx,RHCAH1,RHCAH2
!    3,SP,SPX,Ca_2p_activity**5*H0PO4_3e_activity**3*OH_1e_activity,SPCAH,SHCAH1,SHCAH2
!    3,H0PO4_3e_conc,H1PO4_2e_aque_mole_conc,H2PO4_1e_aque_mole_conc,XOH_mole_conc,XROH1_mole_conc,XROH2_mole_conc,XHPO4_mole_conc,XH2PO4_mole_conc
!    4,RHA0P1,RHA1P1,RHA2P1,RHA3P1
!    2,RHA4P1,RHF0P1,RHF1P1,RHF2P1
!    3,RHF3P1,RHF4P1,RPCAD1,3.0_r8*RHCAH1
!    4,H1PO4_to_XHPO4_ROH_flx,RH1P,H2PO4_e_to_HPO4_2e_flx,RF1P,RC1P,RM1P
!    5,RHA0P2,RHA1P2,RHA2P2,RHA3P2
!    2,RHA4P2,RHF0P2,RHF1P2,RHF2P2
!    3,RHF3P2,RHF4P2,RHCAD2,3.0_r8*RHCAH2
!    4,H2PO4_1e_to_XH2PO4_ROH2_flx,H2PO4_to_XH2PO4_ROH_flx,H2PO4_e_to_HPO4_2e_flx,RH3P,RF2P,RC2P,RH3P
!     ENDIF
!
!     MONOCALCIUM PHOSPHATE
!
    P1                           = Ca_2p_activity
    P2                           = H2PO4_1e_activity
    P1                           = AMAX1(ZERO,P1)
    P2                           = AMAX1(ZERO,P2)
    ! Ca(H2PO4)2(s) <-> Ca(2+)+2H2PO4(-): K=a_Ca*a_H2PO4**2.
    ! Relax phosphate toward equilibrium at the current calcium activity.
    ! Positive flux precipitates one formula unit; transfers use one Ca and two P.
    SPX                          = SQRT(SPCAM/P1)
    H2PO4_1e_CaH4P2O8_dissol_flx = AMAX1(-Precp_CaH4P2O8_mole_conc,TPDX*(P2-SPX))
  ELSE
    H2PO4_1e_AlPO4_dissol_flx    = 0._r8
    H2PO4_1e_FePO4_dissol_flx    = 0._r8
    H2PO4_1e_CaHPO4_dissol_flx   = 0._r8
    H2PO4_1e_apatite_dissol_flx  = 0._r8
    RHA0P1                       = 0._r8
    RHA1P1                       = 0._r8
    RHA2P1                       = 0._r8
    RHA3P1                       = 0._r8
    RHA4P1                       = 0._r8
    RHA0P2                       = 0._r8
    RHA1P2                       = 0._r8
    RHA2P2                       = 0._r8
    RHA3P2                       = 0._r8
    RHA4P2                       = 0._r8
    RHF0P1                       = 0._r8
    RHF1P1                       = 0._r8
    RHF2P1                       = 0._r8
    RHF3P1                       = 0._r8
    RHF4P1                       = 0._r8
    RHF0P2                       = 0._r8
    RHF1P2                       = 0._r8
    RHF2P2                       = 0._r8
    RHF3P2                       = 0._r8
    RHF4P2                       = 0._r8
    RPCAD1                       = 0._r8
    RHCAD2                       = 0._r8
    RHCAH1                       = 0._r8
    RHCAH2                       = 0._r8
    H2PO4_1e_CaH4P2O8_dissol_flx = 0._r8
  ENDIF
  end subroutine PhospPrecipDissolNonBand

!----------------------------------------------------------------------------
  subroutine SummarizeIonFluxes

  implicit none
  !     begin_execution
  !     CONVERT TOTAL ION FLUXES FROM CHANGES IN CONCENTRATION
  !     TO CHANGES IN MASS PER UNIT AREA FOR USE IN 'REDIST'
  !
  TRChem_NH4_soil_mole      = TRChem_NH4_soil_mole*VLWatMicPNH
  TRChem_NH4_band_soil_mole = TRChem_NH4_band_soil_mole*VLWatMicPNB
  TRChem_NH3_soil_mole      = TRChem_NH3_soil_mole*VLWatMicPNH
  TRChem_NH3_band_soil      = TRChem_NH3_band_soil*VLWatMicPNB
  TRChem_Al_3p_soil         = TRChem_Al_3p_soil*VLWatMicPM
  TRChem_Fe_3p_soil         = TRChem_Fe_3p_soil*VLWatMicPM
  TRChem_H_p_soil           = TRChem_H_p_soil*VLWatMicPM
  TRChem_Ca_2p_soil         = TRChem_Ca_2p_soil*VLWatMicPM
  TRChem_Mg_2p_soil         = TRChem_Mg_2p_soil*VLWatMicPM
  TRChem_Na_p_soil          = TRChem_Na_p_soil*VLWatMicPM
  TRChem_K_1p_soil          = TRChem_K_1p_soil*VLWatMicPM
  TRChem_OH_1e_soil         = TRChem_OH_1e_soil*VLWatMicPM
  TRChem_SO4_2e_soil        = TRChem_SO4_2e_soil*VLWatMicPM
  TRChem_CO3_2e_soil        = TRChem_CO3_2e_soil*VLWatMicPM
  TRChem_HCO3_soil          = TRChem_HCO3_soil*VLWatMicPM
  TRChem_CO2_gchem_soil     = TRChem_CO2_gchem_soil*VLWatMicPM
  TRChem_AlOH_soil          = TRChem_AlOH_soil*VLWatMicPM
  TRChem_AlO2H2_soil        = TRChem_AlO2H2_soil*VLWatMicPM
  TRChem_AlO3H3_soil        = TRChem_AlO3H3_soil*VLWatMicPM
  TRChem_AlO4H4_soil               = TRChem_AlO4H4_soil*VLWatMicPM
  TRChem_AlSO4_soil                = TRChem_AlSO4_soil*VLWatMicPM
  TRChem_FeOH_soil                 = TRChem_FeOH_soil*VLWatMicPM
  TRChem_FeO2H2_soil               = TRChem_FeO2H2_soil*VLWatMicPM
  TRChem_FeO3H3_soil               = TRChem_FeO3H3_soil*VLWatMicPM
  TRChem_FeO4H4_soil               = TRChem_FeO4H4_soil*VLWatMicPM
  TRChem_FeSO4_soil                = TRChem_FeSO4_soil*VLWatMicPM
  TRChem_CaOH_soil                 = TRChem_CaOH_soil*VLWatMicPM
  TRChem_CaCO3_soil                = TRChem_CaCO3_soil*VLWatMicPM
  TRChem_CaHCO3_soil               = TRChem_CaHCO3_soil*VLWatMicPM
  TRChem_CaSO4_soil                = TRChem_CaSO4_soil*VLWatMicPM
  TRChem_MgOH_soil                 = TRChem_MgOH_soil*VLWatMicPM
  TRChem_MgCO3_soil                = TRChem_MgCO3_soil*VLWatMicPM
  TRChem_MgHCO3_soil               = TRChem_MgHCO3_soil*VLWatMicPM
  TRChem_MgSO4_soil                = TRChem_MgSO4_soil*VLWatMicPM
  TRChem_NaCO3_soil                = TRChem_NaCO3_soil*VLWatMicPM
  TRChem_NaSO4_soil                = TRChem_NaSO4_soil*VLWatMicPM
  TRChem_KSO4_soil                 = TRChem_KSO4_soil*VLWatMicPM
  TRChem_PO4_soil                  = TRChem_PO4_soil*VLWatMicPPO
  TRChem_H1PO4_soil                = TRChem_H1PO4_soil*VLWatMicPPO
  TRChem_H2PO4_soil                = TRChem_H2PO4_soil*VLWatMicPPO
  TRChem_H3PO4_sorbed_soil         = TRChem_H3PO4_sorbed_soil*VLWatMicPPO
  TRChem_FeHPO4_soil               = TRChem_FeHPO4_soil*VLWatMicPPO
  TRChem_FeH2PO4_soil              = TRChem_FeH2PO4_soil*VLWatMicPPO
  TRChem_CaPO4_soil                = TRChem_CaPO4_soil*VLWatMicPPO
  TRChem_CaHPO4_soil               = TRChem_CaHPO4_soil*VLWatMicPPO
  TRChem_CaH4P2O8_soil             = TRChem_CaH4P2O8_soil*VLWatMicPPO
  TRChem_MgHPO4_soil               = TRChem_MgHPO4_soil*VLWatMicPPO
  TRChem_PO4_band_soil             = TRChem_PO4_band_soil*VLWatMicPPB
  TRChem_H1PO4_band_soil           = TRChem_H1PO4_band_soil*VLWatMicPPB
  TRChem_H2PO4_band_soil           = TRChem_H2PO4_band_soil*VLWatMicPPB
  TRChem_H3PO4_band_soil           = TRChem_H3PO4_band_soil*VLWatMicPPB
  TRChem_FeHPO4_band_soil          = TRChem_FeHPO4_band_soil*VLWatMicPPB
  TRChem_FeH2PO4_band_soil         = TRChem_FeH2PO4_band_soil*VLWatMicPPB
  TRChem_CaPO4_band_soil           = TRChem_CaPO4_band_soil*VLWatMicPPB
  TRChem_CaHPO4_band_soil          = TRChem_CaHPO4_band_soil*VLWatMicPPB
  TRChem_CaH4P2O8_band_soil        = TRChem_CaH4P2O8_band_soil*VLWatMicPPB
  TRChem_MgHPO4_band_soil          = TRChem_MgHPO4_band_soil*VLWatMicPPB
  TRChem_NH4_sorbed_soil           = TRChem_NH4_sorbed_soil*VLWatMicPNH
  TRChem_NH4_sorbed_band_soil      = TRChem_NH4_sorbed_band_soil*VLWatMicPNB
  TRChem_H_p_sorbed_soil           = TRChem_H_p_sorbed_soil*VLWatMicPM
  TRChem_Al_sorbed_soil            = TRChem_Al_sorbed_soil*VLWatMicPM
  TRChem_Fe_sorbed_soil            = TRChem_Fe_sorbed_soil*VLWatMicPM
  TRChem_Ca_sorbed_soil            = TRChem_Ca_sorbed_soil*VLWatMicPM
  TRChem_Mg_sorbed_soil            = TRChem_Mg_sorbed_soil*VLWatMicPM
  TRChem_Na_sorbed_soil            = TRChem_Na_sorbed_soil*VLWatMicPM
  TRChem_K_sorbed_soil             = TRChem_K_sorbed_soil*VLWatMicPM
  TRChem_HCO3_sorbed_soil          = TRChem_HCO3_sorbed_soil*VLWatMicPM
  TRChem_AlO2H2_sorbed_soil        = TRChem_AlO2H2_sorbed_soil*VLWatMicPM
  TRChem_FeO2H2_sorbed_soil        = TRChem_FeO2H2_sorbed_soil*VLWatMicPM
  TRChem_RO_sorbed_soil            = TRChem_RO_sorbed_soil*VLWatMicPPO
  TRChem_ROH_sorbed_soil           = TRChem_ROH_sorbed_soil*VLWatMicPPO
  TRChem_ROH2_sorbed_soil          = TRChem_ROH2_sorbed_soil*VLWatMicPPO
  TRChem_RHPO4_sorbed_soil         = TRChem_RHPO4_sorbed_soil*VLWatMicPPO
  TRChem_RH2PO4_sorbed_soil        = TRChem_RH2PO4_sorbed_soil*VLWatMicPPO
  TRChem_RO_sorbed_band_soil       = TRChem_RO_sorbed_band_soil*VLWatMicPPB
  TRChem_ROH_sorbed_band_soil      = TRChem_ROH_sorbed_band_soil*VLWatMicPPB
  TRChem_ROH2_sorbed_band_soil     = TRChem_ROH2_sorbed_band_soil*VLWatMicPPB
  TRChem_RHPO4_sorbed_band_soil    = TRChem_RHPO4_sorbed_band_soil*VLWatMicPPB
  TRChem_RH2PO4_sorbed_band_soil   = TRChem_RH2PO4_sorbed_band_soil*VLWatMicPPB
  TRChem_AlOH3_precip_soil         = TRChem_AlOH3_precip_soil*VLWatMicPM
  TRChem_FeOH3_precip_soil         = TRChem_FeOH3_precip_soil*VLWatMicPM
  TRChem_CaCO3_precip_soil         = TRChem_CaCO3_precip_soil*VLWatMicPM
  TRChem_CaSO4_precip_soil         = TRChem_CaSO4_precip_soil*VLWatMicPM
  TRChem_AlPO4_precip_soil         = TRChem_AlPO4_precip_soil*VLWatMicPPO
  TRChem_FePO4_precip_soil         = TRChem_FePO4_precip_soil*VLWatMicPPO
  TRChem_CaHPO4_precip_soil        = TRChem_CaHPO4_precip_soil*VLWatMicPPO
  TRChem_apatite_precip_soil       = TRChem_apatite_precip_soil*VLWatMicPPO
  TRChem_CaH4P2O8_precip_soil      = TRChem_CaH4P2O8_precip_soil*VLWatMicPPO
  TRChem_AlPO4_precip_band_soil    = TRChem_AlPO4_precip_band_soil*VLWatMicPPB
  TRChem_FePO4_precip_band_soil    = TRChem_FePO4_precip_band_soil*VLWatMicPPB
  TRChem_CaHPO4_precip_band_soil   = TRChem_CaHPO4_precip_band_soil*VLWatMicPPB
  TRChem_apatite_precip_band_soil  = TRChem_apatite_precip_band_soil*VLWatMicPPB
  TRChem_CaH4P2O8_precip_band_soil = TRChem_CaH4P2O8_precip_band_soil*VLWatMicPPB
!
!     BOUNDARY SALT FLUXES FOR C, H, OH, P, AL+FE, CA, OH
!     USED TO CHECK MATERIAL BALANCES IN REDIST.F
!
!     TRAquaChem_CO2_soil=CO2 net change from all solute equilibria
!     TRH2O_soil=H2O net change from all solute equilibria
!     TBION_soil=total solute net change from all solute equilibria
!
  TRAquaChem_CO2_soil=TRChem_CO3_2e_soil+TRChem_CaCO3_soil+TRChem_MgCO3_soil+TRChem_NaCO3_soil+TRChem_CaCO3_precip_soil &
    +2.0_r8*(TRChem_HCO3_soil+TRChem_CaHCO3_soil+TRChem_MgHCO3_soil)
  TRH2O_soil=TRChem_H_p_soil+TRChem_OH_1e_soil+TRChem_H_p_sorbed_soil+TRChem_HCO3_sorbed_soil
  TBION_soil=4.0_r8*(TRChem_H3PO4_sorbed_soil+TRChem_H3PO4_band_soil) &
    +3.0_r8*(TRChem_FeH2PO4_soil+TRChem_CaH4P2O8_soil &
    +TRChem_FeH2PO4_band_soil+TRChem_CaH4P2O8_band_soil) &
    +2.0_r8*(TRChem_FeHPO4_soil+TRChem_CaHPO4_soil+TRChem_MgHPO4_soil &
    +TRChem_FeHPO4_band_soil+TRChem_CaHPO4_band_soil+TRChem_MgHPO4_band_soil) &
    +TRChem_PO4_soil+TRChem_CaPO4_soil+TRChem_PO4_band_soil+TRChem_CaPO4_band_soil &
    -(TRChem_AlPO4_precip_soil+TRChem_FePO4_precip_soil &
    +TRChem_AlPO4_precip_band_soil+TRChem_FePO4_precip_band_soil) &
    -(TRChem_CaHPO4_precip_soil+TRChem_CaH4P2O8_precip_soil &
    +TRChem_CaHPO4_precip_band_soil+TRChem_CaH4P2O8_precip_band_soil &
    +5.0_r8*(TRChem_apatite_precip_soil+TRChem_apatite_precip_band_soil)) &
    +TRChem_AlOH_soil+TRChem_FeOH_soil &
    +2.0_r8*(TRChem_AlO2H2_soil+TRChem_FeO2H2_soil) &
    +3.0_r8*(TRChem_AlO3H3_soil+TRChem_FeO3H3_soil) &
    +4.0_r8*(TRChem_AlO4H4_soil+TRChem_FeO4H4_soil) &
    +TRChem_CaOH_soil+TRChem_MgOH_soil &
    +3.0_r8*(TRChem_AlOH3_precip_soil+TRChem_FeOH3_precip_soil)
!     IF(L.EQ.11)THEN
!     WRITE(*,1111)'TRChem_CO2_gchem_soil',I,J,L,M,TRChem_CO2_gchem_soil,TRChem_CO3_2e_soil
!    2,TRChem_HCO3,TRChem_CaCO3_soil,TRChem_MgCO3_soil
!    2,TRChem_NaCO3_soil,TRChem_CaHCO3_soil
!    2,TRChem_MgHCO3_soil,TRChem_CaCO3_precip_soil,VLWatMicPM,RCO2
!    3,RHCO,RHCACH,RCO2Q,RCAH,RMGH,RHCO3,H_1p_activity,HCO3_e_activity,H2CO3_activity,DPCO2
!     WRITE(*,1111)'TBION_soil',I,J,L,M,TBION_soil
!     ENDIF
  end subroutine SummarizeIonFluxes
!------------------------------------------------------------------------------------------

  subroutine GetSoluteConcentrations
  implicit none
  character(len=*), parameter :: subname='GetSoluteConcentrations'
!     begin_execution
  call PrintInfo('beg '//subname)
  NH4_1p_aque_mole_conc      = AMAX1(ZERO,NH4_1p_aque_mole_conc)
  NH4_1p_band_mole_conc      = AMAX1(ZERO,NH4_1p_band_mole_conc)
  NH3_aque_mole_conc         = AMAX1(ZERO,NH3_aque_mole_conc)
  NH3_aqu_band_mole_conc     = AMAX1(ZERO,NH3_aqu_band_mole_conc)
  H_1p_aque_mole_conc        = AMAX1(ZERO,10.0_r8**(-(PH-3.0)))
  ! Retain incoming OH(-) until its activity correction can enter the flux ledger.
  OH_1e_aque_mole_conc       = AMAX1(ZERO,OH_1e_aque_mole_conc)
  CO3_2e_aque_mole_conc      = AMAX1(ZERO,CO3_2e_aque_mole_conc)
  Al_3p_aque_mole_conc       = AMAX1(ZERO,Al_3p_aque_mole_conc)
  Fe_3p_aque_mole_conc       = AMAX1(ZERO,Fe_3p_aque_mole_conc)
  Ca_2p_aque_mole_conc       = AMAX1(ZERO,AMIN1(CCAMX,Ca_2p_aque_mole_conc))
  Mg_2p_aque_mole_conc       = AMAX1(ZERO,Mg_2p_aque_mole_conc)
  Na_1p_aque_mole_conc       = AMAX1(ZERO,Na_1p_aque_mole_conc)
  K_1p_aque_mole_conc        = AMAX1(ZERO,K_1p_aque_mole_conc)
  SO4_2e_aque_mole_conc      = AMAX1(ZERO,SO4_2e_aque_mole_conc)
  HCO3_e_mole_conc           = AMAX1(ZERO,HCO3_e_mole_conc)
  H2CO3_aque_mole_conc       = AMAX1(ZERO,H2CO3_aque_mole_conc)
  AlOH_2p_aque_mole_conc     = AMAX1(ZERO,AlOH_2p_aque_mole_conc)
  AlO2H2_1p_aque_mole_conc   = AMAX1(ZERO,AlO2H2_1p_aque_mole_conc)
  AlO3H3_mole_conc                = AMAX1(ZERO,AlO3H3_mole_conc)
  AlO4H4_1e_aque_mole_conc   = AMAX1(ZERO,AlO4H4_1e_aque_mole_conc)
  AlSO4_1p_aque_mole_conc    = AMAX1(ZERO,AlSO4_1p_aque_mole_conc)
  FeOH_2p_aque_mole_conc     = AMAX1(ZERO,FeOH_2p_aque_mole_conc)
  FeO2H2_p_mole_conc              = AMAX1(ZERO,FeO2H2_p_mole_conc)
  FeO3H3_mole_conc                = AMAX1(ZERO,FeO3H3_mole_conc)
  FeO4H4_1e_aque_mole_conc   = AMAX1(ZERO,FeO4H4_1e_aque_mole_conc)
  FeSO4_1p_aque_mole_conc    = AMAX1(ZERO,FeSO4_1p_aque_mole_conc)
  CaO2H2_mole_conc                = AMAX1(ZERO,CaO2H2_mole_conc)
  CaCO3_mole_conc                 = AMAX1(ZERO,CaCO3_mole_conc)
  CaHCO3_1p_aque_mole_conc   = AMAX1(ZERO,CaHCO3_1p_aque_mole_conc)
  CaSO4_mole_conc                 = AMAX1(ZERO,CaSO4_mole_conc)
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
  H0PO4_3e_band_conc         = AMAX1(ZERO,H0PO4_3e_band_conc)
  H1PO4_2e_band_mole_conc    = AMAX1(ZERO,H1PO4_2e_band_mole_conc)
  H2PO4_1e_band_mole_conc    = AMAX1(ZERO,H2PO4_1e_band_mole_conc)
  H3PO4_band_conc            = AMAX1(ZERO,H3PO4_band_conc)
  FeHPO4_1p_band_conc        = AMAX1(ZERO,FeHPO4_1p_band_conc)
  FeH2PO4_2p_band_conc       = AMAX1(ZERO,FeH2PO4_2p_band_conc)
  CaPO4_1e_band_conc         = AMAX1(ZERO,CaPO4_1e_band_conc)
  CaHPO4_band_conc           = AMAX1(ZERO,CaHPO4_band_conc)
  CaH4P2O8_1p_band_conc      = AMAX1(ZERO,CaH4P2O8_1p_band_conc)
  MgHPO4_band_conc           = AMAX1(ZERO,MgHPO4_band_conc)
  XCOO                       = AZMAX1(XCOOH_mole_conc-XHC1-XAlO2H2_conc-XFeO2H2_conc)
  call PrintInfo('end '//subname)
  end subroutine GetSoluteConcentrations
!------------------------------------------------------------------------------------------

  subroutine IonStrengthActivity
  implicit none
  character(len=*), parameter :: subname='IonStrengthActivity'
  real(r8) :: cation_3p_aque_mole_conc,anion_3e_conc,cation_2p_aque_mole_conc
  real(r8) :: anion_2e_aque_mole_conc,cation_1p_aque_mole_conc,anion_1e_aque_mole_conc,CSTR1,CSTR2
  real(r8) :: A1,A3,FSTR2
  !     begin_execution
  call PrintInfo('beg '//subname)
  !     IONIC STRENGTH FROM SUMS OF ION CONCENTRATIONS
  !
  !     cation_3p_aque_mole_conc,anion_3e_conc,cation_2p_aque_mole_conc,anion_2e_aque_mole_conc,cation_1p_aque_mole_conc,anion_1e_aque_mole_conc=total tri-,di-,univalent cations C,anions A
  !     CSTR1=ion strength
  !
  cation_3p_aque_mole_conc=Al_3p_aque_mole_conc+Fe_3p_aque_mole_conc
  anion_3e_conc=H0PO4_3e_conc*VLPO4+H0PO4_3e_band_conc*VLPOB
  cation_2p_aque_mole_conc=Ca_2p_aque_mole_conc+Mg_2p_aque_mole_conc+AlOH_2p_aque_mole_conc+FeOH_2p_aque_mole_conc+FeH2PO4_2p_aque_mole_conc*VLPO4+FeH2PO4_2p_band_conc*VLPOB
  anion_2e_aque_mole_conc=SO4_2e_aque_mole_conc+CO3_2e_aque_mole_conc+H1PO4_2e_aque_mole_conc*VLPO4+H1PO4_2e_band_mole_conc*VLPOB
  cation_1p_aque_mole_conc=NH4_1p_aque_mole_conc*VLNH4+NH4_1p_band_mole_conc*VLNHB+H_1p_aque_mole_conc+Na_1p_aque_mole_conc+K_1p_aque_mole_conc &
    +AlO2H2_1p_aque_mole_conc+FeO2H2_p_mole_conc+AlSO4_1p_aque_mole_conc+FeSO4_1p_aque_mole_conc+CaO2H2_mole_conc+CaHCO3_1p_aque_mole_conc+MgOH_1p_aque_mole_conc+MgHCO3_1p_aque_mole_conc &
    +(FeHPO4_p_conc+CaH4P2O8_1p_aque_mole_conc)*VLPO4+(FeHPO4_1p_band_conc+CaH4P2O8_1p_band_conc)*VLPOB
  anion_1e_aque_mole_conc=NO3_1e_aque_mole_conc*VLNO3+NO3_1e_band_conc*VLNOB+OH_1e_aque_mole_conc+HCO3_e_mole_conc+Cl_e_conc &
    +AlO4H4_1e_aque_mole_conc+FeO4H4_1e_aque_mole_conc+NaCO3_1e_aque_mole_conc+NaSO4_1e_aque_mole_conc+KSO4_1e_aque_mole_conc+(H2PO4_1e_aque_mole_conc+CaPO4_1e_con)*VLPO4 &
    +(H2PO4_1e_band_mole_conc+CaPO4_1e_band_conc)*VLPOB
  CSTR1=AZMAX1(0.5E-03_r8*(9.0_r8*(cation_3p_aque_mole_conc+anion_3e_conc) &
    +4.0_r8*(cation_2p_aque_mole_conc+anion_2e_aque_mole_conc)+cation_1p_aque_mole_conc+anion_1e_aque_mole_conc))

  CSTR2=SQRT(CSTR1)
  FSTR2=CSTR2/(1.0_r8+CSTR2)
!
!     ACTIVITY COEFFICIENTS CALCULATED FROM ION STRENGTH
!
  A1 = AMIN1(1.0_r8,10.0_r8**(-0.509_r8*1.0_r8*FSTR2+0.20_r8*CSTR2))
  A2 = AMIN1(1.0_r8,10.0_r8**(-0.509_r8*4.0_r8*FSTR2+0.20_r8*CSTR2))
  A3 = AMIN1(1.0_r8,10.0_r8**(-0.509_r8*9.0_r8*FSTR2+0.20_r8*CSTR2))
!
!     PRECIPITATION-DISSOLUTION REACTIONS IN NON-BAND, BAND
!     CALCULATED FROM ACTIVITIES OF REACTANTS AND PRODUCTS THROUGH
!     CONVERGENCE SOLUTIONS FOR THEIR EQUILIBRIUM CONCENTRATIONS USING
!     SOLUTE FORMS CURRENTLY AT HIGHEST CONCENTRATIONS
!
!     for all reactions:
!     A*=ion activity
!     PX,PY=solute forms with greatest activity
!     R*,P*=reactant,product
!     NR*,NP*=reactant,product stoichiometry
!     SP=solubility product of PX from parameters above
!     SPX=equilibrium product concentration
!     R*X=precipitation-dissolution rate
!
  ! Prescribed [H(+)] is unchanged; impose a_H*a_OH=DPH2O using A1 for both.
  ! Record the initial correction separately from the subsequent reaction fluxes.
  OH_1e_equil_conc=DPH2O/(A1*A1*H_1p_aque_mole_conc)
  RHOH_equil=OH_1e_equil_conc-OH_1e_aque_mole_conc
  OH_1e_aque_mole_conc=OH_1e_equil_conc
  H_1p_activity             = H_1p_aque_mole_conc*A1
  OH_1e_activity            = OH_1e_aque_mole_conc*A1
  Al_3p_activity            = Al_3p_aque_mole_conc*A3
  AlOH_2p_activity          = AlOH_2p_aque_mole_conc*A2
  AlO2H2_1p_activity        = AlO2H2_1p_aque_mole_conc*A1
  AlO3H3_activity           = AlO3H3_mole_conc*A0
  AlO4H4_1e_activity        = AlO4H4_1e_aque_mole_conc*A1
  Fe_3p_activity            = Fe_3p_aque_mole_conc*A3
  FeOH_2p_activity          = FeOH_2p_aque_mole_conc*A2
  FeO2H2_p_activity         = FeO2H2_p_mole_conc*A1
  FeO3H3_activity           = FeO3H3_mole_conc*A0
  FeO4H4_1e_activity        = FeO4H4_1e_aque_mole_conc*A1
  Ca_2p_activity            = Ca_2p_aque_mole_conc*A2
  CO3_2e_activity           = CO3_2e_aque_mole_conc*A2
  HCO3_e_activity           = HCO3_e_mole_conc*A1
  H2CO3_activity            = H2CO3_aque_mole_conc*A0
  SO4_2e_activity           = SO4_2e_aque_mole_conc*A2
  H0PO4_3e_activity         = H0PO4_3e_conc*A3
  H1PO4_2e_activity         = H1PO4_2e_aque_mole_conc*A2
  H2PO4_1e_activity         = H2PO4_1e_aque_mole_conc*A1
  H3PO4_activity            = H3PO4_conc*A0
  FeHPO4_p_activity         = FeHPO4_p_conc*A1
  FeH2PO4_2p_activity       = FeH2PO4_2p_aque_mole_conc*A2
  CaPO4_1e_activity         = CaPO4_1e_con*A1
  CaHPO4_activity           = CaHPO4_conc*A0
  CaH4P2O8_1p_activity      = CaH4P2O8_1p_aque_mole_conc*A1
  MgHPO4_activity           = MgHPO4_conc*A0
  H0PO4_3e_band_activity    = H0PO4_3e_band_conc*A3
  H1PO4_2e_band_activity    = H1PO4_2e_band_mole_conc*A2
  H2PO4_1e_band_activity    = H2PO4_1e_band_mole_conc*A1
  H3PO4_band_activity       = H3PO4_band_conc*A0
  FeHPO4_1p_band_activity   = FeHPO4_1p_band_conc*A1
  FeH2PO4_2p_band_activity  = FeH2PO4_2p_band_conc*A2
  CaPO4_1e_band_activity    = CaPO4_1e_band_conc*A1
  CaHPO4_band_activity      = CaHPO4_band_conc*A0
  CaH4P2O8_1p_band_activity = CaH4P2O8_1p_band_conc*A1
  MgHPO4_band_activity      = MgHPO4_band_conc*A0
  NH4_1p_activity           = NH4_1p_aque_mole_conc*A1
  NH4_1p_band_activity      = NH4_1p_band_mole_conc*A1
  NH3_activity              = NH3_aque_mole_conc*A0
  NH3_band_activity         = NH3_aqu_band_mole_conc*A0
  Mg_2p_activity            = Mg_2p_aque_mole_conc*A2
  Na_1p_activity            = Na_1p_aque_mole_conc*A1
  K_1p_activity             = K_1p_aque_mole_conc*A1
  AlSO4_1p_activity         = AlSO4_1p_aque_mole_conc*A1
  FeSO4_1p_activity         = FeSO4_1p_aque_mole_conc*A1
  CaO2H2_activity           = CaO2H2_mole_conc*A1
  CaCO3_activity            = CaCO3_mole_conc*A0
  CaSO4_activity            = CaSO4_mole_conc*A0
  CaHCO3_1p_activity        = CaHCO3_1p_aque_mole_conc*A1
  MgOH_1p_activity          = MgOH_1p_aque_mole_conc*A1
  AMGC1                     = MgCO3_conc*A0
  MgHCO3_1p_activity        = MgHCO3_1p_aque_mole_conc*A1
  MgSO4_activity            = MgSO4_conc*A0
  NaCO3_1e_activity         = NaCO3_1e_aque_mole_conc*A1
  NaSO4_1e_activity         = NaSO4_1e_aque_mole_conc*A1
  AKAS1                     = KSO4_1e_aque_mole_conc*A1
  call PrintInfo('end '//subname)
  end subroutine IonStrengthActivity
!------------------------------------------------------------------------------------------

  subroutine PhospPrecipDissolBand

  implicit none
  real(r8) :: SP,SPX,S0,S1,P1,P2,P3,PX,PY
  integer :: NR1,NP3
!     begin_execution

  IF(VLWatMicPPB.GT.ZEROS2)THEN
!
!     ALUMINUM PHOSPHATE (VARISCITE)
!
    ! PX/PY copy an AMAX1 operand, so exact equality identifies the selected
    ! species even at trace activity; the first exact tie keeps precedence.
    PX = AMAX1(Al_3p_activity,AlOH_2p_activity,AlO2H2_1p_activity,AlO3H3_activity,AlO4H4_1e_activity)
    PY = AMAX1(H1PO4_2e_band_activity,H2PO4_1e_band_activity)
    R1 = H_1p_activity
    P3 = H_1p_activity
    IF(PY.EQ.H1PO4_2e_band_activity)THEN
      P2=H1PO4_2e_band_activity
      IF(PX.EQ.Al_3p_activity)THEN
        P1  = Al_3p_activity
        NR1 = 1
        NP3 = 0
        SP  = SHA0P1
      ELSEIF(PX.EQ.AlOH_2p_activity)THEN
        P1  = AlOH_2p_activity
        NR1 = 0
        NP3 = 0
        SP  = SPA1P1
      ELSEIF(PX.EQ.AlO2H2_1p_activity)THEN
        P1  = AlO2H2_1p_activity
        NR1 = 0
        NP3 = 1
        SP  = SHA2P1
      ELSEIF(PX.EQ.AlO3H3_activity)THEN
        P1  = AlO3H3_activity
        NR1 = 0
        NP3 = 2
        SP  = SHA3P1
      ELSEIF(PX.EQ.AlO4H4_1e_activity)THEN
        P1  = AlO4H4_1e_activity
        NR1 = 0
        NP3 = 3
        SP  = SHA4P1
      ENDIF
    ELSE
      P2=H2PO4_1e_band_activity
      IF(PX.EQ.Al_3p_activity)THEN
        P1  = Al_3p_activity
        NR1 = 2
        NP3 = 0
        SP  = SHA0P2
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
    RHA0B1                     = 0._r8
    RHA1B1                     = 0._r8
    RHA2B1                     = 0._r8
    RHA3B1                     = 0._r8
    RHA4B1                     = 0._r8
    RHA0B2                     = 0._r8
    RHA1B2                     = 0._r8
    RHA2B2                     = 0._r8
    RHA3B2                     = 0._r8
    RHA4B2                     = 0._r8
    R1                         = AMAX1(ZERO,R1)
    P1                         = AMAX1(ZERO,P1)
    P2                         = AMAX1(ZERO,P2)
    P3                         = AMAX1(ZERO,P3)
    SPX                        = SP*R1**NR1/P3**NP3
    S0                         = P1+P2
    S1                         = MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    H2PO4_1e_AlPO4_dissolB_flx = AMAX1(-PrecpB_AlPO4_mole_conc,TPDX*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))
    IF(PY.EQ.H1PO4_2e_band_activity)THEN
      IF(PX.EQ.Al_3p_activity)THEN
        RHA0B1=H2PO4_1e_AlPO4_dissolB_flx
      ELSEIF(PX.EQ.AlOH_2p_activity)THEN
        RHA1B1=H2PO4_1e_AlPO4_dissolB_flx
      ELSEIF(PX.EQ.AlO2H2_1p_activity)THEN
        RHA2B1=H2PO4_1e_AlPO4_dissolB_flx
      ELSEIF(PX.EQ.AlO3H3_activity)THEN
        RHA3B1=H2PO4_1e_AlPO4_dissolB_flx
      ELSEIF(PX.EQ.AlO4H4_1e_activity)THEN
        RHA4B1=H2PO4_1e_AlPO4_dissolB_flx
      ENDIF
    ELSE
      IF(PX.EQ.Al_3p_activity)THEN
        RHA0B2=H2PO4_1e_AlPO4_dissolB_flx
      ELSEIF(PX.EQ.AlOH_2p_activity)THEN
        RHA1B2=H2PO4_1e_AlPO4_dissolB_flx
      ELSEIF(PX.EQ.AlO2H2_1p_activity)THEN
        RHA2B2=H2PO4_1e_AlPO4_dissolB_flx
      ELSEIF(PX.EQ.AlO3H3_activity)THEN
        RHA3B2=H2PO4_1e_AlPO4_dissolB_flx
      ELSEIF(PX.EQ.AlO4H4_1e_activity)THEN
        RHA4B2=H2PO4_1e_AlPO4_dissolB_flx
      ENDIF
    ENDIF
!
!     IRON PHOSPHATE (STRENGITE)
!
    PX=AMAX1(Fe_3p_activity,FeOH_2p_activity,FeO2H2_p_activity,FeO3H3_activity,FeO4H4_1e_activity)
    PY=AMAX1(H1PO4_2e_band_activity,H2PO4_1e_band_activity)
    R1=H_1p_activity
    P3=H_1p_activity
    IF(PY.EQ.H1PO4_2e_band_activity)THEN
      P2=H1PO4_2e_band_activity
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
      P2=H2PO4_1e_band_activity
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
    RHF0B1 = 0._r8
    RHF1B1 = 0._r8
    RHF2B1 = 0._r8
    RHF3B1 = 0._r8
    RHF4B1 = 0._r8
    RHF0B2 = 0._r8
    RHF1B2 = 0._r8
    RHF2B2 = 0._r8
    RHF3B2 = 0._r8
    RHF4B2 = 0._r8
    R1     = AMAX1(ZERO,R1)
    P1     = AMAX1(ZERO,P1)
    P2     = AMAX1(ZERO,P2)
    P3     = AMAX1(ZERO,P3)
    SPX    = SP*R1**NR1/P3**NP3
    S0     = P1+P2
    S1     = MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    H2PO4_1e_FePO4_dissolB_flx=AMAX1(-PrecpB_FePO4_mole_conc,TPDX*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))
    IF(PY.EQ.H1PO4_2e_band_activity)THEN
      IF(PX.EQ.Fe_3p_activity)THEN
        RHF0B1=H2PO4_1e_FePO4_dissolB_flx
      ELSEIF(PX.EQ.FeOH_2p_activity)THEN
        RHF1B1=H2PO4_1e_FePO4_dissolB_flx
      ELSEIF(PX.EQ.FeO2H2_p_activity)THEN
        RHF2B1=H2PO4_1e_FePO4_dissolB_flx
      ELSEIF(PX.EQ.FeO3H3_activity)THEN
        RHF3B1=H2PO4_1e_FePO4_dissolB_flx
      ELSEIF(PX.EQ.FeO4H4_1e_activity)THEN
        RHF4B1=H2PO4_1e_FePO4_dissolB_flx
      ENDIF
    ELSE
      IF(PX.EQ.Fe_3p_activity)THEN
        RHF0B2=H2PO4_1e_FePO4_dissolB_flx
      ELSEIF(PX.EQ.FeOH_2p_activity)THEN
        RHF1B2=H2PO4_1e_FePO4_dissolB_flx
      ELSEIF(PX.EQ.FeO2H2_p_activity)THEN
        RHF2B2=H2PO4_1e_FePO4_dissolB_flx
      ELSEIF(PX.EQ.FeO3H3_activity)THEN
        RHF3B2=H2PO4_1e_FePO4_dissolB_flx
      ELSEIF(PX.EQ.FeO4H4_1e_activity)THEN
        RHF4B2=H2PO4_1e_FePO4_dissolB_flx
      ENDIF
    ENDIF
!
!     DICALCIUM PHOSPHATE
!
    PX=AMAX1(H1PO4_2e_band_activity,H2PO4_1e_band_activity)
    R1=H_1p_activity
    P1=Ca_2p_activity
    IF(PX.EQ.H1PO4_2e_band_activity)THEN
      P2=H1PO4_2e_band_activity
      NR1=0
      SP=SPCAD
    ELSEIF(PX.EQ.H2PO4_1e_band_activity)THEN
      P2=H2PO4_1e_band_activity
      NR1=1
      SP=SHCAD2
    ENDIF
    RPCDB1=0._r8
    RHCDB2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1
    S0=P1+P2
    S1=MAX(0._r8,(P1-P2)**2+4._r8*SPX)
    H2PO4_1e_CaHPO4_dissolB_flx=AMAX1(-PrecpB_CaHPO4_mole_conc,TPDX*(4._r8*((P1*P2-SPX)/(S0+SQRT(S1)))))
    IF(PX.EQ.H1PO4_2e_band_activity)THEN
      RPCDB1=H2PO4_1e_CaHPO4_dissolB_flx
    ELSEIF(PX.EQ.H2PO4_1e_band_activity)THEN
      RHCDB2=H2PO4_1e_CaHPO4_dissolB_flx
    ENDIF
!
!     HYDROXYAPATITE
!
    PX=AMAX1(H1PO4_2e_band_activity,H2PO4_1e_band_activity)
    R1=H_1p_activity
    P1=Ca_2p_activity
    IF(PX.EQ.H1PO4_2e_band_activity)THEN
      P2=H1PO4_2e_band_activity
      NR1=4
      SP=SHCAH1
    ELSEIF(PX.EQ.H2PO4_1e_band_activity)THEN
      P2=H2PO4_1e_band_activity
      NR1=7
      SP=SHCAH2
    ENDIF
    RHCHB1=0._r8
    RHCHB2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=(SP*R1**NR1/P1**5._r8)**0.333_r8
    H2PO4_1e_apatite_dissolB_flx=AMAX1(-PrecpB_Ca5P3O12O3H3_mole_conc,TPDX*(P2-SPX))
    IF(PX.EQ.H1PO4_2e_band_activity)THEN
      RHCHB1=H2PO4_1e_apatite_dissolB_flx
    ELSEIF(PX.EQ.H2PO4_1e_band_activity)THEN
      RHCHB2=H2PO4_1e_apatite_dissolB_flx
    ENDIF
!
!     MONOCALCIUM PHOSPHATE
!
    P1=Ca_2p_activity
    P2=H2PO4_1e_band_activity
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    ! Ca(H2PO4)2(s) <-> Ca(2+)+2H2PO4(-): K=a_Ca*a_H2PO4**2.
    ! Use the same fixed-calcium relaxation and formula-unit flux as non-band soil.
    SPX=SQRT(SPCAM/P1)
    H2PO4_1e_CaH4P2O8_dissolB_flx=AMAX1(-PrecpB_CaH4P2O8_mole_conc,TPDX*(P2-SPX))
  ELSE
    H2PO4_1e_AlPO4_dissolB_flx=0._r8
    H2PO4_1e_FePO4_dissolB_flx=0._r8
    H2PO4_1e_CaHPO4_dissolB_flx=0._r8
    H2PO4_1e_apatite_dissolB_flx=0._r8
    H2PO4_1e_CaH4P2O8_dissolB_flx=0._r8
    RHA0B1=0._r8
    RHA1B1=0._r8
    RHA2B1=0._r8
    RHA3B1=0._r8
    RHA4B1=0._r8
    RHA0B2=0._r8
    RHA1B2=0._r8
    RHA2B2=0._r8
    RHA3B2=0._r8
    RHA4B2=0._r8
    RHF0B1=0._r8
    RHF1B1=0._r8
    RHF2B1=0._r8
    RHF3B1=0._r8
    RHF4B1=0._r8
    RHF0B2=0._r8
    RHF1B2=0._r8
    RHF2B2=0._r8
    RHF3B2=0._r8
    RHF4B2=0._r8
    RPCDB1=0._r8
    RHCDB2=0._r8
    RHCHB1=0._r8
    RHCHB2=0._r8
  ENDIF
  end subroutine PhospPrecipDissolBand
!------------------------------------------------------------------------------------------

  subroutine PhospAnionExchNoBand
  implicit none
  real(r8) :: SPH1P,SPH2P
!     begin_execution
!     PHOSPHORUS ANION EXCHANGE IN NON-BAND SOIL ZONE
!     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
!     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
!     EXCHANGE SITES
!
!     SoilMicPMassLayer=soil mass
!     VLWatMicPM=soil water volume
!     XAEC=anion exchange capacity
!     VLWatMicPPO=soil water volume in non-band
!     TADAX=adsorption rate constant
!     RXOH2,RXOH1=OH2,OH exchange with R-OH2,R-OH in non-band
!     SXOH2,SXOH1=equilibrium constant for OH2,OH exchange with R-OH2,R-OH
!
  IF(VLWatMicPM.GT.ZEROS2)THEN
    VLWatMicPBK=AMIN1(1.0,SoilMicPMassLayer/VLWatMicPM)
  ELSE
    VLWatMicPBK=1._r8
  ENDIF
  IF(VLWatMicPPO.GT.ZEROS2.AND.XAEC.GT.ZEROS)THEN
    RXOH2=TADAX*(XROH1_mole_conc*H_1p_activity-SXOH2*XROH2_mole_conc)/(XROH1_mole_conc+SXOH2)*VLWatMicPBK
    RXOH1=TADAX*(XOH_mole_conc*H_1p_activity-SXOH1*XROH1_mole_conc)/(XOH_mole_conc+SXOH1)*VLWatMicPBK
!
!     H2PO4 EXCHANGE IN NON-BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG H2PO4-, H+, OH-, X-OH
!     AND X-H2PO4
!
!     SPH2P,SXH2P=equilibrium constant for H2PO4 exchange with R-OH2,R-OH
!     H2PO4_1e_to_XH2PO4_ROH2_flx,H2PO4_to_XH2PO4_ROH_flx=H2PO4 exchange with R-OH2,R-OH in non-band
!
    ! R-H2PO4 + H2O <=> R-OH2(+) + H2PO4(-).
    ! Combine neutral-site exchange (SXH2P), water dissociation (DPH2O),
    ! and R-OH2(+) <=> R-OH + H(+) (SXOH2) to close the equilibrium cycle.
    SPH2P = SXH2P*DPH2O/SXOH2
    H2PO4_1e_to_XH2PO4_ROH2_flx = TADAX*(XROH2_mole_conc*H2PO4_1e_activity-SPH2P*XH2PO4_mole_conc)/(XROH2_mole_conc+SPH2P)*VLWatMicPBK
    H2PO4_to_XH2PO4_ROH_flx     = TADAX*(XROH1_mole_conc*H2PO4_1e_activity-SXH2P*XH2PO4_mole_conc*OH_1e_activity)/(XROH1_mole_conc+SXH2P)*VLWatMicPBK
!
!     HPO4 EXCHANGE IN NON-BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG HPO4--, H+, OH-, X-OH
!     AND X-HPO4
!
!     SPH1P=equilibrium constant for HPO4 exchange with R-OH
!     H1PO4_to_XHPO4_ROH_flx=HPO4 exchange with R-OH in non-band
!
    SPH1P=SXH1P*DPH2O/DPH2P
    H1PO4_to_XHPO4_ROH_flx=TADAX*(XROH1_mole_conc*H1PO4_2e_activity-SPH1P*XHPO4_mole_conc)/(XROH1_mole_conc+SPH1P)*VLWatMicPBK
  ELSE
    RXOH2                       = 0._r8
    RXOH1                       = 0._r8
    H2PO4_1e_to_XH2PO4_ROH2_flx = 0._r8
    H2PO4_to_XH2PO4_ROH_flx     = 0._r8
    H1PO4_to_XHPO4_ROH_flx      = 0._r8
  ENDIF
  end subroutine PhospAnionExchNoBand
!------------------------------------------------------------------------------------------

  subroutine PhospAnionExchBand
  implicit none
  real(r8) :: SPH1P,SPH2P
!     begin_execution
!     PHOSPHORUS ANION EXCHANGE IN BAND SOIL ZONE
!     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
!     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
!     EXCHANGE SITES
!
  IF(VLWatMicPPB.GT.ZEROS2.AND.XAEC.GT.ZEROS)THEN
!
!     RXO2B,RXO1B=OH2,OH exchange with R-OH2,R-OH in band
!     SXOH2,SXOH1=equilibrium constant for OH2,OH exchange with R-OH2,R-OH
!
    RXO2B=TADAX*(XROH_band_mole_conc*H_1p_activity-SXOH2*XROH2_band_mole_conc)/(XROH_band_mole_conc+SXOH2)*VLWatMicPBK
    RXO1B=TADAX*(XROH1_band_mole_conc*H_1p_activity-SXOH1*XROH_band_mole_conc)/(XROH1_band_mole_conc+SXOH1)*VLWatMicPBK
!
!     H2PO4 EXCHANGE IN BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG H2PO4-, H+, OH-, X-OH
!     AND X-H2PO4
!
!     SPH2P,SXH2P=equilibrium constant for H2PO4 exchange with R-OH2,R-OH
!     H2PO4_1e_to_XH2PO4_ROH2_Bflx,H2PO4_to_XH2PO4_ROH_Bflx=H2PO4 exchange with R-OH2,R-OH in band
!
    ! R-H2PO4 + H2O <=> R-OH2(+) + H2PO4(-).
    ! Combine neutral-site exchange (SXH2P), water dissociation (DPH2O),
    ! and R-OH2(+) <=> R-OH + H(+) (SXOH2) to close the equilibrium cycle.
    SPH2P = SXH2P*DPH2O/SXOH2
    H2PO4_1e_to_XH2PO4_ROH2_Bflx=TADAX*(XROH2_band_mole_conc*H2PO4_1e_band_activity-SPH2P*XH2PO4_band_mole_conc)/(XROH2_band_mole_conc+SPH2P)*VLWatMicPBK
    H2PO4_to_XH2PO4_ROH_Bflx=TADAX*(XROH_band_mole_conc*H2PO4_1e_band_activity-SXH2P*XH2PO4_band_mole_conc*OH_1e_activity)/(XROH_band_mole_conc+SXH2P)*VLWatMicPBK
!
!     HPO4 EXCHANGE IN BAND SOIL ZONE FROM CONVERGENCE
!     SOLUTION FOR EQUILIBRIUM AMONG HPO4--, H+, OH-, X-OH
!     AND X-HPO4
!
!     SPH1P=equilibrium constant for HPO4 exchange with R-OH
!     RXH1B=HPO4 exchange with R-OH in band
!
    SPH1P=SXH1P*DPH2O/DPH2P
    RXH1B=TADAX*(XROH_band_mole_conc*H1PO4_2e_band_activity-SPH1P*XHPO4_band_mole_conc)/(XROH_band_mole_conc+SPH1P)*VLWatMicPBK
!     WRITE(*,2226)'RXH1B',I,J,L,M,RXH1B,XROH_band_mole_conc,XROH2_band_mole_conc,H1PO4_2e_band_mole_conc
!    2,SPH1P,XHPO4_band_mole_conc,H2PO4_to_XH2PO4_ROH_Bflx,H2PO4_1e_band_mole_conc,SXH2P,XH2PO4_band_mole_conc,OH_1e_aque_mole_conc,H_1p_aque_mole_conc,ROH
!2226  FORMAT(A8,4I4,20E12.4)
  ELSE
    RXO2B=0._r8
    RXO1B=0._r8
    H2PO4_1e_to_XH2PO4_ROH2_Bflx=0._r8
    H2PO4_to_XH2PO4_ROH_Bflx=0._r8
    RXH1B=0._r8
  ENDIF
  end subroutine PhospAnionExchBand
!------------------------------------------------------------------------------------------

  subroutine CationExchange
  implicit none
  real(r8) :: FX,XALQ,XCAQ,XCAXl,XKAQ,XCAX,XFEQ
  real(r8) :: XHYQ,XMGQ,XN4Q,XNBQ,XTLQ,XNAQ
!     begin_execution

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
!     Ca-NH4,Ca-H,Ca-Al,Ca-Mg,Ca-Na,Ca-K
!     X*Q=equilibrium exchangeable concentrations
!     XTLQ=total equilibrium exchangeable concentration
!
    AALX=Al_3p_activity**0.333_r8
    AFEX=Fe_3p_activity**0.333_r8
    ACAX=Ca_2p_activity**0.500_r8
    AMGX=Mg_2p_activity**0.500_r8
    XCAX=CEC_conc/(1.0+GKC4*NH4_1p_activity/ACAX*VLNH4 &
      +GKC4*NH4_1p_band_activity/ACAX*VLNHB &
      +GKCH*H_1p_activity/ACAX+GKCA*AALX/ACAX &
      +GKCA*AFEX/ACAX+GKCM*AMGX/ACAX &
      +GKCN*Na_1p_activity/ACAX+GKCK*K_1p_activity/ACAX)
    XN4Q=XCAX*NH4_1p_activity*GKC4
    XNBQ=XCAX*NH4_1p_band_activity*GKC4
    XHYQ=XCAX*H_1p_activity*GKCH
    XALQ=XCAX*AALX*GKCA
    XFEQ=XCAX*AFEX*GKCA
    XCAQ=XCAX*ACAX
    XMGQ=XCAX*AMGX*GKCM
    XNAQ=XCAX*Na_1p_activity*GKCN
    XKAQ=XCAX*K_1p_activity*GKCK
    XTLQ=XN4Q*VLNH4+XNBQ*VLNHB+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ
    IF(XTLQ.GT.ZERO)THEN
      FX=CEC_conc/XTLQ
    ELSE
      FX=0._r8
    ENDIF
    XN4Q=FX*XN4Q
    XNBQ=FX*XNBQ
    XHYQ=FX*XHYQ
    XALQ=FX*XALQ/3.0_r8
    XFEQ=FX*XFEQ/3.0_r8
    XCAQ=FX*XCAQ/2.0_r8
    XMGQ=FX*XMGQ/2.0_r8
    XNAQ=FX*XNAQ
    XKAQ=FX*XKAQ
!
!     NH4 EXCHANGE IN NON-BAND AND BAND SOIL ZONES
!
!     RXN4,RXNB=NH4 adsorption in non-band,band
!     TADCX=adsorption rate constant
!
    RXN4=TADCX*AMAX1(AMIN1((XN4Q-XNH4_mole_conc)*NH4_1p_activity/XN4Q,NH4_1p_aque_mole_conc),-XNH4_mole_conc)
    RXNB=TADCX*AMAX1(AMIN1((XNBQ-XNH4_band_mole_conc)*NH4_1p_band_activity/XNBQ,NH4_1p_band_mole_conc),-XNH4_band_mole_conc)
!
!     H,AL,FE,CA,MG,NA,K EXCHANGE
!
!     RX*=ion adsorption
!
    ! TODO: Desorption below can overdraw exchangeable donors for small targets.
    ! NH4 above has a donor limit, but reactions sharing aqueous donors are not
    ! limited jointly. Retain these equations pending a geochemistry solver that
    ! handles all reactions concurrently, preserves nonnegative pools, and
    ! conserves elemental totals (also applies to initialization chemistry).
    RXHY=TADCX*AMIN1((XHYQ-XHY1)*H_1p_activity/XHYQ,H_1p_aque_mole_conc)
    RXAL=TADCX*AMIN1((XALQ-XAl_conc)*AALX/XALQ,Al_3p_aque_mole_conc)
    RXFE=TADCX*AMIN1((XFEQ-XFe_conc)*AFEX/XFEQ,Fe_3p_aque_mole_conc)
    RXCA=TADCX*AMIN1((XCAQ-XCa_conc)*ACAX/XCAQ,Ca_2p_aque_mole_conc)
    RXMG=TADCX*AMIN1((XMGQ-XMg_conc)*AMGX/XMGQ,Mg_2p_aque_mole_conc)
    RXNA=TADCX*AMIN1((XNAQ-XNa_conc)*Na_1p_activity/XNAQ,Na_1p_aque_mole_conc)
    RXKA=TADCX*AMIN1((XKAQ-XK_conc)*K_1p_activity/XKAQ,K_1p_aque_mole_conc)
!     IF(I.EQ.256.AND.L.EQ.1)THEN
!     WRITE(*,1112)'RXAL',I,J,L,M,RXAL,TADCX,XALQ,XAl_conc,Al_3p_aque_mole_conc
!    2,Al_3p_activity,AALX,CEC_conc,XCAX,GKCA,FX
!     ENDIF
  ELSE
    RXN4=0._r8
    RXNB=0._r8
    RXHY=0._r8
    RXAL=0._r8
    RXFE=0._r8
    RXCA=0._r8
    RXMG=0._r8
    RXNA=0._r8
    RXKA=0._r8
  ENDIF
  end subroutine CationExchange
!------------------------------------------------------------------------------------------

  subroutine SoluteDissociation
  implicit none
  real(r8) :: S0,S1
!     begin_execution
!     for all reactions:
!     S0,S1=equilibrium constant,equilibrium solute concentration**2
!     TSLX=dissociation rate constant
!
!     DISSOCIATION OF CARBOXYL RADICALS
!     AND ADSORPTION OF AL AND FE OH2
!
!     RXHC=COOH-COO+H dissociation
!     RXALO2,RXFLO2=Al(OH2)2+,FeOH2+ adsorption
! COOH <-> COO(-) + H(+)
  S0     = H_1p_activity+XCOO+DPCOH
  S1     = MAX(0._r8,(H_1p_activity-XCOO)**2+2._r8*DPCOH*(H_1p_activity+XCOO)+DPCOH**2+4._r8*DPCOH*XHC1)
  RXHC   = TADCX*(4._r8*((H_1p_activity*XCOO-DPCOH*XHC1)/(S0+SQRT(S1))))

  S0     = AlO2H2_1p_activity+XCOO+DPALO
  S1     = MAX(0._r8,(AlO2H2_1p_activity-XCOO)**2+2._r8*DPALO*(AlO2H2_1p_activity+XCOO)+DPALO**2+4._r8*DPALO*XAlO2H2_conc)
  RXALO2 = TADAX*(4._r8*((AlO2H2_1p_activity*XCOO-DPALO*XAlO2H2_conc)/(S0+SQRT(S1))))
  
  S0     = FeO2H2_p_activity+XCOO+DPFEO
  S1     = MAX(0._r8,(FeO2H2_p_activity-XCOO)**2+2._r8*DPFEO*(FeO2H2_p_activity+XCOO)+DPFEO**2+4._r8*DPFEO*XFeO2H2_conc)
  RXFEO2 = TADAX*(4._r8*((FeO2H2_p_activity*XCOO-DPFEO*XFeO2H2_conc)/(S0+SQRT(S1))))

!
!     RNH4,RNHB-NH4-NH3+H dissociation in non-band,band
!     DPN4=NH4 dissociation constant
! NH4(+) <-> NH3 + H(+), RNH4>0 for increasing NH4(+)
  IF(VLWatMicPNH.GT.ZEROS2)THEN
    RNH4=TSLX*(H_1p_activity*NH3_activity-DPN4*NH4_1p_activity)/(DPN4+H_1p_activity)
  ELSE
    RNH4=0._r8
  ENDIF

  IF(VLWatMicPNB.GT.ZEROS2)THEN
    RNHB=TSLX*(H_1p_activity*NH3_band_activity-DPN4*NH4_1p_band_activity)/(DPN4+H_1p_activity)
  ELSE
    RNHB=0._r8
  ENDIF

!
! RCO2Q=CO2-HCO3+H dissociation, > 0 for increasing H2CO3
! H2CO3 <-> HCO3(-) + H(+)

  S0    = H_1p_activity+HCO3_e_activity+DPCO2
  S1    = MAX(0._r8,(H_1p_activity-HCO3_e_activity)**2+2._r8*DPCO2*(H_1p_activity+HCO3_e_activity)+DPCO2**2+4._r8*DPCO2*H2CO3_activity)
  RCO2Q = TSLX*(4._r8*((H_1p_activity*HCO3_e_activity-DPCO2*H2CO3_activity)/(S0+SQRT(S1))))
!
!     RHCO3=HCO3-CO3+H dissociation, > 0 for increasing HCO3
! HCO3(-) <-> CO3(--)+H(+)
  S0    = H_1p_activity+CO3_2e_activity+DPHCO
  S1    = MAX(0._r8,(H_1p_activity-CO3_2e_activity)**2+2._r8*DPHCO*(H_1p_activity+CO3_2e_activity)+DPHCO**2+4._r8*DPHCO*HCO3_e_activity)
  RHCO3 = TSLX*(4._r8*((H_1p_activity*CO3_2e_activity-DPHCO*HCO3_e_activity)/(S0+SQRT(S1))))
!
!     RALO1=ALOH(++) <-> AL(+++)+OH(-) dissociation
!
  RALO1=TSLX*(Al_3p_activity*OH_1e_activity-DPAL1*AlOH_2p_activity)/(OH_1e_activity+DPAL1)
!
!     RALO2=ALOH2-ALOH+OH dissociation
! Al(OH)2(+) <-> Al(OH)(++)+OH(-)
  RALO2=TSLX*(AlOH_2p_activity*OH_1e_activity-DPAL2*AlO2H2_1p_activity)/(OH_1e_activity+DPAL2)
!
!     RALO3=ALOH3 <-> ALOH2+OH dissociation
! Al(OH)3 <-> Al(OH)2(+)+OH(-)
  RALO3=TSLX*(AlO2H2_1p_activity*OH_1e_activity-DPAL3*AlO3H3_activity)/(OH_1e_activity+DPAL3)
!
!     RALO4=ALOH4 <-> ALOH3+OH dissociation
! Al(OH)4(-) <-> Al(OH)3+OH(-)
  RALO4=TSLX*(AlO3H3_activity*OH_1e_activity-DPAL4*AlO4H4_1e_activity)/(OH_1e_activity+DPAL4)
!
!     RALS=ALSO4 <-> AL+SO4 dissociation
! AlSO4(+) <-> Al(3+) + SO4(2-)
  S0=Al_3p_activity+SO4_2e_activity+DPALS
  S1=MAX(0._r8,(Al_3p_activity-SO4_2e_activity)**2+2._r8*DPALS*(Al_3p_activity+SO4_2e_activity)+DPALS**2+4._r8*DPALS*AlSO4_1p_activity)
  RALS=TSLX*(4._r8*((Al_3p_activity*SO4_2e_activity-DPALS*AlSO4_1p_activity)/(S0+SQRT(S1))))
!
!     RFEO1=FEOH <-> FE+OH dissociation
! Fe(OH)(++) <-> Fe(+++)+OH(-)
  RFEO1=TSLX*(Fe_3p_activity*OH_1e_activity-DPFE1*FeOH_2p_activity)/(OH_1e_activity+DPFE1)
!
!     RFEO2=FEOH2 <-> FEOH+OH dissociation
! Fe(OH)2(+) <-> Fe(OH)(++)+OH(-)
  RFEO2=TSLX*(FeOH_2p_activity*OH_1e_activity-DPFE2*FeO2H2_p_activity)/(OH_1e_activity+DPFE2)
!
!     RFEO3=FEOH3 <-> FEOH2+OH dissociation
! Fe(OH)3 <-> Fe(OH)2(+)+OH(-)
  RFEO3=TSLX*(FeO2H2_p_activity*OH_1e_activity-DPFE3*FeO3H3_activity)/(OH_1e_activity+DPFE3)
!
!     RFE04=FEOH4 <-> FEOH3+OH dissociation
! Fe(OH)4(-) <-> Fe(OH)3 + OH(-)
  RFEO4=TSLX*(FeO3H3_activity*OH_1e_activity-DPFE4*FeO4H4_1e_activity)/(OH_1e_activity+DPFE4)
!
!     RFES <-> FE+SO4 dissociation
! FeSO4(+) <-> Fe(+++)+SO4(--)
  S0=Fe_3p_activity+SO4_2e_activity+DPFES
  S1=MAX(0._r8,(Fe_3p_activity-SO4_2e_activity)**2+2._r8*DPFES*(Fe_3p_activity+SO4_2e_activity)+DPFES**2+4._r8*DPFES*FeSO4_1p_activity)
  RFES=TSLX*(4._r8*((Fe_3p_activity*SO4_2e_activity-DPFES*FeSO4_1p_activity)/(S0+SQRT(S1))))
!
!     RCAO=CAOH <-> CA+OH dissociation
! Ca(OH)(+) <-> Ca(++)+OH(-)
  RCAO=TSLX*(Ca_2p_activity*OH_1e_activity-DPCAO*CaO2H2_activity)/(OH_1e_activity+DPCAO)
!
!     RCAC=CACO3 <-> CA+CO3 dissociation
! CaCO3 <-> Ca(++)+CO3(--)
  S0=Ca_2p_activity+CO3_2e_activity+DPCAC
  S1=MAX(0._r8,(Ca_2p_activity-CO3_2e_activity)**2+2._r8*DPCAC*(Ca_2p_activity+CO3_2e_activity)+DPCAC**2+4._r8*DPCAC*CaCO3_activity)
  RCAC=TSLX*(4._r8*((Ca_2p_activity*CO3_2e_activity-DPCAC*CaCO3_activity)/(S0+SQRT(S1))))
!
!     RCAH=CAHCO3<->CA+HCO3 dissociation
! CaHCO3(+) <-> Ca(++) + HCO3(-)
  S0=Ca_2p_activity+HCO3_e_activity+DPCAH
  S1=MAX(0._r8,(Ca_2p_activity-HCO3_e_activity)**2+2._r8*DPCAH*(Ca_2p_activity+HCO3_e_activity)+DPCAH**2+4._r8*DPCAH*CaHCO3_1p_activity)
  RCAH=TSLX*(4._r8*((Ca_2p_activity*HCO3_e_activity-DPCAH*CaHCO3_1p_activity)/(S0+SQRT(S1))))
!
!     RCAS=CASO4<->CA+SO4 dissociation
! CaSO4 <-> Ca(++) + SO4(--)
  S0=Ca_2p_activity+SO4_2e_activity+DPCAS
  S1=MAX(0._r8,(Ca_2p_activity-SO4_2e_activity)**2+2._r8*DPCAS*(Ca_2p_activity+SO4_2e_activity)+DPCAS**2+4._r8*DPCAS*CaSO4_activity)
  RCAS=TSLX*(4._r8*((Ca_2p_activity*SO4_2e_activity-DPCAS*CaSO4_activity)/(S0+SQRT(S1))))
!
!     RMGO=MGOH<->MG+OH dissociation
! MgOH(+) <-> Mg(++)+OH(-)
  RMGO=TSLX*(Mg_2p_activity*OH_1e_activity-DPMGO*MgOH_1p_activity)/(OH_1e_activity+DPMGO)
!
!     RMGC=MGCO3<->MG+CO3 dissociation
! MgCO3 <-> Mg(++)+CO3(--)
  S0=Mg_2p_activity+CO3_2e_activity+DPMGC
  S1=MAX(0._r8,(Mg_2p_activity-CO3_2e_activity)**2+2._r8*DPMGC*(Mg_2p_activity+CO3_2e_activity)+DPMGC**2+4._r8*DPMGC*AMGC1)
  RMGC=TSLX*(4._r8*((Mg_2p_activity*CO3_2e_activity-DPMGC*AMGC1)/(S0+SQRT(S1))))
!
!     RMGH=MGHCO3<->MG+HCO3 dissociation
! MgHCO3(+) <-> Mg(++) + HCO3(-)
  S0=Mg_2p_activity+HCO3_e_activity+DPMGH
  S1=MAX(0._r8,(Mg_2p_activity-HCO3_e_activity)**2+2._r8*DPMGH*(Mg_2p_activity+HCO3_e_activity)+DPMGH**2+4._r8*DPMGH*MgHCO3_1p_activity)
  RMGH=TSLX*(4._r8*((Mg_2p_activity*HCO3_e_activity-DPMGH*MgHCO3_1p_activity)/(S0+SQRT(S1))))
!
!     RMGS=MGSO4<->MG+SO4 dissociation
! MgSO4 <-> Mg(++) + SO4(--)
  S0=Mg_2p_activity+SO4_2e_activity+DPMGS
  S1=MAX(0._r8,(Mg_2p_activity-SO4_2e_activity)**2+2._r8*DPMGS*(Mg_2p_activity+SO4_2e_activity)+DPMGS**2+4._r8*DPMGS*MgSO4_activity)
  RMGS=TSLX*(4._r8*((Mg_2p_activity*SO4_2e_activity-DPMGS*MgSO4_activity)/(S0+SQRT(S1))))
!
!     RNAC=NACO3<->NA+CO3 dissociation
! NaCO3(-) <-> Na(+) +CO3(--)
  S0=Na_1p_activity+CO3_2e_activity+DPNAC
  S1=MAX(0._r8,(Na_1p_activity-CO3_2e_activity)**2+2._r8*DPNAC*(Na_1p_activity+CO3_2e_activity)+DPNAC**2+4._r8*DPNAC*NaCO3_1e_activity)
  RNAC=TSLX*(4._r8*((Na_1p_activity*CO3_2e_activity-DPNAC*NaCO3_1e_activity)/(S0+SQRT(S1))))
!
!     RNAS=NASO4<->NA+SO4 dissociation
! NaSO4(-) <-> Na(+)+SO4(--)
  S0=Na_1p_activity+SO4_2e_activity+DPNAS
  S1=MAX(0._r8,(Na_1p_activity-SO4_2e_activity)**2+2._r8*DPNAS*(Na_1p_activity+SO4_2e_activity)+DPNAS**2+4._r8*DPNAS*NaSO4_1e_activity)
  RNAS=TSLX*(4._r8*((Na_1p_activity*SO4_2e_activity-DPNAS*NaSO4_1e_activity)/(S0+SQRT(S1))))
!
!     RKAS=KSO4<->K+SO4 dissociation, RKAS>0 for increasing KSO4
! KSO4(-) <-> K(+)+SO4(--)
  S0=K_1p_activity+SO4_2e_activity+DPKAS
  S1=MAX(0._r8,(K_1p_activity-SO4_2e_activity)**2+2._r8*DPKAS*(K_1p_activity+SO4_2e_activity)+DPKAS**2+4._r8*DPKAS*AKAS1)
  RKAS=TSLX*(4._r8*((K_1p_activity*SO4_2e_activity-DPKAS*AKAS1)/(S0+SQRT(S1))))
!
!     PHOSPHORUS IN NON-BAND SOIL ZONE
!
  IF(VLWatMicPPO.GT.ZEROS2)THEN
!
!     RH1P=HPO4<->H+PO4 dissociation in non-band
! HPO4(--) <-> H(+)+PO4(---)
    RH1P=TSLX*(H0PO4_3e_activity*H_1p_activity-DPH1P*H1PO4_2e_activity)/(DPH1P+H_1p_activity)
!
!     H2PO4_e_to_HPO4_2e_flx=H2PO4 <-> H+HPO4 dissociation in non-band
! H2PO4(-) <-> H(+) + HPO4(--)
    H2PO4_e_to_HPO4_2e_flx=TSLX*(H1PO4_2e_activity*H_1p_activity-DPH2P*H2PO4_1e_activity)/(DPH2P+H_1p_activity)
!
!     RH3P=H3PO4 <-> H+H2PO4 dissociation in non-band
! H3PO4 <-> H(+)+H2PO4(-)
    RH3P=TSLX*(H2PO4_1e_activity*H_1p_activity-DPH3P*H3PO4_activity)/(DPH3P+H_1p_activity)
!
!     RF1P=FEHPO4 <-> FE+HPO4 dissociation in non-band
! FeHPO4(+) <-> Fe(+++)+HPO4(--)
    S0=Fe_3p_activity+H1PO4_2e_activity+DPF1P
    S1=MAX(0._r8,(Fe_3p_activity-H1PO4_2e_activity)**2+2._r8*DPF1P*(Fe_3p_activity+H1PO4_2e_activity)+DPF1P**2+4._r8*DPF1P*FeHPO4_p_activity)
    RF1P=TSLX*(4._r8*((Fe_3p_activity*H1PO4_2e_activity-DPF1P*FeHPO4_p_activity)/(S0+SQRT(S1))))
!
!     RF2P=FEH2PO4 <-> FE+H2PO4 dissociation in non-band
! FeH2PO4(++) <-> Fe(+++)+H2PO4(-)
    S0=Fe_3p_activity+H2PO4_1e_activity+DPF2P
    S1=MAX(0._r8,(Fe_3p_activity-H2PO4_1e_activity)**2+2._r8*DPF2P*(Fe_3p_activity+H2PO4_1e_activity)+DPF2P**2+4._r8*DPF2P*FeH2PO4_2p_activity)
    RF2P=TSLX*(4._r8*((Fe_3p_activity*H2PO4_1e_activity-DPF2P*FeH2PO4_2p_activity)/(S0+SQRT(S1))))
!
!     RC0P=CAPO4 <-> CA+PO4 dissociation in non-band
! CaPO4(-) <-> Ca(++)+PO4(---)
    S0=Ca_2p_activity+H0PO4_3e_activity+DPC0P
    S1=MAX(0._r8,(Ca_2p_activity-H0PO4_3e_activity)**2+2._r8*DPC0P*(Ca_2p_activity+H0PO4_3e_activity)+DPC0P**2+4._r8*DPC0P*CaPO4_1e_activity)
    RC0P=TSLX*(4._r8*((Ca_2p_activity*H0PO4_3e_activity-DPC0P*CaPO4_1e_activity)/(S0+SQRT(S1))))
!
!     RC1P=CAHPO4-CA+HPO4 dissociation in non-band
! CaHPO4 <-> Ca(++)+HPO4(--)

    S0=Ca_2p_activity+H1PO4_2e_activity+DPC1P
    S1=MAX(0._r8,(Ca_2p_activity-H1PO4_2e_activity)**2+2._r8*DPC1P*(Ca_2p_activity+H1PO4_2e_activity)+DPC1P**2+4._r8*DPC1P*CaHPO4_activity)
    RC1P=TSLX*(4._r8*((Ca_2p_activity*H1PO4_2e_activity-DPC1P*CaHPO4_activity)/(S0+SQRT(S1))))
!
!     RC2P=CaH2PO4(+)-Ca+H2PO4 dissociation in non-band
! CaH2PO4(+) <-> Ca(++)+H2PO4(-); CaH4P2O8 is the legacy dissolved variable name.
    S0=Ca_2p_activity+H2PO4_1e_activity+DPC2P
    S1=MAX(0._r8,(Ca_2p_activity-H2PO4_1e_activity)**2+2._r8*DPC2P*(Ca_2p_activity+H2PO4_1e_activity)+DPC2P**2+4._r8*DPC2P*CaH4P2O8_1p_activity)
    RC2P=TSLX*(4._r8*((Ca_2p_activity*H2PO4_1e_activity-DPC2P*CaH4P2O8_1p_activity)/(S0+SQRT(S1))))
!
!     RM1P=MGHPO4-MG+HPO4 dissociation in non-band
! MgHPO4 <-> Mg(++)+HPO4(--)
    S0   = Mg_2p_activity+H1PO4_2e_activity+DPM1P
    S1   = MAX(0._r8,(Mg_2p_activity-H1PO4_2e_activity)**2+2._r8*DPM1P*(Mg_2p_activity+H1PO4_2e_activity)+DPM1P**2+4._r8*DPM1P*MgHPO4_activity)
    RM1P = TSLX*(4._r8*((Mg_2p_activity*H1PO4_2e_activity-DPM1P*MgHPO4_activity)/(S0+SQRT(S1))))
  ELSE
    RH1P=0._r8
    H2PO4_e_to_HPO4_2e_flx=0._r8
    RH3P=0._r8
    RF1P=0._r8
    RF2P=0._r8
    RC0P=0._r8
    RC1P=0._r8
    RC2P=0._r8
    RM1P=0._r8
  ENDIF
!
!     PHOSPHORUS IN BAND SOIL ZONE
!
  IF(VLWatMicPPB.GT.ZEROS2)THEN
!
!     RH1B=HPO4-H+PO4 dissociation in band
! HPO4(--) <-> H(+)+PO4(---)
    RH1B=TSLX*(H0PO4_3e_band_activity*H_1p_activity-DPH1P*H1PO4_2e_band_activity)/(H_1p_activity+DPH1P)
!
!     RH2B=H2PO4-H+HPO4 dissociation in band
! H2PO4(-) <-> H(+)+HPO4(--)
    RH2B=TSLX*(H1PO4_2e_band_activity*H_1p_activity-DPH2P*H2PO4_1e_band_activity)/(H_1p_activity+DPH2P)
!
!     RH3B=H3PO4-H+H2PO4 dissociation in band
! H3PO4 <-> H(+) + H2PO4(-)
    RH3B=TSLX*(H2PO4_1e_band_activity*H_1p_activity-DPH3P*H3PO4_band_activity)/(H_1p_activity+DPH3P)
!
!     RF1B=FEHPO4-FE+HPO4 dissociation in band
! FeHPO4(+) <-> Fe(+++)+HPO4(--)
    S0=Fe_3p_activity+H1PO4_2e_band_activity+DPF1P
    S1=MAX(0._r8,(Fe_3p_activity-H1PO4_2e_band_activity)**2+2._r8*DPF1P*(Fe_3p_activity+H1PO4_2e_band_activity)+DPF1P**2+4._r8*DPF1P*FeHPO4_1p_band_activity)
    RF1B=TSLX*(4._r8*((Fe_3p_activity*H1PO4_2e_band_activity-DPF1P*FeHPO4_1p_band_activity)/(S0+SQRT(S1))))
!
!     RF2B=FEH2PO4-FE+H2PO4 dissociation in band
! FeH2PO4(+) <-> Fe(+++)+H2PO4(-)
    S0=Fe_3p_activity+H2PO4_1e_band_activity+DPF2P
    S1=MAX(0._r8,(Fe_3p_activity-H2PO4_1e_band_activity)**2+2._r8*DPF2P*(Fe_3p_activity+H2PO4_1e_band_activity)+DPF2P**2+4._r8*DPF2P*FeH2PO4_2p_band_activity)
    RF2B=TSLX*(4._r8*((Fe_3p_activity*H2PO4_1e_band_activity-DPF2P*FeH2PO4_2p_band_activity)/(S0+SQRT(S1))))
!
!     RC0B=CAPO4-CA+PO4 dissociation in band
! CaPO4(-) <-> Ca(++)+PO4(---)
    S0=Ca_2p_activity+H0PO4_3e_band_activity+DPC0P
    S1=MAX(0._r8,(Ca_2p_activity-H0PO4_3e_band_activity)**2+2._r8*DPC0P*(Ca_2p_activity+H0PO4_3e_band_activity)+DPC0P**2+4._r8*DPC0P*CaPO4_1e_band_activity)
    RC0B=TSLX*(4._r8*((Ca_2p_activity*H0PO4_3e_band_activity-DPC0P*CaPO4_1e_band_activity)/(S0+SQRT(S1))))
!
!     RC1B=CAHPO4-CA+HPO4 dissociation in band
! CaHPO4 <-> Ca(++)+HPO4(--)
    S0=Ca_2p_activity+H1PO4_2e_band_activity+DPC1P
    S1=MAX(0._r8,(Ca_2p_activity-H1PO4_2e_band_activity)**2+2._r8*DPC1P*(Ca_2p_activity+H1PO4_2e_band_activity)+DPC1P**2+4._r8*DPC1P*CaHPO4_band_activity)
    RC1B=TSLX*(4._r8*((Ca_2p_activity*H1PO4_2e_band_activity-DPC1P*CaHPO4_band_activity)/(S0+SQRT(S1))))
!
!     RC2B=CaH2PO4(+)-Ca+H2PO4 dissociation in band
! CaH2PO4(+) <-> Ca(++)+H2PO4(-); the two-P solid Ca(H2PO4)2 is a separate pool.
    S0=Ca_2p_activity+H2PO4_1e_band_activity+DPC2P
    S1=MAX(0._r8,(Ca_2p_activity-H2PO4_1e_band_activity)**2+2._r8*DPC2P*(Ca_2p_activity+H2PO4_1e_band_activity)+DPC2P**2+4._r8*DPC2P*CaH4P2O8_1p_band_activity)
    RC2B=TSLX*(4._r8*((Ca_2p_activity*H2PO4_1e_band_activity-DPC2P*CaH4P2O8_1p_band_activity)/(S0+SQRT(S1))))
!
!     RM1B=MGHPO4-MG+HPO4 dissociation in band
! MgHPO4 <-> Mg(++)+HPO4(--)
    S0=Mg_2p_activity+H1PO4_2e_band_activity+DPM1P
    S1=MAX(0._r8,(Mg_2p_activity-H1PO4_2e_band_activity)**2+2._r8*DPM1P*(Mg_2p_activity+H1PO4_2e_band_activity)+DPM1P**2+4._r8*DPM1P*MgHPO4_band_activity)
    RM1B=TSLX*(4._r8*((Mg_2p_activity*H1PO4_2e_band_activity-DPM1P*MgHPO4_band_activity)/(S0+SQRT(S1))))
  ELSE
    RH1B=0._r8
    RH2B=0._r8
    RH3B=0._r8
    RF1B=0._r8
    RF2B=0._r8
    RC0B=0._r8
    RC1B=0._r8
    RC2B=0._r8
    RM1B=0._r8
  ENDIF
  end subroutine SoluteDissociation
!------------------------------------------------------------------------------------------

  subroutine UpdateIonFluxCurentIter
  implicit none
!     begin_execution
!     RN4S,RN4B=net NH4 flux in non-band,band
!     RN3S,RN3B=net NH3 flux in non-band,band
!     RAL,RFE,RHY,RCA,RMG,RNA,RKA,ROH=net Al,Fe,H,Ca,Mg,Na,K,OH flux
!     RSO4,RCO3,RHCO,RCO2=net SO4,CO3,HCO3,CO2 flux
!     RAL1,RAL2,RAL3,RAL4,RALS=net AlOH,AlOH2,AlOH3,AlOH4,AlSO4
!     RFE1,RFE2,RFE3,RFE4,RFES=net FeOH,FeOH2,FeOH3,FeOH4,FeSO4
!     RHP0,RHP1,RHP2,RHP3=net PO4,HPO4,H2PO4,H3PO4 flux in non-band
!     RXH0,RXH1,RXH2,RX1P,RX2P=net R-O,R-OH,R-OH2,R-HPO4,R-H2PO4 in non-band
!     RHB0,RHB1,RHB2,RHB3=net PO4,HPO4,H2PO4,H3PO4 flux in band
!     RBH0,RBH1,RBH2,RB1P,RB2P=net R-O,R-OH,R-OH2,R-HPO4,R-H2PO4 in band
!
  RN4S=RNH4-RXN4
  RN4B=RNHB-RXNB
  RN3S=-RNH4
  RN3B=-RNHB
  RAL=-RHAL1-RXAL-RALO1-RALS-(RHA0P1+RHA0P2)*VLPO4 &
    -(RHA0B1+RHA0B2)*VLPOB
  RFE=-RHFE1-RXFE-RFEO1-RFES-(RHF0P1+RHF0P2+RF1P+RF2P)*VLPO4 &
    -(RHF0B1+RHF0B2+RF1B+RF2B)*VLPOB
  RHY=-RNH4*VLNH4-RNHB*VLNHB &
    -RXHY-RXHC+2.0_r8*(RHALO1+RHFEO1+RHCACO &
    +(RHA0P2+RHF0P2-RHA3P1-RHA4P2-RHF3P1-RHF4P2)*VLPO4 &
    +(RHA0B2+RHF0B2-RHA3B1-RHA4B2-RHF3B1-RHF4B2)*VLPOB) &
    +3.0_r8*(RHAL1+RHFE1-(RHA4P1+RHF4P1)*VLPO4 &
    -(RHF4B1+RHA4B1)*VLPOB)+4.0_r8*(RHCAH1*VLPO4+RHCHB1*VLPOB) &
    +7.0_r8*(RHCAH2*VLPO4+RHCHB2*VLPOB) &
    +RHALO2+RHFEO2-RHALO4-RHFEO4+RHCACH-RCO2Q-RHCO3 &
    +(RHA0P1-RHA2P1+RHA1P2-RHA3P2+RHF0P1-RHF2P1+RHF1P2-RHF3P2 &
    +RHCAD2-RXOH2-RXOH1-RH1P-H2PO4_e_to_HPO4_2e_flx-RH3P)*VLPO4 &
    +(RHA0B1-RHA2B1+RHA1B2-RHA3B2+RHF0B1-RHF2B1+RHF1B2-RHF3B2 &
    +RHCDB2-RXO2B-RXO1B-RH1B-RH2B-RH3B)*VLPOB
  RCA=-RPCACX-RPCASO-RXCA-RCAO-RCAC-RCAH-RCAS &
    -(H2PO4_1e_CaHPO4_dissol_flx+H2PO4_1e_CaH4P2O8_dissol_flx+RC0P+RC1P+RC2P)*VLPO4 &
    -(H2PO4_1e_CaHPO4_dissolB_flx+H2PO4_1e_CaH4P2O8_dissolB_flx+RC0B+RC1B+RC2B)*VLPOB &
    -5.0_r8*(H2PO4_1e_apatite_dissol_flx*VLPO4+H2PO4_1e_apatite_dissolB_flx*VLPOB)
  RMG=-RXMG-RMGO-RMGC-RMGH-RMGS-RM1P*VLPO4-RM1B*VLPOB
  RNA=-RXNA-RNAC-RNAS
  RKA=-RXKA-RKAS
  ROH=-RCAO-RMGO-RALO1-RALO2-RALO3-RALO4-RFEO1-RFEO2-RFEO3-RFEO4 &
    -(-H2PO4_to_XH2PO4_ROH_flx-H1PO4_to_XHPO4_ROH_flx)*VLPO4-(-H2PO4_to_XH2PO4_ROH_Bflx-RXH1B)*VLPOB
  RSO4=-RPCASO-RALS-RFES-RCAS-RMGS-RNAS-RKAS
  RCO3=-RHCAC3-RHCO3-RCAC-RMGC-RNAC
  RHCO=-RHCACH-RCO2Q-RCAH-RMGH+RHCO3
  RCO2=-RHCACO+RCO2Q
  RAL1=-RHALO1+RALO1-RALO2-(RHA1P1+RHA1P2)*VLPO4-(RHA1B1+RHA1B2)*VLPOB
  RAL2=-RHALO2+RALO2-RALO3-(RHA2P1+RHA2P2)*VLPO4-(RHA2B1+RHA2B2)*VLPOB
  RAL3=-RHALO3+RALO3-RALO4-(RHA3P1+RHA3P2)*VLPO4-(RHA3B1+RHA3B2)*VLPOB
  RAL4=-RHALO4+RALO4-(RHA4P1+RHA4P2)*VLPO4-(RHA4B1+RHA4B2)*VLPOB
  RFE1=-RHFEO1+RFEO1-RFEO2-(RHF1P1+RHF1P2)*VLPO4-(RHF1B1+RHF1B2)*VLPOB
  RFE2=-RHFEO2+RFEO2-RFEO3-(RHF2P1+RHF2P2)*VLPO4-(RHF2B1+RHF2B2)*VLPOB
  RFE3=-RHFEO3+RFEO3-RFEO4-(RHF3P1+RHF3P2)*VLPO4-(RHF3B1+RHF3B2)*VLPOB
  RFE4=-RHFEO4+RFEO4-(RHF4P1+RHF4P2)*VLPO4-(RHF4B1+RHF4B2)*VLPOB
  RHP0=-RH1P-RC0P
  RHP1=-RHA0P1-RHA1P1-RHA2P1-RHA3P1 &
    -RHA4P1-RHF0P1-RHF1P1-RHF2P1 &
    -RHF3P1-RHF4P1-RPCAD1-3.0_r8*RHCAH1-H1PO4_to_XHPO4_ROH_flx &
    +RH1P-H2PO4_e_to_HPO4_2e_flx-RF1P-RC1P-RM1P
  RHP2=-RHA0P2-RHA1P2-RHA2P2-RHA3P2 &
    -RHA4P2-RHF0P2-RHF1P2-RHF2P2 &
    -RHF3P2-RHF4P2-RHCAD2-3.0_r8*RHCAH2 &
    -2.0_r8*H2PO4_1e_CaH4P2O8_dissol_flx-H2PO4_1e_to_XH2PO4_ROH2_flx-H2PO4_to_XH2PO4_ROH_flx+H2PO4_e_to_HPO4_2e_flx-RH3P-RF2P-RC2P
  RHP3=RH3P
  RXH0=-RXOH1
  RXH1=RXOH1-RXOH2-H2PO4_to_XH2PO4_ROH_flx-H1PO4_to_XHPO4_ROH_flx
  RXH2=RXOH2-H2PO4_1e_to_XH2PO4_ROH2_flx
  RX1P=H1PO4_to_XHPO4_ROH_flx
  RX2P=H2PO4_1e_to_XH2PO4_ROH2_flx+H2PO4_to_XH2PO4_ROH_flx
  RHB0=-RH1B-RC0B
  RHB1=-RHA0B1-RHA1B1-RHA2B1-RHA3B1 &
    -RHA4B1-RHF0B1-RHF1B1-RHF2B1 &
    -RHF3B1-RHF4B1-RPCDB1-3.0_r8*RHCHB1-RXH1B &
    +RH1B-RH2B-RF1B-RC1B-RM1B
  RHB2=-RHA0B2-RHA1B2-RHA2B2-RHA3B2 &
    -RHA4B2-RHF0B2-RHF1B2-RHF2B2 &
    -RHF3B2-RHF4B2-RHCDB2-3.0_r8*RHCHB2 &
    -2.0_r8*H2PO4_1e_CaH4P2O8_dissolB_flx-H2PO4_1e_to_XH2PO4_ROH2_Bflx-H2PO4_to_XH2PO4_ROH_Bflx+RH2B-RH3B-RF2B-RC2B
  RHB3=RH3B
  RBH0=-RXO1B
  RBH1=RXO1B-RXO2B-H2PO4_to_XH2PO4_ROH_Bflx-RXH1B
  RBH2=RXO2B-H2PO4_1e_to_XH2PO4_ROH2_Bflx
  RB1P=RXH1B
  RB2P=H2PO4_1e_to_XH2PO4_ROH2_Bflx+H2PO4_to_XH2PO4_ROH_Bflx
!
  end subroutine UpdateIonFluxCurentIter
!------------------------------------------------------------------------------------------

  subroutine UpdateIonConcCurrentIter
  implicit none
  real(r8) :: CHY2_conc   !molar concentration of H(+), mol/m3
  real(r8) :: COH2_conc   !molar concentration of OH(-), mol/m3
!     begin_execution
!     UPDATE ION CONCENTRATIONS FOR CURRENT ITERATION
!     FROM TOTAL ION FLUXES
!
  NH4_1p_aque_mole_conc      = NH4_1p_aque_mole_conc+RN4S
  NH4_1p_band_mole_conc      = NH4_1p_band_mole_conc+RN4B
  NH3_aque_mole_conc         = NH3_aque_mole_conc+RN3S
  NH3_aqu_band_mole_conc     = NH3_aqu_band_mole_conc+RN3B
  Al_3p_aque_mole_conc       = Al_3p_aque_mole_conc+RAL
  Fe_3p_aque_mole_conc       = Fe_3p_aque_mole_conc+RFE
  H_1p_aque_mole_conc        = H_1p_aque_mole_conc+RHY
  Ca_2p_aque_mole_conc       = Ca_2p_aque_mole_conc+RCA
  Mg_2p_aque_mole_conc       = Mg_2p_aque_mole_conc+RMG
  Na_1p_aque_mole_conc       = Na_1p_aque_mole_conc+RNA
  K_1p_aque_mole_conc        = K_1p_aque_mole_conc+RKA
  OH_1e_aque_mole_conc       = OH_1e_aque_mole_conc+ROH
  SO4_2e_aque_mole_conc      = SO4_2e_aque_mole_conc+RSO4
  CO3_2e_aque_mole_conc      = CO3_2e_aque_mole_conc+RCO3
  HCO3_e_mole_conc           = HCO3_e_mole_conc+RHCO
  H2CO3_aque_mole_conc       = H2CO3_aque_mole_conc+RCO2
  AlOH_2p_aque_mole_conc     = AlOH_2p_aque_mole_conc+RAL1
  AlO2H2_1p_aque_mole_conc   = AlO2H2_1p_aque_mole_conc+RAL2
  AlO3H3_mole_conc           = AlO3H3_mole_conc+RAL3
  AlO4H4_1e_aque_mole_conc   = AlO4H4_1e_aque_mole_conc+RAL4
  AlSO4_1p_aque_mole_conc    = AlSO4_1p_aque_mole_conc+RALS
  FeOH_2p_aque_mole_conc     = FeOH_2p_aque_mole_conc+RFE1
  FeO2H2_p_mole_conc         = FeO2H2_p_mole_conc+RFE2
  FeO3H3_mole_conc           = FeO3H3_mole_conc+RFE3
  FeO4H4_1e_aque_mole_conc   = FeO4H4_1e_aque_mole_conc+RFE4
  FeSO4_1p_aque_mole_conc    = FeSO4_1p_aque_mole_conc+RFES
  CaO2H2_mole_conc           = CaO2H2_mole_conc+RCAO
  CaCO3_mole_conc            = CaCO3_mole_conc+RCAC
  CaHCO3_1p_aque_mole_conc   = CaHCO3_1p_aque_mole_conc+RCAH
  CaSO4_mole_conc            = CaSO4_mole_conc+RCAS
  MgOH_1p_aque_mole_conc     = MgOH_1p_aque_mole_conc+RMGO
  MgCO3_conc                 = MgCO3_conc+RMGC
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
  H0PO4_3e_band_conc         = H0PO4_3e_band_conc+RHB0
  H1PO4_2e_band_mole_conc    = H1PO4_2e_band_mole_conc+RHB1
  H2PO4_1e_band_mole_conc    = H2PO4_1e_band_mole_conc+RHB2
  H3PO4_band_conc            = H3PO4_band_conc+RHB3
  FeHPO4_1p_band_conc        = FeHPO4_1p_band_conc+RF1B
  FeH2PO4_2p_band_conc       = FeH2PO4_2p_band_conc+RF2B
  CaPO4_1e_band_conc         = CaPO4_1e_band_conc+RC0B
  CaHPO4_band_conc           = CaHPO4_band_conc+RC1B
  CaH4P2O8_1p_band_conc      = CaH4P2O8_1p_band_conc+RC2B
  MgHPO4_band_conc           = MgHPO4_band_conc+RM1B
!
!     RHHY,RHOH=H2O-H+OH equilibration
!  
  CHY2_conc  = 10.0_r8**(-PH)*1.0E+03
  ! Reuse the activity-consistent target held fixed for this reaction iteration.
  COH2_conc  = OH_1e_equil_conc
  RHHY       = CHY2_conc-H_1p_aque_mole_conc
  RHOH       = COH2_conc-OH_1e_aque_mole_conc
  H_1p_aque_mole_conc  = H_1p_aque_mole_conc+RHHY
  OH_1e_aque_mole_conc = OH_1e_aque_mole_conc+RHOH
  ! Include the pre-reaction correction once in the accumulated OH(-) ledger.
  ! Add it after updating the state, where it has already been applied above.
  RHOH=RHOH+RHOH_equil

!
!     UPDATE EXCHANGEABLE ION CONCENTRATIONS IN CURRENT
!     ITERATION FROM TOTAL ION FLUXES
!
  XNH4_mole_conc        = XNH4_mole_conc+RXN4
  XNH4_band_mole_conc   = XNH4_band_mole_conc+RXNB
  XHY1             = XHY1+RXHY
  XAl_conc         = XAl_conc+RXAL
  XFe_conc         = XFe_conc+RXFE
  XCa_conc         = XCa_conc+RXCA
  XMg_conc         = XMg_conc+RXMG
  XNa_conc         = XNa_conc+RXNA
  XK_conc          = XK_conc+RXKA
  XHC1             = XHC1+RXHC
  XAlO2H2_conc     = XAlO2H2_conc+RXALO2
  XFeO2H2_conc     = XFeO2H2_conc+RXFEO2
  XOH_mole_conc         = XOH_mole_conc+RXH0
  XROH1_mole_conc       = XROH1_mole_conc+RXH1
  XROH2_mole_conc       = XROH2_mole_conc+RXH2
  XHPO4_mole_conc       = XHPO4_mole_conc+RX1P
  XH2PO4_mole_conc      = XH2PO4_mole_conc+RX2P
  XROH1_band_mole_conc  = XROH1_band_mole_conc+RBH0
  XROH_band_mole_conc   = XROH_band_mole_conc+RBH1
  XROH2_band_mole_conc  = XROH2_band_mole_conc+RBH2
  XHPO4_band_mole_conc  = XHPO4_band_mole_conc+RB1P
  XH2PO4_band_mole_conc = XH2PO4_band_mole_conc+RB2P
!
!     UPDATE PRECIPITATE CONCENTRATIONS IN CURRENT
!     ITERATION FROM TOTAL ION FLUXES
!
  Precp_AlO3H3_conc        = Precp_AlO3H3_conc+RPALOX       !amorphous Al(OH)3 precipitation
  Precp_FeO3H3_conc        = Precp_FeO3H3_conc+RPFEOX       !Fe(OH)3 precipitation
  Precp_CaCO3_conc         = Precp_CaCO3_conc+RPCACX       !Calcite CaCO3 precipitation
  Precp_CaSO4_conc         = Precp_CaSO4_conc+RPCASO       !Gypsum CaSO4 precipitation
  Precp_AlPO4_mole_conc         = Precp_AlPO4_mole_conc+H2PO4_1e_AlPO4_dissol_flx       !Variscite AlPO4 precipitation non-band
  Precp_FePO4_mole_conc         = Precp_FePO4_mole_conc+H2PO4_1e_FePO4_dissol_flx       !FePO4 precipitation
  Precp_CaHPO4_mole_conc        = Precp_CaHPO4_mole_conc+H2PO4_1e_CaHPO4_dissol_flx       !CaHPO4 precpitation
  Precp_Ca5P3O12O3H3_mole_conc  = Precp_Ca5P3O12O3H3_mole_conc+H2PO4_1e_apatite_dissol_flx       !Ca5(PO4)3OH (hydroxyapatite) precipitation
  Precp_CaH4P2O8_mole_conc      = Precp_CaH4P2O8_mole_conc+H2PO4_1e_CaH4P2O8_dissol_flx       !Ca(H2PO4)2 precipitation non-band
  PrecpB_AlPO4_mole_conc        = PrecpB_AlPO4_mole_conc+H2PO4_1e_AlPO4_dissolB_flx       !AlPO4 precpitation band
  PrecpB_FePO4_mole_conc         = PrecpB_FePO4_mole_conc+H2PO4_1e_FePO4_dissolB_flx       !FePO4 precpitation band
  PrecpB_CaHPO4_mole_conc       = PrecpB_CaHPO4_mole_conc+H2PO4_1e_CaHPO4_dissolB_flx       !CaHPO4 precipitation band
  PrecpB_Ca5P3O12O3H3_mole_conc = PrecpB_Ca5P3O12O3H3_mole_conc+H2PO4_1e_apatite_dissolB_flx       !Ca5(PO4)3OH hydroxyapatite precpitation band
  PrecpB_CaH4P2O8_mole_conc     = PrecpB_CaH4P2O8_mole_conc+H2PO4_1e_CaH4P2O8_dissolB_flx       !Ca(H2PO4)2 precipitation band
  end subroutine UpdateIonConcCurrentIter
!------------------------------------------------------------------------------------------

  subroutine AccumulateIonFlux
  implicit none
!     begin_execution
!     ACCUMULATE TOTAL ION FLUXES FOR ALL ITERATIONS
!
!
  TRChem_NH4_soil_mole             = TRChem_NH4_soil_mole+RN4S    !net NH4 flux in non-band
  TRChem_NH4_band_soil_mole        = TRChem_NH4_band_soil_mole+RN4B    !net NH4 flux in band
  TRChem_NH3_soil_mole             = TRChem_NH3_soil_mole+RN3S    !net NH3 flux in non-band
  TRChem_NH3_band_soil             = TRChem_NH3_band_soil+RN3B    !net NH3 flux in band
  TRChem_Al_3p_soil                = TRChem_Al_3p_soil+RAL       !total Al flux
  TRChem_Fe_3p_soil                = TRChem_Fe_3p_soil+RFE       !total Fe(3+) flux
  TRChem_H_p_soil                  = TRChem_H_p_soil+RHY+RHHY  !total H(+) flux
  TRChem_Ca_2p_soil                = TRChem_Ca_2p_soil+RCA       !total Ca(2+) flux
  TRChem_Mg_2p_soil                = TRChem_Mg_2p_soil+RMG       !total Mg(2+) flux
  TRChem_Na_p_soil                 = TRChem_Na_p_soil+RNA       !total Na(+) flux
  TRChem_K_1p_soil                 = TRChem_K_1p_soil+RKA       !total K(+) flux
  TRChem_OH_1e_soil                = TRChem_OH_1e_soil+ROH+RHOH  !total OH(-) flux
  TRChem_SO4_2e_soil               = TRChem_SO4_2e_soil+RSO4     !total SO4(2-) production flux
  TRChem_CO3_2e_soil               = TRChem_CO3_2e_soil+RCO3     !total CO3(2-) production flux
  TRChem_HCO3_soil                 = TRChem_HCO3_soil+RHCO       !total HCO3(-) production flux
  TRChem_CO2_gchem_soil            = TRChem_CO2_gchem_soil+RCO2    !total CO2 flux due to dissociation
  TRChem_AlOH_soil                 = TRChem_AlOH_soil+RAL1    !total Al(OH)(2+) flux
  TRChem_AlO2H2_soil               = TRChem_AlO2H2_soil+RAL2    !total Al(OH)2(+) flux
  TRChem_AlO3H3_soil               = TRChem_AlO3H3_soil+RAL3    !total Al(OH)3 flux
  TRChem_AlO4H4_soil               = TRChem_AlO4H4_soil+RAL4    !total Al(OH)4(-) flux
  TRChem_AlSO4_soil                = TRChem_AlSO4_soil+RALS    !total Al(SO4)(+) flux
  TRChem_FeOH_soil                 = TRChem_FeOH_soil+RFE1    !total Fe(OH)(2+) flux
  TRChem_FeO2H2_soil               = TRChem_FeO2H2_soil+RFE2    !total Fe(OH)2(-) flux
  TRChem_FeO3H3_soil               = TRChem_FeO3H3_soil+RFE3    !total Fe(OH)3 flux
  TRChem_FeO4H4_soil               = TRChem_FeO4H4_soil+RFE4    !total Fe(OH4) flux
  TRChem_FeSO4_soil                = TRChem_FeSO4_soil+RFES    !total FeSO4(+) flux
  TRChem_CaOH_soil                 = TRChem_CaOH_soil+RCAO    !total Ca(OH)(+) flux
  TRChem_CaCO3_soil                = TRChem_CaCO3_soil+RCAC    !total CaCO3 flux
  TRChem_CaHCO3_soil               = TRChem_CaHCO3_soil+RCAH    !CaHCO3 flux
  TRChem_CaSO4_soil                = TRChem_CaSO4_soil+RCAS    !CaSO4 flux
  TRChem_MgOH_soil                 = TRChem_MgOH_soil+RMGO    !MgOH  (+)
  TRChem_MgCO3_soil                = TRChem_MgCO3_soil+RMGC    !MgCO3
  TRChem_MgHCO3_soil               = TRChem_MgHCO3_soil+RMGH    !MgHCO3 (+)
  TRChem_MgSO4_soil                = TRChem_MgSO4_soil+RMGS    !MgSO4
  TRChem_NaCO3_soil                = TRChem_NaCO3_soil+RNAC    !NaCO3(-)
  TRChem_NaSO4_soil                = TRChem_NaSO4_soil+RNAS    !NaSO4(-)
  TRChem_KSO4_soil                 = TRChem_KSO4_soil+RKAS    !KSO4(-)
  TRChem_PO4_soil                  = TRChem_PO4_soil+RHP0    !PO4 (3-)
  TRChem_H1PO4_soil                = TRChem_H1PO4_soil+RHP1    !HPO4 (2-)
  TRChem_H2PO4_soil                = TRChem_H2PO4_soil+RHP2    !H2PO4 (-)
  TRChem_H3PO4_sorbed_soil         = TRChem_H3PO4_sorbed_soil+RHP3    !H3PO4
  TRChem_FeHPO4_soil               = TRChem_FeHPO4_soil+RF1P    !FeHPO4
  TRChem_FeH2PO4_soil              = TRChem_FeH2PO4_soil+RF2P    !FeH2PO4
  TRChem_CaPO4_soil                = TRChem_CaPO4_soil+RC0P
  TRChem_CaHPO4_soil               = TRChem_CaHPO4_soil+RC1P
  TRChem_CaH4P2O8_soil             = TRChem_CaH4P2O8_soil+RC2P
  TRChem_MgHPO4_soil               = TRChem_MgHPO4_soil+RM1P
  TRChem_PO4_band_soil             = TRChem_PO4_band_soil+RHB0
  TRChem_H1PO4_band_soil           = TRChem_H1PO4_band_soil+RHB1
  TRChem_H2PO4_band_soil           = TRChem_H2PO4_band_soil+RHB2
  TRChem_H3PO4_band_soil           = TRChem_H3PO4_band_soil+RHB3
  TRChem_FeHPO4_band_soil          = TRChem_FeHPO4_band_soil+RF1B
  TRChem_FeH2PO4_band_soil         = TRChem_FeH2PO4_band_soil+RF2B
  TRChem_CaPO4_band_soil           = TRChem_CaPO4_band_soil+RC0B
  TRChem_CaHPO4_band_soil          = TRChem_CaHPO4_band_soil+RC1B
  TRChem_CaH4P2O8_band_soil        = TRChem_CaH4P2O8_band_soil+RC2B
  TRChem_MgHPO4_band_soil          = TRChem_MgHPO4_band_soil+RM1B
  TRChem_NH4_sorbed_soil           = TRChem_NH4_sorbed_soil+RXN4
  TRChem_NH4_sorbed_band_soil      = TRChem_NH4_sorbed_band_soil+RXNB
  TRChem_H_p_sorbed_soil           = TRChem_H_p_sorbed_soil+RXHY
  TRChem_Al_sorbed_soil            = TRChem_Al_sorbed_soil+RXAL
  TRChem_Fe_sorbed_soil            = TRChem_Fe_sorbed_soil+RXFE
  TRChem_Ca_sorbed_soil            = TRChem_Ca_sorbed_soil+RXCA
  TRChem_Mg_sorbed_soil            = TRChem_Mg_sorbed_soil+RXMG
  TRChem_Na_sorbed_soil            = TRChem_Na_sorbed_soil+RXNA
  TRChem_K_sorbed_soil             = TRChem_K_sorbed_soil+RXKA
  TRChem_HCO3_sorbed_soil          = TRChem_HCO3_sorbed_soil+RXHC
  TRChem_AlO2H2_sorbed_soil        = TRChem_AlO2H2_sorbed_soil+RXALO2
  TRChem_FeO2H2_sorbed_soil        = TRChem_FeO2H2_sorbed_soil+RXFEO2
  TRChem_RO_sorbed_soil            = TRChem_RO_sorbed_soil+RXH0
  TRChem_ROH_sorbed_soil           = TRChem_ROH_sorbed_soil+RXH1
  TRChem_ROH2_sorbed_soil          = TRChem_ROH2_sorbed_soil+RXH2
  TRChem_RHPO4_sorbed_soil         = TRChem_RHPO4_sorbed_soil+RX1P
  TRChem_RH2PO4_sorbed_soil        = TRChem_RH2PO4_sorbed_soil+RX2P
  TRChem_RO_sorbed_band_soil       = TRChem_RO_sorbed_band_soil+RBH0
  TRChem_ROH_sorbed_band_soil      = TRChem_ROH_sorbed_band_soil+RBH1
  TRChem_ROH2_sorbed_band_soil     = TRChem_ROH2_sorbed_band_soil+RBH2
  TRChem_RHPO4_sorbed_band_soil    = TRChem_RHPO4_sorbed_band_soil+RB1P
  TRChem_RH2PO4_sorbed_band_soil   = TRChem_RH2PO4_sorbed_band_soil+RB2P
  TRChem_AlOH3_precip_soil         = TRChem_AlOH3_precip_soil+RPALOX
  TRChem_FeOH3_precip_soil         = TRChem_FeOH3_precip_soil+RPFEOX
  TRChem_CaCO3_precip_soil         = TRChem_CaCO3_precip_soil+RPCACX
  TRChem_CaSO4_precip_soil         = TRChem_CaSO4_precip_soil+RPCASO
  TRChem_AlPO4_precip_soil         = TRChem_AlPO4_precip_soil+H2PO4_1e_AlPO4_dissol_flx
  TRChem_FePO4_precip_soil         = TRChem_FePO4_precip_soil+H2PO4_1e_FePO4_dissol_flx
  TRChem_CaHPO4_precip_soil        = TRChem_CaHPO4_precip_soil+H2PO4_1e_CaHPO4_dissol_flx
  TRChem_apatite_precip_soil       = TRChem_apatite_precip_soil+H2PO4_1e_apatite_dissol_flx
  TRChem_CaH4P2O8_precip_soil      = TRChem_CaH4P2O8_precip_soil+H2PO4_1e_CaH4P2O8_dissol_flx
  TRChem_AlPO4_precip_band_soil    = TRChem_AlPO4_precip_band_soil+H2PO4_1e_AlPO4_dissolB_flx
  TRChem_FePO4_precip_band_soil    = TRChem_FePO4_precip_band_soil+H2PO4_1e_FePO4_dissolB_flx
  TRChem_CaHPO4_precip_band_soil   = TRChem_CaHPO4_precip_band_soil+H2PO4_1e_CaHPO4_dissolB_flx
  TRChem_apatite_precip_band_soil  = TRChem_apatite_precip_band_soil+H2PO4_1e_apatite_dissolB_flx
  TRChem_CaH4P2O8_precip_band_soil = TRChem_CaH4P2O8_precip_band_soil+H2PO4_1e_CaH4P2O8_dissolB_flx

  end subroutine AccumulateIonFlux

end module SaltChemEquilibriaMod
