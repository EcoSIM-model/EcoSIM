module GeochemAPI

  use data_kind_mod,      only: r8 => DAT_KIND_R8
  use SoluteChemDataType, only: solute_flx_type, chem_var_type
  use EcoSimConst,        only: catomw
  use DebugToolMod
  use SoluteMod
  use AqueChemDatatype
  use SoilPropertyDataType
  use SoilBGCDataType
  use EcoSIMCtrlDataType
  use SOMDataType
  use EcoSIMSolverPar
  use SoilWaterDataType
  use GridDataType
  implicit none
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: soluteModel
  contains


  subroutine soluteModel(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES ALL SOLUTE TRANSFORMATIONS
!     FROM THERMODYNAMIC EQUILIBRIA
!
  implicit none

  integer, intent(in) :: I, J
  integer, intent(in) :: NHW, NHE, NVN, NVS

  ! declaration of local variables
  integer :: L,NY,NX,NPI
  real(r8) :: RHP1,RHP2,RN3S,RN4S
  type(solute_flx_type) :: solflx
  type(chem_var_type) :: chemvar


    chemvar%ZMG    = 0._r8
    chemvar%ZNA    = 0._r8
    chemvar%ZKA    = 0._r8
    chemvar%ZMG    = 0._r8
    chemvar%ZNA    = 0._r8
    chemvar%ZKA    = 0._r8
    chemvar%ZHY    = 0._r8
    chemvar%ZOH    = 0._r8
    chemvar%ZAL    = 0._r8
    chemvar%ZFE    = 0._r8
    chemvar%ZCA    = 0._r8
    chemvar%ZSO4   = 0._r8
    chemvar%ZCL    = 0._r8
    chemvar%ZCO3   = 0._r8
    chemvar%ZHCO3  = 0._r8
    chemvar%ZALOH1 = 0._r8
    chemvar%ZALOH2 = 0._r8
    chemvar%ZALOH3 = 0._r8
    chemvar%ZALOH4 = 0._r8
    chemvar%ZALS   = 0._r8
    chemvar%ZFEOH1 = 0._r8
    chemvar%ZFEOH2 = 0._r8
    chemvar%ZFEOH3 = 0._r8
    chemvar%ZFEOH4 = 0._r8
    chemvar%ZFES   = 0._r8
    chemvar%ZCAO   = 0._r8
    chemvar%ZCAC   = 0._r8
    chemvar%ZCAH   = 0._r8
    chemvar%ZCAS   = 0._r8
    chemvar%ZMGO   = 0._r8
    chemvar%ZMGC   = 0._r8
    chemvar%ZMGH   = 0._r8
    chemvar%ZMGS   = 0._r8
    chemvar%ZNAC   = 0._r8
    chemvar%ZNAS   = 0._r8
    chemvar%ZKAS   = 0._r8
    chemvar%H0PO4  = 0._r8
    chemvar%H3PO4  = 0._r8
    chemvar%ZFE1P  = 0._r8
    chemvar%ZFE2P  = 0._r8
    chemvar%ZCA0P  = 0._r8
    chemvar%ZCA1P  = 0._r8
    chemvar%ZCA2P  = 0._r8
    chemvar%ZMG1P  = 0._r8
    chemvar%H0POB  = 0._r8
    chemvar%H3POB  = 0._r8
    chemvar%ZFE1PB = 0._r8
    chemvar%ZFE2PB = 0._r8
    chemvar%ZCA0PB = 0._r8
    chemvar%ZCA1PB = 0._r8
    chemvar%ZCA2PB = 0._r8
    chemvar%ZMG1PB = 0._r8

!     begin_execution
  NPI=INT(NPH/2)
  DO   NX=NHW,NHE
    DO   NY=NVN,NVS
      DO   L=NU_col(NY,NX),NL_col(NY,NX)
        IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX).AND.VLWatMicPM_vr(NPH,L,NY,NX).GT.ZEROS2(NY,NX))THEN
!
!     WATER VOLUME IN NON-BAND AND BAND SOIL ZONES
!
!     VLWatMicPM=soil water volume
!     VLNH4,VLNHB=fractions of soil volume in NH4 non-band,band
!     VLNO3,VLNOB=fractions of soil volume in N03 non-band,band
!     VLPO4,VLPOB=fractions of soil volume in H2PO4 non-band,band
!     SoilMicPMassLayer=soil mass
!
          chemvar%VLWatMicPNH = VLWatMicPM_vr(NPH,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
          chemvar%VLWatMicPNB = VLWatMicPM_vr(NPH,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
          chemvar%VLWatMicPNO = VLWatMicPM_vr(NPH,L,NY,NX)*trcs_VLN_vr(ids_NO3,L,NY,NX)
          chemvar%VLWatMicPNZ = VLWatMicPM_vr(NPH,L,NY,NX)*trcs_VLN_vr(ids_NO3B,L,NY,NX)
          chemvar%VLWatMicPPO = VLWatMicPM_vr(NPH,L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          chemvar%VLWatMicPPB = VLWatMicPM_vr(NPH,L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
            chemvar%SoilMicPMassLayerX = VLSoilMicPMass_vr(L,NY,NX)
            chemvar%BKVLNH             = VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
            chemvar%BKVLNB             = VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
            chemvar%BKVLNO             = VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_NO3,L,NY,NX)
            chemvar%BKVLNZ             = VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_NO3B,L,NY,NX)
            chemvar%BKVLPO             = VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
            chemvar%BKVLPB             = VLSoilMicPMass_vr(L,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          ELSE
            chemvar%SoilMicPMassLayerX = VLMicP_vr(L,NY,NX)
            chemvar%BKVLNH             = chemvar%VLWatMicPNH
            chemvar%BKVLNB             = chemvar%VLWatMicPNB
            chemvar%BKVLNO             = chemvar%VLWatMicPNO
            chemvar%BKVLNZ             = chemvar%VLWatMicPNZ
            chemvar%BKVLPO             = chemvar%VLWatMicPPO
            chemvar%BKVLPB             = chemvar%VLWatMicPPB
          ENDIF
          !set other variables of chemvar
          call UpdateSoilFertlizer(I,J,L,NY,NX,chemvar)

          call GeochemAPISend(L,NY,NX,chemvar,solflx)

          call GeoChemEquilibria(I,J,L,NY,NX,chemvar, solflx)

          call GeochemAPIRecv(L,NY,NX,solflx)

          call UpdateFertilizerBand(L,NY,NX)
        ENDIF
      ENDDO
!
!     SURFACE RESIDUE
!
      call UpdateSurfResidueSolute(NX,NY)
!
    ENDDO
  ENDDO
  end subroutine soluteModel

!------------------------------------------------------------------------

  subroutine GeochemAPISend(L,NY,NX,chemvar,solflx)

  implicit none
  integer, intent(in) :: L,NY,NX
  type(solute_flx_type), intent(inout) :: solflx
  type(chem_var_type), intent(inout) :: chemvar

  chemvar%XCEC       = trcx_solml_vr(idx_CEC,L,NY,NX)
  chemvar%PH         = PH_vr(L,NY,NX)
  chemvar%CAL        = CAL_vr(L,NY,NX)
  chemvar%CFE        = CFE_vr(L,NY,NX)
  chemvar%VLWatMicPM = VLWatMicPM_vr(NPH,L,NY,NX)
  if(salt_model)then
    chemvar%ZMG    = trcSalt_solml_vr(idsalt_Mg,L,NY,NX)
    chemvar%ZNA    = trcSalt_solml_vr(idsalt_Na,L,NY,NX)
    chemvar%ZKA    = trcSalt_solml_vr(idsalt_K,L,NY,NX)
    chemvar%ZMG    = trcSalt_solml_vr(idsalt_Mg,L,NY,NX)
    chemvar%ZNA    = trcSalt_solml_vr(idsalt_Na,L,NY,NX)
    chemvar%ZKA    = trcSalt_solml_vr(idsalt_K,L,NY,NX)
    chemvar%ZHY    = trcSalt_solml_vr(idsalt_Hp,L,NY,NX)
    chemvar%ZOH    = trcSalt_solml_vr(idsalt_OH,L,NY,NX)
    chemvar%ZAL    = trcSalt_solml_vr(idsalt_Al,L,NY,NX)
    chemvar%ZFE    = trcSalt_solml_vr(idsalt_Fe,L,NY,NX)
    chemvar%ZCA    = trcSalt_solml_vr(idsalt_Ca,L,NY,NX)
    chemvar%ZSO4   = trcSalt_solml_vr(idsalt_SO4,L,NY,NX)
    chemvar%ZCL    = trcSalt_solml_vr(idsalt_Cl,L,NY,NX)
    chemvar%ZCO3   = trcSalt_solml_vr(idsalt_CO3,L,NY,NX)
    chemvar%ZHCO3  = trcSalt_solml_vr(idsalt_HCO3,L,NY,NX)
    chemvar%ZALOH1 = trcSalt_solml_vr(idsalt_AlOH,L,NY,NX)
    chemvar%ZALOH2 = trcSalt_solml_vr(idsalt_AlOH2,L,NY,NX)
    chemvar%ZALOH3 = trcSalt_solml_vr(idsalt_AlOH3,L,NY,NX)
    chemvar%ZALOH4 = trcSalt_solml_vr(idsalt_AlOH4,L,NY,NX)
    chemvar%ZALS   = trcSalt_solml_vr(idsalt_AlSO4,L,NY,NX)
    chemvar%ZFEOH1 = trcSalt_solml_vr(idsalt_FeOH,L,NY,NX)
    chemvar%ZFEOH2 = trcSalt_solml_vr(idsalt_FeOH2,L,NY,NX)
    chemvar%ZFEOH3 = trcSalt_solml_vr(idsalt_FeOH3,L,NY,NX)
    chemvar%ZFEOH4 = trcSalt_solml_vr(idsalt_FeOH4,L,NY,NX)
    chemvar%ZFES   = trcSalt_solml_vr(idsalt_FeSO4,L,NY,NX)
    chemvar%ZCAO   = trcSalt_solml_vr(idsalt_CaOH,L,NY,NX)
    chemvar%ZCAC   = trcSalt_solml_vr(idsalt_CaCO3,L,NY,NX)
    chemvar%ZCAH   = trcSalt_solml_vr(idsalt_CaHCO3,L,NY,NX)
    chemvar%ZCAS   = trcSalt_solml_vr(idsalt_CaSO4,L,NY,NX)
    chemvar%ZMGO   = trcSalt_solml_vr(idsalt_MgOH2,L,NY,NX)
    chemvar%ZMGC   = trcSalt_solml_vr(idsalt_MgCO3,L,NY,NX)
    chemvar%ZMGH   = trcSalt_solml_vr(idsalt_MgHCO3,L,NY,NX)
    chemvar%ZMGS   = trcSalt_solml_vr(idsalt_MgSO4,L,NY,NX)
    chemvar%ZNAC   = trcSalt_solml_vr(idsalt_NaCO3,L,NY,NX)
    chemvar%ZNAS   = trcSalt_solml_vr(idsalt_NaSO4,L,NY,NX)
    chemvar%ZKAS   = trcSalt_solml_vr(idsalt_KSO4,L,NY,NX)
    chemvar%H0PO4  = trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)
    chemvar%H3PO4  = trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX)
    chemvar%ZFE1P  = trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)
    chemvar%ZFE2P  = trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)
    chemvar%ZCA0P  = trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX)
    chemvar%ZCA1P  = trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX)
    chemvar%ZCA2P  = trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)
    chemvar%ZMG1P  = trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX)
    chemvar%H0POB  = trcSalt_solml_vr(idsalt_H0PO4B,L,NY,NX)
    chemvar%H3POB  = trcSalt_solml_vr(idsalt_H3PO4B,L,NY,NX)
    chemvar%ZFE1PB = trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX)
    chemvar%ZFE2PB = trcSalt_solml_vr(idsalt_FeH2PO4B,L,NY,NX)
    chemvar%ZCA0PB = trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX)
    chemvar%ZCA1PB = trcSalt_solml_vr(idsalt_CaHPO4B,L,NY,NX)
    chemvar%ZCA2PB = trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX)
    chemvar%ZMG1PB = trcSalt_solml_vr(idsalt_MgHPO4B,L,NY,NX)
  endif
  chemvar%CCA               = CCA_vr(L,NY,NX)
  chemvar%SoilMicPMassLayer = VLSoilMicPMass_vr(L,NY,NX)
  chemvar%XAEC              = trcx_solml_vr(idx_AEC,L,NY,NX)
  chemvar%VLNH4             = trcs_VLN_vr(ids_NH4,L,NY,NX)
  chemvar%GKC4              = GKC4_vr(L,NY,NX)
  chemvar%VLNHB             = trcs_VLN_vr(ids_NH4B,L,NY,NX)
  chemvar%GKCA              = GKCA_vr(L,NY,NX)
  chemvar%GKCH              = GKCH_vr(L,NY,NX)
  chemvar%GKCM              = GKCM_vr(L,NY,NX)
  chemvar%GKCK              = GKCK_vr(L,NY,NX)
  chemvar%GKCN              = GKCN_vr(L,NY,NX)
  chemvar%ZNO3S             = trcs_solml_vr(ids_NO3,L,NY,NX)
  chemvar%ZNO3B             = trcs_solml_vr(ids_NO3B,L,NY,NX)
  chemvar%RProd_Hp          = RProd_Hp_vr(L,NY,NX)
  chemvar%CO2S              = trcs_solml_vr(idg_CO2,L,NY,NX)
  
  chemvar%XHY   = trcx_solml_vr(idx_Hp,L,NY,NX)
  chemvar%XAL   = trcx_solml_vr(idx_Al,L,NY,NX)
  chemvar%XFE   = trcx_solml_vr(idx_Fe,L,NY,NX)
  chemvar%XCA   = trcx_solml_vr(idx_Ca,L,NY,NX)
  chemvar%XMG   = trcx_solml_vr(idx_Mg,L,NY,NX)
  chemvar%XNA   = trcx_solml_vr(idx_Na,L,NY,NX)
  chemvar%XKA   = trcx_solml_vr(idx_K,L,NY,NX)
  chemvar%XHC   = trcx_solml_vr(idx_COOH,L,NY,NX)
  chemvar%XALO2 = trcx_solml_vr(idx_AlOH2,L,NY,NX)
  chemvar%XFEO2 = trcx_solml_vr(idx_FeOH2,L,NY,NX)
  chemvar%ORGC  = SoilOrgM_vr(ielmc,L,NY,NX)
  chemvar%PALOH = trcp_saltpml_vr(idsp_AlOH3,L,NY,NX)
  chemvar%PFEOH = trcp_saltpml_vr(idsp_FeOH3,L,NY,NX)
  chemvar%PCACO = trcp_saltpml_vr(idsp_CaCO3,L,NY,NX)
  chemvar%PCASO = trcp_saltpml_vr(idsp_CaSO4,L,NY,NX)
  chemvar%VLPOB = trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
  chemvar%VLPO4 = trcs_VLN_vr(ids_H1PO4,L,NY,NX)
  chemvar%VLNOB = trcs_VLN_vr(ids_NO3B,L,NY,NX)
  chemvar%VLNO3 = trcs_VLN_vr(ids_NO3,L,NY,NX)

  end subroutine GeochemAPISend

!-----------------------------------------------------------------------

  subroutine GeochemAPIRecv(L,NY,NX,solflx)

  implicit none
  integer, intent(in) :: L,NY,NX
  type(solute_flx_type), intent(in) :: solflx
  character(len=*), parameter :: subname='GeochemAPIRecv'

  call PrintInfo('beg '//subname)
  !units are now of [mol d-2], which will be converted into [g d-2] later.
  trcn_RprodChem_soil_vr(ids_NH4,L,NY,NX)         = trcn_RprodChem_soil_vr(ids_NH4,L,NY,NX)       +solflx%TRChem_NH4_soil_mole
  trcn_RProdChem_band_soil_vr(ids_NH4B,L,NY,NX)   = trcn_RProdChem_band_soil_vr(ids_NH4B,L,NY,NX)+solflx%TRChem_NH4_band_soil_mole
  TRProd_chem_sol_NH3_soil_vr(L,NY,NX)            = TRProd_chem_sol_NH3_soil_vr(L,NY,NX)+solflx%TRChem_NH3_soil_mole
  trcn_RProdChem_band_soil_vr(idg_NH3B,L,NY,NX)   = trcn_RProdChem_band_soil_vr(idg_NH3B,L,NY,NX)+solflx%TRChem_NH3_band_soil
  trcn_RprodChem_soil_vr(ids_H1PO4,L,NY,NX)       = solflx%TRChem_H1PO4_soil
  trcn_RprodChem_soil_vr(ids_H2PO4,L,NY,NX)       = solflx%TRChem_H2PO4_soil
  trcn_RProdChem_band_soil_vr(ids_H1PO4B,L,NY,NX) = solflx%TRChem_H1PO4_band_soil
  trcn_RProdChem_band_soil_vr(ids_H2PO4B,L,NY,NX) = solflx%TRChem_H2PO4_band_soil
  trcx_TRSoilChem_vr(idx_NH4,L,NY,NX)             = solflx%TRChem_NH4_sorbed_soil
  trcx_TRSoilChem_vr(idx_NH4B,L,NY,NX)            = solflx%TRChem_NH4_sorbed_band_soil
  trcx_TRSoilChem_vr(idx_OH,L,NY,NX)              = solflx%TRChem_ROH_sorbed_soil
  trcx_TRSoilChem_vr(idx_OHp,L,NY,NX)             = solflx%TRChem_ROH2_sorbed_soil
  trcx_TRSoilChem_vr(idx_HPO4,L,NY,NX)            = solflx%TRChem_RHPO4_sorbed_soil
  trcx_TRSoilChem_vr(idx_H2PO4,L,NY,NX)           = solflx%TRChem_RH2PO4_sorbed_soil
  trcx_TRSoilChem_vr(idx_OHB,L,NY,NX)             = solflx%TRChem_ROH_sorbed_band_soil
  trcx_TRSoilChem_vr(idx_OHpB,L,NY,NX)            = solflx%TRChem_ROH2_sorbed_band_soil
  trcx_TRSoilChem_vr(idx_HPO4B,L,NY,NX)           = solflx%TRChem_RHPO4_sorbed_band_soil
  trcx_TRSoilChem_vr(idx_H2PO4B,L,NY,NX)          = solflx%TRChem_RH2PO4_sorbed_band_soil
  trcp_RChem_soil_vr(idsp_AlPO4,L,NY,NX)          = solflx%TRChem_AlPO4_precip_soil
  trcp_RChem_soil_vr(idsp_FePO4,L,NY,NX)          = solflx%TRChem_FePO4_precip_soil
  trcp_RChem_soil_vr(idsp_CaHPO4,L,NY,NX)         = solflx%TRChem_CaHPO4_precip_soil
  trcp_RChem_soil_vr(idsp_Apatite,L,NY,NX)             = solflx%TRChem_apatite_precip_soil
  trcp_RChem_soil_vr(idsp_CaH4P2O8,L,NY,NX)       = solflx%TRChem_CaH4P2O8_precip_soil
  trcp_RChem_soil_vr(idsp_AlPO4B,L,NY,NX)         = solflx%TRChem_AlPO4_precip_band_soil
  trcp_RChem_soil_vr(idsp_FePO4B,L,NY,NX)         = solflx%TRChem_FePO4_precip_band_soil
  trcp_RChem_soil_vr(idsp_CaHPO4B,L,NY,NX)        = solflx%TRChem_CaHPO4_precip_band_soil
  trcp_RChem_soil_vr(idsp_ApatiteBand,L,NY,NX)    = solflx%TRChem_apatite_precip_band_soil
  trcp_RChem_soil_vr(idsp_CaH4P2O8B,L,NY,NX)      = solflx%TRChem_CaH4P2O8_precip_band_soil
  TProd_CO2_geochem_soil_vr(L,NY,NX)              = solflx%TRChem_CO2_gchem_soil
  if(salt_model)then
    trcSalt_RGeoChem_flx_vr(idsalt_Al,L,NY,NX)        = solflx%TRChem_Al_3p_soil
    trcSalt_RGeoChem_flx_vr(idsalt_Fe,L,NY,NX)        = solflx%TRChem_Fe_3p_soil
    trcSalt_RGeoChem_flx_vr(idsalt_Hp,L,NY,NX)        = solflx%TRChem_H_p_soil
    trcSalt_RGeoChem_flx_vr(idsalt_Ca,L,NY,NX)        = solflx%TRChem_Ca_2p_soil
    trcSalt_RGeoChem_flx_vr(idsalt_Mg,L,NY,NX)        = solflx%TRChem_Mg_2p_soil
    trcSalt_RGeoChem_flx_vr(idsalt_Na,L,NY,NX)        = solflx%TRChem_Na_p_soil
    trcSalt_RGeoChem_flx_vr(idsalt_K,L,NY,NX)         = solflx%TRChem_K_1p_soil
    trcSalt_RGeoChem_flx_vr(idsalt_OH,L,NY,NX)        = solflx%TRChem_OH_1e_soil
    trcSalt_RGeoChem_flx_vr(idsalt_SO4,L,NY,NX)       = solflx%TRChem_SO4_2e_soil
    trcSalt_RGeoChem_flx_vr(idsalt_CO3,L,NY,NX)       = solflx%TRChem_CO3_2e_soil
    trcSalt_RGeoChem_flx_vr(idsalt_HCO3,L,NY,NX)      = solflx%TRChem_HCO3_soil
    trcSalt_RGeoChem_flx_vr(idsalt_AlOH,L,NY,NX)      = solflx%TRChem_AlOH_soil
    trcSalt_RGeoChem_flx_vr(idsalt_AlOH2,L,NY,NX)     = solflx%TRChem_AlO2H2_soil
    trcSalt_RGeoChem_flx_vr(idsalt_AlOH3,L,NY,NX)     = solflx%TRChem_AlO3H3_soil
    trcSalt_RGeoChem_flx_vr(idsalt_AlOH4,L,NY,NX)     = solflx%TRChem_AlO4H4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_AlSO4,L,NY,NX)     = solflx%TRChem_AlSO4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_FeOH,L,NY,NX)      = solflx%TRChem_FeOH_soil
    trcSalt_RGeoChem_flx_vr(idsalt_FeOH2,L,NY,NX)     = solflx%TRChem_FeO2H2_soil
    trcSalt_RGeoChem_flx_vr(idsalt_FeOH3,L,NY,NX)     = solflx%TRChem_FeO3H3_soil
    trcSalt_RGeoChem_flx_vr(idsalt_FeOH4,L,NY,NX)     = solflx%TRChem_FeO4H4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_FeSO4,L,NY,NX)     = solflx%TRChem_FeSO4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_CaOH,L,NY,NX)      = solflx%TRChem_CaOH_soil
    trcSalt_RGeoChem_flx_vr(idsalt_CaCO3,L,NY,NX)     = solflx%TRChem_CaCO3_soil
    trcSalt_RGeoChem_flx_vr(idsalt_CaHCO3,L,NY,NX)    = solflx%TRChem_CaHCO3_soil
    trcSalt_RGeoChem_flx_vr(idsalt_CaSO4,L,NY,NX)     = solflx%TRChem_CaSO4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_MgOH2,L,NY,NX)     = solflx%TRChem_MgOH_soil
    trcSalt_RGeoChem_flx_vr(idsalt_MgCO3,L,NY,NX)     = solflx%TRChem_MgCO3_soil
    trcSalt_RGeoChem_flx_vr(idsalt_MgHCO3,L,NY,NX)    = solflx%TRChem_MgHCO3_soil
    trcSalt_RGeoChem_flx_vr(idsalt_MgSO4,L,NY,NX)     = solflx%TRChem_MgSO4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_NaCO3,L,NY,NX)     = solflx%TRChem_NaCO3_soil
    trcSalt_RGeoChem_flx_vr(idsalt_NaSO4,L,NY,NX)     = solflx%TRChem_NaSO4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_KSO4,L,NY,NX)      = solflx%TRChem_KSO4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_H0PO4,L,NY,NX)     = solflx%TRChem_PO4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_H3PO4,L,NY,NX)     = solflx%TRChem_H3PO4_sorbed_soil
    trcSalt_RGeoChem_flx_vr(idsalt_FeHPO4,L,NY,NX)    = solflx%TRChem_FeHPO4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_FeH2PO4,L,NY,NX)   = solflx%TRChem_FeH2PO4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_CaPO4,L,NY,NX)     = solflx%TRChem_CaPO4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_CaHPO4,L,NY,NX)    = solflx%TRChem_CaHPO4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_CaH4P2O8,L,NY,NX)  = solflx%TRChem_CaH4P2O8_soil
    trcSalt_RGeoChem_flx_vr(idsalt_MgHPO4,L,NY,NX)    = solflx%TRChem_MgHPO4_soil
    trcSalt_RGeoChem_flx_vr(idsalt_H0PO4B,L,NY,NX)    = solflx%TRChem_PO4_band_soil
    trcSalt_RGeoChem_flx_vr(idsalt_H3PO4B,L,NY,NX)    = solflx%TRChem_H3PO4_band_soil
    trcSalt_RGeoChem_flx_vr(idsalt_FeHPO4B,L,NY,NX)   = solflx%TRChem_FeHPO4_band_soil
    trcSalt_RGeoChem_flx_vr(idsalt_FeH2PO4B,L,NY,NX)  = solflx%TRChem_FeH2PO4_band_soil
    trcSalt_RGeoChem_flx_vr(idsalt_CaPO4B,L,NY,NX)    = solflx%TRChem_CaPO4_band_soil
    trcSalt_RGeoChem_flx_vr(idsalt_CaHPO4B,L,NY,NX)   = solflx%TRChem_CaHPO4_band_soil
    trcSalt_RGeoChem_flx_vr(idsalt_CaH4P2O8B,L,NY,NX) = solflx%TRChem_CaH4P2O8_band_soil
    trcSalt_RGeoChem_flx_vr(idsalt_MgHPO4B,L,NY,NX)   = solflx%TRChem_MgHPO4_band_soil
  endif
  TRChem_H_p_sorbed_soil_vr(L,NY,NX)     = solflx%TRChem_H_p_sorbed_soil
  TRChem_Al_sorbed_soil_vr(L,NY,NX)      = solflx%TRChem_Al_sorbed_soil
  TRChem_Fe_sorbed_soil_vr(L,NY,NX)      = solflx%TRChem_Fe_sorbed_soil
  TRChem_Ca_sorbed_soil_vr(L,NY,NX)      = solflx%TRChem_Ca_sorbed_soil
  TRChem_Mg_sorbed_soil_vr(L,NY,NX)      = solflx%TRChem_Mg_sorbed_soil
  TRChem_Na_sorbed_soil_vr(L,NY,NX)      = solflx%TRChem_Na_sorbed_soil
  TRChem_K_sorbed_soil_vr(L,NY,NX)       = solflx%TRChem_K_sorbed_soil
  TRChem_HCO3_sorbed_soil_vr(L,NY,NX)    = solflx%TRChem_HCO3_sorbed_soil
  TRChem_AlO2H2_sorbed_soil_vr(L,NY,NX)  = solflx%TRChem_AlO2H2_sorbed_soil
  TRChem_FeO2H2_sorbed_soil_vr(L,NY,NX)  = solflx%TRChem_FeO2H2_sorbed_soil
  trcx_TRSoilChem_vr(idx_OHe,L,NY,NX)    = solflx%TRChem_RO_sorbed_soil
  trcx_TRSoilChem_vr(idx_OHeB,L,NY,NX)   = solflx%TRChem_RO_sorbed_band_soil
  trcp_RChem_soil_vr(idsp_AlOH3,L,NY,NX) = solflx%TRChem_AlOH3_precip_soil
  trcp_RChem_soil_vr(idsp_FeOH3,L,NY,NX) = solflx%TRChem_FeOH3_precip_soil
  trcp_RChem_soil_vr(idsp_CaCO3,L,NY,NX) = solflx%TRChem_CaCO3_precip_soil
  trcp_RChem_soil_vr(idsp_CaSO4,L,NY,NX) = solflx%TRChem_CaSO4_precip_soil
  TRChem_H2O_vr(L,NY,NX)                 = solflx%TRH2O_soil
  TBION_vr(L,NY,NX)                      = solflx%TBION_soil
  Txchem_CO2_vr(L,NY,NX)                 = solflx%TRAquaChem_CO2_soil

  call PrintInfo('end '//subname)
  end subroutine GeochemAPIRecv

end module GeochemAPI
