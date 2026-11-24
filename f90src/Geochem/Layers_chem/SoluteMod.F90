module SoluteMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod  , only : AZMAX1,AZMIN1
  use DebugToolMod
  use SoluteParMod
  use SOMDataType
  use EcosimConst
  use ChemTranspDataType
  use FertilizerDataType
  use GridConsts
  use FlagDataType
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SoilBGCDataType
  use SurfLitterDataType
  use RootDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use PlantDataRateType
  use GridDataType
  use EcoSIMSolverPar
  use ChemEquilibriaMod
  use SaltChemEquilibriaMod
  use SoluteChemDataType
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  public :: UreaHydrolysis
  public :: UpdateSoilFertlizer
  public :: UpdateFertilizerBand
  public :: GeoChemEquilibria
  public :: UpdateSurfResidueSolute
  contains


!------------------------------------------------------------------------

  subroutine GeoChemEquilibria(I,J,L,NY,NX,chemvar, solflx)

  implicit none
  type(chem_var_type), intent(inout) :: chemvar
  type(solute_flx_type), intent(inout) :: solflx
  integer, intent(in) :: I,J,L,NY,NX

  IF(salt_model)THEN
    call SaltChemEquilibria(chemvar,solflx)
  ELSE
    call NoSaltChemEquilibria(I,J,L,NY,NX,chemvar,solflx)
  ENDIF

  end subroutine GeoChemEquilibria


!------------------------------------------------------------------------

  subroutine UpdateFertilizerBand(L,NY,NX)

  implicit none
  integer, intent(in) :: L,NY,NX
  character(len=*), parameter :: subname='UpdateFertilizerBand'
  real(r8) :: FLWD
!

  call PrintInfo('beg '//subname)
!     CHANGE IN WIDTHS AND DEPTHS OF FERTILIZER BANDS FROM
!     VERTICAL AND HORIZONTAL DIFFUSION DRIVEN BY CONCENTRATION
!     DIFFERENCES BETWEEN BAND AND NON-BAND SOIL ZONES
!
!     ROWI=band row width
!     FLWD=net vertical flow relative to area
!
!     IF(ROWI(I,NY,NX).GT.0.0)THEN
  FLWD=0.5_r8*(WaterFlowSoiMicP_3D(3,L,NY,NX)+WaterFlowSoiMicP_3D(3,L+1,NY,NX))/AREA_3D(3,L,NY,NX)
!
!     NH4 FERTILIZER BAND
!
  call UpdateNH3FertilizerBandinfo(L,NY,NX,FLWD)
!
!     NO3 FERTILIZER BAND
!
  call UpdateNO3FertilizerBandinfo(L,NY,NX,FLWD)
!
!     PO4 FERTILIZER BAND
!
  call UpdatePO4FertilizerBandinfo(L,NY,NX,FLWD)
!     ENDIF
!
!     SUBTRACT FERTILIZER DISSOLUTION FROM FERTILIZER POOLS
!
!     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=broadcast NH4,NH3,urea,NO3 fertilizer
!     ZNH4FB,ZNH3FB,ZNHUFB,ZNO3FB=banded NH4,NH3,urea,NO3 fertilizer
!     RFertReleaz_NH4,RFertReleaz_NH4Soil2Band=rate of broadcast NH4 fertilizer dissoln in non-band,band
!     RFertReleaz_NH3,RFertReleaz_NH3Soil2Band=rate of broadcast NH3 fertilizer dissoln in non-band,band
!     RFertReleaz_Urea2SoilNH3,RFertReleaz_Urea2BandNH3=rate of broadcast urea fertr dissoln in non-band,band
!     RFertReleaz_NO3,RFertReleaz_NO3Soil2Band=rate of broadcast NO3 fertilizer dissoln in non-band,band
!     RFertReleaz_NH4Band2Band=rate of banded NH4 fertilizer dissolution in band
!     RFertReleaz_NH3Band2Band=rate of banded NH3 fertilizer dissolution in band
!     RFertReleaz_UreaBand2NH3Band=rate of banded urea fertilizer dissolution in band
!     RFertReleaz_NO3Band2Band=rate of banded NO3 fertilizer dissolution in band
!
  
!
!     ADD FERTILIZER DISSOLUTION TO ION FLUXES AND CONVERT TO MASS
!
!     TRN3G=NH3 dissolution from NH3
!     TRChem_NH4_soil,TRChem_NH4_band_soil_mole=NH4 dissolution in non-band,band
!     TRChem_NH3_soil_mole,TRChem_NH3_band_soil=NH3 dissolution from urea in non-band,band
!     TRNO3,TRNOB=NO3 dissolution in non-band,band
!
! convert the unit from [mol d-2 h-1] to [g d-2 h-1]
!------------------------------------------------------------------------
  TProd_CO2_geochem_soil_vr(L,NY,NX)        = TProd_CO2_geochem_soil_vr(L,NY,NX)*catomw
  TProd_gas_NH3_geochem_vr(L,NY,NX)         = TProd_gas_NH3_geochem_vr(L,NY,NX)*natomw
  TRProd_chem_sol_NH3_soil_vr(L,NY,NX)      = TRProd_chem_sol_NH3_soil_vr(L,NY,NX)*natomw

  trcn_RprodChem_soil_vr(ids_NH4,L,NY,NX)   = trcn_RprodChem_soil_vr(ids_NH4,L,NY,NX)*natomw
  trcn_RprodChem_soil_vr(ids_NO3,L,NY,NX)   = trcn_RprodChem_soil_vr(ids_NO3,L,NY,NX)*natomw
  trcn_RprodChem_soil_vr(ids_NO2,L,NY,NX)   = trcn_RprodChem_soil_vr(ids_NO2,L,NY,NX)*natomw
  trcn_RprodChem_soil_vr(ids_H1PO4,L,NY,NX) = trcn_RprodChem_soil_vr(ids_H1PO4,L,NY,NX)*patomw
  trcn_RprodChem_soil_vr(ids_H2PO4,L,NY,NX) = trcn_RprodChem_soil_vr(ids_H2PO4,L,NY,NX)*patomw

  trcn_RProdChem_band_soil_vr(ids_NH4B,L,NY,NX)   = trcn_RProdChem_band_soil_vr(ids_NH4B,L,NY,NX)*natomw
  trcn_RProdChem_band_soil_vr(idg_NH3B,L,NY,NX)   = trcn_RProdChem_band_soil_vr(idg_NH3B,L,NY,NX)*natomw
  trcn_RProdChem_band_soil_vr(ids_NO3B,L,NY,NX)   = trcn_RProdChem_band_soil_vr(ids_NO3B,L,NY,NX)*natomw
  trcn_RProdChem_band_soil_vr(ids_NO2B,L,NY,NX)   = trcn_RProdChem_band_soil_vr(ids_NO2B,L,NY,NX)*natomw
  trcn_RProdChem_band_soil_vr(ids_H1PO4B,L,NY,NX) = trcn_RProdChem_band_soil_vr(ids_H1PO4B,L,NY,NX)*patomw
  trcn_RProdChem_band_soil_vr(ids_H2PO4B,L,NY,NX) = trcn_RProdChem_band_soil_vr(ids_H2PO4B,L,NY,NX)*patomw

  call PrintInfo('end '//subname)
  end subroutine UpdateFertilizerBand

!------------------------------------------------------------------------

  subroutine UreaHydrolysis(L,NY,NX,RFertReleaz_Urea,RFertReleaz_UreaBand)

  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(out) :: RFertReleaz_Urea,RFertReleaz_UreaBand
  real(r8) :: CNHUB,COMA
  real(r8) :: DFNSB,Urea_mole_conc,DFNSA,DUKD
!     begin_execution
!
!     UREA HYDROLYSIS IN BAND AND NON-BAND SOIL ZONES
!
!     VWatMicrobAct_vr=biologically active soil water volume from nitro.f
!     COMA=concentration of active biomass
!     TMicHeterActivity_vr=total microbial activity from nitro.f
!     DUKD=Km for urea hydrolysis
!
  IF(VWatMicrobAct_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    COMA=AMIN1(0.1E+06_r8,TMicHeterActivity_vr(L,NY,NX)/VWatMicrobAct_vr(L,NY,NX))
  ELSE
    COMA=0.1E+06_r8
  ENDIF
  DUKD=DUKM*(1.0_r8+COMA/DUKI)
!
!     UREA HYDROLYSIS INHIBITION
!
!     ZNHU0,ZNHUI=initial,current inhibition activity
!     RUreaInhibtorConst=rate constants for decline in urea hydrolysis inhibition
!
  IF(ZNHU0_vr(L,NY,NX).GT.ZEROS(NY,NX).AND.ZNHUI_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
    ZNHUI_vr(L,NY,NX)=ZNHUI_vr(L,NY,NX)-RUreaInhibtorConst(iUreaHydInhibitorType_col(NY,NX))*ZNHUI_vr(L,NY,NX) &
      *AMAX1(RUreaInhibtorConst(iUreaHydInhibitorType_col(NY,NX)),1.0-ZNHUI_vr(L,NY,NX)/ZNHU0_vr(L,NY,NX))
  ELSE
    ZNHUI_vr(L,NY,NX)=0._r8
  ENDIF
!
!     UREA CONCENTRATION AND HYDROLYSIS IN NON-BAND
!
!     ZNHUFA=urea fertilizer in non-band
!     Urea_mole_conc=concentration of urea fertilizer in non-band
!     DFNSA=effect of microbial concentration on urea hydrolysis in non-band
!     RFertReleaz_Urea=rate of urea hydrolysis in non-band
!     RFertUreaSpecHydrol=specific rate constant for urea hydrolysis
!     TSens4MicbGrwoth_vr=temperature effect on microbial activity from nitro.f
!     ZNHUI=current inhibition activity
!
  IF(FertN_mole_soil_vr(ifert_N_urea,L,NY,NX).GT.ZEROS(NY,NX).AND.VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
    Urea_mole_conc=FertN_mole_soil_vr(ifert_N_urea,L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX)
  ELSEIF(VLWatMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    Urea_mole_conc=FertN_mole_soil_vr(ifert_N_urea,L,NY,NX)/VLWatMicP_vr(L,NY,NX)
  ELSE
    Urea_mole_conc=0._r8
  ENDIF
  DFNSA=Urea_mole_conc/(Urea_mole_conc+DUKD)
  RFertReleaz_Urea=AMIN1(FertN_mole_soil_vr(ifert_N_urea,L,NY,NX),RFertUreaSpecHydrol*TMicHeterActivity_vr(L,NY,NX) &
    *DFNSA*TSens4MicbGrwoth_vr(L,NY,NX))*(1.0-ZNHUI_vr(L,NY,NX))
!
!     UREA CONCENTRATION AND HYDROLYSIS IN BAND
!
!     ZNHUFB=urea fertilizer in band
!     CNHUB=concentration of urea fertilizer in band
!     DFNSB=effect of microbial concentration on urea hydrolysis in band
!     RFertReleaz_UreaBand=rate of urea hydrolysis in non-band
!     RFertUreaSpecHydrol=specific rate constant for urea hydrolysis
!     TSens4MicbGrwoth_vr=temperature effect on microbial activity from nitro.f
!
  IF(FertN_mole_Band_vr(ifert_N_urea_band,L,NY,NX).GT.ZEROS(NY,NX).AND.VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
    CNHUB=FertN_mole_Band_vr(ifert_N_urea_band,L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX)
  ELSEIF(VLWatMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    CNHUB=FertN_mole_Band_vr(ifert_N_urea_band,L,NY,NX)/VLWatMicP_vr(L,NY,NX)
  ELSE
    CNHUB=0._r8
  ENDIF
  DFNSB=CNHUB/(CNHUB+DUKD)
  RFertReleaz_UreaBand=AMIN1(FertN_mole_Band_vr(ifert_N_urea_band,L,NY,NX),RFertUreaSpecHydrol*TMicHeterActivity_vr(L,NY,NX) &
    *DFNSB*TSens4MicbGrwoth_vr(L,NY,NX))*(1.0-ZNHUI_vr(L,NY,NX))

  end subroutine UreaHydrolysis
!------------------------------------------------------------------------

  subroutine UpdateSoilFertlizer(I,J,L,NY,NX,chemvar)

  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  type(chem_var_type), target, intent(inout) :: chemvar
  real(r8) :: REcoReleaz_H2PO4BX,REcoReleaz_NH4BX,REcoReleaz_NH3BX,REcoReleaz_H1PO4BX
  real(r8) :: VLWatMicPPX,VLWatMicPNX
  real(r8) :: REcoReleaz_H2PO4X,REcoReleaz_H1PO4X,REcoReleaz_NH3X,REcoReleaz_NH4X
  real(r8) :: RFertReleaz_NH3,RFertReleaz_NH4,RFertReleaz_Urea2SoilNH3 
  real(r8) :: RFertReleaz_Urea2BandNH3,RFertReleaz_NH4Soil2Band 
  real(r8) :: RFertReleaz_NH4Band2Band,RFertReleaz_NH3Soil2Band
  real(r8) :: RFertReleaz_NH3Band2Band,RFertReleaz_NO3,RFertReleaz_NO3Soil2Band
  real(r8) :: RFertReleaz_UreaBand2NH3Band,RFertReleaz_NO3Band2Band
  real(r8) :: RFertReleaz_Urea,RFertReleaz_UreaBand
  real(r8), pointer :: VLWatMicPNH
  real(r8), pointer :: BKVLNH
  real(r8), pointer :: BKVLNB
  real(r8), pointer :: BKVLPO
  real(r8), pointer :: BKVLPB
  real(r8), pointer :: VLWatMicPNB
  real(r8), pointer :: VLWatMicPPO
  real(r8), pointer :: H1PO4_2e_aque_mole_conc
  real(r8), pointer :: H1PO4_2e_band_mole_conc
  real(r8), pointer :: H2PO4_1e_aque_mole_conc
  real(r8), pointer :: H2PO4_1e_band_mole_conc
  real(r8), pointer :: NH3_aque_mole_conc
  real(r8), pointer :: NH3_aqu_band_mole_conc
  real(r8), pointer :: NH4_1p_aque_mole_conc
  real(r8), pointer :: NH4_1p_band_mole_conc
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
  real(r8), pointer :: XHPO4_band_mole_conc
  real(r8), pointer :: XH2PO4_band_mole_conc
  real(r8), pointer :: XROH_band_mole_conc
  real(r8), pointer :: XHPO4_mole_conc
  real(r8), pointer :: XROH2_band_mole_conc
  real(r8), pointer :: XH2PO4_mole_conc
  real(r8), pointer :: XNH4_mole_conc
  real(r8), pointer :: XNH4_band_mole_conc
  real(r8), pointer :: XROH1_mole_conc
  real(r8), pointer :: XROH2_mole_conc
  real(r8), pointer :: XROH1_band_mole_conc
  real(r8), pointer :: XOH_mole_conc
  real(r8), pointer :: VLWatMicPPB

  XOH_mole_conc                 => chemvar%XOH_mole_conc
  XROH1_band_mole_conc          => chemvar%XROH1_band_mole_conc
  XROH2_mole_conc               => chemvar%XROH2_mole_conc
  XROH1_mole_conc               => chemvar%XROH1_mole_conc
  XROH2_band_mole_conc          => chemvar%XROH2_band_mole_conc
  XHPO4_mole_conc               => chemvar%XHPO4_mole_conc
  XROH_band_mole_conc           => chemvar%XROH_band_mole_conc
  XHPO4_band_mole_conc          => chemvar%XHPO4_band_mole_conc
  XH2PO4_band_mole_conc         => chemvar%XH2PO4_band_mole_conc
  H2PO4_1e_band_mole_conc       => chemvar%H2PO4_1e_band_mole_conc
  XNH4_mole_conc                => chemvar%XNH4_mole_conc
  XNH4_band_mole_conc           => chemvar%XNH4_band_mole_conc
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
  H2PO4_1e_aque_mole_conc       => chemvar%H2PO4_1e_aque_mole_conc
  H1PO4_2e_band_mole_conc       => chemvar%H1PO4_2e_band_mole_conc
  H1PO4_2e_aque_mole_conc       => chemvar%H1PO4_2e_aque_mole_conc
  XH2PO4_mole_conc              => chemvar%XH2PO4_mole_conc
  VLWatMicPNH                   => chemvar%VLWatMicPNH
  BKVLNH                        => chemvar%BKVLNH
  BKVLNB                        => chemvar%BKVLNB
  BKVLPO                        => chemvar%BKVLPO
  BKVLPB                        => chemvar%BKVLPB
  VLWatMicPNB                   => chemvar%VLWatMicPNB
  VLWatMicPPB                   => chemvar%VLWatMicPPB
  VLWatMicPPO                   => chemvar%VLWatMicPPO
!     begin_execution

  call UreaHydrolysis(L,NY,NX,RFertReleaz_Urea,RFertReleaz_UreaBand)
!
!     NH4, NH3, UREA, NO3 DISSOLUTION IN BAND AND NON-BAND
!     SOIL ZONES FROM FIRST-ORDER FUNCTIONS OF REMAINING
!     FERTILIZER (NOTE: SUPERPHOSPHATE AND ROCK PHOSPHATE
!     ARE REPRESENTED AS MONOCALCIUM PHOSPHATE AND HYDROXYAPATITE
!     MODELLED IN PHOSPHORUS REACTIONS BELOW)
!
!     RFertReleaz_NH4,RFertReleaz_NH4Soil2Band=rate of broadcast NH4 fertilizer dissoln in non-band,band
!     RFertReleaz_NH3,RFertReleaz_NH3Soil2Band=rate of broadcast NH3 fertilizer dissoln in non-band,band
!     RFertReleaz_Urea2SoilNH3,RFertReleaz_Urea2BandNH3=rate of broadcast urea fertr dissoln in non-band,band
!     RFertReleaz_NO3,RFertReleaz_NO3Soil2Band=rate of broadcast NO3 fertilizer dissoln in non-band,band
!     RFertReleaz_NH4Band2Band= rate of banded NH4 fertilizer dissolution in band
!     RFertReleaz_NH3Band2Band= rate of banded NH3 fertilizer dissolution in band
!     RFertReleaz_UreaBand2NH3Band= rate of banded urea fertilizer dissolution in band
!     RFertReleaz_NO3Band2Band= rate of banded NO3 fertilizer dissolution in band
!     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=broadcast NH4,NH3,urea,NO3 fertilizer
!     ZNH4FB,ZNH3FB,ZNHUFB,ZNO3FB=banded NH4,NH3,urea,NO3 fertilizer
!     VLNH4,VLNHB=fractions of soil volume in NH4 non-band,band
!     VLNO3,VLNOB=fractions of soil volume in N03 non-band,band
!     VLPO4,VLPOB=fractions of soil volume in H2PO4 non-band,band
!     THETW=soil water concentration
!
  RFertReleaz_NH4          = RFertNH4SpecReleaz*FertN_mole_soil_vr(ifert_N_nh4,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)*THETW_vr(L,NY,NX)
  RFertReleaz_NH4Band2Band = RFertNH4SpecReleaz*FertN_mole_Band_vr(ifert_N_nh4_band,L,NY,NX)*THETW_vr(L,NY,NX)
  RFertReleaz_NH4Soil2Band = RFertNH4SpecReleaz*FertN_mole_soil_vr(ifert_N_nh4,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)*THETW_vr(L,NY,NX)

  RFertReleaz_NH3          = RFertNH3SpecReleaz*FertN_mole_soil_vr(ifert_N_nh3,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
  RFertReleaz_NH3Soil2Band = RFertNH3SpecReleaz*FertN_mole_soil_vr(ifert_N_nh3,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
  RFertReleaz_NH3Band2Band = RFertNH3SpecReleaz*FertN_mole_Band_vr(ifert_N_nh3_band,L,NY,NX)

  RFertReleaz_Urea2SoilNH3     = RFertReleaz_Urea*trcs_VLN_vr(ids_NH4,L,NY,NX)*THETW_vr(L,NY,NX)
  RFertReleaz_Urea2BandNH3     = RFertReleaz_Urea*trcs_VLN_vr(ids_NH4B,L,NY,NX)*THETW_vr(L,NY,NX)
  RFertReleaz_UreaBand2NH3Band = RFertReleaz_UreaBand*trcs_VLN_vr(ids_NH4B,L,NY,NX)*THETW_vr(L,NY,NX)

  RFertReleaz_NO3          = RFertNO3SpecReleaz*FertN_mole_soil_vr(ifert_N_no3,L,NY,NX)*trcs_VLN_vr(ids_NO3,L,NY,NX)*THETW_vr(L,NY,NX)
  RFertReleaz_NO3Soil2Band = RFertNO3SpecReleaz*FertN_mole_soil_vr(ifert_N_no3,L,NY,NX)*trcs_VLN_vr(ids_NO3B,L,NY,NX)*THETW_vr(L,NY,NX)
  RFertReleaz_NO3Band2Band = RFertNO3SpecReleaz*FertN_mole_Band_vr(ifert_N_no3_band,L,NY,NX)*THETW_vr(L,NY,NX)
!
!     SOLUBLE AND EXCHANGEABLE NH4 CONCENTRATIONS
!     IN NON-BAND AND BAND SOIL ZONES
!
!     VLWatMicPNH,VLWatMicPNB=water volume in NH4 non-band,band
!     REcoReleaz_NH4X,REcoReleaz_NH3X=NH4,NH3 input from uptake, mineraln, dissoln in non-band
!     REcoReleaz_NH4BX,REcoReleaz_NH3BX=NH4,NH3 input from uptake, mineraln, dissoln in band
!     TRootNH4Uptake_pft,TUPNH3=soil-root exchange of NH4,NH3 in non-band from uptake.f
!     TUPNHB,TUPN3B=soil-root exchange of NH4,NH3 in band from uptake.f
!     RNH4MicbReliz2Soil_vr,RNH4MicbReliz2Band_vr=net change in NH4 in band,non-band from nitro.f
!     NH4_1p_aque_mole_conc,NH3_aque_mole_conc,NH4_1p_band_mole_conc,NH3_aqu_band_mole_conc=total NH4,NH3 concentration in non-band,band
!     XNH4_mole_conc,XNH4_band_mole_conc=adsorbed NH4 concentration in non-band,band
!
  IF(VLWatMicPNH.GT.ZEROS2(NY,NX))THEN
    VLWatMicPNX=natomw*VLWatMicPNH
    REcoReleaz_NH4X       = (RNut_MicbRelease_vr(ids_NH4,L,NY,NX)-trcs_Soil2plant_uptake_vr(ids_NH4,L,NY,NX)+natomw*RFertReleaz_NH4)/VLWatMicPNX
    REcoReleaz_NH3X       = (natomw*RFertReleaz_Urea2SoilNH3)/VLWatMicPNX
    NH4_1p_aque_mole_conc = AZMAX1(trcs_solml_vr(ids_NH4,L,NY,NX)/VLWatMicPNX+REcoReleaz_NH4X)
    NH3_aque_mole_conc    = AZMAX1(trcs_solml_vr(idg_NH3,L,NY,NX)/VLWatMicPNX+REcoReleaz_NH3X)

    XNH4_mole_conc        = AZMAX1(trcx_solml_vr(idx_NH4,L,NY,NX)/BKVLNH)
  ELSE
    REcoReleaz_NH4X       = 0._r8
    REcoReleaz_NH3X       = 0._r8
    NH4_1p_aque_mole_conc = 0._r8
    NH3_aque_mole_conc    = 0._r8
    XNH4_mole_conc        = 0._r8
  ENDIF
  
  IF(VLWatMicPNB.GT.ZEROS2(NY,NX))THEN
    VLWatMicPNX            = natomw*VLWatMicPNB
    REcoReleaz_NH4BX       = (RNut_MicbRelease_vr(ids_NH4B,L,NY,NX)-trcs_Soil2plant_uptake_vr(ids_NH4B,L,NY,NX) &
      +natomw*(RFertReleaz_NH4Soil2Band+RFertReleaz_NH4Band2Band))/VLWatMicPNX
    REcoReleaz_NH3BX       = (natomw*(RFertReleaz_Urea2BandNH3+RFertReleaz_UreaBand2NH3Band))/VLWatMicPNX
    NH4_1p_band_mole_conc  = AZMAX1(trcs_solml_vr(ids_NH4B,L,NY,NX)/VLWatMicPNX+REcoReleaz_NH4BX)
    NH3_aqu_band_mole_conc = AZMAX1(trcs_solml_vr(idg_NH3B,L,NY,NX)/VLWatMicPNX+REcoReleaz_NH3BX)

    XNH4_band_mole_conc    = AZMAX1(trcx_solml_vr(idx_NH4B,L,NY,NX)/BKVLNB)
  ELSE
    REcoReleaz_NH4BX       = 0._r8
    REcoReleaz_NH3BX       = 0._r8
    NH4_1p_band_mole_conc  = 0._r8
    NH3_aqu_band_mole_conc = 0._r8
    XNH4_band_mole_conc    = 0._r8
  ENDIF
!
!     SOLUBLE, EXCHANGEABLE AND PRECIPITATED PO4 CONCENTRATIONS IN
!     NON-BAND AND BAND SOIL ZONES
!
!     VLWatMicPPO,VLWatMicPPB=water volume in H2PO4 non-band,band
!     REcoReleaz_H1PO4X,REcoReleaz_H2PO4X=HPO4,H2PO4 inputs from mineraln, uptake in non-band
!     REcoReleaz_H1PO4BX,REcoReleaz_H2PO4BX=HPO4,H2PO4 inputs from mineraln, uptake in band
!     RH1PO4MicbReliz2Soil_vr,RH1PO4MicbReliz2Band_vr=net change in HPO4 in band,non-band from nitro.f
!     TRootHPO4Uptake_pft,TRootH2PO4Uptake_pft=soil-root exch of HPO4,H2PO4 in non-band from uptake.f
!     TUPH1B,TUPH2B=soil-root exch of HPO4,H2PO4 in band from uptake.f
!     H1PO4_2e_aque_mole_conc,H2PO4_1e_aque_mole_conc=HPO4,H2PO4 concentrations in non-band
!     H1PO4_2e_band_mole_conc,H2PO4_1e_band_mole_conc=HPO4,H2PO4 concentrations in band
!     XROH1_mole_conc,XROH1_mole_conc,XROH2_mole_conc=concn of adsorption sites R-,R-OH,R-OH2 in non-band
!     XOH1B,XROH_band_mole_conc,XROH2_band_mole_conc=concn of adsorption sites R-,R-OH,R-OH2 in band
!     XHPO4_mole_conc,XH2PO4_mole_conc=concentration of adsorbed HPO4,H2PO4 in non-band
!     XROH_band_mole_conc,XH2PO4_band_mole_conc=concentration of adsorbed HPO4,H2PO4 in band
!     Precp_AlPO4_mole_conc,Precp_FePO4_mole_conc=concn of precip AlPO4,FEPO4 in non-band
!     PrecpB_AlPO4_mole_conc,PrecpB_FePO4_mole_conc=concn of precip AlPO4,FEPO4 in band
!     Precp_CaH4P2O8_mole_conc,Precp_CaHPO4_mole_conc,Precp_Ca5P3O12O3H3_mole_conc=concn of precip CaH4P2O8,CaHPO4,apatite in non-band
!     PrecpB_CaH4P2O8_con,PrecpB_CaHPO4_mole_conc,PrecpB_Ca5P3O12O3H3_mole_conc=concn of precip CaH4P2O8,CaHPO4,apatite in band
!
  IF(VLWatMicPPO.GT.ZEROS2(NY,NX))THEN
    VLWatMicPPX             = patomw*VLWatMicPPO
    !compute biochemical release of H1PO4 and H2PO4
    REcoReleaz_H1PO4X            = (RNut_MicbRelease_vr(ids_H1PO4,L,NY,NX)-trcs_Soil2plant_uptake_vr(ids_H1PO4,L,NY,NX))/VLWatMicPPX
    REcoReleaz_H2PO4X            = (RNut_MicbRelease_vr(ids_H2PO4,L,NY,NX)-trcs_Soil2plant_uptake_vr(ids_H2PO4,L,NY,NX))/VLWatMicPPX
    !update aquoues concentration 
    H1PO4_2e_aque_mole_conc      = AZMAX1(trcs_solml_vr(ids_H1PO4,L,NY,NX)/VLWatMicPPX+REcoReleaz_H1PO4X)
    H2PO4_1e_aque_mole_conc      = AZMAX1(trcs_solml_vr(ids_H2PO4,L,NY,NX)/VLWatMicPPX+REcoReleaz_H2PO4X)
    
    XOH_mole_conc                = AZMAX1(trcx_solml_vr(idx_OHe,L,NY,NX))/BKVLPO
    XROH1_mole_conc              = AZMAX1(trcx_solml_vr(idx_OH,L,NY,NX))/BKVLPO
    XROH2_mole_conc              = AZMAX1(trcx_solml_vr(idx_OHp,L,NY,NX))/BKVLPO
    XHPO4_mole_conc              = AZMAX1(trcx_solml_vr(idx_HPO4,L,NY,NX))/BKVLPO
    XH2PO4_mole_conc             = AZMAX1(trcx_solml_vr(idx_H2PO4,L,NY,NX))/BKVLPO

    Precp_AlPO4_mole_conc        = AZMAX1(trcp_saltpml_vr(idsp_AlPO4,L,NY,NX))/BKVLPO
    Precp_FePO4_mole_conc        = AZMAX1(trcp_saltpml_vr(idsp_FePO4,L,NY,NX))/BKVLPO
    Precp_CaH4P2O8_mole_conc     = AZMAX1(trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX))/BKVLPO
    Precp_CaHPO4_mole_conc       = AZMAX1(trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX))/BKVLPO
    Precp_Ca5P3O12O3H3_mole_conc = AZMAX1(trcp_saltpml_vr(idsp_Apatite,L,NY,NX))/BKVLPO
  ELSE
    REcoReleaz_H1PO4X            = 0._r8
    REcoReleaz_H2PO4X            = 0._r8
    H1PO4_2e_aque_mole_conc      = 0._r8
    H2PO4_1e_aque_mole_conc      = 0._r8
    XOH_mole_conc                = 0._r8
    XROH1_mole_conc              = 0._r8
    XROH2_mole_conc              = 0._r8
    XHPO4_mole_conc              = 0._r8
    XH2PO4_mole_conc             = 0._r8
    Precp_AlPO4_mole_conc        = 0._r8
    Precp_FePO4_mole_conc        = 0._r8
    Precp_CaH4P2O8_mole_conc     = 0._r8
    Precp_CaHPO4_mole_conc       = 0._r8
    Precp_Ca5P3O12O3H3_mole_conc = 0._r8
  ENDIF

  IF(VLWatMicPPB.GT.ZEROS2(NY,NX))THEN
    VLWatMicPPX              = patomw*VLWatMicPPB
    !biochemical release
    REcoReleaz_H1PO4BX = (RNut_MicbRelease_vr(ids_H1PO4B,L,NY,NX)-trcs_Soil2plant_uptake_vr(ids_H1PO4B,L,NY,NX))/VLWatMicPPX
    REcoReleaz_H2PO4BX = (RNut_MicbRelease_vr(ids_H2PO4B,L,NY,NX)-trcs_Soil2plant_uptake_vr(ids_H2PO4B,L,NY,NX))/VLWatMicPPX
    !aqueous update
    H1PO4_2e_band_mole_conc       = AZMAX1(trcs_solml_vr(ids_H1PO4B,L,NY,NX)/VLWatMicPPX+REcoReleaz_H1PO4BX)
    H2PO4_1e_band_mole_conc       = AZMAX1(trcs_solml_vr(ids_H2PO4B,L,NY,NX)/VLWatMicPPX+REcoReleaz_H2PO4BX)

    XROH1_band_mole_conc          = AZMAX1(trcx_solml_vr(idx_OHeB,L,NY,NX))/BKVLPB
    XROH_band_mole_conc           = AZMAX1(trcx_solml_vr(idx_OHB,L,NY,NX))/BKVLPB
    XROH2_band_mole_conc          = AZMAX1(trcx_solml_vr(idx_OHpB,L,NY,NX))/BKVLPB
    XHPO4_band_mole_conc          = AZMAX1(trcx_solml_vr(idx_HPO4B,L,NY,NX))/BKVLPB
    XH2PO4_band_mole_conc         = AZMAX1(trcx_solml_vr(idx_H2PO4B,L,NY,NX))/BKVLPB
    
    PrecpB_AlPO4_mole_conc        = AZMAX1(trcp_saltpml_vr(idsp_AlPO4B,L,NY,NX))/BKVLPB
    PrecpB_FePO4_mole_conc        = AZMAX1(trcp_saltpml_vr(idsp_FePO4B,L,NY,NX))/BKVLPB
    PrecpB_CaH4P2O8_mole_conc     = AZMAX1(trcp_saltpml_vr(idsp_CaH4P2O8B,L,NY,NX))/BKVLPB
    PrecpB_CaHPO4_mole_conc       = AZMAX1(trcp_saltpml_vr(idsp_CaHPO4B,L,NY,NX))/BKVLPB
    PrecpB_Ca5P3O12O3H3_mole_conc = AZMAX1(trcp_saltpml_vr(idsp_ApatiteBand,L,NY,NX))/BKVLPB
  ELSE
    REcoReleaz_H1PO4BX        = 0._r8
    REcoReleaz_H2PO4BX        = 0._r8
    H1PO4_2e_band_mole_conc   = 0._r8
    H2PO4_1e_band_mole_conc   = 0._r8
    XROH1_band_mole_conc      = 0._r8
    XROH_band_mole_conc       = 0._r8
    XROH2_band_mole_conc      = 0._r8
    XHPO4_band_mole_conc      = 0._r8
    XH2PO4_band_mole_conc     = 0._r8
    PrecpB_AlPO4_mole_conc    = 0._r8
    PrecpB_FePO4_mole_conc    = 0._r8
    PrecpB_CaH4P2O8_mole_conc = 0._r8
    PrecpB_CaHPO4_mole_conc   = 0._r8
    PrecpB_Ca5P3O12O3H3_mole_conc  = 0._r8
  ENDIF

  FertN_mole_soil_vr(ifert_N_nh4,L,NY,NX)  = AZMAX1(FertN_mole_soil_vr(ifert_N_nh4,L,NY,NX)-RFertReleaz_NH4-RFertReleaz_NH4Soil2Band)
  FertN_mole_soil_vr(ifert_N_nh3,L,NY,NX)  = AZMAX1(FertN_mole_soil_vr(ifert_N_nh3,L,NY,NX)-RFertReleaz_NH3-RFertReleaz_NH3Soil2Band)
  FertN_mole_soil_vr(ifert_N_urea,L,NY,NX) = AZMAX1(FertN_mole_soil_vr(ifert_N_urea,L,NY,NX)-RFertReleaz_Urea2SoilNH3-RFertReleaz_Urea2BandNH3)
  FertN_mole_soil_vr(ifert_N_no3,L,NY,NX)  = AZMAX1(FertN_mole_soil_vr(ifert_N_no3,L,NY,NX)-RFertReleaz_NO3-RFertReleaz_NO3Soil2Band)

  FertN_mole_Band_vr(ifert_N_nh4_band,L,NY,NX)  = AZMAX1(FertN_mole_Band_vr(ifert_N_nh4_band,L,NY,NX)-RFertReleaz_NH4Band2Band)
  FertN_mole_Band_vr(ifert_N_nh3_band,L,NY,NX)  = AZMAX1(FertN_mole_Band_vr(ifert_N_nh3_band,L,NY,NX)-RFertReleaz_NH3Band2Band)
  FertN_mole_Band_vr(ifert_N_urea_band,L,NY,NX) = AZMAX1(FertN_mole_Band_vr(ifert_N_urea_band,L,NY,NX)-RFertReleaz_UreaBand2NH3Band)
  FertN_mole_Band_vr(ifert_N_no3_band,L,NY,NX)  = AZMAX1(FertN_mole_Band_vr(ifert_N_no3_band,L,NY,NX)-RFertReleaz_NO3Band2Band)

  TProd_gas_NH3_geochem_vr(L,NY,NX)       = TProd_gas_NH3_geochem_vr(L,NY,NX)+RFertReleaz_NH3+RFertReleaz_NH3Soil2Band+RFertReleaz_NH3Band2Band
  trcn_RprodChem_soil_vr(ids_NH4,L,NY,NX) = trcn_RprodChem_soil_vr(ids_NH4,L,NY,NX)+RFertReleaz_NH4
  TRProd_chem_sol_NH3_soil_vr(L,NY,NX)    = TRProd_chem_sol_NH3_soil_vr(L,NY,NX)+RFertReleaz_Urea2SoilNH3   !add to dissolved NH3
  trcn_RprodChem_soil_vr(ids_NO3,L,NY,NX) = trcn_RprodChem_soil_vr(ids_NO3,L,NY,NX)+RFertReleaz_NO3

  trcn_RProdChem_band_soil_vr(ids_NH4B,L,NY,NX)   = trcn_RProdChem_band_soil_vr(ids_NH4B,L,NY,NX)+RFertReleaz_NH4Soil2Band+RFertReleaz_NH4Band2Band
  trcn_RProdChem_band_soil_vr(idg_NH3B,L,NY,NX)   = trcn_RProdChem_band_soil_vr(idg_NH3B,L,NY,NX)+RFertReleaz_Urea2BandNH3+RFertReleaz_UreaBand2NH3Band
  trcn_RProdChem_band_soil_vr(ids_NO3B,L,NY,NX)   = trcn_RProdChem_band_soil_vr(ids_NO3B,L,NY,NX)+RFertReleaz_NO3Soil2Band+RFertReleaz_NO3Band2Band

  end subroutine UpdateSoilFertlizer

!------------------------------------------------------------------------------------------

  subroutine UpdateFertilizerBandDepth(L,NY,NX,DepthChange,BandDepth,BandThickness,BandWidth)
  ! Update band geometry in meters; this routine does not move or remove nutrient mass.
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in) :: DepthChange
  real(r8), intent(inout) :: BandDepth,BandThickness(:),BandWidth(:)
  character(len=*), parameter :: subname='UpdateFertilizerBandDepth'
  real(r8) :: MinChange,MaxChange,AppliedChange,Excess

  call PrintInfo('beg '//subname)
  !
  ! The redistribution stencil reaches only one adjacent layer. Limit movement
  ! to its available thickness/capacity so large steps cannot overfill or empty
  ! that neighbor. NU:NL are active soil bounds; indices outside them are not
  ! destinations for band geometry, even when allocated in the arrays.
  MinChange=-BandThickness(L)
  if(L>NU_col(NY,NX))MinChange=MinChange-BandThickness(L-1)
  MaxChange=DLYR_3D(3,L,NY,NX)-BandThickness(L)
  if(L<NL_col(NY,NX))MaxChange=MaxChange &
    +AMAX1(0._r8,DLYR_3D(3,L+1,NY,NX)-BandThickness(L+1))
  AppliedChange=AMAX1(MinChange,AMIN1(MaxChange,DepthChange))
  !
  ! Use the same limited change for total depth and layer thickness.
  BandDepth        = BandDepth+AppliedChange
  BandThickness(L) = BandThickness(L)+AppliedChange
  if(BandThickness(L)>DLYR_3D(3,L,NY,NX))then
    Excess=BandThickness(L)-DLYR_3D(3,L,NY,NX)
    if(L<NL_col(NY,NX))then
      BandThickness(L+1) = BandThickness(L+1)+Excess
      BandWidth(L+1)     = BandWidth(L)
    endif
    BandThickness(L)=DLYR_3D(3,L,NY,NX)
  elseif(BandThickness(L)<=0._r8)then
    if(L>NU_col(NY,NX))BandThickness(L-1)=BandThickness(L-1)+BandThickness(L)
    BandThickness(L) = 0._r8
    BandWidth(L)     = 0._r8
  endif
  call PrintInfo('end '//subname)
  end subroutine UpdateFertilizerBandDepth

!------------------------------------------------------------------------------------------
  subroutine UpdateNH3FertilizerBandinfo(L,NY,NX,FLWD)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in) :: FLWD

  real(r8) :: XVLNH4,DWNH4,DXNH4
  real(r8) :: DPFLW,DNH4S,DNH3S
  real(r8) :: FVLNH4
!     begin_execution
!     iFertNH4Band_col=banded NH4 fertilizer flag
!     ROWN=NH4 fertilizer band row width
!     DPNH4=NH4 fertilizer band depth
!
  IF(iFertNH4Band_col(NY,NX).EQ.1 .AND. ROWSpaceNH4_col(NY,NX).GT.0.0)THEN
    IF(L.EQ.NU_col(NY,NX).OR.CumDepz2LayBottom_vr(L-1,NY,NX).LT.BandDepthNH4_col(NY,NX))THEN
!
!     NH4 BAND WIDTH
!
!     DWNH4=change in NH4 fertilizer band width
!     WDNHB=layer NH4 fertilizer band width
!     ZNSGL=NH4 diffusivity
!     TortMicPM_vr=tortuosity
!
      DWNH4                    = 0.5*SQRT(SoluteDifusvtyT_vr(idg_NH3,L,NY,NX))*TortMicPM_vr(NPH,L,NY,NX)
      BandWidthNH4_vr(L,NY,NX) = AMIN1(ROWSpaceNH4_col(NY,NX),AMAX1(0.025,BandWidthNH4_vr(L,NY,NX))+DWNH4)
!
!     NH4 BAND DEPTH
!
!     DPFLW=change in NH4 fertilizer band depth
!     DPNH4,DPNHB=total,layer NH4 fertilizer band depth
!
      IF(CumDepz2LayBottom_vr(L,NY,NX).GE.BandDepthNH4_col(NY,NX))THEN
        DPFLW=FLWD+DWNH4
        call UpdateFertilizerBandDepth(L,NY,NX,DPFLW,BandDepthNH4_col(NY,NX), &
          BandThicknessNH4_vr(:,NY,NX),BandWidthNH4_vr(:,NY,NX))
      ENDIF
!
!     FRACTION OF SOIL LAYER OCCUPIED BY NH4 BAND
!     FROM BAND WIDTH X DEPTH
!
!     VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
!     DLYR=soil layer thickness
!     FVLNH4=relative change in VLNH4
!
      XVLNH4=trcs_VLN_vr(ids_NH4,L,NY,NX)
      IF(DLYR_3D(3,L,NY,NX).GT.ZERO)THEN
        trcs_VLN_vr(ids_NH4B,L,NY,NX)=AZMAX1(AMIN1(0.999_r8,BandWidthNH4_vr(L,NY,NX) &
          /ROWSpaceNH4_col(NY,NX)*BandThicknessNH4_vr(L,NY,NX)/DLYR_3D(3,L,NY,NX)))
      ELSE
        trcs_VLN_vr(ids_NH4B,L,NY,NX)=0._r8
      ENDIF
      trcs_VLN_vr(ids_NH4,L,NY,NX)  = 1._r8-trcs_VLN_vr(ids_NH4B,L,NY,NX)
      trcs_VLN_vr(idg_NH3,L,NY,NX)  = trcs_VLN_vr(ids_NH4,L,NY,NX)
      trcs_VLN_vr(idg_NH3B,L,NY,NX) = trcs_VLN_vr(ids_NH4B,L,NY,NX)
      FVLNH4                        = AZMIN1((trcs_VLN_vr(ids_NH4,L,NY,NX)-XVLNH4)/XVLNH4)
!
!     TRANSFER NH4, NH3 FROM NON-BAND TO BAND
!     DURING BAND GROWTH
!
!     DNH4S,DNH3S,DXNH4=transfer of NH4,NH3,exchangeable NH4
!
      DNH4S                                     = FVLNH4*trcs_solml_vr(ids_NH4,L,NY,NX)/natomw
      DNH3S                                     = FVLNH4*trcs_solml_vr(idg_NH3,L,NY,NX)/natomw
      DXNH4                                     = FVLNH4*trcx_solml_vr(idx_NH4,L,NY,NX)
      trcn_RprodChem_soil_vr(ids_NH4,L,NY,NX)     = trcn_RprodChem_soil_vr(ids_NH4,L,NY,NX)+DNH4S
      trcn_RProdChem_band_soil_vr(ids_NH4B,L,NY,NX) = trcn_RProdChem_band_soil_vr(ids_NH4B,L,NY,NX)-DNH4S
      TRProd_chem_sol_NH3_soil_vr(L,NY,NX)           = TRProd_chem_sol_NH3_soil_vr(L,NY,NX)+DNH3S
      trcn_RProdChem_band_soil_vr(idg_NH3B,L,NY,NX) = trcn_RProdChem_band_soil_vr(idg_NH3B,L,NY,NX)-DNH3S
      trcx_TRSoilChem_vr(idx_NH4,L,NY,NX)       = trcx_TRSoilChem_vr(idx_NH4,L,NY,NX)+DXNH4
      trcx_TRSoilChem_vr(idx_NH4B,L,NY,NX)      = trcx_TRSoilChem_vr(idx_NH4B,L,NY,NX)-DXNH4
    ELSE
!
!     AMALGAMATE NH4 BAND WITH NON-BAND IF BAND NO LONGER EXISTS
!
      BandThicknessNH4_vr(L,NY,NX)  = 0._r8
      BandWidthNH4_vr(L,NY,NX)      = 0._r8
      trcs_VLN_vr(ids_NH4,L,NY,NX)  = 1._r8
      trcs_VLN_vr(ids_NH4B,L,NY,NX) = 0._r8
      trcs_VLN_vr(idg_NH3,L,NY,NX)  = trcs_VLN_vr(ids_NH4,L,NY,NX)
      trcs_VLN_vr(idg_NH3B,L,NY,NX) = trcs_VLN_vr(ids_NH4B,L,NY,NX)

      trcs_solml_vr(ids_NH4,L,NY,NX)   = trcs_solml_vr(ids_NH4,L,NY,NX)+trcs_solml_vr(ids_NH4B,L,NY,NX)
      trcs_solml_vr(idg_NH3,L,NY,NX)   = trcs_solml_vr(idg_NH3,L,NY,NX)+trcs_solml_vr(idg_NH3B,L,NY,NX)
      trcs_solml_vr(ids_NH4B,L,NY,NX)  = 0._r8
      trcs_solml_vr(idg_NH3B,L,NY,NX)  = 0._r8
      
      trcx_solml_vr(idx_NH4,L,NY,NX)  = trcx_solml_vr(idx_NH4,L,NY,NX)+trcx_solml_vr(idx_NH4B,L,NY,NX)
      trcx_solml_vr(idx_NH4B,L,NY,NX) = 0._r8
    ENDIF
  ENDIF
  
  end subroutine UpdateNH3FertilizerBandinfo
!------------------------------------------------------------------------------------------

  subroutine UpdatePO4FertilizerBandinfo(L,NY,NX,FLWD)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in) :: FLWD
  real(r8) :: XVLPO4,DPFLW
  real(r8) :: DWPO4,DZH0P,DZH1P
  real(r8) :: DZH2P,DZH3P,DZF1P
  real(r8) :: DZF2P,DZC0P
  real(r8) :: DZC1P,DZC2P,DZM1P,DXOH0,DXOH1,DXOH2,DXH1P,DXH2P
  real(r8) :: DPALP,DPFEP,DPCDP,DPCHP,DPCMP,FVLPO4
!     begin_execution
!     iFertPO4Band_col=banded H2PO4 fertilizer flag
!     ROWP=H2PO4 fertilizer band row width
!     DPPO4=H2PO4 fertilizer band depth
!
  IF(iFertPO4Band_col(NY,NX).EQ.ifert_on .AND. ROWSpacePO4_col(NY,NX).GT.0.0_r8)THEN
    IF(L.EQ.NU_col(NY,NX) .OR. CumDepz2LayBottom_vr(L-1,NY,NX).LT.BandDepthPO4_col(NY,NX))THEN
!
!     PO4 BAND WIDTH
!     DWPO4=change in H2PO4 fertilizer band width
!     WDPO4=layer H2PO4 fertilizer band width
!     POSGL=H2PO4 diffusivity
!     TortMicPM_vr=tortuosity
!
      DWPO4=0.5_r8*SQRT(SoluteDifusvtyT_vr(ids_H1PO4,L,NY,NX))*TortMicPM_vr(NPH,L,NY,NX)
      BandWidthPO4_vr(L,NY,NX)=AMIN1(ROWSpacePO4_col(NY,NX),BandWidthPO4_vr(L,NY,NX)+DWPO4)
!
!     PO4 BAND DEPTH
!
!     DPFLW=change in H2PO4 fertilizer band depth
!     DPPO4,DPPOB=total,layer H2PO4 fertilizer band depth
!
      IF(CumDepz2LayBottom_vr(L,NY,NX).GE.BandDepthPO4_col(NY,NX))THEN
        DPFLW=FLWD+DWPO4
        call UpdateFertilizerBandDepth(L,NY,NX,DPFLW,BandDepthPO4_col(NY,NX), &
          BandThicknessPO4_vr(:,NY,NX),BandWidthPO4_vr(:,NY,NX))
      ENDIF
!
!     FRACTION OF SOIL LAYER OCCUPIED BY PO4 BAND
!     FROM BAND WIDTH X DEPTH
!
!     VLPO43,VLPOB=fraction of soil volume in H2PO4 non-band,band
!     DLYR=soil layer thickness
!     FVLPO4=relative change in VLPO4
!
      XVLPO4=trcs_VLN_vr(ids_H1PO4,L,NY,NX)
      IF(DLYR_3D(3,L,NY,NX).GT.ZERO)THEN
        trcs_VLN_vr(ids_H1PO4B,L,NY,NX)=AZMAX1(AMIN1(0.999,BandWidthPO4_vr(L,NY,NX) &
          /ROWSpacePO4_col(NY,NX)*BandThicknessPO4_vr(L,NY,NX)/DLYR_3D(3,L,NY,NX)))
      ELSE
        trcs_VLN_vr(ids_H1PO4B,L,NY,NX)=0._r8
      ENDIF
      trcs_VLN_vr(ids_H1PO4,L,NY,NX) = 1._r8-trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
      trcs_VLN_vr(ids_H2PO4B,L,NY,NX) = trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
      trcs_VLN_vr(ids_H2PO4,L,NY,NX)  = trcs_VLN_vr(ids_H1PO4,L,NY,NX)
      FVLPO4                          = AZMIN1((trcs_VLN_vr(ids_H1PO4,L,NY,NX)-XVLPO4)/XVLPO4)
!
!     TRANSFER HPO4,H2PO4 FROM NON-BAND TO BAND
!     DURING BAND GROWTH DEPENDING ON SALT
!     VS. NON-SALT OPTION
!
!     DZ*,DX*,DP*=transfer of solute,adsorbed,precipitated HPO4,H2PO4
!
      IF(salt_model)THEN
        DZH1P = FVLPO4*trcs_solml_vr(ids_H1PO4,L,NY,NX)/patomw
        DZH2P = FVLPO4*trcs_solml_vr(ids_H2PO4,L,NY,NX)/patomw

        DZH0P = FVLPO4*trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)
        DZH3P = FVLPO4*trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX)
        DZF1P = FVLPO4*trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)
        DZF2P = FVLPO4*trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)
        DZC0P = FVLPO4*trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX)
        DZC1P = FVLPO4*trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX)
        DZC2P = FVLPO4*trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)
        DZM1P = FVLPO4*trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX)
        DXOH0 = FVLPO4*trcx_solml_vr(idx_OHe,L,NY,NX)
        DXOH1 = FVLPO4*trcx_solml_vr(idx_OH,L,NY,NX)
        DXOH2 = FVLPO4*trcx_solml_vr(idx_OHp,L,NY,NX)
        DXH1P = FVLPO4*trcx_solml_vr(idx_HPO4,L,NY,NX)
        DXH2P = FVLPO4*trcx_solml_vr(idx_H2PO4,L,NY,NX)
        DPALP = FVLPO4*trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)
        DPFEP = FVLPO4*trcp_saltpml_vr(idsp_FePO4,L,NY,NX)
        DPCDP = FVLPO4*trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)
        DPCHP = FVLPO4*trcp_saltpml_vr(idsp_Apatite,L,NY,NX)
        DPCMP = FVLPO4*trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)

        trcn_RProdChem_band_soil_vr(ids_H1PO4B,L,NY,NX) = trcn_RProdChem_band_soil_vr(ids_H1PO4B,L,NY,NX)-DZH1P
        trcn_RProdChem_band_soil_vr(ids_H2PO4B,L,NY,NX) = trcn_RProdChem_band_soil_vr(ids_H2PO4B,L,NY,NX)-DZH2P

        trcSalt_RGeoChem_flx_vr(idsalt_H0PO4,L,NY,NX)    = trcSalt_RGeoChem_flx_vr(idsalt_H0PO4,L,NY,NX)+DZH0P
        trcn_RprodChem_soil_vr(ids_H1PO4,L,NY,NX)        = trcn_RprodChem_soil_vr(ids_H1PO4,L,NY,NX)+DZH1P
        trcn_RprodChem_soil_vr(ids_H2PO4,L,NY,NX)        = trcn_RprodChem_soil_vr(ids_H2PO4,L,NY,NX)+DZH2P
        trcSalt_RGeoChem_flx_vr(idsalt_H3PO4,L,NY,NX)    = trcSalt_RGeoChem_flx_vr(idsalt_H3PO4,L,NY,NX)+DZH3P
        trcSalt_RGeoChem_flx_vr(idsalt_FeHPO4,L,NY,NX)   = trcSalt_RGeoChem_flx_vr(idsalt_FeHPO4,L,NY,NX)+DZF1P
        trcSalt_RGeoChem_flx_vr(idsalt_FeH2PO4,L,NY,NX)  = trcSalt_RGeoChem_flx_vr(idsalt_FeH2PO4,L,NY,NX)+DZF2P
        trcSalt_RGeoChem_flx_vr(idsalt_CaPO4,L,NY,NX)    = trcSalt_RGeoChem_flx_vr(idsalt_CaPO4,L,NY,NX)+DZC0P
        trcSalt_RGeoChem_flx_vr(idsalt_CaHPO4,L,NY,NX)   = trcSalt_RGeoChem_flx_vr(idsalt_CaHPO4,L,NY,NX)+DZC1P
        trcSalt_RGeoChem_flx_vr(idsalt_CaH4P2O8,L,NY,NX) = trcSalt_RGeoChem_flx_vr(idsalt_CaH4P2O8,L,NY,NX)+DZC2P
        trcSalt_RGeoChem_flx_vr(idsalt_MgHPO4,L,NY,NX)   = trcSalt_RGeoChem_flx_vr(idsalt_MgHPO4,L,NY,NX)+DZM1P
        trcSalt_RGeoChem_flx_vr(idsalt_H0PO4B,L,NY,NX)   = trcSalt_RGeoChem_flx_vr(idsalt_H0PO4B,L,NY,NX)-DZH0P

        trcSalt_RGeoChem_flx_vr(idsalt_H3PO4B,L,NY,NX)    = trcSalt_RGeoChem_flx_vr(idsalt_H3PO4B,L,NY,NX)-DZH3P
        trcSalt_RGeoChem_flx_vr(idsalt_FeHPO4B,L,NY,NX)   = trcSalt_RGeoChem_flx_vr(idsalt_FeHPO4B,L,NY,NX)-DZF1P
        trcSalt_RGeoChem_flx_vr(idsalt_FeH2PO4B,L,NY,NX)  = trcSalt_RGeoChem_flx_vr(idsalt_FeH2PO4B,L,NY,NX)-DZF2P
        trcSalt_RGeoChem_flx_vr(idsalt_CaPO4B,L,NY,NX)    = trcSalt_RGeoChem_flx_vr(idsalt_CaPO4B,L,NY,NX)-DZC0P
        trcSalt_RGeoChem_flx_vr(idsalt_CaHPO4B,L,NY,NX)   = trcSalt_RGeoChem_flx_vr(idsalt_CaHPO4B,L,NY,NX)-DZC1P
        trcSalt_RGeoChem_flx_vr(idsalt_CaH4P2O8B,L,NY,NX) = trcSalt_RGeoChem_flx_vr(idsalt_CaH4P2O8B,L,NY,NX)-DZC2P
        trcSalt_RGeoChem_flx_vr(idsalt_MgHPO4B,L,NY,NX)   = trcSalt_RGeoChem_flx_vr(idsalt_MgHPO4B,L,NY,NX)-DZM1P
        
        trcx_TRSoilChem_vr(idx_OHe,L,NY,NX)        = trcx_TRSoilChem_vr(idx_OHe,L,NY,NX)+DXOH0
        trcx_TRSoilChem_vr(idx_OH,L,NY,NX)         = trcx_TRSoilChem_vr(idx_OH,L,NY,NX)+DXOH1
        trcx_TRSoilChem_vr(idx_OHp,L,NY,NX)        = trcx_TRSoilChem_vr(idx_OHp,L,NY,NX)+DXOH2
        trcx_TRSoilChem_vr(idx_HPO4,L,NY,NX)       = trcx_TRSoilChem_vr(idx_HPO4,L,NY,NX)+DXH1P
        trcx_TRSoilChem_vr(idx_H2PO4,L,NY,NX)      = trcx_TRSoilChem_vr(idx_H2PO4,L,NY,NX)+DXH2P
        trcx_TRSoilChem_vr(idx_OHeB,L,NY,NX)       = trcx_TRSoilChem_vr(idx_OHeB,L,NY,NX)-DXOH0
        trcx_TRSoilChem_vr(idx_OHB,L,NY,NX)        = trcx_TRSoilChem_vr(idx_OHB,L,NY,NX)-DXOH1
        trcx_TRSoilChem_vr(idx_OHpB,L,NY,NX)       = trcx_TRSoilChem_vr(idx_OHpB,L,NY,NX)-DXOH2
        trcx_TRSoilChem_vr(idx_HPO4B,L,NY,NX)      = trcx_TRSoilChem_vr(idx_HPO4B,L,NY,NX)-DXH1P
        trcx_TRSoilChem_vr(idx_H2PO4B,L,NY,NX)     = trcx_TRSoilChem_vr(idx_H2PO4B,L,NY,NX)-DXH2P
        trcp_RChem_soil_vr(idsp_AlPO4,L,NY,NX)     = trcp_RChem_soil_vr(idsp_AlPO4,L,NY,NX)+DPALP
        trcp_RChem_soil_vr(idsp_FePO4,L,NY,NX)     = trcp_RChem_soil_vr(idsp_FePO4,L,NY,NX)+DPFEP
        trcp_RChem_soil_vr(idsp_CaHPO4,L,NY,NX)    = trcp_RChem_soil_vr(idsp_CaHPO4,L,NY,NX)+DPCDP
        trcp_RChem_soil_vr(idsp_Apatite,L,NY,NX)        = trcp_RChem_soil_vr(idsp_Apatite,L,NY,NX)+DPCHP
        trcp_RChem_soil_vr(idsp_CaH4P2O8,L,NY,NX)  = trcp_RChem_soil_vr(idsp_CaH4P2O8,L,NY,NX)+DPCMP
        trcp_RChem_soil_vr(idsp_AlPO4B,L,NY,NX)    = trcp_RChem_soil_vr(idsp_AlPO4B,L,NY,NX)-DPALP
        trcp_RChem_soil_vr(idsp_FePO4B,L,NY,NX)    = trcp_RChem_soil_vr(idsp_FePO4B,L,NY,NX)-DPFEP
        trcp_RChem_soil_vr(idsp_CaHPO4B,L,NY,NX)   = trcp_RChem_soil_vr(idsp_CaHPO4B,L,NY,NX)-DPCDP
        trcp_RChem_soil_vr(idsp_ApatiteBand,L,NY,NX)       = trcp_RChem_soil_vr(idsp_ApatiteBand,L,NY,NX)-DPCHP
        trcp_RChem_soil_vr(idsp_CaH4P2O8B,L,NY,NX) = trcp_RChem_soil_vr(idsp_CaH4P2O8B,L,NY,NX)-DPCMP
      ELSE
        DZH1P                                       = FVLPO4*trcs_solml_vr(ids_H1PO4,L,NY,NX)/patomw
        DZH2P                                       = FVLPO4*trcs_solml_vr(ids_H2PO4,L,NY,NX)/patomw
        DXOH1                                       = FVLPO4*trcx_solml_vr(idx_OH,L,NY,NX)
        DXOH2                                       = FVLPO4*trcx_solml_vr(idx_OHp,L,NY,NX)
        DXH2P                                       = FVLPO4*trcx_solml_vr(idx_H2PO4,L,NY,NX)
        DPALP                                       = FVLPO4*trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)
        DPFEP                                       = FVLPO4*trcp_saltpml_vr(idsp_FePO4,L,NY,NX)
        DPCDP                                       = FVLPO4*trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)
        DPCHP                                       = FVLPO4*trcp_saltpml_vr(idsp_Apatite,L,NY,NX)
        DPCMP                                       = FVLPO4*trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)
        trcn_RprodChem_soil_vr(ids_H1PO4,L,NY,NX)     = trcn_RprodChem_soil_vr(ids_H1PO4,L,NY,NX)+DZH1P
        trcn_RprodChem_soil_vr(ids_H2PO4,L,NY,NX)     = trcn_RprodChem_soil_vr(ids_H2PO4,L,NY,NX)+DZH2P
        trcx_TRSoilChem_vr(idx_OH,L,NY,NX)          = trcx_TRSoilChem_vr(idx_OH,L,NY,NX)+DXOH1
        trcx_TRSoilChem_vr(idx_OHp,L,NY,NX)         = trcx_TRSoilChem_vr(idx_OHp,L,NY,NX)+DXOH2
        trcx_TRSoilChem_vr(idx_H2PO4,L,NY,NX)       = trcx_TRSoilChem_vr(idx_H2PO4,L,NY,NX)+DXH2P
        trcn_RProdChem_band_soil_vr(ids_H1PO4B,L,NY,NX) = trcn_RProdChem_band_soil_vr(ids_H1PO4B,L,NY,NX)-DZH1P
        trcn_RProdChem_band_soil_vr(ids_H2PO4B,L,NY,NX) = trcn_RProdChem_band_soil_vr(ids_H2PO4B,L,NY,NX)-DZH2P
        trcx_TRSoilChem_vr(idx_OHB,L,NY,NX)         = trcx_TRSoilChem_vr(idx_OHB,L,NY,NX)-DXOH1
        trcx_TRSoilChem_vr(idx_OHpB,L,NY,NX)        = trcx_TRSoilChem_vr(idx_OHpB,L,NY,NX)-DXOH2
        trcx_TRSoilChem_vr(idx_H2PO4B,L,NY,NX)      = trcx_TRSoilChem_vr(idx_H2PO4B,L,NY,NX)-DXH2P
        trcp_RChem_soil_vr(idsp_AlPO4,L,NY,NX)      = trcp_RChem_soil_vr(idsp_AlPO4,L,NY,NX)+DPALP
        trcp_RChem_soil_vr(idsp_FePO4,L,NY,NX)      = trcp_RChem_soil_vr(idsp_FePO4,L,NY,NX)+DPFEP
        trcp_RChem_soil_vr(idsp_CaHPO4,L,NY,NX)     = trcp_RChem_soil_vr(idsp_CaHPO4,L,NY,NX)+DPCDP
        trcp_RChem_soil_vr(idsp_Apatite,L,NY,NX)         = trcp_RChem_soil_vr(idsp_Apatite,L,NY,NX)+DPCHP
        trcp_RChem_soil_vr(idsp_CaH4P2O8,L,NY,NX)   = trcp_RChem_soil_vr(idsp_CaH4P2O8,L,NY,NX)+DPCMP
        trcp_RChem_soil_vr(idsp_AlPO4B,L,NY,NX)     = trcp_RChem_soil_vr(idsp_AlPO4B,L,NY,NX)-DPALP
        trcp_RChem_soil_vr(idsp_FePO4B,L,NY,NX)     = trcp_RChem_soil_vr(idsp_FePO4B,L,NY,NX)-DPFEP
        trcp_RChem_soil_vr(idsp_CaHPO4B,L,NY,NX)    = trcp_RChem_soil_vr(idsp_CaHPO4B,L,NY,NX)-DPCDP
        trcp_RChem_soil_vr(idsp_ApatiteBand,L,NY,NX)        = trcp_RChem_soil_vr(idsp_ApatiteBand,L,NY,NX)-DPCHP
        trcp_RChem_soil_vr(idsp_CaH4P2O8B,L,NY,NX)  = trcp_RChem_soil_vr(idsp_CaH4P2O8B,L,NY,NX)-DPCMP
      ENDIF
    ELSE
!
!     AMALGAMATE PO4 BAND WITH NON-BAND IF BAND NO LONGER EXISTS
!
      BandThicknessPO4_vr(L,NY,NX)    = 0._r8
      BandWidthPO4_vr(L,NY,NX)        = 0._r8
      trcs_VLN_vr(ids_H1PO4B,L,NY,NX) = 0._r8
      trcs_VLN_vr(ids_H1PO4,L,NY,NX)  = 1._r8
      trcs_VLN_vr(ids_H2PO4B,L,NY,NX) = trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
      trcs_VLN_vr(ids_H2PO4,L,NY,NX)  = trcs_VLN_vr(ids_H1PO4,L,NY,NX)

      trcs_solml_vr(ids_H1PO4,L,NY,NX)  = trcs_solml_vr(ids_H1PO4,L,NY,NX)+trcs_solml_vr(ids_H1PO4B,L,NY,NX)
      trcs_solml_vr(ids_H2PO4,L,NY,NX)  = trcs_solml_vr(ids_H2PO4,L,NY,NX)+trcs_solml_vr(ids_H2PO4B,L,NY,NX)
      trcs_solml_vr(ids_H1PO4B,L,NY,NX) = 0._r8
      trcs_solml_vr(ids_H2PO4B,L,NY,NX) = 0._r8
      IF(salt_model)THEN
        trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)     = trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H0PO4B,L,NY,NX)
        trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX)     = trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H3PO4B,L,NY,NX)
        trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)    = trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX)
        trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)   = trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeH2PO4B,L,NY,NX)
        trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX)     = trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX)
        trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX)    = trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaHPO4B,L,NY,NX)
        trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)  = trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)+trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX)
        trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX)    = trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_MgHPO4B,L,NY,NX)
        trcSalt_solml_vr(idsalt_H0PO4B,L,NY,NX)    = 0._r8
        trcSalt_solml_vr(idsalt_H3PO4B,L,NY,NX)    = 0._r8
        trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX)   = 0._r8
        trcSalt_solml_vr(idsalt_FeH2PO4B,L,NY,NX)  = 0._r8
        trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX)    = 0._r8
        trcSalt_solml_vr(idsalt_CaHPO4B,L,NY,NX)   = 0._r8
        trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX) = 0._r8
        trcSalt_solml_vr(idsalt_MgHPO4B,L,NY,NX)   = 0._r8
      endif
      trcx_solml_vr(idx_OHe,L,NY,NX)    = trcx_solml_vr(idx_OHe,L,NY,NX)+trcx_solml_vr(idx_OHeB,L,NY,NX)
      trcx_solml_vr(idx_OH,L,NY,NX)     = trcx_solml_vr(idx_OH,L,NY,NX)+trcx_solml_vr(idx_OHB,L,NY,NX)
      trcx_solml_vr(idx_OHp,L,NY,NX)    = trcx_solml_vr(idx_OHp,L,NY,NX)+trcx_solml_vr(idx_OHpB,L,NY,NX)
      trcx_solml_vr(idx_HPO4,L,NY,NX)   = trcx_solml_vr(idx_HPO4,L,NY,NX)+trcx_solml_vr(idx_HPO4B,L,NY,NX)
      trcx_solml_vr(idx_H2PO4,L,NY,NX)  = trcx_solml_vr(idx_H2PO4,L,NY,NX)+trcx_solml_vr(idx_H2PO4B,L,NY,NX)
      trcx_solml_vr(idx_OHeB,L,NY,NX)   = 0._r8
      trcx_solml_vr(idx_OHB,L,NY,NX)    = 0._r8
      trcx_solml_vr(idx_OHpB,L,NY,NX)   = 0._r8
      trcx_solml_vr(idx_HPO4B,L,NY,NX)  = 0._r8
      trcx_solml_vr(idx_H2PO4B,L,NY,NX) = 0._r8

      trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)     = trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)+trcp_saltpml_vr(idsp_AlPO4B,L,NY,NX)
      trcp_saltpml_vr(idsp_FePO4,L,NY,NX)     = trcp_saltpml_vr(idsp_FePO4,L,NY,NX)+trcp_saltpml_vr(idsp_FePO4B,L,NY,NX)
      trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)    = trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)+trcp_saltpml_vr(idsp_CaHPO4B,L,NY,NX)
      trcp_saltpml_vr(idsp_Apatite,L,NY,NX)        = trcp_saltpml_vr(idsp_Apatite,L,NY,NX)+trcp_saltpml_vr(idsp_ApatiteBand,L,NY,NX)
      trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)  = trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)+trcp_saltpml_vr(idsp_CaH4P2O8B,L,NY,NX)
      trcp_saltpml_vr(idsp_AlPO4B,L,NY,NX)    = 0._r8
      trcp_saltpml_vr(idsp_FePO4B,L,NY,NX)    = 0._r8
      trcp_saltpml_vr(idsp_CaHPO4B,L,NY,NX)   = 0._r8
      trcp_saltpml_vr(idsp_ApatiteBand,L,NY,NX)       = 0._r8
      trcp_saltpml_vr(idsp_CaH4P2O8B,L,NY,NX) = 0._r8
    ENDIF
  ENDIF
  end subroutine UpdatePO4FertilizerBandinfo
!------------------------------------------------------------------------------------------

  subroutine UpdateNO3FertilizerBandinfo(L,NY,NX,FLWD)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in) :: FLWD
  real(r8) :: XVLNO3,DWNO3,DNO3S
  real(r8) :: DNO2S,FVLNO3,DPFLW

!     begin_execution
!     iFertNO3Band_col=banded NO3 fertilizer flag
!     ROWO=NO3 fertilizer band row width
!     DPNO3=NO3 fertilizer band depth
!
  IF(iFertNO3Band_col(NY,NX).EQ.ifert_on.AND.ROWSpaceNO3Band_col(NY,NX).GT.0.0)THEN

    IF(L.EQ.NU_col(NY,NX).OR.CumDepz2LayBottom_vr(L-1,NY,NX).LT.BandDepthNO3_col(NY,NX))THEN
!
!     NO3 BAND WIDTH
!
!     DWNO3=change in NO3 fertilizer band width
!     WDNOB=layer NO3 fertilizer band width
!     ZOSGL=NO3 diffusivity
!     TortMicPM_vr=tortuosity
!
      DWNO3=0.5_r8*SQRT(SoluteDifusvtyT_vr(ids_NO3,L,NY,NX))*TortMicPM_vr(NPH,L,NY,NX)
      BandWidthNO3_vr(L,NY,NX)=AMIN1(ROWSpaceNO3Band_col(NY,NX),BandWidthNO3_vr(L,NY,NX)+DWNO3)
!
!     NO3 BAND DEPTH
!
!     DPFLW=change in NO3 fertilizer band depth
!     BandDepthNO3_col,DPNOB=total,layer NO3 fertilizer band depth
!
      IF(CumDepz2LayBottom_vr(L,NY,NX).GE.BandDepthNO3_col(NY,NX))THEN
        DPFLW=FLWD+DWNO3
        call UpdateFertilizerBandDepth(L,NY,NX,DPFLW,BandDepthNO3_col(NY,NX), &
          BandThicknessNO3_vr(:,NY,NX),BandWidthNO3_vr(:,NY,NX))
      ENDIF
!
!     FRACTION OF SOIL LAYER OCCUPIED BY NO3 BAND
!     FROM BAND WIDTH X DEPTH
!
!     VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
!     DLYR=soil layer thickness
!     FVLNO3=relative change in VLNO3
!
      XVLNO3=trcs_VLN_vr(ids_NO3,L,NY,NX)
      IF(DLYR_3D(3,L,NY,NX).GT.ZERO)THEN
        trcs_VLN_vr(ids_NO3B,L,NY,NX)=AZMAX1(AMIN1(0.999,BandWidthNO3_vr(L,NY,NX) &
          /ROWSpaceNO3Band_col(NY,NX)*BandThicknessNO3_vr(L,NY,NX)/DLYR_3D(3,L,NY,NX)))
      ELSE
        trcs_VLN_vr(ids_NO3B,L,NY,NX)=0._r8
      ENDIF
      trcs_VLN_vr(ids_NO3,L,NY,NX)  = 1._r8-trcs_VLN_vr(ids_NO3B,L,NY,NX)
      trcs_VLN_vr(ids_NO2B,L,NY,NX) = trcs_VLN_vr(ids_NO3B,L,NY,NX)
      trcs_VLN_vr(ids_NO2,L,NY,NX)  = trcs_VLN_vr(ids_NO3,L,NY,NX)
      FVLNO3                        = AZMIN1((trcs_VLN_vr(ids_NO3,L,NY,NX)-XVLNO3)/XVLNO3)
!
!     TRANSFER NO3 FROM NON-BAND TO BAND
!     DURING BAND GROWTH
!
!     DNO3S,DNO2S=transfer of NO3,NO2
!
      DNO3S                                         = FVLNO3*trcs_solml_vr(ids_NO3,L,NY,NX)/natomw
      DNO2S                                         = FVLNO3*trcs_solml_vr(ids_NO2,L,NY,NX)/natomw
      trcn_RprodChem_soil_vr(ids_NO3,L,NY,NX)       = trcn_RprodChem_soil_vr(ids_NO3,L,NY,NX)+DNO3S
      trcn_RprodChem_soil_vr(ids_NO2,L,NY,NX)       = trcn_RprodChem_soil_vr(ids_NO2,L,NY,NX)+DNO2S
      trcn_RProdChem_band_soil_vr(ids_NO3B,L,NY,NX) = trcn_RProdChem_band_soil_vr(ids_NO3B,L,NY,NX)-DNO3S
      trcn_RProdChem_band_soil_vr(ids_NO2B,L,NY,NX) = trcn_RProdChem_band_soil_vr(ids_NO2B,L,NY,NX)-DNO2S
    ELSE
!
!     AMALGAMATE NO3 BAND WITH NON-BAND IF BAND NO LONGER EXISTS
!
      BandThicknessNO3_vr(L,NY,NX)  = 0._r8
      BandWidthNO3_vr(L,NY,NX)      = 0._r8
      trcs_VLN_vr(ids_NO3,L,NY,NX)  = 1._r8
      trcs_VLN_vr(ids_NO3B,L,NY,NX) = 0._r8
      trcs_VLN_vr(ids_NO2,L,NY,NX)  = trcs_VLN_vr(ids_NO3,L,NY,NX)
      trcs_VLN_vr(ids_NO2B,L,NY,NX) = trcs_VLN_vr(ids_NO3B,L,NY,NX)

      trcs_solml_vr(ids_NO3,L,NY,NX)  = trcs_solml_vr(ids_NO3,L,NY,NX)+trcs_solml_vr(ids_NO3B,L,NY,NX)
      trcs_solml_vr(ids_NO2,L,NY,NX)  = trcs_solml_vr(ids_NO2,L,NY,NX)+trcs_solml_vr(ids_NO2B,L,NY,NX)
      trcs_solml_vr(ids_NO3B,L,NY,NX) = 0._r8
      trcs_solml_vr(ids_NO2B,L,NY,NX) = 0._r8
    ENDIF
  ENDIF
  end subroutine UpdateNO3FertilizerBandinfo
!------------------------------------------------------------------------------------------

  subroutine UpdateSurfResidueSolute(NX,NY)
  implicit none
  integer, intent(in) :: NY,NX
  character(len=*), parameter :: subname='UpdateSurfResidueSolute'

  real(r8) :: Al_3p_aque_mole_conc,CALX,Ca_2p_aque_mole_conc,CaX_conc
  real(r8) :: CCEC0,CO3_2e_aque_mole_conc
  real(r8) :: Fe_3p_aque_mole_conc,CFEX
  real(r8) :: H2PO4_1e_AlPO4_eqv,H2PO4_1e_CaHPO4_eqv,H2PO4_1e_FePO4_eqv,H2PO4_1e_apatite_eqv,H2PO4_1e_CaH4P2O8_eqv
  real(r8) :: H2PO4_e_to_HPO4_2e_flx
  real(r8) :: H_1p_aque_mole_conc,OH_1e_aque_mole_conc,DP
  real(r8) :: FX,S0,S1,XALQ,XCAQ,XCAX
  real(r8) :: XFEQ,XHYQ,XN4Q_mole_conc,XTLQ
  real(r8) :: VLWatMicPMX,VLWatMicPMP,BulkSoilMass
  real(r8) :: ThetaWLitR,COMA
  real(r8) :: Urea_mole_conc
  real(r8) :: CCO20
  real(r8) :: REcoReleaz_H1PO4X,REcoReleaz_H2PO4X,DUKD,DFNSA,REcoReleaz_NH3X,REcoReleaz_NH4X
  real(r8) :: H1PO4_2e_aque_mole_conc,H2PO4_1e_aque_mole_conc,NH3_aque_mole_conc,NH4_1p_aque_mole_conc
  real(r8) :: Precp_AlPO4_mole_conc,Precp_CaHPO4_mole_conc
  real(r8) :: Precp_Ca5P3O12O3H3_mole_conc,Precp_CaH4P2O8_mole_conc,Precp_FePO4_mole_conc,RNH4,H2PO4_1e_AlPO4_dissol_flx
  real(r8) :: H2PO4_1e_CaHPO4_dissol_flx,H2PO4_1e_apatite_dissol_flx,H2PO4_1e_CaH4P2O8_dissol_flx
  real(r8) :: H2PO4_1e_FePO4_dissol_flx,RXN4,XNH4_mole_conc
  real(r8) :: RHP1,RHP2,RN3S,RN4S,RFertReleaz_NH4,RFertReleaz_Urea
  real(r8) :: RFertReleaz_NH3,RFertReleaz_Urea2SoilNH3,RFertReleaz_NO3

  !     begin_execution
  !     BKVL=litter mass
  !
  call PrintInfo('beg '//subname)
  IF(VLWatMicPM_vr(NPH,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    BulkSoilMass=VLSoilMicPMass_vr(0,NY,NX)
!
!     UREA HYDROLYSIS IN SURFACE RESIDUE
!
!     VWatMicrobAct_vr=biologically active litter water volume from nitro.f
!     COMA=concentration of active biomass
!     TMicHeterActivity_vr=total microbial activity from nitro.f
!     DUKD=Km for urea hydrolysis
!
    IF(VWatMicrobAct_vr(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      COMA=AMIN1(0.1E+06_r8,TMicHeterActivity_vr(0,NY,NX)/VWatMicrobAct_vr(0,NY,NX))
    ELSE
      COMA=0.1E+06_r8
    ENDIF
    DUKD=DUKM*(1.0_r8+COMA/DUKI)
!
!     UREA HYDROLYSIS INHIBITION
!
!     ZNHU0,ZNHUI=initial,current inhibition activity
!     RUreaInhibtorConst=rate constant for decline in urea hydrolysis inhibition
!
    IF(ZNHU0_vr(0,NY,NX).GT.ZEROS(NY,NX) .AND. ZNHUI_vr(0,NY,NX).GT.ZEROS(NY,NX))THEN
      ZNHUI_vr(0,NY,NX)=ZNHUI_vr(0,NY,NX)-TSens4MicbGrwoth_vr(0,NY,NX)**0.25_r8 &
        *RUreaInhibtorConst(iUreaHydInhibitorType_col(NY,NX))*ZNHUI_vr(0,NY,NX) &
        *AMAX1(RUreaInhibtorConst(iUreaHydInhibitorType_col(NY,NX)),1.0_r8-ZNHUI_vr(0,NY,NX)/ZNHU0_vr(0,NY,NX))
    ELSE
      ZNHUI_vr(0,NY,NX)=0._r8
    ENDIF
    !
    !     UREA CONCENTRATION AND HYDROLYSIS IN SURFACE RESIDUE
    !
    !     ZNHUFA=urea fertilizer
    !     Urea_mole_conc=concentration of urea fertilizer
    !     DFNSA=effect of microbial concentration on urea hydrolysis
    !     RFertReleaz_Urea=rate of urea hydrolysis
    !     RFertUreaSpecHydrol=specific rate constant for urea hydrolysis
    !     TSens4MicbGrwoth_vr=temperature effect on microbial activity from nitro.f
    !
    IF(FertN_mole_soil_vr(ifert_N_urea,0,NY,NX).GT.ZEROS(NY,NX) &
      .AND. VLSoilMicPMass_vr(0,NY,NX).GT.ZEROS(NY,NX))THEN
      Urea_mole_conc = FertN_mole_soil_vr(ifert_N_urea,0,NY,NX)/VLSoilMicPMass_vr(0,NY,NX)
      DFNSA = Urea_mole_conc/(Urea_mole_conc+DUKD)
      RFertReleaz_Urea = AMIN1(FertN_mole_soil_vr(ifert_N_urea,0,NY,NX),RFertUreaSpecHydrol*TMicHeterActivity_vr(0,NY,NX) &
        *DFNSA*TSens4MicbGrwoth_vr(0,NY,NX))*(1.0_r8-ZNHUI_vr(0,NY,NX))
    ELSE
      RFertReleaz_Urea=0._r8
    ENDIF
!
!     NH4, NH3, UREA, NO3 DISSOLUTION IN SURFACE RESIDUE
!     FROM FIRST-ORDER FUNCTIONS OF REMAINING
!     FERTILIZER (NOTE: SUPERPHOSPHATE AND ROCK PHOSPHATE
!     ARE REPRESENTED AS MONOCALCIUM PHOSPHATE AND HYDROXYAPATITE
!     MODELLED IN PHOSPHORUS REACTIONS BELOW)
!
!     RFertReleaz_NH4=rate of broadcast NH4 fertilizer dissoln
!     RFertReleaz_NH3=rate of broadcast NH3 fertilizer dissoln
!     RFertReleaz_Urea2SoilNH3=rate of broadcast urea fertr dissoln
!     RFertReleaz_NO3=rate of broadcast NO3 fertilizer dissoln
!
    IF(VWatLitRHoldCapcity_col(NY,NX).GT.ZEROS(NY,NX))THEN
      ThetaWLitR=AMIN1(1.0,VLWatMicP_vr(0,NY,NX)/VWatLitRHoldCapcity_col(NY,NX))
    ELSE
      ThetaWLitR=1._r8
    ENDIF
    RFertReleaz_NH4          = RFertNH4SpecReleaz*FertN_mole_soil_vr(ifert_N_nh4,0,NY,NX)*ThetaWLitR
    RFertReleaz_NH3          = RFertNH3SpecReleaz*FertN_mole_soil_vr(ifert_N_nh3,0,NY,NX)
    RFertReleaz_Urea2SoilNH3 = RFertReleaz_Urea*ThetaWLitR
    RFertReleaz_NO3          = RFertNO3SpecReleaz*FertN_mole_soil_vr(ifert_N_no3,0,NY,NX)*ThetaWLitR
!
!     SOLUBLE AND EXCHANGEABLE NH4 CONCENTRATIONS
!
!     VLWatMicPNH=water volume
!     REcoReleaz_NH4X,REcoReleaz_NH3X=NH4,NH3 input from uptake, mineraln, dissoln
!     RNH4MicbReliz2Soil_vr=net change in NH4 from nitro.f
!     NH4_1p_aque_mole_conc,NH3_aque_mole_conc=total NH4,NH3 concentration
!     XNH4_mole_conc=adsorbed NH4 concentration
!
    IF(VLWatMicPM_vr(NPH,0,NY,NX).GT.ZEROS2(NY,NX))THEN
      VLWatMicPMX           = natomw*VLWatMicPM_vr(NPH,0,NY,NX)
      REcoReleaz_NH4X       = (RNut_MicbRelease_vr(ids_NH4,0,NY,NX)+natomw*RFertReleaz_NH4)/VLWatMicPMX
      REcoReleaz_NH3X       = natomw*RFertReleaz_Urea2SoilNH3/VLWatMicPMX
      NH4_1p_aque_mole_conc = AMAX1(ZERO,trcs_solml_vr(ids_NH4,0,NY,NX)/VLWatMicPMX+REcoReleaz_NH4X)
      NH3_aque_mole_conc    = AMAX1(ZERO,trcs_solml_vr(idg_NH3,0,NY,NX)/VLWatMicPMX+REcoReleaz_NH3X)
      IF(BulkSoilMass.GT.ZEROS(NY,NX))THEN
        XNH4_mole_conc=AMAX1(ZERO,trcx_solml_vr(idx_NH4,0,NY,NX)/BulkSoilMass)
      ELSE
        XNH4_mole_conc=0._r8
      ENDIF
!
!     SOLUBLE, EXCHANGEABLE AND PRECIPITATED PO4 CONCENTRATIONS
!
!     VLWatMicPMP=water volume
!     REcoReleaz_H1PO4X,REcoReleaz_H2PO4X=HPO4,H2PO4 inputs from mineraln, uptake
!     RH1PO4MicbReliz2Soil_vr=net change in HPO4 from nitro.f
!     H1PO4_2e_aque_mole_conc,H2PO4_1e_aque_mole_conc=HPO4,H2PO4 concentrations
!
      VLWatMicPMP             = patomw*VLWatMicPM_vr(NPH,0,NY,NX)
      REcoReleaz_H1PO4X       = RNut_MicbRelease_vr(ids_H1PO4,0,NY,NX)/VLWatMicPMP
      REcoReleaz_H2PO4X       = RNut_MicbRelease_vr(ids_H2PO4,0,NY,NX)/VLWatMicPMP
      H1PO4_2e_aque_mole_conc = AZMAX1(trcs_solml_vr(ids_H1PO4,0,NY,NX)/VLWatMicPMP+REcoReleaz_H1PO4X)
      H2PO4_1e_aque_mole_conc = AZMAX1(trcs_solml_vr(ids_H2PO4,0,NY,NX)/VLWatMicPMP+REcoReleaz_H2PO4X)
    ELSE
      REcoReleaz_NH4X         = 0._r8
      REcoReleaz_NH3X         = 0._r8
      NH4_1p_aque_mole_conc   = 0._r8
      NH3_aque_mole_conc      = 0._r8
      XNH4_mole_conc          = 0._r8
      REcoReleaz_H1PO4X       = 0._r8
      REcoReleaz_H2PO4X       = 0._r8
      H1PO4_2e_aque_mole_conc = 0._r8
      H2PO4_1e_aque_mole_conc = 0._r8
    ENDIF
!
!     PHOSPHORUS TRANSFORMATIONS IN SURFACE RESIDUE
!
!     Precp_AlPO4_mole_conc,Precp_FePO4_mole_conc=concn of precip AlPO4,FEPO4
!     Precp_CaH4P2O8_mole_conc,Precp_CaHPO4_mole_conc,Precp_Ca5P3O12O3H3_mole_conc=concn of precip CaH4P2O8,CaHPO4,apatite
!
    IF(BulkSoilMass.GT.ZEROS(NY,NX))THEN
      Precp_AlPO4_mole_conc        = AZMAX1(trcp_saltpml_vr(idsp_AlPO4,0,NY,NX)/BulkSoilMass)
      Precp_FePO4_mole_conc        = AZMAX1(trcp_saltpml_vr(idsp_FePO4,0,NY,NX)/BulkSoilMass)
      Precp_CaH4P2O8_mole_conc     = AZMAX1(trcp_saltpml_vr(idsp_CaH4P2O8,0,NY,NX)/BulkSoilMass)
      Precp_CaHPO4_mole_conc       = AZMAX1(trcp_saltpml_vr(idsp_CaHPO4,0,NY,NX)/BulkSoilMass)
      Precp_Ca5P3O12O3H3_mole_conc = AZMAX1(trcp_saltpml_vr(idsp_Apatite,0,NY,NX)/BulkSoilMass)
    ELSE
      Precp_AlPO4_mole_conc        = 0._r8
      Precp_FePO4_mole_conc        = 0._r8
      Precp_CaH4P2O8_mole_conc     = 0._r8
      Precp_CaHPO4_mole_conc       = 0._r8
      Precp_Ca5P3O12O3H3_mole_conc = 0._r8
    ENDIF
!
!     CALCULATE H2PO4 COPRECIPITATES FRPM LITTER PH
!
    H_1p_aque_mole_conc   = AMAX1(ZERO,10.0_r8**(-(PH_vr(0,NY,NX)-3.0_r8)))
    OH_1e_aque_mole_conc  = AMAX1(ZERO,DPH2O/H_1p_aque_mole_conc)
    Al_3p_aque_mole_conc  = AMAX1(ZERO,SPALO/OH_1e_aque_mole_conc**3._r8)
    Fe_3p_aque_mole_conc  = AMAX1(ZERO,SPFEO/OH_1e_aque_mole_conc**3_r8)
    CCO20                 = AMAX1(ZERO,trc_solcl_vr(idg_CO2,0,NY,NX)/catomw)
    CO3_2e_aque_mole_conc = AMAX1(ZERO,CCO20*DPCO3/H_1p_aque_mole_conc**2._r8)
    Ca_2p_aque_mole_conc  = AMAX1(ZERO,AMIN1(CCAMX,SPCAC/CO3_2e_aque_mole_conc))
!
!     ALUMINUM PHOSPHATE (VARISCITE)
!
!     H2PO4_1e_AlPO4_eqv,H2PO4_1e_aque_mole_conc=equilibrium,current H2PO4 concentration in litter
!     SYA0P2=solubility product derived from SPALO
!     H2PO4_1e_AlPO4_dissol_flx=H2PO4 dissolution from AlPO4 in litter
!
    H2PO4_1e_AlPO4_eqv        = SYA0P2/(Al_3p_aque_mole_conc*OH_1e_aque_mole_conc**2)
    H2PO4_1e_AlPO4_dissol_flx = AMIN1(AZMAX1(4.0E-08*SoilOrgM_vr(ielmc,0,NY,NX)-Precp_AlPO4_mole_conc),&
      AMAX1(-Precp_AlPO4_mole_conc,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_AlPO4_eqv)))
!
!     IRON PHOSPHATE (STRENGITE)
!
!     H2PO4_1e_FePO4_eqv,H2PO4_1e_aque_mole_conc=equilibrium,current H2PO4 concentration in litter
!     SYF0P2=solubility product derived from SPALO
!     H2PO4_1e_FePO4_dissol_flx=H2PO4 dissolution from FePO4 in litter
!
    H2PO4_1e_FePO4_eqv        = SYF0P2/(Fe_3p_aque_mole_conc*OH_1e_aque_mole_conc**2)
    H2PO4_1e_FePO4_dissol_flx = AMIN1(AZMAX1(2.0E-06*SoilOrgM_vr(ielmc,0,NY,NX)-Precp_FePO4_mole_conc),&
      AMAX1(-Precp_FePO4_mole_conc,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_FePO4_eqv)))
!
!     DICALCIUM PHOSPHATE
!
!     H2PO4_1e_CaHPO4_eqv,H2PO4_1e_aque_mole_conc=equilibrium,current H2PO4 concentration in litter
!     SYCAD2=solubility product derived from SPALO
!     H2PO4_1e_CaHPO4_dissol_flx=H2PO4 dissolution from CaHPO4 in litter
!   CaHPO4.H2O <->HPO4(-)+Ca(2+)+OH(-) 
    H2PO4_1e_CaHPO4_eqv        = SYCAD2/(Ca_2p_aque_mole_conc*OH_1e_aque_mole_conc)
    H2PO4_1e_CaHPO4_dissol_flx = AMIN1(AMAX1(0-.0,5.0E-05*SoilOrgM_vr(ielmc,0,NY,NX)-Precp_CaHPO4_mole_conc),&
      AMAX1(-Precp_CaHPO4_mole_conc,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_CaHPO4_eqv)))
!
!     HYDROXYAPATITE
!
!     H2PO4_1e_apatite_eqv,H2PO4_1e_aque_mole_conc=equilibrium,current H2PO4 concentration in litter
!     SYCAH2=solubility product derived from SPALO
!     H2PO4_1e_apatite_dissol_flx=H2PO4 dissolution from apatite in litter
!
    H2PO4_1e_apatite_eqv        = (SYCAH2/(Ca_2p_aque_mole_conc**5_r8*OH_1e_aque_mole_conc**7))**0.333_r8
    H2PO4_1e_apatite_dissol_flx = AMIN1(AZMAX1(5.0E-05_r8*SoilOrgM_vr(ielmc,0,NY,NX)-Precp_Ca5P3O12O3H3_mole_conc),&
      AMAX1(-Precp_Ca5P3O12O3H3_mole_conc,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_apatite_eqv)))
!
!     MONOCALCIUM PHOSPHATE
!
!     H2PO4_1e_CaH4P2O8_eqv,H2PO4_1e_aque_mole_conc=equilibrium,current H2PO4 concentration in litter
!     SPCAM=solubility product for Ca(H2PO4)2
!     H2PO4_1e_CaH4P2O8_dissol_flx=H2PO4 dissolution from Ca(H2PO4)2 in litter
!
    H2PO4_1e_CaH4P2O8_eqv        = SQRT(SPCAM/Ca_2p_aque_mole_conc)
    H2PO4_1e_CaH4P2O8_dissol_flx = AMIN1(AZMAX1(5.0E-05_r8*SoilOrgM_vr(ielmc,0,NY,NX)-Precp_CaH4P2O8_mole_conc),&
      AMAX1(-Precp_CaH4P2O8_mole_conc*SPPO4,TPD*(H2PO4_1e_aque_mole_conc-H2PO4_1e_CaH4P2O8_eqv)))

!
!     PHOSPHORUS ANION EXCHANGE IN SURFACE REDISUE
!     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
!     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
!     EXCHANGE SITES (NOT CALCULATED)
!
!     H2PO4-H+HPO4
!
!     DPH2P=dissociation constant
!     S1=equilibrium concentration in litter
!     H2PO4_e_to_HPO4_2e_flx=H2PO4-H+HPO4 dissociation in litter
!
    DP                     = DPH2P
    S0                     = H1PO4_2e_aque_mole_conc+H_1p_aque_mole_conc+DP
    ! Keep small positive discriminants: concentration cutoffs do not apply to squares.
    ! Use a sum of nonnegative terms and 4*(B*C-K*A)/(S0+sqrt(D)) to avoid cancellation.
    S1                     = MAX(0._r8,(H1PO4_2e_aque_mole_conc-H_1p_aque_mole_conc)**2+2._r8*DP*(H1PO4_2e_aque_mole_conc+H_1p_aque_mole_conc)+DP**2+4._r8*DP*H2PO4_1e_aque_mole_conc)
    H2PO4_e_to_HPO4_2e_flx = TSL*(4._r8*((H1PO4_2e_aque_mole_conc*H_1p_aque_mole_conc-DP*H2PO4_1e_aque_mole_conc)/(S0+SQRT(S1))))
!
!     EQUILIBRIUM X-CA CONCENTRATION FROM CEC, GAPON COEFFICIENTS
!     AND CATION CONCENTRATIONS
!
!     XCEC=cation exchange concentration,capacity
!     XCAX=equilibrium R-Ca concentration
!     GKC4,GKCH,GKCA,GKCM,GKCN,GKCK=Gapon selectivity coefficients for
!     CA-NH4,CA-H,CA-AL,CA-MG,CA-NA,CA-K
!     X*Q=equilibrium exchangeable concentrations
!     XTLQ=total equilibrium exchangeable concentration
!     BulkSoilMass=Mg/d2
!     ORGC_vr=gC/d-2
!     COOH=mol/gC, CEC of organic matter
    IF(BulkSoilMass.GT.ZEROS(NY,NX))THEN
      CCEC0=AMAX1(ZERO,COOH*SoilOrgM_vr(ielmc,0,NY,NX)/BulkSoilMass)
    ELSE
      CCEC0=ZERO
    ENDIF
    CALX     = AMAX1(ZERO,Al_3p_aque_mole_conc)**0.333_r8
    CFEX     = AMAX1(ZERO,Fe_3p_aque_mole_conc)**0.333_r8
    CaX_conc = AMAX1(ZERO,Ca_2p_aque_mole_conc)**0.500_r8
    !
    !     EQUILIBRIUM X-CA CONCENTRATION FROM CEC AND CATION
    !     CONCENTRATIONS
    !
    XCAX=CCEC0/(1.0_r8+GKC4_vr(NU_col(NY,NX),NY,NX)*NH4_1p_aque_mole_conc/CaX_conc &
      +GKCH_vr(NU_col(NY,NX),NY,NX)*H_1p_aque_mole_conc/CaX_conc &
      +GKCA_vr(NU_col(NY,NX),NY,NX)*CALX/CaX_conc &
      +GKCA_vr(NU_col(NY,NX),NY,NX)*CFEX/CaX_conc)
    XN4Q_mole_conc = XCAX*NH4_1p_aque_mole_conc*GKC4_vr(NU_col(NY,NX),NY,NX)
    XHYQ = XCAX*H_1p_aque_mole_conc*GKCH_vr(NU_col(NY,NX),NY,NX)
    XALQ = XCAX*CALX*GKCA_vr(NU_col(NY,NX),NY,NX)
    XFEQ = XCAX*CFEX*GKCA_vr(NU_col(NY,NX),NY,NX)
    XCAQ = XCAX*CaX_conc
    XTLQ = XN4Q_mole_conc+XHYQ+XALQ+XFEQ+XCAQ
    IF(XTLQ.GT.ZERO)THEN
      FX=CCEC0/XTLQ
    ELSE
      FX=0._r8
    ENDIF
    XN4Q_mole_conc=FX*XN4Q_mole_conc
!
!     NH4 AND NH3 EXCHANGE IN SURFACE RESIDUE
!
!     RXN4=NH4 adsorption in litter
!     TADC0=adsorption rate constant
!     RNH4=NH4-NH3+H dissociation in litter
!     DPN4=NH4 dissociation constant
!
    RXN4 = AMAX1(TADC0*(XN4Q_mole_conc-XNH4_mole_conc),-trcx_solml_vr(idx_NH4,0,NY,NX)/VLWatMicPM_vr(NPH,0,NY,NX)*0.9999_r8)
    RNH4 = (H_1p_aque_mole_conc*NH3_aque_mole_conc-DPN4*NH4_1p_aque_mole_conc)/(DPN4+H_1p_aque_mole_conc)

  ELSE
    RFertReleaz_NH4              = 0._r8
    RFertReleaz_NH3              = 0._r8
    RFertReleaz_Urea2SoilNH3     = 0._r8
    RFertReleaz_NO3              = 0._r8
    H2PO4_1e_AlPO4_dissol_flx    = 0._r8
    H2PO4_1e_FePO4_dissol_flx    = 0._r8
    H2PO4_1e_CaHPO4_dissol_flx   = 0._r8
    H2PO4_1e_apatite_dissol_flx  = 0._r8
    H2PO4_1e_CaH4P2O8_dissol_flx = 0._r8
    RXN4                         = 0._r8
    RNH4                         = 0._r8
    H2PO4_e_to_HPO4_2e_flx       = 0._r8
  ENDIF
!     TOTAL ION FLUXES FOR ALL REACTIONS ABOVE
!
!     RN4S=net NH4 flux in litter
!     RN3S=net NH3 flux in litter
!     RHP1,RHP2=net HPO4,H2PO4 flux in litter
!
      RN4S = RNH4-RXN4
      RN3S = -RNH4
      RHP1 = -H2PO4_e_to_HPO4_2e_flx
      RHP2 = H2PO4_e_to_HPO4_2e_flx-H2PO4_1e_AlPO4_dissol_flx-H2PO4_1e_FePO4_dissol_flx &
        -H2PO4_1e_CaHPO4_dissol_flx-2.0*H2PO4_1e_CaH4P2O8_dissol_flx-3.0_r8*H2PO4_1e_apatite_dissol_flx
!
!     CONVERT TOTAL ION FLUXES FROM CHANGES IN CONCENTRATION
!     TO CHANGES IN MASS PER UNIT AREA FOR USE IN 'REDIST'
!
!     TRChem_NH4_soil=total NH4 flux
!     TRChem_NH3_soil_mole=total NH3 flux
!     TRChem_H1PO4_soil,TRChem_H2PO4_soil=net HPO4,H2PO4 flux
!     TRNX4=total NH4 adsorption
!     TRChem_AlPO4_precip_soil,TRChem_FePO4_precip_soil,TRChem_CaHPO4_precip_soil,TRChem_apatite_precip_soil,TRChem_CaH4P2O8_precip_soil
!     =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation
!
      trcn_RprodChem_soil_vr(ids_NH4,0,NY,NX)     = trcn_RprodChem_soil_vr(ids_NH4,0,NY,NX)+RN4S*VLWatMicPM_vr(NPH,0,NY,NX)
      TRProd_chem_sol_NH3_soil_vr(0,NY,NX)           = TRProd_chem_sol_NH3_soil_vr(0,NY,NX)+RN3S*VLWatMicPM_vr(NPH,0,NY,NX)
      trcn_RprodChem_soil_vr(ids_H1PO4,0,NY,NX)   = trcn_RprodChem_soil_vr(ids_H1PO4,0,NY,NX)+RHP1*VLWatMicPM_vr(NPH,0,NY,NX)
      trcn_RprodChem_soil_vr(ids_H2PO4,0,NY,NX)   = trcn_RprodChem_soil_vr(ids_H2PO4,0,NY,NX)+RHP2*VLWatMicPM_vr(NPH,0,NY,NX)
      trcx_TRSoilChem_vr(idx_NH4,0,NY,NX)       = trcx_TRSoilChem_vr(idx_NH4,0,NY,NX)+RXN4*VLWatMicPM_vr(NPH,0,NY,NX)
      
      trcp_RChem_soil_vr(idsp_AlPO4,0,NY,NX)    = trcp_RChem_soil_vr(idsp_AlPO4,0,NY,NX)+H2PO4_1e_AlPO4_dissol_flx*VLWatMicPM_vr(NPH,0,NY,NX)
      trcp_RChem_soil_vr(idsp_FePO4,0,NY,NX)    = trcp_RChem_soil_vr(idsp_FePO4,0,NY,NX)+H2PO4_1e_FePO4_dissol_flx*VLWatMicPM_vr(NPH,0,NY,NX)
      trcp_RChem_soil_vr(idsp_CaHPO4,0,NY,NX)   = trcp_RChem_soil_vr(idsp_CaHPO4,0,NY,NX)+H2PO4_1e_CaHPO4_dissol_flx*VLWatMicPM_vr(NPH,0,NY,NX)
      trcp_RChem_soil_vr(idsp_Apatite,0,NY,NX)       = trcp_RChem_soil_vr(idsp_Apatite,0,NY,NX)+H2PO4_1e_apatite_dissol_flx*VLWatMicPM_vr(NPH,0,NY,NX)
      trcp_RChem_soil_vr(idsp_CaH4P2O8,0,NY,NX) = trcp_RChem_soil_vr(idsp_CaH4P2O8,0,NY,NX)+H2PO4_1e_CaH4P2O8_dissol_flx*VLWatMicPM_vr(NPH,0,NY,NX)

      FertN_mole_soil_vr(ifert_N_nh4,0,NY,NX)     = FertN_mole_soil_vr(ifert_N_nh4,0,NY,NX)-RFertReleaz_NH4
      FertN_mole_soil_vr(ifert_N_nh3,0,NY,NX)     = FertN_mole_soil_vr(ifert_N_nh3,0,NY,NX)-RFertReleaz_NH3
      FertN_mole_soil_vr(ifert_N_urea,0,NY,NX)    = FertN_mole_soil_vr(ifert_N_urea,0,NY,NX)-RFertReleaz_Urea2SoilNH3
      FertN_mole_soil_vr(ifert_N_no3,0,NY,NX)     = FertN_mole_soil_vr(ifert_N_no3,0,NY,NX)-RFertReleaz_NO3
      trcn_RprodChem_soil_vr(ids_NH4,0,NY,NX)     = trcn_RprodChem_soil_vr(ids_NH4,0,NY,NX)+RFertReleaz_NH4
      TRProd_chem_sol_NH3_soil_vr(0,NY,NX)           = TRProd_chem_sol_NH3_soil_vr(0,NY,NX)+RFertReleaz_NH3+RFertReleaz_Urea2SoilNH3
      trcn_RprodChem_soil_vr(ids_NO3,0,NY,NX)     = trcn_RprodChem_soil_vr(ids_NO3,0,NY,NX)+RFertReleaz_NO3

      trcn_RprodChem_soil_vr(ids_NH4,0,NY,NX)     = trcn_RprodChem_soil_vr(ids_NH4,0,NY,NX)*natomw
      TRProd_chem_sol_NH3_soil_vr(0,NY,NX)           = TRProd_chem_sol_NH3_soil_vr(0,NY,NX)*natomw
      trcn_RprodChem_soil_vr(ids_NO3,0,NY,NX)     = trcn_RprodChem_soil_vr(ids_NO3,0,NY,NX)*natomw
      trcn_RprodChem_soil_vr(ids_H1PO4,0,NY,NX)   = trcn_RprodChem_soil_vr(ids_H1PO4,0,NY,NX)*patomw
      trcn_RprodChem_soil_vr(ids_H2PO4,0,NY,NX)   = trcn_RprodChem_soil_vr(ids_H2PO4,0,NY,NX)*patomw

!9989  FORMAT(A8,2I4,12E12.4)
  call PrintInfo('end '//subname)
  end subroutine UpdateSurfResidueSolute

END module SoluteMod
