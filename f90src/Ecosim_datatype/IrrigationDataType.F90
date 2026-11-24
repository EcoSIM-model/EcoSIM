module IrrigationDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use EcoSIMConfig, only : jcplx => jcplxc
  use TracerIDMod
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8) :: DIRRI(12)                                              !change factor for irrigation, [-]
  real(r8) :: CCOU                                                   !subsurface irrigation  CO2 concentration,	[g m-3]
  real(r8) :: CCHU                                                   !subsurface irrigation  CH4 concentration,	[g m-3]
  real(r8) :: COXU                                                   !subsurface irrigation  O2 concentration,	[g m-3]
  real(r8) :: CNNU                                                   !subsurface irrigation  N2 concentration,	[g m-3]
  real(r8) :: CN2U                                                   !subsurface irrigation  N2O concentration	[g m-3]
  ! Annual management forcing, reread from inputs on startup and restart.
  logical, allocatable :: Lirri_auto_col(:,:)                         !automatic irrigation enabled by management, [-]
  integer ,target,allocatable ::  IIRRA(:,:,:)                       !start and end dates of automated irrigation, [-]
  real(r8),target,allocatable ::  RRIG(:,:,:,:)                      !irrigation application depth, [m h-1]
  real(r8),target,allocatable ::  WDPTH(:,:,:)                       !depth of irrigation application, [m]
  real(r8),target,allocatable ::  IrrigSubsurf_col(:,:)              !underground irrigation, [m3 d-2 h-1]
  real(r8),target,allocatable ::  IrrigSurface_col(:,:)              !surface irrigation, [m3 d-2 h-1]
  real(r8),target,allocatable ::  FIRRA_col(:,:)                     !fraction of FC-WP below which automatic irrigation applied, [-]
  real(r8),target,allocatable ::  CIRRA_col(:,:)                     !fraction of FC-WP to which automatic irrigation applied, [-]
  real(r8),target,allocatable ::  DIRRA(:,:,:)                       !depth to which automatic irrigation applied, [m]
  real(r8),target,allocatable ::  TDIRI(:,:,:)                       !accumulated relative change  for irrigation, [-]
  real(r8),target,allocatable ::  PHQ(:,:,:)                         !surface irrigation  pH, [-]
  real(r8),target,allocatable ::  NH4_irrig_mole_conc(:,:,:)         !surface irrigation  NH4 concentration, [g m-3]
  real(r8),target,allocatable ::  NO3_irrig_mole_conc(:,:,:)         !surface irrigation  NO3 concentration, [g m-3]
  real(r8),target,allocatable ::  H2PO4_irrig_mole_conc(:,:,:)       !surface irrigation  H2PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CALQ(:,:,:)                       !surface irrigation  Al concentration, [g m-3]
  real(r8),target,allocatable ::  CFEQ(:,:,:)                       !surface irrigation  Fe concentration, [g m-3]
  real(r8),target,allocatable ::  CHYQ(:,:,:)                       !surface irrigation  H concentration, [g m-3]
  real(r8),target,allocatable ::  CCAQ(:,:,:)                       !surface irrigation  Ca concentration, [g m-3]
  real(r8),target,allocatable ::  CMGQ(:,:,:)                       !surface irrigation  Mg concentration, [g m-3]
  real(r8),target,allocatable ::  CNAQ(:,:,:)                       !surface irrigation  Na concentration, [g m-3]
  real(r8),target,allocatable ::  CKAQ(:,:,:)                       !surface irrigation  K concentration, [g m-3]
  real(r8),target,allocatable ::  COHQ(:,:,:)                       !surface irrigation  OH concentration, [g m-3]
  real(r8),target,allocatable ::  CSOQ(:,:,:)                       !surface irrigation  SO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CCLQ(:,:,:)                       !surface irrigation  Cl concentration, [g m-3]
  real(r8),target,allocatable ::  CC3Q(:,:,:)                       !surface irrigation  CO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CHCQ(:,:,:)                       !surface irrigation  HCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CAL1Q(:,:,:)                      !surface irrigation  AlOH concentration, [g m-3]
  real(r8),target,allocatable ::  CAL2Q(:,:,:)                      !surface irrigation  AlOH2 concentration, [g m-3]
  real(r8),target,allocatable ::  CAL3Q(:,:,:)                      !surface irrigation  AlOH3 concentration, [g m-3]
  real(r8),target,allocatable ::  CAL4Q(:,:,:)                      !surface irrigation  AlOH4 concentration, [g m-3]
  real(r8),target,allocatable ::  CALSQ(:,:,:)                      !surface irrigation  AlSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CFE1Q(:,:,:)                      !surface irrigation  FeOH concentration, [g m-3]
  real(r8),target,allocatable ::  CFE2Q(:,:,:)                      !surface irrigation  FeOH2 concentration, [g m-3]
  real(r8),target,allocatable ::  CFE3Q(:,:,:)                      !surface irrigation  FeOH3 concentration, [g m-3]
  real(r8),target,allocatable ::  CFE4Q(:,:,:)                      !surface irrigation  FeOH4 concentration, [g m-3]
  real(r8),target,allocatable ::  CFESQ(:,:,:)                      !surface irrigation  FeSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CCAOQ(:,:,:)                      !surface irrigation  CaOH concentration, [g m-3]
  real(r8),target,allocatable ::  CCACQ(:,:,:)                      !surface irrigation  CaCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CCAHQ(:,:,:)                      !surface irrigation  CaHCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CCASQ(:,:,:)                      !surface irrigation  CaSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CMGOQ(:,:,:)                      !surface irrigation  MgOH concentration, [g m-3]
  real(r8),target,allocatable ::  CMGCQ(:,:,:)                      !surface irrigation  MgCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CMGHQ(:,:,:)                      !surface irrigation  MgHCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CMGSQ(:,:,:)                      !surface irrigation  MgSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CNACQ(:,:,:)                      !surface irrigation  NaCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CNASQ(:,:,:)                      !surface irrigation  NaSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CKASQ(:,:,:)                      !surface irrigation  K concentration, [g m-3]
  real(r8),target,allocatable ::  CH0PQ(:,:,:)                      !surface irrigation  PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  HPO4_irrig_mole_conc(:,:,:)       !surface irrigation  HPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CH3PQ(:,:,:)                      !surface irrigation  H3PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CF1PQ(:,:,:)                      !surface irrigation  FeHPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CF2PQ(:,:,:)                      !surface irrigation  FeH2PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CC0PQ(:,:,:)                      !surface irrigation  CaPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CC1PQ(:,:,:)                      !surface irrigation  CaHPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CC2PQ(:,:,:)                      !surface irrigation  CaH4P2O8 concentration, [g m-3]
  real(r8),target,allocatable ::  CM1PQ(:,:,:)                      !surface irrigation  MgHPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CSTRQ(:,:,:)                      !surface irrigation ion strength, [g m-3]
  real(r8),target,allocatable ::  SurfIrrig_IonStrenth_col(:,:)     !surface irrigation ion strength, [g m-3]
  real(r8),target,allocatable ::  trcg_irrig_mole_conc_col(:,:,:)   !surface irrigation  volatile concentration, [mol m-3]
  real(r8),target,allocatable ::  trcn_irrig_vr(:,:,:,:)            !subsurface irrigation  nutrient concentration, [g m-3]
  real(r8),target,allocatable ::  trcSalt_Irrig_vr(:,:,:,:)         !prescribed external-water salt concentration, unbanded [mol m-3]
  real(r8),target,allocatable ::  FWatIrrigate2MicP_vr(:,:,:)       !underground irrigation, [m3 d-2 h-1]
  real(r8),target,allocatable ::  HeatIrrigation_vr(:,:,:)          !convective heat due to underground irrigation, [MJ d-2 h-1]
  real(r8),target,allocatable ::  trcs_Irrig_flx_vr(:,:,:,:)            !aqueous non-salt solutes in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  trcsalt_irrig_mole_conc_col(:,:,:,:)        !salt tracer concentration in irrigation, [g m-3]
  real(r8),target,allocatable ::  trcs_irrig_flx_col(:,:,:)         !tracer flux through irrigation, [g d-2 h-1]
  private :: InitAllocate
  contains

  subroutine InitIrrigation

  implicit none

  call InitAllocate
  end subroutine InitIrrigation

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none

  allocate(PHQ(366,JY,JX));     PHQ=0._r8
  allocate(trcs_irrig_flx_col(ids_beg:ids_end,JY,JX)); trcs_irrig_flx_col=0._r8
  allocate(NH4_irrig_mole_conc(366,JY,JX));    NH4_irrig_mole_conc=0._r8
  allocate(trcg_irrig_mole_conc_col(idg_beg:idg_NH3,JY,JX));   trcg_irrig_mole_conc_col=0._r8
  allocate(NO3_irrig_mole_conc(366,JY,JX));    NO3_irrig_mole_conc=0._r8
  allocate(H2PO4_irrig_mole_conc(366,JY,JX));    H2PO4_irrig_mole_conc=0._r8
  allocate(CALQ(366,JY,JX));    CALQ=0._r8
  allocate(CFEQ(366,JY,JX));    CFEQ=0._r8
  allocate(CHYQ(366,JY,JX));    CHYQ=0._r8
  allocate(CCAQ(366,JY,JX));    CCAQ=0._r8
  allocate(CMGQ(366,JY,JX));    CMGQ=0._r8
  allocate(CNAQ(366,JY,JX));    CNAQ=0._r8
  allocate(CKAQ(366,JY,JX));    CKAQ=0._r8
  allocate(COHQ(366,JY,JX));    COHQ=0._r8
  allocate(CSOQ(366,JY,JX));    CSOQ=0._r8
  allocate(CCLQ(366,JY,JX));    CCLQ=0._r8
  allocate(CC3Q(366,JY,JX));    CC3Q=0._r8
  allocate(CHCQ(366,JY,JX));    CHCQ=0._r8
  allocate(CAL1Q(366,JY,JX));   CAL1Q=0._r8
  allocate(CAL2Q(366,JY,JX));   CAL2Q=0._r8
  allocate(CAL3Q(366,JY,JX));   CAL3Q=0._r8
  allocate(CAL4Q(366,JY,JX));   CAL4Q=0._r8
  allocate(CALSQ(366,JY,JX));   CALSQ=0._r8
  allocate(CFE1Q(366,JY,JX));   CFE1Q=0._r8
  allocate(CFE2Q(366,JY,JX));   CFE2Q=0._r8
  allocate(CFE3Q(366,JY,JX));   CFE3Q=0._r8
  allocate(CFE4Q(366,JY,JX));   CFE4Q=0._r8
  allocate(CFESQ(366,JY,JX));   CFESQ=0._r8
  allocate(CCAOQ(366,JY,JX));   CCAOQ=0._r8
  allocate(CCACQ(366,JY,JX));   CCACQ=0._r8
  allocate(CCAHQ(366,JY,JX));   CCAHQ=0._r8
  allocate(CCASQ(366,JY,JX));   CCASQ=0._r8
  allocate(CMGOQ(366,JY,JX));   CMGOQ=0._r8
  allocate(CMGCQ(366,JY,JX));   CMGCQ=0._r8
  allocate(CMGHQ(366,JY,JX));   CMGHQ=0._r8
  allocate(CMGSQ(366,JY,JX));   CMGSQ=0._r8
  allocate(CNACQ(366,JY,JX));   CNACQ=0._r8
  allocate(CNASQ(366,JY,JX));   CNASQ=0._r8
  allocate(CKASQ(366,JY,JX));   CKASQ=0._r8
  allocate(CH0PQ(366,JY,JX));   CH0PQ=0._r8
  allocate(HPO4_irrig_mole_conc(366,JY,JX));   HPO4_irrig_mole_conc=0._r8
  allocate(CH3PQ(366,JY,JX));   CH3PQ=0._r8
  allocate(CF1PQ(366,JY,JX));   CF1PQ=0._r8
  allocate(CF2PQ(366,JY,JX));   CF2PQ=0._r8
  allocate(CC0PQ(366,JY,JX));   CC0PQ=0._r8
  allocate(CC1PQ(366,JY,JX));   CC1PQ=0._r8
  allocate(CC2PQ(366,JY,JX));   CC2PQ=0._r8
  allocate(CM1PQ(366,JY,JX));   CM1PQ=0._r8
  allocate(CSTRQ(366,JY,JX));   CSTRQ=0._r8
  allocate(RRIG(24,366,JY,JX)); RRIG=0._r8
  allocate(WDPTH(366,JY,JX));   WDPTH=0._r8

  allocate(Lirri_auto_col(JY,JX)); Lirri_auto_col=.false.
  allocate(IIRRA(4,JY,JX));     IIRRA=0
  allocate(IrrigSubsurf_col(JY,JX));       IrrigSubsurf_col=0._r8
  allocate(IrrigSurface_col(JY,JX));       IrrigSurface_col=0._r8
  allocate(FIRRA_col(JY,JX));       FIRRA_col=0._r8
  allocate(CIRRA_col(JY,JX));       CIRRA_col=0._r8
  allocate(DIRRA(2,JY,JX));     DIRRA=0._r8
  allocate(TDIRI(12,JY,JX));    TDIRI=0._r8
  allocate(SurfIrrig_IonStrenth_col(JY,JX));       SurfIrrig_IonStrenth_col=0._r8

  allocate(trcn_irrig_vr(ids_nuts_beg:ids_nuts_end,JZ,JY,JX)); trcn_irrig_vr=0._r8
  allocate(trcSalt_Irrig_vr(idsalt_beg:idsalt_end,JZ,JY,JX));     trcSalt_irrig_vr=0._r8
  allocate(FWatIrrigate2MicP_vr(JZ,JY,JX));      FWatIrrigate2MicP_vr=0._r8
  allocate(HeatIrrigation_vr(JZ,JY,JX));    HeatIrrigation_vr=0._r8
  allocate(trcs_Irrig_flx_vr(ids_beg:ids_end,JZ,JY,JX));   trcs_Irrig_flx_vr=0._r8
  allocate(trcsalt_irrig_mole_conc_col(idsalt_beg:idsaltb_end,366,JY,JX));trcsalt_irrig_mole_conc_col=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructIrrigation
  use abortutils, only : destroy
  implicit none

  call destroy(trcsalt_irrig_mole_conc_col)
  call destroy(trcs_irrig_flx_col)
  call destroy(trcn_irrig_vr)
  call destroy(trcSalt_Irrig_vr)
  call destroy(trcs_Irrig_flx_vr)
  call destroy(PHQ)
  call destroy(NH4_irrig_mole_conc)
  call destroy(NO3_irrig_mole_conc)
  call destroy(H2PO4_irrig_mole_conc)
  call destroy(CALQ)
  call destroy(CFEQ)
  call destroy(CHYQ)
  call destroy(CCAQ)
  call destroy(CMGQ)
  call destroy(CNAQ)
  call destroy(CKAQ)
  call destroy(COHQ)
  call destroy(CSOQ)
  call destroy(CCLQ)
  call destroy(CC3Q)
  call destroy(CHCQ)
  call destroy(CAL1Q)
  call destroy(CAL2Q)
  call destroy(CAL3Q)
  call destroy(CAL4Q)
  call destroy(CALSQ)
  call destroy(CFE1Q)
  call destroy(CFE2Q)
  call destroy(CFE3Q)
  call destroy(CFE4Q)
  call destroy(CFESQ)
  call destroy(CCAOQ)
  call destroy(CCACQ)
  call destroy(CCAHQ)
  call destroy(CCASQ)
  call destroy(CMGOQ)
  call destroy(CMGCQ)
  call destroy(CMGHQ)
  call destroy(CMGSQ)
  call destroy(CNACQ)
  call destroy(CNASQ)
  call destroy(CKASQ)
  call destroy(CH0PQ)
  call destroy(HPO4_irrig_mole_conc)
  call destroy(CH3PQ)
  call destroy(CF1PQ)
  call destroy(CF2PQ)
  call destroy(CC0PQ)
  call destroy(CC1PQ)
  call destroy(CC2PQ)
  call destroy(CM1PQ)
  call destroy(CSTRQ)
  call destroy(RRIG)
  call destroy(WDPTH)
  call destroy(Lirri_auto_col)
  call destroy(IIRRA)
  call destroy(IrrigSubsurf_col)
  call destroy(IrrigSurface_col)
  call destroy(FIRRA_col)
  call destroy(CIRRA_col)
  call destroy(DIRRA)
  call destroy(TDIRI)
  call destroy(SurfIrrig_IonStrenth_col)
  call destroy(trcg_irrig_mole_conc_col)
  call destroy(trcSalt_irrig_vr)

  call destroy(FWatIrrigate2MicP_vr)
  call destroy(HeatIrrigation_vr)

  end subroutine DestructIrrigation

end module IrrigationDataType
