module SnowPhysData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod

  implicit none

  character(len=*), private, parameter :: mod_filename=&
  __FILE__

  real(r8),allocatable ::  trcg_AquaAdv_NetFlx_snvr(:,:,:,:)          !dissolved volatile advection in snow [g d-2]
  real(r8),allocatable ::  trcn_AquaAdv_NetFlx_snvr(:,:,:,:)          !dissolved nutrient advection in snow [g d-2]
  real(r8),allocatable ::  trcSalt_AquaAdv_NetFlx_snvr(:,:,:,:)       !dissolved salt advection in snow [mol d-2]
  real(r8),allocatable ::  TDrysnoByRedist_col(:,:)                   !total dry snow moved by redistribution in grid [m3 d-2]
  real(r8),allocatable ::  trcSalt_LossXSnowRedist_col(:,:,:)         
  real(r8),allocatable ::  trcg_LossXSnowRedist_col(:,:,:)
  real(r8),allocatable ::  trcn_LossXSnowRedist_col(:,:,:)
  real(r8),allocatable ::  XSnowThawMassL_snvr(:,:,:)                !hourly convective heat flux from snow transfer
  real(r8),allocatable ::  XIceThawMassL_snvr(:,:,:)                 !hourly convective heat flux from ice transfer  
  real(r8),allocatable ::  THeatSnoByRedist_col(:,:)                 !total heat moved by snow redistribution in grid [MJ d-2] 
  real(r8),allocatable ::  TIceSnoByRedist_col(:,:)                  !total ice volume moved by snow redistribution in grid [m3 d-2]
  real(r8),allocatable ::  TWatSnoByRedist_col(:,:)                  !total water volume moved by snow redistribution in grid [m3 -d2]
  real(r8),allocatable ::  cumWatFlx2LitRByRunoff_col(:,:)           !
  real(r8),allocatable ::  cumHeatFlx2LitRByRunoff_col(:,:)          !
  real(r8),allocatable ::  cumDrySnoByRedistM_col(:,:)               !cumulative dry snow added to grid by redistribution [m3 d-2]
  real(r8),allocatable ::  cumWatSnoByRedistM_col(:,:)               !cumulative water in snow added to grid by redistribution [m3 d-2]
  real(r8),allocatable ::  cumIceSnoByRedistM_col(:,:)               !cumulative ice in snow added to grid by redistribution [m3 d-2]
  real(r8),allocatable ::  cumHeatSnoByRedistM_col(:,:)              !cumulative heat in snow added to grid by redistribution [MJ d-2]
  real(r8),allocatable ::  TKSnow1_snvr(:,:,:)                       !temporary snow layer temperature [K]
  real(r8),allocatable ::  VLHeatCapSnowM1_snvr(:,:,:)               !temporary snow layer volumetric heat capacity [MJ/K]
  real(r8),allocatable ::  VLIceSnow0M_snvr(:,:,:)                   !
  real(r8),allocatable ::  VLWatSnow0M_snvr(:,:,:)                   !snow-held water volume during iteration
  real(r8),allocatable ::  CumSno2SnowLM_snvr(:,:,:)
  real(r8),allocatable ::  CumWat2SnowLM_snvr(:,:,:)
  real(r8),allocatable ::  CumIce2SnowLM_snvr(:,:,:)
  real(r8),allocatable ::  CumHeat2SnowLM_snvr(:,:,:)
  real(r8),allocatable ::  XSnowThawMassLM_snvr(:,:,:)
  real(r8),allocatable ::  XIceThawMassLM_snvr(:,:,:)
  real(r8),allocatable ::  tEnGYM_snvr(:,:,:)  
  real(r8),allocatable ::  CumSno2SnowL_snvr(:,:,:)                !
  real(r8),allocatable ::  CumWat2SnowL_snvr(:,:,:)                !
  real(r8),allocatable ::  CumIce2SnowL_snvr(:,:,:)                !
  real(r8),allocatable ::  CumHeat2SnowL_snvr(:,:,:)               !cumulative heat into snow layer [MJ d-2]
  real(r8),allocatable ::  XPhaseChangeHeatLM_snvr(:,:,:)          !heat released due to ice/liquid phase change [MJ d-2]
  real(r8),allocatable ::  WatX2SnoLay_snvr(:,:,:)                 !liquid water flux into snow layer [m3 d-2]
  real(r8),allocatable ::  SnoX2SnoLay_snvr(:,:,:)                 !dry snow flux into snow layer [m3 d-2]
  real(r8),allocatable ::  IceX2SnoLay_snvr(:,:,:)                 !ice flux into snow layer  [m3 d-2]
  real(r8),allocatable ::  HeatX2SnoLay_snvr(:,:,:)                !heat flux into snow layer [MJ d-2]
  real(r8),allocatable ::  VLDrySnoWE0M_snvr(:,:,:)                !dry snow layer volume during iteration [m3 d-2]
  real(r8),allocatable ::  VLDrySnoWE0_snvr(:,:,:)                 !dry snow layer volume before update [m3 d-2]
  real(r8),allocatable ::  VLIceSnow0_snvr(:,:,:)                  !snow-held ice layer volume before update [m3 d-2]
  real(r8),allocatable ::  VLWatSnow0_snvr(:,:,:)                  !snow-held water layer volume before update [m3 d-2]
  real(r8),allocatable ::  VLSnoDWI1_snvr(:,:,:)                   !total snow volume during iteration [m3 d-2] 
  real(r8),allocatable ::  SnowThickL0_snvr(:,:,:)                 !snow thickness of each layer [m]
  real(r8),allocatable ::  WatSnoFlxByRedistM_2DH(:,:,:,:,:)       !water snow due to redistribution at iteration M [m3 d-2]
  real(r8),allocatable ::  IceSnoFlxByRedistM_2DH(:,:,:,:,:)       !Ice snow due to redistribution at iteration M [m3 d-2]
  real(r8),allocatable ::  HeatSnoFlxByRedistM_2DH(:,:,:,:,:)      !Heat snow due to redistribution at iteration M [m3 d-2]
  
  real(r8),allocatable :: HeatFlowSno2SoiByWat_col(:,:)
  real(r8),allocatable :: HAltFlowSno2SoiByWat_col(:,:)
  real(r8),allocatable :: HeatConvFlxSno2Soi_col(:,:)
  real(r8),allocatable :: HeatCndFlxSno2Soi_col(:,:)
  real(r8),allocatable :: HeatConvFlxLitr2Soi_col(:,:)
  real(r8),allocatable :: HeatCndFlxLitr2Soi_col(:,:)

  public :: InitSnowPhysData
  public :: DestructSnowPhysData
  contains
!----------------------------------------------------------------------  
  subroutine  InitSnowPhysData
  implicit none

  allocate(WatSnoFlxByRedistM_2DH(60,2,2,JY,JX)); WatSnoFlxByRedistM_2DH=0._r8
  allocate(IceSnoFlxByRedistM_2DH(60,2,2,JY,JX)); IceSnoFlxByRedistM_2DH=0._r8
  allocate(HeatSnoFlxByRedistM_2DH(60,2,2,JY,JX)); HeatSnoFlxByRedistM_2DH=0._r8

  allocate(trcg_AquaAdv_NetFlx_snvr(idg_beg:idg_NH3,JS,JY,JX));        trcg_AquaAdv_NetFlx_snvr=0._r8
  allocate(trcn_AquaAdv_NetFlx_snvr(ids_nut_beg:ids_nuts_end,JS,JY,JX)); trcn_AquaAdv_NetFlx_snvr=0._r8
  allocate(trcSalt_AquaAdv_NetFlx_snvr(idsalt_beg:idsalt_end,JS,JY,JX));       trcSalt_AquaAdv_NetFlx_snvr=0._r8  
  allocate(trcSalt_LossXSnowRedist_col(idsalt_beg:idsalt_end,JY,JX));           trcSalt_LossXSnowRedist_col=0._r8  

  allocate(trcg_LossXSnowRedist_col(idg_beg:idg_NH3,JY,JX));            trcg_LossXSnowRedist_col=0._r8
  allocate(trcn_LossXSnowRedist_col(ids_nut_beg:ids_nuts_end,JY,JX));trcn_LossXSnowRedist_col=0._r8

  allocate(TIceSnoByRedist_col(JY,JX));         TIceSnoByRedist_col=0._r8
  allocate(TWatSnoByRedist_col(JY,JX));         TWatSnoByRedist_col=0._r8
  allocate(TDrysnoByRedist_col(JY,JX));         TDrysnoByRedist_col=0._r8
  allocate(THeatSnoByRedist_col(JY,JX));        THeatSnoByRedist_col=0._r8
  allocate(TKSnow1_snvr(JS,JY,JX));     TKSnow1_snvr=0._r8
  allocate(cumHeatFlx2LitRByRunoff_col(JY,JX));       cumHeatFlx2LitRByRunoff_col=0._r8
  allocate(cumHeatSnoByRedistM_col(JY,JX));       cumHeatSnoByRedistM_col=0._r8
  allocate(cumWatFlx2LitRByRunoff_col(JY,JX));        cumWatFlx2LitRByRunoff_col=0._r8
  allocate(cumDrySnoByRedistM_col(JY,JX));        cumDrySnoByRedistM_col=0._r8
  allocate(cumWatSnoByRedistM_col(JY,JX));        cumWatSnoByRedistM_col=0._r8
  allocate(cumIceSnoByRedistM_col(JY,JX));        cumIceSnoByRedistM_col=0._r8
  allocate(VLHeatCapSnowM1_snvr(JS,JY,JX));  VLHeatCapSnowM1_snvr=0._r8
  allocate(VLIceSnow0M_snvr(JS,JY,JX));   VLIceSnow0M_snvr=0._r8
  allocate(VLWatSnow0M_snvr(JS,JY,JX));   VLWatSnow0M_snvr=0._r8
  allocate(CumSno2SnowL_snvr(JS,JY,JX));    CumSno2SnowL_snvr=0._r8
  allocate(tEnGYM_snvr(JS,JY,JX));  tEnGYM_snvr=0._r8
  allocate(CumSno2SnowLM_snvr(JS,JY,JX)); CumSno2SnowLM_snvr=0._r8
  allocate(CumWat2SnowLM_snvr(JS,JY,JX)); CumWat2SnowLM_snvr=0._r8
  allocate(CumIce2SnowLM_snvr(JS,JY,JX)); CumIce2SnowLM_snvr=0._r8
  allocate(CumHeat2SnowLM_snvr(JS,JY,JX)); CumHeat2SnowLM_snvr=0._r8
  allocate(XSnowThawMassLM_snvr(JS,JY,JX)); XSnowThawMassLM_snvr=0._r8
  allocate(XIceThawMassLM_snvr(JS,JY,JX)); XIceThawMassLM_snvr=0._r8

  allocate(CumWat2SnowL_snvr(JS,JY,JX));    CumWat2SnowL_snvr=0._r8
  allocate(CumIce2SnowL_snvr(JS,JY,JX));    CumIce2SnowL_snvr=0._r8
  allocate(CumHeat2SnowL_snvr(JS,JY,JX));   CumHeat2SnowL_snvr=0._r8
  allocate(XPhaseChangeHeatLM_snvr(JS,JY,JX));    XPhaseChangeHeatLM_snvr=0._r8
  allocate(WatX2SnoLay_snvr(JS,JY,JX));    WatX2SnoLay_snvr=0._r8
  allocate(SnoX2SnoLay_snvr(JS,JY,JX));    SnoX2SnoLay_snvr=0._r8
  allocate(IceX2SnoLay_snvr(JS,JY,JX));    IceX2SnoLay_snvr=0._r8
  allocate(HeatX2SnoLay_snvr(JS,JY,JX));   HeatX2SnoLay_snvr=0._r8
  allocate(VLDrySnoWE0M_snvr(JS,JY,JX));   VLDrySnoWE0M_snvr=0._r8
  allocate(VLDrySnoWE0_snvr(JS,JY,JX));    VLDrySnoWE0_snvr=0._r8
  allocate(VLIceSnow0_snvr(JS,JY,JX));    VLIceSnow0_snvr=0._r8
  allocate(VLWatSnow0_snvr(JS,JY,JX));    VLWatSnow0_snvr=0._r8
  allocate(VLSnoDWI1_snvr(JS,JY,JX));    VLSnoDWI1_snvr=0._r8
  allocate(SnowThickL0_snvr(JS,JY,JX));   SnowThickL0_snvr=0._r8
  allocate(XSnowThawMassL_snvr(JS,JY,JX));   XSnowThawMassL_snvr=0._r8
  allocate(XIceThawMassL_snvr(JS,JY,JX));   XIceThawMassL_snvr=0._r8  
  
  allocate(HeatFlowSno2SoiByWat_col(JY,JX))
  allocate(HAltFlowSno2SoiByWat_col(JY,JX))
  allocate(HeatConvFlxSno2Soi_col(JY,JX))
  allocate(HeatCndFlxSno2Soi_col(JY,JX))
  allocate(HeatConvFlxLitr2Soi_col(JY,JX))
  allocate(HeatCndFlxLitr2Soi_col(JY,JX))
  
  end subroutine  InitSnowPhysData
!----------------------------------------------------------------------  
  subroutine DestructSnowPhysData
  use abortutils, only : destroy
  implicit none

  call destroy(trcg_AquaAdv_NetFlx_snvr)
  call destroy(trcn_AquaAdv_NetFlx_snvr)
  call destroy(trcSalt_AquaAdv_NetFlx_snvr)
  call destroy(tEnGYM_snvr)
  call destroy(trcg_LossXSnowRedist_col)
  call destroy(trcn_LossXSnowRedist_col)  
  call destroy(WatSnoFlxByRedistM_2DH)
  call destroy(IceSnoFlxByRedistM_2DH)
  call destroy(HeatSnoFlxByRedistM_2DH)
  call destroy(TWatSnoByRedist_col)  
  call destroy(cumWatFlx2LitRByRunoff_col)
  call destroy(cumDrySnoByRedistM_col)
  call destroy(cumWatSnoByRedistM_col)  
  call destroy(cumIceSnoByRedistM_col)  
  call destroy(cumHeatSnoByRedistM_col)  
  call destroy(TIceSnoByRedist_col)  
  call destroy(THeatSnoByRedist_col)  
  call destroy(TKSnow1_snvr)  
  call destroy(VLHeatCapSnowM1_snvr)
  call destroy(VLIceSnow0M_snvr)
  call destroy(VLWatSnow0M_snvr)
  call destroy(CumSno2SnowL_snvr)
  call destroy(CumSno2SnowLM_snvr)
  call destroy(CumWat2SnowLM_snvr)
  call destroy(CumIce2SnowLM_snvr)
  call destroy(CumHeat2SnowLM_snvr)
  call destroy(XSnowThawMassLM_snvr)
  call destroy(XIceThawMassLM_snvr)
  call destroy(CumWat2SnowL_snvr)
  call destroy(CumIce2SnowL_snvr)
  call destroy(CumHeat2SnowL_snvr)
  call destroy(XPhaseChangeHeatLM_snvr)
  call destroy(WatX2SnoLay_snvr)
  call destroy(SnoX2SnoLay_snvr)
  call destroy(IceX2SnoLay_snvr)
  call destroy(HeatX2SnoLay_snvr)
  call destroy(VLDrySnoWE0M_snvr)
  call destroy(VLDrySnoWE0_snvr)
  call destroy(VLIceSnow0_snvr)
  call destroy(VLWatSnow0_snvr)
  call destroy(VLSnoDWI1_snvr)
  call destroy(SnowThickL0_snvr)
  call destroy(XSnowThawMassL_snvr)  
  call destroy(trcSalt_LossXSnowRedist_col)
  call destroy(cumHeatFlx2LitRByRunoff_col)
  call destroy(TDrysnoByRedist_col)
  call destroy(XIceThawMassL_snvr)  
  
  call destroy(HeatFlowSno2SoiByWat_col)
  call destroy(HAltFlowSno2SoiByWat_col)
  call destroy(HeatConvFlxSno2Soi_col)
  call destroy(HeatCndFlxSno2Soi_col)
  call destroy(HeatConvFlxLitr2Soi_col)
  call destroy(HeatCndFlxLitr2Soi_col)
  
  end subroutine DestructSnowPhysData
end module SnowPhysData
