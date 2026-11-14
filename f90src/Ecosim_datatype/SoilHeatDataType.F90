module SoilHeatDatatype
  ! Description:
  ! Store native soil temperatures, capacities and accepted heat transfers for
  ! physical updates and column/landscape energy accounting.

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  TKS_vr(:,:,:)                      !soil temperature, [K]
  real(r8),target,allocatable ::  TLIceThawMicP_vr(:,:,:)            !hourly accumulated freeze-thaw flux in micropores, [m2 d-2 h-1]
  real(r8),target,allocatable ::  TLPhaseChangeHeat2Soi_vr(:,:,:)    !hourly accumulated freeze-thaw latent heat flux from soil, [MJ d-2 h-1]
  real(r8),target,allocatable ::  TLIceThawMacP_vr(:,:,:)            !hourly accumulated freeze-thaw flux in macropores, [m2 d-2 h-1]
  real(r8),target,allocatable ::  XPhaseChangeHeatL_snvr(:,:,:)      !hourly accumulated freeze-thaw latent heat flux from snow,[MJ d-2 h-1]
  real(r8),target,allocatable ::  VHeatCapacity_vr(:,:,:)            !soil heat capacity, [MJ m-3 K-1]
  real(r8),target,allocatable ::  TCS_vr(:,:,:)                      !soil temperature, [oC]
  real(r8),target,allocatable ::  HeatStore_col(:,:)                 !heat stored over the grid, including soil, litter and canopy, [MJ d-2]
  real(r8),target,allocatable ::  HeatSource_vr(:,:,:)               !heat source for warming,[MJ d-2 h-1]
  real(r8),target,allocatable ::  NumerSolidThermCond_vr(:,:,:)      !numerator for soil solid thermal conductivity, [MJ m h-1 K-1]
  real(r8),target,allocatable ::  DenomSolidThermCond_vr(:,:,:)      !denominator for soil solid thermal conductivity,[MJ K-1]
  real(r8),target,allocatable ::  HeatFlx2Grnd_col(:,:)              !heat flux into ground, computed from surface energy balance model, [MJ d-2 h-1]
  real(r8),target,allocatable ::  THeatFlowCellSoil_vr(:,:,:)        !hourly heat flux into soil layer,  [MJ m-3]
  real(r8),target,allocatable ::  HeatDrain_col(:,:)                 !heat loss through drainage [MJ d-2 h-1]
  real(r8),target,allocatable ::  HeatRunSurf_col(:,:)               !heat loss through surface runoff, [MJ d-2 h-1]
  real(r8),target,allocatable ::  HeatDischar_col(:,:)               !Heat loss through discharge, [MJ d-2 h-1]
  real(r8),target,allocatable ::  THeatFlow2Soil_col(:,:)            !Heat flow into colum, [MJ d-2 h-1]
  real(r8), allocatable :: HeatLateral2Cell_col(:,:) ! hour-integrated internal lateral heat, positive into column [MJ/cell]
  real(r8),target,allocatable ::  HeatSource_col(:,:)                !Heat source from heating, [MJ d-2 h-1]
  real(r8),target,allocatable ::  THeatSoiThaw_col(:,:)              !Heat associated with freeze-thaw, [MJ d-2 h-1]
  real(r8),target,allocatable ::  QSnoHeatXfer2Soil_col(:,:)         !Heat flux from snow into soil, [MJ d-2 h-1]
  real(r8),target,allocatable ::  QIceInflx_vr(:,:,:)                !Ice influx to layer, essential for pond/lake,  [m3 H2O d-2 h-1]

!----------------------------------------------------------------------

contains
  subroutine InitSoilHeatData

  implicit none

  allocate(QIceInflx_vr(JZ,JY,JX));  QIceInflx_vr = 0._r8
  allocate(QSnoHeatXfer2Soil_col(JY,JX)); QSnoHeatXfer2Soil_col = 0._r8
  allocate(HeatLateral2Cell_col(JY,JX)); HeatLateral2Cell_col=0._r8
  allocate(HeatSource_col(JY,JX)); HeatSource_col = 0._r8
  allocate(TKS_vr(0:JZ,JY,JX));    TKS_vr                              = 0._r8
  allocate(HeatFlx2Grnd_col(JY,JX));   HeatFlx2Grnd_col                = 0._r8
  allocate(HeatStore_col(JY,JX));  HeatStore_col                       = 0._r8
  allocate(TLIceThawMicP_vr(JZ,JY,JX));     TLIceThawMicP_vr                 = 0._r8
  allocate(TLPhaseChangeHeat2Soi_vr(JZ,JY,JX));    TLPhaseChangeHeat2Soi_vr  = 0._r8
  allocate(TLIceThawMacP_vr(JZ,JY,JX));    TLIceThawMacP_vr                  = 0._r8
  allocate(XPhaseChangeHeatL_snvr(JS,JY,JX));   XPhaseChangeHeatL_snvr = 0._r8
  allocate(VHeatCapacity_vr(0:JZ,JY,JX));   VHeatCapacity_vr           = 0._r8
  allocate(TCS_vr(0:JZ,JY,JX));    TCS_vr                         = 0._r8
  allocate(NumerSolidThermCond_vr(JZ,JY,JX));      NumerSolidThermCond_vr    = 0._r8
  allocate(DenomSolidThermCond_vr(JZ,JY,JX));      DenomSolidThermCond_vr    = 0._r8
  allocate(THeatFlowCellSoil_vr(JZ,JY,JX));    THeatFlowCellSoil_vr          = 0._r8
  allocate(HeatSource_vr(JZ,JY,JX)); HeatSource_vr=0._r8
  allocate(HeatDrain_col(JY,JX));  HeatDrain_col                       = 0._r8
  allocate(HeatRunSurf_col(JY,JX)); HeatRunSurf_col                    = 0._r8
  allocate(HeatDischar_col(JY,JX)); HeatDischar_col                    = 0._r8
  allocate(THeatFlow2Soil_col(JY,JX));THeatFlow2Soil_col               = 0._r8
  allocate(THeatSoiThaw_col(JY,JX)); THeatSoiThaw_col = 0._r8
  end subroutine InitSoilHeatData

!----------------------------------------------------------------------
  subroutine DestructSoilHeatData
  use abortutils, only : destroy
  implicit none

  call destroy(QIceInflx_vr)
  call destroy(QSnoHeatXfer2Soil_col)
  call destroy(THeatSoiThaw_col)
  call destroy(HeatSource_col)
  call destroy(TKS_vr)
  call destroy(HeatSource_vr)
  call destroy(TLIceThawMicP_vr)
  call destroy(TLPhaseChangeHeat2Soi_vr)
  call destroy(TLIceThawMacP_vr)
  call destroy(HeatStore_col)
  call destroy(XPhaseChangeHeatL_snvr)
  call destroy(VHeatCapacity_vr)
  call destroy(TCS_vr)
  call destroy(NumerSolidThermCond_vr)
  call destroy(DenomSolidThermCond_vr)
  call destroy(THeatFlowCellSoil_vr)  
  call destroy(HeatDrain_col)
  call destroy(HeatRunSurf_col)
  call destroy(HeatDischar_col)
  call destroy(HeatLateral2Cell_col)
  call destroy(THeatFlow2Soil_col)
  end subroutine DestructSoilHeatData

end module SoilHeatDatatype
