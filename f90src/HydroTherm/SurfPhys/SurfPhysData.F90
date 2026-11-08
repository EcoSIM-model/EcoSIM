module SurfPhysData
  ! Own surface-exchange work arrays, local resistances and area/time prefactors.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  save
  character(len=*), private, parameter :: mod_filename=&
  __FILE__

  real(r8),allocatable ::  XVLMobileWaterLitR_col(:,:)                   !water exceeds holding capacity of surface litter that is available for store on surface [m3 H2O/d2]
  real(r8),allocatable ::  XVLMobileWatMicP(:,:)                         !
  real(r8),allocatable ::  XVLiceMicP_col(:,:)                           !
  real(r8),allocatable ::  VLPoreLitR_col(:,:)                        !
  real(r8),allocatable ::  m3OM_col(:,:)                           !surface litter OM [m3 d-2]
  real(r8),allocatable ::  RawH2OLitr2CanAir_col(:,:)          !litter surface to canopy-air resistance [h/m]
  real(r8),allocatable ::  AScaledCdWOverSoil_col(:,:)             !area scaled conductance for latent heat flux over exposed soil [m^2 h]
  real(r8),allocatable ::  AScaledCdWOverLitr_col(:,:)             !area scaled conductance for latent heat flux over litter [m^2 h]
  real(r8),allocatable ::  RawIsoTLitr2CanAir_col(:,:) !neutral litter surface to canopy-air resistance [h/m]
  real(r8),allocatable ::  AScaledCdHOverLitr_col(:,:)             !area scaled conductance for sensible heat flux over litter [MJ h /(m K)]
  real(r8),allocatable ::  VapDiffusResistanceLitR(:,:)            !
  real(r8),allocatable ::  WatFLo2LitrM_col(:,:)                   !
  real(r8),allocatable ::  HeatFLoByWat2LitRM_col(:,:)             !  
  real(r8),allocatable ::  AScaledCdHOverSoil_col(:,:)             !area scaled conductance for sensible heat flux over exposed soil [MJ h /(m K)]
!  real(r8),allocatable :: VapXAir2TopLay(:,:)                     !water vapor flux from canopy air to top layer of soi/lake
  real(r8),allocatable ::  LitrIceFlxThaw_col(:,:)                 !Water flux from ice thaw in litter  (>0 thaw)
  real(r8),allocatable ::  LitrIceHeatFlxFrez_col(:,:)             !Heat associated with ice freeze in litter (<0 thaw) 
  real(r8),allocatable ::  ResistGrndGas2CanAir_col(:,:)                 !H2O-reference diagnostic of the ground-to-canopy-air gas network [h/m]
  real(r8),allocatable ::  FracEffAsLitR_col(:,:)                  !fraction of surface effectively coverd by litter [none]
  real(r8),allocatable ::  Rain2SoiMacP1_col(:,:)                      !
  real(r8),allocatable ::  PRECM_col(:,:)                          !
  real(r8),allocatable ::  Rain2SoiMicP1_col(:,:)                          !
  real(r8),allocatable ::  RainHeat2SoilP1_col(:,:)             !heat to soil from precipitation (+irrigation) [MJ/d2]
  real(r8),allocatable ::  RainHeat2LitR1_col(:,:)                        !
  real(r8),allocatable ::  Rain2LitR1_col(:,:)                          !
  real(r8),allocatable ::  FracAsExposedSoil_col(:,:)              !bare fraction of ground excluding litter/free-water cover, before snow partition [-]
!  real(r8),allocatable ::  HCNDLitr_col(:,:)                         !
  real(r8),allocatable :: TEvapXAir2Toplay_col(:,:)
  real(r8),allocatable :: TEvapXAir2LitR_col(:,:)
  real(r8),allocatable :: TEvapXAir2Snow_col(:,:)
  
  !real(r8),allocatable :: HeatFlowSno2SoiByWat_col(:,:)
  !real(r8),allocatable :: HeatConvFlxSno2Soi_col(:,:)
  !real(r8),allocatable :: HeatCndFlxSno2Soi_col(:,:)
  !real(r8),allocatable :: HeatConvFlxLitr2Soi_col(:,:)
  !real(r8),allocatable :: HeatCndFlxLitr2Soi_col(:,:)

  real(r8),allocatable :: HeatSensAir2Grnd_col(:,:)
  real(r8),allocatable :: Radnet2Grnd_col(:,:)
  real(r8),allocatable :: LatentHeatEvapAir2Grnd_col(:,:)
  real(r8),allocatable :: HeatSensVapAir2Grnd_col(:,:)
  
  real(r8),allocatable :: CdLitRHSens_col(:,:)            !conductance of sensible heat [MJ/K]
  real(r8),allocatable :: watflw(:,:)
  real(r8),allocatable :: waticefl(:,:)
  real(r8),allocatable :: CdSoiHSens_col(:,:)     ! scaled conductance for sensible heat [MJ/K]
  real(r8),allocatable :: HeatSnow2LitrM_col(:,:) !heat flux from snow to litter [MJ/d-2]
  real(r8),allocatable :: HeatAir2LitrM_col(:,:)  !heat flux from air to litter [MJ/d-2]
  real(r8),allocatable :: HeatSoil2LitrM_col(:,:) !heat flux from soil to litter [MJ/d-2]
  public :: InitSurfPhysData,DestructSurfPhysData

  contains
!------------------------------------------------------------------------------------------  

  subroutine InitSurfPhysData
  ! Allocate and initialize surface-exchange work arrays.
  implicit none

  allocate(HeatSoil2LitrM_col(JY,JX)); HeatSoil2LitrM_col=0._r8
  allocate(watflw(JY, JX))
  allocate(waticefl(JY,JX))
  allocate(m3OM_col(JY,JX))
  allocate(XVLMobileWaterLitR_col(JY,JX));       XVLMobileWaterLitR_col=0._r8
  allocate(XVLMobileWatMicP(JY,JX));       XVLMobileWatMicP=0._r8
  allocate(XVLiceMicP_col(JY,JX));       XVLiceMicP_col=0._r8
  allocate(VLPoreLitR_col(JY,JX));      VLPoreLitR_col=0._r8
  allocate(CdLitRHSens_col(JY,JX)); CdLitRHSens_col=0._r8
  allocate(RawH2OLitr2CanAir_col(JY,JX));        RawH2OLitr2CanAir_col=0._r8
  allocate(AScaledCdWOverSoil_col(JY,JX));       AScaledCdWOverSoil_col=0._r8
  allocate(AScaledCdWOverLitr_col(JY,JX));       AScaledCdWOverLitr_col=0._r8
  allocate(RawIsoTLitr2CanAir_col(JY,JX));        RawIsoTLitr2CanAir_col=0._r8
  allocate(AScaledCdHOverLitr_col(JY,JX));       AScaledCdHOverLitr_col=0._r8
  allocate(VapDiffusResistanceLitR(JY,JX));         VapDiffusResistanceLitR=0._r8  
  allocate(WatFLo2LitrM_col(JY,JX));       WatFLo2LitrM_col=0._r8
  allocate(HeatFLoByWat2LitRM_col(JY,JX));      HeatFLoByWat2LitRM_col=0._r8  
  allocate(AScaledCdHOverSoil_col(JY,JX));       AScaledCdHOverSoil_col=0._r8  
!  allocate(VapXAir2TopLay(JY,JX));       VapXAir2TopLay=0._r8  
  allocate(LitrIceFlxThaw_col(JY,JX));       LitrIceFlxThaw_col         = 0._r8
  allocate(ResistGrndGas2CanAir_col(JY,JX));        ResistGrndGas2CanAir_col                            = 0._r8
  allocate(LitrIceHeatFlxFrez_col(JY,JX));       LitrIceHeatFlxFrez_col = 0._r8
  allocate(FracEffAsLitR_col(JY,JX));       FracEffAsLitR_col                           = 0._r8
  allocate(Rain2SoiMacP1_col(JY,JX));        Rain2SoiMacP1_col          = 0._r8
  allocate(PRECM_col(JY,JX));       PRECM_col                   = 0._r8
  allocate(Rain2SoiMicP1_col(JY,JX));        Rain2SoiMicP1_col          = 0._r8
  allocate(RainHeat2SoilP1_col(JY,JX));      RainHeat2SoilP1_col    = 0._r8
  allocate(RainHeat2LitR1_col(JY,JX));      RainHeat2LitR1_col          = 0._r8
  allocate(Rain2LitR1_col(JY,JX));        Rain2LitR1_col                = 0._r8
  allocate(FracAsExposedSoil_col(JY,JX));       FracAsExposedSoil_col                           = 0._r8
  allocate(TEvapXAir2Toplay_col(JY,JX)); TEvapXAir2Toplay_col   = 0._r8
  allocate(TEvapXAir2LitR_col(JY,JX)); TEvapXAir2LitR_col       = 0._r8
  allocate(TEvapXAir2Snow_col(JY,JX)); TEvapXAir2Snow_col       = 0._r8
  allocate(HeatSensAir2Grnd_col(JY,JX));    HeatSensAir2Grnd_col=0._r8
  allocate(Radnet2Grnd_col(JY,JX));         Radnet2Grnd_col=0._r8
  allocate(LatentHeatEvapAir2Grnd_col(JY,JX)); LatentHeatEvapAir2Grnd_col=0._r8
  allocate(HeatSensVapAir2Grnd_col(JY,JX)); HeatSensVapAir2Grnd_col=0._r8
  allocate(HeatSnow2LitrM_col(JY,JX));     HeatSnow2LitrM_col=0._r8
  allocate(HeatAir2LitrM_col(JY,JX));      HeatAir2LitrM_col=0._r8
  allocate(CdSoiHSens_col(JY,JX));         CdSoiHSens_col=0._r8
  end subroutine InitSurfPhysData  
!------------------------------------------------------------------------------------------  

  subroutine DestructSurfPhysData
  ! Release surface-exchange work arrays.
  use abortutils, only : destroy
  implicit none

  call destroy(HeatSoil2LitrM_col)
  call destroy(HeatAir2LitrM_col)
  call destroy(HeatSnow2LitrM_col)
  call destroy(CdSoiHSens_col)
  call destroy(CdLitRHSens_col)
  call destroy(m3OM_col)
  call destroy(TEvapXAir2Toplay_col)
  call destroy(TEvapXAir2LitR_col) 
  call destroy(TEvapXAir2Snow_col) 
  call destroy(XVLMobileWaterLitR_col)
  call destroy(XVLMobileWatMicP)
  call destroy(XVLiceMicP_col)
  call destroy(VLPoreLitR_col)  
  call destroy(RawH2OLitr2CanAir_col)
  call destroy(AScaledCdWOverSoil_col)
  call destroy(AScaledCdWOverLitr_col)
  call destroy(RawIsoTLitr2CanAir_col)
  call destroy(AScaledCdHOverLitr_col)
  call destroy(VapDiffusResistanceLitR)  
  call destroy(WatFLo2LitrM_col)
  call destroy(HeatFLoByWat2LitRM_col)  
  call destroy(AScaledCdHOverSoil_col)  
!  call destroy(VapXAir2TopLay)
  call destroy(LitrIceFlxThaw_col)  
  call destroy(ResistGrndGas2CanAir_col)
  call destroy(LitrIceHeatFlxFrez_col)
  call destroy(FracEffAsLitR_col)
  call destroy(Rain2SoiMacP1_col)
  call destroy(PRECM_col)
  call destroy(Rain2SoiMicP1_col)
  call destroy(RainHeat2SoilP1_col)
  call destroy(RainHeat2LitR1_col)
  call destroy(Rain2LitR1_col)
  call destroy(FracAsExposedSoil_col)
!  call destroy(HCNDLitr_col)  
  call destroy(HeatSensAir2Grnd_col)
  call destroy(Radnet2Grnd_col)
  call destroy(LatentHeatEvapAir2Grnd_col)
  call destroy(HeatSensVapAir2Grnd_col)
  end subroutine DestructSurfPhysData

end module SurfPhysData
