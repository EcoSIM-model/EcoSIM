module HydroThermData
  ! Own substep hydrological/thermal work arrays and local surface resistances.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8),allocatable ::  LWEmscefLitR_col(:,:)                         !  
  real(r8),allocatable ::  LWRad2LitR_col(:,:)                         !
  real(r8),allocatable ::  LWEmscefSoil_col(:,:)                         !
  real(r8),allocatable ::  LWRad2Soil_col(:,:)                         !
  real(r8),allocatable ::  RadSW2LitR_col(:,:)                         !
  real(r8),allocatable ::  RadSW2Soil_col(:,:)                         !
  real(r8),allocatable ::  LWEmscefSnow_col(:,:)                 !  
  real(r8),allocatable ::  LWRad2Snow_col(:,:)                   !    
  real(r8),allocatable ::  FracSoiPAsWat_vr(:,:,:)               !
  real(r8),allocatable ::  PSISM1_vr(:,:,:)                      !intermediate time step soil matric pressure, [MPa]  
  real(r8),allocatable ::  TKSoil1_vr(:,:,:)                     !intermeediate time step soil temeprature, [K]  
  real(r8),allocatable ::  DLYRR_COL(:,:)                        !  
  real(r8),allocatable ::  FracSoiPAsIce_vr(:,:,:)               !  
  real(r8),allocatable ::  FracAirFilledSoilPore_vr(:,:,:)               !  
  real(r8),allocatable ::  VHeatCapacity1_vr(:,:,:)              !whole layer heat capacity (not snow)
  real(r8),allocatable ::  FracSoilAsAirt(:,:,:)                 !fraction of soil volume as air, normalized using current pore volume
  real(r8),allocatable ::  VLWatMicP1_vr(:,:,:)                  !    
  real(r8),allocatable ::  VLiceMicP1_vr(:,:,:)                  !  
  real(r8),allocatable ::  SnowFallt_col(:,:)                    !  snowfall to snow, [m3 d-2]
  real(r8),allocatable ::  Rain2Snowt_col(:,:)                   !  rainfall to snow, [m3 d-2]
  real(r8),allocatable ::  VLairMacP1_vr(:,:,:)                  !  positively corrected macropore volume, [m^3/d^2]
  real(r8),allocatable ::  EVAPW_col(:,:)                        !  evaporation between atmosphere and surface snow [m3 H2O d-2/h]
  real(r8),allocatable ::  EVAPS_col(:,:)                        !  sublimation from surface snow into atmosphere, [m3 H2O d-2/h]
  real(r8),allocatable ::  VapXAir2Sno_col(:,:)                  ! air-snow exchange of water vapor through evaporation/condensation, sublimation/deposition
  real(r8),allocatable ::  SoilFracAsMacP1_vr(:,:,:)             !  
  real(r8),allocatable ::  PrecHeat2Snowt_col(:,:)               ! total heat flux to snow due to precipitation, [m3 d-2]

  real(r8),allocatable ::  AScaledCdHOverSnow_col(:,:)           ! area scaled sensible heat flux conductance over snow  [MJ h /(m K)]
  real(r8),allocatable ::  Ice2Snowt_col(:,:)                    ! icefall to snow, [m3 d-2]

  real(r8),allocatable :: RawH2OSnow2CanAir_col(:,:) !snow to canopy-air resistance [h/m]
  real(r8),allocatable :: RawH2OSoil2CanAir_col(:,:) !bare soil to canopy-air resistance [h/m]
  real(r8),allocatable :: SurfSensHeatConduct_col(:,:) !sum of active soil/litter/snow sensible conductances over this hour [MJ/K]
  real(r8),allocatable :: SurfLatHeatVapDeriv_col(:,:) !sum of active, water-limited latent-energy derivatives over this hour [MJ m3 ton-1]
  real(r8),allocatable ::  VLairMicP1_vr(:,:,:)                  ! corrected air-filled micropore volume [m3/d2]
  real(r8),allocatable ::  TKSnow0_snvr(:,:,:)                   !  
  real(r8),allocatable ::  VLWatMacP1_vr(:,:,:)                  !
  real(r8),allocatable ::  SoilFracAsMicP_vr(:,:,:)              !
  real(r8),allocatable ::  VLHeatCapacityA_vr(:,:,:)             !
  real(r8),allocatable ::  VLHeatCapacityB_vr(:,:,:)             !
  real(r8),allocatable ::  RawIsoTSoil2CanAir_col(:,:)        !neutral soil-to-canopy-air resistance [h/m]
  real(r8),allocatable ::  AScaledCdWOverSnow_col(:,:)           ! area scaled latent heat flux conductance over snow, [m^2 h]
  real(r8),allocatable ::  Altitude_col(:,:)                    ! grid altitude [m]
  real(r8),allocatable ::  VLiceMacP1_vr(:,:,:)                      !
  real(r8),allocatable ::  RadSW2Sno_col(:,:)                   !short wave radiation on snow [MJ]  
  real(r8),allocatable ::  WatFlx2LitRByRunoff_2DH(:,:,:,:)      !surface runoff flux  
  real(r8),allocatable ::  HeatFlx2LitRByRunoff_2DH(:,:,:,:)                      !  
  real(r8),allocatable ::  HeatFlow2Soili_3D(:,:,:,:)                     !  
  real(r8),allocatable ::  WaterFlow2Micpt_3D(:,:,:,:)                    !Micropore flow into the grid cell along direction N at iteration M [m3 H2O dt-1]
  real(r8),allocatable ::  WaterFlow2Macpt_3D(:,:,:,:)                    !Macropore flow into the grid cell along direction N at iteration M [m3 H2O dt-1] 
  real(r8),allocatable ::  WaterFlow2Micptl_3D(:,:,:,:)                   !temporarily copy of Micropore flow into the grid cell along direction N  [m3 H2O dt-1]    
  real(r8),allocatable ::  WaterFlow2Macptl_3D(:,:,:,:)                   !temporarily copy of Macropore flow into the grid cell along direction N [m3 H2O dt-1]     
  real(r8),allocatable ::  WaterFlow2MicptX_3D(:,:,:,:)                   !before saturation correction Micropore flow into the grid cell along direction N  [m3 H2O dt-1]  
  real(r8),allocatable ::  VLWatMicP2_vr(:,:,:)                           !liquid water volume 
  real(r8),allocatable ::  VLairMicP_vr(:,:,:)                            !gas volume content
  real(r8),allocatable ::  VLWatMicPX1_vr(:,:,:)                          !micropore water volume behind wetting front
  real(r8),allocatable ::  Qinflx2SoilM_col(:,:)                          !infiltration flux at iteration M
  real(r8),save :: HeatAdv_scal = 1._r8
  public :: InitHydroThermData
  public :: DestructHydroThermData
  contains

!------------------------------------------------------------------------------------------
  subroutine InitHydroThermData
  ! Allocate and initialize substep state and surface-resistance work arrays.
  implicit none


  allocate(FracSoiPAsWat_vr(0:JZ,JY,JX)); FracSoiPAsWat_vr=0._r8
  allocate(PSISM1_vr(0:JZ,JY,JX)); PSISM1_vr=0._r8
  allocate(TKSoil1_vr(0:JZ,JY,JX));    TKSoil1_vr=0._r8    
  allocate(DLYRR_COL(JY,JX));       DLYRR_col=0._r8  
  allocate(FracSoiPAsIce_vr(0:JZ,JY,JX)); FracSoiPAsIce_vr=0._r8 
  allocate(FracAirFilledSoilPore_vr(0:JZ,JY,JX)); FracAirFilledSoilPore_vr=0._r8   
  allocate(VHeatCapacity1_vr(0:JZ,JY,JX));  VHeatCapacity1_vr=0._r8  
  allocate(FracSoilAsAirt(0:JZ,JY,JX)); FracSoilAsAirt=0._r8  
  allocate(VLWatMicP1_vr(0:JZ,JY,JX));  VLWatMicP1_vr=0._r8  
  allocate(VLiceMicP1_vr(0:JZ,JY,JX));  VLiceMicP1_vr=0._r8
  allocate(SnowFallt_col(JY,JX));       SnowFallt_col=0._r8
  allocate(Rain2Snowt_col(JY,JX));       Rain2Snowt_col=0._r8
  allocate(VLairMacP1_vr(JZ,JY,JX));   VLairMacP1_vr=0._r8 
  allocate(VapXAir2Sno_col(JY,JX));      VapXAir2Sno_col=0._r8  
  allocate(EVAPW_col(JY,JX));       EVAPW_col=0._r8    
  allocate(EVAPS_col(JY,JX));       EVAPS_col=0._r8  
  allocate(SoilFracAsMacP1_vr(JZ,JY,JX));     SoilFracAsMacP1_vr=0._r8  
  allocate(PrecHeat2Snowt_col(JY,JX));      PrecHeat2Snowt_col=0._r8  

  allocate(AScaledCdHOverSnow_col(JY,JX));       AScaledCdHOverSnow_col=0._r8
  allocate(Ice2Snowt_col(JY,JX));       Ice2Snowt_col=0._r8    

  allocate(LWEmscefSnow_col(JY,JX));       LWEmscefSnow_col=0._r8  
  allocate(RawH2OSoil2CanAir_col(JY,JX)); RawH2OSoil2CanAir_col=0._r8
  allocate(RawH2OSnow2CanAir_col(JY,JX));        RawH2OSnow2CanAir_col=0._r8
  allocate(SurfSensHeatConduct_col(JY,JX)); SurfSensHeatConduct_col=0._r8
  allocate(SurfLatHeatVapDeriv_col(JY,JX)); SurfLatHeatVapDeriv_col=0._r8
  allocate(LWRad2Snow_col(JY,JX));       LWRad2Snow_col=0._r8  
  allocate(VLairMicP1_vr(0:JZ,JY,JX));  VLairMicP1_vr=0._r8  
  allocate(TKSnow0_snvr(JS,JY,JX));      TKSnow0_snvr=0._r8  
  allocate(VLWatMacP1_vr(JZ,JY,JX));   VLWatMacP1_vr=0._r8
  allocate(SoilFracAsMicP_vr(JZ,JY,JX));     SoilFracAsMicP_vr=0._r8
  allocate(VLHeatCapacityA_vr(JZ,JY,JX));   VLHeatCapacityA_vr=0._r8
  allocate(VLHeatCapacityB_vr(JZ,JY,JX));   VLHeatCapacityB_vr=0._r8  
  allocate(RawIsoTSoil2CanAir_col(JY,JX));         RawIsoTSoil2CanAir_col=0._r8
  allocate(AScaledCdWOverSnow_col(JY,JX));       AScaledCdWOverSnow_col=0._r8
  allocate(Altitude_col(JY,JX));        Altitude_col=0._r8  
  allocate(VLiceMacP1_vr(JZ,JY,JX));   VLiceMacP1_vr=0._r8  
  allocate(RadSW2Sno_col(JY,JX));       RadSW2Sno_col=0._r8  
  allocate(WatFlx2LitRByRunoff_2DH(2,2,JV,JH));     WatFlx2LitRByRunoff_2DH=0._r8  
  allocate(HeatFlx2LitRByRunoff_2DH(2,2,JV,JH));    HeatFlx2LitRByRunoff_2DH=0._r8

  allocate(HeatFlow2Soili_3D(3,JD,JV,JH));  HeatFlow2Soili_3D=0._r8
  allocate(WaterFlow2Micpt_3D(3,JD,JV,JH));   WaterFlow2Micpt_3D=0._r8  
  allocate(WaterFlow2Macpt_3D(3,JD,JV,JH));  WaterFlow2Macpt_3D=0._r8
  allocate(WaterFlow2Micptl_3D(3,JD,JV,JH));   WaterFlow2Micptl_3D=0._r8  
  allocate(WaterFlow2Macptl_3D(3,JD,JV,JH));  WaterFlow2Macptl_3D=0._r8
  allocate(WaterFlow2MicptX_3D(3,JD,JV,JH));  WaterFlow2MicptX_3D=0._r8
  allocate(VLWatMicP2_vr(JZ,JY,JX));    VLWatMicP2_vr=0._r8
  allocate(VLairMicP_vr(JZ,JY,JX));   VLairMicP_vr=0._r8
  allocate(VLWatMicPX1_vr(JZ,JY,JX));   VLWatMicPX1_vr=0._r8

  allocate(LWEmscefLitR_col(JY,JX));       LWEmscefLitR_col=0._r8
  allocate(LWRad2LitR_col(JY,JX));       LWRad2LitR_col=0._r8
  allocate(LWEmscefSoil_col(JY,JX));       LWEmscefSoil_col=0._r8
  allocate(RadSW2LitR_col(JY,JX));       RadSW2LitR_col=0._r8    
  allocate(LWRad2Soil_col(JY,JX));       LWRad2Soil_col=0._r8
  allocate(RadSW2Soil_col(JY,JX));       RadSW2Soil_col=0._r8
  allocate(Qinflx2SoilM_col(JY,JX)); Qinflx2SoilM_col=0._r8
  end subroutine InitHydroThermData

!------------------------------------------------------------------------------------------
  subroutine DestructHydroThermData
  ! Release substep hydrological/thermal work arrays.
  use abortutils, only : destroy  
  implicit none

  call destroy(WatFlx2LitRByRunoff_2DH)
  call destroy(FracSoiPAsWat_vr)
  call destroy(PSISM1_vr)  
  call destroy(TKSoil1_vr)  
  call destroy(DLYRR_COL)  
  call destroy(FracSoiPAsIce_vr)  
  call destroy(FracAirFilledSoilPore_vr)  
  call destroy(VHeatCapacity1_vr)  
  call destroy(FracSoilAsAirt)  
  call destroy(VLWatMicP1_vr)  
  call destroy(VLiceMicP1_vr)
  call destroy(SnowFallt_col)
  call destroy(Rain2Snowt_col)
  call destroy(VLairMacP1_vr)    
  call destroy(EVAPW_col)  
  call destroy(EVAPS_col)  
  call destroy(VapXAir2Sno_col)  
  call destroy(SoilFracAsMacP1_vr)  
  call destroy(PrecHeat2Snowt_col)

  call destroy(AScaledCdHOverSnow_col)  
  call destroy(Ice2Snowt_col)

  call destroy(LWEmscefSnow_col)
  call destroy(RawH2OSoil2CanAir_col)
  call destroy(RawH2OSnow2CanAir_col)
  call destroy(SurfSensHeatConduct_col)
  call destroy(SurfLatHeatVapDeriv_col)
  call destroy(LWRad2Snow_col)
  call destroy(VLairMicP1_vr)
  call destroy(TKSnow0_snvr)  
  call destroy(VLWatMacP1_vr)
  call destroy(SoilFracAsMicP_vr)
  call destroy(VLHeatCapacityA_vr)
  call destroy(VLHeatCapacityB_vr)  
  call destroy(RawIsoTSoil2CanAir_col)
  call destroy(AScaledCdWOverSnow_col)
  call destroy(Altitude_col)
  call destroy(VLiceMacP1_vr)
  call destroy(RadSW2Sno_col)
  call destroy(HeatFlx2LitRByRunoff_2DH)    

  call destroy(HeatFlow2Soili_3D)  
  call destroy(WaterFlow2Micpt_3D)  
  call destroy(WaterFlow2Macpt_3D)
  call destroy(WaterFlow2Micptl_3D)  
  call destroy(WaterFlow2Macptl_3D)
  call destroy(WaterFlow2MicptX_3D)
  call destroy(VLWatMicP2_vr)
  call destroy(VLairMicP_vr)
  call destroy(VLWatMicPX1_vr)
  call destroy(LWEmscefLitR_col)
  call destroy(LWRad2LitR_col)
  call destroy(LWEmscefSoil_col)
  call destroy(RadSW2LitR_col)  
  call destroy(LWRad2Soil_col)
  call destroy(RadSW2Soil_col)
  call destroy(Qinflx2SoilM_col)

  end subroutine DestructHydroThermData
end module HydroThermData