module WatsubDataMod
! Description:
! Own temporary soil hydrology and heat workspaces shared by surface staging
! and the native soil solver. Allocation remains with InitWatsubData.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),allocatable ::  TMLiceThawedMicP_vr(:,:,:)              !micropore layer integrated ice mass loss due to thaw,
  real(r8),allocatable ::  TMLiceThawedMacP_vr(:,:,:)              !macropore layer integrated ice mass loss due to thaw,

  real(r8),allocatable ::  FracLayVolBelowExtWTBL_vr(:,:,:)          !fraction of layer volume below external water table [0-1]
  real(r8),allocatable ::  FracLayVolBelowTileWTBL_vr(:,:,:)                !
  real(r8),allocatable ::  VLairMacP_vr(:,:,:)                     !
  real(r8),allocatable ::  TLPhaseChangeHeat2Soi1_vr(:,:,:)         !total soil layer latent heat release from melting
  real(r8),allocatable ::  QLaterFlow2CellM_col(:,:)               !lateral water flow from neighbors into grid
  real(r8),allocatable ::  FWatExMacP2MicPiM_vr(:,:,:)             !pressure-driven water flow from macpore to micpore
  real(r8),allocatable ::  dWaterPlantRoot2SoilPrev_vr(:,:,:)      !water flux taken/release by plant roots [m3 H2O d-2]
  real(r8),allocatable ::  dHeatPlantRoot2SoilPrev_vr(:,:,:)       !heat flux taken/released by plant roots [MJ d-2]
  real(r8),allocatable ::  WatNetFlow2MicP_3DM_vr(:,:,:)           !total water charge to micropore in iteration M
  real(r8),allocatable ::  WatNetFlow2Macpt_3DM_vr(:,:,:)          !total water charge to macropore in iteration M
  real(r8),allocatable ::  THeatFlow2Soil_3DM_vr(:,:,:)            !total heat flux add to soil layer in iteration M
  real(r8),allocatable ::  FIceThawedMicP_vr(:,:,:)                !Ice mass thawed in micropore during an iteration time step [ton H2O/d2]
  real(r8),allocatable ::  SoilWatFrezHeatRelease_vr(:,:,:)        !latent energy released during freeze-thaw
  real(r8),allocatable ::  AVCNHL_3D(:,:,:,:)                      !Macropore hydraulic conductivity
  real(r8),allocatable :: TWaterPlantRoot2SoilX_vr(:,:,:)          !water flux from plant root to soil at iteration M [m3 H2O d-2]
  real(r8),allocatable ::  WatNetFlow2MicptX_3DM_vr(:,:,:)         !water flux in iteration M [m3 H2O/d2]
  real(r8),allocatable ::  FWatIrrigate2MicP1_vr(:,:,:)            !water flux due to irrigation [m3 H2O/d2]
  real(r8),allocatable ::  HeatIrrigation1_vr(:,:,:)               !heat flux due to irrigation, [MJ/d2]

  real(r8),allocatable :: FIceThawedMacP_vr(:,:,:)                 !Ice mass thawed in macropre during an iteration time step, [ton H2O/d2]
  real(r8),allocatable :: TWaterPlantRoot2SoilX_col(:,:)            !cumulative water flux from roots to soil [m3 H2O d-2]
  real(r8),allocatable :: TWaterPlantRoot2SoilXM_col(:,:)
  real(r8),allocatable ::  HydroCondMacP1_vr(:,:,:)                 !
  real(r8),allocatable ::  VLMicP1_vr(:,:,:)                        !
  real(r8),allocatable ::  VLMacroPores1_vr(:,:,:)                        !

  real(r8),allocatable ::  WatFlow2TopSoiMicP_col(:,:)              !
  real(r8),allocatable ::  WatFlow2TopSoiMicPX_col(:,:)             !
  real(r8),allocatable ::  WatFlow2TopSoiMacP_col(:,:)              !
  real(r8),allocatable ::  HeatFlow2TopSoi_col(:,:)                 !
  real(r8),allocatable ::  PSISoilMatricPtmp_vr(:,:,:)              !
  real(r8),allocatable :: QDischarM_col(:,:)
  real(r8),allocatable :: QDrainM_col(:,:)

  integer, allocatable ::  N6X(:,:)

!----------------------------------------------------------------------
  public :: InitWatsubData
contains
  subroutine InitWatSubData

  implicit none
  allocate(TWaterPlantRoot2SoilXM_col(JY,JX)); TWaterPlantRoot2SoilXM_col=0._r8
  allocate(TWaterPlantRoot2SoilX_col(JY,JX)); TWaterPlantRoot2SoilX_col=0._r8
  allocate(TWaterPlantRoot2SoilX_vr(JZ,JY,JX)); TWaterPlantRoot2SoilX_vr=0._r8
  allocate(QLaterFlow2CellM_col(JY,JX)); QLaterFlow2CellM_col=0._r8
  allocate(N6X(JY,JX));         N6X=0
  allocate(TMLiceThawedMicP_vr(JZ,JY,JX));   TMLiceThawedMicP_vr=0._r8
  allocate(TMLiceThawedMacP_vr(JZ,JY,JX));   TMLiceThawedMacP_vr=0._r8
  allocate(FracLayVolBelowExtWTBL_vr(JZ,JY,JX));    FracLayVolBelowExtWTBL_vr=0._r8
  allocate(FracLayVolBelowTileWTBL_vr(JZ,JY,JX));   FracLayVolBelowTileWTBL_vr=0._r8
  allocate(QDischarM_col(JY,JX)); QDischarM_col(:,:)=0._r8
  allocate(QDrainM_col(JY,JX)); QDrainM_col(:,:)=0._r8
  allocate(VLairMacP_vr(JZ,JY,JX));  VLairMacP_vr=0._r8
  allocate(TLPhaseChangeHeat2Soi1_vr(JZ,JY,JX));   TLPhaseChangeHeat2Soi1_vr=0._r8
  allocate(FWatExMacP2MicPiM_vr(JZ,JY,JX));    FWatExMacP2MicPiM_vr=0._r8

  allocate(WatNetFlow2MicP_3DM_vr(JZ,JY,JX));    WatNetFlow2MicP_3DM_vr=0._r8
  allocate(WatNetFlow2Macpt_3DM_vr(JZ,JY,JX));   WatNetFlow2Macpt_3DM_vr=0._r8
  allocate(THeatFlow2Soil_3DM_vr(JZ,JY,JX));   THeatFlow2Soil_3DM_vr=0._r8
  allocate(FIceThawedMicP_vr(JZ,JY,JX));    FIceThawedMicP_vr=0._r8
  allocate(SoilWatFrezHeatRelease_vr(JZ,JY,JX));    SoilWatFrezHeatRelease_vr=0._r8
  allocate(AVCNHL_3D(3,JD,JV,JH)); AVCNHL_3D=0._r8

  allocate(WatNetFlow2MicptX_3DM_vr(JZ,JY,JX));   WatNetFlow2MicptX_3DM_vr=0._r8
  allocate(FWatIrrigate2MicP1_vr(JZ,JY,JX));     FWatIrrigate2MicP1_vr=0._r8
  allocate(HeatIrrigation1_vr(JZ,JY,JX));   HeatIrrigation1_vr=0._r8

  allocate(FIceThawedMacP_vr(JZ,JY,JX));   FIceThawedMacP_vr=0._r8
  allocate(dWaterPlantRoot2SoilPrev_vr(JZ,JY,JX)); dWaterPlantRoot2SoilPrev_vr=0._r8
  allocate(dHeatPlantRoot2SoilPrev_vr(JZ,JY,JX)); dHeatPlantRoot2SoilPrev_vr=0._r8
  allocate(HydroCondMacP1_vr(JZ,JY,JX));    HydroCondMacP1_vr=0._r8
  allocate(VLMicP1_vr(0:JZ,JY,JX));  VLMicP1_vr=0._r8
  allocate(VLMacroPores1_vr(JZ,JY,JX));   VLMacroPores1_vr=0._r8

  allocate(WatFlow2TopSoiMicP_col(JY,JX));       WatFlow2TopSoiMicP_col=0._r8
  allocate(WatFlow2TopSoiMicPX_col(JY,JX));      WatFlow2TopSoiMicPX_col=0._r8
  allocate(WatFlow2TopSoiMacP_col(JY,JX));      WatFlow2TopSoiMacP_col=0._r8
  allocate(HeatFlow2TopSoi_col(JY,JX));      HeatFlow2TopSoi_col=0._r8
  allocate(PSISoilMatricPtmp_vr(JZ,JY,JX));   PSISoilMatricPtmp_vr=0._r8

  end subroutine InitWatSubData

!----------------------------------------------------------------------
  subroutine DestructWatSubData
  use abortutils, only : destroy
  implicit none

  call destroy(TWaterPlantRoot2SoilXM_col)
  call destroy(TWaterPlantRoot2SoilX_col)
  call destroy(TWaterPlantRoot2SoilX_vr)
  call destroy(QLaterFlow2CellM_col)
  call destroy(N6X)
  call destroy(dWaterPlantRoot2SoilPrev_vr)
  call destroy(dHeatPlantRoot2SoilPrev_vr)
  call destroy(TMLiceThawedMicP_vr)
  call destroy(TMLiceThawedMacP_vr)
  call destroy(FracLayVolBelowExtWTBL_vr)
  call destroy(FracLayVolBelowTileWTBL_vr)
  call destroy(VLairMacP_vr)
  call destroy(TLPhaseChangeHeat2Soi1_vr)
  call destroy(FWatExMacP2MicPiM_vr)
  call destroy(WatNetFlow2MicP_3DM_vr)
  call destroy(WatNetFlow2Macpt_3DM_vr)
  call destroy(THeatFlow2Soil_3DM_vr)
  call destroy(FIceThawedMicP_vr)
  call destroy(SoilWatFrezHeatRelease_vr)
  call destroy(AVCNHL_3D)
  call destroy(WatNetFlow2MicptX_3DM_vr)
  call destroy(FWatIrrigate2MicP1_vr)
  call destroy(HeatIrrigation1_vr)
  call destroy(FIceThawedMacP_vr)
  call destroy(HydroCondMacP1_vr)
  call destroy(VLMicP1_vr)
  call destroy(VLMacroPores1_vr)
  call destroy(WatFlow2TopSoiMicP_col)
  call destroy(WatFlow2TopSoiMicPX_col)
  call destroy(WatFlow2TopSoiMacP_col)
  call destroy(HeatFlow2TopSoi_col)
  call destroy(PSISoilMatricPtmp_vr)
  call destroy(QDischarM_col)
  call destroy(QDrainM_col)

  end subroutine DestructWatSubData

end module WatsubDataMod
