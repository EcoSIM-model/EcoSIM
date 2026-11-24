module RedistDataMod

  use GridDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : destroy
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc
  use EcoSIMConfig, only : nlbiomcp=>NumLiveMicrbCompts,ndbiomcp=>NumDeadMicrbCompts
implicit none

  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public

  real(r8),allocatable ::  trcs_TransptMicP_vr(:,:,:,:)             !net tracer flux into the grid cell through micropore flow (<0 loss) [g d-2 h-1]
  real(r8),allocatable ::  trcs_TransptMacP_vr(:,:,:,:)             !net tracer flux into the grid cell through macropore flow (<0 loss) [g d-2 h-1] 
  real(r8),allocatable ::  trcSalt_Flo2MicP_vr(:,:,:,:)             !net salt flux into the grid cell through micropore flow (<0 loss) [mol d-2 h-1]
  real(r8),allocatable ::  trcSalt_Flo2MacP_vr(:,:,:,:)             !net salt flux into the grid cell through macropore flow (<0 loss) [mol d-2 h-1]
  real(r8),allocatable ::  trcSalt_SurfRunoff_flx(:,:,:)            !salt tracer loss through surface runoff [mol d-2 h-1]

  real(r8),allocatable ::  TSandErosed_col(:,:)                       !column sand loss through erosion (<0 loss) [g d-2 h-1]
  real(r8),allocatable ::  TSiltErosed_col(:,:)                       !column silt loss through erosion (<0 loss) [g d-2 h-1]
  real(r8),allocatable ::  TCLAYErosed_col(:,:)                       !column clay loss through erosion (<0 loss) [g d-2 h-1]
  real(r8),allocatable ::  TNH4Erosed_molN_col(:,:)                        !column NH4 loss through erosion (<0 loss) [mol d-2 h-1]
  real(r8),allocatable ::  TNH3Erosed_molN_col(:,:)                        !column NH3 loss through erosion (<0 loss) [mol d-2 h-1]
  real(r8),allocatable ::  TNUreaErosed_molN_col(:,:)                      !column Urea loss through erosion (<0 loss) [mol d-2 h-1]
  real(r8),allocatable ::  TNO3Erosed_molN_col(:,:)                        !column NO3 loss through erosion (<0 loss) [mol d-2 h-1]
  real(r8),allocatable ::  TNH4ErosBand_molN_col(:,:)                    !column band NH4 loss through erosion (<0 loss) [molN d-2 h-1]
  real(r8),allocatable ::  TNH3ErosBand_molN_col(:,:)                    !column band NH3 loss through erosion (<0 loss) [molN d-2 h-1]
  real(r8),allocatable ::  TNUreaErosBand_molN_col(:,:)                  !column band urea loss through erosion (<0 loss) [molN d-2 h-1]
  real(r8),allocatable ::  TNO3ErosBand_molN_col(:,:)                    !column band NO3 loss through erosion (<0 loss) [molN d-2 h-1]
  real(r8),allocatable ::  TPO4Erosed_molP_col(:,:)                      !column soil PO4 los through erosion (<0 loss) [molP d-2 h-1]
  real(r8),allocatable ::  TPO4ErosBand_molP_col(:,:)                    !column band PO4 loss through erosion (<0 loss) [molP d-2 h-1]
  real(r8),allocatable ::  trcx_TER_col(:,:,:)                      !column exchangeable cation loss through erosion (<0 loss)[mol d-2 h-1]
  real(r8),allocatable ::  trcp_TER_col(:,:,:)                      !column precipitated salt loss through erosion (<0 loss) [mol d-2 h-1]
  real(r8),allocatable ::  tErosionSedmLoss_col(:,:)                !column sediment loss through erosion (<0 loss) [g d-2 h-1]

  real(r8),allocatable ::  TWatFlowCellMicPX_vr(:,:,:)              !

  real(r8),allocatable ::  trcg_SurfRunoff_flx(:,:,:)               !
  real(r8),allocatable ::  Gas_AdvDif_Flx_vr(:,:,:,:)               !gas flux into the grid through advection and diffusion [g d-2 h-1]

  real(r8),allocatable ::  WatIceThawMicP_vr(:,:,:)                 !water added due to ice thaw to the grid mairopore (>0 thaw) [m3 H2O d-2 h-1]
  real(r8),allocatable ::  THeatSoiThaw_vr(:,:,:)                   !heat released due to water freezing in the grid (>0 release) [MJ d-2 h-1]
  real(r8),allocatable ::  WatIceThawMacP_vr(:,:,:)                 !water added due to ice thaw to the grid macropore (>0 thaw) [m3 H2O d-2 h-1]
  real(r8),allocatable ::  VLWatMicP2_vr(:,:,:)                     !local copy of soil water in micropore
  real(r8),allocatable ::  VLiceMicP2_vr(:,:,:)                     !local copy of soil ice in micropore
  real(r8),allocatable ::  VLWatMacP2_vr(:,:,:)                     !local copy of soil water in macropore
  real(r8),allocatable ::  VLiceMacP2_vr(:,:,:)                     !local copy of soil ice in macropore

  real(r8),allocatable :: TOMEERhetr_col(:,:,:,:,:)                 !hetetroph biomass loss due to erosion [g d-2 h-1]
  real(r8),allocatable :: TOMEERauto_col(:,:,:,:)                   !autotroph biomass loss due to erosion [g d-2 h-1]
 
  real(r8),allocatable ::  DOM_Transp2Micp_vr(:,:,:,:,:)            !DOM added to the grid cell due to micropore transport [g d-2 h-1]
  real(r8),allocatable ::  DOM_Transp2Macp_flx(:,:,:,:,:)           !DOM added to the grid cell due to macropore transport [g d-2 h-1]
  real(r8),allocatable ::  DOM_SurfRunoff_flx(:,:,:,:)                      !DOM added (>0) to the litter layer due to surface runoff [g d-2 h-1]
  real(r8),allocatable ::  TORMER_col(:,:,:,:,:)                    !Microbial residue loss due to erosion (<0 loss) [g d-2 h-1]
  real(r8),allocatable ::  TOHMER_col(:,:,:,:)                      !sorbed OM loss due to erosion (<0 loss) [g d-2 h-1]
  real(r8),allocatable ::  TOSMER_col(:,:,:,:,:)                    !total solid OM loss due to erosion (<0 loss) [g d-2 h-1]
  real(r8),allocatable ::  TOSAER_col(:,:,:,:)                      !active solid OM loss due to erosion (<0 loss) [g d-2 h-1]


  contains

!------------------------------------------------------------------------------------------

  subroutine InitTflxType()

  implicit none

  allocate(trcs_TransptMacP_vr(ids_beg:ids_end,JZ,JY,JX));   trcs_TransptMacP_vr=0._r8
  allocate(Gas_AdvDif_Flx_vr(idg_beg:idg_NH3,JZ,JY,JX));   Gas_AdvDif_Flx_vr=0._r8

  allocate(trcSalt_Flo2MicP_vr(idsalt_beg:idsaltb_end,JZ,JY,JX)); trcSalt_Flo2MicP_vr=0._r8
  allocate(trcSalt_Flo2MacP_vr(idsalt_beg:idsaltb_end,JZ,JY,JX)); trcSalt_Flo2MacP_vr=0._r8
  allocate(trcSalt_SurfRunoff_flx(idsalt_beg:idsalt_end,JY,JX));           trcSalt_SurfRunoff_flx=0._r8  
  allocate(TSandErosed_col(JY,JX));      TSandErosed_col=0._r8
  allocate(TSiltErosed_col(JY,JX));      TSiltErosed_col=0._r8
  allocate(TCLAYErosed_col(JY,JX));      TCLAYErosed_col=0._r8
  allocate(TPO4Erosed_molP_col(JY,JX));  TPO4Erosed_molP_col=0._r8
  allocate(TPO4ErosBand_molP_col(JY,JX)); TPO4ErosBand_molP_col=0._r8
  allocate(TNH4Erosed_molN_col(JY,JX));      TNH4Erosed_molN_col=0._r8
  allocate(TNH3Erosed_molN_col(JY,JX));      TNH3Erosed_molN_col=0._r8
  allocate(TNUreaErosed_molN_col(JY,JX));      TNUreaErosed_molN_col=0._r8
  allocate(TNO3Erosed_molN_col(JY,JX));      TNO3Erosed_molN_col=0._r8
  allocate(TNH4ErosBand_molN_col(JY,JX));      TNH4ErosBand_molN_col=0._r8
  allocate(TNH3ErosBand_molN_col(JY,JX));      TNH3ErosBand_molN_col=0._r8
  allocate(TNUreaErosBand_molN_col(JY,JX));      TNUreaErosBand_molN_col=0._r8
  allocate(TNO3ErosBand_molN_col(JY,JX));      TNO3ErosBand_molN_col=0._r8
  allocate(trcg_SurfRunoff_flx(idg_beg:idg_NH3,JY,JX));      trcg_SurfRunoff_flx=0._r8

  allocate(trcx_TER_col(idx_beg:idx_end,JY,JX));    trcx_TER_col=0._r8
  allocate(trcp_TER_col(idsp_beg:idsp_end,JY,JX));      trcp_TER_col=0._r8
  allocate(trcs_TransptMicP_vr(ids_beg:ids_end,JZ,JY,JX));   trcs_TransptMicP_vr=0._r8

  allocate(tErosionSedmLoss_col(JY,JX));      tErosionSedmLoss_col=0._r8

  allocate(TWatFlowCellMicPX_vr(JZ,JY,JX));    TWatFlowCellMicPX_vr=0._r8

  allocate(WatIceThawMicP_vr(JZ,JY,JX));    WatIceThawMicP_vr=0._r8
  allocate(THeatSoiThaw_vr(JZ,JY,JX));   THeatSoiThaw_vr=0._r8
  allocate(WatIceThawMacP_vr(JZ,JY,JX));   WatIceThawMacP_vr=0._r8

!  allocate(VLWatMicP2_vr(JZ,JY,JX));    VLWatMicP2_vr=0._r8
!  allocate(VLiceMicP2_vr(JZ,JY,JX));    VLiceMicP2_vr=0._r8
!  allocate(VLWatMacP2_vr(JZ,JY,JX));   VLWatMacP2_vr=0._r8
!  allocate(VLiceMacP2_vr(JZ,JY,JX));   VLiceMacP2_vr=0._r8
  allocate(TOMEERhetr_col(NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,JY,JX)); TOMEERhetr_col=0._r8

  allocate(TOMEERauto_col(NumPlantChemElms,1:NumLiveAutoBioms,JY,JX));TOMEERauto_col=0._r8

  allocate(DOM_Transp2Micp_vr(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Micp_vr=0._r8
  allocate(DOM_Transp2Macp_flx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Macp_flx=0._r8
  allocate(DOM_SurfRunoff_flx(idom_beg:idom_end,1:jcplx,JY,JX));DOM_SurfRunoff_flx=0._r8
  allocate(TORMER_col(NumPlantChemElms,ndbiomcp,1:jcplx,JY,JX));TORMER_col=0._r8
  allocate(TOHMER_col(idom_beg:idom_end,1:jcplx,JY,JX));TOHMER_col=0._r8
  allocate(TOSMER_col(NumPlantChemElms,jsken,1:jcplx,JY,JX));TOSMER_col=0._r8
  allocate(TOSAER_col(jsken,1:jcplx,JY,JX));TOSAER_col=0._r8

  end subroutine InitTflxType


!------------------------------------------------------------------------------------------

  subroutine DestructTflxType

  implicit none

  call destroy(trcx_TER_col)
  call destroy(TOMEERhetr_col)
  call destroy(TOMEERauto_col)
  call destroy(DOM_SurfRunoff_flx)
  call destroy(TORMER_col)
  call destroy(TOHMER_col)
  call destroy(TOSMER_col)
  call destroy(TOSAER_col)
  call destroy(trcSalt_Flo2MicP_vr)
  call destroy(trcs_TransptMicP_vr)
  call destroy(trcSalt_Flo2MacP_vr)
!  call destroy(TFLWS)
!  call destroy(TFLWW)
!  call destroy(TFLWI)
!  call destroy(THFLWW)

  call destroy(trcs_TransptMacP_vr)
  call destroy(TPO4ErosBand_molP_col)
  call destroy(TPO4Erosed_molP_col)
  call destroy(TSandErosed_col)
  call destroy(TSiltErosed_col)
  call destroy(TCLAYErosed_col)
  call destroy(TNH4Erosed_molN_col)
  call destroy(TNH3Erosed_molN_col)
  call destroy(TNUreaErosed_molN_col)
  call destroy(TNO3Erosed_molN_col)
  call destroy(TNH4ErosBand_molN_col)
  call destroy(TNH3ErosBand_molN_col)
  call destroy(TNUreaErosBand_molN_col)
  call destroy(TNO3ErosBand_molN_col)
  call destroy(tErosionSedmLoss_col)
  call destroy(TWatFlowCellMicPX_vr)
  call destroy(WatIceThawMicP_vr)
  call destroy(THeatSoiThaw_vr)
  call destroy(WatIceThawMacP_vr)

!  call destroy(VLWatMicP2_vr)
!  call destroy(VLiceMicP2_vr)
!  call destroy(VLWatMacP2_vr)
!  call destroy(VLiceMacP2_vr)

  call destroy(DOM_Transp2Micp_vr)
  call destroy(DOM_Transp2Macp_flx)
  call destroy(trcp_TER_col)
  call destroy(Gas_AdvDif_Flx_vr)
  call destroy(trcg_SurfRunoff_flx)
  call destroy(trcSalt_SurfRunoff_flx)

  end subroutine DestructTflxType
end module RedistDataMod
