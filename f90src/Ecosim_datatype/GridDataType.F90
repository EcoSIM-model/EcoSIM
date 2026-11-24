module GridDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8) :: TAREA                                                  !total area of landscape,	[m2]
  real(r8),target,allocatable ::  CumDepz2LayBottom_vr(:,:,:)        !depth to bottom of soil layer, [m]
  real(r8),target,allocatable ::  DLYR_3D(:,:,:,:)                   !thickness of soil layer, [m]
  real(r8),target,allocatable ::  DLYRI_3D(:,:,:,:)                  !thickness of soil layer in 3 directions, [m]
  real(r8),target,allocatable ::  SoilDepthMidLay_vr(:,:,:)          !depth to middle of soil layer, [m]
  real(r8),target,allocatable ::  CumSoilThickness_vr(:,:,:)         !depth to bottom of soil layer from  surface of grid cell, [m]
  real(r8),target,allocatable ::  CumSoilThickMidL_vr(:,:,:)         !depth to middle of soil layer from  surface of grid cell, [m]
  real(r8),target,allocatable ::  AREA_3D(:,:,:,:)                   !cross-sectional area,  [m2 d-2]
  integer,target,allocatable ::  NU_col(:,:)                         !soil surface layer number,[-]
  integer,target,allocatable ::  NUI_col(:,:)                        !initial soil surface layer number,[-]
  integer,target,allocatable ::  NLF_col(:,:)                        !id of the frozen layer, [-]
  integer,target,allocatable ::  MaxNumRootLays_col(:,:)             !maximum root layer number,[-]
  integer,target,allocatable ::  NK_col(:,:)                         !additional soil lower boundary layers,[-]
  integer,target,allocatable ::  NLI_col(:,:)                        !initial lowest soil layer number,[-]
  integer,target,allocatable ::  NL_col(:,:)                         !lowest soil layer number,[-]
  integer,target,allocatable ::  NUM_col(:,:)                        !new surface layer number,[-]
  real(r8),target,allocatable ::  CumLitRDepzInit_col(:,:)           !initial position of the bottom of liter layer, [m]
  REAL(R8),target,allocatable ::  ALAT_col(:,:)                      !latitude,	[degrees north]
  real(r8),target,allocatable ::  DH_col(:,:)                        !Eeast-West width of the grid cells, [m]
  real(r8),target,allocatable ::  DV_col(:,:)                        !North-South width of the grid cells, [m]
  integer,target,allocatable ::  FlowDirIndicator_col(:,:)           !dimension of low,[-]
!----------------------------------------------------------------------

contains

  subroutine InitGridData

  implicit none
  allocate(NLF_col(JY,JX)); NLF_col=0
  allocate(CumDepz2LayBottom_vr(0:JZ,JY,JX));  CumDepz2LayBottom_vr=0._r8
  allocate(DLYR_3D(3,0:JZ,JY,JX)); DLYR_3D=0._r8
  allocate(DLYRI_3D(3,0:JZ,JY,JX));DLYRI_3D=0._r8
  allocate(SoilDepthMidLay_vr(JZ,JY,JX));     SoilDepthMidLay_vr=0._r8
  allocate(CumSoilThickness_vr(0:JZ,JY,JX)); CumSoilThickness_vr=0._r8
  allocate(CumSoilThickMidL_vr(JZ,JY,JX));    CumSoilThickMidL_vr=0._r8
  allocate(AREA_3D(3,0:JZ,JY,JX)); AREA_3D=0._r8
  allocate(NU_col(JY,JX));          NU_col=0
  allocate(NUI_col(JY,JX));         NUI_col=0
  allocate(MaxNumRootLays_col(JY,JX));          MaxNumRootLays_col=0
  allocate(NK_col(JY,JX));          NK_col=0
  allocate(NLI_col(JV,JH));         NLI_col=0
  allocate(NL_col(JV,JH));          NL_col=0
  allocate(NUM_col(JY,JX));         NUM_col=0
  allocate(CumLitRDepzInit_col(JY,JX));      CumLitRDepzInit_col=0._r8
  allocate(ALAT_col(JY,JX));        ALAT_col=0._r8
  allocate(DH_col(JY,JX));          DH_col=0._r8
  allocate(DV_col(JY,JX));          DV_col=0._r8
  allocate(FlowDirIndicator_col(JY,JX));         FlowDirIndicator_col=3   !vertical by default

  end subroutine InitGridData

!----------------------------------------------------------------------
  subroutine DestructGridData
  use abortutils, only : destroy
  implicit none

  call destroy(NLF_col)
  call destroy(CumDepz2LayBottom_vr)
  call destroy(DLYR_3D)
  call destroy(DLYRI_3D)
  call destroy(SoilDepthMidLay_vr)
  call destroy(CumSoilThickness_vr)
  call destroy(CumSoilThickMidL_vr)
  call destroy(AREA_3D)
  call destroy(NU_col)
  call destroy(NUI_col)
  call destroy(MaxNumRootLays_col)
  call destroy(NK_col)
  call destroy(NLI_col)
  call destroy(NL_col)
  call destroy(NUM_col)
  call destroy(CumLitRDepzInit_col)
  call destroy(ALAT_col)
  call destroy(DH_col)
  call destroy(DV_col)
  call destroy(FlowDirIndicator_col)
  end subroutine DestructGridData

end module GridDataType
