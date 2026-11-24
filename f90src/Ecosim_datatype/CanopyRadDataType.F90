module CanopyRadDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use GridConsts
  implicit none
  public
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  real(r8),target,allocatable :: SineLeafAngle(:)                        !sine of leaf angle,	[-]
  real(r8),target,allocatable :: CosineLeafAngle(:)                      !cosine of leaf angle,	[-]
  real(r8),target,allocatable :: OMEGA2Leaf(:,:,:)                       !sine of indirect sky radiation on leaf surface, [-]
  real(r8),target,allocatable :: OMEGX(:,:,:)                            !sine of indirect sky radiation on leaf surface/sine of indirect sky radiation, [-]
  integer,target,allocatable :: iScatteringDiffus(:,:,:)                 !flag for calculating backscattering of radiation in canopy, [-]
  real(r8),target,allocatable :: RadDifPARAbsorption_zsec(:,:,:,:,:,:)   !shade incoming PAR, [umol m-2 s-1]
  real(r8),target,allocatable :: RadTotPARAbsorption_zsec(:,:,:,:,:,:)   !sunlit incoming PAR, [umol m-2 s-1]
  real(r8),target,allocatable :: LeafAngleClass_pft(:,:,:,:)             !fractionction of leaves in different angle classes, [-]
  real(r8),target,allocatable :: LeafAreaZsec_brch(:,:,:,:,:,:,:)        !leaf surface area, [m2 d-2]
  real(r8),target,allocatable :: LeafEffArea_zsec(:,:,:,:,:,:,:)         !leaf irradiated surface area, [m2 d-2]
  real(r8),target,allocatable :: StemAreaZsec_brch(:,:,:,:,:,:)          !stem surface area, [m2 d-2]
  real(r8),target,allocatable :: RadSW_Canopy_col(:,:)                   !canopy intercepted shortwave radiation, [MJ d-2]
  real(r8),target,allocatable :: RadSWCanopyLAbsroption_pft(:,:,:,:)     !profile of canopy absorbed shortwave radiation, [MJ d-2 h-1]     
  real(r8) :: TotSineSkyAngles_grd
  real(r8) :: dangle
  private :: InitAllocate

  contains

  subroutine  InitCanopyRad

  implicit none

  real(r8) :: da
  real(r8) :: aa
  integer :: N
! NumLeafInclinationClasses: number of leaf inclination groups, 90 deg into NumLeafInclinationClasses groups

  call InitAllocate

  dangle=PICON2h/real(NumLeafInclinationClasses,r8)         !the angle section width

  DO N = 1, NumLeafInclinationClasses
    aa                 = real(N-0.5,r8)*dangle !leaf angle start from flat (0) to vertical (pi/2)
    SineLeafAngle(N)   = sin(aa)
    CosineLeafAngle(N) = cos(aa)
  ENDDO
  TotSineSkyAngles_grd = 0._r8
  end subroutine InitCanopyRad

!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none

  !JLS: number of leaf sectors divide the whole horizontal circle
  integer :: ncols
  ncols = bounds%ncols

  allocate(SineLeafAngle(NumLeafInclinationClasses))
  allocate(CosineLeafAngle(NumLeafInclinationClasses))
  allocate(OMEGA2Leaf(NumOfSkyAzimuthSects,NumLeafInclinationClasses,NumOfLeafAzimuthSectors));OMEGA2Leaf=0._r8
  allocate(OMEGX(NumOfSkyAzimuthSects,NumLeafInclinationClasses,NumOfLeafAzimuthSectors));OMEGX=0._r8
  allocate(iScatteringDiffus(NumOfSkyAzimuthSects,NumLeafInclinationClasses,NumOfLeafAzimuthSectors))
  allocate(LeafAngleClass_pft(NumLeafInclinationClasses,JP,JY,JX));LeafAngleClass_pft=0._r8
  allocate(LeafAreaZsec_brch(NumLeafInclinationClasses,NumCanopyLayers,MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafAreaZsec_brch=0._r8
  allocate(LeafEffArea_zsec(NumLeafInclinationClasses,NumCanopyLayers,MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafEffArea_zsec=0._r8
  allocate(RadTotPARAbsorption_zsec(NumLeafInclinationClasses,NumOfSkyAzimuthSects,NumCanopyLayers,JP,JY,JX));RadTotPARAbsorption_zsec=0._r8
  allocate(RadDifPARAbsorption_zsec(NumLeafInclinationClasses,NumOfSkyAzimuthSects,NumCanopyLayers,JP,JY,JX));RadDifPARAbsorption_zsec=0._r8
  allocate(StemAreaZsec_brch(NumLeafInclinationClasses,NumCanopyLayers,MaxNumBranches,JP,JY,JX));StemAreaZsec_brch=0._r8
  allocate(RadSWCanopyLAbsroption_pft(NumCanopyLayers,JP,JY,JX)); RadSWCanopyLAbsroption_pft=0._r8
  allocate(RadSW_Canopy_col(JY,JX)); RadSW_Canopy_col=0._r8
  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructCanopyRad
  use abortutils, only : destroy
  implicit none
  call destroy(SineLeafAngle)
  call destroy(CosineLeafAngle)
  call destroy(OMEGA2Leaf)
  call destroy(OMEGX)
  call destroy(iScatteringDiffus)
  call destroy(LeafAngleClass_pft)
  call destroy(LeafAreaZsec_brch)
  call destroy(LeafEffArea_zsec)
  call destroy(RadSWCanopyLAbsroption_pft)
  call destroy(RadTotPARAbsorption_zsec)
  call destroy(RadDifPARAbsorption_zsec)
  call destroy(StemAreaZsec_brch)
  call destroy(RadSW_Canopy_col)
  end subroutine DestructCanopyRad
end module CanopyRadDataType
