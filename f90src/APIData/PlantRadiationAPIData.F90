module PlantRadiationAPIData
  ! Description:
  ! Owns the radiation API type, instance and lifecycle procedures.
  use DebugToolMod, only: PrintInfo
  use PlantAPICommonData
  implicit none
  private
  public :: plt_rad_init, plt_rad_destroy

  type, public ::  plant_radiation_type
  real(r8) :: TotSineSkyAngles_grd           !sine of sky angles,[-]
  real(r8) :: SoilAlbedo                     !soil albedo,[-]
  real(r8) :: SurfAlbedo_col                 !Surface albedo,[-]
  real(r8) :: RadPARSolarBeam_col            !PAR radiation in solar beam, [umol m-2 s-1]
  real(r8) :: RadSWDiffus_col                !diffuse shortwave radiation, [W m-2]
  real(r8) :: RadSWDirect_col                !direct shortwave radiation, [W m-2]
  real(r8) :: RadPARDiffus_col               !diffuse PAR, [umol m-2 s-1]
  real(r8) :: RadDirectPAR_col               !direct PAR, [umol m-2 s-1]
  real(r8) :: Eco_NetRad_col                 !ecosystem net radiation, [MJ d-2 h-1]
  real(r8) :: RadSWSolarBeam_col             !shortwave radiation in solar beam, [MJ m-2 h-1]
  real(r8) :: FracSWRad2Grnd_col             !fraction of radiation intercepted by ground surface, [-]
  real(r8) :: RadSWGrnd_col                  !radiation intercepted by ground surface, [MJ m-2 h-1]
  real(r8) :: RadPARGrnd_col                 !PAR radiation reaching the ground, [umol m-2 s-1]
  real(r8) :: SineGrndSlope_col              !sine of slope, [-]
  real(r8) :: GroundSurfaceAzimuth_col          !azimuth of slope, [-]
  real(r8) :: CosineGrndSlope_col            !cosine of slope, [-]
  real(r8) :: LWRadGrnd_col                      !longwave radiation emitted by ground surface, [MJ m-2 h-1]
  real(r8) :: LWRadSky_col                   !sky longwave radiation , [MJ d-2 h-1]
  real(r8) :: SineSunInclAnglNxtHour_col     !sine of solar angle next hour, [-]
  real(r8) :: SineSunInclinationAngle_col           !sine of solar angle, [-]

  integer,  pointer :: iScatteringDiffus(:,:,:)  => null() !flag for calculating backscattering of radiation in canopy,[-]
  real(r8), pointer :: RadSWLeafAlbedo_pft(:)    => null() !canopy shortwave albedo,                           [-]
  real(r8), pointer :: RadPARLeafAlbedo_pft(:)    => null() !canopy PAR albedo,                                 [-]
  real(r8), pointer :: TAU_DirectSunLit(:)    => null() !fraction of radiation intercepted by canopy layer, [-]
  real(r8), pointer :: TAU_DirectSunSha(:)            => null() !fraction of radiation transmitted by canopy layer, [-]
  real(r8), pointer :: LWRadCanopy_pft(:)        => null() !canopy longwave radiation,                         [MJ d-2 h-1]
  real(r8), pointer :: RadSWCanopyAbsorption_pft(:)      => null() !canopy absorbed shortwave radiation,               [MJ d-2 h-1]
  real(r8), pointer :: RadSWCanopyLAbsroption_pft(:,:)   => null() !profile of canopy absorbed shortwave radiation, [MJ d-2 h-1]
  real(r8), pointer :: RadPARCanopyLAbsorption_pft(:,:)  => null() !profile of canopy absorbed PAR, [MJ d-2 h-1]
  real(r8), pointer :: OMEGX(:,:,:)              => null() !sine of indirect sky radiation on leaf surface/sine of indirect sky radiation,[-]
  real(r8), pointer :: OMEGA2Ground(:)           => null() !sine of solar beam on ground surface,                [-]
  real(r8), pointer :: OMEGA2Leaf(:,:,:)              => null() !sine of indirect sky radiation on leaf surface,[-]
  real(r8), pointer :: SineLeafAngle(:)          => null() !sine of leaf angle,[-]
  real(r8), pointer :: CosineLeafAngle(:)        => null() !cosine of leaf angle,[-]
  real(r8), pointer :: RadNet2Canopy_pft(:)      => null() !canopy net radiation,                              [MJ d-2 h-1]
  real(r8), pointer :: LeafSWabsorptivity_pft(:)     => null() !canopy shortwave absorptivity,                     [-]
  real(r8), pointer :: LeafPARabsorptivity_pft(:)    => null() !canopy PAR absorptivity,[-]
  real(r8), pointer :: RadSWLeafTransmitance_pft(:)  => null() !canopy shortwave transmissivity,                   [-]
  real(r8), pointer :: RadPARLeafTransmitance_pft(:) => null() !canopy PAR transmissivity,                         [-]
  real(r8), pointer :: RadPARCanopyAbsorption_pft(:)     => null() !canopy absorbed PAR,                               [umol m-2 s-1]
  real(r8), pointer :: FracPARads2Canopy_pft(:)     => null() !fraction of incoming PAR absorbed by total canopy, [-]
  real(r8), pointer :: FracPARads2LiveCanopy_pft(:) => null() !fraction of incoming PAR absorbed by live canopy,  [-]
  real(r8), pointer :: RadTotPARAbsorption_zsec(:,:,:,:)      => null()     !direct incoming PAR,                           [umol m-2 s-1]
  real(r8), pointer :: RadDifPARAbsorption_zsec(:,:,:,:)   => null()  !diffuse incoming PAR,                             [umol m-2 s-1]
  contains
    procedure, public :: Init    => plt_rad_init
    procedure, public :: Destroy => plt_rad_destroy
  end type plant_radiation_type


  type(plant_radiation_type), public, target :: plt_rad       !plant radiation type

contains

  subroutine plt_rad_init(this)
! DESCRIPTION
! initialize data type for plant_radiation_type
  implicit none
  class(plant_radiation_type) :: this

  allocate(this%RadTotPARAbsorption_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,NumCanopyLayers1,JP1));this%RadTotPARAbsorption_zsec=0._r8
  allocate(this%RadDifPARAbsorption_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,NumCanopyLayers1,JP1));this%RadDifPARAbsorption_zsec=0._r8
  allocate(this%RadSWLeafAlbedo_pft(JP1))
  allocate(this%RadPARLeafAlbedo_pft(JP1))
  allocate(this%TAU_DirectSunLit(NumCanopyLayers1+1));this%TAU_DirectSunLit=0._r8
  allocate(this%TAU_DirectSunSha(NumCanopyLayers1+1));this%TAU_DirectSunSha=0._r8
  allocate(this%LWRadCanopy_pft(JP1))
  allocate(this%RadSWCanopyLAbsroption_pft(NumCanopyLayers1,JP1)); this%RadSWCanopyLAbsroption_pft=0._r8
  allocate(this%RadPARCanopyLAbsorption_pft(NumCanopyLayers1,JP1)); this%RadPARCanopyLAbsorption_pft=0._r8
  allocate(this%RadSWCanopyAbsorption_pft(JP1))
  allocate(this%OMEGX(NumOfSkyAzimuthSects1,NumLeafInclinationClasses1,NumOfLeafAzimuthSectors1));this%OMEGX=0._r8
  allocate(this%OMEGA2Ground(NumOfSkyAzimuthSects1));this%OMEGA2Ground=0._r8
  allocate(this%OMEGA2Leaf(NumOfSkyAzimuthSects1,NumLeafInclinationClasses1,NumOfLeafAzimuthSectors1));this%OMEGA2Leaf=0._r8
  allocate(this%SineLeafAngle(NumLeafInclinationClasses1));this%SineLeafAngle=0._r8
  allocate(this%CosineLeafAngle(NumLeafInclinationClasses1));this%CosineLeafAngle=0._r8
  allocate(this%iScatteringDiffus(NumOfSkyAzimuthSects1,NumLeafInclinationClasses1,NumOfLeafAzimuthSectors1))
  allocate(this%RadNet2Canopy_pft(JP1))
  allocate(this%LeafSWabsorptivity_pft(JP1))
  allocate(this%LeafPARabsorptivity_pft(JP1))
  allocate(this%RadPARLeafTransmitance_pft(JP1))
  allocate(this%RadSWLeafTransmitance_pft(JP1))
  allocate(this%RadPARCanopyAbsorption_pft(JP1))
  allocate(this%FracPARads2Canopy_pft(JP1));     this%FracPARads2Canopy_pft=0._r8
  allocate(this%FracPARads2LiveCanopy_pft(JP1)); this%FracPARads2LiveCanopy_pft=0._r8
  end subroutine plt_rad_init

  subroutine plt_rad_destroy(this)
  ! Description:
  ! Release every array owned by this PlantRadiationAPIData workspace.
  class(plant_radiation_type), intent(inout) :: this ! owned API state; component units and dimensions as declared above
  character(len=*), parameter :: subname='plt_rad_destroy'

  call PrintInfo('beg '//subname)
  if(associated(this%iScatteringDiffus))deallocate(this%iScatteringDiffus);nullify(this%iScatteringDiffus)
  if(associated(this%RadSWLeafAlbedo_pft))deallocate(this%RadSWLeafAlbedo_pft);nullify(this%RadSWLeafAlbedo_pft)
  if(associated(this%RadPARLeafAlbedo_pft))deallocate(this%RadPARLeafAlbedo_pft);nullify(this%RadPARLeafAlbedo_pft)
  if(associated(this%TAU_DirectSunLit))deallocate(this%TAU_DirectSunLit);nullify(this%TAU_DirectSunLit)
  if(associated(this%TAU_DirectSunSha))deallocate(this%TAU_DirectSunSha);nullify(this%TAU_DirectSunSha)
  if(associated(this%LWRadCanopy_pft))deallocate(this%LWRadCanopy_pft);nullify(this%LWRadCanopy_pft)
  if(associated(this%RadSWCanopyAbsorption_pft))deallocate(this%RadSWCanopyAbsorption_pft);nullify(this%RadSWCanopyAbsorption_pft)
  if(associated(this%RadSWCanopyLAbsroption_pft))deallocate(this%RadSWCanopyLAbsroption_pft);nullify(this%RadSWCanopyLAbsroption_pft)
  if(associated(this%RadPARCanopyLAbsorption_pft))deallocate(this%RadPARCanopyLAbsorption_pft);nullify(this%RadPARCanopyLAbsorption_pft)
  if(associated(this%OMEGX))deallocate(this%OMEGX);nullify(this%OMEGX)
  if(associated(this%OMEGA2Ground))deallocate(this%OMEGA2Ground);nullify(this%OMEGA2Ground)
  if(associated(this%OMEGA2Leaf))deallocate(this%OMEGA2Leaf);nullify(this%OMEGA2Leaf)
  if(associated(this%SineLeafAngle))deallocate(this%SineLeafAngle);nullify(this%SineLeafAngle)
  if(associated(this%CosineLeafAngle))deallocate(this%CosineLeafAngle);nullify(this%CosineLeafAngle)
  if(associated(this%RadNet2Canopy_pft))deallocate(this%RadNet2Canopy_pft);nullify(this%RadNet2Canopy_pft)
  if(associated(this%LeafSWabsorptivity_pft))deallocate(this%LeafSWabsorptivity_pft);nullify(this%LeafSWabsorptivity_pft)
  if(associated(this%LeafPARabsorptivity_pft))deallocate(this%LeafPARabsorptivity_pft);nullify(this%LeafPARabsorptivity_pft)
  if(associated(this%RadSWLeafTransmitance_pft))deallocate(this%RadSWLeafTransmitance_pft);nullify(this%RadSWLeafTransmitance_pft)
  if(associated(this%RadPARLeafTransmitance_pft))deallocate(this%RadPARLeafTransmitance_pft);nullify(this%RadPARLeafTransmitance_pft)
  if(associated(this%RadPARCanopyAbsorption_pft))deallocate(this%RadPARCanopyAbsorption_pft);nullify(this%RadPARCanopyAbsorption_pft)
  if(associated(this%FracPARads2Canopy_pft))deallocate(this%FracPARads2Canopy_pft);nullify(this%FracPARads2Canopy_pft)
  if(associated(this%FracPARads2LiveCanopy_pft))deallocate(this%FracPARads2LiveCanopy_pft);nullify(this%FracPARads2LiveCanopy_pft)
  if(associated(this%RadTotPARAbsorption_zsec))deallocate(this%RadTotPARAbsorption_zsec);nullify(this%RadTotPARAbsorption_zsec)
  if(associated(this%RadDifPARAbsorption_zsec))deallocate(this%RadDifPARAbsorption_zsec);nullify(this%RadDifPARAbsorption_zsec)
  call PrintInfo('end '//subname)
  end subroutine plt_rad_destroy
end module PlantRadiationAPIData
