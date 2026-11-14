module PlantPhenologyAPIData
  ! Description:
  ! Own plant phenology API state and trait identifiers used by plant and canopy processes.
  ! Owns the phenology API type, instance and lifecycle procedures.
  use DebugToolMod, only: PrintInfo
  use PlantAPICommonData
  implicit none
  private
  public :: plt_pheno_init, plt_pheno_destroy

  type, public :: plant_pheno_type
  real(r8), pointer :: fTgrowRootP_vr(:,:)                => null()     !root layer temperature growth functiom,                              [-]
  real(r8), pointer :: ShootRootNonstElmConduts_pft(:)   => null()     !shoot-root rate constant for nonstructural C exchange,               [h-1]
  real(r8), pointer :: GrainFillRate25C_pft(:)            => null()     !maximum rate of fill per grain,                                      [g h-1]
  real(r8), pointer :: TempOffset_pft(:)                  => null()     !adjustment of Arhhenius curves for plant thermal acclimation,        [oC]
  real(r8), pointer :: CanPhenoMoistStress_pft(:)         => null()     !moisture stress for plant phenology development,[-]
  real(r8), pointer :: CanPhenoTempStress_pft(:)          => null()     !temperature stress for plant phenology development,[-]
  real(r8), pointer :: PlantO2Stress_pft(:)               => null()     !plant O2 stress indicator,                                           [-]
  real(r8), pointer :: NonstCMinConc2InitBranch_pft(:)    => null()     !branch nonstructural C content required for new branch,              [gC gC-1]
  real(r8), pointer :: NonstCMinCon2InitRoot_pft(:)       => null()     !threshold root nonstructural C content for initiating new root axis, [gC gC-1]
  real(r8), pointer :: LeafElmntRemobFlx_brch(:,:,:)      => null()    !cached remobilizable leaf element mass, [g d-2]
  real(r8), pointer :: LeafSenescInitialElms_brch(:,:,:) => null() !initial senescing leaf C/N/P mass; -1 means unset, [g d-2]
  real(r8), pointer :: PetolShethChemElmRemobFlx_brch(:,:,:) => null()    !cached remobilizable sheath element mass, [g d-2]
  real(r8), pointer :: PetolSenescInitialElms_brch(:,:,:) => null() !initial senescing sheath/petiole C/N/P mass; -1 means unset, [g d-2]
  real(r8), pointer :: TC4LeafOut_pft(:)                  => null()     !threshold temperature for spring leafout/dehardening,                [oC]
  real(r8), pointer :: TCGroth_pft(:)                     => null()     !canopy growth temperature,                                           [oC]
  real(r8), pointer :: TC4LeafOff_pft(:)                  => null()     !threshold temperature for autumn leafoff/hardening,                  [oC]
  real(r8), pointer :: TKGroth_pft(:)                     => null()     !canopy growth temperature,                                           [K]
  real(r8), pointer :: fTCanopyGroth_pft(:)               => null()     !canopy temperature growth function,                                  [-]
  real(r8), pointer :: MorphogenBase_pft(:)               => null()     !!baseline morphogen signal strength, [-]
  real(r8), pointer :: HoursTooLowPsiCan_pft(:)           => null()     !canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY), [h]
  real(r8), pointer :: TCChill4Seed_pft(:)           => null()     !temperature below which seed set is adversely affected, [oC]
  real(r8), pointer :: rPlantThermoAdaptZone_pft(:)  => null()     !plant thermal adaptation zone,                          [-]
  real(r8), pointer :: PlantInitThermoAdaptZone_pft(:)   => null()     !initial plant thermal adaptation zone,                  [-]
  real(r8), pointer :: HighTempLimitSeed_pft(:)      => null()     !temperature above which seed set is adversely affected, [oC]
  real(r8), pointer :: SeedTempSens_pft(:)           => null()     !sensitivity to HTC (seeds oC-1 above HTC),[oC-1]
  real(r8), pointer :: NetCumElmntFlx2Plant_pft(:,:) => null()     !effect of canopy element status on seed set,            [-]
  real(r8), pointer :: RefNodeInitRate_pft(:)        => null()     !rate of node initiation,                                [h-1 at 25 oC]
  real(r8), pointer :: RateRefLeafAppearance_pft(:)      => null()     !rate of leaf initiation,                                [h-1 at 25 oC]
  real(r8), pointer :: CriticPhotoPeriod_pft(:)      => null()     !critical daylength for phenological progress,           [h]
  real(r8), pointer :: PhotoPeriodSens_pft(:)        => null()     !difference between current and critical daylengths used to calculate  phenological progress, [h]
  integer,  pointer :: iEmbryophyteType_pft(:)       => null()      !plant embrophyte
  integer,  pointer :: iPlantStateLive_pft(:)            => null()  !flag for species death, [-]
  integer,  pointer :: FireReSet_pft(:)              => null()      !flag to skill startq for fire rejuvenation, [-]
  integer,  pointer :: iMaintPlantTrait_pft(:)          => null()   !flag for maintaining or reset plant traits, [-]
  integer,  pointer :: IsPlantActive_pft(:)              => null()  !flag for living pft, [-]
  integer,  pointer :: isPlantBranchAlive_brch(:,:)       => null()  !flag to detect branch death,                      [-]
  integer,  pointer :: doRemobilization_brch(:,:)        => null()  !branch phenology flag,                            [-]
  integer,  pointer :: doPlantLeaveOff_brch(:,:)         => null()  !branch phenology flag,                            [-]
  integer,  pointer :: EnablePlantLeafOut_brch(:,:)          => null()  !branch phenology flag,                            [-]
  integer,  pointer :: doInitLeafOut_brch(:,:)           => null()  !branch phenology flag,                            [-]
  logical,  pointer :: doReSeed_pft(:)                   => null()  !flag to do annual plant reseeding, [-]
  integer,  pointer :: doSenescence_brch(:,:)            => null()  !branch phenology flag,                            [-]
  integer,  pointer :: Prep4Literfall_brch(:,:)          => null()  !branch phenology flag,                            [-]
  integer,  pointer :: Hours4LiterfalAftMature_brch(:,:) => null()  !branch phenology flag,                            [h]
  integer,  pointer :: KHiestGroLeafNode_brch(:,:)       => null()  !leaf growth stage counter,                        [-]
  integer,  pointer :: iPlantPhenolType_pft(:)           => null()  !climate signal for phenological progress: none,   temperature, water stress,[-]
  integer,  pointer :: Days4FalseBreak_pft(:)            => null()  !accumulated days to singifying false break, [day]
  integer,  pointer :: iPlantTurnoverPattern_pft(:)      => null()  !phenologically-driven above-ground turnover: all, foliar only, none,[-]
  integer,  pointer :: iPlant2ndGrothPattern_pft(:)       => null() !does the plant express secondary growth, [-]
  integer,  pointer :: isPlantShootAlive_pft(:)           => null()  !flag to detect canopy death,[-]
  integer,  pointer :: iPlantPhenolPattern_pft(:)        => null()  !plant growth habit: annual or perennial,[-]
  integer,  pointer :: isPlantRootAlive_pft(:)            => null()  !flag to detect root system death,[-]
  integer,  pointer :: iPlantDevelopPattern_pft(:)       => null()  !plant growth habit (determinate or indeterminate),[-]
  integer,  pointer :: iPlantPhotoperiodType_pft(:)      => null()  !photoperiod type (neutral, long day, short day),[-]
  integer,  pointer :: iPlantRootProfile_pft(:)          => null()  !plant growth type (vascular, non-vascular),[-]
  integer, pointer :: iPlantSnowIntercepType_pft(:) => null() ! snow interception class [0:4, PFT]
  integer,  pointer :: doInitPlant_pft(:)                => null()  !PFT initialization flag:0=no,1=yes,[-]
  integer,  pointer :: KLowestGroLeafNode_brch(:,:)      => null()  !leaf growth stage counter,                        [-]
  integer,  pointer :: iPlantCalendar_brch(:,:,:)        => null()  !plant growth stage,                               [-]

  real(r8), pointer :: fRootGrowPSISense_pvr(:,:,:)           => null()  !water stress to plant root growth, [-]
  real(r8), pointer :: TotalNodeNumNormByMatgrp_brch(:,:)     => null()  !normalized node number during vegetative growth stages,                      [-]
  real(r8), pointer :: TotReproNodeNumNormByMatrgrp_brch(:,:) => null()  !normalized node number during reproductive growth stages,                    [-]
  real(r8), pointer :: LeafNumberAtFloralInit_brch(:,:)       => null()  !leaf number at floral initiation,                                            [-]
  real(r8), pointer :: Hours4LenthenPhotoPeriod_brch(:,:)     => null()  !initial heat requirement for spring leafout/dehardening,                     [h]
  real(r8), pointer :: Hours4ShortenPhotoPeriod_brch(:,:)     => null()  !initial cold requirement for autumn leafoff/hardening,                       [h]
  real(r8), pointer :: Hours4Leafout_brch(:,:)                => null()  !heat requirement for spring leafout/dehardening,                             [h]
  real(r8), pointer :: HourReq4LeafOut_brch(:,:)              => null()  !hours above threshold temperature required for spring leafout/dehardening,   [-]
  real(r8), pointer :: Hours4LeafOff_brch(:,:)                => null()  !cold requirement for autumn leafoff/hardening,                               [h]
  real(r8), pointer :: HourReq4LeafOff_brch(:,:)              => null()  !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8), pointer :: Hours2LeafOut_brch(:,:)                => null()  !counter for mobilizing nonstructural C during spring leafout/dehardening,    [h]
  real(r8), pointer :: HoursDoingRemob_brch(:,:)              => null()  !counter for mobilizing nonstructural C during autumn leafoff/hardening,      [h]
  real(r8), pointer :: HourlyNodeNumNormByMatgrp_brch(:,:)    => null()  !gain in normalized node number during vegetative growth stages,              [h-1]
  real(r8), pointer :: dReproNodeNumNormByMatG_brch(:,:)      => null()  !gain in normalized node number during reproductive growth stages,            [h-1]
  real(r8), pointer :: MatureGroup_brch(:,:)                  => null()  !plant maturity group,                                                        [-]
  real(r8), pointer :: NodeNumNormByMatgrp_brch(:,:)          => null()  !normalized node number during vegetative growth stages,                      [-]
  real(r8), pointer :: ReprodNodeNumNormByMatrgrp_brch(:,:)   => null()  !normalized node number during reproductive growth stages,                    [-]
  real(r8), pointer :: HourFailGrainFill_brch(:,:)            => null()  !flag to detect physiological maturity from  grain fill,                      [-]
  real(r8), pointer :: MatureGroup_pft(:)                     => null()  !acclimated plant maturity group,                                             [-]
  real(r8), pointer :: fNCLFW_brch(:,:)                       => null()  !NC ratio of growing leaf on branch, [gN/gC]
  real(r8), pointer :: fPCLFW_brch(:,:)                       => null()  !PC ratio of growing leaf on branch, [gP/gC]
  real(r8), pointer :: fNCLFW_pft(:)                          => null()  !NC ratio of growing leaf, [gN/gC]
  real(r8), pointer :: fPCLFW_pft(:)                          => null()  !PC ratio of growing leaf, [gP/gC]

  contains
    procedure, public :: Init    =>  plt_pheno_init
    procedure, public :: Destroy =>  plt_pheno_destroy
  end type plant_pheno_type


  type(plant_pheno_type)    , public, target :: plt_pheno     !plant phenology

contains

  subroutine plt_pheno_init(this)
  ! Description:
  ! Allocate phenology state and PFT trait identifiers for the plant API.
  implicit none
  class(plant_pheno_type) :: this


  allocate(this%TCChill4Seed_pft(JP1));this%TCChill4Seed_pft=spval
  allocate(this%TempOffset_pft(JP1));this%TempOffset_pft=spval
  allocate(this%MatureGroup_pft(JP1));this%MatureGroup_pft=spval
  allocate(this%CanPhenoMoistStress_pft(JP1));this%CanPhenoMoistStress_pft=1._r8
  allocate(this%CanPhenoTempStress_pft(JP1));this%CanPhenoTempStress_pft=1._r8
  allocate(this%PlantO2Stress_pft(JP1));this%PlantO2Stress_pft=spval
  allocate(this%NonstCMinConc2InitBranch_pft(JP1));this%NonstCMinConc2InitBranch_pft=spval
  allocate(this%NonstCMinCon2InitRoot_pft(JP1));this%NonstCMinCon2InitRoot_pft=spval
  allocate(this%fTCanopyGroth_pft(JP1));this%fTCanopyGroth_pft=spval
  allocate(this%MorphogenBase_pft(JP1));this%MorphogenBase_pft=spval
  allocate(this%TC4LeafOut_pft(JP1));this%TC4LeafOut_pft=spval
  allocate(this%TCGroth_pft(JP1));this%TCGroth_pft=spval
  allocate(this%TKGroth_pft(JP1));this%TKGroth_pft=spval
  allocate(this%TC4LeafOff_pft(JP1));this%TC4LeafOff_pft=spval
  allocate(this%HoursTooLowPsiCan_pft(JP1));this%HoursTooLowPsiCan_pft=spval
  allocate(this%LeafElmntRemobFlx_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LeafElmntRemobFlx_brch=spval
  allocate(this%LeafSenescInitialElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LeafSenescInitialElms_brch=-1._r8
  allocate(this%PetolShethChemElmRemobFlx_brch(NumPlantChemElms,MaxNumBranches,JP1));this%PetolShethChemElmRemobFlx_brch=spval
  allocate(this%PetolSenescInitialElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%PetolSenescInitialElms_brch=-1._r8
  allocate(this%fNCLFW_pft(JP1)); this%fNCLFW_pft=0._r8
  allocate(this%fPCLFW_pft(JP1)); this%fPCLFW_pft=0._r8
  allocate(this%fTgrowRootP_vr(JZ1,JP1));this%fTgrowRootP_vr=spval
  allocate(this%GrainFillRate25C_pft(JP1));this%GrainFillRate25C_pft=spval
  allocate(this%ShootRootNonstElmConduts_pft(JP1));this%ShootRootNonstElmConduts_pft=spval
  allocate(this%SeedTempSens_pft(JP1));this%SeedTempSens_pft=spval
  allocate(this%HighTempLimitSeed_pft(JP1));this%HighTempLimitSeed_pft=spval
  allocate(this%PlantInitThermoAdaptZone_pft(JP1));this%PlantInitThermoAdaptZone_pft=spval
  allocate(this%rPlantThermoAdaptZone_pft(JP1));this%rPlantThermoAdaptZone_pft=0
  allocate(this%IsPlantActive_pft(JP1));this%IsPlantActive_pft=0
  allocate(this%iMaintPlantTrait_pft(JP1)); this%iMaintPlantTrait_pft=0
  allocate(this%iPlantStateLive_pft(JP1));this%iPlantStateLive_pft=0
  allocate(this%FireReSet_pft(JP1)); this%FireReSet_pft=ifalse
  allocate(this%NetCumElmntFlx2Plant_pft(NumPlantChemElms,JP1));this%NetCumElmntFlx2Plant_pft=spval
  allocate(this%MatureGroup_brch(MaxNumBranches,JP1));this%MatureGroup_brch=spval
  allocate(this%isPlantShootAlive_pft(JP1));this%isPlantShootAlive_pft=0
  allocate(this%Hours2LeafOut_brch(MaxNumBranches,JP1));this%Hours2LeafOut_brch=spval
  allocate(this%HoursDoingRemob_brch(MaxNumBranches,JP1));this%HoursDoingRemob_brch=spval
  allocate(this%HourlyNodeNumNormByMatgrp_brch(MaxNumBranches,JP1));this%HourlyNodeNumNormByMatgrp_brch=spval
  allocate(this%dReproNodeNumNormByMatG_brch(MaxNumBranches,JP1));this%dReproNodeNumNormByMatG_brch=spval
  allocate(this%NodeNumNormByMatgrp_brch(MaxNumBranches,JP1));this%NodeNumNormByMatgrp_brch=spval
  allocate(this%ReprodNodeNumNormByMatrgrp_brch(MaxNumBranches,JP1));this%ReprodNodeNumNormByMatrgrp_brch=spval
  allocate(this%HourFailGrainFill_brch(MaxNumBranches,JP1));this%HourFailGrainFill_brch=spval
  allocate(this%fNCLFW_brch(MaxNumBranches,JP1)); this%fNCLFW_brch=spval
  allocate(this%fPCLFW_brch(MaxNumBranches,JP1)); this%fPCLFW_brch=spval
  allocate(this%iPlantPhenolType_pft(JP1));this%iPlantPhenolType_pft=0
  allocate(this%Days4FalseBreak_pft(JP1)); this%Days4FalseBreak_pft = 0
  allocate(this%iEmbryophyteType_pft(JP1)); this%iEmbryophyteType_pft=0
  allocate(this%iPlantPhenolPattern_pft(JP1));this%iPlantPhenolPattern_pft=0
  allocate(this%iPlantTurnoverPattern_pft(JP1));this%iPlantTurnoverPattern_pft=0
  allocate(this%iPlant2ndGrothPattern_pft(JP1));this%iPlant2ndGrothPattern_pft=0
  allocate(this%isPlantRootAlive_pft(JP1));this%isPlantRootAlive_pft=0
  allocate(this%iPlantDevelopPattern_pft(JP1));this%iPlantDevelopPattern_pft=0
  allocate(this%iPlantPhotoperiodType_pft(JP1));this%iPlantPhotoperiodType_pft=0
  allocate(this%doInitPlant_pft(JP1));this%doInitPlant_pft=0
  allocate(this%doReSeed_pft(JP1)); this%doReSeed_pft=.false.
  allocate(this%iPlantRootProfile_pft(JP1));this%iPlantRootProfile_pft=0
  allocate(this%iPlantSnowIntercepType_pft(JP1));this%iPlantSnowIntercepType_pft=0
  allocate(this%KHiestGroLeafNode_brch(MaxNumBranches,JP1));this%KHiestGroLeafNode_brch=0
  allocate(this%KLowestGroLeafNode_brch(MaxNumBranches,JP1));this%KLowestGroLeafNode_brch=0
  allocate(this%fRootGrowPSISense_pvr(jroots,JZ1,JP1)); this%fRootGrowPSISense_pvr=spval
  allocate(this%iPlantCalendar_brch(NumGrowthStages,MaxNumBranches,JP1));this%iPlantCalendar_brch=0
  allocate(this%TotalNodeNumNormByMatgrp_brch(MaxNumBranches,JP1));this%TotalNodeNumNormByMatgrp_brch=spval
  allocate(this%TotReproNodeNumNormByMatrgrp_brch(MaxNumBranches,JP1));this%TotReproNodeNumNormByMatrgrp_brch=spval
  allocate(this%LeafNumberAtFloralInit_brch(MaxNumBranches,JP1));this%LeafNumberAtFloralInit_brch=spval
  allocate(this%RateRefLeafAppearance_pft(JP1));this%RateRefLeafAppearance_pft=spval
  allocate(this%RefNodeInitRate_pft(JP1));this%RefNodeInitRate_pft=spval
  allocate(this%CriticPhotoPeriod_pft(JP1));this%CriticPhotoPeriod_pft=spval
  allocate(this%PhotoPeriodSens_pft(JP1));this%PhotoPeriodSens_pft=spval
  allocate(this%isPlantBranchAlive_brch(MaxNumBranches,JP1));this%isPlantBranchAlive_brch=iFalse
  allocate(this%doRemobilization_brch(MaxNumBranches,JP1));this%doRemobilization_brch=0
  allocate(this%doPlantLeaveOff_brch(MaxNumBranches,JP1));this%doPlantLeaveOff_brch=0
  allocate(this%EnablePlantLeafOut_brch(MaxNumBranches,JP1));this%EnablePlantLeafOut_brch=0
  allocate(this%doInitLeafOut_brch(MaxNumBranches,JP1));this%doInitLeafOut_brch=0
  allocate(this%doSenescence_brch(MaxNumBranches,JP1));this%doSenescence_brch=0
  allocate(this%Prep4Literfall_brch(MaxNumBranches,JP1));this%Prep4Literfall_brch=0
  allocate(this%Hours4LiterfalAftMature_brch(MaxNumBranches,JP1));this%Hours4LiterfalAftMature_brch=0
  allocate(this%Hours4LenthenPhotoPeriod_brch(MaxNumBranches,JP1));this%Hours4LenthenPhotoPeriod_brch=spval
  allocate(this%Hours4ShortenPhotoPeriod_brch(MaxNumBranches,JP1));this%Hours4ShortenPhotoPeriod_brch=spval
  allocate(this%Hours4Leafout_brch(MaxNumBranches,JP1));this%Hours4Leafout_brch=spval
  allocate(this%HourReq4LeafOut_brch(NumCanopyLayers1,JP1));this%HourReq4LeafOut_brch=spval
  allocate(this%Hours4LeafOff_brch(MaxNumBranches,JP1));this%Hours4LeafOff_brch=spval
  allocate(this%HourReq4LeafOff_brch(NumCanopyLayers1,JP1));this%HourReq4LeafOff_brch=spval

  end subroutine plt_pheno_init

  subroutine plt_pheno_destroy(this)
  ! Description:
  ! Release every array owned by this PlantPhenologyAPIData workspace.
  class(plant_pheno_type), intent(inout) :: this ! owned API state; component units and dimensions as declared above
  character(len=*), parameter :: subname='plt_pheno_destroy'

  call PrintInfo('beg '//subname)
  if(associated(this%fTgrowRootP_vr))deallocate(this%fTgrowRootP_vr);nullify(this%fTgrowRootP_vr)
  if(associated(this%ShootRootNonstElmConduts_pft))deallocate(this%ShootRootNonstElmConduts_pft);nullify(this%ShootRootNonstElmConduts_pft)
  if(associated(this%GrainFillRate25C_pft))deallocate(this%GrainFillRate25C_pft);nullify(this%GrainFillRate25C_pft)
  if(associated(this%TempOffset_pft))deallocate(this%TempOffset_pft);nullify(this%TempOffset_pft)
  if(associated(this%CanPhenoMoistStress_pft))deallocate(this%CanPhenoMoistStress_pft);nullify(this%CanPhenoMoistStress_pft)
  if(associated(this%CanPhenoTempStress_pft))deallocate(this%CanPhenoTempStress_pft);nullify(this%CanPhenoTempStress_pft)
  if(associated(this%PlantO2Stress_pft))deallocate(this%PlantO2Stress_pft);nullify(this%PlantO2Stress_pft)
  if(associated(this%NonstCMinConc2InitBranch_pft))deallocate(this%NonstCMinConc2InitBranch_pft);nullify(this%NonstCMinConc2InitBranch_pft)
  if(associated(this%NonstCMinCon2InitRoot_pft))deallocate(this%NonstCMinCon2InitRoot_pft);nullify(this%NonstCMinCon2InitRoot_pft)
  if(associated(this%LeafElmntRemobFlx_brch))deallocate(this%LeafElmntRemobFlx_brch);nullify(this%LeafElmntRemobFlx_brch)
  if(associated(this%LeafSenescInitialElms_brch))deallocate(this%LeafSenescInitialElms_brch);nullify(this%LeafSenescInitialElms_brch)
  if(associated(this%PetolShethChemElmRemobFlx_brch))deallocate(this%PetolShethChemElmRemobFlx_brch);nullify(this%PetolShethChemElmRemobFlx_brch)
  if(associated(this%PetolSenescInitialElms_brch))deallocate(this%PetolSenescInitialElms_brch);nullify(this%PetolSenescInitialElms_brch)
  if(associated(this%TC4LeafOut_pft))deallocate(this%TC4LeafOut_pft);nullify(this%TC4LeafOut_pft)
  if(associated(this%TCGroth_pft))deallocate(this%TCGroth_pft);nullify(this%TCGroth_pft)
  if(associated(this%TC4LeafOff_pft))deallocate(this%TC4LeafOff_pft);nullify(this%TC4LeafOff_pft)
  if(associated(this%TKGroth_pft))deallocate(this%TKGroth_pft);nullify(this%TKGroth_pft)
  if(associated(this%fTCanopyGroth_pft))deallocate(this%fTCanopyGroth_pft);nullify(this%fTCanopyGroth_pft)
  if(associated(this%MorphogenBase_pft))deallocate(this%MorphogenBase_pft);nullify(this%MorphogenBase_pft)
  if(associated(this%HoursTooLowPsiCan_pft))deallocate(this%HoursTooLowPsiCan_pft);nullify(this%HoursTooLowPsiCan_pft)
  if(associated(this%TCChill4Seed_pft))deallocate(this%TCChill4Seed_pft);nullify(this%TCChill4Seed_pft)
  if(associated(this%rPlantThermoAdaptZone_pft))deallocate(this%rPlantThermoAdaptZone_pft);nullify(this%rPlantThermoAdaptZone_pft)
  if(associated(this%PlantInitThermoAdaptZone_pft))deallocate(this%PlantInitThermoAdaptZone_pft);nullify(this%PlantInitThermoAdaptZone_pft)
  if(associated(this%HighTempLimitSeed_pft))deallocate(this%HighTempLimitSeed_pft);nullify(this%HighTempLimitSeed_pft)
  if(associated(this%SeedTempSens_pft))deallocate(this%SeedTempSens_pft);nullify(this%SeedTempSens_pft)
  if(associated(this%NetCumElmntFlx2Plant_pft))deallocate(this%NetCumElmntFlx2Plant_pft);nullify(this%NetCumElmntFlx2Plant_pft)
  if(associated(this%RefNodeInitRate_pft))deallocate(this%RefNodeInitRate_pft);nullify(this%RefNodeInitRate_pft)
  if(associated(this%RateRefLeafAppearance_pft))deallocate(this%RateRefLeafAppearance_pft);nullify(this%RateRefLeafAppearance_pft)
  if(associated(this%CriticPhotoPeriod_pft))deallocate(this%CriticPhotoPeriod_pft);nullify(this%CriticPhotoPeriod_pft)
  if(associated(this%PhotoPeriodSens_pft))deallocate(this%PhotoPeriodSens_pft);nullify(this%PhotoPeriodSens_pft)
  if(associated(this%iEmbryophyteType_pft))deallocate(this%iEmbryophyteType_pft);nullify(this%iEmbryophyteType_pft)
  if(associated(this%iPlantStateLive_pft))deallocate(this%iPlantStateLive_pft);nullify(this%iPlantStateLive_pft)
  if(associated(this%FireReSet_pft))deallocate(this%FireReSet_pft);nullify(this%FireReSet_pft)
  if(associated(this%iMaintPlantTrait_pft))deallocate(this%iMaintPlantTrait_pft);nullify(this%iMaintPlantTrait_pft)
  if(associated(this%IsPlantActive_pft))deallocate(this%IsPlantActive_pft);nullify(this%IsPlantActive_pft)
  if(associated(this%isPlantBranchAlive_brch))deallocate(this%isPlantBranchAlive_brch);nullify(this%isPlantBranchAlive_brch)
  if(associated(this%doRemobilization_brch))deallocate(this%doRemobilization_brch);nullify(this%doRemobilization_brch)
  if(associated(this%doPlantLeaveOff_brch))deallocate(this%doPlantLeaveOff_brch);nullify(this%doPlantLeaveOff_brch)
  if(associated(this%EnablePlantLeafOut_brch))deallocate(this%EnablePlantLeafOut_brch);nullify(this%EnablePlantLeafOut_brch)
  if(associated(this%doInitLeafOut_brch))deallocate(this%doInitLeafOut_brch);nullify(this%doInitLeafOut_brch)
  if(associated(this%doReSeed_pft))deallocate(this%doReSeed_pft);nullify(this%doReSeed_pft)
  if(associated(this%doSenescence_brch))deallocate(this%doSenescence_brch);nullify(this%doSenescence_brch)
  if(associated(this%Prep4Literfall_brch))deallocate(this%Prep4Literfall_brch);nullify(this%Prep4Literfall_brch)
  if(associated(this%Hours4LiterfalAftMature_brch))deallocate(this%Hours4LiterfalAftMature_brch);nullify(this%Hours4LiterfalAftMature_brch)
  if(associated(this%KHiestGroLeafNode_brch))deallocate(this%KHiestGroLeafNode_brch);nullify(this%KHiestGroLeafNode_brch)
  if(associated(this%iPlantPhenolType_pft))deallocate(this%iPlantPhenolType_pft);nullify(this%iPlantPhenolType_pft)
  if(associated(this%Days4FalseBreak_pft))deallocate(this%Days4FalseBreak_pft);nullify(this%Days4FalseBreak_pft)
  if(associated(this%iPlantTurnoverPattern_pft))deallocate(this%iPlantTurnoverPattern_pft);nullify(this%iPlantTurnoverPattern_pft)
  if(associated(this%iPlant2ndGrothPattern_pft))deallocate(this%iPlant2ndGrothPattern_pft);nullify(this%iPlant2ndGrothPattern_pft)
  if(associated(this%isPlantShootAlive_pft))deallocate(this%isPlantShootAlive_pft);nullify(this%isPlantShootAlive_pft)
  if(associated(this%iPlantPhenolPattern_pft))deallocate(this%iPlantPhenolPattern_pft);nullify(this%iPlantPhenolPattern_pft)
  if(associated(this%isPlantRootAlive_pft))deallocate(this%isPlantRootAlive_pft);nullify(this%isPlantRootAlive_pft)
  if(associated(this%iPlantDevelopPattern_pft))deallocate(this%iPlantDevelopPattern_pft);nullify(this%iPlantDevelopPattern_pft)
  if(associated(this%iPlantPhotoperiodType_pft))deallocate(this%iPlantPhotoperiodType_pft);nullify(this%iPlantPhotoperiodType_pft)
  if(associated(this%iPlantRootProfile_pft))deallocate(this%iPlantRootProfile_pft);nullify(this%iPlantRootProfile_pft)
  if(associated(this%iPlantSnowIntercepType_pft))deallocate(this%iPlantSnowIntercepType_pft);nullify(this%iPlantSnowIntercepType_pft)
  if(associated(this%doInitPlant_pft))deallocate(this%doInitPlant_pft);nullify(this%doInitPlant_pft)
  if(associated(this%KLowestGroLeafNode_brch))deallocate(this%KLowestGroLeafNode_brch);nullify(this%KLowestGroLeafNode_brch)
  if(associated(this%iPlantCalendar_brch))deallocate(this%iPlantCalendar_brch);nullify(this%iPlantCalendar_brch)
  if(associated(this%fRootGrowPSISense_pvr))deallocate(this%fRootGrowPSISense_pvr);nullify(this%fRootGrowPSISense_pvr)
  if(associated(this%TotalNodeNumNormByMatgrp_brch))deallocate(this%TotalNodeNumNormByMatgrp_brch);nullify(this%TotalNodeNumNormByMatgrp_brch)
  if(associated(this%TotReproNodeNumNormByMatrgrp_brch))deallocate(this%TotReproNodeNumNormByMatrgrp_brch);nullify(this%TotReproNodeNumNormByMatrgrp_brch)
  if(associated(this%LeafNumberAtFloralInit_brch))deallocate(this%LeafNumberAtFloralInit_brch);nullify(this%LeafNumberAtFloralInit_brch)
  if(associated(this%Hours4LenthenPhotoPeriod_brch))deallocate(this%Hours4LenthenPhotoPeriod_brch);nullify(this%Hours4LenthenPhotoPeriod_brch)
  if(associated(this%Hours4ShortenPhotoPeriod_brch))deallocate(this%Hours4ShortenPhotoPeriod_brch);nullify(this%Hours4ShortenPhotoPeriod_brch)
  if(associated(this%Hours4Leafout_brch))deallocate(this%Hours4Leafout_brch);nullify(this%Hours4Leafout_brch)
  if(associated(this%HourReq4LeafOut_brch))deallocate(this%HourReq4LeafOut_brch);nullify(this%HourReq4LeafOut_brch)
  if(associated(this%Hours4LeafOff_brch))deallocate(this%Hours4LeafOff_brch);nullify(this%Hours4LeafOff_brch)
  if(associated(this%HourReq4LeafOff_brch))deallocate(this%HourReq4LeafOff_brch);nullify(this%HourReq4LeafOff_brch)
  if(associated(this%Hours2LeafOut_brch))deallocate(this%Hours2LeafOut_brch);nullify(this%Hours2LeafOut_brch)
  if(associated(this%HoursDoingRemob_brch))deallocate(this%HoursDoingRemob_brch);nullify(this%HoursDoingRemob_brch)
  if(associated(this%HourlyNodeNumNormByMatgrp_brch))deallocate(this%HourlyNodeNumNormByMatgrp_brch);nullify(this%HourlyNodeNumNormByMatgrp_brch)
  if(associated(this%dReproNodeNumNormByMatG_brch))deallocate(this%dReproNodeNumNormByMatG_brch);nullify(this%dReproNodeNumNormByMatG_brch)
  if(associated(this%MatureGroup_brch))deallocate(this%MatureGroup_brch);nullify(this%MatureGroup_brch)
  if(associated(this%NodeNumNormByMatgrp_brch))deallocate(this%NodeNumNormByMatgrp_brch);nullify(this%NodeNumNormByMatgrp_brch)
  if(associated(this%ReprodNodeNumNormByMatrgrp_brch))deallocate(this%ReprodNodeNumNormByMatrgrp_brch);nullify(this%ReprodNodeNumNormByMatrgrp_brch)
  if(associated(this%HourFailGrainFill_brch))deallocate(this%HourFailGrainFill_brch);nullify(this%HourFailGrainFill_brch)
  if(associated(this%MatureGroup_pft))deallocate(this%MatureGroup_pft);nullify(this%MatureGroup_pft)
  if(associated(this%fNCLFW_brch))deallocate(this%fNCLFW_brch);nullify(this%fNCLFW_brch)
  if(associated(this%fPCLFW_brch))deallocate(this%fPCLFW_brch);nullify(this%fPCLFW_brch)
  if(associated(this%fNCLFW_pft))deallocate(this%fNCLFW_pft);nullify(this%fNCLFW_pft)
  if(associated(this%fPCLFW_pft))deallocate(this%fPCLFW_pft);nullify(this%fPCLFW_pft)
  call PrintInfo('end '//subname)
  end subroutine plt_pheno_destroy
end module PlantPhenologyAPIData
