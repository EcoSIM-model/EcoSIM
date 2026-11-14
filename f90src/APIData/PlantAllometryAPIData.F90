module PlantAllometryAPIData
  ! Description:
  ! Owns the allometry API type, instance and lifecycle procedures.
  use DebugToolMod, only: PrintInfo
  use PlantAPICommonData
  implicit none
  private
  public :: plt_allom_init, plt_allom_destroy

  type, public :: plant_allometry_type
  real(r8), pointer :: CPRTS_pft(:)                     => null()  !root P:C ratio x root growth yield,                      [-]
  real(r8), pointer :: CNRTS_pft(:)                     => null()  !root N:C ratio x root growth yield,                      [-]
  real(r8), pointer :: rNCNodule_pft(:)                 => null()  !nodule N:C ratio,                                        [gN gC-1]
  real(r8), pointer :: rPCNoduler_pft(:)                => null()  !nodule P:C ratio,                                        [gP gC-1]
  real(r8), pointer :: rNCRoot_pft(:)                   => null()  !root N:C ratio,                                          [gN gC-1]
  real(r8), pointer :: rPCRootr_pft(:)                  => null()  !root P:C ratio,                                          [gP gC-1]
  real(r8), pointer :: rProteinC2RootN_pft(:)           => null()  !Protein C to root N ratio in remobilizable nonstructural biomass,        [-]
  real(r8), pointer :: rProteinC2RootP_pft(:)           => null()  !Protein C to root P ratio in remobilizable nonstructural biomass,        [-]
  real(r8), pointer :: rProteinC2LeafN_pft(:)          => null()  !Protein C to leaf N ratio in remobilizable nonstructural biomass,        [-]
  real(r8), pointer :: rProteinC2LeafP_pft(:)          => null()  !Protein C to leaf P ratio in remobilizable nonstructural biomass,        [-]
  real(r8), pointer :: NoduGrowthYield_pft(:)           => null()  !nodule growth yield,                                     [g g-1]
  real(r8), pointer :: RootBiomGrosYld_pft(:)           => null()  !root growth yield,                                       [g g-1]
  real(r8), pointer :: rPCEar_pft(:)                    => null()  !ear P:C ratio,                                           [gP gC-1]
  real(r8), pointer :: PetolShethBiomGrowthYld_pft(:)      => null()  !sheath growth yield,                                     [g g-1]
  real(r8), pointer :: rPCHusk_pft(:)                   => null()  !husk P:C ratio,                                          [gP gC-1]
  real(r8), pointer :: StalkBiomGrowthYld_pft(:)        => null()  !stalk growth yield,                                      [gC gC-1]
  real(r8), pointer :: HuskBiomGrowthYld_pft(:)         => null()  !husk growth yield,                                       [gC gC-1]
  real(r8), pointer :: ReserveBiomGrowthYld_pft(:)      => null()  !reserve growth yield,                                    [gC gC-1]
  real(r8), pointer :: GrainBiomGrowthYld_pft(:)        => null()  !grain growth yield,                                      [gC gC-1]
  real(r8), pointer :: EarBiomGrowthYld_pft(:)          => null()  !ear growth yield,                                        [gC gC-1]
  real(r8), pointer :: rNCHusk_pft(:)                   => null()  !husk N:C ratio,                                          [gN gC-1]
  real(r8), pointer :: rNCReserve_pft(:)                => null()  !reserve N:C ratio,                                       [gN gC-1]
  real(r8), pointer :: rNCEar_pft(:)                    => null()  !ear N:C ratio,                                           [gN gC-1]
  real(r8), pointer :: rPCReserve_pft(:)                => null()  !reserve P:C ratio,                                       [gP gC-1]
  real(r8), pointer :: rPCGrain_pft(:)                      => null()  !grain P:C ratio,                                         [gP gP-1]
  real(r8), pointer :: rNCStalk_pft(:)                  => null()  !stalk N:C ratio,                                         [gN gC-1]
  real(r8), pointer :: rECLiveCRoot_pft(:,:)           => null()   ! element:C ratio of live coarse root, [gE gC-1]
  real(r8), pointer :: rECDeadCRoot_pft(:,:)           => null()   ! element:C ratio of dead coarse root, [gE gC-1]
  real(r8), pointer :: rNCLigRoot_pft(:)               => null()  !NC ratio of lignified root, [gN gC-1]
  real(r8), pointer :: rPCLigRoot_pft(:)               => null()  !PC ratio of lignified root, [gP gC-1]
  real(r8), pointer :: FracLeafShethElmAlloc2Litr(:,:)  => null()  !woody element allocation, [-]
  real(r8), pointer :: FracPetolShethAlloc2Litr(:,:) => null()  !leaf element allocation,[-]
  real(r8), pointer :: FracRootElmAllocm(:,:)       => null()  !C woody fraction in root,[-]
  real(r8), pointer :: FracWoodStalkElmAlloc2Litr(:,:)  => null()  !woody element allocation,[-]
  real(r8), pointer :: LeafBiomGrowthYld_pft(:)         => null()  !leaf growth yield,                                       [g g-1]
  real(r8), pointer :: rNCGrain_pft(:)                      => null()  !grain N:C ratio,                                         [g g-1]
  real(r8), pointer :: rPCLeaf_pft(:)                      => null()  !maximum leaf P:C ratio,                                  [g g-1]
  real(r8), pointer :: rPCSheath_pft(:)                     => null()  !sheath P:C ratio,                                        [g g-1]
  real(r8), pointer :: rNCSheath_pft(:)                     => null()  !sheath N:C ratio,                                        [g g-1]
  real(r8), pointer :: rPCStalk_pft(:)                  => null()  !stalk P:C ratio,                                         [g g-1]
  real(r8), pointer :: rNCLeaf_pft(:)                      => null()  !maximum leaf N:C ratio,                                  [g g-1]
  real(r8), pointer :: SingleGrainMeanBiomC_brch(:,:)     => null()  !potential carbon mass per grain, [gC seed-1]
  real(r8), pointer :: FracGroth2Node_pft(:)            => null()  !parameter for allocation of growth to nodes,             [-]
  real(r8), pointer ::RootProteinCMax_pft(:)      => null()  !reference root protein N, [gN g-1]

  contains
    procedure, public :: Init => plt_allom_init
    procedure, public :: Destroy => plt_allom_destroy
  end type plant_allometry_type


  type(plant_allometry_type), public, target :: plt_allom     !plant allometric parameters

contains

  subroutine plt_allom_init(this)
  implicit none
  class(plant_allometry_type) :: this


  allocate(this%RootProteinCMax_pft(JP1));this%RootProteinCMax_pft=spval
  allocate(this%FracGroth2Node_pft(JP1));this%FracGroth2Node_pft=spval
  allocate(this%SingleGrainMeanBiomC_brch(MaxNumBranches,JP1));this%SingleGrainMeanBiomC_brch=spval
  allocate(this%NoduGrowthYield_pft(JP1));this%NoduGrowthYield_pft=spval
  allocate(this%RootBiomGrosYld_pft(JP1));this%RootBiomGrosYld_pft=spval
  allocate(this%rPCRootr_pft(JP1));this%rPCRootr_pft=spval
  allocate(this%rProteinC2RootN_pft(JP1));this%rProteinC2RootN_pft=spval
  allocate(this%rProteinC2RootP_pft(JP1));this%rProteinC2RootP_pft=spval
  allocate(this%rProteinC2LeafN_pft(JP1));this%rProteinC2LeafN_pft=spval
  allocate(this%rProteinC2LeafP_pft(JP1));this%rProteinC2LeafP_pft=spval
  allocate(this%CPRTS_pft(JP1));this%CPRTS_pft=spval
  allocate(this%CNRTS_pft(JP1));this%CNRTS_pft=spval
  allocate(this%rNCNodule_pft(JP1));this%rNCNodule_pft=spval
  allocate(this%rPCNoduler_pft(JP1));this%rPCNoduler_pft=spval
  allocate(this%rNCRoot_pft(JP1));this%rNCRoot_pft=spval
  allocate(this%rPCLeaf_pft(JP1));this%rPCLeaf_pft=spval
  allocate(this%rPCSheath_pft(JP1));this%rPCSheath_pft=spval
  allocate(this%rNCLeaf_pft(JP1));this%rNCLeaf_pft=spval
  allocate(this%rNCSheath_pft(JP1));this%rNCSheath_pft=spval
  allocate(this%rNCGrain_pft(JP1));this%rNCGrain_pft=spval
  allocate(this%rPCStalk_pft(JP1));this%rPCStalk_pft=spval
  allocate(this%rNCStalk_pft(JP1));this%rNCStalk_pft=spval
  allocate(this%rECDeadCRoot_pft(NumPlantChemElms,JP1));    this%rECDeadCRoot_pft=spval
  allocate(this%rECLiveCRoot_pft(NumPlantChemElms,JP1)); this%rECLiveCRoot_pft=spval
  allocate(this%rNCLigRoot_pft(JP1));this%rNCLigRoot_pft=spval
  allocate(this%rPCLigRoot_pft(JP1)); this%rPCLigRoot_pft=spval
  allocate(this%rPCGrain_pft(JP1));this%rPCGrain_pft=spval
  allocate(this%rPCEar_pft(JP1));this%rPCEar_pft=spval
  allocate(this%rPCReserve_pft(JP1));this%rPCReserve_pft=spval
  allocate(this%rNCReserve_pft(JP1));this%rNCReserve_pft=spval
  allocate(this%rPCHusk_pft(JP1));this%rPCHusk_pft=spval
  allocate(this%FracWoodStalkElmAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracWoodStalkElmAlloc2Litr=spval
  allocate(this%FracPetolShethAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracPetolShethAlloc2Litr=spval
  allocate(this%FracRootElmAllocm(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracRootElmAllocm=spval
  allocate(this%FracLeafShethElmAlloc2Litr(NumPlantChemElms,NumOfPlantLitrCmplxs));this%FracLeafShethElmAlloc2Litr=spval

  allocate(this%PetolShethBiomGrowthYld_pft(JP1));this%PetolShethBiomGrowthYld_pft=spval
  allocate(this%HuskBiomGrowthYld_pft(JP1));this%HuskBiomGrowthYld_pft=spval
  allocate(this%StalkBiomGrowthYld_pft(JP1));this%StalkBiomGrowthYld_pft=spval
  allocate(this%ReserveBiomGrowthYld_pft(JP1));this%ReserveBiomGrowthYld_pft=spval
  allocate(this%EarBiomGrowthYld_pft(JP1));this%EarBiomGrowthYld_pft=spval
  allocate(this%GrainBiomGrowthYld_pft(JP1));this%GrainBiomGrowthYld_pft=spval
  allocate(this%rNCHusk_pft(JP1));this%rNCHusk_pft=spval
  allocate(this%rNCEar_pft(JP1));this%rNCEar_pft=spval
  allocate(this%LeafBiomGrowthYld_pft(JP1));this%LeafBiomGrowthYld_pft=spval

  end subroutine plt_allom_init

  subroutine plt_allom_destroy(this)
  ! Description:
  ! Release every array owned by this PlantAllometryAPIData workspace.
  class(plant_allometry_type), intent(inout) :: this ! owned API state; component units and dimensions as declared above
  character(len=*), parameter :: subname='plt_allom_destroy'

  call PrintInfo('beg '//subname)
  if(associated(this%CPRTS_pft))deallocate(this%CPRTS_pft);nullify(this%CPRTS_pft)
  if(associated(this%CNRTS_pft))deallocate(this%CNRTS_pft);nullify(this%CNRTS_pft)
  if(associated(this%rNCNodule_pft))deallocate(this%rNCNodule_pft);nullify(this%rNCNodule_pft)
  if(associated(this%rPCNoduler_pft))deallocate(this%rPCNoduler_pft);nullify(this%rPCNoduler_pft)
  if(associated(this%rNCRoot_pft))deallocate(this%rNCRoot_pft);nullify(this%rNCRoot_pft)
  if(associated(this%rPCRootr_pft))deallocate(this%rPCRootr_pft);nullify(this%rPCRootr_pft)
  if(associated(this%rProteinC2RootN_pft))deallocate(this%rProteinC2RootN_pft);nullify(this%rProteinC2RootN_pft)
  if(associated(this%rProteinC2RootP_pft))deallocate(this%rProteinC2RootP_pft);nullify(this%rProteinC2RootP_pft)
  if(associated(this%rProteinC2LeafN_pft))deallocate(this%rProteinC2LeafN_pft);nullify(this%rProteinC2LeafN_pft)
  if(associated(this%rProteinC2LeafP_pft))deallocate(this%rProteinC2LeafP_pft);nullify(this%rProteinC2LeafP_pft)
  if(associated(this%NoduGrowthYield_pft))deallocate(this%NoduGrowthYield_pft);nullify(this%NoduGrowthYield_pft)
  if(associated(this%RootBiomGrosYld_pft))deallocate(this%RootBiomGrosYld_pft);nullify(this%RootBiomGrosYld_pft)
  if(associated(this%rPCEar_pft))deallocate(this%rPCEar_pft);nullify(this%rPCEar_pft)
  if(associated(this%PetolShethBiomGrowthYld_pft))deallocate(this%PetolShethBiomGrowthYld_pft);nullify(this%PetolShethBiomGrowthYld_pft)
  if(associated(this%rPCHusk_pft))deallocate(this%rPCHusk_pft);nullify(this%rPCHusk_pft)
  if(associated(this%StalkBiomGrowthYld_pft))deallocate(this%StalkBiomGrowthYld_pft);nullify(this%StalkBiomGrowthYld_pft)
  if(associated(this%HuskBiomGrowthYld_pft))deallocate(this%HuskBiomGrowthYld_pft);nullify(this%HuskBiomGrowthYld_pft)
  if(associated(this%ReserveBiomGrowthYld_pft))deallocate(this%ReserveBiomGrowthYld_pft);nullify(this%ReserveBiomGrowthYld_pft)
  if(associated(this%GrainBiomGrowthYld_pft))deallocate(this%GrainBiomGrowthYld_pft);nullify(this%GrainBiomGrowthYld_pft)
  if(associated(this%EarBiomGrowthYld_pft))deallocate(this%EarBiomGrowthYld_pft);nullify(this%EarBiomGrowthYld_pft)
  if(associated(this%rNCHusk_pft))deallocate(this%rNCHusk_pft);nullify(this%rNCHusk_pft)
  if(associated(this%rNCReserve_pft))deallocate(this%rNCReserve_pft);nullify(this%rNCReserve_pft)
  if(associated(this%rNCEar_pft))deallocate(this%rNCEar_pft);nullify(this%rNCEar_pft)
  if(associated(this%rPCReserve_pft))deallocate(this%rPCReserve_pft);nullify(this%rPCReserve_pft)
  if(associated(this%rPCGrain_pft))deallocate(this%rPCGrain_pft);nullify(this%rPCGrain_pft)
  if(associated(this%rNCStalk_pft))deallocate(this%rNCStalk_pft);nullify(this%rNCStalk_pft)
  if(associated(this%rECLiveCRoot_pft))deallocate(this%rECLiveCRoot_pft);nullify(this%rECLiveCRoot_pft)
  if(associated(this%rECDeadCRoot_pft))deallocate(this%rECDeadCRoot_pft);nullify(this%rECDeadCRoot_pft)
  if(associated(this%rNCLigRoot_pft))deallocate(this%rNCLigRoot_pft);nullify(this%rNCLigRoot_pft)
  if(associated(this%rPCLigRoot_pft))deallocate(this%rPCLigRoot_pft);nullify(this%rPCLigRoot_pft)
  if(associated(this%FracLeafShethElmAlloc2Litr))deallocate(this%FracLeafShethElmAlloc2Litr);nullify(this%FracLeafShethElmAlloc2Litr)
  if(associated(this%FracPetolShethAlloc2Litr))deallocate(this%FracPetolShethAlloc2Litr);nullify(this%FracPetolShethAlloc2Litr)
  if(associated(this%FracRootElmAllocm))deallocate(this%FracRootElmAllocm);nullify(this%FracRootElmAllocm)
  if(associated(this%FracWoodStalkElmAlloc2Litr))deallocate(this%FracWoodStalkElmAlloc2Litr);nullify(this%FracWoodStalkElmAlloc2Litr)
  if(associated(this%LeafBiomGrowthYld_pft))deallocate(this%LeafBiomGrowthYld_pft);nullify(this%LeafBiomGrowthYld_pft)
  if(associated(this%rNCGrain_pft))deallocate(this%rNCGrain_pft);nullify(this%rNCGrain_pft)
  if(associated(this%rPCLeaf_pft))deallocate(this%rPCLeaf_pft);nullify(this%rPCLeaf_pft)
  if(associated(this%rPCSheath_pft))deallocate(this%rPCSheath_pft);nullify(this%rPCSheath_pft)
  if(associated(this%rNCSheath_pft))deallocate(this%rNCSheath_pft);nullify(this%rNCSheath_pft)
  if(associated(this%rPCStalk_pft))deallocate(this%rPCStalk_pft);nullify(this%rPCStalk_pft)
  if(associated(this%rNCLeaf_pft))deallocate(this%rNCLeaf_pft);nullify(this%rNCLeaf_pft)
  if(associated(this%SingleGrainMeanBiomC_brch))deallocate(this%SingleGrainMeanBiomC_brch);nullify(this%SingleGrainMeanBiomC_brch)
  if(associated(this%FracGroth2Node_pft))deallocate(this%FracGroth2Node_pft);nullify(this%FracGroth2Node_pft)
  if(associated(this%RootProteinCMax_pft))deallocate(this%RootProteinCMax_pft);nullify(this%RootProteinCMax_pft)
  call PrintInfo('end '//subname)
  end subroutine plt_allom_destroy
end module PlantAllometryAPIData
