module PlantBiomassAPIData
  ! Owns the biomass API type, instance and lifecycle procedures.
  use PlantAPICommonData
  implicit none
  private
  public :: plt_biom_init, plt_biom_destroy

  type, public :: plant_biom_type
  real(r8), pointer :: TotBegVegE_pft(:,:)                  => null()    !total vegetation biomass at the beginning of time step, [g d-2]
  real(r8), pointer :: TotEndVegE_pft(:,:)                  => null()    !total vegetation biomass at the end of time step, [g d-2]
  real(r8), pointer :: StomatalStress_pft(:)                => null()    !stomatal stress from root turgor (0-1),             [-]
  real(r8), pointer :: RootMycoMassElm_pvr(:,:,:,:)         => null()    !root biomass in chemical elements,                  [g d-2]
  real(r8), pointer :: StandingDeadStrutElms_col(:)         => null()    !total standing dead biomass chemical element,       [g d-2]
  real(r8), pointer :: ZERO4LeafVar_pft(:)                  => null()    !threshold zero for leaf calculation,                [-]
  real(r8), pointer :: LeafProteinC_brch(:,:)               => null()    !Protein C for the branches, [gC protein d-2]
  real(r8), pointer :: LeafProteinCperm2LA_pft(:)           => null()    !Protein C for the plant, [gC protein m-2 leaf area]
  real(r8), pointer :: ZERO4Groth_pft(:)                    => null()    !threshold zero for plang growth calculation,        [-]
  real(r8), pointer :: RootNodulNonstElms_rpvr(:,:,:)       => null()    !root  layer nonstructural element,                  [g d-2]
  real(r8), pointer :: RootNodulStrutElms_rpvr(:,:,:)       => null()    !root layer nodule element,                          [g d-2]
  real(r8), pointer :: CanopyLeafCLyr_pft(:,:)              => null()    !canopy layer leaf C,                                [g d-2]
  real(r8), pointer :: RootMycoNonstElms_pft(:,:,:)         => null()    !nonstructural root-myco chemical element,           [g d-2]
  real(r8), pointer :: RootMyco2ndStrutElms_rpvr(:,:,:,:,:) => null()    !root layer element secondary axes,                  [g d-2]
  real(r8), pointer :: RootMyco2ndStrutElms_pvr(:,:,:,:)   => null()    !secondary structural elements: NR sum or prescribed SP root C, [g element cell-1]
  real(r8), pointer :: RootMyco1stStrutElms_rpvr(:,:,:,:)   => null()    !root layer element primary axes,                    [g d-2]
  real(r8), pointer :: Root1stActStructElms_rpvr(:,:,:,:)   => null()    !root layer active zone element in primary axes, [g d-2]
  real(r8), pointer :: Root1stLigStructElms_rpvr(:,:,:,:)   => null()    !root layer lignified zone element in primary axes, [g d-2]
  real(r8), pointer :: KLigMax_pft(:)                       => null()    !Maximum lignification rate [h-1]
  real(r8), pointer :: KLigMM_pft(:)                        => null()    !Half saturation parameter for coarse root lignification, [h-1]
  real(r8), pointer :: RootMyco1stElm_raxs(:,:,:)           => null()    !root C primary axes,                                [g d-2]
  real(r8), pointer :: StandDeadCompKElms_pft(:,:,:)        => null()    !standing dead element fraction,                     [g d-2]
  real(r8), pointer :: CanopyNonstElmConc_pft(:,:)          => null()    !canopy nonstructural element concentration,         [g d-2]
  real(r8), pointer :: CanopyNonstOMConc_pft(:)             => null()    !canopy C+N+P concentration: derived or prescribed CCPOLT, [g gC-1]
  real(r8), pointer :: CanopyNonstElms_pft(:,:)             => null()    !canopy nonstructural element concentration,         [g d-2]
  real(r8), pointer :: CanopyNodulNonstElms_pft(:,:)        => null()    !canopy nodule nonstructural element,                [g d-2]
  real(r8), pointer :: CanopyNoduleNonstCConc_pft(:)        => null()    !nodule nonstructural C,                             [gC d-2]
  real(r8), pointer :: RootMycoActiveBiomC_pvr(:,:,:)       => null()    !root layer structural C,                            [gC d-2]
  real(r8), pointer :: RootMediumStructElms_rpvr(:,:,:,:)   => null()    !root layer medium size root structrual elements,    [g d-2]
  real(r8), pointer :: RootMedStruct_pvr(:,:,:)             => null()    !root layer element biomass for medium size roots, [g d-2]
  real(r8), pointer :: PopuRootMycoC_pvr(:,:,:)             => null()    !root layer C,                                       [gC d-2]
  real(r8), pointer :: RootProteinC_pvr(:,:,:)              => null()    !root layer protein C,                               [gC d-2]
  real(r8), pointer :: RootProteinConc_rpvr(:,:,:)          => null()    !root layer protein C concentration,                 [g g-1]
  real(r8), pointer :: RootMycoNonstElms_rpvr(:,:,:,:)      => null()    !root  layer nonstructural element,                  [g d-2]
  real(r8), pointer :: RootNonstructElmConc_pvr(:,:,:,:)   => null()    !root  layer nonstructural C concentration,          [g g-1]
  real(r8), pointer :: RootNonstOMConc_pvr(:,:,:)         => null()    !sum of nonstructural C/N/P concentrations by root group/layer/PFT, [g gC-1]
  real(r8), pointer :: LeafPetoNonstElmConc_brch(:,:,:)     => null()    !branch nonstructural C concentration,               [g d-2]
  real(r8), pointer :: StructInternodeElms_brch(:,:,:,:)    => null()    !internode C,                                        [g d-2]
  real(r8), pointer :: LeafElmntNode_brch(:,:,:,:)          => null()    !leaf element,                                       [g d-2]
  real(r8), pointer :: LeafProteinC_node(:,:,:)             => null()    !layer leaf protein N,                               [g d-2]
  real(r8), pointer :: PetolShethElmntNode_brch(:,:,:,:)       => null()    !sheath chemical element,                            [g d-2]
  real(r8), pointer :: PetoleProteinC_node(:,:,:)       => null()    !layer sheath protein C,                             [g d-2]
  real(r8), pointer :: LeafLayerElms_node(:,:,:,:,:)        => null()    !layer leaf element,                                 [g d-2]
  real(r8), pointer :: tCanLeafC_clyr(:)                    => null()    !total leaf carbon mass in canopy layers,            [gC d-2]
  real(r8), pointer :: StandingDeadInitC_pft(:)             => null()    !initial standing dead C,                            [g C m-2]
  real(r8), pointer :: RootElms_pft(:,:)                    => null()    !plant root element mass,                            [g d-2]
  real(r8), pointer :: RootActiveC_pft(:)                   => null()    !derived root C excluding inactive primary C and nodules, [gC population-1]
  real(r8), pointer :: RootNoduleElms_pft(:,:)              => null()    !plant root nodule element mass,                     [g d-2]
  real(r8), pointer :: RootNoduleElmsBeg_pft(:,:)           => null()    !previous time step plant root nodule element mass,  [g d-2]
  real(r8), pointer :: RootElmsBeg_pft(:,:)                 => null()    !plant root element at previous time step,           [g d-2]
  real(r8), pointer :: StandDeadStrutElmsBeg_pft(:,:)       => null()    !standing dead element at previous time step,        [g d-2]
  real(r8), pointer :: Root1stActStruct_pvr(:,:,:)          => null()    !active zone of primary roots, [g d-2]
  real(r8), pointer :: Root1stLigStruct_pvr(:,:,:)          => null()    !lignifed primary root biomass, [g d-2]
  real(r8), pointer :: RootStrutElms_pft(:,:)               => null()    !plant root structural element mass,                 [g d-2]
  real(r8), pointer :: SeedPlantedElm_pft(:,:)              => null()    !plant stored nonstructural chemical elements at planting,           [gC d-2]
  real(r8), pointer :: SeasonalNonstCDayAve_pft(:)          => null()    !daily average seasonal storage C for annual plant death check, [g d-2]
  real(r8), pointer :: SeasonalNonstElms_pft(:,:)           => null()    !plant stored nonstructural element at current step, [g d-2]
  real(r8), pointer :: SeasonalNonstElmsbeg_pft(:,:)        => null()    !plant stored nonstructural element at prev step,    [g d-2]
  real(r8), pointer :: CanopyLeafSheathC_pft(:)              => null()    !canopy leaf + sheath C,                             [g d-2]
  real(r8), pointer :: AvgCanopyBiomC2Graze_pft(:)          => null()    !section mean shoot C density times local cell area, [g d-2]
  real(r8), pointer :: StandDeadStrutElms_pft(:,:)          => null()    !standing dead element,                              [g d-2]
  real(r8), pointer :: CanopyNonstElms_brch(:,:,:)          => null()    !branch nonstructural element,                       [g d-2]
  real(r8), pointer :: C4PhotoShootNonstC_brch(:,:)         => null()    !branch shoot nonstrucal elelment,                   [g d-2]
  real(r8), pointer :: CanopyNodulNonstElms_brch(:,:,:)     => null()    !branch nodule nonstructural element,                [g d-2]
  real(r8), pointer :: CanopyLeafSheathC_brch(:,:)          => null()    !plant branch leaf + sheath C,                       [g d-2]
  real(r8), pointer :: StalkRsrvElms_brch(:,:,:)            => null()    !branch reserve element mass,                        [g d-2]
  real(r8), pointer :: LeafStrutElms_brch(:,:,:)            => null()    !branch leaf structural element mass,                [g d-2]
  real(r8), pointer :: CanopyNodulStrutElms_brch(:,:,:)     => null()   !branch nodule structural element,                    [g d-2]
  real(r8), pointer :: PetolShethStrutElms_brch(:,:,:)          => null()   !branch sheath structural element,                    [g d-2]
  real(r8), pointer :: EarStrutElms_brch(:,:,:)             => null()   !branch ear structural chemical element mass,         [g d-2]
  real(r8), pointer :: HuskStrutElms_brch(:,:,:)            => null()   !branch husk structural element mass,                 [g d-2]
  real(r8), pointer :: GrainStrutElms_brch(:,:,:)           => null()   !branch grain structural element mass,                [g d-2]
  real(r8), pointer :: StalkStrutElms_brch(:,:,:)           => null()   !branch stalk structural element mass,                [g d-2]
  real(r8), pointer :: ShootElms_brch(:,:,:)                => null()   !branch shoot structural element mass,                [g d-2]
  real(r8), pointer :: SenecStalkStrutElms_brch(:,:,:)      => null()   !branch stalk structural element,                     [g d-2]
  real(r8), pointer :: SapwoodBiomassC_brch(:,:)            => null()   !branch live stalk C,                                 [gC d-2]
  real(r8), pointer :: StalkStrutElms_pft(:,:)              => null()   !canopy stalk structural element mass,                [g d-2]
  real(r8), pointer :: ShootElms_pft(:,:)                   => null()   !current time whole plant shoot element mass,         [g d-2]
  real(r8), pointer :: ShootElmsBeg_pft(:,:)                => null()   !previous whole plant shoot element mass,             [g d-2]
  real(r8), pointer :: CanopySapwoodC_pft(:)                => null()   !canopy active stalk C,                               [g d-2]
  real(r8), pointer :: LeafStrutElms_pft(:,:)               => null()   !canopy leaf structural element mass,                 [g d-2]
  real(r8), pointer :: PetolShethStrutElms_pft(:,:)             => null()   !canopy sheath structural element mass,               [g d-2]
  real(r8), pointer :: StalkRsrvElms_pft(:,:)               => null()   !canopy reserve element mass,                         [g d-2]
  real(r8), pointer :: HuskStrutElms_pft(:,:)               => null()   !canopy husk structural element mass,                 [g d-2]
  real(r8), pointer :: RootBiomCPerPlant_pft(:)             => null()   !root C biomass per plant,                            [g p-1]
  real(r8), pointer :: GrainStrutElms_pft(:,:)              => null()   !canopy grain structural element,                     [g d-2]
  real(r8), pointer :: EarStrutElms_pft(:,:)                => null()   !canopy ear structural element,                       [g d-2]
  real(r8), pointer :: CanopyMassC_pft(:)                   => null()   !Canopy biomass C,                                    [g d-2]
  real(r8), pointer :: ROOTNLim_rpvr(:,:,:)                 => null()   !root N-limitation, 0->1 weaker limitation, [-]
  real(r8), pointer :: ROOTPLim_rpvr(:,:,:)                 => null()   !root P-limitation, 0->1 weaker limitation, [-]
  real(r8), pointer :: LeafC3ChlC_brch(:,:)                 => null()   !Bundle sheath C4/mesophyll C3 chlorophyll C for the branches, [gC chlorophyll d-2]
  real(r8), pointer :: LeafC4ChlC_brch(:,:)                 => null()   !Mesophyll chlorophyll C for the branches, [gC chlorophyll d-2]
  real(r8), pointer :: LeafRubiscoC_brch(:,:)               => null()   !Bundle sheath C4/mesophyll C3 Rubisco C for the branches, [gC Rubisco d-2]
  real(r8), pointer :: LeafPEPC_brch(:,:)                   => null()   !PEP C for the branches, [gC PEP d-2]
  real(r8), pointer :: LeafC3ChlCperm2LA_pft(:)             => null()   !Bundle sheath C4/mesophyll C3 chlorophyll C for the branches, [gC chlorophyll d-2]
  real(r8), pointer :: LeafC4ChlCperm2LA_pft(:)             => null()   !Mesophyll chlorophyll C for the branches, [gC chlorophyll d-2]
  real(r8), pointer :: LeafRubiscoCperm2LA_pft(:)           => null()   !Bundle sheath C4/mesophyll C3 Rubisco C for the branches, [gC Rubisco d-2]
  real(r8), pointer :: LeafPEPCperm2LA_pft(:)               => null()   !PEP C for the branches, [gC PEP d-2]
  real(r8), pointer :: SpecificLeafArea_pft(:)              => null()   !specifc leaf area per g C of leaf mass, [m2 leaf area (gC leaf C)-1]
  real(r8), pointer :: ShootNoduleElms_pft(:,:)             => null()   !current time canopy nodule element mass, [g d-2]
  real(r8), pointer :: ShootNoduleElmsBeg_pft(:,:)          => null()   !previous time canopy nodule element mass, [g d-2]
  real(r8), pointer :: ShootFineNonLeafElms_pft(:,:)        => null()   !canopy non-leaf fine structural element mass, [g d-2]
  real(r8), pointer :: ShootNonstElms_pft(:,:)              => null()   !shoot nonstructural element mass, [g d-2]
  real(r8), pointer :: ShootWoodyElms_pft(:,:)              => null()   !canopy woody element mass, [g d-2]
  real(r8), pointer :: ShootLeafElms_pft(:,:)               => null()   !shoot leaf element mass, [g d-2]
  contains
    procedure, public :: Init => plt_biom_init
    procedure, public :: Destroy => plt_biom_destroy
  end type plant_biom_type


  type(plant_biom_type)     , public, target :: plt_biom      !plant biomass variables

contains

  subroutine plt_biom_init(this)
  ! Allocate biomass pools and initialize derived root aggregates.
  implicit none
  class(plant_biom_type) :: this

  allocate(this%StomatalStress_pft(JP1));  this%StomatalStress_pft=spval
  allocate(this%ZERO4LeafVar_pft(JP1));this%ZERO4LeafVar_pft=spval
  allocate(this%ZERO4Groth_pft(JP1));this%ZERO4Groth_pft=spval
  allocate(this%StandingDeadStrutElms_col(NumPlantChemElms));this%StandingDeadStrutElms_col=spval
  allocate(this%RootNodulStrutElms_rpvr(NumPlantChemElms,JZ1,JP1));this%RootNodulStrutElms_rpvr=spval
  allocate(this%CanopyLeafCLyr_pft(NumCanopyLayers1,JP1));this%CanopyLeafCLyr_pft=spval
  allocate(this%RootNodulNonstElms_rpvr(NumPlantChemElms,JZ1,JP1));this%RootNodulNonstElms_rpvr=spval
  allocate(this%StandDeadCompKElms_pft(NumPlantChemElms,jsken,JP1));this%StandDeadCompKElms_pft=0._r8
  allocate(this%RootMyco2ndStrutElms_rpvr(NumPlantChemElms,jroots,JZ1,MaxNumRootAxes,JP1))
  this%RootMyco2ndStrutElms_rpvr=spval
  allocate(this%RootMyco2ndStrutElms_pvr(NumPlantChemElms,jroots,JZ1,JP1))
  this%RootMyco2ndStrutElms_pvr=0._r8
  allocate(this%RootMycoNonstElms_pft(NumPlantChemElms,jroots,JP1));this%RootMycoNonstElms_pft=spval
  allocate(this%RootMyco1stStrutElms_rpvr(NumPlantChemElms,JZ1,MaxNumRootAxes,JP1));this%RootMyco1stStrutElms_rpvr=0._r8
  allocate(this%Root1stActStructElms_rpvr(NumPlantChemElms,JZ1,MaxNumRootAxes,JP1)); this%Root1stActStructElms_rpvr=0._r8
  allocate(this%Root1stLigStructElms_rpvr(NumPlantChemElms,JZ1,MaxNumRootAxes,JP1)); this%Root1stLigStructElms_rpvr=0._r8
  allocate(this%KLigMax_pft(JP1));this%KLigMax_pft=0.0_r8
  allocate(this%KLigMM_pft(JP1));this%KLigMM_pft=0._r8
  allocate(this%CanopyNonstElmConc_pft(NumPlantChemElms,JP1));this%CanopyNonstElmConc_pft=spval
  allocate(this%CanopyNonstOMConc_pft(JP1));this%CanopyNonstOMConc_pft=0._r8
  allocate(this%CanopyNonstElms_pft(NumPlantChemElms,JP1));this%CanopyNonstElms_pft=spval
  allocate(this%CanopyNodulNonstElms_pft(NumPlantChemElms,JP1));this%CanopyNodulNonstElms_pft=spval
  allocate(this%CanopyNoduleNonstCConc_pft(JP1));this%CanopyNoduleNonstCConc_pft=spval
  allocate(this%RootProteinConc_rpvr(jroots,JZ1,JP1));this%RootProteinConc_rpvr=spval
  allocate(this%RootProteinC_pvr(jroots,JZ1,JP1));this%RootProteinC_pvr=0._r8
  allocate(this%RootMediumStructElms_rpvr(NumPlantChemElms,JZ1,MaxNumRootAxes,JP1)); this%RootMediumStructElms_rpvr=spval
  allocate(this%RootMycoActiveBiomC_pvr(jroots,JZ1,JP1));this%RootMycoActiveBiomC_pvr=spval
  allocate(this%RootMycoMassElm_pvr(NumPlantChemElms,jroots,JZ1,JP1)); this%RootMycoMassElm_pvr = 0._r8
  allocate(this%PopuRootMycoC_pvr(jroots,JZ1,JP1));this%PopuRootMycoC_pvr=spval
  allocate(this%RootMycoNonstElms_rpvr(NumPlantChemElms,jroots,JZ1,JP1));this%RootMycoNonstElms_rpvr=spval
  allocate(this%RootNonstructElmConc_pvr(NumPlantChemElms,jroots,JZ1,JP1));this%RootNonstructElmConc_pvr=spval
  allocate(this%RootNonstOMConc_pvr(jroots,JZ1,JP1));this%RootNonstOMConc_pvr=0._r8
  allocate(this%ROOTNLim_rpvr(jroots,JZ1,JP1)); this%ROOTNLim_rpvr=0._r8
  allocate(this%ROOTPLim_rpvr(jroots,JZ1,JP1)); this%ROOTPLim_rpvr=0._r8
  allocate(this%CanopyNonstElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%CanopyNonstElms_brch=spval
  allocate(this%C4PhotoShootNonstC_brch(MaxNumBranches,JP1));this%C4PhotoShootNonstC_brch=spval
  allocate(this%CanopySapwoodC_pft(JP1));this%CanopySapwoodC_pft=spval
  allocate(this%ShootElms_pft(NumPlantChemElms,JP1));this%ShootElms_pft=spval
  allocate(this%ShootElmsBeg_pft(NumPlantChemElms,JP1));this%ShootElmsBeg_pft=spval
  allocate(this%LeafProteinC_node(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%LeafProteinC_node=spval
  allocate(this%PetoleProteinC_node(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%PetoleProteinC_node=spval
  allocate(this%StructInternodeElms_brch(NumPlantChemElms,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%StructInternodeElms_brch=spval
  allocate(this%LeafElmntNode_brch(NumPlantChemElms,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%LeafElmntNode_brch=spval
  allocate(this%PetolShethElmntNode_brch(NumPlantChemElms,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%PetolShethElmntNode_brch=spval
  allocate(this%LeafLayerElms_node(NumPlantChemElms,NumCanopyLayers1,0:MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%LeafLayerElms_node=spval
  allocate(this%LeafProteinC_brch(MaxNumBranches,JP1)); this%LeafProteinC_brch=spval
  allocate(this%LeafProteinCperm2LA_pft(JP1)); this%LeafProteinCperm2LA_pft=spval
  allocate(this%SapwoodBiomassC_brch(MaxNumBranches,JP1));this%SapwoodBiomassC_brch=spval
  allocate(this%CanopyNodulNonstElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%CanopyNodulNonstElms_brch=spval
  allocate(this%LeafPetoNonstElmConc_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LeafPetoNonstElmConc_brch=spval
  allocate(this%RootStrutElms_pft(NumPlantChemElms,JP1));this%RootStrutElms_pft=spval
  allocate(this%RootMedStruct_pvr(NumPlantChemElms,JZ1,JP1));this%RootMedStruct_pvr=0._r8
  allocate(this%Root1stActStruct_pvr(NumPlantChemElms,JZ1,JP1)); this%Root1stActStruct_pvr=0._r8
  allocate(this%Root1stLigStruct_pvr(NumPlantChemElms,JZ1,JP1)); this%Root1stLigStruct_pvr=0._r8
  allocate(this%StandDeadStrutElmsBeg_pft(NumPlantChemElms,JP1));this%StandDeadStrutElmsBeg_pft=spval
  allocate(this%tCanLeafC_clyr(NumCanopyLayers1));this%tCanLeafC_clyr=spval
  allocate(this%RootElms_pft(NumPlantChemElms,JP1));this%RootElms_pft=spval
  allocate(this%RootActiveC_pft(JP1));this%RootActiveC_pft=0._r8
  allocate(this%RootNoduleElms_pft(NumPlantChemElms,JP1));this%RootNoduleElms_pft=spval
  allocate(this%RootNoduleElmsBeg_pft(NumPlantChemElms,JP1));this%RootNoduleElmsBeg_pft=spval
  allocate(this%RootElmsBeg_pft(NumPlantChemElms,JP1));this%RootElmsBeg_pft=spval
  allocate(this%SeedPlantedElm_pft(NumPlantChemElms,JP1));this%SeedPlantedElm_pft=spval
  allocate(this%SeasonalNonstElms_pft(NumPlantChemElms,JP1));this%SeasonalNonstElms_pft=spval
  allocate(this%SeasonalNonstCDayAve_pft(JP1)); this%SeasonalNonstCDayAve_pft=spval
  allocate(this%SeasonalNonstElmsbeg_pft(NumPlantChemElms,JP1));this%SeasonalNonstElmsbeg_pft=spval
  allocate(this%TotBegVegE_pft(NumPlantChemElms,JP1)); this%TotBegVegE_pft=spval
  allocate(this%TotEndVegE_pft(NumPlantChemElms,JP1)); this%TotEndVegE_pft=spval
  allocate(this%CanopyLeafSheathC_pft(JP1));this%CanopyLeafSheathC_pft=spval
  allocate(this%StandDeadStrutElms_pft(NumPlantChemElms,JP1));this%StandDeadStrutElms_pft=spval

  allocate(this%PetolShethStrutElms_pft(NumPlantChemElms,JP1));this%PetolShethStrutElms_pft=spval
  allocate(this%StalkStrutElms_pft(NumPlantChemElms,JP1));this%StalkStrutElms_pft=spval
  allocate(this%StalkRsrvElms_pft(NumPlantChemElms,JP1));this%StalkRsrvElms_pft=spval
  allocate(this%GrainStrutElms_pft(NumPlantChemElms,JP1));this%GrainStrutElms_pft=0._r8
  allocate(this%HuskStrutElms_pft(NumPlantChemElms,JP1));this%HuskStrutElms_pft=spval
  allocate(this%RootBiomCPerPlant_pft(JP1));this%RootBiomCPerPlant_pft=spval
  allocate(this%EarStrutElms_pft(NumPlantChemElms,JP1));this%EarStrutElms_pft=spval
  allocate(this%CanopyMassC_pft(JP1)); this%CanopyMassC_pft=0._r8
  allocate(this%ShootElms_pft(NumPlantChemElms,JP1));this%ShootElms_pft=spval
  allocate(this%AvgCanopyBiomC2Graze_pft(JP1));this%AvgCanopyBiomC2Graze_pft=spval
  allocate(this%LeafStrutElms_pft(NumPlantChemElms,JP1));this%LeafStrutElms_pft=spval
  allocate(this%StandingDeadInitC_pft(JP1));this%StandingDeadInitC_pft=spval
  allocate(this%CanopyLeafSheathC_brch(MaxNumBranches,JP1));this%CanopyLeafSheathC_brch=spval
  allocate(this%StalkRsrvElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%StalkRsrvElms_brch=spval
  allocate(this%LeafStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%LeafStrutElms_brch=spval
  allocate(this%CanopyNodulStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%CanopyNodulStrutElms_brch=spval
  allocate(this%PetolShethStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%PetolShethStrutElms_brch=spval
  allocate(this%EarStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%EarStrutElms_brch=spval
  allocate(this%HuskStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%HuskStrutElms_brch=spval
  allocate(this%GrainStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%GrainStrutElms_brch=0._r8
  allocate(this%StalkStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%StalkStrutElms_brch=spval
  allocate(this%ShootElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%ShootElms_brch=spval
  allocate(this%SenecStalkStrutElms_brch(NumPlantChemElms,MaxNumBranches,JP1));this%SenecStalkStrutElms_brch=spval
  allocate(this%RootMyco1stElm_raxs(NumPlantChemElms,MaxNumRootAxes,JP1));this%RootMyco1stElm_raxs=spval

  allocate(this%LeafC3ChlC_brch(MaxNumBranches,JP1));this%LeafC3ChlC_brch=0._r8
  allocate(this%LeafC4ChlC_brch(MaxNumBranches,JP1));this%LeafC4ChlC_brch=0._r8
  allocate(this%LeafRubiscoC_brch(MaxNumBranches,JP1));this%LeafRubiscoC_brch=0._r8
  allocate(this%LeafPEPC_brch(MaxNumBranches,JP1));this%LeafPEPC_brch=0._r8
  allocate(this%LeafC3ChlCperm2LA_pft(JP1))     ;this%LeafC3ChlCperm2LA_pft=0._r8
  allocate(this%LeafC4ChlCperm2LA_pft(JP1))     ;this%LeafC4ChlCperm2LA_pft=0._r8
  allocate(this%LeafRubiscoCperm2LA_pft(JP1))   ;this%LeafRubiscoCperm2LA_pft=0._r8
  allocate(this%LeafPEPCperm2LA_pft(JP1))       ;this%LeafPEPCperm2LA_pft=0._r8
  allocate(this%SpecificLeafArea_pft(JP1)) ; this%SpecificLeafArea_pft=0._r8
  allocate(this%ShootNoduleElms_pft(NumPlantChemElms,JP1)); this%ShootNoduleElms_pft=0._r8
  allocate(this%ShootNoduleElmsBeg_pft(NumPlantChemElms,JP1));this%ShootNoduleElmsBeg_pft=0._r8

  allocate(this%ShootFineNonLeafElms_pft(NumPlantChemElms,JP1));this%ShootFineNonLeafElms_pft=0._r8
  allocate(this%ShootNonstElms_pft(NumPlantChemElms,JP1));       this%ShootNonstElms_pft=0._r8
  allocate(this%ShootWoodyElms_pft(NumPlantChemElms,JP1));       this%ShootWoodyElms_pft=0._r8
  allocate(this%ShootLeafElms_pft(NumPlantChemElms,JP1)); this%ShootLeafElms_pft=0._r8
  end subroutine plt_biom_init

  subroutine plt_biom_destroy(this)
  ! Release biomass API storage handled by this lifecycle procedure.

  implicit none
  class(plant_biom_type) :: this

  if(associated(this%RootMyco2ndStrutElms_pvr))deallocate(this%RootMyco2ndStrutElms_pvr)
  if(associated(this%CanopyNonstOMConc_pft))deallocate(this%CanopyNonstOMConc_pft)
  if(associated(this%RootNonstOMConc_pvr))deallocate(this%RootNonstOMConc_pvr)
  if(associated(this%RootActiveC_pft))deallocate(this%RootActiveC_pft)

!  if(allocated(ZERO4LeafVar_pft))deallocate(ZERO4LeafVar_pft)
!  if(allocated(ZERO4Groth_pft))deallocate(ZERO4Groth_pft)
!  if(allocated(RootNodulNonstElms_rpvr))deallocate(RootNodulNonstElms_rpvr)
!  if(allocated(RootNodulStrutElms_rpvr))deallocate(RootNodulStrutElms_rpvr)
!  if(allocated(CanopyLeafCLyr_pft))deallocate(CanopyLeafCLyr_pft)
!  call destroy(RootMyco1stElm_raxs)
!  if(allocated(RootMyco1stStrutElms_rpvr))deallocate(RootMyco1stStrutElms_rpvr)
!  if(allocated(StandDeadCompKElms_pft))deallocate(StandDeadCompKElms_pft)
!  if(allocated(RootMyco2ndStrutElms_rpvr))deallocate(RootMyco2ndStrutElms_rpvr)
!  if(allocated(CanopyNonstElmConc_pft))deallocate(CanopyNonstElmConc_pft)
!  if(allocated(CanopyNonstElms_pft))deallocate(CanopyNonstElms_pft)
!  if(allocated(CanopyNodulNonstElms_pft))deallocate(CanopyNodulNonstElms_pft)
!  if(allocated(CanopyNoduleNonstCConc_pft))deallocate(CanopyNoduleNonstCConc_pft)
!  if(allocated(RootProteinConc_rpvr))deallocate(RootProteinConc_rpvr)
!  if(allocated(RootProteinC_pvr))deallocate(RootProteinC_pvr)
!  if(allocated(RootMycoActiveBiomC_pvr))deallocate(RootMycoActiveBiomC_pvr)
!  if(allocated( PopuRootMycoC_pvr))deallocate( PopuRootMycoC_pvr)
!  if(allocated(RootMycoNonstElms_rpvr))deallocate(RootMycoNonstElms_rpvr)
!  if(allocated(WVSTK))deallocate(WVSTK)
!  if(allocated(LeafElmntNode_brch))deallocate(LeafElmntNode_brch)
!  if(allocated(LeafProteinC_node))deallocate(LeafProteinC_node)
!  if(allocated(PetoleProteinC_node))deallocate(PetoleProteinC_node)
!  if(allocated(LeafLayerElms_node))deallocate(LeafLayerElms_node)
!  if(allocated(StructInternodeElms_brch))deallocate(StructInternodeElms_brch)
!  if(allocated(PetolShethElmntNode_brch))deallocate(PetolShethElmntNode_brch)
!  if(allocated(SapwoodBiomassC_brch))deallocate(SapwoodBiomassC_brch)
!  if(allocated(PPOOL))deallocate(PPOOL)
!  if(allocated(CanopyNodulNonstElms_brch))deallocate(CanopyNodulNonstElms_brch)
!  if(allocated(LeafPetoNonstElmConc_brch))deallocate(LeafPetoNonstElmConc_brch)
!  if(allocated(RootStrutElms_pft))deallocate(RootStrutElms_pft)
!  if(allocated(tCanLeafC_clyr))deallocate(tCanLeafC_clyr)
!  if(allocated(CanopyNonstElms_brch))deallocate(CanopyNonstElms_brch)
!  if(allocated(RootBiomCPerPlant_pft))deallocate(RootBiomCPerPlant_pft)
!  if(allocated(WTLF))deallocate(WTLF)
!  if(allocated(WTSHE))deallocate(WTSHE)
!  if(allocated(WTRSV))deallocate(WTRSV)
!  if(allocated(WTSTK))deallocate(WTSTK)
!  if(allocated(WTLFP))deallocate(WTLFP)
!  if(allocated(WTGR))deallocate(WTGR)
!  if(allocated(WTLFN))deallocate(WTLFN)
!  if(allocated(WTEAR))deallocate(WTEAR)
!  if(allocated(HuskStrutElms_pft))deallocate(HuskStrutElms_pft)
!  if(allocated(LeafChemElmRemob_brch))deallocate(LeafChemElmRemob_brch)
!  if(allocated(PetolShethChemElmRemob_brch))deallocate(PetolShethChemElmRemob_brch)
!  if(allocated(CanopyLeafSheathC_brch))deallocate(CanopyLeafSheathC_brch)
!  if(allocated(StalkRsrvElms_brch))deallocate(StalkRsrvElms_brch)
!  if(allocated(LeafStrutElms_brch))deallocate(LeafStrutElms_brch)
!  if(allocated(CanopyNodulStrutElms_brch))deallocate(CanopyNodulStrutElms_brch)
!  if(allocated(PetolShethStrutElms_brch))deallocate(PetolShethStrutElms_brch)
!  if(allocated(EarStrutElms_brch))deallocate(EarStrutElms_brch)
!  if(allocated(HuskStrutElms_brch))deallocate(HuskStrutElms_brch)
!  if(allocated(GrainStrutElms_brch))deallocate(GrainStrutElms_brch)
!  if(allocated(StalkStrutElms_brch))deallocate(StalkStrutElms_brch)
!  if(allocated(ShootElms_brch))deallocate(ShootElms_brch)
!  if(allocated(SenecStalkStrutElms_brch))deallocate(SenecStalkStrutElms_brch)
!  if(allocated(WTSTDI))deallocate(WTSTDI)
!  if(allocated(SeedPlantedElm_pft))deallocate(SeedPlantedElm_pft)
!  if(allocated(WTLS))deallocate(WTLS)
!  if(allocated(ShootElms_pft))deallocate(ShootElms_pft)
!  if(allocated(AvgCanopyBiomC2Graze_pft))deallocate(AvgCanopyBiomC2Graze_pft)
!  if(allocated(StandDeadStrutElms_pft)deallocate(StandDeadStrutElms_pft)
!  if(allocated(NodulStrutElms_pft))deallocate(NodulStrutElms_pft)
  end subroutine plt_biom_destroy
end module PlantBiomassAPIData
