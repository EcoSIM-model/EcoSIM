module PlantMorphologyAPIData
  ! Owns the morphology API type, instance and lifecycle procedures.
  use PlantAPICommonData
  implicit none
  private
  public :: plt_morph_init, plt_morph_destroy

  type, public :: plant_morph_type
  real(r8) :: LeafStalkAreaAll_col                       !stalk area of combined, each PFT canopy,[m^2 d-2]
  real(r8) :: CanopyLeafArea_col                      !grid canopy leaf area, [m2 d-2]
  real(r8) :: StemArea_col                            !grid canopy stem area, [m2 d-2]
  real(r8) :: StandDeadSurfArea_col                   !grid canopy standing-dead surface area, [m2 d-2]
  real(r8) :: CanopyHeight_col                        !canopy height , [m]
  real(r8), pointer :: StalkAxialResist_pft(:)         => null() !stalk axial resistance per m for water transport, [MPa h m-4]
  real(r8), pointer :: RootSingleVesselArea_pft(:)               => null() !
  REAL(R8), POINTER :: enh_cyto_pft(:)                 => null() !cytokinin sensitivity of corase root thickening, [-]
  real(r8), pointer :: tlai_day_pft(:)                 => null() !prescribed leaf area, [m2 cell-1]
  real(r8), pointer :: tsai_day_pft(:)                 => null() !prescribed stem area, [m2 cell-1]
  real(r8), pointer :: PARTS_brch(:,:,:)               => null() !fraction of C allocated to each morph unit,                                 [-]
  real(r8), pointer :: FineRootVolPerMassC_pft(:,:)    => null() !Fine root volume:mass ratio,                                                     [m3 g-1]
  real(r8), pointer :: CRootActVolPerMassC_pft(:)    => null()   !Coarse root active zone volume:mass ratio,  [m3 g-1]
  real(r8), pointer :: CRootLigVolPerMassC_pft(:)     => null()  !Coarse root inactive zone volume:mass ratio,  [m3 g-1]
  real(r8), pointer :: RootPorosity_pft(:,:)           => null() !root porosity,                                                              [m3 m-3]
  real(r8), pointer :: Root2ndXSecArea_pft(:,:)        => null() !root  cross-sectional area  secondary axes,                                 [m2]
  real(r8), pointer :: Root1stXSecArea_pft(:,:)        => null() !root cross-sectional area primary axes,                                     [m2]
  real(r8), pointer :: Root1stMaxRadius1_pft(:,:)      => null() !root diameter primary axes,                                                 [m]
  real(r8), pointer :: Root2ndMaxRadius1_pft(:,:)      => null() !root diameter secondary axes,                                               [m]
  real(r8), pointer :: SeedCMass_pft(:)                => null() !grain size at seeding,                                                      [gC/seed]
  real(r8), pointer :: SeedWidth2LenRatio_pft(:)       => null() !Seed width to length ratio, assuming prolate spheroid
  real(r8), pointer :: RootPoreTortu4Gas_pft(:,:)      => null() !power function of root porosity used to calculate root gaseous diffusivity, [-]
  logical,  pointer :: flag2ndGrowth_pvr(:,:,:)        => null() !flag for secondary growth of primary roots, [-]
  real(r8), pointer :: Root1stLenPP_rpvr(:,:,:)        => null() !primary root axis length in soil layer,                                             [m d-2]
  real(r8), pointer :: Root2ndLen_rpvr(:,:,:,:)        => null() !root layer length secondary axes,                                           [m d-2]
  real(r8), pointer :: RootAge_rpvr(:,:,:)             => null() !root age, [h]
  real(r8), pointer :: RootTotLenPerPlant_pvr(:,:,:)   => null() !total root length per plant,                                                [m p-1]
  real(r8), pointer :: RootAbsorbLenPerPlant_pvr(:,:,:)=> null() !total absorptive root length per plant in layer, [m p-1]
  real(r8), pointer :: RootLenPerPlant_pvr(:,:,:)      => null() !fine root length per plant, [m p-1]
  real(r8), pointer :: Root2ndEffLen4uptk_rpvr(:,:,:)  => null() !Layer effective root length four resource uptake, [m]
  real(r8), pointer :: DistRootEffDepz_pvr(:,:)        => null() !Effective shoot-root transport depth, [m]
  real(r8), pointer :: RootSinkScalar_pvr(:,:,:)       => null() !root sink scalar to account for the missing of intermediate size roots, [0-1]
  real(r8), pointer :: Root1stSpecLen_pft(:,:)         => null() !specific root length primary axes,                                          [m g-1]
  real(r8), pointer :: Root2ndSpecLen_pft(:,:)         => null() !specific root length secondary axes,                                        [m g-1]
  real(r8), pointer :: Root2ndXNum_rpvr(:,:,:,:)       => null() !root layer number secondary axes,                                           [d-2]
  real(r8), pointer :: RootMediumXNum_rpvr(:,:,:)      => null() !living population medium-root axes by layer and structural-root group, [# d-2]
  real(r8), pointer :: RootFineFrac2Med_rpvr(:,:,:)       => null() !fine-axis-count-weighted fraction attached to medium roots, by root category, [-]
  real(r8), pointer :: CRootLumenArea_rpvr(:,:,:)      => null() !coarse roots lumen area for root axes, [m2]
  real(r8), pointer :: CRootLumenArea_pvr(:,:)         => null() !coarse roots lumen area, [m2]
  real(r8), pointer :: MRootLumenArea_pvr(:,:)         => null() !medium roots lumen area, [m2]
  real(r8), pointer :: MRootLumenArea_rpvr(:,:,:)         => null() !population medium-root lumen area per structural axis, [m2]
  real(r8), pointer :: Root1stDepz_raxes(:,:)           => null() !root layer depth,                                                           [m]
  real(r8), pointer :: Root1stAxesTipDepz2Surf_pft(:,:)   => null() !plant primary depth relative to column surface, [m]
  real(r8), pointer :: Root1stLenLoc_rpvr(:,:,:)       => null() !local structrual root length in layer, [m]
  real(r8), pointer :: ClumpFactorInit_pft(:)          => null() !initial clumping factor for self-shading in canopy layer,                   [-]
  real(r8), pointer :: ClumpFactorNow_pft(:)           => null() !clumping factor for self-shading in canopy layer at current LAI,            [-]
  real(r8), pointer :: FineRootBranchFreq_pft(:)           => null() !Fine root brancing frequency,                                                    [m-1]
  real(r8), pointer :: MediumRootBranchFreq_pft(:)     => null() !Medium root branch frequency, [m-1]
  real(r8), pointer :: HypocotHeight_pft(:)            => null() !cotyledon height,                                                           [m]
  real(r8), pointer :: CanopyHeight4WatUptake_pft(:)   => null() !canopy height,                                                              [m]
  real(r8), pointer :: CanopyHeightZ_col(:)            => null() !canopy layer height,                                                        [m]
  real(r8), pointer :: CanopyStemSurfArea_pft(:)           => null() !plant stem area,                                                            [m2 d-2]
  real(r8), pointer :: CanopyLeafArea_pft(:)           => null() !plant canopy leaf area,                                                     [m2 d-2]
  real(r8), pointer :: ShootNodeNum_brch(:,:)          => null() !shoot node number,                                                          [-]
  real(r8), pointer :: ShootNodeNumAtInitFloral_brch(:,:)    => null() !shoot node number at floral initiation,                                     [-]
  real(r8), pointer :: ShootNodeNumAtAnthesis_brch(:,:)  => null() !shoot node number at anthesis,                                              [-]
  real(r8), pointer :: SineBranchAngle_pft(:)          => null() !branching angle,                                                            [degree from horizontal]
  real(r8), pointer :: LeafAreaZsec_brch(:,:,:,:,:)    => null() !leaf surface area,                                                          [m2 d-2]
  real(r8), pointer :: PotentialSeedSites_brch(:,:)    => null() !branch potential grain number,                                              [d-2]
  real(r8), pointer :: SetNumberSeeds_brch(:,:)            => null() !branch grain number,                                                        [d-2]
  real(r8), pointer :: CanPBranchHeight(:,:)           => null() !branch height,                                                              [m]
  real(r8), pointer :: LeafAreaDying_brch(:,:)         => null() !branch leaf area,                                                           [m2 d-2]
  real(r8), pointer :: LeafAreaLive_brch(:,:)          => null() !branch leaf area,                                                           [m2 d-2]
  real(r8), pointer :: SinePetolShethAngle_pft(:)         => null() !sheath angle,                                                               [degree from horizontal]
  real(r8), pointer :: LeafAngleClass_pft(:,:)         => null() !fractionction of leaves in different angle classes,                         [-]
  real(r8), pointer :: CanopyStemSurfAreaZ_pft(:,:)        => null() !plant canopy layer stem area,                                               [m2 d-2]
  real(r8), pointer :: CanopyLeafAreaZ_pft(:,:)        => null() !canopy layer leaf area,                                                     [m2 d-2]
  real(r8), pointer :: LeafArea_node(:,:,:)            => null() !leaf area,                                                                  [m2 d-2]
  real(r8), pointer :: CanopySeedNum_pft(:)            => null() !canopy grain number,                                                        [d-2]
  real(r8), pointer :: CanopySeedNumX_pft(:)           => null() !last nonzero canopy grain number,                                           [d-2]
  real(r8), pointer :: SeedDepth_pft(:)                => null() !seeding depth,                                                              [m]
  real(r8), pointer :: PlantinDepz_pft(:)              => null() !planting depth,                                                             [m]
  real(r8), pointer :: SeedMeanLen_pft(:)              => null() !seed length,                                                                [m]
  real(r8), pointer :: SeedVolumeMean_pft(:)           => null() !seed volume,                                                                [m3 ]
  real(r8), pointer :: SeedAreaMean_pft(:)             => null() !seed surface area,                                                          [m2]
  real(r8), pointer :: CanopyStemAareZ_col(:)          => null() !total stem area,                                                            [m2 d-2]
  real(r8), pointer :: PetolShethLen2Mass_pft(:)             => null() !PetolSheth length:mass during growth,                                          [m gC-1]
  real(r8), pointer :: NodeLenPergC_pft(:)             => null() !internode length:mass during growth,                                        [m gC-1]
  real(r8), pointer :: SLA1_pft(:)                     => null() !leaf area:mass during growth,                                               [m2 gC-1]
  real(r8), pointer :: CanopyLeafAareZ_col(:)          => null() !total leaf area,                                                            [m2 d-2]
  real(r8), pointer :: LeafStalkAreaAct_pft(:)         => null() !radiation-active plant leaf+stem/stalk area,                                                 [m2 d-2]
  real(r8), pointer :: StandDeadSurfAreaAct_pft(:)     => null() !radiation-active standing dead surface area, [m2 d-2]
  real(r8), pointer :: StalkNodeVertLength_brch(:,:,:) => null() !internode height,                                                           [m]
  real(r8), pointer :: PetoleLength_node(:,:,:)      => null() !sheath height,                                                              [m]
  real(r8), pointer :: StalkNodeHeight_brch(:,:,:)  => null() !internode height,                                                           [m]
  real(r8), pointer :: StemAreaZsec_brch(:,:,:,:)      => null() !stem surface area,                                                          [m2 d-2]
  real(r8), pointer :: CanopyLeafArea_lnode(:,:,:,:)   => null() !layer/node/branch leaf area,                                                [m2 d-2]
  real(r8), pointer :: CanopyStalkSurfArea_lbrch(:,:,:)    => null() !plant canopy layer branch stem area,                                        [m2 d-2]
  real(r8), pointer :: CanopySurfAreaProfDead_pft(:,:) => null() !standing dead plant canopy surface area, [m2 d-2]
  real(r8), pointer :: StandDeadSurfArea_pft(:)        => null() !standing dead canopy surface area profile, [m2 d-2]
  real(r8), pointer :: ClumpFactor_pft(:)              => null() !clumping factor for self-shading in canopy layer,                           [-]
  real(r8), pointer :: ShootNodeNumAtPlanting_pft(:)   => null() !number of nodes in seed,                                                    [-]
  real(r8), pointer :: CanopyHeightLive_pft(:)             => null() !live plant canopy height,                                                              [m]
  real(r8), pointer :: CanopyHeightDead_pft(:)        => null() !canopy height for standing dead, [m]
  real(r8), pointer :: StalkHeight_pft(:)              => null() !stalk height/length, [m]
  real(r8), pointer :: StemSpecVolume_pft(:)               => null()  !stalk specific volume, [m3 gC-1]
  real(r8), pointer :: CanopyLeafAreaMAX_pft(:)       => null() !running maximum leaf area, [m2 d-2]
  logical, pointer  :: lreset_laimax_pft(:)   => null() !toggle to reset laimax for woody vascular plants, [-]
  real(r8), pointer :: StalkAveRadius_pft(:)        => null() !main stalk radius,[m]
  integer , pointer :: iPlantGrainType_pft(:)        => null() !grain type (below or above-ground),[-]
  integer,  pointer :: iPlantNfixType_pft(:)          => null() !N2 fixation type,[-]
  integer,  pointer :: Myco_pft(:)                    => null() !mycorrhizal type (no or yes),[-]
  integer,  pointer :: MainBranchNum_pft(:)           => null() !number of main branch,[-]
  integer,  pointer :: MaxSoilLays4Root_pft(:)            => null() !maximum soil layer number for all root axes,[-]
  integer,  pointer :: NMaxRootBotLayer_pft(:)         => null() !maximum soil layer number for all root axes, [-]
  integer,  pointer :: NumStructuralRootAxes_pft(:)             => null() !number of structural root axes,[-]
  real(r8), pointer :: RootMatureAge_pft(:)           => null() !Root age to trigger secondary growth, [h]
  integer,  pointer :: NumCogrowthNode_pft(:)         => null() !number of concurrently growing nodes,[-]
  integer,  pointer :: BranchNumber_pft(:)            => null() !main branch numeric id,[-]
  integer,  pointer :: NumOfBranches_pft(:)           => null() !number of branches,[-]
  integer,  pointer :: NRoot1stTipLay_raxes(:,:)      => null() !maximum soil layer number for root axes,     [-]
  integer,  pointer :: KMinNumLeaf4GroAlloc_brch(:,:) => null() !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION,[-]
  integer,  pointer :: BranchNumerID_brch(:,:)         => null() !branch meric id,                             [-]
  integer,  pointer :: NGTopRootLayer_pft(:)          => null() !soil layer at planting depth,                [-]
  real(r8), pointer :: RootSinkWeight_pvr(:,:)        => null() !Root nonst element sink profile, [d-2]
  real(r8), pointer :: Root2ndSinkWeight_pvr(:,:,:)   => null() !Secondary root nonst element sink profile, [d-2]
  real(r8), pointer :: Root1stSinkWeight_pvr(:,:)     => null() !primary root nonst element sink profile, [d-2]
  real(r8), pointer :: RootMSinkWeight_pvr(:,:)       => null() !medium size roots nonst element sink profile, [d-2]
  real(r8), pointer :: Root1stTipSinkWeight_pft(:)     => null() !primary root tip nonst element sink, [d-2]
  real(r8), pointer :: Root1stTransptArea_pvr(:,:,:)        => null()    !root cross section area for water/gas transport,    [g d-2]
  real(r8), pointer :: RootMedTransptArea_pvr(:,:,:)  => null()  !transport area by medisum size roots, [-]
  integer,  pointer :: KLeafNumber_brch(:,:)          => null() !leaf number,                                 [-]
  real(r8), pointer :: RootSegAges_raxes(:,:,:)        => null()   !age of different active root segments, [h]
  integer , pointer :: NActiveRootSegs_raxes(:,:)       => null()   !number of active root segments, [-]
  real(r8), pointer :: RootSegBaseDepth_raxes(:,:)      => null()   !base depth of different root axes, [m]
  real(r8), pointer :: RootSeglengths_raxes(:,:,:)      => null()   !root length in each segment, [m]
  integer , pointer :: IndRootSegBase_raxes(:,:)        => null()   !Index for base segment under tracking, [-]
  integer , pointer :: IndRootSegTip_raxes(:,:)       => null()   !Index for tip segment under tracking, [-]
  real(r8), pointer :: NumOfLeaves_brch(:,:)          => null() !leaf number,                          [-]
  real(r8), pointer :: GrothStalkMaxSeedSites_pft(:)  => null() !maximum grain node number per branch, [-]
  real(r8), pointer :: MaxSeedNumPerSite_pft(:)       => null() !maximum grain number per node,        [-]
  real(r8), pointer :: rLen2WidthLeaf_pft(:)          => null() !leaf length:width ratio,              [-]
  real(r8), pointer :: SeedCMassMax_pft(:)            => null() !maximum grain size,                   [g]
  real(r8), pointer :: Root1stRadius_pvr(:,:,:)       => null() !root layer diameter primary axes,     [m]
  real(r8), pointer :: Root1stRadius_rpvr(:,:,:)      => null() !root layer diameter for each primary axes,  [m]
  real(r8), pointer :: RootCRRadius0_rpvr(:,:,:)      => null() !initial radius of roots that may undergo secondary growth, [m]
  real(r8), pointer :: Root2ndRadius_rpvr(:,:,:)      => null() !root layer diameter secondary axes,   [m]
  real(r8), pointer :: fRootTube_rpvr(:,:,:)          => null() !fraction of root for transport,[-]
  real(r8), pointer :: RootRaidus_rpft(:,:)           => null() !root internal radius,                 [m]
  real(r8), pointer :: Root1stMaxRadius_pft(:,:)      => null() !maximum radius of primary roots,      [m]
  real(r8), pointer :: Root2ndMaxRadius_pft(:,:)      => null() !maximum radius of secondary roots,    [m]
  real(r8), pointer :: RootSingleVesselRstaxial_pft(:)       => null() !axial resistance for a single 1 m water transport vessel [MPa h m-4]
  real(r8), pointer :: RootVesselRadius_pft(:)        => null() !typical radius of the water transport vessel in primary roots, [m]
  real(r8), pointer :: RootRadialResist_pft(:,:)      => null() !root radial resistivity,              [MPa h m-2]
  real(r8), pointer :: Root2ndAxialResist_pft(:,:)       => null() !root axial resistivity,               [MPa h m-4]
  real(r8), pointer :: totRootLenDens_vr(:)           => null() !total root length density,            [m m-3]
  real(r8), pointer :: Root1stXNumL_pvr(:,:)          => null() !root layer number primary axes,       [d-2]
  real(r8), pointer :: fctyok_scalar_rpvr(:,:,:)      => null() !cytokinin scalar for corase root sink, [-]
  REAL(R8), POINTER :: Num1stAxesPerStructRootX_pft(:)      => null() !primary root axes number per structural root axis, [d-2]
  REAL(R8), POINTER :: NumMediumRootAxes_rpvr(:,:,:)         => null() ! Number of medium size root axes in layer for structrual axes, [d-2]
  REAL(R8), POINTER :: RootMediumXNum_pvr(:,:)          => null() ! Number of medium size root axes in layer, [d-2]
  real(r8), pointer :: RootMediumLength_rpvr(:,:,:)        => null()      !root layer length for medium size axes, [m d-2]
  real(r8), pointer :: RootMediumMeanLength_rpvr(:,:,:) => null() !derived mean medium-root length per layer/structural axis, [m]; internal only
  real(r8), pointer :: RootMediumLength_pvr(:,:)      => null() !root layer mean length of individual medium roots, [m]
  real(r8), pointer :: RootMediumRadius_rpvr(:,:,:)    => null() !root layer radius for medium size axes, [m]
  real(r8), pointer :: Num1stAxesPerStructRootXPOP_pft(:)   =>null() !population primary root axes number on one structrual axis, [d-2]
  real(r8), pointer :: Num1stRootAxesPP_pft(:)   => null() !number of primary root axesr per plant, [d-2]
  real(r8), pointer :: Radius95pctMature_pft(:)       => null() !Critical radius where the woody radius is considered 95% mature, [m]
  real(r8), pointer :: Root2ndXNumL_rpvr(:,:,:)       => null() !root layer number axes,               [d-2]
  real(r8), pointer :: Root2ndVH2O_rpvr(:,:,:,:)      => null()  !water-occupied 2nd root volume, [m3 m-3]
  real(r8), pointer :: RootMediumVH2O_rpvr(:,:,:)     => null()  !population medium-root lumen volume per structural axis, [m3 H2O]
  real(r8), pointer :: Root1stVH2O_rpvr(:,:,:)        => null()  !water-occupied xylem volume in corase roots, [m3 m-3]
  real(r8), pointer :: xylemPhi_min_pft(:)            => null()  !the fraction found in the youngest xylem that as lumen for tree, [m2/m2]
  real(r8), pointer :: xylemPhi_max_pft(:)            => null()  !asymptotic limit fraction of the xyxlem area as lumen for tree, [m2/m2]
  real(r8), pointer :: xylemPhi_mean_pft(:)           => null()  !the mean fraction found in the seminal root as lumen for non-tree roots, [m2/m2]

  real(r8), pointer :: RootLenDensPerPlant_pvr(:,:,:) => null() !layer root length density,            [m m-3]
  real(r8), pointer :: RootPoreVol_pvr(:,:,:)        => null() !root layer volume air,                [m2 d-2]
  real(r8), pointer :: RootVH2O_pvr(:,:,:)            => null() !root layer volume water,              [m2 d-2]
  real(r8), pointer :: RootSAreaPerPlant_pvr(:,:,:)    => null() !layer root area per plant,            [m2 p-1]
  real(r8), pointer :: RootArea1stPP_pvr(:,:,:)       => null() !layer 1st root area per plant, [m2 plant-1]
  real(r8), pointer :: RootArea2ndPP_pvr(:,:,:)       => null() !layer 2nd root area per plant, [m2 plant-1]

  contains
    procedure, public :: Init    => plt_morph_init
    procedure, public :: Destroy => plt_morph_destroy
    procedure, public :: RefreshMediumRootMeanLength => plt_morph_refresh_medium_root_mean_length
  end type plant_morph_type


  type(plant_morph_type)    , public, target :: plt_morph     !plant morphology

contains

  subroutine plt_morph_init(this)
  implicit none
  class(plant_morph_type) :: this

  this%StandDeadSurfArea_col = 0._r8

  allocate(this%RootMedTransptArea_pvr(jroots,JZ1,JP1)); this%RootMedTransptArea_pvr=spval
  allocate(this%Root1stTransptArea_pvr(jroots,JZ1,JP1)); this%Root1stTransptArea_pvr=spval
  allocate(this%RootSAreaPerPlant_pvr(jroots,JZ1,JP1));this%RootSAreaPerPlant_pvr=0._r8
  allocate(this%RootLenDensPerPlant_pvr(jroots,JZ1,JP1));this%RootLenDensPerPlant_pvr=spval
  allocate(this%RootArea1stPP_pvr(jroots,JZ1,JP1)); this%RootArea1stPP_pvr=0._r8
  allocate(this%RootArea2ndPP_pvr(jroots,JZ1,JP1)); this%RootArea2ndPP_pvr=0._r8
  allocate(this%RootPoreVol_pvr(jroots,JZ1,JP1));this%RootPoreVol_pvr=spval
  allocate(this%RootVH2O_pvr(jroots,JZ1,JP1));this%RootVH2O_pvr=spval
  allocate(this%Root1stXNumL_pvr(JZ1,JP1));this%Root1stXNumL_pvr=spval
  allocate(this%fctyok_scalar_rpvr(JZ1,MaxNumRootAxes,JP1)); this%fctyok_scalar_rpvr=spval
  allocate(this%Root2ndXNumL_rpvr(jroots,JZ1,JP1));this%Root2ndXNumL_rpvr=spval
  allocate(this%RootMediumVH2O_rpvr(JZ1,MaxNumRootAxes,JP1)); this%RootMediumVH2O_rpvr=spval
  allocate(this%Root2ndVH2O_rpvr(jroots,JZ1,MaxNumRootAxes,JP1)); this%Root2ndVH2O_rpvr=0._r8
  allocate(this%Root1stVH2O_rpvr(JZ1,MaxNumRootAxes,JP1)); this%Root1stVH2O_rpvr=0._r8
  allocate(this%xylemPhi_min_pft(JP1)); this%xylemPhi_min_pft=0._r8
  allocate(this%xylemPhi_max_pft(JP1)); this%xylemPhi_max_pft=0._r8
  allocate(this%xylemPhi_mean_pft(JP1)); this%xylemPhi_mean_pft=0._r8
  allocate(this%SeedCMass_pft(JP1));this%SeedCMass_pft=spval
  allocate(this%SeedWidth2LenRatio_pft(JP1));this%SeedWidth2LenRatio_pft=spval
  allocate(this%totRootLenDens_vr(JZ1));this%totRootLenDens_vr=spval
  allocate(this%FineRootBranchFreq_pft(JP1));this%FineRootBranchFreq_pft=spval
  allocate(this%MediumRootBranchFreq_pft(JP1)); this%MediumRootBranchFreq_pft=spval
  allocate(this%ClumpFactorInit_pft(JP1));this%ClumpFactorInit_pft=spval
  allocate(this%ClumpFactorNow_pft(JP1));this%ClumpFactorNow_pft=spval
  allocate(this%HypocotHeight_pft(JP1));this%HypocotHeight_pft=spval
  allocate(this%RootPoreTortu4Gas_pft(jroots,JP1));this%RootPoreTortu4Gas_pft=spval
  allocate(this%rLen2WidthLeaf_pft(JP1));this%rLen2WidthLeaf_pft=spval
  allocate(this%MaxSeedNumPerSite_pft(JP1));this%MaxSeedNumPerSite_pft=spval
  allocate(this%GrothStalkMaxSeedSites_pft(JP1));this%GrothStalkMaxSeedSites_pft=spval
  allocate(this%SeedCMassMax_pft(JP1));this%SeedCMassMax_pft=spval
  allocate(this%Root1stMaxRadius1_pft(jroots,JP1));this%Root1stMaxRadius1_pft=spval
  allocate(this%Root2ndMaxRadius1_pft(jroots,JP1));this%Root2ndMaxRadius1_pft=spval
  allocate(this%RootRaidus_rpft(jroots,JP1));this%RootRaidus_rpft=spval
  allocate(this%Root1stRadius_pvr(jroots,JZ1,JP1));this%Root1stRadius_pvr=0._r8
  allocate(this%Root1stRadius_rpvr(JZ1,MaxNumRootAxes,JP1));this%Root1stRadius_rpvr=0._r8
  allocate(this%RootCRRadius0_rpvr(JZ1,MaxNumRootAxes,JP1)); this%RootCRRadius0_rpvr=0._r8
  allocate(this%Root2ndRadius_rpvr(jroots,JZ1,JP1));this%Root2ndRadius_rpvr=spval
  allocate(this%fRootTube_rpvr(JZ1,MaxNumRootAxes,JP1));this%fRootTube_rpvr=spval
  allocate(this%Root1stMaxRadius_pft(jroots,JP1));this%Root1stMaxRadius_pft=spval
  allocate(this%Root2ndMaxRadius_pft(jroots,JP1));this%Root2ndMaxRadius_pft=spval

  allocate(this%Root1stDepz_raxes(MaxNumRootAxes,JP1));this%Root1stDepz_raxes=spval
  allocate(this%Root1stAxesTipDepz2Surf_pft(MaxNumRootAxes,JP1)); this%Root1stAxesTipDepz2Surf_pft=spval
  allocate(this%Root1stLenLoc_rpvr(JZ1,MaxNumRootAxes,JP1)); this%Root1stLenLoc_rpvr=spval
  allocate(this%RootTotLenPerPlant_pvr(jroots,JZ1,JP1));this%RootTotLenPerPlant_pvr=spval
  allocate(this%RootAbsorbLenPerPlant_pvr(jroots,JZ1,JP1));this%RootAbsorbLenPerPlant_pvr=0._r8
  allocate(this%RootLenPerPlant_pvr(jroots,JZ1,JP1));this%RootLenPerPlant_pvr=0._r8
  allocate(this%Root2ndEffLen4uptk_rpvr(jroots,JZ1,JP1));this%Root2ndEffLen4uptk_rpvr=spval
  allocate(this%DistRootEffDepz_pvr(JZ1,JP1)); this%DistRootEffDepz_pvr=spval
  allocate(this%RootSinkScalar_pvr(JZ1,MaxNumRootAxes,JP1)); this%RootSinkScalar_pvr=spval
  allocate(this%Root1stSpecLen_pft(jroots,JP1));this%Root1stSpecLen_pft=spval
  allocate(this%Root2ndSpecLen_pft(jroots,JP1));this%Root2ndSpecLen_pft=spval
  allocate(this%Root1stLenPP_rpvr(JZ1,MaxNumRootAxes,JP1));this%Root1stLenPP_rpvr=spval
  allocate(this%flag2ndGrowth_pvr(JZ1,MaxNumRootAxes,JP1));this%flag2ndGrowth_pvr=.false.
  allocate(this%RootAge_rpvr(JZ1,MaxNumRootAxes,JP1)); this%RootAge_rpvr=spval
  allocate(this%Root2ndLen_rpvr(jroots,JZ1,MaxNumRootAxes,JP1));this%Root2ndLen_rpvr=spval
  allocate(this%CRootLumenArea_pvr(JZ1,JP1)); this%CRootLumenArea_pvr=0._r8
  allocate(this%MRootLumenArea_pvr(JZ1,JP1)); this%MRootLumenArea_pvr=0._r8
  allocate(this%MRootLumenArea_rpvr(JZ1,MaxNumRootAxes,JP1)); this%MRootLumenArea_rpvr=0._r8
  allocate(this%CRootLumenArea_rpvr(JZ1,MaxNumRootAxes,JP1)); this%CRootLumenArea_rpvr=0._r8
  allocate(this%Root2ndXNum_rpvr(jroots,JZ1,MaxNumRootAxes,JP1));this%Root2ndXNum_rpvr=0._r8
  allocate(this%RootMediumXNum_rpvr(JZ1,MaxNumRootAxes,JP1)); this%RootMediumXNum_rpvr=0._r8
  allocate(this%RootFineFrac2Med_rpvr(jroots,JZ1,JP1)); this%RootFineFrac2Med_rpvr=0._r8
  allocate(this%iPlantNfixType_pft(JP1));this%iPlantNfixType_pft=0
  allocate(this%Myco_pft(JP1));this%Myco_pft=0
  allocate(this%CanopyHeight4WatUptake_pft(JP1));this%CanopyHeight4WatUptake_pft=spval
  allocate(this%KLeafNumber_brch(MaxNumBranches,JP1));this%KLeafNumber_brch=0
  allocate(this%NumOfLeaves_brch(MaxNumBranches,JP1));this%NumOfLeaves_brch=spval
  allocate(this%NGTopRootLayer_pft(JP1));this%NGTopRootLayer_pft=0;
  allocate(this%RootSinkWeight_pvr(JZ1,JP1)); this%RootSinkWeight_pvr=0._r8
  allocate(this%Root2ndSinkWeight_pvr(JZ1,jroots,JP1));this%Root2ndSinkWeight_pvr=0._r8
  allocate(this%RootMSinkWeight_pvr(JZ1,JP1));this%RootMSinkWeight_pvr=0._r8
  allocate(this%Root1stSinkWeight_pvr(JZ1,JP1));this%Root1stSinkWeight_pvr=0._r8
  allocate(this%Root1stTipSinkWeight_pft(JP1)); this%Root1stTipSinkWeight_pft=0._r8
  allocate(this%RootMediumRadius_rpvr(JZ1,MaxNumRootAxes,JP1));this%RootMediumRadius_rpvr=0._r8
  allocate(this%NumMediumRootAxes_rpvr(JZ1,MaxNumRootAxes,JP1));this%NumMediumRootAxes_rpvr=0._r8
  allocate(this%RootMediumXNum_pvr(JZ1,JP1)); this%RootMediumXNum_pvr=0._r8
  allocate(this%RootMediumLength_pvr(JZ1,JP1)); this%RootMediumLength_pvr=0._r8
  allocate(this%RootMediumLength_rpvr(JZ1,MaxNumRootAxes,JP1));this%RootMediumLength_rpvr=0._r8
  allocate(this%RootMediumMeanLength_rpvr(JZ1,MaxNumRootAxes,JP1));this%RootMediumMeanLength_rpvr=0._r8
  allocate(this%Num1stAxesPerStructRootX_pft(JP1)); this%Num1stAxesPerStructRootX_pft=0._r8
  allocate(this%Num1stAxesPerStructRootXPOP_pft(JP1)); this%Num1stAxesPerStructRootXPOP_pft=0._r8
  allocate(this%Num1stRootAxesPP_pft(JP1)); this%Num1stRootAxesPP_pft=0._r8
  allocate(this%Radius95pctMature_pft(JP1)); this%Radius95pctMature_pft=0._r8
  allocate(this%CanopyHeightLive_pft(JP1));this%CanopyHeightLive_pft=spval
  allocate(this%CanopyHeightDead_pft(JP1)); this%CanopyHeightDead_pft=spval
  allocate(this%StalkHeight_pft(JP1)); this%StalkHeight_pft=spval
  allocate(this%StemSpecVolume_pft(JP1)); this%StemSpecVolume_pft=0._r8
  allocate(this%ShootNodeNumAtPlanting_pft(JP1));this%ShootNodeNumAtPlanting_pft=spval
  allocate(this%CanopyHeightZ_col(0:NumCanopyLayers1));this%CanopyHeightZ_col=spval
  allocate(this%CanopyStemSurfArea_pft(JP1));this%CanopyStemSurfArea_pft=spval
  allocate(this%CanopyLeafArea_pft(JP1));this%CanopyLeafArea_pft=spval
  allocate(this%MainBranchNum_pft(JP1));this%MainBranchNum_pft=0
  allocate(this%NMaxRootBotLayer_pft(JP1));this%NMaxRootBotLayer_pft=0
  allocate(this%NumStructuralRootAxes_pft(JP1));this%NumStructuralRootAxes_pft=0
  allocate(this%RootMatureAge_pft(JP1)); this%RootMatureAge_pft=0._r8
  allocate(this%RootSegAges_raxes(1:pltpar%NMaxRootSegs,1:MaxNumRootAxes,JP1));this%RootSegAges_raxes=0._r8
  allocate(this%NActiveRootSegs_raxes(1:MaxNumRootAxes,JP1));this%NActiveRootSegs_raxes=0
  allocate(this%RootSegBaseDepth_raxes(1:MaxNumRootAxes,JP1));this%RootSegBaseDepth_raxes=0._r8
  allocate(this%RootSeglengths_raxes(1:pltpar%NMaxRootSegs,1:MaxNumRootAxes,JP1));this%RootSeglengths_raxes=0._r8
  allocate(this%IndRootSegBase_raxes(1:MaxNumRootAxes,JP1));this%IndRootSegBase_raxes=0
  allocate(this%IndRootSegTip_raxes(1:MaxNumRootAxes,JP1));this%IndRootSegTip_raxes=0
  allocate(this%NumCogrowthNode_pft(JP1));this%NumCogrowthNode_pft=0
  allocate(this%BranchNumber_pft(JP1));this%BranchNumber_pft=0
  allocate(this%NumOfBranches_pft(JP1));this%NumOfBranches_pft=0
  allocate(this%NRoot1stTipLay_raxes(MaxNumRootAxes,JP1));this%NRoot1stTipLay_raxes=0
  allocate(this%PARTS_brch(NumOfPlantMorphUnits,MaxNumBranches,JP1));this%PARTS_brch=spval
  allocate(this%tlai_day_pft(JP1)); this%tlai_day_pft=spval
  allocate(this%enh_cyto_pft(JP1)); this%enh_cyto_pft=spval
  allocate(this%RootSingleVesselArea_pft(JP1));this%RootSingleVesselArea_pft=spval
  allocate(this%StalkAxialResist_pft(JP1));this%StalkAxialResist_pft=spval
  allocate(this%tsai_day_pft(JP1)); this%tsai_day_pft=spval
  allocate(this%ShootNodeNum_brch(MaxNumBranches,JP1));this%ShootNodeNum_brch=spval
  allocate(this%ShootNodeNumAtInitFloral_brch(MaxNumBranches,JP1));this%ShootNodeNumAtInitFloral_brch=spval
  allocate(this%ShootNodeNumAtAnthesis_brch(MaxNumBranches,JP1));this%ShootNodeNumAtAnthesis_brch=spval
  allocate(this%SineBranchAngle_pft(JP1));this%SineBranchAngle_pft=spval
  allocate(this%LeafAreaZsec_brch(NumLeafInclinationClasses1,NumCanopyLayers1,MaxNodesPerBranch1,MaxNumBranches,JP1))
  this%LeafAreaZsec_brch=spval
  allocate(this%KMinNumLeaf4GroAlloc_brch(MaxNumBranches,JP1));this%KMinNumLeaf4GroAlloc_brch=0
  allocate(this%BranchNumerID_brch(MaxNumBranches,JP1));this%BranchNumerID_brch=0
  allocate(this%PotentialSeedSites_brch(MaxNumBranches,JP1));this%PotentialSeedSites_brch=spval
  allocate(this%CanPBranchHeight(MaxNumBranches,JP1));this%CanPBranchHeight=spval
  allocate(this%LeafAreaDying_brch(MaxNumBranches,JP1));this%LeafAreaDying_brch=spval
  allocate(this%LeafAreaLive_brch(MaxNumBranches,JP1));this%LeafAreaLive_brch=spval
  allocate(this%SinePetolShethAngle_pft(JP1));this%SinePetolShethAngle_pft=spval
  allocate(this%LeafAngleClass_pft(NumLeafInclinationClasses1,JP1));this%LeafAngleClass_pft=spval
  allocate(this%CanopyStemSurfAreaZ_pft(NumCanopyLayers1,JP1));this%CanopyStemSurfAreaZ_pft=spval
  allocate(this%CanopyLeafAreaZ_pft(NumCanopyLayers1,JP1));this%CanopyLeafAreaZ_pft=spval
  allocate(this%LeafArea_node(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%LeafArea_node=spval
  allocate(this%CanopySeedNum_pft(JP1));this%CanopySeedNum_pft=spval
  allocate(this%CanopySeedNumX_pft(JP1));this%CanopySeedNumX_pft=spval
  allocate(this%SeedDepth_pft(JP1));this%SeedDepth_pft=spval
  allocate(this%PlantinDepz_pft(JP1));this%PlantinDepz_pft=spval
  allocate(this%SeedMeanLen_pft(JP1));this%SeedMeanLen_pft=spval
  allocate(this%SeedVolumeMean_pft(JP1));this%SeedVolumeMean_pft=spval
  allocate(this%SeedAreaMean_pft(JP1));this%SeedAreaMean_pft=spval
  allocate(this%CanopyStemAareZ_col(NumCanopyLayers1));this%CanopyStemAareZ_col=spval
  allocate(this%PetolShethLen2Mass_pft(JP1));this%PetolShethLen2Mass_pft=spval
  allocate(this%NodeLenPergC_pft(JP1));this%NodeLenPergC_pft=spval
  allocate(this%SLA1_pft(JP1));this%SLA1_pft=spval
  allocate(this%CanopyLeafAareZ_col(NumCanopyLayers1));this%CanopyLeafAareZ_col=spval
  allocate(this%LeafStalkAreaAct_pft(JP1));this%LeafStalkAreaAct_pft=spval
  allocate(this%StandDeadSurfAreaAct_pft(JP1)); this%StandDeadSurfAreaAct_pft=spval
  allocate(this%StalkNodeVertLength_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%StalkNodeVertLength_brch=spval
  allocate(this%PetoleLength_node(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%PetoleLength_node=spval
  allocate(this%StalkNodeHeight_brch(0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%StalkNodeHeight_brch=spval
  allocate(this%StemAreaZsec_brch(NumLeafInclinationClasses1,NumCanopyLayers1,MaxNumBranches,JP1));this%StemAreaZsec_brch=0._r8
  allocate(this%CanopyLeafArea_lnode(NumCanopyLayers1,0:MaxNodesPerBranch1,MaxNumBranches,JP1));this%CanopyLeafArea_lnode=0._r8
  allocate(this%CanopyStalkSurfArea_lbrch(NumCanopyLayers1,MaxNumBranches,JP1));this%CanopyStalkSurfArea_lbrch=spval
  allocate(this%CanopySurfAreaProfDead_pft(NumCanopyLayers1,JP1)); this%CanopySurfAreaProfDead_pft=spval
  allocate(this%StandDeadSurfArea_pft(JP1)); this%StandDeadSurfArea_pft=spval
  allocate(this%StalkAveRadius_pft(JP1));this%StalkAveRadius_pft=spval
  allocate(this%CanopyLeafAreaMAX_pft(JP1)); this%CanopyLeafAreaMAX_pft=spval
  allocate(this%lreset_laimax_pft(JP1)); this%lreset_laimax_pft=.false.
  allocate(this%MaxSoilLays4Root_pft(JP1));this%MaxSoilLays4Root_pft=0
  allocate(this%SetNumberSeeds_brch(MaxNumBranches,JP1));this%SetNumberSeeds_brch=spval
  allocate(this%ClumpFactor_pft(JP1));this%ClumpFactor_pft=spval
  allocate(this%FineRootVolPerMassC_pft(jroots,JP1));this%FineRootVolPerMassC_pft=spval
  allocate(this%CRootActVolPerMassC_pft(JP1)); this%CRootActVolPerMassC_pft=spval
  allocate(this%CRootLigVolPerMassC_pft(JP1)); this%CRootLigVolPerMassC_pft=spval
  allocate(this%RootPorosity_pft(jroots,JP1));this%RootPorosity_pft=spval
  allocate(this%Root2ndXSecArea_pft(jroots,JP1));this%Root2ndXSecArea_pft=spval
  allocate(this%Root1stXSecArea_pft(jroots,JP1));this%Root1stXSecArea_pft=spval
  allocate(this%RootRadialResist_pft(jroots,JP1));this%RootRadialResist_pft=spval
  allocate(this%RootSingleVesselRstaxial_pft(JP1));      this%RootSingleVesselRstaxial_pft=spval
  allocate(this%RootVesselRadius_pft(JP1));      this%RootVesselRadius_pft=spval
  allocate(this%Root2ndAxialResist_pft(jroots,JP1));this%Root2ndAxialResist_pft=spval
  allocate(this%iPlantGrainType_pft(JP1));this%iPlantGrainType_pft=0
  end subroutine plt_morph_init

  subroutine plt_morph_refresh_medium_root_mean_length(this,NZ)
  implicit none
  class(plant_morph_type), intent(inout) :: this
  integer, intent(in) :: NZ
  integer :: L,NR

  ! Derived from population totals; recompute before use, including after restart.
  this%RootMediumMeanLength_rpvr(:,:,NZ) = 0._r8
  DO NR=1,this%NumStructuralRootAxes_pft(NZ)
    DO L=1,SIZE(this%RootMediumLength_rpvr,1)
      IF(this%RootMediumXNum_rpvr(L,NR,NZ).GT.0._r8 .and. &
         this%RootMediumLength_rpvr(L,NR,NZ).GT.0._r8)THEN
        this%RootMediumMeanLength_rpvr(L,NR,NZ) = &
          this%RootMediumLength_rpvr(L,NR,NZ)/this%RootMediumXNum_rpvr(L,NR,NZ)
      ENDIF
    ENDDO
  ENDDO
  end subroutine plt_morph_refresh_medium_root_mean_length

  subroutine plt_morph_destroy(this)
  implicit none
  class(plant_morph_type) :: this

  if(associated(this%RootMediumMeanLength_rpvr))deallocate(this%RootMediumMeanLength_rpvr)
  if(associated(this%RootFineFrac2Med_rpvr))deallocate(this%RootFineFrac2Med_rpvr)
  end subroutine plt_morph_destroy
end module PlantMorphologyAPIData
