module PlantTraitDataType
  ! Description:
  ! Store plant traits and canopy geometry, including hourly exposed live and
  ! standing-dead areas used by radiation and atmospheric exchange.

!
!!
! data types of plant trait characteristics that cannot be grouped into canopy
! or roots
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__


!allocation parameter

  REAL(R8),target,allocatable :: FracLeafShethElmAlloc2Litr(:,:)             !fraction of shoot leaf element allocation to woody/fine litter,[-]
  real(r8),target,allocatable :: FracPetolShethAlloc2Litr(:,:)            !fraction of shoot stalk element allocation to woody/fine litter,[-]
  real(r8),target,allocatable :: FracRootElmAllocm(:,:)                  !fraction of root element allocation to woody/fine litter,[-]
  real(r8),target,allocatable :: FracWoodStalkElmAlloc2Litr(:,:)             !fraction of root stalk element allocation to woody/fine litter,[-]
  real(r8),target,allocatable :: PARTS_brch(:,:,:,:,:)                       !C partitioning coefficient in a branch, [-]
  real(r8),target,allocatable ::  CanopyStalkSurfArea_lbrch(:,:,:,:,:)       !Canopy stem layer area, [m2 d-2]
  real(r8),target,allocatable ::  CanopySurfAreaProfDead_pft(:,:,:,:)        !standing dead canopy surface area profile, [m2 d-2]
  real(r8),target,allocatable ::  CanopyLeafArea_pft(:,:,:)                  !Canopy leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyLeafAreaMAX_pft(:,:,:)               !running maximum leaf area index, [m2 d-2]
  real(r8),target,allocatable ::  LeafStalkAreaAct_pft(:,:,:)                   !plant canopy leaf+stem/stalk area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyStemSurfArea_pft(:,:,:)                  !plant stem area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyHeightLive_pft(:,:,:)                !live pft canopy height, [m]
  real(r8),target,allocatable ::  StandDeadSurfArea_pft(:,:,:)               !surface area of standing dead, [m2 d-2]
  real(r8),target,allocatable :: StandDeadSurfAreaAct_pft(:,:,:) ! hourly snow-exposed dead surface area [PFT,NY,NX], [m2]
  real(r8),target,allocatable ::  CanopyHeightDead_pft(:,:,:)                !standing dead plant canopy height, [m]
  real(r8),target,allocatable ::  StalkHeight_pft(:,:,:)                     !pft stalk height, [m]
  real(r8),target,allocatable ::  StalkAveRadius_pft(:,:,:)               !pft tree ring mean radius, [m]
  real(r8),target,allocatable ::  CanopyLeafAareZ_col(:,:,:)                 !total leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyStemAareZ_col(:,:,:)                 !total stem area, [m2 d-2]
  real(r8),target,allocatable ::  CanopyLeafArea_col(:,:)                    !grid level plant canopy leaf area, [m2 d-2]
  real(r8),target,allocatable ::  StemArea_col(:,:)                          !total canopy stem area, [m2 d-2]
  real(r8),target,allocatable ::  LeafStalkAreaAll_col(:,:)                     !canopy area of combined over the grid [m2 d-2]
  real(r8),target,allocatable ::  BulkFactor4Snow_col(:,:)                   !grid bulking factor for canopy snow interception effect on radiation, [m2 (kg SWE)-1]
  real(r8),target,allocatable ::  BulkFactor4Snow_pft(:,:,:)                 !pft bulking factor for canopy snow interception effect on radiation, [m2 (kg SWE)-1]
  integer ,target,allocatable ::  NGTopRootLayer_pft(:,:,:)                  !soil layer at planting depth, [-]
  real(r8),target,allocatable ::  PlantinDepz_pft(:,:,:)                     !planting depth, [m]
  real(r8),target,allocatable ::  SeedDepth_pft(:,:,:)                       !seeding depth, [m]
  real(r8),target,allocatable ::  SeedVolumeMean_pft(:,:,:)                  !seed volume, [m3 ]
  real(r8),target,allocatable ::  SeedMeanLen_pft(:,:,:)                     !seed length, [m]
  real(r8),target,allocatable ::  SeedAreaMean_pft(:,:,:)                    !seed surface area, [m2]
  real(r8),target,allocatable ::  HypocotHeight_pft(:,:,:)                   !cotyledon height, [m]
  real(r8),target,allocatable ::  CanopyHeight_col(:,:)                      !canopy height over grid, [m]
  real(r8),target,allocatable ::  CanopyHeightZ_col(:,:,:)                   !canopy layer height , [m]
  real(r8),target,allocatable ::  BranchAngle_pft(:,:,:)                     !branching angle, [degree from horizontal]
  real(r8),target,allocatable ::  PetolShethAngle_pft(:,:,:)                    !sheath angle, [degree from horizontal]
  real(r8),target,allocatable ::  SineBranchAngle_pft(:,:,:)                 !branching angle, [degree from horizontal]
  real(r8),target,allocatable ::  SinePetolShethAngle_pft(:,:,:)                !sheath angle, [degree from horizontal]
  real(r8),target,allocatable ::  RawIsoTCanopy2Atm_pft(:,:,:)           !ccanopy isothermal boundary later resistance, [h m-1]
  real(r8),target,allocatable ::  CanopyHeight4WatUptake_pft(:,:,:)          !effecive canopy height for water uptake, [m]
  real(r8),target,allocatable ::  LeafArea_node(:,:,:,:,:)               !leaf area, [m2 d-2]
  real(r8),target,allocatable ::  PetoleLength_node(:,:,:,:,:)             !sheath height, [m]
  real(r8),target,allocatable ::  StalkNodeHeight_brch(:,:,:,:,:)         !Live internode height, [m]
  real(r8),target,allocatable ::  LeafAreaLive_brch(:,:,:,:)                 !branch leaf area, [m2 d-2]
  real(r8),target,allocatable ::  LeafAreaDying_brch(:,:,:,:)                !branch leaf area, [m2 d-2]
  real(r8),target,allocatable ::  CanPBranchHeight(:,:,:,:)                  !branch height, [m]
  real(r8),target,allocatable ::  SetNumberSeeds_brch(:,:,:,:)                   !branch grain number, [d-2]
  real(r8),target,allocatable ::  PotentialSeedSites_brch(:,:,:,:)           !branch potential grain number, [d-2]
  real(r8),target,allocatable ::  CanopySeedNum_pft(:,:,:)                   !canopy grain number, [d-2]
  real(r8),target,allocatable ::  CanopySeedNumX_pft(:,:,:)                   !last nonzero canopy grain number, [d-2]  
  real(r8),target,allocatable ::  PlantPopuLive_pft(:,:,:)                   !live plant population, [d-2]
  real(r8),target,allocatable ::  PlantPopuDead_pft(:,:,:)                   !standing dead plant population, [d-2]
  real(r8),target,allocatable ::  StalkNodeVertLength_brch(:,:,:,:,:)        !Dead internode height, [m]
  REAL(r8),target,allocatable ::  rECDeadCRoot_pft(:,:,:,:)                  ! element:C ratio of dead coarse root, [gE gC-1]
  real(r8),target,allocatable ::  rECLiveCRoot_pft(:,:,:,:)                  ! element:C ratio of live coarse root, [gE gC-1]
  real(r8),target,allocatable ::  rNCLeaf_pft(:,:,:)                            !maximum leaf N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rPCLeaf_pft(:,:,:)                            !maximum leaf P:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCSheath_pft(:,:,:)                           !sheath N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCLigRoot_pft(:,:,:)         !Root lignified zone N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rPCLigRoot_pft(:,:,:)         !Root lignified zone P:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCStalk_pft(:,:,:)                        !stalk N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCReserve_pft(:,:,:)                      !reserve N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCHusk_pft(:,:,:)                         !husk N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCEar_pft(:,:,:)                          !ear N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCGrain_pft(:,:,:)                            !grain N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rNCNodule_pft(:,:,:)                       !nodule N:C ratio, [g g-1]
  real(r8),target,allocatable ::  rPCSheath_pft(:,:,:)                           !sheath P:C ratio, [g g-1]
  real(r8),target,allocatable ::  rPCStalk_pft(:,:,:)                        !stalk P:C ratio, [g g-1]
  real(r8),target,allocatable ::  KLigMax_pft(:,:,:)                         !Maximum lignification rate [h-1]
  real(r8),target,allocatable ::  KLigMM_pft(:,:,:)                          !Half saturation parameter for coarse root lignification, [h-1]
  real(r8),target,allocatable ::  rPCReserve_pft(:,:,:)                      !reserve P:C ratio, [g g-1]
  real(r8),target,allocatable ::  rPCHusk_pft(:,:,:)                         !husk P:C ratio, [g g-1]
  real(r8),target,allocatable ::  StalkAxialResist_pft(:,:,:)                !stalk axial resistance per m for water transport, [MPa h m-4]
  REAL(r8),target,allocatable ::  enh_cyto_pft(:,:,:)                        !cytokinin sensitivity of coarse root thickening
  real(r8),target,allocatable ::  rPCEar_pft(:,:,:)                          !ear P:C ratio, [g g-1]
  real(r8),target,allocatable ::  rPCGrain_pft(:,:,:)                            !grain P:C ratio, [g g-1]
  real(r8),target,allocatable ::  rPCNoduler_pft(:,:,:)                       !nodule P:C ratio, [g g-1]
  real(r8),target,allocatable ::  rProteinC2RootN_pft(:,:,:)                 !Protein C to root N ratio in remobilizable nonstructural biomass, [-]
  real(r8),target,allocatable ::  rProteinC2RootP_pft(:,:,:)                 !Protein C to root P ratio in remobilizable nonstructural biomass, [-]
  real(r8),target,allocatable ::  rProteinC2LeafN_pft(:,:,:)                 !Protein C to leaf N ratio in remobilizable nonstructural biomass, [-]
  real(r8),target,allocatable ::  rProteinC2LeafP_pft(:,:,:)                 !Protein C to leaf P ratio in remobilizable nonstructural biomass, [-]
  real(r8),target,allocatable ::  OrganOsmoPsi0pt_pft(:,:,:)                   !Organ osmotic potential when canopy water potential = 0 MPa, [MPa]
  real(r8),target,allocatable ::  TC4LeafOff_pft(:,:,:)                      !threshold temperature for autumn leafoff/hardening, [oC]
  real(r8),target,allocatable ::  PlantInitThermoAdaptZone_pft(:,:,:)            !initial plant thermal adaptation zone, [-]
  real(r8),target,allocatable ::  rPlantThermoAdaptZone_pft(:,:,:)           !plant thermal adaptation zone, [-]
  real(r8),target,allocatable ::  MatureGroup_brch(:,:,:,:)                  !plant maturity group, [-]
  real(r8),target,allocatable ::  MatureGroup_pft(:,:,:)                     !acclimated plant maturity group, [-]
  real(r8),target,allocatable ::  GROUPX_pft(:,:,:)                          !initial plant maturity group, [-]
  real(r8),target,allocatable ::  PPI_pft(:,:,:)                             !initial plant population, [# m-2]
  real(r8),target,allocatable ::  PPmax_pft(:,:,:)                            !maximum plant population, [# m-2]
  real(r8),target,allocatable ::  StandingDeadInitC_pft(:,:,:)               !initial standing dead C, [g C m-2]
  real(r8),target,allocatable ::  PPX_pft(:,:,:)                             !plant population, [# m-2]
  real(r8), target, allocatable :: RootSingleVesselRstaxial_pft(:,:,:)          !axial resistance for a single 1 m water transport vessel [MPa h m-4]
  real(r8), target, allocatable :: AlphaVesselAxialResist_pft(:,:,:)         !axial resistance corrector for a single 1 m water transport vessel [-]  
  real(r8), target, allocatable :: RootVesselRadius_pft(:,:,:)               !typical radius of the water transport vessel in primary roots, [m]
  real(r8), target, allocatable :: RootSingleVesselArea_pft(:,:,:)                  !single root vessel cross section area, [m2]
  integer,target,allocatable ::   NumActivePlants_col(:,:)                        !number of active PFT
  real(r8),target,allocatable ::  PlantPopu_col(:,:)                         !total plant population, [d-2]
  real(r8),target,allocatable ::  PPatSeeding_pft(:,:,:)                     !plant population at seeding, [m-2]
  real(r8),target,allocatable ::  HoursTooLowPsiCan_pft(:,:,:)               !canopy plant water stress indicator, number of hours PSILT < PSILY, []
  real(r8),target,allocatable ::  CanPhenoMoistStress_pft(:,:,:)              !moisture stress for plant phenology development,[-]
  real(r8),target,allocatable ::  CanPhenoTempStress_pft(:,:,:)               !temperature stress for plant phenology development,[-]  
  real(r8),target,allocatable ::  PlantO2Stress_pft(:,:,:)                   !plant O2 stress indicator, []
  real(r8),target,allocatable ::  MorphogenBase_pft(:,:,:)                   !base line morphogen concentration for secondary growth, [-]
  real(r8),target,allocatable ::  fTCanopyGroth_pft(:,:,:)                   !canopy temperature growth function, [-]
  real(r8),target,allocatable ::  TCGroth_pft(:,:,:)                         !canopy growth temperature, [oC]
  real(r8),target,allocatable ::  TKGroth_pft(:,:,:)                         !canopy growth temperature, [K]
  real(r8),target,allocatable ::  PetolShethBiomGrowthYld_pft(:,:,:)            !sheath growth yield, [g g-1]
  real(r8),target,allocatable ::  StalkBiomGrowthYld_pft(:,:,:)              !stalk growth yield, [g g-1]
  real(r8),target,allocatable ::  ReserveBiomGrowthYld_pft(:,:,:)            !reserve growth yield, [g g-1]
  real(r8),target,allocatable ::  HuskBiomGrowthYld_pft(:,:,:)               !husk growth yield, [g g-1]
  real(r8),target,allocatable ::  EarBiomGrowthYld_pft(:,:,:)                !ear growth yield, [g g-1]
  real(r8),target,allocatable ::  GrainBiomGrowthYld_pft(:,:,:)              !grain growth yield, [g g-1]
  real(r8),target,allocatable ::  NoduGrowthYield_pft(:,:,:)                 !nodule growth yield, [g g-1]
  real(r8),target,allocatable ::  LeafBiomGrowthYld_pft(:,:,:)               !leaf growth yield, [g g-1]
  real(r8),target,allocatable ::  Hours4LenthenPhotoPeriod_brch(:,:,:,:)     !initial heat requirement for spring leafout/dehardening, [h]
  real(r8),target,allocatable ::  Hours4ShortenPhotoPeriod_brch(:,:,:,:)     !initial cold requirement for autumn leafoff/hardening, [h]
  real(r8),target,allocatable ::  NumOfLeaves_brch(:,:,:,:)                  !leaf number, [-]
  real(r8),target,allocatable ::  LeafNumberAtFloralInit_brch(:,:,:,:)       !leaf number at floral initiation, [-]
  real(r8),target,allocatable ::  Hours4Leafout_brch(:,:,:,:)                !heat requirement for spring leafout/dehardening, [h]
  real(r8),target,allocatable ::  Hours4LeafOff_brch(:,:,:,:)                !cold requirement for autumn leafoff/hardening, [h]
  integer,target,allocatable ::  KLeafNumber_brch(:,:,:,:)                   !leaf number, [-]
  integer,target,allocatable ::  KLowestGroLeafNode_brch(:,:,:,:)            !leaf growth stage counter, [-]
  integer,target,allocatable ::  KMinNumLeaf4GroAlloc_brch(:,:,:,:)          !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION,[-]
  integer,target,allocatable ::  KHiestGroLeafNode_brch(:,:,:,:)             !leaf growth stage counter, [-]
  real(r8),target,allocatable ::  RateRefLeafAppearance_pft(:,:,:)               !rate of leaf initiation at 25 oC, [h-1]
  real(r8),target,allocatable ::  rLen2WidthLeaf_pft(:,:,:)                  !leaf length:width ratio, [-]
  real(r8),target,allocatable ::  SLA1_pft(:,:,:)                            !leaf area:mass during growth, [m2 g-1]
  real(r8),target,allocatable ::  TC4LeafOut_pft(:,:,:)                      !threshold temperature for spring leafout/dehardening, [oC]
  real(r8),target,allocatable ::  PetolShethLen2Mass_pft(:,:,:)                    !PetolSheth length:mass during growth, [m g-1]
  real(r8),target,allocatable ::  HourReq4LeafOut_brch(:,:,:,:)              !hours above threshold temperature required for spring leafout/dehardening, [h]
  integer,target,allocatable ::  NumOfBranches_pft(:,:,:)                    !number of branches of the plant, [-]
  integer,target,allocatable ::  BranchNumber_pft(:,:,:)                     !main branch number, [-]
  integer,target,allocatable ::  BranchNumerID_brch(:,:,:,:)                  !branch number id, [-]
  integer,target,allocatable ::  MainBranchNum_pft(:,:,:)                    !number of main branch, [-]
  integer,target,allocatable ::  Prep4Literfall_brch(:,:,:,:)                !branch phenology flag for senescence, [-]
  integer,target,allocatable ::  Hours4LiterfalAftMature_brch(:,:,:,:)       !hour counter for phenological senescence of a branch , [h]
  integer,target,allocatable ::  doSenescence_brch(:,:,:,:)                  !branch phenological senescence flag, [-]
  integer,target,allocatable ::  doRemobilization_brch(:,:,:,:)              !branch phenological remobilization flag, [-]
  integer,target,allocatable ::  doInitLeafOut_brch(:,:,:,:)                 !branch phenological flag for leafout initialization, [-]
  integer,target,allocatable ::  EnablePlantLeafOut_brch(:,:,:,:)                !branch phenological flag for leafout, [-]
  integer,target,allocatable ::  doPlantLeaveOff_brch(:,:,:,:)               !branch phenological flag for leaf off, [-]
  integer,target,allocatable ::  isPlantBranchAlive_brch(:,:,:,:)             !flag to detect branch death , [-]
  real(r8),target,allocatable ::  NonstCMinConc2InitBranch_pft(:,:,:)        !branch nonstructural C content required for new branch, [g g-1]
  real(r8),target,allocatable ::  NodeNumNormByMatgrp_brch(:,:,:,:)          !normalized node number during vegetative growth stages , [-]
  real(r8),target,allocatable ::  HourlyNodeNumNormByMatgrp_brch(:,:,:,:)    !gain in normalized node number during vegetative growth stages , [h-1]
  real(r8),target,allocatable ::  dReproNodeNumNormByMatG_brch(:,:,:,:)      !gain in normalized node number during reproductive growth stages, [h-1]
  real(r8),target,allocatable ::  ShootNodeNum_brch(:,:,:,:)                 !shoot node number, [-]
  real(r8),target,allocatable ::  ShootNodeNumAtInitFloral_brch(:,:,:,:)           !node number at floral initiation, [-]
  real(r8),target,allocatable ::  ReprodNodeNumNormByMatrgrp_brch(:,:,:,:)   !normalized node number during reproductive growth stages, [-]
  real(r8),target,allocatable ::  ShootNodeNumAtAnthesis_brch(:,:,:,:)         !node number at anthesis, [-]
  real(r8),target,allocatable ::  TotalNodeNumNormByMatgrp_brch(:,:,:,:)     !normalized node number during vegetative growth stages , [-]
  real(r8),target,allocatable ::  TotReproNodeNumNormByMatrgrp_brch(:,:,:,:) !normalized node number during reproductive growth stages , [-]
  real(r8),target,allocatable ::  RefNodeInitRate_pft(:,:,:)                 !rate of node initiation  at 25 oC, [h-1]
  real(r8),target,allocatable ::  NodeLenPergC_pft(:,:,:)                    !internode length:mass during growth, [m g-1]
  real(r8),target,allocatable ::  FracGroth2Node_pft(:,:,:)                  !parameter for allocation of growth to nodes, [-]
  integer,target,allocatable ::  NumCogrowthNode_pft(:,:,:)                  !number of concurrently growing nodes
  real(r8),target,allocatable ::  PSICanPDailyMin_pft(:,:,:)                     !minimum daily canopy water potential, [MPa]
  real(r8),target,allocatable ::  ClumpFactorNow_pft(:,:,:)                  !clumping factor for self-shading in canopy layer at current LAI, [-]
  real(r8),target,allocatable ::  ClumpFactor_pft(:,:,:)                     !clumping factor for self-shading in canopy layer, [-]
  integer,target,allocatable ::  isPlantShootAlive_pft(:,:,:)                 !flag to detect canopy death , [-]
  real(r8),target,allocatable ::  GrothStalkMaxSeedSites_pft(:,:,:)             !maximum grain node number per branch, [-]
  real(r8),target,allocatable ::  MaxSeedNumPerSite_pft(:,:,:)               !maximum grain number per node , [-]
  real(r8),target,allocatable ::  SeedCMassMax_pft(:,:,:)                    !maximum grain size   , [g]
  real(r8),target,allocatable ::  ShootNodeNumAtPlanting_pft(:,:,:)          !number of nodes in seed, [-]
  real(r8),target,allocatable ::  SeedCMass_pft(:,:,:)                       !grain size at seeding, [gC/seed]
  real(r8),target,allocatable ::  SeedWidth2LenRatio_pft(:,:,:)               !Seed width to length ratio, assuming prolate spheroid  
  real(r8),target,allocatable ::  GrainFillRate25C_pft(:,:,:)                !maximum rate of fill per grain, [g h-1]
  real(r8),target,allocatable ::  HourFailGrainFill_brch(:,:,:,:)            !flag to detect physiological maturity from  grain fill , [-]
  real(r8),target,allocatable ::  Hours2LeafOut_brch(:,:,:,:)                !counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
  real(r8),target,allocatable ::  HoursDoingRemob_brch(:,:,:,:)              !counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
  integer,target,allocatable ::  iPlantCalendar_brch(:,:,:,:,:)              !plant branch growth stage, [-]
  real(r8),target,allocatable ::  TCChill4Seed_pft(:,:,:)                    !temperature below which seed set is adversely affected, [oC]
  real(r8),target,allocatable ::  HighTempLimitSeed_pft(:,:,:)               !temperature above which seed set is adversely affected, [oC]
  real(r8),target,allocatable ::  SeedTempSens_pft(:,:,:)                    !sensitivity to canopy temperature, [oC-1]
  real(r8),target,allocatable ::  CriticPhotoPeriod_pft(:,:,:)               !critical daylength for phenological progress, [h]
  real(r8),target,allocatable ::  PhotoPeriodSens_pft(:,:,:)                 !difference between current and critical daylengths used to calculate  phenological progress, [h]
  real(r8),target,allocatable ::  ClumpFactorInit_pft(:,:,:)                 !initial clumping factor for self-shading in canopy layer, [-]
  real(r8),target,allocatable ::  HourReq4LeafOff_brch(:,:,:,:)              !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8),target,allocatable ::  StemSpecVolume_pft(:,:,:)                  !stalk specific volume, [m3 gC-1]
  real(r8),target,allocatable ::  TempOffset_pft(:,:,:)                      !adjustment of Arhhenius curves for plant thermal acclimation, [oC]
  integer,target,allocatable ::  iEmbryophyteType_pft(:,:,:)                 !plant embroytype [Bryophytes,Pteridophytes,Gymnosperms,Monocots and Eudicots]
  integer,target,allocatable ::  iPlantPhotosynsType_pft(:,:,:)             !plant photosynthetic type (C3 or C4),[-]
  integer,target,allocatable ::  iPlantRootProfile_pft(:,:,:)                !plant growth type (vascular, non-vascular),[-]
  integer,target,allocatable ::  iPlantSnowIntercepType_pft(:,:,:)           !plant snow interception type, (bryophyte, grasses, shrubs, deciduous trees, confierous trees)
  integer,target,allocatable ::  iPlantPhenolPattern_pft(:,:,:)              !plant growth habit (annual or perennial),[-]
  integer,target,allocatable ::  iPlantDevelopPattern_pft(:,:,:)             !plant growth habit (determinate or indeterminate),[-]  
  integer,target,allocatable ::  iPlantNfixType_pft(:,:,:)                   !N2 fixation type,[-]
  integer,target,allocatable ::  Days4FalseBreak_pft(:,:,:)                  !accumulated days to singifying false break, [day]    
  integer,target,allocatable ::  iPlantPhenolType_pft(:,:,:)                 !climate signal for phenological progress none, temperature, water stress,[-]
  integer,target,allocatable ::  iPlantPhotoperiodType_pft(:,:,:)            !photoperiod type (neutral, long day, short day),[-]
  integer,target,allocatable ::  iPlantTurnoverPattern_pft(:,:,:)            !phenologically-driven above-ground turnover (all, foliar only, none),[-]
  integer,target,allocatable :: iPlant2ndGrothPattern_pft(:,:,:)             !does the plant express secondary growth, [-]
  integer,target,allocatable ::  iPlantGrainType_pft(:,:,:)                  !grain type (below or above-ground), e.g. potato and onion are below,[-]
  integer,target,allocatable ::  Myco_pft(:,:,:)                               !mycorrhizal type, 1, 2 ,[-]
  

!----------------------------------------------------------------------

contains
  subroutine InitPlantTraits(NumOfPlantLitrCmplxs)

  implicit none
  integer, intent(in) :: NumOfPlantLitrCmplxs

  allocate(FracPetolShethAlloc2Litr(NumPlantChemElms,1:NumOfPlantLitrCmplxs));  FracPetolShethAlloc2Litr=0._r8
  allocate(FracLeafShethElmAlloc2Litr(NumPlantChemElms,1:NumOfPlantLitrCmplxs));  FracLeafShethElmAlloc2Litr=0._r8
  allocate(FracRootElmAllocm(NumPlantChemElms,1:NumOfPlantLitrCmplxs));  FracRootElmAllocm=0._r8         !
  allocate(FracWoodStalkElmAlloc2Litr(NumPlantChemElms,1:NumOfPlantLitrCmplxs));  FracWoodStalkElmAlloc2Litr=0._r8         !woody element allocation
  allocate(CanopyStalkSurfArea_lbrch(NumCanopyLayers,MaxNumBranches,JP,JY,JX));CanopyStalkSurfArea_lbrch=0._r8
  allocate(CanopySurfAreaProfDead_pft(NumCanopyLayers,JP,JY,JX));CanopySurfAreaProfDead_pft=0._r8
  allocate(iPlantSnowIntercepType_pft(JP,JY,JX)); iPlantSnowIntercepType_pft=0
  allocate(CanopyLeafArea_pft(JP,JY,JX));    CanopyLeafArea_pft=0._r8
  allocate(CanopyLeafAreaMAX_pft(JP,JY,JX)); CanopyLeafAreaMAX_pft=0._r8
  allocate(LeafStalkAreaAct_pft(JP,JY,JX));    LeafStalkAreaAct_pft=0._r8
  allocate(CanopyStemSurfArea_pft(JP,JY,JX));    CanopyStemSurfArea_pft=0._r8
  allocate(CanopyHeightLive_pft(JP,JY,JX));       CanopyHeightLive_pft=0._r8
  allocate(CanopyHeightDead_pft(JP,JY,JX)); CanopyHeightDead_pft=0._r8
  allocate(StandDeadSurfArea_pft(JP,JY,JX)); StandDeadSurfArea_pft=0._r8
  allocate(StandDeadSurfAreaAct_pft(JP,JY,JX)); StandDeadSurfAreaAct_pft=0._r8
  allocate(StalkHeight_pft(JP,JY,JX)); StalkHeight_pft=0._r8
  allocate(StalkAveRadius_pft(JP,JY,JX)); StalkAveRadius_pft=0._r8
  allocate(CanopyLeafAareZ_col(NumCanopyLayers,JY,JX));    CanopyLeafAareZ_col=0._r8
  allocate(CanopyStemAareZ_col(NumCanopyLayers,JY,JX));    CanopyStemAareZ_col=0._r8
  allocate(CanopyLeafArea_col(JY,JX));       CanopyLeafArea_col=0._r8
  allocate(StemArea_col(JY,JX));       StemArea_col=0._r8
  allocate(LeafStalkAreaAll_col(JY,JX));       LeafStalkAreaAll_col=0._r8
  allocate(BulkFactor4Snow_col(JY,JX)); BulkFactor4Snow_col=0._r8
  allocate(BulkFactor4Snow_pft(JP,JY,JX)); BulkFactor4Snow_pft=0._r8
  allocate(NGTopRootLayer_pft(JP,JY,JX));       NGTopRootLayer_pft=0
  allocate(PlantinDepz_pft(JP,JY,JX));   PlantinDepz_pft=0._r8
  allocate(SeedDepth_pft(JP,JY,JX));    SeedDepth_pft=0._r8
  allocate(SeedVolumeMean_pft(JP,JY,JX));     SeedVolumeMean_pft=0._r8
  allocate(SeedMeanLen_pft(JP,JY,JX));     SeedMeanLen_pft=0._r8
  allocate(SeedAreaMean_pft(JP,JY,JX));     SeedAreaMean_pft=0._r8
  allocate(HypocotHeight_pft(JP,JY,JX));    HypocotHeight_pft=0._r8
  allocate(CanopyHeight_col(JY,JX));          CanopyHeight_col=0._r8
  allocate(CanopyHeightZ_col(0:NumCanopyLayers,JY,JX));     CanopyHeightZ_col=0._r8
  allocate(BranchAngle_pft(JP,JY,JX));    BranchAngle_pft=0._r8
  allocate(PetolShethAngle_pft(JP,JY,JX));    PetolShethAngle_pft=0._r8
  allocate(SineBranchAngle_pft(JP,JY,JX));    SineBranchAngle_pft=0._r8
  allocate(SinePetolShethAngle_pft(JP,JY,JX));    SinePetolShethAngle_pft=0._r8
  allocate(RawIsoTCanopy2Atm_pft(JP,JY,JX));      RawIsoTCanopy2Atm_pft=0._r8
  allocate(CanopyHeight4WatUptake_pft(JP,JY,JX));    CanopyHeight4WatUptake_pft=0._r8
  allocate(PARTS_brch(NumOfPlantMorphUnits,MaxNumBranches,JP,JY,JX));PARTS_brch=0._r8
  allocate(LeafArea_node(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));LeafArea_node=0._r8
  allocate(PetoleLength_node(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));PetoleLength_node=0._r8
  allocate(StalkNodeHeight_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));StalkNodeHeight_brch=0._r8
  allocate(LeafAreaLive_brch(MaxNumBranches,JP,JY,JX)); LeafAreaLive_brch=0._r8
  allocate(LeafAreaDying_brch(MaxNumBranches,JP,JY,JX)); LeafAreaDying_brch=0._r8
  allocate(CanPBranchHeight(MaxNumBranches,JP,JY,JX));CanPBranchHeight=0._r8
  allocate(SetNumberSeeds_brch(MaxNumBranches,JP,JY,JX)); SetNumberSeeds_brch=0._r8
  allocate(PotentialSeedSites_brch(MaxNumBranches,JP,JY,JX)); PotentialSeedSites_brch=0._r8
  allocate(CanopySeedNum_pft(JP,JY,JX));     CanopySeedNum_pft=0._r8
  allocate(CanopySeedNumX_pft(JP,JY,JX));     CanopySeedNumX_pft=0._r8  
  allocate(PlantPopuLive_pft(JP,JY,JX));       PlantPopuLive_pft=0._r8
  allocate(PlantPopuDead_pft(JP,JY,JX)); PlantPopuDead_pft=0._r8
  allocate(rECDeadCRoot_pft(NumPlantChemElms,JP,JY,JX)); rECDeadCRoot_pft=0._r8
  allocate(rECLiveCRoot_pft(NumPlantChemElms,JP,JY,JX)); rECLiveCRoot_pft=0._r8
  allocate(StalkNodeVertLength_brch(0:MaxNodesPerBranch,MaxNumBranches,JP,JY,JX));StalkNodeVertLength_brch=0._r8
  allocate(rNCLeaf_pft(JP,JY,JX));     rNCLeaf_pft=0._r8
  allocate(rPCLeaf_pft(JP,JY,JX));     rPCLeaf_pft=0._r8
  allocate(rNCSheath_pft(JP,JY,JX));    rNCSheath_pft=0._r8
  allocate(rNCStalk_pft(JP,JY,JX));    rNCStalk_pft=0._r8
  allocate(rNCLigRoot_pft(JP,JY,JX)); rNCLigRoot_pft=0._r8
  allocate(rPCLigRoot_pft(JP,JY,JX)); rPCLigRoot_pft=0._r8
  allocate(rNCReserve_pft(JP,JY,JX));    rNCReserve_pft=0._r8
  allocate(rNCHusk_pft(JP,JY,JX));    rNCHusk_pft=0._r8
  allocate(rNCEar_pft(JP,JY,JX));    rNCEar_pft=0._r8
  allocate(rNCGrain_pft(JP,JY,JX));     rNCGrain_pft=0._r8
  allocate(rNCNodule_pft(JP,JY,JX));     rNCNodule_pft=0._r8
  allocate(rPCSheath_pft(JP,JY,JX));    rPCSheath_pft=0._r8
  allocate(rPCStalk_pft(JP,JY,JX));    rPCStalk_pft=0._r8
  allocate(KLigMax_pft(JP,JY,JX)); KLigMax_pft=0._r8
  allocate(KLigMM_pft(JP,JY,JX)); KLigMM_pft=0._r8
  allocate(rPCReserve_pft(JP,JY,JX));    rPCReserve_pft=0._r8
  allocate(rPCHusk_pft(JP,JY,JX));    rPCHusk_pft=0._r8
  allocate(rPCEar_pft(JP,JY,JX));    rPCEar_pft=0._r8
  allocate(StalkAxialResist_pft(JP,JY,JX)); StalkAxialResist_pft=0._r8
  allocate(enh_cyto_pft(JP,JY,JX));  enh_cyto_pft=0._r8
  allocate(rPCGrain_pft(JP,JY,JX));     rPCGrain_pft=0._r8
  allocate(rPCNoduler_pft(JP,JY,JX));     rPCNoduler_pft=0._r8
  allocate(rProteinC2RootN_pft(JP,JY,JX)); rProteinC2RootN_pft=0._r8
  allocate(rProteinC2RootP_pft(JP,JY,JX)); rProteinC2RootP_pft=0._r8  
  allocate(rProteinC2LeafN_pft(JP,JY,JX));     rProteinC2LeafN_pft=0._r8
  allocate(rProteinC2LeafP_pft(JP,JY,JX));     rProteinC2LeafP_pft=0._r8
  allocate(OrganOsmoPsi0pt_pft(JP,JY,JX));     OrganOsmoPsi0pt_pft=0._r8
  allocate(TC4LeafOff_pft(JP,JY,JX));      TC4LeafOff_pft=0._r8
  allocate(PlantInitThermoAdaptZone_pft(JP,JY,JX));    PlantInitThermoAdaptZone_pft=0._r8
  allocate(rPlantThermoAdaptZone_pft(JP,JY,JX));     rPlantThermoAdaptZone_pft=0._r8
  allocate(MatureGroup_brch(MaxNumBranches,JP,JY,JX)); MatureGroup_brch=0._r8
  allocate(MatureGroup_pft(JP,JY,JX));   MatureGroup_pft=0._r8
  allocate(GROUPX_pft(JP,JY,JX));   GROUPX_pft=0._r8
  allocate(PPI_pft(JP,JY,JX));      PPI_pft=0._r8
  allocate(PPmax_pft(JP,JY,JX));    PPmax_pft=0._r8
  allocate(StandingDeadInitC_pft(JP,JY,JX));   StandingDeadInitC_pft=0._r8
  allocate(PPX_pft(JP,JY,JX));      PPX_pft=0._r8
  allocate(RootSingleVesselArea_pft(JP,JY,JX)); RootSingleVesselArea_pft=0._r8
  allocate(RootSingleVesselRstaxial_pft(JP,JY,JX)); RootSingleVesselRstaxial_pft=0._R8
  allocate(AlphaVesselAxialResist_pft(JP,JY,JX)); AlphaVesselAxialResist_pft=0._R8  
  allocate(RootVesselRadius_pft(JP,JY,JX)); RootVesselRadius_pft=0._r8
  allocate(NumActivePlants_col(JY,JX));       NumActivePlants_col=0
  allocate(PlantPopu_col(JY,JX));         PlantPopu_col=0._r8
  allocate(PPatSeeding_pft(JP,JY,JX));      PPatSeeding_pft=0._r8
  allocate(HoursTooLowPsiCan_pft(JP,JY,JX));     HoursTooLowPsiCan_pft=0._r8
  allocate(CanPhenoMoistStress_pft(JP,JY,JX)); CanPhenoMoistStress_pft=1._R8
  allocate(CanPhenoTempStress_pft(JP,JY,JX)); CanPhenoTempStress_pft=1._r8
  allocate(PlantO2Stress_pft(JP,JY,JX));     PlantO2Stress_pft=0._r8
  allocate(fTCanopyGroth_pft(JP,JY,JX));     fTCanopyGroth_pft=0._r8
  allocate(MorphogenBase_pft(JP,JY,JX)); MorphogenBase_pft=0._r8
  allocate(TCGroth_pft(JP,JY,JX));      TCGroth_pft=0._r8
  allocate(TKGroth_pft(JP,JY,JX));      TKGroth_pft=0._r8
  allocate(PetolShethBiomGrowthYld_pft(JP,JY,JX));    PetolShethBiomGrowthYld_pft=0._r8
  allocate(StalkBiomGrowthYld_pft(JP,JY,JX));    StalkBiomGrowthYld_pft=0._r8
  allocate(ReserveBiomGrowthYld_pft(JP,JY,JX));    ReserveBiomGrowthYld_pft=0._r8
  allocate(HuskBiomGrowthYld_pft(JP,JY,JX));    HuskBiomGrowthYld_pft=0._r8
  allocate(EarBiomGrowthYld_pft(JP,JY,JX));    EarBiomGrowthYld_pft=0._r8
  allocate(GrainBiomGrowthYld_pft(JP,JY,JX));     GrainBiomGrowthYld_pft=0._r8
  allocate(NoduGrowthYield_pft(JP,JY,JX));     NoduGrowthYield_pft=0._r8
  allocate(LeafBiomGrowthYld_pft(JP,JY,JX));     LeafBiomGrowthYld_pft=0._r8
  allocate(Hours4LenthenPhotoPeriod_brch(MaxNumBranches,JP,JY,JX));  Hours4LenthenPhotoPeriod_brch=0._r8
  allocate(Hours4ShortenPhotoPeriod_brch(MaxNumBranches,JP,JY,JX));  Hours4ShortenPhotoPeriod_brch=0._r8
  allocate(NumOfLeaves_brch(MaxNumBranches,JP,JY,JX));  NumOfLeaves_brch=0._r8
  allocate(LeafNumberAtFloralInit_brch(MaxNumBranches,JP,JY,JX)); LeafNumberAtFloralInit_brch=0._r8
  allocate(Hours4Leafout_brch(MaxNumBranches,JP,JY,JX));  Hours4Leafout_brch=0._r8
  allocate(Hours4LeafOff_brch(MaxNumBranches,JP,JY,JX));  Hours4LeafOff_brch=0._r8
  allocate(KLeafNumber_brch(MaxNumBranches,JP,JY,JX)); KLeafNumber_brch=0
  allocate(KLowestGroLeafNode_brch(MaxNumBranches,JP,JY,JX));KLowestGroLeafNode_brch=0
  allocate(KMinNumLeaf4GroAlloc_brch(MaxNumBranches,JP,JY,JX));KMinNumLeaf4GroAlloc_brch=0
  allocate(KHiestGroLeafNode_brch(MaxNumBranches,JP,JY,JX)); KHiestGroLeafNode_brch=0
  allocate(RateRefLeafAppearance_pft(JP,JY,JX));     RateRefLeafAppearance_pft=0._r8
  allocate(rLen2WidthLeaf_pft(JP,JY,JX));     rLen2WidthLeaf_pft=0._r8
  allocate(SLA1_pft(JP,JY,JX));     SLA1_pft=0._r8
  allocate(TC4LeafOut_pft(JP,JY,JX));      TC4LeafOut_pft=0._r8
  allocate(PetolShethLen2Mass_pft(JP,JY,JX));     PetolShethLen2Mass_pft=0._r8
  allocate(HourReq4LeafOut_brch(NumCanopyLayers,JP,JY,JX));  HourReq4LeafOut_brch=0._r8
  allocate(NumOfBranches_pft(JP,JY,JX));      NumOfBranches_pft=0
  allocate(BranchNumber_pft(JP,JY,JX));      BranchNumber_pft=0
  allocate(BranchNumerID_brch(MaxNumBranches,JP,JY,JX));  BranchNumerID_brch=0
  allocate(MainBranchNum_pft(JP,JY,JX));      MainBranchNum_pft=0
  allocate(Prep4Literfall_brch(MaxNumBranches,JP,JY,JX)); Prep4Literfall_brch=ifalse
  allocate(Hours4LiterfalAftMature_brch(MaxNumBranches,JP,JY,JX)); Hours4LiterfalAftMature_brch=0
  allocate(doSenescence_brch(MaxNumBranches,JP,JY,JX)); doSenescence_brch=0
  allocate(doRemobilization_brch(MaxNumBranches,JP,JY,JX)); doRemobilization_brch=0
  allocate(doInitLeafOut_brch(MaxNumBranches,JP,JY,JX)); doInitLeafOut_brch=0
  allocate(EnablePlantLeafOut_brch(MaxNumBranches,JP,JY,JX)); EnablePlantLeafOut_brch=0
  allocate(doPlantLeaveOff_brch(MaxNumBranches,JP,JY,JX)); doPlantLeaveOff_brch=0
  allocate(isPlantBranchAlive_brch(MaxNumBranches,JP,JY,JX)); isPlantBranchAlive_brch=iFalse
  allocate(NonstCMinConc2InitBranch_pft(JP,JY,JX));       NonstCMinConc2InitBranch_pft=0._r8
  allocate(NodeNumNormByMatgrp_brch(MaxNumBranches,JP,JY,JX)); NodeNumNormByMatgrp_brch=0._r8
  allocate(HourlyNodeNumNormByMatgrp_brch(MaxNumBranches,JP,JY,JX));HourlyNodeNumNormByMatgrp_brch=0._r8
  allocate(dReproNodeNumNormByMatG_brch(MaxNumBranches,JP,JY,JX));dReproNodeNumNormByMatG_brch=0._r8
  allocate(ShootNodeNum_brch(MaxNumBranches,JP,JY,JX));  ShootNodeNum_brch=0._r8
  allocate(ShootNodeNumAtInitFloral_brch(MaxNumBranches,JP,JY,JX)); ShootNodeNumAtInitFloral_brch=0._r8
  allocate(ReprodNodeNumNormByMatrgrp_brch(MaxNumBranches,JP,JY,JX)); ReprodNodeNumNormByMatrgrp_brch=0._r8
  allocate(ShootNodeNumAtAnthesis_brch(MaxNumBranches,JP,JY,JX)); ShootNodeNumAtAnthesis_brch=0._r8
  allocate(TotalNodeNumNormByMatgrp_brch(MaxNumBranches,JP,JY,JX));TotalNodeNumNormByMatgrp_brch=0._r8
  allocate(TotReproNodeNumNormByMatrgrp_brch(MaxNumBranches,JP,JY,JX));TotReproNodeNumNormByMatrgrp_brch=0._r8
  allocate(RefNodeInitRate_pft(JP,JY,JX));     RefNodeInitRate_pft=0._r8
  allocate(NodeLenPergC_pft(JP,JY,JX));     NodeLenPergC_pft=0._r8
  allocate(FracGroth2Node_pft(JP,JY,JX));     FracGroth2Node_pft=0._r8
  allocate(NumCogrowthNode_pft(JP,JY,JX));     NumCogrowthNode_pft=0
  allocate(PSICanPDailyMin_pft(JP,JY,JX));    PSICanPDailyMin_pft=0._r8
  allocate(ClumpFactorNow_pft(JP,JY,JX));      ClumpFactorNow_pft=0._r8
  allocate(ClumpFactor_pft(JP,JY,JX));       ClumpFactor_pft=0._r8
  allocate(isPlantShootAlive_pft(JP,JY,JX));    isPlantShootAlive_pft=0
  allocate(GrothStalkMaxSeedSites_pft(JP,JY,JX));     GrothStalkMaxSeedSites_pft=0._r8
  allocate(MaxSeedNumPerSite_pft(JP,JY,JX));     MaxSeedNumPerSite_pft=0._r8
  allocate(SeedCMassMax_pft(JP,JY,JX));     SeedCMassMax_pft=0._r8
  allocate(ShootNodeNumAtPlanting_pft(JP,JY,JX));     ShootNodeNumAtPlanting_pft=0._r8
  allocate(SeedCMass_pft(JP,JY,JX));     SeedCMass_pft=0._r8
  allocate(SeedWidth2LenRatio_pft(JP,JY,JX)); SeedWidth2LenRatio_pft=0._r8
  allocate(GrainFillRate25C_pft(JP,JY,JX));    GrainFillRate25C_pft=0._r8
  allocate(HourFailGrainFill_brch(MaxNumBranches,JP,JY,JX));  HourFailGrainFill_brch=0._r8
  allocate(Hours2LeafOut_brch(MaxNumBranches,JP,JY,JX));  Hours2LeafOut_brch=0._r8
  allocate(HoursDoingRemob_brch(MaxNumBranches,JP,JY,JX));  HoursDoingRemob_brch=0._r8
  allocate(iPlantCalendar_brch(NumGrowthStages,MaxNumBranches,JP,JY,JX));iPlantCalendar_brch=0
  allocate(TCChill4Seed_pft(JP,JY,JX));      TCChill4Seed_pft=0._r8
  allocate(HighTempLimitSeed_pft(JP,JY,JX));      HighTempLimitSeed_pft=0._r8
  allocate(SeedTempSens_pft(JP,JY,JX));     SeedTempSens_pft=0._r8
  allocate(CriticPhotoPeriod_pft(JP,JY,JX));      CriticPhotoPeriod_pft=0._r8
  allocate(PhotoPeriodSens_pft(JP,JY,JX));     PhotoPeriodSens_pft=0._r8
  allocate(ClumpFactorInit_pft(JP,JY,JX));      ClumpFactorInit_pft=0._r8
  allocate(HourReq4LeafOff_brch(NumCanopyLayers,JP,JY,JX));  HourReq4LeafOff_brch=0._r8
  allocate(TempOffset_pft(JP,JY,JX));    TempOffset_pft=0._r8
  allocate(StemSpecVolume_pft(JP,JY,JX)); StemSpecVolume_pft =0._r8
  allocate(iPlantPhotosynsType_pft(JP,JY,JX));    iPlantPhotosynsType_pft=0
  allocate(iEmbryophyteType_pft(JP,JY,JX)); iEmbryophyteType_pft=0
  allocate(iPlantRootProfile_pft(JP,JY,JX));    iPlantRootProfile_pft=0
  allocate(iPlantPhenolPattern_pft(JP,JY,JX));    iPlantPhenolPattern_pft=0
  allocate(iPlantDevelopPattern_pft(JP,JY,JX));    iPlantDevelopPattern_pft=0
  allocate(iPlantNfixType_pft(JP,JY,JX));    iPlantNfixType_pft=0
  allocate(iPlantPhenolType_pft(JP,JY,JX));    iPlantPhenolType_pft=0
  allocate(Days4FalseBreak_pft(JP,JY,JX)); Days4FalseBreak_pft=0
  allocate(iPlantPhotoperiodType_pft(JP,JY,JX));    iPlantPhotoperiodType_pft=0
  allocate(iPlantTurnoverPattern_pft(JP,JY,JX));    iPlantTurnoverPattern_pft=0
  allocate(iPlant2ndGrothPattern_pft(JP,JY,JX)); iPlant2ndGrothPattern_pft=0
  allocate(iPlantGrainType_pft(JP,JY,JX));    iPlantGrainType_pft=0
  allocate(Myco_pft(JP,JY,JX));       Myco_pft=0

  end subroutine InitPlantTraits

!----------------------------------------------------------------------
  subroutine DestructPlantTraits
  use abortutils, only : destroy
  implicit none

  call destroy(StemSpecVolume_pft)
  call destroy(iPlant2ndGrothPattern_pft)
  call destroy(iPlantSnowIntercepType_pft)
  call destroy(Days4FalseBreak_pft)
  call destroy(FracLeafShethElmAlloc2Litr)
  call destroy(FracRootElmAllocm)
  call destroy(FracWoodStalkElmAlloc2Litr)
  call destroy(CanopySurfAreaProfDead_pft)
  call destroy(CanopyStalkSurfArea_lbrch)
  call destroy(CanopyLeafAreaMAX_pft)
  call destroy(CanopyLeafArea_pft)
  call destroy(LeafStalkAreaAct_pft)
  call destroy(CanopyStemSurfArea_pft)
  call destroy(CanopyHeightDead_pft)
  call destroy(enh_cyto_pft)
  call destroy(StalkAxialResist_pft)
  call destroy(StandDeadSurfArea_pft)
  call destroy(StandDeadSurfAreaAct_pft)
  call destroy(CanopyHeightLive_pft)
  call destroy(StalkHeight_pft)
  call destroy(CanopyLeafAareZ_col)
  call destroy(CanopyStemAareZ_col)
  call destroy(CanopyLeafArea_col)
  call destroy(StemArea_col)
  call destroy(LeafStalkAreaAll_col)
  call destroy(NGTopRootLayer_pft)
  call destroy(BulkFactor4Snow_pft)
  call destroy(BulkFactor4Snow_col)
  call destroy(PlantinDepz_pft)
  call destroy(SeedDepth_pft)
  call destroy(PlantPopuDead_pft)
  call destroy(SeedVolumeMean_pft)
  call destroy(SeedMeanLen_pft)
  call destroy(SeedAreaMean_pft)
  call destroy(HypocotHeight_pft)
  call destroy(CanopyHeight_col)
  call destroy(CanopyHeightZ_col)
  call destroy(BranchAngle_pft)
  call destroy(PetolShethAngle_pft)
  call destroy(SineBranchAngle_pft)
  call destroy(SinePetolShethAngle_pft)
  call destroy(RawIsoTCanopy2Atm_pft)
  call destroy(CanopyHeight4WatUptake_pft)
  call destroy(LeafArea_node)
  call destroy(PetoleLength_node)
  call destroy(StalkNodeHeight_brch)
  call destroy(LeafAreaLive_brch)
  call destroy(LeafAreaDying_brch)
  call destroy(CanPBranchHeight)
  call destroy(SetNumberSeeds_brch)
  call destroy(PotentialSeedSites_brch)
  call destroy(CanopySeedNum_pft)
  call destroy(CanopySeedNumX_pft)  
  call destroy(PlantPopuLive_pft)
  call destroy(rECDeadCRoot_pft)
  call destroy(rECLiveCRoot_pft)
  call destroy(StalkNodeVertLength_brch)
  call destroy(PARTS_brch)
  call destroy(rNCLeaf_pft)
  call destroy(rPCLeaf_pft)
  call destroy(rNCSheath_pft)
  call destroy(rNCStalk_pft)
  call destroy(rNCLigRoot_pft)
  call destroy(rPCLigRoot_pft)
  call destroy(rNCReserve_pft)
  call destroy(rNCHusk_pft)
  call destroy(rNCEar_pft)
  call destroy(rNCGrain_pft)
  call destroy(rNCNodule_pft)
  call destroy(rPCSheath_pft)
  call destroy(rPCStalk_pft)
  call destroy(KLigMax_pft)
  call destroy(KLigMM_pft)
  call destroy(rPCReserve_pft)
  call destroy(rPCHusk_pft)
  call destroy(rPCEar_pft)
  call destroy(rPCGrain_pft)
  call destroy(rPCNoduler_pft)
  call destroy(rProteinC2RootP_pft)  
  call destroy(rProteinC2RootN_pft)
  call destroy(rProteinC2LeafN_pft)
  call destroy(rProteinC2LeafP_pft)
  call destroy(OrganOsmoPsi0pt_pft)
  call destroy(TC4LeafOff_pft)
  call destroy(PlantInitThermoAdaptZone_pft)
  call destroy(rPlantThermoAdaptZone_pft)
  call destroy(MatureGroup_brch)
  call destroy(MatureGroup_pft)
  call destroy(GROUPX_pft)
  call destroy(PPI_pft)
  call destroy(PPmax_pft)
  call destroy(StandingDeadInitC_pft)
  call destroy(PPX_pft)
  call destroy(RootSingleVesselArea_pft)
  call destroy(RootSingleVesselRstaxial_pft)
  call destroy(AlphaVesselAxialResist_pft)  
  call destroy(RootVesselRadius_pft)
  call destroy(NumActivePlants_col)
  call destroy(PlantPopu_col)
  call destroy(PPatSeeding_pft)
  call destroy(HoursTooLowPsiCan_pft)
  call destroy(CanPhenoTempStress_pft)
  call destroy(CanPhenoMoistStress_pft)
  call destroy(PlantO2Stress_pft)
  call destroy(fTCanopyGroth_pft)
  call destroy(MorphogenBase_pft)
  call destroy(TCGroth_pft)
  call destroy(TKGroth_pft)
  call destroy(PetolShethBiomGrowthYld_pft)
  call destroy(StalkBiomGrowthYld_pft)
  call destroy(ReserveBiomGrowthYld_pft)
  call destroy(HuskBiomGrowthYld_pft)
  call destroy(EarBiomGrowthYld_pft)
  call destroy(GrainBiomGrowthYld_pft)
  call destroy(NoduGrowthYield_pft)
  call destroy(LeafBiomGrowthYld_pft)
  call destroy(Hours4LenthenPhotoPeriod_brch)
  call destroy(Hours4ShortenPhotoPeriod_brch)
  call destroy(NumOfLeaves_brch)
  call destroy(LeafNumberAtFloralInit_brch)
  call destroy(Hours4Leafout_brch)
  call destroy(Hours4LeafOff_brch)
  call destroy(KLeafNumber_brch)
  call destroy(KLowestGroLeafNode_brch)
  call destroy(KMinNumLeaf4GroAlloc_brch)
  call destroy(KHiestGroLeafNode_brch)
  call destroy(RateRefLeafAppearance_pft)
  call destroy(rLen2WidthLeaf_pft)
  call destroy(SLA1_pft)
  call destroy(TC4LeafOut_pft)
  call destroy(PetolShethLen2Mass_pft)
  call destroy(HourReq4LeafOut_brch)
  call destroy(NumOfBranches_pft)
  call destroy(BranchNumber_pft)
  call destroy(BranchNumerID_brch)
  call destroy(MainBranchNum_pft)
  call destroy(Prep4Literfall_brch)
  call destroy(Hours4LiterfalAftMature_brch)
  call destroy(doSenescence_brch)
  call destroy(doRemobilization_brch)
  call destroy(doInitLeafOut_brch)
  call destroy(EnablePlantLeafOut_brch)
  call destroy(doPlantLeaveOff_brch)
  call destroy(isPlantBranchAlive_brch)
  call destroy(NonstCMinConc2InitBranch_pft)
  call destroy(NodeNumNormByMatgrp_brch)
  call destroy(HourlyNodeNumNormByMatgrp_brch)
  call destroy(dReproNodeNumNormByMatG_brch)
  call destroy(ShootNodeNum_brch)
  call destroy(ShootNodeNumAtInitFloral_brch)
  call destroy(ReprodNodeNumNormByMatrgrp_brch)
  call destroy(ShootNodeNumAtAnthesis_brch)
  call destroy(TotalNodeNumNormByMatgrp_brch)
  call destroy(TotReproNodeNumNormByMatrgrp_brch)
  call destroy(RefNodeInitRate_pft)
  call destroy(NodeLenPergC_pft)
  call destroy(FracGroth2Node_pft)
  call destroy(NumCogrowthNode_pft)
  call destroy(PSICanPDailyMin_pft)
  call destroy(ClumpFactorNow_pft)
  call destroy(ClumpFactor_pft)
  call destroy(isPlantShootAlive_pft)
  call destroy(GrothStalkMaxSeedSites_pft)
  call destroy(MaxSeedNumPerSite_pft)
  call destroy(SeedCMassMax_pft)
  call destroy(ShootNodeNumAtPlanting_pft)
  call destroy(SeedCMass_pft)
  call destroy(SeedWidth2LenRatio_pft)
  call destroy(GrainFillRate25C_pft)
  call destroy(HourFailGrainFill_brch)
  call destroy(Hours2LeafOut_brch)
  call destroy(HoursDoingRemob_brch)
  call destroy(iPlantCalendar_brch)
  call destroy(TCChill4Seed_pft)
  call destroy(HighTempLimitSeed_pft)
  call destroy(SeedTempSens_pft)
  call destroy(CriticPhotoPeriod_pft)
  call destroy(PhotoPeriodSens_pft)
  call destroy(ClumpFactorInit_pft)
  call destroy(HourReq4LeafOff_brch)
  call destroy(TempOffset_pft)
  call destroy(iEmbryophyteType_pft)
  call destroy(iPlantPhotosynsType_pft)
  call destroy(iPlantRootProfile_pft)
  call destroy(iPlantPhenolPattern_pft)
  call destroy(iPlantDevelopPattern_pft)
  call destroy(iPlantNfixType_pft)
  call destroy(iPlantPhenolType_pft)
  call destroy(iPlantPhotoperiodType_pft)
  call destroy(iPlantTurnoverPattern_pft)
  call destroy(iPlantGrainType_pft)
  call destroy(Myco_pft)

  end subroutine DestructPlantTraits

end module PlantTraitDataType
