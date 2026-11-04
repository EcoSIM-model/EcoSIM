module InitPlantMod
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use UnitMod,          only: units
  use minimathmod,      only: AZMAX1,isclose,isALeftinBC
  use EcoSiMParDataMod, only: pltpar
  use EcosimConst
  use EcoSIMConfig
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantPhotosynthesisAPIData, only : plt_photo
  use PlantRadiationAPIData, only : plt_rad
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantSoilChemistryAPIData, only : plt_soilchem
  use PlantAllometryAPIData, only : plt_allom
  use PlantBiomassAPIData, only : plt_biom
  use PlantEnergyWaterAPIData, only : plt_ew
  use PlantDisturbanceAPIData, only : plt_distb
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantRootBGCAPIData, only : plt_rbgc
  use TracerIDMod
  use ElmIDMod
  use PlantMathFuncMod
  use PlantBGCPars
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: StartPlants
  public :: InitPlantPhenoMorphoBio
  public :: InitRootMychorMorphoBio
  public :: ApplySeedBank
  public :: SeedDeposition
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  SUBROUTINE StartPlants(NZ1Q,NZ2Q)
!
!     THIS SUBROUTINE INITIALIZES ALL PLANT VARIABLES
!
  implicit none
  integer, intent(in) :: NZ1Q,NZ2Q

  integer :: K,L,M,NZ,NZ2X
  character(len=*), parameter :: subname='StartPlants'
!     begin_execution

  associate(                                                   &
    iYearPlanting_pft     => plt_distb%iYearPlanting_pft      ,& !input  :year of planting,[-]
    iYearPlantHarvest_pft => plt_distb%iYearPlantHarvest_pft  ,& !input  :year of harvest,[-]
    iDayPlanting_pft      => plt_distb%iDayPlanting_pft       ,& !input  :day of planting,[-]
    iDayPlantHarvest_pft  => plt_distb%iDayPlantHarvest_pft   ,& !input  :day of harvest,[-]
    AREA3                 => plt_site%AREA3                   ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    iYearCurrent          => plt_site%iYearCurrent            ,& !input  :current year,[-]    
    iHarvestDay_pft       => plt_distb%iHarvestDay_pft        ,& !input  :day of harvest, [-]
    iHarvestYear_pft      => plt_distb%iHarvestYear_pft       ,& !input  :year of harvest,[-]
    iPlantingDay_pft      => plt_distb%iPlantingDay_pft       ,& !input  :day of planting,[-]
    iPlantingYear_pft     => plt_distb%iPlantingYear_pft      ,& !input  :year of planting,[-]    
    IsPlantActive_pft     => plt_pheno%IsPlantActive_pft      ,& !input  :flag for living pft, [-]
    PPatSeeding_pft       => plt_site%PPatSeeding_pft         ,& !input  :plant population at seeding, [plants d-2]    
    PPI_pft               => plt_site%PPI_pft                 ,& !output :initial plant population, [plants d-2]    
    PPX_pft               => plt_site%PPX_pft                 ,& !output :plant population, [plants m-2]    
    NL                    => plt_site%NL                      ,& !input  :lowest soil layer number,[-]
    NP                    => plt_site%NP                      ,& !input  :current number of plant species,[-]
    NP0                   => plt_site%NP0                     ,& !input  :number of tracked slots, including fallow residue,[-]
    NU                    => plt_site%NU                      ,& !input  :current soil surface layer number, [-]
    PlantPopuLive_pft     => plt_site%PlantPopuLive_pft       ,& !input  :plant population, [d-2]
    PlantPopuDead_pft     => plt_site%PlantPopuDead_pft       ,& !input  :live+standing dead plant population, [d-2]
    ZERO                  => plt_site%ZERO                    ,& !input  :threshold zero for numerical stability, [-]
    iMaintPlantTrait_pft  => plt_pheno%iMaintPlantTrait_pft   ,& !inoput :flag for maintain or reset plant traits, [-]    
    ZERO4Groth_pft        => plt_biom%ZERO4Groth_pft          ,& !output :threshold zero for plang growth calculation, [-]
    ZERO4LeafVar_pft      => plt_biom%ZERO4LeafVar_pft        ,& !output :threshold zero for leaf calculation, [-]
    ZERO4Uptk_pft         => plt_rbgc%ZERO4Uptk_pft            & !output :threshold zero for uptake calculation, [-]
  )

  call PrintInfo('beg '//subname)
!
!     INITIALIZE SHOOT GROWTH VARIABLES
!
!     IsPlantActive_pft=PFT flag:0=not active,1=active
!     PPI,PPX=initial,current population (m-2)
!     CF,ClumpFactorInit_pft=current,initial clumping factor
!     H2OCuticleResist_pft=cuticular resistance to water (h m-1)
!     CO2CuticleResist_pft=cuticular resistance to CO2 (s m-1)
!     CNWS,rProteinC2LeafP_pft=protein:N,protein:P ratios
!     RootProteinCMax_pft=maximum root protein concentration (g g-1)
!     O2I=intercellular O2 concentration in C3,C4 PFT (umol mol-1)
!

  NZ2X=MIN(NZ2Q,NP)
  D9985: DO NZ=NZ1Q,NZ2X

    IF(iMaintPlantTrait_pft(NZ).EQ.iFalse)THEN
      iYearPlanting_pft(NZ)     = iPlantingYear_pft(NZ)
      iDayPlanting_pft(NZ)      = iPlantingDay_pft(NZ)
      iYearPlantHarvest_pft(NZ) = iHarvestYear_pft(NZ)
      iDayPlantHarvest_pft(NZ)  = iHarvestDay_pft(NZ)
    ENDIF

    IF(iYearPlanting_pft(NZ).EQ.iYearCurrent)THEN
      IF(PPI_pft(NZ).GT.0.0)THEN
        PlantPopuLive_pft(NZ) = PPI_pft(NZ)*AREA3(NU)
        PPatSeeding_pft(NZ)   = PPI_pft(NZ)
        PPX_pft(NZ)           = PPatSeeding_pft(NZ)
      ELSE
        PlantPopuLive_pft(NZ) = PPatSeeding_pft(NZ)*AREA3(NU)
        PPX_pft(NZ)           = PPatSeeding_pft(NZ)
      ENDIF
    ELSE
      PlantPopuLive_pft(NZ)=PPX_pft(NZ)*AREA3(NU)
    ENDIF
    PlantPopuDead_pft(NZ)=PlantPopuLive_pft(NZ)

    IF(IsPlantActive_pft(NZ).EQ.iFalse)THEN

      call InitShootGrowth(NZ)

      call PlantLitterFraction(NZ)

      call PFTThermalAcclimation(NZ)

      call InitDimensionsandUptake(NZ)

      call InitPlantPhenoMorphoBio(NZ)

      call InitMassBalance(NZ)

      call InitRootMychorMorphoBio(NZ)

      call InitSeedMorphoBio(NZ)

      call InitPlantHeatWater(NZ)

    ENDIF
  ENDDO D9985

  DO NZ=NZ1Q,NZ2X
    ZERO4Groth_pft(NZ)   = ZERO*PlantPopuLive_pft(NZ)
    ZERO4Uptk_pft(NZ)    = ZERO*PlantPopuLive_pft(NZ)/AREA3(NU)
    ZERO4LeafVar_pft(NZ) = AMIN1(ZERO*PlantPopuLive_pft(NZ)*1.0E+06_r8,1.e-8_r8)    
  ENDDO  
!
!     FILL OUT UNUSED ARRAYS
!
  D9986: DO NZ=NP0+1,JP1
    plt_bgcr%SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ) = 0._r8
    plt_bgcr%LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ)     = 0._r8
    plt_biom%StandDeadStrutElms_pft(1:NumPlantChemElms,NZ)         = 0._r8
    D6401: DO L=1,NL
      DO  K=1,pltpar%NumOfPlantLitrCmplxs
        plt_bgcr%LitrfallElms_pvr(1:NumPlantChemElms,1:jsken,K,L,NZ)=0._r8
      enddo
    ENDDO D6401
  ENDDO D9986
  call PrintInfo('end '//subname)    
  RETURN
  end associate
  END subroutine StartPlants

!----------------------------------------------------------------------------------------------------
  subroutine InitShootGrowth(NZ)

  implicit none
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='InitShootGrowth'

  associate(                                                                   &
    ClumpFactorInit_pft           => plt_morph%ClumpFactorInit_pft            ,& !input  :initial clumping factor for self-shading in canopy layer, [-]
    CuticleResist_pft             => plt_photo%CuticleResist_pft              ,& !input  :maximum stomatal resistance to vapor, [s m-1]
    rNCRoot_pft                   => plt_allom%rNCRoot_pft                    ,& !input  :root N:C ratio, [gN gC-1]
    rPCRootr_pft                  => plt_allom%rPCRootr_pft                   ,& !input  :root P:C ratio, [gP gC-1]
    iPlantPhotosynsType_pft       => plt_photo%iPlantPhotosynsType_pft        ,& !input  :plant photosynthetic type (C3 or C4),[-]    
    PPI_pft                       => plt_site%PPI_pft                         ,& !output :initial plant population, [plants d-2]
    rProteinC2RootP_pft           => plt_allom%rProteinC2RootP_pft            ,& !output :Protein C to root P ratio in remobilizable nonstructural biomass, [-]    
    rProteinC2RootN_pft           => plt_allom%rProteinC2RootN_pft            ,& !output :Protein C to root N ratio in remobilizable nonstructural biomass, [-]
    rProteinC2LeafN_pft           => plt_allom%rProteinC2LeafN_pft            ,& !output :Protein C to leaf N ratio in remobilizable nonstructural biomass, [-]
    rProteinC2LeafP_pft           => plt_allom%rProteinC2LeafP_pft            ,& !output :Protein C to leaf P ratio in remobilizable nonstructural biomass, [-]
    CO2CuticleResist_pft          => plt_photo%CO2CuticleResist_pft           ,& !output :maximum stomatal resistance to CO2, [s h-1]
    ClumpFactor_pft               => plt_morph%ClumpFactor_pft                ,& !output :clumping factor for self-shading in canopy layer, [-]
    H2OCuticleResist_pft          => plt_photo%H2OCuticleResist_pft           ,& !output :maximum stomatal resistance to vapor, [s h-1]
    O2I_pft                       => plt_photo%O2I_pft                        ,& !output :leaf gaseous O2 concentration, [umol m-3]
    PPX_pft                       => plt_site%PPX_pft                         ,& !output :plant population, [plants m-2]
    RootProteinCMax_pft           => plt_allom%RootProteinCMax_pft             & !output :maximum root protein concentration, [gC g root C-1]
  )
  call PrintInfo('beg '//subname)
  PPX_pft(NZ)               = PPI_pft(NZ)
  ClumpFactor_pft(NZ)       = ClumpFactorInit_pft(NZ)

  H2OCuticleResist_pft(NZ)          = CuticleResist_pft(NZ)/3600.0_r8
  CO2CuticleResist_pft(NZ)          = CuticleResist_pft(NZ)*1.56_r8    !1.56 = sqrt(44./18.)

  RootProteinCMax_pft(NZ) = AMIN1(rNCRoot_pft(NZ)*rProteinC2RootN_pft(NZ),rPCRootr_pft(NZ)*rProteinC2RootP_pft(NZ))  !groot N/g rootC * gprotein C/g root N = g protein C/g root C
  IF(iPlantPhotosynsType_pft(NZ).EQ.ic3_photo)THEN
    O2I_pft(NZ)=2.10E+05_r8
  ELSE
    O2I_pft(NZ)=3.96E+05_r8
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine InitShootGrowth

!----------------------------------------------------------------------------------------------------
  subroutine PlantLitterFraction(NZ)
  implicit none
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='PlantLitterFraction'
  integer :: N,M
  real(r8) :: CNOPC(4),CPOPC(4)
  REAL(R8) :: CNOPCT,CPOPCT,litter_sum
  character(len=256) :: msg

  associate(                                                             &
    MatureGroup_pft             => plt_pheno%MatureGroup_pft            ,& !input  :acclimated plant maturity group, [-]
    NumLitterGroups             => pltpar%NumLitterGroups               ,& !input  :number of litter groups nonstructural(0,*)
    RateRefLeafAppearance_pft   => plt_pheno%RateRefLeafAppearance_pft  ,& !input  :rate of leaf initiation, [h-1 at 25 oC]
    iPlantNfixType_pft          => plt_morph%iPlantNfixType_pft         ,& !input  :N2 fixation type,[-]
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft      ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlantTurnoverPattern_pft   => plt_pheno%iPlantTurnoverPattern_pft  ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    icarbhyro                   => pltpar%icarbhyro                     ,& !input  :kinetic id of litter component as carbonhydrate
    icellulos                   => pltpar%icellulos                     ,& !input  :kinetic id of litter component as cellulose
    iPlant2ndGrothPattern_pft   => plt_pheno%iPlant2ndGrothPattern_pft  ,& !input  :plant expression of secondary growth, [-]                
    icwood                      => pltpar%icwood                        ,& !input  :group id of coarse woody litter
    ifoliar                     => pltpar%ifoliar                       ,& !input  :group id of plant foliar litter
    ilignin                     => pltpar%ilignin                       ,& !input  :kinetic id of litter component as lignin
    inonfoliar                  => pltpar%inonfoliar                    ,& !input  :group id of plant non-foliar litter group
    inonstruct                  => pltpar%inonstruct                    ,& !input  :group id of plant nonstructural litter
    iprotein                    => pltpar%iprotein                      ,& !input  :kinetic id of litter component as protein
    iroot                       => pltpar%iroot                         ,& !input  :group id of plant root litter
    istalk                      => pltpar%istalk                        ,& !input  :group id of plant stalk litter group
    PlantElmAllocMat4Litr       => plt_soilchem%PlantElmAllocMat4Litr   ,& !inoput :litter kinetic fraction, [-]
    FracGroth2Node_pft          => plt_allom%FracGroth2Node_pft         ,& !output :parameter for allocation of growth to nodes, [-]
    NumCogrowthNode_pft         => plt_morph%NumCogrowthNode_pft         & !output :number of concurrently growing nodes,[-]
  )
  call PrintInfo('beg '//subname)
  !
  !     FRACTIONS OF PLANT LITTER ALLOCATED TO KINETIC COMPONENTS
  !     PROTEIN(*,1),CH2O(*,2),CELLULOSE(*,3),LIGNIN(*,4) IN SOIL LITTER
  !
  !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
  !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
  !
  !     NONSTRUCTURAL
  !
  PlantElmAllocMat4Litr(ielmc,inonstruct,iprotein,NZ)  = 0.0_r8
  PlantElmAllocMat4Litr(ielmc,inonstruct,icarbhyro,NZ) = 0.67_r8
  PlantElmAllocMat4Litr(ielmc,inonstruct,icellulos,NZ) = 0.33_r8
  PlantElmAllocMat4Litr(ielmc,inonstruct,ilignin,NZ)   = 0.0_r8
  !
  !     NON-VASCULAR (E.G. MOSSES)
  !
  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ)))THEN
    PlantElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ) = 0.25_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ) = 0.30_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ)   = 0.38_r8

    PlantElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ) = 0.25_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ) = 0.30_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ)   = 0.38_r8
    !
    !     LEGUMES
    !
  ELSEIF(is_plant_N2fix(iPlantNfixType_pft(NZ)))THEN
    PlantElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ)  = 0.16_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ) = 0.38_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ) = 0.34_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ)   = 0.12_r8

    PlantElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ) = 0.41_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ) = 0.37_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ)   = 0.15_r8
    !
    !     ANNUALS, GRASSES, SHRUBS
    !
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ))))THEN
    PlantElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ)  = 0.08_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ) = 0.41_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ) = 0.36_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ)   = 0.15_r8

    PlantElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ) = 0.41_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ) = 0.36_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ)   = 0.16_r8
    !
    !     DECIDUOUS TREES
    !
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.1 &
    .OR.iPlantTurnoverPattern_pft(NZ).GE.3)THEN
    PlantElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ) = 0.34_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ) = 0.36_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ)   = 0.23_r8

    PlantElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ)  = 0.0_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ) = 0.045_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ) = 0.660_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ)   = 0.295_r8
    !
    !     CONIFEROUS TREES
    !
  ELSE
    PlantElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ) = 0.25_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ) = 0.38_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ)   = 0.30_r8

    PlantElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ)  = 0.0_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ) = 0.045_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ) = 0.660_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ)   = 0.295_r8
  ENDIF
!
!     FRACTIONS OF WOODY LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     NON-VASCULAR
!
  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ)))THEN
    PlantElmAllocMat4Litr(ielmc,istalk,iprotein,NZ)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,istalk,icarbhyro,NZ) = 0.25_r8
    PlantElmAllocMat4Litr(ielmc,istalk,icellulos,NZ) = 0.30_r8
    PlantElmAllocMat4Litr(ielmc,istalk,ilignin,NZ)   = 0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ))))THEN
    PlantElmAllocMat4Litr(ielmc,istalk,iprotein,NZ)  = 0.03_r8
    PlantElmAllocMat4Litr(ielmc,istalk,icarbhyro,NZ) = 0.25_r8
    PlantElmAllocMat4Litr(ielmc,istalk,icellulos,NZ) = 0.57_r8
    PlantElmAllocMat4Litr(ielmc,istalk,ilignin,NZ)   = 0.15_r8
!
!     DECIDUOUS AND CONIFEROUS TREES
!
  ELSE
    PlantElmAllocMat4Litr(ielmc,istalk,iprotein,NZ)  = 0.0_r8
    PlantElmAllocMat4Litr(ielmc,istalk,icarbhyro,NZ) = 0.045_r8
    PlantElmAllocMat4Litr(ielmc,istalk,icellulos,NZ) = 0.660_r8
    PlantElmAllocMat4Litr(ielmc,istalk,ilignin,NZ)   = 0.295_r8
  ENDIF
!
!     FRACTIONS OF FINE ROOT LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN PC&E 25:601-608
!
!     NON-VASCULAR
!
  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ)))THEN
    PlantElmAllocMat4Litr(ielmc,iroot,iprotein,NZ)=0.07_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icarbhyro,NZ)=0.25_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icellulos,NZ)=0.30_r8
    PlantElmAllocMat4Litr(ielmc,iroot,ilignin,NZ)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ))))THEN
    PlantElmAllocMat4Litr(ielmc,iroot,iprotein,NZ)  = 0.057_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icarbhyro,NZ) = 0.263_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icellulos,NZ) = 0.542_r8
    PlantElmAllocMat4Litr(ielmc,iroot,ilignin,NZ)   = 0.138_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.1&
    .OR.iPlantTurnoverPattern_pft(NZ).GE.3)THEN
    PlantElmAllocMat4Litr(ielmc,iroot,iprotein,NZ)  = 0.059_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icarbhyro,NZ) = 0.308_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icellulos,NZ) = 0.464_r8
    PlantElmAllocMat4Litr(ielmc,iroot,ilignin,NZ)   = 0.169_r8
!
!     CONIFEROUS TREES
!
  ELSE
    PlantElmAllocMat4Litr(ielmc,iroot,iprotein,NZ)  = 0.059_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icarbhyro,NZ) = 0.308_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icellulos,NZ) = 0.464_r8
    PlantElmAllocMat4Litr(ielmc,iroot,ilignin,NZ)   = 0.169_r8
  ENDIF
  !
  !     COARSE WOODY LITTER FROM BOLES AND ROOTS
  !
  PlantElmAllocMat4Litr(ielmc,icwood,iprotein,NZ)  = 0.00_r8
  PlantElmAllocMat4Litr(ielmc,icwood,icarbhyro,NZ) = 0.045_r8
  PlantElmAllocMat4Litr(ielmc,icwood,icellulos,NZ) = 0.660_r8
  PlantElmAllocMat4Litr(ielmc,icwood,ilignin,NZ)   = 0.295_r8
  !
  !     Validate that carbon litter fractions across kinetic components sum to one.
  !
  D95: DO N=0,NumLitterGroups
    litter_sum = SUM(PlantElmAllocMat4Litr(ielmc,N,1:jsken,NZ))
    IF(.not.isclose(litter_sum,1.0_r8))THEN
      write(msg,'(A,I0,A,I0,A,ES16.8)') &
        'Error: PlantElmAllocMat4Litr C fractions for PFT ',NZ, &
        ', litter group ',N,' sum to ',litter_sum
      call endrun(trim(msg)//' in '//trim(mod_filename),__LINE__)
    ENDIF
  ENDDO D95
  !
  !     INITIALIZE C-N AND C-P RATIOS IN PLANT LITTER
  !
  !     CNOPC,CPOPC=fractions to allocate N,P to kinetic components
  !     CFOPN,CFOPP=distribution of litter N,P to kinetic components
  !
  CNOPC(iprotein)  = 0.020_r8
  CNOPC(icarbhyro) = 0.010_r8
  CNOPC(icellulos) = 0.010_r8
  CNOPC(ilignin)   = 0.020_r8

  CPOPC(iprotein)  = 0.0020_r8
  CPOPC(icarbhyro) = 0.0010_r8
  CPOPC(icellulos) = 0.0010_r8
  CPOPC(ilignin)   = 0.0020_r8

  D110: DO N=0,NumLitterGroups
    CNOPCT=0.0_r8
    CPOPCT=0.0_r8
    D100: DO M=1,jsken
      CNOPCT=CNOPCT+PlantElmAllocMat4Litr(ielmc,N,M,NZ)*CNOPC(M)
      CPOPCT=CPOPCT+PlantElmAllocMat4Litr(ielmc,N,M,NZ)*CPOPC(M)
    ENDDO D100
    D105: DO M=1,jsken
      PlantElmAllocMat4Litr(ielmn,N,M,NZ)=PlantElmAllocMat4Litr(ielmc,N,M,NZ)*CNOPC(M)/CNOPCT
      PlantElmAllocMat4Litr(ielmp,N,M,NZ)=PlantElmAllocMat4Litr(ielmc,N,M,NZ)*CPOPC(M)/CPOPCT
    ENDDO D105
  ENDDO D110
!
!     CONCURRENT NODE GROWTH
!
!     FracGroth2Node_pft=scales node number for perennial vegetation (e.g. trees)
!     NumCogrowthNode_pft=number of concurrently growing nodes
!     iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+PetolSheth,2=none,3=between 1,2!
  IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. (.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ))))THEN
    !Annual plant (e.g.crops) or grass or bryophyte
    !The following may be revised for maize or soybean, or crops in general
    FracGroth2Node_pft(NZ) = 1.0_r8
    IF(MatureGroup_pft(NZ).LE.6)THEN
      NumCogrowthNode_pft(NZ) = 2
    ELSEIF(MatureGroup_pft(NZ).LE.10)THEN
      NumCogrowthNode_pft(NZ) = 3
    ELSEIF(MatureGroup_pft(NZ).LE.15)THEN
      NumCogrowthNode_pft(NZ) = 4
    ELSE
      NumCogrowthNode_pft(NZ) = 5
    ENDIF
  ELSE
    !perrenial tree-like
    FracGroth2Node_pft(NZ)  = AMAX1(1.0_r8,0.04_r8/RateRefLeafAppearance_pft(NZ))
    NumCogrowthNode_pft(NZ) = MaxNodesPerBranch1-1
  ENDIF
  end associate
  end subroutine PlantLitterFraction

!----------------------------------------------------------------------------------------------------
  subroutine PFTThermalAcclimation(NZ)

  implicit none
  integer, intent(in) :: NZ
  real(r8), parameter :: TCZD = 5.0_r8        !basal value for threshold temperature for spring leafout/dehardening	oC
  real(r8), parameter :: TCXD = 12.0_r8       !basal value for threshold temperature for autumn leafoff/hardening	oC

  associate(                                                                   &
    DATAP                         => plt_site%DATAP                           ,& !input  :parameter file name,[-]
    PlantInitThermoAdaptZone_pft  => plt_pheno%PlantInitThermoAdaptZone_pft   ,& !input  :initial plant thermal adaptation zone, [-]
    iPlantPhotosynsType_pft       => plt_photo%iPlantPhotosynsType_pft        ,& !input  :plant photosynthetic type (C3 or C4),[-]
    TempOffset_pft                => plt_pheno%TempOffset_pft                 ,& !output :adjustment of Arhhenius curves for plant thermal acclimation, [oC]
    rPlantThermoAdaptZone_pft     => plt_pheno%rPlantThermoAdaptZone_pft      ,& !output :plant thermal adaptation zone, [-]
    HighTempLimitSeed_pft         => plt_pheno%HighTempLimitSeed_pft          ,& !output :temperature above which seed set is adversely affected, [oC]
    SeedTempSens_pft              => plt_pheno%SeedTempSens_pft               ,& !output :sensitivity to HTC (seeds oC-1 above HTC),[oC-1]
    TC4LeafOff_pft                => plt_pheno%TC4LeafOff_pft                 ,& !output :threshold temperature for autumn leafoff/hardening, [oC]
    TC4LeafOut_pft                => plt_pheno%TC4LeafOut_pft                  & !output :threshold temperature for spring leafout/dehardening, [oC]
  )
!
!     PFT THERMAL ACCLIMATION
!
!     ZTYP,PlantInitThermoAdaptZone_pft=dynamic,initial thermal adaptation zone from PFT file
!     TempOffset_pft=shift in Arrhenius curve for thermal adaptation (oC)
!     TCZ,TC4LeafOff_pft=threshold temperature for leafout,leafoff
!     HTC=high temperature threshold for grain number loss (oC)
!     SeedTempSens_pft=sensitivity to HTC (seeds oC-1 above HTC)
!
  rPlantThermoAdaptZone_pft(NZ) = PlantInitThermoAdaptZone_pft(NZ)
  !greater thermal adaptation zone, greater TC4LeafOut_pft and TC4LeafOff_pft values
  TempOffset_pft(NZ)            = 2.667_r8*(2.5_r8-rPlantThermoAdaptZone_pft(NZ))
  TC4LeafOut_pft(NZ)            = TCZD-TempOffset_pft(NZ)
  TC4LeafOff_pft(NZ)            = AMIN1(15.0_r8,TCXD-TempOffset_pft(NZ))
  IF(iPlantPhotosynsType_pft(NZ).EQ.ic3_photo)THEN
    IF(DATAP(NZ)(1:4).EQ.'soyb')THEN
      HighTempLimitSeed_pft(NZ)=30.0_r8+3.0_r8*rPlantThermoAdaptZone_pft(NZ)
      SeedTempSens_pft(NZ)=0.002_r8
    ELSE
      HighTempLimitSeed_pft(NZ)=27.0_r8+3.0_r8*rPlantThermoAdaptZone_pft(NZ)
      SeedTempSens_pft(NZ)=0.002_r8
    ENDIF
  ELSE
    HighTempLimitSeed_pft(NZ)=27.0_r8+3.0_r8*rPlantThermoAdaptZone_pft(NZ)
    SeedTempSens_pft(NZ)=0.005_r8
  ENDIF
  end associate
  end subroutine PFTThermalAcclimation

!----------------------------------------------------------------------------------------------------
  subroutine InitDimensionsandUptake(NZ)

  implicit none
  integer, intent(in) :: NZ
  INTEGER :: L,N,NR
  associate(                                                          &
    CumSoilThickness_vr       => plt_site%CumSoilThickness_vr        ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    NL                        => plt_site%NL                         ,& !input  :lowest soil layer number,[-]
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    PlantinDepz_pft           => plt_morph%PlantinDepz_pft           ,& !input  :planting depth, [m]
    RootBiomGrosYld_pft       => plt_allom%RootBiomGrosYld_pft       ,& !input  :root growth yield, [g g-1]
    rNCRoot_pft               => plt_allom%rNCRoot_pft               ,& !input  :root N:C ratio, [gN gC-1]
    rPCRootr_pft              => plt_allom%rPCRootr_pft              ,& !input  :root P:C ratio, [gP gC-1]
    SeedAreaMean_pft          => plt_morph%SeedAreaMean_pft          ,& !input  :seed surface area, [m2]
    SeedCMass_pft             => plt_morph%SeedCMass_pft             ,& !input  :grain size at seeding, [gC/per seed]
    SeedWidth2LenRatio_pft    => plt_morph%SeedWidth2LenRatio_pft    ,& !input  :Seed width to length ratio, assuming prolate spheroid
    SeedMeanLen_pft           => plt_morph%SeedMeanLen_pft           ,& !input  :seed length, [m]
    SeedVolumeMean_pft        => plt_morph%SeedVolumeMean_pft        ,& !input  :seed volume, [m3 ]
    CMinNH4Root_pft           => plt_rbgc%CMinNH4Root_pft            ,& !inoput :minimum NH4 concentration for root NH4 uptake, [g m-3]
    CMinPO4Root_pft           => plt_rbgc%CMinPO4Root_pft            ,& !inoput :minimum PO4 concentration for root NH4 uptake, [g m-3]
    CminNO3Root_pft           => plt_rbgc%CminNO3Root_pft            ,& !inoput :minimum NO3 concentration for root NH4 uptake, [g m-3]
    KmNH4Root_pft             => plt_rbgc%KmNH4Root_pft              ,& !inoput :Km for root NH4 uptake, [g m-3]
    KmNO3Root_pft             => plt_rbgc%KmNO3Root_pft              ,& !inoput :Km for root NO3 uptake, [g m-3]
    KmPO4Root_pft             => plt_rbgc%KmPO4Root_pft              ,& !inoput :Km for root PO4 uptake, [g m-3]
    RootPorosity_pft          => plt_morph%RootPorosity_pft          ,& !inoput :root porosity, [m3 m-3]
    VmaxNH4Root_pft           => plt_rbgc%VmaxNH4Root_pft            ,& !inoput :maximum root NH4 uptake rate, [g m-2 h-1]
    VmaxNO3Root_pft           => plt_rbgc%VmaxNO3Root_pft            ,& !inoput :maximum root NO3 uptake rate, [g m-2 h-1]
    VmaxPO4Root_pft           => plt_rbgc%VmaxPO4Root_pft            ,& !inoput :maximum root PO4 uptake rate, [g m-2 h-1]
    Root1stMaxRadius1_pft     => plt_morph%Root1stMaxRadius1_pft     ,& !output :root diameter primary axes, [m]
    Root1stMaxRadius_pft      => plt_morph%Root1stMaxRadius_pft      ,& !output :maximum radius of primary roots, [m]
    Root2ndMaxRadius1_pft     => plt_morph%Root2ndMaxRadius1_pft     ,& !output :root diameter secondary axes, [m]
    Root2ndMaxRadius_pft      => plt_morph%Root2ndMaxRadius_pft      ,& !output :maximum radius of secondary roots, [m]
    FineRootVolPerMassC_pft   => plt_morph%FineRootVolPerMassC_pft   ,& !output :fine root volume:mass ratio, [m3 g-1]
    SeedDepth_pft             => plt_morph%SeedDepth_pft             ,& !output :seeding depth, [m]
    SoilSurfDepZ_col          => plt_site%SoilSurfDepZ_col           ,& !input  :soil surface depth, [m]
    CNRTS_pft                 => plt_allom%CNRTS_pft                 ,& !output :root N:C ratio x root growth yield, [-]
    CPRTS_pft                 => plt_allom%CPRTS_pft                 ,& !output :root P:C ratio x root growth yield, [-]
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft        ,& !output :soil layer at planting depth, [-]
    NMaxRootBotLayer_pft      => plt_morph%NMaxRootBotLayer_pft      ,& !output :maximum soil layer number for all root axes, [-]
    NRoot1stTipLay_raxes      => plt_morph%NRoot1stTipLay_raxes      ,& !output :maximum soil layer number for root axes, [-]
    Root1stSpecLen_pft        => plt_morph%Root1stSpecLen_pft        ,& !output :specific root length primary axes, [m g-1]
    Root1stXSecArea_pft       => plt_morph%Root1stXSecArea_pft       ,& !output :root cross-sectional area primary axes, [m2]
    Root2ndSpecLen_pft        => plt_morph%Root2ndSpecLen_pft        ,& !output :specific root length secondary axes, [m g-1]
    Root2ndXSecArea_pft       => plt_morph%Root2ndXSecArea_pft       ,& !output :root cross-sectional area secondary axes, [m2]
    Root2ndAxialResist_pft       => plt_morph%Root2ndAxialResist_pft       ,& !output :root axial resistivity, [MPa h m-4]
    RootPoreTortu4Gas_pft     => plt_morph%RootPoreTortu4Gas_pft     ,& !output :power function of root porosity used to calculate root gaseous diffusivity, [-]
    RootRadialResist_pft      => plt_morph%RootRadialResist_pft      ,& !output :root radial resistivity, [MPa h m-2]
    RootRaidus_rpft           => plt_morph%RootRaidus_rpft            & !output :root internal radius, [m]
  )
!
!     SEED CHARACTERISTICS
!
!     SeedVolumeMean_pft,SeedMeanLen_pft,SeedAreaMean_pft=seed volume(m3),length(m),AREA3(NU)(m2)
!     SeedCMass=seed C mass (g) from PFT file
!
  call calc_seed_geometry(SeedCMass_pft(NZ),SeedWidth2LenRatio_pft(NZ),SeedVolumeMean_pft(NZ),SeedMeanLen_pft(NZ),SeedAreaMean_pft(NZ))
 
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) DIMENSIONS, UPTAKE PARAMETERS
!
!     SeedDepth_pft=seeding depth(m) from PFT management file
!     CumSoilThickness_vr=depth to soil layer bottom from surface(m)
!     NG,NIX,NRoot1stTipLay_raxes=seeding,upper,lower rooting layer
!     CNRTS_pft,CPRTS_pft=N,P root growth yield
!     Root1stMaxRadius_pft,Root2ndMaxRadius_pft=maximum primary,secondary mycorrhizal radius (m)
!     PORT=mycorrhizal porosity
!     VmaxNH4Root_pft,KmNH4Root_pft,CMinNH4Root_pft=NH4 max uptake(g m-2 h-1),Km(uM),min concn (uM)
!     VmaxNO3Root_pft,KmNO3Root_pft,CminNO3Root_pft=NO3 max uptake(g m-2 h-1),Km(uM), min concn (uM)
!     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake(g m-2 h-1),Km(uM),min concn (uM)
!     RootRadialResist_pft,Root2ndAxialResist_pft=radial,axial root resistivity (m2 MPa-1 h-1)
!
  SeedDepth_pft(NZ)=PlantinDepz_pft(NZ)+SoilSurfDepZ_col
  DO L=NU,NL
    if(isclose(SeedDepth_pft(NZ),CumSoilThickness_vr(L)))then
      SeedDepth_pft(NZ)=AMAX1(CumSoilThickness_vr(L)-ppmc,ppmc)
      exit
    endif
  ENDDO

  D9795: DO L=NU,NL
    IF(isALeftinBC(SeedDepth_pft(NZ),CumSoilThickness_vr(L-1),CumSoilThickness_vr(L)))THEN
      NGTopRootLayer_pft(NZ)   = L
      NMaxRootBotLayer_pft(NZ) = L
      D9790: DO NR=1,pltpar%MaxNumRootAxes
        NRoot1stTipLay_raxes(NR,NZ)=L
      ENDDO D9790
    ENDIF
  ENDDO D9795
  CNRTS_pft(NZ)              = rNCRoot_pft(NZ)*RootBiomGrosYld_pft(NZ)
  CPRTS_pft(NZ)              = rPCRootr_pft(NZ)*RootBiomGrosYld_pft(NZ)
  Root1stMaxRadius_pft(2,NZ) = 5.0E-06_r8      !mycorrhizae
  Root2ndMaxRadius_pft(2,NZ) = 5.0E-06_r8      !mycorrhizae
  RootPorosity_pft(2,NZ)     = RootPorosity_pft(ipltroot,NZ)
  VmaxNH4Root_pft(2,NZ)      = VmaxNH4Root_pft(ipltroot,NZ)
  KmNH4Root_pft(2,NZ)        = KmNH4Root_pft(ipltroot,NZ)
  CMinNH4Root_pft(2,NZ)      = CMinNH4Root_pft(ipltroot,NZ)
  VmaxNO3Root_pft(2,NZ)      = VmaxNO3Root_pft(ipltroot,NZ)
  KmNO3Root_pft(2,NZ)        = KmNO3Root_pft(ipltroot,NZ)
  CminNO3Root_pft(2,NZ)      = CminNO3Root_pft(ipltroot,NZ)
  VmaxPO4Root_pft(2,NZ)      = VmaxPO4Root_pft(ipltroot,NZ)
  KmPO4Root_pft(2,NZ)        = KmPO4Root_pft(ipltroot,NZ)
  CMinPO4Root_pft(2,NZ)      = CMinPO4Root_pft(ipltroot,NZ)
  RootRadialResist_pft(2,NZ) = 1.0E+04_r8            !mycorrhizae
  Root2ndAxialResist_pft(2,NZ)  = 1.0E+12_r8            !mycorrhizae
!
!     RootPoreTortu4Gas_pft=tortuosity for gas transport
!     RootRaidus_rpft=path length for radial diffusion within root (m)
!     FineRootVolPerMassC_pft=volume:C ratio (m3 g-1)
!     Root1stSpecLen_pft,Root2ndSpecLen_pft=specific primary,secondary root length (m g-1)
!     Root1stXSecArea_pft,Root2ndXSecArea_pft=specific primary,secondary root area (m2 g-1)
!

  D500: DO N=1,2
    RootPoreTortu4Gas_pft(N,NZ)   = RootPorosity_pft(N,NZ)**1.33_r8
    !Lemon, 1962, Soil Aeration and Plant Root Relations I. Theory, Figure 1.
    RootRaidus_rpft(N,NZ)         = LOG(1.0_r8/SQRT(AMAX1(0.01_r8,RootPorosity_pft(N,NZ))))
    FineRootVolPerMassC_pft(N,NZ) = 1.e-6_r8/(BlkDensFineRoots*(1.0_r8-RootPorosity_pft(N,NZ)))    !0.05gC cm-3, ~ 0.10g cm-3.
    Root1stSpecLen_pft(N,NZ)      = FineRootVolPerMassC_pft(N,NZ)/(PICON*Root1stMaxRadius_pft(N,NZ)**2)
    Root2ndSpecLen_pft(N,NZ)      = FineRootVolPerMassC_pft(N,NZ)/(PICON*Root2ndMaxRadius_pft(N,NZ)**2)
    Root1stMaxRadius1_pft(N,NZ)   = Root1stMaxRadius_pft(N,NZ)
!    2*SQRT(0.25*(1.0-RootPorosity_pft(N,NZ)))
    Root2ndMaxRadius1_pft(N,NZ)=Root2ndMaxRadius_pft(N,NZ)
!    2*SQRT(0.25*(1.0-RootPorosity_pft(N,NZ)))
    Root1stXSecArea_pft(N,NZ)=PICON*Root1stMaxRadius1_pft(N,NZ)**2
    Root2ndXSecArea_pft(N,NZ)=PICON*Root2ndMaxRadius1_pft(N,NZ)**2
  ENDDO D500
  end associate
  end subroutine InitDimensionsandUptake

!----------------------------------------------------------------------------------------------------
  subroutine InitPlantPhenoMorphoBio(NZ)

  implicit none
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='InitPlantPhenoMorphoBio'

  integer :: K,L,M,N,NB
  real(r8) :: NumVessels
  associate(                                                                           &
    AREA3                             => plt_site%AREA3                               ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    AlphaVesselAxialResist_pft        => plt_morph%AlphaVesselAxialResist_pft         ,& !input  :axial resistance corrector for a single 1 m water transport vessel, [-]
    RootVesselRadius_pft              => plt_morph%RootVesselRadius_pft               ,& !input  :typical radius of the water transport vessel in primary roots, [m]
    Root2ndMaxRadius_pft              => plt_morph%Root2ndMaxRadius_pft               ,& !input  :maximum radius of secondary roots, [m]
    RootSingleVesselRstaxial_pft      => plt_morph%RootSingleVesselRstaxial_pft       ,& !output :axial resistance for a single 1 m water transport vessel, [MPa h m-4]
    RootSingleVesselArea_pft          => plt_morph%RootSingleVesselArea_pft           ,& !output :cross-sectional area of a single water transport vessel, [m2]
    Root2ndAxialResist_pft            => plt_morph%Root2ndAxialResist_pft             ,& !output :root axial resistivity, [MPa h m-4]
    MatureGroup_pft                   => plt_pheno%MatureGroup_pft                    ,& !input  :acclimated plant maturity group, [-]
    NU                                => plt_site%NU                                  ,& !input  :current soil surface layer number, [-]
    PPX_pft                           => plt_site%PPX_pft                             ,& !input  :plant population, [plants m-2]
    PetolShethChemElmRemobFlx_brch    => plt_pheno%PetolShethChemElmRemobFlx_brch     ,& !input  :cached remobilizable sheath element mass, [g d-2]
    ShootNodeNumAtPlanting_pft        => plt_morph%ShootNodeNumAtPlanting_pft         ,& !input  :number of nodes in seed, [-]
    isPlantBranchAlive_brch           => plt_pheno%isPlantBranchAlive_brch            ,& !input  :flag to detect branch death, [-]
    Hours4LenthenPhotoPeriod_brch     => plt_pheno%Hours4LenthenPhotoPeriod_brch      ,& !output :initial heat requirement for spring leafout/dehardening, [h]
    Hours4ShortenPhotoPeriod_brch     => plt_pheno%Hours4ShortenPhotoPeriod_brch      ,& !output :initial cold requirement for autumn leafoff/hardening, [h]
    ShootNodeNum_brch                 => plt_morph%ShootNodeNum_brch                  ,& !output :shoot node number, [-]
    BranchNumerID_brch                => plt_morph%BranchNumerID_brch                 ,& !output :branch meric id, [-]
    BranchNumber_pft                  => plt_morph%BranchNumber_pft                   ,& !output :main branch numeric id,[-]
    GrainFillDowreg_brch              => plt_photo%GrainFillDowreg_brch               ,& !output :down-regulation of C4 photosynthesis, [-]
    CMassCO2BundleSheath_node         => plt_photo%CMassCO2BundleSheath_node          ,& !output :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CMassHCO3BundleSheath_node        => plt_photo%CMassHCO3BundleSheath_node         ,& !output :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CPOOL3_node                       => plt_photo%CPOOL3_node                        ,& !output :minimum sink strength for nonstructural C transfer, [g d-2]
    CPOOL4_node                       => plt_photo%CPOOL4_node                        ,& !output :leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
    CanPBranchHeight                  => plt_morph%CanPBranchHeight                   ,& !output :branch height, [m]
    CanopyLeafAreaMAX_pft             => plt_morph%CanopyLeafAreaMAX_pft              ,& !output :running maximum leaf area, [m2 d-2]    
    CanopyHeightLive_pft              => plt_morph%CanopyHeightLive_pft               ,& !output :canopy height, [m]
    CanopyLeafAreaZ_pft               => plt_morph%CanopyLeafAreaZ_pft                ,& !output :canopy layer leaf area, [m2 d-2]
    CanopyLeafArea_lnode              => plt_morph%CanopyLeafArea_lnode               ,& !output :layer/node/branch leaf area, [m2 d-2]
    CanopyLeafArea_pft                => plt_morph%CanopyLeafArea_pft                 ,& !output :plant canopy leaf area, [m2 d-2]
    CanopyStalkSurfArea_lbrch         => plt_morph%CanopyStalkSurfArea_lbrch          ,& !output :plant canopy layer branch stem area, [m2 d-2]
    CanopyStemSurfAreaZ_pft           => plt_morph%CanopyStemSurfAreaZ_pft            ,& !output :plant canopy layer stem area, [m2 d-2]
    CanopyStemSurfArea_pft            => plt_morph%CanopyStemSurfArea_pft             ,& !output :plant stem area, [m2 d-2]
    ChillHours_pft                    => plt_photo%ChillHours_pft                     ,& !output :chilling effect on CO2 fixation, [-]
    HourFailGrainFill_brch            => plt_pheno%HourFailGrainFill_brch             ,& !output :flag to detect physiological maturity from grain fill, [-]
    Hours2LeafOut_brch                => plt_pheno%Hours2LeafOut_brch                 ,& !output :counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
    Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch                 ,& !output :cold requirement for autumn leafoff/hardening, [h]
    Hours4Leafout_brch                => plt_pheno%Hours4Leafout_brch                 ,& !output :heat requirement for spring leafout/dehardening, [h]
    HoursDoingRemob_brch              => plt_pheno%HoursDoingRemob_brch               ,& !output :counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
    HoursTooLowPsiCan_pft             => plt_pheno%HoursTooLowPsiCan_pft              ,& !output :canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY), [h]
    HypocotHeight_pft                 => plt_morph%HypocotHeight_pft                  ,& !output :cotyledon height, [m]
    KHiestGroLeafNode_brch            => plt_pheno%KHiestGroLeafNode_brch             ,& !output :leaf growth stage counter, [-]
    KLeafNumber_brch                  => plt_morph%KLeafNumber_brch                   ,& !output :leaf number, [-]
    KLowestGroLeafNode_brch           => plt_pheno%KLowestGroLeafNode_brch            ,& !output :leaf growth stage counter, [-]
    KMinNumLeaf4GroAlloc_brch         => plt_morph%KMinNumLeaf4GroAlloc_brch          ,& !output :NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION,[-]
    LeafAreaDying_brch                => plt_morph%LeafAreaDying_brch                 ,& !output :branch leaf area, [m2 d-2]
    LeafAreaLive_brch                 => plt_morph%LeafAreaLive_brch                  ,& !output :branch leaf area, [m2 d-2]
    LeafAreaZsec_brch                 => plt_morph%LeafAreaZsec_brch                  ,& !output :leaf surface area, [m2 d-2]
    LeafArea_node                     => plt_morph%LeafArea_node                      ,& !output :leaf area, [m2 d-2]
    LeafNumberAtFloralInit_brch       => plt_pheno%LeafNumberAtFloralInit_brch        ,& !output :leaf number at floral initiation, [-]
    StalkNodeHeight_brch              => plt_morph%StalkNodeHeight_brch               ,& !output :internode height, [m]
    MatureGroup_brch                  => plt_pheno%MatureGroup_brch                   ,& !output :plant maturity group, [-]
    ShootNodeNumAtInitFloral_brch     => plt_morph%ShootNodeNumAtInitFloral_brch      ,& !output :shoot node number at floral initiation, [-]
    NodeNumNormByMatgrp_brch          => plt_pheno%NodeNumNormByMatgrp_brch           ,& !output :normalized node number during vegetative growth stages, [-]
    ShootNodeNumAtAnthesis_brch       => plt_morph%ShootNodeNumAtAnthesis_brch        ,& !output :shoot node number at anthesis, [-]
    NumOfBranches_pft                 => plt_morph%NumOfBranches_pft                  ,& !output :number of branches,[-]
    NumOfLeaves_brch                  => plt_morph%NumOfLeaves_brch                   ,& !output :leaf number, [-]
    PlantPopuLive_pft                 => plt_site%PlantPopuLive_pft                   ,& !output :plant population, [d-2]
    PotentialSeedSites_brch           => plt_morph%PotentialSeedSites_brch            ,& !output :branch potential grain number, [d-2]
    ReprodNodeNumNormByMatrgrp_brch   => plt_pheno%ReprodNodeNumNormByMatrgrp_brch    ,& !output :normalized node number during reproductive growth stages, [-]
    RubiscoActivity_brch              => plt_photo%RubiscoActivity_brch               ,& !output :branch down-regulation of CO2 fixation, [-]
    SetNumberSeeds_brch               => plt_morph%SetNumberSeeds_brch                ,& !output :branch grain number, [d-2]
    StemAreaZsec_brch                 => plt_morph%StemAreaZsec_brch                  ,& !output :stem surface area, [m2 d-2]
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch  ,& !output :normalized node number during reproductive growth stages, [-]
    dReproNodeNumNormByMatG_brch      => plt_pheno%dReproNodeNumNormByMatG_brch ,& !output :current hourly reproductive development increment, [h-1]
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch      ,& !output :normalized node number during vegetative growth stages, [-]
    Cytokinin2ndConc_rpvr             => plt_rbgc%Cytokinin2ndConc_rpvr               ,& !output :cytokinin concentration in fine roots, [gC m-3 H2O]    
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch                 & !output :plant growth stage, [-]
  )
  !
  !     INITIALIZE PLANT PHENOLOGY
  !
  !     PP=population (grid cell-1)
  !
  call PrintInfo('beg '//subname)
  CanopyLeafAreaMAX_pft(NZ) = 0._R8
  plt_pheno%Days4FalseBreak_pft(NZ)  = 0
  plt_pheno%doInitPlant_pft(NZ)      = ifalse
  plt_pheno%isPlantShootAlive_pft(NZ) = iTrue
  plt_pheno%isPlantRootAlive_pft(NZ)  = iTrue
  BranchNumber_pft(NZ)               = 0
  NumOfBranches_pft(NZ)              = 0
  HypocotHeight_pft(NZ)              = 0._r8
  CanopyHeightLive_pft(NZ)               = 0._r8
  Cytokinin2ndConc_rpvr(:,:,:,NZ)    = 0._r8
  !recompute vessel axial resistances from the current PFT traits, which may change between plantings
  RootSingleVesselRstaxial_pft(NZ)     = AlphaVesselAxialResist_pft(NZ)*Rax_ref*(1.e-6_r8/RootVesselRadius_pft(NZ))**4
  RootSingleVesselArea_pft(NZ)         = PICON*RootVesselRadius_pft(NZ)**2
  NumVessels                           = 0.2_r8*PICON*Root2ndMaxRadius_pft(ipltroot,NZ)**2/RootSingleVesselArea_pft(NZ)
  Root2ndAxialResist_pft(ipltroot,NZ)  = RootSingleVesselRstaxial_pft(NZ)/NumVessels
  D10: DO NB=1,MaxNumBranches
    plt_pheno%doInitLeafOut_brch(NB,NZ)           = iTrue
    plt_pheno%EnablePlantLeafOut_brch(NB,NZ)      = iTrue
    plt_pheno%doPlantLeaveOff_brch(NB,NZ)         = iTrue
    plt_pheno%Prep4Literfall_brch(NB,NZ)          = ifalse
    plt_pheno%Hours4LiterfalAftMature_brch(NB,NZ) = 0
    MatureGroup_brch(NB,NZ)                       = MatureGroup_pft(NZ)
    ShootNodeNum_brch(NB,NZ)                      = ShootNodeNumAtPlanting_pft(NZ)
    ShootNodeNumAtInitFloral_brch(NB,NZ)          = ShootNodeNum_brch(NB,NZ)
    ShootNodeNumAtAnthesis_brch(NB,NZ)              = 0._r8
    NumOfLeaves_brch(NB,NZ)                       = 0._r8
    LeafNumberAtFloralInit_brch(NB,NZ)            = 0._r8
    KLeafNumber_brch(NB,NZ)                       = 1
    KMinNumLeaf4GroAlloc_brch(NB,NZ)              = 1
    KHiestGroLeafNode_brch(NB,NZ)                 = 1
    KLowestGroLeafNode_brch(NB,NZ)                = 0
    NodeNumNormByMatgrp_brch(NB,NZ)               = 0._r8
    ReprodNodeNumNormByMatrgrp_brch(NB,NZ)        = 0._r8
    TotalNodeNumNormByMatgrp_brch(NB,NZ)          = 0._r8
    TotReproNodeNumNormByMatrgrp_brch(NB,NZ)      = 0._r8
    dReproNodeNumNormByMatG_brch(NB,NZ)           = 0._r8
    Hours4LenthenPhotoPeriod_brch(NB,NZ)          = 0._r8
    Hours4ShortenPhotoPeriod_brch(NB,NZ)          = 0._r8
    Hours4Leafout_brch(NB,NZ)                     = Hours4LenthenPhotoPeriod_brch(NB,NZ)
    Hours4LeafOff_brch(NB,NZ)                     = Hours4ShortenPhotoPeriod_brch(NB,NZ)
    Hours2LeafOut_brch(NB,NZ)               = 0._r8
    RubiscoActivity_brch(NB,NZ)             = 1.0_r8
    GrainFillDowreg_brch(NB,NZ)             = 1.0_r8
    HourFailGrainFill_brch(NB,NZ)           = 0
    HoursDoingRemob_brch(NB,NZ)             = 0
    BranchNumerID_brch(NB,NZ)               = 0
    plt_pheno%isPlantBranchAlive_brch(NB,NZ) = iFalse
    D15: DO M=1,pltpar%NumGrowthStages
      iPlantCalendar_brch(M,NB,NZ)=0
    ENDDO D15
  ENDDO D10
  !
  !     INITIALIZE PLANT MORPHOLOGY AND BIOMASS
  !
  HoursTooLowPsiCan_pft(NZ)=0._r8
  ChillHours_pft(NZ)=0._r8
  plt_biom%CanopyNonstElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)            = 0._r8
  plt_biom%CanopyNodulNonstElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)       = 0._r8
  plt_biom%PetolShethStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)        = 0._r8
  plt_biom%ShootElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)                  = 0._r8
  plt_biom%StalkStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)             = 0._r8
  plt_biom%LeafStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)              = 0._r8
  plt_biom%StalkRsrvElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)              = 0._r8
  plt_biom%HuskStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)              = 0._r8
  plt_biom%GrainStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)             = 0._r8
  plt_biom%EarStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)               = 0._r8
  plt_biom%CanopyNodulStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)       = 0._r8
  plt_pheno%LeafElmntRemobFlx_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)         = 0._r8
  plt_pheno%LeafSenescInitialElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)         = -1._r8
  plt_pheno%PetolShethChemElmRemobFlx_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ) = 0._r8
  plt_pheno%PetolSenescInitialElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ) = -1._r8
  plt_biom%SenecStalkStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)        = 0._r8
  plt_morph%NActiveRootSegs_raxes(:,NZ) = 0
  plt_morph%IndRootSegBase_raxes(:,NZ)   = 1
  plt_morph%IndRootSegTip_raxes(:,NZ)    = 1
  plt_morph%RootSeglengths_raxes(:,:,NZ) = 0.0_r8
  plt_morph%RootSegAges_raxes(:,:,NZ)    = 0.0_r8  
  D25: DO NB=1,MaxNumBranches
    plt_biom%SapwoodBiomassC_brch(NB,NZ)     = 0._r8
    plt_biom%CanopyLeafSheathC_brch(NB,NZ)   = 0._r8
    PotentialSeedSites_brch(NB,NZ)           = 0._r8
    SetNumberSeeds_brch(NB,NZ)                 = 0._r8
    plt_allom%SingleGrainMeanBiomC_brch(NB,NZ) = 0._r8
    LeafAreaLive_brch(NB,NZ)                 = 0._r8
    plt_rbgc%NH3Dep2Can_brch(NB,NZ)          = 0._r8
    LeafAreaDying_brch(NB,NZ)                = 0._r8
    CanPBranchHeight(NB,NZ)                  = 0._r8
    D5: DO L=1,NumCanopyLayers1
      CanopyStalkSurfArea_lbrch(L,NB,NZ)=0._r8
      DO N=1,NumLeafInclinationClasses1
        StemAreaZsec_brch(N,L,NB,NZ)=0._r8
      enddo
    ENDDO D5

    DO K=0,MaxNodesPerBranch1
      LeafArea_node(K,NB,NZ)                                        = 0._r8
      StalkNodeHeight_brch(K,NB,NZ)                                 = 0._r8
      plt_morph%StalkNodeVertLength_brch(K,NB,NZ)                   = 0._r8
      plt_morph%PetoleLength_node(K,NB,NZ)                          = 0._r8
      plt_biom%LeafProteinC_node(K,NB,NZ)                           = 0._r8
      plt_biom%PetoleProteinC_node(K,NB,NZ)                         = 0._r8
      plt_biom%LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)       = 0._r8
      plt_biom%StructInternodeElms_brch(1:NumPlantChemElms,K,NB,NZ) = 0._r8
      plt_biom%PetolShethElmntNode_brch(1:NumPlantChemElms,K,NB,NZ) = 0._r8

      D55: DO L=1,NumCanopyLayers1
        CanopyLeafArea_lnode(L,K,NB,NZ)=0._r8
        plt_biom%LeafLayerElms_node(1:NumPlantChemElms,L,K,NB,NZ)=0._r8
      ENDDO D55

      IF(K.NE.0)THEN
        CPOOL3_node(K,NB,NZ)                = 0._r8
        CMassCO2BundleSheath_node(K,NB,NZ)  = 0._r8
        CMassHCO3BundleSheath_node(K,NB,NZ) = 0._r8
        CPOOL4_node(K,NB,NZ)                = 0._r8
        D45: DO L=1,NumCanopyLayers1
          DO N=1,NumLeafInclinationClasses1
            LeafAreaZsec_brch(N,L,K,NB,NZ)=0._r8
          enddo
        ENDDO D45
      ENDIF
    enddo
  ENDDO D25
  
  D35: DO L=1,NumCanopyLayers1
    CanopyLeafAreaZ_pft(L,NZ)         = 0._r8
    plt_biom%CanopyLeafCLyr_pft(L,NZ) = 0._r8
    CanopyStemSurfAreaZ_pft(L,NZ)         = 0._r8
  ENDDO D35
  plt_biom%CanopyNonstElms_pft(1:NumPlantChemElms,NZ)     = 0._r8
  plt_biom%CanopyNonstElmConc_pft(1:NumPlantChemElms,NZ)  = 0._r8
  plt_biom%CanopyNoduleNonstCConc_pft(NZ)                 = 0._r8
  plt_biom%ShootElms_pft(1:NumPlantChemElms,NZ)           = 0._r8
  plt_biom%LeafStrutElms_pft(1:NumPlantChemElms,NZ)       = 0._r8
  plt_biom%PetolShethStrutElms_pft(1:NumPlantChemElms,NZ) = 0._r8
  plt_biom%StalkStrutElms_pft(1:NumPlantChemElms,NZ)      = 0._r8
  plt_biom%CanopySapwoodC_pft(NZ)                         = 0._r8
  plt_biom%StalkRsrvElms_pft(1:NumPlantChemElms,NZ)       = 0._r8
  plt_biom%HuskStrutElms_pft(1:NumPlantChemElms,NZ)       = 0._r8
  plt_biom%EarStrutElms_pft(1:NumPlantChemElms,NZ)        = 0._r8
  plt_biom%GrainStrutElms_pft(1:NumPlantChemElms,NZ)      = 0._r8
  plt_biom%RootStrutElms_pft(1:NumPlantChemElms,NZ)       = 0._r8
  plt_biom%CanopyLeafSheathC_pft(NZ)                      = 0._r8
  plt_biom%RootNoduleElms_pft(:,NZ)                       = 0._r8
  plt_biom%ShootNoduleElms_pft(:,NZ)                      = 0._r8
  CanopyLeafArea_pft(NZ)                                  = 0._r8
  plt_biom%RootBiomCPerPlant_pft(NZ)                      = 0._r8
  CanopyStemSurfArea_pft(NZ)                              = 0._r8
  call PrintInfo('end '//subname)
  end associate
  end subroutine InitPlantPhenoMorphoBio

!----------------------------------------------------------------------------------------------------
  subroutine InitMassBalance(NZ)

  implicit none
  integer, intent(in) :: NZ
  integer :: M,NE
  real(r8) :: StandDeadC_pft

  associate(                                                                    &
    AREA3                          => plt_site%AREA3                           ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    PlantElmAllocMat4Litr          => plt_soilchem%PlantElmAllocMat4Litr       ,& !input  :litter kinetic fraction, [-]
    NU                             => plt_site%NU                              ,& !input  :current soil surface layer number, [-]
    StandingDeadInitC_pft          => plt_biom%StandingDeadInitC_pft           ,& !input  :initial standing dead C, [g C m-2]
    icwood                         => pltpar%icwood                            ,& !input  :group id of coarse woody litter
    rNCStalk_pft                   => plt_allom%rNCStalk_pft                   ,& !input  :stalk N:C ratio, [gN gC-1]
    rPCStalk_pft                   => plt_allom%rPCStalk_pft                   ,& !input  :stalk P:C ratio, [g g-1]
    StandDeadStrutElms_pft         => plt_biom%StandDeadStrutElms_pft          ,& !inoput :standing dead element, [g d-2]
    StandDeadCompKElms_pft         => plt_biom%StandDeadCompKElms_pft          ,& !output :standing dead element fraction, [g d-2]
    CanopyRespC_CumYr_pft          => plt_bgcr%CanopyRespC_CumYr_pft           ,& !output :total autotrophic respiration, [gC d-2 ]
    ETCanopy_CumYr_pft             => plt_ew%ETCanopy_CumYr_pft                ,& !output :total transpiration, [m H2O d-2]
    GrossCO2Fix_pft                => plt_bgcr%GrossCO2Fix_pft                 ,& !output :total gross CO2 fixation, [gC d-2 ]
    GrossResp_pft                  => plt_bgcr%GrossResp_pft                   ,& !output :total plant respiration, [gC d-2 ]
    LitrfalStrutElms_CumYr_pft     => plt_bgcr%LitrfalStrutElms_CumYr_pft      ,& !output :total plant element LitrFall, [g d-2 ]
    NH3Dep2Can_pft                 => plt_bgcr%NH3Dep2Can_pft                  ,& !output :canopy NH3 flux, [g d-2 h-1]
    NH3Emis_CumYr_pft              => plt_bgcr%NH3Emis_CumYr_pft               ,& !output :total canopy NH3 flux, [gN d-2 ]
    NetCumElmntFlx2Plant_pft       => plt_pheno%NetCumElmntFlx2Plant_pft       ,& !output :effect of canopy element status on seed set, [-]
    PlantExudElm_CumYr_pft         => plt_rbgc%PlantExudElm_CumYr_pft          ,& !output :total net root element uptake (+ve) - exudation (-ve), [gC d-2 ]
    PlantN2Fix_CumYr_pft           => plt_bgcr%PlantN2Fix_CumYr_pft            ,& !output :total plant N2 fixation, [g d-2 ]
    RootUptk_N_CumYr_pft           => plt_rbgc%RootUptk_N_CumYr_pft            ,& !output :cumulative plant N uptake, [gN d-2]
    RootUptk_P_CumYr_pft           => plt_rbgc%RootUptk_P_CumYr_pft            ,& !output :cumulative plant P uptake, [gP d-2]
    SurfLitrfalStrutElms_CumYr_pft => plt_bgcr%SurfLitrfalStrutElms_CumYr_pft   & !output :total surface LitrFall element, [g d-2]
  )
!
!     INITIALIZE MASS BALANCE CHECKS
!  
  IF(is_cold_run().AND.is_first_year)THEN
    GrossCO2Fix_pft(NZ)                                   = 0._r8
    SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ) = 0._r8
    GrossResp_pft(NZ)                                     = 0._r8
    CanopyRespC_CumYr_pft(NZ)                             = 0._r8
    PlantExudElm_CumYr_pft(1:NumPlantChemElms,NZ)         = 0._r8
    RootUptk_N_CumYr_pft(NZ)                              = 0._r8
    RootUptk_P_CumYr_pft(NZ)                              = 0._r8
    LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ)     = 0._r8

    PlantN2Fix_CumYr_pft(NZ)                                 = 0._r8
    NH3Dep2Can_pft(NZ)                                       = 0._r8
    NH3Emis_CumYr_pft(NZ)                                    = 0._r8
    plt_distb%CO2ByFire_CumYr_pft(NZ)                        = 0._r8
    plt_distb%CH4ByFire_CumYr_pft(NZ)                        = 0._r8
    plt_distb%O2ByFire_CumYr_pft(NZ)                         = 0._r8
    plt_distb%NH3byFire_CumYr_pft(NZ)                        = 0._r8
    plt_distb%N2ObyFire_CumYr_pft(NZ)                        = 0._r8
    plt_distb%PO4byFire_CumYr_pft(NZ)                        = 0._r8
    plt_distb%EcoHavstElmntCum_pft(1:NumPlantChemElms,NZ)    = 0._r8
    plt_distb%EcoHavstElmnt_CumYr_pft(1:NumPlantChemElms,NZ) = 0._r8
    NetCumElmntFlx2Plant_pft(1:NumPlantChemElms,NZ)          = 0._r8
    ETCanopy_CumYr_pft(NZ)                                   = 0._r8
  ENDIF
  end associate
  end subroutine InitMassBalance

!----------------------------------------------------------------------------------------------------
  subroutine InitPlantHeatWater(NZ)

  implicit none
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='InitPlantHeatWater'
  REAL(R8) :: FDM

  associate(                                                    &
    ATCA                    => plt_site%ATCA                   ,& !input  :mean annual air temperature, [oC]
    OrganOsmoPsi0pt_pft     => plt_ew%OrganOsmoPsi0pt_pft      ,& !input  :Organ osmotic potential when canopy water potential = 0 MPa, [MPa]
    CanopyLeafSheathC_pft   => plt_biom%CanopyLeafSheathC_pft  ,& !input  :canopy leaf + sheath C, [g d-2]
    ShootElms_pft           => plt_biom%ShootElms_pft          ,& !input  :canopy shoot structural chemical element mass, [g d-2]
    StemSpecVolume_pft      => plt_morph%StemSpecVolume_pft    ,& !input  :stalk specific volume, [m3 gC-1]        
    CanopyHeatLoss2Dist_col     => plt_ew%CanopyHeatLoss2Dist_col      ,& !inoput :canopy energy +/- due to disturbance, [MJ /d2]
    QCanopyWatLoss2Dist_col     => plt_ew%QCanopyWatLoss2Dist_col      ,& !inoput :canopy water +/- due to disturbance, [m3 H2O/d2]
    CanopyBiomWater_pft     => plt_ew%CanopyBiomWater_pft      ,& !output :canopy water content, [m3 d-2]
    PSICanopyOsmo_pft       => plt_ew%PSICanopyOsmo_pft        ,& !output :canopy osmotic water potential, [Mpa]
    PSICanopy_pft           => plt_ew%PSICanopy_pft            ,& !output :canopy total water potential, [Mpa]
    TCGroth_pft             => plt_pheno%TCGroth_pft           ,& !output :canopy growth temperature, [oC]
    TKC_pft                 => plt_ew%TKC_pft                  ,& !output :canopy temperature, [K]
    TdegCCanopy_pft         => plt_ew%TdegCCanopy_pft          ,& !output :canopy temperature, [oC]
    VHeatCapCanopy_pft      => plt_ew%VHeatCapCanopy_pft       ,& !output :canopy heat capacity, [MJ d-2 K-1]
    DeltaTKC_pft            => plt_ew%DeltaTKC_pft             ,& !output :change in canopy temperature, [K]
    ENGYX_pft               => plt_ew%ENGYX_pft                ,& !output :canopy heat storage from previous time step, [MJ d-2]
    FracPARads2Canopy_pft     => plt_rad%FracPARads2Canopy_pft       ,& !output :fraction of incoming PAR absorbed by total canopy, [-]
    FracPARads2LiveCanopy_pft => plt_rad%FracPARads2LiveCanopy_pft   ,& !output :fraction of incoming PAR absorbed by live canopy, [-]
    PSICanopyTurg_pft       => plt_ew%PSICanopyTurg_pft        ,& !output :plant canopy turgor water potential, [MPa]
    TKGroth_pft             => plt_pheno%TKGroth_pft           ,& !output :canopy growth temperature, [K]
    Transpiration_pft       => plt_ew%Transpiration_pft        ,& !output :canopy transpiration, [m2 d-2 h-1]
    fTCanopyGroth_pft       => plt_pheno%fTCanopyGroth_pft      & !output :canopy temperature growth function, [-]
  )
  call PrintInfo('beg '//subname)
  !
  !     INITIALIZE PLANT HEAT AND WATER STATUS
  !
  !     VHeatCapCanopy_pft=canopy heat capacity (MJ m-3 K-1)
  !     TdegCCanopy_pft,TKC=canopy temperature for growth (oC,K)
  !     TCGroth_pft,TKGroth_pft=canopy temperature for phenology (oC,K)
  !     PSICanopy_pft,PSICanopyOsmo_pft,PSICanopyTurg_pft=canopy total,osmotic,turgor water potl(MPa)
  !
  VHeatCapCanopy_pft(NZ)    = cpw*ShootElms_pft(ielmc,NZ)*10.0E-06
  ENGYX_pft(NZ)             = 0._r8
  DeltaTKC_pft(NZ)          = 0._r8
  TdegCCanopy_pft(NZ)       = ATCA
  TKC_pft(NZ)               = units%Celcius2Kelvin(TdegCCanopy_pft(NZ))
  TCGroth_pft(NZ)           = TdegCCanopy_pft(NZ)
  TKGroth_pft(NZ)           = units%Celcius2Kelvin(TCGroth_pft(NZ))
  fTCanopyGroth_pft(NZ)     = 1.0
  PSICanopy_pft(NZ)         = -1.0E-03
  PSICanopyOsmo_pft(NZ)     = OrganOsmoPsi0pt_pft(NZ)+PSICanopy_pft(NZ)
  PSICanopyTurg_pft(NZ)     = AZMAX1(PSICanopy_pft(NZ)-PSICanopyOsmo_pft(NZ))
  Transpiration_pft(NZ)     = 0._r8
  FracPARads2Canopy_pft(NZ)     = 0._r8
  FracPARads2LiveCanopy_pft(NZ) = 0._r8
  FDM                       = get_FDM(PSICanopy_pft(NZ))
  CanopyBiomWater_pft(NZ)   = ppmc*CanopyLeafSheathC_pft(NZ)/FDM
  VHeatCapCanopy_pft(NZ)    = cpw*(ShootElms_pft(ielmc,NZ)*StemSpecVolume_pft(NZ)+CanopyBiomWater_pft(NZ))
  QCanopyWatLoss2Dist_col       = QCanopyWatLoss2Dist_col-CanopyBiomWater_pft(NZ)
  CanopyHeatLoss2Dist_col       = CanopyHeatLoss2Dist_col-VHeatCapCanopy_pft(NZ)*TKC_pft(NZ)
  call PrintInfo('end '//subname)
  end associate
  end subroutine InitPlantHeatWater

!----------------------------------------------------------------------------------------------------
  subroutine InitRootMychorMorphoBio(NZ)
  implicit none
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='InitRootMychorMorphoBio'

  integer :: K,L,M,N,NR,idg
  REAL(R8) :: CCO2A
  REAL(R8) :: CCO2P
  REAL(R8) :: COXYA
  REAL(R8) :: COXYP
  associate(                                                                   &
    ATCA                          => plt_site%ATCA                            ,& !input  :mean annual air temperature, [oC]
    CCO2EI_gperm3                 => plt_site%CCO2EI_gperm3                   ,& !input  :initial atmospheric CO2 concentration, [g m-3]
    CO2EI                         => plt_site%CO2EI                           ,& !input  :initial atmospheric CO2 concentration, [umol mol-1]
    COXYE                         => plt_site%COXYE                           ,& !input  :current atmospheric O2 concentration, [g m-3]
    OrganOsmoPsi0pt_pft           => plt_ew%OrganOsmoPsi0pt_pft               ,& !input  :Organ osmotic potential when canopy water potential = 0 MPa, [MPa]
    NL                            => plt_site%NL                              ,& !input  :lowest soil layer number,[-]
    OXYE                          => plt_site%OXYE                            ,& !input  :atmospheric O2 concentration, [umol mol-1]
    Root1stMaxRadius_pft          => plt_morph%Root1stMaxRadius_pft           ,& !input  :maximum radius of primary roots, [m]
    Root2ndMaxRadius_pft          => plt_morph%Root2ndMaxRadius_pft           ,& !input  :maximum radius of secondary roots, [m]
    RootProteinCMax_pft           => plt_allom%RootProteinCMax_pft            ,& !input  :reference root protein N, [gN g-1]
    SeedDepth_pft                 => plt_morph%SeedDepth_pft                  ,& !input  :seeding depth, [m]
    trcg_rootml_pvr               => plt_rbgc%trcg_rootml_pvr                 ,& !inoput :root gas content, [g d-2]
    trcs_deadroot2soil_pvr        => plt_rbgc%trcs_deadroot2soil_pvr          ,& !inoput :gases released to soil upong dying roots, [g d-2 h-1]    
    trcs_rootml_pvr               => plt_rbgc%trcs_rootml_pvr                 ,& !inoput :root aqueous content, [g d-2]
    PSIRootOSMO_vr                => plt_ew%PSIRootOSMO_vr                    ,& !output :root osmotic water potential, [Mpa]
    PSIRoot_pvr                   => plt_ew%PSIRoot_pvr                       ,& !output :root total water potential, [Mpa]
    RootPoreVol_pvr               => plt_morph%RootPoreVol_pvr                ,& !output :root layer volume air, [m2 d-2]
    RootVH2O_pvr                  => plt_morph%RootVH2O_pvr                   ,& !output :root space volume occupied by water in each layer, [m2 d-2]   
    NumStructuralRootAxes_pft     => plt_morph%NumStructuralRootAxes_pft      ,& !output :number of structural root axes,[-]
    PSIRootTurg_vr                => plt_ew%PSIRootTurg_vr                    ,& !output :root turgor water potential, [Mpa]
    Root1stTransptArea_pvr        => plt_morph%Root1stTransptArea_pvr         ,& !output :transport area by 1st order root, [m2 d-2]         
    RootMedTransptArea_pvr        => plt_morph%RootMedTransptArea_pvr         ,& !output  :transport area by medisum size roots, [-]        
    Root1stRadius_pvr             => plt_morph%Root1stRadius_pvr              ,& !output :root layer diameter primary axes, [m]
    Root2ndRadius_rpvr            => plt_morph%Root2ndRadius_rpvr             ,& !output :root layer diameter secondary axes, [m]
    RootProteinC_pvr              => plt_biom%RootProteinC_pvr                ,& !output :root layer protein C, [gC d-2]
    RootProteinConc_rpvr          => plt_biom%RootProteinConc_rpvr             & !output :root layer protein C concentration, [g g-1]
  )
  call PrintInfo('beg '//subname)
  !
  !     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) MORPHOLOGY AND BIOMASS
  !
  !     PSIRoot_pvr,PSIRootOSMO_vr,PSIRootTurg_vr=root,myco total,osmotic,turgor water potl(MPa)
  !     CO2A,CO2P=root,myco gaseous,aqueous CO2 content (g)
  !     OXYA,OXYP=root,myco gaseous,aqueous O2 content (g)
  !

  NumStructuralRootAxes_pft(NZ)=0
  plt_rbgc%RootNH4Uptake_pft(NZ)   = 0._r8
  plt_rbgc%RootNO3Uptake_pft(NZ)   = 0._r8
  plt_rbgc%RootH2PO4Uptake_pft(NZ) = 0._r8
  plt_rbgc%RootHPO4Uptake_pft(NZ)  = 0._r8
  plt_rbgc%RootN2Fix_pft(NZ)       = 0._r8
  DO NR=1,MaxNumRootAxes
    plt_morph%RootSegBaseDepth_raxes(NR,NZ)                = SeedDepth_pft(NZ)
    plt_morph%Root1stDepz_raxes(NR,NZ)                     = SeedDepth_pft(NZ)
    plt_biom%RootMyco1stElm_raxs(1:NumPlantChemElms,NR,NZ) = 0._r8    
  ENDDO      

  D20: DO L=1,NL
    D40: DO N=1,pltpar%jroots
      plt_ew%RootH2OUptkStress_pvr(N,L,NZ)                        = 0._r8
      plt_ew%RPlantRootH2OUptk_pvr(N,L,NZ)                        = 0._r8
      PSIRoot_pvr(N,L,NZ)                                          = -0.01_r8
      PSIRootOSMO_vr(N,L,NZ)                                       = OrganOsmoPsi0pt_pft(NZ)+PSIRoot_pvr(N,L,NZ)
      PSIRootTurg_vr(N,L,NZ)                                       = AZMAX1(PSIRoot_pvr(N,L,NZ)-PSIRootOSMO_vr(N,L,NZ))
      plt_biom%RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ)   = 0._r8
      plt_biom%RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ) = 0._r8
      RootProteinConc_rpvr(N,L,NZ)                                  = RootProteinCMax_pft(NZ)
      plt_biom%RootMycoActiveBiomC_pvr(N,L,NZ)                     = 0._r8
      plt_biom%PopuRootMycoC_pvr(N,L,NZ)=0._r8
      RootProteinC_pvr(N,L,NZ)                                 = 0._r8

      plt_morph%Root2ndXNumL_rpvr(N,L,NZ)                      = 0._r8
      plt_morph%RootTotLenPerPlant_pvr(N,L,NZ)                 = 0._r8
      plt_morph%RootLenDensPerPlant_pvr(N,L,NZ)                = 0._r8
      RootPoreVol_pvr(N,L,NZ)                                 = 0._r8
      RootVH2O_pvr(N,L,NZ)                                     = 0._r8
      Root1stRadius_pvr(N,L,NZ)                                = Root1stMaxRadius_pft(N,NZ)
      Root2ndRadius_rpvr(N,L,NZ)                               = Root2ndMaxRadius_pft(N,NZ)
      Root1stTransptArea_pvr(N,L,NZ)                           = PICON*Root1stMaxRadius_pft(N,NZ)**2  
      RootMedTransptArea_pvr(N,L,NZ)                           = 0._r8    
      plt_morph%RootSAreaPerPlant_pvr(N,L,NZ)                  = 0._r8
      plt_morph%Root2ndEffLen4uptk_rpvr(N,L,NZ)                = 1.0E-03

      plt_rbgc%RootNutUptake_pvr(ids_NH4B:ids_nuts_end,N,L,NZ) = 0._r8
      plt_rbgc%RootO2Dmnd4Resp_pvr(N,L,NZ)                     = 0._r8
      plt_rbgc%RootNH4DmndSoil_pvr(N,L,NZ)                     = 0._r8
      plt_rbgc%RootNH4DmndBand_pvr(N,L,NZ)                     = 0._r8
      plt_rbgc%RootNO3DmndSoil_pvr(N,L,NZ)                     = 0._r8
      plt_rbgc%RootNO3DmndBand_pvr(N,L,NZ)                     = 0._r8
      plt_rbgc%RootH2PO4DmndSoil_pvr(N,L,NZ)                   = 0._r8
      plt_rbgc%RootH1PO4DmndSoil_pvr(N,L,NZ)                   = 0._r8
      plt_rbgc%RootH2PO4DmndBand_pvr(N,L,NZ)                   = 0._r8
      plt_rbgc%RootH1PO4DmndBand_pvr(N,L,NZ)                   = 0._r8
      
      DO idg=idg_beg,idg_NH3
        trcs_deadroot2soil_pvr(idg,L,NZ) = trcs_deadroot2soil_pvr(idg,L,NZ) + trcg_rootml_pvr(idg,N,L,NZ)
        trcs_deadroot2soil_pvr(idg,L,NZ) = trcs_deadroot2soil_pvr(idg,L,NZ) + trcs_rootml_pvr(idg,N,L,NZ)
        plt_rbgc%trcg_rootml_pvr(idg,N,L,NZ)         = 0._r8
        plt_rbgc%trcs_rootml_pvr(idg,N,L,NZ)         = 0._r8
      ENDDO
      CCO2A                                             = CCO2EI_gperm3
      CCO2P                                             = 0.030*EXP(-2.621_r8-0.0317_r8*ATCA)*CO2EI
      plt_rbgc%trcg_air2root_flx_pvr(idg_CO2,N,L,NZ)    = 0._r8
      plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_CO2,N,L,NZ) = 0._r8
      
      COXYA                                             = COXYE
      COXYP                                             = 0.032_r8*EXP(-6.175_r8-0.0211_r8*ATCA)*OXYE
      plt_rbgc%RAutoRootO2Limter_rpvr(N,L,NZ)=1.0
      D30: DO NR=1,MaxNumRootAxes
        plt_morph%Root2ndXNum_rpvr(N,L,NR,NZ)                            = 0._r8
        plt_morph%Root2ndLen_rpvr(N,L,NR,NZ)                              = 0._r8
        plt_biom%RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = 0._r8
      ENDDO D30
    ENDDO D40

    plt_morph%Root1stXNumL_pvr(L,NZ)                      = 0._r8  
    DO NR=1,MaxNumRootAxes
      plt_biom%RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,L,NR,NZ) = 0._r8
      plt_biom%Root1stActStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ) = 0._r8
      plt_biom%Root1stLigStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ) = 0._r8
      plt_morph%Root1stLenPP_rpvr(L,NR,NZ)                           = 0._r8
      plt_morph%RootAge_rpvr(L,NR,NZ)                                = 0._r8
    ENDDO  

    plt_biom%RootNodulNonstElms_rpvr(1:NumPlantChemElms,L,NZ) = 0._r8
    plt_biom%RootNodulStrutElms_rpvr(1:NumPlantChemElms,L,NZ) = 0._r8

  ENDDO D20

  plt_morph%RootLenDensPerPlant_pvr(1:2,NL+1:JZ1,NZ)=0._r8
  call PrintInfo('end '//subname)
  end associate
  end subroutine InitRootMychorMorphoBio
!----------------------------------------------------------------------------------------------------
  subroutine ApplySeedBank(NZ,I,RootC)
  implicit none
  integer, intent(in) :: NZ
  integer, intent(in) :: I
  real(r8),intent(in) :: RootC

  character(len=*), parameter :: subname='ApplySeedBank'
  integer :: NE
  associate(                                                       &
    rNCGrain_pft             => plt_allom%rNCGrain_pft            ,& !input  :grain N:C ratio, [g g-1]
    rPCGrain_pft             => plt_allom%rPCGrain_pft            ,& !input  :grain P:C ratio, [gP gP-1]
    AREA3                    => plt_site%AREA3                    ,& !input  :soil cross section area, [m2]
    NU                       => plt_site%NU                       ,& !input  :current soil surface layer number, [-]
    PPatSeeding_pft          => plt_site%PPatSeeding_pft          ,& !input  :plant population at seeding, [plants d-2]
    iYearCurrent             => plt_site%iYearCurrent             ,& !input  :current year,[-]
    iPlantPhenolPattern_pft  => plt_pheno%iPlantPhenolPattern_pft ,& !input  :plant growth habit: annual or perennial,[-]    
    PPX_pft                  => plt_site%PPX_pft                  ,& !output :plant population, [plants m-2]
    FireReSet_pft            => plt_pheno%FireReSet_pft           ,& !output :flag to skill startq for fire rejuvenation    
    PlantElmDistLoss_pft     => plt_distb%PlantElmDistLoss_pft    ,& !inoput :plant loss to disturbance,    [g d-2 h-1]       
    PlantPopuLive_pft        => plt_site%PlantPopuLive_pft        ,& !output :plant population, [d-2]
    PlantPopuDead_pft        => plt_site%PlantPopuDead_pft        ,& !output :live+standing dead plant population, [d-2]
    iYearPlanting_pft        => plt_distb%iYearPlanting_pft       ,& !input  :year of planting,[-]    
    iDayPlanting_pft         => plt_distb%iDayPlanting_pft        ,& !input  :day of planting,[-]    
    SeedCMass_pft            => plt_morph%SeedCMass_pft           ,& !input  :grain size at seeding, [g]
    Eco_NBP_CumYr_col        => plt_bgcr%Eco_NBP_CumYr_col        ,& !inoput :total NBP, [g d-2]    
    doInitPlant_pft          => plt_pheno%doInitPlant_pft         ,& !inoput :PFT initialization flag:0=no,1=yes,[-]    
    SeasonalNonstElms_pft    => plt_biom%SeasonalNonstElms_pft    ,& !output :plant stored nonstructural element at current step, [g d-2]
    IsPlantActive_pft        => plt_pheno%IsPlantActive_pft       ,& !output :flag for living pft, [-]    
    SeedPlantedElm_pft       => plt_biom%SeedPlantedElm_pft        & !output :plant stored nonstructural C at planting, [gC d-2]
  )
  call PrintInfo('beg '//subname)
  PPX_pft(NZ)                = PPatSeeding_pft(NZ)
  PlantPopuLive_pft(NZ)      = PPatSeeding_pft(NZ)*AREA3(NU)
  PlantPopuDead_pft(NZ)      = PlantPopuLive_pft(NZ)

  IsPlantActive_pft(NZ)      = iTrue
  SeedPlantedElm_pft(:,NZ)   = 0._r8
  SeedPlantedElm_pft(ielmc,NZ) = SeedCMass_pft(NZ)*PlantPopuLive_pft(NZ)

  if(SeedPlantedElm_pft(ielmc,NZ).LE.SeasonalNonstElms_pft(ielmc,NZ))THEN
    !reseed 
    SeedPlantedElm_pft(ielmn,NZ) = rNCGrain_pft(NZ)*SeedPlantedElm_pft(ielmc,NZ)
    SeedPlantedElm_pft(ielmp,NZ) = rPCGrain_pft(NZ)*SeedPlantedElm_pft(ielmc,NZ)
    DO NE=1,NumPlantChemElms
      PlantElmDistLoss_pft(NE,NZ)  = PlantElmDistLoss_pft(NE,NZ)+SeasonalNonstElms_pft(NE,NZ)
    ENDDO
    SeasonalNonstElms_pft(:,NZ)  = SeedPlantedElm_pft(:,NZ)
    
    Eco_NBP_CumYr_col = Eco_NBP_CumYr_col+SeedPlantedElm_pft(ielmc,NZ)
  else
    SeedPlantedElm_pft(:,NZ)=0._r8
  endif  
  !the following is only for natural plants
  IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)then  
    iDayPlanting_pft(NZ)                                  = -1E+06
    iYearPlanting_pft(NZ)                                 = -1E+06
    doInitPlant_pft(NZ)                                   = itrue  
  else
    iYearPlanting_pft(NZ)        = iYearCurrent+1
    iDayPlanting_pft(NZ)         = 1 
  endif 
  FireReSet_pft(NZ)= iTrue
  call InitPlantPhenoMorphoBio(NZ)
  
  call PrintInfo('end '//subname)
  end associate
  end subroutine ApplySeedBank
!----------------------------------------------------------------------------------------------------
  subroutine InitSeedMorphoBio(NZ)

  implicit none
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='InitSeedMorphoBio'

  associate(                                                                   &
    rNCGrain_pft                  => plt_allom%rNCGrain_pft                   ,& !input  :grain N:C ratio, [g g-1]
    rPCGrain_pft                  => plt_allom%rPCGrain_pft                   ,& !input  :grain P:C ratio, [gP gP-1]
    jHarvstType_pft               => plt_distb%jHarvstType_pft                ,& !input  :flag for stand replacing disturbance,[-]    
    NGTopRootLayer_pft            => plt_morph%NGTopRootLayer_pft             ,& !input  :soil layer at planting depth, [-]
    PetolShethStrutElms_brch      => plt_biom%PetolShethStrutElms_brch        ,& !input  :branch sheath structural element, [g d-2]
    PlantPopuLive_pft             => plt_site%PlantPopuLive_pft               ,& !input  :plant population, [d-2]
    PopuRootMycoC_pvr             => plt_biom% PopuRootMycoC_pvr              ,& !input  :root layer C, [gC d-2]
    RootProteinCMax_pft           => plt_allom%RootProteinCMax_pft            ,& !input  :reference root protein N, [gN g-1]
    RootMyco1stStrutElms_rpvr     => plt_biom%RootMyco1stStrutElms_rpvr       ,& !input  :root layer element primary axes, [g d-2]
    Root1stActStructElms_rpvr     => plt_biom%Root1stActStructElms_rpvr       ,& !inoput :root layer active zone element in primary axes, [g d-2]
    RootMycoActiveBiomC_pvr       => plt_biom%RootMycoActiveBiomC_pvr         ,& !input  :root layer structural C, [gC d-2]
    RootMycoNonstElms_rpvr        => plt_biom%RootMycoNonstElms_rpvr          ,& !input  :root layer nonstructural element, [g d-2]
    RootProteinC_pvr              => plt_biom%RootProteinC_pvr                ,& !input  :root layer protein C, [gC d-2]
    SeedCMass_pft                 => plt_morph%SeedCMass_pft                  ,& !input  :grain size at seeding, [g]
    CanopyLeafSheathC_pft         => plt_biom%CanopyLeafSheathC_pft           ,& !inoput :canopy leaf + sheath C, [g d-2]
    CanopyNonstElms_brch          => plt_biom%CanopyNonstElms_brch            ,& !inoput :branch nonstructural element, [g d-2]
    LeafStrutElms_brch            => plt_biom%LeafStrutElms_brch              ,& !inoput :branch leaf structural element mass, [g d-2]
    RootMyco1stElm_raxs           => plt_biom%RootMyco1stElm_raxs             ,& !inoput :root C primary axes, [g d-2]
    SeasonalNonstElms_pft         => plt_biom%SeasonalNonstElms_pft           ,& !inoput :plant stored nonstructural element at current step, [g d-2]
    CanopyLeafSheathC_brch        => plt_biom%CanopyLeafSheathC_brch          ,& !output :plant branch leaf + sheath C, [g d-2]
    SeedPlantedElm_pft            => plt_biom%SeedPlantedElm_pft              ,& !output :plant stored nonstructural C at planting, [gC d-2]
    WatHeldOnCanopy_pft           => plt_ew%WatHeldOnCanopy_pft                & !output :canopy surface water content, [m3 d-2]
  )
  call PrintInfo('beg '//subname)
  !
  !     INITIALIZE SEED MORPHOLOGY AND BIOMASS
  !
  !     WTRVC,WTRVN,WTRVP=C,N,P in storage reserves (g)
  !     WTLFB,WTLFBN,WTLFBP=C,N,P in leaves (g)
  !     CanopyLeafSheathC_brch=C in leaves+PetolSheths (g)
  !     FDM-dry matter fraction (g DM C g FM C-1)
  !     CanopyBiomWater_pft,WatHeldOnCanopy_pft=water volume in,on canopy (m3)
  !     CPOOL,ZPOOL,PPOOL=C,N,P in canopy nonstructural pools (g)
  !     WTRT1,WTRT1N,WTRT1P=C,N,P in primary root layer (g)
  !     RTWT1,RTWT1N,RTWT1P=total C,N,P in primary root (g)
  !     RootMycoActiveBiomC_pvr, PopuRootMycoC_pvr=total root C mass (g)
  !     RootProteinC_pvr=total root protein C mass (g)
  !     CPOOLR,ZPOOLR,PPOOLR=C,N,P in root,myco nonstructural pools (g)
  !
  IF(jHarvstType_pft(NZ).EQ.jharvtyp_tmareseed .and. SeasonalNonstElms_pft(ielmc,NZ)>1.e-10_r8)then
    SeedPlantedElm_pft(:,NZ)=0._r8
  else
    SeedPlantedElm_pft(ielmc,NZ)    = SeedCMass_pft(NZ)*PlantPopuLive_pft(NZ)
    SeasonalNonstElms_pft(ielmc,NZ) = SeedPlantedElm_pft(ielmc,NZ)
    SeasonalNonstElms_pft(ielmn,NZ) = rNCGrain_pft(NZ)*SeasonalNonstElms_pft(ielmc,NZ)
    SeasonalNonstElms_pft(ielmp,NZ) = rPCGrain_pft(NZ)*SeasonalNonstElms_pft(ielmc,NZ)
    SeedPlantedElm_pft(2:NumPlantChemElms,NZ)=SeasonalNonstElms_pft(2:NumPlantChemElms,NZ)
  endif
  
  LeafStrutElms_brch(ielmn,1,NZ)  = rNCGrain_pft(NZ)*LeafStrutElms_brch(ielmc,1,NZ)
  LeafStrutElms_brch(ielmp,1,NZ)  = rPCGrain_pft(NZ)*LeafStrutElms_brch(ielmc,1,NZ)
  CanopyLeafSheathC_brch(1,NZ)    = LeafStrutElms_brch(ielmc,1,NZ)+PetolShethStrutElms_brch(ielmc,1,NZ)
  CanopyLeafSheathC_pft(NZ)        = CanopyLeafSheathC_pft(NZ)+CanopyLeafSheathC_brch(1,NZ)

  WatHeldOnCanopy_pft(NZ)                                      = 0._r8
  CanopyNonstElms_brch(ielmn,1,NZ)                             = rNCGrain_pft(NZ)*CanopyNonstElms_brch(ielmc,1,NZ)
  CanopyNonstElms_brch(ielmp,1,NZ)                             = rPCGrain_pft(NZ)*CanopyNonstElms_brch(ielmc,1,NZ)
  RootMyco1stStrutElms_rpvr(ielmn,NGTopRootLayer_pft(NZ),1,NZ) = rNCGrain_pft(NZ)*RootMyco1stStrutElms_rpvr(ielmc,NGTopRootLayer_pft(NZ),1,NZ)
  RootMyco1stStrutElms_rpvr(ielmp,NGTopRootLayer_pft(NZ),1,NZ) = rPCGrain_pft(NZ)*RootMyco1stStrutElms_rpvr(ielmc,NGTopRootLayer_pft(NZ),1,NZ)
  Root1stActStructElms_rpvr(:,NGTopRootLayer_pft(NZ),1,NZ)     = RootMyco1stStrutElms_rpvr(:,NGTopRootLayer_pft(NZ),1,NZ)
  RootMyco1stElm_raxs(ielmn,1,NZ)                              = rNCGrain_pft(NZ)*RootMyco1stElm_raxs(ielmc,1,NZ)
  RootMyco1stElm_raxs(ielmp,1,NZ)                              = rPCGrain_pft(NZ)*RootMyco1stElm_raxs(ielmc,1,NZ)
  RootMycoActiveBiomC_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)  = RootMyco1stStrutElms_rpvr(ielmc,NGTopRootLayer_pft(NZ),1,NZ)
  PopuRootMycoC_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)        = RootMyco1stStrutElms_rpvr(ielmc,NGTopRootLayer_pft(NZ),1,NZ)
  RootProteinC_pvr(1,NGTopRootLayer_pft(NZ),NZ)                = RootMycoActiveBiomC_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)*RootProteinCMax_pft(NZ)
  RootMycoNonstElms_rpvr(ielmn,1,NGTopRootLayer_pft(NZ),NZ)=rNCGrain_pft(NZ)&
    *RootMycoNonstElms_rpvr(ielmc,1,NGTopRootLayer_pft(NZ),NZ)
  RootMycoNonstElms_rpvr(ielmp,1,NGTopRootLayer_pft(NZ),NZ)=rPCGrain_pft(NZ) &
    *RootMycoNonstElms_rpvr(ielmc,1,NGTopRootLayer_pft(NZ),NZ)
  call PrintInfo('end '//subname)
  end associate
  end subroutine InitSeedMorphoBio
!----------------------------------------------------------------------------------------------------
  subroutine SeedDeposition(I,SeedCDeposition_pft_1d,PPmax_pft_1d,NP_col)
  implicit none
  INTEGER, INTENT(IN) :: I
  real(r8), intent(in) :: SeedCDeposition_pft_1d(:)  !input: seeding C array for all PFTs at current grid, [gC m-2]
  real(r8), intent(in) :: PPmax_pft_1d(:)            !input: maximum plant population, [# m-2]
  integer, intent(in) :: NP_col                      !input: number of pts
  integer :: NZ                !PFT index
  integer :: NE                                     !element index
  REAL(R8) :: dSeedPlantedElm_pft(NumPlantChemElms)
  
  associate(                                                  &
    AREA3                  => plt_site%AREA3                 ,& !soil area
    NU                     => plt_site%NU                    ,& !current soil surface layer
    SeedPlantedElm_pft     => plt_biom%SeedPlantedElm_pft    ,& !seeded plant elements
    SeasonalNonstElms_pft  => plt_biom%SeasonalNonstElms_pft ,& !seasonal non-structural elements
    PPatSeeding_pft        => plt_site%PPatSeeding_pft       ,& !plant population at seeding
    PPX_pft                => plt_site%PPX_pft               ,& !output :plant population, [plants m-2]        
    PlantPopuLive_pft      => plt_site%PlantPopuLive_pft     ,& !live plant population
    SeedCMass_pft          => plt_morph%SeedCMass_pft        ,& !seed carbon mass
    rNCGrain_pft           => plt_allom%rNCGrain_pft         ,& !grain N:C ratio
    rPCGrain_pft           => plt_allom%rPCGrain_pft          &  !grain P:C ratio
    )

    DO NZ=1,NP_col
      IF(SeedCDeposition_pft_1d(NZ) > 0._r8 .and. PPmax_pft_1d(NZ)>0._r8)then
        !update seeding population, update seeding planting mass
        
        dSeedPlantedElm_pft(ielmc) = AREA3(NU)*SeedCDeposition_pft_1d(NZ)
        dSeedPlantedElm_pft(ielmn) = rNCGrain_pft(NZ)*dSeedPlantedElm_pft(ielmc)
        dSeedPlantedElm_pft(ielmp) = rPCGrain_pft(NZ)*dSeedPlantedElm_pft(ielmc)

        DO NE=1,NumPlantChemElms
          SeedPlantedElm_pft(NE,NZ)=SeedPlantedElm_pft(NE,NZ)+dSeedPlantedElm_pft(NE)
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+dSeedPlantedElm_pft(NE)
        ENDDO
        !update population
        PPatSeeding_pft(NZ)   = AMIN1(PPatSeeding_pft(NZ)+SeedCDeposition_pft_1d(NZ)/SeedCMass_pft(NZ),PPmax_pft_1d(NZ))
        PlantPopuLive_pft(NZ) = PPatSeeding_pft(NZ)*AREA3(NU)
        PPX_pft(NZ)           = PPatSeeding_pft(NZ)
      endif
    ENDDO
  end associate
  end subroutine SeedDeposition  
  ![tail]
  end module InitPlantMod
