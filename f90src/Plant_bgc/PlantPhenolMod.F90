module PlantPhenolMod
!!
! Description:
! code to do plant phenology

  use data_kind_mod, only : r8 => DAT_KIND_R8,yearIJ_type
  use minimathmod, only : AZMAX1
  use InitPlantMod, only : StartPlants
  use EcoSIMCtrlMod, only : etimer
  use PlantDebugMod
  use DebugToolMod
  use EcosimConst
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantAllometryAPIData, only : plt_allom
  use PlantBiomassAPIData, only : plt_biom
  use PlantEnergyWaterAPIData, only : plt_ew
  use PlantDisturbanceAPIData, only : plt_distb
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantRootBGCAPIData, only : plt_rbgc
  use PlantMathFuncMod

  implicit none

  private


  character(len=*), parameter :: mod_filename = &
  __FILE__
!
  real(r8), PARAMETER :: PSIMin4LeafExpansion=0.1_r8                                   !minimum canopy turgor potential for leaf expansion, [MPa]
  real(r8), parameter :: PSIMin4LeafOut(2:4)=(/-0.5_r8,-1.5_r8,-0.5_r8/)               !minimum canopy water potential for leafout of drought-deciduous PFT, [MPa]
  real(r8) ,PARAMETER :: GrothStageNorm4VegetaPheno= 2.0_r8                            !normalized growth stage durations for vegetative phenology,[-]
  real(r8), PARAMETER :: GrothStageNorm4ReprodPheno= 0.667_r8                          !normalized growth stage durations for reproductive phenology,[-]
  real(r8), PARAMETER :: MaxPhotoHour4LeafOutOff   = 3600.0_r8                         !maximum light exposure hours for leaf off [h]
  !PSIMin4LeafOff should be revisited for ferns
  real(r8), parameter :: PSIMin4LeafOff(0:4)=real((/-200.0,-2.0,-2.0,-2.0,-2.0/),r8)   !minimum leaf water potential of leave off, [h]
  integer , parameter :: BranchNumMax(0:3)=(/5,1,1,1/)                            !maximum branch number

  public :: PhenologyUpdate
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine PhenologyUpdate(yearIJ)
!
!     THIS subroutine CALCULATES PLANT PHENOLOGY
!
  use PlantBalMod, only : SumPlantBiomStates
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer :: I,J
  INTEGER :: NB, NZ
  integer :: NE
  character(len=*), parameter :: subname='PhenologyUpdate'
! begin_execution
  associate(                                                                   &
    dReproNodeNumNormByMatG_brch  => plt_pheno%dReproNodeNumNormByMatG_brch ,& !output :current hourly reproductive development increment, [h-1]
    doInitPlant_pft               => plt_pheno%doInitPlant_pft                ,& !input  :PFT initialization flag:0=no,1=yes,[-]
    IsPlantActive_pft             => plt_pheno%IsPlantActive_pft              ,& !input  :flag for living pft, [-]
    iPlantCalendar_brch           => plt_pheno%iPlantCalendar_brch            ,& !input  :plant growth stage, [-]
    DATAP                         => plt_site%DATAP                           ,& !input  :parameter file name,[-]
    NP                            => plt_site%NP                              ,& !input  :current number of plant species,[-]
    PlantPopuLive_pft             => plt_site%PlantPopuLive_pft               ,& !input  :plant population, [d-2]
    MainBranchNum_pft             => plt_morph%MainBranchNum_pft              ,& !input  :id number of main branch,[-]
    PlantPopu_col                 => plt_site%PlantPopu_col                    & !inoput :total plant population, [plants d-2]
  )
  call PrintInfo('beg '//subname)
  I=yearIJ%I; J=yearIJ%J
  ! Clear the hourly increment even when a PFT or branch skips development.
  dReproNodeNumNormByMatG_brch(:,:) = 0._r8
  D9985: DO NZ=1,NP

    IF(DATAP(NZ).NE.'NO')THEN
      !
      !     PlantPopu_col=total biome population
      !
      PlantPopu_col=PlantPopu_col+PlantPopuLive_pft(NZ)
      !
      !         SET CROP FLAG ACCORDINGTopRootLayer_pftTO PLANTING, HARVEST DATES, DEATH,
      !
      call set_plant_flags(I,J,NZ)
      !
      !         INITIALIZE VARIABLES IN ACTIVE PFT
      !      
      IF(IsPlantActive_pft(NZ).EQ.iTrue)THEN

        call FindMainBranchNumber(NZ)

        call StagePlantPhenology(I,J,NZ)

        call TestPlantEmergence(I,J,NZ)

        call root_shoot_branching(I,J,NZ)
        
        IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0 .OR. doInitPlant_pft(NZ).EQ.itrue)THEN          
          call Emerged_plant_Phenology(I,J,NZ)
        ENDIF        
      ENDIF
      
    ENDIF
    
  ENDDO D9985

  call PrintInfo('end '//subname)
  RETURN
  end associate
  END subroutine PhenologyUpdate

!----------------------------------------------------------------------------------------------------
  subroutine Emerged_plant_Phenology(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), parameter :: subname='Emerged_plant_Phenology'
  integer :: NB
  REAL(R8) :: TFNP,WFNG,OFNG
  integer :: LeafNumberGrowing
  associate(                                                               &
    isPlantBranchAlive_brch     => plt_pheno%isPlantBranchAlive_brch      ,& !input  :flag to detect branch death, [-]
    LeafNumberAtFloralInit_brch => plt_pheno%LeafNumberAtFloralInit_brch  ,& !input  :leaf number at floral initiation, [-]
    PSICanopy_pft               => plt_ew%PSICanopy_pft                   ,& !input  :canopy total water potential, [Mpa]
    NumOfLeaves_brch            => plt_morph%NumOfLeaves_brch             ,& !input  :leaf number, [-]
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft            ,& !input  :number of branches,[-]
    iEmbryophyteType_pft        => plt_pheno%iEmbryophyteType_pft         ,& !input  :plant embrophyte
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft        ,& !input  :plant growth type (vascular, non-vascular),[-]
    HoursTooLowPsiCan_pft       => plt_pheno%HoursTooLowPsiCan_pft        ,& !inoput :canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY), [h]
    KHiestGroLeafNode_brch      => plt_pheno%KHiestGroLeafNode_brch       ,& !inoput :leaf growth stage counter, [-]
    doRemobilization_brch       => plt_pheno%doRemobilization_brch        ,& !output :branch phenology flag, [-]
    KLeafNumber_brch            => plt_morph%KLeafNumber_brch              & !output :leaf number, [-]
  )
  call PrintInfo('beg '//subname)
    
  call CalcPhenolEnvfactor(I,J,NZ,TFNP,WFNG,OFNG)

  D2010: DO NB=1,NumOfBranches_pft(NZ)

    IF(isPlantBranchAlive_brch(NB,NZ).EQ.iTrue)THEN
      call live_branch_phenology(I,J,NB,NZ,TFNP,WFNG,OFNG)
    ENDIF
    !
    !  KHiestGroLeafNode_brch=integer of most recent leaf number currently growing
    !
    LeafNumberGrowing=KHiestGroLeafNode_brch(NB,NZ)
    IF(LeafNumberAtFloralInit_brch(NB,NZ).LE.ppmc)THEN
      KHiestGroLeafNode_brch(NB,NZ)=INT(NumOfLeaves_brch(NB,NZ))+1
    ELSE
      KHiestGroLeafNode_brch(NB,NZ)=INT(AMIN1(NumOfLeaves_brch(NB,NZ),LeafNumberAtFloralInit_brch(NB,NZ)))+1
    ENDIF

    KLeafNumber_brch(NB,NZ)=MIN(MaxNodesPerBranch1-1,KHiestGroLeafNode_brch(NB,NZ))
    IF(KHiestGroLeafNode_brch(NB,NZ).GT.LeafNumberGrowing)THEN
      doRemobilization_brch(NB,NZ)=itrue
    ELSE
      doRemobilization_brch(NB,NZ)=ifalse
    ENDIF
    !
    call branch_specific_phenology(I,J,NB,NZ)

  ENDDO D2010
  !
  !             WATER STRESS INDICATOR
  !
  !             PSICanopy_pft=canopy total water potential
  !             PSIMin4LeafOff=minimum canopy water potential for leafoff
  !
  IF(PSICanopy_pft(NZ).LT.PSIMin4LeafOff(iEmbryophyteType_pft(NZ)))THEN
    HoursTooLowPsiCan_pft(NZ)=HoursTooLowPsiCan_pft(NZ)+1.0_r8
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine Emerged_plant_Phenology          

!----------------------------------------------------------------------------------------------------
  subroutine set_plant_flags(I,J,NZ)
!
! Description
!
  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), parameter :: subname='set_plant_flags'
  INTEGER :: L
  logical :: HarvestChk,PlantingChk,postharvst_check

! begin_execution
  associate(                                                   &
    iYearPlanting_pft     => plt_distb%iYearPlanting_pft      ,& !input  :year of planting,[-]
    iYearPlantHarvest_pft => plt_distb%iYearPlantHarvest_pft  ,& !input  :year of harvest,[-]
    iDayPlanting_pft      => plt_distb%iDayPlanting_pft       ,& !input  :day of planting,[-]
    iDayPlantHarvest_pft  => plt_distb%iDayPlantHarvest_pft   ,& !input  :day of harvest,[-]
    DATAP                 => plt_site%DATAP                   ,& !input  :parameter file name,[-]
    FireReSet_pft         => plt_pheno%FireReSet_pft          ,& !input  :flag to skill startq for fire rejuvenation
    iYearCurrent          => plt_site%iYearCurrent            ,& !input  :current year,[-]
    SeedPlantedElm_pft    => plt_biom%SeedPlantedElm_pft      ,& !input  :plant stored nonstructural C at planting, [gC d-2]
    NumActivePlants       => plt_site%NumActivePlants         ,& !inoput :number of active PFT in the grid, [-]
    Eco_NBP_CumYr_col     => plt_bgcr%Eco_NBP_CumYr_col       ,& !inoput :total NBP, [g d-2]
    iPlantStateLive_pft   => plt_pheno%iPlantStateLive_pft    ,& !inoput :flag for species death, [-]
    iMaintPlantTrait_pft  => plt_pheno%iMaintPlantTrait_pft   ,& !inoput :flag for maintain or reset plant traits, [-]
    IsPlantActive_pft     => plt_pheno%IsPlantActive_pft       & !output :flag for living pft, [-]
  )
  call PrintInfo('beg '//subname)
  !first hour of day
  IF(J.EQ.1)THEN
    HarvestChk=iDayPlanting_pft(NZ).LE.iDayPlantHarvest_pft(NZ) .OR. iYearPlanting_pft(NZ).LT.iYearPlantHarvest_pft(NZ)
    
    !Before harvest  
    IF(HarvestChk)THEN
      !planting is feasible
      IF(I.GE.iDayPlanting_pft(NZ) .OR. iYearCurrent.GT.iYearPlanting_pft(NZ))THEN
        !planted 
        postharvst_check=I.GT.iDayPlantHarvest_pft(NZ) .AND. iYearCurrent.GE.iYearPlantHarvest_pft(NZ) .AND. iPlantStateLive_pft(NZ).EQ.iFalse

        IF(postharvst_check)THEN
          !post harvest
          IsPlantActive_pft(NZ)=iFalse
        ELSE
          
          IF(I.EQ.iDayPlanting_pft(NZ) .AND. iYearCurrent.EQ.iYearPlanting_pft(NZ))THEN            
            !planting day of year
            IsPlantActive_pft(NZ)   = iFalse
            iPlantStateLive_pft(NZ) = iTrue
            iMaintPlantTrait_pft(NZ)= iTrue
            !
            if(FireReSet_pft(NZ).EQ.iFalse)THEN
              CALL StartPlants(NZ,NZ)
              Eco_NBP_CumYr_col=Eco_NBP_CumYr_col+SeedPlantedElm_pft(ielmc,NZ)
            endif  
            FireReSet_pft(NZ)=ifalse
          ENDIF

          !the living plant has actual properties set
          IF(DATAP(NZ).NE.'NO' .AND. iPlantStateLive_pft(NZ).EQ.iTrue)then
            IsPlantActive_pft(NZ)=iTrue
          endif
        ENDIF
      ELSE
        IsPlantActive_pft(NZ)=iFalse
      ENDIF
      !After harvest
    ELSE
      HarvestChk  = I.LT.iDayPlanting_pft(NZ) .AND. I.GT.iDayPlantHarvest_pft(NZ) .AND. iYearCurrent.GE.iYearPlantHarvest_pft(NZ)
      PlantingChk = I.LT.iDayPlanting_pft(NZ) .AND. iYearPlanting_pft(NZ).GT.iYearPlantHarvest_pft(NZ)

      !not planted

      IF((HarvestChk .AND. iPlantStateLive_pft(NZ).EQ.iFalse) .OR. PlantingChk)THEN
        IsPlantActive_pft(NZ)=iFalse
      ELSE
        !planting
        IF(I.EQ.iDayPlanting_pft(NZ) .AND. iYearCurrent.EQ.iYearPlanting_pft(NZ))THEN
          IsPlantActive_pft(NZ)    = iFalse
          iPlantStateLive_pft(NZ)  = iTrue
          iMaintPlantTrait_pft(NZ) = iTrue

          if(FireReSet_pft(NZ).eq.iFalse)THEN
            CALL StartPlants(NZ,NZ)
            Eco_NBP_CumYr_col = Eco_NBP_CumYr_col+SeedPlantedElm_pft(ielmc,NZ)
          ENDIF
          FireReSet_pft(NZ)=ifalse
        ENDIF
        
        IF(DATAP(NZ).NE.'NO' .AND. iPlantStateLive_pft(NZ).EQ.iTrue)then
          IsPlantActive_pft(NZ)=iTrue
        endif
      ENDIF
    ENDIF
    NumActivePlants=NumActivePlants+IsPlantActive_pft(NZ)
  ENDIF
  
  call PrintInfo('end '//subname)
  end associate
  end subroutine set_plant_flags

!----------------------------------------------------------------------------------------------------
  subroutine root_shoot_branching(I,J,NZ)

  implicit none
  integer, intent(in) :: I,J,NZ
  integer :: NB
  character(len=*), parameter :: subname='root_shoot_branching'
  integer :: BranchNumber_new
  logical :: checkRootInitializer,shootCheck
! begin_execution
  associate(                                                                 &
    CanopyNonstElmConc_pft       => plt_biom%CanopyNonstElmConc_pft         ,& !input  :canopy nonstructural element concentration, [g d-2]
    SeasonalNonstElms_pft        => plt_biom%SeasonalNonstElms_pft          ,& !input  :plant stored nonstructural element at current step, [g d-2]
    iPlantCalendar_brch          => plt_pheno%iPlantCalendar_brch           ,& !input  :plant growth stage, [-]
    doInitPlant_pft              => plt_pheno%doInitPlant_pft               ,& !input  :PFT initialization flag:0=no,1=yes,[-]
    iPlantPhenolPattern_pft      => plt_pheno%iPlantPhenolPattern_pft       ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantTurnoverPattern_pft    => plt_pheno%iPlantTurnoverPattern_pft     ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    NonstCMinCon2InitRoot_pft    => plt_pheno%NonstCMinCon2InitRoot_pft     ,& !input  :threshold root nonstructural C content for initiating new root axis, [gC gC-1]
    MatureGroup_pft              => plt_pheno%MatureGroup_pft               ,& !input  :acclimated plant maturity group, [-]
    NonstCMinConc2InitBranch_pft => plt_pheno%NonstCMinConc2InitBranch_pft  ,& !input  :branch nonstructural C content required for new branch, [gC gC-1]
    PlantPopuLive_pft            => plt_site%PlantPopuLive_pft              ,& !input  :plant population, [d-2]
    PSIRootTurg_vr               => plt_ew%PSIRootTurg_vr                   ,& !input  :root turgor water potential, [Mpa]
    FracGroth2Node_pft           => plt_allom%FracGroth2Node_pft            ,& !input  :parameter for allocation of growth to nodes, [-]
    MainBranchNum_pft            => plt_morph%MainBranchNum_pft             ,& !input  :id number of main branch,[-]
    NumCogrowthNode_pft          => plt_morph%NumCogrowthNode_pft           ,& !input  :number of concurrently growing nodes,[-]
    NGTopRootLayer_pft           => plt_morph%NGTopRootLayer_pft            ,& !input  :soil layer at planting depth, [-]
    ShootNodeNumAtPlanting_pft   => plt_morph%ShootNodeNumAtPlanting_pft    ,& !input  :number of nodes in seed, [-]
    ShootNodeNum_brch            => plt_morph%ShootNodeNum_brch             ,& !input  :shoot node number, [-]
    isPlantBranchAlive_brch      => plt_pheno%isPlantBranchAlive_brch       ,& !inoput :flag to detect branch death, [-]
    NumStructuralRootAxes_pft    => plt_morph%NumStructuralRootAxes_pft     ,& !inoput :number of structural root axes per plant,[-]
    NumOfBranches_pft            => plt_morph%NumOfBranches_pft             ,& !inoput :number of branches,[-]
    BranchNumber_pft             => plt_morph%BranchNumber_pft              ,& !inoput :main branch numeric id,[-]
    MatureGroup_brch             => plt_pheno%MatureGroup_brch              ,& !output :branch level plant maturity group, [-]
    isPlantRootAlive_pft         => plt_pheno%isPlantRootAlive_pft          ,& !output :flag to detect root system death,[-]
    isPlantShootAlive_pft        => plt_pheno%isPlantShootAlive_pft         ,& !output :flag to detect canopy death,[-]
    Hours4Leafout_brch           => plt_pheno%Hours4Leafout_brch            ,& !output :heat requirement for spring leafout/dehardening, [h]
    BranchNumerID_brch           => plt_morph%BranchNumerID_brch             & !output :branch meric id, [-]
  )
  call PrintInfo('beg '//subname)
  !
  ! ADD BRANCH TO SHOOT IF PLANT GROWTH STAGE, SHOOT NON-STRUCTURAL
  ! CONCENTRATION PERMIT
  !  
  IF(doInitPlant_pft(NZ).EQ.ifalse)THEN
    !plant initialized

    IF(J.EQ.1 .AND. PlantPopuLive_pft(NZ).GT.0.0_r8)THEN
      !first hour of the day, population > 0

      IF(PSIRootTurg_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ).GT.PSIMin4LeafExpansion)THEN
        IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .OR. iPlantCalendar_brch(ipltcal_InitFloral,MainBranchNum_pft(NZ),NZ).EQ.0)THEN
          !perennial plant or flower not initiated for annual plant 

          IF((NumOfBranches_pft(NZ).EQ.0 .AND. SeasonalNonstElms_pft(ielmc,NZ).GT.0.0_r8) &
            .OR. (CanopyNonstElmConc_pft(ielmc,NZ).GT.NonstCMinConc2InitBranch_pft(NZ) &
            .AND. NonstCMinConc2InitBranch_pft(NZ).GT.0.0_r8))THEN

            D120: DO NB=1,MaxNumBranches

              IF(isPlantBranchAlive_brch(NB,NZ).EQ.iFalse)THEN
                BranchNumber_new=BranchNumber_pft(NZ)+NumCogrowthNode_pft(NZ)/FracGroth2Node_pft(NZ)+ShootNodeNumAtPlanting_pft(NZ)
                
                IF(NB.EQ.MainBranchNum_pft(NZ) .OR. ShootNodeNum_brch(MainBranchNum_pft(NZ),NZ).GT.BranchNumber_new)THEN
                  !initiate a new branch
                  BranchNumber_pft(NZ)          = BranchNumber_pft(NZ)+1
                  NumOfBranches_pft(NZ)         = MIN(BranchNumMax(iPlantTurnoverPattern_pft(NZ)),MAX(NB,NumOfBranches_pft(NZ)))
                  BranchNumerID_brch(NB,NZ)     = BranchNumber_pft(NZ)-1
                  isPlantShootAlive_pft(NZ)      = iTrue
                  isPlantBranchAlive_brch(NB,NZ) = iTrue
                  Hours4Leafout_brch(NB,NZ)     = 0.0_r8

                  IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
                    MatureGroup_brch(NB,NZ)=AZMAX1(MatureGroup_pft(NZ)-BranchNumerID_brch(NB,NZ))
                  ELSE
                    MatureGroup_brch(NB,NZ)=MatureGroup_pft(NZ)
                  ENDIF
                  exit
                ENDIF
              ENDIF
            ENDDO D120            
          ENDIF
        ENDIF
      ENDIF
      !
      !     ADD STRUCTURAL AXIS TO ROOT IF PLANT GROWTH STAGE, ROOT NON-STRUCTURAL C
      !     CONCENTRATION PERMIT
      !
      !     PR=nonstructural C concentration needed for root branching
      !     FracGroth2Node_pft: parameter for allocation of growth to nodes
      !     XLI: number of nodes in seed
      !     PSTG: node number
      !     MainBranchNum_pft: number of main branch
      !     CanopyNonstElmConc_pft: canopy nonstructural element concentration
      !     PSIRootTurg_vr: root turgor pressure
      !     SeasonalNonstElms_pft: non-structural carbon
      !     root axis initialization      
      IF(PSIRootTurg_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ).GT.PSIMin4LeafExpansion)THEN !water condition met
        IF(NumStructuralRootAxes_pft(NZ).EQ.0 .OR. ShootNodeNum_brch(MainBranchNum_pft(NZ),NZ) &
          .GT.NumStructuralRootAxes_pft(NZ)/FracGroth2Node_pft(NZ)+ShootNodeNumAtPlanting_pft(NZ))THEN

          checkRootInitializer= (NumStructuralRootAxes_pft(NZ).EQ.0 .AND. SeasonalNonstElms_pft(ielmc,NZ).GT.0.0_r8) &              !storage/seed 
            .OR. (CanopyNonstElmConc_pft(ielmc,NZ).GT.NonstCMinCon2InitRoot_pft(NZ) .AND. NonstCMinCon2InitRoot_pft(NZ).GT.0.0_r8)  !plant status

          IF(checkRootInitializer)THEN
            NumStructuralRootAxes_pft(NZ) = MIN(MaxNumRootAxes,NumStructuralRootAxes_pft(NZ)+1)
            isPlantRootAlive_pft(NZ)  = iTrue
          ENDIF
        ENDIF
      ENDIF
    ENDIF

  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine root_shoot_branching

!----------------------------------------------------------------------------------------------------
  subroutine FindMainBranchNumber(NZ)
  implicit none
  integer, intent(in) :: NZ
  integer :: NB,BranchNumberX_pft
  character(len=*), parameter :: subname='FindMainBranchNumber'
  associate(                                                      &
    isPlantBranchAlive_brch => plt_pheno%isPlantBranchAlive_brch  ,& !input  :flag to detect branch death, [-]
    NumOfBranches_pft       => plt_morph%NumOfBranches_pft        ,& !input  :number of branches,[-]
    BranchNumerID_brch      => plt_morph%BranchNumerID_brch       ,& !input  :branch meric id, [-]
    MainBranchNum_pft       => plt_morph%MainBranchNum_pft         & !output :id number of main branch,[-]
  )
  call PrintInfo('beg '//subname)
  MainBranchNum_pft(NZ) = 1
  BranchNumberX_pft     = 1.0E+06_r8

  DD140: DO NB=1,NumOfBranches_pft(NZ)
    IF(isPlantBranchAlive_brch(NB,NZ).EQ.iTrue)THEN
      !find main branch number, which is the most recent live branch      
      IF(BranchNumerID_brch(NB,NZ).LT.BranchNumberX_pft)THEN
        MainBranchNum_pft(NZ) = NB
        BranchNumberX_pft     = BranchNumerID_brch(NB,NZ)
      ENDIF      
    ENDIF
  ENDDO DD140
  call PrintInfo('end '//subname)
  end associate
  end subroutine FindMainBranchNumber

!----------------------------------------------------------------------------------------------------
  subroutine StagePlantPhenology(I,J,NZ)
  ! Refresh nonstructural element concentrations and their canopy C+N+P sum.
  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), parameter :: subname='StagePlantPhenology'
  integer :: NB,N,L,NE

  associate(                                                            &
    CanopyLeafSheathC_pft       => plt_biom%CanopyLeafSheathC_pft      ,& !input  :canopy leaf + sheath C, [g d-2]
    CanopyLeafSheathC_brch      => plt_biom%CanopyLeafSheathC_brch     ,& !input  :plant branch leaf + sheath C, [g d-2]
    CanopyNodulNonstElms_brch   => plt_biom%CanopyNodulNonstElms_brch  ,& !input  :branch nodule nonstructural element, [g d-2]
    RootMycoNonstElms_rpvr      => plt_biom%RootMycoNonstElms_rpvr     ,& !input  :root layer nonstructural element, [g d-2]
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch       ,& !input  :branch nonstructural element, [g d-2]
    ZERO4LeafVar_pft            => plt_biom%ZERO4LeafVar_pft           ,& !input  :threshold zero for leaf calculation, [-]
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft             ,& !input  :threshold zero for plang growth calculation, [-]
    RootMycoActiveBiomC_pvr     => plt_biom%RootMycoActiveBiomC_pvr    ,& !input  :root layer structural C, [gC d-2]
    NU                          => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    isPlantBranchAlive_brch     => plt_pheno%isPlantBranchAlive_brch   ,& !input  :flag to detect branch death, [-]
    Myco_pft                    => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    NMaxRootBotLayer_pft        => plt_morph%NMaxRootBotLayer_pft      ,& !input  :maximum soil layer number for all root axes, [-]
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft         ,& !input  :number of branches,[-]
    CanopyNonstElms_pft         => plt_biom%CanopyNonstElms_pft        ,& !inoput :canopy nonstructural element concentration, [g d-2]
    CanopyNodulNonstElms_pft    => plt_biom%CanopyNodulNonstElms_pft   ,& !inoput :canopy nodule nonstructural element, [g d-2]
    NGTopRootLayer_pft          => plt_morph%NGTopRootLayer_pft        ,& !inoput :soil layer at planting depth, [-]
    CanopyNoduleNonstCConc_pft  => plt_biom%CanopyNoduleNonstCConc_pft ,& !output :canopy nodule nonstructural C concentration, [gC (gC)-1]
    LeafPetoNonstElmConc_brch   => plt_biom%LeafPetoNonstElmConc_brch  ,& !output :branch nonstructural C concentration, [g d-2]
    CanopyNonstElmConc_pft      => plt_biom%CanopyNonstElmConc_pft     ,& !output :canopy nonstructural element concentration, [g d-2]
    CanopyNonstOMConc_pft       => plt_biom%CanopyNonstOMConc_pft      ,& !output :canopy nonstructural C+N+P concentration, [g gC-1]
    RootNonstructElmConc_pvr   => plt_biom%RootNonstructElmConc_pvr  ,& !output :root layer nonstructural C concentration, [g g-1]
    RootNonstOMConc_pvr        => plt_biom%RootNonstOMConc_pvr       ,& !output :root/mycorrhizal nonstructural C+N+P concentration, [g gC-1]
    MaxSoilLays4Root_pft        => plt_morph%MaxSoilLays4Root_pft       & !output :maximum soil layer number for all root axes,[-]
  )
  call PrintInfo('beg '//subname)
  plt_bgcr%RootGasLossDisturb_pft(idg_beg:idg_NH3,NZ)=0.0_r8
  CanopyNodulNonstElms_pft(1:NumPlantChemElms,NZ)=0._r8
  CanopyNonstElms_pft(1:NumPlantChemElms,NZ)=0.0_r8
  MaxSoilLays4Root_pft(NZ)   = NMaxRootBotLayer_pft(NZ)
  NGTopRootLayer_pft(NZ) = MIN(MaxSoilLays4Root_pft(NZ),MAX(NGTopRootLayer_pft(NZ),NU))
  
  !
  ! TOTAL PLANT NON-STRUCTURAL C, N, P
  !
  ! CPOOL*,ZPOOL*,PPOOL*=non-structl C,N,P in branch(NB),canopy(g)
  ! CPOLN*,ZPOLN*,PPOLN*=non-structl C,N,P in branch,canopy nodules (g)
  ! MainBranchNum_pft=main branch number
  !
  D140: DO NB=1,NumOfBranches_pft(NZ)
    IF(isPlantBranchAlive_brch(NB,NZ).EQ.iTrue)THEN
      DO NE=1,NumPlantChemElms      
        CanopyNonstElms_pft(NE,NZ)=CanopyNonstElms_pft(NE,NZ)+CanopyNonstElms_brch(NE,NB,NZ)
        CanopyNodulNonstElms_pft(NE,NZ)=CanopyNodulNonstElms_pft(NE,NZ)+CanopyNodulNonstElms_brch(NE,NB,NZ)    
      ENDDO     
    ENDIF
  ENDDO D140
  !
  ! NON-STRUCTURAL C, N, P CONCENTRATIONS IN ROOT
  !
  ! WTRTL=root mas(g)
  ! CPOOLR,ZPOOLR,PPOOLR=non-structl C,N,P in root(1),myco(2)(g)
  ! CCPOLR,CZPOLR,CPPOLR=non-structl C,N,P concn in root(1),myco(2)(g g-1)
  !
  D180: DO N=1,Myco_pft(NZ)
    D160: DO L=NU,MaxSoilLays4Root_pft(NZ)
      IF(RootMycoActiveBiomC_pvr(N,L,NZ).GT.ZERO4LeafVar_pft(NZ))THEN
        DO NE=1,NumPlantChemElms
          RootNonstructElmConc_pvr(NE,N,L,NZ)=AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)/RootMycoActiveBiomC_pvr(N,L,NZ))
        ENDDO
      ELSE
        DO NE=1,NumPlantChemElms
          RootNonstructElmConc_pvr(NE,N,L,NZ)=1.0_r8
        ENDDO
      ENDIF
      RootNonstOMConc_pvr(N,L,NZ)=sum(RootNonstructElmConc_pvr(1:NumPlantChemElms,N,L,NZ))
    ENDDO D160
  ENDDO D180
  !
  ! NON-STRUCTURAL C, N, P CONCENTRATIONS IN SHOOT
  !
  ! CCPOLP,CZPOLP,CPPOLP=nonstructural C,N,P concn in canopy(g g-1)
  ! CanopyNoduleNonstCConc_pft=nonstructural C concentration in canopy nodules
  ! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
  !
  IF(CanopyLeafSheathC_pft(NZ).GT.ZERO4LeafVar_pft(NZ))THEN
    DO NE=1,NumPlantChemElms
      CanopyNonstElmConc_pft(NE,NZ)=AZMAX1(AMIN1(1.0_r8,CanopyNonstElms_pft(NE,NZ)/CanopyLeafSheathC_pft(NZ)))
    ENDDO
    CanopyNoduleNonstCConc_pft(NZ)=AZMAX1(AMIN1(1.0_r8,CanopyNodulNonstElms_pft(ielmc,NZ)/CanopyLeafSheathC_pft(NZ)))
  ELSE
    CanopyNonstElmConc_pft(1:NumPlantChemElms,NZ) = 1.0_r8
    CanopyNoduleNonstCConc_pft(NZ)                = 1.0_r8
  ENDIF
  CanopyNonstOMConc_pft(NZ)=CanopyNonstElmConc_pft(ielmc,NZ) &
    +CanopyNonstElmConc_pft(ielmn,NZ)+CanopyNonstElmConc_pft(ielmp,NZ)
  
  D190: DO NB=1,NumOfBranches_pft(NZ)
    IF(CanopyLeafSheathC_brch(NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      DO NE=1,NumPlantChemElms      
        LeafPetoNonstElmConc_brch(NE,NB,NZ)=AZMAX1(CanopyNonstElms_brch(NE,NB,NZ)/CanopyLeafSheathC_brch(NB,NZ))
      ENDDO  
    ELSE
      LeafPetoNonstElmConc_brch(1:NumPlantChemElms,NB,NZ)=1.0_r8
    ENDIF    
  ENDDO D190
  call PrintInfo('end '//subname)
  end associate
  end subroutine StagePlantPhenology

!----------------------------------------------------------------------------------------------------
  subroutine TestPlantEmergence(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), parameter :: subname='TestPlantEmergence'
  real(r8):: ShootArea
  logical :: CanopyChk,RootChk

  associate(                                                       &
    MainBranchNum_pft       => plt_morph%MainBranchNum_pft        ,& !input  :number of main branch,[-]
    CanopyLeafArea_pft      => plt_morph%CanopyLeafArea_pft       ,& !input  :plant canopy leaf area, [m2 d-2]
    ShootElms_pft           => plt_biom%ShootElms_pft             ,& !input  :canopy shoot structural chemical element mass, [g d-2]
    HypocotHeight_pft       => plt_morph%HypocotHeight_pft        ,& !input  :cotyledon height, [m]
    Root1stDepz_raxes       => plt_morph%Root1stDepz_raxes        ,& !input  :root layer depth, [m]
    ZERO4LeafVar_pft        => plt_biom%ZERO4LeafVar_pft          ,& !input  :threshold zero for leaf calculation, [-]
    WatHeldOnCanopy_pft     => plt_ew%WatHeldOnCanopy_pft         ,& !input  :canopy surface water content, [m3 d-2]
    CanopyStemSurfArea_pft  => plt_morph%CanopyStemSurfArea_pft   ,& !input  :plant stem area, [m2 d-2]
    SeedDepth_pft           => plt_morph%SeedDepth_pft            ,& !input  :seeding depth, [m]
    iPlantCalendar_brch     => plt_pheno%iPlantCalendar_brch      ,& !input  :plant growth stage, [-]
    VHeatCapCanopy_pft      => plt_ew%VHeatCapCanopy_pft           & !output :canopy heat capacity, [MJ d-2 K-1]
  )
  call PrintInfo('beg '//subname)
  !
  ! EMERGENCE DATE FROM COTYLEDON HEIGHT, LEAF AREA, ROOT DEPTH
  !
  ! iPlantCalendar_brch(ipltcal_Emerge,=emergence date
  ! CanopyLeafArea_pft,CanopyStemSurfArea_pft=leaf,stalk areas
  ! HypocotHeight_pft=hypocotyledon height
  ! SeedDepth_pft=seeding depth
  ! Root1stDepz_raxes=primary root depth,avoid plant the seed at the grid interface
  ! VHeatCapCanopy_pft,WTSHT,WatHeldOnCanopy_pft=canopy heat capacity,mass,water content
  !
  
  IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).EQ.0)THEN
    ShootArea = CanopyLeafArea_pft(NZ)+CanopyStemSurfArea_pft(NZ)
    CanopyChk = (HypocotHeight_pft(NZ)+1.e-3_r8.GT.SeedDepth_pft(NZ)).AND.(ShootArea.GT.ZERO4LeafVar_pft(NZ))
    RootChk   = Root1stDepz_raxes(1,NZ).GT.(SeedDepth_pft(NZ)+1.e-8_r8)
    
    IF(CanopyChk .AND. RootChk)THEN
      iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ)=I
      VHeatCapCanopy_pft(NZ)=cpw*(ShootElms_pft(ielmc,NZ)*10.0E-06_r8+WatHeldOnCanopy_pft(NZ))
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine TestPlantEmergence

!----------------------------------------------------------------------------------------------------
  subroutine branch_specific_phenology(I,J,NB,NZ)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  character(len=*), parameter :: subname='branch_specific_phenology'

  associate(                                                                   &
    PSICanopy_pft                 => plt_ew%PSICanopy_pft                     ,& !input  :canopy total water potential, [Mpa]
    DayLenthCurrent               => plt_site%DayLenthCurrent                 ,& !input  :current daylength of the grid, [h]
    DayLenthPrev                  => plt_site%DayLenthPrev                    ,& !input  :daylength of previous day, [h]
    iPlantCalendar_brch           => plt_pheno%iPlantCalendar_brch            ,& !input  :plant growth stage, [-]
    iPlantRootProfile_pft         => plt_pheno%iPlantRootProfile_pft          ,& !input  :plant growth type (vascular, non-vascular),[-]
    HourReq4LeafOut_brch          => plt_pheno%HourReq4LeafOut_brch           ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    iEmbryophyteType_pft          => plt_pheno%iEmbryophyteType_pft           ,& !input  :plant embrophyte
    iPlantPhenolType_pft          => plt_pheno%iPlantPhenolType_pft           ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    doInitPlant_pft               => plt_pheno%doInitPlant_pft                ,& !input  :PFT initialization flag:0=no,1=yes,[-]
    isPlantBranchAlive_brch       => plt_pheno%isPlantBranchAlive_brch        ,& !input  :flag to detect branch death, [-]
    HourReq4LeafOff_brch          => plt_pheno%HourReq4LeafOff_brch           ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    doPlantLeaveOff_brch          => plt_pheno%doPlantLeaveOff_brch           ,& !input  :branch phenology flag, [-]
    EnablePlantLeafOut_brch       => plt_pheno%EnablePlantLeafOut_brch        ,& !inoput :branch phenology flag, [-]
    Hours4ShortenPhotoPeriod_brch => plt_pheno%Hours4ShortenPhotoPeriod_brch  ,& !inoput :initial cold requirement for autumn leafoff/hardening, [h]
    Hours4Leafout_brch            => plt_pheno%Hours4Leafout_brch             ,& !inoput :heat requirement for spring leafout/dehardening, [h]
    Hours4LenthenPhotoPeriod_brch => plt_pheno%Hours4LenthenPhotoPeriod_brch  ,& !inoput :initial heat requirement for spring leafout/dehardening, [h]
    Hours4LeafOff_brch            => plt_pheno%Hours4LeafOff_brch              & !inoput :cold requirement for autumn leafoff/hardening, [h]
  )
  call PrintInfo('beg '//subname)
  !               PHENOLOGY
  !
  !           DayLenthPrev,DLYN=daylength of previous,current day
  !           Hours4LenthenPhotoPeriod_brch,Hours4ShortenPhotoPeriod_brch
  !            =hourly counter for lengthening,shortening photoperiods
  !  
  IF(isPlantBranchAlive_brch(NB,NZ).EQ.iTrue .OR. doInitPlant_pft(NZ).EQ.itrue)THEN
    IF(DayLenthCurrent.GE.DayLenthPrev)THEN
      Hours4LenthenPhotoPeriod_brch(NB,NZ)=Hours4LenthenPhotoPeriod_brch(NB,NZ)+1.0_r8
      Hours4ShortenPhotoPeriod_brch(NB,NZ)=0.0_r8
    ELSE
      Hours4LenthenPhotoPeriod_brch(NB,NZ)=0.0_r8
      Hours4ShortenPhotoPeriod_brch(NB,NZ)=Hours4ShortenPhotoPeriod_brch(NB,NZ)+1.0_r8
    ENDIF
    !
    IF(iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen)THEN

      call CropEvergreenPhenology(I,J,NB,NZ)
      !
    ELSEIF(iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid)THEN

      call ColdDeciduousBranchPhenology(I,J,NB,NZ)
      !
      !     CALCULATE DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayer_pftHOURS
      !     ABOVE SPECIFIED WATER POTENTIAL DURINGTopRootLayer_pftDORMANCY
      !
    ELSEIF(iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu &           !drought deciduos
      .OR. iPlantPhenolType_pft(NZ).EQ.4 &
      .OR. iPlantPhenolType_pft(NZ).EQ.5)THEN
      !
      !type 4 and 5 are place holders for subtropical/tropical evergreen
      !
      IF(EnablePlantLeafOut_brch(NB,NZ).EQ.iTrue)THEN
        IF(PSICanopy_pft(NZ).GE.PSIMin4LeafOut(iEmbryophyteType_pft(NZ)))THEN
          Hours4Leafout_brch(NB,NZ)=Hours4Leafout_brch(NB,NZ)+1.0_r8
        ENDIF

        IF(Hours4Leafout_brch(NB,NZ).LT.HourReq4LeafOut_brch(NB,NZ))THEN
          IF(PSICanopy_pft(NZ).LT.PSIMin4LeafOff(iEmbryophyteType_pft(NZ)))THEN
            Hours4Leafout_brch(NB,NZ)=AZMAX1(Hours4Leafout_brch(NB,NZ)-12.0_r8)
          ENDIF
        ENDIF

        IF(Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN
          Hours4LeafOff_brch(NB,NZ)=0.0_r8
          IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff_brch(NB,NZ)=iTrue
        ENDIF
      ENDIF

      IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff_brch(NB,NZ)=iTrue
      !
      !     CALCULATE DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayer_pftHOURS
      !     BELOW SPECIFIED WATER POTENTIAL DURINGTopRootLayer_pftGROWINGTopRootLayer_pftSEASON
      !
      !     Hours4Leafout_brch=leafout hours,hours required for leafout
      !     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
      !     EnablePlantLeafOut_brch,doPlantLeaveOff_brch=flag for enabling leafout,leafoff:0=enable,1=disable
      !     PSICanopy_pft=canopy total water potential
      !     ALAT=latitude
      !     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
      !     Hours4LenthenPhotoPeriod_brch,Hours4ShortenPhotoPeriod_brch=hourly counter for lengthening,shortening photoperiods
      !
      IF(EnablePlantLeafOut_brch(NB,NZ).EQ.iFalse .AND. doPlantLeaveOff_brch(NB,NZ).EQ.iTrue)THEN
        IF(PSICanopy_pft(NZ).LT.PSIMin4LeafOff(iEmbryophyteType_pft(NZ)))THEN
          Hours4LeafOff_brch(NB,NZ)=Hours4LeafOff_brch(NB,NZ)+1.0_r8
        ENDIF

        IF(iPlantPhenolType_pft(NZ).EQ.4)THEN
          IF(Hours4ShortenPhotoPeriod_brch(NB,NZ).GT.MaxPhotoHour4LeafOutOff)THEN
            Hours4LeafOff_brch(NB,NZ)=Hours4ShortenPhotoPeriod_brch(NB,NZ)
          ENDIF
        ELSEIF(iPlantPhenolType_pft(NZ).EQ.5)THEN
          IF(Hours4LenthenPhotoPeriod_brch(NB,NZ).GT.MaxPhotoHour4LeafOutOff)THEN
            Hours4LeafOff_brch(NB,NZ)=Hours4LenthenPhotoPeriod_brch(NB,NZ)
          ENDIF
        ENDIF

        IF(Hours4LeafOff_brch(NB,NZ).GE.HourReq4LeafOff_brch(NB,NZ) .AND. EnablePlantLeafOut_brch(NB,NZ).EQ.iFalse)THEN
          Hours4Leafout_brch(NB,NZ)  = 0.0_r8
          EnablePlantLeafOut_brch(NB,NZ) = iTrue
        ENDIF
      ENDIF
  
    ELSEIF(iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid)THEN
      call ColdDroughtDeciduPhenology(I,J,NB,NZ)

    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine branch_specific_phenology

!----------------------------------------------------------------------------------------------------
  subroutine ColdDroughtDeciduPhenology(I,J,NB,NZ)
  ! 
  !     CALCULATE WINTER AND DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING
  !     HOURS ABOVE SPECIFIED TEMPERATURE OR WATER POTENTIAL DURING
  !     LENGTHENINGTopRootLayer_pftPHOTOPERIODS
  !
  implicit none
  integer, intent(in) :: I   !day
  integer, intent(in) :: J   !hour
  integer, intent(in) :: NB  !branch id
  integer, intent(in) :: NZ  !pft id
  logical :: daylight_check
  character(len=*), parameter :: subname='ColdDroughtDeciduPhenology'
  associate(                                                           &
    DayLenthCurrent           => plt_site%DayLenthCurrent             ,& !input  :current daylength of the grid, [h]
    DayLenthMax_col           => plt_site%DayLenthMax_col             ,& !input  :maximum daylength of the grid, [h]
    PSICanopyTurg_pft         => plt_ew%PSICanopyTurg_pft             ,& !input  :plant canopy turgor water potential, [MPa]
    DayLenthPrev              => plt_site%DayLenthPrev                ,& !input  :daylength of previous day, [h]
    HourReq4LeafOut_brch      => plt_pheno%HourReq4LeafOut_brch       ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch        ,& !input  :plant growth stage, [-]
    doPlantLeaveOff_brch      => plt_pheno%doPlantLeaveOff_brch       ,& !input  :branch phenology flag, [-]
    iEmbryophyteType_pft      => plt_pheno%iEmbryophyteType_pft       ,& !input  :plant embrophyte    
    TCGroth_pft               => plt_pheno%TCGroth_pft                ,& !input  :canopy growth temperature, [oC]
    TC4LeafOff_pft            => plt_pheno%TC4LeafOff_pft             ,& !input  :threshold temperature for autumn leafoff/hardening, [oC]
    TCChill4Seed_pft          => plt_pheno%TCChill4Seed_pft           ,& !input  :temperature below which seed set is adversely affected, [oC]
    TC4LeafOut_pft            => plt_pheno%TC4LeafOut_pft             ,& !input  :threshold temperature for spring leafout/dehardening, [oC]
    PSICanopy_pft             => plt_ew%PSICanopy_pft                 ,& !input  :canopy total water potential, [Mpa]
    HourReq4LeafOff_brch      => plt_pheno%HourReq4LeafOff_brch       ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft      ,& !input  :plant growth type (vascular, non-vascular),[-]
    EnablePlantLeafOut_brch   => plt_pheno%EnablePlantLeafOut_brch    ,& !inoput :branch phenology flag, [-]
    Hours4LeafOff_brch        => plt_pheno%Hours4LeafOff_brch         ,& !inoput :cold requirement for autumn leafoff/hardening, [h]
    Hours4Leafout_brch        => plt_pheno%Hours4Leafout_brch          & !inoput :heat requirement for spring leafout/dehardening, [h]
  )
  call PrintInfo('beg '//subname)
  daylight_check=(DayLenthCurrent.GE.DayLenthPrev .OR. DayLenthCurrent.GE.DayLenthMax_col-2.0_r8)
  IF( daylight_check .AND. EnablePlantLeafOut_brch(NB,NZ).EQ.iTrue)THEN
    !
    IF(TCGroth_pft(NZ).GE.TC4LeafOut_pft(NZ) .AND. PSICanopyTurg_pft(NZ).GT.PSIMin4LeafExpansion)THEN !temeprature and moisture condition
      Hours4Leafout_brch(NB,NZ)=Hours4Leafout_brch(NB,NZ)+1.0_r8
    ENDIF

    !chilling effect
    IF(Hours4Leafout_brch(NB,NZ).LT.HourReq4LeafOut_brch(NB,NZ))THEN
      IF(TCGroth_pft(NZ).LT.TCChill4Seed_pft(NZ) .OR. PSICanopyTurg_pft(NZ).LT.PSIMin4LeafExpansion)THEN
        Hours4Leafout_brch(NB,NZ)=AZMAX1(Hours4Leafout_brch(NB,NZ)-1.5_r8)
      ENDIF
    ENDIF

    IF(Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN
      Hours4LeafOff_brch(NB,NZ)=0.0_r8
      IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff_brch(NB,NZ)=iTrue
    ENDIF
  ENDIF

  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff_brch(NB,NZ)=iTrue
  !
  !     CALCULATE WINTER AND DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING
  !     HOURS BELOW SPECIFIED TEMPERATURE OR WATER POTENTIAL DURING
  !     SHORTENINGTopRootLayer_pftPHOTOPERIODS
  !
  !     DayLenthPrev,DLYN=daylength of previous,current day
  !     Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
  !     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
  !     EnablePlantLeafOut_brch,doPlantLeaveOff_brch=flag for enabling leafout,leafoff:0=enable,1=disable
  !     TCGroth_pft,TCZ,TCChill4Seed_pft=canopy temp,leafout threshold temp,chilling temp
  !     PSICanopy_pft=canopy total water potential
  !     ALAT=latitude
  !     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
  !
  IF((DayLenthCurrent.LT.DayLenthPrev .OR. DayLenthCurrent.LT.24.0_r8-DayLenthMax_col+2.0_r8) &
    .AND. doPlantLeaveOff_brch(NB,NZ).EQ.iTrue)THEN
    IF(TCGroth_pft(NZ).LE.TC4LeafOff_pft(NZ) .OR. &
      PSICanopy_pft(NZ).LT.PSIMin4LeafOff(iEmbryophyteType_pft(NZ)))THEN
      Hours4LeafOff_brch(NB,NZ)=Hours4LeafOff_brch(NB,NZ)+1.0_r8
    ENDIF
    IF(Hours4LeafOff_brch(NB,NZ).GE.HourReq4LeafOff_brch(NB,NZ) &
      .AND. EnablePlantLeafOut_brch(NB,NZ).EQ.iFalse)THEN
      Hours4Leafout_brch(NB,NZ)  = 0.0_r8
      EnablePlantLeafOut_brch(NB,NZ) = iTrue
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine ColdDroughtDeciduPhenology

!----------------------------------------------------------------------------------------------------
  subroutine CropEvergreenPhenology(I,J,NB,NZ)

  ! CALCULATE EVERGREEN PHENOLOGY DURING TopRootLayer_pft LENGTHENING TopRootLayer_pft PHOTOPERIODS
  implicit none
  integer, intent(in) :: I   !day
  integer, intent(in) :: J   !hour
  integer, intent(in) :: NB  !branch id
  integer, intent(in) :: NZ  !pft id
  character(len=*), parameter :: subname='CropEvergreenPhenology'

  associate(                                                                   &
    DayLenthCurrent               => plt_site%DayLenthCurrent                 ,& !input  :current daylength of the grid, [h]
    DayLenthPrev                  => plt_site%DayLenthPrev                    ,& !input  :daylength of previous day, [h]
    HourReq4LeafOff_brch          => plt_pheno%HourReq4LeafOff_brch           ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    Hours4LenthenPhotoPeriod_brch => plt_pheno%Hours4LenthenPhotoPeriod_brch  ,& !input  :initial heat requirement for spring leafout/dehardening, [h]
    HourReq4LeafOut_brch          => plt_pheno%HourReq4LeafOut_brch           ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    Hours4ShortenPhotoPeriod_brch => plt_pheno%Hours4ShortenPhotoPeriod_brch  ,& !input  :initial cold requirement for autumn leafoff/hardening, [h]
    ALAT                          => plt_site%ALAT                            ,& !input  :latitude, [degrees north]
    Hours4Leafout_brch            => plt_pheno%Hours4Leafout_brch             ,& !output :heat requirement for spring leafout/dehardening, [h]
    Hours4LeafOff_brch            => plt_pheno%Hours4LeafOff_brch             ,& !output :cold requirement for autumn leafoff/hardening, [h]
    EnablePlantLeafOut_brch       => plt_pheno%EnablePlantLeafOut_brch        ,& !output :branch phenology flag, [-]
    doPlantLeaveOff_brch          => plt_pheno%doPlantLeaveOff_brch            & !output :branch phenology flag, [-]
  )
  call PrintInfo('beg '//subname)
  IF(DayLenthCurrent.GE.DayLenthPrev)THEN
    Hours4Leafout_brch(NB,NZ)=Hours4LenthenPhotoPeriod_brch(NB,NZ)

    IF(Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ) &
      .OR.(ALAT.GT.0.0_r8 .AND. I.EQ.173) &     !northern hemisphere
      .OR.(ALAT.LT.0.0_r8 .AND. I.EQ.355))THEN  !southern hemisphere
      Hours4LeafOff_brch(NB,NZ)   = 0.0_r8
      doPlantLeaveOff_brch(NB,NZ) = iTrue
    ENDIF
  ENDIF
  !
  !   CALCULATE EVERGREEN PHENOLOGY DURINGTopRootLayer_pftSHORTENINGTopRootLayer_pftPHOTOPERIODS
  !
  !   Hours4Leafout_brch,Hours4LeafOff_brch=leafout,leafoff hours
  !   Hours4ShortenPhotoPeriod_brch=hourly counter for shortening photoperiods
  !   EnablePlantLeafOut_brch=flag for enabling leafout:0=enable,1=disable
  !   ALAT=latitude
  !
  IF(DayLenthCurrent.LT.DayLenthPrev)THEN
    Hours4LeafOff_brch(NB,NZ)=Hours4ShortenPhotoPeriod_brch(NB,NZ)
    IF(Hours4LeafOff_brch(NB,NZ).GE.HourReq4LeafOff_brch(NB,NZ) &
      .OR.(ALAT.GT.0.0_r8 .AND. I.EQ.355) &     !northern hemisphere
      .OR.(ALAT.LT.0.0_r8 .AND. I.EQ.173))THEN  !southern hemisphere
      Hours4Leafout_brch(NB,NZ)  = 0.0_r8
      EnablePlantLeafOut_brch(NB,NZ) = iTrue
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate  
  end subroutine CropEvergreenPhenology

!----------------------------------------------------------------------------------------------------
  subroutine ColdDeciduousBranchPhenology(I,J,NB,NZ)
  !   CALCULATE WINTER DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayer_pftHOURS ABOVE
  !   SPECIFIED TEMPERATURE DURINGTopRootLayer_pftLENGTHENINGTopRootLayer_pftPHOTOPERIODS

  implicit none
  integer, intent(in) :: I   !day
  integer, intent(in) :: J   !hour
  integer, intent(in) :: NB  !branch id
  integer, intent(in) :: NZ  !pft id
  character(len=*), parameter :: subname='ColdDeciduousBranchPhenology'

  associate(                                                         &
    DayLenthPrev             => plt_site%DayLenthPrev               ,& !input  :daylength of previous day, [h]
    DayLenthCurrent          => plt_site%DayLenthCurrent            ,& !input  :current daylength of the grid, [h]
    HourReq4LeafOff_brch     => plt_pheno%HourReq4LeafOff_brch      ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantCalendar_brch      => plt_pheno%iPlantCalendar_brch       ,& !input  :plant growth stage, [-]
    TCGroth_pft              => plt_pheno%TCGroth_pft               ,& !input  :canopy growth temperature, [oC]
    TC4LeafOff_pft           => plt_pheno%TC4LeafOff_pft            ,& !input  :threshold temperature for autumn leafoff/hardening, [oC]
    TC4LeafOut_pft           => plt_pheno%TC4LeafOut_pft            ,& !input  :threshold temperature for spring leafout/dehardening, [oC]
    TCChill4Seed_pft         => plt_pheno%TCChill4Seed_pft          ,& !input  :temperature below which seed set is adversely affected, [oC]
    HourReq4LeafOut_brch     => plt_pheno%HourReq4LeafOut_brch      ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    ALAT                     => plt_site%ALAT                       ,& !input  :latitude, [degrees north]
    EnablePlantLeafOut_brch  => plt_pheno%EnablePlantLeafOut_brch   ,& !inoput :branch phenology flag, [-]
    Hours4LeafOff_brch       => plt_pheno%Hours4LeafOff_brch        ,& !inoput :cold requirement for autumn leafoff/hardening, [h]
    Hours4Leafout_brch       => plt_pheno%Hours4Leafout_brch        ,& !inoput :heat requirement for spring leafout/dehardening, [h]
    doPlantLeaveOff_brch     => plt_pheno%doPlantLeaveOff_brch       & !output :branch phenology flag, [-]
  )

  call PrintInfo('beg '//subname)
  
  IF((DayLenthCurrent.GE.DayLenthPrev .OR. (DayLenthCurrent.LT.DayLenthPrev .AND. &
    Hours4LeafOff_brch(NB,NZ).LT.HourReq4LeafOff_brch(NB,NZ))) &
    .AND.EnablePlantLeafOut_brch(NB,NZ).EQ.iTrue)THEN

    !favorable temperature for leafout
    IF(TCGroth_pft(NZ).GE.TC4LeafOut_pft(NZ))THEN
      Hours4Leafout_brch(NB,NZ)=Hours4Leafout_brch(NB,NZ)+1.0_r8
    ENDIF
    !not sufficient leafout hour
    IF(Hours4Leafout_brch(NB,NZ).LT.HourReq4LeafOut_brch(NB,NZ))THEN
      !unfavorable cold temperature
      IF(TCGroth_pft(NZ).LT.TCChill4Seed_pft(NZ))THEN
        Hours4Leafout_brch(NB,NZ)=AZMAX1(Hours4Leafout_brch(NB,NZ)-1.0_r8)
      ENDIF
    ENDIF
    IF(Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ) &
      .OR.(ALAT.GT.0.0_r8 .AND.I.EQ.173) &     !northern hemisphere
      .OR.(ALAT.LT.0.0_r8 .AND.I.EQ.355))THEN  !southern hemisphere
      Hours4LeafOff_brch(NB,NZ)=0.0_r8
    ENDIF
  ENDIF

  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0 .OR. &
    (DayLenthCurrent.LT.DayLenthPrev .AND. DayLenthCurrent.LT.12.0_r8))THEN
    !has flower or (day length is decreasing, and current daytime is short)
    doPlantLeaveOff_brch(NB,NZ)=iTrue
  ENDIF
  !
  !     CALCULATE WINTER DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayer_pftHOURS BELOW
  !     SPECIFIED TEMPERATURE DURINGTopRootLayer_pftSHORTENINGTopRootLayer_pftPHOTOPERIODS
  !
  !     DayLenthPrev,DLYN=daylength of previous,current day
  !     Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
  !     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
  !     EnablePlantLeafOut_brch,doPlantLeaveOff_brch=flag for enabling leafout,leafoff:0=enable,1=disable
  !     TCGroth_pft,TCZ,TCChill4Seed_pft=canopy temp,leafout threshold temp,chilling temp
  !     ALAT=latitude
  !     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
  !
  IF(DayLenthCurrent.LT.DayLenthPrev .AND. doPlantLeaveOff_brch(NB,NZ).EQ.iTrue &
    .AND.iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)THEN
    IF(TCGroth_pft(NZ).LE.TC4LeafOff_pft(NZ))THEN
      Hours4LeafOff_brch(NB,NZ)=Hours4LeafOff_brch(NB,NZ)+1.0_r8
    ENDIF

    IF(Hours4LeafOff_brch(NB,NZ).GE.HourReq4LeafOff_brch(NB,NZ).AND.&
      EnablePlantLeafOut_brch(NB,NZ).EQ.iFalse)THEN
      Hours4Leafout_brch(NB,NZ)  = 0.0_r8
      EnablePlantLeafOut_brch(NB,NZ) = iTrue
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine ColdDeciduousBranchPhenology

!----------------------------------------------------------------------------------------------------
  subroutine live_branch_phenology(I,J,NB,NZ,TFNP,WFNG,OFNG)
  implicit none
  integer,  intent(in) :: I,J
  integer,  intent(in) :: NB  !plant branch id
  integer,  intent(in) :: NZ  !plant species id
  real(r8), intent(in) :: TFNP
  real(r8), intent(in) :: WFNG
  real(r8), intent(in) :: OFNG

  character(len=*), parameter :: subname='live_branch_phenology'

  integer :: kk
  logical :: DayLensShortenChk     !true when day length is decreasing
  
! begin_execution
  associate(                                                                           &
    DayLenthCurrent                   => plt_site%DayLenthCurrent                     ,& !input  :current daylength of the grid, [h]
    DayLenthPrev                      => plt_site%DayLenthPrev                        ,& !input  :daylength of previous day, [h]
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft               ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch                 ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    HourReq4LeafOff_brch              => plt_pheno%HourReq4LeafOff_brch               ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    MainBranchNum_pft                 => plt_morph%MainBranchNum_pft                  ,& !input  :id number of main branch,[-]
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch                ,& !inoput :plant growth stage, [-]
    Hours4Leafout_brch                => plt_pheno%Hours4Leafout_brch                 ,& !inoput :heat requirement for spring leafout/dehardening, [h]
    EnablePlantLeafOut_brch           => plt_pheno%EnablePlantLeafOut_brch            ,& !output :branch phenology flag, [-]
    doSenescence_brch                 => plt_pheno%doSenescence_brch                  ,& !output :branch phenology flag, [-]
    doInitLeafOut_brch                => plt_pheno%doInitLeafOut_brch                  & !output :branch phenology flag, [-]
  )

  call PrintInfo('beg '//subname)
  IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).EQ.0)THEN
    !plant emergence
    iPlantCalendar_brch(ipltcal_Emerge,NB,NZ) = I
    doInitLeafOut_brch(NB,NZ)                 = iFalse
    EnablePlantLeafOut_brch(NB,NZ)            = iTrue
    Hours4Leafout_brch(NB,NZ)                 = 0.5_r8*Hours4Leafout_brch(MainBranchNum_pft(NZ),NZ)
  ENDIF
  !
  ! CALCULATE NODE INITIATION AND LEAF APPEARANCE RATES
  ! FROM TEMPERATURE FUNCTION CALCULATED IN 'UPTAKE' AND
  ! RATES AT 25C ENTERED IN 'READQ' EXCEPT WHEN DORMANT
  !  
  IF(iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen .OR. Hours4LeafOff_brch(NB,NZ).LT.HourReq4LeafOff_brch(NB,NZ))THEN
    call UpdateBranchNodeNumber(I,J,NB,NZ,TFNP,WFNG,OFNG)
    doSenescence_brch(NB,NZ)=itrue
  ELSE
    doSenescence_brch(NB,NZ)=ifalse
  ENDIF

  ! REPRODUCTIVE GROWTH STAGES ADVANCE WHEN THRESHOLD NUMBER
  ! OF NODES HAVE BEEN INITIATED. FIRST DETERMINE PHOTOPERIOD
  ! AND TEMPERATURE EFFECTS ON FINAL VEG NODE NUMBER FROM
  ! NUMBER OF INITIATED NODES

  !
  DayLensShortenChk=DayLenthCurrent.LT.DayLenthPrev

  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).EQ.0)THEN    !->2
    !Main shoot or parent shoot leaf production
    call InitiateBranchFloral(I,J,NB,NZ,DayLensShortenChk)
    !
  ELSEIF(iPlantCalendar_brch(ipltcal_Jointing,NB,NZ).EQ.0)THEN  !->3
    !Stem jointing in phenological development marks the onset of rapid stem elongation in cereals 
    !like maize or wheat, when internodes begin to lengthen and the first node becomes visible or 
    !palpable above the soil surface.
    !
    !Tiller production
    call BranchStemJointing(I,J,NB,NZ,DayLensShortenChk)
    !
    !   ipltcal_Elongation,=mid stem elongation
    !
  ELSEIF(iPlantCalendar_brch(ipltcal_Elongation,NB,NZ).EQ.0)THEN !->4
    !Stem elongation and booting
    !Stem elongation in phenological development is the phase when cereal stems 
    !rapidly lengthen due to internode expansion, transitioning from tillering to booting/heading
!    if(NZ==1)write(683,*)(plt_site%iYearCurrent*1000+I)*10+J/24.,NB,'elongation'    
    call BranchStemElongation(I,J,NB,NZ,DayLensShortenChk)

  ELSEIF(iPlantCalendar_brch(ipltcal_heading,NB,NZ).EQ.0)THEN !->5
    !“heading” is the stage when the plant’s inflorescence (flower head) becomes visible and emerges from the enclosing leaves or sheath.
    !
    call BranchFlowerHeading(I,J,NB,NZ,DayLensShortenChk) 
    !!
  ELSEIF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).EQ.0)THEN !->6
    !Anthesis in phenological development is the stage when flowers are fully open and functional, allowing pollination to occur.
    !Anthesis stage or flowering marks the beginning of grain setting and filling

    call BranchAnthesis(I,J,NB,NZ,DayLensShortenChk) 
!
  ELSEIF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN !->7
    !
!    if(NZ==1)write(683,*)(plt_site%iYearCurrent*1000+I)*10+J/24.,NB,'bgseedfill'    
    call BeginBranchSeedFill(I,J,NB,NZ,DayLensShortenChk) 
!
!   END SEED NUMBER SET PERIOD
!
!   end date setting for final seed number
!
  ELSEIF(iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).EQ.0)THEN !->8

!    if(NZ==1)write(683,*)(plt_site%iYearCurrent*1000+I)*10+J/24.,NB,'seedmass'      
    IF(CheckBranchNodeState(ipltcal_SetSeedNumber,NB,NZ,I,J))THEN
      iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ)=I
    ENDIF
  !
  !   END SEED SIZE SET PERIOD
  !
  !   end of setting max seed size
  !
  ELSEIF(iPlantCalendar_brch(ipltcal_SetSeedMass,NB,NZ).EQ.0)THEN  !->9->10
!    if(NZ==1)write(683,*)(plt_site%iYearCurrent*1000+I)*10+J/24.,NB,'endseedfill'    
    IF(CheckBranchNodeState(ipltcal_SetSeedMass,NB,NZ,I,J))THEN
      iPlantCalendar_brch(ipltcal_SetSeedMass,NB,NZ)=I
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine live_branch_phenology
!----------------------------------------------------------------------------------------------------
  subroutine CalcPhenolEnvfactor(I,J,NZ,TFNP,WFNG,OFNG)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8),intent(out) :: TFNP              !temperature function for phenology (25 oC =1 )
  real(r8),intent(out) :: WFNG              !water stress effect on phenology
  real(r8),intent(out) :: OFNG              !oxygen stress
  character(len=*), parameter :: subname='CalcPhenolEnvfactor'
  real(r8) :: ACTV
  real(r8) :: RTK
  real(r8) :: STK,TKCO

  associate(                                                                           &
    iPlantPhenolPattern_pft           => plt_pheno%iPlantPhenolPattern_pft            ,& !input  :plant growth habit: annual or perennial,[-]  
    TKGroth_pft                       => plt_pheno%TKGroth_pft                        ,& !input  :canopy growth temperature, [K]
    PSICanopy_pft                     => plt_ew%PSICanopy_pft                         ,& !input  :canopy total water potential, [Mpa]        
    PlantO2Stress_pft                 => plt_pheno%PlantO2Stress_pft                  ,& !input  :plant O2 stress indicator, [-]       
    CanPhenoMoistStress_pft           => plt_pheno%CanPhenoMoistStress_pft            ,& !output :moisture stress for plant phenology development,[-]
    CanPhenoTempStress_pft            => plt_pheno%CanPhenoTempStress_pft             ,& !output :temperature stress for plant phenology development,[-]
    TempOffset_pft                    => plt_pheno%TempOffset_pft                      & !input  :adjustment of Arhhenius curves for plant thermal acclimation, [oC]    
  )
  call PrintInfo('beg '//subname)
  !
  TKCO = TKGroth_pft(NZ)+TempOffset_pft(NZ)
  TFNP = calc_leave_grow_tempf(TKCO)
  IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
    WFNG = EXP(0.025_r8*AMAX1(PSICanopy_pft(NZ),-1000._r8))
    OFNG = SQRT(PlantO2Stress_pft(NZ))
  ELSE
    WFNG=1._r8
    OFNG=1._r8
  ENDIF
  CanPhenoMoistStress_pft(NZ) = WFNG*OFNG
  CanPhenoTempStress_pft(NZ)  = TFNP
  call PrintInfo('end '//subname)
  end associate
  end subroutine CalcPhenolEnvfactor

!----------------------------------------------------------------------------------------------------
  subroutine UpdateBranchNodeNumber(I,J,NB,NZ,TFNP,WFNG,OFNG)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8),intent(in) :: TFNP              !temperature function for phenology (25 oC =1 )
  real(r8),intent(in) :: WFNG              !water stress effect on phenology
  real(r8),intent(in) :: OFNG              !oxygen stress
  character(len=*), parameter :: subname='UpdateBranchNodeNumber'  
  real(r8) :: NodeInitRate      !rates of node initiation
  real(r8) :: LeafAppearRate    !leaf appearance

  real(r8) :: HourlyNodeNumNormByMatgrp_brch

  associate(                                                                           &
    MatureGroup_pft                   => plt_pheno%MatureGroup_pft                    ,& !input  :acclimated plant maturity group, [-]    
    NumOfLeaves_brch                  => plt_morph%NumOfLeaves_brch                   ,& !inoput :leaf number, [-]    
    ShootNodeNumAtAnthesis_brch       => plt_morph%ShootNodeNumAtAnthesis_brch        ,& !input  :shoot node number at anthesis, [-]    
    RefNodeInitRate_pft               => plt_pheno%RefNodeInitRate_pft                ,& !input  :rate of node initiation at 25 oC, [h-1]    
    ShootNodeNumAtInitFloral_brch     => plt_morph%ShootNodeNumAtInitFloral_brch      ,& !input  :shoot node number at floral initiation, [-]    
    RateRefLeafAppearance_pft         => plt_pheno%RateRefLeafAppearance_pft          ,& !input  :rate of leaf initiation at 25 oC, [h-1]
    iPlantPhenolPattern_pft           => plt_pheno%iPlantPhenolPattern_pft            ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch                ,& !inoput :plant growth stage, [-]    
    ShootNodeNum_brch                 => plt_morph%ShootNodeNum_brch                  ,& !inoput :shoot node number, [-]    
    RNodeInitiate_pft                 => plt_rbgc%RNodeInitiate_pft                   ,& !inoput :node initiation rate, [h-1]
    RLeafAppear_pft                   => plt_rbgc%RLeafAppear_pft                     ,& !inoput :leaf appearing rate, [h-1]
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch  ,& !inoput :normalized node number during reproductive growth stages, [-]
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch      ,& !inoput :normalized node number during vegetative growth stages, [-]
    NodeNumNormByMatgrp_brch          => plt_pheno%NodeNumNormByMatgrp_brch           ,& !output :normalized node number during vegetative growth stages, [-]
    dReproNodeNumNormByMatG_brch      => plt_pheno%dReproNodeNumNormByMatG_brch       ,& !output :gain in normalized node number during reproductive growth stages, [h-1]
    ReprodNodeNumNormByMatrgrp_brch   => plt_pheno%ReprodNodeNumNormByMatrgrp_brch     & !output :normalized node number during reproductive growth stages, [-]    
  )
  call PrintInfo('beg '//subname)

  NodeInitRate   = AZMAX1(RefNodeInitRate_pft(NZ)*TFNP)
  LeafAppearRate = AZMAX1(RateRefLeafAppearance_pft(NZ)*TFNP)
  !
  !   NODE INITIATION AND LEAF APPEARANCE RATES SLOWED BY LOW TURGOR
  !
  !for leaf/node growth, only annual plants depends on moisture
  IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
    IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).EQ.0)THEN
      NodeInitRate   = NodeInitRate*OFNG
      LeafAppearRate = LeafAppearRate*OFNG
    ENDIF

    IF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).EQ.0)THEN
      NodeInitRate   = NodeInitRate*WFNG
      LeafAppearRate = LeafAppearRate*WFNG
    ENDIF

  ELSE 
    !do nothing for perrennial plants at the moment, Mar 25, 2026
    !perennial plants tend to maintain stronger antioxidant defenses and show delayed senescence compared to annuals    
    !      IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).EQ.0)THEN
    !        OFNG           = PlantO2Stress_pft(NZ)**0.5_r8
    !        NodeInitRate   = NodeInitRate*OFNG
    !        LeafAppearRate = LeafAppearRate*OFNG
    !      ENDIF      
  ENDIF
  !
  !   ACCUMULATE NODE INITIATION AND LEAF APPEARANCE RATES
  !   INTO TOTAL NUMBER OF NODES AND LEAVES
  !
  !   PSTG,NumOfLeaves_brch=number of nodes initiated,leaves appeared
  RNodeInitiate_pft(NZ)    = RNodeInitiate_pft(NZ)+NodeInitRate
  RLeafAppear_pft(NZ)      = RLeafAppear_pft(NZ)+LeafAppearRate
  ShootNodeNum_brch(NB,NZ) = ShootNodeNum_brch(NB,NZ)+NodeInitRate
  NumOfLeaves_brch(NB,NZ)  = NumOfLeaves_brch(NB,NZ)+LeafAppearRate
  !
  !   USE TOTAL NUMBER OF NODES TO CALCULATE PROGRESSION THROUGH
  !   VEGETATIVE AND REPRODUCTIVE GROWTH STAGES. THIS PROGRESSION
  !   IS USED TO SET START AND END DATES FOR GROWTH STAGES BELOW
  !
  !   NodeNumNormByMatgrp_brch=vegetative node number normalized for maturity group
  !   MatureGroup_pft=node number required for floral initiation
  
  
  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)THEN
    NodeNumNormByMatgrp_brch(NB,NZ)      = (ShootNodeNum_brch(NB,NZ)-ShootNodeNumAtInitFloral_brch(NB,NZ))/MatureGroup_pft(NZ)    
    HourlyNodeNumNormByMatgrp_brch       = NodeInitRate/(MatureGroup_pft(NZ)*GrothStageNorm4VegetaPheno)
    TotalNodeNumNormByMatgrp_brch(NB,NZ) = TotalNodeNumNormByMatgrp_brch(NB,NZ)+HourlyNodeNumNormByMatgrp_brch
  ENDIF

  IF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).NE.0)THEN
    ReprodNodeNumNormByMatrgrp_brch(NB,NZ)   = (ShootNodeNum_brch(NB,NZ)-ShootNodeNumAtAnthesis_brch(NB,NZ))/MatureGroup_pft(NZ)
    dReproNodeNumNormByMatG_brch(NB,NZ)      = NodeInitRate/(MatureGroup_pft(NZ)*GrothStageNorm4ReprodPheno)
    TotReproNodeNumNormByMatrgrp_brch(NB,NZ) = TotReproNodeNumNormByMatrgrp_brch(NB,NZ)+dReproNodeNumNormByMatG_brch(NB,NZ)
  ENDIF    
  call PrintInfo('end '//subname)
  end associate
  END SUBROUTINE UpdateBranchNodeNumber

!----------------------------------------------------------------------------------------------------
  function CheckBranchNodeState(ipltcal,NB,NZ,I,J)result(NodeNumChk)
  !
  !Description:
  !
  implicit none
  integer, intent(in) :: ipltcal   !phenological stage
  integer, intent(in) :: NB,NZ,I,J
  character(len=*), parameter :: subname='CheckBranchNodeState'
  logical :: NodeNumChk

  associate(                                                                        &
    MatureGroup_pft                 => plt_pheno%MatureGroup_pft                  , & !input  :acclimated plant maturity group, [-]  
    MatureGroup_brch                => plt_pheno%MatureGroup_brch                 , & !input  :plant maturity group, [-]  
    ShootNodeNum_brch               => plt_morph%ShootNodeNum_brch                , & !input  :shoot node number, [-]  
    NumOfLeaves_brch                => plt_morph%NumOfLeaves_brch                 , & !input  :leaf number, [-]  
    ShootNodeNumAtInitFloral_brch   => plt_morph%ShootNodeNumAtInitFloral_brch    , & !input  :shoot node number at floral initiation, [-]    
    ReprodNodeNumNormByMatrgrp_brch => plt_pheno%ReprodNodeNumNormByMatrgrp_brch  , & !input  :normalized node number during reproductive growth stages, [-]  
    NodeNumNormByMatgrp_brch        => plt_pheno%NodeNumNormByMatgrp_brch           & !input  :normalized node number during vegetative growth stages, [-]
  )
  call PrintInfo('beg '//subname)

  select case(ipltcal)
    case(ipltcal_InitFloral) !2
    !ensures the plant is tall enough 
    !the transition of the shoot apical meristem from vegetative to reproductive growth
    NodeNumChk  = ShootNodeNum_brch(NB,NZ).GT.MatureGroup_pft(NZ)+ShootNodeNumAtInitFloral_brch(NB,NZ)    

    case(ipltcal_Jointing) !3
    NodeNumChk = NodeNumNormByMatgrp_brch(NB,NZ).GT.0.25_r8*GrothStageNorm4VegetaPheno

    case(ipltcal_Elongation) !4
    NodeNumChk = NodeNumNormByMatgrp_brch(NB,NZ).GT.0.50_r8*GrothStageNorm4VegetaPheno

    case(ipltcal_heading) !5
    NodeNumChk = NodeNumNormByMatgrp_brch(NB,NZ).GT.1.00_r8*GrothStageNorm4VegetaPheno

    case(ipltcal_Anthesis)   !6 
    NodeNumChk = NumOfLeaves_brch(NB,NZ).GT.ShootNodeNumAtInitFloral_brch(NB,NZ)

    case(ipltcal_BeginSeedFill) !7
    NodeNumChk = ReprodNodeNumNormByMatrgrp_brch(NB,NZ).GT.0.50_r8*GrothStageNorm4ReprodPheno

    case(ipltcal_SetSeedNumber) !8
    NodeNumChk=ReprodNodeNumNormByMatrgrp_brch(NB,NZ).GT.1.00_r8*GrothStageNorm4ReprodPheno

    case(ipltcal_SetSeedMass) !9
    NodeNumChk=ReprodNodeNumNormByMatrgrp_brch(NB,NZ).GT.1.50_r8*GrothStageNorm4ReprodPheno

    case default
    NodeNumChk=.False.
  end select

  call PrintInfo('end '//subname)
  end associate
  end function CheckBranchNodeState
  
!----------------------------------------------------------------------------------------------------
  subroutine BeginBranchSeedFill(I,J,NB,NZ,DayLensShortenChk) 
  implicit none
  integer, intent(in) :: I
  integer, intent(in) :: J
  integer, intent(in) :: NB
  integer, intent(in) :: NZ
  logical, intent(in) :: DayLensShortenChk      !true when day length is decreasing
  character(len=*), parameter :: subname='BeginBranchSeedFill'
  logical :: LeafOffChk,PhenoChk1,PhenoChk2

  associate(                                                                       &
    Hours4LeafOff_brch              => plt_pheno%Hours4LeafOff_brch               ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    iPlantPhenolPattern_pft         => plt_pheno%iPlantPhenolPattern_pft          ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantPhenolType_pft            => plt_pheno%iPlantPhenolType_pft             ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    HourReq4LeafOff_brch            => plt_pheno%HourReq4LeafOff_brch             ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    EnablePlantLeafOut_brch         => plt_pheno%EnablePlantLeafOut_brch          ,& !input  :branch phenology flag, [-]
    iPlantPhotoperiodType_pft       => plt_pheno%iPlantPhotoperiodType_pft        ,& !input  :photoperiod type (neutral, long day, short day),[-]
    iPlantCalendar_brch             => plt_pheno%iPlantCalendar_brch               & !output :plant growth stage, [-]
  )
  call PrintInfo('beg '//subname)

  LeafOffChk = Hours4LeafOff_brch(NB,NZ).GT.HourReq4LeafOff_brch(NB,NZ)

  PhenoChk1  = is_cold_deciduos(iPlantPhenolType_pft(NZ)) .AND. DayLensShortenChk &
    .AND.iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .AND. iPlantPhotoperiodType_pft(NZ).NE.iphotop_short 

  PhenoChk2 = iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_annual

  IF(CheckBranchNodeState(ipltcal_BeginSeedFill,NB,NZ,I,J) .OR.(PhenoChk1.AND.EnablePlantLeafOut_brch(NB,NZ).EQ.iFalse.AND.LeafOffChk) &
                .OR. PhenoChk2.AND.EnablePlantLeafOut_brch(NB,NZ).EQ.iFalse.AND.LeafOffChk)THEN
      iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ)=I
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine BeginBranchSeedFill

!----------------------------------------------------------------------------------------------------
  subroutine BranchAnthesis(I,J,NB,NZ,DayLensShortenChk) 
  !
  !Description:
  !   ANTHESIS OCCURS WHEN THE NUMBER OF LEAVES THAT HAVE APPEARED
  !   EQUALS THE NUMBER OF NODES INITIATED WHEN THE FINAL VEGETATIVE
  !   NODE NUMBER WAS SET ABOVE
  !
  implicit none
  integer, intent(in) :: I
  integer, intent(in) :: J
  integer, intent(in) :: NB
  integer, intent(in) :: NZ
  logical, intent(in) :: DayLensShortenChk
  character(len=*), parameter :: subname='BranchAnthesis'
  logical :: LeafOffChk,PerennialPhenoChk,AnnualPhenoChk,CalChk

  associate(                                                           &
    Hours4LeafOff_brch        => plt_pheno%Hours4LeafOff_brch         ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    iPlantPhenolPattern_pft   => plt_pheno%iPlantPhenolPattern_pft    ,& !input  :plant growth habit: annual or perennial,[-]
    HourReq4LeafOff_brch      => plt_pheno%HourReq4LeafOff_brch       ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft       ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iPlantPhotoperiodType_pft => plt_pheno%iPlantPhotoperiodType_pft  ,& !input  :photoperiod type (neutral, long day, short day),[-]
    EnablePlantLeafOut_brch       => plt_pheno%EnablePlantLeafOut_brch        ,& !input  :branch phenology flag, [-]
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft          ,& !input  :id number of main branch,[-]
    ShootNodeNum_brch         => plt_morph%ShootNodeNum_brch          ,& !input  :shoot node number, [-]
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch        ,& !inoput :plant growth stage, [-]
    ShootNodeNumAtAnthesis_brch => plt_morph%ShootNodeNumAtAnthesis_brch   & !output :shoot node number at anthesis, [-]
  )
  call PrintInfo('beg '//subname)

  LeafOffChk = Hours4LeafOff_brch(NB,NZ).GT.HourReq4LeafOff_brch(NB,NZ) .AND. EnablePlantLeafOut_brch(NB,NZ).EQ.iFalse

  PerennialPhenoChk  = is_cold_deciduos(iPlantPhenolType_pft(NZ)) .AND. DayLensShortenChk &
    .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .AND. iPlantPhotoperiodType_pft(NZ).NE.iphotop_short 

  AnnualPhenoChk = iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_annual
  CalChk         = iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .AND. iPlantCalendar_brch(ipltcal_heading,NB,NZ).NE.0

  IF(CheckBranchNodeState(ipltcal_Anthesis,NB,NZ,I,J) .OR. CalChk &
    .OR. (PerennialPhenoChk  .AND. LeafOffChk) .OR. (AnnualPhenoChk .AND. LeafOffChk))THEN

    IF(NB.EQ.MainBranchNum_pft(NZ) .OR. iPlantCalendar_brch(ipltcal_Anthesis,MainBranchNum_pft(NZ),NZ).NE.0)THEN
      iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ) = I
      ShootNodeNumAtAnthesis_brch(NB,NZ)            = ShootNodeNum_brch(NB,NZ)
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate  
  end subroutine BranchAnthesis

!----------------------------------------------------------------------------------------------------
  subroutine BranchFlowerHeading(I,J,NB,NZ,DayLensShortenChk)
  !
  !Description:
  !For flowering plants, the flowering heading stage (often simply called heading) is the critical 
  !transition from vegetative growth to reproductive development. In this stage:
  !Stalk is characterized by a process called "bolting" or rapid elongation.
  !Leave senescence is triggered, and nutrients are redirected from the leaves to the reproductive organs.
  !The plant's reproductive organs become visible
  implicit none
  integer, intent(in) :: I
  integer, intent(in) :: J
  integer, intent(in) :: NB
  integer, intent(in) :: NZ
  logical, intent(in) :: DayLensShortenChk   !true when day length is decreasing
  character(len=*), parameter :: subname='BranchFlowerHeading'
  logical :: LeafOffChk,PerennialPhenoChk,AnnualPhenoChk

  associate(                                                           &
    Hours4LeafOff_brch        => plt_pheno%Hours4LeafOff_brch         ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    iPlantPhenolPattern_pft   => plt_pheno%iPlantPhenolPattern_pft    ,& !input  :plant growth habit: annual or perennial,[-]
    EnablePlantLeafOut_brch       => plt_pheno%EnablePlantLeafOut_brch        ,& !input  :branch phenology flag, [-]
    HourReq4LeafOff_brch      => plt_pheno%HourReq4LeafOff_brch       ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantPhotoperiodType_pft => plt_pheno%iPlantPhotoperiodType_pft  ,& !input  :photoperiod type (neutral, long day, short day),[-]
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft       ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch         & !output :plant growth stage, [-]
  )

  call PrintInfo('beg '//subname)

  !check for leave senescence
  LeafOffChk = EnablePlantLeafOut_brch(NB,NZ).EQ.iFalse .AND. Hours4LeafOff_brch(NB,NZ).GT.HourReq4LeafOff_brch(NB,NZ) 

  !Perennial (Cold deciduous, or cold drought deciduous) .and. not short photoperiod and day length is decreasing
  PerennialPhenoChk  = is_cold_deciduos(iPlantPhenolType_pft(NZ)) .AND. DayLensShortenChk &
    .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .AND. iPlantPhotoperiodType_pft(NZ).NE.iphotop_short 

  !Annual drought deciduous  
  AnnualPhenoChk=iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_annual  

  IF(CheckBranchNodeState(ipltcal_heading,NB,NZ,I,J) .OR. ((PerennialPhenoChk .OR. AnnualPhenoChk)  .AND.  LeafOffChk))THEN
    iPlantCalendar_brch(ipltcal_heading,NB,NZ)=I
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine BranchFlowerHeading

!----------------------------------------------------------------------------------------------------
  subroutine BranchStemElongation(I,J,NB,NZ,DayLensShortenChk)
  !
  !During stem elongation stage (often called "jointing" in grasses), 
  !the plant shifts its energy from making a bushy base to building height.
  !The Flag (last) Leaf Emergence, the plant feels less "leafy" and more "stiff" or "reedy."
  implicit none
  integer, intent(in) :: I
  integer, intent(in) :: J
  integer, intent(in) :: NB
  integer, intent(in) :: NZ
  logical, intent(in) :: DayLensShortenChk   !is day length decreasing
  character(len=*), parameter :: subname='BranchStemElongation'
  logical :: NodeNumChk,LeafOffChk,PerennialPhenoChk,AnnualPhenoChk

  associate(                                                               &
    NodeNumNormByMatgrp_brch    => plt_pheno%NodeNumNormByMatgrp_brch     ,& !input  :normalized node number during vegetative growth stages, [-]
    iPlantPhenolType_pft        => plt_pheno%iPlantPhenolType_pft         ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    Hours4LeafOff_brch          => plt_pheno%Hours4LeafOff_brch           ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    iPlantPhenolPattern_pft     => plt_pheno%iPlantPhenolPattern_pft      ,& !input  :plant growth habit: annual or perennial,[-]
    EnablePlantLeafOut_brch     => plt_pheno%EnablePlantLeafOut_brch      ,& !input  :branch phenology flag, [-]
    ShootNodeNum_brch           => plt_morph%ShootNodeNum_brch            ,& !input  :shoot node number, [-]
    iPlantDevelopPattern_pft    => plt_pheno%iPlantDevelopPattern_pft     ,& !input  :plant growth habit (determinate or indeterminate),[-]
    iPlantPhotoperiodType_pft   => plt_pheno%iPlantPhotoperiodType_pft    ,& !input  :photoperiod type (neutral, long day, short day),[-]
    HourReq4LeafOff_brch        => plt_pheno%HourReq4LeafOff_brch         ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantCalendar_brch         => plt_pheno%iPlantCalendar_brch          ,& !output :plant growth stage, [-]
    LeafNumberAtFloralInit_brch => plt_pheno%LeafNumberAtFloralInit_brch   & !output :leaf number at floral initiation, [-]
  )
  call PrintInfo('beg '//subname)

  LeafOffChk = EnablePlantLeafOut_brch(NB,NZ).EQ.iFalse .AND. Hours4LeafOff_brch(NB,NZ).GT.HourReq4LeafOff_brch(NB,NZ)

  !perrennial (cold deciduous or drought deciduous)
  PerennialPhenoChk  = is_cold_deciduos(iPlantPhenolType_pft(NZ)) .AND. DayLensShortenChk & 
    .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .AND. iPlantPhotoperiodType_pft(NZ).NE.iphotop_short 

  !annual drought deciduous  
  AnnualPhenoChk= iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu 

  IF(CheckBranchNodeState(ipltcal_Elongation,NB,NZ,I,J) .OR. ((PerennialPhenoChk .OR. AnnualPhenoChk) .AND.  LeafOffChk))THEN
    iPlantCalendar_brch(ipltcal_Elongation,NB,NZ)=I

    !indeterminate
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantDevelopPattern_pft(NZ).NE.ideterminate)THEN
      LeafNumberAtFloralInit_brch(NB,NZ)=ShootNodeNum_brch(NB,NZ)
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate  
  end subroutine BranchStemElongation

!----------------------------------------------------------------------------------------------------
  subroutine BranchStemJointing(I,J,NB,NZ,DayLensShortenChk)
!   STEM ELONGATION
  implicit none
  integer, intent(in) :: I
  integer, intent(in) :: J
  integer, intent(in) :: NB
  integer, intent(in) :: NZ
  logical, intent(in) :: DayLensShortenChk
  character(len=*), parameter :: subname='BranchStemJointing'
  logical :: PerennialPhenoChk,AnnualPhenoChk,LeafOffChk

  associate(                                                           &
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft       ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iPlantPhotoperiodType_pft => plt_pheno%iPlantPhotoperiodType_pft  ,& !input  :photoperiod type (neutral, long day, short day),[-]
    iPlantPhenolPattern_pft   => plt_pheno%iPlantPhenolPattern_pft    ,& !input  :plant growth habit: annual or perennial,[-]
    EnablePlantLeafOut_brch   => plt_pheno%EnablePlantLeafOut_brch    ,& !input  :branch phenology flag, [-]
    Hours4LeafOff_brch        => plt_pheno%Hours4LeafOff_brch         ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    HourReq4LeafOff_brch      => plt_pheno%HourReq4LeafOff_brch       ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch         & !output :plant growth stage, [-]
  )
  call PrintInfo('beg '//subname)

  PerennialPhenoChk  = is_cold_deciduos(iPlantPhenolType_pft(NZ)) .AND. DayLensShortenChk &
    .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .AND. iPlantPhotoperiodType_pft(NZ).NE.iphotop_short 

  AnnualPhenoChk  = iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_annual

  LeafOffChk = EnablePlantLeafOut_brch(NB,NZ).EQ.iFalse .AND. Hours4LeafOff_brch(NB,NZ).GT.HourReq4LeafOff_brch(NB,NZ)

  IF(CheckBranchNodeState(ipltcal_Jointing,NB,NZ,I,J) &
    .OR. (PerennialPhenoChk .AND. LeafOffChk) .OR. (AnnualPhenoChk .AND. LeafOffChk))THEN
    
    iPlantCalendar_brch(ipltcal_Jointing,NB,NZ)=I
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine BranchStemJointing

!----------------------------------------------------------------------------------------------------
  subroutine InitiateBranchFloral(I,J,NB,NZ,DayLensShortenChk)
  !
  !Description:
  implicit none
  integer, intent(in) :: I
  integer, intent(in) :: J
  integer, intent(in) :: NB
  integer, intent(in) :: NZ
  logical, intent(in) :: DayLensShortenChk  !true when day length is decreasing
  character(len=*), parameter :: subname='InitiateBranchFloral'
  logical :: NodeNumChk,LeafOutChk,PlantDayChk,CanopyHeightChk,PerennialPhenoDecidChk
  logical :: CropPhenoChk,PhotoPrdChk,PhenolCheck,PhotoLongCheck,PhotoShortCheck,AnnualFloralCheck
  real(r8) :: PPD

  associate(                                                                     &
    ShootNodeNum_brch               => plt_morph%ShootNodeNum_brch              ,& !input  :shoot node number, [-]
    ZERO                            => plt_site%ZERO                            ,& !input  :threshold zero for numerical stability, [-]
    iPlantPhenolPattern_pft         => plt_pheno%iPlantPhenolPattern_pft        ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantDevelopPattern_pft        => plt_pheno%iPlantDevelopPattern_pft       ,& !input  :plant growth habit (determinate or indeterminate),[-]
    iYearCurrent                    => plt_site%iYearCurrent                    ,& !input  :current year,[-]
    SnowDepth                       => plt_ew%SnowDepth                         ,& !input  :snowpack depth, [m]
    CanopyHeightLive_pft                => plt_morph%CanopyHeightLive_pft               ,& !input  :canopy height, [m]
    iDayPlanting_pft                => plt_distb%iDayPlanting_pft               ,& !input  :day of planting,[-]
    PhotoPeriodSens_pft             => plt_pheno%PhotoPeriodSens_pft            ,& !input  :difference between current and critical daylengths used to calculate phenological progress, [h]
    HourReq4LeafOut_brch            => plt_pheno%HourReq4LeafOut_brch           ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    DayLenthCurrent                 => plt_site%DayLenthCurrent                 ,& !input  :current daylength of the grid, [h]
    DayLenthPrev                    => plt_site%DayLenthPrev                    ,& !input  :daylength of previous day, [h]
    iPlantPhenolType_pft            => plt_pheno%iPlantPhenolType_pft           ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iYearPlanting_pft               => plt_distb%iYearPlanting_pft              ,& !input  :year of planting,[-]
    Hours4Leafout_brch              => plt_pheno%Hours4Leafout_brch             ,& !input  :heat requirement for spring leafout/dehardening, [h]
    CriticPhotoPeriod_pft           => plt_pheno%CriticPhotoPeriod_pft          ,& !input  :critical daylength for phenological progress, [h]
    iPlantPhotoperiodType_pft       => plt_pheno%iPlantPhotoperiodType_pft      ,& !input  :photoperiod type (neutral, long day, short day),[-]
    ShootNodeNumAtInitFloral_brch   => plt_morph%ShootNodeNumAtInitFloral_brch  ,& !inoput :shoot node number at floral initiation, [-]
    LeafNumberAtFloralInit_brch     => plt_pheno%LeafNumberAtFloralInit_brch    ,& !output :leaf number at floral initiation, [-]
    iPlantCalendar_brch             => plt_pheno%iPlantCalendar_brch             & !output :plant growth stage, [-]
  )

  call PrintInfo('beg '//subname)
  
  LeafOutChk      = Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ)
  PlantDayChk     = I.GE.iDayPlanting_pft(NZ) .AND. iYearCurrent.EQ.iYearPlanting_pft(NZ) .AND. DayLenthCurrent.GT.DayLenthPrev
  CanopyHeightChk = CanopyHeightLive_pft(NZ).GE.SnowDepth-ZERO
  PerennialPhenoDecidChk   = iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .AND. is_cold_deciduos(iPlantPhenolType_pft(NZ))

  !Annual crop plants (i.e. those seeded by human) are set as evergreen, if it is self-seeding, then 
  !iPlantPhenolType_pft should be set according to koppen climate zone.  

  CropPhenoChk      = iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen
  PhenolCheck       = ((PerennialPhenoDecidChk .OR. CropPhenoChk) .AND. CanopyHeightChk .AND. DayLensShortenChk)
  AnnualFloralCheck = (CheckBranchNodeState(ipltcal_InitFloral,NB,NZ,I,J) .AND. (LeafOutChk .OR. PlantDayChk))
  !
  IF( AnnualFloralCheck .OR. PhenolCheck)THEN
    !
    !     FINAL VEGETATIVE NODE NUMBER DEPENDS ON PHOTOPERIOD FROM 'DAY'
    !     AND ON MATURITY GROUP, CRITICAL PHOTOPERIOD AND PHOTOPERIOD
    !     SENSITIVITY ENTERED IN 'READQ'
    !
    !     iPlantPhotoperiodType_pft=photoperiod type from PFT file
    !     PPD=photoperiod sensitivity
    !     CriticPhotoPeriod_pft=critical photoperiod from PFT file
    !     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
    !     VSTGX=node number on date of floral initiation
    !
    IF(iPlantPhotoperiodType_pft(NZ).EQ.iphotop_neutral)THEN
      PPD=0.0_r8
    ELSE
      PPD=AZMAX1(CriticPhotoPeriod_pft(NZ)-DayLenthCurrent)
      IF(iPlantPhotoperiodType_pft(NZ).EQ.iphotop_short .AND. DayLenthCurrent.GE.DayLenthPrev)PPD=0.0_r8
    ENDIF

    PhotoShortCheck = (iPlantPhotoperiodType_pft(NZ).EQ.iphotop_short.AND. PPD.GT.PhotoPeriodSens_pft(NZ))
    PhotoLongCheck  = (iPlantPhotoperiodType_pft(NZ).EQ.iphotop_long .AND. PPD.LT.PhotoPeriodSens_pft(NZ))
    PhotoPrdChk     = iPlantPhotoperiodType_pft(NZ).EQ.iphotop_neutral .OR. PhotoShortCheck .OR. PHotoLongCheck
    
    IF( PhotoPrdChk .OR. PhenolCheck)THEN
      iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ) = I
      ShootNodeNumAtInitFloral_brch(NB,NZ)          = ShootNodeNum_brch(NB,NZ)

      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantDevelopPattern_pft(NZ).EQ.ideterminate)THEN
        LeafNumberAtFloralInit_brch(NB,NZ)=ShootNodeNum_brch(NB,NZ)
      ENDIF
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine InitiateBranchFloral
  ![tail]
end module PlantPhenolMod

