module PlantBranchMod

! Description:
! module for plant biological transformations
  use minimathmod,   only: isclose, safe_adb, AZMAX1
  use data_kind_mod, only: r8 => DAT_KIND_R8,yearIJ_type
  use EcoSIMCtrlMod, only: etimer, guard_floor_dying_stand, CPOOL_DIV_FLOOR_dying
  use PlantBalMod
  use ElmIDMod
  use DebugToolMod
  use EcosimConst
  use PlantBGCPars
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
  use PhotoSynsMod
  use PlantMathFuncMod
  use NoduleBGCMod  
  use LitterFallMod
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: GrowOneBranch

  integer, parameter :: ibrch_leaf   = 1
  integer, parameter :: ibrch_petole = 2
  integer, parameter :: ibrch_stalk  = 3
  integer, parameter :: ibrch_resrv  = 4
  integer, parameter :: ibrch_husk   = 5
  integer, parameter :: ibrch_ear    = 6
  integer, parameter :: ibrch_grain  = 7
  logical :: stop_flag
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine GrowOneBranch(yearIJ,NB,NZ,TFN6_vr,CanopyHeight_copy,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,&
    CPRTW,TFN5,WaterStress4Groth,Stomata_Stress,TurgEff4LeafPetolExpansion,TurgEff4CanopyResp,&
    GrothPART2LeafPetole,BegRemoblize)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in)  :: NB,NZ
  REAL(R8), INTENT(IN) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: CanopyHeight_copy(JP1)
  real(r8), intent(in) :: CNLFW,CPLFW   !leaf
  real(r8), intent(in) :: CNSHW,CPSHW   !petole
  real(r8), intent(in) :: CNRTW,CPRTW   !root/stalk
  real(r8), intent(in) :: TFN5,WaterStress4Groth
  real(r8), intent(in) :: Stomata_Stress
  real(r8), intent(in) :: TurgEff4LeafPetolExpansion
  real(r8), intent(in) :: TurgEff4CanopyResp
  real(r8), intent(out) :: GrothPART2LeafPetole !new growth allocated to leaf and PetolSheth
  integer, intent(out) :: BegRemoblize !phenological eligibility of this branch
  character(len=*), parameter :: subname='GrowOneBranch'
  integer :: I,J,K
  real(r8) :: YCO2Gro_brch                     !respiraiton ratio of canopy growth, 1./YCO2Gro_brch, total biomass required to give 1 unit of CO2 respiration
  integer  :: MinNodeID,NE  
  REAL(R8) :: CH2O3(pltpar%MaxNodesPerBranch1),CH2O4(pltpar%MaxNodesPerBranch1)
  integer  :: LRemob_brch
  REAL(R8) :: PART(pltpar%NumOfPlantMorphUnits)
  real(r8) :: ALLOCL
  REAL(R8) :: CNPG
  real(r8) :: CCE
  real(r8) :: CCC,CNC,CPC
  real(r8) :: EtiolationCoeff
  real(r8) :: Growth_brch(NumPlantChemElms,pltpar%NumOfPlantMorphUnits)
  real(r8) :: RCO2NonstC_brch
  REAL(R8) :: RCO2Maint_brch
  real(r8) :: RMxess_brch              !excess maintenance respiration, drives remobilization & senescence
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: SSL
  real(r8) :: CNLFM,CPLFM
  real(r8) :: CNLFX,CPLFX,CNSHX,CPSHX
  real(r8) :: RNonstC4Groth_brch           !nonstructural carbon flux for growth
  real(r8) :: RCO2NonstC4Nassim_brch
  real(r8) :: ShootStructE_brch(1:NumPlantChemElms)
  real(r8) :: DMLFB
  real(r8) :: DMSHB
  real(r8) :: CNLFB
  real(r8) :: CPLFB  
  real(r8) :: CNSHB,CPSHB,arr(12)
  real(r8) :: RCCE(NumPlantChemElms)
  character(len=32)::header
! begin_execution
  associate(                                                            &
    LeafPetoNonstElmConc_brch  => plt_biom%LeafPetoNonstElmConc_brch   ,& !input  :branch nonstructural C concentration, [g d-2]
    LeafStrutElms_brch         => plt_biom%LeafStrutElms_brch          ,& !input  :branch leaf structural element mass, [g d-2]
    LeafElmntNode_brch         => plt_biom%LeafElmntNode_brch          ,& !input :leaf element, [g d-2]    
    PetolShethStrutElms_brch   => plt_biom%PetolShethStrutElms_brch    ,& !input  :branch sheath structural element, [g d-2]
    iPlantPhotosynsType_pft    => plt_photo%iPlantPhotosynsType_pft    ,& !input  :plant photosynthetic type (C3 or C4),[-]
    iPlantRootProfile_pft      => plt_pheno%iPlantRootProfile_pft      ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlantCalendar_brch        => plt_pheno%iPlantCalendar_brch        ,& !input  :plant growth stage, [-]
    iPlantTurnoverPattern_pft  => plt_pheno%iPlantTurnoverPattern_pft  ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    iPlant2ndGrothPattern_pft  => plt_pheno%iPlant2ndGrothPattern_pft  ,& !input  :plant expression of secondary growth, [-]                
    ZERO                       => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]
    MainBranchNum_pft          => plt_morph%MainBranchNum_pft          ,& !input  :number of main branch,[-]
    isPlantBranchAlive_brch    => plt_pheno%isPlantBranchAlive_brch    ,& !inoput :flag to detect branch death, [-]
    canopy_growth_pft          => plt_rbgc%canopy_growth_pft           ,& !inoput :canopy structural C growth rate, [gC d-2 h-1]
    KHiestGroLeafNode_brch     => plt_pheno%KHiestGroLeafNode_brch     ,& !inoput :leaf growth stage counter, [-]    
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft          ,& !input  :number of branches, [-]
    SapwoodBiomassC_brch       => plt_biom%SapwoodBiomassC_brch        ,& !input  :branch sapwood C, [gC d-2]
    CanopySapwoodC_pft         => plt_biom%CanopySapwoodC_pft          ,& !output :canopy sapwood C, [gC d-2]
    CanopyLeafSheathC_pft      => plt_biom%CanopyLeafSheathC_pft       ,& !output :canopy leaf + sheath C, [gC d-2]
    CanopyLeafSheathC_brch     => plt_biom%CanopyLeafSheathC_brch       & !output :plant branch leaf + sheath C, [g d-2]
  )
  call PrintInfo('beg '//subname)
  stop_flag=.false.
  I=yearIJ%I;J=yearIJ%J

  ! Define outputs even for dead branches; these are consumed by GrowOnePlant.
  BegRemoblize = ifalse
  GrothPART2LeafPetole = 0._r8

  IF(isPlantBranchAlive_brch(NB,NZ).EQ.iTrue)THEN
    !correct the higest growing leaf node
    DO K=MaxNodesPerBranch1,0,-1
      if(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.1.e-6_r8)then
        KHiestGroLeafNode_brch(NB,NZ)=K
        exit
      endif
    ENDDO
    KHiestGroLeafNode_brch(NB,NZ)=MAX(KHiestGroLeafNode_brch(NB,NZ),1)
    
    call CalcPartitionCoeff(I,J,NB,NZ,PART, &
      BegRemoblize=BegRemoblize,LRemob_brch=LRemob_brch)

    IF(NB.EQ.MainBranchNum_pft(NZ))THEN
      GrothPART2LeafPetole=PART(ibrch_leaf)+PART(ibrch_petole)
    ENDIF
    !
    call UpdateBranchAllometry(I,J,NZ,NB,PART,CNLFW,CNRTW,CNSHW,CPLFW,&
      CPRTW,CPSHW,ShootStructE_brch,YCO2Gro_brch,CNLFM,CPLFM,&
      CNSHX,CPSHX,CNLFX,CPLFX,DMLFB,DMSHB,CNLFB,CPLFB,CNSHB,CPSHB)
    !
    !   GROSS PRIMARY PRODUCTIVITY
    !
    call UpdatePhotosynthates(I,J,NB,NZ,TFN6_vr,YCO2Gro_brch,CNLFM,CPLFM,CNSHX,CPSHX &
      ,CNLFX,CPLFX,ShootStructE_brch(ielmn),TFN5,WaterStress4Groth,Stomata_Stress,TurgEff4CanopyResp &
      ,CH2O3,CH2O4,CNPG,RMxess_brch,RNonstC4Groth_brch)
    !
    !   TRANSFER OF C4 FIXATION PRODUCTS FROM NON-STRUCTURAL POOLS
    !   IN MESOPHYLL TO THOSE IN BUNDLE SHEATH, DECARBOXYLATION
    !   OF C4 FIXATION PRODUCTS IN BUNDLE SHEATH, LEAKAGE OF DECARBOXYLATION
    !   PRODUCTS BACK TO MESOPHYLL IN C4 PLANTS
    !
    !   iPlantPhotosynsType_pft=photosynthesis type:3=C3,4=C4
    !
    IF(iPlantPhotosynsType_pft(NZ).EQ.ic4_photo)THEN
      call C4PhotoProductTransfer(I,J,NZ,NB,CH2O3,CH2O4)
    ENDIF
    !
    canopy_growth_pft(NZ) = canopy_growth_pft(NZ)+RNonstC4Groth_brch
    !
    call BranchBiomAllocate(I,J,NB,NZ,PART,RNonstC4Groth_brch,&
      DMLFB,DMSHB,CNLFB,CPLFB,CNSHB,CPSHB,ZPLFD,CNPG,Growth_brch,EtiolationCoeff,MinNodeID)
    !
    CALL GrowLeavesOnBranch(I,J,NZ,NB,MinNodeID,Growth_brch(:,ibrch_leaf),EtiolationCoeff,TurgEff4LeafPetolExpansion,ALLOCL)
    !
    !     DISTRIBUTE SHEATH OR PetolSheth GROWTH AMONG CURRENTLY GROWING NODES
    !
    CALL GrowPetolShethOnBranch(yearIJ,NZ,NB,MinNodeID,Growth_brch(:,ibrch_petole),EtiolationCoeff,TurgEff4LeafPetolExpansion,ALLOCL)
    !
    !   DISTRIBUTE STALK GROWTH AMONG CURRENTLY GROWING NODES
    !
    call GrowStalkOnBranch(I,J,NZ,NB,Growth_brch(:,ibrch_stalk),EtiolationCoeff)
    !
    !   RECOVERY OF REMOBILIZABLE N,P DURING REMOBILIZATION DEPENDS
    !   ON SHOOT NON-STRUCTURAL C:N:P
    !
    !   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
    !   CNKI,CPKI=nonstruct N,P inhibition constant on growth (g N,P g-1 C)
    !   RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !   RCCZ,RCCY=min,max fractions for shoot C recycling
    !   RCCX,RCCQ=max fractions for shoot N,P recycling
    !   iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
    !
    !plant already emerged and has meaningful leaves for the branch
    IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).NE.0 .AND. LeafPetoNonstElmConc_brch(ielmc,NB,NZ).GT.ZERO)THEN
      CCC=AZMAX1(AMIN1(1.0_r8 &
        ,LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmn,NB,NZ)+LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CNKI) &
        ,LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmp,NB,NZ)+LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CPKI)))
      CNC=AZMAX1(AMIN1(1.0_r8 &
        ,LeafPetoNonstElmConc_brch(ielmc,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmc,NB,NZ)+LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/CNKI)))
      CPC=AZMAX1(AMIN1(1.0_r8 &
        ,LeafPetoNonstElmConc_brch(ielmc,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmc,NB,NZ)+LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/CPKI)))
    ELSE
      CCC=0._r8
      CNC=0._r8
      CPC=0._r8
    ENDIF
    RCCC = RCCZ(iPlantRootProfile_pft(NZ))+CCC*RCCY(iPlantRootProfile_pft(NZ))
    RCCN = CNC*RCCX(iPlantRootProfile_pft(NZ))
    RCCP = CPC*RCCQ(iPlantRootProfile_pft(NZ))

    RCCE(ielmc)=RCCC
    RCCE(ielmn)=(RCCN+(1.0_r8-RCCN)*RCCC)
    RCCE(ielmp)=(RCCP+(1.0_r8-RCCP)*RCCC)

    !
    !       WITHDRAW REMOBILIZABLE C,N,P FROM LOWEST NODE AFTER
    !       MAXIMUM NODE NUMBER OF 25 IS REACHED
    !
    !       doSenescence_brch=PFT branch senescence flag
    !       KHiestGroLeafNode_brch=integer of most recent leaf number
    !       fTCanopyGroth_pft=temperature function for canopy growth
    !       XRLA=rate of leaf appearance at 25 oC (h-1)
    !
    call SenescenceBranch(NZ,NB,RCCE)
    !
    call RemobilizeBranch(I,J,NZ,NB,BegRemoblize,LRemob_brch,RCCC,RCCN,RCCP,RMxess_brch)
    !
    !   DEATH IF MAIN STALK OF TREE DIES
    !
    !   iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+PetolSheth,2=none,3=between 1,2
    !   iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
    !   KMinGroingLeafNodeNum,KHiestGroLeafNode_brch=integer of lowest,highest leaf number currently growing
    !
    IF(iPlantTurnoverPattern_pft(NZ).NE.0 .AND. is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)) &
      .AND. isPlantBranchAlive_brch(MainBranchNum_pft(NZ),NZ).EQ.iFalse)then
      isPlantBranchAlive_brch(NB,NZ)=iFalse
    endif
    !     REMOBILIZE EXCESS LEAF STRUCTURAL N,P
    call WithdrawNutBranchLeaves(I,J,NB,NZ,CNLFB,CPLFB)
    !
    call GrainFillOnBranch(I,J,NB,NZ,Growth_brch(:,ibrch_grain),Growth_brch(ielmc,ibrch_stalk))
    call ResetBranchPhenology(I,J,NB,NZ)
    !
    ! Seasonal turnover can clear leaves or remove stalk biomass. Build geometry
    ! once, after these removals, so layer totals and heights use surviving organs.
    call AllocateLeaf2CanopyLayers(I,J,NB,NZ,CanopyHeight_copy)
    ! Keep angular profiles consistent even when a seasonal reset occurs at night.
    call LeafClassAllocation(NB,NZ)
    ! Refresh the biomass denominators used by the following reserve transfers.
    CanopyLeafSheathC_brch(NB,NZ)=AZMAX1(LeafStrutElms_brch(ielmc,NB,NZ)+PetolShethStrutElms_brch(ielmc,NB,NZ))
    CanopyLeafSheathC_pft(NZ)=SUM(CanopyLeafSheathC_brch(1:NumOfBranches_pft(NZ),NZ))
    CanopySapwoodC_pft(NZ)=SUM(SapwoodBiomassC_brch(1:NumOfBranches_pft(NZ),NZ))
    !
    call BranchElmntTransfer(I,J,NB,NZ,BegRemoblize,WaterStress4Groth,TurgEff4CanopyResp)
    !
    !   CANOPY N2 FIXATION (CYANOBACTERIA)
    !
    call CanopyNoduleBiochemistry(I,J,NZ,NB,TFN5,WaterStress4Groth)

  ENDIF
  
  call PrintInfo('end '//subname)
  end associate
  end subroutine GrowOneBranch

!----------------------------------------------------------------------------------------------------
  subroutine WithdrawNutBranchLeaves(I,J,NB,NZ,CNLFB,CPLFB)
  !
  !Description:
  !withdraw excessive nutrients from leaves
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in):: CNLFB,CPLFB
  integer :: KK,NE,K
  integer :: KMinGroingLeafNodeNum
  real(r8) :: ZPOOLD,XFRN1,XFRP1
  real(r8) :: CPOOLT
! --- Jing Tao --- floored copy of CPOOLT used ONLY as the division denominator below.
  real(r8) :: CPOOLT_div
  real(r8) :: XFRE(1:NumPlantChemElms)
  real(r8) :: PPOOLD
  character(len=*), parameter :: subname='WithdrawNutBranchLeaves'

  associate(                                                     &
    ZERO4Groth_pft         => plt_biom%ZERO4Groth_pft           ,& !input  :threshold zero for plang growth calculation, [-]
    KHiestGroLeafNode_brch => plt_pheno%KHiestGroLeafNode_brch  ,& !input  :leaf growth stage counter, [-]
    rProteinC2LeafN_pft    => plt_allom%rProteinC2LeafN_pft     ,& !input  :Plant leaf protein C to N ratio [g protein C (g leaf N)-1]
    rProteinC2LeafP_pft    => plt_allom%rProteinC2LeafP_pft     ,& !input  :Plant leaf protein C to P ratio [g protein C (g leaf P)-1]
    LeafElmntNode_brch     => plt_biom%LeafElmntNode_brch       ,& !inoput :leaf element, [g d-2]
    LeafStrutElms_brch     => plt_biom%LeafStrutElms_brch       ,& !inoput :branch leaf structural element mass, [g d-2]
    LeafProteinC_node      => plt_biom%LeafProteinC_node        ,& !inoput :layer leaf protein C, [g d-2]
    CanopyNonstElms_brch   => plt_biom%CanopyNonstElms_brch      & !inoput :branch nonstructural element, [g d-2]
  )
  call DebugPrint('beg '//subname//' NZ',NZ)

  KMinGroingLeafNodeNum=MAX(0,KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+1)
  
  D495: DO KK=KMinGroingLeafNodeNum,KHiestGroLeafNode_brch(NB,NZ)
    K=pMOD(KK,MaxNodesPerBranch1)
    IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.0.0_r8)THEN
      CPOOLT=LeafElmntNode_brch(ielmc,K,NB,NZ)+CanopyNonstElms_brch(ielmc,NB,NZ)
      IF(CPOOLT.GT.ZERO4Groth_pft(NZ))THEN
        ZPOOLD=LeafElmntNode_brch(ielmn,K,NB,NZ)*CanopyNonstElms_brch(ielmc,NB,NZ) &
          -CanopyNonstElms_brch(ielmn,NB,NZ)*LeafElmntNode_brch(ielmc,K,NB,NZ)
        PPOOLD=LeafElmntNode_brch(ielmp,K,NB,NZ)*CanopyNonstElms_brch(ielmc,NB,NZ) &
          -CanopyNonstElms_brch(ielmp,NB,NZ)*LeafElmntNode_brch(ielmc,K,NB,NZ)
          
! --- Jing Tao ---
! the guard above (CPOOLT > ZERO4Groth_pft) does NOT bound this division once the
! stand is dying, because ZERO4Groth_pft = ZERO*PlantPopuLive_pft collapses with the
! population and the test degenerates to "CPOOLT > 0". Floor the DENOMINATOR ONLY,
! leaving CPOOLT itself untouched for every other use. Default-off: when the switch
! is .false. CPOOLT_div IS CPOOLT and the arithmetic is bit-for-bit unchanged.
        CPOOLT_div = CPOOLT
        if(guard_floor_dying_stand) CPOOLT_div = AMAX1(CPOOLT, CPOOL_DIV_FLOOR_dying)
        XFRN1=AZMAX1(AMIN1(1.0E-03_r8*ZPOOLD/CPOOLT_div,LeafElmntNode_brch(ielmn,K,NB,NZ) &
          -ZPLFM*CNLFB*LeafElmntNode_brch(ielmc,K,NB,NZ)))            
        XFRP1=AZMAX1(AMIN1(1.0E-03_r8*PPOOLD/CPOOLT_div,LeafElmntNode_brch(ielmp,K,NB,NZ) &
          -ZPLFM*CPLFB*LeafElmntNode_brch(ielmc,K,NB,NZ)))
          
        XFRE(ielmn)=AMAX1(XFRN1,10.0_r8*XFRP1)
        XFRE(ielmp)=AMAX1(XFRP1,0.10_r8*XFRN1)
        DO NE=2,NumPlantChemElms
          call ExchFluxLimiter(LeafElmntNode_brch(NE,K,NB,NZ) ,CanopyNonstElms_brch(NE,NB,NZ),XFRE(NE))
          LeafElmntNode_brch(NE,K,NB,NZ) = LeafElmntNode_brch(NE,K,NB,NZ)-XFRE(NE)
          LeafStrutElms_brch(NE,NB,NZ)   = LeafStrutElms_brch(NE,NB,NZ)-XFRE(NE)
          CanopyNonstElms_brch(NE,NB,NZ) = CanopyNonstElms_brch(NE,NB,NZ)+XFRE(NE)
        ENDDO
        LeafProteinC_node(K,NB,NZ)=AZMAX1(LeafProteinC_node(K,NB,NZ)-&
          AMAX1(XFRE(ielmn)*rProteinC2LeafN_pft(NZ),XFRE(ielmp)*rProteinC2LeafP_pft(NZ)))
      ENDIF
    ENDIF
  ENDDO D495
  call PrintInfo('end '//subname)
  end associate
  end subroutine WithdrawNutBranchLeaves

!----------------------------------------------------------------------------------------------------
  subroutine BranchBiomAllocate(I,J,NB,NZ,PART,RNonstC4Groth_brch,&
    DMLFB,DMSHB,CNLFB,CPLFB,CNSHB,CPSHB,ZPLFD,CNPG,Growth_brch, &
    EtiolationCoeff,MinNodeID)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: PART(pltpar%NumOfPlantMorphUnits)
  real(r8), intent(in) :: RNonstC4Groth_brch                  !non-structural C used for branch growth, [gC h-1 d-2]
  real(r8), intent(in) :: DMLFB,DMSHB,CNLFB
  real(r8), intent(in) :: CPLFB,CNSHB,CPSHB,ZPLFD,CNPG
  real(r8), intent(out) :: Growth_brch(NumPlantChemElms,pltpar%NumOfPlantMorphUnits)
  integer , intent(out) :: MinNodeID
  real(r8), intent(out) :: EtiolationCoeff       !nutrient effect on etiolation

  character(len=*), parameter :: subname='BranchBiomAllocate'
  real(r8) :: CCE
  integer  :: NE
  
  associate(                                                          &
    LeafPetoNonstElmConc_brch => plt_biom%LeafPetoNonstElmConc_brch  ,& !input  :branch nonstructural C concentration, [g d-2]
    EarBiomGrowthYld_pft      => plt_allom%EarBiomGrowthYld_pft      ,& !input  :ear growth yield, [gC gC-1]
    rNCEar_pft                => plt_allom%rNCEar_pft                ,& !input  :ear N:C ratio, [gN gC-1]
    rPCReserve_pft            => plt_allom%rPCReserve_pft            ,& !input  :reserve P:C ratio, [gP gC-1]
    rPCHusk_pft               => plt_allom%rPCHusk_pft               ,& !input  :husk P:C ratio, [gP gC-1]
    rPCStalk_pft              => plt_allom%rPCStalk_pft              ,& !input  :stalk P:C ratio, [g g-1]
    rNCStalk_pft              => plt_allom%rNCStalk_pft              ,& !input  :stalk N:C ratio, [gN gC-1]
    StalkBiomGrowthYld_pft    => plt_allom%StalkBiomGrowthYld_pft    ,& !input  :stalk growth yield, [gC gC-1]
    HuskBiomGrowthYld_pft     => plt_allom%HuskBiomGrowthYld_pft     ,& !input  :husk growth yield, [gC gC-1]
    rNCHusk_pft               => plt_allom%rNCHusk_pft               ,& !input  :husk N:C ratio, [gN gC-1]
    rPCEar_pft                => plt_allom%rPCEar_pft                ,& !input  :ear P:C ratio, [gP gC-1]
    ReserveBiomGrowthYld_pft  => plt_allom%ReserveBiomGrowthYld_pft  ,& !input  :reserve growth yield, [gC gC-1]
    GrainBiomGrowthYld_pft    => plt_allom%GrainBiomGrowthYld_pft    ,& !input  :grain growth yield, [gC gC-1]
    SeedDepth_pft             => plt_morph%SeedDepth_pft             ,& !input  :seeding depth, [m]
    HypocotHeight_pft         => plt_morph%HypocotHeight_pft         ,& !input  :cotyledon height, [m]
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft         ,& !input  :number of main branch,[-]
    rNCReserve_pft            => plt_allom%rNCReserve_pft            ,& !input  :reserve N:C ratio, [gN gC-1]
    StalkStrutElms_brch       => plt_biom%StalkStrutElms_brch        ,& !inoput :branch stalk structural element mass, [g d-2]
    StalkRsrvElms_brch        => plt_biom%StalkRsrvElms_brch         ,& !inoput :branch reserve element mass, [g d-2]
    HuskStrutElms_brch        => plt_biom%HuskStrutElms_brch         ,& !inoput :branch husk structural element mass, [g d-2]
    EarStrutElms_brch         => plt_biom%EarStrutElms_brch          ,& !inoput :branch ear structural chemical element mass, [g d-2]
    PetolShethStrutElms_brch  => plt_biom%PetolShethStrutElms_brch   ,& !inoput :branch sheath structural element, [g d-2]
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch         ,& !inoput :branch leaf structural element mass, [g d-2]
    fNCLFW_brch               => plt_pheno%fNCLFW_brch               ,& !inoput :NC ratio of growing leaf, [gN/gC]
    fPCLFW_brch               => plt_pheno%fPCLFW_brch                & !inoput :PC ratio of growing leaf, [gP/gC]
  )
!
!   C,N,P GROWTH OF LEAF, SHEATH OR PetolSheth, STALK,
!   STALK RESERVES, REPRODUCTIVE ORGANS, GRAIN
!

  call DebugPrint('beg '//subname//' NZ',NZ)

  Growth_brch(ielmc,ibrch_leaf)   =RNonstC4Groth_brch*PART(ibrch_leaf)  *DMLFB       !leaf
  Growth_brch(ielmc,ibrch_petole) =RNonstC4Groth_brch*PART(ibrch_petole)*DMSHB       !PetolSheth
  Growth_brch(ielmc,ibrch_stalk)  =RNonstC4Groth_brch*PART(ibrch_stalk) *StalkBiomGrowthYld_pft(NZ)  !stalk
  Growth_brch(ielmc,ibrch_resrv)  =RNonstC4Groth_brch*PART(ibrch_resrv) *ReserveBiomGrowthYld_pft(NZ)  !reserve
  Growth_brch(ielmc,ibrch_husk)   =RNonstC4Groth_brch*PART(ibrch_husk)  *HuskBiomGrowthYld_pft(NZ)  !husk
  Growth_brch(ielmc,ibrch_ear)    =RNonstC4Groth_brch*PART(ibrch_ear)   *EarBiomGrowthYld_pft(NZ)  !ear
  Growth_brch(ielmc,ibrch_grain)  =RNonstC4Groth_brch*PART(ibrch_grain) *GrainBiomGrowthYld_pft(NZ)    !grain

  fNCLFW_brch(NB,NZ)              = CNLFB*(ZPLFM+ZPLFD*CNPG)
  Growth_brch(ielmn,ibrch_leaf)   = Growth_brch(ielmc,ibrch_leaf)  *fNCLFW_brch(NB,NZ)
  Growth_brch(ielmn,ibrch_petole) = Growth_brch(ielmc,ibrch_petole)*CNSHB
  Growth_brch(ielmn,ibrch_stalk)  = Growth_brch(ielmc,ibrch_stalk) *rNCStalk_pft(NZ)
  Growth_brch(ielmn,ibrch_resrv)  = Growth_brch(ielmc,ibrch_resrv) *rNCReserve_pft(NZ)
  Growth_brch(ielmn,ibrch_husk)   = Growth_brch(ielmc,ibrch_husk)  *rNCHusk_pft(NZ)
  Growth_brch(ielmn,ibrch_ear)    = Growth_brch(ielmc,ibrch_ear)   *rNCEar_pft(NZ)
  Growth_brch(ielmn,ibrch_grain)  = Growth_brch(ielmc,ibrch_grain) *rNCReserve_pft(NZ)

  fPCLFW_brch(NB,NZ)              = CPLFB*(ZPLFM+ZPLFD*CNPG)
  Growth_brch(ielmp,ibrch_leaf)   = Growth_brch(ielmc,ibrch_leaf)  *fPCLFW_brch(NB,NZ)
  Growth_brch(ielmp,ibrch_petole) = Growth_brch(ielmc,ibrch_petole)*CPSHB
  Growth_brch(ielmp,ibrch_stalk)  = Growth_brch(ielmc,ibrch_stalk) *rPCStalk_pft(NZ)
  Growth_brch(ielmp,ibrch_resrv)  = Growth_brch(ielmc,ibrch_resrv) *rPCReserve_pft(NZ)
  Growth_brch(ielmp,ibrch_husk)   = Growth_brch(ielmc,ibrch_husk)  *rPCHusk_pft(NZ)
  Growth_brch(ielmp,ibrch_ear)    = Growth_brch(ielmc,ibrch_ear)   *rPCEar_pft(NZ)
  Growth_brch(ielmp,ibrch_grain)  = Growth_brch(ielmc,ibrch_grain) *rPCReserve_pft(NZ)
  
  DO NE=1,NumPlantChemElms
    LeafStrutElms_brch(NE,NB,NZ)   = LeafStrutElms_brch(NE,NB,NZ)+Growth_brch(NE,ibrch_leaf)
    PetolShethStrutElms_brch(NE,NB,NZ) = PetolShethStrutElms_brch(NE,NB,NZ)+Growth_brch(NE,ibrch_petole)
    StalkStrutElms_brch(NE,NB,NZ)  = StalkStrutElms_brch(NE,NB,NZ)+Growth_brch(NE,ibrch_stalk)
    StalkRsrvElms_brch(NE,NB,NZ)   = StalkRsrvElms_brch(NE,NB,NZ)+Growth_brch(NE,ibrch_resrv)
    HuskStrutElms_brch(NE,NB,NZ)   = HuskStrutElms_brch(NE,NB,NZ)+Growth_brch(NE,ibrch_husk)
    EarStrutElms_brch(NE,NB,NZ)    = EarStrutElms_brch(NE,NB,NZ)+Growth_brch(NE,ibrch_ear)
  ENDDO

  !
  !   ETOLIATION
  !
  !   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
  !   CNKI,CPKI=nonstruct N,P inhibition constant on growth (g N,P g-1 C)
  !   ETOL=coefficient for etoliation effects on expansion,extension
  !
  CCE=AMIN1(safe_adb(LeafPetoNonstElmConc_brch(ielmn,NB,NZ),LeafPetoNonstElmConc_brch(ielmn,NB,NZ)+LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CNKI) &
    ,safe_adb(LeafPetoNonstElmConc_brch(ielmp,NB,NZ),LeafPetoNonstElmConc_brch(ielmp,NB,NZ)+LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CPKI))

  EtiolationCoeff = AMAX1(1.0_r8+CCE,1._r8)
  !
  !   DISTRIBUTE LEAF GROWTH AMONG CURRENTLY GROWING NODES
  !
  IF(NB.EQ.MainBranchNum_pft(NZ) .AND. HypocotHeight_pft(NZ).LE.SeedDepth_pft(NZ))THEN
    MinNodeID=0
  ELSE
    MinNodeID=1
  ENDIF
  call PrintInfo('end '//subname)
  end associate  
  end subroutine BranchBiomAllocate

!----------------------------------------------------------------------------------------------------
  subroutine UpdateBranchAllometry(I,J,NZ,NB,PART,CNLFW,CNRTW,CNSHW,CPLFW,&
    CPRTW,CPSHW,ShootStructE_brch,YCO2Gro_brch,CNLFM,CPLFM,&
    CNSHX,CPSHX,CNLFX,CPLFX,DMLFB,DMSHB,CNLFB,CPLFB,CNSHB,CPSHB)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: CNLFW,CPLFW   !leaf
  real(r8), intent(in) :: CNSHW,CPSHW   !petole
  real(r8), intent(in) :: CNRTW,CPRTW   !root/stalk
  real(r8), intent(in) :: PART(pltpar%NumOfPlantMorphUnits)
  real(r8), intent(out) :: ShootStructE_brch(NumPlantChemElms)    
  real(r8), intent(out) :: YCO2Gro_brch     !respiraiton ratio of canopy growth, 1./YCO2Gro_brch, total biomass required to give 1 unit of CO2 respiration
  real(r8), intent(out) :: CNLFM,CPLFM
  real(r8), intent(out) :: CNSHX,CPSHX
  real(r8), intent(out) :: CNLFX,CPLFX
  real(r8), intent(out) :: DMLFB,DMSHB,CNLFB,CPLFB,CNSHB,CPSHB

  character(len=*), parameter :: subname='UpdateBranchAllometry'
  real(r8)   :: DMSHT  !shoot C biomass yield, averaged across all canopy organs

  associate(                                                               &
    LeafBiomGrowthYld_pft       => plt_allom%LeafBiomGrowthYld_pft        ,& !input  :leaf growth yield, [g g-1]
    PetolShethBiomGrowthYld_pft => plt_allom%PetolShethBiomGrowthYld_pft  ,& !input  :sheath growth yield, [g g-1]
    EarBiomGrowthYld_pft        => plt_allom%EarBiomGrowthYld_pft         ,& !input  :ear growth yield, [gC gC-1]
    rNCEar_pft                  => plt_allom%rNCEar_pft                   ,& !input  :ear N:C ratio, [gN gC-1]
    rPCReserve_pft              => plt_allom%rPCReserve_pft               ,& !input  :reserve P:C ratio, [gP gC-1]
    rPCHusk_pft                 => plt_allom%rPCHusk_pft                  ,& !input  :husk P:C ratio, [gP gC-1]
    rPCStalk_pft                => plt_allom%rPCStalk_pft                 ,& !input  :stalk P:C ratio, [g g-1]
    rNCStalk_pft                => plt_allom%rNCStalk_pft                 ,& !input  :stalk N:C ratio, [gN gC-1]
    rNCReserve_pft              => plt_allom%rNCReserve_pft               ,& !input  :reserve N:C ratio, [gN gC-1]
    StalkBiomGrowthYld_pft      => plt_allom%StalkBiomGrowthYld_pft       ,& !input  :stalk growth yield, [gC gC-1]
    HuskBiomGrowthYld_pft       => plt_allom%HuskBiomGrowthYld_pft        ,& !input  :husk growth yield, [gC gC-1]
    rNCHusk_pft                 => plt_allom%rNCHusk_pft                  ,& !input  :husk N:C ratio, [gN gC-1]
    rPCEar_pft                  => plt_allom%rPCEar_pft                   ,& !input  :ear P:C ratio, [gP gC-1]
    ReserveBiomGrowthYld_pft    => plt_allom%ReserveBiomGrowthYld_pft     ,& !input  :reserve growth yield, [gC gC-1]
    GrainBiomGrowthYld_pft      => plt_allom%GrainBiomGrowthYld_pft       ,& !input  :grain growth yield, [gC gC-1]
    iPlantCalendar_brch         => plt_pheno%iPlantCalendar_brch          ,& !input  :plant growth stage, [-]
    PetolShethStrutElms_brch    => plt_biom%PetolShethStrutElms_brch      ,& !input  :branch sheath structural element, [g d-2]
    GrainStrutElms_brch         => plt_biom%GrainStrutElms_brch           ,& !input  :branch grain structural element mass, [g d-2]
    EarStrutElms_brch           => plt_biom%EarStrutElms_brch             ,& !input  :branch ear structural chemical element mass, [g d-2]
    HuskStrutElms_brch          => plt_biom%HuskStrutElms_brch            ,& !input  :branch husk structural element mass, [g d-2]
    LeafStrutElms_brch          => plt_biom%LeafStrutElms_brch            ,& !input  :branch leaf structural element mass, [g d-2]
    SapwoodBiomassC_brch        => plt_biom%SapwoodBiomassC_brch          ,& !input  :branch sapwood C, [gC d-2]
    RootBiomGrosYld_pft         => plt_allom%RootBiomGrosYld_pft           & !input  :root growth yield, [g g-1]
  )

  call DebugPrint('beg '//subname//' NZ',NZ)
!   SHOOT COEFFICIENTS FOR GROWTH RESPIRATION AND N,P CONTENTS
!   FROM GROWTH YIELDS ENTERED IN 'READQ', AND FROM PARTITIONING
!   COEFFICIENTS ABOVE
!
!   ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!   ZPLFD=1.0_r8-ZPLFM
!
  IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).NE.0)THEN
    DMLFB = LeafBiomGrowthYld_pft(NZ)
    DMSHB = PetolShethBiomGrowthYld_pft(NZ)
    CNLFB = CNLFW
    CNSHB = CNSHW
    CPLFB = CPLFW
    CPSHB = CPSHW
  ELSE
    DMLFB = RootBiomGrosYld_pft(NZ)
    DMSHB = RootBiomGrosYld_pft(NZ)
    CNLFB = CNRTW
    CNSHB = CNRTW
    CPLFB = CPRTW
    CPSHB = CPRTW
  ENDIF
    !part 1(leaf), (2) PetolSheth, (3) stalk, (4) reserve, (5) husk, (6) ear, (7) grain
  DMSHT = PART(ibrch_leaf)*DMLFB+PART(ibrch_petole)*DMSHB+PART(ibrch_stalk)*StalkBiomGrowthYld_pft(NZ) &
    +PART(ibrch_resrv)*ReserveBiomGrowthYld_pft(NZ)+PART(ibrch_husk)*HuskBiomGrowthYld_pft(NZ) &
    +PART(ibrch_ear)*EarBiomGrowthYld_pft(NZ)+PART(ibrch_grain)*GrainBiomGrowthYld_pft(NZ)

  YCO2Gro_brch = 1.0_r8-DMSHT
  CNLFM = PART(ibrch_leaf)*DMLFB*ZPLFM*CNLFB
  CPLFM = PART(ibrch_leaf)*DMLFB*ZPLFM*CPLFB
  CNLFX = PART(ibrch_leaf)*DMLFB*ZPLFD*CNLFB
  CPLFX = PART(ibrch_leaf)*DMLFB*ZPLFD*CPLFB
  CNSHX = PART(ibrch_petole)*DMSHB*CNSHB                               &
    +PART(ibrch_stalk)*StalkBiomGrowthYld_pft(NZ)*rNCStalk_pft(NZ)     &
    +PART(ibrch_resrv)*ReserveBiomGrowthYld_pft(NZ)*rNCReserve_pft(NZ) &
    +PART(ibrch_husk)*HuskBiomGrowthYld_pft(NZ)*rNCHusk_pft(NZ)        &
    +PART(ibrch_ear)*EarBiomGrowthYld_pft(NZ)*rNCEar_pft(NZ)           &
    +PART(ibrch_grain)*GrainBiomGrowthYld_pft(NZ)*rNCReserve_pft(NZ)

  CPSHX = PART(ibrch_petole)*DMSHB*CPSHB                               &
    +PART(ibrch_stalk)*StalkBiomGrowthYld_pft(NZ)*rPCStalk_pft(NZ)     &
    +PART(ibrch_resrv)*ReserveBiomGrowthYld_pft(NZ)*rPCReserve_pft(NZ) &
    +PART(ibrch_husk)*HuskBiomGrowthYld_pft(NZ)*rPCHusk_pft(NZ)        &
    +PART(ibrch_ear)*EarBiomGrowthYld_pft(NZ)*rPCEar_pft(NZ)           &
    +PART(ibrch_grain)*GrainBiomGrowthYld_pft(NZ)*rPCReserve_pft(NZ)
!
!   TOTAL SHOOT STRUCTURAL N MASS FOR MAINTENANCE RESPIRATION
!
!   ShootStructE_brch(ielmn)=shoot structural N mass
!   WTLFBN,WTSHBN,WTHSBN,WTEARN,WTFRBN=leaf,PetolSheth,husk,ear,grain N mass
!   rNCStalk_pft,SapwoodBiomassC_brch=stalk N:C,sapwood mass
!   iPlantCalendar_brch(10=date of physiological maturity
!
  ShootStructE_brch(ielmc) = AZMAX1(LeafStrutElms_brch(ielmc,NB,NZ)+PetolShethStrutElms_brch(ielmc,NB,NZ)+SapwoodBiomassC_brch(NB,NZ))
  ShootStructE_brch(ielmn) = AZMAX1(LeafStrutElms_brch(ielmn,NB,NZ)+PetolShethStrutElms_brch(ielmn,NB,NZ)+SapwoodBiomassC_brch(NB,NZ)*rNCStalk_pft(NZ))
  ShootStructE_brch(ielmp) = AZMAX1(LeafStrutElms_brch(ielmp,NB,NZ)+PetolShethStrutElms_brch(ielmp,NB,NZ)+SapwoodBiomassC_brch(NB,NZ)*rPCStalk_pft(NZ))
  
  IF(iPlantCalendar_brch(ipltcal_EndSeedFill,NB,NZ).EQ.0)THEN
    ShootStructE_brch(ielmn)=ShootStructE_brch(ielmn)+AZMAX1(HuskStrutElms_brch(ielmn,NB,NZ) &
      +EarStrutElms_brch(ielmn,NB,NZ)+GrainStrutElms_brch(ielmn,NB,NZ))

    ShootStructE_brch(ielmp)=ShootStructE_brch(ielmp)+AZMAX1(HuskStrutElms_brch(ielmp,NB,NZ) &
      +EarStrutElms_brch(ielmp,NB,NZ)+GrainStrutElms_brch(ielmp,NB,NZ))      
  ENDIF

  call PrintInfo('end '//subname)
  end associate
  end subroutine UpdateBranchAllometry  

!----------------------------------------------------------------------------------------------------
  subroutine CalcPartitionCoeff(I,J,NB,NZ,PART,BegRemoblize,LRemob_brch)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  integer, intent(out) :: BegRemoblize,LRemob_brch
  real(r8), intent(out):: PART(pltpar%NumOfPlantMorphUnits)

  character(len=*), parameter :: subname='CalcPartitionCoeff'
  REAL(R8) :: ARLFI
  integer :: N,LP
  real(r8) :: FPARTL
  real(r8) :: PARTS
  real(r8) :: PARTX
  real(r8) :: TOTAL
  logical :: check_perennial,check_annual,check_cold_deciduous,check_drought_deciduous
  real(r8) :: PSILY(0:3)
  real(r8), parameter :: FPART1=1.00_r8
  real(r8), parameter :: FPART2=0.40_r8

  associate(                                                                           &
    TdegCCanopy_pft                   => plt_ew%TdegCCanopy_pft                       ,& !input  :canopy temperature, [oC]
    PSICanopy_pft                     => plt_ew%PSICanopy_pft                         ,& !input  :canopy total water potential, [Mpa]
    SapwoodBiomassC_brch              => plt_biom%SapwoodBiomassC_brch                ,& !input  :branch sapwood C, [gC d-2]
    StalkRsrvElms_brch                => plt_biom%StalkRsrvElms_brch                  ,& !input  :branch reserve element mass, [g d-2]
    ZERO                              => plt_site%ZERO                                ,& !input  :threshold zero for numerical stability, [-]
    NU                                => plt_site%NU                                  ,& !input  :current soil surface layer number, [-]
    AREA3                             => plt_site%AREA3                               ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    HourReq4LeafOff_brch              => plt_pheno%HourReq4LeafOff_brch               ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch                ,& !input  :plant growth stage, [-]
    iPlantRootProfile_pft             => plt_pheno%iPlantRootProfile_pft              ,& !input  :plant growth type (vascular, non-vascular),[-]
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch  ,& !input  :normalized node number during reproductive growth stages, [-]
    iPlantDevelopPattern_pft          => plt_pheno%iPlantDevelopPattern_pft           ,& !input  :plant growth habit (determinate or indeterminate),[-]
    iPlantTurnoverPattern_pft         => plt_pheno%iPlantTurnoverPattern_pft          ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch                 ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch      ,& !input  :normalized node number during vegetative growth stages, [-]
    iPlantPhenolPattern_pft           => plt_pheno%iPlantPhenolPattern_pft            ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft               ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    TCChill4Seed_pft                  => plt_pheno%TCChill4Seed_pft                   ,& !input  :temperature below which seed set is adversely affected, [oC]
    NodeLenPergC_pft                  => plt_morph%NodeLenPergC_pft                   ,& !input  :internode length:mass during growth, [m gC-1]
    CanopyLeafArea_pft                => plt_morph%CanopyLeafArea_pft                 ,& !input  :plant canopy leaf area, [m2 d-2]
    MainBranchNum_pft                 => plt_morph%MainBranchNum_pft                  ,& !input  :id number of main branch,[-]
    HoursDoingRemob_brch              => plt_pheno%HoursDoingRemob_brch               ,& !inoput :counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
    PARTS_brch                        => plt_morph%PARTS_brch                          & !output :fraction of C allocated to each morph unit, [-]
  )

  call DebugPrint('beg '//subname//' NZ',NZ)
  PSILY=real((/-200.0_r8,-2.0,-2.0,-2.0/),r8)

  !     begin_execution

  !
  !     PARTITION GROWTH WITHIN EACH BRANCH FROM GROWTH STAGE
  !     1=LEAF,2=SHEATH OR PetolSheth,3=STALK,4=RESERVE,
  !     5,6=REPRODUCTIVE ORGANS,7=GRAIN
  !
  !     PART=organ partitioning fraction
  !

  TOTAL=0._r8
  PART(1:NumOfPlantMorphUnits)=0._r8
  !
  !     IF BEFORE FLORAL INDUCTION
  !
  !plant calendar includes:
  !  planting, emergence, init floral, jointing, elongation, heading, anthesis, begin seed fill, set seed number, set seed mass, end seed fill.
  !
  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).EQ.0)THEN
    PART(ibrch_leaf)   = 0.725_r8
    PART(ibrch_petole) = 0.275_r8
    !
    !     IF BEFORE ANTHESIS
    !
    !  iPlantCalendar_brch(ipltcal_Anthesis,=start of anthesis and setting final seed number
    !  greater maturity group, lower TotalNodeNumNormByMatgrp_brch, more allocation to leaf, later grain-filling
  ELSEIF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).EQ.0)THEN
    PART(ibrch_leaf)   = AMAX1(PART2LEAF_MIN,0.725_r8-FPART1*TotalNodeNumNormByMatgrp_brch(NB,NZ))
    PART(ibrch_petole) = AMAX1(PART2PETOL_MIN,0.275_r8-FPART2*TotalNodeNumNormByMatgrp_brch(NB,NZ))
    PARTS              = 1.0_r8-PART(ibrch_leaf)-PART(ibrch_petole)
    PART(ibrch_stalk)  = 0.60_r8*PARTS
    PART(ibrch_resrv)  = 0.30_r8*PARTS
    PARTX              = PARTS-PART(ibrch_stalk)-PART(ibrch_resrv)
    PART(ibrch_husk)   = 0.5_r8*PARTX
    PART(ibrch_ear)    = 0.5_r8*PARTX
    !
    !     IF BEFORE GRAIN FILLING, DETERMINATE OR INDETERMINATE
    !
    !     iPlantCalendar_brch(ipltcal_BeginSeedFill,=start of grain filling and setting max seed size
    !     iPlantDevelopPattern_pft=growth habit:0=determinate,1=indetermimate from PFT file
    !     TotReproNodeNumNormByMatrgrp_brch=total change in reprv node number normalized for maturity group
    !
  ELSEIF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN
    IF(iPlantDevelopPattern_pft(NZ).EQ.ideterminate)THEN
      PART(ibrch_leaf)   = 0._r8
      PART(ibrch_petole) = 0._r8
    ELSE
      PART(ibrch_leaf)   = AMAX1(PART2LEAF_MIN,(0.725_r8-FPART1)*(1.0_r8-TotReproNodeNumNormByMatrgrp_brch(NB,NZ)))
      PART(ibrch_petole) = AMAX1(PART2PETOL_MIN,(0.275_r8-FPART2)*(1.0_r8-TotReproNodeNumNormByMatrgrp_brch(NB,NZ)))
    ENDIF
    PARTS             = 1.0_r8-PART(ibrch_leaf)-PART(ibrch_petole)
    PART(ibrch_stalk) = 0.60_r8*AZMAX1(PARTS*(1.0_r8-TotReproNodeNumNormByMatrgrp_brch(NB,NZ)))
    PART(ibrch_resrv) = 0.30_r8*AZMAX1(PARTS*(1.0_r8-TotReproNodeNumNormByMatrgrp_brch(NB,NZ)))
    PARTX             = PARTS-PART(ibrch_stalk)-PART(ibrch_resrv)
    PART(ibrch_husk)  = 0.5_r8*PARTX
    PART(ibrch_ear)   = 0.5_r8*PARTX
    !
    !     DURING GRAIN FILLING, DETERMINATE OR INDETERMINATE
    !
    !     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
    !     iPlantDevelopPattern_pft=growth habit:0=determinate,1=indetermimate
    !
  ELSE
    !for determinate growth pattern, stops growing and all put to grain
    IF(iPlantDevelopPattern_pft(NZ).EQ.ideterminate)THEN
      PART(ibrch_grain)=1.0_r8
    ELSE
      PART(ibrch_leaf)   = PART2LEAF_MIN
      PART(ibrch_petole) = PART2PETOL_MIN
      PARTS              = 1.0_r8-PART(ibrch_leaf)-PART(ibrch_petole)

      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
        PART(ibrch_stalk) = 0.125_r8*PARTS
        PART(ibrch_husk)  = 0.125_r8*PARTS
        PART(ibrch_ear)   = 0.125_r8*PARTS
        PART(ibrch_grain) = 0.625_r8*PARTS
      ELSE
        PART(ibrch_stalk) = 0.75_r8*PARTS
        PART(ibrch_grain) = 0.25_r8*PARTS
      ENDIF
    ENDIF
  ENDIF
  !
  !     IF AFTER GRAIN FILLING
  !
  !     iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+PetolSheth,2=none,3=between 1,2
  !
  IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .AND. iPlantCalendar_brch(ipltcal_EndSeedFill,NB,NZ).NE.0)THEN
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
      PART(ibrch_resrv) = 0._r8
      PART(ibrch_stalk) = 0._r8
      PART(ibrch_grain) = 0._r8
    ELSE
      !keep the reserve for perrennial
      PART(ibrch_resrv) = PART(ibrch_resrv)+PART(ibrch_stalk)
      PART(ibrch_stalk) = 0._r8
      PART(ibrch_grain) = 0._r8
    ENDIF
  ENDIF
  !
  !     REDIRECT FROM STALK TO STALK RESERVES IF RESERVES BECOME LOW
  !
  !     WTRSVB,SapwoodBiomassC_brch=stalk reserve,sapwood mass
  !     XFRX=maximum storage C content for remobiln from stalk,root reserves
  !
  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)THEN
    IF(StalkRsrvElms_brch(ielmc,NB,NZ).LT.XFRX*SapwoodBiomassC_brch(NB,NZ))THEN
      D1020: DO N=1,NumOfPlantMorphUnits
        IF(N.NE.ibrch_resrv)THEN
          PART(ibrch_resrv) = PART(ibrch_resrv)+0.10_r8*PART(N)
          PART(N)           = (1._r8-0.10_r8)*PART(N)
        ENDIF
      ENDDO D1020
      !
      !     REDIRECT FROM STALK RESERVES TO STALK IF RESERVES BECOME TOO LARGE
      !
    ELSEIF(StalkRsrvElms_brch(ielmc,NB,NZ).GT.1.0_r8*SapwoodBiomassC_brch(NB,NZ))THEN
      PART(ibrch_stalk) = PART(ibrch_stalk)+PART(ibrch_resrv)+PART(ibrch_grain)
      PART(ibrch_resrv) = 0._r8
      PART(ibrch_grain) = 0._r8
    ENDIF
  ENDIF
  !
  !     REDIRECT FROM LEAVES TO STALK IF LAI BECOMES TOO LARGE
  !
  !     CanopyLeafArea_pft=PFT leaf area
  !
  ARLFI=CanopyLeafArea_pft(NZ)/AREA3(NU)
  IF(ARLFI.GT.5.0_r8)THEN
    FPARTL             = AZMAX1((10.0_r8-ARLFI)/5.0_r8)
    PART(ibrch_stalk)  = PART(ibrch_stalk)+(1.0_r8-FPARTL)*(PART(ibrch_leaf)+PART(ibrch_petole))
    PART(ibrch_leaf)   = FPARTL*PART(ibrch_leaf)
    PART(ibrch_petole) = FPARTL*PART(ibrch_petole)
  ENDIF

  !
  !     DECIDUOUS LEAF FALL AFTER GRAIN FILL IN DETERMINATES,
  !     AFTER AUTUMNIZATION IN INDETERMINATES, OR AFTER SUSTAINED
  !     WATER STRESS
  !
  !     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
  !     FracHour4LeafoffRemob=fraction of hours required for leafoff to initiate remobilization
  !     iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date for setting final seed number
  !     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
  !     LRemob_brch,BegRemoblize=remobilization flags
  !     FLGZ=control rate of remobilization
  !
  ! IFLGZ and IFLGY in grosub.f: phenological eligibility and hourly activation.
  ! Define both outputs even when an eligible deciduous branch has no trigger.
  BegRemoblize = ifalse
  LRemob_brch  = ifalse
  check_annual    = (iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).NE.0)
  check_perennial = (iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .AND. &
    Hours4LeafOff_brch(NB,NZ).GE.FracHour4LeafoffRemob(iPlantPhenolType_pft(NZ))*HourReq4LeafOff_brch(NB,NZ))
  check_cold_deciduous = iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid .OR. &
    iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid
  check_drought_deciduous = iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu .OR. &
    iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid

  IF(check_perennial .OR. check_annual)THEN 
    !set remobilization true
    BegRemoblize = itrue
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen)THEN
      !annual plant or evergreen perennial
      LRemob_brch                 = itrue
      HoursDoingRemob_brch(NB,NZ) = HoursDoingRemob_brch(NB,NZ)+1.0_r8
    ELSEIF(check_cold_deciduous .AND. TdegCCanopy_pft(NZ).LT.TCChill4Seed_pft(NZ))THEN           !chill temperature in effect
      LRemob_brch                 = itrue
      HoursDoingRemob_brch(NB,NZ) = HoursDoingRemob_brch(NB,NZ)+1.0_r8
    ELSEIF(check_drought_deciduous .AND. PSICanopy_pft(NZ).LT.PSILY(iPlantRootProfile_pft(NZ)))THEN !drought stress in effect
      LRemob_brch                 = itrue
      HoursDoingRemob_brch(NB,NZ) = HoursDoingRemob_brch(NB,NZ)+1.0_r8
    ENDIF

    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .AND. iPlantPhenolType_pft(NZ).NE.iphenotyp_evgreen)THEN
      PART(ibrch_stalk)  = PART(ibrch_stalk)+0.5_r8*(PART(ibrch_leaf)+PART(ibrch_petole))
      PART(ibrch_resrv)  = PART(ibrch_resrv)+0.5_r8*(PART(ibrch_leaf)+PART(ibrch_petole))
      PART(ibrch_leaf)   = 0._r8
      PART(ibrch_petole) = 0._r8
    ENDIF
  ELSE
    HoursDoingRemob_brch(NB,NZ) = 0._r8
  ENDIF
  !
  !     CHECK PARTITIONING COEFFICIENTS
  !
  D1000: DO N=1,NumOfPlantMorphUnits
    IF(N.EQ.ibrch_stalk.AND.isclose(NodeLenPergC_pft(NZ),0._r8))THEN
      PART(N)=0._r8
    ELSE
      PART(N)=AZMAX1(PART(N))
    ENDIF
    TOTAL=TOTAL+PART(N)
  ENDDO D1000
  
  IF(TOTAL.GT.ZERO)THEN
    D1010: DO N=1,NumOfPlantMorphUnits
      PART(N)=PART(N)/TOTAL
    ENDDO D1010
  ELSE
    D1015: DO N=1,NumOfPlantMorphUnits
      PART(N)=0._r8
    ENDDO D1015
  ENDIF

  PARTS_brch(:,NB,NZ)=PART  
  call PrintInfo('end '//subname)
  end associate
  end subroutine CalcPartitionCoeff

!----------------------------------------------------------------------------------------------------
  subroutine UpdatePhotosynthates(I,J,NB,NZ,TFN6_vr,YCO2Gro_brch,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,&
    ShootStructN_brch,TFN5,WaterStress4Groth,Stomata_Stress,TurgEff4CanopyResp,CH2O3,CH2O4,CNPG,&
    RMxess_brch,RNonstC4Groth_brch)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: YCO2Gro_brch
  real(r8), intent(in) :: CNLFM,CPLFM
  real(r8), intent(in) :: CNSHX,CPSHX
  real(r8), intent(in) :: CNLFX,CPLFX
  real(r8), intent(in) :: ShootStructN_brch
  real(r8), intent(in) :: TFN5,WaterStress4Groth
  real(r8), intent(in) :: Stomata_Stress
  real(r8), intent(in) :: TurgEff4CanopyResp
  real(r8), intent(out) :: CH2O3(pltpar%MaxNodesPerBranch1)
  real(r8), intent(out) :: CH2O4(pltpar%MaxNodesPerBranch1)
  REAL(R8), INTENT(OUT) :: CNPG
  real(r8), intent(out) :: RMxess_brch                 !excess maintenance respiration, drives remobilization & senescence
  real(r8), intent(out) :: RNonstC4Groth_brch          !nonstructural C for growth, [gC d-2 h-1]

  character(len=*), parameter :: subname='UpdatePhotosynthates'
  real(r8) :: CH2O,CO2F        !CO2 fixation due to photosynthesis
  real(r8) :: CanopyNonstElm4Gros(NumPlantChemElms)      !nonstrucal C/N/P used for growth and respiraiton
  real(r8) :: CH2OClmt,CH2OLlmt
  integer :: K
  associate(                                                 &
    iPlantCalendar_brch  => plt_pheno%iPlantCalendar_brch   ,& !input  :plant growth stage, [-]
    CO2FixCL_pft         => plt_rbgc%CO2FixCL_pft           ,& !inoput :Rubisco-limited CO2 fixation, [gC d-2 h-1]
    CO2FixLL_pft         => plt_rbgc%CO2FixLL_pft           ,& !inoput :Light-limited CO2 fixation, [gC d-h2 h-1]
    GPP_brch             => plt_rbgc%GPP_brch               ,& !inoput :dGPP (C4-C3 product) over branch, [gC d-2 h-1]
    CanopyNonstElms_brch => plt_biom%CanopyNonstElms_brch    & !inoput :branch nonstructural element, [g d-2]
  )
! begin_execution

  call DebugPrint('beg '//subname//' NZ',NZ)
! FDBK=N,P feedback inhibition on C3 CO2 fixation
!

  IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).NE.0)THEN
    !  
    call ComputeGPP(I,J,NB,NZ,WaterStress4Groth,Stomata_Stress,CH2O3,CH2O4,CH2O,CO2F,CH2OClmt,CH2OLlmt)
    
    GPP_brch(NB,NZ)  = GPP_brch(NB,NZ)+CO2F
!   SHOOT AUTOTROPHIC RESPIRATION AFTER EMERGENCE
!   
    call ComputRAutoAfEmergence(I,J,NB,NZ,YCO2Gro_brch,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,&
      CO2F,CH2O,TFN5,WaterStress4Groth,TurgEff4CanopyResp,ShootStructN_brch,CanopyNonstElm4Gros,CNPG,&
      RMxess_brch,RNonstC4Groth_brch)

    CO2FixCL_pft(NZ) = CO2FixCL_pft(NZ)+CH2OClmt   !carbon-dependent photosynthesis
    CO2FixLL_pft(NZ) = CO2FixLL_pft(NZ)+CH2OLlmt   !light-dependent photosynthesis
!   SHOOT AUTOTROPHIC RESPIRATION BEFORE EMERGENCE
!
  ELSE
    CH2O  = 0._r8
    CH2O3 = 0._r8
    CH2O4 = 0._r8    
    call ComputRAutoB4Emergence(I,J,NB,NZ,TFN6_vr,YCO2Gro_brch,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,ShootStructN_brch,&
      WaterStress4Groth,TurgEff4CanopyResp,CanopyNonstElm4Gros,CNPG,RMxess_brch,RNonstC4Groth_brch)
  ENDIF

  !   REMOVE C,N,P USED IN MAINTENANCE + GROWTH REPIRATION AND GROWTH
  !   FROM NON-STRUCTURAL POOLS
  !
  !   CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
  !   CH2O=total CH2O production
  !   RCO2Maint_brch=maintenance respiration
  !   RCO2NonstC_brch=respiration from non-structural C
  !   RCO2NonstC4Nassim_brch=respiration for N assimilation
  !
  CanopyNonstElms_brch(ielmc,NB,NZ) = CanopyNonstElms_brch(ielmc,NB,NZ)-CanopyNonstElm4Gros(ielmc)+CH2O
  CanopyNonstElms_brch(ielmn,NB,NZ) = CanopyNonstElms_brch(ielmn,NB,NZ)-CanopyNonstElm4Gros(ielmn)
  CanopyNonstElms_brch(ielmp,NB,NZ) = CanopyNonstElms_brch(ielmp,NB,NZ)-CanopyNonstElm4Gros(ielmp)

  call PrintInfo('end '//subname)
  end associate
  end subroutine UpdatePhotosynthates

!----------------------------------------------------------------------------------------------------
  subroutine C4PhotoProductTransfer(I,J,NZ,NB,CH2O3,CH2O4)
  implicit none
  integer, intent(in) :: I,J,NZ,NB
  real(r8), intent(in) :: CH2O3(pltpar%MaxNodesPerBranch1)
  real(r8), intent(in) :: CH2O4(pltpar%MaxNodesPerBranch1)
  character(len=*), parameter :: subname='C4PhotoProductTransfer'
  integer :: K
  real(r8) :: CPL4M,CCBS,CPL3K,BundleSheatCO2Leak

!     begin_execution
  associate(                                                             &
    LeafElmntNode_brch         => plt_biom%LeafElmntNode_brch           ,& !input  :leaf element, [g d-2]
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft               ,& !input  :threshold zero for plang growth calculation, [-]
    aquCO2Intraleaf_pft        => plt_photo%aquCO2Intraleaf_pft         ,& !input  :leaf aqueous CO2 concentration, [uM]
    CO2NetFix_pft              => plt_bgcr%CO2NetFix_pft                ,& !inoput :canopy net CO2 exchange, [gC d-2 h-1]
    CanopyResp_brch            => plt_bgcr%CanopyResp_brch              ,& !inoput :canopy respiration for a branch, [gC d-2 h-1]    
    CanopyGrosRCO2_pft         => plt_bgcr%CanopyGrosRCO2_pft           ,& !inoput :canopy plant+nodule autotrophic respiraiton, [gC d-2]
    CanopyRespC_CumYr_pft      => plt_bgcr%CanopyRespC_CumYr_pft        ,& !inoput :total autotrophic respiration, [gC d-2 ]
    Eco_AutoR_CumYr_col        => plt_bgcr%Eco_AutoR_CumYr_col          ,& !inoput :ecosystem autotrophic respiration, [g d-2 h-1]
    ECO_ER_col                 => plt_bgcr%ECO_ER_col                   ,& !inoput :ecosystem respiration, [g d-2 h-1]
    CMassHCO3BundleSheath_node => plt_photo%CMassHCO3BundleSheath_node  ,& !inoput :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CMassCO2BundleSheath_node  => plt_photo%CMassCO2BundleSheath_node   ,& !inoput :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CPOOL3_node                => plt_photo%CPOOL3_node                 ,& !inoput :bundle sheath C4 carbon product to support C3 photosynthesis, [gC d-2]
    CPOOL4_node                => plt_photo%CPOOL4_node                  & !inoput :mesophyll C4 carbon, [g d-2]
  )
  call DebugPrint('beg '//subname//' NZ',NZ)
  D170: DO K=1,MaxNodesPerBranch1
    IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
!
!     MESOPHYLL TO BUNDLE SHEATH TRANSFER
!
!     WGLF=node leaf C mass
!     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
!     CH2O3,CH2O4=total CO2 fixation in bundle sheath,mesophyll
!     CPL4M=mesophyll to bundle sheath transfer of nonstructural C4 that fuels C3 photosynthesis in bundle sheath cells. 
!     FWCBundlSheath,FWCMesophyll=leaf water content in bundle sheath, mesophyll
!
      CPOOL3_node(K,NB,NZ) = CPOOL3_node(K,NB,NZ)-CH2O3(K)
      CPOOL4_node(K,NB,NZ) = CPOOL4_node(K,NB,NZ)+CH2O4(K)
      CPL4M                = 1.0_r8*(CPOOL4_node(K,NB,NZ)*FWCBundlSheath-CPOOL3_node(K,NB,NZ)*FWCMesophyll)/(FWCBundlSheath+FWCMesophyll)
      CPOOL4_node(K,NB,NZ)=CPOOL4_node(K,NB,NZ)-CPL4M
      CPOOL3_node(K,NB,NZ)=CPOOL3_node(K,NB,NZ)+CPL4M
      !
      !     BUNDLE SHEATH CO2 DECARBOXYLATION
      !
      !     CCBS=CO2 concn in bundle sheath (uM)
      !     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
      !     WGLF=node leaf C mass
      !     FWCBundlSheath,FWCMesophyll=leaf water content in bundle sheath, mesophyll
      !     CPL3K=bundle sheath CO2 decarboxylation
      !     CO2KI=Ki for C3 leakage from bundle sheath to mesophyll in C4 (uM)
      !     FCMassCO2BundleSheath_node,FCMassHCO3BundleSheath_node=partition decarboxylation to CO2,HCO3
      !     CPOOL3_node=C4 nonstructural C mass in bundle sheath
      !
      CCBS                 = AZMAX1(0.083E+09_r8*CMassCO2BundleSheath_node(K,NB,NZ)/(LeafElmntNode_brch(ielmc,K,NB,NZ)*FWCBundlSheath))
      CPL3K                = 2.5E-02_r8*CPOOL3_node(K,NB,NZ)/(1.0_r8+CCBS/CO2KI)
      CPOOL3_node(K,NB,NZ) = CPOOL3_node(K,NB,NZ)-CPL3K
      CMassCO2BundleSheath_node(K,NB,NZ)  = CMassCO2BundleSheath_node(K,NB,NZ)+FCMassCO2BundleSheath_node*CPL3K
      CMassHCO3BundleSheath_node(K,NB,NZ) = CMassHCO3BundleSheath_node(K,NB,NZ)+FCMassHCO3BundleSheath_node*CPL3K
      !
      !     BUNDLE SHEATH LEAKAGE
      !
      !     BundleSheatCO2Leak=bundle sheath CO2 leakage, >0, leak
      !     CCBS=CO2 concn in bundle sheath (uM)
      !     aquCO2Intraleaf_pft=intercellular CO2 concentration (uM)
      !     WGLF=node leaf C mass
      !     FWCBundlSheath=leaf water content in bundle sheath
      !     FCMassCO2BundleSheath_node,FCMassHCO3BundleSheath_node=partition decarboxylation to CO2,HCO3
      !     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
      !
      BundleSheatCO2Leak                  = 5.0E-07_r8*(CCBS-aquCO2Intraleaf_pft(NZ))*LeafElmntNode_brch(ielmc,K,NB,NZ)*FWCBundlSheath
      CMassCO2BundleSheath_node(K,NB,NZ)  = CMassCO2BundleSheath_node(K,NB,NZ)-FCMassCO2BundleSheath_node*BundleSheatCO2Leak
      CMassHCO3BundleSheath_node(K,NB,NZ) = CMassHCO3BundleSheath_node(K,NB,NZ)-FCMassHCO3BundleSheath_node*BundleSheatCO2Leak

      !
      !     TOTAL C EXCHANGE
      !
      !     CanopyGrosRCO2_pft,CanopyRespC_CumYr_pft=total,above-ground PFT respiration
      !     CO2NetFix_pft=PFT net CO2 fixation
      !     ECO_ER_col=ecosystem respiration
      !     Eco_AutoR_CumYr_col=total autotrophic respiration
      !     BundleSheatCO2Leak=bundle sheath CO2 leakage
      !
      CanopyResp_brch(NB,NZ)    = CanopyResp_brch(NB,NZ)-BundleSheatCO2Leak
      CanopyGrosRCO2_pft(NZ)    = CanopyGrosRCO2_pft(NZ)-BundleSheatCO2Leak
      CanopyRespC_CumYr_pft(NZ) = CanopyRespC_CumYr_pft(NZ)-BundleSheatCO2Leak
      CO2NetFix_pft(NZ)         = CO2NetFix_pft(NZ)-BundleSheatCO2Leak
      ECO_ER_col                = ECO_ER_col-BundleSheatCO2Leak
      Eco_AutoR_CumYr_col       = Eco_AutoR_CumYr_col-BundleSheatCO2Leak
    ENDIF
  ENDDO D170
  call PrintInfo('end '//subname)
  end associate
  end subroutine C4PhotoProductTransfer

!----------------------------------------------------------------------------------------------------
  subroutine RemobizeLeafNodes(I,J,NB,NZ,KN,RCCE,FracRecyled,RespSenes_node,lgoto565,RCO2e)
  implicit none
  integer,intent(in) :: I,J,NB,NZ,KN
  real(r8), intent(in) :: RCCE(NumPlantChemElms)
  real(r8), intent(in):: FracRecyled            !recyled fraction
  real(r8), intent(inout) :: RespSenes_node     !potential senescence C flux per leaf node, [gC h-1]  
  logical,intent(out):: lgoto565
  real(r8), intent(inout) :: RCO2e

  character(len=*), parameter :: subname='RemobizeLeafNodes'
  real(r8) :: RLeaf2Senes   !leaf senescence flux [gC h-1]
  real(r8) :: FSNAS,FracNodeAsLeaf
  real(r8) :: RSpecKillPetol  !specific petole turnover rate due to remobilization [h-1]
  real(r8) :: SensResp2Use,Frac2fall,dfall
  integer  :: K,KK,NE,M
  real(r8) :: RCES(NumPlantChemElms),dleaf,dPetole,dmassE  
  real(r8) :: RemobledLeafElm(NumPlantChemElms)

  !begin_execution
  associate(                                                               &
    KHiestGroLeafNode_brch      => plt_pheno%KHiestGroLeafNode_brch       ,& !input  :leaf growth stage counter, [-]
    inonfoliar                  => pltpar%inonfoliar                      ,& !input  :group id of plant non-foliar litter group
    icarbhyro                   => pltpar%icarbhyro                       ,& !input  :kinetic id of litter component as carbonhydrate
    icwood                      => pltpar%icwood                          ,& !input  :group id of coarse woody litter
    ifoliar                     => pltpar%ifoliar                         ,& !input  :group id of plant foliar litter
    FracLeafShethElmAlloc2Litr  => plt_allom%FracLeafShethElmAlloc2Litr   ,& !input  :woody element allocation, [-]
    FracPetolShethAlloc2Litr    => plt_allom%FracPetolShethAlloc2Litr     ,& !input  :leaf element allocation,[-]
    PlantElmAllocMat4Litr       => plt_soilchem%PlantElmAllocMat4Litr     ,& !input  :litter kinetic fraction, [-]
    rProteinC2LeafN_pft         => plt_allom%rProteinC2LeafN_pft          ,& !input  :Plant leaf protein C to N ratio [g protein C (g leaf N)-1]
    rProteinC2LeafP_pft         => plt_allom%rProteinC2LeafP_pft          ,& !input  :Plant leaf protein C to P ratio [g protein C (g leaf P)-1]
    k_fine_comp                 => pltpar%k_fine_comp                     ,& !input  :fine litter complex id
    k_woody_comp                => pltpar%k_woody_comp                    ,& !input  :woody litter complex id
    ZERO4LeafVar_pft            => plt_biom%ZERO4LeafVar_pft              ,& !input  :threshold zero for leaf calculation, [-]
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft                ,& !input  :threshold zero for plang growth calculation, [-]
    PetolShethElmntNode_brch    => plt_biom%PetolShethElmntNode_brch      ,& !inoput :Petiole/sheath chemical element, [g d-2]
    LitrfallElms_pvr            => plt_bgcr%LitrfallElms_pvr              ,& !inoput :plant LitrFall element, [g d-2 h-1]
    LitrFallElms_brch           => plt_bgcr%LitrFallElms_brch             ,& !inoput :litterfall from the branch, [g d-2 h-1]    
    CPOOL3_node                 => plt_photo%CPOOL3_node                  ,& !inoput :minimum sink strength for nonstructural C transfer, [g d-2]
    CPOOL4_node                 => plt_photo%CPOOL4_node                  ,& !inoput :leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
    PetoleProteinC_node         => plt_biom%PetoleProteinC_node           ,& !inoput :layer sheath protein C, [g d-2]
    LeafProteinC_node           => plt_biom%LeafProteinC_node             ,& !inoput :layer leaf protein C, [g d-2]
    LeafAreaLive_brch           => plt_morph%LeafAreaLive_brch            ,& !inoput :branch leaf area, [m2 d-2]
    LeafArea_node               => plt_morph%LeafArea_node                ,& !inoput :leaf area, [m2 d-2]
    LeafStrutElms_brch          => plt_biom%LeafStrutElms_brch            ,& !inoput :branch leaf structural element mass, [g d-2]
    LeafElmntNode_brch          => plt_biom%LeafElmntNode_brch            ,& !inoput :leaf element, [g d-2]
    PetolShethStrutElms_brch    => plt_biom%PetolShethStrutElms_brch      ,& !inoput :branch sheath structural element, [g d-2]
    PetoleLength_node           => plt_morph%PetoleLength_node            ,& !inoput :sheath height, [m]
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch           & !inoput :branch nonstructural element, [g d-2]
  )

  call DebugPrint('beg '//subname//' NZ',NZ)
  lgoto565=.false.
  D650: DO KK=KN,KHiestGroLeafNode_brch(NB,NZ)
    RLeaf2Senes  = 0._r8
    SensResp2Use = 0._r8

    K  = pMOD(KK,MaxNodesPerBranch1)
!
!       REMOBILIZATION OF LEAF C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
!
!       WGLF,PetolShethElmntNode_brch=node leaf,PetolSheth C mass
  !       SCNF,SCNSH=leaf,PetolSheth senescence respiration
  !       RemobledLeafElm(:)=remobilization of C,N,P from senescing leaf
  !       RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
  !       RCCZ,RCCY=min,max fractions for shoot C recycling
  !
    IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FracNodeAsLeaf = LeafElmntNode_brch(ielmc,K,NB,NZ)/(LeafElmntNode_brch(ielmc,K,NB,NZ)+PetolShethElmntNode_brch(ielmc,K,NB,NZ))
      RLeaf2Senes    = FracNodeAsLeaf*RespSenes_node
      SensResp2Use   = RespSenes_node-RLeaf2Senes
      DO NE=1,NumPlantChemElms
        RemobledLeafElm(NE) = LeafElmntNode_brch(NE,K,NB,NZ)*RCCE(NE)
      ENDDO
      !
      !         FRACTION OF CURRENT LEAF TO BE REMOBILIZED
      !
      !         Frac2fall=fraction of current leaf C,area to be remobilized
      !senescence flux = recycled + litter fall

      IF(RemobledLeafElm(ielmc).GT.ZERO4Groth_pft(NZ))THEN
        Frac2fall=AZMAX1(AMIN1(1.0_r8,RLeaf2Senes/RemobledLeafElm(ielmc)))
      ELSE
        Frac2fall=1.0_r8
      ENDIF
      !
      !     NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
      !     TO FRACTIONS SET IN 'STARTQ'
      !
      !     CSNC,ZSNC,PSNC=literfall C,N,P
      !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
      !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
      !     Frac2fall=fraction of current leaf to be remobilized
      !     WGLF,WGLFN,WGLFP=node leaf C,N,P mass
      !     RemobledLeafElm(:)=remobilization of C,N,P from senescing leaf
      !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
      !     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
      !     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
      ! the senescing leaf = remobilized + falling litter
      !
      D6310: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          dleaf = AZMAX1(LeafElmntNode_brch(NE,K,NB,NZ)-RemobledLeafElm(NE))*Frac2fall
          dFall = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dleaf*FracLeafShethElmAlloc2Litr(NE,k_woody_comp)
          LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)              = LitrFallElms_brch(NE,NB,NZ)+dFall

          dFall=PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)*dleaf*FracLeafShethElmAlloc2Litr(NE,k_fine_comp)
          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) +dFall
          LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ)+dFall
        ENDDO
      ENDDO D6310

      IF(K.NE.0)THEN
        dFall=Frac2fall*AZMAX1(CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ))
        LitrfallElms_pvr(ielmc,icarbhyro,k_fine_comp,0,NZ) = LitrfallElms_pvr(ielmc,icarbhyro,k_fine_comp,0,NZ)+dFall
        LitrFallElms_brch(ielmc,NB,NZ)                     = LitrFallElms_brch(ielmc,NB,NZ)+dFall

        CPOOL3_node(K,NB,NZ)=CPOOL3_node(K,NB,NZ)*AZMAX1(1._r8-Frac2fall)
        CPOOL4_node(K,NB,NZ)=CPOOL4_node(K,NB,NZ)*AZMAX1(1._r8-Frac2fall)
      ENDIF
      !
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
      !
      !     LeafAreaLive_brch=leaf area
      !     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
      !     WGLF,WGLFN,WGLFP,LeafProteinC_node=node leaf C,N,P,protein mass
      !     Frac2fall=fraction of current leaf to be remobilized
      !     CNWS,rProteinC2LeafP_pft=protein:N,protein:P ratios from startq.f
      !
      LeafAreaLive_brch(NB,NZ) = AZMAX1(LeafAreaLive_brch(NB,NZ)-Frac2fall*LeafArea_node(K,NB,NZ))
      LeafArea_node(K,NB,NZ)   = AZMAX1(LeafArea_node(K,NB,NZ)-Frac2fall*LeafArea_node(K,NB,NZ))

      LeafProteinC_node(K,NB,NZ)=AZMAX1(LeafProteinC_node(K,NB,NZ)-Frac2fall* &
        AMAX1(LeafElmntNode_brch(ielmn,K,NB,NZ)*rProteinC2LeafN_pft(NZ),&
          LeafElmntNode_brch(ielmp,K,NB,NZ)*rProteinC2LeafP_pft(NZ)))
      DO NE=1,NumPlantChemElms
        LeafStrutElms_brch(NE,NB,NZ)   = AZMAX1(LeafStrutElms_brch(NE,NB,NZ)-Frac2fall*LeafElmntNode_brch(NE,K,NB,NZ))
        LeafElmntNode_brch(NE,K,NB,NZ) = AZMAX1(LeafElmntNode_brch(NE,K,NB,NZ)-Frac2fall*LeafElmntNode_brch(NE,K,NB,NZ))
      ENDDO
      !
      !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
      !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
      !
      !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
      !     Frac2fall=fraction of current leaf C to be remobilized
      !     RemobledLeafElm(:)=remobilization of C,N,P from senescing leaf
      !     RLeaf2Senes,RespSenes_node=remaining senescence respiration carried to next node
      !
      RCO2e=RCO2e+Frac2fall*RemobledLeafElm(ielmc)*(1._r8-FracRecyled)
      CanopyNonstElms_brch(ielmc,NB,NZ)=CanopyNonstElms_brch(ielmc,NB,NZ)+Frac2fall*RemobledLeafElm(ielmc)*FracRecyled

      DO NE=2,NumPlantChemElms
        CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+Frac2fall*RemobledLeafElm(NE)
      ENDDO
      RLeaf2Senes    = RLeaf2Senes-Frac2fall*RemobledLeafElm(ielmc)
      RespSenes_node = RespSenes_node-Frac2fall*RemobledLeafElm(ielmc)
      IF(LeafStrutElms_brch(ielmc,NB,NZ).LE.ZERO4LeafVar_pft(NZ))THEN
        LeafStrutElms_brch(ielmc,NB,NZ) = 0._r8
        LeafAreaLive_brch(NB,NZ)        = 0._r8
      ENDIF
      !
      !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
      !
!        IF(RLeaf2Senes.LE.ZERO4Groth_pft(NZ))GO TO 564
      !
      !     OTHERWISE REMAINING C,N,P IN LEAF GOES TO LitrFall
      !
      !     CSNC,ZSNC,PSNC=literfall C,N,P
      !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
      !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
      !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
      !     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
      !     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
      !     LeafAreaLive_brch=leaf area
      !     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
      !     WGLF,WGLFN,WGLFP,LeafProteinC_node=node leaf C,N,P,protein mass
      !
    ELSE        
      D6315: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          dFall=PlantElmAllocMat4Litr(NE,icwood,M,NZ)*AZMAX1(LeafElmntNode_brch(NE,K,NB,NZ))*FracLeafShethElmAlloc2Litr(NE,k_woody_comp)
          LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)  =LitrFallElms_brch(NE,NB,NZ)+dFall  

          dFall=PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)*AZMAX1(LeafElmntNode_brch(NE,K,NB,NZ))*FracLeafShethElmAlloc2Litr(NE,k_fine_comp)
          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)  =LitrFallElms_brch(NE,NB,NZ)+dFall  
            
        ENDDO
      ENDDO D6315

      IF(K.NE.0)THEN
        dFall=AZMAX1(CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ))
        LitrfallElms_pvr(ielmc,icarbhyro,k_fine_comp,0,NZ)=LitrfallElms_pvr(ielmc,icarbhyro,k_fine_comp,0,NZ)+dFall
        LitrFallElms_brch(ielmc,NB,NZ)  =LitrFallElms_brch(ielmc,NB,NZ)+dFall  

        CPOOL3_node(K,NB,NZ) = 0._r8
        CPOOL4_node(K,NB,NZ) = 0._r8
      ENDIF

      LeafAreaLive_brch(NB,NZ)=AZMAX1(LeafAreaLive_brch(NB,NZ)-LeafArea_node(K,NB,NZ))
      DO NE=1,NumPlantChemElms
        LeafStrutElms_brch(NE,NB,NZ)   = AZMAX1(LeafStrutElms_brch(NE,NB,NZ)-LeafElmntNode_brch(NE,K,NB,NZ))
        LeafElmntNode_brch(NE,K,NB,NZ) = 0._r8
      ENDDO
      LeafArea_node(K,NB,NZ)     = 0._r8
      LeafProteinC_node(K,NB,NZ) = 0._r8
      IF(LeafStrutElms_brch(ielmc,NB,NZ).LE.ZERO4LeafVar_pft(NZ))THEN
        LeafStrutElms_brch(ielmc,NB,NZ) = 0._r8
        LeafAreaLive_brch(NB,NZ)        = 0._r8
      ENDIF
    ENDIF
    !
    !     REMOBILIZATION OF SHEATHS OR PetolSheth C,N,P DEPENDS ON
    !     NON-STRUCTURAL C:N:P
    !
    !     PetolShethElmntNode_brch,WGSHN,WGSHP=node PetolSheth C,N,P mass
    !     RCES(ielmc),RCES(ielmn),RCES(ielmp)=remobilization of C,N,P from senescing PetolSheth
    !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !
    IF(PetolShethElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      DO NE=1,NumPlantChemElms
        RCES(NE)=PetolShethElmntNode_brch(NE,K,NB,NZ)*RCCE(NE)
      ENDDO
      !
      !     FRACTION OF REMOBILIZATION THAT CAN BE MET FROM CURRENT SHEATH
      !     OR PetolSheth
      !
      IF(RCES(ielmc).GT.ZERO4Groth_pft(NZ))THEN
        RSpecKillPetol=AZMAX1(AMIN1(1.0_r8,SensResp2Use/RCES(ielmc)))
      ELSE
        RSpecKillPetol=1.0_r8
      ENDIF
      
      !
      !     NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
      !     TO FRACTIONS SET IN 'STARTQ'
      !
      !     CSNC,ZSNC,PSNC=literfall C,N,P
      !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
      !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
      !     RSpecKillPetol=fraction of current PetolSheth to be remobilized
      !     PetolShethElmntNode_brch,WGSHN,WGSHP=node PetolSheth C,N,P mass
      !     RCES(ielmc),RCES(ielmn),RCES(ielmp)=remobilization of C,N,P from senescing PetolSheth
      !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
      !     FWODSN,FWODSP=N,P woody fraction in PetolSheth:0=woody,1=non-woody
      !
      
      D6320: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          dPetole                                  = RSpecKillPetol*AZMAX1(PetolShethElmntNode_brch(NE,K,NB,NZ)-RCES(NE))
          dFall                                    = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dPetole*FracPetolShethAlloc2Litr(NE,k_woody_comp)
          LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)              = LitrFallElms_brch(NE,NB,NZ)+dFall
            
          dFall                                   = PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*dPetole*FracPetolShethAlloc2Litr(NE,k_fine_comp)
          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ)+dFall

        ENDDO    
      ENDDO D6320
      PetoleProteinC_node(K,NB,NZ) = AZMAX1(PetoleProteinC_node(K,NB,NZ)-RSpecKillPetol*&
        AMAX1(PetolShethElmntNode_brch(ielmn,K,NB,NZ)*rProteinC2LeafN_pft(NZ),&
          PetolShethElmntNode_brch(ielmp,K,NB,NZ)*rProteinC2LeafP_pft(NZ)))

      DO NE=1,NumPlantChemElms
        dmassE=RSpecKillPetol*PetolShethElmntNode_brch(NE,K,NB,NZ)
        PetolShethStrutElms_brch(NE,NB,NZ)    = AZMAX1(PetolShethStrutElms_brch(NE,NB,NZ)-dmassE)
        PetolShethElmntNode_brch(NE,K,NB,NZ) = AZMAX1(PetolShethElmntNode_brch(NE,K,NB,NZ)-dmassE)
      ENDDO
      !
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
      !
      !     PetoleLength_node=PetolSheth length
      !     WTSHEB,WTLFBN,WTSHBP=branch PetolSheth C,N,P mass
      !     PetolShethElmntNode_brch,WGSHN,WGSHP,PetoleProteinC_node=node PetolSheth C,N,P,protein mass
      !     RSpecKillPetol=fraction of current PetolSheth to be remobilized
      !     CNWS,rProteinC2LeafP_pft=protein:N,protein:P ratios from startq.f
      !
      PetoleLength_node(K,NB,NZ)     = PetoleLength_node(K,NB,NZ)-RSpecKillPetol*PetoleLength_node(K,NB,NZ)
      !
      !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
      !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
      !
      !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
      !     RSpecKillPetol=fraction of current PetolSheth C to be remobilized
      !     RCES(ielmc),RCES(ielmn),RCES(ielmp)=remobilization of C,N,P from senescing PetolSheth
      !     SensResp2Use,RespSenes_node=remaining senescence respiration carried to next node
      !
      RCO2e=RCO2e+RSpecKillPetol*RCES(ielmc)*(1._r8-FracRecyled)
      CanopyNonstElms_brch(ielmc,NB,NZ)=CanopyNonstElms_brch(ielmc,NB,NZ)+RSpecKillPetol*RCES(ielmc)*FracRecyled
      DO NE=2,NumPlantChemElms
        CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+RSpecKillPetol*RCES(NE)
      ENDDO

      SensResp2Use   = SensResp2Use-RSpecKillPetol*RCES(ielmc)
      RespSenes_node = RespSenes_node-RSpecKillPetol*RCES(ielmc)
      IF(PetolShethStrutElms_brch(ielmc,NB,NZ).LE.ZERO4LeafVar_pft(NZ))THEN
        PetolShethStrutElms_brch(ielmc,NB,NZ)=0._r8
      ENDIF
!
!     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
!
      lgoto565=(SensResp2Use.LE.ZERO4Groth_pft(NZ))
!
    !     OTHERWISE REMAINING C,N,P IN SHEATH OR PetolSheth GOES TO LitrFall
    !
    !     CSNC,ZSNC,PSNC=literfall C,N,P
    !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
    !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
    !     FWODB=C woody fraction in branch:0=woody,1=non-woody
    !     FWODSN,FWODSP=N,P woody fraction in PetolSheth:0=woody,1=non-woody
    !     PetoleLength_node=PetolSheth length
    !     WTSHEB,WTSHBN,WTSHBP=branch PetolSheth C,N,P mass
    !     PetolShethElmntNode_brch,WGSHN,WGSHP,PetoleProteinC_node=node PetolSheth C,N,P,protein mass
    !
    ELSE      
      D6325: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          dFall                                    = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*AZMAX1(PetolShethElmntNode_brch(NE,K,NB,NZ))*FracPetolShethAlloc2Litr(NE,k_woody_comp)
          LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)              = LitrFallElms_brch(NE,NB,NZ)+dFall

          dFall                                   = PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(PetolShethElmntNode_brch(NE,K,NB,NZ))*FracPetolShethAlloc2Litr(NE,k_fine_comp)
          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ)+dFall

        ENDDO    
      ENDDO D6325
      DO NE=1,NumPlantChemElms
        PetolShethStrutElms_brch(NE,NB,NZ)    = AZMAX1(PetolShethStrutElms_brch(NE,NB,NZ)-PetolShethElmntNode_brch(NE,K,NB,NZ))
        PetolShethElmntNode_brch(NE,K,NB,NZ) = 0._r8
      ENDDO
      PetoleLength_node(K,NB,NZ)     = 0._r8
      PetoleProteinC_node(K,NB,NZ) = 0._r8
      IF(PetolShethStrutElms_brch(ielmc,NB,NZ).LE.ZERO4LeafVar_pft(NZ))THEN
        PetolShethStrutElms_brch(ielmc,NB,NZ)=0._r8
      ENDIF
    ENDIF
  ENDDO D650
  call PrintInfo('end '//subname)
  end associate
  END subroutine RemobizeLeafNodes

!----------------------------------------------------------------------------------------------------
  subroutine RemobilizeLeafLayers(I,J,NumRemobLeafNodes,NB,NZ,RespSenesTot_brch,RCCC,RCCN,RCCP,FracRecyled,RCO2e)
  implicit none
  integer, intent(in)  :: I,J
  INTEGER, intent(in)  :: NB,NZ,NumRemobLeafNodes
  real(r8), intent(in) :: RespSenesTot_brch   !total C flux needs to filled through remobilization or litterfall
  real(r8), intent(in) :: RCCC,RCCN,RCCP
  real(r8), intent(inout):: FracRecyled       !fraction of biomass that is recycled into reserve
  real(r8), intent(inout) :: RCO2e            !CO2 flux associated with senescence, [gC h-1]

  character(len=*), parameter :: subname='RemobilizeLeafLayers'
  integer :: N,M,K,KK,MXNOD,MNNOD,NE
  real(r8) :: FRCC
  real(r8) :: FSNCK
  real(r8) :: FSNCR
  real(r8) :: HTNODZ
  real(r8) :: FStalkSenes
  real(r8) :: RCSC,RCSN,RCSP
  real(r8) :: RCEK(NumPlantChemElms)
  REAL(R8)  :: RCCE(NumPlantChemElms)  
  real(r8) :: RMxess_brch
  real(r8) :: RespPhenolSenes_brch
  real(r8) :: RespSenes_node,dIntNode,dStalk,dNodeH,dFall
  integer :: KN  
  logical :: lgoto565
! begin_execution
  associate(                                                             &
    istalk                     => pltpar%istalk                         ,& !input  :group id of plant stalk litter group
    k_fine_comp                => pltpar%k_fine_comp                    ,& !input  :fine litter complex id
    k_woody_comp               => pltpar%k_woody_comp                   ,& !input  :woody litter complex id
    SapwoodBiomassC_brch       => plt_biom%SapwoodBiomassC_brch         ,& !input  :branch sapwood C, [gC d-2]
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft               ,& !input  :threshold zero for plang growth calculation, [-]
    FracWoodStalkElmAlloc2Litr => plt_allom%FracWoodStalkElmAlloc2Litr  ,& !input  :woody element allocation,[-]
    PlantElmAllocMat4Litr      => plt_soilchem%PlantElmAllocMat4Litr    ,& !input  :litter kinetic fraction, [-]
    KLowestGroLeafNode_brch    => plt_pheno%KLowestGroLeafNode_brch     ,& !input  :leaf growth stage counter, [-]
    KHiestGroLeafNode_brch     => plt_pheno%KHiestGroLeafNode_brch      ,& !input  :leaf growth stage counter, [-]
    iPlantPhenolPattern_pft    => plt_pheno%iPlantPhenolPattern_pft     ,& !input  :plant growth habit: annual or perennial,[-]
    NumCogrowthNode_pft        => plt_morph%NumCogrowthNode_pft         ,& !input  :number of concurrently growing nodes,[-]
    SSXferElms_pft             => plt_bgcr%SSXferElms_pft               ,& !inoput :flux export from seasonal storage, [g h-1 d-2]
    SSXfer2ShootElms_pft       => plt_bgcr%SSXfer2ShootElms_pft         ,& !inoput :flux export from seasonal storage to shoot, [g h-1 d-2]        
    StalkRsrvElms_brch         => plt_biom%StalkRsrvElms_brch           ,& !inoput :branch reserve element mass, [g d-2]
    StalkStrutElms_brch        => plt_biom%StalkStrutElms_brch          ,& !inoput :branch stalk structural element mass, [g d-2]
    StructInternodeElms_brch   => plt_biom%StructInternodeElms_brch     ,& !inoput :internode structural elelments, [g d-2]
    SeasonalNonstElms_pft      => plt_biom%SeasonalNonstElms_pft        ,& !inoput :plant stored nonstructural element at current step, [g d-2]
    SenecStalkStrutElms_brch   => plt_biom%SenecStalkStrutElms_brch     ,& !inoput :senescing branch stalk structural element, [g d-2]
    LitrfallElms_pvr           => plt_bgcr%LitrfallElms_pvr             ,& !inoput :plant LitrFall element, [g d-2 h-1]
    LitrFallElms_brch          => plt_bgcr%LitrFallElms_brch            ,& !inoput :litterfall from the branch, [g d-2 h-1]    
    StalkNodeHeight_brch       => plt_morph%StalkNodeHeight_brch        ,& !inoput :live internode height, [m]
    StalkNodeVertLength_brch   => plt_morph%StalkNodeVertLength_brch    ,& !inoput :dead internode height, [m]
    isPlantBranchAlive_brch    => plt_pheno%isPlantBranchAlive_brch      & !output :flag to detect branch death, [-]
  )
  !     REMOBILIZATION AND LitrFall WHEN GROWTH RESPIRATION < 0
  !     STARTING FROM LOWEST LEAFED NODE AND PROCEEDING UPWARDS
  !
  !     RespSenesTot_brch,RespSenes_node=branch,node senescence respiration
  !     KSNC=number of nodes undergoing remobilization
  !
  call DebugPrint('beg '//subname//' NZ',NZ)
  RCCE(ielmc)=RCCC
  RCCE(ielmn)=(RCCN+(1.0_r8-RCCN)*RCCC)
  RCCE(ielmp)=(RCCP+(1.0_r8-RCCP)*RCCC)

  KN=MAX(0,KLowestGroLeafNode_brch(NB,NZ)-1)
  D575: DO N=1,NumRemobLeafNodes
    RespSenes_node=RespSenesTot_brch/real(NumRemobLeafNodes,kind=r8)

    call RemobizeLeafNodes(I,J,NB,NZ,KN,RCCE,FracRecyled,RespSenes_node,lgoto565,RCO2e)
    if(lgoto565)cycle
    
    KN=KN+1
    RMxess_brch=RespSenes_node*(1.0_r8-FracRecyled)
    !
    !     REMOBILIZATION OF RESERVE C
    !
    !     WTRSVB=stalk reserve C mass
    !     RMxess_brch=excess maintenance respiration
    !
    IF(StalkRsrvElms_brch(ielmc,NB,NZ).GT.RMxess_brch)THEN
      StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)-RMxess_brch
      RCO2e       = RCO2e+RMxess_brch
      RMxess_brch = 0._r8
      cycle
    ENDIF
    !
    !     REMOBILIZATION OF STALK C,N,P
    !
    !     FXFS=rate constant for remobilization of stalk C,N,P (h-1)
    !     RespPhenolSenes_brch=phenologically-driven respiration senescence during late-season
    !     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
    !     WTSTKB,SapwoodBiomassC_brch=stalk,sapwood C mass
    !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !     MXNOD,MNNOD=max,min node number currently growing
    !     NumCogrowthNode_pft=number of concurrently growing nodes
    !     KHiestGroLeafNode_brch=integer of most recent leaf number
    !

    RespPhenolSenes_brch = FXFS*RMxess_brch
    RespSenes_node       = RMxess_brch+RespPhenolSenes_brch
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .AND. RespSenes_node.GT.ZERO4Groth_pft(NZ) &
      .AND. StalkStrutElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FracRecyled = RespPhenolSenes_brch/RespSenes_node
      FRCC        = SapwoodBiomassC_brch(NB,NZ)/StalkStrutElms_brch(ielmc,NB,NZ)
      RCSC        = RCCC*FRCC
      RCSN        = RCCN*FRCC
      RCSP        = RCCP*FRCC
      MXNOD       = KHiestGroLeafNode_brch(NB,NZ)
      MNNOD       = MAX(MIN(0,MAX(0,MXNOD-NumCogrowthNode_pft(NZ))),KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+2)
      MXNOD       = MAX(MXNOD,MNNOD)

      D1650: DO KK = MXNOD, MNNOD, -1
        K=pMOD(KK,MaxNodesPerBranch1)
        !
        !     REMOBILIZATION OF STALK C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
        !
        !     RCEK(:)=remobilization of C,N,P from senescing internode
        !     StructInternodeElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
        !
        IF(StructInternodeElms_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
          RCEK(ielmc)=StructInternodeElms_brch(ielmc,K,NB,NZ)*RCSC
          RCEK(ielmn)=StructInternodeElms_brch(ielmn,K,NB,NZ)*(RCSN+(1.0_r8-RCSN)*RCSC)
          RCEK(ielmp)=StructInternodeElms_brch(ielmp,K,NB,NZ)*(RCSP+(1.0_r8-RCSP)*RCSC)
          !
          !     FRACTION OF CURRENT NODE TO BE REMOBILIZED
          !
          !
          IF(RCEK(ielmc).GT.ZERO4Groth_pft(NZ))THEN
            FSNCK=AZMAX1(AMIN1(1.0,RespSenes_node/RCEK(ielmc)))
          ELSE
            FSNCK=1.0_r8
          ENDIF
          !
          !     NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
          !     TO FRACTIONS SET IN 'STARTQ'
          !
          !     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
          !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
          !     FSNCK=fraction of lowest internode to be remobilized
          !     RCEK(:)=remobilization of C,N,P from senescing internode
          !     StructInternodeElms_brch,WGNODN,WGNODP=senescing internode C,N,P mass
          !
          
          D7310: DO M=1,jsken
            DO NE=1,NumPlantChemElms
              dIntNode                                 = AZMAX1(StructInternodeElms_brch(NE,K,NB,NZ)-RCEK(NE))*FSNCK
              dFall                                    = PlantElmAllocMat4Litr(NE,istalk,M,NZ)*dIntNode*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp)
              LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)+dFall
              LitrFallElms_brch(NE,NB,NZ)              = LitrFallElms_brch(NE,NB,NZ) +dFall

              dFall                                   = PlantElmAllocMat4Litr(NE,istalk,M,NZ)*dIntNode*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp)
              LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+dFall
              LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ) +dFall

            ENDDO    
          ENDDO D7310
          !
          !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
          !
          !     FSNCK=fraction of lowest internode to be remobilized
          !     StalkNodeHeight_brch,StalkNodeVertLength_brch=living,senescing internode length
          !     WTSTKB,WTSTBN,WTSTBP,StructInternodeElms_brch,WGNODN,WGNODP=C,N,P mass in senescing internode
          !
          DO NE=1,NumPlantChemElms
            dNodeH                               = FSNCK*StructInternodeElms_brch(NE,K,NB,NZ)
            StalkStrutElms_brch(NE,NB,NZ)        = AZMAX1(StalkStrutElms_brch(NE,NB,NZ)-dNodeH)
            StructInternodeElms_brch(NE,K,NB,NZ) = AZMAX1(StructInternodeElms_brch(NE,K,NB,NZ)-dNodeH)
          ENDDO
          !
          dNodeH                            = FSNCK*StalkNodeVertLength_brch(K,NB,NZ)
          StalkNodeHeight_brch(K,NB,NZ)     = AZMAX1(StalkNodeHeight_brch(K,NB,NZ)-dNodeH)
          StalkNodeVertLength_brch(K,NB,NZ) = AZMAX1(StalkNodeVertLength_brch(K,NB,NZ)-dNodeH)

          !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
          !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
          !
          !     RCEK(:)=remobilization of C,N,P from senescing internode
          !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
          !     FSNCK=fraction of lowest internode to be remobilized
          !     RespSenes_node=remaining node senescence respiration
          !
          RCO2e=RCO2e+FSNCK*RCEK(ielmc)*(1._r8-FracRecyled)
          StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)+FSNCK*RCEK(ielmc)*FracRecyled
          DO NE=2,NumPlantChemElms
            StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+FSNCK*RCEK(NE)
          ENDDO
          
          RespSenes_node=RespSenes_node-FSNCK*RCEK(ielmc)
          !
          !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
          !
          IF(RespSenes_node.LE.ZERO4Groth_pft(NZ))then
            lgoto565=.true.
            exit
          else
            lgoto565=.false.
          endif
          !
          !       OTHERWISE REMAINING C,N,P IN NODE GOES TO LitrFall
          !
          !       CSNC,ZSNC,PSNC=literfall C,N,P
          !       CFOPC=fraction of plant litter allocated in nonstructural(0,*),
          !       WTSTKB,WTSTBN,WTSTBP,StructInternodeElms_brch,WGNODN,WGNODP=C,N,P mass in senescing internode
          !       StalkNodeHeight_brch,StalkNodeVertLength_brch=living,senescing internode length
          !
        ELSE
          
          D7315: DO M=1,jsken
            DO NE=1,NumPlantChemElms
              dFall                                    = PlantElmAllocMat4Litr(NE,istalk,M,NZ)*AZMAX1(StructInternodeElms_brch(NE,K,NB,NZ))*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp)
              LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)+dFall
              LitrFallElms_brch(NE,NB,NZ)              = LitrFallElms_brch(NE,NB,NZ) +dFall
                
              dFall                                   = PlantElmAllocMat4Litr(NE,istalk,M,NZ)*AZMAX1(StructInternodeElms_brch(NE,K,NB,NZ))*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp)
              LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) +dFall
              LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ) +dFall

            ENDDO    
          ENDDO D7315

          DO NE=1,NumPlantChemElms
            StalkStrutElms_brch(NE,NB,NZ)        = AZMAX1(StalkStrutElms_brch(NE,NB,NZ)-StructInternodeElms_brch(NE,K,NB,NZ))
            StructInternodeElms_brch(NE,K,NB,NZ) = 0._r8
          ENDDO
          StalkNodeHeight_brch(K,NB,NZ)     = StalkNodeHeight_brch(K,NB,NZ)-StalkNodeVertLength_brch(K,NB,NZ)
          StalkNodeVertLength_brch(K,NB,NZ) = 0._r8

        ENDIF
      ENDDO D1650

      if(lgoto565)cycle
      !
      !   RESIDUAL STALK
      !
      !   RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
      !   WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
      !
      IF(SenecStalkStrutElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
        RCEK(ielmc)=SenecStalkStrutElms_brch(ielmc,NB,NZ)*RCSC
        RCEK(ielmn)=SenecStalkStrutElms_brch(ielmn,NB,NZ)*(RCSN+(1.0_r8-RCSN)*RCSC)
        RCEK(ielmp)=SenecStalkStrutElms_brch(ielmp,NB,NZ)*(RCSP+(1.0_r8-RCSP)*RCSC)
        !
        !     FRACTION OF RESIDUAL STALK TO BE REMOBILIZED
        !
        !     FSNCR=fraction of residual stalk to be remobilized
        !     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
        !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
        !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
        !
        IF(RCEK(ielmc).GT.ZERO4Groth_pft(NZ))THEN
          FSNCR=AZMAX1(AMIN1(1.0,RespSenes_node/RCEK(ielmc)))
        ELSE
          FSNCR=1.0_r8
        ENDIF
        !
        !     NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
        !     TO FRACTIONS SET IN 'STARTQ'
        !
        D8310: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            dStalk                                   = FSNCR*AZMAX1(SenecStalkStrutElms_brch(NE,NB,NZ)-RCEK(NE))
            dFall                                    = PlantElmAllocMat4Litr(NE,istalk,M,NZ)*dStalk*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp)
            LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)+dFall
            LitrFallElms_brch(NE,NB,NZ)              = LitrFallElms_brch(NE,NB,NZ) +dFall

            dFall                                   = PlantElmAllocMat4Litr(NE,istalk,M,NZ)*dStalk*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp)
            LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+dFall
            LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ) +dFall
          ENDDO    
        ENDDO D8310

        !
        !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
        !
        !     FSNCR=fraction of residual stalk to be remobilized
        !     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in senescing stalk
        !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
        !     StalkNodeHeight_brch,StalkNodeVertLength_brch=living,senescing internode length
        !
        DO NE=1,NumPlantChemElms
          FStalkSenes                        = FSNCR*SenecStalkStrutElms_brch(NE,NB,NZ)
          StalkStrutElms_brch(NE,NB,NZ)      = AZMAX1(StalkStrutElms_brch(NE,NB,NZ)-FStalkSenes)
          SenecStalkStrutElms_brch(NE,NB,NZ) = AZMAX1(SenecStalkStrutElms_brch(NE,NB,NZ)-FStalkSenes)
        ENDDO
        !update 
        HTNODZ=0._r8
        D8320: DO K=0,MaxNodesPerBranch1
          HTNODZ=AMAX1(HTNODZ,StalkNodeHeight_brch(K,NB,NZ))
        ENDDO D8320

        HTNODZ=AZMAX1(1._r8-FSNCR)*HTNODZ
        D8325: DO K=0,MaxNodesPerBranch1
          StalkNodeHeight_brch(K,NB,NZ)=AMIN1(HTNODZ,StalkNodeHeight_brch(K,NB,NZ))
        ENDDO D8325

        !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
        !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
        !
        !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
        !     RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
        !     FSNCR=fraction of residual stalk to be remobilized
        !     RespSenes_node=remaining node senescence respiration
        !
        RCO2e=RCO2e+FSNCR*RCEK(ielmc)*(1._r8-FracRecyled)
        StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)+FSNCR*RCEK(ielmc)*FracRecyled
        DO NE=2,NumPlantChemElms
          StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+FSNCR*RCEK(NE)
        ENDDO
        RespSenes_node=RespSenes_node-FSNCR*RCEK(ielmc)
      ENDIF
      !
      !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
      !
      IF(RespSenes_node.LE.ZERO4Groth_pft(NZ))cycle
      !
      !     OTHERWISE REMAINING C,N,P IN NODE GOES TO LitrFall
      !
      !     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
      !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
      !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
      !
    ELSE      

      D8315: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          dFall  = PlantElmAllocMat4Litr(NE,istalk,M,NZ)*SenecStalkStrutElms_brch(NE,NB,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp)

          LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)              = LitrFallElms_brch(NE,NB,NZ) +dFall

          dFall  = PlantElmAllocMat4Litr(NE,istalk,M,NZ)*SenecStalkStrutElms_brch(NE,NB,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp)

          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) +dFall
          LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ) +dFall
        ENDDO    
      ENDDO D8315

      DO NE=1,NumPlantChemElms
        StalkStrutElms_brch(NE,NB,NZ)      = AZMAX1(StalkStrutElms_brch(NE,NB,NZ)-SenecStalkStrutElms_brch(NE,NB,NZ))
        SenecStalkStrutElms_brch(NE,NB,NZ) = 0._r8
      ENDDO
      !
      !     REMOBILIZATION OF STORAGE C,N,P
      !
      !     WTRVC=storage C
      !     RMxess_brch=remaining excess maintenance respiration
      !
      RMxess_brch=RespSenes_node/(1.0_r8+FXFS)
      IF(SeasonalNonstElms_pft(ielmc,NZ).GT.RMxess_brch)THEN
        SeasonalNonstElms_pft(ielmc,NZ) = SeasonalNonstElms_pft(ielmc,NZ)-RMxess_brch
        SSXferElms_pft(ielmc,NZ)        = SSXferElms_pft(ielmc,NZ)-RMxess_brch
        SSXfer2ShootElms_pft(ielmc,NZ)  = SSXfer2ShootElms_pft(ielmc,NZ)-RMxess_brch
        RCO2e                           = RCO2e+RMxess_brch
        RMxess_brch                     = 0._r8
      ELSEIF(iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial)THEN
        isPlantBranchAlive_brch(NB,NZ)=iFalse
      ENDIF      
    ENDIF

  ENDDO D575
  call PrintInfo('end '//subname)
  end associate
  end subroutine RemobilizeLeafLayers

!----------------------------------------------------------------------------------------------------
  subroutine AllocateLeaf2CanopyLayers(I,J,NB,NZ,CanopyHeight_copy)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NB,NZ
  real(r8), intent(in) :: CanopyHeight_copy(JP1)
  character(len=*), parameter :: subname='AllocateLeaf2CanopyLayers'

  integer  :: LL,LU,L,K,k1,k2,KK,NE
  integer  :: KMinGroingLeafNodeNum,KLeafNumHighestGrowing
  integer  :: LNumHeightLeafTip,LNumHeightLeafBase
  integer  :: LNumHeightBranchTip,LNumHeightBranchBase,N
  real(r8) :: DTransportTube        !thicknees of stalk out layer where xylem and pholem are located, [m]
  real(r8) :: YLeafElmnt_node(NumPlantChemElms)
  real(r8) :: YLeafArea_node,LeafElevation,LeafLength
  real(r8) :: StalkSurfArea,StalkSectionArea
  real(r8) :: FRACL,FRACLS(NumCanopyLayers1)
  real(r8) :: HeightBranchBase
  real(r8) :: HeightStalk
  real(r8) :: HeightLeafNode,HeightLeafLow,HeightLeafTip
  real(r8) :: HeightLeafBase
  real(r8) :: StalkRadius
  real(r8) :: TotLeafElevation
! begin_execution
  associate(                                                           &
    LeafElmntNode_brch        => plt_biom%LeafElmntNode_brch          ,& !input  :leaf element, [g d-2]
    StalkStrutElms_brch       => plt_biom%StalkStrutElms_brch         ,& !input  :branch stalk structural element mass, [g d-2]
    PetoleProteinC_node       => plt_biom%PetoleProteinC_node         ,& !input  :layer sheath protein C, [g d-2]
    FracGroth2Node_pft        => plt_allom%FracGroth2Node_pft         ,& !input  :parameter for allocation of growth to nodes, [-]
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft      ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlantPhenolPattern_pft   => plt_pheno%iPlantPhenolPattern_pft    ,& !input  :plant growth habit: annual or perennial,[-]
    KHiestGroLeafNode_brch    => plt_pheno%KHiestGroLeafNode_brch     ,& !input  :leaf growth stage counter, [-]
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft  ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    rLen2WidthLeaf_pft        => plt_morph%rLen2WidthLeaf_pft         ,& !input  :leaf length:width ratio, [-]
    SeedDepth_pft             => plt_morph%SeedDepth_pft              ,& !input  :seeding depth, [m]
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft          ,& !input  :id number of main branch,[-]
    PetoleLength_node         => plt_morph%PetoleLength_node          ,& !input  :sheath height, [m]
    BranchNumerID_brch        => plt_morph%BranchNumerID_brch         ,& !input  :branch meric id, [-]
    CanopyHeightZ_col         => plt_morph%CanopyHeightZ_col          ,& !input  :canopy layer height, [m]
    LeafAngleClass_pft        => plt_morph%LeafAngleClass_pft         ,& !input  :fractionction of leaves in different angle classes, [-]
    LeafArea_node             => plt_morph%LeafArea_node              ,& !input  :leaf area, [m2 d-2]
    iPlant2ndGrothPattern_pft => plt_pheno%iPlant2ndGrothPattern_pft  ,& !input  :plant expression of secondary growth, [-]                
    PlantPopuLive_pft         => plt_site%PlantPopuLive_pft           ,& !input  :plant population, [d-2]
    StemSpecVolume_pft        => plt_morph%StemSpecVolume_pft         ,& !input  :stalk specific volume, [m3 gC-1]    
    ZERO                      => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch          ,& !input  :branch leaf structural element mass, [g d-2]    
    SineLeafAngle             => plt_rad%SineLeafAngle                ,& !input  :sine of leaf angle,[-]
    LeafLayerElms_node        => plt_biom%LeafLayerElms_node          ,& !inoput :layer leaf element, [g d-2]
    CanopyLeafCLyr_pft        => plt_biom%CanopyLeafCLyr_pft          ,& !inoput :canopy layer leaf C, [g d-2]
    CanopyLeafArea_lnode      => plt_morph%CanopyLeafArea_lnode       ,& !inoput :layer/node/branch leaf area, [m2 d-2]
    StalkNodeHeight_brch      => plt_morph%StalkNodeHeight_brch       ,& !inoput :internode height, [m]
    CanopyHeightLive_pft      => plt_morph%CanopyHeightLive_pft       ,& !inoput :canopy height, [m]
    StalkHeight_pft           => plt_morph%StalkHeight_pft            ,& !inoput :stalk Height/length, [m]
    HypocotHeight_pft         => plt_morph%HypocotHeight_pft          ,& !inoput :cotyledon height (>0), [m]
    CanopyLeafAreaZ_pft       => plt_morph%CanopyLeafAreaZ_pft        ,& !inoput :canopy layer leaf area, [m2 d-2]
    SapwoodBiomassC_brch      => plt_biom%SapwoodBiomassC_brch        ,& !output :branch sapwood C, [gC d-2]
    KLowestGroLeafNode_brch   => plt_pheno%KLowestGroLeafNode_brch    ,& !output :leaf growth stage counter, [-]
    StalkAveRadius_pft     => plt_morph%StalkAveRadius_pft      ,& !output :main stalk radius,[m]
    CanopyStalkSurfArea_lbrch => plt_morph%CanopyStalkSurfArea_lbrch   & !output :plant canopy layer branch stem area, [m2 d-2]
  )
  call DebugPrint('beg '//subname//' NZ',NZ)
  !   ALLOCATION OF LEAF AREA TO CANOPY LAYERS
  !
  !   HypocotHeight_pft=hypocotyledon height, it emerges from the seeding depth, gradually upwards, and eventually getting out of the soil, i.e. HypocotHeight_pft(NZ).GE.SeedDepth_pft(NZ).
  !   SeedDepth_pft=seeding depth
  !   LeafArea_node=node leaf area
  !   PetoleLength_node=PetolSheth length
  !
  KLowestGroLeafNode_brch(NB,NZ)=0;LeafLength=0._r8
  IF(NB.EQ.MainBranchNum_pft(NZ))StalkAveRadius_pft(NZ)=0._r8
  
  IF(HypocotHeight_pft(NZ).LE.SeedDepth_pft(NZ) .AND. LeafArea_node(0,MainBranchNum_pft(NZ),NZ).GT.0.0_r8)THEN
    !plant not emerged yet
    LeafLength            = SQRT(1.0E+02_r8*LeafArea_node(0,MainBranchNum_pft(NZ),NZ)/PlantPopuLive_pft(NZ))
    !pushed upward, 
    HypocotHeight_pft(NZ) = LeafLength+PetoleLength_node(0,MainBranchNum_pft(NZ),NZ)+StalkNodeHeight_brch(0,MainBranchNum_pft(NZ),NZ)
  ENDIF
  !
  ! IF CANOPY HAS EMERGED
  !  
  IF(HypocotHeight_pft(NZ).GT.SeedDepth_pft(NZ))THEN
    D540: DO K=0,MaxNodesPerBranch1
      DO  L=1,NumCanopyLayers1
        CanopyLeafArea_lnode(L,K,NB,NZ)                  = 0._r8
        LeafLayerElms_node(1:NumPlantChemElms,L,K,NB,NZ) = 0._r8
      enddo
    ENDDO D540

    D535: DO L=1,NumCanopyLayers1
      CanopyStalkSurfArea_lbrch(L,NB,NZ)=0._r8
    ENDDO D535
    !
    !   BRANCH HEIGHT
    !
    !   iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+PetolSheth,2=none,3=between 1,2
    !   iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
    !   KLeafNumHighestGrowing,KLowestGroLeafNode_brch=integer of highest,lowest leaf number currently growing
    !   StalkNodeHeight_brch=internode length
    !   HeightBranchBase=branch base height
    !
    IF(iPlantTurnoverPattern_pft(NZ).NE.0 .AND. is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))THEN
      IF(NB.NE.MainBranchNum_pft(NZ))THEN
        KLeafNumHighestGrowing=MAX(1,KHiestGroLeafNode_brch(MainBranchNum_pft(NZ),NZ)-MaxNodesPerBranch1+1)
        IF(BranchNumerID_brch(NB,NZ).GE.KLeafNumHighestGrowing)THEN
          K                = pMOD(BranchNumerID_brch(NB,NZ),MaxNodesPerBranch1)
          HeightBranchBase = StalkNodeHeight_brch(K,MainBranchNum_pft(NZ),NZ)
        ELSE
          HeightBranchBase=0._r8
        ENDIF
      ELSE
        HeightBranchBase=0._r8
      ENDIF
    ELSE
      HeightBranchBase=0._r8
    ENDIF
    KMinGroingLeafNodeNum=MAX(0,KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+1)
    !
    !   FOR ALL LEAFED NODES
    !
    D560: DO KK=KMinGroingLeafNodeNum,KHiestGroLeafNode_brch(NB,NZ)
      K=pMOD(KK,MaxNodesPerBranch1)
      !
      !     HEIGHT OF STALK INTERNODE + SHEATH OR PetolSheth
      !     AND LENGTH OF LEAF
      !
      !     HeightStalk=stalk height
      !     StalkNodeHeight_brch=internode length
      !     HeightLeafNode=leaf node height
      !     LeafArea_node=leaf node area
      !     PP=plant population
      !     FracGroth2Node_pft=scales node number for perennial vegetation (e.g. trees)
      !     LeafLength=leaf length
      !
      HeightStalk         = HeightBranchBase+StalkNodeHeight_brch(K,NB,NZ)
      StalkHeight_pft(NZ) = AMAX1(StalkHeight_pft(NZ),HeightStalk)
      HeightLeafNode      = HeightStalk+PetoleLength_node(K,NB,NZ)
      LeafLength          = AZMAX1(SQRT(rLen2WidthLeaf_pft(NZ)*AZMAX1(LeafArea_node(K,NB,NZ))/(PlantPopuLive_pft(NZ)*FracGroth2Node_pft(NZ))))
      !
      !   ALLOCATE FRACTIONS OF LEAF IN EACH INCLINATION CLASS
      !     FROM HIGHEST TO LOWEST TO CANOPY LAYER
      !
      !     LeafElevation=leaf elevation
      !     CLASS=leaf inclination class
      !     LeafLength=leaf length
      !     ZC,CanopyHeight_copy=canopy height
      !     HeightLeafBase,HeightLeafTip=height of leaf base,tip
      !     ZL=height to bottom of each canopy layer
      !     LNumHeightLeafBase,LNumHeightLeafTip=layer number of leaf base,tip
      !     FRACL=leaf fraction in each layer
      !
      TotLeafElevation = 0._r8
      D555: DO N=NumLeafInclinationClasses1,1,-1
        !assuming the maximum growth rate in 30 mm/hour, note bamboo can growth as fast as 38 mm per hour.
        LeafElevation = SineLeafAngle(N)*LeafAngleClass_pft(N,NZ)*LeafLength
        HeightLeafLow = AMIN1(CanopyHeight_copy(NZ)+0.01_r8-LeafElevation,HeightLeafNode+TotLeafElevation)
        HeightLeafTip = AMIN1(CanopyHeight_copy(NZ)+0.01_r8,HeightLeafLow+LeafElevation)
        !
        LU = 0;LL = 0
        D550: DO L = NumCanopyLayers1, 1, -1
          IF(LU.EQ.1 .AND. LL.EQ.1)exit
          IF((HeightLeafTip.GT.CanopyHeightZ_col(L-1) .OR. CanopyHeightZ_col(L-1).LE.ZERO) .AND. LU.EQ.0)THEN
            LNumHeightLeafTip = MAX(1,L)
            LU                = 1
          ENDIF
          IF((HeightLeafLow.GT.CanopyHeightZ_col(L-1) .OR. CanopyHeightZ_col(L-1).LE.ZERO) .AND. LL.EQ.0)THEN
            LNumHeightLeafBase = MAX(1,L)
            LL                 = 1
          ENDIF
        ENDDO D550
        !
        FRACLS=0._r8
        D570: DO L=LNumHeightLeafBase,LNumHeightLeafTip
          IF(LNumHeightLeafTip.EQ.LNumHeightLeafBase)THEN
            FRACL=LeafAngleClass_pft(N,NZ)
          ELSEIF(HeightLeafTip.GT.HeightLeafLow .AND. CanopyHeightZ_col(L).GT.HeightLeafLow)THEN
            FRACLS(L)=(AMIN1(HeightLeafTip,CanopyHeightZ_col(L))-AMAX1(HeightLeafLow,CanopyHeightZ_col(L-1)))/(HeightLeafTip-HeightLeafLow)
            FRACL=LeafAngleClass_pft(N,NZ)*FRACLS(L)
          ELSE
            FRACL=LeafAngleClass_pft(N,NZ)
          ENDIF
          !
          YLeafArea_node=FRACL*LeafArea_node(K,NB,NZ)
          DO NE=1,NumPlantChemElms
            YLeafElmnt_node(NE)=FRACL*LeafElmntNode_brch(NE,K,NB,NZ)
          ENDDO

          !     ACCUMULATE LAYER LEAF AREAS, C, N AND P CONTENTS
          !
          !     CanopyLeafArea_lnode=leaf node area in canopy layer
          !     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
          !     CanopyLeafAreaZ_pft,CanopyLeafCLyr_pft=total leaf area,C in canopy layer
          !     StalkNodeHeight_brch=internode length
          !
          CanopyLeafArea_lnode(L,K,NB,NZ)=CanopyLeafArea_lnode(L,K,NB,NZ)+YLeafArea_node
          !        
          DO NE=1,NumPlantChemElms
            LeafLayerElms_node(NE,L,K,NB,NZ)=LeafLayerElms_node(NE,L,K,NB,NZ)+YLeafElmnt_node(NE)
          ENDDO
          CanopyLeafAreaZ_pft(L,NZ) = CanopyLeafAreaZ_pft(L,NZ)+YLeafArea_node
          CanopyLeafCLyr_pft(L,NZ)  = CanopyLeafCLyr_pft(L,NZ)+YLeafElmnt_node(ielmc)
        ENDDO D570
        TotLeafElevation         = TotLeafElevation+LeafElevation
        CanopyHeightLive_pft(NZ) = AMAX1(CanopyHeightLive_pft(NZ),HeightLeafTip)
      ENDDO D555

      IF(PetoleProteinC_node(K,NB,NZ).GT.0.0_r8)THEN
        IF(KLowestGroLeafNode_brch(NB,NZ).EQ.0)KLowestGroLeafNode_brch(NB,NZ)=MIN(KK,KHiestGroLeafNode_brch(NB,NZ))
      ENDIF
    ENDDO D560

    IF(KLowestGroLeafNode_brch(NB,NZ).EQ.0)KLowestGroLeafNode_brch(NB,NZ)=KHiestGroLeafNode_brch(NB,NZ)
    K1 = pMOD(KHiestGroLeafNode_brch(NB,NZ),MaxNodesPerBranch1)
    K2 = pMOD(KHiestGroLeafNode_brch(NB,NZ)-1,MaxNodesPerBranch1)
    !
    IF(isclose(StalkNodeHeight_brch(K1,NB,NZ),0._r8))THEN
      StalkNodeHeight_brch(K1,NB,NZ)=StalkNodeHeight_brch(K2,NB,NZ)
    ENDIF

    HeightLeafBase=HeightBranchBase+AZMAX1(StalkNodeHeight_brch(K1,NB,NZ))
    !
    !     ALLOCATE STALK SURFACE AREA TO CANOPY LAYERS
    !
    !     StalkNodeHeight_brch=internode length
    !     HeightLeafBase=leaf base height
    !     ZL=height to bottom of each canopy layer
    !     LNumHeightBranchBase,LNumHeightBranchTip=layer number of branch base,tip
    !     WTSTKB,StalkSurfArea=branch stalk mass,surface area
    !     FSTK=fraction of stalk area contributing to water,heat flow
    !     SapwoodBiomassC_brch=stalk sapwood mass, the living, outermost portion of a tree trunk or branch. 
    !          It lies just beneath the bark and cambium layer but surrounds the darker, inner "heartwood
    !     FRACL=stalk fraction in each layer
    !     CanopyStalkSurfArea_lbrch=total branch stalk surface area in each layer
    !
    IF(StalkNodeHeight_brch(K1,NB,NZ).GT.0.0_r8)THEN
      LU=0;LL=0
      D545: DO L=NumCanopyLayers1,1,-1
        IF(LU.EQ.1 .AND. LL.EQ.1)exit
        IF((HeightLeafBase.GT.CanopyHeightZ_col(L-1) .OR. CanopyHeightZ_col(L-1).LE.ZERO) .AND. LU.EQ.0)THEN
          LNumHeightBranchTip = MAX(1,L)
          LU                  = 1
        ENDIF
        IF((HeightBranchBase.GT.CanopyHeightZ_col(L-1).OR.CanopyHeightZ_col(L-1).LE.ZERO) .AND. LL.EQ.0)THEN
          LNumHeightBranchBase = MAX(1,L)
          LL                   = 1
        ENDIF
      ENDDO D545      
      !
      StalkRadius=SQRT(StemSpecVolume_pft(NZ)*(AZMAX1(StalkStrutElms_brch(ielmc,NB,NZ))/PlantPopuLive_pft(NZ))/(PICON*StalkNodeHeight_brch(K1,NB,NZ)))
      
      IF(NB.EQ.MainBranchNum_pft(NZ))StalkAveRadius_pft(NZ)=StalkRadius
      
      !assume stalk is cylindrical
      StalkSurfArea = PICON*StalkNodeHeight_brch(K1,NB,NZ)*StalkRadius*PlantPopuLive_pft(NZ)

      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
        SapwoodBiomassC_brch(NB,NZ)=StalkStrutElms_brch(ielmc,NB,NZ)
      ELSE        
        DTransportTube     = AMIN1(ZSTX,FSTK*StalkRadius)
        !the following is based on the observation that xylem and pholem are limited to a thicknees of DTranpTube on the outside of stalk.
        !the definition of sapwood here is not exact. For a tree, from the core to bark, they are heartwood, sapwood, vascular cambium, 
        !living phloem, cork cambium and cork. 
        StalkSectionArea            = PICON*(2._r8*StalkRadius*DTransportTube-DTransportTube**2)
        SapwoodBiomassC_brch(NB,NZ) = StalkSectionArea/StemSpecVolume_pft(NZ)*StalkNodeHeight_brch(K1,NB,NZ)*PlantPopuLive_pft(NZ)
      ENDIF

      D445: DO L=LNumHeightBranchBase,LNumHeightBranchTip
        IF(HeightLeafBase.GT.HeightBranchBase)THEN
          IF(HeightLeafBase.GT.CanopyHeightZ_col(L-1))THEN
            FRACL=(AMIN1(HeightLeafBase,CanopyHeightZ_col(L))-AMAX1(HeightBranchBase,CanopyHeightZ_col(L-1))) &
              /(HeightLeafBase-HeightBranchBase)
          ELSE
            FRACL=0._r8
          ENDIF
        ELSE
          FRACL=1.0_r8
        ENDIF
        CanopyStalkSurfArea_lbrch(L,NB,NZ)=FRACL*StalkSurfArea
      ENDDO D445
    ELSE
      SapwoodBiomassC_brch(NB,NZ)=0._r8
      D450: DO L=1,NumCanopyLayers1
        CanopyStalkSurfArea_lbrch(L,NB,NZ)=0._r8
      ENDDO D450
    ENDIF
  ELSE
    L=1;K=0
    SapwoodBiomassC_brch(NB,NZ)                     = 0._r8
    CanopyStalkSurfArea_lbrch(1:NumCanopyLayers1,NB,NZ) = 0._r8
    CanopyLeafArea_lnode(L,K,NB,NZ)                 = LeafArea_node(0,NB,NZ)
    LeafLayerElms_node(:,L,K,NB,NZ)                 = LeafStrutElms_brch(:,NB,NZ)
  ENDIF

  call PrintInfo('end '//subname)
  end associate
  end subroutine AllocateLeaf2CanopyLayers

!----------------------------------------------------------------------------------------------------
  subroutine LeafClassAllocation(NB,NZ)
  implicit none
  integer, intent(in) :: NB,NZ
  character(len=*), parameter :: subname='LeafClassAllocation'
  real(r8) :: dangle
  integer :: L,K,N
  ! begin_execution
  associate(                                                           &
    SineBranchAngle_pft       => plt_morph%SineBranchAngle_pft        ,& !input  :branching angle, [degree from horizontal]
    LeafArea_node             => plt_morph%LeafArea_node              ,& !input  :leaf area, [m2 d-2]
    CanopyStalkSurfArea_lbrch => plt_morph%CanopyStalkSurfArea_lbrch  ,& !input  :plant canopy layer branch stem area, [m2 d-2]
    LeafAngleClass_pft        => plt_morph%LeafAngleClass_pft         ,& !input  :fractionction of leaves in different angle classes, [-]
    CanopyLeafArea_lnode      => plt_morph%CanopyLeafArea_lnode       ,& !input  :layer/node/branch leaf area, [m2 d-2]
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft          ,& !input  :number of main branch,[-]
    StemAreaZsec_brch         => plt_morph%StemAreaZsec_brch          ,& !output :stem surface area, [m2 d-2]
    LeafAreaZsec_brch         => plt_morph%LeafAreaZsec_brch           & !output :leaf surface area, [m2 d-2]
  )
  call PrintInfo('beg '//subname)
  D900: DO K=1,MaxNodesPerBranch1
    DO  L=1,NumCanopyLayers1
      DO  N=1,NumLeafInclinationClasses1
        LeafAreaZsec_brch(N,L,K,NB,NZ)=0._r8
      enddo
    enddo
  ENDDO D900

  D500: DO K=1,MaxNodesPerBranch1
    IF(LeafArea_node(K,NB,NZ).GT.0.0_r8)THEN
      D700: DO L=NumCanopyLayers1,1,-1
        D800: DO N=1,NumLeafInclinationClasses1
          LeafAreaZsec_brch(N,L,K,NB,NZ)=AZMAX1(LeafAngleClass_pft(N,NZ)*CanopyLeafArea_lnode(L,K,NB,NZ))/REAL(NumOfLeafAzimuthSectors1,R8)
        ENDDO D800
      ENDDO D700
    ENDIF
  ENDDO D500
  !
  ! ALLOCATE STALK AREA TO INCLINATION CLASSES ACCORDING TO
  ! BRANCH ANGLE ENTERED IN 'READQ' ASSUMING AZIMUTH IS UNIFORM
  !
  ! StemAreaZsec_brch=stalk surface area in canopy layer
  ! SineBranchAngle_pft=stem angle from horizontal
  ! CanopyStalkSurfArea_lbrch=total branch stalk surface area in each layer
  !
  D910: DO L=1,NumCanopyLayers1
    DO  N=1,NumLeafInclinationClasses1
      StemAreaZsec_brch(N,L,NB,NZ)=0._r8
    enddo
  ENDDO D910

  !locate the stem angle class id for the branch, 1 means horizontal, NumLeafInclinationClasses1 means vertical
  !the main branch always vertical
  IF(NB.EQ.MainBranchNum_pft(NZ))THEN
    N=NumLeafInclinationClasses1
  ELSE
    dangle = PICON2h/real(NumLeafInclinationClasses1,r8)       !pi/2=PICON2h 
    N      = MIN(NumLeafInclinationClasses1,INT(ASIN(SineBranchAngle_pft(NZ))/dangle)+1)
  ENDIF
  D710: DO L=NumCanopyLayers1,1,-1
    StemAreaZsec_brch(N,L,NB,NZ)=CanopyStalkSurfArea_lbrch(L,NB,NZ)/real(NumOfLeafAzimuthSectors1,r8)
  ENDDO D710
  call PrintInfo('end '//subname)
  end associate
  end subroutine LeafClassAllocation

!----------------------------------------------------------------------------------------------------
  subroutine GrainFillOnBranch(I,J,NB,NZ,GrainGroth,GROSTKC)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: GrainGroth(NumPlantChemElms)
  real(r8), intent(in) :: GROSTKC                       !stalk growth rate, [gC/hr]
  character(len=*), parameter :: subname='GrainFillOnBranch'
  real(r8) :: ZPGRP,ZPGRN,ZNPGP,ZNPGN
  real(r8) :: Reserve2GrainCMax,Reserve2GrainE(NumPlantChemElms)
  real(r8) :: FGRNX
  real(r8) :: GRMXB
  real(r8) :: GrainFillRateMax
  real(r8) :: GrainFillRate
  real(r8) :: GrainCRemaining                   !unfilled grain carbon capacity, [gC d-2]
  real(r8) :: SeedSET
  integer :: NE
! begin_execution
  associate(                                                                 &
    TdegCCanopy_pft              => plt_ew%TdegCCanopy_pft                  ,& !input  :canopy temperature, [oC]
    LeafPetoNonstElmConc_brch    => plt_biom%LeafPetoNonstElmConc_brch      ,& !input  :branch nonstructural C concentration, [g d-2]
    ZERO4Groth_pft               => plt_biom%ZERO4Groth_pft                 ,& !input  :threshold zero for plang growth calculation, [-]
    rNCGrain_pft                 => plt_allom%rNCGrain_pft                  ,& !input  :grain N:C ratio, [g g-1]
    rPCGrain_pft                 => plt_allom%rPCGrain_pft                  ,& !input  :grain P:C ratio, [gP gP-1]
    fTgrowRootP_vr               => plt_pheno%fTgrowRootP_vr                ,& !input  :root layer temperature growth functiom, [-]
    fTCanopyGroth_pft            => plt_pheno%fTCanopyGroth_pft             ,& !input  :canopy temperature growth function, [-]
    iPlantPhenolType_pft         => plt_pheno%iPlantPhenolType_pft          ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    HourReq4LeafOff_brch         => plt_pheno%HourReq4LeafOff_brch          ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    GrainFillRate25C_pft         => plt_pheno%GrainFillRate25C_pft          ,& !input  :maximum rate of fill per grain, [g h-1]
    SeedTempSens_pft             => plt_pheno%SeedTempSens_pft              ,& !input  :sensitivity to HTC (seeds oC-1 above HTC),[oC-1]
    iPlantPhenolPattern_pft      => plt_pheno%iPlantPhenolPattern_pft       ,& !input  :plant growth habit: annual or perennial,[-]
    HighTempLimitSeed_pft        => plt_pheno%HighTempLimitSeed_pft         ,& !input  :temperature above which seed set is adversely affected, [oC]
    dReproNodeNumNormByMatG_brch => plt_pheno%dReproNodeNumNormByMatG_brch  ,& !input  :gain in normalized node number during reproductive growth stages, [h-1]
    TCChill4Seed_pft             => plt_pheno%TCChill4Seed_pft              ,& !input  :temperature below which seed set is adversely affected, [oC]
    PlantPopuLive_pft            => plt_site%PlantPopuLive_pft              ,& !input  :plant population, [d-2]
    MaxSeedNumPerSite_pft        => plt_morph%MaxSeedNumPerSite_pft         ,& !input  :maximum grain number per node, [-]
    SeedCMassMax_pft             => plt_morph%SeedCMassMax_pft              ,& !input  :maximum grain size, [g]
    GrothStalkMaxSeedSites_pft   => plt_morph%GrothStalkMaxSeedSites_pft    ,& !input  :maximum grain sites per unit of stalk/branch growth, [-]
    NGTopRootLayer_pft           => plt_morph%NGTopRootLayer_pft            ,& !input  :soil layer at planting depth, [-]
    iPlantGrainType_pft          => plt_morph%iPlantGrainType_pft           ,& !input  :grain type (below or above-ground),[-]
    GrainStrutElms_brch          => plt_biom%GrainStrutElms_brch            ,& !inoput :branch grain structural element mass, [g d-2]
    StalkRsrvElms_brch           => plt_biom%StalkRsrvElms_brch             ,& !inoput :branch reserve element mass, [g d-2]
    SingleGrainMeanBiomC_brch      => plt_allom%SingleGrainMeanBiomC_brch       ,& !inoput :potential carbon mass per grain, [gC seed-1]
    HourFailGrainFill_brch       => plt_pheno%HourFailGrainFill_brch        ,& !inoput :flag to detect physiological maturity from grain fill, [-]
    iPlantCalendar_brch          => plt_pheno%iPlantCalendar_brch           ,& !inoput :plant growth stage, [-]
    PotentialSeedSites_brch      => plt_morph%PotentialSeedSites_brch       ,& !inoput :branch potential grain number, [d-2]
    Hours4LeafOff_brch           => plt_pheno%Hours4LeafOff_brch            ,& !output :cold requirement for autumn leafoff/hardening, [h]
    SetNumberSeeds_brch          => plt_morph%SetNumberSeeds_brch            & !output :branch grain sites number, [d-2]
  )
  call DebugPrint('beg '//subname//' NZ',NZ)
  !
  !   SET MAXIMUM GRAIN NUMBER FROM SHOOT MASS BEFORE ANTHESIS
  !
  !   iPlantCalendar_brch(ipltcal_Jointing,=start of stem elongation and setting max seed number
  !   iPlantCalendar_brch(ipltcal_Anthesis,=start of anthesis and setting final seed number
  !   PotentialSeedSites_brch=potential number of seed set sites
  !   GrothStalkMaxSeedSites_pft=maximum potential seed number from PFT file
  !     GROSTKC=stalk growth rate
  !
  IF(iPlantCalendar_brch(ipltcal_Jointing,NB,NZ).NE.0 .AND. iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).EQ.0)THEN
    PotentialSeedSites_brch(NB,NZ)=PotentialSeedSites_brch(NB,NZ)+GrothStalkMaxSeedSites_pft(NZ)*AZMAX1(GROSTKC)
  ENDIF
  !
  !   SET FINAL GRAIN NUMBER FROM C,N,P NON-STRUCTURAL POOLS AFTER ANTHESIS
  !
  !   iPlantCalendar_brch(ipltcal_Anthesis,=start of anthesis and setting final seed number
  !   iPlantCalendar_brch(ipltcal_BeginSeedFill,=start of grain filling and setting max seed size
  !   iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date setting for final seed number
  !   iPlantCalendar_brch(ipltcal_SetSeedMass,=end of setting max seed size
  !   SET=seed set limited by nonstructural C,N,P
  !   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
  !   TdegCCanopy_pft=canopy temperature
  !   TCChill4Seed_pft=chilling temperature for CO2 fixation, seed loss (oC)
  !   HTC=high temperature threshold for grain number loss
  !   FGRNX=loss of seed set
  !   SeedTempSens_pft=sensitivity to TdegCCanopy_pft> HTC,TdegCCanopy_pft< TCChill4Seed_pftfrom startq.f (seeds oC-1)
  !   SetNumberSeeds_brch=seed set number
  !   PotentialSeedSites_brch=potential number of seed set sites
  !   MaxSeedNumPerSite_pft=maximum seed number per GrothStalkMaxSeedSites_pft from PFT file
  !   dReproNodeNumNormByMatG_brch=change in reproductive node number normalized for maturity group
  !
  IF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).NE.0 .AND. iPlantCalendar_brch(ipltcal_SetSeedMass,NB,NZ).EQ.0)THEN !anthesis stage
    SeedSET=AMIN1(LeafPetoNonstElmConc_brch(ielmc,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmc,NB,NZ)+SETC) &
      ,LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmn,NB,NZ)+SETN) &
      ,LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmp,NB,NZ)+SETP))

    IF(TdegCCanopy_pft(NZ).LT.TCChill4Seed_pft(NZ))THEN
      IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN
        FGRNX=SeedTempSens_pft(NZ)*(TCChill4Seed_pft(NZ)-TdegCCanopy_pft(NZ))
      ELSEIF(iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).EQ.0)THEN
        FGRNX=SeedTempSens_pft(NZ)*(TCChill4Seed_pft(NZ)-TdegCCanopy_pft(NZ))
      ELSE
        FGRNX=0._r8
      ENDIF
    ELSEIF(TdegCCanopy_pft(NZ).GT.HighTempLimitSeed_pft(NZ))THEN
      IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN
        FGRNX=SeedTempSens_pft(NZ)*(TdegCCanopy_pft(NZ)-HighTempLimitSeed_pft(NZ))
      ELSEIF(iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).EQ.0)THEN
        FGRNX=SeedTempSens_pft(NZ)*(TdegCCanopy_pft(NZ)-HighTempLimitSeed_pft(NZ))
      ELSE
        FGRNX=0._r8
      ENDIF
    ELSE
      FGRNX=0._r8
    ENDIF

    IF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).NE.0 .AND. iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).EQ.0)THEN
!     PotentialSeedSites_brch(NB,NZ)=PotentialSeedSites_brch(NB,NZ)*FGRNX
      SetNumberSeeds_brch(NB,NZ)=AMIN1(MaxSeedNumPerSite_pft(NZ)*PotentialSeedSites_brch(NB,NZ) &
        ,SetNumberSeeds_brch(NB,NZ)+MaxSeedNumPerSite_pft(NZ)*PotentialSeedSites_brch(NB,NZ) &
        *SeedSET*dReproNodeNumNormByMatG_brch(NB,NZ)-FGRNX*SetNumberSeeds_brch(NB,NZ))
    ENDIF
    !
    !     SET MAXIMUM GRAIN SIZE FROM C,N,P NON-STRUCTURAL POOLS AFTER ANTHESIS
    !
    !     GRMX=maximum individual seed size from PFT file (g)
    !     dReproNodeNumNormByMatG_brch=change in reproductive node number normalized for maturity group
    !     SingleGrainMeanBiomC_brch=individual seed size
    !     SET=seed set limited by nonstructural C,N,P
    !
    IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).NE.0 &
      .AND.iPlantCalendar_brch(ipltcal_SetSeedMass,NB,NZ).EQ.0)THEN
      GRMXB                          = SeedCMassMax_pft(NZ)
      SingleGrainMeanBiomC_brch(NB,NZ) = AMIN1(SeedCMassMax_pft(NZ),SingleGrainMeanBiomC_brch(NB,NZ) &
        +GRMXB*AMAX1(0.50_r8,SeedSET**0.25_r8)*dReproNodeNumNormByMatG_brch(NB,NZ))
    ENDIF
  ENDIF
  !
  !   GRAIN FILL BY TRANSLOCATION FROM STALK RESERVES
  !   UNTIL GRAIN SINK (=FINAL GRAIN NUMBER X MAXIMUM
  !   GRAIN SIZE) IS FILLED OR RESERVES ARE EXHAUSTED
  !
  !   iPlantCalendar_brch(ipltcal_BeginSeedFill,=start of grain filling and setting max seed size
  !   WTGRB=total seed C mass
  !   SingleGrainMeanBiomC_brch=individual seed size
  !   SetNumberSeeds_brch=seed set number
  !   GrainFillRate25C_pft=grain filling rate at 25 oC from PFT file
  !   fTCanopyGroth_pft=temperature function for canopy growth
  !   fTgrowRootP_vr=temperature function for root growth
  !
  IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).NE.0)THEN
    GrainCRemaining = AZMAX1(SingleGrainMeanBiomC_brch(NB,NZ)*SetNumberSeeds_brch(NB,NZ) &
      -GrainStrutElms_brch(ielmc,NB,NZ))
    IF(GrainCRemaining.LE.0._r8)THEN
      GrainFillRateMax=0._r8
    ELSEIF(iPlantGrainType_pft(NZ).EQ.igraintyp_abvgrnd)THEN
      GrainFillRateMax=AZMAX1(GrainFillRate25C_pft(NZ)*SetNumberSeeds_brch(NB,NZ)*SQRT(fTCanopyGroth_pft(NZ)))
    ELSE
      !seed/storage belowground 
      GrainFillRateMax=AZMAX1(GrainFillRate25C_pft(NZ)*SetNumberSeeds_brch(NB,NZ)*SQRT(fTgrowRootP_vr(NGTopRootLayer_pft(NZ),NZ)))
    ENDIF
    ! Limit this hourly transfer to the remaining sink before applying C/N/P constraints.
    GrainFillRateMax = AMIN1(GrainFillRateMax,GrainCRemaining)
!
!     GRAIN FILL RATE MAY BE CONSTRAINED BY HIGH GRAIN C:N OR C:P
!
!     WTGRB,WTGRBN,WTGRBP=total seed C,N,P mass
!     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
!     CNGR,CPGR=maximum N:C,P:C ratios in grain from PFT file
!     Reserve2GrainCMax,Reserve2GrainE(ielmc)=maximum,actual C translocation rate from reserve to grain
!
    IF(GrainStrutElms_brch(ielmn,NB,NZ).LT.ZPGRM*rNCGrain_pft(NZ)*GrainStrutElms_brch(ielmc,NB,NZ) &
      .OR. GrainStrutElms_brch(ielmp,NB,NZ).LT.ZPGRM*rPCGrain_pft(NZ)*GrainStrutElms_brch(ielmc,NB,NZ))THEN
      GrainFillRate=0._r8
    ELSE
      GrainFillRate=GrainFillRateMax
    ENDIF
    Reserve2GrainCMax     = AMIN1(GrainFillRateMax,StalkRsrvElms_brch(ielmc,NB,NZ))
    Reserve2GrainE(ielmc) = AMIN1(GrainFillRate,StalkRsrvElms_brch(ielmc,NB,NZ))
    !
    !     GRAIN N OR P FILL RATE MAY BE LIMITED BY C:N OR C:P RATIOS
    !     OF STALK RESERVES
    !
    !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
    !     ZNPGN,ZNPGP=effect of reserve N:C,P:C on grain fill N:C,P:C
    !     SETN,SETP=Km for nonstructural N,P concn on seed set (g g-1)
    !     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
    !     ZPGRD=1.0_r8-ZPGRM
    !     ZPGRN,ZPGRP=N:C,P:C ratios during grain fill
    !     Reserve2GrainCMax,Reserve2GrainE(ielmc)=maximum,actual C translocation rate from reserve to grain
    !     CNGR,CPGR=maximum N:C,P:C ratios in grain from PFT file
    !     Reserve2GrainE(ielmn:ielmp)=N,P translocation rate from reserve to grain
    !
    IF(StalkRsrvElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      ZNPGN = StalkRsrvElms_brch(ielmn,NB,NZ)/(StalkRsrvElms_brch(ielmn,NB,NZ)+SETN*StalkRsrvElms_brch(ielmc,NB,NZ))
      ZNPGP = StalkRsrvElms_brch(ielmp,NB,NZ)/(StalkRsrvElms_brch(ielmp,NB,NZ)+SETP*StalkRsrvElms_brch(ielmc,NB,NZ))
      ZPGRN = ZPGRM+ZPGRD*AZMAX1(AMIN1(1.0_r8,ZNPGN))
      ZPGRP = ZPGRM+ZPGRD*AZMAX1(AMIN1(1.0_r8,ZNPGP))
      Reserve2GrainE(ielmn) = AMIN1(Reserve2GrainCMax*rNCGrain_pft(NZ) &
        ,AZMAX1(StalkRsrvElms_brch(ielmn,NB,NZ)*ZPGRN) &
        ,(GrainStrutElms_brch(ielmc,NB,NZ)+Reserve2GrainE(ielmc))*rNCGrain_pft(NZ) &
        -GrainStrutElms_brch(ielmn,NB,NZ))
      Reserve2GrainE(ielmp)=AMIN1(Reserve2GrainCMax*rPCGrain_pft(NZ) &
        ,AZMAX1(StalkRsrvElms_brch(ielmp,NB,NZ)*ZPGRP) &
        ,(GrainStrutElms_brch(ielmc,NB,NZ)+Reserve2GrainE(ielmc))*rPCGrain_pft(NZ) &
        -GrainStrutElms_brch(ielmp,NB,NZ))
    ELSE
      Reserve2GrainE(ielmn)=0._r8
      Reserve2GrainE(ielmp)=0._r8
    ENDIF
    !
    !     TRANSLOCATE C,N,P FROM STALK RESERVES TO GRAIN
    !
    !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
    !     GROGR=grain growth rate
    !     Reserve2GrainE=C,N,P translocation rate from reserve to grain
    !
    DO NE=1,NumPlantChemElms
      StalkRsrvElms_brch(NE,NB,NZ)  = StalkRsrvElms_brch(NE,NB,NZ)+GrainGroth(NE)-Reserve2GrainE(NE)
      GrainStrutElms_brch(NE,NB,NZ) = GrainStrutElms_brch(NE,NB,NZ)+Reserve2GrainE(NE)
    ENDDO
  ELSE
    Reserve2GrainE(1:NumPlantChemElms)=0._r8
  ENDIF
  !
  !   SET DATE OF PHYSIOLOGICAL MATURITY WHEN GRAIN FILL
  !   HAS STOPPED FOR SET PERIOD OF TIME
  !
  !   iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date setting for final seed number
  !   Reserve2GrainE(ielmc)=C translocation rate from reserve to grain
  !   PP=PFT population
  !   HourFailGrainFill_brch=number of hours with no grain fill
  !   Hours4PhyslMature=number of hours with no grain filling until physl maturity
  !   iPlantCalendar_brch(ipltcal_EndSeedFill,=date of physiological maturity
  !
  IF(iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).NE.0)THEN
    IF(Reserve2GrainE(ielmc).LE.1.0E-09_r8*PlantPopuLive_pft(NZ))THEN
      HourFailGrainFill_brch(NB,NZ)=HourFailGrainFill_brch(NB,NZ)+1.0
    ELSE
      HourFailGrainFill_brch(NB,NZ)=0._r8
    ENDIF
    IF(HourFailGrainFill_brch(NB,NZ).GE.Hours4PhyslMature)THEN
      IF(iPlantCalendar_brch(ipltcal_EndSeedFill,NB,NZ).EQ.0)THEN
        iPlantCalendar_brch(ipltcal_EndSeedFill,NB,NZ)=I
      ENDIF
    ENDIF
    !
    !     TERMINATE ANNUALS AFTER GRAIN FILL
    !
    !     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial
    !     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
    !     HourFailGrainFill_brch=number of hours with no grain fill
    !     Hours4PhyslMature=number of hours with no grain filling until physiological maturity
    !     Hours4SenesAftMature=number of hours after physiol maturity required for senescence
    !     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
    !
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantPhenolType_pft(NZ).NE.iphenotyp_evgreen)THEN
      !real annual
      IF(HourFailGrainFill_brch(NB,NZ).GT.Hours4PhyslMature+Hours4SenesAftMature(iPlantPhenolType_pft(NZ)))THEN
        Hours4LeafOff_brch(NB,NZ)=HourReq4LeafOff_brch(NB,NZ)+0.5_r8
      ENDIF
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine GrainFillOnBranch

!----------------------------------------------------------------------------------------------------
  subroutine ResetNonAnnualBranch(I,J,NB,NZ)      
!  !
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NB,NZ
  integer :: K,M,NE
  real(r8) :: dFall
  character(len=*), parameter :: subname='ResetNonAnnualBranch'

  associate(                                                                           &
    MatureGroup_pft                   => plt_pheno%MatureGroup_pft                    ,& !input  :acclimated plant maturity group, [-]
    icwood                            => pltpar%icwood                                ,& !input  :group id of coarse woody litter
    istalk                            => pltpar%istalk                                ,& !input  :group id of plant stalk litter group
    ifoliar                           => pltpar%ifoliar                               ,& !input  :group id of plant foliar litter
    inonfoliar                        => pltpar%inonfoliar                            ,& !input  :group id of plant non-foliar litter group
    k_fine_comp                       => pltpar%k_fine_comp                           ,& !input  :fine litter complex id
    k_woody_comp                      => pltpar%k_woody_comp                          ,& !input  :woody litter complex id
    MainBranchNum_pft                 => plt_morph%MainBranchNum_pft                  ,& !input  :id number of main branch,[-]
    iPlantPhenolPattern_pft           => plt_pheno%iPlantPhenolPattern_pft            ,& !input  :plant growth habit: annual or perennial,[-]
    Hours4Leafout_brch                => plt_pheno%Hours4Leafout_brch                 ,& !input  :heat requirement for spring leafout/dehardening, [h]
    iPlantTurnoverPattern_pft         => plt_pheno%iPlantTurnoverPattern_pft          ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    PlantElmAllocMat4Litr             => plt_soilchem%PlantElmAllocMat4Litr           ,& !input  :litter kinetic fraction, [-]
    FracPetolShethAlloc2Litr          => plt_allom%FracPetolShethAlloc2Litr           ,& !input  :leaf element allocation,[-]
    FracLeafShethElmAlloc2Litr        => plt_allom%FracLeafShethElmAlloc2Litr         ,& !input  :woody element allocation, [-]
    HourReq4LeafOut_brch              => plt_pheno%HourReq4LeafOut_brch               ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    ShootNodeNumAtPlanting_pft        => plt_morph%ShootNodeNumAtPlanting_pft         ,& !input  :number of nodes in seed, [-]
    iPlantRootProfile_pft             => plt_pheno%iPlantRootProfile_pft              ,& !input  :plant growth type (vascular, non-vascular),[-]
    EnablePlantLeafOut_brch           => plt_pheno%EnablePlantLeafOut_brch            ,& !inoput :branch phenology flag, [-]
    ShootNodeNum_brch                 => plt_morph%ShootNodeNum_brch                  ,& !inoput :shoot node number, [-]
    LitrfallElms_pvr                  => plt_bgcr%LitrfallElms_pvr                    ,& !inoput :plant LitrFall element, [g d-2 h-1]
    LitrFallElms_brch                 => plt_bgcr%LitrFallElms_brch                   ,& !inoput :litterfall from the branch, [g d-2 h-1]    
    LeafStrutElms_brch                => plt_biom%LeafStrutElms_brch                  ,& !inoput :branch leaf structural element mass, [g d-2]
    PetolShethStrutElms_brch          => plt_biom%PetolShethStrutElms_brch            ,& !inoput :branch sheath structural element, [g d-2]
    HuskStrutElms_brch                => plt_biom%HuskStrutElms_brch                  ,& !inoput :branch husk structural element mass, [g d-2]
    EarStrutElms_brch                 => plt_biom%EarStrutElms_brch                   ,& !inoput :branch ear structural chemical element mass, [g d-2]
    GrainStrutElms_brch               => plt_biom%GrainStrutElms_brch                 ,& !inoput :branch grain structural element mass, [g d-2]
    StalkStrutElms_brch               => plt_biom%StalkStrutElms_brch                 ,& !inoput :branch stalk structural element mass, [g d-2]
    iPlant2ndGrothPattern_pft         => plt_pheno%iPlant2ndGrothPattern_pft          ,& !input  :plant expression of secondary growth, [-]                
    ShootNodeNumAtInitFloral_brch     => plt_morph%ShootNodeNumAtInitFloral_brch      ,& !output :shoot node number at plant initiation, [-]
    ShootNodeNumAtAnthesis_brch       => plt_morph%ShootNodeNumAtAnthesis_brch        ,& !output :shoot node number at anthesis, [-]
    LeafNumberAtFloralInit_brch       => plt_pheno%LeafNumberAtFloralInit_brch        ,& !output :leaf number at floral initiation, [-]
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch      ,& !output :normalized node number during vegetative growth stages, [-]
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch  ,& !output :normalized node number during reproductive growth stages, [-]
    dReproNodeNumNormByMatG_brch      => plt_pheno%dReproNodeNumNormByMatG_brch ,& !output :current hourly reproductive development increment, [h-1]
    MatureGroup_brch                  => plt_pheno%MatureGroup_brch                   ,& !output :plant maturity group, [-]
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch                ,& !output :plant growth stage, [-]
    HoursTooLowPsiCan_pft             => plt_pheno%HoursTooLowPsiCan_pft              ,& !output :canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY), [h]
    NumOfLeaves_brch                  => plt_morph%NumOfLeaves_brch                   ,& !output :leaf number, [-]
    KLeafNumber_brch                  => plt_morph%KLeafNumber_brch                   ,& !output :leaf number, [-]
    KHiestGroLeafNode_brch            => plt_pheno%KHiestGroLeafNode_brch             ,& !output :leaf growth stage counter, [-]
    HourFailGrainFill_brch            => plt_pheno%HourFailGrainFill_brch             ,& !output :flag to detect physiological maturity from grain fill, [-]
    LeafAreaLive_brch                 => plt_morph%LeafAreaLive_brch                  ,& !output :branch leaf area, [m2 d-2]
    LeafElmntNode_brch                => plt_biom%LeafElmntNode_brch                  ,& !output :leaf element, [g d-2]
    LeafArea_node                     => plt_morph%LeafArea_node                      ,& !output :leaf area, [m2 d-2]
    PetolShethElmntNode_brch          => plt_biom%PetolShethElmntNode_brch            ,& !output :sheath chemical element, [g d-2]
    PetoleProteinC_node               => plt_biom%PetoleProteinC_node                 ,& !output :layer sheath protein C, [g d-2]
    PetoleLength_node                 => plt_morph%PetoleLength_node                  ,& !output :sheath height, [m]
    LeafProteinC_node                 => plt_biom%LeafProteinC_node                   ,& !output :layer leaf protein C, [g d-2]
    PotentialSeedSites_brch           => plt_morph%PotentialSeedSites_brch            ,& !output :branch potential grain number, [d-2]
    doInitLeafOut_brch                => plt_pheno%doInitLeafOut_brch                 ,& !output :branch phenology flag, [-]
    SingleGrainMeanBiomC_brch           => plt_allom%SingleGrainMeanBiomC_brch            ,& !output :potential carbon mass per grain, [gC seed-1]
    SenecStalkStrutElms_brch          => plt_biom%SenecStalkStrutElms_brch            ,& !output :branch stalk structural element, [g d-2]
    StructInternodeElms_brch          => plt_biom%StructInternodeElms_brch            ,& !output :internode chemical element, [g d-2]
    StalkNodeHeight_brch              => plt_morph%StalkNodeHeight_brch               ,& !output :internode height, [m]
    StalkNodeVertLength_brch          => plt_morph%StalkNodeVertLength_brch           ,& !output :internode height, [m]
    doPlantLeaveOff_brch              => plt_pheno%doPlantLeaveOff_brch               ,& !output :branch phenology flag, [-]
    Prep4Literfall_brch               => plt_pheno%Prep4Literfall_brch                ,& !output :branch phenology flag, [-]
    Hours4LiterfalAftMature_brch      => plt_pheno%Hours4LiterfalAftMature_brch       ,& !output :branch phenology flag, [h]
    SetNumberSeeds_brch               => plt_morph%SetNumberSeeds_brch                 & !output :branch grain number, [d-2]
  )
!    GROUP,MatureGroup_pft=node number required for floral initiation
!    ShootNodeNumAtInitFloral_brch=node number at floral initiation
!    ShootNodeNumAtAnthesis_brch=node number at flowering
!    VSTGX=leaf number on date of floral initiation
!    TotalNodeNumNormByMatgrp_brch=total change in vegve node number normalized for maturity group
!    TotReproNodeNumNormByMatrgrp_brch=total change in reprve node number normalized for maturity group
!    iPlantCalendar_brch(ipltcal_Emerge)=emergence date
!
  call DebugPrint('beg '//subname//' NZ',NZ)
  IF((EnablePlantLeafOut_brch(NB,NZ).EQ.iTrue .AND. &
    iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial) & !perrenial
    .AND. (Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ)))THEN
    ! branch is ready to do leaf out
    MatureGroup_brch(NB,NZ)                      = MatureGroup_pft(NZ)
    ShootNodeNumAtAnthesis_brch(NB,NZ)           = 0._r8
    LeafNumberAtFloralInit_brch(NB,NZ)           = 0._r8
    TotalNodeNumNormByMatgrp_brch(NB,NZ)         = 0._r8
    TotReproNodeNumNormByMatrgrp_brch(NB,NZ)     = 0._r8
    dReproNodeNumNormByMatG_brch(NB,NZ)          = 0._r8
    iPlantCalendar_brch(ipltcal_Emerge,NB,NZ)    = I
    iPlantCalendar_brch(2:NumGrowthStages,NB,NZ) = 0
    ! Start each perennial reproductive cycle with no failed grain-fill hours.
    HourFailGrainFill_brch(NB,NZ)               = 0._r8

    IF(NB.EQ.MainBranchNum_pft(NZ))THEN
      HoursTooLowPsiCan_pft(NZ)=0._r8
    ENDIF
    !
    !   SPRING LEAF AND SHEATH RESET
    !
    !   doInitLeafOut_brch,EnablePlantLeafOut_brch=flag for initializing,enabling leafout
    !   Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
    !   PSTG=node number
    !   NumOfLeaves_brch=number of leaves appeared
    !     KHiestGroLeafNode_brch=integer of most recent leaf number currently growing
    !     HourFailGrainFill_brch=number of hours with no grain fill
    !     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
    !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
    !     WT*,WT*N,WT*P=branch organ C,N,P mass
    !     WG,WG*N,WG*P=node organ C,N,P mass
    !     organ key:LF=leaf,SHE=PetolSheth,STK=stalk,RSV=reserve
    !     HSK=husk,EAR=ear,GR=grain,SHT=shoot
    !     LeafAreaLive_brch,LeafArea_node=branch,node leaf area
    !     StructInternodeElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
    !     StalkNodeVertLength_brch,StalkNodeHeight_brch=stalk height,stalk internode length
    !     SetNumberSeeds_brch=seed set number
    !     PotentialSeedSites_brch=potential number of seed set sites
    !     SingleGrainMeanBiomC_brch=individual seed size
    !
    IF(EnablePlantLeafOut_brch(NB,NZ).EQ.iTrue .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial &
      .AND. Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN
      IF(iPlantTurnoverPattern_pft(NZ).EQ.0)THEN
        ShootNodeNum_brch(NB,NZ)      = ShootNodeNumAtPlanting_pft(NZ)
        NumOfLeaves_brch(NB,NZ)       = 0._r8
        KLeafNumber_brch(NB,NZ)       = 1
        KHiestGroLeafNode_brch(NB,NZ) = 1
        D5330: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            dFall= PlantElmAllocMat4Litr(NE,icwood,M,NZ)*AZMAX1(LeafStrutElms_brch(NE,NB,NZ))*FracLeafShethElmAlloc2Litr(NE,k_woody_comp) &
              +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*AZMAX1(PetolShethStrutElms_brch(NE,NB,NZ))*FracPetolShethAlloc2Litr(NE,k_woody_comp)

            LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) +dFall
            LitrFallElms_brch(NE,NB,NZ)              = LitrFallElms_brch(NE,NB,NZ)+dFall

            dFall= PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)*AZMAX1(LeafStrutElms_brch(NE,NB,NZ))*FracLeafShethElmAlloc2Litr(NE,k_fine_comp) &
              +PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(PetolShethStrutElms_brch(NE,NB,NZ))*FracPetolShethAlloc2Litr(NE,k_fine_comp)
            LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) +dFall
            LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ)+dFall
          ENDDO  
        ENDDO D5330
        LeafAreaLive_brch(NB,NZ)                       = 0._r8
        LeafStrutElms_brch(1:NumPlantChemElms,NB,NZ)   = 0._r8
        PetolShethStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
        D5335: DO K=0,MaxNodesPerBranch1
          LeafArea_node(K,NB,NZ)                            = 0._r8
          PetoleLength_node(K,NB,NZ)                        = 0._r8
          LeafProteinC_node(K,NB,NZ)                        = 0._r8
          LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)    = 0._r8
          PetolShethElmntNode_brch(1:NumPlantChemElms,K,NB,NZ) = 0._r8
          PetoleProteinC_node(K,NB,NZ)                      = 0._r8
        ENDDO D5335
      ENDIF
    ENDIF
    !
    !     RESIDUAL STALKS BECOME LitrFall IN GRASSES, SHRUBS AT
    !     START OF SEASON
    !
    IF((EnablePlantLeafOut_brch(NB,NZ).EQ.iTrue.AND.iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial)&
      .AND.Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN

      D6232: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          dFall= PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(HuskStrutElms_brch(NE,NB,NZ) &
            +EarStrutElms_brch(NE,NB,NZ)+GrainStrutElms_brch(NE,NB,NZ))

          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ)+dFall

        ENDDO
      ENDDO D6232
        
      HuskStrutElms_brch(1:NumPlantChemElms,NB,NZ)  = 0._r8
      EarStrutElms_brch(1:NumPlantChemElms,NB,NZ)   = 0._r8
      GrainStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
      
      PotentialSeedSites_brch(NB,NZ) = 0._r8
      SetNumberSeeds_brch(NB,NZ)     = 0._r8
      SingleGrainMeanBiomC_brch(NB,NZ) = 0._r8
      IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. (.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ))))THEN
        D6345: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            dFall=PlantElmAllocMat4Litr(NE,istalk,M,NZ)*AZMAX1(StalkStrutElms_brch(NE,NB,NZ))
            LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+dFall
            LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ)+dFall
          ENDDO
        ENDDO D6345
        StalkStrutElms_brch(1:NumPlantChemElms,NB,NZ)      = 0._r8
        SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
        DO K=0,MaxNodesPerBranch1
          DO NE=1,NumPlantChemElms
            StructInternodeElms_brch(NE,K,NB,NZ)=0._r8
          ENDDO
        ENDDO
        D6340: DO K=0,MaxNodesPerBranch1
          StalkNodeHeight_brch(K,NB,NZ)     = 0._r8
          StalkNodeVertLength_brch(K,NB,NZ) = 0._r8
        ENDDO D6340
      ENDIF
    ENDIF
    ! Use the new season's node count after any turnover-specific reset.
    ShootNodeNumAtInitFloral_brch(NB,NZ) = ShootNodeNum_brch(NB,NZ)
  ENDIF
  !
  !   SPRING OR FALL FLAG RESET
  !
  !
  IF(EnablePlantLeafOut_brch(NB,NZ).EQ.iTrue .AND. Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN
    EnablePlantLeafOut_brch(NB,NZ)          = iFalse
    doPlantLeaveOff_brch(NB,NZ)         = iTrue
    Prep4Literfall_brch(NB,NZ)          = ifalse
    Hours4LiterfalAftMature_brch(NB,NZ) = 0
    !doing leave off    
  ELSE
    EnablePlantLeafOut_brch(NB,NZ)          = iTrue
    doPlantLeaveOff_brch(NB,NZ)         = iFalse
    Prep4Literfall_brch(NB,NZ)          = itrue
    Hours4LiterfalAftMature_brch(NB,NZ) = 0
    doInitLeafOut_brch(NB,NZ)           = iTrue
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine ResetNonAnnualBranch

!----------------------------------------------------------------------------------------------------
  subroutine ResetBranchPhenology(I,J,NB,NZ)

  implicit none
  integer, intent(in) :: I,J,nb,nz
  character(len=*), parameter :: subname='ResetBranchPhenology'
  integer :: K,M,NE
  real(r8) :: RSpecLiterFall1,dFall
  logical :: LeafOffTest,LeafOutTest
! begin_execution
  associate(                                                                 &
    istalk                       => pltpar%istalk                           ,& !input  :group id of plant stalk litter group
    k_fine_comp                  => pltpar%k_fine_comp                      ,& !input  :fine litter complex id
    inonfoliar                   => pltpar%inonfoliar                       ,& !input  :group id of plant non-foliar litter group
    HourReq4LeafOff_brch         => plt_pheno%HourReq4LeafOff_brch          ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    iPlantPhenolType_pft         => plt_pheno%iPlantPhenolType_pft          ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iPlantRootProfile_pft        => plt_pheno%iPlantRootProfile_pft         ,& !input  :plant growth type (vascular, non-vascular),[-]
    Hours4LeafOff_brch           => plt_pheno%Hours4LeafOff_brch            ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    HourReq4LeafOut_brch         => plt_pheno%HourReq4LeafOut_brch          ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    iPlantPhenolPattern_pft      => plt_pheno%iPlantPhenolPattern_pft       ,& !input  :plant growth habit: annual or perennial,[-]
    doPlantLeaveOff_brch         => plt_pheno%doPlantLeaveOff_brch          ,& !input  :branch phenology flag, [-]
    EnablePlantLeafOut_brch      => plt_pheno%EnablePlantLeafOut_brch       ,& !input  :branch phenology flag, [-]
    iPlantTurnoverPattern_pft    => plt_pheno%iPlantTurnoverPattern_pft     ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    Hours4Leafout_brch           => plt_pheno%Hours4Leafout_brch            ,& !input  :heat requirement for spring leafout/dehardening, [h]
    PlantElmAllocMat4Litr        => plt_soilchem%PlantElmAllocMat4Litr      ,& !input  :litter kinetic fraction, [-]
    iYearCurrent                 => plt_site%iYearCurrent                   ,& !input  :current year,[-]
    iPlant2ndGrothPattern_pft    => plt_pheno%iPlant2ndGrothPattern_pft     ,& !input  :plant expression of secondary growth, [-]                
    MainBranchNum_pft            => plt_morph%MainBranchNum_pft             ,& !input  :id number of main branch,[-]
    GrainStrutElms_brch          => plt_biom%GrainStrutElms_brch            ,& !inoput :branch grain structural element mass, [g d-2]
    SenecStalkStrutElms_brch     => plt_biom%SenecStalkStrutElms_brch       ,& !inoput :branch stalk structural element, [g d-2]
    StalkStrutElms_brch          => plt_biom%StalkStrutElms_brch            ,& !inoput :branch stalk structural element mass, [g d-2]
    EarStrutElms_brch            => plt_biom%EarStrutElms_brch              ,& !inoput :branch ear structural chemical element mass, [g d-2]
    HuskStrutElms_brch           => plt_biom%HuskStrutElms_brch             ,& !inoput :branch husk structural element mass, [g d-2]
    SeasonalNonstElms_pft        => plt_biom%SeasonalNonstElms_pft          ,& !inoput :plant stored nonstructural element at current step, [g d-2]
    StructInternodeElms_brch     => plt_biom%StructInternodeElms_brch       ,& !inoput :internode chemical element, [g d-2]
    SingleGrainMeanBiomC_brch      => plt_allom%SingleGrainMeanBiomC_brch       ,& !inoput :potential carbon mass per grain, [gC seed-1]
    Prep4Literfall_brch          => plt_pheno%Prep4Literfall_brch           ,& !inoput :branch phenology flag, [-]
    Hours4LiterfalAftMature_brch => plt_pheno%Hours4LiterfalAftMature_brch  ,& !inoput :branch phenology flag, [h]
    LitrfallElms_pvr             => plt_bgcr%LitrfallElms_pvr               ,& !inoput :plant LitrFall element, [g d-2 h-1]
    SSXferElms_pft               => plt_bgcr%SSXferElms_pft                 ,& !inoput :export flux from the seasonal storage, [g h-1 d-2]    
    SSXfer2ShootElms_pft         => plt_bgcr%SSXfer2ShootElms_pft           ,& !inoput :flux export from seasonal storage to shoot, [g h-1 d-2]        
    LitrFallElms_brch            => plt_bgcr%LitrFallElms_brch              ,& !inoput :litterfall from the branch, [g d-2 h-1]    
    StalkNodeHeight_brch         => plt_morph%StalkNodeHeight_brch          ,& !inoput :cumulative node height, [m]
    StalkNodeVertLength_brch     => plt_morph%StalkNodeVertLength_brch      ,& !inoput :internode height, [m]
    PotentialSeedSites_brch      => plt_morph%PotentialSeedSites_brch       ,& !inoput :branch potential grain number, [d-2]
    SetNumberSeeds_brch          => plt_morph%SetNumberSeeds_brch           ,& !inoput :branch grain number, [d-2]
    FracBiomHarvsted             => plt_distb%FracBiomHarvsted              ,& !output :harvest efficiency, [-]
    iYearPlantHarvest_pft        => plt_distb%iYearPlantHarvest_pft         ,& !output :year of harvest,[-]
    THIN_pft                     => plt_distb%THIN_pft                      ,& !output :thinning of plant population, [-]
    CanopyCutProxy_pft           => plt_distb%CanopyCutProxy_pft            ,& !output :harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
    iHarvstType_pft              => plt_distb%iHarvstType_pft               ,& !output :type of harvest,[-]
    jHarvstType_pft              => plt_distb%jHarvstType_pft               ,& !output :flag for stand replacing disturbance,[-]
    iDayPlanting_pft             => plt_distb%iDayPlanting_pft              ,& !output :day of planting,[-]
    iDayPlantHarvest_pft         => plt_distb%iDayPlantHarvest_pft          ,& !output :day of harvest,[-]
    iYearPlanting_pft            => plt_distb%iYearPlanting_pft             ,& !output :year of planting,[-]
    doInitPlant_pft              => plt_pheno%doInitPlant_pft                & !output :PFT initialization flag:0=no,1=yes,[-]
  )
  call DebugPrint('beg '//subname//' NZ',NZ)
  !   RESET PHENOLOGY AT EMERGENCE ('Hours4Leafout_brch' > 'VRNL')
  !   AND END OF SEASON ('Hours4LeafOff_brch' > 'VRNX')
  !
  !   iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
  !   iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
  !   EnablePlantLeafOut_brch=flag for enabling leafout:0=enable,1=disable
  !   Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
  !   doPlantLeaveOff_brch=flag for enabling leafoff:0=enable,1=disable
  !   Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
  !
  IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .OR. &                                    !perrenial
    (iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. is_drought_deciduos(iPlantPhenolType_pft(NZ))))THEN 

    LeafOutTest=(EnablePlantLeafOut_brch(NB,NZ).EQ.iTrue .AND. Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))
    LeafOffTest=(doPlantLeaveOff_brch(NB,NZ).EQ.iTrue .AND. Hours4LeafOff_brch(NB,NZ).GE.HourReq4LeafOff_brch(NB,NZ))
    IF(LeafOutTest .OR. LeafOffTest)THEN     
      !may include annual crop, which is set to evergreen
      CALL ResetNonAnnualBranch(I,J,NB,NZ)      
    ENDIF
  ENDIF
  !
  !   REPRODUCTIVE MATERIAL BECOMES LitrFall AT END OF SEASON
  !
  IF(Prep4Literfall_brch(NB,NZ).EQ.itrue)THEN
    Hours4LiterfalAftMature_brch(NB,NZ)=Hours4LiterfalAftMature_brch(NB,NZ)+1
    IF(Hours4LiterfalAftMature_brch(NB,NZ).EQ.HoursReq4LiterfalAftMature)THEN
      Prep4Literfall_brch(NB,NZ)          = ifalse
      Hours4LiterfalAftMature_brch(NB,NZ) = 0
    ENDIF
    RSpecLiterFall1=1.0_r8-RSpecLiterFall

    D6330: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        dFall=RSpecLiterFall*PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ))
        LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+ dFall
        LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ)+dFall

        IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantPhenolType_pft(NZ).NE.iphenotyp_evgreen)THEN
          dFall=RSpecLiterFall*PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(GrainStrutElms_brch(NE,NB,NZ))
          SeasonalNonstElms_pft(NE,NZ) = SeasonalNonstElms_pft(NE,NZ)+dFall
          SSXferElms_pft(NE,NZ)        = SSXferElms_pft(NE,NZ)+ dFall
          SSXfer2ShootElms_pft(NE,NZ)  = SSXfer2ShootElms_pft(NE,NZ)+dFall
        ELSE
          dFall=RSpecLiterFall*PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(GrainStrutElms_brch(NE,NB,NZ))        
          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ)+dFall
        ENDIF
      ENDDO
    ENDDO D6330

    DO NE=1,NumPlantChemElms  
      HuskStrutElms_brch(NE,NB,NZ)  = RSpecLiterFall1*AZMAX1(HuskStrutElms_brch(NE,NB,NZ))
      EarStrutElms_brch(NE,NB,NZ)   = RSpecLiterFall1*AZMAX1(EarStrutElms_brch(NE,NB,NZ))
      GrainStrutElms_brch(NE,NB,NZ) = RSpecLiterFall1*AZMAX1(GrainStrutElms_brch(NE,NB,NZ))
    ENDDO
    PotentialSeedSites_brch(NB,NZ) = RSpecLiterFall1*PotentialSeedSites_brch(NB,NZ)
    SetNumberSeeds_brch(NB,NZ)     = RSpecLiterFall1*SetNumberSeeds_brch(NB,NZ)
    ! Partial removal changes grain number, not the size potential of surviving grains.
    IF(RSpecLiterFall1.LE.0._r8) SingleGrainMeanBiomC_brch(NB,NZ) = 0._r8
    !
    !     STALKS BECOME LitrFall IN GRASSES AT END OF SEASON
    !
    IF((iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. (.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))) &
      .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial)THEN

      D6335: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          dFall=RSpecLiterFall*PlantElmAllocMat4Litr(NE,istalk,M,NZ)*AZMAX1(StalkStrutElms_brch(NE,NB,NZ))
          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ)+dFall
        ENDDO
      ENDDO D6335
      DO NE=1,NumPlantChemElms  
        StalkStrutElms_brch(NE,NB,NZ)      = RSpecLiterFall1*AZMAX1(StalkStrutElms_brch(NE,NB,NZ))
        SenecStalkStrutElms_brch(NE,NB,NZ) = RSpecLiterFall1*SenecStalkStrutElms_brch(NE,NB,NZ)
      ENDDO
      DO K=0,MaxNodesPerBranch1
        DO NE=1,NumPlantChemElms
          StructInternodeElms_brch(NE,K,NB,NZ)=RSpecLiterFall1*StructInternodeElms_brch(NE,K,NB,NZ)
        ENDDO
      ENDDO
      D2010: DO K=0,MaxNodesPerBranch1  
        StalkNodeVertLength_brch(K,NB,NZ)=RSpecLiterFall1*StalkNodeVertLength_brch(K,NB,NZ)
        StalkNodeHeight_brch(K,NB,NZ)=RSpecLiterFall1*StalkNodeHeight_brch(K,NB,NZ)
      ENDDO D2010
    ENDIF

    !
    !   SELF-SEEDING ANNUALS IF COLD OR DROUGHT DECIDUOUS
    !
    !     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial
    !     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
    !     iDayPlantHarvest_pft,iYearPlantHarvest_pft=day,year of harvesting
    !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
    !                       ,3=pruning,4=grazing,5=fire,6=herbivory
    !     jHarvstType_pft=terminate PFT:0=no,1=yes,2=yes,but reseed
    !     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
    !          iHarvstType_pft=3:reduction of clumping factor
    !          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
    !     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
    !          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
    !     FracBiomHarvsted(1,1:4)=fraction of leaf,non-foliar,woody, standing dead removed from PFT
    !     FracBiomHarvsted(2,1:4)=fraction of leaf,non-foliar,woody, standing dead removed from ecosystem
    !     iDayPlanting_pft,iYearPlanting_pft=day,year of planting
    !     doInitPlant_pft=PFT initialization flag:0=no,1=yes
    !
    IF(NB.EQ.MainBranchNum_pft(NZ))THEN
      !deciduous annual plant      

      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantPhenolType_pft(NZ).NE.iphenotyp_evgreen)THEN
        iDayPlantHarvest_pft(NZ)                              = I
        iYearPlantHarvest_pft(NZ)                             = iYearCurrent
        iHarvstType_pft(NZ)                                   = iharvtyp_grain
        jHarvstType_pft(NZ)                                   = jharvtyp_tmareseed !terminate and reseed for annual plants
        CanopyCutProxy_pft(NZ)                                = 0._r8
        THIN_pft(NZ)                                          = 0._r8   !do not 
        FracBiomHarvsted(iHarvst_pft,:,NZ)                    = 1.0_r8  !remove everything from the aboveground plant
        FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ)        = 0._r8   !drop leaf as litter
        FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ) = 1.0_r8  !grain is removed
        FracBiomHarvsted(iHarvst_col,iplthvst_stalk,NZ)       = 0._r8   !drop stalk as litter
        FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ)      = 0._r8   !drop standing dead as litter
        iDayPlanting_pft(NZ)                                  = -1E+06
        iYearPlanting_pft(NZ)                                 = -1E+06
        doInitPlant_pft(NZ)                                   = itrue
        Hours4Leafout_brch(NB,NZ)                             = 0
        Hours4LeafOff_brch(NB,NZ)                             = 0
      ENDIF
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine ResetBranchPhenology

!----------------------------------------------------------------------------------------------------
  subroutine BranchElmntTransfer(I,J,NB,NZ,BegRemoblize,WaterStress4Groth,TurgEff4CanopyResp)
  use PlantNonstElmDynMod
  implicit none
  integer, intent(in) :: I,J,NB,NZ,BegRemoblize
  real(r8), intent(in) :: WaterStress4Groth
  real(r8), intent(in) :: TurgEff4CanopyResp
  integer :: NE
  real(r8) :: XFRE(1:NumPlantChemElms)
  logical :: AnnualPlantChk,PerennialPlantChk  
  logical :: PlantingChk,RemobChk,LeafOutChk,LeafMobil4OutChk
  character(len=*), parameter :: subname='BranchElmntTransfer'

  ! begin_execution
  associate(                                                       &
    iDayPlanting_pft        => plt_distb%iDayPlanting_pft         ,& !input  :day of planting,[-]
    iYearPlanting_pft       => plt_distb%iYearPlanting_pft        ,& !input  :year of planting,[-]
    iYearCurrent            => plt_site%iYearCurrent              ,& !input  :current year,[-]
    iPlantCalendar_brch     => plt_pheno%iPlantCalendar_brch      ,& !input  :plant growth stage, [-]
    fTCanopyGroth_pft       => plt_pheno%fTCanopyGroth_pft        ,& !input  :canopy temperature growth function, [-]
    HourReq4LeafOff_brch    => plt_pheno%HourReq4LeafOff_brch     ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    Hours4Leafout_brch      => plt_pheno%Hours4Leafout_brch       ,& !input  :heat requirement for spring leafout/dehardening, [h]
    HourReq4LeafOut_brch    => plt_pheno%HourReq4LeafOut_brch     ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    doInitPlant_pft         => plt_pheno%doInitPlant_pft          ,& !input  :PFT initialization flag:0=no,1=yes,[-]
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft  ,& !input  :plant growth habit: annual or perennial,[-]
    Hours4LeafOff_brch      => plt_pheno%Hours4LeafOff_brch       ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    iPlantPhenolType_pft    => plt_pheno%iPlantPhenolType_pft     ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    MainBranchNum_pft       => plt_morph%MainBranchNum_pft        ,& !input  :number of main branch,[-]
    CanopyNonstElms_brch    => plt_biom%CanopyNonstElms_brch      ,& !inoput :branch nonstructural element, [g d-2]
    doInitLeafOut_brch      => plt_pheno%doInitLeafOut_brch       ,& !inoput :branch phenology flag, [-]
    Hours2LeafOut_brch      => plt_pheno%Hours2LeafOut_brch        & !inoput :counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
  )
  call DebugPrint('beg '//subname//' NZ',NZ)
!   TRANSFER C,N,P FROM SEASONAL STORAGE TO SHOOT AND ROOT
!   NON-STRUCTURAL C DURING SEED GERMINATION OR LEAFOUT
!
  PlantingChk    = I.GE.iDayPlanting_pft(NZ) .AND. iYearCurrent.EQ.iYearPlanting_pft(NZ)
  RemobChk       = Hours4LeafOff_brch(NB,NZ).LT.FracHour4LeafoffRemob(iPlantPhenolType_pft(NZ))*HourReq4LeafOff_brch(NB,NZ)
  LeafOutChk     = Hours4Leafout_brch(MainBranchNum_pft(NZ),NZ).GE.HourReq4LeafOut_brch(NB,NZ)
  AnnualPlantChk = iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. doInitPlant_pft(NZ).EQ.ifalse

  IF(AnnualPlantChk .OR. (PlantingChk.AND.RemobChk) .OR. (LeafOutChk.AND.RemobChk))THEN
    !
    ! RESET TIME COUNTER
    !
    ! Hours2LeafOut_brch=hourly leafout counter
    ! doInitLeafOut_brch=flag for initializing leafout
    !
    IF(doInitLeafOut_brch(NB,NZ).EQ.iTrue)THEN
      Hours2LeafOut_brch(NB,NZ) = 0._r8
      doInitLeafOut_brch(NB,NZ) = iFalse
    ENDIF
    !
    ! INCREMENT TIME COUNTER
    !
    ! iPlantPhotoperiodType_pft=photoperiod type:0=day neutral,1=short day,2=long day
    ! iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
    ! CriticPhotoPeriod_pft=critical photoperiod (h):<0=maximum daylength from site file
    ! PhotoPeriodSens_pft=photoperiod sensitivity (node h-1)
    ! DayLenthCurrent=daylength
    ! TurgEff4CanopyResp=expansion,extension function of canopy water potential
    ! fTCanopyGroth_pft=temperature function for canopy growth
    ! HourReq2InitSStor4LeafOut=number of hours required to initiate remobilization of storage C for leafout
    ! main branch leaf out
    
    IF(NB.EQ.MainBranchNum_pft(NZ))THEN      
      call DevelopMainBranch(I,J,NB,NZ,TurgEff4CanopyResp)
    ENDIF
    !
    ! REDISTRIBUTE TRANFERRED C FROM MAIN STEM TO OTHER BRANCHES
    !
    ! ATRP=hourly leafout counter
    ! fTCanopyGroth_pft=temperature function for canopy growth
    ! HourReq2InitSStor4LeafOut=number of hours required for remobilization of storage C during leafout
    ! WaterStress4Groth=growth function of canopy water potential
    ! 
    LeafMobil4OutChk=Hours2LeafOut_brch(NB,NZ).LE.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ))
    !not main branch
    IF(NB.NE.MainBranchNum_pft(NZ) .AND. LeafMobil4OutChk)THEN
      Hours2LeafOut_brch(NB,NZ)=Hours2LeafOut_brch(NB,NZ)+fTCanopyGroth_pft(NZ)*WaterStress4Groth
      !
      DO NE=1,NumPlantChemElms
        XFRE(NE)=AZMAX1(0.05_r8*fTCanopyGroth_pft(NZ) &
          *(0.5_r8*(CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)+CanopyNonstElms_brch(NE,NB,NZ)) &
          -CanopyNonstElms_brch(NE,NB,NZ)))
        call ExchFluxLimiter(CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ),CanopyNonstElms_brch(NE,NB,NZ),XFRE(NE))  
        CanopyNonstElms_brch(NE,NB,NZ)                    = CanopyNonstElms_brch(NE,NB,NZ)+XFRE(NE)
        CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ) = CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)-XFRE(NE)
      ENDDO
    ENDIF
  ENDIF
  !
  ! TRANSFER LEAF AND STALK NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
  ! IN PERENNIALS AFTER GRAIN FILL IN DETERMINATES, AFTER AUTUMNIZ'N
  ! IN INDETERMINATES, OR AFTER SUSTAINED WATER STRESS
  !
  IF(BegRemoblize.EQ.itrue .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial)THEN
    call SeasonStoreShootTransfer(I,J,NB,NZ)
  ENDIF
  !
  !   TRANSFER NON-STRUCTURAL C,N,P FROM LEAVES AND ROOTS TO RESERVES
  !   IN STALKS DURING GRAIN FILL IN ANNUALS OR BETWEEN STALK RESERVES
  !   AND LEAVES IN PERENNIALS ACCORDING TO CONCENTRATION DIFFERENCES
  !
  !   iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
  !   iPlantCalendar_brch(ipltcal_Jointing,=start of stem elongation and setting max seed number
  !   iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date setting for final seed number
  !   CanopyLeafSheathC_brch=leaf+PetolSheth mass
  !   SapwoodBiomassC_brch=stalk sapwood mass
  !
  AnnualPlantChk    = iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).NE.0
  PerennialPlantChk = iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial .AND. iPlantCalendar_brch(ipltcal_Jointing,NB,NZ).NE.0
  IF(AnnualPlantChk .OR. PerennialPlantChk)THEN
    !
    call StalkRsrvShootNonstTransfer(I,J,NB,NZ)
    !
    IF(AnnualPlantChk)THEN
      !stalk-root transfer
      call StalkRsrvRootNonstTransfer(I,J,NB,NZ)
    ENDIF
  ENDIF
  !
  !   REPLENISH BRANCH NON-STRUCTURAL POOL FROM
  !   SEASONAL STORAGE POOL 
  !
  call RepleteSeaStoreByStalk(I,J,NB,NZ)

  !
  call PrintInfo('end '//subname)
  end associate
  end subroutine BranchElmntTransfer

!----------------------------------------------------------------------------------------------------
  subroutine DevelopMainBranch(I,J,NB,NZ,TurgEff4CanopyResp)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: TurgEff4CanopyResp
  real(r8)  :: TotPopuPlantRootC  
  real(r8) :: CH2OH
  real(r8) :: FXFC,FXFN  
  REAL(R8) :: cpoolt  
  real(r8) :: ATRPPD
  real(r8) :: NonstGradt  
  real(r8) :: PPDX,WFNSP
  real(r8) :: GFNX  
  real(r8) :: NonstElm2RootMyco(NumPlantChemElms)  
  real(r8) :: TotalRootNonstElms(1:NumPlantChemElms)  
  real(r8) :: ElmXferStore2Shoot(1:NumPlantChemElms)  
  real(r8) :: dTimeIncLeafout,CNonst2Root
  integer :: L,NE
  logical ::  PlantChk,checkDeciduos
  character(len=*), parameter :: subname='DevelopMainBranch'

  associate(                                                           &
    PhotoPeriodSens_pft       => plt_pheno%PhotoPeriodSens_pft        ,& !input  :difference between current and critical daylengths used to calculate phenological progress, [h]
    NU                        => plt_site%NU                          ,& !input  :current soil surface layer number, [-]
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft              ,& !input  :threshold zero for plang growth calculation, [-]
    CriticPhotoPeriod_pft     => plt_pheno%CriticPhotoPeriod_pft      ,& !input  :critical daylength for phenological progress, [h]
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft      ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft       ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr          ,& !input  :root layer C, [gC d-2]
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft         ,& !input  :soil layer at planting depth, [-]
    fTCanopyGroth_pft         => plt_pheno%fTCanopyGroth_pft          ,& !input  :canopy temperature growth function, [-]
    iPlantPhenolPattern_pft   => plt_pheno%iPlantPhenolPattern_pft    ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft  ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    iPlantPhotoperiodType_pft => plt_pheno%iPlantPhotoperiodType_pft  ,& !input  :photoperiod type (neutral, long day, short day),[-]
    MaxSoilLays4Root_pft      => plt_morph%MaxSoilLays4Root_pft       ,& !input  :maximum soil layer number for all root axes,[-]
    DayLenthCurrent           => plt_site%DayLenthCurrent             ,& !input  :current daylength of the grid, [h]
    Xfer2RootsC_pft           => plt_bgcr%Xfer2RootsC_pft             ,& !inoput :carbon transfer from other nonstructural source to roots, [gC d-2 h-1]
    SSXferElms_pft            => plt_bgcr%SSXferElms_pft              ,& !inoput :export flux from the seasonal storage, [g h-1 d-2]
    SSXfer2ShootElms_pft      => plt_bgcr%SSXfer2ShootElms_pft        ,& !inoput :flux export from seasonal storage to shoot, [g h-1 d-2]    
    Hours2LeafOut_brch        => plt_pheno%Hours2LeafOut_brch         ,& !inoput :counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr      ,& !inoput :root layer nonstructural element, [g d-2]
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft       ,& !inoput :plant stored nonstructural element at current step, [g d-2]
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch         & !inoput :branch nonstructural element, [g d-2]
  )
  call DebugPrint('beg '//subname//' NZ',NZ)
  !write(1114,*)I*1000+J/24.,1,2,0.0,plt_biom%RootMycoNonstElms_rpvr(ielmc,1,2,NZ),'bf'//subname  
  TotPopuPlantRootC         = 0._r8
  TotalRootNonstElms(ielmc) = 0._r8
  D4: DO L=NU,MaxSoilLays4Root_pft(NZ)
    TotPopuPlantRootC         = TotPopuPlantRootC+AZMAX1(PopuRootMycoC_pvr(ipltroot,L,NZ))
    TotalRootNonstElms(ielmc) = TotalRootNonstElms(ielmc)+AZMAX1(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))
  ENDDO D4
  !
  IF(iPlantPhotoperiodType_pft(NZ).EQ.iphotop_long .AND. is_drought_deciduos(iPlantPhenolType_pft(NZ)))THEN
    PPDX   = AZMAX1(CriticPhotoPeriod_pft(NZ)-DayLenthCurrent-PhotoPeriodSens_pft(NZ))
    ATRPPD = EXP(-0.0_r8*PPDX)
  ELSE
    ATRPPD=1.0_r8
  ENDIF
  !
  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ)))THEN
    WFNSP=1.0_r8
  ELSE
    WFNSP=TurgEff4CanopyResp    
  ENDIF

  dTimeIncLeafout           = ATRPPD*fTCanopyGroth_pft(NZ)*WFNSP
  Hours2LeafOut_brch(NB,NZ) = Hours2LeafOut_brch(NB,NZ)+dTimeIncLeafout
  PlantChk                  = iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen
  CH2OH                     = 0._r8

  IF(Hours2LeafOut_brch(NB,NZ).LE.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ)) .OR. plantChk)THEN
    IF(SeasonalNonstElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      !CPOOLT:=total root nonst + branch nonst
      CPOOLT=TotalRootNonstElms(ielmc)+CanopyNonstElms_brch(ielmc,NB,NZ)
      !
      !       REMOBILIZE C FROM SEASONAL STORAGE AT FIRST-ORDER RATE
      !       MODIFIED BY SOIL TEMPERATURE AT SEED DEPTH
      !
      !     GVMX=specific oxidation rate of storage C during leafout at 25 C
      !     WTRVC=storage C
      !     CH2OH=storage C oxidation rate during leafout
      !     CPOOL,CPOOLR=non-structural C mass in branch,root
      !     FXSH,FXRT=shoot-root partitioning of storage C during leafout
      !     WTRTD=root C mass
      !
      GFNX                              = GVMX(iPlantPhenolPattern_pft(NZ))*dTimeIncLeafout
      CH2OH                             = AZMAX1(GFNX*SeasonalNonstElms_pft(ielmc,NZ))
      SeasonalNonstElms_pft(ielmc,NZ)   = SeasonalNonstElms_pft(ielmc,NZ)-CH2OH
      SSXferElms_pft(ielmc,NZ)          = SSXferElms_pft(ielmc,NZ)-CH2OH
      SSXfer2ShootElms_pft(ielmc,NZ)    = SSXfer2ShootElms_pft(ielmc,NZ)-CH2OH*FX2SH(iPlantPhenolPattern_pft(NZ))
      CanopyNonstElms_brch(ielmc,NB,NZ) = CanopyNonstElms_brch(ielmc,NB,NZ)+CH2OH*FX2SH(iPlantPhenolPattern_pft(NZ))
      !
      ! if root condition met
      CNonst2Root=CH2OH*FX2RT(iPlantPhenolPattern_pft(NZ))
      IF(TotPopuPlantRootC.GT.ZERO4Groth_pft(NZ) .AND. TotalRootNonstElms(ielmc).GT.ZERO4Groth_pft(NZ))THEN
        D50: DO L=NU,MaxSoilLays4Root_pft(NZ)
          FXFC                                        = AZMAX1(PopuRootMycoC_pvr(ipltroot,L,NZ))/TotPopuPlantRootC
          RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)+FXFC*CNonst2Root
          Xfer2RootsC_pft(NZ) = Xfer2RootsC_pft(NZ) + FXFC*CNonst2Root
        ENDDO D50
      ELSE
          RootMycoNonstElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),NZ)= &
            RootMycoNonstElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),NZ)+CNonst2Root
          Xfer2RootsC_pft(NZ) = Xfer2RootsC_pft(NZ) + CNonst2Root
      ENDIF
    ENDIF
  ENDIF
  !
  !     REMOBILIZE N,P FROM SEASONAL STORAGE AT FIRST-ORDER RATE
  !     MODIFIED BY SOIL TEMPERATURE AT SEED DEPTH
  !
  !     WTRVC,WTRVN,WTRVP=storage C,N,P
  !     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
  !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
  !     NXferStore2Shoot,PXferStore2Shoot=N,P transfer from storage to shoot
  !     CH2OH=storage C oxidation rate during leafout
  !     FRSV=rate constant for remobiln of storage C,N,P during leafout C
  !     FXSH=shoot partitioning of storage C during leafout
  !
  IF(SeasonalNonstElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial)THEN
      !perennial plant
      CPOOLT=AZMAX1(SeasonalNonstElms_pft(ielmc,NZ)+CanopyNonstElms_brch(ielmc,NB,NZ))          
! --- Jing Tao ---
! DIVIDE-BY-ZERO, measured not inferred. The IF above guards
! SeasonalNonstElms_pft(ielmc), but the divisor is CPOOLT -- a DIFFERENT quantity. CPOOLT
! comes from AZMAX1, whose single-argument form (MiniMathMod::AZMAX1_s) returns EXACTLY
! 0._r8 when its argument is below tiny_val=1.e-12. So when both nonstructural C pools sit
! near 1e-13 their sum clamps to zero, the guard still passes, and this divides by zero.
! Evidence, EcoSIM_BioCON R3 case 7143 at simulated 2022, caught under gdb: MXCSR = 0x1926
! (ZE set, IE and OE clear), divisor xmm1 = exactly 0.0, dividend xmm0 = -2.955e-24 --
! which is the product of two ~1e-12 pools, corroborating this site rather than another.
! 42 of the 50 cases still failing after the R3 relaunch die this way.
!
! FIX: skip, do not floor. A zero CPOOLT means both pools are empty, so there is no
! gradient and no transfer; flooring the denominator with an epsilon would manufacture a
! large spurious flux out of nothing.
!
! NOT switch-gated, deliberately, and this is the reason: on every path where the old code
! RETURNED a value (CPOOLT > 0) the arithmetic is bit-for-bit unchanged, and the only path
! it alters is the one that previously aborted the run. There is no baseline behaviour to
! preserve in that branch, so a default-off switch would only preserve the crash.
! safe_adb() was considered and rejected: its success path calls AZERO(), which clamps any
! result below 1.e-10 to zero, so it would change values where the denominator is healthy
! and break V0-at-equality.
!
! The adjacent root-transfer block below already floors its CPOOLT with
! AMAX1(ZERO4Groth_pft,...). One of the two was protected and this one was not.
      DO NE=2,NumPlantChemElms
        if(CPOOLT>0._r8)then
          NonstGradt=(SeasonalNonstElms_pft(NE,NZ)*CanopyNonstElms_brch(ielmc,NB,NZ)- &
            CanopyNonstElms_brch(NE,NB,NZ)*SeasonalNonstElms_pft(ielmc,NZ))/CPOOLT            
        else
          NonstGradt=0._r8
        endif
        ElmXferStore2Shoot(NE)=AZMAX1(FRSV(iPlantTurnoverPattern_pft(NZ))*NonstGradt)
      ENDDO
      !annual plant
    ELSE
      DO NE=2,NumPlantChemElms
        ElmXferStore2Shoot(NE)=AZMAX1(FX2SH(iPlantPhenolPattern_pft(NZ))*CH2OH*SeasonalNonstElms_pft(NE,NZ) &
          /SeasonalNonstElms_pft(ielmc,NZ))
      ENDDO  
    ENDIF
  ELSE
    !when there is no seasonal storage C
    DO NE=2,NumPlantChemElms
      ElmXferStore2Shoot(NE)=AZMAX1(FX2SH(iPlantPhenolPattern_pft(NZ))*SeasonalNonstElms_pft(NE,NZ))
    ENDDO
  ENDIF
  !
  ! ADD TO NON-STRUCTURAL POOLS IN ROOT
  !
  ! CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
  ! WTRVC,WTRVN,WTRVP=storage C,N,P
  ! iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
  ! UPNH4R,UPPO4R=N,P transfer from storage to root
  ! FRSV=rate constant for remobiln of storage C,N,P during leafout
  ! FXRT=root partitioning of storage C during leafout
  !

  DO NE=1,NumPlantChemElms
    TotalRootNonstElms(NE)=0._r8
    D3: DO L=NU,MaxSoilLays4Root_pft(NZ)
      TotalRootNonstElms(NE)=TotalRootNonstElms(NE)+AZMAX1(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ))
    ENDDO D3
  ENDDO

  IF(SeasonalNonstElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial)THEN
      CPOOLT=AMAX1(ZERO4Groth_pft(NZ),SeasonalNonstElms_pft(ielmc,NZ)+TotalRootNonstElms(ielmc))
      DO NE=2,NumPlantChemElms
        NonstGradt=(SeasonalNonstElms_pft(NE,NZ)*TotalRootNonstElms(ielmc)-&
          TotalRootNonstElms(NE)*SeasonalNonstElms_pft(ielmc,NZ))/CPOOLT            
        NonstElm2RootMyco(NE)=AZMAX1(FRSV(iPlantTurnoverPattern_pft(NZ))*NonstGradt)
      ENDDO
    ELSE
      DO NE=2,NumPlantChemElms
        NonstElm2RootMyco(NE)=AZMAX1(FX2RT(iPlantPhenolPattern_pft(NZ))*CH2OH*&
          SeasonalNonstElms_pft(NE,NZ)/SeasonalNonstElms_pft(ielmc,NZ))
      ENDDO    
    ENDIF
  ELSE
    DO NE=2,NumPlantChemElms
      NonstElm2RootMyco(NE)=AZMAX1(FX2RT(iPlantPhenolPattern_pft(NZ))*SeasonalNonstElms_pft(NE,NZ))
    ENDDO  
  ENDIF
  !
  !     TRANSFER STORAGE FLUXES
  !
  !     WTRVC,WTRVN,WTRVP=storage C,N,P
  !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
  !     NXferStore2Shoot,PXferStore2Shoot=N,P transfer from storage to shoot
  !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
  !     UPNH4R,UPPO4R=N,P transfer from storage to root
  !     FXFN=root layer allocation
  !
  DO NE=2,NumPlantChemElms
    SeasonalNonstElms_pft(NE,NZ)   = SeasonalNonstElms_pft(NE,NZ)-ElmXferStore2Shoot(NE)-NonstElm2RootMyco(NE)
    CanopyNonstElms_brch(NE,NB,NZ) = CanopyNonstElms_brch(NE,NB,NZ)+ElmXferStore2Shoot(NE)
    SSXferElms_pft(NE,NZ)          = SSXferElms_pft(NE,NZ)+ElmXferStore2Shoot(NE)+NonstElm2RootMyco(NE)
    SSXfer2ShootElms_pft(NE,NZ)    = SSXfer2ShootElms_pft(NE,NZ)+ElmXferStore2Shoot(NE)
  ENDDO  
  !
  IF(TotPopuPlantRootC.GT.ZERO4Groth_pft(NZ) .AND. TotalRootNonstElms(ielmc).GT.ZERO4Groth_pft(NZ))THEN
    D51: DO L=NU,MaxSoilLays4Root_pft(NZ)
      FXFN=AZMAX1(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))/TotalRootNonstElms(ielmc)
      DO NE=2,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ) &
          +FXFN*NonstElm2RootMyco(NE)
      ENDDO  
    ENDDO D51
  ELSE
    DO NE=2,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,ipltroot,NGTopRootLayer_pft(NZ),NZ)=&
        RootMycoNonstElms_rpvr(NE,ipltroot,NGTopRootLayer_pft(NZ),NZ)+NonstElm2RootMyco(NE)   
    ENDDO            
  ENDIF
  !write(1114,*)I*1000+J/24.,1,2,0.0,plt_biom%RootMycoNonstElms_rpvr(ielmc,1,2,NZ),'af'//subname 
  call PrintInfo('end '//subname)
  end associate
  end subroutine DevelopMainBranch

!----------------------------------------------------------------------------------------------------
  subroutine ComputRAutoAfEmergence(I,J,NB,NZ,YCO2Gro_brch,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,CO2F,&
    CH2O,TFN5,WaterStress4Groth,TurgEff4CanopyResp,ShootStructN_brch,CanopyNonstElm4Gros,CNPG,RMxess_brch,&
    RNonstC4Groth_brch)

  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: YCO2Gro_brch
  real(r8), intent(in) :: CNLFM,CPLFM
  real(r8), intent(in) :: CNSHX,CPSHX
  real(r8), intent(in) :: CNLFX,CPLFX
  real(r8), intent(in) :: ShootStructN_brch
  real(r8), intent(in) :: CO2F      !CO2 fixation rate, [gC h-1 d-2]
  real(r8), intent(in) :: CH2O      !total CH2O production
  real(r8), intent(in) :: TFN5      !temperature function for canopy maintenance respiration
  real(r8), intent(in) :: WaterStress4Groth
  real(r8), intent(in) :: TurgEff4CanopyResp
  real(r8), intent(out) :: CanopyNonstElm4Gros(NumPlantChemElms)   !nonstructural C/N/P metabolized for canopy growth
  real(r8), INTENT(OUT) :: CNPG
  real(r8), intent(out) :: RMxess_brch          !excess maintenance respiration, drives remobilization & senescence
  real(r8), intent(out) :: RNonstC4Groth_brch   !nonstrucal C used for canopy growth [gC d-2 h-1]
  real(r8) :: RCO2NonstC_brch
  real(r8) :: RCO2Maint_brch
  real(r8) :: RCO2NonstC4Nassim_brch
  real(r8) :: ZPOOLB
  real(r8) :: PPOOLB
  real(r8) :: RCO2X,RCO2Y,RGFNP
  real(r8) :: RgroCO2_ltd  !Nutient limited growth respiration
  real(r8) :: Rauto_brch   !branch autotrophic respiraiton
  real(r8) :: CNG,CPG,RFN1,RFP1
  character(len=*), parameter :: subname='ComputRAutoAfEmergence'
! begin_execution
  associate(                                                          &
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch       ,& !input  :branch nonstructural element, [g d-2]
    LeafPetoNonstElmConc_brch => plt_biom%LeafPetoNonstElmConc_brch  ,& !input  :branch nonstructural C concentration, [g d-2]
    ZERO                      => plt_site%ZERO                       ,& !input  :threshold zero for numerical stability, [-]
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft     ,& !input  :plant growth type (vascular, non-vascular),[-]
    fTCanopyGroth_pft         => plt_pheno%fTCanopyGroth_pft         ,& !input  :canopy temperature growth function, [-]
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft      ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    GrainFillDowreg_brch      => plt_photo%GrainFillDowreg_brch      ,& !input  :down-regulation of C4 photosynthesis, [-]
    CO2NetFix_pft             => plt_bgcr%CO2NetFix_pft              ,& !inoput :canopy net CO2 exchange, [gC d-2 h-1]
    CanopyGrosRCO2_pft        => plt_bgcr%CanopyGrosRCO2_pft         ,& !inoput :canopy plant+nodule autotrophic respiraiton, [gC d-2]
    CanopyResp_brch           => plt_bgcr%CanopyResp_brch            ,& !inoput :canopy respiration for a branch, [gC d-2 h-1]
    Eco_AutoR_CumYr_col       => plt_bgcr%Eco_AutoR_CumYr_col        ,& !inoput :ecosystem autotrophic respiration, [g d-2 h-1]
    ECO_ER_col                => plt_bgcr%ECO_ER_col                 ,& !inoput :ecosystem respiration, [g d-2 h-1]
    GrossCO2Fix_pft           => plt_bgcr%GrossCO2Fix_pft            ,& !inoput :total gross CO2 fixation, [gC d-2 ]
    RCanMaintDef_CO2_pft      => plt_bgcr%RCanMaintDef_CO2_pft       ,& !inoput :canopy maintenance respiraiton deficit as CO2, [gC d-2 h-1]
    CanopyRespC_CumYr_pft     => plt_bgcr%CanopyRespC_CumYr_pft      ,& !inoput :total autotrophic respiration, [gC d-2 ]
    CanopyNLimFactor_brch     => plt_bgcr%CanopyNLimFactor_brch      ,& !output :Canopy N-limitation factor, [0->1] stronger limitation
    CanopyPLimFactor_brch     => plt_bgcr%CanopyPLimFactor_brch       & !output :Canopy P-limitation factor, [0->1] stronger limitation    
  )

  call DebugPrint('beg '//subname//' NZ',NZ)
  ! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
  !
  ! CNPG=N,P constraint on growth respiration
  ! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
  ! CNKI,CPKI=nonstruct N,P inhibn constant on growth(g N,P g-1 C)
  !
  IF(LeafPetoNonstElmConc_brch(ielmc,NB,NZ).GT.ZERO)THEN
    CNG=LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmn,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CNKI)
    CPG=LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmp,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CPKI)
    CNPG=AMIN1(CNG,CPG)
  ELSE
    CNG=1._r8
    CPG=1._r8
    CNPG=1.0_r8
  ENDIF
  CanopyNLimFactor_brch(NB,NZ)=CNG
  CanopyPLimFactor_brch(NB,NZ)=CPG

  !
  ! RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
  ! NON-STRUCTURAL C:N:P
  !
  ! RCO2NonstC_brch=respiration from non-structural C
  ! VMXC=rate constant for nonstructural C oxidation in respiration (h-1)
  ! CPOOL=non-structural C mass
  ! fTCanopyGroth_pft=temperature function for canopy growth
  ! WaterStress4Groth=growth function of canopy water potential
  ! CNPG=N,P constraint on respiration
  ! GrainFillDowreg_brch=termination feedback inhibition on C3 CO2
  !
  RCO2NonstC_brch=AZMAX1(VMXC*CanopyNonstElms_brch(ielmc,NB,NZ) &
    *fTCanopyGroth_pft(NZ))*CNPG*GrainFillDowreg_brch(NB,NZ)*WaterStress4Groth

  !
  ! MAINTENANCE RESPIRATION FROM TEMPERATURE, PLANT STRUCTURAL N
  !
  ! RCO2Maint_brch=maintenance respiration
  ! TFN5=
  ! ShootStructE_brch(ielmn)=shoot structural N mass
  ! iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
  ! iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
  ! WaterStress4Groth=growth function of canopy water potential
  !
  !maintenance respiraiton is computed based on whole branch
  RCO2Maint_brch=AZMAX1(RmSpecPlant*TFN5*ShootStructN_brch)     
  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ)) .OR. &
    iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    RCO2Maint_brch=RCO2Maint_brch*WaterStress4Groth
  ENDIF
  
!
! GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
! IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
! RCO2X=difference between non-structural C respn and mntc respn
! RCO2Y=growth respiration unlimited by N,P
! TurgEff4CanopyResp=expansion,extension function of canopy water potential

!
  RCO2X       = RCO2NonstC_brch-RCO2Maint_brch
  RCO2Y       = AZMAX1(RCO2X)*TurgEff4CanopyResp  
  RMxess_brch = AZMAX1(-RCO2X)

!
! GROWTH RESPIRATION MAY BE LIMITED BY NON-STRUCTURAL N,P
! AVAILABLE FOR GROWTH
!
! RCO2Y,RgroCO2_ltd=growth respiration unlimited,limited by N,P
! CNLFX=diff between min and max leaf N prodn vs nonstruct C consumption
! CNSHX=N production vs nonstructural C consumption in rest of shoot
! ZPOOL,PPOOL=nonstructural N,P mass
! YCO2Gro_brch=branch C respiration vs nonstructural C consumption
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
!
  IF(RCO2Y.GT.0.0_r8 .AND. (CNSHX.GT.0.0_r8 .OR. CNLFX.GT.0.0_r8))THEN
    ZPOOLB      = AZMAX1(CanopyNonstElms_brch(ielmn,NB,NZ))
    PPOOLB      = AZMAX1(CanopyNonstElms_brch(ielmp,NB,NZ))
    RFN1        = ZPOOLB*YCO2Gro_brch/(CNSHX+CNLFM+CNLFX*CNPG)
    RFP1        = PPOOLB*YCO2Gro_brch/(CPSHX+CPLFM+CPLFX*CNPG)
    RGFNP       = AMIN1(RFN1,RFP1)
    RgroCO2_ltd = AMIN1(RCO2Y,RGFNP)
  ELSE
    RgroCO2_ltd=0._r8
  ENDIF

!
! TOTAL NON-STRUCTURAL C,N,P USED IN GROWTH
! AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELDS
! ENTERED IN 'READQ'
!
! RNonstC4Groth_brch=total non-structural C used in growth and growth respiration
! RgroCO2_ltd=growth respiration limited by N,P
! YCO2Gro_brch=branch C respiration vs nonstructural C consumption, growth respiration quotient
! CanopyNonstElm4Gros(ielmn),CanopyNonstElm4Gros(ielmp)=nonstructural N,P used in growth
! ZPOOL,PPOOL=nonstructural N,P mass
! CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
! RCO2NonstC4Nassim_brch=respiration for N assimilation
! CH2O= 1/4 newly fixed C is used for N-assimilation (where is this number from?)
!
  
  RNonstC4Groth_brch         = RgroCO2_ltd/YCO2Gro_brch
  CanopyNonstElm4Gros(ielmn) = AZMAX1(AMIN1(CanopyNonstElms_brch(ielmn,NB,NZ),RNonstC4Groth_brch*(CNSHX+CNLFM+CNLFX*CNPG)))
  CanopyNonstElm4Gros(ielmp) = AZMAX1(AMIN1(CanopyNonstElms_brch(ielmp,NB,NZ),RNonstC4Groth_brch*(CPSHX+CPLFM+CPLFX*CNPG)))
  RCO2NonstC4Nassim_brch     = AZMAX1(1.70_r8*CanopyNonstElm4Gros(ielmn)-0.025_r8*CH2O)
  CanopyNonstElm4Gros(ielmc) = AMIN1(RCO2Maint_brch,RCO2NonstC_brch)+RNonstC4Groth_brch+RCO2NonstC4Nassim_brch
!
! TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
! ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
!
! Rauto_brch=total C respiration
! RCO2Maint_brch=maintenance respiration
! RCO2NonstC_brch=respiration from non-structural C
! RgroCO2_ltd=growth respiration limited by N,P
! RMxess_brch=excess maintenance respiration
! RCO2NonstC4Nassim_brch=respiration for N assimilation
! GrossCO2Fix_pft=total PFT CO2 fixation
! CO2F=total CO2 fixation
! CanopyGrosRCO2_pft,CanopyRespC_CumYr_pft=total,above-ground PFT respiration
! CO2NetFix_pft=PFT net CO2 fixation
! ECO_ER_col=ecosystem respiration
! Eco_AutoR_CumYr_col=total autotrophic respiration
!
  Rauto_brch                = AMIN1(RCO2Maint_brch,RCO2NonstC_brch)+RgroCO2_ltd+RCO2NonstC4Nassim_brch
  GrossCO2Fix_pft(NZ)       = GrossCO2Fix_pft(NZ)+CO2F
  CanopyGrosRCO2_pft(NZ)    = CanopyGrosRCO2_pft(NZ)-Rauto_brch
  CanopyResp_brch(NB,NZ)    = CanopyResp_brch(NB,NZ)-Rauto_brch
  CanopyRespC_CumYr_pft(NZ) = CanopyRespC_CumYr_pft(NZ)-Rauto_brch
  CO2NetFix_pft(NZ)         = CO2NetFix_pft(NZ)+CO2F-Rauto_brch
  ECO_ER_col                = ECO_ER_col-Rauto_brch
  Eco_AutoR_CumYr_col       = Eco_AutoR_CumYr_col-Rauto_brch
  RCanMaintDef_CO2_pft(NZ)  = RCanMaintDef_CO2_pft(NZ)-RMxess_brch
  
  call PrintInfo('end '//subname)
  end associate
  end subroutine ComputRAutoAfEmergence

!----------------------------------------------------------------------------------------------------
  subroutine ComputRAutoB4Emergence(I,J,NB,NZ,TFN6_vr,YCO2Gro_brch,CNLFM,CPLFM,CNSHX,CPSHX,&
    CNLFX,CPLFX,ShootStructN_brch,WaterStress4Groth,TurgEff4CanopyResp,CanopyNonstElm4Gros,CNPG,&
    RMxess_brch,RNonstC4Groth_brch)

  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8),intent(in) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: YCO2Gro_brch !respiration quotient
  real(r8), intent(in) :: CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX
  real(r8), intent(in) :: ShootStructN_brch
  real(r8), intent(in) :: WaterStress4Groth
  real(r8), intent(in) :: TurgEff4CanopyResp
  real(r8), intent(out) :: CanopyNonstElm4Gros(NumPlantChemElms)
  real(r8), INTENT(OUT) :: CNPG
  real(r8), intent(out) :: RMxess_brch
  real(r8), intent(out) :: RNonstC4Groth_brch    !nonstrucal C to drive structural biomass growth
  real(r8) :: RCO2NonstC_brch  
  real(r8) :: RCO2Maint_brch  
  real(r8) :: RCO2NonstC4Nassim_brch

  character(len=*), parameter :: subname='ComputRAutoB4Emergence'
  real(r8) :: ZPOOLB,ZADDBM,CGROSM
  real(r8) :: RCO2NonstC4NassimOUltd                !O2 unlimited respiration for N assimilation
  real(r8) :: FNP,RFN1,RFP1
  real(r8) :: PPOOLB
  real(r8) :: RCO2X_O2ld,RCO2Y
  real(r8) :: RgroCO2_ltd   !oxygen & nutrient-limited growth respiration
  real(r8) :: Rauto_brch,RCO2NonstC_O2ulm
  real(r8) :: RCO2X_O2ulm,RCO2YM
  real(r8) :: RgroCO2_ulm
  real(r8) :: RCO2TM  !total CO2 respiration unlimited by O2
  real(r8) :: RMxess_Ulm,CNG,CPG
  
! begin_execution
  associate(                                                          &
    LeafPetoNonstElmConc_brch => plt_biom%LeafPetoNonstElmConc_brch  ,& !input  :branch nonstructural C concentration, [g d-2]
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch       ,& !input  :branch nonstructural element, [g d-2]
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft     ,& !input  :plant growth type (vascular, non-vascular),[-]
    fTgrowRootP_vr            => plt_pheno%fTgrowRootP_vr            ,& !input  :root layer temperature growth functiom, [-]
    RAutoRootO2Limter_rpvr    => plt_rbgc%RAutoRootO2Limter_rpvr     ,& !input  :O2 constraint to root respiration (0-1), [-]
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft      ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    RootRespPotent_pvr        => plt_rbgc%RootRespPotent_pvr         ,& !input  :root respiration unconstrained by O2, [g d-2 h-1]
    RootCO2EmisPot_pvr        => plt_rbgc%RootCO2EmisPot_pvr         ,& !input  :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootCO2Autor_pvr          => plt_rbgc%RootCO2Autor_pvr           ,& !input  :root respiration constrained by O2, [g d-2 h-1]
    ZERO                      => plt_site%ZERO                       ,& !input  :threshold zero for numerical stability, [-]
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft        ,& !input  :soil layer at planting depth, [-]
    GrainFillDowreg_brch      => plt_photo%GrainFillDowreg_brch       & !input  :down-regulation of C4 photosynthesis, [-]
  )  
  !
  ! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
  !
  ! CNPG=N,P constraint on growth respiration
  ! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
  ! CNKI,CPKI=nonstructural N,P inhibition constant on growth
  !
  call DebugPrint('beg '//subname//' NZ',NZ)
  IF(LeafPetoNonstElmConc_brch(ielmc,NB,NZ).GT.ZERO)THEN
    CNG=LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmn,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CNKI)
    CPG=LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmp,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CPKI) 
    CNPG=AMIN1(CNG,CPG)
  ELSE
    CNPG=1.0_r8
  ENDIF
  !
  ! RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
  ! NON-STRUCTURAL C:N:P, O2 UPTAKE
  !
  ! RCO2NonstC_O2ulm,RCO2NonstC_brch=respiration from non-structural C unlimited,limited by O2
  ! VMXC=rate constant for nonstructural C oxidation in respiration (h-1)
  ! CPOOL=non-structural C mass
  ! fTgrowRootP_vr=temperature function for root growth
  ! WaterStress4Groth=growth function of canopy water potential
  ! CNPG=N,P constraint on respiration
  ! GrainFillDowreg_brch=termination feedback inhibition on C3 CO2
  ! RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
  !
  RCO2NonstC_O2ulm=AZMAX1(VMXC*CanopyNonstElms_brch(ielmc,NB,NZ) &
    *fTgrowRootP_vr(NGTopRootLayer_pft(NZ),NZ))*WaterStress4Groth*CNPG*GrainFillDowreg_brch(NB,NZ)

  RCO2NonstC_brch=RCO2NonstC_O2ulm*RAutoRootO2Limter_rpvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)  
  !
  ! MAINTENANCE RESPIRATION FROM TEMPERATURE, PLANT STRUCTURAL N
  !
  ! RCO2Maint_brch=maintenance respiration
  ! TFN6_vr=temperature function for root maintenance respiration
  ! ShootStructE_brch(ielmn)=shoot structural N mass
  ! iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
  ! iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
  ! WaterStress4Groth=growth function of canopy water potential
  !
  RCO2Maint_brch=AZMAX1(RmSpecPlant*TFN6_vr(NGTopRootLayer_pft(NZ))*ShootStructN_brch)
  !
  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ)) .OR. &
    iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    RCO2Maint_brch=RCO2Maint_brch*WaterStress4Groth
  ENDIF
  !
  ! GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
  ! IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
  !
  ! RCO2X_O2ulm,RCO2X_O2ld=diff between C respn unltd,ltd by O2 and mntc respn
  ! RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
  ! TurgEff4CanopyResp=expansion,extension function of canopy water potential
  ! RMxess_Ulm,RMxess_brch=excess maintenance respiration unltd,ltd by O2
  !
  RCO2X_O2ulm = RCO2NonstC_O2ulm-RCO2Maint_brch
  RCO2X_O2ld  = RCO2NonstC_brch-RCO2Maint_brch
  RCO2YM      = AZMAX1(RCO2X_O2ulm)*TurgEff4CanopyResp
  RCO2Y       = AZMAX1(RCO2X_O2ld)*TurgEff4CanopyResp
  RMxess_Ulm  = AZMAX1(-RCO2X_O2ulm)
  RMxess_brch = AZMAX1(-RCO2X_O2ld)

  !
  ! GROWTH RESPIRATION MAY BE LIMITED BY NON-STRUCTURAL N,P
  ! AVAILABLE FOR GROWTH
  !
  ! RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
  ! CNLFX=diff between min and max leaf N prodn vs nonstruct C consumption
  ! CNSHX=N production vs nonstructural C consumption in rest of shoot
  ! ZPOOL,PPOOL=nonstructural N,P mass
  ! FNP=growth respiration limited by O2 and N,P
  ! YCO2Gro_brch=branch C respiration vs nonstructural C consumption
  ! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
  ! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
  ! CNPG=N,P constraint on growth respiration
  ! RgroCO2_ulm,RgroCO2_ltd=growth respiration unltd,ltd by O2 and limited by N,P
  ! RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
  !
  FNP=1._r8
  IF(CNSHX.GT.0.0_r8.OR.CNLFX.GT.0.0_r8)THEN
    ZPOOLB = AZMAX1(CanopyNonstElms_brch(ielmn,NB,NZ))
    PPOOLB = AZMAX1(CanopyNonstElms_brch(ielmp,NB,NZ))
    RFN1   = ZPOOLB*YCO2Gro_brch/(CNSHX+CNLFM+CNLFX*CNPG)
    RFP1   = PPOOLB*YCO2Gro_brch/(CPSHX+CPLFM+CPLFX*CNPG)
    FNP    = AMIN1(RFN1,RFP1)
    IF(RCO2YM.GT.0.0_r8)THEN
      RgroCO2_ulm=AMIN1(RCO2YM,FNP)
    ELSE
      RgroCO2_ulm=0._r8
    ENDIF

    IF(RCO2Y.GT.0.0_r8)THEN
      RgroCO2_ltd=AMIN1(RCO2Y,FNP*RAutoRootO2Limter_rpvr(ipltroot,NGTopRootLayer_pft(NZ),NZ))
    ELSE
      RgroCO2_ltd=0._r8
    ENDIF
  ELSE
    RgroCO2_ulm=0._r8
    RgroCO2_ltd=0._r8
  ENDIF
  !
  ! TOTAL NON-STRUCTURAL C,N,P USED IN GROWTH
  ! AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELDS
  ! ENTERED IN 'READQ'
  !
  ! CGROSM,RNonstC4Groth_brch=total non-structural C used in growth and respn unltd,ltd by O2
  ! RgroCO2_ulm,RgroCO2_ltd=growth respiration unltd,ltd by O2 and limited by N,P
  ! YCO2Gro_brch=branch C respiration vs nonstructural C consumption
  ! ZADDBM,CanopyNonstElm4Gros(ielmn),CanopyNonstElm4Gros(ielmp)=nonstructural N,P unltd,ltd by O2 used in growth
  ! ZPOOL,PPOOL=nonstructural N,P mass
  ! CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
  ! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
  ! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
  ! CNPG=N,P constraint on growth respiration
  ! RCO2NonstC4NassimOUltd,RCO2NonstC4Nassim_brch=respiration for N assimilation unltd,ltd by O2
  !
  CGROSM                     = RgroCO2_ulm/YCO2Gro_brch
  RNonstC4Groth_brch         = RgroCO2_ltd/YCO2Gro_brch
  ZADDBM                     = AZMAX1(CGROSM*(CNSHX+CNLFM+CNLFX*CNPG))

  CanopyNonstElm4Gros(ielmn) = AZMAX1(RNonstC4Groth_brch*(CNSHX+CNLFM+CNLFX*CNPG))
  CanopyNonstElm4Gros(ielmp) = AZMAX1(RNonstC4Groth_brch*(CPSHX+CPLFM+CPLFX*CNPG))
  RCO2NonstC4NassimOUltd     = AZMAX1(1.70_r8*ZADDBM)
  RCO2NonstC4Nassim_brch     = AZMAX1(1.70_r8*CanopyNonstElm4Gros(ielmn))
  CanopyNonstElm4Gros(ielmc) = AMIN1(RCO2Maint_brch,RCO2NonstC_brch)+RNonstC4Groth_brch+RCO2NonstC4Nassim_brch
  !
  ! TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
  ! ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
  !
  ! RCO2TM,Rauto_brch=total C respiration unltd,ltd by O2
  ! RCO2Maint_brch=maintenance respiration
  ! RgroCO2_ulm,RgroCO2_ltd=growth respiration limited by N,P unltd,ltd by O2
  ! RMxess_Ulm,RMxess_brch=excess maintenance respiration unltd,ltd by O2
  ! RCO2NonstC4NassimOUltd,RCO2NonstC4Nassim_brch=respiration for N assimilation unltd,ltd by O2
  ! RootCO2Autor_pvr=total root respiration
  ! RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
  !
  RCO2TM     = AMIN1(RCO2Maint_brch,RCO2NonstC_O2ulm)+RgroCO2_ulm+RCO2NonstC4NassimOUltd
  Rauto_brch = AMIN1(RCO2Maint_brch,RCO2NonstC_brch)+RgroCO2_ltd+RCO2NonstC4Nassim_brch

  RootRespPotent_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ) = RootRespPotent_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+RCO2TM
  RootCO2EmisPot_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ) = RootCO2EmisPot_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+Rauto_brch
  RootCO2Autor_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)   = RootCO2Autor_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)-Rauto_brch

  call PrintInfo('end '//subname)
  end associate
  end subroutine ComputRAutoB4Emergence

!----------------------------------------------------------------------------------------------------
  subroutine GrowLeavesOnBranch(I,J,NZ,NB,MinNodeID,GrowthLeaf,EtiolationCoeff,TurgEff4LeafPetolExpansion,ALLOCL)
  implicit none
  integer , intent(in) :: I,J  
  integer , intent(in) :: NB,NZ
  integer , intent(in) :: MinNodeID
  real(r8), intent(in) :: GrowthLeaf(NumPlantChemElms)
  real(r8), intent(in) :: EtiolationCoeff
  real(r8), intent(in) :: TurgEff4LeafPetolExpansion
  REAL(R8), INTENT(OUT):: ALLOCL
  integer :: K,NE,KK
  integer :: MXNOD,MNNOD,KNOD
  real(r8) :: GNOD,dMC
  REAL(R8) :: GrozNodesPerPlant,LeafAreaGrowth,SpecAreaLeafGrowth
  REAL(R8) :: GrowthElms(NumPlantChemElms)
  character(len=*), parameter :: subname='GrowLeavesOnBranch'

  associate(                                                     &
    NumCogrowthNode_pft    => plt_morph%NumCogrowthNode_pft     ,& !input  :number of concurrently growing nodes,[-]
    SLA1_pft               => plt_morph%SLA1_pft                ,& !input  :leaf area:mass during growth, [m2 gC-1]
    KHiestGroLeafNode_brch => plt_pheno%KHiestGroLeafNode_brch  ,& !input  :leaf growth stage counter, [-]
    ZERO4LeafVar_pft       => plt_biom%ZERO4LeafVar_pft         ,& !input  :threshold zero for leaf calculation, [-]
    rProteinC2LeafN_pft    => plt_allom%rProteinC2LeafN_pft     ,& !input  :Plant leaf protein C to N ratio [g protein C (g leaf N)-1]
    rProteinC2LeafP_pft    => plt_allom%rProteinC2LeafP_pft     ,& !input  :Plant leaf protein C to P ratio [g protein C (g leaf P)-1]
    FracGroth2Node_pft     => plt_allom%FracGroth2Node_pft      ,& !input  :parameter for allocation of growth to nodes, [-]
    PlantPopuLive_pft      => plt_site%PlantPopuLive_pft        ,& !input  :plant population, [d-2]
    LeafArea_node          => plt_morph%LeafArea_node           ,& !inoput :leaf area, [m2 d-2]
    LeafAreaLive_brch      => plt_morph%LeafAreaLive_brch       ,& !inoput :branch leaf area, [m2 d-2]
    LeafElmntNode_brch     => plt_biom%LeafElmntNode_brch       ,& !inoput :leaf element, [g d-2]
    LeafProteinC_node      => plt_biom%LeafProteinC_node         & !inoput :leaf node protein C, [g d-2]
  )
  call DebugPrint('beg '//subname//' NZ',NZ)
  
  IF(GrowthLeaf(ielmc).GT.0.0_r8)THEN
    MXNOD  = KHiestGroLeafNode_brch(NB,NZ)
    MNNOD  = MAX(MinNodeID,MXNOD-NumCogrowthNode_pft(NZ)+1)
    MXNOD  = MAX(MXNOD,MNNOD)
    KNOD   = MXNOD-MNNOD+1
    GNOD   = KNOD           !number of growing nodes
    ALLOCL = 1.0_r8/GNOD    !allocation to each growing nodes
    DO NE  = 1, NumPlantChemElms
      GrowthElms(NE)=ALLOCL*GrowthLeaf(NE)
    ENDDO

    GrozNodesPerPlant=ALLOCL*FracGroth2Node_pft(NZ)*NumCogrowthNode_pft(NZ)
    !
    !     GROWTH AT EACH CURRENT NODE
    !
    !     WGLF,WGLFN,WGLFP,LeafProteinC_node=node leaf C,N,P,protein mass
    !     CNWS,rProteinC2LeafP_pft=protein:N,protein:P ratios from startq.f
    !
    D490: DO KK=MNNOD,MXNOD
      K=MOD(KK,MaxNodesPerBranch1)
      IF(K.EQ.0.AND.KK.NE.0)K=MaxNodesPerBranch1

      DO NE=1,NumPlantChemElms
        LeafElmntNode_brch(NE,K,NB,NZ)=LeafElmntNode_brch(NE,K,NB,NZ)+GrowthElms(NE)
      ENDDO
      LeafProteinC_node(K,NB,NZ)=LeafProteinC_node(K,NB,NZ)+&
        AMIN1(GrowthElms(ielmn)*rProteinC2LeafN_pft(NZ),GrowthElms(ielmp)*rProteinC2LeafP_pft(NZ))
      !
      !         SPECIFIC LEAF AREA FUNCTION OF CURRENT LEAF MASS
      !         AT EACH NODE
      !
      !         SpecAreaLeafGrowth=specific area of leaf growth
      !         ETOL=coefficient for etoliation effects on expansion,extension
      !         SLA1_pft=growth in leaf area vs mass from PFT file
      !         SLA2=parameter for calculating leaf area expansion
      !         WGLF=leaf C mass
      !         PP=PFT population
      !         GrozNodesPerPlant=number of growing nodes per plant
      !         TurgEff4LeafPetolExpansion=turgor expansion,extension function
      !         LeafAreaGrowth,GRO=leaf area,mass growth
      !         LeafAreaLive_brch,LeafArea_node=branch,node leaf area
      !Eq.(3) in Grant and Hesketh, 1992, with updates for small biomass
      dMC                      = AZMAX1(1.e-12_r8,LeafElmntNode_brch(ielmc,K,NB,NZ)/(PlantPopuLive_pft(NZ)*GrozNodesPerPlant))
      SpecAreaLeafGrowth       = EtiolationCoeff*SLA1_pft(NZ)*dMC**SLA2*TurgEff4LeafPetolExpansion
      LeafAreaGrowth           = GrowthElms(ielmc)*SpecAreaLeafGrowth
      LeafAreaLive_brch(NB,NZ) = LeafAreaLive_brch(NB,NZ)+LeafAreaGrowth
      LeafArea_node(K,NB,NZ)   = LeafArea_node(K,NB,NZ)+LeafAreaGrowth
    ENDDO D490
  ENDIF
  call PrintInfo('end '//subname)
  end associate  
  END subroutine GrowLeavesOnBranch

!----------------------------------------------------------------------------------------------------
  subroutine GrowPetolShethOnBranch(yearIJ,NZ,NB,MinNodeID,GrowthPetolSheth,EtiolationCoeff,TurgEff4LeafPetolExpansion,ALLOCL)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ,MinNodeID,NB
  real(r8), intent(in) :: GrowthPetolSheth(NumPlantChemElms)
  REAL(R8), INTENT(IN) :: EtiolationCoeff
  real(r8), intent(in) :: TurgEff4LeafPetolExpansion
  REAL(R8), INTENT(IN) :: ALLOCL
  integer :: MXNOD,MNNOD,NE,KK,K
  real(r8) :: GNOD,ALLOCS,GSSL,GROS
  REAL(R8) :: GrowthElms(NumPlantChemElms)
  REAL(R8) :: SSL,dMC
  character(len=*), parameter :: subname='GrowPetolShethOnBranch'

  associate(                                                        &
    NumCogrowthNode_pft       => plt_morph%NumCogrowthNode_pft     ,& !input  :number of concurrently growing nodes,[-]
    PetolShethLen2Mass_pft    => plt_morph%PetolShethLen2Mass_pft  ,& !input  :PetolSheth length:mass during growth, [m gC-1]
    SinePetolShethAngle_pft   => plt_morph%SinePetolShethAngle_pft ,& !input  :sheath angle, [degree from horizontal]
    LeafElmntNode_brch        => plt_biom%LeafElmntNode_brch       ,& !input  :leaf element, [g d-2]
    ZERO4LeafVar_pft          => plt_biom%ZERO4LeafVar_pft         ,& !input  :threshold zero for leaf calculation, [-]
    rProteinC2LeafN_pft       => plt_allom%rProteinC2LeafN_pft     ,& !input  :Plant leaf protein C to N ratio [g protein C (g leaf N)-1]
    rProteinC2LeafP_pft       => plt_allom%rProteinC2LeafP_pft     ,& !input  :Plant leaf protein C to P ratio [g protein C (g leaf P)-1]
    FracGroth2Node_pft        => plt_allom%FracGroth2Node_pft      ,& !input  :parameter for allocation of growth to nodes, [-]
    PlantPopuLive_pft         => plt_site%PlantPopuLive_pft        ,& !input  :plant population, [d-2]
    KHiestGroLeafNode_brch    => plt_pheno%KHiestGroLeafNode_brch  ,& !input  :leaf growth stage counter, [-]
    PetoleLength_node         => plt_morph%PetoleLength_node       ,& !inoput :sheath height, [m]
    PetolShethElmntNode_brch  => plt_biom%PetolShethElmntNode_brch ,& !inoput :sheath chemical element, [g d-2]
    PetoleProteinC_node       => plt_biom%PetoleProteinC_node       & !inoput :layer sheath protein C, [g d-2]
  )
  call DebugPrint('beg '//subname//' NZ',NZ)

  IF(GrowthPetolSheth(ielmc).GT.0.0_r8)THEN
    MXNOD  = KHiestGroLeafNode_brch(NB,NZ)
    MNNOD  = MAX(MinNodeID,MXNOD-NumCogrowthNode_pft(NZ)+1)
    MXNOD  = MAX(MXNOD,MNNOD)
    GNOD   = MXNOD-MNNOD+1
    ALLOCS = 1.0_r8/GNOD
    DO NE=1,NumPlantChemElms
      GrowthElms(NE)=ALLOCS*GrowthPetolSheth(NE)
    ENDDO
    GSSL=ALLOCL*FracGroth2Node_pft(NZ)*NumCogrowthNode_pft(NZ)
    !
    !       GROWTH AT EACH CURRENT NODE
    !
    !       PetolShethElmntNode_brch,WGSHN,WGSHP,PetoleProteinC_node=node PetolSheth C,N,P,protein mass
    !       GRO,GrowthElms(ielmn),GrowthElms(ielmp)=PetolSheth C,N,P growth at each node
    !       CNWS,rProteinC2LeafP_pft=protein:N,protein:P ratios from startq.f
    !
    D505: DO KK=MNNOD,MXNOD
      K=pMOD(KK,MaxNodesPerBranch1)

      DO NE=1,NumPlantChemElms
        PetolShethElmntNode_brch(NE,K,NB,NZ)=PetolShethElmntNode_brch(NE,K,NB,NZ)+GrowthElms(NE)
      ENDDO
      PetoleProteinC_node(K,NB,NZ)=PetoleProteinC_node(K,NB,NZ)+&
        AMIN1(GrowthElms(ielmn)*rProteinC2LeafN_pft(NZ),GrowthElms(ielmp)*rProteinC2LeafP_pft(NZ))
      !
      !           SPECIFIC SHEATH OR PetolSheth LENGTH FUNCTION OF CURRENT MASS
      !           AT EACH NODE
      !
      !   SSL=specific length of PetolSheth growth
      !   ETOL=coefficient for etoliation effects on expansion,extension
      !   PetolShethLen2Mass_pft=growth in PetolSheth length vs mass from PFT file
      !   SSL2=parameter for calculating PetolSheth extension
      !   PetolShethElmntNode_brch=PetolSheth C mass
      !   PP=PFT population
      !   GSSL=allocation of PetolSheth length growth to each node
      !   TurgEff4LeafPetolExpansion=turgor expansion,extension function
      !   GROS,GRO=PetolSheth length,mass growth
      !   PetoleLength_node=PetolSheth length
      !
      IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.0.0_r8)THEN
        dMC                        = AZMAX1(1.e-12_r8,PetolShethElmntNode_brch(ielmc,K,NB,NZ)/(PlantPopuLive_pft(NZ)*GSSL))
        SSL                        = EtiolationCoeff*PetolShethLen2Mass_pft(NZ)*dMC**SSL2*TurgEff4LeafPetolExpansion
        GROS                       = GrowthElms(ielmc)/PlantPopuLive_pft(NZ)*SSL
        PetoleLength_node(K,NB,NZ) = PetoleLength_node(K,NB,NZ)+GROS*SinePetolShethAngle_pft(NZ)
      ENDIF
    ENDDO D505
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine GrowPetolShethOnBranch  

!----------------------------------------------------------------------------------------------------
  subroutine GrowStalkOnBranch(I,J,NZ,NB,GrowthStalk,EtiolationCoeff)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ,NB
  real(r8), INTENT(IN) :: GrowthStalk(NumPlantChemElms)    !gross growth rate on stalk, [g h-1]
  REAL(R8), INTENT(IN) :: EtiolationCoeff                  !etiolation coefficient
  
  character(len=*), parameter :: subname='GrowStalkOnBranch'  
  REAL(R8) :: GrowthElms(NumPlantChemElms)
  integer :: NN,MXNOD,MNNOD,K1,K2,NE,KK
  REAL(R8) :: ALLOCN,GNOD,SpecLenStalkGrowth
  real(r8) :: StalkLenGrowth,StalkAveStrutC_brch

  associate(                                                            &
    MainBranchNum_pft          => plt_morph%MainBranchNum_pft          ,& !input  :number of main branch,[-]  
    NumCogrowthNode_pft        => plt_morph%NumCogrowthNode_pft        ,& !input  :number of concurrently growing nodes,[-]
    NodeLenPergC_pft           => plt_morph%NodeLenPergC_pft           ,& !input  :internode length:mass during growth, [m gC-1]
    SineBranchAngle_pft        => plt_morph%SineBranchAngle_pft        ,& !input  :branching angle, [degree from horizontal]
    StalkStrutElms_brch        => plt_biom%StalkStrutElms_brch         ,& !input  :branch stalk structural element mass, [g d-2]
    PlantPopuLive_pft          => plt_site%PlantPopuLive_pft           ,& !input  :plant population, [d-2]
    iPlantCalendar_brch        => plt_pheno%iPlantCalendar_brch        ,& !input  :plant growth stage, [-]
    KHiestGroLeafNode_brch     => plt_pheno%KHiestGroLeafNode_brch     ,& !input  :leaf growth stage counter, [-]
    StalkNodeVertLength_brch   => plt_morph%StalkNodeVertLength_brch   ,& !inoput :internode height, [m]
    StalkNodeHeight_brch       => plt_morph%StalkNodeHeight_brch       ,& !inoput :internode height, [m]
    StructInternodeElms_brch   => plt_biom%StructInternodeElms_brch     & !inoput :internode chemical elelment, [g d-2]
  )
  call DebugPrint('beg '//subname//' NZ',NZ)
  !plant hasnot emerged
  IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).EQ.0)THEN  
    NN=0
    !plant emerged  
  ELSE
    NN=1
  ENDIF

  MXNOD = KHiestGroLeafNode_brch(NB,NZ)
  MNNOD = MAX(MIN(NN,MAX(NN,MXNOD-NumCogrowthNode_pft(NZ))),KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+2)
  MXNOD = MAX(MXNOD,MNNOD)

  IF(GrowthStalk(ielmc).GT.0.0_r8)THEN
    GNOD=MXNOD-MNNOD+1
    ALLOCN=1.0_r8/GNOD   !allocate uniformly to the different nodes
    DO NE=1,NumPlantChemElms
      GrowthElms(NE)=ALLOCN*GrowthStalk(NE)
    ENDDO
    !
    !     SPECIFIC INTERNODE LENGTH FUNCTION OF CURRENT STALK MASS
    !     AT EACH NODE
    !
    !     SpecLenStalkGrowth=specific length of stalk growth
    !     EtiolationCoeff=coefficient for etioliation effects on expansion,extension
    !     NodeLenPergC_pft=growth in stalk length vs mass from PFT file
    !     SNL2=parameter for calculating stalk extension
    !     WTSKB=stalk C mass
    !     PP=PFT population
    !     StalkLenGrowth,GRO=stalk length,mass growth
    !
    StalkAveStrutC_brch = AZMAX1(1.e-12_r8,StalkStrutElms_brch(ielmc,NB,NZ)/PlantPopuLive_pft(NZ))
    SpecLenStalkGrowth  = EtiolationCoeff*NodeLenPergC_pft(NZ)*StalkAveStrutC_brch**SNL2    
    StalkLenGrowth      = GrowthElms(ielmc)/PlantPopuLive_pft(NZ)*SpecLenStalkGrowth
    !write(9081,*)I*1000+J/24.,NB,GrowthElms(ielmc),PlantPopuLive_pft(NZ),StalkLenGrowth,SpecLenStalkGrowth,StalkAveStrutC_brch,EtiolationCoeff

    !
    !     GROWTH AT EACH CURRENT NODE
    !
    !     StructInternodeElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
    !     GRO,GrowthElms(ielmn),GrowthElms(ielmp)=stalk C,N,P growth at each node
    !     StalkNodeVertLength_brch,StalkNodeHeight_brch=stalk node length project in the vertical direction,stalk node height
    !     SineBranchAngle_pft=sine of stalk angle from horizontal from PFT file
    !
    D510: DO KK=MNNOD,MXNOD
      K1 = pMOD(KK,MaxNodesPerBranch1)
      K2 = pMOD(KK-1,MaxNodesPerBranch1)
      !grow stalk biomass 
      DO NE = 1, NumPlantChemElms
        StructInternodeElms_brch(NE,K1,NB,NZ)=StructInternodeElms_brch(NE,K1,NB,NZ)+GrowthElms(NE)
      ENDDO
      !
      !grow stalk length      
      StalkNodeVertLength_brch(K1,NB,NZ)=StalkNodeVertLength_brch(K1,NB,NZ)+StalkLenGrowth*SineBranchAngle_pft(NZ)
      IF(K1.NE.0)THEN
        StalkNodeHeight_brch(K1,NB,NZ)=StalkNodeHeight_brch(K2,NB,NZ)+StalkNodeVertLength_brch(K1,NB,NZ)
      ELSE
        StalkNodeHeight_brch(K1,NB,NZ)=StalkNodeVertLength_brch(K1,NB,NZ)
      ENDIF
    ENDDO D510
  ENDIF
  call PrintInfo('end '//subname)
  END associate  
  end subroutine GrowStalkOnBranch

!----------------------------------------------------------------------------------------------------
  subroutine RemobilizeBranch(I,J,NZ,NB,BegRemoblize,LRemob_brch,RCCC,RCCN,RCCP,RMxess_brch)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ,NB
  integer, intent(in) :: LRemob_brch,BegRemoblize
  real(r8), intent(in) :: RCCC,RCCN,RCCP
  real(r8), intent(inout) :: RMxess_brch
  integer :: NumRemobLeafNodes  
  real(r8) :: RespSenesTot_brch
  real(r8) :: FracRecyled  !ratio of phenologically-driven vs total senescence respiration
  integer  :: NBZ(MaxNumBranches)
  real(r8) :: XFRE(1:NumPlantChemElms)
  REAL(R8) :: RCO2V,RCO2e
  real(r8) :: RespPhenolSenes_brch  !respiration associated with senescence [gC h-1]
  integer :: NBY,NBX,NBL
  integer  :: NBK,NE
  character(len=*), parameter :: subname='RemobilizeBranch'

  associate(                                                           &
    fTCanopyGroth_pft         => plt_pheno%fTCanopyGroth_pft          ,& !input  :canopy temperature growth function, [-]
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft  ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    HoursDoingRemob_brch      => plt_pheno%HoursDoingRemob_brch       ,& !input  :counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
    KHiestGroLeafNode_brch    => plt_pheno%KHiestGroLeafNode_brch     ,& !input  :leaf growth stage counter, [-]
    iPlantPhenolPattern_pft   => plt_pheno%iPlantPhenolPattern_pft    ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft      ,& !input  :plant growth type (vascular, non-vascular),[-]
    isPlantBranchAlive_brch   => plt_pheno%isPlantBranchAlive_brch    ,& !input  :flag to detect branch death, [-]
    KLowestGroLeafNode_brch   => plt_pheno%KLowestGroLeafNode_brch    ,& !input  :leaf growth stage counter, [-]
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft          ,& !input  :number of branches,[-]
    iPlant2ndGrothPattern_pft => plt_pheno%iPlant2ndGrothPattern_pft  ,& !input  :plant expression of secondary growth, [-]                
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft          ,& !input  :number of main branch,[-]
    BranchNumerID_brch        => plt_morph%BranchNumerID_brch         ,& !input  :branch meric id, [-]
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch        ,& !input  :branch nonstructural element, [g d-2]
    CanopyLeafSheathC_brch    => plt_biom%CanopyLeafSheathC_brch      ,& !input  :plant branch leaf + sheath C, [g d-2]
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft              ,& !input  :threshold zero for plang growth calculation, [-]
    CanopyGrosRCO2_pft        => plt_bgcr%CanopyGrosRCO2_pft          ,& !inoput :canopy plant+nodule autotrophic respiraiton, [gC d-2]    
    CanopyResp_brch           => plt_bgcr%CanopyResp_brch             ,& !inoput :canopy respiration for a branch, [gC d-2 h-1]        
    CanopyRespC_CumYr_pft     => plt_bgcr%CanopyRespC_CumYr_pft       ,& !inoput :total autotrophic respiration, [gC d-2 ]    
    CO2NetFix_pft             => plt_bgcr%CO2NetFix_pft               ,& !inoput :canopy net CO2 exchange, [gC d-2 h-1]    
    ECO_ER_col                => plt_bgcr%ECO_ER_col                  ,& !inoput :ecosystem respiration, [g d-2 h-1]    
    Eco_AutoR_CumYr_col       => plt_bgcr%Eco_AutoR_CumYr_col         ,& !inoput :ecosystem autotrophic respiration, [g d-2 h-1]    
    StalkRsrvElms_brch        => plt_biom%StalkRsrvElms_brch           & !inoput :branch reserve element mass, [g d-2]
  )
  !     REMOBILIZATION OF STALK RESERVE C,N,P IF GROWTH RESPIRATION < 0
  call DebugPrint('beg '//subname//' NZ',NZ)
  !
  !     RMxess_brch=excess maintenance respiration
  !     WTRSVB=stalk reserve C mass
  !     RCO2V=remobilization of stalk reserve C
  !     VMXC=rate constant for nonstructural C oxidation in respiration
  !     fTCanopyGroth_pft=temperature function for canopy growth
  !

  IF(BegRemoblize.EQ.ifalse)THEN
    IF(RMxess_brch.GT.0.0_r8 .AND. StalkRsrvElms_brch(ielmc,NB,NZ).GT.0.0_r8)THEN
      RCO2V                           = AMIN1(RMxess_brch,VMXC*StalkRsrvElms_brch(ielmc,NB,NZ)*fTCanopyGroth_pft(NZ))
      StalkRsrvElms_brch(ielmc,NB,NZ) = StalkRsrvElms_brch(ielmc,NB,NZ)-RCO2V
      RMxess_brch                     = RMxess_brch-RCO2V

      CanopyGrosRCO2_pft(NZ)    = CanopyGrosRCO2_pft(NZ)-RCO2V
      CanopyResp_brch(NB,NZ)    = CanopyResp_brch(NB,NZ)-RCO2V
      CanopyRespC_CumYr_pft(NZ) = CanopyRespC_CumYr_pft(NZ)-RCO2V
      CO2NetFix_pft(NZ)         = CO2NetFix_pft(NZ)-RCO2V
      ECO_ER_col                = ECO_ER_col-RCO2V
      Eco_AutoR_CumYr_col       = Eco_AutoR_CumYr_col-RCO2V
    ENDIF
  ENDIF
!
!       TOTAL REMOBILIZATION = GROWTH RESPIRATION < 0 + DECIDUOUS LEAF
!       FALL DURING AUTUMN + REMOBILZATION DURING GRAIN FILL IN ANNUALS
!
!       iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!       LRemob_brch,BegRemoblize=remobilization flags
!       RespPhenolSenes_brch=phenologically-driven respiration senescence during late-season
!       RateK4ShootSeaStoreNonstEXfer=rate constant for plant-storage nonstructural C,N,P exchange
!       iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+PetolSheth,2=none,3=between 1,2
!       CanopyLeafSheathC_brch=leaf+PetolSheth mass
!       FLGZ=control rate of remobilization
!       Hours4FullSenes=number of hours until full senescence after physl maturity
!       RespSenesTot_brch=total senescence respiration
!       KHiestGroLeafNode_brch,KLowestGroLeafNode_brch=integer of highest,lowest leaf number currently growing
!       KSNC=number of nodes undergoing remobilization
!
  IF(BegRemoblize.EQ.itrue .AND. LRemob_brch.EQ.itrue .AND. iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial)THEN
    RespPhenolSenes_brch=RateK4ShootSeaStoreNonstEXfer(iPlantTurnoverPattern_pft(NZ))*CanopyLeafSheathC_brch(NB,NZ) &
      *AMIN1(1.0_r8,HoursDoingRemob_brch(NB,NZ)/Hours4FullSenes)
  ELSE
    RespPhenolSenes_brch=0._r8
  ENDIF
  RespSenesTot_brch=RMxess_brch+RespPhenolSenes_brch

  IF(RespSenesTot_brch.GT.ZERO4Groth_pft(NZ))THEN
    RCO2e=0._r8
    
    FracRecyled         = RespPhenolSenes_brch/RespSenesTot_brch
    NumRemobLeafNodes = INT(0.5_r8*(KHiestGroLeafNode_brch(NB,NZ)-KLowestGroLeafNode_brch(NB,NZ)))+1
    !
    !     TRANSFER NON-STRUCTURAL C,N,P FROM BRANCHES TO MAIN STEM
    !     IF MAIN STEM POOLS ARE DEPLETED
    !
    !     iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+PetolSheth,2=none,3=between 1,2
    !     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
    !
    IF(iPlantTurnoverPattern_pft(NZ).NE.0 .AND. is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)) &
      .AND. NB.EQ.MainBranchNum_pft(NZ) .AND. isclose(FracRecyled,0._r8))THEN
      NBY=0
      D584: DO NBL=1,NumOfBranches_pft(NZ)
        NBZ(NBL)=0
      ENDDO D584

      D586: DO NBL=1,NumOfBranches_pft(NZ)
        NBX=KHiestGroLeafNode_brch(NB,NZ)
        D585: DO NBK=1,NumOfBranches_pft(NZ)
          IF(isPlantBranchAlive_brch(NBK,NZ).EQ.iTrue .AND. NBK.NE.MainBranchNum_pft(NZ) &
            .AND. BranchNumerID_brch(NBK,NZ).LT.NBX .AND. BranchNumerID_brch(NBK,NZ).GT.NBY)THEN
            NBZ(NBL) = NBK
            NBX      = BranchNumerID_brch(NBK,NZ)
          ENDIF
        ENDDO D585

        IF(NBZ(NBL).NE.0)THEN
          NBY=BranchNumerID_brch(NBZ(NBL),NZ)
        ENDIF
      ENDDO D586

      D580: DO NBL=1,NumOfBranches_pft(NZ)
        IF(NBZ(NBL).NE.0)THEN          
          IF(BranchNumerID_brch(NBZ(NBL),NZ).NE.MainBranchNum_pft(NZ))THEN
            IF(CanopyNonstElms_brch(ielmc,NBZ(NBL),NZ).GT.0._r8)THEN
              XFRE(ielmc)=1.0E-02_r8*AMIN1(RespSenesTot_brch,CanopyNonstElms_brch(ielmc,NBZ(NBL),NZ))
              DO NE=2,NumPlantChemElms
                XFRE(NE)=XFRE(ielmc)*CanopyNonstElms_brch(NE,NBZ(NBL),NZ)/CanopyNonstElms_brch(ielmc,NBZ(NBL),NZ)
              ENDDO
            ELSE
              XFRE(ielmc)=0._r8
              DO NE=2,NumPlantChemElms
                XFRE(NE)=1.0E-02_r8*CanopyNonstElms_brch(NE,NBZ(NBL),NZ)
              ENDDO
            ENDIF

            DO NE=1,NumPlantChemElms
              CanopyNonstElms_brch(NE,NBZ(NBL),NZ)=CanopyNonstElms_brch(NE,NBZ(NBL),NZ)-XFRE(NE)
            ENDDO

            CanopyNonstElms_brch(ielmc,MainBranchNum_pft(NZ),NZ)=CanopyNonstElms_brch(ielmc,MainBranchNum_pft(NZ),NZ)+XFRE(ielmc)*FracRecyled
            DO NE=2,NumPlantChemElms
              CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)=CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)+XFRE(NE)
            ENDDO
            RCO2e             = RCO2e+XFRE(ielmc)*(1._r8-FracRecyled)
            RespSenesTot_brch = RespSenesTot_brch-XFRE(ielmc)
            IF(RespSenesTot_brch.LE.0.0_r8)exit
          ENDIF
        ENDIF
      ENDDO D580
    ENDIF

    IF(RespSenesTot_brch.GT.0.0_r8)then
      call RemobilizeLeafLayers(I,J,NumRemobLeafNodes,NB,nz,RespSenesTot_brch,RCCC,RCCN,RCCP,FracRecyled,RCO2e)
    endif  

    if(RCO2e>0._r8)then
      CanopyGrosRCO2_pft(NZ)    = CanopyGrosRCO2_pft(NZ)-RCO2e
      CanopyResp_brch(NB,NZ)    = CanopyResp_brch(NB,NZ)-RCO2e
      CanopyRespC_CumYr_pft(NZ) = CanopyRespC_CumYr_pft(NZ)-RCO2e
      CO2NetFix_pft(NZ)         = CO2NetFix_pft(NZ)-RCO2e
      ECO_ER_col                = ECO_ER_col-RCO2e
      Eco_AutoR_CumYr_col       = Eco_AutoR_CumYr_col-RCO2e
    endif      
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine RemobilizeBranch

!----------------------------------------------------------------------------------------------------
  subroutine SenescenceBranch(NZ,NB,RCCE)
  implicit none
  integer, intent(in) :: NZ,NB
  real(r8), intent(in) :: RCCE(NumPlantChemElms)

  character(len=*), parameter :: subname='SenescenceBranch'
  INTEGER :: K,KMinGroingLeafNodeNum,M,NE
  real(r8) :: RSpecKillLeafPetol,dFall
  real(r8) :: dRemoblE
  real(r8) :: RSpecKillLeaf               !specific leaf turnover rate
  real(r8) :: RSpecKillPetol              !specific petiole turnover rate [h-1]
  real(r8) :: PetolECycl(NumPlantChemElms),PetolEKill(NumPlantChemElms)
  real(r8) :: LeafEcycl(NumPlantChemElms),LeafEKill(NumPlantChemElms)
  
  associate(                                                                     &
    icwood                         => pltpar%icwood                             ,& !input  :group id of coarse woody litter
    ifoliar                        => pltpar%ifoliar                            ,& !input  :group id of plant foliar litter
    inonfoliar                     => pltpar%inonfoliar                         ,& !input  :group id of plant non-foliar litter group
    k_woody_comp                   => pltpar%k_woody_comp                       ,& !input  :woody litter complex id
    k_fine_comp                    => pltpar%k_fine_comp                        ,& !input  :fine litter complex id
    PlantElmAllocMat4Litr          => plt_soilchem%PlantElmAllocMat4Litr        ,& !input  :litter kinetic fraction, [-]
    FracPetolShethAlloc2Litr       => plt_allom%FracPetolShethAlloc2Litr        ,& !input  :leaf element allocation,[-]
    FracLeafShethElmAlloc2Litr     => plt_allom%FracLeafShethElmAlloc2Litr      ,& !input  :woody element allocation, [-]
    rProteinC2LeafN_pft            => plt_allom%rProteinC2LeafN_pft             ,& !input  :Plant leaf protein C to N ratio [g protein C (g leaf N)-1]
    rProteinC2LeafP_pft            => plt_allom%rProteinC2LeafP_pft             ,& !input  :Plant leaf protein C to P ratio [g protein C (g leaf P)-1]
    ZERO4Groth_pft                 => plt_biom%ZERO4Groth_pft                   ,& !input  :threshold zero for plang growth calculation, [-]
    doSenescence_brch              => plt_pheno%doSenescence_brch               ,& !input  :branch phenology flag, [-]
    doRemobilization_brch          => plt_pheno%doRemobilization_brch           ,& !input  :branch phenology flag, [-]
    RateRefLeafAppearance_pft      => plt_pheno%RateRefLeafAppearance_pft       ,& !input  :rate of leaf initiation, [h-1 at 25 oC]
    KHiestGroLeafNode_brch         => plt_pheno%KHiestGroLeafNode_brch          ,& !input  :leaf growth stage counter, [-]
    fTCanopyGroth_pft              => plt_pheno%fTCanopyGroth_pft               ,& !input  :canopy temperature growth function, [-]
    LitrfallElms_pvr               => plt_bgcr%LitrfallElms_pvr                 ,& !inoput :plant LitrFall element, [g d-2 h-1]
    LitrFallElms_brch              => plt_bgcr%LitrFallElms_brch                ,& !inoput :litterfall from the branch, [g d-2 h-1]
    SenecStalkStrutElms_brch       => plt_biom%SenecStalkStrutElms_brch         ,& !inoput :branch stalk structural element, [g d-2]
    PetolShethElmntNode_brch       => plt_biom%PetolShethElmntNode_brch         ,& !inoput :sheath chemical element, [g d-2]
    LeafProteinC_node              => plt_biom%LeafProteinC_node                ,& !inoput :layer leaf protein C, [g d-2]
    CanopyNonstElms_brch           => plt_biom%CanopyNonstElms_brch             ,& !inoput :branch nonstructural element, [g d-2]
    LeafStrutElms_brch             => plt_biom%LeafStrutElms_brch               ,& !inoput :branch leaf structural element mass, [g d-2]
    PetoleProteinC_node            => plt_biom%PetoleProteinC_node              ,& !inoput :layer sheath protein C, [g d-2]
    StructInternodeElms_brch       => plt_biom%StructInternodeElms_brch         ,& !inoput :internode C, [g d-2]
    LeafElmntNode_brch             => plt_biom%LeafElmntNode_brch               ,& !inoput :leaf element, [g d-2]
    PetolShethStrutElms_brch       => plt_biom%PetolShethStrutElms_brch         ,& !inoput :branch sheath structural element, [g d-2]
    LeafArea_node                  => plt_morph%LeafArea_node                   ,& !inoput :leaf area, [m2 d-2]
    LeafAreaLive_brch              => plt_morph%LeafAreaLive_brch               ,& !inoput :branch leaf area, [m2 d-2]
    PetoleLength_node              => plt_morph%PetoleLength_node               ,& !inoput :sheath height, [m]
    CanPBranchHeight               => plt_morph%CanPBranchHeight                ,& !output :branch height, [m]
    StalkNodeVertLength_brch       => plt_morph%StalkNodeVertLength_brch        ,& !output :internode height, [m]
    LeafAreaDying_brch             => plt_morph%LeafAreaDying_brch              ,& !output :branch leaf area, [m2 d-2]
    LeafSenescInitialElms_brch     => plt_pheno%LeafSenescInitialElms_brch      ,& !inoput :initial senescing leaf C/N/P mass, [g d-2]
    PetolSenescInitialElms_brch    => plt_pheno%PetolSenescInitialElms_brch     ,& !inoput :initial senescing sheath/petiole C/N/P mass, [g d-2]
    PetolShethChemElmRemobFlx_brch => plt_pheno%PetolShethChemElmRemobFlx_brch  ,& !output :cached remobilizable sheath element mass, [g d-2]
    LeafElmntRemobFlx_brch         => plt_pheno%LeafElmntRemobFlx_brch           & !output :cached remobilizable leaf element mass, [g d-2]
  )

  call DebugPrint('beg '//subname//' NZ',NZ)
  IF(doSenescence_brch(NB,NZ).EQ.itrue)THEN
    KMinGroingLeafNodeNum=MAX(0,KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+1)

    IF(KMinGroingLeafNodeNum.GT.0)THEN
      K    = pMOD(KMinGroingLeafNodeNum,MaxNodesPerBranch1)
      RSpecKillLeafPetol = fTCanopyGroth_pft(NZ)*RateRefLeafAppearance_pft(NZ)
      !
      !   REMOBILIZATION OF LEAF C,N,P ALSO DEPENDS ON STRUCTURAL C:N:P
      !
      !   doRemobilization_brch=flag for remobilization
      !   WGLF,WGLFN,WGLFP=node leaf C,N,P mass
      !   LeafArea_node=node leaf area
      !
      ! Cache the complete donor and its recyclable portion together. A negative snapshot
      ! marks initialization or a legacy restart, whose original snapshot cannot be recovered.
      IF(doRemobilization_brch(NB,NZ).EQ.itrue .OR. ANY(LeafSenescInitialElms_brch(:,NB,NZ).LT.0._r8))THEN

        LeafAreaDying_brch(NB,NZ)=AZMAX1(LeafArea_node(K,NB,NZ))
        LeafSenescInitialElms_brch(:,NB,NZ)=MAX(0._r8,LeafElmntNode_brch(:,K,NB,NZ))

        IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
          DO NE=1,NumPlantChemElms
            LeafElmntRemobFlx_brch(NE,NB,NZ) = LeafSenescInitialElms_brch(NE,NB,NZ)*MAX(0._r8,MIN(1._r8,RCCE(NE)))
          ENDDO
        ELSE
          LeafElmntRemobFlx_brch(1:NumPlantChemElms,NB,NZ)=0._r8
        ENDIF
      ENDIF
      !
      !       FRACTION OF CURRENT LEAF TO BE REMOBILIZED
      !
      !       RSpecKillLeafPetol,FSNCL=fraction of lowest leaf to be remobilized
      !
      ! Apply one fraction of the initial tissue to C/N/P, recycling, litter, and geometry.
      ! Bound it by every remaining donor pool; other processes may change tissue stoichiometry.
      RSpecKillLeaf=MAX(0._r8,MIN(1._r8,RSpecKillLeafPetol))
      IF(LeafSenescInitialElms_brch(ielmc,NB,NZ).LE.ZERO4Groth_pft(NZ))RSpecKillLeaf=0._r8
      DO NE=1,NumPlantChemElms
        IF(LeafSenescInitialElms_brch(NE,NB,NZ).GT.0._r8)THEN
          RSpecKillLeaf=MIN(RSpecKillLeaf,MAX(0._r8,LeafElmntNode_brch(NE,K,NB,NZ))/LeafSenescInitialElms_brch(NE,NB,NZ), &
            MAX(0._r8,LeafStrutElms_brch(NE,NB,NZ))/LeafSenescInitialElms_brch(NE,NB,NZ))
        ENDIF
      ENDDO
      IF(LeafAreaDying_brch(NB,NZ).GT.0._r8) &
        RSpecKillLeaf=MIN(RSpecKillLeaf,MAX(0._r8,LeafArea_node(K,NB,NZ))/LeafAreaDying_brch(NB,NZ))
      IF(LeafAreaDying_brch(NB,NZ).GT.0._r8) &
        RSpecKillLeaf=MIN(RSpecKillLeaf,MAX(0._r8,LeafAreaLive_brch(NB,NZ))/LeafAreaDying_brch(NB,NZ))
      !
      !       NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
      !       TO FRACTIONS SET IN 'STARTQ'
      !
      !       CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
      !       CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
      !       WGLFX,WGLFNX,WGLFPX=senescing leaf C,N,P mass
      !       FWODB=C woody fraction in other organs:0=woody,1=non-woody
      !       FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
      !
      DO NE=1,NumPlantChemElms
        LeafEKill(NE)=RSpecKillLeaf*LeafSenescInitialElms_brch(NE,NB,NZ)
        LeafEcycl(NE)=RSpecKillLeaf*LeafElmntRemobFlx_brch(NE,NB,NZ)
      ENDDO  

      D6300: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          dRemoblE = AZMAX1(LeafEKill(NE)-LeafEcycl(NE))
          dFall    = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dRemoblE*FracLeafShethElmAlloc2Litr(NE,k_woody_comp)

          LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) +dFall
          LitrFallElms_brch(NE,NB,NZ)              = LitrFallElms_brch(NE,NB,NZ)+dFall

          dFall=PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)*dRemoblE*FracLeafShethElmAlloc2Litr(NE,k_fine_comp)
          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ)+dFall
        ENDDO
      ENDDO D6300      
      !
      !       UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
      !
      !       LeafAreaLive_brch,LeafAreaDying_brch=branch living,senescing leaf area
      !       WTLFB,WTLFBN,WTLFBP,WGLFX,WGLFNX,WGLFPX=C,N,P mass in living,senescing leaf
      !       LeafProteinC_node=leaf protein mass
      !       CNWS,rProteinC2LeafP_pft=protein:N,protein:P ratios from startq.f
      !       CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
      !     canopy nonstructural Carbon/nutrient is not stored in leaves
      LeafAreaLive_brch(NB,NZ) = MAX(0._r8,LeafAreaLive_brch(NB,NZ)-RSpecKillLeaf*LeafAreaDying_brch(NB,NZ))
      LeafArea_node(K,NB,NZ)   = MAX(0._r8,LeafArea_node(K,NB,NZ)-RSpecKillLeaf*LeafAreaDying_brch(NB,NZ))
      LeafProteinC_node(K,NB,NZ) = AZMAX1(LeafProteinC_node(K,NB,NZ)-AMAX1(LeafEKill(ielmn)*rProteinC2LeafN_pft(NZ),LeafEKill(ielmp)*rProteinC2LeafP_pft(NZ)))

      DO NE=1,NumPlantChemElms
        LeafStrutElms_brch(NE,NB,NZ)   = MAX(0._r8,LeafStrutElms_brch(NE,NB,NZ)-LeafEKill(NE))
        LeafElmntNode_brch(NE,K,NB,NZ) = MAX(0._r8,LeafElmntNode_brch(NE,K,NB,NZ)-LeafEKill(NE))
        CanopyNonstElms_brch(NE,NB,NZ) = CanopyNonstElms_brch(NE,NB,NZ)+LeafEcycl(NE)
      ENDDO
      !
      !       REMOBILIZATION OF SHEATHS OR PetolSheth C,N,P ALSO DEPENDS ON
      !       STRUCTURAL C:N:P
      !
      !       doRemobilization_brch=flag for remobilization
      !       PetolShethElmntNode_brch,WGSHN,WGSHP=node PetolSheth C,N,P mass
      !       CanPBranchHeight=PetolSheth length
      !       RCES(ielmc)X,RCES(ielmn)X,RCES(ielmp)X=remobilization of C,N,P from senescing PetolSheth
      !

      ! Cache the complete donor and its recyclable portion together. A negative snapshot
      ! marks initialization or a legacy restart, whose original snapshot cannot be recovered.
      IF(doRemobilization_brch(NB,NZ).EQ.itrue .OR. ANY(PetolSenescInitialElms_brch(:,NB,NZ).LT.0._r8))THEN
        CanPBranchHeight(NB,NZ)=AZMAX1(PetoleLength_node(K,NB,NZ))
        PetolSenescInitialElms_brch(:,NB,NZ)=MAX(0._r8,PetolShethElmntNode_brch(:,K,NB,NZ))
        !
        IF(PetolShethElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
          DO NE=1,NumPlantChemElms
            PetolShethChemElmRemobFlx_brch(NE,NB,NZ) = PetolSenescInitialElms_brch(NE,NB,NZ)*MAX(0._r8,MIN(1._r8,RCCE(NE)))
          ENDDO
        ELSE
          PetolShethChemElmRemobFlx_brch(1:NumPlantChemElms,NB,NZ)=0._r8
        ENDIF
        
      ENDIF
      IF(doRemobilization_brch(NB,NZ).EQ.itrue)THEN
        DO NE=1,NumPlantChemElms
          SenecStalkStrutElms_brch(NE,NB,NZ)=SenecStalkStrutElms_brch(NE,NB,NZ)+StructInternodeElms_brch(NE,K,NB,NZ)
        ENDDO
        StructInternodeElms_brch(1:NumPlantChemElms,K,NB,NZ) = 0._r8
        StalkNodeVertLength_brch(K,NB,NZ)                    = 0._r8
      ENDIF
      !
      !       FRACTION OF CURRENT SHEATH TO BE REMOBILIZED
      !
      !       RSpecKillPetol=fraction of lowest PetolSheth to be remobilized
      !
      ! Apply one fraction of the initial tissue to C/N/P, recycling, litter, and geometry.
      ! Bound it by every remaining donor pool; other processes may change tissue stoichiometry.
      RSpecKillPetol=MAX(0._r8,MIN(1._r8,RSpecKillLeafPetol))
      IF(PetolSenescInitialElms_brch(ielmc,NB,NZ).LE.ZERO4Groth_pft(NZ))RSpecKillPetol=0._r8
      DO NE=1,NumPlantChemElms
        IF(PetolSenescInitialElms_brch(NE,NB,NZ).GT.0._r8)THEN
          RSpecKillPetol=MIN(RSpecKillPetol,MAX(0._r8,PetolShethElmntNode_brch(NE,K,NB,NZ))/PetolSenescInitialElms_brch(NE,NB,NZ), &
            MAX(0._r8,PetolShethStrutElms_brch(NE,NB,NZ))/PetolSenescInitialElms_brch(NE,NB,NZ))
        ENDIF
      ENDDO
      IF(CanPBranchHeight(NB,NZ).GT.0._r8) &
        RSpecKillPetol=MIN(RSpecKillPetol,MAX(0._r8,PetoleLength_node(K,NB,NZ))/CanPBranchHeight(NB,NZ))
      !
      !       NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
      !       TO FRACTIONS SET IN 'STARTQ'
      !
      !       CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
      !       CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
      !       RSpecKillPetol=fraction of lowest PetolSheth to be remobilized
      !       RCES(ielmc)X,RCES(ielmn)X,RCES(ielmp)X=remobilization of C,N,P from senescing PetolSheth
      !       WGSHX,WGSHNX,WGSHPX=senescing PetolSheth C,N,P mass
      !       FWODB=C woody fraction in other organs:0=woody,1=non-woody
      !       FWODSN,FWODSP=N,P woody fraction in PetolSheth:0=woody,1=non-woody
      !
      DO NE=1,NumPlantChemElms
        PetolEKill(NE) = RSpecKillPetol*PetolSenescInitialElms_brch(NE,NB,NZ)
        PetolECycl(NE) = RSpecKillPetol*PetolShethChemElmRemobFlx_brch(NE,NB,NZ)
      ENDDO
      D6305: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          dRemoblE                                 = AZMAX1(PetolEKill(NE)-PetolECycl(NE))
          dFall                                    = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dRemoblE*FracPetolShethAlloc2Litr(NE,k_woody_comp)
          LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)              = LitrFallElms_brch(NE,NB,NZ)+dFall

          dFall                                   = PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*dRemoblE*FracPetolShethAlloc2Litr(NE,k_fine_comp)
          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+dFall
          LitrFallElms_brch(NE,NB,NZ)             = LitrFallElms_brch(NE,NB,NZ)+dFall
        ENDDO    
      ENDDO D6305      
      !
      !       UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
      !
      !       RSpecKillPetol=fraction of lowest PetolSheth to be remobilized
      !       PetoleLength_node,CanPBranchHeight=living,senescing PetolSheth length
      !       WTSHB,WTSHBN,WTSHBP,WGSHEX,WGSHNX,WGSHPX=C,N,P mass in living,senescing PetolSheth
      !       PetoleProteinC_node=PetolSheth protein mass
      !       CNWS,rProteinC2LeafP_pft=protein:N,protein:P ratios from startq.f
      !       CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
      !       RCES(ielmc)X,RCES(ielmn)X,RCES(ielmp)X=remobilization of C,N,P from senescing PetolSheth
      !
      DO NE=1,NumPlantChemElms
        PetolShethStrutElms_brch(NE,NB,NZ)   = AZMAX1(PetolShethStrutElms_brch(NE,NB,NZ)-PetolEKill(NE))
        PetolShethElmntNode_brch(NE,K,NB,NZ) = AZMAX1(PetolShethElmntNode_brch(NE,K,NB,NZ)-PetolEKill(NE))
        CanopyNonstElms_brch(NE,NB,NZ)       = CanopyNonstElms_brch(NE,NB,NZ)+PetolECycl(NE)
      ENDDO
      PetoleLength_node(K,NB,NZ)   = MAX(0._r8,PetoleLength_node(K,NB,NZ)-RSpecKillPetol*CanPBranchHeight(NB,NZ))
      PetoleProteinC_node(K,NB,NZ) = AZMAX1(PetoleProteinC_node(K,NB,NZ)-AMAX1(PetolEKill(ielmn)*rProteinC2LeafN_pft(NZ),PetolEKill(ielmp)*rProteinC2LeafP_pft(NZ)))
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  END associate
  end subroutine SenescenceBranch  
  ![tail]
end module PlantBranchMod
