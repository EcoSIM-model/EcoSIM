module PlantDisturbsMod
!
!! Description:
! code to apply distance to plants
  use data_kind_mod,      only: r8 => DAT_KIND_R8, yearIJ_type
  use minimathmod,        only: isclose, AZMAX1
  use EcoSIMCtrlDataType, only: DazCurrYear
  use PlantBalMod,        only: SumPlantBranchBiome,SumPlantBiomStates,SumRootBiome,CheckPlantBalanceZ,SumPlantBiome,SumLitfallBlg
  use LitterFallMod, only : SetDeadPlant
  use InitPlantMod,  only : InitPlantPhenoMorphoBio,InitRootMychorMorphoBio  
  use PlantDebugMod
  use ElmIDMod
  use EcosimConst
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantPhotosynthesisAPIData, only : plt_photo
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantSoilChemistryAPIData, only : plt_soilchem
  use PlantAllometryAPIData, only : plt_allom
  use PlantBiomassAPIData, only : plt_biom
  use PlantEnergyWaterAPIData, only : plt_ew
  use PlantDisturbanceAPIData, only : plt_distb
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantRootBGCAPIData, only : plt_rbgc
  use PlantBGCPars
  use PlantMathFuncMod
  use PlantDisturbByFireMod
  use PlantDisturbByTillageMod
  use PlantDisturbByGrazingMod
  implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
! end_include_section

! disturbance variables
  real(r8) :: CanopyNonstElmRemoval(NumPlantChemElms)
  real(r8) :: LeafElmntRemoval(NumPlantChemElms)
  real(r8) :: PetoleShethElmntRemoval(NumPlantChemElms)
  real(r8) :: WoodyElmntRemoval(NumPlantChemElms)
  real(r8) :: StandeadElmntRemoval(NumPlantChemElms)
  real(r8) :: LeafElmnt2Litr(NumPlantChemElms)
  real(r8) :: FineNonleafElmnt2Litr(NumPlantChemElms)
  real(r8) :: WoodyElmnt2Litr(NumPlantChemElms)
  real(r8) :: StandeadElmnt2Litr(NumPlantChemElms)
  real(r8) :: LeafElmntHarv2Litr(NumPlantChemElms)
  real(r8) :: PetolShethElmntHarv2Litr(NumPlantChemElms)
  real(r8) :: WoodyElmntHarv2Litr(NumPlantChemElms)      !harvested woody element goes to litter
  real(r8) :: StandeadElmntHarv2Litr(NumPlantChemElms)   !harvested standing dead that turns into litter
  real(r8) :: GrainHarvst(NumPlantChemElms)

  public :: RemoveBiomassByDisturbance
  public :: InitPlantDisturbance
  public :: StageDisturbances
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine InitPlantDisturbance
  implicit none
  call InitPlantFireMod
  end subroutine InitPlantDisturbance

!----------------------------------------------------------------------------------------------------
  subroutine StageDisturbances(yearIJ,NZ)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ    
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='StageDisturbances'

  call PrintInfo('beg '//subname)
  StandeadElmntRemoval(1:NumPlantChemElms)     = 0._r8
  CanopyNonstElmRemoval(1:NumPlantChemElms)    = 0._r8
  LeafElmntRemoval(1:NumPlantChemElms)         = 0._r8
  PetoleShethElmntRemoval(1:NumPlantChemElms)  = 0._r8
  WoodyElmntRemoval(1:NumPlantChemElms)        = 0._r8
  LeafElmnt2Litr(1:NumPlantChemElms)           = 0._r8
  FineNonleafElmnt2Litr(1:NumPlantChemElms)    = 0._r8
  WoodyElmnt2Litr(1:NumPlantChemElms)          = 0._r8
  StandeadElmnt2Litr(1:NumPlantChemElms)       = 0._r8
  LeafElmntHarv2Litr(1:NumPlantChemElms)       = 0._r8
  PetolShethElmntHarv2Litr(1:NumPlantChemElms) = 0._r8
  WoodyElmntHarv2Litr(1:NumPlantChemElms)      = 0._r8
  GrainHarvst(1:NumPlantChemElms)              = 0._r8
  call PrintInfo('end '//subname)
  end subroutine StageDisturbances
!----------------------------------------------------------------------------------------------------
  subroutine RemoveBiomByMgmt(yearIJ,NZ)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  integer, intent(in) :: NZ
  integer :: K
  real(r8) ::  massroot1(NumPlantChemElms)
  real(r8) :: tvegE(NumPlantChemElms),tvegE1(NumPlantChemElms)
  character(len=*), parameter :: subname='RemoveBiomByMgmt'

  !     TRANSFER ABOVE-GROUND C,N,P AT HARVEST OR DISTURBANCE
  !
  call PrintInfo('beg '//subname)

  !prepare for disturbance
  CALL SumPlantBiomStates(yearIJ,NZ,subname)
  !
  call RemoveBiomByHarvest(yearIJ,NZ)
  !
  !     REDUCE OR REMOVE PLANT POPULATIONS DURING TILLAGE
  !
  call RemoveBiomByTillage(yearIJ,NZ)
!  IF(NZ.EQ.2)call MassBalCheck(yearIJ,NZ,opt='exit',checktype='c',iut=992,info=subname//'1')  
  !
  call PrintInfo('end '//subname)
  end subroutine RemoveBiomByMgmt

!----------------------------------------------------------------------------------------------------
  subroutine RemoveDeadAnnual(yearIJ,NZ)
  !
  !remove annual plant biomass 
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='RemoveDeadAnnual'

  logical :: checkDroughtDeciduos
  real(r8) :: dshootLoss(NumPlantChemElms)
  real(r8) :: drootLoss(NumPlantChemElms)
  integer :: NE,M,NB,L,N
  associate(                                                               &
    ifoliar                     => pltpar%ifoliar                         ,& !input  :group id of plant foliar litter    
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft            ,& !input  :number of branches,[-]    
    LeafStrutElms_brch          => plt_biom%LeafStrutElms_brch            ,& !input  :branch leaf structural element mass, [g d-2]    
    PlantElmAllocMat4Litr       => plt_soilchem%PlantElmAllocMat4Litr     ,& !input  :litter kinetic fraction, [-]
    PetolShethStrutElms_brch    => plt_biom%PetolShethStrutElms_brch      ,& !input  :branch sheath structural element, [g d-2]    
    StalkStrutElms_brch         => plt_biom%StalkStrutElms_brch           ,& !input  :branch stalk structural element mass, [g d-2]
    StalkRsrvElms_brch          => plt_biom%StalkRsrvElms_brch            ,& !input  :branch reserve element mass, [g d-2]
    HuskStrutElms_brch          => plt_biom%HuskStrutElms_brch            ,& !input  :branch husk structural element mass, [g d-2]
    EarStrutElms_brch           => plt_biom%EarStrutElms_brch             ,& !input  :branch ear structural chemical element mass, [g d-2]
    GrainStrutElms_brch         => plt_biom%GrainStrutElms_brch           ,& !input  :branch grain structural element mass, [g d-2]    
    k_fine_comp                 => pltpar%k_fine_comp                     ,& !input  :fine litter complex id    
    MaxNumRootLays              => plt_site%MaxNumRootLays                ,& !input  :maximum root layer number,[-]        
    NU                          => plt_site%NU                            ,& !input  :current soil surface layer number, [-]    
    LitrfallElms_pvr            => plt_bgcr%LitrfallElms_pvr              ,& !inoput :plant LitrFall element, [g d-2 h-1]    
    SeasonalNonstElms_pft       => plt_biom%SeasonalNonstElms_pft         ,& !input: plant stored nonstructural element at current step, [g d-2]  
    MainBranchNum_pft           => plt_morph%MainBranchNum_pft            ,& !input  :id number of main branch,[-]    
    Myco_pft                    => plt_morph%Myco_pft                     ,& !input  :mycorrhizal type (no or yes),[-]    
    iPlantPhenolPattern_pft     => plt_pheno%iPlantPhenolPattern_pft      ,& !input  :plant growth habit: annual or perennial,[-]    
    iPlantPhenolType_pft        => plt_pheno%iPlantPhenolType_pft         ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]    
    RootNodulStrutElms_rpvr     => plt_biom%RootNodulStrutElms_rpvr       ,& !input  :root layer nodule element, [g d-2]
    RootNodulNonstElms_rpvr     => plt_biom%RootNodulNonstElms_rpvr       ,& !input  :root layer nonstructural element, [g d-2]
    RootMyco2ndStrutElms_rpvr   => plt_biom%RootMyco2ndStrutElms_rpvr     ,& !input  :root layer element secondary axes, [g d-2]    
    RootMyco1stStrutElms_rpvr   => plt_biom%RootMyco1stStrutElms_rpvr     ,& !input  :root layer element primary axes, [g d-2]    
    NGTopRootLayer_pft          => plt_morph%NGTopRootLayer_pft           ,& !input  :soil layer at planting depth, [-]    
    NumStructuralRootAxes_pft        => plt_morph%NumStructuralRootAxes_pft         ,& !input  :number of structural root axes,[-]    
    RootMycoNonstElms_rpvr      => plt_biom%RootMycoNonstElms_rpvr        ,& !inoput :root layer nonstructural element, [g d-2]    
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch          ,& !inoput :branch nonstructural element, [g d-2]    
    Days4FalseBreak_pft         => plt_pheno%Days4FalseBreak_pft          ,& !inoput :accumulated days to singifying false break
    SeasonalNonstCDayAve_pft    => plt_biom%SeasonalNonstCDayAve_pft       & !inoput: daily average seasonal storage C for annual plant death check, [g d-2]   
  )
  call PrintInfo('beg '//subname)
  checkDroughtDeciduos=iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid

  IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .and. checkDroughtDeciduos)THEN
    SeasonalNonstCDayAve_pft(NZ)=SeasonalNonstCDayAve_pft(NZ)+SeasonalNonstElms_pft(ielmc,NZ)/24._r8
    !Evaluate the completed daily mean and count consecutive qualifying days.
    if(yearIJ%J.EQ.24)then
      if(eval_annual_false_break_death(yearIJ,NZ))then
        Days4FalseBreak_pft(NZ)=Days4FalseBreak_pft(NZ)+1
      else
        Days4FalseBreak_pft(NZ)=0
      endif
    endif
    
    if(yearIJ%J.EQ.24 .and. Days4FalseBreak_pft(NZ).GE.Days2CallFalseBreak)then
      !add CanopyNonstElms_brch to below ground litter at layer NGTopRootLayer_pft(NZ)

      DO M=1,jsken
        DO  NB=1,NumOfBranches_pft(NZ)
          DO NE=1,NumPlantChemElms
            LitrfallElms_pvr(NE,M,k_fine_comp,NGTopRootLayer_pft(NZ),NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,NGTopRootLayer_pft(NZ),NZ)+&
              (CanopyNonstElms_brch(NE,NB,NZ)+LeafStrutElms_brch(NE,NB,NZ)+PetolShethStrutElms_brch(NE,NB,NZ)+ &
              StalkStrutElms_brch(NE,NB,NZ)+StalkRsrvElms_brch(NE,NB,NZ)+HuskStrutElms_brch(NE,NB,NZ)+&
              EarStrutElms_brch(NE,NB,NZ)+GrainStrutElms_brch(NE,NB,NZ))*PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)
            dshootLoss(NE) =dshootLoss(NE) +(CanopyNonstElms_brch(NE,NB,NZ)+LeafStrutElms_brch(NE,NB,NZ)+PetolShethStrutElms_brch(NE,NB,NZ)+ &
              StalkStrutElms_brch(NE,NB,NZ)+StalkRsrvElms_brch(NE,NB,NZ) + &
              HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ)+GrainStrutElms_brch(NE,NB,NZ))*PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)
          ENDDO
        ENDDO

        DO L=NU,MaxNumRootLays
          DO NE=1,NumPlantChemElms        
            DO N=1,Myco_pft(NZ)
              LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+&
              (sum(RootMyco2ndStrutElms_rpvr(NE,N,L,1:NumStructuralRootAxes_pft(NZ),NZ)) + &              
                RootMycoNonstElms_rpvr(NE,N,L,NZ))*PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)
              drootLoss(NE) =drootLoss(NE)+(sum(RootMyco2ndStrutElms_rpvr(NE,N,L,1:NumStructuralRootAxes_pft(NZ),NZ)) + &              
                RootMycoNonstElms_rpvr(NE,N,L,NZ))*PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)
            ENDDO  
            LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+&
              (sum(RootMyco1stStrutElms_rpvr(NE,L,1:NumStructuralRootAxes_pft(NZ),NZ))  &
              +RootNodulStrutElms_rpvr(NE,L,NZ)+ RootNodulNonstElms_rpvr(NE,L,NZ))*PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)            
            drootLoss(NE) =drootLoss(NE)+(sum(RootMyco1stStrutElms_rpvr(NE,L,1:NumStructuralRootAxes_pft(NZ),NZ))  &
              +RootNodulStrutElms_rpvr(NE,L,NZ)+ RootNodulNonstElms_rpvr(NE,L,NZ))*PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)              
          ENDDO  
        ENDDO  
      ENDDO 
      !      
      call SetDeadPlant(yearIJ,NZ,NumOfBranches_pft(NZ))
      !
      plt_pheno%doReSeed_pft(NZ)=.FALSE.

      call InitPlantPhenoMorphoBio(NZ)
      !
      call InitRootMychorMorphoBio(NZ)
    endif
  endif  

  call PrintInfo('end '//subname)
  end associate    
  END subroutine RemoveDeadAnnual
!----------------------------------------------------------------------------------------------------
  function eval_annual_false_break_death(yearIJ,NZ)result(tokill)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  logical :: tokill
  associate(                                                           &
    CanopyHeightLive_pft       => plt_morph%CanopyHeightLive_pft     , & !inoput :canopy height, [m]    
    SeasonalNonstCDayAve_pft   => plt_biom%SeasonalNonstCDayAve_pft  , & !input : daily average seasonal storage C for annual plant death check, [g d-2]
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft            , & !input :threshold zero for plang growth calculation, [-]    
    Hours2LeafOut_brch         => plt_pheno%Hours2LeafOut_brch       , & !inoput:counter for mobilizing nonstructural C during spring leafout/dehardening, [h]    
    NumOfBranches_pft          => plt_morph%NumOfBranches_pft        , & !input :number of branches,[-]    
    iPlantPhenolPattern_pft    => plt_pheno%iPlantPhenolPattern_pft  , & !input :plant growth habit: annual or perennial,[-]        
    SeasonalNonstElms_pft      => plt_biom%SeasonalNonstElms_pft     , & !inoput:plant stored nonstructural element at current step, [g d-2]  
    CanopyNonstElms_brch       => plt_biom%CanopyNonstElms_brch        & !input :branch nonstructural element, [g d-2]
  )
  
  tokill= any(CanopyNonstElms_brch(ielmc,1:NumOfBranches_pft(NZ),NZ)>0._r8) &
    .and. isclose(SeasonalNonstCDayAve_pft(NZ),SeasonalNonstElms_pft(ielmc,NZ)) &
    .and. any(Hours2LeafOut_brch(1:NumOfBranches_pft(NZ),NZ).GT.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ))) &
    .and. SeasonalNonstElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ) .and. isclose(CanopyHeightLive_pft(NZ),0._r8)

  end associate
  end function eval_annual_false_break_death
!----------------------------------------------------------------------------------------------------
  subroutine RemoveStandingDead(yearIJ,NZ)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='RemoveStandingDead'

  real(r8) :: FracStdeadLeft   !fraction of standing dead left after removal
  real(r8) :: FHVSH            !fraction of standing dead left after removal
  integer :: M,NE
  associate(                                                    &
    SolarNoonHour_col      => plt_site%SolarNoonHour_col       ,& !input  :time of solar noon, [h]
    FracBiomHarvsted       => plt_distb%FracBiomHarvsted       ,& !input  :harvest efficiency, [-]
    THIN_pft               => plt_distb%THIN_pft               ,& !input  :thinning of plant population, [-]
    iHarvstType_pft        => plt_distb%iHarvstType_pft        ,& !input  :type of harvest,[-]
    StandDeadCompKElms_pft => plt_biom%StandDeadCompKElms_pft   & !inoput :standing dead element fraction, [g d-2]
  )
  call PrintInfo('beg '//subname)
  StandeadElmntRemoval(1:NumPlantChemElms)=0._r8
  StandeadElmntHarv2Litr(1:NumPlantChemElms)=0._r8

  IF(yearIJ%J.EQ.INT(SolarNoonHour_col) .AND. iHarvstType_pft(NZ).NE.iharvtyp_grazing &
    .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(isclose(THIN_pft(NZ),0._r8))THEN
      FracStdeadLeft = AZMAX1(1._r8-FracBiomHarvsted(iHarvst_pft,iplthvst_stdead,NZ))
      FHVSH = FracStdeadLeft
    ELSE
      FracStdeadLeft=AZMAX1(1._r8-THIN_pft(NZ))
      IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
        FHVSH=AZMAX1(1._r8-FracBiomHarvsted(iHarvst_pft,iplthvst_stdead,NZ)*THIN_pft(NZ))
      ELSE
        FHVSH=FracStdeadLeft
      ENDIF
    ENDIF
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
    call RemoveStandDeadByGrazing(yearIJ,NZ,FracStdeadLeft,FHVSH)
  ELSE
    FracStdeadLeft = 1.0_r8
    FHVSH          = 1.0_r8
  ENDIF

  D6475: DO M=1,jsken
    DO NE=1,NumPlantChemElms
      StandeadElmntRemoval(NE)        = StandeadElmntRemoval(NE)+(1._r8-FHVSH)*StandDeadCompKElms_pft(NE,M,NZ)
      StandeadElmntHarv2Litr(NE)      = StandeadElmntHarv2Litr(NE)+(FHVSH-FracStdeadLeft)*StandDeadCompKElms_pft(NE,M,NZ)
      StandDeadCompKElms_pft(NE,M,NZ) = FracStdeadLeft*StandDeadCompKElms_pft(NE,M,NZ)
    ENDDO
  ENDDO D6475
  call PrintInfo('end '//subname)
  end associate
  end subroutine RemoveStandingDead  

!----------------------------------------------------------------------------------------------------
  subroutine RemoveBiomassByDisturbance(yearIJ,NZ)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  integer , intent(in) :: NZ
  character(len=*), parameter :: subname='RemoveBiomassByDisturbance'
  real(r8) :: tvegE(NumPlantChemElms),canopyE(NumPlantChemElms),RootE(NumPlantChemElms)
  real(r8) :: tvegE1(NumPlantChemElms),canopyE1(NumPlantChemElms),RootE1(NumPlantChemElms)

!     begin_execution
  associate(                                                      &
    NU                      => plt_site%NU                       ,& !input  :current soil surface layer number, [-]
    ZERO                    => plt_site%ZERO                     ,& !input  :threshold zero for numerical stability, [-]
    AREA3                   => plt_site%AREA3                    ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    iHarvstType_pft         => plt_distb%iHarvstType_pft         ,& !input  :type of harvest,[-]
    PlantPopuLive_pft       => plt_site%PlantPopuLive_pft        ,& !input  :plant population, [d-2]
    ZERO4Uptk_pft           => plt_rbgc%ZERO4Uptk_pft            ,& !output :threshold zero for uptake calculation, [-]
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft           ,& !output :threshold zero for plang growth calculation, [-]
    ZERO4LeafVar_pft        => plt_biom%ZERO4LeafVar_pft          & !output :threshold zero for leaf calculation, [-]
  )
  call PrintInfo('beg '//subname)

  call StageDisturbances(yearIJ,NZ)

  CALL RemoveDeadAnnual(yearIJ,NZ)
  !
  !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
  !                       ,3=pruning,4=grazing,5=fire,6=herbivory
  !    
  call RemoveBiomByMgmt(yearIJ,NZ)  
  
  IF(iHarvstType_pft(NZ).GE.iharvtyp_none)THEN
    !
    CALL RemoveStandingDead(yearIJ,NZ)
    !
    call PlantShootDisturbance(yearIJ,NZ)
  
    ZERO4Groth_pft(NZ)   = ZERO*PlantPopuLive_pft(NZ)
    ZERO4Uptk_pft(NZ)    = ZERO*PlantPopuLive_pft(NZ)/AREA3(NU)
    ZERO4LeafVar_pft(NZ) = ZERO*PlantPopuLive_pft(NZ)*1.0E+06_r8
    
  ENDIF
  call PrintInfo('end '//subname) 
  end associate
  end subroutine RemoveBiomassByDisturbance

!----------------------------------------------------------------------------------------------------
  subroutine PlantShootDisturbance(yearIJ,NZ)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  integer, intent(in) :: NZ

  character(len=*), parameter :: subname='PlantShootDisturbance'

  real(r8) :: NonstructElmntOffEcosystem(NumPlantChemElms)
  real(r8) :: LeafElmntOffEcosystem(NumPlantChemElms)
  real(r8) :: FineNonleafElmOffEcosystem(NumPlantChemElms)
  real(r8) :: WoodyElmntOffEcosystem(NumPlantChemElms)
  real(r8) :: StandeadElmntOffEcosystem(NumPlantChemElms)
  real(r8) :: CanopyNonstElm2Litr(NumPlantChemElms)
  real(r8) :: TotalAbgElmnt2Litr(NumPlantChemElms)
  real(r8) :: HarvestAbgElmnt2Litr(NumPlantChemElms)
  real(r8) :: tvegE(NumPlantChemElms)
  real(r8) :: tvegE1(NumPlantChemElms)

  call PrintInfo('beg '//subname)
  NonstructElmntOffEcosystem(1:NumPlantChemElms) = 0._r8
  LeafElmntOffEcosystem(1:NumPlantChemElms)      = 0._r8
  FineNonleafElmOffEcosystem(1:NumPlantChemElms) = 0._r8
  WoodyElmntOffEcosystem(1:NumPlantChemElms)     = 0._r8
  StandeadElmntOffEcosystem(1:NumPlantChemElms)  = 0._r8
  CanopyNonstElm2Litr(1:NumPlantChemElms)        = 0._r8

  call ApplyDisturbanceBiomRemoval(yearIJ,NZ,CanopyNonstElm2Litr,NonstructElmntOffEcosystem,&
    LeafElmntOffEcosystem,FineNonleafElmOffEcosystem,WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)

  !     TOTAL C,N,P REMOVAL FROM DISTURBANCE
  call AbvgBiomRemovalByDisturb(yearIJ,NZ,CanopyNonstElm2Litr,HarvestAbgElmnt2Litr,TotalAbgElmnt2Litr)
  !
  !     ABOVE-GROUND LitrFall FROM HARVESTING
  !
  call LiterfallByDisturbance(yearIJ,NZ,HarvestAbgElmnt2Litr,TotalAbgElmnt2Litr,CanopyNonstElm2Litr,&
    NonstructElmntOffEcosystem,LeafElmntOffEcosystem,FineNonleafElmOffEcosystem,&
    WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)

  call PrintInfo('end '//subname)
  end subroutine PlantShootDisturbance

!----------------------------------------------------------------------------------------------------
  subroutine LiterfallByDisturbance(yearIJ,NZ,HarvestAbgElmnt2Litr,TotalAbgElmnt2Litr,CanopyNonstElm2Litr,&
    NonstructElmntOffEcosystem,LeafElmntOffEcosystem,FineNonleafElmOffEcosystem,&
    WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)

  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  integer , intent(in) :: NZ
  real(r8), intent(in) :: HarvestAbgElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: TotalAbgElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: CanopyNonstElm2Litr(NumPlantChemElms)
  real(r8), intent(in) :: NonstructElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(in) :: LeafElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(in) :: FineNonleafElmOffEcosystem(NumPlantChemElms)
  real(r8), intent(in) :: WoodyElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(in) :: StandeadElmntOffEcosystem(NumPlantChemElms)
  character(len=*), parameter :: subname='LiterfallByDisturbance'
  integer :: M,NE
  real(r8) :: dWoody
  real(r8) :: LeafElm2Litr(NumPlantChemElms),WoodElm2Litr(NumPlantChemElms)
  real(r8) :: FineElm2Litr(NumPlantChemElms),StdeadElm2Litr(NumPlantChemElms)  
!     begin_execution
  associate(                                                                    &
    inonstruct                     => pltpar%inonstruct                        ,& !input  :group id of plant nonstructural litter
    k_fine_comp                    => pltpar%k_fine_comp                       ,& !input  :fine litter complex id
    k_woody_comp                   => pltpar%k_woody_comp                      ,& !input  :woody litter complex id
    ifoliar                        => pltpar%ifoliar                           ,& !input  :group id of plant foliar litter
    inonfoliar                     => pltpar%inonfoliar                        ,& !input  :group id of plant non-foliar litter group
    istalk                         => pltpar%istalk                            ,& !input  :group id of plant stalk litter group
    icwood                         => pltpar%icwood                            ,& !input  :group id of coarse woody litter
    iHarvstType_pft                => plt_distb%iHarvstType_pft                ,& !input  :type of harvest,[-]
    FracWoodStalkElmAlloc2Litr     => plt_allom%FracWoodStalkElmAlloc2Litr     ,& !input  :woody element allocation,[-]
    PlantElmAllocMat4Litr          => plt_soilchem%PlantElmAllocMat4Litr       ,& !input  :litter kinetic fraction, [-]
    iPlantTurnoverPattern_pft      => plt_pheno%iPlantTurnoverPattern_pft      ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    iPlant2ndGrothPattern_pft      => plt_pheno%iPlant2ndGrothPattern_pft      ,& !input  :plant expression of secondary growth, [-]                
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft          ,& !input  :plant growth type (vascular, non-vascular),[-]
    StandDeadCompKElms_pft         => plt_biom%StandDeadCompKElms_pft          ,& !inoput :standing dead element fraction, [g d-2]
    LitrfallElms_pvr               => plt_bgcr%LitrfallElms_pvr                ,& !inoput :plant LitrFall element, [g d-2 h-1]
    LitrfalStrutElms_CumYr_pft     => plt_bgcr%LitrfalStrutElms_CumYr_pft      ,& !inoput :total plant element LitrFall, [g d-2 ]
    SurfLitrfalStrutElms_CumYr_pft => plt_bgcr%SurfLitrfalStrutElms_CumYr_pft   & !inoput :total surface LitrFall element, [g d-2]
  )

  call PrintInfo('beg '//subname)
  !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
  !                       ,3=pruning,4=grazing,5=fire,6=herbivory
  !     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
  !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
  !     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+PetolSheth,2=none,3=between 1,2
  !     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
  !
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    !not by fire
    IF(iHarvstType_pft(NZ).NE.iharvtyp_fire)THEN
      LeafElm2Litr   = LeafElmnt2Litr(:)+LeafElmntHarv2Litr(:)
      WoodElm2Litr   = WoodyElmnt2Litr(:)+WoodyElmntHarv2Litr(:)
      FineElm2Litr   = FineNonleafElmnt2Litr(:)+PetolShethElmntHarv2Litr(:)
      StdeadElm2Litr = StandeadElmnt2Litr(:)+StandeadElmntHarv2Litr(:)
      D6375: DO M  = 1, jsken
        DO NE=1,NumPlantChemElms        
          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) &
            +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*AZMAX1(CanopyNonstElm2Litr(NE)) &
            +PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)*LeafElm2Litr(NE) &
            +PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*FineElm2Litr(NE)

          IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. (.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ))))THEN
            LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) &
              +PlantElmAllocMat4Litr(NE,istalk,M,NZ)*AZMAX1(WoodElm2Litr(NE)+StdeadElm2Litr(NE))
          ELSE
            StandDeadCompKElms_pft(NE,M,NZ)=StandDeadCompKElms_pft(NE,M,NZ) &
              +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*AZMAX1(WoodyElmntHarv2Litr(NE)+StandeadElmntHarv2Litr(NE))

            dWoody=AZMAX1(WoodyElmnt2Litr(NE)+StandeadElmnt2Litr(NE))

            LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) &
              +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dWoody*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp)

            LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) &
              +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dWoody*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp)
          ENDIF
        ENDDO
      ENDDO D6375
      !
      !     ABOVE-GROUND LitrFall FROM FIRE
      !
      !     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+PetolSheth,2=none,3=between 1,2
      !     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
      !
    ELSE
      !
      call AbvGrndLiterFallByFire(yearIJ,NZ,CanopyNonstElm2Litr,StandeadElmntOffEcosystem, &
        FineNonleafElmOffEcosystem,LeafElmnt2Litr,LeafElmntOffEcosystem,NonstructElmntOffEcosystem,&
        WoodyElmntOffEcosystem,WoodyElmnt2Litr,StandeadElmnt2Litr,PetolShethElmntHarv2Litr,&
        FineNonleafElmnt2Litr,LeafElmntHarv2Litr,StandeadElmntHarv2Litr,WoodyElmntHarv2Litr)
    ENDIF
    !by grazing or herbivory
  ELSE
    !
    !     ABOVE-GROUND LitrFall FROM GRAZING
    !
    !     TCSNC,TZSNC,TPSNC=cumulative C,N,P LitrFall
    !     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P LitrFall
    !
    DO NE=1,NumPlantChemElms
      LitrfalStrutElms_CumYr_pft(NE,NZ)     = LitrfalStrutElms_CumYr_pft(NE,NZ)+TotalAbgElmnt2Litr(NE)+HarvestAbgElmnt2Litr(NE)
      SurfLitrfalStrutElms_CumYr_pft(NE,NZ) = SurfLitrfalStrutElms_CumYr_pft(NE,NZ)+TotalAbgElmnt2Litr(NE)+HarvestAbgElmnt2Litr(NE)
    ENDDO
  ENDIF
  call PrintInfo('beg '//subname)
  end associate
  end subroutine LiterfallByDisturbance

!----------------------------------------------------------------------------------------------------
  subroutine AbvgBiomRemovalByDisturb(yearIJ,NZ,CanopyNonstElm2Litr,HarvestAbgElmnt2Litr,TotalAbgElmnt2Litr)
  !
  !Description:
  !Loss of above ground biomass due to disturbance
  !
  use InitPlantMod, only : ApplySeedBank
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  

  integer , intent(in)  :: NZ
  real(r8), intent(in)  :: CanopyNonstElm2Litr(NumPlantChemElms)
  real(r8), intent(out) :: HarvestAbgElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: TotalAbgElmnt2Litr(NumPlantChemElms)
  character(len=*), parameter :: subname='AbvgBiomRemovalByDisturb'
  real(r8) :: TotalAbgElmntRemoval(NumPlantChemElms) !total aboveground biomass removal
  real(r8) :: RooTE1(NumPlantChemElms)
  integer :: NE,NR
!     begin_execution
  associate(                                                       &  
    iHarvstType_pft         => plt_distb%iHarvstType_pft          ,& !input  :type of harvest,[-]    
    jHarvstType_pft         => plt_distb%jHarvstType_pft          ,& !input  :flag for stand replacing disturbance,[-]
    MaxNumRootAxes          => pltpar%MaxNumRootAxes              ,& !input  : maximum number root axes,[-]
    PlantPopuLive_pft       => plt_site%PlantPopuLive_pft         ,& !input  :plant population, [d-2]    
    NU                      => plt_site%NU                        ,& !input  :current soil surface layer number, [-]
    NumOfBranches_pft       => plt_morph%NumOfBranches_pft        ,& !input  :number of branches,[-]    
    MaxNumRootLays          => plt_site%MaxNumRootLays            ,& !input  :maximum root layer number,[-]
    NumStructuralRootAxes_pft    => plt_morph%NumStructuralRootAxes_pft     ,& !input: number of structural root axes,[-]  
    RootCRRadius0_rpvr      => plt_morph%RootCRRadius0_rpvr       ,& !inoput: initial radius of roots that may undergo secondary growth, [m]
    RootAge_rpvr            => plt_morph%RootAge_rpvr             ,& !inoput :root age,[h]
    SeasonalNonstElms_pft   => plt_biom%SeasonalNonstElms_pft     ,& !inoput :plant stored nonstructural element at current step, [g d-2]
    EcoHavstElmnt_CumYr_pft => plt_distb%EcoHavstElmnt_CumYr_pft  ,& !inoput :plant element harvest, [g d-2 ]
    NMaxRootBotLayer_pft    => plt_morph%NMaxRootBotLayer_pft     ,& !inoput :maximum soil layer number for all root axes, [-]    
    NGTopRootLayer_pft      => plt_morph%NGTopRootLayer_pft       ,& !inoput  :soil layer at planting depth, [-]    
    EcoHavstElmnt_CumYr_col => plt_distb%EcoHavstElmnt_CumYr_col  ,& !inoput :ecosystem harvest element, [gC d-2]
    PlantElmDistLoss_pft    => plt_distb%PlantElmDistLoss_pft     ,& !inoput :plant loss to disturbance,    [g d-2 h-1]   
    RootSegBaseDepth_raxes  => plt_morph%RootSegBaseDepth_raxes   ,& !inoput : base depth of different root axes, [m]
    Root1stDepz_raxes       => plt_morph%Root1stDepz_raxes        ,& !inoput : root layer depth, [m]    
    Eco_NBP_CumYr_col       => plt_bgcr%Eco_NBP_CumYr_col          & !inoput :total NBP, [g d-2]
  )
  call PrintInfo('beg '//subname)
  !
  !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
  !                       ,3=pruning,4=grazing,5=fire,6=herbivory
  !     jHarvstType_pft=terminate PFT:0=no,1=yes,2=yes,but reseed
  !     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
  !     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
  !     WTRVC,WTRVN,WTRVP=storage C,N,P
  ! harvest/disturbance cause actual removal (which is taken away) and addition to litter.
  ! need double check
  !
  DO NE=1,NumPlantChemElms
    TotalAbgElmntRemoval(NE) = CanopyNonstElmRemoval(NE)+LeafElmntRemoval(NE)+PetoleShethElmntRemoval(NE)+WoodyElmntRemoval(NE)+StandeadElmntRemoval(NE)
    TotalAbgElmnt2Litr(NE)   = CanopyNonstElm2Litr(NE)+LeafElmnt2Litr(NE)+FineNonleafElmnt2Litr(NE)+WoodyElmnt2Litr(NE)+StandeadElmnt2Litr(NE)
    HarvestAbgElmnt2Litr(NE) = LeafElmntHarv2Litr(NE)+PetolShethElmntHarv2Litr(NE)+WoodyElmntHarv2Litr(NE)+StandeadElmntHarv2Litr(NE)
  ENDDO

  IF(jHarvstType_pft(NZ).NE.jharvtyp_tmareseed .and. iHarvstType_pft(NZ).NE.iharvtyp_fire)THEN
   !not do harvest and reseed
    DO NE=1,NumPlantChemElms  
      PlantElmDistLoss_pft(NE,NZ) = PlantElmDistLoss_pft(NE,NZ)+TotalAbgElmntRemoval(NE)-TotalAbgElmnt2Litr(NE)
    ENDDO  
  ENDIF
  !
  IF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
    !
    !     C,N,P REMOVED FROM GRAZING
    !  
    CALL AbvgBiomRemovalByGrazing(yearIJ,NZ,TotalAbgElmnt2Litr,TotalAbgElmntRemoval)
    !
  ELSE
    ! 
    IF(iHarvstType_pft(NZ).EQ.iharvtyp_fire)THEN
      !
      !     C,N,P LOST AS GAS IF FIRE
      !      
      call AbvgBiomRemovalByFire(yearIJ,NZ,TotalAbgElmnt2Litr,TotalAbgElmntRemoval)   

      call SumRootBiome(yearIJ,NZ,RootE1)
      
      !dealing with rejuvenation
      if(PlantPopuLive_pft(NZ).LE.0._r8)then

        call ApplySeedBank(NZ,yearIJ%I,RootE1(ielmc))
        !                
        !if root biomass is significant, allow it to rejuvenate
        if(RootE1(ielmc).LE.0._R8)call InitRootMychorMorphoBio(NZ)
        !        
      endif

    ELSE
      
      IF(jHarvstType_pft(NZ).EQ.jharvtyp_tmareseed)THEN
        !terminate and reseed
        
        DO NE=1,NumPlantChemElms
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+TotalAbgElmntRemoval(NE)-TotalAbgElmnt2Litr(NE)
        ENDDO
        !other
      ELSE
        !harvested
        if(PlantPopuLive_pft(NZ).LE.0._r8)then
          DO NE=1,NumPlantChemElms  
            if(SeasonalNonstElms_pft(NE,NZ).GT.0._r8)then
              PlantElmDistLoss_pft(NE,NZ)  = PlantElmDistLoss_pft(NE,NZ)+SeasonalNonstElms_pft(NE,NZ)
              SeasonalNonstElms_pft(NE,NZ) = 0._r8
            endif
          ENDDO
          NMaxRootBotLayer_pft(NZ) = 0
          NGTopRootLayer_pft(NZ)   = 0
          NumStructuralRootAxes_pft(NZ) = 0
          DO NR=1,MaxNumRootAxes
            RootSegBaseDepth_raxes(NR,NZ)         = 0._r8
            Root1stDepz_raxes(NR,NZ)              = 0._r8
            RootAge_rpvr(NU:MaxNumRootLays,NR,NZ) = 0._r8
            RootCRRadius0_rpvr(:,NR,NZ)           = 0._r8
          ENDDO            
        ENDIF

        DO NE=1,NumPlantChemElms
          EcoHavstElmnt_CumYr_pft(NE,NZ) = EcoHavstElmnt_CumYr_pft(NE,NZ)+TotalAbgElmntRemoval(NE)-TotalAbgElmnt2Litr(NE)
          EcoHavstElmnt_CumYr_col(NE)    = EcoHavstElmnt_CumYr_col(NE)+TotalAbgElmntRemoval(NE)-TotalAbgElmnt2Litr(NE)
        ENDDO
        Eco_NBP_CumYr_col=Eco_NBP_CumYr_col+TotalAbgElmnt2Litr(ielmc)-TotalAbgElmntRemoval(ielmc)
      ENDIF
    ENDIF
  ENDIF
  
  call PrintInfo('end '//subname)
  end associate
  end subroutine AbvgBiomRemovalByDisturb

!----------------------------------------------------------------------------------------------------
  subroutine ApplyDisturbanceBiomRemoval(yearIJ,NZ,CanopyNonstElm2Litr,&
    NonstructElmntOffEcosystem,LeafElmntOffEcosystem,FineNonleafElmOffEcosystem,&
    WoodyElmntOffEcosystem,StandeadElmntOffEcosystem)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), intent(out) :: CanopyNonstElm2Litr(NumPlantChemElms)
  real(r8), intent(out) :: NonstructElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: LeafElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: FineNonleafElmOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: WoodyElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: StandeadElmntOffEcosystem(NumPlantChemElms)
  character(len=*), parameter :: subname='ApplyDisturbanceBiomRemoval'

  real(r8) :: FracAftHVST21,FracAftHVST22,FracAftHVST23,FracAftHVST24

  integer  :: NE
  !   begin_execution
  associate(                                         &
    FracBiomHarvsted => plt_distb%FracBiomHarvsted  ,& !input  :harvest efficiency, [-]
    iHarvstType_pft  => plt_distb%iHarvstType_pft    & !input  :type of harvest,[-]
  )

  call PrintInfo('beg '//subname)
  !     IF NO PLANT C,N,P REMOVED AT HARVEST (ALL RESIDUE RETURNED)
  !
  !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
  !                       ,3=pruning,4=grazing,5=fire,6=herbivory
  !
  FracAftHVST21=1._r8-FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ)
  FracAftHVST22=1._r8-FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ)
  FracAftHVST23=1._r8-FracBiomHarvsted(iHarvst_col,iplthvst_stalk,NZ)
  FracAftHVST24=1._r8-FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ)

  IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
    DO NE=1,NumPlantChemElms
      CanopyNonstElm2Litr(NE)   = CanopyNonstElmRemoval(NE)*FracAftHVST21     !non-structural
      LeafElmnt2Litr(NE)        = LeafElmntRemoval(NE)*FracAftHVST21          !leaf
      FineNonleafElmnt2Litr(NE) = PetoleShethElmntRemoval(NE)*FracAftHVST22   !fine, non-woody, grain
      WoodyElmnt2Litr(NE)       = WoodyElmntRemoval(NE)*FracAftHVST23         !woody
      StandeadElmnt2Litr(NE)    = StandeadElmntRemoval(NE)*FracAftHVST24      !standing dead
    ENDDO
    !
    !     IF ONLY GRAIN C,N,P REMOVED AT HARVEST
    !
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_grain)THEN
    DO NE=1,NumPlantChemElms
      CanopyNonstElm2Litr(NE)   = CanopyNonstElmRemoval(NE)
      LeafElmnt2Litr(NE)        = LeafElmntRemoval(NE)
      FineNonleafElmnt2Litr(NE) = PetoleShethElmntRemoval(NE)-GrainHarvst(NE)*FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ)
      WoodyElmnt2Litr(NE)       = WoodyElmntRemoval(NE)
      StandeadElmnt2Litr(NE)    = StandeadElmntRemoval(NE)
    ENDDO
    !
    !     IF ONLY WOOD C,N,P REMOVED AT HARVEST
    !
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_allabvg)THEN
    DO NE=1,NumPlantChemElms
      CanopyNonstElm2Litr(NE)   = CanopyNonstElmRemoval(NE)*FracAftHVST21
      LeafElmnt2Litr(NE)        = LeafElmntRemoval(NE)*FracAftHVST21
      FineNonleafElmnt2Litr(NE) = PetoleShethElmntRemoval(NE)*FracAftHVST22
      WoodyElmnt2Litr(NE)       = WoodyElmntRemoval(NE)*FracAftHVST23
      StandeadElmnt2Litr(NE)    = StandeadElmntRemoval(NE)*FracAftHVST24
    ENDDO
    !
    !     IF ALL PLANT C,N,P REMOVED AT HARVEST (NO RESIDUE RETURNED)
    !
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
    DO NE=1,NumPlantChemElms
      CanopyNonstElm2Litr(NE)   = CanopyNonstElmRemoval(NE)*FracAftHVST21
      LeafElmnt2Litr(NE)        = LeafElmntRemoval(NE)*FracAftHVST21
      FineNonleafElmnt2Litr(NE) = PetoleShethElmntRemoval(NE)*FracAftHVST22
      WoodyElmnt2Litr(NE)       = WoodyElmntRemoval(NE)*FracAftHVST23
      StandeadElmnt2Litr(NE)    = StandeadElmntRemoval(NE)*FracAftHVST24
    ENDDO
    !
    !     IF PLANT C,N,P REMOVED BY GRAZING
    !
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN

    call ApplyBiomRemovalByGrazing(yearIJ%I,yearIJ%J,NZ,FracAftHVST21,FracAftHVST22,FracAftHVST23,FracAftHVST24,&
      CanopyNonstElmRemoval,LeafElmntRemoval,PetoleShethElmntRemoval,WoodyElmntRemoval,StandeadElmntRemoval,&
      CanopyNonstElm2Litr,LeafElmnt2Litr,FineNonleafElmnt2Litr,WoodyElmnt2Litr,StandeadElmnt2Litr)
    !
    !     REMOVALS BY FIRE
    !
    !     EFIRE=combustion  of N,P relative to C
    !
  ELSEIF(iHarvstType_pft(NZ).EQ.iharvtyp_fire)THEN

    call ShootBiomRemovalByFire(yearIJ%I,yearIJ%J,NZ,&
      FracAftHVST21,FracAftHVST22, FracAftHVST23, FracAftHVST24,&
      StandeadElmntRemoval,CanopyNonstElmRemoval,LeafElmntRemoval,WoodyElmntRemoval,&
      PetoleShethElmntRemoval,CanopyNonstElm2Litr,NonstructElmntOffEcosystem,&
      LeafElmntOffEcosystem,FineNonleafElmOffEcosystem,WoodyElmntOffEcosystem,&
      StandeadElmntOffEcosystem,LeafElmnt2Litr,FineNonleafElmnt2Litr,&
      WoodyElmnt2Litr,StandeadElmnt2Litr)

  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine ApplyDisturbanceBiomRemoval

!----------------------------------------------------------------------------------------------------
  subroutine RemoveBiomByHarvest(yearIJ,NZ)

  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='RemoveBiomByHarvest'
  integer :: L,K,M,NR,N,NB,NBX,NE
  real(r8) :: FracLeftThin
  real(r8) :: XHVST1
  REAL(R8) :: LeafLayerC_brch(NumCanopyLayers1,JP1,JP1)
  real(r8) :: LeafArea2Remove !leaf area left after removal
  real(r8) :: LeafAreaRemoved  
  real(r8) :: APSILT
  real(r8) :: FHVSH,fracL2rm
  real(r8) :: HarvestedLeafC,HarvestedShethC,HarvestedEarC,HarvestedGrainC,GrazedCanopyNonstC
  real(r8) :: HarvestedStalkC,HarvestedStalkRsrvC
  real(r8) :: HarvestedPetoleC
  real(r8) :: GrazedCanopyNoduleC
  real(r8) :: tvegE(NumPlantChemElms),tvegE1(NumPlantChemElms)
  real(r8) :: dLitR,dLitR0,massroot1(NumPlantChemElms),rootE(NumPlantChemElms),rootE1(NumPlantChemElms)

!     begin_execution
  associate(                                                             &
    THIN_pft                   => plt_distb%THIN_pft                    ,& !input  :thinning of plant population, [-]
    iHarvstType_pft            => plt_distb%iHarvstType_pft             ,& !input  :type of harvest,[-]
    jHarvstType_pft            => plt_distb%jHarvstType_pft             ,& !input  :flag for stand replacing disturbance,[-]
    PPI_pft                    => plt_site%PPI_pft                      ,& !input  :initial plant population, [plants d-2]
    NU                         => plt_site%NU                           ,& !input  :current soil surface layer number, [-]
    MaxNumRootLays             => plt_site%MaxNumRootLays               ,& !input  :maximum root layer number,[-]
    SolarNoonHour_col          => plt_site%SolarNoonHour_col            ,& !input  :time of solar noon, [h]
    ZEROS                      => plt_site%ZEROS                        ,& !input  :threshold zero for numerical stability,[-]
    AREA3                      => plt_site%AREA3                        ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    NumStructuralRootAxes_pft       => plt_morph%NumStructuralRootAxes_pft        ,& !input  :number of structural root axes,[-]
    SapwoodBiomassC_brch       => plt_biom%SapwoodBiomassC_brch         ,& !input  :branch live stalk C, [gC d-2]
    StalkStrutElms_brch        => plt_biom%StalkStrutElms_brch          ,& !input  :branch stalk structural element mass, [g d-2]
    CanopyLeafSheathC_brch     => plt_biom%CanopyLeafSheathC_brch       ,& !input  :plant branch leaf + sheath C, [g d-2]
    FracWoodStalkElmAlloc2Litr => plt_allom%FracWoodStalkElmAlloc2Litr  ,& !input  :woody element allocation,[-]
    iPlantPhenolPattern_pft    => plt_pheno%iPlantPhenolPattern_pft     ,& !input  :plant growth habit: annual or perennial,[-]
    PlantElmAllocMat4Litr      => plt_soilchem%PlantElmAllocMat4Litr    ,& !input  :litter kinetic fraction, [-]
    inonstruct                 => pltpar%inonstruct                     ,& !input  :group id of plant nonstructural litter
    k_fine_comp                => pltpar%k_fine_comp                    ,& !input  :fine litter complex id
    k_woody_comp               => pltpar%k_woody_comp                   ,& !input  :woody litter complex id
    LitrfallElms_pvr           => plt_bgcr%LitrfallElms_pvr             ,& !input  :plant LitrFall element, [g d-2 h-1]
    NGTopRootLayer_pft         => plt_morph%NGTopRootLayer_pft          ,& !input  :soil layer at planting depth, [-]
    Myco_pft                   => plt_morph%Myco_pft                    ,& !input  :mycorrhizal type (no or yes),[-]
    CanopyLeafAareZ_col        => plt_morph%CanopyLeafAareZ_col         ,& !input  :total leaf area, [m2 d-2]
    CanopyHeightZ_col          => plt_morph%CanopyHeightZ_col           ,& !input  :canopy layer height, [m]
    NumOfBranches_pft          => plt_morph%NumOfBranches_pft           ,& !input  :number of branches,[-]
    CanopyStalkSurfArea_lbrch  => plt_morph%CanopyStalkSurfArea_lbrch   ,& !input  :plant canopy layer branch stem area, [m2 d-2]
!    SeedBankSize_pft          => plt_morph%SeedBankSize_pft            ,& !output :seed bank size, in terms of number of seeds [d-2]
    CanopyLeafArea_col         => plt_morph%CanopyLeafArea_col          ,& !input  :grid canopy leaf area, [m2 d-2]
    CanopyCutProxy_pft         => plt_distb%CanopyCutProxy_pft          ,& !inoput :harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
    PlantPopuLive_pft          => plt_site%PlantPopuLive_pft            ,& !inoput :plant population, [d-2]
    PlantPopuDead_pft          => plt_site%PlantPopuDead_pft            ,& !inoput :live+standing dead plant population, [d-2]    
    CanopyHeightLive_pft       => plt_morph%CanopyHeightLive_pft        ,& !inoput  :canopy height, [m]        
    PPatSeeding_pft            => plt_site%PPatSeeding_pft              ,& !input  :plant population at seeding, [plants d-2]        
    PPX_pft                    => plt_site%PPX_pft                      ,& !inoput :plant population, [plants m-2]
    SeasonalNonstElms_pft      => plt_biom%SeasonalNonstElms_pft        ,& !inoput :plant stored nonstructural element at current step, [g d-2]
    RootMyco1stElm_raxs        => plt_biom%RootMyco1stElm_raxs          ,& !inoput :root C primary axes, [g d-2]
    CanopySapwoodC_pft         => plt_biom%CanopySapwoodC_pft           ,& !inoput :canopy active stalk C, [g d-2]
    CanopyLeafSheathC_pft      => plt_biom%CanopyLeafSheathC_pft        ,& !inoput :canopy leaf + sheath C, [g d-2]
    StalkStrutElms_pft         => plt_biom%StalkStrutElms_pft           ,& !inoput :canopy stalk structural element mass, [g d-2]
    CanopyStemSurfArea_pft     => plt_morph%CanopyStemSurfArea_pft      ,& !inoput :plant stem area, [m2 d-2]
    ClumpFactor_pft            => plt_morph%ClumpFactor_pft              & !inoput :clumping factor for self-shading in canopy layer, [-]
  )
  !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
  !                       ,3=pruning,4=grazing,5=fire,6=herbivory
  !
  call PrintInfo('beg '//subname)

  IF((iHarvstType_pft(NZ).GE.iharvtyp_none .AND. yearIJ%J.EQ.INT(SolarNoonHour_col)              &     !time to do harvest
    .AND. iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo) &     !neither grazing nor herbivory
    .OR. (iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo))THEN   !or (grazing or herbivory, not required at noon)

    !
    !     ACCUMULATE ALL HARVESTED MATERIAL ABOVE CUTTING HEIGHT
    !     ACCOUNTING FOR HARVEST EFFICIENCY ENTERED IN 'READQ'
    !
    !     jHarvstType_pft=terminate PFT:0=no,1=yes,2=yes,and reseed
    !     PPX,PP=PFT population per m2,grid cell
    !     THIN_pft=thinning:fraction of population removed
    !     CF=clumping factor
    !     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
    !          iHarvstType_pft=3:reduction of clumping factor
    !          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
    !     CanopyLeafArea_col,CanopyLeafAareZ_col=leaf area of combined canopy, canopy layer
    !     LeafAreaRemoved,LeafArea2Remove=leaf area harvested,remaining
    !     ZL=height to bottom of each canopy layer
    !
    !neither grazing nor herbivory
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN

      !not terminate and reseed, including fire
      IF(jHarvstType_pft(NZ).NE.jharvtyp_tmareseed)THEN                
        PPX_pft(NZ)           = PPX_pft(NZ)*(1._r8-THIN_pft(NZ))
        PlantPopuLive_pft(NZ) = PlantPopuLive_pft(NZ)*(1._r8-THIN_pft(NZ))
        PlantPopuDead_pft(NZ)  = PlantPopuDead_pft(NZ)*(1._r8-THIN_pft(NZ))
        !terminate and reseed        
      ELSE
        !modeling seed bank effect
        !SeedBankSize_pft(NZ) = 0.5_r8*PPI_pft(NZ)
        !PPI_pft(NZ)          = AMAX1(1.0_r8,0.5_r8*(PPI_pft(NZ)+CanopySeedNum_pft(NZ)/AREA3(NU)))
        PPX_pft(NZ)           = PPatSeeding_pft(NZ)
        PlantPopuLive_pft(NZ) = PPX_pft(NZ)*AREA3(NU)
        PlantPopuDead_pft(NZ)  = PlantPopuLive_pft(NZ)
      ENDIF

      IF(iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
        ClumpFactor_pft(NZ)=ClumpFactor_pft(NZ)*CanopyCutProxy_pft(NZ)
      ENDIF

      !remove leaf area
      IF(iHarvstType_pft(NZ).LE.iharvtyp_allabvg .AND. CanopyCutProxy_pft(NZ).LT.0.0_r8)THEN
        !leaf area left in column, now CanopyCutProxy_pft is the fraction removed
        LeafArea2Remove = ABS(CanopyCutProxy_pft(NZ))*CanopyLeafArea_col
        LeafAreaRemoved = 0._r8
        CanopyCutProxy_pft(NZ)=0._r8
        !find the cut height due to leaf area removal
        D9875: DO L=NumCanopyLayers1,1,-1
          IF(CanopyHeightZ_col(L).GT.CanopyHeightZ_col(L-1) & !canopy height L is meaningful
            .AND. CanopyLeafAareZ_col(L).GT.ZEROS           & !leaf area in L is meaningful
            .AND. LeafAreaRemoved.LT.LeafArea2Remove)THEN     !removal goal not met
             
            IF(LeafAreaRemoved+CanopyLeafAareZ_col(L).GE.LeafArea2Remove)THEN    
              !reached the layer layer for removal
              !fraction of layer L to be removed
              fracL2rm= (LeafArea2Remove-LeafAreaRemoved)/CanopyLeafAareZ_col(L)
              !cut criterion is converted into height
              CanopyCutProxy_pft(NZ)=CanopyHeightZ_col(L)-fracL2rm*(CanopyHeightZ_col(L)-CanopyHeightZ_col(L-1))
              EXIT
            ENDIF
          ENDIF
          LeafAreaRemoved=LeafAreaRemoved+CanopyLeafAareZ_col(L)
        ENDDO D9875
      ENDIF
      HarvestedLeafC      = 0._r8
      HarvestedShethC     = 0._r8
      HarvestedEarC       = 0._r8
      HarvestedGrainC     = 0._r8
      GrazedCanopyNonstC  = 0._r8
      HarvestedStalkC     = 0._r8
      HarvestedStalkRsrvC = 0._r8
      LeafLayerC_brch     = 0._r8          !it is a filler
    ELSE
      !
      !     GRAZING REMOVAL
      call GrazingPlant(yearIJ%I,yearIJ%J,NZ,HarvestedLeafC,HarvestedShethC,HarvestedEarC,HarvestedGrainC,&
        GrazedCanopyNonstC,HarvestedStalkC,HarvestedStalkRsrvC,HarvestedPetoleC,GrazedCanopyNoduleC,LeafLayerC_brch)
      !
    ENDIF

    !
    !     HARVEST REMOVAL FROM TOP TO BOTTOM OF CANOPY
    call HarvestCanopy(yearIJ,NZ,HarvestedLeafC,LeafLayerC_brch)

    CALL CutPlant(yearIJ,NZ,HarvestedPetoleC,GrazedCanopyNonstC,GrazedCanopyNoduleC,HarvestedShethC,HarvestedGrainC,HarvestedEarC,&
      HarvestedStalkRsrvC,HarvestedStalkC)

    CanopyLeafSheathC_pft(NZ)     = 0._r8
    StalkStrutElms_pft(:,NZ)     = 0._r8
    CanopySapwoodC_pft(NZ)       = 0._r8
    CanopyStemSurfArea_pft(NZ)       = 0._r8

    D9840: DO NB=1,NumOfBranches_pft(NZ)
      CanopyLeafSheathC_pft(NZ)    = CanopyLeafSheathC_pft(NZ)+CanopyLeafSheathC_brch(NB,NZ)
      CanopySapwoodC_pft(NZ)       = CanopySapwoodC_pft(NZ)+SapwoodBiomassC_brch(NB,NZ)
      DO NE=1,NumPlantChemElms
        StalkStrutElms_pft(NE,NZ) = StalkStrutElms_pft(NE,NZ)+StalkStrutElms_brch(NE,NB,NZ)
      ENDDO
      D9830: DO L=1,NumCanopyLayers1
        CanopyStemSurfArea_pft(NZ)=CanopyStemSurfArea_pft(NZ)+CanopyStalkSurfArea_lbrch(L,NB,NZ)
      ENDDO D9830
    ENDDO D9840
    !
    !     ROOT LitrFall FROM HARVESTING OR FIRE
    !
    !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
    !                       ,3=pruning,4=grazing,5=fire,6=herbivory
    !     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
    !     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
    !          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
    !     THETW=soil water concentration
    !     CORGC=SOC concentration
    !     EFIRE=combustion  of N,P relative to C
    !
    !
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo & !not grazing/herbiory
      .AND. iHarvstType_pft(NZ).NE.iharvtyp_allabvg)THEN !including belowgorund
      !
      FracLeftThin=1.0_r8-THIN_pft(NZ)
      DO NR=1,NumStructuralRootAxes_pft(NZ)        
        DO NE=1,NumPlantChemElms
          RootMyco1stElm_raxs(NE,NR,NZ)=RootMyco1stElm_raxs(NE,NR,NZ)*FracLeftThin
        ENDDO        
      ENDDO
      !

      D3985: DO N=1,Myco_pft(NZ)
        D3980: DO L=NU,MaxNumRootLays
          !          
          CALL RootRemovalLbyFire(yearIJ,N,L,NZ,FracLeftThin,XHVST1)

          call HarvstUpdateRootStateL(yearIJ,N,L,NZ,FracLeftThin,XHVST1, &
            PopulationThinning=iHarvstType_pft(NZ).NE.iharvtyp_fire .and. &
              jHarvstType_pft(NZ).NE.jharvtyp_tmareseed)

        ENDDO D3980
      ENDDO D3985
      !
      !     STORAGE LitrFall AND STATE VARIABLES DURING HARVESTING
      !
      !     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
      !     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
      !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
      !     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
      !     WTRVC,WTRVN,WTRVP=storage C,N,P
      !
      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial)THEN
        XHVST1=1._r8-FracLeftThin    
        DO NE=1,NumPlantChemElms
          D3400: DO M=1,jsken
            !
            dLitR=XHVST1*PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*AZMAX1(SeasonalNonstElms_pft(NE,NZ))
            !
            LitrfallElms_pvr(NE,M,k_woody_comp,NGTopRootLayer_pft(NZ),NZ)=&
              LitrfallElms_pvr(NE,M,k_woody_comp,NGTopRootLayer_pft(NZ),NZ) &
              +dLitR*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp)

            LitrfallElms_pvr(NE,M,k_fine_comp,NGTopRootLayer_pft(NZ),NZ)=&
              LitrfallElms_pvr(NE,M,k_fine_comp,NGTopRootLayer_pft(NZ),NZ) &
              +dLitR*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp)
          ENDDO D3400
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)*FracLeftThin
        ENDDO
      ENDIF
    ENDIF
  ENDIF
  
  call PrintInfo('end '//subname)
  end associate
  end subroutine RemoveBiomByHarvest

!----------------------------------------------------------------------------------------------------
  subroutine ResetCutBranch(I,J,NZ,NB)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ,NB
  character(len=*), parameter :: subname='ResetCutBranch'

  logical :: nonEvergreenChk
  integer :: NBX,M
 
  associate(                                                                           &
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft               ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    HourReq4LeafOff_brch              => plt_pheno%HourReq4LeafOff_brch               ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    MatureGroup_pft                   => plt_pheno%MatureGroup_pft                    ,& !input  :acclimated plant maturity group, [-]
    ShootNodeNum_brch                 => plt_morph%ShootNodeNum_brch                  ,& !input  :shoot node number, [-]
    NumOfBranches_pft                 => plt_morph%NumOfBranches_pft                  ,& !input  :number of branches,[-]
    MainBranchNum_pft                 => plt_morph%MainBranchNum_pft                  ,& !input  :number of main branch,[-]
    Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch                 ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch                ,& !inoput :plant growth stage, [-]
    MatureGroup_brch                  => plt_pheno%MatureGroup_brch                   ,& !output :plant maturity group, [-]
    ShootNodeNumAtAnthesis_brch       => plt_morph%ShootNodeNumAtAnthesis_brch        ,& !output :shoot node number at anthesis, [-]
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch      ,& !output :normalized node number during vegetative growth stages, [-]
    ShootNodeNumAtInitFloral_brch     => plt_morph%ShootNodeNumAtInitFloral_brch      ,& !output :shoot node number at floral initiation, [-]
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch  ,& !output :normalized node number during reproductive growth stages, [-]
    dReproNodeNumNormByMatG_brch      => plt_pheno%dReproNodeNumNormByMatG_brch ,& !output :current hourly reproductive development increment, [h-1]
    doInitLeafOut_brch                => plt_pheno%doInitLeafOut_brch                 ,& !output :branch phenology flag, [-]
    HourFailGrainFill_brch            => plt_pheno%HourFailGrainFill_brch             ,& !output :flag to detect physiological maturity from grain fill, [-]
    LeafNumberAtFloralInit_brch       => plt_pheno%LeafNumberAtFloralInit_brch         & !output :leaf number at floral initiation, [-]
  )
  call PrintInfo('beg '//subname)
  !     ZC=canopy height
  !     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
  !     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
  !     iPlantCalendar_brch(ipltcal_Emerge,=emergence date
  !     GROUP=node number required for floral initiation
  !     ShootNodeNumAtInitFloral_brch=node number at floral initiation
  !     ShootNodeNumAtAnthesis_brch=node number at flowering
  !     VSTGX=leaf number on date of floral initiation
  !     TotalNodeNumNormByMatgrp_brch=total change in vegve node number normalized for maturity group
  !     TotReproNodeNumNormByMatrgrp_brch=total change in reprve node number normalized for maturity group
  !     HourFailGrainFill_brch=number of hours with no grain fill
  !     doInitLeafOut_brch=flag for initializing leafout
  !

  nonEvergreenChk=iPlantPhenolType_pft(NZ).NE.iphenotyp_evgreen &
    .AND. Hours4LeafOff_brch(NB,NZ).LE.FracHour4LeafoffRemob(iPlantPhenolType_pft(NZ))*HourReq4LeafOff_brch(NB,NZ)

  IF(nonEvergreenChk .OR. (iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen .AND. iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).NE.0))THEN
    MatureGroup_brch(NB,NZ)                      = MatureGroup_pft(NZ)
    ShootNodeNumAtInitFloral_brch(NB,NZ)         = ShootNodeNum_brch(NB,NZ)
    ShootNodeNumAtAnthesis_brch(NB,NZ)             = 0._r8
    LeafNumberAtFloralInit_brch(NB,NZ)           = 0._r8
    TotalNodeNumNormByMatgrp_brch(NB,NZ)         = 0._r8
    TotReproNodeNumNormByMatrgrp_brch(NB,NZ)     = 0._r8
    dReproNodeNumNormByMatG_brch(NB,NZ)          = 0._r8
    HourFailGrainFill_brch(NB,NZ)                = 0._r8
    iPlantCalendar_brch(ipltcal_Emerge,NB,NZ)    = I
    iPlantCalendar_brch(2:NumGrowthStages,NB,NZ) = 0
    doInitLeafOut_brch(NB,NZ)                    = iTrue

    IF(NB.EQ.MainBranchNum_pft(NZ))THEN
      D3010: DO NBX=1,NumOfBranches_pft(NZ)
        IF(NBX.NE.MainBranchNum_pft(NZ))THEN
          MatureGroup_brch(NBX,NZ)                   = MatureGroup_pft(NZ)
          ShootNodeNumAtInitFloral_brch(NBX,NZ)      = ShootNodeNum_brch(NBX,NZ)
          ShootNodeNumAtAnthesis_brch(NBX,NZ)          = 0._r8
          LeafNumberAtFloralInit_brch(NBX,NZ)        = 0._r8
          TotalNodeNumNormByMatgrp_brch(NBX,NZ)      = 0._r8
          TotReproNodeNumNormByMatrgrp_brch(NBX,NZ)  = 0._r8
          dReproNodeNumNormByMatG_brch(NBX,NZ)       = 0._r8
          HourFailGrainFill_brch(NBX,NZ)             = 0._r8
          iPlantCalendar_brch(ipltcal_Emerge,NBX,NZ) = I
          D3015: DO M=2,NumGrowthStages
            iPlantCalendar_brch(M,NBX,NZ)=0
          ENDDO D3015
          doInitLeafOut_brch(NBX,NZ)=0
        ENDIF
      ENDDO D3010
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine ResetCutBranch

!----------------------------------------------------------------------------------------------------
  subroutine CutBranchStalk(yearIJ,NB,NZ,BranchLength,HarvestedStalkC,HarvestedStalkRsrvC)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NB,NZ
  real(r8), intent(in) :: BranchLength     !branch stalk length subject to removal
  real(r8), intent(in) :: HarvestedStalkC  !harvested stalk C, [gC h-1]
  real(r8), intent(in) :: HarvestedStalkRsrvC   !harvested reserve C, [gC h-1]

  character(len=*), parameter :: subname='CutBranchStalk'
  integer :: NE,K  
  real(r8) :: FractionStalkMassLeft,FracHeightLeft
  real(r8) :: FracNodeLen4Cut       !fraction of node length subject to cut
  real(r8) :: FHVSH,FHVSETS         !fraction of organ left after harvest
  

  associate(                                                         &
    CanopyCutProxy_pft       => plt_distb%CanopyCutProxy_pft        ,& !input  :harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
    THIN_pft                 => plt_distb%THIN_pft                  ,& !input  :thinning of plant population, [-]
    FracBiomHarvsted         => plt_distb%FracBiomHarvsted          ,& !input  :harvest efficiency, [-]
    StalkRsrvElms_pft        => plt_biom%StalkRsrvElms_pft          ,& !input  :canopy reserve element mass, [g d-2]
    ZERO4Groth_pft           => plt_biom%ZERO4Groth_pft             ,& !input  :threshold zero for plang growth calculation, [-]
    ZERO4LeafVar_pft         => plt_biom%ZERO4LeafVar_pft           ,& !input  :threshold zero for leaf calculation, [-]
    ZERO                     => plt_site%ZERO                       ,& !input  :threshold zero for numerical stability, [-]
    StalkStrutElms_pft       => plt_biom%StalkStrutElms_pft         ,& !input  :canopy stalk structural element mass, [g d-2]
    iHarvstType_pft          => plt_distb%iHarvstType_pft           ,& !input  :type of harvest,[-]
    StalkRsrvElms_brch       => plt_biom%StalkRsrvElms_brch         ,& !inoput :branch reserve element mass, [g d-2]
    StalkNodeVertLength_brch => plt_morph%StalkNodeVertLength_brch  ,& !inoput :internode height, [m]
    SapwoodBiomassC_brch     => plt_biom%SapwoodBiomassC_brch       ,& !inoput :branch live stalk C, [gC d-2]
    StalkStrutElms_brch      => plt_biom%StalkStrutElms_brch        ,& !inoput :branch stalk structural element mass, [g d-2]
    StalkNodeHeight_brch     => plt_morph%StalkNodeHeight_brch      ,& !inoput :internode height, [m]
    SenecStalkStrutElms_brch => plt_biom%SenecStalkStrutElms_brch   ,& !inoput :branch stalk structural element, [g d-2]
    StructInternodeElms_brch => plt_biom%StructInternodeElms_brch    & !inoput :internode C, [g d-2]
  )
  call PrintInfo('beg '//subname)
  !
  !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
  !                       ,3=pruning,4=grazing,5=fire,6=herbivory
  !     BranchLength=internode length subject to removal
  !     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
  !          iHarvstType_pft=3:reduction of clumping factor
  !          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
  !     FracHeightLeft=fraction of canopy layer height not harvested
  !     FractionStalkMassLeft=fraction of canopy layer mass not harvested
  !     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
  !          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
  !
  !neither grazing nor herbivory
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(BranchLength.GT.ZERO)THEN
      !not pruning
      IF(iHarvstType_pft(NZ).NE.iharvtyp_pruning)THEN
        FracHeightLeft=AZMAX1(AMIN1(1.0_r8,CanopyCutProxy_pft(NZ)/BranchLength)) !above CanopyCutProxy_pft is removed
      ELSE
        FracHeightLeft=0._r8
      ENDIF

      IF(isclose(THIN_pft(NZ),0._r8))THEN
        !fleft_mass:=1-fharvested_height*fharvested_wood
        FractionStalkMassLeft = AZMAX1(1._r8-(1._r8-FracHeightLeft)*FracBiomHarvsted(iHarvst_pft,iplthvst_stalk,NZ))
        FHVSH                 = FractionStalkMassLeft
      ELSE
        FractionStalkMassLeft=AZMAX1(1._r8-THIN_pft(NZ))
        IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
          FHVSH=1.0_r8-(1._r8-FracHeightLeft)*FracBiomHarvsted(iHarvst_pft,iplthvst_stalk,NZ)*THIN_pft(NZ)
        ELSE
          FHVSH=FractionStalkMassLeft
        ENDIF
      ENDIF
    ELSE
      FractionStalkMassLeft = 1.0_r8
      FHVSH                 = 1.0_r8
    ENDIF
  ELSE
    !grazing or herbivory
    IF(StalkStrutElms_pft(ielmc,NZ).GT.ZERO4LeafVar_pft(NZ))THEN
      FractionStalkMassLeft = AZMAX1(AMIN1(1.0_r8,1._r8-HarvestedStalkC/StalkStrutElms_pft(ielmc,NZ)))
      FHVSH                 = FractionStalkMassLeft
    ELSE
      FractionStalkMassLeft = 1.0_r8
      FHVSH                 = 1.0_r8
    ENDIF
  ENDIF

  !
  !     HARVESTED STALK C,N,P
  !
  !     WoodyElmntRemoval=harvested stalk C,N,P
  !     WoodyElmntHarv2Litr=harvested stalk C,N,P to litter
  !     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in harvested stalk
  !
  DO NE=1,NumPlantChemElms
    WoodyElmntRemoval(NE)   = WoodyElmntRemoval(NE)+(1._r8-FHVSH)*StalkStrutElms_brch(NE,NB,NZ)
    WoodyElmntHarv2Litr(NE) = WoodyElmntHarv2Litr(NE)+(FHVSH-FractionStalkMassLeft)*StalkStrutElms_brch(NE,NB,NZ)
    !
    !     REMAINING STALK C,N,P
    !
    !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in harvested stalk
    !
    StalkStrutElms_brch(NE,NB,NZ)      = FractionStalkMassLeft*StalkStrutElms_brch(NE,NB,NZ)
    SenecStalkStrutElms_brch(NE,NB,NZ) = FractionStalkMassLeft*SenecStalkStrutElms_brch(NE,NB,NZ)
  ENDDO

  SapwoodBiomassC_brch(NB,NZ)=FractionStalkMassLeft*SapwoodBiomassC_brch(NB,NZ)
  !
  !     CUT STALK NODES
  !
  !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
  !                       ,3=pruning,4=grazing,5=fire,6=herbivory
  !     StalkNodeVertLength_brch,StalkNodeHeight_brch=stalk height,stalk internode length
  !     FracNodeLen4Cut=fraction of internode length harvested
  !     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed, 
  !          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
  !     FracBiomHarvsted(iHarvst_pft,1:4=fraction of leaf,non-foliar,woody, standing dead removed from PFT
  !     WTSTK=stalk C mass
  !     StructInternodeElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
  !
  D9820: DO K=MaxNodesPerBranch1,0,-1
    !neither grazing nor herbivory
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      IF(StalkNodeVertLength_brch(K,NB,NZ).GT.ZERO)THEN      
        IF(iHarvstType_pft(NZ).NE.iharvtyp_pruning)THEN
          FracNodeLen4Cut=AZMAX1(AMIN1(1.0_r8,(StalkNodeHeight_brch(K,NB,NZ)-CanopyCutProxy_pft(NZ))/StalkNodeVertLength_brch(K,NB,NZ)))
        ELSE
          FracNodeLen4Cut=0._r8
        ENDIF
        
        IF(isclose(THIN_pft(NZ),0._r8))THEN
          FHVSETS=AZMAX1(1._r8-FracNodeLen4Cut*FracBiomHarvsted(iHarvst_pft,iplthvst_stalk,NZ))
        ELSE
          FHVSETS=AZMAX1(1._r8-THIN_pft(NZ))
        ENDIF
      ELSE
        FHVSETS=1.0_r8
      ENDIF
    ELSE
      IF(StalkStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
        FHVSETS=AZMAX1(AMIN1(1.0_r8,1._r8-HarvestedStalkC/StalkStrutElms_pft(ielmc,NZ)))
      ELSE
        FHVSETS=1.0_r8
      ENDIF
    ENDIF

    DO NE=1,NumPlantChemElms
      StructInternodeElms_brch(NE,K,NB,NZ)=FHVSETS*StructInternodeElms_brch(NE,K,NB,NZ)
    ENDDO
    IF(iHarvstType_pft(NZ).LE.iharvtyp_allabvg.AND.isclose(THIN_pft(NZ),0._r8))THEN
      StalkNodeVertLength_brch(K,NB,NZ) = FHVSETS*StalkNodeVertLength_brch(K,NB,NZ)
      StalkNodeHeight_brch(K,NB,NZ)     = AMIN1(StalkNodeHeight_brch(K,NB,NZ),CanopyCutProxy_pft(NZ))
    ENDIF

  ENDDO D9820
  !
  !     CUT STALK RESERVES
  !
  !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
  !                       ,3=pruning,4=grazing,5=fire,6=herbivory
  !     WTSTKB=C mass remaining in harvested stalk
  !     WTRSV=stalk reserve C mass
  !     HarvestedStalkRsrvC=remaining stalk reserve C mass
  !     FractionStalkMassLeft=fraction of reserve mass not harvested
  ! double check below
  !neither grazing nor herbivory
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
    IF(StalkStrutElms_brch(ielmc,NB,NZ).LE.ZERO4Groth_pft(NZ))THEN
      FractionStalkMassLeft = 0._r8
      FHVSH           = 0._r8
    ENDIF
    !grazing or herbivory
  ELSE
    !
    IF(StalkRsrvElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FractionStalkMassLeft = AZMAX1(AMIN1(1.0_r8,1._r8-HarvestedStalkRsrvC/StalkRsrvElms_pft(ielmc,NZ)))
      FHVSH                 = FractionStalkMassLeft
    ELSE
      FractionStalkMassLeft = 0._r8
      FHVSH                 = 0._r8
    ENDIF
  ENDIF
  !
  !     HARVESTED STALK RESERVE C,N,P
  !
  !     WoodyElmntRemoval=harvested stalk C,N,P
  !     WoodyElmntHarv2Litr=harvested stalk C,N,P to litter
  !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
  !
  DO NE=1,NumPlantChemElms
    WoodyElmntRemoval(NE)   = WoodyElmntRemoval(NE)+(1._r8-FHVSH)*StalkRsrvElms_brch(NE,NB,NZ)
    WoodyElmntHarv2Litr(NE) = WoodyElmntHarv2Litr(NE)+(FHVSH-FractionStalkMassLeft)*StalkRsrvElms_brch(NE,NB,NZ)
    !
    !     REMAINING STALK RESERVE C,N,P IF STALK REMAINING
    !
    StalkRsrvElms_brch(NE,NB,NZ)=FractionStalkMassLeft*StalkRsrvElms_brch(NE,NB,NZ)
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine CutBranchStalk

!----------------------------------------------------------------------------------------------------
  subroutine CutBranchReprodOrgans(I,J,NB,NZ,BranchLength,HarvestedShethC,HarvestedGrainC,HarvestedEarC)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: BranchLength
  real(r8), intent(in) :: HarvestedShethC
  real(r8), intent(in) :: HarvestedGrainC
  real(r8), intent(in) :: HarvestedEarC

  character(len=*), parameter :: subname='CutBranchReprodOrgans'
  integer :: NE
  real(r8) :: FracGrainNotHvsted,FracShethGrainNotHvsted
  real(r8) :: FracHuskNotHvsted,FracEarNotHvsted,FracShethHuskNotHvsted,FracShethNotHvsted
  associate(                                                       &
    CanopyCutProxy_pft      => plt_distb%CanopyCutProxy_pft       ,& !input  :harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
    FracBiomHarvsted        => plt_distb%FracBiomHarvsted         ,& !input  :harvest efficiency, [-]
    GrainStrutElms_pft      => plt_biom%GrainStrutElms_pft        ,& !input  :canopy grain structural element, [g d-2]
    THIN_pft                => plt_distb%THIN_pft                 ,& !input  :thinning of plant population, [-]
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft            ,& !input  :threshold zero for plang growth calculation, [-]
    HuskStrutElms_pft       => plt_biom%HuskStrutElms_pft         ,& !input  :canopy husk structural element mass, [g d-2]
    EarStrutElms_pft        => plt_biom%EarStrutElms_pft          ,& !input  :canopy ear structural element, [g d-2]
    iHarvstType_pft         => plt_distb%iHarvstType_pft          ,& !input  :type of harvest,[-]
    GrainStrutElms_brch     => plt_biom%GrainStrutElms_brch       ,& !inoput :branch grain structural element mass, [g d-2]
    EarStrutElms_brch       => plt_biom%EarStrutElms_brch         ,& !inoput :branch ear structural chemical element mass, [g d-2]
    PotentialSeedSites_brch => plt_morph%PotentialSeedSites_brch  ,& !inoput :branch potential grain number, [d-2]
    SetNumberSeeds_brch       => plt_morph%SetNumberSeeds_brch        ,& !inoput :branch grain number, [d-2]
    SingleGrainMeanBiomC_brch => plt_allom%SingleGrainMeanBiomC_brch  ,& !inoput :potential carbon mass per grain, [gC seed-1]
    HuskStrutElms_brch      => plt_biom%HuskStrutElms_brch         & !inoput :branch husk structural element mass, [g d-2]
  )

  call PrintInfo('beg '//subname)
  ! 
  !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
  !                       ,3=pruning,4=grazing,5=fire,6=herbivory
  !     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
  !          iHarvstType_pft=3:reduction of clumping factor
  !          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
  !     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
  !          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
  !     FracGrainNotHvsted,FracHuskNotHvsted,FracEarNotHvsted=fraction of grain,husk,ear mass not harvested
  !     FracBiomHarvsted(iHarvst_pft,1:4=fraction of leaf,non-foliar,woody, standing dead removed from PFT
  !     WTHSK,WTEAR,WTGR=PFT husk,ear,grain C mass
  !
  IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN  !neither grazing nor herbivory
    IF(CanopyCutProxy_pft(NZ).LT.BranchLength .OR. iHarvstType_pft(NZ).EQ.iharvtyp_grain &  !
      .OR. iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
      IF(isclose(THIN_pft(NZ),0._r8))THEN
        FracGrainNotHvsted      = 1.0_r8-FracBiomHarvsted(iHarvst_pft,iplthvst_finenonleaf,NZ)
        FracShethGrainNotHvsted = FracGrainNotHvsted
      ELSE
        FracGrainNotHvsted      = 1.0_r8-THIN_pft(NZ)
        FracShethGrainNotHvsted = 1.0_r8-FracBiomHarvsted(iHarvst_pft,iplthvst_finenonleaf,NZ)*THIN_pft(NZ)
      ENDIF
    ELSE
      FracGrainNotHvsted      = 1.0_r8-THIN_pft(NZ)
      FracShethGrainNotHvsted = FracGrainNotHvsted
    ENDIF
    FracHuskNotHvsted      = FracGrainNotHvsted
    FracEarNotHvsted       = FracGrainNotHvsted
    FracShethHuskNotHvsted = FracShethGrainNotHvsted
    FracShethNotHvsted     = FracShethGrainNotHvsted

    !grazing or herbivory
  ELSE

    IF(HuskStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FracHuskNotHvsted      = AZMAX1(AMIN1(1.0_r8,1._r8-HarvestedShethC/HuskStrutElms_pft(ielmc,NZ)))
      FracShethHuskNotHvsted = FracHuskNotHvsted
    ELSE
      FracHuskNotHvsted      = 1.0_r8
      FracShethHuskNotHvsted = 1.0_r8
    ENDIF

    IF(EarStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FracEarNotHvsted   = AZMAX1(AMIN1(1.0_r8,1._r8-HarvestedEarC/EarStrutElms_pft(ielmc,NZ)))
      FracShethNotHvsted = FracEarNotHvsted
    ELSE
      FracEarNotHvsted   = 1.0_r8
      FracShethNotHvsted = 1.0_r8
    ENDIF

    IF(GrainStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FracGrainNotHvsted      = AZMAX1(AMIN1(1.0_r8,1._r8-HarvestedGrainC/GrainStrutElms_pft(ielmc,NZ)))
      FracShethGrainNotHvsted = FracGrainNotHvsted
    ELSE
      FracGrainNotHvsted      = 1.0_r8
      FracShethGrainNotHvsted = 1.0_r8
    ENDIF
  ENDIF
  !
  !     HARVESTED REPRODUCTIVE C,N,P
  !
  !     PetoleShethElmntRemoval=reproductive C,N,P removed
  !     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
  !     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
  !     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
  !     GrainHarvst()=grain harvested
  !
  DO NE=1,NumPlantChemElms
    PetoleShethElmntRemoval(NE)=PetoleShethElmntRemoval(NE)+(1._r8-FracShethHuskNotHvsted)*HuskStrutElms_brch(NE,NB,NZ)&
      +(1._r8-FracShethNotHvsted)*EarStrutElms_brch(NE,NB,NZ)+(1._r8-FracShethGrainNotHvsted)*GrainStrutElms_brch(NE,NB,NZ)

    PetolShethElmntHarv2Litr(NE)=PetolShethElmntHarv2Litr(NE)+(FracShethHuskNotHvsted-FracHuskNotHvsted)*HuskStrutElms_brch(NE,NB,NZ) &
      +(FracShethNotHvsted-FracEarNotHvsted)*EarStrutElms_brch(NE,NB,NZ) &
      +(FracShethGrainNotHvsted-FracGrainNotHvsted)*GrainStrutElms_brch(NE,NB,NZ)

    GrainHarvst(NE)=GrainHarvst(NE)+(1._r8-FracGrainNotHvsted)*GrainStrutElms_brch(NE,NB,NZ)

    !
    !     REMAINING REPRODUCTIVE C,N,P
    !
    HuskStrutElms_brch(NE,NB,NZ)  = FracHuskNotHvsted*HuskStrutElms_brch(NE,NB,NZ)
    EarStrutElms_brch(NE,NB,NZ)   = FracEarNotHvsted*EarStrutElms_brch(NE,NB,NZ)
    GrainStrutElms_brch(NE,NB,NZ) = FracGrainNotHvsted*GrainStrutElms_brch(NE,NB,NZ)
  ENDDO

  PotentialSeedSites_brch(NB,NZ) = FracGrainNotHvsted*PotentialSeedSites_brch(NB,NZ)
  SetNumberSeeds_brch(NB,NZ)         = FracGrainNotHvsted*SetNumberSeeds_brch(NB,NZ)
  ! Partial removal changes grain number, not the size potential of surviving grains.
  IF(FracGrainNotHvsted.LE.0._r8) SingleGrainMeanBiomC_brch(NB,NZ) = 0._r8
  call PrintInfo('end '//subname)
  end associate
  END subroutine CutBranchReprodOrgans

!----------------------------------------------------------------------------------------------------
  subroutine CutBranchNonstructural(I,J,NB,NZ,LeafCafCut_brch,PetolShethCAfHvst_brch,LeafCB4Cut_brch,&
    PetolShethCB4Hvst_brch,GrazedCanopyNonstC,GrazedCanopyNoduleC)

  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: LeafCafCut_brch,PetolShethCAfHvst_brch,LeafCB4Cut_brch,PetolShethCB4Hvst_brch
  real(r8), intent(in) :: GrazedCanopyNoduleC
  real(r8), intent(in) :: GrazedCanopyNonstC

  character(len=*), parameter :: subname='CutBranchNonstructural'
  real(r8) :: CanopyNonstElmCopy_brch(NumPlantChemElms)  
  real(r8) :: CanopyNodulNonstElmCopy_brch(NumPlantChemElms)
  real(r8) :: CanopyNonstElmAfhvst_brch(NumPlantChemElms)
  real(r8) :: CanopyNodulNonstElmAfhvst(NumPlantChemElms)  
  real(r8) :: CanopyNodulStrutElmAfhvst(NumPlantChemElms)
  real(r8) :: FracNonstLeft,FracNonstRemoved
  integer  :: K,NE
  real(r8) :: FrcLeafMassLeft
  
  associate(                                                             &
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft               ,& !input  :threshold zero for plang growth calculation, [-]
    iPlantPhotosynsType_pft    => plt_photo%iPlantPhotosynsType_pft     ,& !input  :plant photosynthetic type (C3 or C4),[-]
    iHarvstType_pft            => plt_distb%iHarvstType_pft             ,& !input  :type of harvest,[-]
    CanopyNonstElms_brch       => plt_biom%CanopyNonstElms_brch         ,& !inoput :branch nonstructural element, [g d-2]
    CanopyNodulStrutElms_brch  => plt_biom%CanopyNodulStrutElms_brch    ,& !inoput :branch nodule structural element, [g d-2]
    CPOOL3_node                => plt_photo%CPOOL3_node                 ,& !inoput :minimum sink strength for nonstructural C transfer, [g d-2]
    CPOOL4_node                => plt_photo%CPOOL4_node                 ,& !inoput :leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
    CMassHCO3BundleSheath_node => plt_photo%CMassHCO3BundleSheath_node  ,& !inoput :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CMassCO2BundleSheath_node  => plt_photo%CMassCO2BundleSheath_node   ,& !inoput :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CanopyNodulNonstElms_brch  => plt_biom%CanopyNodulNonstElms_brch     & !inoput :branch nodule nonstructural element, [g d-2]
  )
  call PrintInfo('beg '//subname)
  !
  !     CUT NON-STRUCTURAL C,N,P IN HARVESTED BRANCHES
  !
  !     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
  !     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
  !     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
  !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
  !                       ,3=pruning,4=grazing,5=fire,6=herbivory
  !     FrcLeafMassLeft=fraction of leaf+PetolSheth node mass not harvested
  !     CanopyNonstElmAfhvst_brch(),=branch non-structural C,N,P mass after harvest
  !     CanopyNodulNonstElmAfhvst=nonstructural C,N,P in bacteria after harvest
  !     CanopyNodulStrutElmAfhvst=bacterial C,N,P mass after harvest
  !     WTLS,CanopyLeafSheathC_brch=total,branch PFT leaf+PetolSheth C mass
  !     WHVSC*=nonstructural C removed
  ! make a copy of canopy variables preparing for harvest
  DO NE=1,NumPlantChemElms
    CanopyNonstElmCopy_brch(NE)      = AZMAX1(CanopyNonstElms_brch(NE,NB,NZ))
    CanopyNodulNonstElmCopy_brch(NE) = AZMAX1(CanopyNodulNonstElms_brch(NE,NB,NZ))
  ENDDO

  IF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
    !Fallback for depleted C pools: remove reserves with the actual host tissue.
    FrcLeafMassLeft=1._r8
    IF(LeafCB4Cut_brch+PetolShethCB4Hvst_brch.GT.ZERO4Groth_pft(NZ))THEN
      FrcLeafMassLeft=AZMAX1(AMIN1(1._r8,(LeafCafCut_brch+PetolShethCAfHvst_brch) &
        /(LeafCB4Cut_brch+PetolShethCB4Hvst_brch)))
    ENDIF
    call CutBranchNonstalByGrazing(I,J,NB,NZ,GrazedCanopyNonstC,CanopyNonstElmCopy_brch,&
      GrazedCanopyNoduleC,CanopyNodulNonstElmCopy_brch,CanopyNonstElmAfhvst_brch,CanopyNodulNonstElmAfhvst,&
      CanopyNodulStrutElmAfhvst,FrcLeafMassLeft)
  ELSE
    IF(LeafCB4Cut_brch+PetolShethCB4Hvst_brch.GT.ZERO4Groth_pft(NZ))THEN
      FrcLeafMassLeft=AZMAX1(AMIN1(1.0_r8,(LeafCafCut_brch+PetolShethCAfHvst_brch)/(LeafCB4Cut_brch+PetolShethCB4Hvst_brch)))
      DO NE=1,NumPlantChemElms
        CanopyNonstElmAfhvst_brch(NE) = CanopyNonstElmCopy_brch(NE)*FrcLeafMassLeft
        CanopyNodulNonstElmAfhvst(NE) = CanopyNodulNonstElmCopy_brch(NE)*FrcLeafMassLeft
        CanopyNodulStrutElmAfhvst(NE) = CanopyNodulStrutElms_brch(NE,NB,NZ)*FrcLeafMassLeft
      ENDDO
    ELSE
      CanopyNonstElmAfhvst_brch(:)      = 0._r8
      CanopyNodulNonstElmAfhvst(:) = 0._r8
      CanopyNodulStrutElmAfhvst(:) = 0._r8
    ENDIF
  ENDIF
  !
  !     HARVESTED NON-STRUCTURAL C, N, P
  !
  !     CanopyNonstElmRemoval()=nonstructural C,N,P removed
  !

  DO NE=1,NumPlantChemElms
    CanopyNonstElmRemoval(NE)=CanopyNonstElmRemoval(NE)+CanopyNonstElmCopy_brch(NE)-CanopyNonstElmAfhvst_brch(NE)+&
      CanopyNodulNonstElmCopy_brch(NE)-CanopyNodulNonstElmAfhvst(NE)
    CanopyNonstElmRemoval(NE)=CanopyNonstElmRemoval(NE)+CanopyNodulStrutElms_brch(NE,NB,NZ)-CanopyNodulStrutElmAfhvst(NE)
    !
    ! REMAINING NON-STRUCTURAL C, N, P
    !
    ! CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
    ! CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
    ! WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
    !
    CanopyNonstElms_brch(NE,NB,NZ)      = CanopyNonstElmAfhvst_brch(NE)
    CanopyNodulNonstElms_brch(NE,NB,NZ) = CanopyNodulNonstElmAfhvst(NE)
    CanopyNodulStrutElms_brch(NE,NB,NZ) = CanopyNodulStrutElmAfhvst(NE)
  ENDDO  
  !
  !     REMOVE C4 NON-STRUCTURAL C
  !
  !     iPlantPhotosynsType_pft=photosynthesis type:3=C3,4=C4 from PFT file
  !     FracNonstLeft=fraction of nonstructural mass not harvested
  !     CanopyNonstElmAfhvst_brch(ielmc)=branch non-structural C mass after harvest
  !     CanopyNonstElmRemoval()=nonstructural C,N,P removed
  !     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
  !     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
  !
  IF(iPlantPhotosynsType_pft(NZ).EQ.ic4_photo .AND. CanopyNonstElmCopy_brch(ielmc).GT.ZERO4Groth_pft(NZ))THEN
    FracNonstLeft    = CanopyNonstElmAfhvst_brch(ielmc)/CanopyNonstElmCopy_brch(ielmc)
    FracNonstRemoved = 1._r8-FracNonstLeft
    D9810: DO K   = 1, MaxNodesPerBranch1
      CanopyNonstElmRemoval(ielmc)=CanopyNonstElmRemoval(ielmc) &
        +FracNonstRemoved*(CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ) &
        +CMassCO2BundleSheath_node(K,NB,NZ)+CMassHCO3BundleSheath_node(K,NB,NZ))
      CPOOL3_node(K,NB,NZ)                = FracNonstLeft*CPOOL3_node(K,NB,NZ)
      CPOOL4_node(K,NB,NZ)                = FracNonstLeft*CPOOL4_node(K,NB,NZ)
      CMassCO2BundleSheath_node(K,NB,NZ)  = FracNonstLeft*CMassCO2BundleSheath_node(K,NB,NZ)
      CMassHCO3BundleSheath_node(K,NB,NZ) = FracNonstLeft*CMassHCO3BundleSheath_node(K,NB,NZ)
    ENDDO D9810
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine CutBranchNonstructural

!----------------------------------------------------------------------------------------------------
  subroutine CutBranchSheathPetole(yearIJ,NB,NZ,HarvestedPetoleC,FracIntnodeNotHvsted,&
    FracNodeNotHvsted,BranchLength,PetolShethCAfHvst_brch,PetolShethCB4Hvst_brch)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NB,NZ
  real(r8), intent(in) :: HarvestedPetoleC
  real(r8), intent(inout) :: FracIntnodeNotHvsted(0:MaxNodesPerBranch1)
  real(r8), intent(inout) :: FracNodeNotHvsted(0:MaxNodesPerBranch1)  
  real(r8), intent(out) :: BranchLength
  real(r8), intent(out) :: PetolShethCAfHvst_brch,PetolShethCB4Hvst_brch

  real(r8) :: FinePetoleMassE,WoodyPetoleMassE
  character(len=*), parameter :: subname='CutBranchSheathPetole'
  real(r8) :: PetoleCRmved_brch,FHGT
  integer :: K,NE

  associate(                                                             &
    PetolShethStrutElms_pft    => plt_biom%PetolShethStrutElms_pft      ,& !input  :canopy sheath structural element mass, [g d-2]
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft               ,& !input  :threshold zero for plang growth calculation, [-]
    CanopyCutProxy_pft         => plt_distb%CanopyCutProxy_pft          ,& !input  :harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
    StalkNodeHeight_brch       => plt_morph%StalkNodeHeight_brch        ,& !input  :internode height, [m]
    k_fine_comp                => pltpar%k_fine_comp                    ,& !input  :fine litter complex id
    k_woody_comp               => pltpar%k_woody_comp                   ,& !input  :woody litter complex id
    FracPetolShethAlloc2Litr   => plt_allom%FracPetolShethAlloc2Litr    ,& !input  :woody element allocation, [-]
    iHarvstType_pft            => plt_distb%iHarvstType_pft             ,& !input  :type of harvest,[-]
    PetolShethStrutElms_brch   => plt_biom%PetolShethStrutElms_brch     ,& !inoput :branch sheath structural element, [g d-2]
    PetolShethElmntNode_brch   => plt_biom%PetolShethElmntNode_brch     ,& !inoput :sheath chemical element, [g d-2]
    PetoleLength_node          => plt_morph%PetoleLength_node           ,& !inoput :sheath height, [m]
    PetoleProteinC_node        => plt_biom%PetoleProteinC_node           & !inoput :layer sheath protein C, [g d-2]
  )
  call PrintInfo('beg '//subname)
  PetolShethCB4Hvst_brch=0._r8;PetolShethCAfHvst_brch=0._r8
  !
  !     CUT SHEATHS OR PetolShethS AND STALKS HARVESTED NODES AND LAYERS
  !
  !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
  !                       ,3=pruning,4=grazing,5=fire,6=herbivory
  !     PetoleCRmved_brch,HarvestedPetoleC=branch, PFT PetolSheth C mass removed
  !     StalkNodeHeight_brch=internode length
  !     BranchLength=internode length removed
  !
  !grazing or herbivory
  IF((iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo) &
    .AND. PetolShethStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
    PetoleCRmved_brch=HarvestedPetoleC*PetolShethStrutElms_brch(ielmc,NB,NZ)/PetolShethStrutElms_pft(ielmc,NZ)
  ELSE
    PetoleCRmved_brch=0._r8
  ENDIF

  BranchLength=0._r8
  D9805: DO K=MaxNodesPerBranch1,0,-1
    IF(StalkNodeHeight_brch(K,NB,NZ).GT.0.0_r8) BranchLength=AMAX1(BranchLength,StalkNodeHeight_brch(K,NB,NZ))
    !
    !     HARVESTED SHEATH OR PetolSheth C,N,P
    !
    !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
    !                       ,3=pruning,4=grazing,5=fire,6=herbivory
    !     PetoleCRmved_brch=branch PetolSheth C mass removed
    !     PetolShethElmntNode_brch,WGSHN,WGSHP,PetoleProteinC_node=node PetolSheth C,N,P,protein mass
    !     FracIntnodeNotHvsted=fraction of internode layer mass not harvested
    !     PetoleShethElmntRemoval=harvested PetolSheth C,N,P
    !     PetolShethElmntHarv2Litr=harvested PetolSheth C,N,P to litter
    !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
    !     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
    !     PetoleLength_node,StalkNodeHeight_brch=PetolSheth,internode length
    !
    !neither grazing nor herbivory
    IF((iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo) &
      .OR. PetoleCRmved_brch.GT.0.0_r8)THEN
      !
      !grazing or herbivory
      IF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
        IF(PetolShethElmntNode_brch(ielmc,K,NB,NZ).GT.PetoleCRmved_brch)THEN
          FracIntnodeNotHvsted(K) = (PetolShethElmntNode_brch(ielmc,K,NB,NZ)-PetoleCRmved_brch)/PetolShethElmntNode_brch(ielmc,K,NB,NZ)
          FracIntnodeNotHvsted(K) = AZMAX1(AMIN1(1.0_r8,FracIntnodeNotHvsted(K)))
          FracNodeNotHvsted(K)    = FracIntnodeNotHvsted(K)
        ELSE
          FracIntnodeNotHvsted(K) = 0._r8
          FracNodeNotHvsted(K)    = 0._r8
        ENDIF
      ENDIF
      !
      !remove PetolSheth mass from node K
      PetoleCRmved_brch=PetoleCRmved_brch-(1._r8-FracIntnodeNotHvsted(K))*PetolShethElmntNode_brch(ielmc,K,NB,NZ)

      DO NE=1,NumPlantChemElms
        FinePetoleMassE              = PetolShethElmntNode_brch(NE,K,NB,NZ)*FracPetolShethAlloc2Litr(NE,k_fine_comp)
        PetoleShethElmntRemoval(NE)  = PetoleShethElmntRemoval(NE)+(1._r8-FracNodeNotHvsted(K))*FinePetoleMassE
        PetolShethElmntHarv2Litr(NE) = PetolShethElmntHarv2Litr(NE)+(FracNodeNotHvsted(K)-FracIntnodeNotHvsted(K))*FinePetoleMassE

        WoodyPetoleMassE        = PetolShethElmntNode_brch(NE,K,NB,NZ)*FracPetolShethAlloc2Litr(NE,k_woody_comp)
        WoodyElmntRemoval(NE)   = WoodyElmntRemoval(NE)+(1._r8-FracNodeNotHvsted(K))*WoodyPetoleMassE
        WoodyElmntHarv2Litr(NE) = WoodyElmntHarv2Litr(NE)+(FracNodeNotHvsted(K)-FracIntnodeNotHvsted(K))*WoodyPetoleMassE
      ENDDO
      !
      !     ACCUMULATE REMAINING SHEATH OR PetolSheth C,N,P AND LENGTH
      !
      !     PetolShethElmntNode_brch=PetolSheth node C mass
      !     WTSHEB,WTSHBN,WTSHBP=branch PetolSheth C,N,P mass
      !     PetoleLength_node=node PetolSheth height
      !     PetoleProteinC_node=PetolSheth protein mass
      !
      PetolShethCB4Hvst_brch       = PetolShethCB4Hvst_brch+PetolShethElmntNode_brch(ielmc,K,NB,NZ)
      PetoleProteinC_node(K,NB,NZ) = FracIntnodeNotHvsted(K)*PetoleProteinC_node(K,NB,NZ)

      DO NE=1,NumPlantChemElms
        PetolShethStrutElms_brch(NE,NB,NZ)=PetolShethStrutElms_brch(NE,NB,NZ) &
          -(1._r8-FracIntnodeNotHvsted(K))*PetolShethElmntNode_brch(NE,K,NB,NZ)
        PetolShethElmntNode_brch(NE,K,NB,NZ)=FracIntnodeNotHvsted(K)*PetolShethElmntNode_brch(NE,K,NB,NZ)
      ENDDO
      !      PetoleProteinC_node(K,NB,NZ)=FracIntnodeNotHvsted(K)*PetoleProteinC_node(K,NB,NZ)
      !
      IF(iHarvstType_pft(NZ).LE.iharvtyp_allabvg .AND. PetoleLength_node(K,NB,NZ).GT.0.0_r8)THEN
        FHGT=AZMAX1(AMIN1(1.0_r8,(StalkNodeHeight_brch(K,NB,NZ) &
          +PetoleLength_node(K,NB,NZ)-CanopyCutProxy_pft(NZ))/PetoleLength_node(K,NB,NZ)))
        PetoleLength_node(K,NB,NZ)=(1._r8-FHGT)*PetoleLength_node(K,NB,NZ)
      ELSE
        PetoleLength_node(K,NB,NZ)=FracIntnodeNotHvsted(K)*PetoleLength_node(K,NB,NZ)
      ENDIF
      PetolShethCAfHvst_brch=PetolShethCAfHvst_brch+PetolShethElmntNode_brch(ielmc,K,NB,NZ)

      !     IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing.AND.iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      !     IF(StalkNodeHeight_brch(K,NB,NZ).GT.CanopyCutProxy_pft(NZ)
      !    2.OR.iHarvstType_pft(NZ).EQ.iharvtyp_pruning)THEN
      !     IF(isclose(FracIntnodeNotHvsted(K),0._r8).AND.K.GT.0)THEN
      !     IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.(.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ)))THEN
      !     NumOfLeaves_brch(NB,NZ)=AZMAX1(NumOfLeaves_brch(NB,NZ)-1.0)
      !     ELSE
      !     NumOfLeaves_brch(NB,NZ)=AZMAX1(NumOfLeaves_brch(NB,NZ)-0.04)
      !     ENDIF
      !     ENDIF
      !     ENDIF
      !     ENDIF
    ENDIF
  ENDDO D9805
  call PrintInfo('end '//subname)
  end associate
  end subroutine CutBranchSheathPetole

!----------------------------------------------------------------------------------------------------
  subroutine HarvestCanopy(yearIJ,NZ,HarvestedLeafC,LeafLayerC_brch)
  !
  !Canopy harvest
  !from layer to branch and to leaf node

  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), intent(in):: HarvestedLeafC
  REAL(R8), intent(in) :: LeafLayerC_brch(NumCanopyLayers1,JP1,JP1)

  character(len=*), parameter :: subname='HarvestCanopy'
  integer :: L,NB,K,NE
  real(r8) :: FracHeightLeft,FHVSH,FracHarvested
  real(r8) :: FrcLeafMassLeft
  real(r8) :: HvestedLeafCLayer_brch
  real(r8) :: FracHvst2Litr

  associate(                                                               &
    CanopyHeightZ_col           => plt_morph%CanopyHeightZ_col            ,& !input  :canopy layer height, [m]
    CanopyCutProxy_pft          => plt_distb%CanopyCutProxy_pft           ,& !input  :harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
    FracBiomHarvsted            => plt_distb%FracBiomHarvsted             ,& !input  :harvest efficiency, [-]
    k_fine_comp                 => pltpar%k_fine_comp                     ,& !input  :fine litter complex id
    k_woody_comp                => pltpar%k_woody_comp                    ,& !input  :woody litter complex id
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft            ,& !input  :number of branches,[-]
    ZERO4LeafVar_pft            => plt_biom%ZERO4LeafVar_pft              ,& !input  :threshold zero for leaf calculation, [-]
    THIN_pft                    => plt_distb%THIN_pft                     ,& !input  :thinning of plant population, [-]
    LeafStrutElms_pft           => plt_biom%LeafStrutElms_pft             ,& !input  :canopy leaf structural element mass, [g d-2]
    FracLeafShethElmAlloc2Litr  => plt_allom%FracLeafShethElmAlloc2Litr   ,& !input  :leaf element allocation,[-]
    iHarvstType_pft             => plt_distb%iHarvstType_pft              ,& !input  :type of harvest,[-]
    LeafLayerElms_node          => plt_biom%LeafLayerElms_node            ,& !inoput :layer leaf element, [g d-2]
    CanopyStalkSurfArea_lbrch   => plt_morph%CanopyStalkSurfArea_lbrch    ,& !inoput :plant canopy layer branch stem area, [m2 d-2]
    CanopyLeafArea_lnode        => plt_morph%CanopyLeafArea_lnode         ,& !inoput :layer/node/branch leaf area, [m2 d-2]
    CanopyStemSurfAreaZ_pft     => plt_morph%CanopyStemSurfAreaZ_pft      ,& !inoput :plant canopy layer stem area, [m2 d-2]
    CanopyLeafAreaZ_pft         => plt_morph%CanopyLeafAreaZ_pft          ,& !output :canopy layer leaf area, [m2 d-2]
    CanopyLeafCLyr_pft          => plt_biom%CanopyLeafCLyr_pft             & !output :canopy layer leaf C, [g d-2]
  )
  call PrintInfo('beg '//subname)
  !
  !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
  !                       ,3=pruning,4=grazing,5=fire,6=herbivory
  !     ZL=height to bottom of each canopy layer
  !     FracHeightLeft=fraction of canopy layer height not harvested
  !     FrcLeafMassLeft=fraction of canopy layer mass not harvested
  !     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
  !          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)  
  !
  D9865: DO L=NumCanopyLayers1,1,-1
    !neither grazing nor herbivory
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      IF(iHarvstType_pft(NZ).NE.iharvtyp_pruning)THEN
        IF(CanopyHeightZ_col(L).GT.CanopyHeightZ_col(L-1))THEN
          FracHeightLeft=AZMAX1(AMIN1(1.0_r8,1._r8-((CanopyHeightZ_col(L))-CanopyCutProxy_pft(NZ))/ &
            (CanopyHeightZ_col(L)-CanopyHeightZ_col(L-1))))
        ELSE
          FracHeightLeft=1.0_r8
        ENDIF
      ELSE
        FracHeightLeft=0._r8
      ENDIF

      IF(isclose(THIN_pft(NZ),0._r8))THEN
        FrcLeafMassLeft = AZMAX1(1._r8-(1._r8-FracHeightLeft)*FracBiomHarvsted(iHarvst_pft,iplthvst_leaf,NZ))
        FHVSH           = FrcLeafMassLeft
      ELSE
        FrcLeafMassLeft=AZMAX1(1._r8-THIN_pft(NZ))
        IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
          FHVSH=1.0_r8-(1._r8-FracHeightLeft)*FracBiomHarvsted(iHarvst_pft,iplthvst_leaf,NZ)*THIN_pft(NZ)
        ELSE
          FHVSH=FrcLeafMassLeft
        ENDIF
      ENDIF
    ELSE
      FrcLeafMassLeft = 0._r8
      FHVSH           = 0._r8
    ENDIF
    !
    !     CUT LEAVES AT HARVESTED NODES AND LAYERS
    !
    !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
    !                       ,3=pruning,4=grazing,5=fire,6=herbivory
    !     WTLF=PFT leaf C mass
    !     LeafLayerC_brch=branch leaf C mass in canopy layer
    !     WHVBSL,HarvestedLeafC=layer,total leaf C mass removed
    !     WGLFL=leaf node C in canopy layer
    !     FrcLeafMassLeft=fraction of leaf node mass not harvested
    ! 
    D9855: DO NB=1,NumOfBranches_pft(NZ)
      !grazing or herbivory
      IF((iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo) &
        .AND. LeafStrutElms_pft(ielmc,NZ).GT.ZERO4LeafVar_pft(NZ))THEN
        HvestedLeafCLayer_brch=HarvestedLeafC*AZMAX1(LeafLayerC_brch(L,NB,NZ))/LeafStrutElms_pft(ielmc,NZ)
      ELSE
        HvestedLeafCLayer_brch=0._r8
      ENDIF

      !distribute the harvest to leaf nodes on branch NB
      D9845: DO K=MaxNodesPerBranch1,0,-1
        !neither grazing nor herbivory
        IF((iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo) &
          .OR. HvestedLeafCLayer_brch.GT.0.0_r8)THEN

          !grazing or herbivory          
          IF(iHarvstType_pft(NZ).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN            
            IF(LeafLayerElms_node(ielmc,L,K,NB,NZ).GT.HvestedLeafCLayer_brch)THEN
              !leaf C in node is greater than harvest demand
              FrcLeafMassLeft=AZMAX1(AMIN1(1.0_r8,(LeafLayerElms_node(ielmc,L,K,NB,NZ)-HvestedLeafCLayer_brch) &
                /LeafLayerElms_node(ielmc,L,K,NB,NZ)))
              FHVSH=FrcLeafMassLeft
            ELSE
              !leaf C in node is harvest demand
              FrcLeafMassLeft = 0.0_r8
              FHVSH           = 0.0_r8
            ENDIF
          ENDIF
          !
          !     HARVESTED LEAF AREA, C, N, P
          !
          !     FrcLeafMassLeft=fraction of leaf node mass not harvested
          !     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
          !
          HvestedLeafCLayer_brch = HvestedLeafCLayer_brch-(1._r8-FrcLeafMassLeft)*LeafLayerElms_node(ielmc,L,K,NB,NZ)
          FracHarvested          = 1._r8-FHVSH
          FracHvst2Litr          = FHVSH-FrcLeafMassLeft
          DO NE = 1, NumPlantChemElms
            LeafElmntRemoval(NE)=LeafElmntRemoval(NE) &
              +FracHarvested*LeafLayerElms_node(NE,L,K,NB,NZ)*FracLeafShethElmAlloc2Litr(NE,k_fine_comp)

            LeafElmntHarv2Litr(NE)=LeafElmntHarv2Litr(NE) &
              +FracHvst2Litr*LeafLayerElms_node(NE,L,K,NB,NZ)*FracLeafShethElmAlloc2Litr(NE,k_fine_comp)

            WoodyElmntRemoval(NE)=WoodyElmntRemoval(NE) &
              +FracHarvested*LeafLayerElms_node(NE,L,K,NB,NZ)*FracLeafShethElmAlloc2Litr(NE,k_woody_comp)

            WoodyElmntHarv2Litr(NE)=WoodyElmntHarv2Litr(NE) &
              +FracHvst2Litr*LeafLayerElms_node(NE,L,K,NB,NZ)*FracLeafShethElmAlloc2Litr(NE,k_woody_comp)

            LeafLayerElms_node(NE,L,K,NB,NZ)=FrcLeafMassLeft*LeafLayerElms_node(NE,L,K,NB,NZ)
          ENDDO
          !
          !     REMAINING LEAF C,N,P AND AREA
          !
          CanopyLeafArea_lnode(L,K,NB,NZ)=FrcLeafMassLeft*CanopyLeafArea_lnode(L,K,NB,NZ)
          IF(K.EQ.1)THEN
            CanopyStalkSurfArea_lbrch(L,NB,NZ)=FrcLeafMassLeft*CanopyStalkSurfArea_lbrch(L,NB,NZ)
          ENDIF
        ENDIF
      ENDDO D9845
    ENDDO D9855
    CanopyLeafAreaZ_pft(L,NZ) = 0._r8
    CanopyLeafCLyr_pft(L,NZ)  = 0._r8
    CanopyStemSurfAreaZ_pft(L,NZ) = CanopyStemSurfAreaZ_pft(L,NZ)*FrcLeafMassLeft
  ENDDO D9865
  call PrintInfo('end '//subname)
  end associate
  end subroutine HarvestCanopy

!----------------------------------------------------------------------------------------------------
  subroutine CutBranchLeaves(yearIJ,NB,NZ,LeafCafCut_brch,LeafCB4Cut_brch,FracIntnodeNotHvsted,FracNodeNotHvsted)

  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NB,NZ
  real(r8), intent(out) :: LeafCafCut_brch,LeafCB4Cut_brch            !leaf C before and after cut
  real(r8), intent(out) :: FracIntnodeNotHvsted(0:MaxNodesPerBranch1)  
  real(r8), intent(out) :: FracNodeNotHvsted(0:MaxNodesPerBranch1)

  character(len=*), parameter :: subname='CutBranchLeaves'
  integer :: K,NE,L
  real(r8) :: ARLFG               !leaf area accumulator
  real(r8) :: FracIntnodeKHvsted
  real(r8) :: LeafElmNodeK_brch(NumPlantChemElms)  

  associate(                                                        &
    CanopyLeafArea_lnode     => plt_morph%CanopyLeafArea_lnode     ,& !input  :layer/node/branch leaf area, [m2 d-2]
    ZERO4Groth_pft           => plt_biom%ZERO4Groth_pft            ,& !input  :threshold zero for plang growth calculation, [-]
    THIN_pft                 => plt_distb%THIN_pft                 ,& !input  :thinning of plant population, [-]
    FracBiomHarvsted         => plt_distb%FracBiomHarvsted         ,& !input  :harvest efficiency, [-]
    LeafLayerElms_node       => plt_biom%LeafLayerElms_node        ,& !input  :layer leaf element, [g d-2]
    iHarvstType_pft          => plt_distb%iHarvstType_pft          ,& !input  :type of harvest,[-]
    LeafElmntNode_brch       => plt_biom%LeafElmntNode_brch        ,& !inoput :leaf element, [g d-2]
    CanopyLeafCLyr_pft       => plt_biom%CanopyLeafCLyr_pft        ,& !inoput :canopy layer leaf C, [g d-2]
    LeafAreaLive_brch        => plt_morph%LeafAreaLive_brch        ,& !inoput :branch leaf area, [m2 d-2]
    LeafArea_node            => plt_morph%LeafArea_node            ,& !inoput :leaf area, [m2 d-2]
    LeafProteinC_node        => plt_biom%LeafProteinC_node         ,& !inoput :layer leaf protein C, [g d-2]
    LeafStrutElms_brch       => plt_biom%LeafStrutElms_brch        ,& !inoput :branch leaf structural element mass, [g d-2]
    CanopyLeafAreaZ_pft      => plt_morph%CanopyLeafAreaZ_pft       & !inoput :canopy layer leaf area, [m2 d-2]
  )

  call PrintInfo('beg '//subname)
  LeafCafCut_brch=0._r8;LeafCB4Cut_brch=0._r8

  D9825: DO K=0,MaxNodesPerBranch1
    ARLFG=0._r8
    LeafElmNodeK_brch(1:NumPlantChemElms)=0._r8
    !
    !     ACCUMULATE REMAINING LEAF AREA, C, N, P
    !
    !     CanopyLeafArea_lnode,CanopyLeafAreaZ_pft=leaf node,total area in canopy layer
    !
    D9815: DO L=1,NumCanopyLayers1
      ARLFG=ARLFG+CanopyLeafArea_lnode(L,K,NB,NZ)
      DO NE=1,NumPlantChemElms
        LeafElmNodeK_brch(NE)=LeafElmNodeK_brch(NE)+LeafLayerElms_node(NE,L,K,NB,NZ)
      ENDDO
      CanopyLeafAreaZ_pft(L,NZ) = CanopyLeafAreaZ_pft(L,NZ)+CanopyLeafArea_lnode(L,K,NB,NZ)
      CanopyLeafCLyr_pft(L,NZ)  = CanopyLeafCLyr_pft(L,NZ)+LeafLayerElms_node(ielmc,L,K,NB,NZ)
    ENDDO D9815
    !
    !     CUT STALK AT HARVESTED NODES AND LAYERS
    !
    !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
    !                       ,3=pruning,4=grazing,5=fire,6=herbivory
    !     WGLF=leaf node C mass
    !     FracBiomHarvsted(iHarvst_pft,1:4=fraction of leaf,non-foliar,woody, standing dead removed from PFT
    !     FracIntnodeNotHvsted=fraction of internode layer mass not harvested
    !     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
    !          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
    !neither grazing nor herbivory
    IF(iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)THEN
      IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ) .AND. FracBiomHarvsted(iHarvst_pft,iplthvst_leaf,NZ).GT.0.0)THEN
        FracIntnodeKHvsted=(1._r8-AZMAX1(LeafElmNodeK_brch(ielmc))/LeafElmntNode_brch(ielmc,K,NB,NZ)) &
          *FracBiomHarvsted(iHarvst_pft,iplthvst_finenonleaf,NZ)/FracBiomHarvsted(iHarvst_pft,iplthvst_leaf,NZ)
        FracIntnodeNotHvsted(K) = AZMAX1(AMIN1(1.0_r8,(1._r8-FracIntnodeKHvsted)))
        FracNodeNotHvsted(K)    = FracIntnodeNotHvsted(K)
      ELSE
        IF(isclose(THIN_pft(NZ),0._r8))THEN
          FracIntnodeNotHvsted(K) = 1.0_r8-FracBiomHarvsted(iHarvst_pft,iplthvst_finenonleaf,NZ)
          FracNodeNotHvsted(K)    = FracIntnodeNotHvsted(K)
        ELSE
          FracIntnodeNotHvsted(K)=1.0_r8-THIN_pft(NZ)
          IF(iHarvstType_pft(NZ).EQ.iharvtyp_none)THEN
            FracNodeNotHvsted(K)=1.0_r8-FracBiomHarvsted(iHarvst_pft,iplthvst_finenonleaf,NZ)*THIN_pft(NZ)
          ELSE
            FracNodeNotHvsted(K)=FracIntnodeNotHvsted(K)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      FracIntnodeNotHvsted(K) = 0._r8
      FracNodeNotHvsted(K)    = 0._r8
    ENDIF
    !
    !     ACCUMULATE REMAINING BRANCH LEAF AREA, C, N, P
    !
    !     WGLF=leaf node C mass
    !     LeafAreaLive_brch,LeafArea_node=branch,node leaf area
    !     LeafProteinC_node=leaf protein mass
    !
    LeafCB4Cut_brch=LeafCB4Cut_brch+LeafElmntNode_brch(ielmc,K,NB,NZ)
    !-total cut + remain
    DO NE=1,NumPlantChemElms
      LeafStrutElms_brch(NE,NB,NZ)=LeafStrutElms_brch(NE,NB,NZ)-LeafElmntNode_brch(NE,K,NB,NZ)+LeafElmNodeK_brch(NE)
    ENDDO
    !
    LeafAreaLive_brch(NB,NZ)=LeafAreaLive_brch(NB,NZ)-LeafArea_node(K,NB,NZ)+ARLFG
    IF(LeafArea_node(K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      LeafProteinC_node(K,NB,NZ)=LeafProteinC_node(K,NB,NZ)*ARLFG/LeafArea_node(K,NB,NZ)
    ELSE
      LeafProteinC_node(K,NB,NZ)=0._r8
    ENDIF
    LeafArea_node(K,NB,NZ)=ARLFG

    DO NE=1,NumPlantChemElms
      LeafElmntNode_brch(NE,K,NB,NZ)=LeafElmNodeK_brch(NE)
    ENDDO
    LeafCafCut_brch=LeafCafCut_brch+LeafElmntNode_brch(ielmc,K,NB,NZ)
  ENDDO D9825
  call PrintInfo('end '//subname)
  end associate
  end subroutine CutBranchLeaves    

!----------------------------------------------------------------------------------------------------
  subroutine CutPlant(yearIJ,NZ,HarvestedPetoleC,GrazedCanopyNonstC,GrazedCanopyNoduleC,HarvestedShethC,&
    HarvestedGrainC,HarvestedEarC,HarvestedStalkRsrvC,HarvestedStalkC)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), intent(in) :: HarvestedPetoleC
  real(r8), intent(in) :: GrazedCanopyNonstC
  real(r8), intent(in) :: GrazedCanopyNoduleC
  real(r8), intent(in) :: HarvestedShethC
  real(r8), intent(in) :: HarvestedGrainC
  real(r8), intent(in) :: HarvestedEarC
  real(r8), intent(in) :: HarvestedStalkRsrvC
  real(r8), intent(in) :: HarvestedStalkC

  character(len=*), parameter :: subname='CutPlant'
  integer :: NB,I,J
  real(r8) :: FracNodeNotHvsted(0:MaxNodesPerBranch1),FracIntnodeNotHvsted(0:MaxNodesPerBranch1)
  real(r8) :: LeafCafCut_brch,PetolShethCAfHvst_brch,LeafCB4Cut_brch,PetolShethCB4Hvst_brch
  real(r8) :: VOLWPX,WVPLT  
  real(r8) :: BranchLength  !Branch node length subject to removal
  real(r8) :: FDM  
  real(r8) :: watflx,VHeatCapCanopyPrev,CanopyMassC,ThermalCarbon,StemArea,DeadCarbon
  real(r8) :: tvegE(NumPlantChemElms),tvegE1(NumPlantChemElms)
  
  call PrintInfo('beg '//subname)
  associate(                                                           &
    WatHeldOnCanopy_pft        => plt_ew%WatHeldOnCanopy_pft, &
    SnowOnCanopy_pft           => plt_ew%SnowOnCanopy_pft, &
    CanopyStalkSurfArea_lbrch  => plt_morph%CanopyStalkSurfArea_lbrch, &
    StandDeadSurfArea_pft      => plt_morph%StandDeadSurfArea_pft, &
    StandDeadCompKElms_pft     => plt_biom%StandDeadCompKElms_pft, &
    TKC_pft                   => plt_ew%TKC_pft                       ,& !input  :canopy temperature, [K]
    PlantPopuLive_pft         => plt_site%PlantPopuLive_pft           ,& !input  :plant population, [d-2]
    jHarvstType_pft           => plt_distb%jHarvstType_pft            ,& !input  :flag for stand replacing disturbance,[-]
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft  ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft      ,& !input  :plant growth type (vascular, non-vascular),[-]
    PetolShethStrutElms_brch  => plt_biom%PetolShethStrutElms_brch    ,& !input  :branch sheath structural element, [g d-2]
    CanopyCutProxy_pft        => plt_distb%CanopyCutProxy_pft         ,& !input  :harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
    PSICanopy_pft             => plt_ew%PSICanopy_pft                 ,& !input  :canopy total water potential, [Mpa]
    SapwoodBiomassC_brch     => plt_biom%SapwoodBiomassC_brch      ,& !input  :remaining branch sapwood C, [gC d-2]
    CanopySapwoodC_pft        => plt_biom%CanopySapwoodC_pft          ,& !input  :canopy active stalk C, [g d-2]
    CanopyLeafSheathC_pft     => plt_biom%CanopyLeafSheathC_pft       ,& !input  :canopy leaf + sheath C, [g d-2]
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch          ,& !input  :branch leaf structural element mass, [g d-2]
    CanopyHeightLive_pft      => plt_morph%CanopyHeightLive_pft       ,& !input  :canopy height, [m]
    iPlant2ndGrothPattern_pft => plt_pheno%iPlant2ndGrothPattern_pft  ,& !input  :plant expression of secondary growth, [-]                
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft          ,& !input  :number of branches,[-]
    StemSpecVolume_pft        => plt_morph%StemSpecVolume_pft         ,& !input  :stalk specific volume, [m3 gC-1]        
    iHarvstType_pft           => plt_distb%iHarvstType_pft            ,& !input  :type of harvest,[-]
    CanopyBiomWater_pft       => plt_ew%CanopyBiomWater_pft           ,& !inoput :canopy water content, [m3 d-2]
    QCanopyWatLoss2Dist_col       => plt_ew%QCanopyWatLoss2Dist_col           ,& !inoput :canopy water +/- due to disturbance, [m3 H2O/d2]
    VHeatCapCanopy_pft        => plt_ew%VHeatCapCanopy_pft            ,& !inoput :canopy heat capacity, [MJ d-2 K-1]
    CanopyHeatLoss2Dist_col   => plt_ew%CanopyHeatLoss2Dist_col       ,& !inoput :canopy energy +/- due to disturbance, [MJ /d2]
    QH2OLoss_lnds             => plt_site%QH2OLoss_lnds               ,& !inoput :total subsurface water loss flux over the landscape, [m3 d-2]
    H2OLoss_CumYr_col         => plt_ew%H2OLoss_CumYr_col             ,& !inoput :total subsurface water flux, [m3 d-2]
    CanopyLeafSheathC_brch    => plt_biom%CanopyLeafSheathC_brch      ,& !output :plant branch leaf + sheath C, [g d-2]
    isPlantBranchAlive_brch   => plt_pheno%isPlantBranchAlive_brch     & !output :flag to detect branch death, [-]
  )
  I=yearIJ%I;J=yearIJ%J
  
  D9835: DO NB=1,NumOfBranches_pft(NZ)

    CALL CutBranchLeaves(yearIJ,NB,NZ,LeafCafCut_brch,LeafCB4Cut_brch,FracIntnodeNotHvsted,FracNodeNotHvsted)

    call CutBranchSheathPetole(yearIJ,NB,NZ,HarvestedPetoleC,FracIntnodeNotHvsted,FracNodeNotHvsted,BranchLength,PetolShethCAfHvst_brch,PetolShethCB4Hvst_brch)
    
    call CutBranchNonstructural(I,J,NB,NZ,LeafCafCut_brch,PetolShethCAfHvst_brch,LeafCB4Cut_brch,PetolShethCB4Hvst_brch,GrazedCanopyNonstC,GrazedCanopyNoduleC)    
    !
    !     CUT STALKS
    call CutBranchStalk(yearIJ,NB,NZ,BranchLength,HarvestedStalkC,HarvestedStalkRsrvC)
    !
    !     CUT REPRODUCTIVE ORGANS FracHuskNotHvsted
    call CutBranchReprodOrgans(I,J,NB,NZ,BranchLength,HarvestedShethC,HarvestedGrainC,HarvestedEarC)
    !
    !     REMAINING TOTAL BRANCH C,N,P AND LEAF, STALK AREA
    !
    !     C4PhotoShootNonstC_brch=total C4 nonstructural C in branch
    !     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
    !     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
    !     SapwoodBiomassC_brch=stalk sapwood mass
    !     PSICanopy_pft=canopy water potential
    !     CanopyBiomWater_pft=water volume in canopy
    !     QH2OLoss_lnds,H2OLoss_CumYr_col=accumulated water loss for water balance calculation
    !
    call SumPlantBranchBiome(NB,NZ)
    !
    CanopyLeafSheathC_brch(NB,NZ)=AZMAX1(LeafStrutElms_brch(ielmc,NB,NZ)+PetolShethStrutElms_brch(ielmc,NB,NZ))

    !
    !     RESET PHENOLOGY, GROWTH STAGE IF STALKS ARE CUT
    !
    !     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+PetolSheth,2=none,3=between 1,2
    !     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
    !     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
    !                       ,3=pruning,4=grazing,5=fire,6=herbivory
    !     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
    !          iHarvstType_pft=3:reduction of clumping factor
    !          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
    !
    IF((iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. (.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))) & !grass-like
      .AND. (iHarvstType_pft(NZ).NE.iharvtyp_grazing .AND. iHarvstType_pft(NZ).NE.iharvtyp_herbivo)        & !neither grazing nor herbivory
      .AND. CanopyHeightLive_pft(NZ).GT.CanopyCutProxy_pft(NZ))THEN
      call ResetCutBranch(I,J,NZ,NB)
    ENDIF
    !
    !     DEATH OF BRANCH IF KILLING HARVEST ENTERED IN 'READQ'
    !
    !     jHarvstType_pft=terminate PFT:0=no,1=yes,2=yes,and reseed
    !     isPlantBranchAlive_brch=branch living flag: 0=alive,1=dead
    !     PP=PFT population
    !     WTLS=total PFT leaf+PetolSheth C mass
    !     WTSTK=total PFT stalk C mass
    !     WVSTK=total PFT sapwood C mass
    !     CanopyStalkSurfArea_lbrch=total PFT stalk surface area
    !    
    IF(jHarvstType_pft(NZ).NE.jharvtyp_noaction .or. PlantPopuLive_pft(NZ).LE.0.0_r8)then
      isPlantBranchAlive_brch(NB,NZ)=iFalse      
    endif
    
  ENDDO D9835  
  !Keep pre-harvest PFT totals fixed while allocating removal among branches.
  !After all branches are cut, update canopy water and heat once per PFT.
  VOLWPX                    = CanopyBiomWater_pft(NZ)
  VHeatCapCanopyPrev        = VHeatCapCanopy_pft(NZ)
  CanopyLeafSheathC_pft(NZ) = SUM(CanopyLeafSheathC_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopySapwoodC_pft(NZ)    = SUM(SapwoodBiomassC_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopyMassC               = AZMAX1(CanopyLeafSheathC_pft(NZ)+CanopySapwoodC_pft(NZ))
  FDM                       = get_FDM(PSICanopy_pft(NZ))    !drymatter/water = fdm

  CanopyBiomWater_pft(NZ) = 1.e-6_r8*CanopyMassC/FDM
  watflx              = VOLWPX-CanopyBiomWater_pft(NZ)
  QH2OLoss_lnds       = QH2OLoss_lnds+VOLWPX-CanopyBiomWater_pft(NZ)
  H2OLoss_CumYr_col   = H2OLoss_CumYr_col+watflx

  StemArea=SUM(CanopyStalkSurfArea_lbrch(:,1:NumOfBranches_pft(NZ),NZ))
  ThermalCarbon=CanopyLeafSheathC_pft(NZ)+CanopyStemThermalCarbon(CanopySapwoodC_pft(NZ), &
    StemArea,StemSpecVolume_pft(NZ))
  DeadCarbon=CanopyStemThermalCarbon(SUM(StandDeadCompKElms_pft(ielmc,:,NZ)), &
    StandDeadSurfArea_pft(NZ),StemSpecVolume_pft(NZ))
  VHeatCapCanopy_pft(NZ)=cpo*gOC_to_m3_OM(ThermalCarbon+DeadCarbon) &
    +cpw*(CanopyBiomWater_pft(NZ)+WatHeldOnCanopy_pft(NZ))+cps*SnowOnCanopy_pft(NZ)
  QCanopyWatLoss2Dist_col = QCanopyWatLoss2Dist_col+watflx
  CanopyHeatLoss2Dist_col = CanopyHeatLoss2Dist_col+(VHeatCapCanopyPrev-VHeatCapCanopy_pft(NZ))*TKC_pft(NZ)
  end associate
  call PrintInfo('end '//subname)
  end subroutine CutPlant

!----------------------------------------------------------------------------------------------------

  subroutine RootRemovalL4Annual(yearIJ,N,L,NZ,FracLeftThin,XHVST1)

  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer , intent(in) :: N,L,NZ
  real(r8), intent(out) :: FracLeftThin
  real(r8), intent(out):: XHVST1
  character(len=*), parameter :: subname='RootRemovalL4Annual'

  real(r8) :: HarvestedBiomass(NumPlantChemElms)
  real(r8) :: FFIRE(NumPlantChemElms)
  integer  :: NR,idg,M,NE

  associate(                                                          &
    PlantElmAllocMat4Litr     => plt_soilchem%PlantElmAllocMat4Litr  ,& !input  :litter kinetic fraction, [-]
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     ,& !input  :root layer nonstructural element, [g d-2]
    THIN_pft                  => plt_distb%THIN_pft                  ,& !input  :thinning of plant population, [-]
    NumStructuralRootAxes_pft      => plt_morph%NumStructuralRootAxes_pft      ,& !input  :number of structural root axes,[-]
    DCORP                     => plt_distb%DCORP                     ,& !input  :soil mixing fraction with tillage, [-]
    FracRootElmAllocm         => plt_allom%FracRootElmAllocm         ,& !input  :C woody fraction in root,[-]
    k_fine_comp               => pltpar%k_fine_comp                  ,& !input  :fine litter complex id
    k_woody_comp              => pltpar%k_woody_comp                 ,& !input  :woody litter complex id
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr  ,& !input  :root layer element secondary axes, [g d-2]
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr  ,& !input  :root layer element primary axes, [g d-2]
    iroot                     => pltpar%iroot                        ,& !input  :group id of plant root litter
    inonstruct                => pltpar%inonstruct                   ,& !input  :group id of plant nonstructural litter
    icwood                    => pltpar%icwood                       ,& !input  :group id of coarse woody litter
    iHarvstType_pft           => plt_distb%iHarvstType_pft           ,& !input  :type of harvest,[-]
    LitrfallElms_pvr          => plt_bgcr%LitrfallElms_pvr           ,& !inoput :plant LitrFall element, [g d-2 h-1]
    RootGasLossDisturb_pft    => plt_bgcr%RootGasLossDisturb_pft     ,& !inoput :gaseous flux fron root disturbance, [g d-2 h-1]
    trcg_rootml_pvr           => plt_rbgc%trcg_rootml_pvr            ,& !inoput :root gas content, [g d-2]
    trcs_rootml_pvr           => plt_rbgc%trcs_rootml_pvr             & !inoput :root aqueous content, [g d-2]
  )
  call PrintInfo('beg '//subname)
  FracLeftThin = 0._r8
  XHVST1       = 1._r8-FracLeftThin
  D3385: DO M=1,jsken
    DO NE=1,NumPlantChemElms
      HarvestedBiomass(NE)=XHVST1*PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ))
      LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+HarvestedBiomass(NE)
    ENDDO

    DO NR=1,NumStructuralRootAxes_pft(NZ)
      if(N==ipltroot)THEN
        DO NE=1,NumPlantChemElms
          HarvestedBiomass(NE)=XHVST1*PlantElmAllocMat4Litr(NE,icwood,M,NZ)*AZMAX1(RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)) &
            *FracRootElmAllocm(NE,k_woody_comp)
          LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+HarvestedBiomass(NE)       
        ENDDO
      ENDIF

      DO NE=1,NumPlantChemElms
        HarvestedBiomass(NE)=XHVST1*PlantElmAllocMat4Litr(NE,icwood,M,NZ)*AZMAX1(RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)) &
          *FracRootElmAllocm(NE,k_woody_comp)
        LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+HarvestedBiomass(NE)       
      ENDDO

      !woody roots
      if(N==ipltroot)THEN
        DO NE=1,NumPlantChemElms
          HarvestedBiomass(NE)=XHVST1*PlantElmAllocMat4Litr(NE,iroot,M,NZ)*AZMAX1(RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)) &
            *FracRootElmAllocm(NE,k_fine_comp)
          LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+HarvestedBiomass(NE)
        ENDDO
      ENDIF

      DO NE=1,NumPlantChemElms
        HarvestedBiomass(NE)=XHVST1*PlantElmAllocMat4Litr(NE,iroot,M,NZ)*AZMAX1(RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)) &
          *FracRootElmAllocm(NE,k_fine_comp)
        LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+HarvestedBiomass(NE)
      ENDDO

    enddo
  ENDDO D3385
  !
  !     RELEASE ROOT GAS CONTENTS DURING HARVESTING
  !
  !     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
  !     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
  !     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
  !
  DO idg=idg_beg,idg_NH3
    RootGasLossDisturb_pft(idg,NZ)=RootGasLossDisturb_pft(idg,NZ)-XHVST1 &
      *(trcg_rootml_pvr(idg,N,L,NZ)+trcs_rootml_pvr(idg,N,L,NZ))
    trcg_rootml_pvr(idg,N,L,NZ)=FracLeftThin*trcg_rootml_pvr(idg,N,L,NZ)
    trcs_rootml_pvr(idg,N,L,NZ)=FracLeftThin*trcs_rootml_pvr(idg,N,L,NZ)
  ENDDO
  call PrintInfo('end '//subname)
  end associate          
  end subroutine RootRemovalL4Annual
!----------------------------------------------------------------------------------------------------
  subroutine TerminateRoots4Annuals(yearIJ,NZ) 
  !
  !!Description
  !terminate roots for annual grasses that terminate and reseed
   
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='TerminateRoots4Annuals'
  real(r8) :: FracLeftThin,XHVST1
  integer :: N,L
  associate(                                                             &
    Myco_pft                   => plt_morph%Myco_pft                    ,& !input  :mycorrhizal type (no or yes),[-]  
    NU                         => plt_site%NU                           ,& !input  :current soil surface layer number, [-]
    MaxNumRootLays             => plt_site%MaxNumRootLays                & !input  :maximum root layer number,[-]
  )
  call PrintInfo('beg '//subname)
  DO N=1,Myco_pft(NZ)
    DO L=NU,MaxNumRootLays
      call RootRemovalL4Annual(yearIJ,N,L,NZ,FracLeftThin,XHVST1)

      call HarvstUpdateRootStateL(yearIJ,N,L,NZ,FracLeftThin,XHVST1,PopulationThinning=.false.)
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine TerminateRoots4Annuals
!----------------------------------------------------------------------------------------------------
  subroutine HarvstUpdateRootStateL(yearIJ,N,L,NZ,FracLeftThin,XHVST1,PopulationThinning)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  

  integer,  intent(in) :: N,L,NZ
  real(r8), intent(in) :: FracLeftThin    !fraction of biomass remaining alive after removal
  real(r8), intent(in) :: XHVST1          !fraction of biomass removed

  logical, intent(in) :: PopulationThinning !whole plants removed, rather than damage to surviving roots
  real(r8) :: FracGeometryLeft              !remaining per-axis/per-plant dimensions

  character(len=*), parameter :: subname='HarvstUpdateRootStateL'
  integer :: NE,NR,M
  associate(                                                          &
    NumStructuralRootAxes_pft => plt_morph%NumStructuralRootAxes_pft ,& !input  :number of structural root axes,[-]
    inonstruct                => pltpar%inonstruct                   ,& !input  :group id of plant nonstructural litter
    iroot                     => pltpar%iroot                        ,& !input  :group id of plant root litter
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft        ,& !input  :N2 fixation type,[-]
    PlantElmAllocMat4Litr     => plt_soilchem%PlantElmAllocMat4Litr  ,& !input  :litter kinetic fraction, [-]
    k_fine_comp               => pltpar%k_fine_comp                  ,& !input  :fine litter complex id
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr  ,& !inoput :root layer element primary axes, [g d-2]
    RootMediumStructElms_rpvr => plt_biom%RootMediumStructElms_rpvr  ,& !inoput :medium root structural C/N/P, [g d-2]
    RootMedStruct_pvr         => plt_biom%RootMedStruct_pvr          ,& !output :layer medium root structural C/N/P, [g d-2]
    RootMediumLength_rpvr     => plt_morph%RootMediumLength_rpvr     ,& !inoput :total medium root length, [m d-2]
    RootMediumXNum_rpvr       => plt_morph%RootMediumXNum_rpvr       ,& !inoput :medium root axis count, [d-2]
    RootMediumXNum_pvr        => plt_morph%RootMediumXNum_pvr        ,& !output :layer medium root axis count, [d-2]
    RootMediumLength_pvr      => plt_morph%RootMediumLength_pvr      ,& !output :layer mean medium root length, [m]
    Root1stActStructElms_rpvr => plt_biom%Root1stActStructElms_rpvr  ,& !inoput :root layer active zone element in primary axes, [g d-2]
    Root1stLigStructElms_rpvr => plt_biom%Root1stLigStructElms_rpvr  ,& !inoput :root layer lignified zone element in primary axes, [g d-2]
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr  ,& !inoput :root layer element secondary axes, [g d-2]
    Root1stLenPP_rpvr         => plt_morph%Root1stLenPP_rpvr         ,& !inoput :primary root length per axis in soil layer, [m]
    Root2ndLen_rpvr           => plt_morph%Root2ndLen_rpvr           ,& !inoput :root layer length secondary axes, [m d-2]
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     ,& !inoput :root layer nonstructural element, [g d-2]
    Root2ndXNum_rpvr          => plt_morph%Root2ndXNum_rpvr          ,& !inoput :root layer number secondary axes, [d-2]
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr           ,& !inoput :root layer protein C, [gC d-2]
    Root1stXNumL_pvr          => plt_morph%Root1stXNumL_pvr          ,& !inoput :root layer number primary axes, [d-2]
    Root2ndXNumL_rpvr         => plt_morph%Root2ndXNumL_rpvr         ,& !inoput :root layer number axes, [d-2]
    RootTotLenPerPlant_pvr    => plt_morph%RootTotLenPerPlant_pvr    ,& !inoput :root layer length per plant, [m p-1]
    RootLenDensPerPlant_pvr   => plt_morph%RootLenDensPerPlant_pvr   ,& !inoput :root layer length density, [m m-3]
    RootPoreVol_pvr           => plt_morph%RootPoreVol_pvr           ,& !inoput :root layer volume air, [m2 d-2]
    RootSAreaPerPlant_pvr     => plt_morph%RootSAreaPerPlant_pvr     ,& !inoput :root layer area per plant, [m p-1]
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr         ,& !inoput :root layer C, [gC d-2]
    RootVH2O_pvr              => plt_morph%RootVH2O_pvr              ,& !inoput :root layer volume water, [m2 d-2]
    RootCO2Autor_pvr          => plt_rbgc%RootCO2Autor_pvr           ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    RootCO2EmisPot_pvr        => plt_rbgc%RootCO2EmisPot_pvr         ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootRespPotent_pvr        => plt_rbgc%RootRespPotent_pvr         ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]
    LitrfallElms_pvr          => plt_bgcr%LitrfallElms_pvr           ,& !inoput :plant LitrFall element, [g d-2 h-1]
    RootNodulStrutElms_rpvr   => plt_biom%RootNodulStrutElms_rpvr    ,& !inoput :root layer nodule element, [g d-2]
    RootNodulNonstElms_rpvr   => plt_biom%RootNodulNonstElms_rpvr    ,& !inoput :root layer nonstructural element, [g d-2]
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr     & !inoput :root layer structural C, [gC d-2]
  )
  !
  !     REDUCE ROOT STATE VARIABLES DURING HARVESTING
  !
  !     RootLenDensPerPlant_pvr,RootTotLenPerPlant_pvr=root length density,root length per plant
  !     RootVH2O_pvr,RootPoreVol_pvr=root or myco aqueous,gaseous volume
  !     RootSAreaPerPlant_pvr=root surface area per plant
  !     RootRespPotent_pvr,RootCO2EmisPot_pvr,RootCO2Autor_pvr unlimited by O2,nonstructural C
  !    
  call PrintInfo('beg '//subname)
  ! Population totals below already account for the removed plants. Do not
  ! shorten the roots of each survivor a second time. Complete removal still
  ! clears geometry; tissue damage retains its layer-specific reduction.
  FracGeometryLeft=FracLeftThin
  IF(PopulationThinning .and. FracLeftThin.GT.0._r8)FracGeometryLeft=1._r8

  if(N.EQ.ipltroot)then
    DO NR=1,NumStructuralRootAxes_pft(NZ)
      DO NE=1,NumPlantChemElms        
        RootMediumStructElms_rpvr(NE,L,NR,NZ) = RootMediumStructElms_rpvr(NE,L,NR,NZ)*FracLeftThin
        Root1stActStructElms_rpvr(NE,L,NR,NZ) = Root1stActStructElms_rpvr(NE,L,NR,NZ)*FracLeftThin
        Root1stLigStructElms_rpvr(NE,L,NR,NZ) = Root1stLigStructElms_rpvr(NE,L,NR,NZ)*FracLeftThin
        RootMyco1stStrutElms_rpvr(NE,L,NR,NZ) = Root1stActStructElms_rpvr(NE,L,NR,NZ)+Root1stLigStructElms_rpvr(NE,L,NR,NZ)
      ENDDO
      !Thin total length and axis count together; surviving root radius is unchanged.
      RootMediumLength_rpvr(L,NR,NZ) = RootMediumLength_rpvr(L,NR,NZ)*FracLeftThin
      RootMediumXNum_rpvr(L,NR,NZ)   = RootMediumXNum_rpvr(L,NR,NZ)*FracLeftThin
      Root1stLenPP_rpvr(L,NR,NZ)  = Root1stLenPP_rpvr(L,NR,NZ)*FracGeometryLeft
    ENDDO
    call plt_morph%RefreshMediumRootMeanLength(NZ)
    !Refresh layer totals once for plant roots, not again for mycorrhizae.
    DO NE=1,NumPlantChemElms
      RootMedStruct_pvr(NE,L,NZ)=SUM(RootMediumStructElms_rpvr(NE,L,1:NumStructuralRootAxes_pft(NZ),NZ))
    ENDDO
    ! Per-group counts can include axes reserved for growth initiation.
    ! Rebuild transport totals from surviving medium-root tissue only.
    RootMediumXNum_pvr(L,NZ)=0._r8
    RootMediumLength_pvr(L,NZ)=0._r8
    DO NR=1,NumStructuralRootAxes_pft(NZ)
      IF(RootMediumStructElms_rpvr(ielmc,L,NR,NZ).GT.0._r8 .and. &
         RootMediumLength_rpvr(L,NR,NZ).GT.0._r8 .and. &
         RootMediumXNum_rpvr(L,NR,NZ).GT.0._r8)THEN
        RootMediumXNum_pvr(L,NZ)=RootMediumXNum_pvr(L,NZ)+RootMediumXNum_rpvr(L,NR,NZ)
        RootMediumLength_pvr(L,NZ)=RootMediumLength_pvr(L,NZ)+RootMediumLength_rpvr(L,NR,NZ)
      ENDIF
    ENDDO
    IF(RootMediumXNum_pvr(L,NZ).GT.0._r8)THEN
      RootMediumLength_pvr(L,NZ)=RootMediumLength_pvr(L,NZ)/RootMediumXNum_pvr(L,NZ)
    ENDIF
    Root1stXNumL_pvr(L,NZ)        = Root1stXNumL_pvr(L,NZ)*FracLeftThin      
  ENDIF

  D3960: DO NR=1,NumStructuralRootAxes_pft(NZ)
    DO NE=1,NumPlantChemElms
      RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ) = RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)*FracLeftThin
    ENDDO
    Root2ndLen_rpvr(N,L,NR,NZ)  = Root2ndLen_rpvr(N,L,NR,NZ)*FracLeftThin
    Root2ndXNum_rpvr(N,L,NR,NZ) = Root2ndXNum_rpvr(N,L,NR,NZ)*FracLeftThin
  ENDDO D3960
  !
  DO NE=1,NumPlantChemElms
    RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)*FracLeftThin
  ENDDO

  RootMycoActiveBiomC_pvr(N,L,NZ) = RootMycoActiveBiomC_pvr(N,L,NZ)*FracLeftThin
  PopuRootMycoC_pvr(N,L,NZ)       = PopuRootMycoC_pvr(N,L,NZ)*FracLeftThin
  RootProteinC_pvr(N,L,NZ)        = RootProteinC_pvr(N,L,NZ)*FracLeftThin
  Root2ndXNumL_rpvr(N,L,NZ)       = Root2ndXNumL_rpvr(N,L,NZ)*FracLeftThin
  RootTotLenPerPlant_pvr(N,L,NZ)  = RootTotLenPerPlant_pvr(N,L,NZ)*FracGeometryLeft
  RootLenDensPerPlant_pvr(N,L,NZ) = RootLenDensPerPlant_pvr(N,L,NZ)*FracGeometryLeft
  RootPoreVol_pvr(N,L,NZ)         = RootPoreVol_pvr(N,L,NZ)*FracLeftThin
  RootVH2O_pvr(N,L,NZ)            = RootVH2O_pvr(N,L,NZ)*FracLeftThin
  RootSAreaPerPlant_pvr(N,L,NZ)   = RootSAreaPerPlant_pvr(N,L,NZ)*FracGeometryLeft
  RootRespPotent_pvr(N,L,NZ)      = RootRespPotent_pvr(N,L,NZ)*FracLeftThin
  RootCO2EmisPot_pvr(N,L,NZ)      = RootCO2EmisPot_pvr(N,L,NZ)*FracLeftThin
  !
  !     NODULE LitrFall AND STATE VARIABLES DURING HARVESTING
  !
  !     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
  !     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
  !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
  !     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
  !     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
  !
  IF(N.EQ.ipltroot .and. is_plant_N2fix(iPlantNfixType_pft(NZ)) .AND. is_root_N2fix(iPlantNfixType_pft(NZ)))THEN
    DO NE=1,NumPlantChemElms
      D3395: DO M=1,jsken
        LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+ &
          XHVST1*AZMAX1(PlantElmAllocMat4Litr(NE,iroot,M,NZ)*RootNodulStrutElms_rpvr(NE,L,NZ) &
          +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*RootNodulNonstElms_rpvr(NE,L,NZ))
      ENDDO D3395
      RootNodulStrutElms_rpvr(NE,L,NZ) = RootNodulStrutElms_rpvr(NE,L,NZ)*FracLeftThin
      RootNodulNonstElms_rpvr(NE,L,NZ) = RootNodulNonstElms_rpvr(NE,L,NZ)*FracLeftThin
    ENDDO
  ENDIF
  call PrintInfo('end '//subname)
  end associate          
  end subroutine HarvstUpdateRootStateL
  ![tail]
end module PlantDisturbsMod
