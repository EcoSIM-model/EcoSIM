module grosubsMod
!!
! Description:
! module for plant biological transformations
  use minimathmod,         only: safe_adb, AZMAX1,real_truncate
  use data_kind_mod,       only: r8 => DAT_KIND_R8,yearIJ_type
  use EcoSIMCtrlMod,       only: lverb,lcoarseroot
  use EcoSiMParDataMod,    only: pltpar
  use RootMod,             only: RootBGCModel
  use PlantNonstElmDynMod, only: PlantNonstElmTransfer
  use PlantDebugMod 
  use DebugToolMod,        only: PrintInfo
  use EcosimConst
  use PlantBGCPars
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantPhotosynthesisAPIData, only : plt_photo
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantAllometryAPIData, only : plt_allom
  use PlantBiomassAPIData, only : plt_biom
  use PlantEnergyWaterAPIData, only : plt_ew
  use PlantDisturbanceAPIData, only : plt_distb
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantRootBGCAPIData, only : plt_rbgc
  use PhotoSynsMod
  use PlantMathFuncMod
  use LitterFallMod
  use PlantBranchMod
  use PlantBalMod
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__
! DIMENSION VCO2(400,366,05)
!
!     RTSK=relative primary root sink strength 0.25=shallow,4.0=deep root profile
!     FXRN=rate constant for plant-bacteria nonstructl C,N,P exchange (h-1)
!     RateK4ShootSeaStoreNonstEXfer=rate constant for leaf-storage nonstructl C,N,P exchange (h-1)
!     RateK4RootSeaStorNonstEXfer=rate constant for root-storage nonstructl C,N,P exchange (h-1)
!     FPART=parameter for allocating nonstructural C to shoot organs
!     FXSH,FXRT=shoot-root partitioning of storage C during leafout
!     FRSV=rate constant for remobiln of storage C,N,P during leafout (h-1)
!     FXFY,FXFZ=rate constant for leaf-reserve nonstructural C,N,P exchange (h-1)
!     EFIRE=combustion  of N,P relative to C
!     PSILY=canopy water potential below which leafoff is induced (MPa)
!     HourFailGrainFill_brchY=number of hours after physiol maturity required for senescence
!     HourReq2InitSStor4LeafOut=number of hours required to initiate remobilization of storage C for leafout
!     GVMX=specific oxidation rate of nonstructural C during leafout at 25 C
!

  integer  :: curday,curhour
  public :: GrowPlants

  contains
  ![header]

!----------------------------------------------------------------------------------------------------
  subroutine GrowPlants(yearIJ)
!
!     THIS subroutine CALCULATES ALL PLANT BIOLOGICAL TRANSFORMATIONS
!
  use PlantDisturbsMod, only : RemoveBiomassByDisturbance
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  character(len=*), parameter :: subname='GrowPlants'
  
  real(r8) :: CanopyHeight_copy(JP1)
  integer :: L,K,M
  integer :: NZ,NE
  real(r8) :: tvegE(NumPlantChemElms),tvegE1(NumPlantChemElms)

  ! begin_execution
  associate(                                                               &
    IsPlantActive_pft            => plt_pheno%IsPlantActive_pft           ,& !input  :flag for living pft, [-]
    NP                           => plt_site%NP                           ,& !input  :current number of plant species,[-]
    NP0                          => plt_site%NP0                          ,& !input  :intitial number of plant species,[-]
    ZEROS                        => plt_site%ZEROS                        ,& !input  :threshold zero for numerical stability,[-]    
    PlantPopuLive_pft            => plt_site%PlantPopuLive_pft            ,& !input  :plant population, [d-2]    
    StalkHeight_pft              => plt_morph%StalkHeight_pft             ,& !inoput :stalk height/length, [m]
    CanopyHeightLive_pft         => plt_morph%CanopyHeightLive_pft         & !inoput :canopy height, [m]
  )
  call PrintInfo('beg '//subname)
  !     TOTAL AGB FOR GRAZING IN LANDSCAPE SECTION
  !
  !
  !     INITIALIZE SENESCENCE ARRAYS
  DO NZ=1,NP0
    CanopyHeight_copy(NZ)                  = CanopyHeightLive_pft(NZ)
    CanopyHeightLive_pft(NZ)               = 0._r8
    StalkHeight_pft(NZ)                    = 0._r8    
    plt_rbgc%canopy_growth_pft(NZ)         = 0._r8
    plt_biom%RootMycoMassElm_pvr(:,:,:,NZ) = 0._r8
  ENDDO  
  !
  !     TRANSFORMATIONS IN LIVING PLANT POPULATIONS
  !
  D9985: DO NZ=1,NP
    !      
    IF(IsPlantActive_pft(NZ).EQ.iTrue .and. PlantPopuLive_pft(NZ).GT.ZEROS)THEN      
      !
      call GrowOnePlant(yearIJ,NZ,CanopyHeight_copy)

      call RemoveBiomassByDisturbance(yearIJ,NZ)

      !   RESET DEAD BRANCHES
      call ResetDeadPlant(yearIJ,NZ)
      
    ENDIF  
    
    call AccumulateStates(yearIJ,NZ)
    !call CheckPlantBalanceZ(yearIJ,NZ,header=subname)
  ENDDO D9985
  !
  call Live2DeadTransformation(yearIJ)     
  !call CheckPlantBalanceZ(yearIJ,NZ=1,header=subname//'L2DEAD')  
  !
  call PrintInfo('end '//subname)
  end associate
  END subroutine GrowPlants

!----------------------------------------------------------------------------------------------------
  subroutine Live2DeadTransformation(yearIJ)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  character(len=*), parameter :: subname='Live2DeadTransformation'
  integer :: L,K,NZ,M,NE,NB
  real(r8) :: XFRC,XFRN,XFRP,XFRE
  real(r8), parameter :: StandingDeadKd_decid=1.5814E-04_r8  !first-order decay rate of deciduous standing dead biomass to litter
  real(r8), parameter :: StandingDeadKd_other=1.5814E-05_r8  !first-order decay rate of non-deciduous standing dead biomass to litter
  real(r8) :: StandingDeadKd
!     begin_execution

  associate(                                                                    &
    CH4ByFire_CumYr_pft            => plt_distb%CH4ByFire_CumYr_pft            ,& !input  :plant CH4 emission from fire, [g d-2 ]
    CO2ByFire_CumYr_pft            => plt_distb%CO2ByFire_CumYr_pft            ,& !input  :plant CO2 emission from fire, [g d-2 ]
    CumSoilThickness_vr            => plt_site%CumSoilThickness_vr             ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    EcoHavstElmntCum_pft           => plt_distb%EcoHavstElmntCum_pft           ,& !input  :total plant element harvest, [gC d-2 ]
    GrossCO2Fix_pft                => plt_bgcr%GrossCO2Fix_pft                 ,& !input  :total gross CO2 fixation, [gC d-2 ]
    GrossResp_pft                  => plt_bgcr%GrossResp_pft                   ,& !input  :total plant respiration, [gC d-2 ]
    HourReq4LeafOut_brch           => plt_pheno%HourReq4LeafOut_brch           ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    Hours4Leafout_brch             => plt_pheno%Hours4Leafout_brch             ,& !input  :heat requirement for spring leafout/dehardening, [h]
    IsPlantActive_pft              => plt_pheno%IsPlantActive_pft              ,& !input  :flag for living pft, [-]
    MaxNumRootLays                 => plt_site%MaxNumRootLays                  ,& !input  :maximum root layer number,[-]
    N2ObyFire_CumYr_pft            => plt_distb%N2ObyFire_CumYr_pft            ,& !input  :plant N2O emission from fire, [g d-2 ]
    NH3Emis_CumYr_pft              => plt_bgcr%NH3Emis_CumYr_pft               ,& !input  :total canopy NH3 flux, [gN d-2 ]
    NH3byFire_CumYr_pft            => plt_distb%NH3byFire_CumYr_pft            ,& !input  :plant NH3 emission from fire, [g d-2 ]
    NP0                            => plt_site%NP0                             ,& !input  :intitial number of plant species,[-]
    NP                             => plt_site%NP                              ,& !input  :current number of plant species,[-]
    NU                             => plt_site%NU                              ,& !input  :current soil surface layer number, [-]
    TKS_vr                         => plt_ew%TKS_vr                            ,& !input  :mean annual soil temperature, [K]
    TempOffset_pft                 => plt_pheno%TempOffset_pft                 ,& !input  :adjustment of Arhhenius curves for plant thermal acclimation, [oC]
    NetCumElmntFlx2Plant_pft       => plt_pheno%NetCumElmntFlx2Plant_pft       ,& !input  :effect of canopy element status on seed set, [-]
    NumOfBranches_pft              => plt_morph%NumOfBranches_pft              ,& !input  :number of branches,[-]
    PO4byFire_CumYr_pft            => plt_distb%PO4byFire_CumYr_pft            ,& !input  :plant PO4 emission from fire, [g d-2 ]
    PlantExudElm_CumYr_pft         => plt_rbgc%PlantExudElm_CumYr_pft          ,& !input  :total net root element uptake (+ve) - exudation (-ve), [gC d-2 ]
    PlantN2Fix_CumYr_pft           => plt_bgcr%PlantN2Fix_CumYr_pft            ,& !input  :total plant N2 fixation, [g d-2 ]
    RootElms_pft                   => plt_biom%RootElms_pft                    ,& !input  :plant root element mass, [g d-2]
    SeasonalNonstElms_pft          => plt_biom%SeasonalNonstElms_pft           ,& !input  :plant stored nonstructural element at current step, [g d-2]
    ShootElms_pft                  => plt_biom%ShootElms_pft                   ,& !input  :canopy shoot structural chemical element mass, [g d-2]
    EnablePlantLeafOut_brch        => plt_pheno%EnablePlantLeafOut_brch        ,& !input  :branch phenology flag, [-]
    iPlant2ndGrothPattern_pft      => plt_pheno%iPlant2ndGrothPattern_pft      ,& !input  :plant expression of secondary growth, [-]                
    fTCanopyGroth_pft              => plt_pheno%fTCanopyGroth_pft              ,& !inoput :canopy temperature growth function, [-]
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft          ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlantTurnoverPattern_pft      => plt_pheno%iPlantTurnoverPattern_pft      ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    iYearCurrent                   => plt_site%iYearCurrent                    ,& !input  :current year,[-]
    k_fine_comp                    => pltpar%k_fine_comp                       ,& !input  :fine litter complex id
    k_woody_comp                   => pltpar%k_woody_comp                      ,& !input  :woody litter complex id
    SurfLitrfallElms_pft           => plt_bgcr%SurfLitrfallElms_pft            ,& !output :surface litterfall
    LitrfalStrutElms_CumYr_pft     => plt_bgcr%LitrfalStrutElms_CumYr_pft      ,& !inoput :total plant element LitrFall, [g d-2 ]
    LitrfallElms_pft               => plt_bgcr%LitrfallElms_pft                ,& !inoput :plant element LitrFall, [g d-2 h-1]
    LitrfallElms_pvr               => plt_bgcr%LitrfallElms_pvr                ,& !inoput :plant LitrFall element, [g d-2 h-1]
    StandDeadCompKElms_pft         => plt_biom%StandDeadCompKElms_pft          ,& !inoput :standing dead element fraction, [g d-2]
    SurfLitrfalStrutElms_CumYr_pft => plt_bgcr%SurfLitrfalStrutElms_CumYr_pft  ,& !inoput :total surface LitrFall element, [g d-2]
    doInitPlant_pft                => plt_pheno%doInitPlant_pft                ,& !inoput :PFT initialization flag:0=no,1=yes,[-]
    NetPrimProduct_pft             => plt_bgcr%NetPrimProduct_pft              ,& !output :total net primary productivity, [gC d-2]
    PlantinDepz_pft                => plt_morph%PlantinDepz_pft                ,& !output :planting depth, [m]
    iDayPlanting_pft               => plt_distb%iDayPlanting_pft               ,& !output :day of planting,[-]
    iYearPlanting_pft              => plt_distb%iYearPlanting_pft               & !output :year of planting,[-]
  )
  call PrintInfo('beg '//subname)
  D9975: DO NZ=1,NP0
    if(NZ.GT.NP)then
      !fallow slot: residue only, so the canopy temperature function is not
      !updated by the uptake model; use the soil surface temperature instead
      fTCanopyGroth_pft(NZ)=calc_canopy_grow_tempf(real_truncate(TKS_vr(NU)+TempOffset_pft(NZ),1.e-3_r8))
    endif
    if(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. .NOT.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))then
      StandingDeadKd=StandingDeadKd_decid
    else
      StandingDeadKd=StandingDeadKd_other
    endif
    !
    !     ACTIVATE DORMANT SEEDS
    !
    IF(doInitPlant_pft(NZ).EQ.itrue)THEN
      D205: DO NB=1,NumOfBranches_pft(NZ)      
        IF(EnablePlantLeafOut_brch(NB,NZ).EQ.iTrue .AND. Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN
          iDayPlanting_pft(NZ)  = yearIJ%I
          iYearPlanting_pft(NZ) = iYearCurrent
          PlantinDepz_pft(NZ)   = 0.005_r8+CumSoilThickness_vr(0)
          doInitPlant_pft(NZ)   = ifalse   !mark plant as initialized
          exit  
        ENDIF  
      ENDDO D205
    ENDIF
    !
    !     LitrFall FROM STANDING DEAD (i.e. loss of standing dead)
    !
    !     XFRC,XFRN,XFRP=LitrFall from standing dead
    !     fTCanopyGroth_pft=temperature function for canopy growth
    !    
    D6235: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        XFRE=StandingDeadKd*fTCanopyGroth_pft(NZ)*StandDeadCompKElms_pft(NE,M,NZ)
        
        IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. .not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))THEN
          !grass/bryophyte
          LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+XFRE
        ELSE
          !trees
          LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)+XFRE
        ENDIF
        StandDeadCompKElms_pft(NE,M,NZ)=StandDeadCompKElms_pft(NE,M,NZ)-XFRE
      ENDDO
    ENDDO D6235
    !
    !     ACCUMULATE TOTAL SURFACE, SUBSURFACE LitrFall
    !
    !   diagnose surface literfall 
    !
    SurfLitrfallElms_pft(:,NZ)=0._r8
    DO K=1,pltpar%NumOfPlantLitrCmplxs      
      D6431: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          SurfLitrfallElms_pft(NE,NZ)           = SurfLitrfallElms_pft(NE,NZ)+LitrfallElms_pvr(NE,M,K,0,NZ)
          SurfLitrfalStrutElms_CumYr_pft(NE,NZ) = SurfLitrfalStrutElms_CumYr_pft(NE,NZ)+LitrfallElms_pvr(NE,M,K,0,NZ)
        ENDDO
      ENDDO D6431  
    ENDDO
    !
    !   the following purposely starts from layer 0.
    !
    LitrfallElms_pft(1:NumPlantChemElms,NZ)=0._r8
    D8955: DO L=0,MaxNumRootLays
      DO K=1,pltpar%NumOfPlantLitrCmplxs      
        D6430: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitrfallElms_pft(NE,NZ)=LitrfallElms_pft(NE,NZ)+LitrfallElms_pvr(NE,M,K,L,NZ)
          enddo         
        ENDDO D6430      
      ENDDO
    ENDDO D8955    

    DO NE=1,NumPlantChemElms
      LitrfalStrutElms_CumYr_pft(NE,NZ)=LitrfalStrutElms_CumYr_pft(NE,NZ)+LitrfallElms_pft(NE,NZ)
    ENDDO

    !   GrossResp_pft < 0 respired into atmosphere
    NetPrimProduct_pft(NZ) = GrossCO2Fix_pft(NZ)+GrossResp_pft(NZ)

  ENDDO D9975
  call PrintInfo('end '//subname)
  end associate
  end subroutine Live2DeadTransformation

!----------------------------------------------------------------------------------------------------
  subroutine GrowOnePlant(yearIJ,NZ,CanopyHeight_copy)
  !
  !Description
  !plant growth
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), intent(in) :: CanopyHeight_copy(JP1)

  character(len=*), parameter :: subname='GrowOnePlant'
  real(r8)  :: CanopyN2Fix_pft(JP1)
  integer  :: NB,NE,KK
  integer  :: BegRemoblizeBranch !phenological eligibility of the current branch
  integer  :: BegRemoblizePlant  !main-branch eligibility controls root seasonal storage
  real(r8) :: TFN6_vr(JZ1)
  real(r8) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW
  real(r8) :: PTRT    !main branch growth allocated to leaf and PetolSheth
  real(r8) :: TFN5
  real(r8) :: WaterStress4Groth
  real(r8) :: Stomata_Stress
  real(r8) :: TurgEff4LeafPetolExpansion,TurgEff4CanopyResp
  real(r8) :: RootSinkC_vr(pltpar%jroots,JZ1),RootSinkC(pltpar%jroots)
  real(r8) :: Root1stSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8) :: tvegE_beg(NumPlantChemElms)
  real(r8) :: tvegE_end(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: canopyE_finale(NumPlantChemElms)
  real(r8) :: err,arr(12),arr1(11),GrothPART2LeafPetole
! begin_execution
  associate(                                                   &
    NumOfBranches_pft     => plt_morph%NumOfBranches_pft      ,& !input  :number of branches,[-]
    PlantPopuLive_pft     => plt_site%PlantPopuLive_pft       ,& !input  :plant population, [d-2]
    ZERO4Groth_pft        => plt_biom%ZERO4Groth_pft          ,& !input  :threshold zero for plang growth calculation, [-]
    isPlantRootAlive_pft  => plt_pheno%isPlantRootAlive_pft   ,& !input  :flag to detect root system death,[-]
    NU                    => plt_site%NU                      ,& !input  :current soil surface layer number, [-]
    MaxSoilLays4Root_pft  => plt_morph%MaxSoilLays4Root_pft   ,& !input  :maximum soil layer number for all root axes,[-]
    MainBranchNum_pft     => plt_morph%MainBranchNum_pft      ,& !input  :number of main branch,[-]    
    isPlantShootAlive_pft => plt_pheno%isPlantShootAlive_pft   & !input  :flag to detect canopy death,[-]
  )

  call PrintInfo('beg '//subname)

  IF(isPlantShootAlive_pft(NZ).EQ.iTrue .OR. isPlantRootAlive_pft(NZ).EQ.iTrue .and. PlantPopuLive_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
    BegRemoblizePlant = ifalse
    PTRT = 0._r8
    
    call StagePlantForGrowth(yearIJ%I,yearIJ%J,NZ,TFN6_vr,CNLFW,CPLFW,&
      CNSHW,CPSHW,CNRTW,CPRTW,TFN5,WaterStress4Groth,Stomata_Stress,TurgEff4LeafPetolExpansion,TurgEff4CanopyResp)
    !
    !     CALCULATE GROWTH OF EACH BRANCH
    DO  NB=1,NumOfBranches_pft(NZ)
      !
      call GrowOneBranch(yearIJ,NB,NZ,TFN6_vr,CanopyHeight_copy,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,&
        TFN5,WaterStress4Groth,Stomata_Stress,TurgEff4LeafPetolExpansion,TurgEff4CanopyResp,&
        GrothPART2LeafPetole,BegRemoblizeBranch)
      !
      ! Use the main branch as the representative plant phenology for root storage.
      ! Other branches retain their own eligibility for branch-local transfers.
      IF(NB.EQ.MainBranchNum_pft(NZ))THEN
        PTRT = GrothPART2LeafPetole
        BegRemoblizePlant = BegRemoblizeBranch
      ENDIF
    ENDDO
 
    call RootBGCModel(yearIJ,NZ,TFN6_vr,CNRTW,CPRTW,RootSinkC_vr,RootSinkC)

    call PlantNonstElmTransfer(yearIJ%I,yearIJ%J,NZ,PTRT,RootSinkC_vr,RootSinkC,BegRemoblizePlant)

  else
    plt_morph%RootSinkWeight_pvr(NU:MaxSoilLays4Root_pft(NZ),NZ)=0._r8   
  ENDIF

  !
  call PrintInfo('end '//subname)  
  end associate
  end subroutine GrowOnePlant

!----------------------------------------------------------------------------------------------------
  subroutine StagePlantForGrowth(I,J,NZ,TFN6_vr,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,&
    TFN5,WaterStress4Groth,Stomata_Stress,TurgEff4LeafPetolExpansion,TurgEff4CanopyResp)
  integer, intent(in) :: I,J,NZ
  REAL(R8), INTENT(OUT) :: TFN6_vr(JZ1)
  REAL(R8), INTENT(OUT) :: CNLFW,CPLFW                 !NC, PC mass ratio of leaf growth
  real(r8), intent(out) :: CNSHW,CPSHW                 !NC, PC mass ratio of shoot growth
  real(r8), intent(out) :: CNRTW,CPRTW                 !NC, PC mass ratio of root growth
  real(r8), intent(out) :: TFN5
  real(r8), intent(out) :: WaterStress4Groth
  real(r8), intent(out) :: Stomata_Stress
  real(r8), intent(out) :: TurgEff4LeafPetolExpansion  !turgor expansion,extension function
  real(r8), intent(out) :: TurgEff4CanopyResp          !expansion,extension function of canopy water potential on respiraiton
  integer :: L,NR,N,NE
  real(r8) :: ACTVM,RTK,STK,TKCM,TKSM,RadiusBulk,RadiusMean
  real(r8), parameter :: dscal=0.99999_r8
  character(len=*), parameter :: subname='StagePlantForGrowth'
!     begin_execution

  associate(                                                                      &
    rNCLeaf_pft                    => plt_allom%rNCLeaf_pft                      ,& !input  :maximum leaf N:C ratio, [g g-1]
    rNCSheath_pft                  => plt_allom%rNCSheath_pft                    ,& !input  :sheath N:C ratio, [g g-1]
    rPCLeaf_pft                    => plt_allom%rPCLeaf_pft                      ,& !input  :maximum leaf P:C ratio, [g g-1]
    rPCSheath_pft                  => plt_allom%rPCSheath_pft                    ,& !input  :sheath P:C ratio, [g g-1]
    CanopySapwoodC_pft             => plt_biom%CanopySapwoodC_pft                ,& !input  :canopy active stalk C, [g d-2]
    Myco_pft                       => plt_morph%Myco_pft                         ,& !input  :mycorrhizal type (no or yes),[-]
    MaxNumRootLays                 => plt_site%MaxNumRootLays                    ,& !input  :maximum root layer number,[-]
    NK                             => plt_site%NK                                ,& !input  :current hydrologically active layer, [-]
    NU                             => plt_site%NU                                ,& !input  :current soil surface layer number, [-]
    PSICanopyTurg_pft              => plt_ew%PSICanopyTurg_pft                   ,& !input  :plant canopy turgor water potential, [MPa]
    PSICanopy_pft                  => plt_ew%PSICanopy_pft                       ,& !input  :canopy total water potential, [Mpa]
    PlantPopuLive_pft              => plt_site%PlantPopuLive_pft                 ,& !input  :plant population, [d-2]
    iPlantPhenolPattern_pft        => plt_pheno%iPlantPhenolPattern_pft          ,& !input  :plant growth habit: annual or perennial,[-]    
    RCS_pft                        => plt_photo%RCS_pft                          ,& !input  :e-folding turgor pressure for stomatal resistance, [MPa]
    RootElms_pft                   => plt_biom%RootElms_pft                      ,& !input  :plant root element mass, [g d-2]
    rNCRoot_pft                    => plt_allom%rNCRoot_pft                      ,& !input  :root N:C ratio, [gN gC-1]
    rPCRootr_pft                   => plt_allom%rPCRootr_pft                     ,& !input  :root P:C ratio, [gP gC-1]
    StalkStrutElms_pft             => plt_biom%StalkStrutElms_pft                ,& !input  :canopy stalk structural element mass, [g d-2]
    TKC_pft                        => plt_ew%TKC_pft                             ,& !input  :canopy temperature, [K]
    TKS_vr                         => plt_ew%TKS_vr                              ,& !input  :mean annual soil temperature, [K]
    TempOffset_pft                 => plt_pheno%TempOffset_pft                   ,& !input  :adjustment of Arhhenius curves for plant thermal acclimation, [oC]
    NumStructuralRootAxes_pft      => plt_morph%NumStructuralRootAxes_pft        ,& !input  :number of structural root axes per plant,[-]    
    ZERO4Groth_pft                 => plt_biom%ZERO4Groth_pft                    ,& !input  :threshold zero for plang growth calculation, [-]
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft            ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlantTurnoverPattern_pft      => plt_pheno%iPlantTurnoverPattern_pft        ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    rNCStalk_pft                   => plt_allom%rNCStalk_pft                     ,& !input  :stalk N:C ratio, [gN gC-1]
    rPCStalk_pft                   => plt_allom%rPCStalk_pft                     ,& !input  :stalk P:C ratio, [g g-1]
    k_fine_comp                    => pltpar%k_fine_comp                         ,& !input  :fine litter complex id
    iPlant2ndGrothPattern_pft      => plt_pheno%iPlant2ndGrothPattern_pft        ,& !input  :plant expression of secondary growth, [-]                
    k_woody_comp                   => pltpar%k_woody_comp                        ,& !input  :woody litter complex id
    Root1stRadius_pvr              => plt_morph%Root1stRadius_pvr                ,& !input  :root layer radius primary axes, [m]      
    FracRootElmAllocm              => plt_allom%FracRootElmAllocm                ,& !inoput :C woody fraction in root,[-]
    FracWoodStalkElmAlloc2Litr     => plt_allom%FracWoodStalkElmAlloc2Litr       ,& !inoput :woody element allocation,[-]
    FracLeafShethElmAlloc2Litr     => plt_allom%FracLeafShethElmAlloc2Litr       ,& !inoput :leave-sheath element allocation, [-]
    FracPetolShethAlloc2Litr       => plt_allom%FracPetolShethAlloc2Litr         ,& !inoput :leaf element allocation,[-]
    RootBiomCPerPlant_pft          => plt_biom%RootBiomCPerPlant_pft             ,& !inoput :root C biomass per plant, [g p-1]
    TurgEff4CanopyResp_pft         => plt_rbgc%TurgEff4CanopyResp_pft            ,& !output :Turg pressure effect on canopy respiration, [-]
    CanopyLeafAreaZ_pft            => plt_morph%CanopyLeafAreaZ_pft              ,& !output :canopy layer leaf area, [m2 d-2]
    CanopyLeafCLyr_pft             => plt_biom%CanopyLeafCLyr_pft                ,& !output :canopy layer leaf C, [g d-2]
    CanopyStemSurfAreaZ_pft        => plt_morph%CanopyStemSurfAreaZ_pft          ,& !output :plant canopy layer stem surface area, [m2 d-2]
    Root1stXNumL_pvr               => plt_morph%Root1stXNumL_pvr                 ,& !output :root layer number primary axes, [d-2]
    RootMediumXNum_pvr             => plt_morph%RootMediumXNum_pvr               ,& !inoput :Number of medium size root axes in layer, [d-2]    
    StalkAveRadius_pft             => plt_morph%StalkAveRadius_pft               ,& !input  :main stalk radius,[m]        
    Num1stAxesPerStructRootX_pft   => plt_morph%Num1stAxesPerStructRootX_pft     ,& !output :number of primary root axes per structural root axis, [d-2]
    Num1stAxesPerStructRootXPOP_pft=> plt_morph%Num1stAxesPerStructRootXPOP_pft  ,& !output :population primary root axes number on one structrual axis, [d-2]    
    Num1stRootAxesPP_pft           => plt_morph%Num1stRootAxesPP_pft             ,& !output :number of primary root axesr per plant, [d-2]    
    Root2ndXNumL_rpvr              => plt_morph%Root2ndXNumL_rpvr                ,& !output :root layer number axes, [d-2]
    RootCO2Autor_pvr               => plt_rbgc%RootCO2Autor_pvr                  ,& !output :root respiration constrained by O2, [g d-2 h-1]
    RootCO2EmisPot_pvr             => plt_rbgc%RootCO2EmisPot_pvr                ,& !output :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootN2Fix_pvr                  => plt_bgcr%RootN2Fix_pvr                     ,& !output :root N2 fixation, [gN d-2 h-1]
    RootProteinC_pvr               => plt_biom%RootProteinC_pvr                  ,& !output :root layer protein C, [gC d-2]
    RootRespPotent_pvr             => plt_rbgc%RootRespPotent_pvr                 & !output :root respiration unconstrained by O2, [g d-2 h-1]
  )
  call PrintInfo('beg '//subname)
  StalkAveRadius_pft(NZ) = 0._r8
  D2: DO L=1,NumCanopyLayers1
    CanopyLeafAreaZ_pft(L,NZ)     = 0._r8
    CanopyLeafCLyr_pft(L,NZ)      = 0._r8
    CanopyStemSurfAreaZ_pft(L,NZ) = 0._r8
  ENDDO D2


  D6: DO L=1,NK
    RootMediumXNum_pvr(L,NZ) = 0._r8
    Root1stXNumL_pvr(L,NZ)      = 0._r8  
    D9: DO N=1,Myco_pft(NZ)    
      RootProteinC_pvr(N,L,NZ)   = 0._r8
      Root2ndXNumL_rpvr(N,L,NZ)  = 0._r8
      RootRespPotent_pvr(N,L,NZ) = 0._r8
      RootCO2EmisPot_pvr(N,L,NZ) = 0._r8
      RootCO2Autor_pvr(N,L,NZ)   = 0._r8
    ENDDO D9
    RootN2Fix_pvr(L,NZ)          = 0._r8
  ENDDO D6

!
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+PetolSheth,2=none,3=between 1,2
!     WTSTK,WVSTK=stalk,sapwood mass
!     FWOOD,FWODB=C woody fraction in stalk,other organs:0=woody,1=non-woody
!     CN*,CP*=N:C,P:C ratios in plant organs from PFT files
!     CN*W,CP*W=N:C,P:C ratios in plant organs weighted for wood content
!     *LF=leaf,*SHE=PetolSheth,*STK=stalk,*RT=root
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     FWODSN,FWODSP=N,P woody fraction in PetolSheth:0=woody,1=non-woody
!     FWOODN,FWOODP=N,P woody fraction in stalk:0=woody,1=non-woody
!
  IF(iPlantTurnoverPattern_pft(NZ).EQ.0                      & !'Rapid, like deciduous tree '
    .OR. (.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ))) &
    .OR. StalkStrutElms_pft(ielmc,NZ).LE.ZERO4Groth_pft(NZ))THEN
    !not tree
    FracLeafShethElmAlloc2Litr(ielmc,k_fine_comp)     = 1.0_r8
    FracWoodStalkElmAlloc2Litr(ielmc,k_fine_comp) = 1.0_r8
    FracRootElmAllocm(ielmc,k_fine_comp)          = 1.0_r8
  ELSE
    !tree
    FracLeafShethElmAlloc2Litr(ielmc,k_fine_comp) = 1.0_r8
    FracWoodStalkElmAlloc2Litr(ielmc,k_fine_comp) = AMIN1(SQRT(CanopySapwoodC_pft(NZ)/StalkStrutElms_pft(ielmc,NZ)),1._r8)
    FracRootElmAllocm(ielmc,k_fine_comp)          = AMIN1((CanopySapwoodC_pft(NZ)/StalkStrutElms_pft(ielmc,NZ))**(1._r8/6),1._R8)
  ENDIF

  FracLeafShethElmAlloc2Litr(ielmc,k_woody_comp) = AZMAX1(1.0_r8-FracLeafShethElmAlloc2Litr(ielmc,k_fine_comp))
  FracWoodStalkElmAlloc2Litr(ielmc,k_woody_comp) = AZMAX1(1.0_r8-FracWoodStalkElmAlloc2Litr(ielmc,k_fine_comp))
  FracRootElmAllocm(ielmc,k_woody_comp)          = AZMAX1(1.0_r8-FracRootElmAllocm(ielmc,k_fine_comp))

  CNLFW=FracLeafShethElmAlloc2Litr(ielmc,k_woody_comp)*rNCStalk_pft(NZ)+FracLeafShethElmAlloc2Litr(ielmc,k_fine_comp)*rNCLeaf_pft(NZ)
  CPLFW=FracLeafShethElmAlloc2Litr(ielmc,k_woody_comp)*rPCStalk_pft(NZ)+FracLeafShethElmAlloc2Litr(ielmc,k_fine_comp)*rPCLeaf_pft(NZ)
  CNSHW=FracLeafShethElmAlloc2Litr(ielmc,k_woody_comp)*rNCStalk_pft(NZ)+FracLeafShethElmAlloc2Litr(ielmc,k_fine_comp)*rNCSheath_pft(NZ)
  CPSHW=FracLeafShethElmAlloc2Litr(ielmc,k_woody_comp)*rPCStalk_pft(NZ)+FracLeafShethElmAlloc2Litr(ielmc,k_fine_comp)*rPCSheath_pft(NZ)
  CNRTW=FracRootElmAllocm(ielmc,k_woody_comp)*rNCStalk_pft(NZ)+FracRootElmAllocm(ielmc,k_fine_comp)*rNCRoot_pft(NZ)
  CPRTW=FracRootElmAllocm(ielmc,k_woody_comp)*rPCStalk_pft(NZ)+FracRootElmAllocm(ielmc,k_fine_comp)*rPCRootr_pft(NZ)

  FracPetolShethAlloc2Litr(ielmc,1:NumOfPlantLitrCmplxs) = FracLeafShethElmAlloc2Litr(ielmc,1:NumOfPlantLitrCmplxs)

  FracLeafShethElmAlloc2Litr(ielmn,k_woody_comp) = AMIN1(FracLeafShethElmAlloc2Litr(ielmc,k_woody_comp)*rNCStalk_pft(NZ)/CNLFW,1._r8)
  FracLeafShethElmAlloc2Litr(ielmp,k_woody_comp) = AMIN1(FracLeafShethElmAlloc2Litr(ielmc,k_woody_comp)*rPCStalk_pft(NZ)/CPLFW,1._r8)
  
  FracPetolShethAlloc2Litr(ielmn,k_woody_comp) = AMIN1(FracLeafShethElmAlloc2Litr(ielmc,k_woody_comp)*rNCStalk_pft(NZ)/CNSHW,1._r8)
  FracPetolShethAlloc2Litr(ielmp,k_woody_comp) = AMIN1(FracLeafShethElmAlloc2Litr(ielmc,k_woody_comp)*rPCStalk_pft(NZ)/CPSHW,1._r8)
  
  FracWoodStalkElmAlloc2Litr(ielmn,k_woody_comp) = AMIN1(FracWoodStalkElmAlloc2Litr(ielmc,k_woody_comp)*rNCStalk_pft(NZ)/CNRTW,1._r8)
  FracWoodStalkElmAlloc2Litr(ielmp,k_woody_comp) = AMIN1(FracWoodStalkElmAlloc2Litr(ielmc,k_woody_comp)*rPCStalk_pft(NZ)/CPRTW,1._r8)

  FracRootElmAllocm(ielmn,k_woody_comp) = AMIN1(FracRootElmAllocm(ielmc,k_woody_comp)*rNCRoot_pft(NZ)/CNRTW,1._r8)
  FracRootElmAllocm(ielmp,k_woody_comp) = AMIN1(FracRootElmAllocm(ielmc,k_woody_comp)*rPCRootr_pft(NZ)/CPRTW,1._r8)

  DO NE=2,NumPlantChemElms
    FracPetolShethAlloc2Litr(NE,k_fine_comp)   = AZMAX1(1.0_r8-FracPetolShethAlloc2Litr(NE,k_woody_comp))
    FracWoodStalkElmAlloc2Litr(NE,k_fine_comp) = AZMAX1(1.0_r8-FracWoodStalkElmAlloc2Litr(NE,k_woody_comp))
    FracLeafShethElmAlloc2Litr(NE,k_fine_comp) = AZMAX1(1.0_r8-FracLeafShethElmAlloc2Litr(NE,k_woody_comp))
    FracRootElmAllocm(NE,k_fine_comp)          = AZMAX1(1.0_r8-FracRootElmAllocm(NE,k_woody_comp))
  ENDDO
  !
  !     SHOOT AND ROOT TEMPERATURE FUNCTIONS FOR MAINTENANCE
  !     RESPIRATION FROM TEMPERATURES WITH OFFSETS FOR THERMAL ADAPTATION
  !
  !     TKC,TKCM=canopy temperature,canopy temp used in Arrhenius eqn
  !     TKS_vr,TKSM=soil temperature,soil temp used in Arrhenius eqn
  !     TempOffset_pft=shift in Arrhenius curve for thermal adaptation
  !     TFN5,TFN6_vr=temperature function for canopy,root mntc respn (25 oC =1)
  !     8.3143,710.0=gas constant,enthalpy
  !     62500,195000,232500=energy of activn,high,low temp inactivn(KJ mol-1)
  !
  TKCM = real_truncate(TKC_pft(NZ)+TempOffset_pft(NZ),1.e-3_r8)
  TFN5 = calc_plant_maint_tempf(TKCM)
  
  D7: DO L=NU,MaxNumRootLays
    TKSM       = real_truncate(TKS_vr(L)+TempOffset_pft(NZ),1.e-3_r8)
    TFN6_vr(L) = calc_plant_maint_tempf(TKSM)
  ENDDO D7
  !
  !     PRIMARY ROOT NUMBER
  !
  !     RootBiomCPerPlant_pft=root mass per plant used to calculate primary root number
  !     WTRT,PP=root mass,PFT population
  !here it tries to enforce the Shinozaki's pipe model (1964), to some extent (Jinyun Tang)
  RootBiomCPerPlant_pft(NZ)  = AMAX1(dscal*RootBiomCPerPlant_pft(NZ),RootElms_pft(ielmc,NZ)/PlantPopuLive_pft(NZ))
  
  IF(NumStructuralRootAxes_pft(NZ).GT.0)THEN
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
      Num1stAxesPerStructRootX_pft(NZ) = AMAX1(1.0_r8,RootBiomCPerPlant_pft(NZ)**0.833_r8)/NumStructuralRootAxes_pft(NZ)
    elseif(is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))then
      if(lcoarseroot)then
        !for woody vascular, structural root axes is counted by one
        Num1stAxesPerStructRootX_pft(NZ) = AMIN1(1._R8,AMAX1(1.0_r8,RootBiomCPerPlant_pft(NZ)**0.7143_r8)/NumStructuralRootAxes_pft(NZ))
      else
        Num1stAxesPerStructRootX_pft(NZ) = AMAX1(1.0_r8,RootBiomCPerPlant_pft(NZ)**0.7143_r8)/NumStructuralRootAxes_pft(NZ)
      endif
    else
      !for herbaceous plant, **(1./0.95) 
      Num1stAxesPerStructRootX_pft(NZ) = AMAX1(1.0_r8,RootBiomCPerPlant_pft(NZ)**1.053_r8)/NumStructuralRootAxes_pft(NZ)
    endif
  else
    Num1stAxesPerStructRootX_pft(NZ) = 1._r8
  ENDIF
  Num1stAxesPerStructRootXPOP_pft(NZ) = Num1stAxesPerStructRootX_pft(NZ)*PlantPopuLive_pft(NZ)
  Num1stRootAxesPP_pft(NZ)            = Num1stAxesPerStructRootX_pft(NZ)*NumStructuralRootAxes_pft(NZ)

  !
  !     WATER STRESS FUNCTIONS FOR EXPANSION AND GROWTH RESPIRATION
  !     FROM CANOPY TURGOR
  !
  !     PSICanopyTurg_pft,TurgPSIMin4OrganExtens=current,minimum canopy turgor potl for expansion,extension
  !     Stomata_Stress=stomatal resistance function of canopy turgor
  !     PSICanopy_pft=canopy water potential
  !     WaterStress4Groth=growth function of canopy water potential
  
  !
  TurgEff4LeafPetolExpansion=real_truncate(AMIN1(1.0_r8,AZMAX1(PSICanopyTurg_pft(NZ)-TurgPSIMin4OrganExtens)),1.e-4_r8)

  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ)))THEN
    !bryophyte, no turgor
    Stomata_Stress    = 1.0_r8
    WaterStress4Groth = EXP(0.05_r8*AMAX1(PSICanopy_pft(NZ),-5000._r8))
  ELSE
    !others
    Stomata_Stress    = EXP(-PSICanopyTurg_pft(NZ)/RCS_pft(NZ))
    WaterStress4Groth = EXP(0.10_r8*AMAX1(PSICanopy_pft(NZ),-5000._r8))
  ENDIF
  WaterStress4Groth               = real_truncate(WaterStress4Groth,1.e-4_r8)
  TurgEff4CanopyResp              = fRespWatSens(TurgEff4LeafPetolExpansion,iPlantRootProfile_pft(NZ))
  plt_biom%StomatalStress_pft(NZ) = Stomata_Stress
  TurgEff4CanopyResp_pft(NZ)      = TurgEff4CanopyResp
  call PrintInfo('end '//subname)
  end associate
  end subroutine StagePlantForGrowth

!----------------------------------------------------------------------------------------------------
  subroutine AccumulateStates(yearIJ,NZ)
  
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='AccumulateStates'
  integer :: L,NR,N,NE,NB

!     begin_execution
  associate(                                                          &
    ZEROS                     => plt_site%ZEROS                      ,& !input  :threshold zero for numerical stability,[-]  
    CanopyNodulNonstElms_brch => plt_biom%CanopyNodulNonstElms_brch  ,& !input  :branch nodule nonstructural element, [g d-2]
    CanopyNodulStrutElms_brch => plt_biom%CanopyNodulStrutElms_brch  ,& !input  :branch nodule structural element, [g d-2]
    CanopyStalkSurfArea_lbrch => plt_morph%CanopyStalkSurfArea_lbrch ,& !input  :plant canopy layer branch stem area, [m2 d-2]
    CanopyLeafSheathC_brch    => plt_biom%CanopyLeafSheathC_brch     ,& !input  :plant branch leaf + sheath C, [g d-2]
    MaxSoilLays4Root_pft      => plt_morph%MaxSoilLays4Root_pft      ,& !input  :maximum soil layer number for all root axes,[-]
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft         ,& !input  :number of branches,[-]
    RootH2PO4Uptake_pft       => plt_rbgc%RootH2PO4Uptake_pft        ,& !input  :total root uptake of PO4, [g d-2 h-1]
    RootHPO4Uptake_pft        => plt_rbgc%RootHPO4Uptake_pft         ,& !input  :total root uptake of HPO4, [g d-2 h-1]
    Soil2RootMycoExudE_pft    => plt_rbgc%Soil2RootMycoExudE_pft     ,& !input  :total root uptake (+ve) - exudation (-ve) of dissolved element, [g d-2 h-1]
    RootN2Fix_pft             => plt_rbgc%RootN2Fix_pft              ,& !input  :total root N2 fixation, [g d-2 h-1]
    RootNH4Uptake_pft         => plt_rbgc%RootNH4Uptake_pft          ,& !input  :total root uptake of NH4, [g d-2 h-1]
    RootNO3Uptake_pft         => plt_rbgc%RootNO3Uptake_pft          ,& !input  :total root uptake of NO3, [g d-2 h-1]
    RootNodulNonstElms_rpvr   => plt_biom%RootNodulNonstElms_rpvr    ,& !input  :root layer nonstructural element, [g d-2]
    RootNodulStrutElms_rpvr   => plt_biom%RootNodulStrutElms_rpvr    ,& !input  :root layer nodule element, [g d-2]
    CanopyN2Fix_pft           => plt_rbgc%CanopyN2Fix_pft            ,& !input :total canopy N2 fixation, [g d-2 h-1]                    
    SetNumberSeeds_brch       => plt_morph%SetNumberSeeds_brch       ,& !input  :branch grain number, [d-2]
    SapwoodBiomassC_brch      => plt_biom%SapwoodBiomassC_brch       ,& !input  :branch live stalk C, [gC d-2]
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft        ,& !input  :N2 fixation type,[-]
    CanopyNodulNonstElms_pft  => plt_biom%CanopyNodulNonstElms_pft   ,& !inoput :canopy nodule nonstructural element, [g d-2]
    CanopyStemSurfAreaZ_pft   => plt_morph%CanopyStemSurfAreaZ_pft   ,& !output :plant canopy layer stem area, [m2 d-2]
    PlantExudElm_CumYr_pft    => plt_rbgc%PlantExudElm_CumYr_pft     ,& !inoput :total net root element uptake (+ve) - exudation (-ve), [gC d-2 ]
    PlantN2Fix_CumYr_pft      => plt_bgcr%PlantN2Fix_CumYr_pft       ,& !inoput :total plant N2 fixation, [g d-2 ]
    PlantRootSoilElmNetX_pft  => plt_rbgc%PlantRootSoilElmNetX_pft   ,& !inoput :net root element uptake (+ve) - exudation (-ve), [gC d-2 h-1]
    RootUptk_N_CumYr_pft      => plt_rbgc%RootUptk_N_CumYr_pft       ,& !inoput :cumulative plant N uptake, [gN d-2]
    RootUptk_P_CumYr_pft      => plt_rbgc%RootUptk_P_CumYr_pft       ,& !inoput :cumulative plant P uptake, [gP d-2]
    VcMaxRubiscoRef_brch      => plt_photo%VcMaxRubiscoRef_brch      ,& !input  :maximum rubisco carboxylation rate at reference temperature, [umol g-1 h-1]    
    VoMaxRubiscoRef_brch      => plt_photo%VoMaxRubiscoRef_brch      ,& !input  :maximum rubisco oxygenation rate at reference temperature, [umol g-1 h-1]    
    VcMaxPEPCarboxyRef_brch   => plt_photo%VcMaxPEPCarboxyRef_brch   ,& !input  :reference maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]        
    fNCLFW_brch               => plt_pheno%fNCLFW_brch               ,& !input : NC ratio of growing leaf on branch, [gN/gC]
    fPCLFW_brch               => plt_pheno%fPCLFW_brch               ,& !input : PC ratio of growing leaf on branch, [gP/gC]
    LeafAreaLive_brch         => plt_morph%LeafAreaLive_brch         ,& !output :branch leaf area, [m2 d-2]
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch         ,& !output :branch leaf structural element mass, [g d-2]
    PetolShethStrutElms_brch  => plt_biom%PetolShethStrutElms_brch   ,& !output :branch sheath structural element, [g d-2]
    CanopyLeafArea_pft        => plt_morph%CanopyLeafArea_pft        ,& !output :plant canopy leaf area, [m2 d-2]
    CanopyLeafSheathC_pft     => plt_biom%CanopyLeafSheathC_pft      ,& !output :canopy leaf + sheath C, [g d-2]
    CanopySeedNum_pft         => plt_morph%CanopySeedNum_pft         ,& !output :canopy grain number, [d-2]
    CanopySeedNumX_pft        => plt_morph%CanopySeedNumX_pft        ,& !output :last nonzero canopy grain number, [d-2]    
    CanopySapwoodC_pft        => plt_biom%CanopySapwoodC_pft         ,& !output :canopy active stalk C, [g d-2]
    CanopyStemSurfArea_pft    => plt_morph%CanopyStemSurfArea_pft    ,& !output :plant stem area, [m2 d-2]
    EarStrutElms_pft          => plt_biom%EarStrutElms_pft           ,& !output :canopy ear structural element, [g d-2]
    GrainStrutElms_pft        => plt_biom%GrainStrutElms_pft         ,& !output :canopy grain structural element, [g d-2]
    HuskStrutElms_pft         => plt_biom%HuskStrutElms_pft          ,& !output :canopy husk structural element mass, [g d-2]
    LeafStrutElms_pft         => plt_biom%LeafStrutElms_pft          ,& !output :canopy leaf structural element mass, [g d-2]
    PetolShethStrutElms_pft   => plt_biom%PetolShethStrutElms_pft    ,& !output :canopy sheath structural element mass, [g d-2]
    StalkRsrvElms_pft         => plt_biom%StalkRsrvElms_pft          ,& !output :canopy reserve element mass, [g d-2]
    StalkStrutElms_pft        => plt_biom%StalkStrutElms_pft         ,& !output :canopy stalk structural element mass, [g d-2]
    fNCLFW_pft                => plt_pheno%fNCLFW_pft                ,& !output : NC ratio of growing leaf, [gN/gC]
    fPCLFW_pft                => plt_pheno%fPCLFW_pft                 & !output : PC ratio of growing leaf, [gP/gC]

  )
  !
  !     ACCUMULATE PFT STATE VARIABLES FROM BRANCH STATE VARIABLES
  !
  !
  call PrintInfo('beg '//subname)
  DO NB=1,NumOfBranches_pft(NZ)    
    if(isclose(CanopyLeafSheathC_brch(NB,NZ),0._r8))then
      LeafStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
      PetolShethStrutElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
      LeafAreaLive_brch(NB,NZ)                     = 0._r8
    endif
 ENDDO  

  CanopySapwoodC_pft(NZ)   = 0._r8
  CanopyLeafSheathC_pft(NZ) = 0._r8
  CanopySeedNum_pft(NZ)    = 0._r8
  CanopyLeafArea_pft(NZ)   = 0._r8
  CanopyStemSurfArea_pft(NZ)   = 0._r8
  ! Rebuild layer areas here so the result does not depend on an earlier growth reset.
  CanopyStemSurfAreaZ_pft(:,NZ) = 0._r8

  DO NB=1,NumOfBranches_pft(NZ)        
    CanopySapwoodC_pft(NZ)     = CanopySapwoodC_pft(NZ)+SapwoodBiomassC_brch(NB,NZ)
    CanopyLeafSheathC_pft(NZ)   = CanopyLeafSheathC_pft(NZ) +CanopyLeafSheathC_brch(NB,NZ)
    CanopySeedNum_pft(NZ)      = CanopySeedNum_pft(NZ)+SetNumberSeeds_brch(NB,NZ)
    CanopyLeafArea_pft(NZ)     = CanopyLeafArea_pft(NZ)+LeafAreaLive_brch(NB,NZ)
    fNCLFW_pft(NZ)             = fNCLFW_pft(NZ)+fNCLFW_brch(NB,NZ)*LeafAreaLive_brch(NB,NZ)
    fPCLFW_pft(NZ)             = fPCLFW_pft(NZ)+fPCLFW_brch(NB,NZ)*LeafAreaLive_brch(NB,NZ)
    DO L=1,NumCanopyLayers1
      CanopyStemSurfAreaZ_pft(L,NZ)=CanopyStemSurfAreaZ_pft(L,NZ)+CanopyStalkSurfArea_lbrch(L,NB,NZ)
    ENDDO
  ENDDO
  
  ! The PFT total and layer profile must describe the same surviving branch geometry.
  CanopyStemSurfArea_pft(NZ) = SUM(CanopyStemSurfAreaZ_pft(1:NumCanopyLayers1,NZ))

  if(CanopySeedNum_pft(NZ)>0._r8)CanopySeedNumX_pft(NZ)=CanopySeedNum_pft(NZ)
  if(CanopyLeafArea_pft(NZ).GT.ZEROs)then
    fNCLFW_pft(NZ)=fNCLFW_pft(NZ)/CanopyLeafArea_pft(NZ)
    fPCLFW_pft(NZ)=fPCLFW_pft(NZ)/CanopyLeafArea_pft(NZ)
  endif
  !
  !     ACCUMULATE NODULE STATE VATIABLES FROM NODULE LAYER VARIABLES
  !
  !     iPlantNfixType_pft=N2 fixation: 1,2,3=rapid to slow root symbiosis
  !     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
  !     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
  !
  !
  !     ACCUMULATE TOTAL SOIL-PLANT C,N,P EXCHANGE
  !
  !     PlantN2Fix_CumYr_pft=cumulative PFT N2 fixation
  !     PlantRootSoilElmNetX_pft= > 0 taking from soil 
  !
  PlantRootSoilElmNetX_pft(1:NumPlantChemElms,NZ)=Soil2RootMycoExudE_pft(1:NumPlantChemElms,NZ)

  PlantRootSoilElmNetX_pft(ielmn,NZ)=PlantRootSoilElmNetX_pft(ielmn,NZ)+&
    RootNH4Uptake_pft(NZ)+RootNO3Uptake_pft(NZ)+RootN2Fix_pft(NZ)
  PlantRootSoilElmNetX_pft(ielmp,NZ)=PlantRootSoilElmNetX_pft(ielmp,NZ)+&
    RootH2PO4Uptake_pft(NZ)+RootHPO4Uptake_pft(NZ)

  PlantExudElm_CumYr_pft(1:NumPlantChemElms,NZ)=PlantExudElm_CumYr_pft(1:NumPlantChemElms,NZ)+&
    Soil2RootMycoExudE_pft(1:NumPlantChemElms,NZ)

  RootUptk_N_CumYr_pft(NZ)=RootUptk_N_CumYr_pft(NZ)+RootNH4Uptake_pft(NZ)+RootNO3Uptake_pft(NZ) + &
    Soil2RootMycoExudE_pft(ielmn,NZ)
  RootUptk_P_CumYr_pft(NZ)=RootUptk_P_CumYr_pft(NZ)+RootH2PO4Uptake_pft(NZ)+RootHPO4Uptake_pft(NZ) + &
    Soil2RootMycoExudE_pft(ielmp,NZ)

  PlantN2Fix_CumYr_pft(NZ)=PlantN2Fix_CumYr_pft(NZ)+RootN2Fix_pft(NZ)+CanopyN2Fix_pft(NZ)
  
  call PrintInfo('end '//subname)
  end associate
  end subroutine AccumulateStates

!----------------------------------------------------------------------------------------------------
  subroutine ZeroPlantStates(I,J,NZ)  
  !
  !for annual plants
  implicit none
  integer, intent(in) :: I,J,NZ

  integer :: NE
  DO NE=1,NumPlantChemElms
    plt_biom%CanopyNonstElms_pft(NE,NZ)=0._r8
    plt_biom%StalkRsrvElms_pft(NE,NZ)=0._r8
    plt_biom%ShootElms_pft(NE,NZ)=0._r8
    plt_biom%PetolShethStrutElms_pft(NE,NZ)=0._r8
    plt_biom%StalkStrutElms_pft(NE,NZ)=0._r8
    plt_biom%LeafStrutElms_pft(NE,NZ)=0._r8
    plt_biom%HuskStrutElms_pft(NE,NZ)=0._r8
    plt_biom%GrainStrutElms_pft(NE,NZ)=0._r8
    plt_biom%EarStrutElms_pft(NE,NZ)=0._r8
    plt_biom%ShootNoduleElms_pft(NE,NZ)=0._r8
    plt_biom%RootNoduleElms_pft(NE,NZ)=0._r8
    plt_rbgc%PlantRootSoilElmNetX_pft(NE,NZ)=0._r8
  ENDDO

  plt_biom%CanopySapwoodC_pft(NZ)=0._r8
  plt_biom%CanopyLeafSheathC_pft(NZ) =0._r8
  plt_morph%CanopySeedNum_pft(NZ) =0._r8
  plt_morph%CanopyLeafArea_pft(NZ)=0._r8
  plt_morph%CanopyStemSurfArea_pft(NZ)=0._r8
  plt_morph%CanopyStemSurfAreaZ_pft(1:NumCanopyLayers1,NZ)=0._r8

  end subroutine ZeroPlantStates

  ![tail]
end module grosubsMod
