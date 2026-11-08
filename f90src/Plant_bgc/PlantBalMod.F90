module PlantBalMod
!
!Description
!code to do mass balance calculation for plant bgc

  use data_kind_mod,    only: r8 => DAT_KIND_R8,yearIJ_type
  use EcoSiMParDataMod, only: pltpar
  use abortutils,       only: endrun
  use EcoSIMCtrlMod, only: guard_floor_dying_stand, VEGC_BAL_FLOOR_dying
  use PlantMathFuncMod
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantPhotosynthesisAPIData, only : plt_photo
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantAllometryAPIData, only : plt_allom
  use PlantBiomassAPIData, only : plt_biom
  use PlantDisturbanceAPIData, only : plt_distb
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantRootBGCAPIData, only : plt_rbgc
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public :: SumPlantBiome
  public :: SumPlantBiomStates
  public :: SumRootBiome
  public :: SumPlantBranchBiome
  public :: EnterPlantBalance
  public :: ExitPlantBalance
  public :: SumPlantRootGas
  public :: SumRootAR
  public :: CheckPlantBalanceZ
  public :: SumLitfallBlg
  public :: SumCanopyBiome
  public :: SumLitfallAbg
  logical,save  :: lfile(2)=.true.
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine SumPlantBiome(yearIJ,NZ,header,tvegE,canopyE,RootE)

  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  character(len=*), intent(in) :: header
  real(r8),optional, intent(out) :: tvegE(NumPlantChemElms)
  real(r8),optional, intent(out) :: canopyE(NumPlantChemElms)
  real(r8),optional, intent(out) :: RootE(NumPlantChemElms)

  real(r8) :: tvegE1(NumPlantChemElms),canopyE1(NumPlantChemElms),RootE1(NumPlantChemElms)
  integer :: L,K,N,NE,NB,M
  real(r8) :: root1st,root2nd

!     begin_execution
  associate(                                                         &
    NU                       => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    NH3Dep2Can_brch          => plt_rbgc%NH3Dep2Can_brch            ,& !input  :gaseous NH3 flux fron root disturbance band, [g d-2 h-1]
    CanopyGrosRCO2_pft       => plt_bgcr%CanopyGrosRCO2_pft         ,& !input  :canopy plant+nodule autotrophic respiraiton, [gC d-2]
    NumOfBranches_pft        => plt_morph%NumOfBranches_pft         ,& !input  :number of branches,[-]
    LitrfallElms_pvr         => plt_bgcr%LitrfallElms_pvr           ,& !input  :plant LitrFall element, [g d-2 h-1]
    Myco_pft                 => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    NMaxRootBotLayer_pft     => plt_morph%NMaxRootBotLayer_pft      ,& !input  :maximum soil layer number for all root axes,[-]
    RootElms_pft             => plt_biom%RootElms_pft               ,& !input  :plant root element mass, [g d-2]
    RootCO2Autor_pvr         => plt_rbgc%RootCO2Autor_pvr           ,& !input  :root respiration constrained by O2, [g d-2 h-1]
    ShootElms_pft            => plt_biom%ShootElms_pft              ,& !input  :current time whole plant shoot element mass, [g d-2]
    NH3Dep2Can_pft           => plt_bgcr%NH3Dep2Can_pft             ,& !inoput :canopy NH3 flux, [gN d-2 h-1]
    RootAutoCO2_pft          => plt_bgcr%RootAutoCO2_pft            ,& !inoput :root autotrophic respiraiton, [gC d-2]
    LitrfallElms_pft         => plt_bgcr%LitrfallElms_pft           ,& !output :plant element LitrFall, [g d-2 h-1]
    LitrfallAbvgElms_pft     => plt_bgcr%LitrfallAbvgElms_pft       ,& !output :aboveground plant element LitrFall, [g d-2 h-1]
    LitrfallBlgrElms_pft     => plt_bgcr%LitrfallBlgrElms_pft       ,& !output :belowground plant element LitrFall, [g d-2 h-1]
    GrossResp_pft            => plt_bgcr%GrossResp_pft               & !output :total plant respiration, [gC d-2 ]
  )

  call SumPlantBiomStates(yearIJ,NZ,header,tvegE=tvegE1,canopyE=canopyE1,RootE=RootE1)

  if(present(tvegE))tvegE=tvegE1
  if(present(canopyE))canopyE=canopyE1
  if(present(RootE))RootE=RootE1

  call SumRootAR(NZ)
  
  call SumLitfallAbg(NZ)

  call SumLitfallBlg(NZ)

  !sum fluxes
  !     NH3Dep2Can_brch,NH3Dep2Can_pft=PFT NH3 flux between atmosphere and branch,canopy
  NH3Dep2Can_pft(NZ)=0._r8
  DO NB=1,NumOfBranches_pft(NZ)
    NH3Dep2Can_pft(NZ)=NH3Dep2Can_pft(NZ)+NH3Dep2Can_brch(NB,NZ)
  ENDDO    

  DO NE=1,NumPlantChemElms
    LitrfallElms_pft(NE,NZ)=LitrfallAbvgElms_pft(NE,NZ)+LitrfallBlgrElms_pft(NE,NZ)
  ENDDO  

  GrossResp_pft(NZ)=RootAutoCO2_pft(NZ)+CanopyGrosRCO2_pft(NZ)

  end associate
  end subroutine SumPlantBiome
!----------------------------------------------------------------------------------------------------
  subroutine SumLitfallAbg(NZ)  
  implicit none
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='SumLitfallAbg'
  integer :: L,K,M,NE

  associate(                                                         &
    LitrfallElms_pvr         => plt_bgcr%LitrfallElms_pvr           ,& !input  :plant LitrFall element, [g d-2 h-1]
    LitrfallAbvgElms_pft     => plt_bgcr%LitrfallAbvgElms_pft        & !output :aboveground plant element LitrFall, [g d-2 h-1]
  )    
  call PrintInfo('beg '//subname)
  LitrfallAbvgElms_pft(1:NumPlantChemElms,NZ) = 0._r8
  L  = 0
  DO K=1,pltpar%NumOfPlantLitrCmplxs      
    DO M=1,jsken
      DO NE=1,NumPlantChemElms
        LitrfallAbvgElms_pft(NE,NZ)=LitrfallAbvgElms_pft(NE,NZ)+LitrfallElms_pvr(NE,M,K,L,NZ)
      ENDDO
    ENDDO
  ENDDO      
  call PrintInfo('end '//subname)
  end associate
  end subroutine SumLitfallAbg
!----------------------------------------------------------------------------------------------------
  subroutine SumLitfallBlg(NZ)
  implicit none
  integer, intent(in) :: NZ
  integer :: L,K,M,NE

  associate(                                                         &
    NU                       => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    MaxNumRootLays           => plt_site%MaxNumRootLays             ,& !input  :maximum root layer number,[-]    
    LitrfallElms_pvr         => plt_bgcr%LitrfallElms_pvr           ,& !input  :plant LitrFall element, [g d-2 h-1]
    LitrfallBlgrElms_pft     => plt_bgcr%LitrfallBlgrElms_pft        & !output :belowground plant element LitrFall, [g d-2 h-1]
  )
  LitrfallBlgrElms_pft(1:NumPlantChemElms,NZ) = 0._r8

  DO L=1,MaxNumRootLays
    DO K=1,pltpar%NumOfPlantLitrCmplxs      
      DO M=1,jsken
        DO NE=1,NumPlantChemElms
          LitrfallBlgrElms_pft(NE,NZ)=LitrfallBlgrElms_pft(NE,NZ)+LitrfallElms_pvr(NE,M,K,L,NZ)
        ENDDO
      ENDDO
    ENDDO      
  ENDDO

  end associate
  end subroutine SumLitfallBlg    
!----------------------------------------------------------------------------------------------------

  subroutine SumRootAR(NZ)
  implicit none
  integer, intent(in) :: NZ
  integer :: N,L
  associate(                                                         &  
    NU                       => plt_site%NU                         ,& !input  :current soil surface layer number, [-]  
    Myco_pft                 => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-] 
    RootCO2Autor_pvr         => plt_rbgc%RootCO2Autor_pvr           ,& !input  :root respiration constrained by O2, [g d-2 h-1]  
    RootAutoCO2_pft          => plt_bgcr%RootAutoCO2_pft            ,& !inoput :root autotrophic respiraiton, [gC d-2]
    NMaxRootBotLayer_pft     => plt_morph%NMaxRootBotLayer_pft       & !input  :maximum soil layer number for all root axes,[-]    
  )
  RootAutoCO2_pft(NZ)=0._r8
  DO N=1,Myco_pft(NZ)
    DO L=1,NMaxRootBotLayer_pft(NZ)
      RootAutoCO2_pft(NZ)=RootAutoCO2_pft(NZ)+RootCO2Autor_pvr(N,L,NZ)
    ENDDO     
  ENDDO  
  end associate
  end subroutine SumRootAR
!----------------------------------------------------------------------------------------------------
  subroutine SumPlantRootGas(I,J,NZ1)

  implicit none
  integer, intent(in) :: I,J
  INTEGER, optional, INTENT(IN) :: NZ1
  integer :: NZ ,L,N,idg,K
  real(r8) :: trcg(idg_beg:idg_NH3)
  associate(                                                 &
    trcg_rootml_pvr      => plt_rbgc%trcg_rootml_pvr        ,& !input  :root gas content, [g d-2]
    trcs_rootml_pvr      => plt_rbgc%trcs_rootml_pvr        ,& !input  :root aqueous content, [g d-2]
    Myco_pft             => plt_morph%Myco_pft              ,& !input  :mycorrhizal type (no or yes),[-]
    NK                   => plt_site%NK                     ,& !input  : total layers with roots,  can be set to NL_col
    MaxSoilLays4Root_pft => plt_morph%MaxSoilLays4Root_pft  ,& !input  :maximum soil layer number for all root axes,[-]
    trcg_root_vr         => plt_rbgc%trcg_root_vr            & !inoput :total root internal gas flux, [g d-2 h-1]
  )
  trcg_root_vr(idg_beg:idg_NH3,:)   = 0._r8

  trcg(:)=0._r8
  IF(present(NZ1))THEN
    NZ=NZ1
    DO L=1,MIN(MaxSoilLays4Root_pft(NZ)+1,NK)
      DO N=1,Myco_pft(NZ)  
        DO idg=idg_beg,idg_NH3
          trcg_root_vr(idg,L)=trcg_root_vr(idg,L)+trcs_rootml_pvr(idg,N,L,NZ)+trcg_rootml_pvr(idg,N,L,NZ)        
        ENDDO        
      ENDDO  
      DO idg=idg_beg,idg_NH3
        trcg(idg)=trcg(idg)+trcg_root_vr(idg,L)
      ENDDO
    ENDDO
  ELSE
    DO NZ=1,plt_site%NP  

      DO L=1,MIN(MaxSoilLays4Root_pft(NZ)+1,NK)
        DO N=1,Myco_pft(NZ)  
          DO idg=idg_beg,idg_NH3
            trcg_root_vr(idg,L)=trcg_root_vr(idg,L)+trcs_rootml_pvr(idg,N,L,NZ)+trcg_rootml_pvr(idg,N,L,NZ)        
          ENDDO        
        ENDDO  
        DO idg=idg_beg,idg_NH3
          trcg(idg)=trcg(idg)+trcg_root_vr(idg,L)
        ENDDO
      ENDDO
    ENDDO
  ENDIF  

  end associate
  end subroutine SumPlantRootGas
!----------------------------------------------------------------------------------------------------
  subroutine SumCanopyBiome(yearIJ,NZ,canopyE)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), optional, intent(out) :: canopyE(NumPlantChemElms)
  character(len=*), parameter :: subname='SumCanopyBiome'
  integer :: NB,NE,jk

  associate(                                                          &
    k_fine_comp               => pltpar%k_fine_comp                  ,& !input  :fine litter complex id
    k_woody_comp              => pltpar%k_woody_comp                 ,& !input  :woody litter complex id
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft        ,& !input  :N2 fixation type,[-]
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft         ,& !input  :number of branches,[-]
    RootElmsBeg_pft           => plt_biom%RootElmsBeg_pft            ,& !input  :plant root element at previous time step, [g d-2]    
    ShootElmsBeg_pft          => plt_biom%ShootElmsBeg_pft           ,& !input  :previous whole plant shoot element mass, [g d-2]    
    CanopyNodulNonstElms_brch => plt_biom%CanopyNodulNonstElms_brch  ,& !input  :branch nodule nonstructural element, [g d-2]
    CanopyNodulStrutElms_brch => plt_biom%CanopyNodulStrutElms_brch  ,& !input  :branch nodule structural element, [g d-2]
    PetolShethStrutElms_brch  => plt_biom%PetolShethStrutElms_brch   ,& !input :branch sheath structural element, [g d-2]    
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch         ,& !input :branch leaf structural element mass, [g d-2]    
    StalkRsrvElms_brch        => plt_biom%StalkRsrvElms_brch         ,& !input  :branch reserve element mass, [g d-2]    
    StalkStrutElms_brch       => plt_biom%StalkStrutElms_brch        ,& !input  :branch stalk structural element mass, [g d-2]    
    HuskStrutElms_brch        => plt_biom%HuskStrutElms_brch         ,& !input  :branch husk structural element mass, [g d-2]    
    EarStrutElms_brch         => plt_biom%EarStrutElms_brch          ,& !input  :branch ear structural chemical element mass, [g d-2]    
    GrainStrutElms_brch       => plt_biom%GrainStrutElms_brch        ,& !input  :branch grain structural element mass, [g d-2]    
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch       ,& !input  :branch nonstructural element, [g d-2]    
    StandDeadCompKElms_pft    => plt_biom%StandDeadCompKElms_pft     ,& !input  :standing dead element fraction, [g d-2]
    ShootElms_brch            => plt_biom%ShootElms_brch             ,& !input  :branch shoot structural element mass, [g d-2]
    ShootElms_pft             => plt_biom%ShootElms_pft              ,& !inoput :current time whole plant shoot element mass, [g d-2]
    CanopyHeightLive_pft      => plt_morph%CanopyHeightLive_pft      ,& !input  :live plant canopy height, [m]
    CanopyHeightDead_pft      => plt_morph%CanopyHeightDead_pft      ,& !input  :canopy height for standing dead, [m]
    LeafStrutElms_pft         => plt_biom%LeafStrutElms_pft          ,& !output :canopy leaf structural element mass, [g d-2]    
    ShootFineNonLeafElms_pft  => plt_biom%ShootFineNonLeafElms_pft   ,& !output :canopy non-leaf fine structural element mass, [g d-2]
    FracLeafShethElmAlloc2Litr=> plt_allom%FracLeafShethElmAlloc2Litr,& !input  :leave-sheath element allocation, [-]    
    ShootNonstElms_pft        => plt_biom%ShootNonstElms_pft         ,& !output :shoot nonstructural element mass, [g d-2]
    ShootLeafElms_pft         => plt_biom%ShootLeafElms_pft          ,& !output :shoot leaf element mass, [g d-2]
    ShootWoodyElms_pft        => plt_biom%ShootWoodyElms_pft         ,& !output :canopy woody element mass, [g d-2]
    FracPetolShethAlloc2Litr  => plt_allom%FracPetolShethAlloc2Litr  ,& !input  :woody element allocation, [-]    
    PetolShethStrutElms_pft   => plt_biom%PetolShethStrutElms_pft    ,& !output :canopy sheath structural element mass, [g d-2]    
    ShootNoduleElms_pft       => plt_biom%ShootNoduleElms_pft        ,& !output :current time canopy nodule element mass, [g d-2]    
    StalkRsrvElms_pft         => plt_biom%StalkRsrvElms_pft          ,& !output :canopy reserve element mass, [g d-2]
    StalkStrutElms_pft        => plt_biom%StalkStrutElms_pft         ,& !output :canopy stalk structural element mass, [g d-2]
    EarStrutElms_pft          => plt_biom%EarStrutElms_pft           ,& !output :canopy ear structural element, [g d-2]
    GrainStrutElms_pft        => plt_biom%GrainStrutElms_pft         ,& !output :canopy grain structural element, [g d-2]
    HuskStrutElms_pft         => plt_biom%HuskStrutElms_pft          ,& !output :canopy husk structural element mass, [g d-2]
    StandDeadStrutElms_pft    => plt_biom%StandDeadStrutElms_pft     ,& !output :standing dead element, [g d-2]
    CanopyNonstElms_pft       => plt_biom%CanopyNonstElms_pft         & !output :canopy nonstructural element concentration, [g d-2]        
  )
  call PrintInfo('beg '//subname)
  !shoots
  DO NB=1,NumOfBranches_pft(NZ)
    CALL SumPlantBranchBiome(NB,NZ)
  ENDDO

  !add nodule, currently, a plant either has canopy or root N-fixing symbiosis, not both
  ShootNoduleElms_pft(:,NZ)=0._r8
  IF(is_plant_N2fix(iPlantNfixType_pft(NZ)))THEN
    IF(is_canopy_N2fix(iPlantNfixType_pft(NZ)))THEN
      DO NE=1,NumPlantChemElms
        ShootNoduleElms_pft(NE,NZ)=ShootNoduleElms_pft(NE,NZ)  &
          +sum(CanopyNodulStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)) &
          +sum(CanopyNodulNonstElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
      ENDDO
    ENDIF
  ENDIF

  DO NE=1,NumPlantChemElms
    ShootElms_pft(NE,NZ)           = AZMAX1(sum(ShootElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    LeafStrutElms_pft(NE,NZ)       = AZMAX1(sum(LeafStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    PetolShethStrutElms_pft(NE,NZ) = AZMAX1(sum(PetolShethStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    StalkStrutElms_pft(NE,NZ)      = AZMAX1(sum(StalkStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    StalkRsrvElms_pft(NE,NZ)       = AZMAX1(sum(StalkRsrvElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    HuskStrutElms_pft(NE,NZ)       = AZMAX1(sum(HuskStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    EarStrutElms_pft(NE,NZ)        = AZMAX1(sum(EarStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    GrainStrutElms_pft(NE,NZ)      = AZMAX1(sum(GrainStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))
    CanopyNonstElms_pft(NE,NZ)     = AZMAX1(sum(CanopyNonstElms_brch(NE,1:NumOfBranches_pft(NZ),NZ)))

    ShootNonstElms_pft(NE,NZ)       = CanopyNonstElms_pft(NE,NZ)    +ShootNoduleElms_pft(NE,NZ)
    ShootleafElms_pft(NE,NZ)        = LeafStrutElms_pft(NE,NZ)*FracLeafShethElmAlloc2Litr(NE,k_fine_comp)
    ShootFineNonLeafElms_pft(NE,NZ) = HuskStrutElms_pft(NE,NZ)+EarStrutElms_pft(NE,NZ)+GrainStrutElms_pft(NE,NZ) &
      +PetolShethStrutElms_pft(NE,NZ)*FracPetolShethAlloc2Litr(NE,k_fine_comp)
    ShootWoodyElms_pft(NE,NZ)       = StalkStrutElms_pft(NE,NZ) +StalkRsrvElms_pft(NE,NZ)+PetolShethStrutElms_pft(NE,NZ) &
      *FracPetolShethAlloc2Litr(NE,k_woody_comp)+LeafStrutElms_pft(NE,NZ)*FracLeafShethElmAlloc2Litr(NE,k_woody_comp)
  ENDDO
  
!
!     TOTAL STANDING DEAD
!
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
! standing dead are made up of different kinetic components
!add standing dead to shoot
  DO NE=1,NumPlantChemElms
    StandDeadStrutElms_pft(NE,NZ)=sum(StandDeadCompKElms_pft(NE,1:jsken,NZ))    
  ENDDO

  if(present(canopyE))then
    canopyE=ShootElms_pft(:,NZ)+ShootNoduleElms_pft(:,NZ)+StandDeadStrutElms_pft(:,NZ)
  endif
  call PrintInfo('end '//subname)
  END associate
  end subroutine SumCanopyBiome
!----------------------------------------------------------------------------------------------------
  subroutine SumPlantBiomStates(yearIJ,NZ,header,tvegE,canopyE,RootE)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  character(len=*),intent(in) :: header
  real(r8), optional, intent(out) :: tvegE(NumPlantChemElms)
  real(r8), optional, intent(out) :: canopyE(NumPlantChemElms)
  real(r8), optional, intent(out) :: RootE(NumPlantChemElms)

  real(r8) :: tvegE1(NumPlantChemElms)
  real(r8) :: canopyE1(NumPlantChemElms)
  real(r8) :: RootE1(NumPlantChemElms)

!     begin_execution
  

  call SumCanopyBiome(yearIJ,NZ,canopyE1)

  call SumRootBiome(yearIJ,NZ,RootE1)

  if(present(tvegE))tvegE=canopyE1+RootE1+plt_biom%SeasonalNonstElms_pft(:,NZ)
  if(present(canopyE))canopyE=canopyE1
  if(present(rootE))RootE=RootE1

  END subroutine SumPlantBiomStates

!----------------------------------------------------------------------------------------------------
  subroutine SumPlantBranchBiome(NB,NZ)
  implicit none
  integer, intent(in) :: NB,NZ
  integer :: NE,K
  associate(                                                             &
    MaxNodesPerBranch1         => pltpar%MaxNodesPerBranch1             ,& !input  :maximum number of canopy nodes, 25
    NumOfBranches_pft          => plt_morph%NumOfBranches_pft           ,& !input  :number of branches,[-]
    LeafStrutElms_brch         => plt_biom%LeafStrutElms_brch           ,& !input  :branch leaf structural element mass, [g d-2]
    PetolShethStrutElms_brch   => plt_biom%PetolShethStrutElms_brch     ,& !input  :branch sheath structural element, [g d-2]
    StalkStrutElms_brch        => plt_biom%StalkStrutElms_brch          ,& !input  :branch stalk structural element mass, [g d-2]
    StalkRsrvElms_brch         => plt_biom%StalkRsrvElms_brch           ,& !input  :branch reserve element mass, [g d-2]
    HuskStrutElms_brch         => plt_biom%HuskStrutElms_brch           ,& !input  :branch husk structural element mass, [g d-2]
    EarStrutElms_brch          => plt_biom%EarStrutElms_brch            ,& !input  :branch ear structural chemical element mass, [g d-2]
    GrainStrutElms_brch        => plt_biom%GrainStrutElms_brch          ,& !input  :branch grain structural element mass, [g d-2]
    CanopyNonstElms_brch       => plt_biom%CanopyNonstElms_brch         ,& !input  :branch nonstructural element, [g d-2]
    CPOOL3_node                => plt_photo%CPOOL3_node                 ,& !input  :minimum sink strength for nonstructural C transfer, [g d-2]
    CPOOL4_node                => plt_photo%CPOOL4_node                 ,& !input  :leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
    iPlantPhotosynsType_pft    => plt_photo%iPlantPhotosynsType_pft     ,& !input  :plant photosynthetic type (C3 or C4),[-]
    CMassHCO3BundleSheath_node => plt_photo%CMassHCO3BundleSheath_node  ,& !input  :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CMassCO2BundleSheath_node  => plt_photo%CMassCO2BundleSheath_node   ,& !input  :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    C4PhotoShootNonstC_brch    => plt_biom%C4PhotoShootNonstC_brch      ,& !inoput :branch shoot nonstrucal elelment, [g d-2]
    ShootElms_brch             => plt_biom%ShootElms_brch                & !output :branch shoot element mass, [g d-2]
  )

  DO NE=1,NumPlantChemElms
    ShootElms_brch(NE,NB,NZ)=LeafStrutElms_brch(NE,NB,NZ) &
      +PetolShethStrutElms_brch(NE,NB,NZ)+StalkStrutElms_brch(NE,NB,NZ)+StalkRsrvElms_brch(NE,NB,NZ) &
      +HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ)+GrainStrutElms_brch(NE,NB,NZ) &
      +CanopyNonstElms_brch(NE,NB,NZ)    
  ENDDO

  !add C4 specific reserve carbon
  C4PhotoShootNonstC_brch(NB,NZ)=0._r8
  IF(iPlantPhotosynsType_pft(NZ).EQ.ic4_photo)THEN      
    D3251: DO K=1,MaxNodesPerBranch1
      C4PhotoShootNonstC_brch(NB,NZ)=C4PhotoShootNonstC_brch(NB,NZ)+CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ) &
        +CMassCO2BundleSheath_node(K,NB,NZ)+CMassHCO3BundleSheath_node(K,NB,NZ)
    ENDDO D3251
    ShootElms_brch(ielmc,NB,NZ)=ShootElms_brch(ielmc,NB,NZ)+C4PhotoShootNonstC_brch(NB,NZ)
  ENDIF

  end associate
  end subroutine SumPlantBranchBiome

!----------------------------------------------------------------------------------------------------
  subroutine ZeroGrosub(NP)
  
  implicit none
  integer, intent(in) :: NP
  integer :: NZ,K,L,M,NE

  associate(                                                   &
    Myco_pft               => plt_morph%Myco_pft              ,& !input  :mycorrhizal type (no or yes),[-]
    NU                     => plt_site%NU                     ,& !input  :current soil surface layer number, [-]
    MaxSoilLays4Root_pft   => plt_morph%MaxSoilLays4Root_pft  ,& !input  :maximum soil layer number for all root axes,[-]
    MaxNumRootLays         => plt_site%MaxNumRootLays         ,& !input  :maximum root layer number,[-]
    Root1stActStruct_pvr   => plt_biom%Root1stActStruct_pvr   ,& !output :derived primary active root element mass
    Root1stLigStruct_pvr   => plt_biom%Root1stLigStruct_pvr   ,& !output :derived primary lignified root element mass
    RootMedStruct_pvr      => plt_biom%RootMedStruct_pvr      ,& !output :derived medium root element mass
    RootMyco2ndStrutElms_pvr => plt_biom%RootMyco2ndStrutElms_pvr ,& !output :derived secondary root/mycorrhizal element mass
    RootSinkWeight_pvr     => plt_morph%RootSinkWeight_pvr    ,& !output :root-layer sink weighting
    RootAtmGasConductance_rpvr => plt_rbgc%RootAtmGasConductance_rpvr ,& !output :root-atmosphere gas conductance
    RootCO2Autor_pvr       => plt_rbgc%RootCO2Autor_pvr       ,& !input  :root respiration constrained by O2, [g d-2 h-1]
    CO2NetFix_pft          => plt_bgcr%CO2NetFix_pft          ,& !output :canopy net CO2 exchange, [gC d-2 h-1]
    RCanMaintDef_CO2_pft   => plt_bgcr%RCanMaintDef_CO2_pft   ,& !output :canopy maintenance respiraiton deficit as CO2, [gC d-2]    
    RootMaintDef_CO2_pvr   => plt_bgcr%RootMaintDef_CO2_pvr   ,& !output :plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]        
    NodulInfectElms_pft    => plt_bgcr%NodulInfectElms_pft    ,& !output :nodule infection chemical element mass, [g d-2]
    Soil2RootMycoExudE_pft => plt_rbgc%Soil2RootMycoExudE_pft ,& !output :total root uptake (+ve) - exudation (-ve) of dissolved element, [g d-2 h-1]
    NH3Dep2Can_pft         => plt_bgcr%NH3Dep2Can_pft         ,& !output :canopy NH3 flux, [g d-2 h-1]
    GrossResp_pft          => plt_bgcr%GrossResp_pft          ,& !output :total plant respiration, [gC d-2 ]
    GrossCO2Fix_pft        => plt_bgcr%GrossCO2Fix_pft        ,& !output :total gross CO2 fixation, [gC d-2 ]
    CO2FixCL_pft           => plt_rbgc%CO2FixCL_pft           ,& !output :Rubisco-limited CO2 fixation, [gC d-2 h-1]
    RootN2Fix_pft          => plt_rbgc%RootN2Fix_pft          ,& !outut  :total root N2 fixation, [g d-2 h-1]    
    CanopyN2Fix_pft        => plt_rbgc%CanopyN2Fix_pft        ,& !outut  :total canopy N2 fixation, [g d-2 h-1]        
    CO2FixLL_pft           => plt_rbgc%CO2FixLL_pft           ,& !output :Light-limited CO2 fixation, [gC d-h2 h-1]
    RootNutUptakeN_pft     => plt_rbgc%RootNutUptakeN_pft     ,& !output :total N uptake by plant roots, [gN d-h2 h-1]
    RootNutUptakeP_pft     => plt_rbgc%RootNutUptakeP_pft     ,& !output :total P uptake by plant roots, [gP d-h2 h-1]
    CanopyGrosRCO2_pft     => plt_bgcr%CanopyGrosRCO2_pft     ,& !output :canopy plant+nodule autotrophic respiraiton, [gC d-2]
    CanopyResp_brch        => plt_bgcr%CanopyResp_brch        ,& !inoput :canopy respiration for a branch, [gC d-2 h-1]       
    SSXferElms_pft         => plt_bgcr%SSXferElms_pft         ,& !inoput :flux export from seasonal storage, [g h-1 d-2] 
    SSXfer2ShootElms_pft   => plt_bgcr%SSXfer2ShootElms_pft   ,& !inoput :flux export from seasonal storage to shoot, [g h-1 d-2]            
    PlantElmDistLoss_pft   => plt_distb%PlantElmDistLoss_pft  ,& !ouput  :plant element loss due to disturbance, [g d-2 h-1]
    RootLost2Fire_pft      => plt_distb%RootLost2Fire_pft     ,& !inoput :plant root biomass lost by fire, [g d-2 h-1]    
    LitrfallElms_pvr       => plt_bgcr%LitrfallElms_pvr       ,& !output :plant LitrFall element, [g d-2 h-1]
    LitrFallElms_brch      => plt_bgcr%LitrFallElms_brch      ,& !inoput :litterfall from the branch, [g d-2 h-1]        
    GPP_brch               => plt_rbgc%GPP_brch               ,& !output :GPP over branch, [gC d-2 h-1]    
    RootHPO4Uptake_pft     => plt_rbgc%RootHPO4Uptake_pft     ,& !output :total root uptake of HPO4, [g d-2 h-1]
    RootH2PO4Uptake_pft    => plt_rbgc%RootH2PO4Uptake_pft    ,& !output :total root uptake of PO4, [g d-2 h-1]
    RootNH4Uptake_pft      => plt_rbgc%RootNH4Uptake_pft      ,& !output :total root uptake of NH4, [g d-2 h-1]
    RootNO3Uptake_pft      => plt_rbgc%RootNO3Uptake_pft      ,& !output :total root uptake of NO3, [g d-2 h-1]
    PARSunlit_pft          => plt_photo%PARSunlit_pft         ,& !output :PAR absorbed by sunlit leaf, [umol m-2 s-1]
    PARSunsha_pft          => plt_photo%PARSunsha_pft         ,& !output :PAR absorbed by sun-shaded leaf, [umol m-2 s-1]
    CH2OSunlit_pft         => plt_photo%CH2OSunlit_pft        ,& !output :carbon fixation by sun-lit leaf, [gC d-2 h-1]
    CH2OSunsha_pft         => plt_photo%CH2OSunsha_pft        ,& !output :carbon fixation by sun-shaded leaf, [gC d-2 h-1]   
    ShootRootXferElm_pft   => plt_bgcr%ShootRootXferElm_pft   ,& !inoput :shoot-root nonstructural element transfer, [ g d-2 h-1]    
    trcs_deadroot2soil_pvr => plt_rbgc%trcs_deadroot2soil_pvr ,& !inoput :gases released to soil upong dying roots, [g d-2 h-1]    
    RNodeInitiate_pft      => plt_rbgc%RNodeInitiate_pft      ,& !inoput :node initiation rate, [h-1]
    RLeafAppear_pft        => plt_rbgc%RLeafAppear_pft        ,& !inoput :leaf appearing rate, [h-1]
    CO2Intra_pft           => plt_photo%CO2Intra_pft          ,& !inoput :leaf-area-weighted intracellular CO2 sum, [umol mol-1 m2 d-2]
    CO2IntraScal_pft       => plt_photo%CO2IntraScal_pft      ,& !inoput :intracellular CO2 leaf-area weight, [m2 d-2]
    fNCLFW_brch            => plt_pheno%fNCLFW_brch           ,& !output : NC ratio of growing leaf on branch, [gN/gC]
    fPCLFW_brch            => plt_pheno%fPCLFW_brch           ,& !output : PC ratio of growing leaf on branch, [gP/gC]
    fNCLFW_pft             => plt_pheno%fNCLFW_pft            ,& !output : NC ratio of growing leaf, [gN/gC]
    RNFixCO2_pft           => plt_bgcr%RNFixCO2_pft           ,& !inoput :CO2 respired by nodules, [gC d-2]    
    fPCLFW_pft             => plt_pheno%fPCLFW_pft             & !output : PC ratio of growing leaf, [gP/gC]
  )
  
  plt_rbgc%trcs_Soil2plant_uptake_vr=0._r8

  D9980: DO NZ=1,NP
    ! Excluded layers can retain legacy mass, but have no active root diagnostics.
    ! Clear shared working arrays so neither another column nor a restart leaves
    ! stale values in history. The underlying structural and reserve pools are kept.
    Root1stActStruct_pvr(:,MaxNumRootLays+1:,NZ) = 0._r8
    Root1stLigStruct_pvr(:,MaxNumRootLays+1:,NZ) = 0._r8
    RootMedStruct_pvr(:,MaxNumRootLays+1:,NZ) = 0._r8
    RootMyco2ndStrutElms_pvr(:,:,MaxNumRootLays+1:,NZ) = 0._r8
    RootSinkWeight_pvr(MaxNumRootLays+1:,NZ) = 0._r8
    RootAtmGasConductance_rpvr(:,:,MaxNumRootLays+1:,NZ) = 0._r8
    CO2Intra_pft(NZ) = 0._R8
    CO2IntraScal_pft(NZ) = 0._R8
    plt_rbgc%NH3Dep2Can_brch(:,NZ) = 0._R8  
    RLeafAppear_pft(NZ) = 0._r8
    RNodeInitiate_pft(NZ) = 0._r8
    SSXfer2ShootElms_pft(:,NZ)                  = 0._r8
    SSXferElms_pft(:,NZ)                        = 0._r8
    trcs_deadroot2soil_pvr(:,:,NZ)              = 0._r8  
    RNFixCO2_pft(NZ)                             = 0._r8
    fNCLFW_pft(NZ)                               = 0._r8
    fPCLFW_pft(NZ)                               = 0._r8
    plt_photo%CanopyVcMaxRubisco25C_pft(NZ)      = 0._r8
    plt_photo%CanopyVoMaxRubisco25C_pft(NZ)      = 0._r8
    plt_photo%CanopyVcMaxPEP25C_pft(NZ)          = 0._r8
    plt_photo%ProteinCperm2LeafArea_node(:,:,NZ) = 0._r8
    plt_photo%ElectronTransptJmax25C_pft(NZ)     = 0._r8
    plt_bgcr%RCO2Nodule_pvr(:,NZ)                = 0._r8
    plt_biom%SeedPlantedElm_pft(:,NZ)            = 0._r8
    ShootRootXferElm_pft(:,NZ)                   = 0._r8
    plt_biom%LeafC3ChlCperm2LA_pft(NZ)           = 0._r8
    plt_biom%LeafC4ChlCperm2LA_pft(NZ)           = 0._r8
    plt_biom%LeafRubiscoCperm2LA_pft(NZ)         = 0._r8
    plt_biom%LeafPEPCperm2LA_pft(NZ)             = 0._r8
    plt_biom%SpecificLeafArea_pft(NZ)            = 0._r8
    plt_biom%LeafProteinCperm2LA_pft(NZ)         = 0._r8
    plt_bgcr%Xfer2RootsC_pft(NZ)                 = 0._R8
    D1: DO L=0,MaxNumRootLays
      DO K=1,pltpar%NumOfPlantLitrCmplxs
        DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitrfallElms_pvr(NE,M,K,L,NZ)=0._r8
          ENDDO
        ENDDO
      ENDDO
    ENDDO D1
    fNCLFW_brch(:,NZ)                                           = 0._r8
    fPCLFW_brch(:,NZ)                                           = 0._r8
    LitrFallElms_brch(:,:,NZ)                                   = 0._r8
    PARSunsha_pft(NZ)                                           = 0._r8
    PARSunlit_pft(NZ)                                           = 0._r8
    CH2OSunlit_pft(NZ)                                          = 0._r8
    CH2OSunsha_pft(NZ)                                          = 0._r8
    RootMaintDef_CO2_pvr(:,:,NZ)                                = 0._r8
    RootCO2Autor_pvr(1:Myco_pft(NZ),NU:MaxSoilLays4Root_pft(NZ),NZ) = 0._r8
    NH3Dep2Can_pft(NZ)                                          = 0._r8
    GrossResp_pft(NZ)                                           = 0._r8
    GrossCO2Fix_pft(NZ)                                         = 0._r8
    CO2NetFix_pft(NZ)                                           = 0._r8
    RCanMaintDef_CO2_pft(NZ)                                    = 0._r8
    CanopyGrosRCO2_pft(NZ)                                      = 0._r8
    CanopyResp_brch(:,NZ)                                       = 0._r8
    PlantElmDistLoss_pft(1:NumPlantChemElms,NZ)                 = 0._r8
    RootLost2Fire_pft(1:NumPlantChemElms,NZ)                    = 0._R8
    Soil2RootMycoExudE_pft(1:NumPlantChemElms,NZ)               = 0._r8
    NodulInfectElms_pft(1:NumPlantChemElms,NZ)                  = 0._r8
    CO2FixCL_pft(NZ)                                            = 0._r8
    CO2FixLL_pft(NZ)                                            = 0._r8
    RootNutUptakeN_pft(NZ)                                      = 0._r8
    RootNutUptakeP_pft(NZ)                                      = 0._r8
    RootN2Fix_pft(NZ)                                           = 0._r8
    CanopyN2Fix_pft(NZ)                                         = 0._R8
    GPP_brch(:,NZ)                                              = 0._R8
    RootNH4Uptake_pft(NZ)                                       = 0._r8
    RootNO3Uptake_pft(NZ)                                       = 0._r8
    RootH2PO4Uptake_pft(NZ)                                     = 0._r8
    RootHPO4Uptake_pft(NZ)                                      = 0._r8
  ENDDO D9980
  end associate
  end subroutine ZeroGrosub

!----------------------------------------------------------------------------------------------------
  subroutine SumRootBiome(yearIJ,NZ,massroot)
  ! Aggregate per-axis structural and reserve pools into root/mycorrhizal totals.

  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), optional, intent(out) :: massroot(NumPlantChemElms)  
  character(len=*), parameter :: subname='SumRootBiome'

  integer :: NE,N,L
  real(r8) :: massr1st1(NumPlantChemElms)
  real(r8) :: massr2nd1(NumPlantChemElms)
  real(r8) :: massnonst1(NumPlantChemElms)
  associate(                                                          &
    PlantPopuLive_pft         => plt_site%PlantPopuLive_pft         , & !inoput :plant population, [d-2]  
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    Myco_pft                  => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    NMaxRootBotLayer_pft      => plt_morph%NMaxRootBotLayer_pft      ,& !input  :maximum soil layer number for all root axes, [-]    
    MaxNumRootLays            => plt_site%MaxNumRootLays             ,& !input  :maximum root layer number,[-]
    MaxSoilLays4Root_pft      => plt_morph%MaxSoilLays4Root_pft      ,& !input  :maximum soil layer number for all root axes,[-]
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft        ,& !input  :N2 fixation type,[-]
    NumStructuralRootAxes_pft => plt_morph%NumStructuralRootAxes_pft ,& !input  :number of structural root axes,[-]
    RootMediumStructElms_rpvr => plt_biom%RootMediumStructElms_rpvr  ,& !inoput :root layer element for medium size root axes, [g d-2]
    RootMedStruct_pvr         => plt_biom%RootMedStruct_pvr          ,& !inoput :root layer element biomass for medium size roots, [g d-2]    
    Root1stActStructElms_rpvr => plt_biom%Root1stActStructElms_rpvr  ,& !inoput :Root layer primary axes Active zone structrual element, [g d-2]    
    Root1stLigStructElms_rpvr => plt_biom%Root1stLigStructElms_rpvr  ,& !inoput :root layer lignified zone element in primary axes, [g d-2]    
    RootNodulStrutElms_rpvr   => plt_biom%RootNodulStrutElms_rpvr    ,& !input  :root layer nodule element, [g d-2]
    RootNodulNonstElms_rpvr   => plt_biom%RootNodulNonstElms_rpvr    ,& !input  :root layer nonstructural element, [g d-2]
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr  ,& !input  :root layer element primary axes, [g d-2]
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr  ,& !input  :root layer element secondary axes, [g d-2]
    RootMyco2ndStrutElms_pvr  => plt_biom%RootMyco2ndStrutElms_pvr   ,& !output :secondary structural elements summed over active NR, [g element cell-1]
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     ,& !input  :root layer nonstructural element, [g d-2]
    RootElms_pft              => plt_biom%RootElms_pft               ,& !output :plant root element mass, [g d-2]    
    RootNoduleElms_pft        => plt_biom%RootNoduleElms_pft         ,& !output :plant root nodule element mass, [g d-2]    
    RootMycoNonstElms_pft     => plt_biom%RootMycoNonstElms_pft      ,& !output :nonstructural root-myco chemical element, [g d-2]
    RootStrutElms_pft         => plt_biom%RootStrutElms_pft          ,& !output :plant root structural element mass, [g d-2]
    Root1stActStruct_pvr      => plt_biom%Root1stActStruct_pvr       ,& !output :active zone of primary roots, [g d-2]
    Root1stLigStruct_pvr      => plt_biom%Root1stLigStruct_pvr       ,& !output :lignifed primary root biomass, [g d-2]    
    RootMycoMassElm_pvr       => plt_biom%RootMycoMassElm_pvr         & !output :root biomass in chemical elements, [g d-2]
  )
  call PrintInfo('beg '//subname)
  massr1st1=0._r8;massr2nd1=0._r8
  RootMyco2ndStrutElms_pvr(:,:,:,NZ)=0._r8
  RootElms_pft(:,NZ)=0._r8
  RootMycoNonstElms_pft(:,:,NZ)=0._r8
  DO NE=1,NumPlantChemElms
    DO L=1,MaxNumRootLays
      DO N=1,Myco_pft(NZ)
        RootMyco2ndStrutElms_pvr(NE,N,L,NZ) = &
          sum(RootMyco2ndStrutElms_rpvr(NE,N,L,1:NumStructuralRootAxes_pft(NZ),NZ))
        RootMycoMassElm_pvr(NE,N,L,NZ) = RootMyco2ndStrutElms_pvr(NE,N,L,NZ)+RootMycoNonstElms_rpvr(NE,N,L,NZ)
        RootMycoNonstElms_pft(NE,N,NZ) = RootMycoNonstElms_pft(NE,N,NZ)+RootMycoNonstElms_rpvr(NE,N,L,NZ)
        massr2nd1(NE)                  = massr2nd1(NE)+RootMyco2ndStrutElms_pvr(NE,N,L,NZ)
      ENDDO  
      RootMedStruct_pvr(NE,L,NZ)    = SUM(RootMediumStructElms_rpvr(NE,L,1:NumStructuralRootAxes_pft(NZ),NZ))
      Root1stActStruct_pvr(NE,L,NZ) = SUM(Root1stActStructElms_rpvr(NE,L,1:NumStructuralRootAxes_pft(NZ),NZ))
      Root1stLigStruct_pvr(NE,L,NZ) = SUM(Root1stLigStructElms_rpvr(NE,L,1:NumStructuralRootAxes_pft(NZ),NZ))
      massr1st1(NE)=massr1st1(NE)+sum(RootMyco1stStrutElms_rpvr(NE,L,1:NumStructuralRootAxes_pft(NZ),NZ))+RootMedStruct_pvr(NE,L,NZ)        
      RootMycoMassElm_pvr(NE,ipltroot,L,NZ)= RootMycoMassElm_pvr(NE,ipltroot,L,NZ)+sum(RootMyco1stStrutElms_rpvr(NE,L,1:NumStructuralRootAxes_pft(NZ),NZ))+RootMedStruct_pvr(NE,L,NZ)        
    ENDDO
    
    RootStrutElms_pft(NE,NZ)=massr1st1(NE)+massr2nd1(NE)
    massnonst1(NE)      = sum(RootMycoNonstElms_pft(NE,1:Myco_pft(NZ),NZ))
    
    RootElms_pft(NE,NZ) = massr1st1(NE)+massr2nd1(NE)+massnonst1(NE)
    if(RootElms_pft(NE,NZ)<0._r8)then
    write(945,*)yearIJ%I*1000+yearIJ%J/24., massr1st1(NE),massr2nd1(NE),massnonst1(NE),'NZ',NZ,NE,NumStructuralRootAxes_pft(NZ)
    endif
    !add reserve to struct
    RootNoduleElms_pft(NE,NZ)=0._r8
    if(is_plant_N2fix(iPlantNfixType_pft(NZ)) .and. is_root_N2fix(iPlantNfixType_pft(NZ)))THEN
      RootNoduleElms_pft(NE,NZ) = RootNoduleElms_pft(NE,NZ)+sum(RootNodulStrutElms_rpvr(NE,NU:NMaxRootBotLayer_pft(NZ),NZ))+sum(RootNodulNonstElms_rpvr(NE,NU:NMaxRootBotLayer_pft(NZ),NZ))
    endif      
  ENDDO
  
  if(present(massroot))massroot=RootElms_pft(:,NZ)+RootNoduleElms_pft(:,NZ)
  call PrintInfo('end '//subname)
  end associate
  end subroutine SumRootBiome

!----------------------------------------------------------------------------------------------------
  subroutine EnterPlantBalance(yearIJ,NP)
  !
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NP
  INTEGER :: NZ,NE
  real(r8) :: mass_finale(NumPlantChemElms)

  DO NZ=1,NP
    
    DO NE=1,NumPlantChemElms
      plt_bgcr%RootShootExch_pvr(NE,:,NZ)=0._r8
      plt_biom%TotEndVegE_pft(NE,NZ) = plt_biom%RootElms_pft(NE,NZ)+plt_biom%ShootElms_pft(NE,NZ)+&
        plt_biom%SeasonalNonstElms_pft(NE,NZ)+plt_biom%StandDeadStrutElms_pft(NE,NZ)+ plt_biom%ShootNoduleElms_pft(NE,NZ)  + &
        plt_biom%RootNoduleElms_pft(NE,NZ) 

!      if(plt_biom%RootElms_pft(NE,NZ)<0._r8)then
!        write(944,*)yearIJ%I*1000+yearIJ%J/24.,NE,plt_biom%RootElms_pft(NE,NZ),plt_biom%ShootElms_pft(NE,NZ),&
!        plt_biom%SeasonalNonstElms_pft(NE,NZ),plt_biom%StandDeadStrutElms_pft(NE,NZ), plt_biom%ShootNoduleElms_pft(NE,NZ), &
!        plt_biom%RootNoduleElms_pft(NE,NZ) 
!      endif
      plt_biom%RootElmsBeg_pft(NE,NZ)           = plt_biom%RootElms_pft(NE,NZ)
      plt_biom%StandDeadStrutElmsBeg_pft(NE,NZ) = plt_biom%StandDeadStrutElms_pft(NE,NZ)
      plt_biom%ShootElmsBeg_pft(NE,NZ)          = plt_biom%ShootElms_pft(NE,NZ)
      plt_biom%SeasonalNonstElmsbeg_pft(NE,NZ)  = plt_biom%SeasonalNonstElms_pft(NE,NZ)
      plt_biom%ShootNoduleElmsBeg_pft(NE,NZ)    = plt_biom%ShootNoduleElms_pft(NE,NZ)
      plt_biom%TotBegVegE_pft(NE,NZ)            = plt_biom%TotEndVegE_pft(NE,NZ)
      plt_biom%RootNoduleElmsBeg_pft(NE,NZ)     = plt_biom%RootNoduleElms_pft(NE,NZ)
      plt_distb%FireLossE_pft(NE,NZ)            = 0._r8
      plt_bgcr%NH3Dep2Can_pft(NZ)               = 0._r8
      plt_bgcr%SurfLitrfallElms_pft(NE,NZ)      = 0._r8
    ENDDO

  ENDDO

  !   UPDATE PLANT PHENOLOGY IN 'HFUNC'
  !     zero out plant hourly fluxes
  call ZeroGrosub(NP)  

  end subroutine EnterPlantBalance


!----------------------------------------------------------------------------------------------------
  subroutine ExitPlantBalance(yearIJ,NP)
  !
  !Description:
  ! mass balance check for plants
  ! dplantE/dt=input-output
  !
  ! plantE = root+shoot+standing dead
  ! C input: photosynthesis, nodule infection, DOM uptake/exudation
  ! C output: respiration, harvest
  
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NP
  integer :: NZ
  character(len=*), parameter :: subname='ExitPlantBalance'

  DO NZ=1,NP      
    call CheckPlantBalanceZ(yearIJ,NZ,'')
  ENDDO

  end subroutine ExitPlantBalance

!----------------------------------------------------------------------------------------------------
  subroutine CheckPlantBalanceZ(yearIJ,NZ,header)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  character(len=*), intent(in) :: header
  character(len=*), parameter :: subname='CheckPlantBalanceZ'
  integer :: NE,I,J,K
  real(r8) :: balE(NumPlantChemElms)
  real(r8) :: err_rel,dGPP  
! --- Jing Tao ---
! set true only when guard_floor_dying_stand is on AND total vegetation carbon has
! collapsed below VEGC_BAL_FLOOR_dying, which makes the RELATIVE balance test
! meaningless. Declared here rather than module-scope so it cannot leak between calls.
  logical :: lveg_too_small_for_relative_test
  logical :: lbalC_abort
  associate(                                                          &
    NY                        => plt_site%NY                         ,& !input  : Y coordinate
    NX                        => plt_site%NX                         ,& !input  : X coordinate
    iYearCurrent              => plt_site%iYearCurrent               ,& !input  :current year,[-]  
    PlantElmBalCum_pft        => plt_site%PlantElmBalCum_pft         ,& !inoput :cumulative plant element balance, [g d-2]    
    NH3Dep2Can_pft            => plt_bgcr%NH3Dep2Can_pft             ,& !input  :canopy NH3 flux, [gN d-2 h-1]    
    GrossCO2Fix_pft           => plt_bgcr%GrossCO2Fix_pft            ,& !input  :total gross CO2 fixation, [gC d-2 ]
    GrossResp_pft             => plt_bgcr%GrossResp_pft              ,& !input  :total plant respiration, [gC d-2 ]
    NodulInfectElms_pft       => plt_bgcr%NodulInfectElms_pft        ,& !input  :nodule infection chemical element mass, [g d-2]
    LitrfallElms_pft          => plt_bgcr%LitrfallElms_pft           ,& !input  :plant element LitrFall, [g d-2 h-1]
    canopy_growth_pft         => plt_rbgc%canopy_growth_pft          ,& !inoput :canopy structural C growth rate, [gC d-2 h-1]    
    RCanMaintDef_CO2_pft      => plt_bgcr%RCanMaintDef_CO2_pft       ,& !inoput :canopy maintenance respiraiton deficit as CO2, [gC d-2 h-1]    
    Soil2RootMycoExudE_pft    => plt_rbgc%Soil2RootMycoExudE_pft     ,& !input  :total root uptake (+ve) - exudation (-ve) of dissolved element, [g d-2 h-1]
    RootElmsBeg_pft           => plt_biom%RootElmsBeg_pft            ,& !input  :plant root element at previous time step, [g d-2]
    RootElms_pft              => plt_biom%RootElms_pft               ,& !input  :plant root element mass, [g d-2]
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft      ,& !input  :plant stored nonstructural element at current step, [g d-2]
    SeasonalNonstElmsbeg_pft  => plt_biom%SeasonalNonstElmsbeg_pft   ,& !input  :plant stored nonstructural element at prev step, [g d-2]
    StandDeadStrutElmsBeg_pft => plt_biom%StandDeadStrutElmsBeg_pft  ,& !input  :standing dead element at previous time step, [g d-2]
    StandDeadStrutElms_pft    => plt_biom%StandDeadStrutElms_pft     ,& !input  :standing dead element, [g d-2]
    ShootElmsBeg_pft          => plt_biom%ShootElmsBeg_pft           ,& !input  :previous whole plant shoot element mass, [g d-2]
    RootNoduleElmsBeg_pft     => plt_biom%RootNoduleElmsBeg_pft      ,& !input  :previous time step plant root nodule element mass, [g d-2]       
    RootNoduleElms_pft        => plt_biom%RootNoduleElms_pft         ,& !input  :plant root nodule element mass, [g d-2]   
    ShootNoduleElmsBeg_pft    => plt_biom%ShootNoduleElmsBeg_pft     ,& !input  :previous time step canopy nodule element mass, [g d-2]             
    ShootNoduleElms_pft       => plt_biom%ShootNoduleElms_pft        ,& !input  :current time canopy nodule element mass, [g d-2]         
    RootNutUptakeN_pft        => plt_rbgc%RootNutUptakeN_pft         ,& !input  :total N uptake by plant roots, [gN d-h2 h-1]
    RootNutUptakeP_pft        => plt_rbgc%RootNutUptakeP_pft         ,& !input  :total P uptake by plant roots, [gP d-h2 h-1]    
    RootAutoCO2_pft           => plt_bgcr%RootAutoCO2_pft            ,& !input  :root autotrophic respiraiton, [gC d-2]    
    CanopyGrosRCO2_pft        => plt_bgcr%CanopyGrosRCO2_pft         ,& !input :canopy plant+nodule autotrophic respiraiton, [gC d-2]    
    ShootElms_pft             => plt_biom%ShootElms_pft              ,& !input  :current time whole plant shoot element mass, [g d-2]    
    PlantElmDistLoss_pft      => plt_distb%PlantElmDistLoss_pft      ,& !input  :plant loss to disturbance,    [g d-2 h-1]        
    RootN2Fix_pft             => plt_rbgc%RootN2Fix_pft              ,& !input  :total root N2 fixation, [gN d-2 h-1]    
    SurfLitrfallElms_pft      => plt_bgcr%SurfLitrfallElms_pft       ,& !input :surface litterfall, [g d-2 h-1]    
    CanopyN2Fix_pft           => plt_rbgc%CanopyN2Fix_pft            ,& !input  :total canopy N2 fixation, [gN d-2 h-1]        
    CanopyCi2CaRatio_pft      => plt_photo%CanopyCi2CaRatio_pft      ,& !input  :prior Ci:Ca ratio, [-]
    CanopyGasCO2_pft          => plt_photo%CanopyGasCO2_pft          ,& !inoput :canopy gaesous CO2 concentration, [umol mol-1]
    DynCi2CaRatio_pft         => plt_photo%DynCi2CaRatio_pft         ,& !output :dynamic Ci:Ca ratio, [-]
    SeedPlantedElm_pft        => plt_biom%SeedPlantedElm_pft         ,& !input  :seed biomass at planting, [g d-2] 
    iDayPlantHarvest_pft      => plt_distb%iDayPlantHarvest_pft      ,& !input  : day of plant harvest,[-]
    LitrfallAbvgElms_pft      => plt_bgcr%LitrfallAbvgElms_pft       ,& !input  :aboveground plant element LitrFall, [g d-2 h-1]
    LitrfallBlgrElms_pft      => plt_bgcr%LitrfallBlgrElms_pft       ,& !input  :belowground plant element LitrFall, [g d-2 h-1]
    CO2Intra_pft              => plt_photo%CO2Intra_pft              ,& !inoput :leaf-area-weighted intracellular CO2 sum, [umol mol-1 m2 d-2]
    CO2IntraScal_pft          => plt_photo%CO2IntraScal_pft          ,& !inoput :intracellular CO2 leaf-area weight, [m2 d-2]
    FireLossE_pft             => plt_distb%FireLossE_pft             ,& !input  :plant element lost by fire, [g d-2 h-1]
    TotBegVegE_pft            => plt_biom%TotBegVegE_pft             ,& !Input  :total vegetation carbon at the beginning of the time step,[g d-2]
    TotEndVegE_pft            => plt_biom%TotEndVegE_pft              & !output :total vegetation carbon at the end of the time step,[g d-2]
  )

  call PrintInfo('beg '//subname)
  
  ! Preserve the weighted sum for the Ci history diagnostic.
  IF(CO2IntraScal_pft(NZ).GT.1.E-5_R8.AND.CanopyGasCO2_pft(NZ).GT.0._r8)THEN
    DynCi2CaRatio_pft(NZ)=(CO2Intra_pft(NZ)/CO2IntraScal_pft(NZ))/CanopyGasCO2_pft(NZ)
  ELSE
    DynCi2CaRatio_pft(NZ)=  CanopyCi2CaRatio_pft(NZ)
  ENDIF
  call SumPlantBiome(yearIJ,NZ,header)

  I=yearIJ%I; J=yearIJ%J  

  DO NE=1,NumPlantChemElms
    TotEndVegE_pft(NE,NZ) = RootElms_pft(NE,NZ)+ShootElms_pft(NE,NZ)+SeasonalNonstElms_pft(NE,NZ)+StandDeadStrutElms_pft(NE,NZ)  + &
      ShootNoduleElms_pft(NE,NZ)  + RootNoduleElms_pft(NE,NZ)   

    balE(NE)=TotEndVegE_pft(NE,NZ)-TotBegVegE_pft(NE,NZ) &
      -NodulInfectElms_pft(NE,NZ)-Soil2RootMycoExudE_pft(NE,NZ) &
      +LitrfallElms_pft(NE,NZ)+PlantElmDistLoss_pft(NE,NZ)-SeedPlantedElm_pft(NE,NZ)
  ENDDO
  
  dGPP=sum(plt_rbgc%GPP_brch(:,NZ))
  balE(ielmc)=balE(ielmc)-GrossCO2Fix_pft(NZ)-GrossResp_pft(NZ)   
  balE(ielmn)=balE(ielmn)-RootNutUptakeN_pft(NZ)-CanopyN2Fix_pft(NZ)-RootN2Fix_pft(NZ)-NH3Dep2Can_pft(NZ)             
  balE(ielmp)=balE(ielmp)-RootNutUptakeP_pft(NZ)      

  NE=ielmc
  if(abs(balE(NE))>1.e-10_r8)then
    err_rel=safe_adb(balE(NE),AMAX1(TotEndVegE_pft(NE,NZ),TotBegVegE_pft(NE,NZ)))
  else
    err_rel=1.e-10_r8
  endif
! --- Jing Tao (2026-10-14, branch a2mc-anchor-0ffe5a06) ---
! err_rel is normalised by total vegetation carbon, which COLLAPSES as the stand
! dies, so a rounding-level absolute residual can read as a large relative error,
! or (via safe_adb's own near-zero-denominator clamp) as a spuriously SMALL one --
! either direction is unreliable below a floor on total vegetation carbon. Compute
! the regime flag and the single abort condition HERE, before the diagnostic-write
! gate below, so the gate and the abort test share ONE condition rather than two
! copies that can drift apart (the previous form let the gate require
! abs(err_rel)>1e-3 while the abort below the floor did not, which could make a
! corrected abort unreachable).
! Default-off: when guard_floor_dying_stand is .false. this reduces exactly to the
! original abs(balE)>1e-6 .and. abs(err_rel)>1e-3.
  lveg_too_small_for_relative_test = .false.
  if(guard_floor_dying_stand) lveg_too_small_for_relative_test = &
    (AMAX1(TotEndVegE_pft(NE,NZ),TotBegVegE_pft(NE,NZ)) < VEGC_BAL_FLOOR_dying)
  lbalC_abort = abs(balE(NE))>1.e-6_r8 .and. (lveg_too_small_for_relative_test .or. abs(err_rel)>1.e-3_r8)
  if(lbalC_abort .or. len_trim(header)>3)then      
    write(888,*)iYearCurrent*1000+I+J/24.,'pft=',NZ,'balC err err_rel',balE(ielmc),err_rel,plt_distb%iDayPlanting_pft(NZ),plt_distb%iDayPlantHarvest_pft(NZ)
    write(888,*)'NY, NX=',NY,NX,header
    write(888,*)'endc, begc        =',TotEndVegE_pft(NE,NZ),TotBegVegE_pft(NE,NZ),TotEndVegE_pft(NE,NZ)-TotBegVegE_pft(NE,NZ)
    write(888,*)'end, beg rootC    =',RootElms_pft(NE,NZ),RootElmsBeg_pft(NE,NZ),RootElms_pft(NE,NZ)-RootElmsBeg_pft(NE,NZ)
    write(888,*)'shootC            =',ShootElms_pft(NE,NZ),ShootElmsBeg_pft(NE,NZ),ShootElms_pft(NE,NZ)-ShootElmsBeg_pft(NE,NZ)
    write(888,*)'sstoreC           =',SeasonalNonstElms_pft(NE,NZ),SeasonalNonstElmsbeg_pft(NE,NZ),SeasonalNonstElms_pft(NE,NZ)-SeasonalNonstElmsbeg_pft(NE,NZ)
    write(888,*)'ssteadC           =',StandDeadStrutElms_pft(NE,NZ),StandDeadStrutElmsBeg_pft(NE,NZ),StandDeadStrutElms_pft(NE,NZ)-StandDeadStrutElmsBeg_pft(NE,NZ)
    write(888,*)'ShootNodulC       =',ShootNoduleElms_pft(NE,NZ),ShootNoduleElmsBeg_pft(NE,NZ),ShootNoduleElms_pft(NE,NZ)-ShootNoduleElmsBeg_pft(NE,NZ)
    write(888,*)'RootNoduleC       =',RootNoduleElms_pft(NE,NZ),RootNoduleElmsBeg_pft(NE,NZ),RootNoduleElms_pft(NE,NZ)-RootNoduleElmsBeg_pft(NE,NZ)    
    write(888,*)'nodulinfectC      =',NodulInfectElms_pft(NE,NZ)
    write(888,*)'rootexudC         =',Soil2RootMycoExudE_pft(NE,NZ)
    write(888,*)'litfallC, abg,blg =',LitrfallElms_pft(NE,NZ),LitrfallAbvgElms_pft(NE,NZ),LitrfallBlgrElms_pft(NE,NZ)
    write(888,*)'disturbC          =',PlantElmDistLoss_pft(NE,NZ)    
    write(888,*)'canMaintDef       =',RCanMaintDef_CO2_pft(NZ)
    write(888,*)'FireLossC         =',FireLossE_pft(NE,NZ)
    write(888,*)'CanoGroth         =',canopy_growth_pft(NZ)
    write(888,*)'GPP               =',GrossCO2Fix_pft(NZ),dGPP,GrossCO2Fix_pft(NZ)-dGPP
    write(888,*)'AR,rootAR,shootAR =',GrossResp_pft(NZ),RootAutoCO2_pft(NZ),CanopyGrosRCO2_pft(NZ)   
    write(888,*)'seed planted      =',plt_biom%SeedPlantedElm_pft(NE,NZ)
    write(888,*)'surf literfall C  =',SurfLitrfallElms_pft(NE,NZ)       
    write(888,*)('-',K=1,100)         
    !if(I/=plt_distb%iDayPlantHarvest_pft(NZ) .and. I/=plt_distb%iDayPlanting_pft(NZ)) &  
    if(lbalC_abort)call endrun('C balance error test failure in '//trim(mod_filename)//' at line',__LINE__)                    
  endif

  NE=ielmn
  if(abs(balE(NE))>1.e-10_r8)then            
    err_rel=safe_adb(balE(NE),AMAX1(TotEndVegE_pft(NE,NZ),TotBegVegE_pft(NE,NZ)))      
  else
    err_rel=1.e-10_r8
  endif

  if(abs(balE(NE))>1.e-6_r8 .and. abs(err_rel)>1.e-3_r8 .or. len_trim(header)>3)then      
    write(888,*)iYearCurrent*1000+I+J/24.,NZ,'balN',balE(NE),err_rel,plt_distb%iDayPlanting_pft(NZ),plt_distb%iDayPlantHarvest_pft(NZ)
    write(888,*)'NY,NX=',plt_site%NY,plt_site%NX,header
    write(888,*)'endN, begN        =',TotEndVegE_pft(NE,NZ),TotBegVegE_pft(NE,NZ),TotEndVegE_pft(NE,NZ)-TotBegVegE_pft(NE,NZ)
    write(888,*)'end/begrootN      =',RootElms_pft(NE,NZ),RootElmsBeg_pft(NE,NZ),RootElms_pft(NE,NZ)-RootElmsBeg_pft(NE,NZ)
    write(888,*)'end/begshootN     =',ShootElms_pft(NE,NZ),ShootElmsBeg_pft(NE,NZ),ShootElms_pft(NE,NZ)-ShootElmsBeg_pft(NE,NZ)
    write(888,*)'end/begsstoreN    =',SeasonalNonstElms_pft(NE,NZ),SeasonalNonstElmsbeg_pft(NE,NZ),SeasonalNonstElms_pft(NE,NZ)-SeasonalNonstElmsbeg_pft(NE,NZ)
    write(888,*)'end/begssteadN    =',StandDeadStrutElms_pft(NE,NZ),StandDeadStrutElmsBeg_pft(NE,NZ),StandDeadStrutElms_pft(NE,NZ)-StandDeadStrutElmsBeg_pft(NE,NZ)
    write(888,*)'end/begShootNodulN=',ShootNoduleElms_pft(NE,NZ),ShootNoduleElmsBeg_pft(NE,NZ),ShootNoduleElms_pft(NE,NZ)-ShootNoduleElmsBeg_pft(NE,NZ)
    write(888,*)'end/begRootNoduleN=',RootNoduleElms_pft(NE,NZ),RootNoduleElmsBeg_pft(NE,NZ),RootNoduleElms_pft(NE,NZ)-RootNoduleElmsBeg_pft(NE,NZ)    
    write(888,*)'nodulinfectN      =',NodulInfectElms_pft(NE,NZ)
    write(888,*)'rootexudN         =',Soil2RootMycoExudE_pft(NE,NZ)
    write(888,*)'litfallN, abg,blg =',LitrfallElms_pft(NE,NZ),LitrfallAbvgElms_pft(NE,NZ),LitrfallBlgrElms_pft(NE,NZ)
    write(888,*)'disturbN          =',PlantElmDistLoss_pft(NE,NZ)    
    write(888,*)'FireLossN         =',FireLossE_pft(NE,NZ)    
    write(888,*)'RootNuptk         =',RootNutUptakeN_pft(NZ)
    write(888,*)'CanopyNFix        =',CanopyN2Fix_pft(NZ)
    write(888,*)'CanopyNH3dep      =',NH3Dep2Can_pft(NZ)
    write(888,*)'RootNFix          =',RootN2Fix_pft(NZ)            
    write(888,*)'seed planted      =',plt_biom%SeedPlantedElm_pft(NE,NZ)     
    write(888,*)'surf literfall N =',SurfLitrfallElms_pft(NE,NZ)          
    write(888,*)('-',K=1,100)             
    if(abs(balE(NE))>1.e-6_r8 .and. abs(err_rel)>1.e-3_r8)call endrun('N balance error test failure in '//trim(mod_filename)//' at line',__LINE__)      
  endif

  NE=ielmp
  if(abs(balE(NE))>1.e-10_r8)then
    err_rel=safe_adb(balE(NE),AMAX1(TotEndVegE_pft(NE,NZ),TotBegVegE_pft(NE,NZ)))
  else
    err_rel=1.e-10_r8
  endif
  if(abs(balE(NE))>1.e-6_r8 .and. abs(err_rel)>1.e-3_r8)then  
    write(888,*)iYearCurrent*1000+I+J/24.,NZ,'P',balE(NE),err_rel,plt_distb%iDayPlanting_pft(NZ),plt_distb%iDayPlantHarvest_pft(NZ)
    write(888,*)'NY,NX=',plt_site%NY,plt_site%NX    
    write(888,*)'endP, begP       =',TotEndVegE_pft(NE,NZ),TotBegVegE_pft(NE,NZ),TotEndVegE_pft(NE,NZ)-TotBegVegE_pft(NE,NZ)
    write(888,*)'rootP            =',RootElms_pft(NE,NZ),RootElmsBeg_pft(NE,NZ),RootElms_pft(NE,NZ)-RootElmsBeg_pft(NE,NZ)
    write(888,*)'shootP           =',ShootElms_pft(NE,NZ),ShootElmsBeg_pft(NE,NZ),ShootElms_pft(NE,NZ)-ShootElmsBeg_pft(NE,NZ)
    write(888,*)'sstoreP          =',SeasonalNonstElms_pft(NE,NZ),SeasonalNonstElmsbeg_pft(NE,NZ),SeasonalNonstElms_pft(NE,NZ)-SeasonalNonstElmsbeg_pft(NE,NZ)
    write(888,*)'ssteadP          =',StandDeadStrutElms_pft(NE,NZ),StandDeadStrutElmsBeg_pft(NE,NZ),StandDeadStrutElms_pft(NE,NZ)-StandDeadStrutElmsBeg_pft(NE,NZ)
    write(888,*)'ShootNodulP      =',ShootNoduleElms_pft(NE,NZ),ShootNoduleElmsBeg_pft(NE,NZ),ShootNoduleElms_pft(NE,NZ)-ShootNoduleElmsBeg_pft(NE,NZ)
    write(888,*)'RootNoduleP      =',RootNoduleElms_pft(NE,NZ),RootNoduleElmsBeg_pft(NE,NZ),RootNoduleElms_pft(NE,NZ)-RootNoduleElmsBeg_pft(NE,NZ)    
    write(888,*)'nodulinfectP     =',NodulInfectElms_pft(NE,NZ)
    write(888,*)'rootexudP        =',Soil2RootMycoExudE_pft(NE,NZ)
    write(888,*)'litfallP,abg,blg =',LitrfallElms_pft(NE,NZ),LitrfallAbvgElms_pft(NE,NZ),LitrfallBlgrElms_pft(NE,NZ)
    write(888,*)'disturbP         =',PlantElmDistLoss_pft(NE,NZ)    
    write(888,*)'FireLossP         =',FireLossE_pft(NE,NZ)    
    write(888,*)'RootPuptk        =',RootNutUptakeP_pft(NZ)      
    write(888,*)'seed planted     =',plt_biom%SeedPlantedElm_pft(NE,NZ)  
    write(888,*)'surf literfall P =',SurfLitrfallElms_pft(NE,NZ)      
    if(I/=plt_distb%iDayPlantHarvest_pft(NZ) .and. I/=plt_distb%iDayPlanting_pft(NZ)) &  
      call endrun('P balance error test failure in '//trim(mod_filename)//' at line',__LINE__)      
  endif
  
  DO NE=1,NumPlantChemElms
    PlantElmBalCum_pft(NE,NZ)=PlantElmBalCum_pft(NE,NZ)+balE(NE)
  ENDDO
  call PrintInfo('end '//subname)  
  end associate

  end subroutine CheckPlantBalanceZ
  ![tail]
end module PlantBalMod
