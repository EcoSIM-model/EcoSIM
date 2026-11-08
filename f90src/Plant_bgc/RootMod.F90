module RootMod
  use data_kind_mod,       only: r8 => DAT_KIND_R8, yearIJ_type
  use minimathmod,         only: safe_adb,  AZMAX1, AZMIN1, AZERO, sfexp
  use EcoSIMCtrlMod,       only: lsoilCompaction,llignification,lcoarseroot
  use DebugToolMod,        only: PrintInfo
  use abortutils,          only: endrun
  use PlantNonstElmDynMod, only: RepleteLowSeaStorByRoot
  use PlantBalMod,         only: SumRootBiome,  SumRootAR, SumLitfallBlg
  use EcosimConst
  use PlantBGCPars
  use ElmIDMod
  use PlantMathFuncMod
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantPhotosynthesisAPIData, only : plt_photo
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantSoilChemistryAPIData, only : plt_soilchem
  use PlantAllometryAPIData, only : plt_allom
  use PlantBiomassAPIData, only : plt_biom
  use PlantEnergyWaterAPIData, only : plt_ew
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantRootBGCAPIData, only : plt_rbgc
  use NoduleBGCMod
implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  public :: RootBGCModel

  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine RootBGCModel(yearIJ,NZ,TFN6_vr,CNRTW,CPRTW,RootSinkC_vr,RootSinkC)
!
! Roots and mycorrhizae both have primary and secondary components
! primary roots are accordant to shoots, seminar root corresponds to main branch.
! not all plants have seminal roots.
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer , intent(in) :: NZ
  real(r8), intent(in) :: TFN6_vr(JZ1)                        !temperature function for root maintenance respiration
  real(r8), intent(in) :: CNRTW,CPRTW                         !NC, PC mass ratio of root growth
  real(r8), intent(out) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(out) :: RootSinkC(pltpar%jroots)

  real(r8) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  real(r8) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8) :: Root1stSinkTip(pltpar%MaxNumRootAxes)
  real(r8) :: RootMSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  real(r8) :: fUSE4MR_pvr(JZ1,pltpar%MaxNumRootAxes)
  character(len=*), parameter :: subname='RootBGCModel'
  integer, parameter  :: NumRootAxes4DeadPlant =0    !
  real(r8) :: TotAbsorbRootVol
  real(r8) :: litrflx(NumPlantChemElms)
  real(r8) :: RCO2flx,err,resp0
  integer :: N,M,NE,L,K
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: tmpval
  associate(                                                            &
    SeasonalNonstElms_pft      => plt_biom%SeasonalNonstElms_pft       ,& !input  :plant stored nonstructural element at current step, [g d-2]
    ZERO4LeafVar_pft           => plt_biom%ZERO4LeafVar_pft            ,& !input  :threshold zero for leaf calculation, [-]
    DLYR3                      => plt_site%DLYR3                       ,& !input  :vertical thickness of soil layer, [m]
    PlantPopuLive_pft          => plt_site%PlantPopuLive_pft           ,& !input  :plant population, [d-2]
    ZERO                       => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]
    NU                         => plt_site%NU                          ,& !input  :current soil surface layer number, [-]    
    MaxSoilLays4Root_pft       => plt_morph%MaxSoilLays4Root_pft       ,& !input  :maximum soil layer number for all root axes,[-]        
    Myco_pft                   => plt_morph%Myco_pft                   ,& !input  :mycorrhizal type (no or yes),[-]    
    iPlantPhenolPattern_pft    => plt_pheno%iPlantPhenolPattern_pft    ,& !input  :plant growth habit: annual or perennial,[-]
    LitrfallElms_pvr           => plt_bgcr%LitrfallElms_pvr            ,& !input  :plant LitrFall element, [g d-2 h-1]        
    RootNutUptakeN_pft         => plt_rbgc%RootNutUptakeN_pft          ,& !inoput :total N uptake by plant roots, [gN d-h2 h-1]
    RootNutUptakeP_pft         => plt_rbgc%RootNutUptakeP_pft          ,& !inoput :total P uptake by plant roots, [gP d-h2 h-1]
    SeedMeanLen_pft            => plt_morph%SeedMeanLen_pft            ,& !input  :seed length, [m]
    RootVH2O_pvr               => plt_morph%RootVH2O_pvr               ,& !input  :root layer volume water, [m2 d-2]
    RootCO2Autor_pvr           => plt_rbgc%RootCO2Autor_pvr            ,& !input  :root respiration constrained by O2, [g d-2 h-1]   
    NMaxRootBotLayer_pft       => plt_morph%NMaxRootBotLayer_pft       ,& !input  :maximum soil layer number for all root axes, [-]         
    RootSAreaPerPlant_pvr      => plt_morph%RootSAreaPerPlant_pvr      ,& !input  :root layer area per plant, [m p-1]
    RootPoreVol_pvr            => plt_morph%RootPoreVol_pvr            ,& !input  :root layer volume air, [m2 d-2]
    RootLenDensPerPlant_pvr    => plt_morph%RootLenDensPerPlant_pvr    ,& !input  :root layer length density, [m m-3]
    RootTotLenPerPlant_pvr     => plt_morph%RootTotLenPerPlant_pvr     ,& !inoput :root layer length per plant, [m p-1]
    RootAbsorbLenPerPlant_pvr  => plt_morph%RootAbsorbLenPerPlant_pvr  ,& !inoput :root layer length per plant, [m p-1]    
    NGTopRootLayer_pft         => plt_morph%NGTopRootLayer_pft         ,& !input  :soil layer at planting depth, [-]
    RootPorosity_pft           => plt_morph%RootPorosity_pft           ,& !input  :root porosity, [m3 m-3]
    SeedVolumeMean_pft         => plt_morph%SeedVolumeMean_pft         ,& !input  :seed volume, [m3 ]
    SeedAreaMean_pft           => plt_morph%SeedAreaMean_pft           ,& !input  :seed surface area, [m2]
    NumStructuralRootAxes_pft       => plt_morph%NumStructuralRootAxes_pft       ,& !input  :number of structural root axes,[-]
    isPlantRootAlive_pft       => plt_pheno%isPlantRootAlive_pft       ,& !output :flag to detect root system death,[-]
    isPlantShootAlive_pft      => plt_pheno%isPlantShootAlive_pft       & !output :flag to detect canopy death,[-]
  )
!     ROOT GROWTH
!
  call PrintInfo('beg '//subname)
  !  call RootCheck(I,J,NZ,'head')
!  call SumRootBiome(yearIJ,NZ,mass_inital)
  !first add newly acquired nutrients
  call SummarizeRootSink(yearIJ,NZ,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,Root1stSinkTip,&
    RootMSink_pvr,RootSinkC,fUSE4MR_pvr)

  call SumRootBiome(yearIJ,NZ,mass_inital)
!  call SumRootBiome(yearIJ,NZ,mass_finale)
  call SumRootAR(NZ);call SumLitfallBlg(NZ)
  tmpval=plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ)

!  write(424,*)NZ,'root1',yearIJ%I*1000+yearIJ%J/24.,mass_finale(ielmn)-mass_inital(ielmn)- &
!    plt_rbgc%RootNutUptakeN_pft(NZ)+plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ),'afnod',mass_finale(ielmn),mass_inital(ielmn),&
!    plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ)

!  write(423,*)NZ,'root',yearIJ%I*1000+yearIJ%J/24.,mass_finale(ielmc)-mass_inital(ielmc)- &
!     plt_bgcr%RootAutoCO2_pft(NZ)+plt_bgcr%LitrfallBlgrElms_pft(ielmc,NZ),'afrsinkc',mass_finale(ielmc),mass_inital(ielmc),&
!     plt_bgcr%RootAutoCO2_pft(NZ),plt_bgcr%LitrfallBlgrElms_pft(ielmc,NZ)

  call RootBiochemistry(yearIJ,NZ,TFN6_vr,CNRTW,CPRTW,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,Root1stSinkTip,&
    RootMSink_pvr,RootSinkC,fUSE4MR_pvr)

  call SumRootBiome(yearIJ,NZ,mass_finale) 
  call SumRootAR(NZ);call SumLitfallBlg(NZ)
!  write(423,*)NZ,'root',yearIJ%I*1000+yearIJ%J/24.,mass_finale(ielmc)-mass_inital(ielmc)- &
!     plt_bgcr%RootAutoCO2_pft(NZ)+plt_bgcr%LitrfallBlgrElms_pft(ielmc,NZ),'afrootbgc',mass_finale(ielmc),mass_inital(ielmc),&
!     plt_bgcr%RootAutoCO2_pft(NZ),plt_bgcr%LitrfallBlgrElms_pft(ielmc,NZ)
  if(.false.)write(426,*)NZ,'root2',yearIJ%I*1000+yearIJ%J/24.,mass_finale(ielmn)-mass_inital(ielmn)+ &
    plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ)-tmpval,'afnod',mass_finale(ielmn),mass_inital(ielmn),plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ),tmpval

  !
  !     ADD SEED DIMENSIONS TO ROOT DIMENSIONS (ONLY IMPORTANT DURING
  !     GERMINATION)
  !  
  RootTotLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)    = RootTotLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+SeedMeanLen_pft(NZ)
  RootAbsorbLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ) = RootAbsorbLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+SeedMeanLen_pft(NZ)
  
  IF(DLYR3(NGTopRootLayer_pft(NZ)).GT.ZERO)THEN
    !by m3
    RootLenDensPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)= &
      RootTotLenPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)/DLYR3(NGTopRootLayer_pft(NZ))
  ELSE
    RootLenDensPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=0._r8
  ENDIF

  TotAbsorbRootVol=RootPoreVol_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+RootVH2O_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+ &
    SeedVolumeMean_pft(NZ)*PlantPopuLive_pft(NZ)
  RootPoreVol_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)       = RootPorosity_pft(ipltroot,NZ)*TotAbsorbRootVol
  RootVH2O_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)          = (1.0_r8-RootPorosity_pft(ipltroot,NZ))*TotAbsorbRootVol
  RootSAreaPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ) = RootSAreaPerPlant_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+&
    SeedAreaMean_pft(NZ)

  IF(NumRootAxes4DeadPlant.EQ.NumStructuralRootAxes_pft(NZ) .OR. (SeasonalNonstElms_pft(ielmc,NZ).LE.ZERO4LeafVar_pft(NZ).AND. &
    iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial))THEN
    isPlantRootAlive_pft(NZ)  = iFalse
    isPlantShootAlive_pft(NZ) = iFalse
  ENDIF
  call SumRootBiome(yearIJ,NZ,mass_inital)
!
!     ROOT N2 FIXATION (RHIZOBIA)
  call RootNodulBiochemistry(yearIJ%I,yearIJ%J,NZ,TFN6_vr)

!  call SumRootBiome(yearIJ,NZ,mass_finale)
!  call SumRootAR(NZ);call SumLitfallBlg(NZ)
!  write(423,*)NZ,'root',yearIJ%I*1000+yearIJ%J/24.,mass_finale(ielmc)-mass_inital(ielmc)- &
!    plt_bgcr%RootAutoCO2_pft(NZ)+plt_bgcr%LitrfallBlgrElms_pft(ielmc,NZ),'afnod',mass_finale(ielmc),mass_inital(ielmc)
!  write(424,*)NZ,'root3',yearIJ%I*1000+yearIJ%J/24.,mass_finale(ielmn)-mass_inital(ielmn)- &
!    plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ),'afnod',mass_finale(ielmn),mass_inital(ielmn),plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ)

  call PrintInfo('end '//subname)  
  end associate
  end subroutine RootBGCModel

!----------------------------------------------------------------------------------------------------
  subroutine RootBiochemistry(yearIJ,NZ,TFN6_vr,CNRTW,CPRTW,RootSinkC_vr,Root1stSink_pvr,&
    Root2ndSink_pvr,Root1stSinkTip,RootMSink_pvr,RootSinkC,fUSE4MR_pvr)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), intent(in) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: CNRTW,CPRTW    !NC, PC mass ratio of root growth
  REAL(R8), INTENT(in) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), INTENT(in) :: RootSinkC(pltpar%jroots)
  real(r8), INTENT(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root1stSinkTip(pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: RootMSink_pvr(JZ1,pltpar%MaxNumRootAxes)   
  real(r8), intent(in) :: fUSE4MR_pvr(JZ1,pltpar%MaxNumRootAxes)  
  integer :: LL,LZ,Lnext,L,N
  logical  :: Root1stTipUpdateFlag(pltpar%MaxNumRootAxes)  !Primary root tip zone update flag
  logical  :: FoundRootAxesTip(pltpar%jroots,pltpar%MaxNumRootAxes)
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: masst_inital(NumPlantChemElms)
  real(r8) :: masst_finale(NumPlantChemElms)  
  real(r8) :: litrflxt(NumPlantChemElms)
  real(r8) :: RCO2flxt,tmpval,dmassN
  character(len=*), parameter :: subname='RootBiochemistry'

!     begin_execution
  associate(                                                                 &
    iPlantRootProfile_pft        => plt_pheno%iPlantRootProfile_pft         ,& !input  :plant growth type (vascular, non-vascular),[-]
    PSIRoot_pvr                  => plt_ew%PSIRoot_pvr                      ,& !input  :root total water potential, [Mpa]
    PSIRootTurg_vr               => plt_ew%PSIRootTurg_vr                   ,& !input  :root turgor water potential, [Mpa]
    VLSoilPoreMicP_vr            => plt_soilchem%VLSoilPoreMicP_vr          ,& !input  :volume of soil layer, [m3 d-2]
    NU                           => plt_site%NU                             ,& !input  :current soil surface layer number, [-]
    MaxNumRootLays               => plt_site%MaxNumRootLays                 ,& !input  :maximum root layer number,[-]
    ZEROS2                       => plt_site%ZEROS2                         ,& !input  :threshold zero for numerical stability,[-]
    NL                           => plt_site%NL                             ,& !input  :lowest soil layer number,[-]
    flag2ndGrowth_pvr            => plt_morph%flag2ndGrowth_pvr             ,& !input  :flag for secondary growth of primary roots, [-]    
    Myco_pft                     => plt_morph%Myco_pft                      ,& !input  :mycorrhizal type (no or yes),[-]
    NGTopRootLayer_pft           => plt_morph%NGTopRootLayer_pft            ,& !input  :soil layer at planting depth, [-]
    fRootGrowPSISense_pvr        => plt_pheno%fRootGrowPSISense_pvr         ,& !output :water stress to plant root growth, [-]
    NMaxRootBotLayer_pft         => plt_morph%NMaxRootBotLayer_pft           & !output :maximum soil layer number for all root axes, [-]
  )
  call PrintInfo('beg '//subname)

  Root1stTipUpdateFlag     = .false. !initialize all primary roots as not updated
  FoundRootAxesTip        = .false.
  NMaxRootBotLayer_pft(NZ) = NGTopRootLayer_pft(NZ) !initialize diagnostic for the deepest root layer
  !
  !     RESPIRATION AND GROWTH OF ROOT, MYCORRHIZAE IN EACH LAYER
  !
  litrflxt=0._r8;RCO2flxt=0._r8
  CALL GetPlantRoot1stDepz(yearIJ,NZ)

!  call SumRootBiome(yearIJ,NZ,masst_inital)
!  call SumRootAR(NZ);call SumLitfallBlg(NZ)
!  tmpval=plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ)
  plt_morph%RootMediumVH2O_rpvr(:,:,NZ)  = 0._r8
  ! Production is a current-hour rate; skipped growth must not reuse the previous rate.
  plt_rbgc%RootMRProdCytok_rpvr(:,:,NZ) = 0._r8
  D5010: DO N=1,Myco_pft(NZ)
    ! Growth and reserve transfers must use the same layers as SumRootBiome.
    ! The additional lower soil layer is not part of the root mass balance.
    D5000: DO L=NU,MaxNumRootLays

      call SumRootBiome(yearIJ,NZ,masst_inital)
      call SumRootAR(NZ);call SumLitfallBlg(NZ)
      tmpval=plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ)


      plt_morph%Root2ndVH2O_rpvr(N,L,:,NZ) = 0._r8      
      !
      !     IDENTIFY NEXT LOWER ROOT LAYER
      !
      !     VLSoilPoreMicP_vr=soil layer volume excluding macropore, rocks
      !
      IF(VLSoilPoreMicP_vr(L).GT.ZEROS2)THEN
        Lnext=L
        D5003: DO LZ=L+1,MaxNumRootLays
          IF(VLSoilPoreMicP_vr(LZ).GT.ZEROS2 .OR. LZ.EQ.NL)THEN
            Lnext=LZ
            EXIT
          ENDIF
        ENDDO D5003

        IF(is_root_bryophyte(iPlantRootProfile_pft(NZ)))THEN
          fRootGrowPSISense_pvr(N,L,NZ)=EXP(0.05_r8*AMAX1(PSIRoot_pvr(N,L,NZ),-5000._r8))
        ELSE
          fRootGrowPSISense_pvr(N,L,NZ)=EXP(0.10_r8*AMAX1(PSIRoot_pvr(N,L,NZ),-5000._r8))
        ENDIF  
        
        call GrowRootMycoAxes(yearIJ,N,L,Lnext,NZ,FoundRootAxesTip,Root1stTipUpdateFlag,TFN6_vr,&
          RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,Root1stSinkTip,RootMSink_pvr,fUSE4MR_pvr, CNRTW,CPRTW,&
          fRootGrowPSISense_pvr(N,L,NZ),litrflxt,RCO2flxt)
        
         call SumRootBiome(yearIJ,NZ,masst_finale)
         call SumRootAR(NZ);call SumLitfallBlg(NZ)
!        if(yearIJ%I>=143)write(431,*)yearIJ%I*1000+yearIJ%J/24.,masst_finale(ielmc)-masst_inital(ielmc)- &
!          plt_bgcr%RootAutoCO2_pft(NZ)+plt_bgcr%LitrfallBlgrElms_pft(ielmc,NZ)-tmpval,'rbgc',N,L
  dmassN=masst_finale(ielmn)-masst_inital(ielmn)+ (plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ)-tmpval)
  if(.false.)write(426,*)'root',yearIJ%I*1000+yearIJ%J/24.,N*10+L,dmassN,'afnod',masst_finale(ielmn),masst_inital(ielmn),plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ)
!  if(abs(dmassN)>1.e-8)stop

!====================================================================================
        !masst_inital=masst_finale
        !     DRAW FROM ROOT NON-STRUCTURAL POOL WHEN
        !     SEASONAL STORAGE POOL IS DEPLETED
        call RepleteLowSeaStorByRoot(yearIJ%I,yearIJ%J,N,L,NZ)
        call SumRootBiome(yearIJ,NZ,masst_finale)
        if(.false.)write(426,*)'root',yearIJ%I*1000+yearIJ%J/24.,N*10+L,dmassN,'repsea',masst_finale(ielmn)
        call DiagRootGeometry(yearIJ%I,yearIJ%J,N,L,NZ)
      !   call SumRootBiome(NZ,masst_finale)

      ENDIF

    ENDDO D5000
  ENDDO D5010
  !write(1114,*)yearIJ%I*1000+yearIJ%J/24.,1,2,0.0,plt_biom%RootMycoNonstElms_rpvr(ielmc,1,2,NZ),'af'//subname

  !Keep the per-layer diagnostics used during growth, then refresh the final state.
  !Withdrawal in layer L can also remove roots/mycorrhizae in the already processed
  !layer L-1. All categories and layers must be consistent before cytokinin transport.
  DO N=1,Myco_pft(NZ)
    DO L=NU,MaxNumRootLays
      IF(VLSoilPoreMicP_vr(L).GT.ZEROS2)THEN
        call DiagRootGeometry(yearIJ%I,yearIJ%J,N,L,NZ)
      ENDIF
    ENDDO
  ENDDO

  call plt_morph%RefreshMediumRootMeanLength(NZ)
  call CytoKininDynamics(yearIJ%I,yearIJ%J,NZ,Root2ndSink_pvr)

  call PrintInfo('end '//subname)
  end associate
  end subroutine RootBiochemistry
!----------------------------------------------------------------------------------------------------

  subroutine DiagRootGeometry(I,J,N,L,NZ)
  !
  !Description:
  !Diagnose geometry properties of roots
  !
  implicit none
  integer , intent(in)  :: I,J,N,L,NZ
  real(r8) :: TotRoot2ndPopuC             !current secondary root C mass in layer L
  real(r8) :: TotPopuRoot2ndLlenAxes      !current population secondary root length
  real(r8) :: CoarseVol                    !coarse root volume, [m3]    
  real(r8) :: TotRoot1stPopuC              !total primary root C mass
  real(r8) :: TotRoot1stLlenAxesPP         !total primary root length  
  real(r8) :: TotPopuRoot1stLen 
  real(r8) :: TotPopuRootLen               !total root length
  real(r8) :: TotFineRootC                 !total root C mass
  real(r8) :: TotAbsorbRootVol             !total volume
  real(r8) :: Root1stSurfAbsorbArea,Root2ndSurfArea
  real(r8) :: Root1stRadiusEst,SCAL,FineRoot1stTotc,TotFineRoot1stLlenAxesPP  
  real(r8) :: TotCoarseRoot1stLlenAxesPP,CoarseRoot1stTotC
  real(r8) :: MediumTransptAxisCount
  real(r8) :: AreaTranspt,DTransptTube,TotAbsorbRootVolPP       
  integer  :: NTG,NE,NR
  logical :: checkRootExt,check2ndGrothRoot
  character(len=*), parameter :: subname='DiagRootGeometry'
  associate(                                                                      &  
    ZERO4Groth_pft                 => plt_biom%ZERO4Groth_pft                    ,& !input  :threshold zero for plang growth calculation, [-]
    PlantPopuLive_pft              => plt_site%PlantPopuLive_pft                 ,& !input  :plant population, [d-2]
    ZERO                           => plt_site%ZERO                              ,& !input  :threshold zero for numerical stability, [-]
    RootPorosity_pft               => plt_morph%RootPorosity_pft                 ,& !input  :root porosity, [m3 m-3]
    DLYR3                          => plt_site%DLYR3                             ,& !input  :vertical thickness of soil layer, [m]
    PSIRoot_pvr                    => plt_ew%PSIRoot_pvr                         ,& !input  :root total water potential, [Mpa]
    PSIRootTurg_vr                 => plt_ew%PSIRootTurg_vr                      ,& !input  :root turgor water potential, [Mpa]
    Root1stMaxRadius1_pft          => plt_morph%Root1stMaxRadius1_pft            ,& !input  :root radius primary axes, [m]
    StalkAveRadius_pft             => plt_morph%StalkAveRadius_pft               ,& !input  :main stalk radius,[m]
    Root2ndXSecArea_pft            => plt_morph%Root2ndXSecArea_pft              ,& !input  :root cross-sectional area of secondary axes, [m2]
    Root2ndMaxRadius1_pft          => plt_morph%Root2ndMaxRadius1_pft            ,& !input  :root radius secondary axes, [m]
    Root2ndMaxRadius_pft           => plt_morph%Root2ndMaxRadius_pft             ,& !input  :maximum radius of secondary roots, [m]
    Root1stMaxRadius_pft           => plt_morph%Root1stMaxRadius_pft             ,& !input  :maximum radius of primary roots, [m]
    k_fine_comp                    => pltpar%k_fine_comp                         ,& !input  :fine litter complex id
    Root1stXSecArea_pft            => plt_morph%Root1stXSecArea_pft              ,& !input  :root cross-sectional area primary axes, [m2]
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft            ,& !input  :plant growth type (vascular, non-vascular),[-]
    RootMatureAge_pft              => plt_morph%RootMatureAge_pft                ,& !input : Root maturation age, [h]
    NumStructuralRootAxes_pft      => plt_morph%NumStructuralRootAxes_pft        ,& !input  :number of structural root axes,[-]
    FracRootElmAllocm              => plt_allom%FracRootElmAllocm                ,& !input  :C woody/soft fraction in root,[-]
    Root1stLenPP_rpvr              => plt_morph%Root1stLenPP_rpvr                ,& !input  :length per primary root axis in soil layer, [m]
    Root2ndXNum_rpvr               => plt_morph%Root2ndXNum_rpvr                 ,& !inoput :secondary root axes per structural axis, [d-2]
    Root2ndXNumL_rpvr              => plt_morph%Root2ndXNumL_rpvr                ,& !output :root layer number axes, [d-2]
    RootMediumXNum_rpvr            => plt_morph%RootMediumXNum_rpvr              ,& !inout :living population medium-root axes per structural group, [d-2]
    RootMediumXNum_pvr             => plt_morph%RootMediumXNum_pvr               ,& !output :surviving population medium-root axes in layer, [-]
    RootMediumLength_rpvr          => plt_morph%RootMediumLength_rpvr            ,& !input :population medium-root length per structural group, [m d-2]
    RootMediumRadius_rpvr          => plt_morph%RootMediumRadius_rpvr            ,& !inoput :root layer radius for medium size axes, [m]    
    NGTopRootLayer_pft             => plt_morph%NGTopRootLayer_pft               ,& !input  :soil layer at planting depth, [-]
    Num1stAxesPerStructRootX_pft           => plt_morph%Num1stAxesPerStructRootX_pft             ,& !input :primary axes per structural-root group per plant, [-]
    Root2ndLen_rpvr                => plt_morph%Root2ndLen_rpvr                  ,& !input :population secondary root length, [m d-2]
    RootMyco2ndStrutElms_rpvr      => plt_biom%RootMyco2ndStrutElms_rpvr         ,& !input :secondary root structural elements, [g d-2]
    Root2ndVH2O_rpvr               => plt_morph%Root2ndVH2O_rpvr                 ,& !output :water-occupied secondary root volume
    RootMyco1stStrutElms_rpvr      => plt_biom%RootMyco1stStrutElms_rpvr         ,& !input  :root layer element primary axes, [g d-2]
    RootMediumStructElms_rpvr      => plt_biom%RootMediumStructElms_rpvr         ,& !inoput :root layer element for medium size root axes, [g d-2]    
    MaxSoilLays4Root_pft           => plt_morph%MaxSoilLays4Root_pft             ,& !input  :maximum soil layer number for all root axes,[-]
    flag2ndGrowth_pvr              => plt_morph%flag2ndGrowth_pvr                ,& !input  :flag for secondary growth of primary roots, [-]
    iPlant2ndGrothPattern_pft      => plt_pheno%iPlant2ndGrothPattern_pft        ,& !input  :plant expression of secondary growth, [-]
    CRootActVolPerMassC_pft        => plt_morph%CRootActVolPerMassC_pft          ,& !input  :coarse root volume:mass ratio, [m3 gC-1]
    FineRootVolPerMassC_pft        => plt_morph%FineRootVolPerMassC_pft          ,& !input  :fine root volume:mass ratio, [m3 gC-1]
    Root1stRadius_rpvr             => plt_morph%Root1stRadius_rpvr               ,& !input  :root layer radius for each primary axes, [m]
    Root1stTransptArea_pvr         => plt_morph%Root1stTransptArea_pvr           ,& !output :mean transport area by 1st order root, [m2 d-2]
    RootMedTransptArea_pvr         => plt_morph%RootMedTransptArea_pvr           ,& !output  :transport area by medisum size roots, [-]    
    Root1stRadius_pvr              => plt_morph%Root1stRadius_pvr                ,& !output :root layer radius primary axes, [m]
    RootVH2O_pvr                   => plt_morph%RootVH2O_pvr                     ,& !output :root space volume occupied by water in each layer, [m2 d-2]
    Root2ndRadius_rpvr             => plt_morph%Root2ndRadius_rpvr               ,& !output :root layer radius secondary axes, [m]
    RootLenDensPerPlant_pvr        => plt_morph%RootLenDensPerPlant_pvr          ,& !output :layer root length density, [m m-3]
    RootArea1stPP_pvr              => plt_morph%RootArea1stPP_pvr                ,& !output :layer 1st root area per plant, [m2 plant-1]
    RootArea2ndPP_pvr              => plt_morph%RootArea2ndPP_pvr                ,& !output :layer 2nd root area per plant, [m2 plant-1]
    RootPoreVol_pvr                => plt_morph%RootPoreVol_pvr                  ,& !output :layer root volume air, [m2 d-2]
    trcg_rootml_pvr                => plt_rbgc%trcg_rootml_pvr                   ,& !inoput :root gas content, [g d-2]
    trcs_rootml_pvr                => plt_rbgc%trcs_rootml_pvr                   ,& !inoput :root aqueous content, [g d-2]
    RootSAreaPerPlant_pvr          => plt_morph%RootSAreaPerPlant_pvr            ,& !output :layer 1st+2nd root area per plant, [m2 plant-1]
    RootTotLenPerPlant_pvr         => plt_morph%RootTotLenPerPlant_pvr           ,& !output :total root length per plant in layer (including root hair), [m p-1]
    RootAbsorbLenPerPlant_pvr      => plt_morph%RootAbsorbLenPerPlant_pvr        ,& !output :total absorptive root length per plant in layer, [m p-1]
    Root2ndEffLen4uptk_rpvr        => plt_morph%Root2ndEffLen4uptk_rpvr          ,& !output :root layer average length, [m]
    RootLenPerPlant_pvr            => plt_morph%RootLenPerPlant_pvr              ,& !output :root length per plant in layer (excluding root hair), [m p-1]
    RootGasLossDisturb_pft         => plt_bgcr%RootGasLossDisturb_pft             & !inoput :gaseous flux fron root disturbance, [g d-2 h-1]
  )
  call PrintInfo('beg '//subname)
  !Derive fine-root geometry from surviving state, including after withdrawal.
  TotRoot2ndPopuC = 0._r8
  TotPopuRoot2ndLlenAxes = 0._r8
  Root2ndVH2O_rpvr(N,L,:,NZ) = 0._r8
  DO NR=1,NumStructuralRootAxes_pft(NZ)
    ! Mortality may shorten roots after their counts were set during growth.
    ! Retain surviving tips only; this refresh must not create new axes.
    IF(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ).LE.ZERO4Groth_pft(NZ) .OR. &
       Root2ndLen_rpvr(N,L,NR,NZ).LE.0._r8)THEN
      Root2ndXNum_rpvr(N,L,NR,NZ) = 0._r8
    ELSE
      Root2ndXNum_rpvr(N,L,NR,NZ) = MIN(AZMAX1(Root2ndXNum_rpvr(N,L,NR,NZ)), &
        Root2ndLen_rpvr(N,L,NR,NZ)/Root2ndTipLen4uptk)
    ENDIF
    TotRoot2ndPopuC = TotRoot2ndPopuC+RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)
    TotPopuRoot2ndLlenAxes = TotPopuRoot2ndLlenAxes+Root2ndLen_rpvr(N,L,NR,NZ)
    Root2ndVH2O_rpvr(N,L,NR,NZ) = AZERO(Root2ndXSecArea_pft(N,NZ)*Root2ndLen_rpvr(N,L,NR,NZ) &
      *(1._r8-RootPorosity_pft(N,NZ)))
  ENDDO
  ! Sum the corrected counts before effective lengths and cytokinin water routing.
  Root2ndXNumL_rpvr(N,L,NZ) = SUM(Root2ndXNum_rpvr(N,L,1:NumStructuralRootAxes_pft(NZ),NZ))
  !
  CoarseRoot1stTotC              = 0._r8
  FineRoot1stTotc                = 0._r8
  TotFineRoot1stLlenAxesPP       = 0._r8
  TotCoarseRoot1stLlenAxesPP     = 0._r8
  TotAbsorbRootVolPP             = 0._r8
  Root1stSurfAbsorbArea          = 0._r8
  Root1stTransptArea_pvr(N,L,NZ) = 0._r8
  RootMedTransptArea_pvr(N,L,NZ) = 0._R8
  MediumTransptAxisCount = 0._r8

  DO NR=1,NumStructuralRootAxes_pft(NZ)
    ! Each NR contributes its own primary-axis count; the NR loop sums groups.
    ! TotRoot1stLlenAxesPP = total primary root length per plant
    ! TotRoot1stPopuC      = total primary root C mass for whole population
    if(N.EQ.ipltroot)then
      DTransptTube     = AMIN1(ZSTX,AMAX1(FSTK*Root1stRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(N,NZ)))
      Root1stRadiusEst = AMAX1(Root1stRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ))
      !
      !transport area for axis NR, 
      AreaTranspt   = PICON*(2._r8*Root1stRadiusEst*DTransptTube-DTransptTube**2)
      !
      !active secondary growth
      if(flag2ndGrowth_pvr(L,NR,NZ) .and. Root1stLenPP_rpvr(L,NR,NZ).GT.0._r8)then
        CoarseRoot1stTotC          = CoarseRoot1stTotC+RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)
        TotCoarseRoot1stLlenAxesPP = TotCoarseRoot1stLlenAxesPP+Root1stLenPP_rpvr(L,NR,NZ)*Num1stAxesPerStructRootX_pft(NZ)
      else
        FineRoot1stTotc          = FineRoot1stTotc+RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)
        TotFineRoot1stLlenAxesPP = TotFineRoot1stLlenAxesPP+Root1stLenPP_rpvr(L,NR,NZ)*Num1stAxesPerStructRootX_pft(NZ)
        !assume cylindrical
        TotAbsorbRootVolPP       = TotAbsorbRootVolPP+AreaTranspt*Root1stLenPP_rpvr(L,NR,NZ)*Num1stAxesPerStructRootX_pft(NZ)
      endif      
      Root1stTransptArea_pvr(N,L,NZ)=Root1stTransptArea_pvr(N,L,NZ)+AreaTranspt
      !Counts describe living axes; establishment capacity is local to the sink calculation.
      if(RootMediumStructElms_rpvr(ielmc,L,NR,NZ).LE.0._r8 .or. &
         RootMediumLength_rpvr(L,NR,NZ).LE.0._r8)RootMediumXNum_rpvr(L,NR,NZ)=0._r8
      !Only existing tissue contributes to the layer count and transport area.
      if(RootMediumStructElms_rpvr(ielmc,L,NR,NZ).GT.0._r8 .and. &
         RootMediumLength_rpvr(L,NR,NZ).GT.0._r8 .and. &
         RootMediumXNum_rpvr(L,NR,NZ).GT.0._r8)then
        DTransptTube                   = AMIN1(ZSTX,AMAX1(FSTK*RootMediumRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(ipltroot,NZ)))
        AreaTranspt                    = PICON*(2._r8*RootMediumRadius_rpvr(L,NR,NZ)*DTransptTube-DTransptTube**2)
        ! Average over medium-root axes, not over their parent structural groups.
        RootMedTransptArea_pvr(N,L,NZ) = RootMedTransptArea_pvr(N,L,NZ) &
          +RootMediumXNum_rpvr(L,NR,NZ)*AreaTranspt
        MediumTransptAxisCount = MediumTransptAxisCount+RootMediumXNum_rpvr(L,NR,NZ)
      endif  
    endif              
  ENDDO
  if(NumStructuralRootAxes_pft(NZ).GT.0)then 
    Root1stTransptArea_pvr(N,L,NZ)=Root1stTransptArea_pvr(N,L,NZ)/NumStructuralRootAxes_pft(NZ)
  endif

  ! Use the same surviving axes for layer count and mean transport area.
  ! The mycorrhizal diagnostic must not overwrite the plant medium-root count.
  if(N.EQ.ipltroot)RootMediumXNum_pvr(L,NZ) = MediumTransptAxisCount
  if(MediumTransptAxisCount.GT.0._r8)then
    RootMedTransptArea_pvr(N,L,NZ) = RootMedTransptArea_pvr(N,L,NZ)/MediumTransptAxisCount
  endif

  ! ROOT AND MYCORRHIZAL LENGTH, DENSITY, VOLUME, RADIUS, AREA
  ! TO CALCULATE WATER AND NUTRIENT UPTAKE IN 'UPTAKE'
  ! account for fine roots only

  IF(N.EQ.ipltroot)THEN
    TotRoot1stLlenAxesPP = TotCoarseRoot1stLlenAxesPP+TotFineRoot1stLlenAxesPP
    CoarseVol            = CoarseRoot1stTotC*CRootActVolPerMassC_pft(NZ)
  else
    TotRoot1stLlenAxesPP = TotCoarseRoot1stLlenAxesPP+TotFineRoot1stLlenAxesPP
    CoarseVol            = 0._r8
  ENDIF

  TotFineRootC      = FineRoot1stTotc+TotRoot2ndPopuC
  TotPopuRoot1stLen = TotRoot1stLlenAxesPP*PlantPopuLive_pft(NZ)
  TotPopuRootLen    = TotPopuRoot1stLen+TotPopuRoot2ndLlenAxes

  checkRootExt      = TotPopuRootLen.GT.ZERO4Groth_pft(NZ) .AND. TotFineRootC.GT.ZERO4Groth_pft(NZ) &
    .AND. PlantPopuLive_pft(NZ).GT.ZERO4Groth_pft(NZ)

  IF(checkRootExt)THEN  
    !there are roots
    RootAbsorbLenPerPlant_pvr(N,L,NZ) = TotFineRoot1stLlenAxesPP+TotPopuRoot2ndLlenAxes/PlantPopuLive_pft(NZ)
    RootTotLenPerPlant_pvr(N,L,NZ)    = TotPopuRootLen/PlantPopuLive_pft(NZ)
    RootLenPerPlant_pvr(N,L,NZ)       = TotRoot1stLlenAxesPP
    IF(DLYR3(L).GT.ZERO)THEN
      !per volume
      RootLenDensPerPlant_pvr(N,L,NZ) = RootTotLenPerPlant_pvr(N,L,NZ)/DLYR3(L)
    ELSE
      RootLenDensPerPlant_pvr(N,L,NZ) = 0._r8
    ENDIF
    !
    !total root volume [m3], with FineRootVolPerMassC_pft including contribution from water, fine roots volume > 80% as water
    !coarse roots, 45-60% volume as water 

    check2ndGrothRoot=N.EQ.ipltroot .and. is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)) &
      .and. lcoarseroot .and. (TotPopuRoot1stLen>0._r8 .and. CoarseVol > 0._r8 .OR. &
        (L.LT.NGTopRootLayer_pft(NZ) .and. MaxSoilLays4Root_pft(NZ).GT.NGTopRootLayer_pft(NZ)))
    
    if(check2ndGrothRoot)then
      !assuming elongation zone is negligible in length, pi*r^2*l=vol/pop       
      if(L.LT.NGTopRootLayer_pft(NZ))then
        Root1stRadiusEst = AZMAX1(StalkAveRadius_pft(NZ),Root1stMaxRadius1_pft(ipltroot,NZ))
      else
        Root1stRadiusEst = sqrt(CoarseVol/(PICON*TotPopuRoot1stLen))
      endif
      !total absorption root volume
      TotAbsorbRootVol = TotAbsorbRootVolPP*PlantPopuLive_pft(NZ) + AMAX1(Root2ndXSecArea_pft(N,NZ)*TotPopuRoot2ndLlenAxes, &
        TotRoot2ndPopuC*FineRootVolPerMassC_pft(N,NZ)*PSIRootTurg_vr(N,L,NZ))
      !primary root radius
      Root1stRadius_pvr(N,L,NZ) = AMAX1(Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ),Root1stRadiusEst)            
    else
      TotAbsorbRootVol=AMAX1(Root1stXSecArea_pft(N,NZ)*TotPopuRoot1stLen+Root2ndXSecArea_pft(N,NZ)*TotPopuRoot2ndLlenAxes &
        ,TotFineRootC*FineRootVolPerMassC_pft(N,NZ)*PSIRootTurg_vr(N,L,NZ))    

      Root1stRadius_pvr(N,L,NZ) = AMAX1(Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ))
      Root1stSurfAbsorbArea     = TwoPiCON*Root1stRadius_pvr(N,L,NZ)*TotPopuRoot1stLen      
    endif  
        
    TotAbsorbRootVol=AZMAX1(TotAbsorbRootVol)
    RootPoreVol_pvr(N,L,NZ) = RootPorosity_pft(N,NZ)*TotAbsorbRootVol
    RootVH2O_pvr(N,L,NZ)    = (1.0_r8-RootPorosity_pft(N,NZ))*TotAbsorbRootVol

    !primary and secondary root radius
    Root2ndRadius_rpvr(N,L,NZ) = AMAX1(Root2ndMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root2ndMaxRadius_pft(N,NZ))
    Root2ndSurfArea            = TwoPiCON*Root2ndRadius_rpvr(N,L,NZ)*TotPopuRoot2ndLlenAxes
    RootArea1stPP_pvr(N,L,NZ)  = Root1stSurfAbsorbArea/PlantPopuLive_pft(NZ)
    RootArea2ndPP_pvr(N,L,NZ)  = Root2ndSurfArea/PlantPopuLive_pft(NZ)

    IF(Root2ndXNumL_rpvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
      Root2ndEffLen4uptk_rpvr(N,L,NZ)=AMAX1(Root2ndTipLen4uptk,TotPopuRoot2ndLlenAxes/Root2ndXNumL_rpvr(N,L,NZ))
    ELSE
      Root2ndEffLen4uptk_rpvr(N,L,NZ)=Root2ndTipLen4uptk
    ENDIF
    RootSAreaPerPlant_pvr(N,L,NZ)=RootArea1stPP_pvr(N,L,NZ)+RootArea2ndPP_pvr(N,L,NZ)
  ELSE
    !no roots
    RootAbsorbLenPerPlant_pvr(N,L,NZ) = 0._r8
    RootArea1stPP_pvr(N,L,NZ)       = 0._r8
    RootArea2ndPP_pvr(N,L,NZ)       = 0._r8
    RootLenPerPlant_pvr(N,L,NZ)     = 0._r8
    RootTotLenPerPlant_pvr(N,L,NZ)  = 0._r8
    RootLenDensPerPlant_pvr(N,L,NZ) = 0._r8
    RootPoreVol_pvr(N,L,NZ)         = 0._r8
    RootVH2O_pvr(N,L,NZ)            = 0._r8
    RootSAreaPerPlant_pvr(N,L,NZ)   = 0._r8
    Root1stRadius_pvr(N,L,NZ)       = Root1stMaxRadius_pft(N,NZ)
    Root2ndRadius_rpvr(N,L,NZ)      = Root2ndMaxRadius_pft(N,NZ)
    Root1stTransptArea_pvr(N,L,NZ)  = PICON*Root1stMaxRadius_pft(N,NZ)**2
    Root2ndEffLen4uptk_rpvr(N,L,NZ) = Root2ndTipLen4uptk
    DO NTG=idg_beg,idg_NH3
      RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-(trcg_rootml_pvr(NTG,N,L,NZ)+trcs_rootml_pvr(NTG,N,L,NZ))
    ENDDO
    trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)=0._r8
    trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)=0._r8
  ENDIF

  call PrintInfo('end '//subname)      
  end associate      
  end subroutine DiagRootGeometry
!----------------------------------------------------------------------------------------------------
  subroutine Grow2ndRootAxes(I,J,L,NZ,N,NR,VMXCSpec, CNRTW,CPRTW,fdLext2nd,&
    RespElongWatSens,DMRespEff,fRootGrowPSISense,TFN6_vr,RootSinkC_vr,Root2ndSink_pvr,&
    litrflx2,RCO2flx2,Root2ndStrutRemob)
  !
  !Description:
  !Grow fine roots and mycorrhizae
  !  
  implicit none
  integer, intent(in) :: I,J,L,NZ
  integer, intent(in) :: N                          !root type id
  integer, intent(in) :: NR                         !root axis id
  real(r8), intent(in) :: VMXCSpec
  real(r8), intent(in) :: CNRTW,CPRTW  
  real(r8), intent(in) :: fdLext2nd                 !soil compaction modification of root extension
  real(r8), intent(in) :: RespElongWatSens          !respiration function of root water potential
  real(r8), intent(in) :: DMRespEff                 !root C respiration vs nonstructural C consumption, i.e. respiration quotient
  real(r8), intent(in) :: fRootGrowPSISense         !growth function of root water potential
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  REAL(R8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)  
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)  

  real(r8), intent(out) :: litrflx2(NumPlantChemElms)      !litterfall flux 
  real(r8), intent(out) :: RCO2flx2                        !CO2 respiraiton flux 
  real(r8), intent(out) :: Root2ndStrutRemob(NumPlantChemElms)      
  character(len=*), parameter :: subname='Grow2ndRootAxes'
  real(r8) :: Rmaint2nd_CO2                   !root maintenance respiration
  real(r8) :: RNonstCO2_OUltd,RNonstCO2_Oltd  !respiration from non-structural C unlimited/limited by O2
  REAL(R8) :: RCO2MaintDef2ndStruct_Oltd      !maintenance deficit met up by structural root biomass remobilization under O2 limitation
  real(r8) :: RCO2MaintDef2ndStruct_OUltd     !maintenance deficit met up by structural root biomass remobilization without O2 limitation
  real(r8) :: RCO2Nonst4Nassim_Oltd           !respiration for N assimilation ltd by O2
  real(r8) :: RCO2Nonst4Nassim_OUltd          !respiration for N assimilation unltd by O2
  real(r8) :: RootMycoNonst4Grow_Oltd(NumPlantChemElms)              !nonst biomass convereted into structural biomass growth
  real(r8) :: RootMycoNonst4GrowC_OUltd, RootMycoNonst4GrowC_Oltd    !Nonst C used for growth respiration + biomass increment unlimited/limited by O2
  real(r8) :: Root2ndPopExtenz  !whole plant population 2ndary root extension
  real(r8) :: RootMycoNonst4Grow_OUltd(NumPlantChemElms)  !nonstrucal biomass converted into biomass 
  real(r8) :: Root2ndNetGrowthElms(NumPlantChemElms)
  real(r8) :: dmass(NumPlantChemElms)  
  real(r8) :: RPotCO2Groth2nd_OUltd, RPotCO2Groth2nd_Oltd !potential growth resp  unltd,ltd by O2 
  real(r8) :: RCO2T2nd_OUltd, RCO2T2nd_Oltd       !total C respiration ltd, unltd by O2
  real(r8) :: RGrowCO2_OUltd,RGrowCO2_Oltd        !growth respiration limited by N,P unltd,ltd by O2
  real(r8) :: DMRTR
  real(r8) :: PPOOLB,ZPOOLB
  real(r8) :: Num2ndRoots,Num3rdRoots
  real(r8) :: Root2ndExtPot  
  real(r8) :: FNP,ratio
  real(r8) :: FracRoot2ndCSinkL,FracRoot2ndSinkN,FracRoot2ndSinkP  
  integer  :: M,NE
  real(r8) :: NumGroFineAxes,MaxNumGroFineAxes

  associate(                                                                      &
    ZERO4Groth_pft                 => plt_biom%ZERO4Groth_pft                    ,& !input  :threshold zero for plang growth calculation, [-]
    RAutoRootO2Limter_rpvr         => plt_rbgc%RAutoRootO2Limter_rpvr            ,& !input  :O2 constraint to root respiration (0-1), [-]
    FracRootElmAllocm              => plt_allom%FracRootElmAllocm                ,& !input  :C woody/soft fraction in root,[-]
    CNRTS_pft                      => plt_allom%CNRTS_pft                        ,& !input  :root N:C ratio x root growth yield, [-]
    CPRTS_pft                      => plt_allom%CPRTS_pft                        ,& !input  :root P:C ratio x root growth yield, [-]
    RootBiomGrosYld_pft            => plt_allom%RootBiomGrosYld_pft              ,& !input  :root growth yield, [g g-1]
    k_fine_comp                    => pltpar%k_fine_comp                         ,& !input  :fine litter complex id
    Num1stAxesPerStructRootXPOP_pft=> plt_morph%Num1stAxesPerStructRootXPOP_pft  ,& !output :population primary root axes number on one structrual axis, [d-2]        
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft            ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlantPhenolType_pft           => plt_pheno%iPlantPhenolType_pft             ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    fTgrowRootP_vr                 => plt_pheno%fTgrowRootP_vr                   ,& !input  :root layer temperature growth functiom, [-]
    PlantPopuLive_pft              => plt_site%PlantPopuLive_pft                 ,& !input  :plant population, [d-2]    
    DLYR3                          => plt_site%DLYR3                             ,& !input  :vertical thickness of soil layer, [m]
    FineRootBranchFreq_pft         => plt_morph%FineRootBranchFreq_pft           ,& !input  :Fine root brancing frequency, [m-1]
    Root2ndSpecLen_pft             => plt_morph%Root2ndSpecLen_pft               ,& !input  :specific root length secondary axes, [m gC-1]
    MainBranchNum_pft              => plt_morph%MainBranchNum_pft                ,& !input  :number of main branch,[-]
    GrainFillDowreg_brch           => plt_photo%GrainFillDowreg_brch             ,& !input  :grain fill down-regulation of annual plants, [-]
    RootMyco2ndStrutElms_rpvr      => plt_biom%RootMyco2ndStrutElms_rpvr         ,& !inoput :root layer element secondary axes, [g d-2]
    RootMycoNonstElms_rpvr         => plt_biom%RootMycoNonstElms_rpvr            ,& !inoput :root layer nonstructural element, [g d-2]
    RootCO2Autor_pvr               => plt_rbgc%RootCO2Autor_pvr                  ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    RootCO2EmisPot_pvr             => plt_rbgc%RootCO2EmisPot_pvr                ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootRespPotent_pvr             => plt_rbgc%RootRespPotent_pvr                ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]
    RootMaintDef_CO2_pvr           => plt_bgcr%RootMaintDef_CO2_pvr              ,& !inoput :plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]
    Root2ndXNumL_rpvr              => plt_morph%Root2ndXNumL_rpvr                ,& !inoput :root layer number axes, [d-2]
    RootMediumLength_rpvr          => plt_morph%RootMediumLength_rpvr            ,& !inoput :root layer length intermediate size axes, [m d-2]        
    Root2ndProdCytok_rpvr          => plt_rbgc%Root2ndProdCytok_rpvr             ,& !output :cytokinin production rate due to fine root/myco elongation, []
    RootMyco2ndSinkC_rpvr          => plt_rbgc%RootMyco2ndSinkC_rpvr             ,& !output :fine root/myco carbon sink, [gC d-2 h-1]
    Root2ndLen_rpvr                => plt_morph%Root2ndLen_rpvr                  ,& !inoput :root layer length secondary axes, [m d-2]
    Num1stAxesPerStructRootX_pft   => plt_morph%Num1stAxesPerStructRootX_pft     ,& !output :primary root axes number, [d-2]
    Root2ndXNum_rpvr               => plt_morph%Root2ndXNum_rpvr                  & !output :root layer number secondary axes, [d-2]
  )
  !     FRACTION OF SECONDARY ROOT SINK IN SOIL LAYER ATTRIBUTED
  !     TO CURRENT AXIS
  !
  !     Root2ndSink_pvr=total secondary root sink strength
  !     RootSinkC_vr=total root sink strength
  !     FracRoot2ndCSinkL=fraction of secondary root sink strength in axis
  !
  call PrintInfo('beg '//subname)
  !return when the root has no structural nitrogen

  DO NE=1,NumPlantChemElms
    dmass(NE)=RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)
  ENDDO
  
  litrflx2=0._r8
  IF(RootSinkC_vr(N,L).GT.ZERO4Groth_pft(NZ))THEN
    FracRoot2ndCSinkL=Root2ndSink_pvr(N,L,NR)/RootSinkC_vr(N,L)
  ELSE
    !used to initialize fine root growth
    FracRoot2ndCSinkL=1.0_r8
  ENDIF  
  !
  !     N,P CONSTRAINT ON SECONDARY ROOT RESPIRATION FROM
  !     NON-STRUCTURAL C:N:P
  !  
  !
  !
  !     SECONDARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
  !     ROOT STRUCTURAL N
  !
  !     RmSpecPlant=specific maintenance respiration rate (g C g-1 N h-1)
  !     WTRT2N=secondary root N mass
  !     TFN6_vr=temperature function for root maintenance respiration
  !     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
  !     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
  
  !
  Rmaint2nd_CO2=AZMAX1(RmSpecPlant*RootMyco2ndStrutElms_rpvr(ielmn,N,L,NR,NZ))*TFN6_vr(L)

  !if herbaceous root or drought deciduous, maintenance is moisture dependent.
  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ)) .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    Rmaint2nd_CO2=Rmaint2nd_CO2*fRootGrowPSISense
  ENDIF
  !
  !     O2-UNLIMITED SECONDARY ROOT RESPIRATION FROM NON-STRUCTURAL C
  !     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
  !
  !     VMXC=rate constant for nonstructural C oxidation in respiration C     
  !  
  RNonstCO2_OUltd=AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ)*AMIN1(VMXCSpec*FracRoot2ndCSinkL,1._R8))
  !
  !     O2-LIMITED SECONDARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
  !
  RNonstCO2_Oltd               = RNonstCO2_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
  RPotCO2Groth2nd_OUltd        = RNonstCO2_OUltd-Rmaint2nd_CO2
  RPotCO2Groth2nd_Oltd         = RNonstCO2_Oltd-Rmaint2nd_CO2
  RGrowCO2_OUltd               = AZMAX1(RPotCO2Groth2nd_OUltd)*RespElongWatSens
  RGrowCO2_Oltd                = AZMAX1(RPotCO2Groth2nd_Oltd)*RespElongWatSens
  RootMaintDef_CO2_pvr(N,L,NZ) = RootMaintDef_CO2_pvr(N,L,NZ)+AMIN1(RPotCO2Groth2nd_Oltd,0._r8)

  !
  !     SECONDARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
  !     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
  !
  !     FracRoot2ndCSinkL=fraction of secondary root sink strength in axis
  !     FNP=growth respiration limited by non-structural N,P
  !
  ZPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))*FracRoot2ndCSinkL
  PPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))*FracRoot2ndCSinkL
  FNP    = AMIN1(ZPOOLB/CNRTS_pft(NZ),PPOOLB/CPRTS_pft(NZ))*DMRespEff
  
  IF(RGrowCO2_OUltd.GT.0.0_r8)THEN
    RGrowCO2_OUltd=AMIN1(RGrowCO2_OUltd,FNP)
  ELSE
    RGrowCO2_OUltd=0._r8
  ENDIF

  !O2-limited active growth
  IF(RGrowCO2_Oltd.GT.0.0_r8)THEN
    RGrowCO2_Oltd=AMIN1(RGrowCO2_Oltd,FNP*RAutoRootO2Limter_rpvr(N,L,NZ))
  ELSE
    RGrowCO2_Oltd=0._r8
  ENDIF

!
!     TOTAL NON-STRUCTURAL C,N,P USED IN SECONDARY ROOT GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD ENTERED IN 'READQ'
!
!     RootBiomGrosYld_pft=root growth yield
!
  Root2ndProdCytok_rpvr(N,L,NR,NZ)= RGrowCO2_Oltd*RCytok(N)  
  RootMycoNonst4GrowC_OUltd       = RGrowCO2_OUltd/DMRespEff
  RootMycoNonst4GrowC_Oltd        = RGrowCO2_Oltd/DMRespEff

  RootMycoNonst4Grow_OUltd(ielmc) = RootMycoNonst4GrowC_OUltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_OUltd(ielmn) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CNRTW)
  RootMycoNonst4Grow_OUltd(ielmp) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CPRTW)
  RCO2Nonst4Nassim_OUltd          = AZMAX1(1.70_r8*RootMycoNonst4Grow_OUltd(ielmn))

  !O2-limited biomass growth
  RootMycoNonst4Grow_Oltd(ielmc) = RootMycoNonst4GrowC_Oltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_Oltd(ielmn) = AZMAX1(AMIN1(FracRoot2ndCSinkL*RootMycoNonstElms_rpvr(ielmn,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*CNRTW))
  RootMycoNonst4Grow_Oltd(ielmp) = AZMAX1(AMIN1(FracRoot2ndCSinkL*RootMycoNonstElms_rpvr(ielmp,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*CPRTW))  
  RootMycoNonst4Grow_Oltd = RootMycoNonst4Grow_Oltd*fdLext2nd
  RCO2Nonst4Nassim_Oltd          = AZMAX1(1.70_r8*RootMycoNonst4Grow_Oltd(ielmn))

  RCO2T2nd_OUltd = AMIN1(Rmaint2nd_CO2,RNonstCO2_OUltd)+RGrowCO2_OUltd+RCO2Nonst4Nassim_OUltd
  RCO2T2nd_Oltd  = AMIN1(Rmaint2nd_CO2,RNonstCO2_Oltd)+RGrowCO2_Oltd+RCO2Nonst4Nassim_Oltd
  RCO2T2nd_Oltd  = AMIN1(RCO2T2nd_Oltd,AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ)))
  !take CO2 respiration from nonst C
  RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-RCO2T2nd_Oltd
  
  !Consume nonstrucal elements for growth
  DO NE=1,NumPlantChemElms  
    RootMycoNonst4Grow_Oltd(NE)       = AMIN1(RootMycoNonstElms_rpvr(NE,N,L,NZ),RootMycoNonst4Grow_Oltd(NE))
    RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)-RootMycoNonst4Grow_Oltd(NE))
  ENDDO

  ! apply turgor pressure limitation on elongation
  Root2ndExtPot    = RootMycoNonst4Grow_Oltd(ielmc)*FracRootElmAllocm(ielmc,k_fine_comp)*Root2ndSpecLen_pft(N,NZ)
  Root2ndPopExtenz = Root2ndExtPot

  ! net growth=positive growth - senescence
  Root2ndNetGrowthElms = RootMycoNonst4Grow_Oltd
  !
  !     SECONDARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
  !     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
  !     SECONDARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LitrFall
  !  
  call DiagSenes2ndRootAxes(I,J,N,NR,L,NZ,RPotCO2Groth2nd_OUltd,RPotCO2Groth2nd_Oltd,Root2ndPopExtenz,Root2ndNetGrowthElms,&
    litrflx2,RCO2MaintDef2ndStruct_OUltd,RCO2MaintDef2ndStruct_Oltd)
  
  !
  !     CONSUMPTION OF NON-STRUCTURAL C,N,P BY SECONDARY ROOT
  !
  !do final summary
  !
  !     TOTAL SECONDARY ROOT RESPIRATION
  !
  !     Rmaint2nd_CO2=root maintenance respiration
  !     RNonstCO2_OUltd,RNonstCO2_Oltd=respiration from non-structural C unltd,ltd by O2
  !     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
  !     RCO2MaintDef2ndStruct_OUltd,RCO2MaintDef2ndStruct_Oltd=excess maintenance respiration unltd,ltd by O2
  !     RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
  !     RootCO2Autor_pvr=total root respiration
  !     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
  !
  RCO2T2nd_OUltd                   = RCO2T2nd_OUltd +RCO2MaintDef2ndStruct_OUltd
  RCO2T2nd_Oltd                    = RCO2T2nd_Oltd+RCO2MaintDef2ndStruct_Oltd
  RootMyco2ndSinkC_rpvr(N,L,NR,NZ) = RCO2T2nd_Oltd+RootMycoNonst4Grow_Oltd(ielmc)

  RootRespPotent_pvr(N,L,NZ) = RootRespPotent_pvr(N,L,NZ)+RCO2T2nd_OUltd
  RootCO2EmisPot_pvr(N,L,NZ) = RootCO2EmisPot_pvr(N,L,NZ)+RCO2T2nd_Oltd
  RootCO2Autor_pvr(N,L,NZ)   = RootCO2Autor_pvr(N,L,NZ)-RCO2T2nd_Oltd
  RCO2flx2                   = -RCO2T2nd_Oltd
  
  !
  !     SECONDARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
  !
  !     Root2ndPopExtenz=secondary root length extension on whole root axis NR
  !     RootMycoNonst4Grow_Oltd(ielmc)=secondary root C growth ltd by O2
  !     Root2ndSpecLen_pft=specific secondary root length from startq.f
  !     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
  !     Frac2Senes2=fraction of secondary root C to be remobilized
  !     Root2ndLen_rpvr=secondary root length
  !
  !     UPDATE STATE VARIABLES FOR SECONDARY ROOT LENGTH, C, N, P
  !     AND AXIS NUMBER
  !
  !     FineRootBranchFreq_pft=root branching frequency from PFT file for generate secondary roots
  !     Root2ndXNum_rpvr,Root2ndXNumL_rpvr=number of secondary root axes
  !
  Root2ndLen_rpvr(N,L,NR,NZ)=AZERO(Root2ndLen_rpvr(N,L,NR,NZ)+Root2ndPopExtenz)
  DO NE=1,NumPlantChemElms
    RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)+Root2ndNetGrowthElms(NE)
  ENDDO
  
  !secondary/fine root axes (and root hair) addition is a quadratic function of branching frequency
  if(lcoarseroot .and. RootMediumLength_rpvr(L,NR,NZ).GT.0._r8)then
    Num2ndRoots    = FineRootBranchFreq_pft(NZ)    !fine root tips
    Num3rdRoots    = FineRootBranchFreq_pft(NZ)*Num2ndRoots                   !root hair tips
    ! Medium length already includes the population; scale primary length by its axis count.
    NumGroFineAxes = (Num2ndRoots+Num3rdRoots)*(RootMediumLength_rpvr(L,NR,NZ) &
      +DLYR3(L)*Num1stAxesPerStructRootXPOP_pft(NZ))
  else
    Num2ndRoots    = FineRootBranchFreq_pft(NZ)*Num1stAxesPerStructRootXPOP_pft(NZ)    !fine root tips
    Num3rdRoots    = FineRootBranchFreq_pft(NZ)*Num2ndRoots                   !root hair tips
    NumGroFineAxes = (Num2ndRoots+Num3rdRoots)*DLYR3(L)    
  endif  
  MaxNumGroFineAxes           = Root2ndLen_rpvr(N,L,NR,NZ)/Root2ndTipLen4uptk
  Root2ndXNum_rpvr(N,L,NR,NZ) = AMIN1(NumGroFineAxes,MaxNumGroFineAxes)
  Root2ndXNumL_rpvr(N,L,NZ)   = Root2ndXNumL_rpvr(N,L,NZ)+Root2ndXNum_rpvr(N,L,NR,NZ)

  DO NE=1,NumPlantChemElms
    dmass(NE)=dmass(NE)-RootMycoNonstElms_rpvr(NE,N,L,NZ)-RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)-litrflx2(NE)
  ENDDO

  dmass(ielmc)=dmass(ielmc)-RCO2T2nd_Oltd
!  if(I>=140)write(433,*)I*1000+J/24.,dmass,RCO2T2nd_Oltd,litrflx2(ielmc),NR,L,'mass',RCO2T2nd_Oltd+litrflx2(ielmc)
!  if(I>=144 .and. abs(dmass(ielmc))>1.e-6)then
!    write(555,*)RCO2MaintDef2ndStruct_Oltd,'def'
!    stop
!  endif
  call PrintInfo('end '//subname)
  end associate  
  end subroutine Grow2ndRootAxes
!----------------------------------------------------------------------------------------------------

  subroutine DiagSenes2ndRootAxes(I,J,N,NR,L,NZ,RPotCO2Groth2nd_OUltd,RPotCO2Groth2nd_Oltd,Root2ndPopExtenz,&
    Root2ndNetGrowthElms,litrflx2,RCO2MaintDef2ndStruct_OUltd,RCO2MaintDef2ndStruct_Oltd)
  !
  !Description:
  !Diagnose secondary root senescence and remolization
  implicit none
  integer, intent(in) :: I,J,N,NR,L,NZ
  REAL(R8), intent(in) :: RPotCO2Groth2nd_OUltd                        !O2-unlimited growth respiraiton, [gC h-1]
  real(r8), intent(in) :: RPotCO2Groth2nd_Oltd                         !O2-limited growth respiraiton, [gC h-1]
  real(r8), intent(out) :: RCO2MaintDef2ndStruct_OUltd                 !O2-unlimited maintenance respiraiton deficit [gC h-1]
  real(r8), intent(out) :: RCO2MaintDef2ndStruct_Oltd                  !O2-limited maintenance respiraiton deficit [gC h-1]
  real(r8), intent(inout) :: litrflx2(NumPlantChemElms)                !total litterfall flux
  real(r8), intent(inout) :: Root2ndNetGrowthElms(NumPlantChemElms)  !structural biomass growth for 2ndary roots
  real(r8), intent(inout) :: Root2ndPopExtenz                        !potential length extension [m]
  real(r8) :: Root2ndStrutRemob(NumPlantChemElms)                   !remobilization of C,N,P from senescing root
  real(r8) :: Remobl2ndcycl(NumPlantChemElms)                       !recycled nutrients by remobilization
  real(r8) :: Remobl2ndKill(NumPlantChemElms)                       !root lost due to maintenance deficit

  real(r8) :: CCC,CNC,CPC
  real(r8) :: Frac2Senes2              !Fraction of secondary root C to be remobilized
  real(r8) :: RCCE(NumPlantChemElms)
  real(r8) :: RCCC,RCCN,RCCP           !remobilization coefficient for C,N,P
  real(r8) :: FracRemobl, remobilC  
  real(r8) :: dsenecE,dfineE,dwoodyE
  integer :: M,NE

  associate(                                                          &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft             ,& !input  :threshold zero for plang growth calculation, [-]  
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft         ,& !input  :number of main branch,[-]  
    FracRootElmAllocm         => plt_allom%FracRootElmAllocm         ,& !input  :C woody/soft fraction in root,[-]    
    RootNonstructElmConc_pvr  => plt_biom%RootNonstructElmConc_pvr  ,& !input  :root layer nonstructural C concentration, [g g-1]
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr  ,& !input  :root layer element secondary axes, [g d-2]    
    ZERO                      => plt_site%ZERO                       ,& !input  :threshold zero for numerical stability, [-]    
    icwood                    => pltpar%icwood                       ,& !input  :group id of coarse woody litter
    iroot                     => pltpar%iroot                        ,& !input  :group id of plant root litter
    Root2ndLen_rpvr           => plt_morph%Root2ndLen_rpvr           ,& !inoput :root layer length secondary axes, [m d-2]    
    PlantElmAllocMat4Litr     => plt_soilchem%PlantElmAllocMat4Litr  ,& !input  :litter kinetic fraction, [-]    
    k_woody_comp              => pltpar%k_woody_comp                 ,& !input  :woody litter complex id
    k_fine_comp               => pltpar%k_fine_comp                  ,& !input  :fine litter complex id    
    RAutoRootO2Limter_rpvr    => plt_rbgc%RAutoRootO2Limter_rpvr     ,& !input  :O2 constraint to root respiration (0-1), [-]    
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     ,& !inoput :root layer nonstructural element, [g d-2]    
    LitrfallElms_pvr          => plt_bgcr%LitrfallElms_pvr           ,& !inoput :plant LitrFall element, [g d-2 h-1]    
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch        & !input  :plant growth stage, [-]
  )
  !     iPlantCalendar_brch(ipltcal_Emerge),=emergence date
  !     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
  !     CNKI,CPKI=nonstructural N,P inhibition constant on growth
  !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
  !     RCCZR,RCCYR=min,max fractions for root C recycling
  !     RCCXR,RCCQR=max fractions for root N,P recycling
  !
  IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0 .AND. &
    RootNonstructElmConc_pvr(ielmc,N,L,NZ).GT.ZERO)THEN
    CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmn,N,L,NZ) &
      ,RootNonstructElmConc_pvr(ielmn,N,L,NZ)+RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CNKI) &
      ,safe_adb(RootNonstructElmConc_pvr(ielmp,N,L,NZ) &
      ,RootNonstructElmConc_pvr(ielmp,N,L,NZ)+RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CPKI)))
    CNC=AZMAX1(AMIN1(1.0_r8 &
      ,safe_adb(RootNonstructElmConc_pvr(ielmc,N,L,NZ) &
      ,RootNonstructElmConc_pvr(ielmc,N,L,NZ)+RootNonstructElmConc_pvr(ielmn,N,L,NZ)/CNKI)))
    CPC=AZMAX1(AMIN1(1.0_r8 &
      ,safe_adb(RootNonstructElmConc_pvr(ielmc,N,L,NZ) &
      ,RootNonstructElmConc_pvr(ielmc,N,L,NZ)+ RootNonstructElmConc_pvr(ielmp,N,L,NZ)/CPKI)))
  ELSE
    CCC = 0._r8
    CNC = 0._r8
    CPC = 0._r8
  ENDIF
  RCCC = RCCZR+CCC*RCCYR
  RCCN = CNC*RCCXR
  RCCP = CPC*RCCQR
  RCCE(ielmc)=RCCC
  RCCE(ielmn)=(RCCN+(1.0_r8-RCCN)*RCCC)
  RCCE(ielmp)=(RCCP+(1.0_r8-RCCP)*RCCC)
  !
  !     RECOVERY OF REMOBILIZABLE N,P FROM SECONDARY ROOT DURING
  !     REMOBILIZATION DEPENDS ON ROOT NON-STRUCTURAL C:N:P
  !
  !     RPotCO2Groth2nd_OUltd,RPotCO2Groth2nd_Oltd=diff between C respn unltd,ltd by O2 and mntc respn
  !     RCO2MaintDef2ndStruct_OUltd,RCO2MaintDef2ndStruct_Oltd=excess maintenance respiration unltd,ltd by O2
  !     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
  !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
  !     Root2ndStrutRemob(:)=remobilization of C,N,P from senescing root

  !
  remobilC=AZMAX1(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)*RCCC)
  
  IF(-RPotCO2Groth2nd_OUltd.GT.0.0_r8)THEN
    !potential growth respiration is negative, i.e. there is maintenance deficit
    IF(-RPotCO2Groth2nd_OUltd.LT.remobilC)THEN
      !maintenance deficit is offset by remobilC
      RCO2MaintDef2ndStruct_OUltd=-RPotCO2Groth2nd_OUltd
    ELSE
      !maintenance deficit can only be Partially met by remobilC
      RCO2MaintDef2ndStruct_OUltd=remobilC
    ENDIF
  ELSE
    RCO2MaintDef2ndStruct_OUltd=0._r8
  ENDIF

  !there is maintenance deficit,
  IF(-RPotCO2Groth2nd_Oltd.GT.0.0_r8)THEN  
    ! Apply oxygen limitation in both carbon-sufficient and carbon-limited cases.
    RCO2MaintDef2ndStruct_Oltd = AMIN1(-RPotCO2Groth2nd_Oltd,remobilC) &
      *RAutoRootO2Limter_rpvr(N,L,NZ)
  ELSE
    RCO2MaintDef2ndStruct_Oltd=0._r8
  ENDIF

  !Maintenance deficit leads to remobilization/degradation of structural biomass to pay off the deficit
  IF(RCO2MaintDef2ndStruct_Oltd.GT.0.0_r8 .AND. RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
    DO NE=1,NumPlantChemElms
      Root2ndStrutRemob(NE) = RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)*RCCE(NE)
    ENDDO
 
    IF(Root2ndStrutRemob(ielmc).GT.ZERO4Groth_pft(NZ))THEN
      Frac2Senes2=AZMAX1(AMIN1(1.0_r8,RCO2MaintDef2ndStruct_Oltd/Root2ndStrutRemob(ielmc)))
    ELSE
      Frac2Senes2=1.0_r8
    ENDIF
  ELSE
    Root2ndStrutRemob(1:NumPlantChemElms) = 0._r8
    Frac2Senes2                           = 0._r8
  ENDIF
  !
  !     SECONDARY ROOT LitrFall CAUSED BY REMOBILIZATION
  !
  !     CSNC,ZSNC,PSNC=literfall C,N,P
  !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
  !
  litrflx2=0._r8
  if(Frac2Senes2.GT.0._r8)then 
    DO NE=1,NumPlantChemElms
      Remobl2ndKill(NE) = RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)*Frac2Senes2
      Remobl2ndcycl(NE) = Root2ndStrutRemob(NE)*Frac2Senes2
    ENDDO  
    
    !Of the killed biomass, a fraction is added to litter, while the other is added back to reserve
    D6350: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        dsenecE = AZMAX1(Remobl2ndKill(NE)-Remobl2ndcycl(NE))
        dwoodyE = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dsenecE*FracRootElmAllocm(NE,k_woody_comp)
        dfineE  = PlantElmAllocMat4Litr(NE,iroot,M,NZ)*dsenecE*FracRootElmAllocm(NE,k_fine_comp)

        LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)+dwoodyE
        LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)  = LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+dfineE
        litrflx2(NE)                             = litrflx2(NE)+dwoodyE+dfineE

      ENDDO    
    ENDDO D6350

    ! Remobilization is added to nonstrucal metabolites       
    Remobl2ndcycl(ielmc)=Remobl2ndcycl(ielmc)-RCO2MaintDef2ndStruct_Oltd
    DO NE=1,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,N,L,NZ) = RootMycoNonstElms_rpvr(NE,N,L,NZ)+Remobl2ndcycl(NE)
      Root2ndNetGrowthElms(NE)          = AZERO(Root2ndNetGrowthElms(NE)-Remobl2ndKill(NE))
    ENDDO

    !remobilization shrinks roots
    Root2ndPopExtenz=Root2ndPopExtenz-Frac2Senes2*Root2ndLen_rpvr(N,L,NR,NZ)

  endif
  end associate
  end subroutine DiagSenes2ndRootAxes
!----------------------------------------------------------------------------------------------------
  subroutine ReserveMediumRootMaintenance(N,L,NZ,TFN6_vr,fRootGrowPSISense,Root1stSinkTip, &
    FoundRootAxesTip,Root1stTipUpdateFlag,MediumRootMaintCAlloc_OUltd,MediumRootMaintCAlloc_Oltd)
  !Demand-weighted maintenance allocation, independent of the growth sinks.
  !All masses are population totals per ground area. Mycorrhizae use another
  !nonstructural pool and therefore do not enter the plant-root denominator.
  implicit none
  integer, intent(in) :: N,L,NZ
  real(r8), intent(in) :: TFN6_vr(JZ1),fRootGrowPSISense
  real(r8), intent(in) :: Root1stSinkTip(pltpar%MaxNumRootAxes)
  logical, intent(in) :: FoundRootAxesTip(pltpar%jroots,pltpar%MaxNumRootAxes)
  logical, intent(in) :: Root1stTipUpdateFlag(pltpar%MaxNumRootAxes)
  real(r8), intent(out) :: MediumRootMaintCAlloc_OUltd(pltpar%MaxNumRootAxes)
  real(r8), intent(out) :: MediumRootMaintCAlloc_Oltd(pltpar%MaxNumRootAxes)
  real(r8) :: MediumDemand(pltpar%MaxNumRootAxes),TotalDemand,MaintRate,PrimaryN
  real(r8) :: FracRootMCMaintL(pltpar%MaxNumRootAxes),AvailableC,MaintCapacity,VMXCMaint
  integer :: NR
  logical :: TipMaintenance,ThickMaintenance

  associate( &
    iPlantRootProfile_pft              => plt_pheno%iPlantRootProfile_pft, &                     !input  :plant root/growth type code (vascular or nonvascular), [-]
    iPlant2ndGrothPattern_pft          => plt_pheno%iPlant2ndGrothPattern_pft, &                 !input  :plant secondary-growth pattern code, [-]
    iPlantPhenolType_pft               => plt_pheno%iPlantPhenolType_pft, &                      !input  :phenological climate-response type code, [-]
    NumStructuralRootAxes_pft          => plt_morph%NumStructuralRootAxes_pft, &                 !input  :number of structural-root groups represented per PFT, [-]
    NRoot1stTipLay_raxes               => plt_morph%NRoot1stTipLay_raxes, &                      !input  :deepest soil-layer index reached by each structural-root group, [-]
    RootMyco2ndStrutElms_rpvr          => plt_biom%RootMyco2ndStrutElms_rpvr, &                  !input  :population fine-root/mycorrhizal structural C/N/P by layer and group, [g element d-2]
    Root1stLenPP_rpvr                  => plt_morph%Root1stLenPP_rpvr, &                         !input  :length per primary-root axis in each soil layer, [m]
    RootMediumStructElms_rpvr          => plt_biom%RootMediumStructElms_rpvr, &                  !input  :population medium-root structural C/N/P by layer and group, [g element d-2]
    Root1stAxesTipDepz2Surf_pft        => plt_morph%Root1stAxesTipDepz2Surf_pft, &               !input  :primary-root tip depth relative to the column surface, [m]
    CumSoilThickness_vr                => plt_site%CumSoilThickness_vr, &                        !input  :depth to the bottom of each soil layer from the column surface, [m]
    Root1stActStructElms_rpvr          => plt_biom%Root1stActStructElms_rpvr, &                  !input  :population primary-root active-zone C/N/P by layer and group, [g element d-2]
    MaxNumRootLays                     => plt_site%MaxNumRootLays, &                             !input  :maximum soil-layer index accessible to roots, [-]
    flag2ndGrowth_pvr                  => plt_morph%flag2ndGrowth_pvr, &                         !input  :secondary-growth activation flag by layer and structural-root group, [-]
    RootMycoNonstElms_rpvr             => plt_biom%RootMycoNonstElms_rpvr, &                     !inout  :population root/mycorrhizal nonstructural C/N/P; C debited for maintenance reservation, [g element d-2]
    fTGrowRootP_vr                     => plt_pheno%fTGrowRootP_vr, &                            !input  :root growth temperature response used for C mobilization, [-]
    RAutoRootO2Limter_rpvr             => plt_rbgc%RAutoRootO2Limter_rpvr &                      !input  :physiological O2 multiplier on root respiration (0-1), [-]
  )
  MediumRootMaintCAlloc_OUltd = 0._r8
  MediumRootMaintCAlloc_Oltd  = 0._r8

  if(N.NE.ipltroot .or. .not.lcoarseroot)return

  if(.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ), &
    iPlant2ndGrothPattern_pft(NZ)))return

  MaintRate=RmSpecPlant*TFN6_vr(L)
  if(iPlantPhenolType_pft(NZ).eq.iphenotyp_drouhtdecidu)MaintRate=MaintRate*fRootGrowPSISense
  MediumDemand = 0._r8
  TotalDemand  = 0._r8

  do NR=1,NumStructuralRootAxes_pft(NZ)
    if(L>NRoot1stTipLay_raxes(NR,NZ) .or. FoundRootAxesTip(N,NR))cycle

    !fine roots
    TotalDemand=TotalDemand+MaintRate*MAX(0._r8,RootMyco2ndStrutElms_rpvr(ielmn,N,L,NR,NZ))
    !medium roots
    if(Root1stLenPP_rpvr(L,NR,NZ)>0._r8)then
      MediumDemand(NR) = MaintRate*MAX(0._r8,RootMediumStructElms_rpvr(ielmn,L,NR,NZ))
      TotalDemand      = TotalDemand+MediumDemand(NR)
    endif

    !Match the primary-root call gates; tip and thickening maintenance are counted once.
    if(Root1stTipUpdateFlag(NR))cycle

    if(Root1stAxesTipDepz2Surf_pft(NR,NZ)<=CumSoilThickness_vr(L-1))cycle

    !primary roots
    PrimaryN       = MAX(0._r8,Root1stActStructElms_rpvr(ielmn,L,NR,NZ))
    TipMaintenance = (Root1stAxesTipDepz2Surf_pft(NR,NZ)<=CumSoilThickness_vr(L) &
      .or. L==MaxNumRootLays) .and. Root1stSinkTip(NR)>0._r8

    ThickMaintenance=flag2ndGrowth_pvr(L,NR,NZ) &
      .and. Root1stLenPP_rpvr(L,NR,NZ)>0._r8 .and. PrimaryN>=1.e-8_r8
    if(TipMaintenance .or. ThickMaintenance)TotalDemand=TotalDemand+MaintRate*PrimaryN
  enddo
  if(TotalDemand<=0._r8)return

  !Maintenance needs stored C, not a nutrient supply for building new tissue.
  !Retain the temperature/water-dependent mobilization capacity, but exclude
  !the nutrient and grain-fill modifiers used for growth allocation.
  AvailableC                  = AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))
  VMXCMaint                   = VMXC*fTGrowRootP_vr(L,NZ)*fRootGrowPSISense
  MaintCapacity               = MIN(TotalDemand,AvailableC*MIN(1._r8,MAX(0._r8,VMXCMaint)))
  FracRootMCMaintL            = MediumDemand/TotalDemand
  MediumRootMaintCAlloc_OUltd = MaintCapacity*FracRootMCMaintL

  !The oxygen limiter is a physiological factor, not a remaining oxygen pool.
  MediumRootMaintCAlloc_Oltd=MediumRootMaintCAlloc_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
  RootMycoNonstElms_rpvr(ielmc,N,L,NZ)= RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-SUM(MediumRootMaintCAlloc_Oltd)
  
  end associate
  end subroutine ReserveMediumRootMaintenance
!----------------------------------------------------------------------------------------------------
  subroutine GrowMediumRootAxes(yearIJ,N,L,NR,NZ,VMXCSpec,RootSinkC_vr,RootMSink_pvr,fUSE4MR_pvr,&
    TFN6_vr,fRootGrowPSISense,DMRespEff,RespElongWatSens,CNRTW,CPRTW,fdLextM,MediumRootMaintCAlloc_OUltd,MediumRootMaintCAlloc_Oltd,litrflxt)
  !
  !Grow intermediate root axes on structual root NR
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: N,L,NR,NZ
  real(r8), intent(in) :: VMXCSpec
  REAL(R8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)  
  real(r8), intent(in) :: RootMSink_pvr(JZ1,pltpar%MaxNumRootAxes)   !C sink for medium size roots
  real(r8), intent(in) :: fRootGrowPSISense  
  real(r8), intent(in) :: fUSE4MR_pvr(JZ1,pltpar%MaxNumRootAxes) !fine-root activation of medium-root branching, [0-1]
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  real(r8), intent(in) :: DMRespEff
  real(r8), intent(in) :: RespElongWatSens
  real(r8), intent(in) :: CNRTW,CPRTW  
  real(r8), intent(in) :: fdLextM
  real(r8), intent(in) :: MediumRootMaintCAlloc_OUltd,MediumRootMaintCAlloc_Oltd !reserved maintenance C, [gC d-2 timestep-1]
  real(r8), intent(out):: litrflxt(NumPlantChemElms)
  character(len=*), parameter :: subname='GrowMediumRootAxes'
  real(r8) :: FracRootMCSinkL,RmaintMR_CO2
  real(r8) :: MaintNonst_OUltd,MaintNonst_Oltd,GrowthNonst_OUltd
  REAL(R8) :: RNonstCO2_Oltd,RPotCO2GrothMR_OUltd,RPotCO2GrothMR_Oltd
  REAL(R8) :: RGrowCO2_Oltd,RGrowCO2_OUltd,ZPOOLB,PPOOLB,FNP
  real(r8) :: RNonstCO2_OUltd
  real(r8) :: RootMycoNonst4GrowC_OUltd,RootMycoNonst4GrowC_Oltd
  real(r8) :: RootMycoNonst4Grow_OUltd(NumPlantChemElms)
  real(r8) :: RootMycoNonst4Grow_Oltd(NumPlantChemElms)
  real(r8) :: RootMRExtPot          !total root extension for axis NR
  real(r8) :: RCO2TMR_OUltd,RCO2TMR_Oltd
  real(r8) :: RootMRNetGrowthElms(NumPlantChemElms)
  real(r8) :: PotentialMediumAxes,TargetMediumAxes  !population branching capacity and activated target, [d-2]
  real(r8) :: C4Elongation,cyto_scal,fctyok,f_elong,f_thick

  integer  :: M,NE

  associate(                                                                      &
    DLYR3                          => plt_site%DLYR3                             ,& !input  :vertical thickness of soil layer, [m]
    ZERO4Groth_pft                 => plt_biom%ZERO4Groth_pft                    ,& !input  :threshold zero for plang growth calculation, [-]    
    iPlantPhenolType_pft           => plt_pheno%iPlantPhenolType_pft             ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft            ,& !input  :plant growth type (vascular, non-vascular),[-]    
    MediumRootBranchFreq_pft       => plt_morph%MediumRootBranchFreq_pft         ,& !input  :Fine root brancing frequency, [m-1]
    k_fine_comp                    => pltpar%k_fine_comp                         ,& !input  :fine litter complex id    
    iroot                          => pltpar%iroot                               ,& !input  :group id of plant root litter
    CNRTS_pft                      => plt_allom%CNRTS_pft                        ,& !input  :root N:C ratio x root growth yield, [-]
    CPRTS_pft                      => plt_allom%CPRTS_pft                        ,& !input  :root P:C ratio x root growth yield, [-]    
    RootBiomGrosYld_pft            => plt_allom%RootBiomGrosYld_pft              ,& !input  :root growth yield, [g g-1]    
    Root1stSpecLen_pft             => plt_morph%Root1stSpecLen_pft               ,& !input  :specific root length primary axes, [m root gC-1]    
    SeedDepth_pft                  => plt_morph%SeedDepth_pft                    ,& !input  :seeding depth, [m]
    CRootActVolPerMassC_pft        => plt_morph%CRootActVolPerMassC_pft          ,& !input  :coarse root active zone volume:mass ratio, [m3 gC-1]        
    enh_cyto_pft                   => plt_morph%enh_cyto_pft                     ,& !input  :cytokinin sensitivity of coarse root thickening
    fTgrowRootP_vr                 => plt_pheno%fTgrowRootP_vr                   ,& !input  :root layer temperature growth functiom, [-]    
    GrainFillDowreg_brch           => plt_photo%GrainFillDowreg_brch             ,& !input  :grain fill down-regulation of annual plants, [-]    
    PlantPopuLive_pft              => plt_site%PlantPopuLive_pft                 ,& !input  :plant population, [d-2]    
    MainBranchNum_pft              => plt_morph%MainBranchNum_pft                ,& !input  :number of main branch,[-]
    RAutoRootO2Limter_rpvr         => plt_rbgc%RAutoRootO2Limter_rpvr            ,& !input  :O2 constraint to root respiration (0-1), [-]
    PlantElmAllocMat4Litr          => plt_soilchem%PlantElmAllocMat4Litr         ,& !input  :litter kinetic fraction, [-]
    RootMaintDef_CO2_pvr           => plt_bgcr%RootMaintDef_CO2_pvr              ,& !inoput :plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]
    Root1stLenPP_rpvr              => plt_morph%Root1stLenPP_rpvr                ,& !input :length per primary-root axis in layer, [m]
    RootMyco2ndStrutElms_rpvr      => plt_biom%RootMyco2ndStrutElms_rpvr         ,& !input :population fine-root/mycorrhizal structural elements, [g d-2]
    Root2ndLen_rpvr                => plt_morph%Root2ndLen_rpvr                  ,& !input :population fine-root/mycorrhizal length, [m d-2]
    RootMRProdCytok_rpvr           => plt_rbgc%RootMRProdCytok_rpvr              ,& !input  :cytokinin production rate due to medium roots metabolism, [gC CK h-1]    
    Num1stAxesPerStructRootXPOP_pft=> plt_morph%Num1stAxesPerStructRootXPOP_pft  ,& !output :population primary root axes number on one structrual axis, [d-2]        
    RootCO2Autor_pvr               => plt_rbgc%RootCO2Autor_pvr                  ,& !inoput :root respiration constrained by O2, [g d-2 h-1]    
    CytokininMRConc_rpvr           => plt_rbgc%CytokininMRConc_rpvr              ,& !output :cytokinin concentration in medium size roots, [gC m-3 H2O]    
    RootMediumStructElms_rpvr      => plt_biom%RootMediumStructElms_rpvr         ,& !inoput :root layer element for medium size root axes, [g d-2]
    RootMyco1stElm_raxs            => plt_biom%RootMyco1stElm_raxs               ,& !inoput :root layer primary axes structrual element, [g d-2]
    RootMycoNonstElms_rpvr         => plt_biom%RootMycoNonstElms_rpvr            ,& !inoput :root layer nonstructural element, [g d-2]    
    RootCO2EmisPot_pvr             => plt_rbgc%RootCO2EmisPot_pvr                ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]    
    RootRespPotent_pvr             => plt_rbgc%RootRespPotent_pvr                ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]    
    RootMediumXNum_rpvr            => plt_morph%RootMediumXNum_rpvr              ,& !inoput :Number of medium size root axes in layer for structrual axes, [d-2]
    RootMediumLength_rpvr          => plt_morph%RootMediumLength_rpvr            ,& !inoput :total length for medium size root for axis NR, [m d-2]    
    RootMediumRadius_rpvr          => plt_morph%RootMediumRadius_rpvr            ,& !inoput :root layer radius for medium size axes, [m]
    LitrfallElms_pvr               => plt_bgcr%LitrfallElms_pvr                  ,& !inoput :plant LitrFall element, [g d-2 h-1]
    FracRootElmAllocm              => plt_allom%FracRootElmAllocm                 & !input  :C woody fraction in root,[-]                                                     
  )
  call PrintInfo('beg '//subname)
  litrflxt = 0._r8

  !Recheck living fine roots after the preceding fine/primary-root updates.
  !A sink computed before withdrawal must not enable growth after fine-root death.
  IF(RootSinkC_vr(N,L).GT.ZERO4Groth_pft(NZ) .and. &
     RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ).GT.ZERO4Groth_pft(NZ) .and. &
     Root2ndLen_rpvr(N,L,NR,NZ).GT.0._r8 .and. fUSE4MR_pvr(L,NR).GT.0._r8)THEN
    FracRootMCSinkL=RootMSink_pvr(L,NR)/RootSinkC_vr(N,L)
  else
    !A zero growth sink does not remove the maintenance demand of living roots.
    FracRootMCSinkL=0._r8
  ENDIF  
  !The branching target is a growth diagnostic, not the surviving axis count.
  PotentialMediumAxes = AZMAX1(MediumRootBranchFreq_pft(NZ)*Root1stLenPP_rpvr(L,NR,NZ) &
    *Num1stAxesPerStructRootXPOP_pft(NZ))
  TargetMediumAxes = MIN(1._r8,AZMAX1(fUSE4MR_pvr(L,NR)))*PotentialMediumAxes
  if(RootMediumStructElms_rpvr(ielmc,L,NR,NZ).LE.0._r8 .or. &
     RootMediumLength_rpvr(L,NR,NZ).LE.0._r8)RootMediumXNum_rpvr(L,NR,NZ)=0._r8
  RmaintMR_CO2                 = AZMAX1(RmSpecPlant*RootMediumStructElms_rpvr(ielmn,L,NR,NZ))*TFN6_vr(L)

  !if herbaceous root or drought deciduous, maintenance is moisture dependent.
  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ)) .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    RmaintMR_CO2=RmaintMR_CO2*fRootGrowPSISense
  ENDIF
  !Restore this group's reservation exactly once. Earlier primary-root withdrawal
  !may have reduced its demand; unused reserved C remains in the shared pool.
  RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)+MediumRootMaintCAlloc_Oltd
  MaintNonst_OUltd = MIN(RmaintMR_CO2,MediumRootMaintCAlloc_OUltd)
  MaintNonst_Oltd  = MaintNonst_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
  !Maintenance uses the growth-based mobilization allowance first. Add only
  !its remainder, so reserving maintenance does not count that allowance twice.
  GrowthNonst_OUltd = AZMAX1(AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ)) &
    *MIN(VMXCSpec*FracRootMCSinkL,1._r8)-MaintNonst_OUltd)
  RNonstCO2_OUltd = MaintNonst_OUltd+GrowthNonst_OUltd
  RNonstCO2_Oltd  = MaintNonst_Oltd+GrowthNonst_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
  RPotCO2GrothMR_OUltd         = RNonstCO2_OUltd-RmaintMR_CO2
  RPotCO2GrothMR_Oltd          = RNonstCO2_Oltd-RmaintMR_CO2
  RGrowCO2_OUltd               = AZMAX1(RPotCO2GrothMR_OUltd)*RespElongWatSens
  RGrowCO2_Oltd                = AZMAX1(RPotCO2GrothMR_Oltd)*RespElongWatSens
  RootMaintDef_CO2_pvr(N,L,NZ) = RootMaintDef_CO2_pvr(N,L,NZ)+AMIN1(RPotCO2GrothMR_Oltd,0._r8)
  !
  !     SECONDARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
  !     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
  !
  !     FracRoot2ndCSinkL=fraction of secondary root sink strength in axis
  !     FNP=growth respiration limited by non-structural N,P
  !
  ZPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))*FracRootMCSinkL
  PPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))*FracRootMCSinkL
  FNP    = AMIN1(ZPOOLB/CNRTS_pft(NZ),PPOOLB/CPRTS_pft(NZ))*DMRespEff

  IF(RGrowCO2_OUltd.GT.0.0_r8)THEN
    RGrowCO2_OUltd=AMIN1(RGrowCO2_OUltd,FNP)
  ELSE
    RGrowCO2_OUltd=0._r8
  ENDIF

  IF(RGrowCO2_Oltd.GT.0.0_r8)THEN
    RGrowCO2_Oltd=AMIN1(RGrowCO2_Oltd,FNP*RAutoRootO2Limter_rpvr(N,L,NZ))
  ELSE
    RGrowCO2_Oltd=0._r8
  ENDIF

  RootMycoNonst4GrowC_OUltd       = RGrowCO2_OUltd/DMRespEff
  RootMycoNonst4GrowC_Oltd        = RGrowCO2_Oltd/DMRespEff

  RootMycoNonst4Grow_OUltd(ielmc) = RootMycoNonst4GrowC_OUltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_OUltd(ielmn) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CNRTW)
  RootMycoNonst4Grow_OUltd(ielmp) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CPRTW)

  !O2-limited biomass growth
  RootMycoNonst4Grow_Oltd(ielmc) = RootMycoNonst4GrowC_Oltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_Oltd(ielmn) = AZMAX1(AMIN1(FracRootMCSinkL*RootMycoNonstElms_rpvr(ielmn,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*CNRTW))
  RootMycoNonst4Grow_Oltd(ielmp) = AZMAX1(AMIN1(FracRootMCSinkL*RootMycoNonstElms_rpvr(ielmp,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*CPRTW))  
  RootMycoNonst4Grow_Oltd = RootMycoNonst4Grow_Oltd*fdLextM

  RCO2TMR_OUltd = AMIN1(RmaintMR_CO2,RNonstCO2_OUltd)+RGrowCO2_OUltd
  RCO2TMR_Oltd  = AMIN1(RmaintMR_CO2,RNonstCO2_Oltd)+RGrowCO2_Oltd
  RCO2TMR_Oltd  = AMIN1(RCO2TMR_Oltd,AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ)))   !>0

  !take CO2 respiration from nonst C
  RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-RCO2TMR_Oltd

  !Consume nonstrucal elements for growth: elongation+thickening
  DO NE=1,NumPlantChemElms  
    RootMycoNonst4Grow_Oltd(NE)       = AMIN1(RootMycoNonstElms_rpvr(NE,N,L,NZ),RootMycoNonst4Grow_Oltd(NE))
    RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)-RootMycoNonst4Grow_Oltd(NE))
  ENDDO
  
  ! net growth=positive growth - senescence
  RootMRNetGrowthElms = RootMycoNonst4Grow_Oltd
  !has maitenance deficit
  if(-RPotCO2GrothMR_Oltd.GT.0._r8)then
    call RemobilizeMediumRoots(yearIJ,N,L,NR,NZ,RPotCO2GrothMR_Oltd,RPotCO2GrothMR_OUltd,RCO2TMR_Oltd,RCO2TMR_OUltd,&
      RootMRNetGrowthElms,litrflxt)
  endif

  RootRespPotent_pvr(N,L,NZ) = RootRespPotent_pvr(N,L,NZ)+RCO2TMR_OUltd
  RootCO2EmisPot_pvr(N,L,NZ) = RootCO2EmisPot_pvr(N,L,NZ)+RCO2TMR_Oltd
  RootCO2Autor_pvr(N,L,NZ)   = RootCO2Autor_pvr(N,L,NZ)-RCO2TMR_Oltd

  IF(RootMRNetGrowthElms(ielmc).LT.0.0_r8)THEN              
    call WithDrawMediumRoots(N,NZ,L,NR,RootMRNetGrowthElms,litrflxt)    
  ELSEIF(RootMRNetGrowthElms(ielmc).GT.0.0_r8)THEN              
    !Only actual growth can establish/add axes. A weaker sink cannot kill them.
    RootMediumXNum_rpvr(L,NR,NZ)=MAX(RootMediumXNum_rpvr(L,NR,NZ),TargetMediumAxes)
    !allocate netgrowth to elongation and thickening  
    cyto_scal = CytokininMRConc_rpvr(L,NR,NZ)
    fctyok    = func_cyto(cyto_scal,enh_cyto_pft(NZ))
    f_thick   = fctyok*AZMAX1(1._r8-RootMediumRadius_rpvr(L,NR,NZ)/dmax)
    f_elong   = 1._r8-f_thick

    ! apply turgor pressure limitation on elongation
    C4Elongation                   = RootMRNetGrowthElms(ielmc)*f_elong
    RootMRExtPot                   = C4Elongation*FracRootElmAllocm(ielmc,k_fine_comp)*Root1stSpecLen_pft(N,NZ)
    RootMediumLength_rpvr(L,NR,NZ) = RootMediumLength_rpvr(L,NR,NZ)+RootMRExtPot    
   
    ! now update the biomass
    RootMediumStructElms_rpvr(:,L,NR,NZ) =  RootMediumStructElms_rpvr(:,L,NR,NZ) +RootMRNetGrowthElms
    ! compute the mean radius 
    RootMediumRadius_rpvr(L,NR,NZ)=sqrt(RootMediumStructElms_rpvr(ielmc,L,NR,NZ)*CRootActVolPerMassC_pft(NZ) &
    /(PiCON*RootMediumLength_rpvr(L,NR,NZ))) 
  ENDIF
  
  RootMRProdCytok_rpvr(L,NR,NZ)=RGrowCO2_Oltd*RCytok(ipltroot)*0.1_r8  
  call PrintInfo('end '//subname)
  end associate
  end subroutine GrowMediumRootAxes
!----------------------------------------------------------------------------------------------------
  subroutine RemobilizeMediumRoots(yearIJ,N,L,NR,NZ,RPotCO2GrothMR_Oltd,RPotCO2GrothMR_OUltd,RCO2TMR_Oltd,RCO2TMR_OUltd,&
    RootMRNetGrowthElms,litrflxt)

  !do remobilization under negative growth  
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in)    :: N,L,NR,NZ
  real(r8),intent(inout) :: RPotCO2GrothMR_Oltd   !potential growth respiration, < 0 meaning in deficiency
  real(r8),intent(inout) :: RPotCO2GrothMR_OUltd
  real(r8),intent(inout)  :: RCO2TMR_Oltd
  real(r8),intent(inout) :: RCO2TMR_OUltd
  real(r8),intent(inout) :: litrflxt(NumPlantChemElms)
  real(r8), intent(inout) :: RootMRNetGrowthElms(NumPlantChemElms)  
  real(r8) :: RootMRcylc(NumPlantChemElms)
  real(r8) :: RootMRKill(NumPlantChemElms)
  real(r8) :: RootMRStrutRemob(NumPlantChemElms)
  integer :: M,NE
  real(r8) :: CCC,CNC,CPC
  real(r8) :: Frac2Senes1
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: RCO2MDefMRStructPayoff_Oltd  !>0
  real(r8) :: RCO2MDefMRStructPayoff_OUltd
  real(r8) :: RCCE(NumPlantChemElms)
  real(r8) :: dsenecE,dlitrfall

  character(len=*), parameter :: subname='RemobilizeMediumRoots'
  associate(                                                                      &
    RootNonstructElmConc_pvr       => plt_biom%RootNonstructElmConc_pvr          ,& !input  :root layer nonstructural C concentration, [g g-1]
    ZERO4Groth_pft                 => plt_biom%ZERO4Groth_pft                    ,& !input  :threshold zero for plang growth calculation, [-]
    ZERO                           => plt_site%ZERO                              ,& !input  :threshold zero for numerical stability, [-]
    icwood                         => pltpar%icwood                              ,& !input  :group id of coarse woody litter
    k_woody_comp                   => pltpar%k_woody_comp                        ,& !input  :woody litter complex id
    k_fine_comp                    => pltpar%k_fine_comp                         ,& !input  :fine litter complex id
    iroot                          => pltpar%iroot                               ,& !input  :group id of plant root litter
    FracRootElmAllocm              => plt_allom%FracRootElmAllocm                ,& !input  :root element woody/fine fraction in biomass,[-]
    PlantElmAllocMat4Litr          => plt_soilchem%PlantElmAllocMat4Litr         ,& !input  :litter kinetic fraction, [-]
    RAutoRootO2Limter_rpvr         => plt_rbgc%RAutoRootO2Limter_rpvr            ,& !input  :O2 constraint to root respiration (0-1), [-]
    iPlantCalendar_brch            => plt_pheno%iPlantCalendar_brch              ,& !input  :plant growth stage, [-]
    MainBranchNum_pft              => plt_morph%MainBranchNum_pft                ,& !input  :number of main branch,[-]
    RootMediumStructElms_rpvr      => plt_biom%RootMediumStructElms_rpvr         ,& !inoput :root layer element for medium size root axes, [g d-2]
    RootMycoNonstElms_rpvr         => plt_biom%RootMycoNonstElms_rpvr            ,& !inoput :root layer nonstructural element, [g d-2]
    LitrfallElms_pvr               => plt_bgcr%LitrfallElms_pvr                   & !inoput :plant LitrFall element, [g d-2 h-1]
  )
  call PrintInfo('beg '//subname)
  RootMRcylc=0._r8

  IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0 &
    .AND. RootNonstructElmConc_pvr(ielmc,N,L,NZ).GT.ZERO)THEN
    CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmn,N,L,NZ), &
      RootNonstructElmConc_pvr(ielmn,N,L,NZ)+RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CNKI), &
      safe_adb(RootNonstructElmConc_pvr(ielmp,N,L,NZ),RootNonstructElmConc_pvr(ielmp,N,L,NZ) &
        +RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CPKI)))
    CNC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmc,N,L,NZ), &
      RootNonstructElmConc_pvr(ielmc,N,L,NZ)+RootNonstructElmConc_pvr(ielmn,N,L,NZ)/CNKI)))
    CPC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmc,N,L,NZ), &
      RootNonstructElmConc_pvr(ielmc,N,L,NZ)+RootNonstructElmConc_pvr(ielmp,N,L,NZ)/CPKI)))
  ELSE
    CCC=0._r8;CNC=0._r8;CPC=0._r8
  ENDIF
  RCCC        = RCCZR+CCC*RCCYR
  RCCN        = CNC*RCCXR
  RCCP        = CPC*RCCQR
  RCCE(ielmc) = AMIN1(RCCC,1._R8)
  RCCE(ielmn) = AMIN1(RCCN+(1.0_r8-RCCN)*RCCC,1._R8)
  RCCE(ielmp) = AMIN1(RCCP+(1.0_r8-RCCP)*RCCC,1._R8)

  !insufficient energy for maintenance
  IF(-RPotCO2GrothMR_OUltd.GT.0.0_r8)THEN
    IF(-RPotCO2GrothMR_OUltd.LT.RootMediumStructElms_rpvr(ielmc,L,NR,NZ)*RCCE(ielmc))THEN
      RCO2MDefMRStructPayoff_OUltd=-RPotCO2GrothMR_OUltd
    ELSE
      RCO2MDefMRStructPayoff_OUltd=AZMAX1(RootMediumStructElms_rpvr(ielmc,L,NR,NZ)*RCCE(ielmc))
    ENDIF
    RPotCO2GrothMR_OUltd = RPotCO2GrothMR_OUltd+RCO2MDefMRStructPayoff_OUltd
  ELSE
    RCO2MDefMRStructPayoff_OUltd=0._r8
  ENDIF 

  !has maintenance deficit
  IF(-RPotCO2GrothMR_Oltd.GT.0.0_r8)THEN
    ! Use this oxygen-limited payment for both respiration and tissue withdrawal.
    RCO2MDefMRStructPayoff_Oltd = AMIN1(-RPotCO2GrothMR_Oltd, &
      AZMAX1(RootMediumStructElms_rpvr(ielmc,L,NR,NZ)*RCCC))*RAutoRootO2Limter_rpvr(N,L,NZ)
    !offset some CO2
    RPotCO2GrothMR_Oltd = RPotCO2GrothMR_Oltd+RCO2MDefMRStructPayoff_Oltd
  ELSE
    RCO2MDefMRStructPayoff_Oltd=0._r8
  ENDIF

  IF(RCO2MDefMRStructPayoff_Oltd.GT.0.0_r8 .AND. RootMediumStructElms_rpvr(ielmc,L,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
    !remobilization upon starvation-induced root retreat
    DO NE=1,NumPlantChemElms
      RootMRStrutRemob(NE) = RootMediumStructElms_rpvr(NE,L,NR,NZ)*RCCE(NE)
    ENDDO

    !maintenance deficit is paid by remobilization
    IF(RootMRStrutRemob(ielmc).GT.ZERO4Groth_pft(NZ))THEN
      Frac2Senes1=AZMAX1(AMIN1(1.0_r8,RCO2MDefMRStructPayoff_Oltd/RootMRStrutRemob(ielmc)))
    ELSE
      Frac2Senes1=1.0_r8
    ENDIF
  ELSE
    RootMRStrutRemob(1:NumPlantChemElms) = 0._r8
    Frac2Senes1                          = 0._r8
  ENDIF

  if(Frac2Senes1.GT.0._r8)then
    DO NE=1,NumPlantChemElms          
      RootMRKill(NE) = RootMediumStructElms_rpvr(NE,L,NR,NZ)*Frac2Senes1  
      RootMRcylc(NE) = RootMRStrutRemob(NE)*Frac2Senes1
    ENDDO

    DO M=1,jsken
      DO NE=1,NumPlantChemElms          
        dsenecE      = AZMAX1(RootMRKill(NE)-RootMRcylc(NE))
        dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,icwood,M,NZ)*FracRootElmAllocm(NE,k_woody_comp)
        litrflxt(NE) = litrflxt(NE)+dlitrfall

        LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) +dlitrfall
        
        dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,iroot,M,NZ)*FracRootElmAllocm(NE,k_fine_comp)
        litrflxt(NE) = litrflxt(NE)+dlitrfall

        LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+dlitrfall
      ENDDO
    ENDDO

    !remobilization reduces growth    
    DO NE=1,NumPlantChemElms
      RootMRNetGrowthElms(NE)=RootMRNetGrowthElms(NE)-RootMRKill(NE)
    ENDDO 

    RootMRcylc(ielmc)=RootMRcylc(ielmc)-RCO2MDefMRStructPayoff_Oltd

    DO NE=1,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootMRcylc(NE))
    ENDDO 
  endif

  ! add CO2 respiraiton from remobilization  
  RCO2TMR_Oltd      = RCO2TMR_Oltd+RCO2MDefMRStructPayoff_Oltd
  RCO2TMR_OUltd      = RCO2TMR_OUltd+RCO2MDefMRStructPayoff_OUltd

  call PrintInfo('end '//subname)
  end associate

  end subroutine RemobilizeMediumRoots
!----------------------------------------------------------------------------------------------------
  subroutine GrowRootMycoAxes(yearIJ,N,L,Lnext,NZ,FoundRootAxesTip,Root1stTipUpdateFlag,TFN6_vr,&
    RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,Root1stSinkTip,RootMSink_pvr,fUSE4MR_pvr,CNRTW,CPRTW,fRootGrowPSISense,&
    litrflx,RCO2flx)
  !
  !Description
  !Grow root axes in laye L  
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  INTEGER , INTENT(IN) :: N         !root type id
  integer , intent(in) :: L         !current root layer
  integer , intent(in) :: Lnext     !next root growable layer, soil porosity >0 
  integer , intent(in) :: NZ     !pft id
  real(r8), intent(in) :: TFN6_vr(JZ1)
  REAL(R8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root1stSinkTip(pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: RootMSink_pvr(JZ1,pltpar%MaxNumRootAxes)   
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: fUSE4MR_pvr(JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: CNRTW,CPRTW
  real(r8), intent(in) :: fRootGrowPSISense          !growth function of root water potential
  logical , intent(inout) :: FoundRootAxesTip(pltpar%jroots,pltpar%MaxNumRootAxes)  
  logical , intent(inout) :: Root1stTipUpdateFlag(pltpar%MaxNumRootAxes)  
  real(r8), intent(inout) :: litrflx(NumPlantChemElms)
  real(r8), intent(inout) :: RCO2flx

  character(len=*), parameter :: subname='GrowRootMycoAxes'

  real(r8) :: Root2ndStrutRemob(NumPlantChemElms)    
  real(r8) :: DMRespEff                          !respiraiton efficiency for new root growth
  real(r8) :: RespElongWatSens                   !respiration function of root water potential
  real(r8) :: CNG,CPG
  real(r8) :: Nutstress4GrossResp                           !nutrient stress on root mobilization of nonstrucal C.
  integer  :: NR
  logical  :: FoundPrimaryRootsLayer
  real(r8) :: litrflxt(NumPlantChemElms),litrflx2(NumPlantChemElms)
  real(r8) :: RCO2flxt,RCO2flx2
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: massr1st(NumPlantChemElms),massr1st1(NumPlantChemElms)
  real(r8) :: massr2nd(NumPlantChemElms),massr2nd1(NumPlantChemElms)
  real(r8) :: massnonst(NumPlantChemElms),massnonst1(NumPlantChemElms)
  real(r8) :: massnodul(NumPlantChemElms),massnodul1(NumPlantChemElms)
  real(r8) :: dlitrfall(NumPlantChemElms)
  real(r8) :: fdLext1st,fdLext2nd          !soil compaction modifier of root elongation 
  real(r8) :: RTN2X,RTN2Y,tmpval,dmassN
  real(r8) :: VMXCSpec                     !specific remobilization rate of nonstrucal carbon [h-1]
  real(r8) :: MediumMediumRootMaintCAlloc_OUltd(pltpar%MaxNumRootAxes)
  real(r8) :: MediumMediumRootMaintCAlloc_Oltd(pltpar%MaxNumRootAxes)

  !begin_execution
  associate(                                                                      &
    FracRootElmAllocm              => plt_allom%FracRootElmAllocm                ,& !input  :C woody fraction in root,[-]
    ZERO4Groth_pft                 => plt_biom%ZERO4Groth_pft                    ,& !input  :threshold zero for plang growth calculation, [-]
    ZERO                           => plt_site%ZERO                              ,& !input  :threshold zero for numerical stability, [-]
    Root1stLenPP_rpvr              => plt_morph%Root1stLenPP_rpvr                ,& !input  :primary root axis length in soil layer, [m d-2]
    Root2ndRadius_rpvr             => plt_morph%Root2ndRadius_rpvr               ,& !input  :root layer radius secondary axes, [m]
    k_fine_comp                    => pltpar%k_fine_comp                         ,& !input  :fine litter complex id
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft            ,& !input  :plant growth type (vascular, non-vascular),[-]
    iPlant2ndGrothPattern_pft      => plt_pheno%iPlant2ndGrothPattern_pft        ,& !input  :plant expression of secondary growth, [-]                    
    CumSoilThickness_vr            => plt_site%CumSoilThickness_vr               ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    NRoot1stTipLay_raxes           => plt_morph%NRoot1stTipLay_raxes             ,& !input  :maximum soil layer number for root axes, [-]
    RootBiomGrosYld_pft            => plt_allom%RootBiomGrosYld_pft              ,& !input  :root growth yield, [g g-1]
    fTGrowRootP_vr                 => plt_pheno%fTgrowRootP_vr                   ,& !input  :root layer temperature growth functiom, [-]    
    Num1stAxesPerStructRootX_pft   => plt_morph%Num1stAxesPerStructRootX_pft     ,& !input :primary root axes number, [d-2]
    PlantPopuLive_pft              => plt_site%PlantPopuLive_pft                 ,& !input  :plant population, [d-2]
    DLYR3                          => plt_site%DLYR3                             ,& !input  :vertical thickness of soil layer, [m]
    GrainFillDowreg_brch           => plt_photo%GrainFillDowreg_brch             ,& !input  :grain fill down-regulation of annual plants, [-]    
    RootNonstructElmConc_pvr       => plt_biom%RootNonstructElmConc_pvr          ,& !input  :root layer nonstructural C concentration, [g g-1]
    Root2ndLen_rpvr                => plt_morph%Root2ndLen_rpvr                  ,& !input  :root layer length secondary axes, [m d-2]
    RootMyco2ndStrutElms_rpvr      => plt_biom%RootMyco2ndStrutElms_rpvr         ,& !input  :root layer element secondary axes, [g d-2]
    RootMyco1stStrutElms_rpvr      => plt_biom%RootMyco1stStrutElms_rpvr         ,& !input  :root layer element primary axes, [g d-2]
    MainBranchNum_pft              => plt_morph%MainBranchNum_pft                ,& !input  :number of main branch,[-]    
    flag2ndGrowth_pvr              => plt_morph%flag2ndGrowth_pvr                ,& !input  :flag for secondary growth of primary roots, [-]
    PSIRootTurg_vr                 => plt_ew%PSIRootTurg_vr                      ,& !input  :root turgor water potential, [Mpa]
    NumStructuralRootAxes_pft      => plt_morph%NumStructuralRootAxes_pft        ,& !input  :number of structural root axes,[-]
    RootPorosity_pft               => plt_morph%RootPorosity_pft                 ,& !input  :root porosity, [m3 m-3]
    SoilBulkModulus4RootPent_vr    => plt_soilchem%SoilBulkModulus4RootPent_vr   ,& !input  :soil hydraulic resistance, [MPa h m-2]
    Root2ndXSecArea_pft            => plt_morph%Root2ndXSecArea_pft              ,& !input  :root cross-sectional area of secondary axes, [m2]
    GroSrcRootStress_pvr           => plt_rbgc%GroSrcRootStress_pvr              ,& !output :root growth stress due to nutrient and water, [-]
    NMaxRootBotLayer_pft           => plt_morph%NMaxRootBotLayer_pft             ,& !inoput :maximum soil layer number for all root axes, [-]
    Root1stAxesTipDepz2Surf_pft    => plt_morph%Root1stAxesTipDepz2Surf_pft      ,& !output :plant primary depth relative to column surface, [m]
    ROOTNLim_rpvr                  => plt_biom%ROOTNLim_rpvr                     ,& !output :root N-limitation, 0->1 weaker limitation, [-]
    ROOTPLim_rpvr                  => plt_biom%ROOTPLim_rpvr                      & !output :root P-limitation, 0->1 weaker limitation, [-]
  )
  call PrintInfo('beg '//subname)
  !     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
  !     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
  !
  !     SoilBulkModulus4RootPent_vr,SoilRest4Root2ndPentrate=soil resistance to secondary root penetration (MPa)
  !     Root2ndRadius_rpvr=secondary root radius
  !     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
  !     PSIRoot_pvr,PSIRootTurg_vr=root total,turgor water potential
  !     DMRT=root growth yield
  !
  !   call SumRootBiome(NZ,mass_inital,massr1st1,massr2nd1,massnonst1,massnodul1)
  !set fine root elongation zone to 2 mm, an intermediate value for fine roots and mycorrhizae

  CALL RootElongationWaterFunc(N,L,NZ,2.e-3_r8,fdLext2nd,RespElongWatSens)
  if(.not.lsoilCompaction)fdLext2nd=1._r8
  
  !respiration fraction
  DMRespEff         = 1.0_r8-RootBiomGrosYld_pft(NZ)     !1-[gC root/gC nonst] = [gC resp/gC nonst]

  !nutrient limitation factor
  IF(RootNonstructElmConc_pvr(ielmc,N,L,NZ).GT.ZERO)THEN
    CNG  = RootNonstructElmConc_pvr(ielmn,N,L,NZ)/(RootNonstructElmConc_pvr(ielmn,N,L,NZ)+RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CNKI)
    CPG  = RootNonstructElmConc_pvr(ielmp,N,L,NZ)/(RootNonstructElmConc_pvr(ielmp,N,L,NZ)+RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CPKI)
    Nutstress4GrossResp = AMIN1(CNG,CPG)
  ELSE
    CNG =0._r8
    CPG =0._r8
    Nutstress4GrossResp=1.0_r8
  ENDIF
  ROOTNLim_rpvr(N,L,NZ)=CNG
  ROOTPLim_rpvr(N,L,NZ)=CPG

  if(N.EQ.ipltroot .and. lsoilCompaction)then
    !10 mm for primary root elongation zone
    CALL RootElongationWaterFunc(N,L,NZ,10.e-3_r8,fdLext1st)    
  else
    fdLext1st=1._r8    
  endif  
  !VMXC=rate constant for nonstructural C oxidation in respiration C     
  VMXCSpec=VMXC*fTGrowRootP_vr(L,NZ)*fRootGrowPSISense*Nutstress4GrossResp*GrainFillDowreg_brch(MainBranchNum_pft(NZ),NZ)
  
  !Reserve medium-root maintenance before other roots consume the shared pool.
  !The demand denominator includes fine and active primary roots in this layer.
  call ReserveMediumRootMaintenance(N,L,NZ,TFN6_vr,fRootGrowPSISense,Root1stSinkTip, &
    FoundRootAxesTip,Root1stTipUpdateFlag,MediumMediumRootMaintCAlloc_OUltd,MediumMediumRootMaintCAlloc_Oltd)

  !     FOR EACH ROOT AXIS
  !=================================================================================
  !loop through all root axes 
  D5050: DO NR=1,NumStructuralRootAxes_pft(NZ)
    
    !     SECONDARY ROOT EXTENSION
    !   make sure current layer is no deeper than the lowest root layer
    IF(L.LE.NRoot1stTipLay_raxes(NR,NZ) .AND. .not.FoundRootAxesTip(N,NR))THEN

      call SumRootBiome(yearIJ,NZ,mass_inital)
      !
      !secondary roots can grow in any root layer that is not the deepest
      call Grow2ndRootAxes(yearIJ%I,yearIJ%J,L,NZ,N,NR,VMXCSpec,CNRTW,CPRTW,fdLext2nd,RespElongWatSens,&
        DMRespEff,fRootGrowPSISense,TFN6_vr,RootSinkC_vr,Root2ndSink_pvr,litrflx2,RCO2flx2,Root2ndStrutRemob)

!      dlitrfall=plt_bgcr%LitrfallBlgrElms_pft(:,NZ)
      call SumRootBiome(yearIJ,NZ,mass_finale)
      !      call SumRootAR(NZ);call SumLitfallBlg(NZ);dlitrfall=dlitrfall-plt_bgcr%LitrfallBlgrElms_pft(:,NZ)
      if(.false.)write(426,*)'2ndroot',yearIJ%I*1000+yearIJ%J/24.,NR*10+L,mass_finale(ielmn)-mass_inital(ielmn) + &
        litrflx2(ielmn),'afnod',mass_finale(ielmn),mass_inital(ielmn),dlitrfall(ielmn),litrflx2(ielmn),plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ)

      !     call SumRootBiome(yearIJ,NZ,mass_inital)
      !     call SumRootAR(NZ);call SumLitfallBlg(NZ)
      !     tmpval=-plt_bgcr%RootAutoCO2_pft(NZ)+plt_bgcr%LitrfallBlgrElms_pft(ielmc,NZ)

      litrflx = litrflx+litrflx2
      RCO2flx = RCO2flx+RCO2flx2
      
      !
      !make sure it is plant root, and root axis NR is below layer L-1.      
      !
      FoundPrimaryRootsLayer=N.EQ.ipltroot .and. Root1stAxesTipDepz2Surf_pft(NR,NZ).GT.CumSoilThickness_vr(L-1)
      if(.false.)write(426,*)'structroootttttttt'
      call SumRootBiome(yearIJ,NZ,mass_inital)      
      IF(FoundPrimaryRootsLayer)THEN 
        if(lcoarseroot .and. is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))THEN
          call GrowWoodyStructuralRootAxes(yearIJ,N,NR,L,Lnext,NZ,VMXCSpec,Root1stSink_pvr,Root2ndSink_pvr,&
            Root1stSinkTip,RootSinkC_vr,fRootGrowPSISense,TFN6_vr,fdLext1st,RespElongWatSens,DMRespEff,CNRTW,CPRTW,&
            Root1stTipUpdateFlag,FoundRootAxesTip,litrflxt)
        else
          ! plant PRIMARY ROOT EXTENSION of fibrous roots
          call GrowFibrous1stRootAxes(yearIJ,N,NR,L,Lnext,NZ,VMXCSpec,Root1stSink_pvr,&
            Root2ndSink_pvr,Root1stSinkTip,RootSinkC_vr,fRootGrowPSISense,TFN6_vr,fdLext1st,RespElongWatSens,&
            DMRespEff,CNRTW,CPRTW,Root1stTipUpdateFlag,FoundRootAxesTip,litrflxt,RCO2flxt)
        endif   
        
        litrflx = litrflx+AZMAX1(litrflxt)
        RCO2flx = RCO2flx+AZMAX1(RCO2flxt)
      ENDIF

      dlitrfall=plt_bgcr%LitrfallBlgrElms_pft(:,NZ)
      call SumRootBiome(yearIJ,NZ,mass_finale)
 !     if(yearIJ%I>=133)write(*,*)yearIJ%I*1000+yearIJ%J,mass_finale(ielmn),mass_inital(ielmn),litrflxt(ielmn)
 !     dmassN=mass_finale(ielmn)-mass_inital(ielmn)+litrflxt(ielmn)
      !      call SumRootAR(NZ);call SumLitfallBlg(NZ); dlitrfall=dlitrfall-plt_bgcr%LitrfallBlgrElms_pft(:,NZ)
      if(.false.)write(426,*)'1stroot',yearIJ%I*1000+yearIJ%J/24.,NR*10+L,mass_finale(ielmn)-mass_inital(ielmn)+ &
        litrflxt(ielmn),'afnod',mass_finale(ielmn),mass_inital(ielmn),dlitrfall(ielmn),litrflxt(ielmn)
!      if(abs(dmassN)>1.e-8)stop  
      call SumRootBiome(yearIJ,NZ,mass_inital)

      if(lcoarseroot .and. is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)) .and. &
        Root1stLenPP_rpvr(L,NR,NZ).GT.0._r8 .and. N.eq.ipltroot)then    
        !
        call GrowMediumRootAxes(yearIJ,N,L,NR,NZ,VMXCSpec,RootSinkC_vr,RootMSink_pvr,fUSE4MR_pvr,&
          TFN6_vr,fRootGrowPSISense,DMRespEff,RespElongWatSens,CNRTW,CPRTW,fdLext1st, &
          MediumMediumRootMaintCAlloc_OUltd(NR),MediumMediumRootMaintCAlloc_Oltd(NR),litrflx2)
        MediumMediumRootMaintCAlloc_Oltd(NR)=0._r8
      endif  
      dlitrfall=plt_bgcr%LitrfallBlgrElms_pft(:,NZ)
      call SumRootBiome(yearIJ,NZ,mass_finale)
      !      call SumRootAR(NZ);call SumLitfallBlg(NZ);dlitrfall=dlitrfall-plt_bgcr%LitrfallBlgrElms_pft(:,NZ)
      if(.false.)write(426,*)'medroot',yearIJ%I*1000+yearIJ%J/24.,NR*10+L,mass_finale(ielmn)-mass_inital(ielmn) + &
        litrflx2(ielmn),'afnod',mass_finale(ielmn),mass_inital(ielmn),dlitrfall(ielmn),litrflx2(ielmn),&
        plt_bgcr%LitrfallBlgrElms_pft(ielmn,NZ)

    ENDIF
    NMaxRootBotLayer_pft(NZ)=MAX(NMaxRootBotLayer_pft(NZ),NRoot1stTipLay_raxes(NR,NZ))
    !   call SumRootBiome(NZ,mass_finale,massr1st,massr2nd,massnonst,massnodul)

    if(plt_biom%RootMycoNonstElms_rpvr(ielmc,N,L,NZ)<-1.e-7_r8)then
      print*,yearIJ%I,yearIJ%J,'=',plt_biom%RootMycoNonstElms_rpvr(ielmc,N,L,NZ),N,L
      call endrun(trim(mod_filename)//' at line',__LINE__)          
    endif  
  ENDDO D5050
  !Return reservations for groups skipped after geometry changed during growth.
  plt_biom%RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = &
    plt_biom%RootMycoNonstElms_rpvr(ielmc,N,L,NZ)+SUM(MediumMediumRootMaintCAlloc_Oltd)
  if(.false.)write(426,*)'-------------------------------'
  call PrintInfo('end '//subname)
  end associate
  end subroutine GrowRootMycoAxes
!----------------------------------------------------------------------------------------------------
  subroutine GetPlantRoot1stDepz(yearIJ,NZ)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='GetPlantRoot1stDepz'
  INTEGER :: N,NR,L
  associate(                                                                &
    ZERO                         => plt_site%ZERO                          ,& !input  :threshold zero for numerical stability, [-]  
    NU                           => plt_site%NU                            ,& !input  :current soil surface layer number, [-]
    DLYR3                        => plt_site%DLYR3                         ,& !input  :vertical thickness of soil layer, [m]    
    MaxNumRootLays               => plt_site%MaxNumRootLays                ,& !input  :maximum root layer number,[-]    
    NumStructuralRootAxes_pft    => plt_morph%NumStructuralRootAxes_pft    ,& !input  :number of structural root axes,[-]    
    Root1stDepz_raxes            => plt_morph%Root1stDepz_raxes            ,& !input  :root layer depth, [m]    
    CumSoilThickness_vr          => plt_site%CumSoilThickness_vr           ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]    
    SoilBulkDensity_vr           => plt_soilchem%SoilBulkDensity_vr        ,& !input  :soil bulk density, [Mg m-3]
    Root1stLenLoc_rpvr           => plt_morph%Root1stLenLoc_rpvr           ,& !inoput :local structrual root length in layer, [m]    
    Root1stAxesTipDepz2Surf_pft  => plt_morph%Root1stAxesTipDepz2Surf_pft   & !output :plant primary depth relative to column surface, [m] 
  )  

  call PrintInfo('beg' //subname)
  N=ipltroot
  Root1stAxesTipDepz2Surf_pft=0._r8
  Root1stLenLoc_rpvr = 0._r8
  DO NR=1,NumStructuralRootAxes_pft(NZ)
    DO L=NU,MaxNumRootLays+1
      !root tip layer
      IF(Root1stDepz_raxes(NR,NZ).GT.CumSoilThickness_vr(L-1) .AND. Root1stDepz_raxes(NR,NZ).LE.CumSoilThickness_vr(L))THEN
        IF(SoilBulkDensity_vr(L).GT.ZERO)THEN
          !in soil
          Root1stAxesTipDepz2Surf_pft(NR,NZ)=Root1stDepz_raxes(NR,NZ)-CumSoilThickness_vr(0)
        ELSE
          !in ponding water
          Root1stAxesTipDepz2Surf_pft(NR,NZ)=Root1stDepz_raxes(NR,NZ)
        ENDIF
        Root1stLenLoc_rpvr(L,NR,NZ)=Root1stDepz_raxes(NR,NZ)-CumSoilThickness_vr(L-1)
        EXIT
      !root full layer  
      elseif(Root1stDepz_raxes(NR,NZ).GT.CumSoilThickness_vr(L-1))then
        Root1stLenLoc_rpvr(L,NR,NZ)=DLYR3(L)
      ENDIF
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine GetPlantRoot1stDepz

!----------------------------------------------------------------------------------------------------
  subroutine Grow1stRootElongZoneLayL(yearIJ,N,NR,L,Lnext,NZ,VMXCSpec,Root1stSinkTip,RootSinkC_vr,&
    fRootGrowPSISense,TFN6_vr,fdLext1st_mech,RespElongWatSens,DMRespEff,CNRTW,CPRTW,&
    litrflxt,RCO2flxt,Root1stExtenzPP)
  !
  !Do elongation of the primary roots in layer L for root axies NR  
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  integer,  intent(in) :: N,NR,NZ
  integer,  intent(in) :: L,Lnext       !current and next layer for root growth
  REAL(R8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)  
  real(r8), intent(in) :: fRootGrowPSISense
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  real(r8), intent(in) :: fdLext1st_mech       !soil compaction modification of elongation
  real(r8), intent(in) :: RespElongWatSens
  real(r8), intent(in) :: Root1stSinkTip(pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: CNRTW,CPRTW  
  real(r8), intent(in) :: DMRespEff                       !respiraiton efficiency for new dry root mass creation from nonstructural C
  real(r8), intent(in) :: VMXCSpec
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)
  real(r8), intent(inout) :: RCO2flxt       !accumulator of CO2 flux
  real(r8), intent(out) :: Root1stExtenzPP
  real(r8) :: litrflx(NumPlantChemElms)
  real(r8) :: FracRoot1stCSinkL
  real(r8) :: GRTWTM  
  real(r8) :: TFRCO2,FRCO2  
  real(r8) :: Rmaint1st_CO2
  real(r8) :: RNonstCO2_OUltd
  real(r8) :: RGrowCO2_Oltd,RGrowCO2_OUltd               !growth respiration unltd,ltd by O2 and unlimited by N,P
  real(r8) :: RNonstCO2_Oltd
  real(r8) :: RCO2PotGroth_OUltd,RCO2PotGroth_Oltd       !O2-unlimited/limited potential growth respiration
  real(r8) :: DMRTR,ZPOOLB,PPOOLB,FNP
  real(r8) :: RootMycoNonst4GrowC_OUltd
  real(r8) :: RootMycoNonst4Grow_Oltd(NumPlantChemElms)  !nonstructural N,P unltd, not ltd by O2 used in growth
  real(r8) :: RootMycoNonst4Grow_OUltd(NumPlantChemElms) !nonstructural N,P unltd, ltd by O2 used in growth
  real(r8) :: RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd
  real(r8) :: RCO2T1st_OUltd,RCO2T1st_Oltd
  real(r8) :: dsenecE
  real(r8) :: RootMycoNonst4GrowC_Oltd
  real(r8) :: Root1stcylc(NumPlantChemElms)              !nonstrucal biomass from remobilization
  real(r8) :: Root1stNetGrowthElms(NumPlantChemElms)
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: RCO2MDef1stStructPayoff_Oltd,RCO2MDef1stStructPayoff_OUltd
  real(r8) :: mass_finale(NumPlantChemElms),RLEN,RootNsink1st
  integer  :: NE,M,LL,LL1,L1,N1
  real(r8) :: RRadius,ratio,rootmy1
  real(r8) :: Root1stPotExtzPP             !per axial pontential/actual root extension due to nonstrucal C mobilization.  
  real(r8), parameter :: RTDPX1=0.005_r8  !the apical growing zone length
  character(len=*), parameter :: subname='Grow1stRootElongZoneLayL'

  associate(                                                                          &
    Myco_pft                          => plt_morph%Myco_pft                          ,& !input  :mycorrhizal type (no or yes),[-]    
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft              ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    RAutoRootO2Limter_rpvr            => plt_rbgc%RAutoRootO2Limter_rpvr             ,& !input  :O2 constraint to root respiration (0-1), [-]
    MaxNumRootLays                    => plt_site%MaxNumRootLays                     ,& !input  :maximum root layer number,[-]
    CumSoilThickness_vr               => plt_site%CumSoilThickness_vr                ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    ZERO4Groth_pft                    => plt_biom%ZERO4Groth_pft                     ,& !input  :threshold zero for plang growth calculation, [-]
    CNRTS_pft                         => plt_allom%CNRTS_pft                         ,& !input  :root N:C ratio x root growth yield, [-]
    CPRTS_pft                         => plt_allom%CPRTS_pft                         ,& !input  :root P:C ratio x root growth yield, [-]
    RootBiomGrosYld_pft               => plt_allom%RootBiomGrosYld_pft               ,& !input  :root growth yield, [g g-1]
    DLYR3                             => plt_site%DLYR3                              ,& !input  :vertical thickness of soil layer, [m]        
    Root1stRadius_pvr                 => plt_morph%Root1stRadius_pvr                 ,& !input  :root layer radius primary axes, [m]  
    PlantElmAllocMat4Litr             => plt_soilchem%PlantElmAllocMat4Litr          ,& !input  :litter kinetic fraction, [-]    
    StalkBiomGrowthYld_pft            => plt_allom%StalkBiomGrowthYld_pft            ,& !input  :stalk growth yield, [gC gC-1]    
    Root1stMaxRadius1_pft             => plt_morph%Root1stMaxRadius1_pft             ,& !input  :root radius primary axes, [m]    
    RootMatureAge_pft                 => plt_morph%RootMatureAge_pft                 ,& !input : Root maturation age, [h]
    PSIRoot_pvr                       => plt_ew%PSIRoot_pvr                          ,& !input :root total water potential, [Mpa]    
    RootAge_rpvr                      => plt_morph%RootAge_rpvr                      ,& !inoput :root age,[h]
    LitrfallElms_pvr                  => plt_bgcr%LitrfallElms_pvr                   ,& !input  :plant LitrFall element, [g d-2 h-1]    
    PlantPopuLive_pft                 => plt_site%PlantPopuLive_pft                  ,& !input  :plant population, [d-2]
    iPlant2ndGrothPattern_pft         => plt_pheno%iPlant2ndGrothPattern_pft         ,& !input  :plant expression of secondary growth, [-]                
    Num1stAxesPerStructRootXPOP_pft   => plt_morph%Num1stAxesPerStructRootXPOP_pft   ,& !output :population primary root axes number on one structrual axis, [d-2]        
    Root1stSpecLen_pft                => plt_morph%Root1stSpecLen_pft                ,& !input  :specific root length primary axes, [m root gC-1]
    k_fine_comp                       => pltpar%k_fine_comp                          ,& !input  :fine litter complex id    
    iroot                             => pltpar%iroot                                ,& !input  :group id of plant root litter    
    Root1stLenPP_rpvr                 => plt_morph%Root1stLenPP_rpvr                 ,& !input  :primary root axis length in soil layer, [m d-2]
    SeedDepth_pft                     => plt_morph%SeedDepth_pft                     ,& !input  :seeding depth, [m]
    flag2ndGrowth_pvr                 => plt_morph%flag2ndGrowth_pvr                 ,& !inoput  :flag for secondary growth of primary roots, [-]        
    NGTopRootLayer_pft                => plt_morph%NGTopRootLayer_pft                ,& !input  :soil layer at planting depth, [-]
    Root1stDepz_raxes                 => plt_morph%Root1stDepz_raxes                 ,& !input  :root layer depth, [m]
    Root1stMaxRadius_pft              => plt_morph%Root1stMaxRadius_pft              ,& !input  :maximum radius of primary roots, [m]    
    FracRootElmAllocm                 => plt_allom%FracRootElmAllocm                 ,& !input  :C woody fraction in root,[-]    
    Root1stAxesTipDepz2Surf_pft       => plt_morph%Root1stAxesTipDepz2Surf_pft       ,& !output :plant primary depth relative to column surface, [m]     
    iPlantRootProfile_pft             => plt_pheno%iPlantRootProfile_pft             ,& !input  :plant growth type (vascular, non-vascular),[-]
    GrainFillDowreg_brch              => plt_photo%GrainFillDowreg_brch              ,& !input  :grain fill down-regulation of annual plants, [-]
    RootMaintDef_CO2_pvr              => plt_bgcr%RootMaintDef_CO2_pvr               ,& !inoput :plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]    
    RootCO2Autor_pvr                  => plt_rbgc%RootCO2Autor_pvr                   ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    RootMyco1stSinkC_rpvr             => plt_rbgc%RootMyco1stSinkC_rpvr              ,& !output :1st root C sink, [gC d-2 h-1]
    RootCO2EmisPot_pvr                => plt_rbgc%RootCO2EmisPot_pvr                 ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootRespPotent_pvr                => plt_rbgc%RootRespPotent_pvr                 ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]
    RootMycoNonstElms_rpvr            => plt_biom%RootMycoNonstElms_rpvr             ,& !inoput :root layer nonstructural element, [g d-2]
    RootMyco1stElm_raxs               => plt_biom%RootMyco1stElm_raxs                ,& !inoput :root layer primary axes structrual element, [g d-2]
    Root1stActStructElms_rpvr         => plt_biom%Root1stActStructElms_rpvr          ,& !inoput :Root layer primary axes Active zone structrual element, [g d-2]
    RootMyco2ndStrutElms_rpvr         => plt_biom%RootMyco2ndStrutElms_rpvr          ,& !inoput :root layer element secondary axes, [g d-2]
    RootMyco1stStrutElms_rpvr         => plt_biom%RootMyco1stStrutElms_rpvr          ,& !inoput :root layer element primary axes, [g d-2]
    NRoot1stTipLay_raxes              => plt_morph%NRoot1stTipLay_raxes               & !inoput :soil layer number for deepest root axes, [-]
  )
  call PrintInfo('beg '//subname)
  !
  !     FRACTION OF PRIMARY ROOT SINK IN SOIL LAYER
  !     ATTRIBUTED TO CURRENT AXIS
  !
  !     Root1stSink_pvr=primary root sink strength
  !     RootSinkC_vr=total root sink strength
  !     FracRoot1stCSinkL=fraction of primary root sink strength in axis
  !
  litrflx     = 0._r8
  !  mass_inital = 0._r8
  !  DO L1  = max(1,L-1), Lnext
  !    DO NE=1,NumPlantChemElms
  !      mass_inital(NE) = mass_inital(NE) +RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
  !    ENDDO  
  !    DO N1=1,Myco_pft(NZ)
  !      DO NE=1,NumPlantChemElms
  !        mass_inital(NE) = mass_inital(NE) +RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
  !      ENDDO
  !    ENDDO
  !  ENDDO
  !
  IF(RootSinkC_vr(N,L).GT.ZERO4Groth_pft(NZ))THEN
    FracRoot1stCSinkL=Root1stSinkTip(NR)/RootSinkC_vr(N,L)
  ELSE
    FracRoot1stCSinkL=1.0_r8
  ENDIF
  !
  !     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
  !     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
  !
  !     N,P CONSTRAINT ON PRIMARY ROOT RESPIRATION FROM
  !     NON-STRUCTURAL C:N:P
  !
  !     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
  !     
  !     CNKI,CPKI=nonstructural N,P inhibition constant on growth
  !
  !     PRIMARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
  !     ROOT STRUCTURAL N
  !
  !     Rmaint1st_CO2=root maintenance respiration
  !     RmSpecPlant=specific maintenance respiration rate (g C g-1 N h-1)
  !     WTRT1N=primary root N mass
  !     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
  !     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
  !     fRootGrowPSISense=growth function of root water potential
  !
  if(lcoarseroot .and. is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))then
    Rmaint1st_CO2 = AZMAX1(RmSpecPlant*Root1stActStructElms_rpvr(ielmn,L,NR,NZ))*TFN6_vr(L)
  else
    Rmaint1st_CO2 = AZMAX1(RmSpecPlant*RootMyco1stElm_raxs(ielmn,NR,NZ))*TFN6_vr(L)
  endif

  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ)) .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    Rmaint1st_CO2=Rmaint1st_CO2*fRootGrowPSISense
  ENDIF
  !
  !     O2-UNLIMITED PRIMARY ROOT RESPIRATION FROM ROOT NON-STRUCTURAL C
  !     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
  !
  !     RNonstCO2_OUltd=respiration from non-structural C unlimited by O2
  !     VMXC=rate constant for nonstructural C oxidation in respiration C     
  !     FracRoot1stCSinkL=fraction of primary root sink strength in axis
  !     CPOOL=non-structural C mass
  !     GrainFillDowreg_brch=termination feedback inhibition on C3 CO2
  !     fRootGrowPSISense=growth function of root water potential
  !
  RNonstCO2_OUltd=AZMAX1(FracRoot1stCSinkL*RootMycoNonstElms_rpvr(ielmc,N,L,NZ))*VMXCSpec 

  !do not let the root grow beyond the maximum allowed layer
  !This is a temporary fix for the lack of deep soil toxicty that prevents root growth
  IF(Root1stAxesTipDepz2Surf_pft(NR,NZ).GE.CumSoilThickness_vr(MaxNumRootLays))THEN
    RNonstCO2_OUltd=AMIN1(Rmaint1st_CO2,RNonstCO2_OUltd)
  ENDIF
  !
  !     O2-LIMITED PRIMARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
  !
  !     RNonstCO2_Oltd=respiration from non-structural C limited by O2
  !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
  !     RespElongWatSens=respiration function of root water potential
  !
  RNonstCO2_Oltd               = RNonstCO2_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
  RCO2PotGroth_OUltd           = RNonstCO2_OUltd-Rmaint1st_CO2
  RCO2PotGroth_Oltd            = RNonstCO2_Oltd-Rmaint1st_CO2
  RGrowCO2_OUltd               = AZMAX1(RCO2PotGroth_OUltd)*RespElongWatSens
  RGrowCO2_Oltd                = AZMAX1(RCO2PotGroth_Oltd)*RespElongWatSens
  RootMaintDef_CO2_pvr(N,L,NZ) = RootMaintDef_CO2_pvr(N,L,NZ)+AMIN1(RCO2PotGroth_Oltd,0._r8)
  !
  !     PRIMARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
  !     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
  !
  !     FracRoot1stCSinkL=fraction of secondary root sink strength in axis
  !     ZPOOLR,PPOOLR=non-structural N,P mass in root 
  !     CNRTS_pft,CPRTS_pft=N,P root growth yield = = [gN,P root/gC root]*[gC root /gC nonst]=[gN,P root/gC nonstC]
  !     FNP=growth respiration limited by non-structural N,P
  !     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
  !
  ZPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))*FracRoot1stCSinkL
  PPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))*FracRoot1stCSinkL
  FNP    = AMIN1(ZPOOLB/CNRTS_pft(NZ),PPOOLB/CPRTS_pft(NZ))*DMRespEff

  IF(RGrowCO2_OUltd.GT.0.0_r8)THEN
    RGrowCO2_OUltd=AMIN1(RGrowCO2_OUltd,FNP)
  ELSE
    RGrowCO2_OUltd=0._r8
  ENDIF

  IF(RGrowCO2_Oltd.GT.0.0_r8)THEN
    RGrowCO2_Oltd=AMIN1(RGrowCO2_Oltd,FNP*RAutoRootO2Limter_rpvr(N,L,NZ))
  ELSE
    RGrowCO2_Oltd=0._r8
  ENDIF
  !
  !     TOTAL NON-STRUCTURAL C,N,P USED IN PRIMARY ROOT GROWTH
  !     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
  !     ENTERED IN 'READQ'
  !
  !     RootMycoNonst4GrowC_OUltd,RootMycoNonst4GrowC_Oltd=total non-structural C used in growth and respn unltd,ltd by O2
  !     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
  !     RootMycoNonst4Grow_OUltd(ielmc),RootMycoNonst4Grow_Oltd(ielmc)=root C growth unltd,ltd by O2
  !     RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
  !
  RootMycoNonst4GrowC_OUltd       = RGrowCO2_OUltd/DMRespEff
  RootMycoNonst4GrowC_Oltd        = RGrowCO2_Oltd/DMRespEff

  RootMycoNonst4Grow_OUltd(ielmc) = RootMycoNonst4GrowC_OUltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_OUltd(ielmn) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CNRTW)
  RootMycoNonst4Grow_OUltd(ielmp) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CPRTW)
  !
  !nonstrucal biomass elements for metabolism
  RootMycoNonst4Grow_Oltd(ielmc) = RootMycoNonst4GrowC_Oltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_Oltd(ielmn) = AZMAX1(AMIN1(FracRoot1stCSinkL*RootMycoNonstElms_rpvr(ielmn,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*CNRTW))
  RootMycoNonst4Grow_Oltd(ielmp) = AZMAX1(AMIN1(FracRoot1stCSinkL*RootMycoNonstElms_rpvr(ielmp,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*CPRTW))
  RootMycoNonst4Grow_Oltd        = RootMycoNonst4Grow_Oltd*fdLext1st_mech
  RCO2Nonst4Nassim_OUltd         = AZMAX1(1.70_r8*RootMycoNonst4Grow_OUltd(ielmn))
  RCO2Nonst4Nassim_Oltd          = AZMAX1(1.70_r8*RootMycoNonst4Grow_Oltd(ielmn))
  !
  !     TOTAL PRIMARY ROOT RESPIRATION
  !
  !     RCO2T1st_OUltd,RCO2T1st_Oltd=total C respiration unltd,ltd by O2
  !     Rmaint1st_CO2=root maintenance respiration
  !     RNonstCO2_OUltd,RNonstCO2_Oltd=respiration from non-structural C unltd,ltd by O2
  !     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
  !     RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
  !     RootCO2Autor_pvr=total root respiration
  !     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
  !
  RCO2T1st_OUltd = AMIN1(Rmaint1st_CO2,RNonstCO2_OUltd)+RGrowCO2_OUltd+RCO2Nonst4Nassim_OUltd
  RCO2T1st_Oltd  = AMIN1(Rmaint1st_CO2,RNonstCO2_Oltd)+RGrowCO2_Oltd+RCO2Nonst4Nassim_Oltd
  !
  !subtract CO2 from reserve 
  RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-RCO2T1st_Oltd

  DO NE=1,NumPlantChemElms
    RootMycoNonst4Grow_Oltd(NE)       = AMIN1(RootMycoNonstElms_rpvr(NE,N,L,NZ),RootMycoNonst4Grow_Oltd(NE))
    RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)-RootMycoNonst4Grow_Oltd(NE))
  ENDDO

  !
  !mean primary root extension based on growth flux, accounting for all root axes
  Root1stPotExtzPP = RootMycoNonst4Grow_Oltd(ielmc)*FracRootElmAllocm(ielmc,k_fine_comp) &
    /Num1stAxesPerStructRootXPOP_pft(NZ)*Root1stSpecLen_pft(N,NZ)
  Root1stExtenzPP  = Root1stPotExtzPP
  
  !dealing with negative growth

  Root1stNetGrowthElms(:)=RootMycoNonst4Grow_Oltd(:) ! > 0

  RootMyco1stSinkC_rpvr(L,NR,NZ)=RCO2T1st_Oltd+RootMycoNonst4Grow_Oltd(ielmc)

!  write(981,*)yearIJ%I*1000+yearIJ%J/24.,RCO2PotGroth_Oltd,Root1stPotExtzPP,NR,L,fdLext1st_mech,&
!    Root1stAxesTipDepz2Surf_pft(NR,NZ).GE.CumSoilThickness_vr(MaxNumRootLays),Root1stDepz_raxes(NR,NZ),&
!    CumSoilThickness_vr(MaxNumRootLays),Root1stAxesTipDepz2Surf_pft(NR,NZ)

!  mass_finale=0._r8
!  DO L1=max(1,L-1),Lnext
!    DO NE=1,NumPlantChemElms
!      mass_finale(NE)=mass_finale(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!    ENDDO    
!    DO N1=1,Myco_pft(NZ)
!      DO NE=1,NumPlantChemElms
!        mass_finale(NE)=mass_finale(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  write(490,*)yearIJ%I*1000+yearIJ%J/24.,'bfremob',mass_finale(ielmn)-mass_inital(ielmn),Root1stNetGrowthElms(ielmn),RCO2PotGroth_Oltd

  if(-RCO2PotGroth_Oltd.GT.0._r8)then !Root1stNetGrowthElm and Root1stPotExtzPP are zero   
    !may turn Root1stNetGrowthElms negative due to remobilzation
    
    call RemobilizePrimeActRoots(N,L,NZ,NR,RCO2PotGroth_Oltd,RCO2PotGroth_OUltd,Root1stcylc,RCO2T1st_OUltd,&
      RCO2T1st_Oltd,Root1stNetGrowthElms,litrflx,RootTipZone=.True.)
    
    DO NE=1,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)+Root1stcylc(NE))
    ENDDO         
  endif

!  mass_finale=0._r8
!  DO L1=max(1,L-1),Lnext
!    DO NE=1,NumPlantChemElms
!      mass_finale(NE)=mass_finale(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!    ENDDO    
!    DO N1=1,Myco_pft(NZ)
!      DO NE=1,NumPlantChemElms
!        mass_finale(NE)=mass_finale(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  write(490,*)yearIJ%I*1000+yearIJ%J/24.,'afremob',mass_finale(ielmn)-mass_inital(ielmn)+litrflx(ielmn)+Root1stNetGrowthElms(ielmn),Root1stNetGrowthElms(ielmn),litrflx(ielmn)

  !
  !     ALLOCATE PRIMARY ROOT TOTAL RESPIRATION TO ALL SOIL LAYERS
  !     THROUGH WHICH PRIMARY ROOTS GROW
  !
  !     Root1stLenPP_rpvr=primary root length, assuming primary roots are of the same radius
  !     SeedDepth_pft=seeding depth
  !     FRCO2=fraction of primary root respiration attributed to layer
  !     RootCO2Autor_pvr=total root respiration
  !     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
  !
  IF(Root1stDepz_raxes(NR,NZ).GT.CumSoilThickness_vr(NGTopRootLayer_pft(NZ)))THEN
    TFRCO2=0._r8;RLEN=0._r8
    LL1=NGTopRootLayer_pft(NZ)

    ! Exclude the mature prefix; retain the tip as a fallback.
    if(is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))then
      DO WHILE(LL1.LT.NRoot1stTipLay_raxes(NR,NZ))
        IF(.not.flag2ndGrowth_pvr(LL1,NR,NZ))exit
        LL1=LL1+1
      ENDDO
    endif

    ! Normalize using exactly the layers that receive this respiration.
    RLEN=SUM(MAX(0._r8,Root1stLenPP_rpvr(LL1:NRoot1stTipLay_raxes(NR,NZ),NR,NZ)))
    D5100: DO LL=LL1,NRoot1stTipLay_raxes(NR,NZ)
      IF(LL.EQ.NRoot1stTipLay_raxes(NR,NZ))THEN
        FRCO2=MAX(0._r8,1._r8-TFRCO2)
      ELSEIF(RLEN.GT.0._r8)THEN
        FRCO2=MIN(MAX(0._r8,1._r8-TFRCO2),MAX(0._r8,Root1stLenPP_rpvr(LL,NR,NZ))/RLEN)
      ELSE
        FRCO2=0._r8
      ENDIF
      TFRCO2                      = TFRCO2+FRCO2
      RootRespPotent_pvr(N,LL,NZ) = RootRespPotent_pvr(N,LL,NZ)+RCO2T1st_OUltd*FRCO2
      RootCO2EmisPot_pvr(N,LL,NZ) = RootCO2EmisPot_pvr(N,LL,NZ)+RCO2T1st_Oltd*FRCO2
      RootCO2Autor_pvr(N,LL,NZ)   = RootCO2Autor_pvr(N,LL,NZ)-RCO2T1st_Oltd*FRCO2
      RCO2flxt                    = RCO2flxt-RCO2T1st_Oltd*FRCO2
    ENDDO D5100
  ELSE
    RootRespPotent_pvr(N,L,NZ) = RootRespPotent_pvr(N,L,NZ)+RCO2T1st_OUltd
    RootCO2EmisPot_pvr(N,L,NZ) = RootCO2EmisPot_pvr(N,L,NZ)+RCO2T1st_Oltd
    RootCO2Autor_pvr(N,L,NZ)   = RootCO2Autor_pvr(N,L,NZ)-RCO2T1st_Oltd
    RCO2flxt                   = RCO2flxt-RCO2T1st_Oltd
  ENDIF
  !
  !     ALLOCATE ANY NEGATIVE PRIMARY ROOT C,N,P GROWTH TO SECONDARY
  !     ROOTS ON THE SAME AXIS IN THE SAME LAYER UNTIL SECONDARY ROOTS
  !     HAVE DISAPPEARED
  !
  !     Frac2Senes1=fraction of primary root C to be remobilized
  !     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
  !     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
  !     Root2ndLen_rpvr=secondary root length
  !
  !if primary root growth is negative, then offset negative growth by killing secondary roots
  IF(Root1stNetGrowthElms(ielmc).LT.0.0_r8)THEN              
    call Withdraw2ndRoots(N,NZ,L,NR,Root1stNetGrowthElms,litrflx)
  ENDIF

!  mass_finale=0._r8
!  DO L1=max(1,L-1),Lnext
!    DO NE=1,NumPlantChemElms
!      mass_finale(NE)=mass_finale(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!    ENDDO    
!    DO N1=1,Myco_pft(NZ)
!      DO NE=1,NumPlantChemElms
!        mass_finale(NE)=mass_finale(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  write(490,*)yearIJ%I*1000+yearIJ%J/24.,'af2ndrot',mass_finale(ielmn)-mass_inital(ielmn)+litrflx(ielmn)+Root1stNetGrowthElms(ielmn),Root1stNetGrowthElms(ielmn),'litr',litrflx(ielmn)

  IF(Root1stNetGrowthElms(ielmc).LT.0.0_r8)THEN
    !negative growth
    IF(RootMyco1stElm_raxs(ielmc,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
      !primary roots withdraw, note that primary root depth was initialized at seedDepth
      !make sure it fits the corase root plants     
      Root1stExtenzPP=(Root1stDepz_raxes(NR,NZ)-SeedDepth_pft(NZ))/RootMyco1stElm_raxs(ielmc,NR,NZ)*Root1stNetGrowthElms(ielmc) 
    ELSE
      !enforce mass balance under a special edge case.
      DO M=1,jsken      
        DO NE=1,NumPlantChemElms
          LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) &
            +PlantElmAllocMat4Litr(NE,iroot,M,NZ)*Root1stNetGrowthElms(NE)
        ENDDO  
      ENDDO
      Root1stNetGrowthElms=0._r8
    ENDIF
  ENDIF
  !
  !the extension should not exceed next soil layer thickness  
  !this means: 
  !1. the extension is vertical, should the extension be per axial?
  !2. the extension may exceed current layer thickness, but not the next layer
  IF(L.LT.MaxNumRootLays)THEN
    Root1stExtenzPP=AMIN1(DLYR3(Lnext),Root1stExtenzPP)
  ENDIF
  
  !the extension can be negative or positive
  call PrimeRootsExtension(yearIJ,L,Lnext,N,NR,NZ,Root1stNetGrowthElms,Root1stExtenzPP)

!  mass_finale=0._r8
!  DO L1=max(1,L-1),Lnext
!    DO NE=1,NumPlantChemElms
!      mass_finale(NE)=mass_finale(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!    ENDDO    
!    DO N1=1,Myco_pft(NZ)
!      DO NE=1,NumPlantChemElms
!        mass_finale(NE)=mass_finale(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
  litrflxt=litrflxt+litrflx
!  if(I>=143)write(433,*)'elong',mass_finale(ielmc)-mass_inital(ielmc)+RCO2T1st_Oltd+litrflxt(ielmc),RCO2T1st_Oltd,litrflxt(ielmc),'RCO2flxt',RCO2flxt
!  write(490,*)'elong',mass_finale(ielmn)-mass_inital(ielmn)+litrflxt(ielmn),mass_finale(ielmn),mass_inital(ielmn),litrflxt(ielmn)
  call PrintInfo('end '//subname)
  end associate
  end subroutine Grow1stRootElongZoneLayL
!----------------------------------------------------------------------------------------------------
  subroutine GrowWoodyStructuralRootAxes(yearIJ,N,NR,L,Lnext,NZ,VMXCSpec,Root1stSink_pvr,Root2ndSink_pvr,&
    Root1stSinkTip,RootSinkC_vr,fRootGrowPSISense,TFN6_vr,fdLext1st,RespElongWatSens,DMRespEff,CNRTW,CPRTW,&
    Root1stTipUpdateFlag,FoundRootAxesTip,litrflxt)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: N,NR
  integer, intent(in) :: L,Lnext
  INTEGER, intent(in) :: NZ
  real(r8), intent(in) :: fRootGrowPSISense
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  real(r8), intent(in) :: fdLext1st,RespElongWatSens
  real(r8), intent(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: Root1stSinkTip(pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: CNRTW,CPRTW  
  real(r8), intent(in) :: DMRespEff              !respiraiton efficiency for fresh root growth
  real(r8), intent(in) :: VMXCSpec
  REAL(R8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)  
  logical,  intent(inout) :: Root1stTipUpdateFlag(pltpar%MaxNumRootAxes)  
  logical,  intent(inout) :: FoundRootAxesTip(pltpar%jroots,pltpar%MaxNumRootAxes)    
  real(r8), intent(out) :: litrflxt(NumPlantChemElms)
  character(len=*), parameter :: subname='GrowWoodyStructuralRootAxes'
  logical  :: FoundRootTipLayer
  logical  :: LDoneTip

  real(r8) :: RCO2flxt,OldRootDepth

  associate(                                                                        &
    MaxNumRootLays                   => plt_site%MaxNumRootLays                    ,& !input  :maximum root layer number,[-]  
    CumSoilThickness_vr              => plt_site%CumSoilThickness_vr               ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]    
    Root1stDepz_raxes               => plt_morph%Root1stDepz_raxes               ,& !input  :primary-root axis depth, [m]
    NRoot1stTipLay_raxes             => plt_morph%NRoot1stTipLay_raxes             ,& !input  :primary-root tip layer
    RootAge_rpvr                     => plt_morph%RootAge_rpvr                     ,& !inoput :root age,[h]    
    flag2ndGrowth_pvr                => plt_morph%flag2ndGrowth_pvr                ,& !input  :flag for secondary growth of primary roots, [-]    
    Root1stLenPP_rpvr                => plt_morph%Root1stLenPP_rpvr                ,& !input  :primary root axis length in soil layer per plant, [m d-2]        
    Num1stAxesPerStructRootXPOP_pft  => plt_morph%Num1stAxesPerStructRootXPOP_pft  ,& !output :population primary root axes number on one structrual axis, [d-2]            
    Root1stXNumL_pvr                 => plt_morph%Root1stXNumL_pvr                 ,& !inoput :soil layer number of primary root axes, [d-2]        
    Root1stAxesTipDepz2Surf_pft      => plt_morph%Root1stAxesTipDepz2Surf_pft      ,& !output :plant primary depth relative to column surface, [m]         
    RootCO2Autor_pvr                 => plt_rbgc%RootCO2Autor_pvr                   & !inoput :root respiration constrained by O2, [g d-2 h-1]
  )
  call PrintInfo('beg '//subname)

  FoundRootTipLayer = Root1stAxesTipDepz2Surf_pft(NR,NZ).LE.CumSoilThickness_vr(L) .OR. L.EQ.MaxNumRootLays
  litrflxt          = 0._r8
  if(.not.Root1stTipUpdateFlag(NR))then !root tip is used for update in the last
    Root1stXNumL_pvr(L,NZ) = Root1stXNumL_pvr(L,NZ)+Num1stAxesPerStructRootXPOP_pft(NZ) 
    RootAge_rpvr(L,NR,NZ)  = RootAge_rpvr(L,NR,NZ)+1._r8

    LDoneTip=.false.
    if(FoundRootTipLayer .and. Root1stSinkTip(NR).GT.0._r8)THEN
      Root1stTipUpdateFlag(NR) = .true.
      LDoneTip                 = .true.
      OldRootDepth=Root1stDepz_raxes(NR,NZ)
      call StructuralRootElongation(yearIJ,N,NR,L,Lnext,NZ,VMXCSpec,Root1stSinkTip,RootSinkC_vr,&
        fRootGrowPSISense,TFN6_vr,fdLext1st,RespElongWatSens,DMRespEff,CNRTW,CPRTW,&
        litrflxt,RCO2flxt)
      if(Root1stDepz_raxes(NR,NZ).LT.OldRootDepth)then
        !L still identifies the old tip layer, even after geometry moves upward.
        call PrimeRootsWithdraw(L,NR,NZ,N,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr)
        FoundRootAxesTip(N,NR)=.true.
      endif
    endif    
      !
      !  
    if(L.LE.NRoot1stTipLay_raxes(NR,NZ) .and. &
       flag2ndGrowth_pvr(L,NR,NZ) .and. Root1stLenPP_rpvr(L,NR,NZ).GT.0._r8)then
      call StructuralRootThickening(yearIJ,N,NR,L,NZ,VMXCSpec,RootSinkC_vr,Root1stSink_pvr, &
        fRootGrowPSISense,TFN6_vr,litrflxt,RCO2flxt,LDoneTip)
    endif
    
  endif

  call PrintInfo('end '//subname)

  end associate
  end subroutine GrowWoodyStructuralRootAxes
!----------------------------------------------------------------------------------------------------
  subroutine StructuralRootThickening(yearIJ,N,NR,L,NZ,VMXCSpec,RootSinkC_vr,Root1stSink_pvr, &
    fRootGrowPSISense,TFN6_vr,litrflxt,RCO2flxt,LDoneTip)
  !
  !Description:
  !Do secondary growth of roots
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: N,NR,L,NZ
  real(r8), intent(in) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  real(r8), intent(in) :: fRootGrowPSISense
  real(r8), intent(in) :: VMXCSpec
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)
  real(r8), intent(inout) :: RCO2flxt
  logical, intent(in) :: LDoneTip
  real(r8) :: DMRespEff
  real(r8) :: FracRoot1stCSinkL
  real(r8) :: Rmaint1st_CO2
  real(r8) :: RNonstCO2_OUltd,RGrowCO2_Oltd,RGrowCO2_OUltd
  real(r8) :: RNonstCO2_Oltd,RCO2PotGroth_OUltd,RCO2PotGroth_Oltd
  real(r8) :: DMRTR,ZPOOLB,PPOOLB,FNP
  real(r8) :: RootMycoNonst4GrowC_Oltd
  real(r8) :: RootMycoNonst4GrowC_OUltd,RootMycoNonst4Lig_OUltd(NumPlantChemElms)
  real(r8) :: RootMycoNonst4Lig_Oltd(NumPlantChemElms)
  real(r8) :: RootMycoNonstC4Thick_Oltd,RootMycoNonstC4Thick_OUltd
  real(r8) :: RootMycoNonstC4Lig_Oltd,RootMycoNonstC4Lig_OUltd     !Nonstructural C used lignification
  real(r8) :: RootMycoNonst4Thick_OUltd(NumPlantChemElms)
  real(r8) :: RootMycoNonst4Thick_Oltd(NumPlantChemElms)   
  real(r8) :: Root1stNetGrowthElms(NumPlantChemElms)
  real(r8) :: Root1stcylc(NumPlantChemElms)              !nonstrucal biomass from remobilization
  real(r8) :: litrflx(NumPlantChemElms)
  real(r8) :: RCO2T1st_OUltd,RCO2T1st_Oltd,RootNsink1st
  real(r8) :: RootThickenWatSens,RespThickenWatSens,dR1stExp,dRPotexp,dRExp,ratio
  real(r8) :: mass_inital(NumPlantChemElms),mass_finale(NumPlantChemElms)
  real(r8) :: RCO2MDef1stStructPayoff_Oltd,RCO2MDef1stStructPayoff_OUltd
  real(r8) :: Rlignif_OUltd, Rlignif_Oltd      !lignification rate of active zone C without/with O limitation
  real(r8) :: RligCO2_OUltd, RligCO2_Oltd      !CO2 respiration for lignification without/with O limitation
  real(r8) :: CO2Qt_lig                        !CO2 quotient of lignification
  real(r8) :: RThickCO2_OUltd, RThickCO2_Oltd  !CO2 for coarse root thickening
  real(r8) :: scal
  real(r8) :: dRCO2MaintDef_Oltd, dRCO2MaintDef_OUltd !CO2 deficit for maintenance after active zone mobilization  
  real(r8) :: alpha_fill                         !C-filling quotient during lignification, [gC reserve /gC lignified-active biomass]
  real(r8) :: mass1(NumPlantChemElms),mass2(NumPlantChemElms),mass3(NumPlantChemElms)
  real(r8) :: RCO2_meta_Oltd,RCO2_meta_OUltd
  real(r8) :: RootMycoDemand,dRThickening
  real(r8) :: f_stress, K_ligmax,errC
  integer  :: NE,L1,N1

  character(len=*), parameter :: subname='SecondaryGrowthZone'

  associate(                                                                          &
    Myco_pft                          => plt_morph%Myco_pft                          ,& !input  :mycorrhizal type (no or yes),[-]    
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft              ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    MainBranchNum_pft                 => plt_morph%MainBranchNum_pft                 ,& !input  :number of main branch,[-]
    RootMyco1stStrutElms_rpvr         => plt_biom%RootMyco1stStrutElms_rpvr          ,& !input  :root layer element primary axes, [g d-2]
    Root1stActStructElms_rpvr         => plt_biom%Root1stActStructElms_rpvr          ,& !inoput :root layer active zone element in primary axes, [g d-2]
    Root1stLigStructElms_rpvr         => plt_biom%Root1stLigStructElms_rpvr          ,& !inoput :root layer lignified zone element in primary axes, [g d-2]
    StalkBiomGrowthYld_pft            => plt_allom%StalkBiomGrowthYld_pft            ,& !input  :stalk growth yield, [gC gC-1]    
    Root1stRadius_pvr                 => plt_morph%Root1stRadius_pvr                 ,& !input  :root layer radius primary axes, [m]  
    rNCRoot_pft                       => plt_allom%rNCRoot_pft                       ,& !input  :root N:C ratio, [gN gC-1]
    Root1stRadius_rpvr                => plt_morph%Root1stRadius_rpvr                ,& !input: root layer radius for each primary axes,  [m]    
    KLigMax_pft                       => plt_biom%KLigMax_pft                        ,& !input  :Maximum lignification rate [h-1]
    KLigMM_pft                        => plt_biom%KLigMM_pft                         ,& !input  :Half saturation parameter for coarse root lignification [h-1]
    CRootActVolPerMassC_pft           => plt_morph%CRootActVolPerMassC_pft           ,& !input  :coarse root active zone volume:mass ratio, [m3 gC-1]    
    CRootLigVolPerMassC_pft           => plt_morph%CRootLigVolPerMassC_pft           ,& !input  :coarse root inactive zone volume:mass ratio, [m3 gC-1]     
    rNCLigRoot_pft                    => plt_allom%rNCLigRoot_pft                    ,& !input  :NC ratio of lignified coarse root, [gN gC-1]
    rPCLigRoot_pft                    => plt_allom%rPCLigRoot_pft                    ,& !input  :PC ratio of lignified coarse root, [gN gC-1]
    Num1stAxesPerStructRootXPOP_pft   => plt_morph%Num1stAxesPerStructRootXPOP_pft   ,& !output :population primary root axes number on one structrual axis, [d-2]        
    RootMyco2ndStrutElms_rpvr         => plt_biom%RootMyco2ndStrutElms_rpvr          ,& !inoput :root layer element secondary axes, [g d-2]    
    GrainFillDowreg_brch              => plt_photo%GrainFillDowreg_brch              ,& !input  :grain fill down-regulation of annual plants, [-]
    PlantPopuLive_pft                 => plt_site%PlantPopuLive_pft                  ,& !input  :plant population, [d-2]    
    Root1stLenPP_rpvr                 => plt_morph%Root1stLenPP_rpvr                 ,& !input :primary root axis length in soil layer per plant, [m d-2]    
    rECLiveCRoot_pft                  => plt_allom%rECLiveCRoot_pft                  ,& !input  element:C ratio of live coarse root, [gE gC-1]    
    ZERO4Groth_pft                    => plt_biom%ZERO4Groth_pft                     ,& !input  :threshold zero for plang growth calculation, [-]
    fTgrowRootP_vr                    => plt_pheno%fTgrowRootP_vr                    ,& !input  :root layer temperature growth functiom, [-]
    RAutoRootO2Limter_rpvr            => plt_rbgc%RAutoRootO2Limter_rpvr             ,& !input  :O2 constraint to root respiration (0-1), [-]
    RootMyco1stSinkC_rpvr             => plt_rbgc%RootMyco1stSinkC_rpvr              ,& !output :primary root C sink, [gC d-2 h-1]
    RootMycoNonstElms_rpvr            => plt_biom%RootMycoNonstElms_rpvr             ,& !inoput :root layer nonstructural element, [g d-2]
    RootCO2EmisPot_pvr                => plt_rbgc%RootCO2EmisPot_pvr                 ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootMaintDef_CO2_pvr              => plt_bgcr%RootMaintDef_CO2_pvr               ,& !inoput :plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]    
    RootCO2Autor_pvr                  => plt_rbgc%RootCO2Autor_pvr                   ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    RootMyco1stElm_raxs               => plt_biom%RootMyco1stElm_raxs                ,& !inoput :root C biomass of primary axes, [g d-2]
    RootCRRadius0_rpvr                => plt_morph%RootCRRadius0_rpvr                ,& !inoput :initial radius of roots that may undergo secondary growth, [m]
    RootRespPotent_pvr                => plt_rbgc%RootRespPotent_pvr                  & !inoput :root respiration unconstrained by O2, [g d-2 h-1]
  )
  call PrintInfo('beg '//subname)

!  mass_inital=0._r8
!  DO L1=max(1,L-1),L
!    DO NE=1,NumPlantChemElms
!      mass_inital(NE) = mass_inital(NE) +RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      DO N1=1,Myco_pft(NZ)
!        mass_inital(NE) = mass_inital(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  write(490,*)yearIJ%I*1000+yearIJ%J/24.,'1'//subname,mass_inital(ielmn)

  !when lignification is not explicitly modeled  
  litrflx=0._r8
  
  if(Root1stActStructElms_rpvr(ielmn,L,NR,NZ).LT.1.e-8_r8)return

!  mass_inital=0._r8
!  DO L1=max(1,L-1),L
!    DO NE=1,NumPlantChemElms
!      mass_inital(NE) = mass_inital(NE) +RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      DO N1=1,Myco_pft(NZ)
!        mass_inital(NE) = mass_inital(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  write(490,*)yearIJ%I*1000+yearIJ%J/24.,subname,mass_inital(ielmn)

  !if the plant has deep roots, aka taproot 
  !grow the maturation zone
  IF(RootSinkC_vr(N,L).GT.ZERO4Groth_pft(NZ))THEN
    FracRoot1stCSinkL=AZMAX1(Root1stSink_pvr(L,NR)/RootSinkC_vr(N,L),1.e-9_r8)
  ELSE
    FracRoot1stCSinkL=1.0_r8
  ENDIF

  call DebugPrint('FracRoot1stCSinkL',FracRoot1stCSinkL)

  !apply nutrient limitation to growth respiraiton
  DMRespEff = 1._r8-StalkBiomGrowthYld_pft(NZ)       ![gC CO2/gC nonst]
  ZPOOLB    = AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))
  PPOOLB    = AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))
  FNP       = AMIN1(ZPOOLB/rECLiveCRoot_pft(ielmn,NZ),PPOOLB/rECLiveCRoot_pft(ielmp,NZ))*DMRespEff/StalkBiomGrowthYld_pft(NZ)*FracRoot1stCSinkL
  !
  !Nutrient depletion prevents structural growth, but living roots still
  !require maintenance and may remobilize biomass when carbon is insufficient.
  if(isclose(FNP,0._r8))FNP=0._r8

  !no explicit lignification, but lump thickening and lignification together
  if(LDoneTip)then
    !Tip maintenance already done
    Rmaint1st_CO2 = 0._r8
  else
    Rmaint1st_CO2 = RmSpecPlant*TFN6_vr(L)*AZMAX1(Root1stActStructElms_rpvr(ielmn,L,NR,NZ))  
  endif

  if(iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    Rmaint1st_CO2=Rmaint1st_CO2*fRootGrowPSISense
  endif
  
  !mobilize storage carbon (in CO2) to support maintenance, thickening and lignification 

  RNonstCO2_OUltd=AZMAX1(AMIN1(FracRoot1stCSinkL*VMXCSpec,1._r8)*RootMycoNonstElms_rpvr(ielmc,N,L,NZ))

  !apply oxygen dependence
  RNonstCO2_Oltd  = AMIN1(RNonstCO2_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ),FracRoot1stCSinkL*RootMycoNonstElms_rpvr(ielmc,N,L,NZ))

!  write(433,*)RootMycoNonstElms_rpvr(ielmc,N,L,NZ),RCO2T1st_Oltd,'xxxx',DMRespEff

  !diagnose growth respiraiton
  RCO2PotGroth_OUltd           = RNonstCO2_OUltd-Rmaint1st_CO2        !excessive CO2-resp O2-unlimited
  RCO2PotGroth_Oltd            = RNonstCO2_Oltd-Rmaint1st_CO2         !excessive CO2-resp O2-limited
  
  RGrowCO2_OUltd               = AZMAX1(RCO2PotGroth_OUltd)                              !whole population respiration for active zone growth+lignification
  RGrowCO2_Oltd                = AZMAX1(RCO2PotGroth_Oltd)                               !total respiration for active zone growth+lignification
  !apply nutrient limitation
  RGrowCO2_Oltd                = AMIN1(RGrowCO2_Oltd,FNP*RAutoRootO2Limter_rpvr(N,L,NZ)) 
  RootMaintDef_CO2_pvr(N,L,NZ) = RootMaintDef_CO2_pvr(N,L,NZ)+AMIN1(RCO2PotGroth_Oltd,0._r8)

  RCO2T1st_OUltd = AMIN1(Rmaint1st_CO2,RNonstCO2_OUltd)
  RCO2T1st_Oltd  = AMIN1(Rmaint1st_CO2,RNonstCO2_Oltd)

  !----------------------------------------
  !positive growth
  if(RCO2PotGroth_Oltd.GT.0._r8)then
    !thickening nonst demand
    call RootThickenWaterFunc(N,L,NZ,Root1stRadius_rpvr(L,NR,NZ),RespThickenWatSens,dR1stExp)

    !carbon required for thickening
    !CHyOz+ O2 -> (Yld_thick)C + (1-Yld_thick)CO2
    !Radius, length, and radial expansion are per root; demand is for the population.
    !The axis count already includes PlantPopuLive_pft, so apply it only once.
    RootMycoNonst4Thick_OUltd(ielmc) = TwoPiCON*Root1stRadius_rpvr(L,NR,NZ) &
      *Root1stLenPP_rpvr(L,NR,NZ)*dR1stExp*Num1stAxesPerStructRootXPOP_pft(NZ) &
      /CRootActVolPerMassC_pft(NZ)
    RootMycoNonst4Thick_OUltd(ielmn) = RootMycoNonst4Thick_OUltd(ielmc)*rECLiveCRoot_pft(ielmn,NZ)
    RootMycoNonst4Thick_OUltd(ielmp) = RootMycoNonst4Thick_OUltd(ielmc)*rECLiveCRoot_pft(ielmp,NZ)  
    !
    !total nonst biomass C required for thickening
    RootMycoNonstC4Thick_OUltd       = RootMycoNonst4Thick_OUltd(ielmc)/StalkBiomGrowthYld_pft(NZ)
    RThickCO2_OUltd                  = RootMycoNonstC4Thick_OUltd*DMRespEff
    !
    !nonst biomass converted to coarse root biomass
    RootMycoNonst4Thick_Oltd(ielmc) = RootMycoNonst4Thick_OUltd(ielmc)*RAutoRootO2Limter_rpvr(N,L,NZ)
    RootMycoNonst4Thick_Oltd(ielmn) = RootMycoNonst4Thick_Oltd(ielmc)*rECLiveCRoot_pft(ielmn,NZ)
    RootMycoNonst4Thick_Oltd(ielmp) = RootMycoNonst4Thick_Oltd(ielmc)*rECLiveCRoot_pft(ielmp,NZ)
    !
    !total O2-limited nonst C required for thickening
    RootMycoNonstC4Thick_Oltd = RootMycoNonstC4Thick_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
    RThickCO2_Oltd            = RThickCO2_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
    
    !O2-limited
    RCO2_meta_Oltd=RThickCO2_Oltd
    if(RCO2_meta_Oltd .GT. RGrowCO2_Oltd)then
      !then scale down
      scal           = RGrowCO2_Oltd/RCO2_meta_Oltd
      RThickCO2_Oltd = RThickCO2_Oltd*scal
      RootMycoNonstC4Thick_Oltd = RootMycoNonstC4Thick_Oltd*scal
      RootMycoNonst4Thick_Oltd=RootMycoNonst4Thick_Oltd*scal
    endif

    !summarize and substract CO2
    RCO2T1st_OUltd = RCO2T1st_OUltd + RThickCO2_OUltd
    RCO2T1st_Oltd  = RCO2T1st_Oltd  + RThickCO2_Oltd

    RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-RCO2T1st_Oltd
    RootMycoNonstC4Lig_Oltd              = 0._r8
    RootMycoNonst4Lig_Oltd               = 0._r8
     
    !summarize NonstC for growth: thickening+lignification
    RootMycoDemand = RootMycoNonst4Thick_Oltd(ielmc) + RootMycoNonstC4Lig_Oltd*Yld_lignif

    !check C limitation
    if(RootMycoDemand.GT.RootMycoNonstElms_rpvr(ielmc,N,L,NZ) .and. RootMycoNonstElms_rpvr(ielmc,N,L,NZ).GT.0._r8)then
      !C-starvation
      scal                     = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)/RootMycoDemand*0.9999_r8
      RootMycoNonstC4Lig_Oltd  = RootMycoNonstC4Lig_Oltd*scal
      RootMycoNonst4Thick_Oltd = RootMycoNonst4Thick_Oltd*scal
      RootMycoNonst4Lig_Oltd   = RootMycoNonst4Lig_Oltd*scal
    endif

    !active lignification
    if(RootMycoNonstC4Lig_Oltd.GT.0._r8)then
      !nonst C lignified into inactive zone 
      RootMycoNonst4Lig_Oltd(ielmc) = RootMycoNonstC4Lig_Oltd*Yld_lignif

      !N and P remobilization from lignification of active zone
      RootMycoNonst4Lig_Oltd(ielmn) = -Rlignif_Oltd*(rECLiveCRoot_pft(ielmn,NZ)-rNCLigRoot_pft(NZ)*CRootActVolPerMassC_pft(NZ)/CRootLigVolPerMassC_pft(NZ))
      RootMycoNonst4Lig_Oltd(ielmp) = -Rlignif_Oltd*(rECLiveCRoot_pft(ielmp,NZ)-rPCLigRoot_pft(NZ)*CRootActVolPerMassC_pft(NZ)/CRootLigVolPerMassC_pft(NZ))

      DO NE=1,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)-RootMycoNonst4Lig_Oltd(NE))
      ENDDO         

      !apply the structrual growth    
      DO NE=1,NumPlantChemElms
        RootMyco1stElm_raxs(NE,NR,NZ)         = RootMyco1stElm_raxs(NE,NR,NZ)+RootMycoNonst4Lig_Oltd(NE)
        RootMyco1stStrutElms_rpvr(NE,L,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)+RootMycoNonst4Lig_Oltd(NE)
        Root1stLigStructElms_rpvr(NE,L,NR,NZ) = Root1stLigStructElms_rpvr(NE,L,NR,NZ)+RootMycoNonst4Lig_Oltd(NE)      
      ENDDO

      !structural biomass transfer from active to lignified zone
      Root1stActStructElms_rpvr(ielmc,L,NR,NZ) = Root1stActStructElms_rpvr(ielmc,L,NR,NZ)-Rlignif_Oltd
      Root1stLigStructElms_rpvr(ielmc,L,NR,NZ) = Root1stLigStructElms_rpvr(ielmc,L,NR,NZ)+Rlignif_Oltd

      Root1stActStructElms_rpvr(ielmn,L,NR,NZ) = Root1stActStructElms_rpvr(ielmn,L,NR,NZ)-Rlignif_Oltd*rECLiveCRoot_pft(ielmn,NZ)
      Root1stLigStructElms_rpvr(ielmn,L,NR,NZ) = Root1stLigStructElms_rpvr(ielmn,L,NR,NZ)+Rlignif_Oltd*rECLiveCRoot_pft(ielmn,NZ)+RootMycoNonst4Lig_Oltd(ielmn)

      Root1stActStructElms_rpvr(ielmp,L,NR,NZ) = Root1stActStructElms_rpvr(ielmp,L,NR,NZ)-Rlignif_Oltd*rECLiveCRoot_pft(ielmp,NZ)
      Root1stLigStructElms_rpvr(ielmp,L,NR,NZ) = Root1stLigStructElms_rpvr(ielmp,L,NR,NZ)+Rlignif_Oltd*rECLiveCRoot_pft(ielmp,NZ)+RootMycoNonst4Lig_Oltd(ielmp)
    endif

    !apply the thickening
    DO NE=1,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)-RootMycoNonst4Thick_Oltd(NE))
    ENDDO

  else
    !potential growth is zero due to maintenance deficit
    RootMycoNonst4Thick_Oltd(:)=0._r8

    RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = &
      RootMycoNonstElms_rpvr(ielmc,N,L,NZ) - RCO2T1st_Oltd
  endif  

!  mass_finale=0._r8;mass1=0._r8;mass2=0._r8;mass3=0._r8
!  DO L1=max(1,L-1),L
!    DO NE=1,NumPlantChemElms
!      mass_finale(NE) = mass_finale(NE) +RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      mass1(NE)=mass1(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      DO N1=1,Myco_pft(NZ)
!        mass_finale(NE) = mass_finale(NE) +RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!        mass2(NE)=mass2(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)
!        mass3(NE)=mass3(NE)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  errC=mass_finale(ielmc)-mass_inital(ielmc)+RCO2T1st_Oltd+RootMycoNonst4Grow_Oltd(ielmc)
!  write(433,*)mass_finale(ielmc),mass_inital(ielmc),'thick111dm',errC,'RCO2T1st_Oltd=',RCO2T1st_Oltd,RCO2PotGroth_Oltd,RootMycoNonst4Grow_Oltd(ielmc),L,NR,yearIJ%I*1000+yearIJ%J/24.
!  if(abs(errC)>1.e-4)stop
!  if(I>=140)write(434,*)'mss1',mass1(ielmc),mass2(ielmc),mass3(ielmc),RootMycoNonst4Thick_Oltd(ielmc)  

  Root1stNetGrowthElms(:)  = RootMycoNonst4Thick_Oltd(:)

  RootMyco1stSinkC_rpvr(L,NR,NZ)=RCO2T1st_Oltd+RootMycoNonst4Thick_Oltd(ielmc)
  
  !update negative growth due to remobilization
  if(-RCO2PotGroth_Oltd.GT.0._r8)then
    !remobilization may turn Root1stNetGrowthElms negative
    call RemobilizePrimeActRoots(N,L,NZ,NR,RCO2PotGroth_Oltd,RCO2PotGroth_OUltd,Root1stcylc,RCO2T1st_OUltd,&
      RCO2T1st_Oltd,Root1stNetGrowthElms,litrflx,RootTipZone=.false.)
    !  
    DO NE=1,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)+Root1stcylc(NE))
    ENDDO 
  endif  
  !
  !if primary root negative growth is negative, offset negative growth by withdrawing from secondary roots, 
  IF(Root1stNetGrowthElms(ielmc).LT.0.0_r8)THEN     
    !negative growth may be offset by 2nd root death   
    call WithDrawMediumRoots(N,NZ,L,NR,Root1stNetGrowthElms,litrflx)      
  ENDIF
  
  !apply the structrual growth
  DO NE=1,NumPlantChemElms
    RootMyco1stElm_raxs(NE,NR,NZ)         = RootMyco1stElm_raxs(NE,NR,NZ)+Root1stNetGrowthElms(NE)
    Root1stActStructElms_rpvr(NE,L,NR,NZ) = Root1stActStructElms_rpvr(NE,L,NR,NZ)+Root1stNetGrowthElms(NE)
    RootMyco1stStrutElms_rpvr(NE,L,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)+Root1stNetGrowthElms(NE)
  ENDDO
  !
  !maintenance deficit is still more than maitenance paid from remobilization  
  !
  RootRespPotent_pvr(N,L,NZ) = RootRespPotent_pvr(N,L,NZ)+RCO2T1st_OUltd
  RootCO2EmisPot_pvr(N,L,NZ) = RootCO2EmisPot_pvr(N,L,NZ)+RCO2T1st_Oltd
  RootCO2Autor_pvr(N,L,NZ)   = RootCO2Autor_pvr(N,L,NZ)-RCO2T1st_Oltd  
  RCO2flxt                   = RCO2flxt-RCO2T1st_Oltd  
  !
  dRThickening=Root1stRadius_rpvr(L,NR,NZ)
  Root1stRadius_rpvr(L,NR,NZ)=sqrt(RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)*CRootActVolPerMassC_pft(NZ) &
    /(PiCON*Root1stLenPP_rpvr(L,NR,NZ)*Num1stAxesPerStructRootXPOP_pft(NZ)))

  !the layer to be pushed into dead zone
  dRThickening=Root1stRadius_rpvr(L,NR,NZ)-dRThickening 

  call CRootsLive2DeadTransformation(yearIJ,NR,L,NZ,Root1stRadius_rpvr(L,NR,NZ),dRThickening)

  litrflxt=litrflxt+litrflx
!  mass_finale=0._r8;mass1=0._r8;mass2=0._r8;mass3=0._r8
!  DO L1=max(1,L-1),L
!    DO NE=1,NumPlantChemElms
!      mass_finale(NE) = mass_finale(NE) +RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      mass1(NE)=mass1(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      DO N1=1,Myco_pft(NZ)
!        mass_finale(NE) = mass_finale(NE) +RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!        mass2(NE)=mass2(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)
!        mass3(NE)=mass3(NE)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  errC=mass_finale(ielmn)-mass_inital(ielmn)+litrflx(ielmn)
!  write(490,*)yearIJ%I*1000+yearIJ%J/24.,subname,errC,mass_finale(ielmn),mass_inital(ielmn),litrflx(ielmn)
!  write(433,*)yearIJ%I*1000+yearIJ%J/24.,RootMycoNonstElms_rpvr(ielmc,N,L,NZ),mass_finale(ielmc),'thick',errC,'co2',RCO2T1st_Oltd,litrflxt(ielmc),&
!    'groth',Root1stNetGrowthElms(ielmc),RootMycoNonst4Grow_Oltd(ielmc),L,NR
!  if(abs(errC)>1.e-4)call endrun(trim(mod_filename)//' at line',__LINE__)           
!  write(433,*)yearIJ%I*1000+yearIJ%J/24.,L,NR,mass_finale(ielmn)-mass_inital(ielmn)+litrflxt(ielmn),mass_finale(ielmn),mass_inital(ielmn)
  call PrintInfo('end '//subname)
  end associate
  end subroutine StructuralRootThickening
!----------------------------------------------------------------------------------------------------
  subroutine WithDrawMediumRoots(N,NZ,L,NR,RootNetGrowthElms,litrflxt)      
  implicit none
  integer, intent(in) :: N        !root type, plant/myco
  integer, intent(in) :: NZ,L,NR
  real(r8), intent(inout) :: RootNetGrowthElms(NumPlantChemElms)
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)          !litterfall due to withdraw of secondary/fine roots
  real(r8) :: litrflx(NumPlantChemElms)
  integer :: LX,LL,NE,M,N1
  real(r8) :: GRTWTM     !negative primary root C growth
  real(r8) :: FSNCM,FSNCP
  real(r8) :: dRootMyco2ndst2Litr(NumPlantChemElms)
  real(r8) :: dRootMycoNonst2Litr(NumPlantChemElms)
  real(r8) :: dsenecE,dwoodyE,dfineE1,dfineE2
  character(len=*), parameter :: subname='WithDrawMediumRoots'
  real(r8) :: mass_inital(NumPlantChemElms)

  associate(                                                          &
    Myco_pft                  => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]  
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft             ,& !input  :threshold zero for plang growth calculation, [-]
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr    ,& !input  :root layer structural C, [gC d-2]
    icwood                    => pltpar%icwood                       ,& !input  :group id of coarse woody litter
    PlantElmAllocMat4Litr     => plt_soilchem%PlantElmAllocMat4Litr  ,& !input  :litter kinetic fraction, [-]
    FracRootElmAllocm         => plt_allom%FracRootElmAllocm         ,& !input  :C woody fraction in root,[-]
    iroot                     => pltpar%iroot                        ,& !input  :group id of plant root litter
    k_woody_comp              => pltpar%k_woody_comp                 ,& !input  :woody litter complex id
    k_fine_comp               => pltpar%k_fine_comp                  ,& !input  :fine litter complex id
    inonstruct                => pltpar%inonstruct                   ,& !input  :group id of plant nonstructural litter
    RootMediumStructElms_rpvr => plt_biom%RootMediumStructElms_rpvr  ,& !inoput :root layer element for medium size root axes, [g d-2]
    LitrfallElms_pvr          => plt_bgcr%LitrfallElms_pvr           ,& !inoput :plant LitrFall element, [g d-2 h-1]
    RootMediumXNum_rpvr       => plt_morph%RootMediumXNum_rpvr       ,& !inout :living population medium-root axes by layer and group, [d-2]
    RootMediumLength_rpvr     => plt_morph%RootMediumLength_rpvr      & !inoput :total length for medium size root for axis NR, [m d-2]    
  )
  call PrintInfo('beg '//subname)

  !withdraw fine roots
  call Withdraw2ndRoots(N,NZ,L,NR,RootNetGrowthElms,litrflxt)

  IF(any(RootNetGrowthElms.LT.0.0_r8))THEN    

    !remove secondary roots and offset some of the negative growth
    DO NE=1,NumPlantChemElms
      !root dying, offset some negative growth, but not CO2, assuming CO2 can only be offset from reserve and primary roots
      IF(RootNetGrowthElms(NE).LT.0.0_r8)THEN    
        !negative growth is less than total fine roots in layer    
        IF(-RootNetGrowthElms(NE).LE.RootMediumStructElms_rpvr(NE,L,NR,NZ))THEN
          !root partially dead, length thrinked
          if(NE.EQ.ielmc)then
            RootMediumLength_rpvr(L,NR,NZ)=RootMediumLength_rpvr(L,NR,NZ)*(1._r8+RootNetGrowthElms(NE)/RootMediumStructElms_rpvr(NE,L,NR,NZ))
          endif   
          RootMediumStructElms_rpvr(NE,L,NR,NZ) = RootMediumStructElms_rpvr(NE,L,NR,NZ)+RootNetGrowthElms(NE)
          RootNetGrowthElms(NE)                    = 0._r8
        ELSE
          !root dead completely
          if(NE.EQ.ielmc)RootMediumLength_rpvr(L,NR,NZ)=0._r8
          RootNetGrowthElms(NE)                    = RootNetGrowthElms(NE)+RootMediumStructElms_rpvr(NE,L,NR,NZ)
          RootMediumStructElms_rpvr(NE,L,NR,NZ) = 0._r8
        ENDIF
      ENDIF
    ENDDO  
  endif
  !Partial withdrawal shortens surviving axes; complete loss removes them.
  if(RootMediumStructElms_rpvr(ielmc,L,NR,NZ).LE.0._r8 .or. &
     RootMediumLength_rpvr(L,NR,NZ).LE.0._r8)RootMediumXNum_rpvr(L,NR,NZ)=0._r8
  call PrintInfo('end '//subname)
  end associate
  end subroutine WithDrawMediumRoots
!----------------------------------------------------------------------------------------------------

  subroutine StructuralRootElongation(yearIJ,N,NR,L,Lnext,NZ,VMXCSpec,Root1stSinkTip,RootSinkC_vr,&
    fRootGrowPSISense,TFN6_vr,fdLext1st_mech,RespElongWatSens,DMRespEff,CNRTW,CPRTW,litrflxt,RCO2flxt)
  !
  !Do elongation of the primary roots  
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  integer,  intent(in) :: N,NR,NZ
  integer,  intent(in) :: L,Lnext       !current and next layer for root growth
  REAL(R8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)  
  real(r8), intent(in) :: fRootGrowPSISense
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  real(r8), intent(in) :: fdLext1st_mech       !soil compaction modification of elongation
  real(r8), intent(in) :: RespElongWatSens
  real(r8), intent(in) :: Root1stSinkTip(pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: CNRTW,CPRTW  
  real(r8), intent(in) :: DMRespEff                       !respiraiton efficiency for new dry root mass creation from nonstructural C
  real(r8), intent(in) :: VMXCSpec
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)
  real(r8), intent(inout) :: RCO2flxt       !accumulator of CO2 flux
  real(r8) :: Root1stExtenzPP
  real(r8) :: litrflx(NumPlantChemElms)
  real(r8) :: FracRoot1stCSinkL
  real(r8) :: GRTWTM  
  real(r8) :: TFRCO2,FRCO2  
  real(r8) :: Rmaint1st_CO2
  real(r8) :: RNonstCO2_OUltd
  real(r8) :: RGrowCO2_Oltd,RGrowCO2_OUltd               !growth respiration unltd,ltd by O2 and unlimited by N,P
  real(r8) :: RNonstCO2_Oltd
  real(r8) :: RCO2PotGroth_OUltd,RCO2PotGroth_Oltd       !O2-unlimited/limited potential growth respiration
  real(r8) :: DMRTR,ZPOOLB,PPOOLB,FNP
  real(r8) :: RootMycoNonst4GrowC_OUltd
  real(r8) :: RootMycoNonst4Grow_Oltd(NumPlantChemElms)  !nonstructural N,P unltd, not ltd by O2 used in growth
  real(r8) :: RootMycoNonst4Grow_OUltd(NumPlantChemElms) !nonstructural N,P unltd, ltd by O2 used in growth
  real(r8) :: RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd
  real(r8) :: RCO2T1st_OUltd,RCO2T1st_Oltd
  real(r8) :: dsenecE
  real(r8) :: RootMycoNonst4GrowC_Oltd
  real(r8) :: Root1stcylc(NumPlantChemElms)              !nonstrucal biomass from remobilization
  real(r8) :: Root1stNetGrowthElms(NumPlantChemElms)
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: RCO2MDef1stStructPayoff_Oltd,RCO2MDef1stStructPayoff_OUltd
  real(r8) :: mass_finale(NumPlantChemElms),RLEN,RootNsink1st
  integer  :: NE,M,LL,LL1,L1,N1
  real(r8) :: RRadius,ratio,rootmy1
  real(r8) :: Root1stPotExtzPP             !per axial pontential/actual root extension due to nonstrucal C mobilization.  
  real(r8), parameter :: RTDPX1=0.005_r8  !the apical growing zone length
  character(len=*), parameter :: subname='Grow1stRootElongZoneLayL'

  associate(                                                                          &
    Myco_pft                          => plt_morph%Myco_pft                          ,& !input  :mycorrhizal type (no or yes),[-]    
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft              ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    fTgrowRootP_vr                    => plt_pheno%fTgrowRootP_vr                    ,& !input  :root layer temperature growth functiom, [-]
    RAutoRootO2Limter_rpvr            => plt_rbgc%RAutoRootO2Limter_rpvr             ,& !input  :O2 constraint to root respiration (0-1), [-]
    MaxNumRootLays                    => plt_site%MaxNumRootLays                     ,& !input  :maximum root layer number,[-]
    CumSoilThickness_vr               => plt_site%CumSoilThickness_vr                ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    ZERO4Groth_pft                    => plt_biom%ZERO4Groth_pft                     ,& !input  :threshold zero for plang growth calculation, [-]
    CNRTS_pft                         => plt_allom%CNRTS_pft                         ,& !input  :root N:C ratio x root growth yield, [-]
    CPRTS_pft                         => plt_allom%CPRTS_pft                         ,& !input  :root P:C ratio x root growth yield, [-]
    RootBiomGrosYld_pft               => plt_allom%RootBiomGrosYld_pft               ,& !input  :root growth yield, [g g-1]
    DLYR3                             => plt_site%DLYR3                              ,& !input  :vertical thickness of soil layer, [m]        
    Root1stRadius_pvr                 => plt_morph%Root1stRadius_pvr                 ,& !input  :root layer radius primary axes, [m]  
    PlantElmAllocMat4Litr             => plt_soilchem%PlantElmAllocMat4Litr          ,& !input  :litter kinetic fraction, [-]    
    StalkBiomGrowthYld_pft            => plt_allom%StalkBiomGrowthYld_pft            ,& !input  :stalk growth yield, [gC gC-1]    
    Root1stMaxRadius1_pft             => plt_morph%Root1stMaxRadius1_pft             ,& !input  :root radius primary axes, [m]    
    RootMatureAge_pft                 => plt_morph%RootMatureAge_pft                 ,& !input : Root maturation age, [h]
    PSIRoot_pvr                       => plt_ew%PSIRoot_pvr                          ,& !input :root total water potential, [Mpa]    
    RootAge_rpvr                      => plt_morph%RootAge_rpvr                      ,& !inoput :root age,[h]
    LitrfallElms_pvr                  => plt_bgcr%LitrfallElms_pvr                   ,& !input  :plant LitrFall element, [g d-2 h-1]    
    PlantPopuLive_pft                 => plt_site%PlantPopuLive_pft                  ,& !input  :plant population, [d-2]
    Num1stAxesPerStructRootXPOP_pft   => plt_morph%Num1stAxesPerStructRootXPOP_pft   ,& !output :population primary root axes number on one structrual axis, [d-2]        
    iPlant2ndGrothPattern_pft         => plt_pheno%iPlant2ndGrothPattern_pft         ,& !input  :plant expression of secondary growth, [-]                
    Root1stSpecLen_pft                => plt_morph%Root1stSpecLen_pft                ,& !input  :specific root length primary axes, [m root gC-1]
    k_fine_comp                       => pltpar%k_fine_comp                          ,& !input  :fine litter complex id    
    iroot                             => pltpar%iroot                                ,& !input  :group id of plant root litter    
    Root1stLenPP_rpvr                 => plt_morph%Root1stLenPP_rpvr                 ,& !input  :primary root axis length in soil layer, [m d-2]
    SeedDepth_pft                     => plt_morph%SeedDepth_pft                     ,& !input  :seeding depth, [m]
    flag2ndGrowth_pvr                 => plt_morph%flag2ndGrowth_pvr                 ,& !inoput  :flag for secondary growth of primary roots, [-]        
    NGTopRootLayer_pft                => plt_morph%NGTopRootLayer_pft                ,& !input  :soil layer at planting depth, [-]
    Root1stDepz_raxes                 => plt_morph%Root1stDepz_raxes                 ,& !input  :root layer depth, [m]
    Root1stMaxRadius_pft              => plt_morph%Root1stMaxRadius_pft              ,& !input  :maximum radius of primary roots, [m]    
    MainBranchNum_pft                 => plt_morph%MainBranchNum_pft                 ,& !input  :number of main branch,[-]
    FracRootElmAllocm                 => plt_allom%FracRootElmAllocm                 ,& !input  :C woody fraction in root,[-]    
    Root1stAxesTipDepz2Surf_pft       => plt_morph%Root1stAxesTipDepz2Surf_pft       ,& !output :plant primary depth relative to column surface, [m]     
    iPlantRootProfile_pft             => plt_pheno%iPlantRootProfile_pft             ,& !input  :plant growth type (vascular, non-vascular),[-]
    GrainFillDowreg_brch              => plt_photo%GrainFillDowreg_brch              ,& !input  :grain fill down-regulation of annual plants, [-]
    RootMaintDef_CO2_pvr              => plt_bgcr%RootMaintDef_CO2_pvr               ,& !inoput :plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]    
    RootCO2Autor_pvr                  => plt_rbgc%RootCO2Autor_pvr                   ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    RootMyco1stSinkC_rpvr             => plt_rbgc%RootMyco1stSinkC_rpvr              ,& !output :1st root C sink, [gC d-2 h-1]
    RootCO2EmisPot_pvr                => plt_rbgc%RootCO2EmisPot_pvr                 ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootRespPotent_pvr                => plt_rbgc%RootRespPotent_pvr                 ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]
    RootMycoNonstElms_rpvr            => plt_biom%RootMycoNonstElms_rpvr             ,& !inoput :root layer nonstructural element, [g d-2]
    RootMyco1stElm_raxs               => plt_biom%RootMyco1stElm_raxs                ,& !inoput :root layer primary axes structrual element, [g d-2]
    Root1stActStructElms_rpvr         => plt_biom%Root1stActStructElms_rpvr          ,& !inoput :Root layer primary axes Active zone structrual element, [g d-2]
    RootMyco2ndStrutElms_rpvr         => plt_biom%RootMyco2ndStrutElms_rpvr          ,& !inoput :root layer element secondary axes, [g d-2]
    RootMyco1stStrutElms_rpvr         => plt_biom%RootMyco1stStrutElms_rpvr          ,& !inoput :root layer element primary axes, [g d-2]
    NRoot1stTipLay_raxes              => plt_morph%NRoot1stTipLay_raxes               & !inoput :soil layer number for deepest root axes, [-]
  )
  call PrintInfo('beg '//subname)
  !
  !     FRACTION OF PRIMARY ROOT SINK IN SOIL LAYER
  !     ATTRIBUTED TO CURRENT AXIS
  !
  !     Root1stSink_pvr=primary root sink strength
  !     RootSinkC_vr=total root sink strength
  !     FracRoot1stCSinkL=fraction of primary root sink strength in axis
  !
  litrflx     = 0._r8
  !  mass_inital = 0._r8
  !  DO L1  = max(1,L-1), Lnext
  !    DO NE=1,NumPlantChemElms
  !      mass_inital(NE) = mass_inital(NE) +RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
  !    ENDDO  
  !    DO N1=1,Myco_pft(NZ)
  !      DO NE=1,NumPlantChemElms
  !        mass_inital(NE) = mass_inital(NE) +RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
  !      ENDDO
  !    ENDDO
  !  ENDDO
  !
  IF(RootSinkC_vr(N,L).GT.ZERO4Groth_pft(NZ))THEN
    FracRoot1stCSinkL=Root1stSinkTip(NR)/RootSinkC_vr(N,L)
  ELSE
    FracRoot1stCSinkL=1.0_r8
  ENDIF
  !
  !     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
  !     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
  !
  !     N,P CONSTRAINT ON PRIMARY ROOT RESPIRATION FROM
  !     NON-STRUCTURAL C:N:P
  !
  !     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
  !     
  !     CNKI,CPKI=nonstructural N,P inhibition constant on growth
  !
  !     PRIMARY ROOT MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
  !     ROOT STRUCTURAL N
  !
  !     Rmaint1st_CO2=root maintenance respiration
  !     RmSpecPlant=specific maintenance respiration rate (g C g-1 N h-1)
  !     WTRT1N=primary root N mass
  !     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
  !     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
  !     fRootGrowPSISense=growth function of root water potential
  !
  if(lcoarseroot .and. is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))then
    Rmaint1st_CO2 = AZMAX1(RmSpecPlant*Root1stActStructElms_rpvr(ielmn,L,NR,NZ))*TFN6_vr(L)
  else
    Rmaint1st_CO2 = AZMAX1(RmSpecPlant*RootMyco1stElm_raxs(ielmn,NR,NZ))*TFN6_vr(L)
  endif

  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ)) .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    Rmaint1st_CO2=Rmaint1st_CO2*fRootGrowPSISense
  ENDIF
  !
  !     O2-UNLIMITED PRIMARY ROOT RESPIRATION FROM ROOT NON-STRUCTURAL C
  !     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
  !
  !     RNonstCO2_OUltd=respiration from non-structural C unlimited by O2
  !     VMXC=rate constant for nonstructural C oxidation in respiration C     
  !     FracRoot1stCSinkL=fraction of primary root sink strength in axis
  !     CPOOL=non-structural C mass
  !     fTgrowRootP_vr=temperature function for root growth
  !     GrainFillDowreg_brch=termination feedback inhibition on C3 CO2
  !     fRootGrowPSISense=growth function of root water potential
  !
  RNonstCO2_OUltd=AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ)*AMIN1(FracRoot1stCSinkL*VMXCSpec,1._r8))

  !do not let the root grow beyond the maximum allowed layer
  !This is a temporary fix for the lack of deep soil toxicty that prevents root growth
  IF(Root1stAxesTipDepz2Surf_pft(NR,NZ).GE.CumSoilThickness_vr(MaxNumRootLays))THEN
    RNonstCO2_OUltd=AMIN1(Rmaint1st_CO2,RNonstCO2_OUltd)
  ENDIF
  !
  !     O2-LIMITED PRIMARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
  !
  !     RNonstCO2_Oltd=respiration from non-structural C limited by O2
  !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
  !     RespElongWatSens=respiration function of root water potential
  !
  RNonstCO2_Oltd               = RNonstCO2_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
  RCO2PotGroth_OUltd           = RNonstCO2_OUltd-Rmaint1st_CO2
  RCO2PotGroth_Oltd            = RNonstCO2_Oltd-Rmaint1st_CO2
  RGrowCO2_OUltd               = AZMAX1(RCO2PotGroth_OUltd)*RespElongWatSens
  RGrowCO2_Oltd                = AZMAX1(RCO2PotGroth_Oltd)*RespElongWatSens
  RootMaintDef_CO2_pvr(N,L,NZ) = RootMaintDef_CO2_pvr(N,L,NZ)+AMIN1(RCO2PotGroth_Oltd,0._r8)
  !
  !     PRIMARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
  !     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
  !
  !     FracRoot1stCSinkL=fraction of secondary root sink strength in axis
  !     ZPOOLR,PPOOLR=non-structural N,P mass in root 
  !     CNRTS_pft,CPRTS_pft=N,P root growth yield = = [gN,P root/gC root]*[gC root /gC nonst]=[gN,P root/gC nonstC]
  !     FNP=growth respiration limited by non-structural N,P
  !     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
  !
  ZPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))*FracRoot1stCSinkL
  PPOOLB = AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))*FracRoot1stCSinkL
  FNP    = AMIN1(ZPOOLB/CNRTS_pft(NZ),PPOOLB/CPRTS_pft(NZ))*DMRespEff

  IF(RGrowCO2_OUltd.GT.0.0_r8)THEN
    RGrowCO2_OUltd=AMIN1(RGrowCO2_OUltd,FNP)
  ELSE
    RGrowCO2_OUltd=0._r8
  ENDIF

  IF(RGrowCO2_Oltd.GT.0.0_r8)THEN
    RGrowCO2_Oltd=AMIN1(RGrowCO2_Oltd,FNP*RAutoRootO2Limter_rpvr(N,L,NZ))
  ELSE
    RGrowCO2_Oltd=0._r8
  ENDIF
  !
  !     TOTAL NON-STRUCTURAL C,N,P USED IN PRIMARY ROOT GROWTH
  !     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
  !     ENTERED IN 'READQ'
  !
  !     RootMycoNonst4GrowC_OUltd,RootMycoNonst4GrowC_Oltd=total non-structural C used in growth and respn unltd,ltd by O2
  !     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
  !     RootMycoNonst4Grow_OUltd(ielmc),RootMycoNonst4Grow_Oltd(ielmc)=root C growth unltd,ltd by O2
  !     RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
  !
  RootMycoNonst4GrowC_OUltd       = RGrowCO2_OUltd/DMRespEff
  RootMycoNonst4GrowC_Oltd        = RGrowCO2_Oltd/DMRespEff

  RootMycoNonst4Grow_OUltd(ielmc) = RootMycoNonst4GrowC_OUltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_OUltd(ielmn) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CNRTW)
  RootMycoNonst4Grow_OUltd(ielmp) = AZMAX1(RootMycoNonst4Grow_OUltd(ielmc)*CPRTW)
  !
  !nonstrucal biomass elements for metabolism
  RootMycoNonst4Grow_Oltd(ielmc) = RootMycoNonst4GrowC_Oltd*RootBiomGrosYld_pft(NZ)
  RootMycoNonst4Grow_Oltd(ielmn) = AZMAX1(AMIN1(FracRoot1stCSinkL*RootMycoNonstElms_rpvr(ielmn,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*CNRTW))
  RootMycoNonst4Grow_Oltd(ielmp) = AZMAX1(AMIN1(FracRoot1stCSinkL*RootMycoNonstElms_rpvr(ielmp,N,L,NZ),RootMycoNonst4Grow_Oltd(ielmc)*CPRTW))
  RootMycoNonst4Grow_Oltd        = RootMycoNonst4Grow_Oltd*fdLext1st_mech
  RCO2Nonst4Nassim_OUltd         = AZMAX1(1.70_r8*RootMycoNonst4Grow_OUltd(ielmn))
  RCO2Nonst4Nassim_Oltd          = AZMAX1(1.70_r8*RootMycoNonst4Grow_Oltd(ielmn))
  !
  !     TOTAL PRIMARY ROOT RESPIRATION
  !
  !     RCO2T1st_OUltd,RCO2T1st_Oltd=total C respiration unltd,ltd by O2
  !     Rmaint1st_CO2=root maintenance respiration
  !     RNonstCO2_OUltd,RNonstCO2_Oltd=respiration from non-structural C unltd,ltd by O2
  !     RGrowCO2_OUltd,RGrowCO2_Oltd=growth respiration limited by N,P unltd,ltd by O2
  !     RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
  !     RootCO2Autor_pvr=total root respiration
  !     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
  !
  RCO2T1st_OUltd = AMIN1(Rmaint1st_CO2,RNonstCO2_OUltd)+RGrowCO2_OUltd+RCO2Nonst4Nassim_OUltd
  RCO2T1st_Oltd  = AMIN1(Rmaint1st_CO2,RNonstCO2_Oltd)+RGrowCO2_Oltd+RCO2Nonst4Nassim_Oltd
  !
  !subtract CO2 from reserve 
  RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-RCO2T1st_Oltd

  DO NE=1,NumPlantChemElms
    RootMycoNonst4Grow_Oltd(NE)       = AMIN1(RootMycoNonstElms_rpvr(NE,N,L,NZ),RootMycoNonst4Grow_Oltd(NE))
    RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)-RootMycoNonst4Grow_Oltd(NE))
  ENDDO

  !
  !primary root extension based on growth flux, accounting for all root axes
  Root1stPotExtzPP = RootMycoNonst4Grow_Oltd(ielmc)*FracRootElmAllocm(ielmc,k_fine_comp)/Num1stAxesPerStructRootXPOP_pft(NZ)*Root1stSpecLen_pft(N,NZ)
  Root1stExtenzPP  = Root1stPotExtzPP
  
  !dealing with negative growth

  Root1stNetGrowthElms(:)=RootMycoNonst4Grow_Oltd(:) ! > 0

  RootMyco1stSinkC_rpvr(L,NR,NZ)=RCO2T1st_Oltd+RootMycoNonst4Grow_Oltd(ielmc)

!  write(981,*)yearIJ%I*1000+yearIJ%J/24.,RCO2PotGroth_Oltd,Root1stPotExtzPP,NR,L,fdLext1st_mech,&
!    Root1stAxesTipDepz2Surf_pft(NR,NZ).GE.CumSoilThickness_vr(MaxNumRootLays),Root1stDepz_raxes(NR,NZ),&
!    CumSoilThickness_vr(MaxNumRootLays),Root1stAxesTipDepz2Surf_pft(NR,NZ)

!  mass_finale=0._r8
!  DO L1=max(1,L-1),Lnext
!    DO NE=1,NumPlantChemElms
!      mass_finale(NE)=mass_finale(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!    ENDDO    
!    DO N1=1,Myco_pft(NZ)
!      DO NE=1,NumPlantChemElms
!        mass_finale(NE)=mass_finale(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  write(490,*)yearIJ%I*1000+yearIJ%J/24.,'bfremob',mass_finale(ielmn)-mass_inital(ielmn),Root1stNetGrowthElms(ielmn),RCO2PotGroth_Oltd

  if(-RCO2PotGroth_Oltd.GT.0._r8)then !Root1stNetGrowthElm and Root1stPotExtzPP are zero   
    !may turn Root1stNetGrowthElms negative due to remobilzation
    !
    call RemobilizePrimeActRoots(N,L,NZ,NR,RCO2PotGroth_Oltd,RCO2PotGroth_OUltd,Root1stcylc,RCO2T1st_OUltd,&
      RCO2T1st_Oltd,Root1stNetGrowthElms,litrflx,RootTipZone=.True.)
    
    DO NE=1,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)+Root1stcylc(NE))
    ENDDO         
  endif

!  mass_finale=0._r8
!  DO L1=max(1,L-1),Lnext
!    DO NE=1,NumPlantChemElms
!      mass_finale(NE)=mass_finale(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!    ENDDO    
!    DO N1=1,Myco_pft(NZ)
!      DO NE=1,NumPlantChemElms
!        mass_finale(NE)=mass_finale(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  write(490,*)yearIJ%I*1000+yearIJ%J/24.,'afremob',mass_finale(ielmn)-mass_inital(ielmn)+litrflx(ielmn)+Root1stNetGrowthElms(ielmn),Root1stNetGrowthElms(ielmn),litrflx(ielmn)

  !
  !     ALLOCATE PRIMARY ROOT TOTAL RESPIRATION TO ALL SOIL LAYERS
  !     THROUGH WHICH PRIMARY ROOTS GROW
  !
  !     Root1stLenPP_rpvr=primary root length, assuming primary roots are of the same radius
  !     SeedDepth_pft=seeding depth
  !     FRCO2=fraction of primary root respiration attributed to layer
  !     RootCO2Autor_pvr=total root respiration
  !     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
  !
  IF(Root1stDepz_raxes(NR,NZ).GT.CumSoilThickness_vr(NGTopRootLayer_pft(NZ)))THEN
    TFRCO2=0._r8;RLEN=0._r8
    LL1=NGTopRootLayer_pft(NZ)

    ! Exclude the mature prefix; retain the tip as a fallback.
    if(is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))then
      DO WHILE(LL1.LT.NRoot1stTipLay_raxes(NR,NZ))
        IF(.not.flag2ndGrowth_pvr(LL1,NR,NZ))exit
        LL1=LL1+1
      ENDDO
    endif

    ! Normalize using exactly the layers that receive this respiration.
    RLEN=SUM(MAX(0._r8,Root1stLenPP_rpvr(LL1:NRoot1stTipLay_raxes(NR,NZ),NR,NZ)))
    D5100: DO LL=LL1,NRoot1stTipLay_raxes(NR,NZ)
      IF(LL.EQ.NRoot1stTipLay_raxes(NR,NZ))THEN
        FRCO2=MAX(0._r8,1._r8-TFRCO2)
      ELSEIF(RLEN.GT.0._r8)THEN
        FRCO2=MIN(MAX(0._r8,1._r8-TFRCO2),MAX(0._r8,Root1stLenPP_rpvr(LL,NR,NZ))/RLEN)
      ELSE
        FRCO2=0._r8
      ENDIF
      TFRCO2                      = TFRCO2+FRCO2
      RootRespPotent_pvr(N,LL,NZ) = RootRespPotent_pvr(N,LL,NZ)+RCO2T1st_OUltd*FRCO2
      RootCO2EmisPot_pvr(N,LL,NZ) = RootCO2EmisPot_pvr(N,LL,NZ)+RCO2T1st_Oltd*FRCO2
      RootCO2Autor_pvr(N,LL,NZ)   = RootCO2Autor_pvr(N,LL,NZ)-RCO2T1st_Oltd*FRCO2
      RCO2flxt                    = RCO2flxt-RCO2T1st_Oltd*FRCO2
    ENDDO D5100
  ELSE
    RootRespPotent_pvr(N,L,NZ) = RootRespPotent_pvr(N,L,NZ)+RCO2T1st_OUltd
    RootCO2EmisPot_pvr(N,L,NZ) = RootCO2EmisPot_pvr(N,L,NZ)+RCO2T1st_Oltd
    RootCO2Autor_pvr(N,L,NZ)   = RootCO2Autor_pvr(N,L,NZ)-RCO2T1st_Oltd
    RCO2flxt                   = RCO2flxt-RCO2T1st_Oltd
  ENDIF
  !
  !     ALLOCATE ANY NEGATIVE PRIMARY ROOT C,N,P GROWTH TO SECONDARY
  !     ROOTS ON THE SAME AXIS IN THE SAME LAYER UNTIL SECONDARY ROOTS
  !     HAVE DISAPPEARED
  !
  !     Frac2Senes1=fraction of primary root C to be remobilized
  !     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
  !     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
  !     Root2ndLen_rpvr=secondary root length
  !
  !if primary root growth is negative, then offset negative growth by killing secondary roots
  IF(Root1stNetGrowthElms(ielmc).LT.0.0_r8)THEN   
    call WithDrawMediumRoots(N,NZ,L,NR,Root1stNetGrowthElms,litrflx)                     
  ENDIF

!  mass_finale=0._r8
!  DO L1=max(1,L-1),Lnext
!    DO NE=1,NumPlantChemElms
!      mass_finale(NE)=mass_finale(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!    ENDDO    
!    DO N1=1,Myco_pft(NZ)
!      DO NE=1,NumPlantChemElms
!        mass_finale(NE)=mass_finale(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  write(490,*)yearIJ%I*1000+yearIJ%J/24.,'af2ndrot',mass_finale(ielmn)-mass_inital(ielmn)+litrflx(ielmn)+Root1stNetGrowthElms(ielmn),Root1stNetGrowthElms(ielmn),'litr',litrflx(ielmn)

  IF(Root1stNetGrowthElms(ielmc).LT.0.0_r8)THEN
    !negative growth
    IF(RootMyco1stElm_raxs(ielmc,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
      !primary roots withdraw, note that primary root depth was initialized at seedDepth
      !make sure it fits the corase root plants     
      Root1stExtenzPP=(Root1stDepz_raxes(NR,NZ)-SeedDepth_pft(NZ))/RootMyco1stElm_raxs(ielmc,NR,NZ)*Root1stNetGrowthElms(ielmc) 
    ELSE
      !enforce mass balance under a special edge case.
      DO M=1,jsken      
        DO NE=1,NumPlantChemElms
          LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) &
            +PlantElmAllocMat4Litr(NE,iroot,M,NZ)*Root1stNetGrowthElms(NE)
        ENDDO  
      ENDDO
      Root1stNetGrowthElms=0._r8
    ENDIF
  ENDIF
  !
  !the extension should not exceed next soil layer thickness  
  !this means: 
  !1. the extension is vertical, should the extension be per axial?
  !2. the extension may exceed current layer thickness, but not the next layer
  IF(L.LT.MaxNumRootLays)THEN
    Root1stExtenzPP=AMIN1(DLYR3(Lnext),Root1stExtenzPP)
  ENDIF
  
  !the extension can be negative or positive
  call PrimeRootsExtension(yearIJ,L,Lnext,N,NR,NZ,Root1stNetGrowthElms,Root1stExtenzPP)

!  mass_finale=0._r8
!  DO L1=max(1,L-1),Lnext
!    DO NE=1,NumPlantChemElms
!      mass_finale(NE)=mass_finale(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!    ENDDO    
!    DO N1=1,Myco_pft(NZ)
!      DO NE=1,NumPlantChemElms
!        mass_finale(NE)=mass_finale(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
  litrflxt=litrflxt+litrflx
!  if(I>=143)write(433,*)'elong',mass_finale(ielmc)-mass_inital(ielmc)+RCO2T1st_Oltd+litrflxt(ielmc),RCO2T1st_Oltd,litrflxt(ielmc),'RCO2flxt',RCO2flxt
!  write(490,*)'elong',mass_finale(ielmn)-mass_inital(ielmn)+litrflxt(ielmn),mass_finale(ielmn),mass_inital(ielmn),litrflxt(ielmn)
  call PrintInfo('end '//subname)
  
  end associate
  end subroutine StructuralRootElongation
!----------------------------------------------------------------------------------------------------
  subroutine GrowFibrous1stRootAxes(yearIJ,N,NR,L,Lnext,NZ,VMXCSpec,Root1stSink_pvr,Root2ndSink_pvr,&
    Root1stSinkTip,RootSinkC_vr,fRootGrowPSISense,TFN6_vr,fdLext1st,RespElongWatSens,&
    DMRespEff,CNRTW,CPRTW,Root1stTipUpdateFlag,FoundRootAxesTip,litrflxt,RCO2flxt)
  !
  !Description:
  !grow primary root axes
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer,  intent(in) :: N,NR,NZ
  integer,  intent(in) :: L,Lnext        !current and next layer for root growth
  REAL(R8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)  
  real(r8), intent(in) :: fRootGrowPSISense
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  real(r8), intent(in) :: fdLext1st,RespElongWatSens
  real(r8), intent(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: Root1stSinkTip(pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)  
  real(r8), intent(in) :: CNRTW,CPRTW  
  real(r8), intent(in) :: DMRespEff              !respiraiton efficiency for fresh root growth
  real(r8), intent(in) :: VMXCSpec
  logical,  intent(inout) :: Root1stTipUpdateFlag(pltpar%MaxNumRootAxes)
  logical,  intent(inout) :: FoundRootAxesTip(pltpar%jroots,pltpar%MaxNumRootAxes)  
  real(r8), intent(out) :: litrflxt(NumPlantChemElms)
  real(r8), intent(out) :: RCO2flxt
  real(r8) :: Root1stExtenzPP
  logical :: RootRetreated
  real(r8) :: Frac2Senes1
  real(r8) :: GRTWTM  
  real(r8) :: TFRCO2,FRCO2  
  real(r8) :: RCO2Nonst4Nassim_OUltd,RCO2Nonst4Nassim_Oltd
  real(r8) :: RCO2T1st_OUltd,RCO2T1st_Oltd
  real(r8) :: FSNCM,FSNCP,CNG,CPG  
  real(r8) :: dsenecE,SoilResit4RootPentration
  real(r8) :: dRootMycoElms(NumPlantChemElms)  
  real(r8) :: Root1stNetGrowthElms(NumPlantChemElms)
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: DMTapRTD,dmass(NumPlantChemElms),dmass0(NumPlantChemElms)
  integer  :: NE,M,L1,N1
  logical  :: FoundRootTipLayer
  logical  :: LDoneTip
  character(len=*), parameter :: subname='GrowFibrous1stRootAxes'

  associate(                                                                          &
    Myco_pft                          => plt_morph%Myco_pft                          ,& !input  :mycorrhizal type (no or yes),[-]    
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft              ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    MaxNumRootLays                    => plt_site%MaxNumRootLays                     ,& !input  :maximum root layer number,[-]
    CumSoilThickness_vr               => plt_site%CumSoilThickness_vr                ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    ZERO4Groth_pft                    => plt_biom%ZERO4Groth_pft                     ,& !input  :threshold zero for plang growth calculation, [-]
    Num1stAxesPerStructRootX_pft      => plt_morph%Num1stAxesPerStructRootX_pft      ,& !input  :number of primary root axes per NR, [d-2]        
    Root1stLenPP_rpvr                 => plt_morph%Root1stLenPP_rpvr                 ,& !input  :primary root axis length in soil layer, [m d-2]
    NGTopRootLayer_pft                => plt_morph%NGTopRootLayer_pft                ,& !input  :soil layer at planting depth, [-]
    iPlantRootProfile_pft             => plt_pheno%iPlantRootProfile_pft             ,& !input  :plant growth type (vascular, non-vascular),[-]
    RootMatureAge_pft                 => plt_morph%RootMatureAge_pft                 ,& !input  : Root maturation age, [h]
    RootAge_rpvr                      => plt_morph%RootAge_rpvr                      ,& !inoput :root age,[h]
    PlantPopuLive_pft                 => plt_site%PlantPopuLive_pft                  ,& !input  :plant population, [d-2]        
    flag2ndGrowth_pvr                 => plt_morph%flag2ndGrowth_pvr                 ,& !input  :flag for secondary growth of primary roots, [-]
    Root1stAxesTipDepz2Surf_pft       => plt_morph%Root1stAxesTipDepz2Surf_pft       ,& !output :plant primary depth relative to column surface, [m]     
    Num1stAxesPerStructRootXPOP_pft   => plt_morph%Num1stAxesPerStructRootXPOP_pft   ,& !output :population primary root axes number on one structrual axis, [d-2]        
    RootCO2Autor_pvr                  => plt_rbgc%RootCO2Autor_pvr                   ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    RootCO2EmisPot_pvr                => plt_rbgc%RootCO2EmisPot_pvr                 ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootRespPotent_pvr                => plt_rbgc%RootRespPotent_pvr                 ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]
    RootMycoNonstElms_rpvr            => plt_biom%RootMycoNonstElms_rpvr             ,& !inoput :root layer nonstructural element, [g d-2]
    RootMyco1stElm_raxs               => plt_biom%RootMyco1stElm_raxs                ,& !inoput :root chemical element for primary axes, [g d-2]
    RootMyco2ndStrutElms_rpvr         => plt_biom%RootMyco2ndStrutElms_rpvr          ,& !inoput :soil layer cheimcal element of secondary roots, [g d-2]
    RootMyco1stStrutElms_rpvr         => plt_biom%RootMyco1stStrutElms_rpvr          ,& !inoput :soil layer chemical element of primary roots, [g d-2]
    Root1stActStructElms_rpvr         => plt_biom%Root1stActStructElms_rpvr          ,& !inoput :root layer active zone element in primary axes, [g d-2]
    Root1stXNumL_pvr                  => plt_morph%Root1stXNumL_pvr                  ,& !inoput :soil layer number of primary root axes, [d-2]
    NRoot1stTipLay_raxes              => plt_morph%NRoot1stTipLay_raxes               & !inoput :soil layer number for deepest root axes, [-]
  )
  !   call SumRootBiome(NZ,mass_inital)

  call PrintInfo('beg '//subname)

  litrflxt=0._r8;RCO2flxt=0._r8
  !identify root tip and do root elongation
  !root pass the top surface of layer L, and the primary root referenced by (N,NR) has not been updated  

  FoundRootTipLayer=Root1stAxesTipDepz2Surf_pft(NR,NZ).LE.CumSoilThickness_vr(L) .OR. L.EQ.MaxNumRootLays
  !  write(490,*)'-----------------------------------------------------------'

  !make sure the last root layer has not been updated
  if(.not.Root1stTipUpdateFlag(NR))then 
    !diagnose primary root number in layer L
    Root1stXNumL_pvr(L,NZ) = Root1stXNumL_pvr(L,NZ)+Num1stAxesPerStructRootXPOP_pft(NZ) 
    RootAge_rpvr(L,NR,NZ)  = RootAge_rpvr(L,NR,NZ)+1._r8

!    dmass0=0._r8
!    DO L1=max(1,L-1),Lnext
!    DO NE=1,NumPlantChemElms
!      dmass0(NE)=dmass0(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ) 
!      DO N1=1,1
!        dmass0(NE)=dmass0(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)
!      ENDDO
!    ENDDO
!    ENDDO
!    write(490,*)L,NR,yearIJ%I*1000+yearIJ%J/24.,'beftip',dmass0(ielmn),litrflxt(ielmn)

    RootRetreated=.false.
    LDoneTip=.false.
    if(FoundRootTipLayer .and. Root1stSinkTip(NR).GT.0._r8)THEN
      !reached last root layer
      Root1stTipUpdateFlag(NR) = .true.
      LDoneTip                 = .true.
      call Grow1stRootElongZoneLayL(yearIJ,N,NR,L,Lnext,NZ,VMXCSpec,Root1stSinkTip,RootSinkC_vr,&
        fRootGrowPSISense,TFN6_vr,fdLext1st,RespElongWatSens,DMRespEff,CNRTW,CPRTW,litrflxt,RCO2flxt,Root1stExtenzPP)
      RootRetreated=Root1stExtenzPP.LT.0._r8

!      dmass=0._r8
!      DO L1=max(1,L-1),Lnext
!      DO NE=1,NumPlantChemElms
!        dmass(NE)=dmass(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!        DO N1=1,Myco_pft(NZ)
!          dmass(NE)=dmass(NE)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)
!        ENDDO
!      ENDDO
!      ENDDO
!      write(490,*)L,NR,yearIJ%I*1000+yearIJ%J/24.,'tip',dmass(ielmn)-dmass0(ielmn)+litrflxt(ielmn), dmass(ielmn),dmass0(ielmn),litrflxt(ielmn)
    endif    

!    dmass0=0._r8
!    DO L1=max(1,L-1),L
!    DO NE=1,NumPlantChemElms
!      dmass0(NE)=dmass0(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ) 
!      DO N1=1,Myco_pft(NZ)
!        dmass0(NE)=dmass0(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)
!      ENDDO
!    ENDDO
!    ENDDO
!    write(490,*)L,NR,yearIJ%I*1000+yearIJ%J/24.,'2ndbf',dmass0(ielmn),litrflxt(ielmn),L
       
!    if(flag2ndGrowth_pvr(L,NR,NZ) .and. Root1stLenPP_rpvr(L,NR,NZ).GT.0._r8)then
!      
!      call SecondaryGrowthZone(yearIJ,N,NR,L,NZ,Nutstress4GrossResp,RootSinkC_vr,Root1stSink_pvr,fRootGrowPSISense,&
!        TFN6_vr,litrflxt,RCO2flxt,LDoneTip)
!    endif

!    dmass=0._r8
!    DO L1=max(1,L-1),L
!    DO NE=1,NumPlantChemElms
!      dmass(NE)=dmass(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      DO N1=1,Myco_pft(NZ)
!        dmass(NE)=dmass(NE)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)
!      ENDDO
!    ENDDO
!    ENDDO
!    dmass(ielmc)=dmass(ielmc)-dmass0(ielmc)-RCO2flxt+litrflxt(ielmc) 

!    write(490,*)L,NR,yearIJ%I*1000+yearIJ%J/24.,'2nd',dmass(ielmn)-dmass0(ielmn)+litrflxt(ielmn), dmass(ielmn),dmass0(ielmn),litrflxt(ielmn)
    !
    !consider root withraw at root tip
    IF(RootRetreated)THEN
      !Use the old tip layer L; the extension helper has updated the tip index.
      call PrimeRootsWithdraw(L,NR,NZ,N,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr)
    ENDIF

    !    dmass=0._r8
    !    DO L1=max(1,L-1),Lnext
    !    DO NE=1,NumPlantChemElms
    !      dmass(NE)=dmass(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
    !      DO N1=1,Myco_pft(NZ)
    !        dmass(NE)=dmass(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)
    !      ENDDO
    !    ENDDO
    !    ENDDO
    !    dmass(ielmc)=dmass(ielmc)-dmass0(ielmc)-RCO2flxt+litrflxt(ielmc) 

    !  REMOVE ANY NEGATIVE ROOT MASS FROM NONSTRUCTURAL C
    IF(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ).LT.0.0_r8)THEN
      DO NE=1,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,N,L,NZ)       = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))
        RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ) = 0._r8
      ENDDO
    ENDIF

    !  what if nonstructural C is insufficient to meet the negative root mass
    IF(RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ).LT.0.0_r8)THEN
      DO NE=1,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,N,L,NZ)     = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootMyco1stStrutElms_rpvr(NE,L,NR,NZ))
        RootMyco1stElm_raxs(NE,NR,NZ)         = RootMyco1stElm_raxs(NE,NR,NZ)+RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)
        RootMyco1stStrutElms_rpvr(NE,L,NR,NZ) = 0._r8 
        Root1stActStructElms_rpvr(NE,L,NR,NZ) = 0._R8    
      ENDDO
    ENDIF

    NRoot1stTipLay_raxes(NR,NZ) = MIN(NRoot1stTipLay_raxes(NR,NZ),MaxNumRootLays)
    FoundRootAxesTip(N,NR)      = RootRetreated .OR. (L.EQ.NRoot1stTipLay_raxes(NR,NZ))

    !    dmass=0._r8
    !    DO L1=max(1,L-1),Lnext
    !    DO NE=1,NumPlantChemElms      
    !      dmass(NE)=dmass(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
    !      DO N1=1,Myco_pft(NZ)
    !        dmass(NE)=dmass(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)
    !      ENDDO
    !    ENDDO
    !    ENDDO
    !    dmass(ielmc)=dmass(ielmc)-dmass0(ielmc)-RCO2flxt+litrflxt(ielmc) 

    !    if(I>=143)write(433,*)I*1000+J/24.,dmass(ielmc),NR,L,'1strot',RCO2flxt,litrflxt(ielmc) 
    !    if(I>=143 .and. abs(dmass(ielmc))>1.e-6_r8)stop
  ENDIF

  call PrintInfo('end '//subname)
  end associate
  end subroutine GrowFibrous1stRootAxes
!----------------------------------------------------------------------------------------------------

  subroutine SecondaryGrowthZone(yearIJ,N,NR,L,NZ,VMXCSpec,RootSinkC_vr,Root1stSink_pvr, &
    fRootGrowPSISense,TFN6_vr,litrflxt,RCO2flxt,LDoneTip)
  !
  !Description:
  !Do secondary growth of roots
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: N,NR,L,NZ
  real(r8), intent(in) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: TFN6_vr(JZ1)  
  real(r8), intent(in) :: fRootGrowPSISense
  real(r8), intent(in) :: VMXCSpec
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)
  real(r8), intent(inout) :: RCO2flxt
  logical, intent(in) :: LDoneTip
  real(r8) :: DMRespEff
  real(r8) :: FracRoot1stCSinkL
  real(r8) :: Rmaint1st_CO2
  real(r8) :: RNonstCO2_OUltd,RGrowCO2_Oltd,RGrowCO2_OUltd
  real(r8) :: RNonstCO2_Oltd,RCO2PotGroth_OUltd,RCO2PotGroth_Oltd
  real(r8) :: DMRTR,ZPOOLB,PPOOLB,FNP
  real(r8) :: RootMycoNonst4GrowC_Oltd
  real(r8) :: RootMycoNonst4GrowC_OUltd,RootMycoNonst4Lig_OUltd(NumPlantChemElms)
  real(r8) :: RootMycoNonst4Lig_Oltd(NumPlantChemElms)
  real(r8) :: RootMycoNonstC4Thick_Oltd,RootMycoNonstC4Thick_OUltd
  real(r8) :: RootMycoNonstC4Lig_Oltd,RootMycoNonstC4Lig_OUltd     !Nonstructural C used lignification
  real(r8) :: RootMycoNonst4Thick_OUltd(NumPlantChemElms)
  real(r8) :: RootMycoNonst4Thick_Oltd(NumPlantChemElms)   
  real(r8) :: RootMycoNonst4Grow_Oltd(NumPlantChemElms)    !nonstrucal biomass converted into structrual biomass with O2 limitation
  real(r8) :: RootMycoNonst4Grow_OUltd(NumPlantChemElms)   !nonstrucal biomass converted into structrual biomass without O2 limitation
  real(r8) :: Root1stNetGrowthElms(NumPlantChemElms)
  real(r8) :: Root1stcylc(NumPlantChemElms)              !nonstrucal biomass from remobilization
  real(r8) :: litrflx(NumPlantChemElms)
  real(r8) :: RCO2T1st_OUltd,RCO2T1st_Oltd,RootNsink1st
  real(r8) :: RootThickenWatSens,RespThickenWatSens,dR1stExp,dRPotexp,dRExp,ratio
  real(r8) :: mass_inital(NumPlantChemElms),mass_finale(NumPlantChemElms)
  real(r8) :: RCO2MDef1stStructPayoff_Oltd,RCO2MDef1stStructPayoff_OUltd
  real(r8) :: Rlignif_OUltd, Rlignif_Oltd      !lignification rate of active zone C without/with O limitation
  real(r8) :: RligCO2_OUltd, RligCO2_Oltd      !CO2 respiration for lignification without/with O limitation
  real(r8) :: CO2Qt_lig                        !CO2 quotient of lignification
  real(r8) :: RThickCO2_OUltd, RThickCO2_Oltd  !CO2 for coarse root thickening
  real(r8) :: scal
  real(r8) :: dRCO2MaintDef_Oltd, dRCO2MaintDef_OUltd !CO2 deficit for maintenance after active zone mobilization  
  real(r8) :: alpha_fill                         !C-filling quotient during lignification, [gC reserve /gC lignified-active biomass]
  real(r8) :: mass1(NumPlantChemElms),mass2(NumPlantChemElms),mass3(NumPlantChemElms)
  real(r8) :: RCO2_meta_Oltd,RCO2_meta_OUltd
  real(r8) :: RootMycoDemand,dRThickening
  real(r8) :: f_stress, K_ligmax,errC
  integer  :: NE,L1,N1

  character(len=*), parameter :: subname='SecondaryGrowthZone'

  associate(                                                                          &
    Myco_pft                          => plt_morph%Myco_pft                          ,& !input  :mycorrhizal type (no or yes),[-]    
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft              ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    MainBranchNum_pft                 => plt_morph%MainBranchNum_pft                 ,& !input  :number of main branch,[-]
    RootMyco1stStrutElms_rpvr         => plt_biom%RootMyco1stStrutElms_rpvr          ,& !input  :root layer element primary axes, [g d-2]
    Root1stActStructElms_rpvr         => plt_biom%Root1stActStructElms_rpvr          ,& !inoput :root layer active zone element in primary axes, [g d-2]
    Root1stLigStructElms_rpvr         => plt_biom%Root1stLigStructElms_rpvr          ,& !inoput :root layer lignified zone element in primary axes, [g d-2]
    StalkBiomGrowthYld_pft            => plt_allom%StalkBiomGrowthYld_pft            ,& !input  :stalk growth yield, [gC gC-1]    
    Root1stRadius_pvr                 => plt_morph%Root1stRadius_pvr                 ,& !input  :root layer radius primary axes, [m]  
    rNCRoot_pft                       => plt_allom%rNCRoot_pft                       ,& !input  :root N:C ratio, [gN gC-1]
    Root1stRadius_rpvr                => plt_morph%Root1stRadius_rpvr                ,& !input: root layer radius for each primary axes,  [m]    
    KLigMax_pft                       => plt_biom%KLigMax_pft                        ,& !input  :Maximum lignification rate [h-1]
    KLigMM_pft                        => plt_biom%KLigMM_pft                         ,& !input  :Half saturation parameter for coarse root lignification [h-1]
    CRootActVolPerMassC_pft           => plt_morph%CRootActVolPerMassC_pft           ,& !input  :coarse root active zone volume:mass ratio, [m3 gC-1]    
    CRootLigVolPerMassC_pft           => plt_morph%CRootLigVolPerMassC_pft           ,& !input  :coarse root inactive zone volume:mass ratio, [m3 gC-1]     
    rNCLigRoot_pft                    => plt_allom%rNCLigRoot_pft                    ,& !input  :NC ratio of lignified coarse root, [gN gC-1]
    rPCLigRoot_pft                    => plt_allom%rPCLigRoot_pft                    ,& !input  :PC ratio of lignified coarse root, [gN gC-1]
    Num1stAxesPerStructRootXPOP_pft   => plt_morph%Num1stAxesPerStructRootXPOP_pft   ,& !output :population primary root axes number on one structrual axis, [d-2]        
    RootMyco2ndStrutElms_rpvr         => plt_biom%RootMyco2ndStrutElms_rpvr          ,& !inoput :root layer element secondary axes, [g d-2]    
    PlantPopuLive_pft                 => plt_site%PlantPopuLive_pft                  ,& !input  :plant population, [d-2]    
    Root1stLenPP_rpvr                 => plt_morph%Root1stLenPP_rpvr                 ,& !input :primary root axis length in soil layer per plant, [m d-2]    
    rECLiveCRoot_pft                  => plt_allom%rECLiveCRoot_pft                  ,& !input  element:C ratio of live coarse root, [gE gC-1]    
    ZERO4Groth_pft                    => plt_biom%ZERO4Groth_pft                     ,& !input  :threshold zero for plang growth calculation, [-]
    fTgrowRootP_vr                    => plt_pheno%fTgrowRootP_vr                    ,& !input  :root layer temperature growth functiom, [-]
    RAutoRootO2Limter_rpvr            => plt_rbgc%RAutoRootO2Limter_rpvr             ,& !input  :O2 constraint to root respiration (0-1), [-]
    RootMyco1stSinkC_rpvr             => plt_rbgc%RootMyco1stSinkC_rpvr              ,& !output :primary root C sink, [gC d-2 h-1]
    RootMycoNonstElms_rpvr            => plt_biom%RootMycoNonstElms_rpvr             ,& !inoput :root layer nonstructural element, [g d-2]
    RootCO2EmisPot_pvr                => plt_rbgc%RootCO2EmisPot_pvr                 ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootMaintDef_CO2_pvr              => plt_bgcr%RootMaintDef_CO2_pvr               ,& !inoput :plant root maintenance respiraiton deficit as CO2, [g d-2 h-1]    
    RootCO2Autor_pvr                  => plt_rbgc%RootCO2Autor_pvr                   ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    RootMyco1stElm_raxs               => plt_biom%RootMyco1stElm_raxs                ,& !inoput :root C biomass of primary axes, [g d-2]
    RootCRRadius0_rpvr                => plt_morph%RootCRRadius0_rpvr                ,& !inoput :initial radius of roots that may undergo secondary growth, [m]
    RootRespPotent_pvr                => plt_rbgc%RootRespPotent_pvr                  & !inoput :root respiration unconstrained by O2, [g d-2 h-1]
  )
  call PrintInfo('beg '//subname)

!  mass_inital=0._r8
!  DO L1=max(1,L-1),L
!    DO NE=1,NumPlantChemElms
!      mass_inital(NE) = mass_inital(NE) +RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      DO N1=1,Myco_pft(NZ)
!        mass_inital(NE) = mass_inital(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  write(490,*)yearIJ%I*1000+yearIJ%J/24.,'1'//subname,mass_inital(ielmn)

  !when lignification is not explicitly modeled  
  litrflx=0._r8
  
  if(Root1stActStructElms_rpvr(ielmn,L,NR,NZ).LT.1.e-8_r8)return

!  mass_inital=0._r8
!  DO L1=max(1,L-1),L
!    DO NE=1,NumPlantChemElms
!      mass_inital(NE) = mass_inital(NE) +RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      DO N1=1,Myco_pft(NZ)
!        mass_inital(NE) = mass_inital(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  write(490,*)yearIJ%I*1000+yearIJ%J/24.,subname,mass_inital(ielmn)

  !if the plant has deep roots, aka taproot 
  !grow the maturation zone
  IF(RootSinkC_vr(N,L).GT.ZERO4Groth_pft(NZ))THEN
    FracRoot1stCSinkL=AZMAX1(Root1stSink_pvr(L,NR)/RootSinkC_vr(N,L),1.e-9_r8)
  ELSE
    FracRoot1stCSinkL=1.0_r8
  ENDIF

  call DebugPrint('FracRoot1stCSinkL',FracRoot1stCSinkL)

  !apply nutrient limitation to growth respiraiton
  DMRespEff = 1._r8-StalkBiomGrowthYld_pft(NZ)       ![gC CO2/gC nonst]
  ZPOOLB    = AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))
  PPOOLB    = AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))
  FNP       = AMIN1(ZPOOLB/rECLiveCRoot_pft(ielmn,NZ),PPOOLB/rECLiveCRoot_pft(ielmp,NZ))*DMRespEff/StalkBiomGrowthYld_pft(NZ)*FracRoot1stCSinkL
  !
  if(isclose(FNP,0._r8))return  

  !no explicit lignification, but lump thickening and lignification together
  if(LDoneTip)then
    !Tip maintenance already done
    Rmaint1st_CO2 = 0._r8
  else
    Rmaint1st_CO2 = RmSpecPlant*TFN6_vr(L)*AZMAX1(Root1stActStructElms_rpvr(ielmn,L,NR,NZ))  
  endif

  if(iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    Rmaint1st_CO2=Rmaint1st_CO2*fRootGrowPSISense
  endif
  
  !mobilize storage carbon (in CO2) to support maintenance, thickening and lignification 

  RNonstCO2_OUltd=AZMAX1(AMIN1(FracRoot1stCSinkL*VMXCSpec,1._r8)*RootMycoNonstElms_rpvr(ielmc,N,L,NZ))

  !apply oxygen dependence
  RNonstCO2_Oltd               = AMIN1(RNonstCO2_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ),FracRoot1stCSinkL*RootMycoNonstElms_rpvr(ielmc,N,L,NZ))

!  write(433,*)RootMycoNonstElms_rpvr(ielmc,N,L,NZ),RCO2T1st_Oltd,'xxxx',DMRespEff

  !diagnose growth respiraiton
  RCO2PotGroth_OUltd           = RNonstCO2_OUltd-Rmaint1st_CO2        !excessive CO2-resp O2-unlimited
  RCO2PotGroth_Oltd            = RNonstCO2_Oltd-Rmaint1st_CO2         !excessive CO2-resp O2-limited
  
  RGrowCO2_OUltd               = AZMAX1(RCO2PotGroth_OUltd)                              !whole population respiration for active zone growth+lignification
  RGrowCO2_Oltd                = AZMAX1(RCO2PotGroth_Oltd)                               !total respiration for active zone growth+lignification
  !apply nutrient limitation
  RGrowCO2_Oltd                = AMIN1(RGrowCO2_Oltd,FNP*RAutoRootO2Limter_rpvr(N,L,NZ)) 
  RootMaintDef_CO2_pvr(N,L,NZ) = RootMaintDef_CO2_pvr(N,L,NZ)+AMIN1(RCO2PotGroth_Oltd,0._r8)

  !----------------------------------------
  !positive growth
  if(RCO2PotGroth_Oltd.GT.0._r8)then
    !thickening nonst demand
    call RootThickenWaterFunc(N,L,NZ,Root1stRadius_rpvr(L,NR,NZ),RespThickenWatSens,dR1stExp)

    !carbon required for thickening
    !CHyOz+ O2 -> (Yld_thick)C + (1-Yld_thick)CO2
    !O2-unlimited biomass production
    RootMycoNonst4Thick_OUltd(ielmc) = (TwoPiCON*Root1stRadius_rpvr(L,NR,NZ)*Root1stLenPP_rpvr(L,NR,NZ)*Num1stAxesPerStructRootXPOP_pft(NZ)) &
      *dR1stExp/CRootActVolPerMassC_pft(NZ)
    RootMycoNonst4Thick_OUltd(ielmn) = RootMycoNonst4Thick_OUltd(ielmc)*rECLiveCRoot_pft(ielmn,NZ)
    RootMycoNonst4Thick_OUltd(ielmp) = RootMycoNonst4Thick_OUltd(ielmc)*rECLiveCRoot_pft(ielmp,NZ)  
    !
    !total nonst biomass C required for thickening
    RootMycoNonstC4Thick_OUltd       = RootMycoNonst4Thick_OUltd(ielmc)/StalkBiomGrowthYld_pft(NZ)
    RThickCO2_OUltd                  = RootMycoNonstC4Thick_OUltd*DMRespEff
    !
    !nonst biomass converted to coarse root biomass
    RootMycoNonst4Thick_Oltd(ielmc) = RootMycoNonst4Thick_OUltd(ielmc)*RAutoRootO2Limter_rpvr(N,L,NZ)
    RootMycoNonst4Thick_Oltd(ielmn) = RootMycoNonst4Thick_Oltd(ielmc)*rECLiveCRoot_pft(ielmn,NZ)
    RootMycoNonst4Thick_Oltd(ielmp) = RootMycoNonst4Thick_Oltd(ielmc)*rECLiveCRoot_pft(ielmp,NZ)
    !
    !total O2-limited nonst C required for thickening
    RootMycoNonstC4Thick_Oltd = RootMycoNonstC4Thick_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
    RThickCO2_Oltd            = RThickCO2_OUltd*RAutoRootO2Limter_rpvr(N,L,NZ)
    
    if(llignification)then
      !computing lignification rate
      !----------------------------------------
      !lignification rate
      f_stress = 1._r8 !stress4lignif(NZ,RAutoRootO2Limter_rpvr(N,L,NZ),Nutstress4GrossResp,fRootGrowPSISense)
      K_ligmax = KLigMax_pft(NZ)*(1._r8+f_stress)
      !respiration quotient of lignification, gC-CO2 per gC lignified
      CO2Qt_lig = (1._r8-Yld_lignif)/Yld_lignif
      !
      !carry-over C mass from active zone for biomass, a function of available energy and active zone C.
      Rlignif_OUltd = K_ligmax*safe_adb(RGrowCO2_OUltd*Root1stActStructElms_rpvr(ielmc,L,NR,NZ),RGrowCO2_OUltd+KLigMM_pft(NZ)*Root1stActStructElms_rpvr(ielmc,L,NR,NZ))
      !relative increase in biomass, (rho_lig-rho_act)/rho_act, extra C biomass needed for converting 1 vol C-act to C-lig
      alpha_fill    = (CRootActVolPerMassC_pft(NZ)-CRootLigVolPerMassC_pft(NZ))/CRootLigVolPerMassC_pft(NZ)
      !
      !CO2 flux required for lignification = condensed C * CO2_quotient_lig
      RligCO2_OUltd = Rlignif_OUltd*alpha_fill*CO2Qt_lig
      !
      !amount of active zone C to be lignified
      Rlignif_Oltd = K_ligmax*safe_adb(RGrowCO2_Oltd*Root1stActStructElms_rpvr(ielmc,L,NR,NZ),RGrowCO2_Oltd+KLigMM_pft(NZ)*Root1stActStructElms_rpvr(ielmc,L,NR,NZ))
      !
      !lignification is a compression and filling process that does not change the bulk volume
      RligCO2_Oltd = Rlignif_Oltd*alpha_fill*CO2Qt_lig
      !
      !downscale when C-limited 
      !O2-unlimited
      RCO2_meta_OUltd= RligCO2_OUltd+RThickCO2_OUltd
      if(RCO2_meta_OUltd .GT. RGrowCO2_OUltd)then
        scal            = RGrowCO2_OUltd/RCO2_meta_OUltd
        RligCO2_OUltd   = RligCO2_OUltd*scal
        RThickCO2_OUltd = RThickCO2_OUltd*scal
        Rlignif_OUltd   = Rlignif_OUltd*scal
      endif
      !
      !O2-limited
      RCO2_meta_Oltd=RligCO2_Oltd + RThickCO2_Oltd
      if(RCO2_meta_Oltd .GT. RGrowCO2_Oltd)then
        !then scale down
        scal           = RGrowCO2_Oltd/RCO2_meta_Oltd
        RligCO2_Oltd   = RligCO2_Oltd*scal
        RThickCO2_Oltd = RThickCO2_Oltd*scal

        !active-C converted into lignifed-C
        Rlignif_Oltd   = Rlignif_Oltd*scal
      endif

      !summarize and substract CO2
      RCO2T1st_OUltd                       = AMIN1(Rmaint1st_CO2,RNonstCO2_OUltd)+RThickCO2_OUltd+RligCO2_OUltd
      RCO2T1st_Oltd                        = AMIN1(Rmaint1st_CO2,RNonstCO2_Oltd)+RThickCO2_Oltd+RligCO2_Oltd
      RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-RCO2T1st_Oltd

      !total nonst C required for lignification
      !CHyOz + O2 -> (Y_lig)C+(1-Y_lig)CO2
      RootMycoNonstC4Lig_OUltd   = RligCO2_OUltd/(1._r8-Yld_lignif)
      RootMycoNonstC4Lig_Oltd    = RligCO2_Oltd/(1._r8-Yld_lignif)
      
      !Nonst converted into lignified C
      RootMycoNonst4Lig_OUltd(ielmc) = RootMycoNonstC4Lig_OUltd*Yld_lignif
      RootMycoNonst4Lig_OUltd(ielmn) = -Rlignif_OUltd*(rECLiveCRoot_pft(ielmn,NZ)-rNCLigRoot_pft(NZ)*CRootActVolPerMassC_pft(NZ)/CRootLigVolPerMassC_pft(NZ))
      RootMycoNonst4Lig_OUltd(ielmp) = -Rlignif_OUltd*(rECLiveCRoot_pft(ielmp,NZ)-rPCLigRoot_pft(NZ)*CRootActVolPerMassC_pft(NZ)/CRootLigVolPerMassC_pft(NZ))

      RootMycoNonst4Lig_Oltd(ielmc) = RootMycoNonstC4Lig_Oltd*Yld_lignif
      RootMycoNonst4Lig_Oltd(ielmn) = -Rlignif_Oltd*(rECLiveCRoot_pft(ielmn,NZ)-rNCLigRoot_pft(NZ)*CRootActVolPerMassC_pft(NZ)/CRootLigVolPerMassC_pft(NZ))
      RootMycoNonst4Lig_Oltd(ielmp) = -Rlignif_Oltd*(rECLiveCRoot_pft(ielmp,NZ)-rPCLigRoot_pft(NZ)*CRootActVolPerMassC_pft(NZ)/CRootLigVolPerMassC_pft(NZ))

    else
      !O2-limited
      RCO2_meta_Oltd=RThickCO2_Oltd
      if(RCO2_meta_Oltd .GT. RGrowCO2_Oltd)then
        !then scale down
        scal           = RGrowCO2_Oltd/RCO2_meta_Oltd
        RligCO2_Oltd   = RligCO2_Oltd*scal
        RThickCO2_Oltd = RThickCO2_Oltd*scal
      endif

      !summarize and substract CO2
      RCO2T1st_OUltd                       = AMIN1(Rmaint1st_CO2,RNonstCO2_OUltd)+RThickCO2_OUltd
      RCO2T1st_Oltd                        = AMIN1(Rmaint1st_CO2,RNonstCO2_Oltd)+RThickCO2_Oltd
!      write(433,*)RootMycoNonstElms_rpvr(ielmc,N,L,NZ),RCO2T1st_Oltd,'xxxx',AMIN1(Rmaint1st_CO2,RNonstCO2_Oltd),RThickCO2_Oltd,RCO2PotGroth_Oltd
      RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-RCO2T1st_Oltd
      RootMycoNonstC4Lig_Oltd              = 0._r8
      RootMycoNonst4Lig_Oltd               = 0._r8
     
    endif

    !summarize NonstC for growth: thickening+lignification
    RootMycoDemand = RootMycoNonst4Thick_Oltd(ielmc) + RootMycoNonstC4Lig_Oltd*Yld_lignif

    !check C limitation
    if(RootMycoDemand.GT.RootMycoNonstElms_rpvr(ielmc,N,L,NZ) .and. RootMycoNonstElms_rpvr(ielmc,N,L,NZ).GT.0._r8)then
      !C-starvation
      scal                     = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)/RootMycoDemand*0.9999_r8
      RootMycoNonstC4Lig_Oltd  = RootMycoNonstC4Lig_Oltd*scal
      RootMycoNonst4Thick_Oltd = RootMycoNonst4Thick_Oltd*scal
      RootMycoNonst4Lig_Oltd   = RootMycoNonst4Lig_Oltd*scal
    endif

    !active lignification
    if(RootMycoNonstC4Lig_Oltd.GT.0._r8)then
      !nonst C lignified into inactive zone 
      RootMycoNonst4Lig_Oltd(ielmc) = RootMycoNonstC4Lig_Oltd*Yld_lignif

      !N and P remobilization from lignification of active zone
      RootMycoNonst4Lig_Oltd(ielmn) = -Rlignif_Oltd*(rECLiveCRoot_pft(ielmn,NZ)-rNCLigRoot_pft(NZ)*CRootActVolPerMassC_pft(NZ)/CRootLigVolPerMassC_pft(NZ))
      RootMycoNonst4Lig_Oltd(ielmp) = -Rlignif_Oltd*(rECLiveCRoot_pft(ielmp,NZ)-rPCLigRoot_pft(NZ)*CRootActVolPerMassC_pft(NZ)/CRootLigVolPerMassC_pft(NZ))

      DO NE=1,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)-RootMycoNonst4Lig_Oltd(NE))
      ENDDO         

      !apply the structrual growth    
      DO NE=1,NumPlantChemElms
        RootMyco1stElm_raxs(NE,NR,NZ)         = RootMyco1stElm_raxs(NE,NR,NZ)+RootMycoNonst4Lig_Oltd(NE)
        RootMyco1stStrutElms_rpvr(NE,L,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)+RootMycoNonst4Lig_Oltd(NE)
        Root1stLigStructElms_rpvr(NE,L,NR,NZ) = Root1stLigStructElms_rpvr(NE,L,NR,NZ)+RootMycoNonst4Lig_Oltd(NE)      
      ENDDO

      !structural biomass transfer from active to lignified zone
      Root1stActStructElms_rpvr(ielmc,L,NR,NZ) = Root1stActStructElms_rpvr(ielmc,L,NR,NZ)-Rlignif_Oltd
      Root1stLigStructElms_rpvr(ielmc,L,NR,NZ) = Root1stLigStructElms_rpvr(ielmc,L,NR,NZ)+Rlignif_Oltd

      Root1stActStructElms_rpvr(ielmn,L,NR,NZ) = Root1stActStructElms_rpvr(ielmn,L,NR,NZ)-Rlignif_Oltd*rECLiveCRoot_pft(ielmn,NZ)
      Root1stLigStructElms_rpvr(ielmn,L,NR,NZ) = Root1stLigStructElms_rpvr(ielmn,L,NR,NZ)+Rlignif_Oltd*rECLiveCRoot_pft(ielmn,NZ)+RootMycoNonst4Lig_Oltd(ielmn)

      Root1stActStructElms_rpvr(ielmp,L,NR,NZ) = Root1stActStructElms_rpvr(ielmp,L,NR,NZ)-Rlignif_Oltd*rECLiveCRoot_pft(ielmp,NZ)
      Root1stLigStructElms_rpvr(ielmp,L,NR,NZ) = Root1stLigStructElms_rpvr(ielmp,L,NR,NZ)+Rlignif_Oltd*rECLiveCRoot_pft(ielmp,NZ)+RootMycoNonst4Lig_Oltd(ielmp)
    endif
!    write(433,*)RootMycoNonstElms_rpvr(ielmc,N,L,NZ),RootMycoNonst4Thick_Oltd(ielmc)
    !apply the thickening
    DO NE=1,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)-RootMycoNonst4Thick_Oltd(NE))
    ENDDO
!    write(433,*)RootMycoNonstElms_rpvr(ielmc,N,L,NZ),RootMycoNonst4Thick_Oltd(ielmc)
  else
    !potential growth is zero due to maintenance deficit
    RootMycoNonst4Thick_Oltd(:)=0._r8
  endif  

!  mass_finale=0._r8;mass1=0._r8;mass2=0._r8;mass3=0._r8
!  DO L1=max(1,L-1),L
!    DO NE=1,NumPlantChemElms
!      mass_finale(NE) = mass_finale(NE) +RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      mass1(NE)=mass1(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      DO N1=1,Myco_pft(NZ)
!        mass_finale(NE) = mass_finale(NE) +RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!        mass2(NE)=mass2(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)
!        mass3(NE)=mass3(NE)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  errC=mass_finale(ielmc)-mass_inital(ielmc)+RCO2T1st_Oltd+RootMycoNonst4Grow_Oltd(ielmc)
!  write(433,*)mass_finale(ielmc),mass_inital(ielmc),'thick111dm',errC,'RCO2T1st_Oltd=',RCO2T1st_Oltd,RCO2PotGroth_Oltd,RootMycoNonst4Grow_Oltd(ielmc),L,NR,yearIJ%I*1000+yearIJ%J/24.
!  if(abs(errC)>1.e-4)stop
!  if(I>=140)write(434,*)'mss1',mass1(ielmc),mass2(ielmc),mass3(ielmc),RootMycoNonst4Thick_Oltd(ielmc)  

  Root1stNetGrowthElms(:)  = RootMycoNonst4Thick_Oltd(:)

  RootMyco1stSinkC_rpvr(L,NR,NZ)=RCO2T1st_Oltd+RootMycoNonst4Grow_Oltd(ielmc)
  
  !update negative growth due to remobilization
  if(-RCO2PotGroth_Oltd.GT.0._r8)then
    !remobilization may turn Root1stNetGrowthElms negative
    call RemobilizePrimeActRoots(N,L,NZ,NR,RCO2PotGroth_Oltd,RCO2PotGroth_OUltd,Root1stcylc,RCO2T1st_OUltd,&
      RCO2T1st_Oltd,Root1stNetGrowthElms,litrflx,RootTipZone=.false.)
    !  
    DO NE=1,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)+Root1stcylc(NE))
    ENDDO 
  endif  
  !
  !if primary root negative growth is negative, offset negative growth by withdrawing from secondary roots, 
  IF(Root1stNetGrowthElms(ielmc).LT.0.0_r8)THEN     
    !negative growth may be offset by 2nd root death     
    call Withdraw2ndRoots(N,NZ,L,NR,Root1stNetGrowthElms,litrflx)    
  ENDIF
  
  !apply the structrual growth
  DO NE=1,NumPlantChemElms
    RootMyco1stElm_raxs(NE,NR,NZ)         = RootMyco1stElm_raxs(NE,NR,NZ)+Root1stNetGrowthElms(NE)
    Root1stActStructElms_rpvr(NE,L,NR,NZ) = Root1stActStructElms_rpvr(NE,L,NR,NZ)+Root1stNetGrowthElms(NE)
    RootMyco1stStrutElms_rpvr(NE,L,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)+Root1stNetGrowthElms(NE)
  ENDDO
  !
  !maintenance deficit is still more than maitenance paid from remobilization
  
  IF(llignification .and. -RCO2PotGroth_Oltd.GT.ZERO4Groth_pft(NZ))then  
    !remobilize from the lignified roots    
    dRCO2MaintDef_Oltd  = -RCO2PotGroth_Oltd
    dRCO2MaintDef_OUltd = -RCO2PotGroth_OUltd
    call RemobilLignifCoarseRoots(L,NZ,NR,dRCO2MaintDef_OUltd,dRCO2MaintDef_Oltd,RCO2T1st_OUltd,RCO2T1st_Oltd,litrflx)
  endif   
  !
  RootRespPotent_pvr(N,L,NZ) = RootRespPotent_pvr(N,L,NZ)+RCO2T1st_OUltd
  RootCO2EmisPot_pvr(N,L,NZ) = RootCO2EmisPot_pvr(N,L,NZ)+RCO2T1st_Oltd
  RootCO2Autor_pvr(N,L,NZ)   = RootCO2Autor_pvr(N,L,NZ)-RCO2T1st_Oltd  
  RCO2flxt                   = RCO2flxt-RCO2T1st_Oltd  
  !
  dRThickening=Root1stRadius_rpvr(L,NR,NZ)
  Root1stRadius_rpvr(L,NR,NZ)=sqrt(RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)*CRootActVolPerMassC_pft(NZ) &
    /(PiCON*Root1stLenPP_rpvr(L,NR,NZ)*Num1stAxesPerStructRootXPOP_pft(NZ)))

  !the layer to be pushed into dead zone
  dRThickening=Root1stRadius_rpvr(L,NR,NZ)-dRThickening 
  if(.not.llignification)call CRootsLive2DeadTransformation(yearIJ,NR,L,NZ,Root1stRadius_rpvr(L,NR,NZ),dRThickening)

  litrflxt=litrflxt+litrflx
!  mass_finale=0._r8;mass1=0._r8;mass2=0._r8;mass3=0._r8
!  DO L1=max(1,L-1),L
!    DO NE=1,NumPlantChemElms
!      mass_finale(NE) = mass_finale(NE) +RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      mass1(NE)=mass1(NE)+RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)
!      DO N1=1,Myco_pft(NZ)
!        mass_finale(NE) = mass_finale(NE) +RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!        mass2(NE)=mass2(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)
!        mass3(NE)=mass3(NE)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO
!  errC=mass_finale(ielmn)-mass_inital(ielmn)+litrflx(ielmn)
!  write(490,*)yearIJ%I*1000+yearIJ%J/24.,subname,errC,mass_finale(ielmn),mass_inital(ielmn),litrflx(ielmn)
!  write(433,*)yearIJ%I*1000+yearIJ%J/24.,RootMycoNonstElms_rpvr(ielmc,N,L,NZ),mass_finale(ielmc),'thick',errC,'co2',RCO2T1st_Oltd,litrflxt(ielmc),&
!    'groth',Root1stNetGrowthElms(ielmc),RootMycoNonst4Grow_Oltd(ielmc),L,NR
!  if(abs(errC)>1.e-4)call endrun(trim(mod_filename)//' at line',__LINE__)           
!  write(433,*)yearIJ%I*1000+yearIJ%J/24.,L,NR,mass_finale(ielmn)-mass_inital(ielmn)+litrflxt(ielmn),mass_finale(ielmn),mass_inital(ielmn)
  call PrintInfo('end '//subname)  
  end associate
  end subroutine SecondaryGrowthZone
!----------------------------------------------------------------------------------------------------
  subroutine CRootsLive2DeadTransformation(yearIJ,NR,L,NZ,Radius,dRThickening)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NR,L,NZ
  real(r8), intent(in):: Radius !root radius [m]
  real(r8),intent(in) :: dRThickening
  REAL(R8) :: Act2DeadElms(NumPlantChemElms)
  real(r8) :: Act2DeadElms1,Act2DeadElms2
  integer :: NE
  real(r8), parameter :: RL2D = 0.7_r8/(365._r8*24._r8) !transformation rate from live to dead coarse roots 1/hr
  CHARACTER(LEN=*), parameter :: subname='CRootsLive2DeadTransformation'
  real(r8) :: RL2DZ
  associate(                                                                &
    RootMyco1stElm_raxs         => plt_biom%RootMyco1stElm_raxs          ,& !inout :axis primary-root structural elements, [g d-2]
    RootMyco1stStrutElms_rpvr    => plt_biom%RootMyco1stStrutElms_rpvr     ,& !input  :root layer element primary axes, [g d-2]
    Root1stActStructElms_rpvr    => plt_biom%Root1stActStructElms_rpvr     ,& !inoput :root layer active zone element in primary axes, [g d-2]
    Root1stLigStructElms_rpvr    => plt_biom%Root1stLigStructElms_rpvr     ,& !inoput :root layer lignified zone element in primary axes, [g d-2]    
    rECLiveCRoot_pft             => plt_allom%rECLiveCRoot_pft             ,& !input  element:C ratio of live coarse root, [gE gC-1]
    rECDeadCRoot_pft             => plt_allom%rECDeadCRoot_pft             ,& !input  element:C ratio of dead coarse root, [gE gC-1]
    RootMycoNonstElms_rpvr       => plt_biom%RootMycoNonstElms_rpvr         & !inoput :root layer nonstructural element, [g d-2]    
  )

  ! do live to dead transformation, the following is an approximation for the lignification process, but does not respire CO2 as actual lignification
  if(dRThickening.LE.0._r8 .or. Radius .LE. dmax)return
  RL2DZ=AMIN1(2._r8*(Radius-dmax)/(dmax*(2._r8*Radius-dmax))*dRThickening,1._R8)

  Act2DeadElms(ielmc)                      = Root1stActStructElms_rpvr(ielmc,L,NR,NZ)*RL2DZ
  Root1stActStructElms_rpvr(ielmc,L,NR,NZ) = AZMAX1(Root1stActStructElms_rpvr(ielmc,L,NR,NZ)-Act2DeadElms(ielmc))
  Root1stLigStructElms_rpvr(ielmc,L,NR,NZ) = Root1stLigStructElms_rpvr(ielmc,L,NR,NZ)+Act2DeadElms(ielmc)
  
  !translocation of N/P from dead to nonstructural
  DO NE=2,NumPlantChemElms
    Act2DeadElms1 = AMIN1(Act2DeadElms(ielmc)*rECLiveCRoot_pft(NE,NZ),Root1stActStructElms_rpvr(NE,L,NR,NZ),&
                    RootMyco1stStrutElms_rpvr(NE,L,NR,NZ))
    Act2DeadElms2 = Act2DeadElms1*rECDeadCRoot_pft(NE,NZ)/rECLiveCRoot_pft(NE,NZ)
    
    Act2DeadElms(NE)                         = (Act2DeadElms1-Act2DeadElms2)
    Root1stActStructElms_rpvr(NE,L,NR,NZ)    = AZMAX1(Root1stActStructElms_rpvr(NE,L,NR,NZ)-Act2DeadElms1)
    Root1stLigStructElms_rpvr(NE,L,NR,NZ)    = AZMAX1(Root1stLigStructElms_rpvr(NE,L,NR,NZ)+Act2DeadElms2)

    RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)    = AZMAX1(RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)-Act2DeadElms(NE))
    RootMyco1stElm_raxs(NE,NR,NZ)          = RootMyco1stElm_raxs(NE,NR,NZ)-Act2DeadElms(NE)
    RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)+Act2DeadElms(NE)
  ENDDO

  end associate
  end subroutine CRootsLive2DeadTransformation

!----------------------------------------------------------------------------------------------------

  subroutine RemobilLignifCoarseRoots(L,NZ,NR,dRCO2MaintDef_OUltd,dRCO2MaintDef_Oltd,RCO2T1st_OUltd,RCO2T1st_Oltd,litrflxt)
  !
  !Description:
  !remobilization of the inactive zone of coarse roots under carbon starvation
  implicit none
  integer, intent(in) :: L,NZ,NR
  real(r8), intent(in) :: dRCO2MaintDef_OUltd,dRCO2MaintDef_Oltd    !>0, there is maintenance deficit
  real(r8), intent(inout) :: RCO2T1st_OUltd,RCO2T1st_Oltd
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)
  real(r8) :: CCC,CNC,CPC
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: Frac2Senes1
  real(r8) :: dsenecE, dlitrfall
  integer  :: NE,M, N
  real(r8) :: Root1stcylc(NumPlantChemElms),Root1stKill(NumPlantChemElms)
  real(r8) :: RCO2MaintDef1stStruct_Oltd,RCO2MDef1stStructPayoff_OUltd
  real(r8) :: RCCE(NumPlantChemElms),Root1stStrutRemob(NumPlantChemElms)

  character(len=*), parameter :: subname = 'RemobilLignifCoarseRoots'

  associate(                                                                &
    RootNonstructElmConc_pvr     => plt_biom%RootNonstructElmConc_pvr      ,& !input  :root layer nonstructural C concentration, [g g-1]
    RootMyco1stStrutElms_rpvr    => plt_biom%RootMyco1stStrutElms_rpvr     ,& !input  :root layer element primary axes, [g d-2]    
    RootMyco1stElm_raxs          => plt_biom%RootMyco1stElm_raxs           ,& !inoput :root C biomass of primary axes, [g d-2]    
    ZERO4Groth_pft               => plt_biom%ZERO4Groth_pft                ,& !input  :threshold zero for plang growth calculation, [-]
    ZERO                         => plt_site%ZERO                          ,& !input  :threshold zero for numerical stability, [-]  
    k_woody_comp                 => pltpar%k_woody_comp                    ,& !input  :woody litter complex id
    k_fine_comp                  => pltpar%k_fine_comp                     ,& !input  :fine litter complex id    
    istalk                       => pltpar%istalk                          ,& !input  :group id of coarse woody litter    
    RootMycoNonstElms_rpvr       => plt_biom%RootMycoNonstElms_rpvr        ,& !inoput :root layer nonstructural element, [g d-2]    
    PlantElmAllocMat4Litr        => plt_soilchem%PlantElmAllocMat4Litr     ,& !input  :litter kinetic fraction, [-]
    FracWoodStalkElmAlloc2Litr   => plt_allom%FracWoodStalkElmAlloc2Litr   ,& !input  :stalk element woody/fine fraction in biomass,[-]
    RAutoRootO2Limter_rpvr       => plt_rbgc%RAutoRootO2Limter_rpvr        ,& !input  :O2 constraint to root respiration (0-1), [-]    
    LitrfallElms_pvr             => plt_bgcr%LitrfallElms_pvr              ,& !inoput  :plant LitrFall element, [g d-2 h-1]     
    Root1stLigStructElms_rpvr    => plt_biom%Root1stLigStructElms_rpvr     ,& !inoput :root layer lignified zone element in primary axes, [g d-2]
    iPlantCalendar_brch          => plt_pheno%iPlantCalendar_brch          ,& !input  :plant growth stage, [-]    
    MainBranchNum_pft            => plt_morph%MainBranchNum_pft             & !input  :number of main branch,[-]    
  )
  call PrintInfo('beg '//subname)
  N=ipltroot

  IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0 &
    .AND. RootNonstructElmConc_pvr(ielmc,N,L,NZ).GT.ZERO)THEN
    CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmn,N,L,NZ), &
      RootNonstructElmConc_pvr(ielmn,N,L,NZ)+RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CNKI), &
      safe_adb(RootNonstructElmConc_pvr(ielmp,N,L,NZ),RootNonstructElmConc_pvr(ielmp,N,L,NZ) &
        +RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CPKI)))
    CNC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmc,N,L,NZ), &
      RootNonstructElmConc_pvr(ielmc,N,L,NZ)+RootNonstructElmConc_pvr(ielmn,N,L,NZ)/CNKI)))
    CPC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmc,N,L,NZ), &
      RootNonstructElmConc_pvr(ielmc,N,L,NZ)+RootNonstructElmConc_pvr(ielmp,N,L,NZ)/CPKI)))
  ELSE
    CCC=0._r8
    CNC=0._r8
    CPC=0._r8
  ENDIF
  RCCC        = RCCZR+CCC*RCCYR
  RCCN        = CNC*RCCXR
  RCCP        = CPC*RCCQR
  RCCE(ielmc) = AMIN1(RCCC,1._R8)
  RCCE(ielmn) = AMIN1(RCCN+(1.0_r8-RCCN)*RCCC,1._r8)
  RCCE(ielmp) = AMIN1(RCCP+(1.0_r8-RCCP)*RCCC,1._r8)

  !insufficient energy for maintenance
  IF(dRCO2MaintDef_OUltd.GT.0.0_r8)THEN
    IF(dRCO2MaintDef_OUltd.LT.Root1stLigStructElms_rpvr(ielmc,L,NR,NZ)*RCCE(ielmc))THEN
      RCO2MDef1stStructPayoff_OUltd=dRCO2MaintDef_OUltd
    ELSE
      RCO2MDef1stStructPayoff_OUltd=AZMAX1(Root1stLigStructElms_rpvr(ielmc,L,NR,NZ)*RCCE(ielmc))
    ENDIF    
  ELSE
    RCO2MDef1stStructPayoff_OUltd=0._r8
  ENDIF

  IF(dRCO2MaintDef_Oltd.GT.0.0_r8)THEN
    !maintenance deficit is less than remobilizable C.
    IF(dRCO2MaintDef_Oltd.LT.Root1stLigStructElms_rpvr(ielmc,L,NR,NZ)*RCCE(ielmc))THEN
      RCO2MaintDef1stStruct_Oltd=dRCO2MaintDef_Oltd
    ELSE
      RCO2MaintDef1stStruct_Oltd=AZMAX1(Root1stLigStructElms_rpvr(ielmc,L,NR,NZ)*RCCE(ielmc))*RAutoRootO2Limter_rpvr(N,L,NZ)
    ENDIF
  ELSE
    RCO2MaintDef1stStruct_Oltd=0._r8
  ENDIF

  IF(RCO2MaintDef1stStruct_Oltd.GT.0.0_r8 .AND. Root1stLigStructElms_rpvr(ielmc,L,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
    !remobilization upon starvation-induced root retreat
    DO NE=1,NumPlantChemElms
      Root1stStrutRemob(NE) = Root1stLigStructElms_rpvr(NE,L,NR,NZ)*RCCE(NE)
    ENDDO

    !maintenance deficit is paid by remobilization
    IF(Root1stStrutRemob(ielmc).GT.ZERO4Groth_pft(NZ))THEN
      Frac2Senes1=AZMAX1(AMIN1(1.0_r8,RCO2MaintDef1stStruct_Oltd/Root1stStrutRemob(ielmc)))
    ELSE
      Frac2Senes1=1.0_r8
    ENDIF
  ELSE
    Root1stStrutRemob(1:NumPlantChemElms)=0._r8
    Frac2Senes1=0._r8
  ENDIF

  Root1stcylc=0._r8
  if(Frac2Senes1>0._r8)then
    DO NE=1,NumPlantChemElms          
      Root1stKill(NE) = Root1stLigStructElms_rpvr(NE,L,NR,NZ)*Frac2Senes1  
      Root1stcylc(NE) = Root1stStrutRemob(NE)*Frac2Senes1
      Root1stLigStructElms_rpvr(NE,L,NR,NZ)=Root1stLigStructElms_rpvr(NE,L,NR,NZ)-Root1stKill(NE)
      RootMyco1stElm_raxs(NE,NR,NZ)         = RootMyco1stElm_raxs(NE,NR,NZ)-Root1stKill(NE)
      RootMyco1stStrutElms_rpvr(NE,L,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)-Root1stKill(NE)
    ENDDO
    !inactive zone
    !first mobilize active zone
    D6356: DO M=1,jsken
      DO NE=1,NumPlantChemElms      
            
        dsenecE      = AZMAX1(Root1stKill(NE)-Root1stcylc(NE))
        dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,istalk,M,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp)
        litrflxt(NE) = litrflxt(NE)+dlitrfall

        LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)+dlitrfall

        dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,istalk,M,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp)
        litrflxt(NE) = litrflxt(NE)+dlitrfall

        LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) +dlitrfall
      ENDDO
    ENDDO D6356
  endif
  RCO2T1st_Oltd      = RCO2T1st_Oltd+RCO2MaintDef1stStruct_Oltd
  RCO2T1st_OUltd     = RCO2T1st_OUltd+RCO2MDef1stStructPayoff_OUltd
  Root1stcylc(ielmc) = Root1stcylc(ielmc)-RCO2MaintDef1stStruct_Oltd

  DO NE=1,NumPlantChemElms
    RootMycoNonstElms_rpvr(NE,N,L,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)+Root1stcylc(NE))
  ENDDO 

  call PrintInfo('end '//subname)
  end associate
  end subroutine RemobilLignifCoarseRoots

!----------------------------------------------------------------------------------------------------

  function stress4lignif(NZ,RAutoRootO2Limter,Nutstress4GrossResp,fRootGrowPSISense)result(f)
  !
  !stress enhancement for coarse root lignification
  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: RAutoRootO2Limter   !hypoxia stress [0-1]
  real(r8), intent(in) :: Nutstress4GrossResp !nutrient limitation [0-1] 
  real(r8), intent(in) :: fRootGrowPSISense   !moisture stress [0-1]
  real(r8) :: f   !by default stress enhancement is set to zero

  f = 0._r8

  end function stress4lignif

!----------------------------------------------------------------------------------------------------

  subroutine KillCoarseRoots(L,NZ,NR,Root1stNetGrowthElms,litrflxt)
  !
  !Description:
  !remobilization of the inactive zone of coarse roots under carbon starvation
  implicit none
  integer, intent(in) :: L,NZ,NR
  real(r8), intent(in) :: Root1stNetGrowthElms(NumPlantChemElms)
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)
  real(r8) :: dsenecE, dlitrfall
  integer  :: NE,M,N,L1
  real(r8) :: Root1stKill(NumPlantChemElms),dRootMycoNonst2Litr(NumPlantChemElms)
  real(r8) :: dRootMyco2ndst2Litr(NumPlantChemElms)
  real(r8) :: dwoodyE,dfineE1,dfineE2

  character(len=*), parameter :: subname = 'KillCoarseRoots'

  associate(                                                                &
    RootNonstructElmConc_pvr    => plt_biom%RootNonstructElmConc_pvr     ,& !input  :root layer nonstructural C concentration, [g g-1]
    NU                           => plt_site%NU                            ,& !input  :current soil surface layer number, [-]     
    ZERO4Groth_pft               => plt_biom%ZERO4Groth_pft                ,& !input  :threshold zero for plang growth calculation, [-]
    ZERO                         => plt_site%ZERO                          ,& !input  :threshold zero for numerical stability, [-]  
    icwood                       => pltpar%icwood                          ,& !input  :group id of coarse woody litter
    iroot                        => pltpar%iroot                           ,& !input  :group id of plant root litter    
    Myco_pft                     => plt_morph%Myco_pft                     ,& !input  :mycorrhizal type (no or yes),[-]      
    NRoot1stTipLay_raxes         => plt_morph%NRoot1stTipLay_raxes         ,& !input  :maximum soil layer number for root axes, [-]    
    Root2ndLen_rpvr              => plt_morph%Root2ndLen_rpvr              ,& !inoput :root layer length secondary axes, [m d-2]    
    inonstruct                   => pltpar%inonstruct                      ,& !input  :group id of plant nonstructural litter    
    FracRootElmAllocm            => plt_allom%FracRootElmAllocm            ,& !input  :C woody/soft fraction in root,[-]     
    k_woody_comp                 => pltpar%k_woody_comp                    ,& !input  :woody litter complex id
    Root1stLenPP_rpvr            => plt_morph%Root1stLenPP_rpvr            ,& !input  :primary root axis length in soil layer, [m d-2]       
    k_fine_comp                  => pltpar%k_fine_comp                     ,& !input  :fine litter complex id    
    istalk                       => pltpar%istalk                          ,& !input  :group id of coarse woody litter    
    FracWoodStalkElmAlloc2Litr   => plt_allom%FracWoodStalkElmAlloc2Litr   ,& !input  :stalk element woody/fine fraction in biomass,[-]    
    PlantElmAllocMat4Litr        => plt_soilchem%PlantElmAllocMat4Litr     ,& !input  :litter kinetic fraction, [-]    
    Root1stActStructElms_rpvr    => plt_biom%Root1stActStructElms_rpvr     ,& !inoput :Root layer primary axes Active zone structrual element, [g d-2]    
    RootMycoNonstElms_rpvr       => plt_biom%RootMycoNonstElms_rpvr        ,& !inoput :root layer nonstructural element, [g d-2]    
    RootMyco1stElm_raxs          => plt_biom%RootMyco1stElm_raxs           ,& !inoput :root layer primary axes structrual element, [g d-2]    
    RootMyco2ndStrutElms_rpvr    => plt_biom%RootMyco2ndStrutElms_rpvr     ,& !inoput :root layer element secondary axes, [g d-2]    
    RootMyco1stStrutElms_rpvr    => plt_biom%RootMyco1stStrutElms_rpvr     ,& !inoput  :root layer element primary axes, [g d-2]     
    LitrfallElms_pvr             => plt_bgcr%LitrfallElms_pvr              ,& !inoput  :plant LitrFall element, [g d-2 h-1]     
    Root1stLigStructElms_rpvr    => plt_biom%Root1stLigStructElms_rpvr      & !inoput :root layer lignified zone element in primary axes, [g d-2]
  )
  call PrintInfo('beg '//subname)

  IF(-Root1stNetGrowthElms(ielmc).GT.ZERO4Groth_pft(NZ))THEN  
    !
    !kill the root in current layer and layer below
    DO L1=L,NRoot1stTipLay_raxes(NR,NZ)+1
      DO NE=1,NumPlantChemElms          
        Root1stKill(NE)                       = Root1stLigStructElms_rpvr(NE,L1,NR,NZ)
        Root1stLigStructElms_rpvr(NE,L1,NR,NZ) = 0._r8
        RootMyco1stElm_raxs(NE,NR,NZ)         = RootMyco1stElm_raxs(NE,NR,NZ)-Root1stKill(NE)
        RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)-Root1stKill(NE)
      ENDDO
      !
      !inactive zone
      D6356: DO M=1,jsken
        DO NE=1,NumPlantChemElms                  
          dsenecE      = AZMAX1(Root1stKill(NE))
          dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,istalk,M,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp)
          litrflxt(NE) = litrflxt(NE)+dlitrfall
          LitrfallElms_pvr(NE,M,k_woody_comp,L1,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,L1,NZ)+dlitrfall

          dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,istalk,M,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp)
          litrflxt(NE) = litrflxt(NE)+dlitrfall
          LitrfallElms_pvr(NE,M,k_fine_comp,L1,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L1,NZ) +dlitrfall
        ENDDO
      ENDDO D6356      
      !
      !active zone
      DO NE=1,NumPlantChemElms          
        Root1stKill(NE)                       = Root1stActStructElms_rpvr(NE,L1,NR,NZ) 
        Root1stActStructElms_rpvr(NE,L1,NR,NZ) = 0._r8
        RootMyco1stElm_raxs(NE,NR,NZ)         = RootMyco1stElm_raxs(NE,NR,NZ)-Root1stKill(NE)
        RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ)-Root1stKill(NE)
      ENDDO
      !
      !active zone
      D6357: DO M=1,jsken
        DO NE=1,NumPlantChemElms                  
          dsenecE      = AZMAX1(Root1stKill(NE))
          dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,istalk,M,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp)          
          litrflxt(NE) = litrflxt(NE)+dlitrfall
          LitrfallElms_pvr(NE,M,k_woody_comp,L1,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,L1,NZ)+dlitrfall

          dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,istalk,M,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp)
          litrflxt(NE) = litrflxt(NE)+dlitrfall
          LitrfallElms_pvr(NE,M,k_fine_comp,L1,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L1,NZ) +dlitrfall
        ENDDO
      ENDDO D6357
      Root1stLenPP_rpvr(L1,NR,NZ) = 0._r8

      !remove secondary roots
      DO N=1,Myco_pft(NZ)
        D6450: DO M=1,jsken
          D64511: DO NE=1,NumPlantChemElms
            dsenecE = AZMAX1(RootMyco2ndStrutElms_rpvr(NE,N,L1,NR,NZ))
            dwoodyE = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dsenecE*FracRootElmAllocm(NE,k_woody_comp)
            dfineE1 = PlantElmAllocMat4Litr(NE,iroot,M,NZ)*dsenecE*FracRootElmAllocm(NE,k_fine_comp)  

            LitrfallElms_pvr(NE,M,k_woody_comp,L1,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,L1,NZ)+dwoodyE
            LitrfallElms_pvr(NE,M,k_fine_comp,L1,NZ)  = LitrfallElms_pvr(NE,M,k_fine_comp,L1,NZ)+dfineE1
            dRootMyco2ndst2Litr(NE)                   = dRootMyco2ndst2Litr(NE)+dwoodyE+dfineE1
            litrflxt(NE)                              = litrflxt(NE)+dwoodyE+dfineE1

            dfineE2=PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*RootMycoNonstElms_rpvr(NE,N,L1,NZ)
            LitrfallElms_pvr(NE,M,k_fine_comp,L1,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,L1,NZ)+dfineE2
            dRootMycoNonst2Litr(NE)                  = dRootMycoNonst2Litr(NE)+dfineE2
            litrflxt(NE)                             = litrflxt(NE)+dfineE2
          ENDDO D64511   
        ENDDO D6450
        
        DO NE=1,NumPlantChemElms  
          RootMyco2ndStrutElms_rpvr(NE,N,L1,NR,NZ) = AZMAX1(RootMyco2ndStrutElms_rpvr(NE,N,L1,NR,NZ)-dRootMyco2ndst2Litr(NE))
          RootMycoNonstElms_rpvr(NE,N,L1,NZ)       = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L1,NZ)-dRootMycoNonst2Litr(NE))
        ENDDO
        Root2ndLen_rpvr(N,L1,NR,NZ)=0._r8
      ENDDO
    ENDDO  
  endif
  NRoot1stTipLay_raxes(NR,NZ)=MAX(L-1,NU)
  call PrintInfo('end '//subname)
  end associate
  end subroutine KillCoarseRoots

!----------------------------------------------------------------------------------------------------
  subroutine RootElongationWaterFunc(N,L,NZ,Lelonz,fdLext_mech,RespElongWatSens)
  !
  !Description:
  !water dependence for root elongation
  implicit none
  integer, intent(in) :: N,L,NZ
  real(r8), intent(in) :: Lelonz                          !length of the elongation zone. [m]       
  real(r8), intent(out) :: fdLext_mech                    !soil compaction scalar on pressure-allowed root penetration/elongation [h m-1]
  real(r8), optional, intent(out) :: RespElongWatSens     !respiration function of root water potential, [-]

  !--- local variables ---
  real(r8), parameter :: phi_lmin = 0.1_r8                !minimum root elastic modulus normalizer [MPa-1 h], i.e. F' in Rickman et al. 1992, Grant, 1992
  real(r8), parameter :: phi_lmax = 1._r8                 !maximum root elastic modulus normalizer [MPa-1 h]
  real(r8), parameter :: Ksens=0.8_r8                     !Haf saturation constant, [MPa]

  real(r8)  :: RootElongWatSens                           !water function for root elongation, [-]
  real(r8) :: V0,V0p                                      !dispaceable soil volume due to root growth [m3]  
  real(r8) :: rv0,Radius_Impact                           !impact radius, [m]
  real(r8), parameter :: r0 =0.5e-3_r8                    !reference radius [m]
  real(r8), parameter :: Vr0=PICON*1.e-2_r8*r0**2         !reference root volume, [m3] 
  real(r8) :: phi_l,SoilResit4RootPentration
  character(len=*), parameter :: subname='RootElongationWaterFunc'

  associate(                                                                    &
    PSIRootTurg_vr               => plt_ew%PSIRootTurg_vr                      ,& !input  :root turgor water potential, [Mpa]
    iPlantRootProfile_pft        => plt_pheno%iPlantRootProfile_pft            ,& !input  :plant growth type (vascular, non-vascular),[-]
    SoilBulkModulus4RootPent_vr  => plt_soilchem%SoilBulkModulus4RootPent_vr    & !input  :elastic modulus of the undisturbed soil, [MPa]
  )
!     SoilBulkModulus4RootPent_vr,SoilResit4PrimRootPentration=soil resistance to primary root penetration (MPa)
!     RRAD1=primary root radius
  !eq.(3) in Rickman et al (1992), Calculating daily root length density profiles by applying elastic theory to agricultural soils
  
  phi_l=phi_lmin+(phi_lmax-phi_lmin)/(1._r8+(SoilBulkModulus4RootPent_vr(L)/Ksens)**2)
  !eq. (10) in Rickman et al. (1992)
  !V0 was obtained by regression, assuming r0=0.05 cm=0.5mm and elongation zone 1cm.
  V0=(0.0327_r8+0.0055_r8*SoilBulkModulus4RootPent_vr(L))*1.e-4_r8

  !!aka SoilResit4RootPentration=E'*dt/V0 in eq. (4) of Rickman et al. (1992), [h m-1]
  SoilResit4RootPentration = phi_l*SoilBulkModulus4RootPent_vr(L)*Vr0/V0

  !The lockhart equation for elongation
  RootElongWatSens = AMIN1(1.0_r8,AZMAX1(PSIRootTurg_vr(N,L,NZ)-TurgPSIMin4OrganExtens))  
  if(present(RespElongWatSens))RespElongWatSens = fRespWatSens(RootElongWatSens,iPlantRootProfile_pft(NZ))
  !dLext_mech          = RootElongWatSens*phi_l*Lelonz
  fdLext_mech = 1._r8/(1._r8+SoilResit4RootPentration)
  end associate    
  END subroutine RootElongationWaterFunc
!----------------------------------------------------------------------------------------------------
  subroutine RootThickenWaterFunc(N,L,NZ,RootRadius,RespThickenWatSens,dR1stExp)

  implicit none
  integer, intent(in) :: N,L,NZ
  real(r8), intent(in)  :: RootRadius                     !root radius [m]
  real(r8), intent(out) :: RespThickenWatSens             !respiration function of root water potential, [-]
  real(r8), intent(out) :: dR1stExp                       !stress allowed radial expansion, [m]

  character(len=*), parameter :: subname='RootThickenWaterFunc'
  real(r8), parameter :: phi_Rmin =0.1_r8      !the baseline extensibility in hydroponics, which keeps root thin [MPa-1 h-1]
  real(r8), parameter :: phi_Rmax =0.5_r8      !The maximum extensibility in hard soil, allows thickening, [MPa-1 h-1]
  real(r8), parameter :: KSensR = 0.1_r8       !The "Sensory Threshold." The pressure at which the plant starts thickening [MPa]
  real(r8) :: phi_R,normSoilStress
  real(r8) :: RootThickenWatSens
  associate(                                                                     &
    PSIRootTurg_vr                => plt_ew%PSIRootTurg_vr                      ,& !input  :root turgor water potential, [Mpa]
    iPlantRootProfile_pft         => plt_pheno%iPlantRootProfile_pft            ,& !input  :plant growth type (vascular, non-vascular),[-]
    SoilModulus4RootRadialexp_vr  => plt_soilchem%SoilModulus4RootRadialexp_vr   & !input  :soil modulus for root radial expansion, [MPa]    
  )

  normSoilStress     = (SoilModulus4RootRadialexp_vr(L)/KSensR)**2
  phi_R              = phi_Rmin+(phi_Rmax-phi_Rmin)*(normSoilStress/(1._r8+normSoilStress))
  RootThickenWatSens = AMIN1(1.0_r8,AZMAX1(PSIRootTurg_vr(N,L,NZ)-TurgPSIMin4OrganExtens-SoilModulus4RootRadialexp_vr(L)))
  RespThickenWatSens = fRespWatSens(RootThickenWatSens,iPlantRootProfile_pft(NZ))
  dR1stExp           = phi_R*RootRadius*RootThickenWatSens
  end associate
  end subroutine RootThickenWaterFunc
!----------------------------------------------------------------------------------------------------
  subroutine Withdraw2ndRoots(N,NZ,L,NR,RootNetGrowthElms,litrflxt)
  !
  !Description
  !Due to negative growth of primary roots, secondary roots are killed to make up the metabolic deficit
  !
  implicit none
  integer, intent(in) :: N        !root type, plant/myco
  integer, intent(in) :: NZ,L,NR
  real(r8), intent(inout) :: RootNetGrowthElms(NumPlantChemElms)
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)          !litterfall due to withdraw of secondary/fine roots
  real(r8) :: litrflx(NumPlantChemElms)
  integer :: LX,LL,NE,M,N1
  real(r8) :: GRTWTM     !negative primary root C growth
  real(r8) :: FSNCM,FSNCP
  real(r8) :: dRootMyco2ndst2Litr(NumPlantChemElms)
  real(r8) :: dRootMycoNonst2Litr(NumPlantChemElms)
  real(r8) :: dsenecE,dwoodyE,dfineE1,dfineE2
  character(len=*), parameter :: subname='Withdraw2ndRoots'
  real(r8) :: mass_inital(NumPlantChemElms)

  associate(                                                          &
    Myco_pft                  => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]  
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft             ,& !input  :threshold zero for plang growth calculation, [-]
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr    ,& !input  :root layer structural C, [gC d-2]
    icwood                    => pltpar%icwood                       ,& !input  :group id of coarse woody litter
    PlantElmAllocMat4Litr     => plt_soilchem%PlantElmAllocMat4Litr  ,& !input  :litter kinetic fraction, [-]
    FracRootElmAllocm         => plt_allom%FracRootElmAllocm         ,& !input  :C woody fraction in root,[-]
    iroot                     => pltpar%iroot                        ,& !input  :group id of plant root litter
    k_woody_comp              => pltpar%k_woody_comp                 ,& !input  :woody litter complex id
    k_fine_comp               => pltpar%k_fine_comp                  ,& !input  :fine litter complex id
    inonstruct                => pltpar%inonstruct                   ,& !input  :group id of plant nonstructural litter
    RootMyco2ndStrutElms_rpvr => plt_biom%RootMyco2ndStrutElms_rpvr  ,& !inoput :root layer element secondary axes, [g d-2]
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr     ,& !inoput :root layer nonstructural element, [g d-2]
    LitrfallElms_pvr          => plt_bgcr%LitrfallElms_pvr           ,& !inoput :plant LitrFall element, [g d-2 h-1]
    Root2ndLen_rpvr           => plt_morph%Root2ndLen_rpvr            & !inoput :root layer length secondary axes, [m d-2]
  )
  call PrintInfo('beg '//subname)
  LX          = MAX(1,L-1)
!  litrflx     = 0._r8
!  mass_inital = 0._r8

!  DO LL = LX, L
!    DO N1=1,Myco_pft(NZ)
!      DO NE=1,NumPlantChemElms
!        mass_inital(NE)=mass_inital(NE)+RootMycoNonstElms_rpvr(NE,N1,L1,NZ)+RootMyco2ndStrutElms_rpvr(NE,N1,L1,NR,NZ)      
!      ENDDO
!    ENDDO
!  ENDDO

  D5105: DO LL=L,LX,-1
    !make a copy of the total negative growth 
    GRTWTM=RootNetGrowthElms(ielmc)
    
    !remove secondary roots and offset some of the negative growth
    D5106: DO NE=1,NumPlantChemElms
      !root dying, offset some negative growth, but not CO2, assuming CO2 can only be offset from reserve and primary roots
      IF(RootNetGrowthElms(NE).LT.0.0_r8)THEN    
        !negative growth is less than total fine roots in layer    
        IF(-RootNetGrowthElms(NE).LE.RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ))THEN
          !root partially dead, length thrinked
          if(NE.EQ.ielmc)then
            Root2ndLen_rpvr(N,LL,NR,NZ)=Root2ndLen_rpvr(N,LL,NR,NZ)*(1._r8+RootNetGrowthElms(NE)/RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ))
          endif   
          RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ) = RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)+RootNetGrowthElms(NE)
          RootNetGrowthElms(NE)                    = 0._r8
        ELSE
          !root dead completely
          if(NE.EQ.ielmc)Root2ndLen_rpvr(N,LL,NR,NZ)=0._r8
          RootNetGrowthElms(NE)                    = RootNetGrowthElms(NE)+RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ)
          RootMyco2ndStrutElms_rpvr(NE,N,LL,NR,NZ) = 0._r8
        ENDIF
      ENDIF
    ENDDO D5106    
    !
    !     CONCURRENT LOSS OF MYCORRHIZAE AND NODULES WITH LOSS
    !     OF SECONDARY ROOTS
    !
    !     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
    !     FSNCM,FSNCP=fraction of mycorrhizal structural,nonstructural C to be remobilized
    !     WTRTL=active root C mass
    !     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
    !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
    !     WTRT2,WTRT2N,WTRT2P=mycorrhizal C,N,P mass
    !     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
    !     Root2ndLen_rpvr=mycorrhizal length
    !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in mycorrhizae
    !mycorrhizae is turned completely into litter, as it is assumed as an extension of fine roots
    IF(GRTWTM.LT.0.0_r8 .and. Myco_pft(NZ).GT.1)THEN

      !fraction of mycorrhizae senescenced
      IF(RootMyco2ndStrutElms_rpvr(ielmc,imycorrhz,LL,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
        FSNCM=AMIN1(1.0_r8,ABS(GRTWTM)/RootMyco2ndStrutElms_rpvr(ielmc,imycorrhz,LL,NR,NZ))
      ELSE
        FSNCM=1.0_r8
      ENDIF

      !fraction of root senescenced
      IF(RootMycoActiveBiomC_pvr(ipltroot,LL,NZ).GT.ZERO4Groth_pft(NZ))THEN
        FSNCP=AMIN1(1.0_r8,ABS(GRTWTM)/RootMycoActiveBiomC_pvr(ipltroot,LL,NZ))
      ELSE
        FSNCP=1.0_r8
      ENDIF

      dRootMyco2ndst2Litr = 0._r8
      dRootMycoNonst2Litr = 0._r8
      
      D6450: DO M=1,jsken
        D64511: DO NE=1,NumPlantChemElms
          !mycorrhizal structrual biomass senescence
          dsenecE = FSNCM*AZMAX1(RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ))
          dwoodyE = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*dsenecE*FracRootElmAllocm(NE,k_woody_comp)
          dfineE1 = PlantElmAllocMat4Litr(NE,iroot,M,NZ)*dsenecE*FracRootElmAllocm(NE,k_fine_comp)  

          LitrfallElms_pvr(NE,M,k_woody_comp,LL,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,LL,NZ)+dwoodyE
          LitrfallElms_pvr(NE,M,k_fine_comp,LL,NZ)  = LitrfallElms_pvr(NE,M,k_fine_comp,LL,NZ)+dfineE1
          dRootMyco2ndst2Litr(NE)                   = dRootMyco2ndst2Litr(NE)+dwoodyE+dfineE1

          !mycorrhizae nonstrucal biomass loss
          litrflxt(NE)=litrflxt(NE)+dwoodyE+dfineE1
          dfineE2=PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*FSNCP*RootMycoNonstElms_rpvr(NE,imycorrhz,LL,NZ)
          LitrfallElms_pvr(NE,M,k_fine_comp,LL,NZ) = LitrfallElms_pvr(NE,M,k_fine_comp,LL,NZ)+dfineE2
          dRootMycoNonst2Litr(NE)                  = dRootMycoNonst2Litr(NE)+dfineE2

          litrflxt(NE)=litrflxt(NE)+dfineE2              
        ENDDO D64511   
      ENDDO D6450
      
      DO NE=1,NumPlantChemElms  
        RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ) = AZMAX1(RootMyco2ndStrutElms_rpvr(NE,imycorrhz,LL,NR,NZ)-dRootMyco2ndst2Litr(NE))
        RootMycoNonstElms_rpvr(NE,imycorrhz,LL,NZ)       = AZMAX1(RootMycoNonstElms_rpvr(NE,imycorrhz,LL,NZ)-dRootMycoNonst2Litr(NE))
      ENDDO
      Root2ndLen_rpvr(imycorrhz,LL,NR,NZ)=AZMAX1(Root2ndLen_rpvr(imycorrhz,LL,NR,NZ))*(1.0_r8-FSNCM)
    ENDIF
  ENDDO D5105
  call PrintInfo('end '//subname)
  end associate
  end subroutine Withdraw2ndRoots

!----------------------------------------------------------------------------------------------------
  subroutine RemobilizePrimeActRoots(N,L,NZ,NR,RCO2PotGroth_Oltd,RCO2PotGroth_OUltd,Root1stcylc,RCO2T1st_OUltd,&
    RCO2T1st_Oltd,Root1stNetGrowthElms,litrflxt,RootTipZone)
  ! 
  !Description:
  !primary roots remobilization
  implicit none
  integer , intent(in)    :: N,L,NZ,NR
  real(r8), intent(inout) :: RCO2PotGroth_Oltd                        !O2-limited potential growth respiration
  real(r8), intent(inout) :: RCO2PotGroth_OUltd                       !O2-unlimited potential growth respiration
  real(r8), intent(inout) :: RCO2T1st_Oltd                            !CO2 flux from 1st roots, with O-limitation  
  real(r8), intent(inout) :: RCO2T1st_OUltd                           !CO2 flux from 1st roots, without O-limitation
  real(r8), intent(inout) :: litrflxt(NumPlantChemElms)
  real(r8), intent(out)   :: Root1stcylc(NumPlantChemElms)            !remobilization flux recycled into reserve
  real(r8), intent(inout) :: Root1stNetGrowthElms(NumPlantChemElms)   !net growth after acouting starvation-kill
  logical,  intent(in)    :: RootTipZone                              !dealing with root tip zone or not?
  integer :: M,NE
  real(r8) :: CCC,CNC,CPC
  real(r8) :: Root1stStrutRemob(NumPlantChemElms)
  real(r8) :: Frac2Senes1         !fraction of biomass under starvation-killing for remobilization
  real(r8) :: dlitrfall
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: RCO2MDef1stStructPayoff_OUltd,RCO2MDef1stStructPayoff_Oltd             !maintenance deficit payoff by remobilization of structural carbon without/with O2 limitation
  real(r8) :: dsenecE,RCCE(NumPlantChemElms)
  real(r8) :: Root1stKill(NumPlantChemElms)
  CHARACTER(LEN=*), parameter :: subname='RemobilizePrimeActRoots'

! begin_execution
  associate(                                                           &
    RootNonstructElmConc_pvr  => plt_biom%RootNonstructElmConc_pvr   ,& !input  :root layer nonstructural C concentration, [g g-1]
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft              ,& !input  :threshold zero for plang growth calculation, [-]
    ZERO                       => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]
    icwood                     => pltpar%icwood                        ,& !input  :group id of coarse woody litter
    istalk                     => pltpar%istalk                        ,& !input  :group id of coarse woody litter
    k_woody_comp               => pltpar%k_woody_comp                  ,& !input  :woody litter complex id
    k_fine_comp                => pltpar%k_fine_comp                   ,& !input  :fine litter complex id
    iroot                      => pltpar%iroot                         ,& !input  :group id of plant root litter
    RootMyco1stElm_raxs        => plt_biom%RootMyco1stElm_raxs         ,& !inoput :root layer primary axes structrual element, [g d-2]    
    Root1stActStructElms_rpvr  => plt_biom%Root1stActStructElms_rpvr   ,& !input  :root layer active zone element in primary axes, [g d-2]
    FracRootElmAllocm          => plt_allom%FracRootElmAllocm          ,& !input  :root element woody/fine fraction in biomass,[-]
    FracWoodStalkElmAlloc2Litr => plt_allom%FracWoodStalkElmAlloc2Litr ,& !input  :stalk element woody/fine fraction in biomass,[-]
    PlantElmAllocMat4Litr      => plt_soilchem%PlantElmAllocMat4Litr   ,& !input  :litter kinetic fraction, [-]
    RAutoRootO2Limter_rpvr     => plt_rbgc%RAutoRootO2Limter_rpvr      ,& !input  :O2 constraint to root respiration (0-1), [-]
    iPlantCalendar_brch        => plt_pheno%iPlantCalendar_brch        ,& !input  :plant growth stage, [-]
    MainBranchNum_pft          => plt_morph%MainBranchNum_pft          ,& !input  :number of main branch,[-]
    LitrfallElms_pvr           => plt_bgcr%LitrfallElms_pvr             & !inoput :plant LitrFall element, [g d-2 h-1]
  )
  !
  !     PRIMARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
  !     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
  !     PRIMARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT LitrFall
  !
  !     iPlantCalendar_brch(ipltcal_Emerge,=emergence date
  !     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
  !     CNKI,CPKI=nonstructural N,P inhibition constant on growth
  !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
  !     RCCZR,RCCYR=min,max fractions for root C recycling
  !     RCCXR,RCCQR=max fractions for root N,P recycling
  !
  call PrintInfo('beg '//subname)
  Root1stcylc=0._r8
  IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0 &
    .AND. RootNonstructElmConc_pvr(ielmc,N,L,NZ).GT.ZERO)THEN
    CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmn,N,L,NZ), &
      RootNonstructElmConc_pvr(ielmn,N,L,NZ)+RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CNKI), &
      safe_adb(RootNonstructElmConc_pvr(ielmp,N,L,NZ),RootNonstructElmConc_pvr(ielmp,N,L,NZ) &
        +RootNonstructElmConc_pvr(ielmc,N,L,NZ)*CPKI)))
    CNC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmc,N,L,NZ), &
      RootNonstructElmConc_pvr(ielmc,N,L,NZ)+RootNonstructElmConc_pvr(ielmn,N,L,NZ)/CNKI)))
    CPC=AZMAX1(AMIN1(1.0_r8,safe_adb(RootNonstructElmConc_pvr(ielmc,N,L,NZ), &
      RootNonstructElmConc_pvr(ielmc,N,L,NZ)+RootNonstructElmConc_pvr(ielmp,N,L,NZ)/CPKI)))
  ELSE
    CCC=0._r8;CNC=0._r8;CPC=0._r8
  ENDIF
  RCCC        = RCCZR+CCC*RCCYR
  RCCN        = CNC*RCCXR
  RCCP        = CPC*RCCQR
  RCCE(ielmc) = AMIN1(RCCC,1._R8)
  RCCE(ielmn) = AMIN1(RCCN+(1.0_r8-RCCN)*RCCC,1._R8)
  RCCE(ielmp) = AMIN1(RCCP+(1.0_r8-RCCP)*RCCC,1._R8)
  !
  !     RECOVERY OF REMOBILIZABLE N,P DURING PRIMARY ROOT REMOBILIZATION
  !     DEPENDS ON ROOT NON-STRUCTURAL C:N:P
  !
  !     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
  !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
  !     Root1stStrutRemob(ielmc),Root1stStrutRemob(ielmn),Root1stStrutRemob(ielmp)=remobilization of C,N,P from senescing root
  !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
  !
  if(.not.lcoarseroot)then
    !insufficient energy for maintenance
    IF(-RCO2PotGroth_OUltd.GT.0.0_r8)THEN
      IF(-RCO2PotGroth_OUltd.LT.RootMyco1stElm_raxs(ielmc,NR,NZ)*RCCE(ielmc))THEN
        RCO2MDef1stStructPayoff_OUltd=-RCO2PotGroth_OUltd
      ELSE
        RCO2MDef1stStructPayoff_OUltd=AZMAX1(RootMyco1stElm_raxs(ielmc,NR,NZ)*RCCE(ielmc))
      ENDIF
      RCO2PotGroth_OUltd = RCO2PotGroth_OUltd+RCO2MDef1stStructPayoff_OUltd
    ELSE
      RCO2MDef1stStructPayoff_OUltd=0._r8
    ENDIF     

    !has maintenance deficit
    IF(-RCO2PotGroth_Oltd.GT.0.0_r8)THEN
      ! Use this oxygen-limited payment for both respiration and tissue withdrawal.
      RCO2MDef1stStructPayoff_Oltd = AMIN1(-RCO2PotGroth_Oltd, &
        AZMAX1(RootMyco1stElm_raxs(ielmc,NR,NZ)*RCCC))*RAutoRootO2Limter_rpvr(N,L,NZ)
      !offset some CO2
      RCO2PotGroth_Oltd=RCO2PotGroth_Oltd+RCO2MDef1stStructPayoff_Oltd
    ELSE
      RCO2MDef1stStructPayoff_Oltd=0._r8
    ENDIF
    
    IF(RCO2MDef1stStructPayoff_Oltd.GT.0.0_r8 .AND. RootMyco1stElm_raxs(ielmc,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
      !remobilization upon starvation-induced root retreat
      DO NE=1,NumPlantChemElms
        Root1stStrutRemob(NE) = RootMyco1stElm_raxs(NE,NR,NZ)*RCCE(NE)
      ENDDO

      !maintenance deficit is paid by remobilization
      IF(Root1stStrutRemob(ielmc).GT.ZERO4Groth_pft(NZ))THEN
        Frac2Senes1=AZMAX1(AMIN1(1.0_r8,RCO2MDef1stStructPayoff_Oltd/Root1stStrutRemob(ielmc)))
      ELSE
        Frac2Senes1=1.0_r8
      ENDIF
    ELSE
      Root1stStrutRemob(1:NumPlantChemElms)=0._r8
      Frac2Senes1=0._r8
    ENDIF
    
    if(Frac2Senes1.GT.0._r8)then
      DO NE=1,NumPlantChemElms          
        Root1stKill(NE) = RootMyco1stElm_raxs(NE,NR,NZ)*Frac2Senes1  
        Root1stcylc(NE) = Root1stStrutRemob(NE)*Frac2Senes1
      ENDDO
      !
      D6351: DO M=1,jsken
        DO NE=1,NumPlantChemElms          
          dsenecE      = AZMAX1(Root1stKill(NE)-Root1stcylc(NE))
          dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,icwood,M,NZ)*FracRootElmAllocm(NE,k_woody_comp)
          litrflxt(NE) = litrflxt(NE)+dlitrfall

          LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) +dlitrfall
          
          dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,iroot,M,NZ)*FracRootElmAllocm(NE,k_fine_comp)
          litrflxt(NE) = litrflxt(NE)+dlitrfall

          LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+dlitrfall
        ENDDO
      ENDDO D6351
      !remobilization reduces growth    
      DO NE=1,NumPlantChemElms
        Root1stNetGrowthElms(NE)=Root1stNetGrowthElms(NE)-Root1stKill(NE)
      ENDDO     
    endif
  else

    !insufficient energy for maintenance
    IF(-RCO2PotGroth_OUltd.GT.0.0_r8)THEN
      IF(-RCO2PotGroth_OUltd.LT.Root1stActStructElms_rpvr(ielmc,L,NR,NZ)*RCCE(ielmc))THEN
        RCO2MDef1stStructPayoff_OUltd=-RCO2PotGroth_OUltd
      ELSE
        RCO2MDef1stStructPayoff_OUltd=AZMAX1(Root1stActStructElms_rpvr(ielmc,L,NR,NZ)*RCCE(ielmc))
      ENDIF
      RCO2PotGroth_OUltd = RCO2PotGroth_OUltd+RCO2MDef1stStructPayoff_OUltd
    ELSE
      RCO2MDef1stStructPayoff_OUltd=0._r8
    ENDIF  

    !has maintenance deficit
    IF(-RCO2PotGroth_Oltd.GT.0.0_r8)THEN
      ! Use this oxygen-limited payment for both respiration and tissue withdrawal.
      RCO2MDef1stStructPayoff_Oltd = AMIN1(-RCO2PotGroth_Oltd, &
        AZMAX1(Root1stActStructElms_rpvr(ielmc,L,NR,NZ)*RCCC))*RAutoRootO2Limter_rpvr(N,L,NZ)
      !offset some CO2
      RCO2PotGroth_Oltd = RCO2PotGroth_Oltd+RCO2MDef1stStructPayoff_Oltd
    ELSE
      RCO2MDef1stStructPayoff_Oltd=0._r8
    ENDIF

    IF(RCO2MDef1stStructPayoff_Oltd.GT.0.0_r8 .AND. Root1stActStructElms_rpvr(ielmc,L,NR,NZ).GT.ZERO4Groth_pft(NZ))THEN
      !remobilization upon starvation-induced root retreat
      DO NE=1,NumPlantChemElms
        Root1stStrutRemob(NE) = Root1stActStructElms_rpvr(NE,L,NR,NZ)*RCCE(NE)
      ENDDO

      !maintenance deficit is paid by remobilization
      IF(Root1stStrutRemob(ielmc).GT.ZERO4Groth_pft(NZ))THEN
        Frac2Senes1=AZMAX1(AMIN1(1.0_r8,RCO2MDef1stStructPayoff_Oltd/Root1stStrutRemob(ielmc)))
      ELSE
        Frac2Senes1=1.0_r8
      ENDIF
    ELSE
      Root1stStrutRemob(1:NumPlantChemElms)=0._r8
      Frac2Senes1=0._r8
    ENDIF
    !
    !     PRIMARY ROOT LitrFall CAUSED BY REMOBILIZATION
    !
    !     CSNC,ZSNC,PSNC=literfall C,N,P
    !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
    !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
    !     Frac2Senes1=fraction of primary root C to be remobilized
    !     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass  
    ! respiration deficit -> some biomass into litterfall (offset negative growth), some biomass recycled (some offset negative CO2)
    if(Frac2Senes1.GT.0._r8)then
      DO NE=1,NumPlantChemElms          
        Root1stKill(NE) = Root1stActStructElms_rpvr(NE,L,NR,NZ)*Frac2Senes1  
        Root1stcylc(NE) = Root1stStrutRemob(NE)*Frac2Senes1
      ENDDO

      if(RootTipZone)then
        D6355: DO M=1,jsken
          DO NE=1,NumPlantChemElms          
            dsenecE      = AZMAX1(Root1stKill(NE)-Root1stcylc(NE))
            dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,icwood,M,NZ)*FracRootElmAllocm(NE,k_woody_comp)
            litrflxt(NE) = litrflxt(NE)+dlitrfall

            LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) +dlitrfall
            
            dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,iroot,M,NZ)*FracRootElmAllocm(NE,k_fine_comp)
            litrflxt(NE) = litrflxt(NE)+dlitrfall

            LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+dlitrfall
          ENDDO
        ENDDO D6355
      else
        !thickening zone
        !first mobilize active zone
        D6356: DO M=1,jsken
          DO NE=1,NumPlantChemElms          
            dsenecE      = AZMAX1(Root1stKill(NE)-Root1stcylc(NE))
            dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,istalk,M,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp)
            litrflxt(NE) = litrflxt(NE)+dlitrfall

            LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)+dlitrfall

            dlitrfall    = dsenecE*PlantElmAllocMat4Litr(NE,istalk,M,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp)
            litrflxt(NE) = litrflxt(NE)+dlitrfall

            LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ) +dlitrfall
          ENDDO
        ENDDO D6356
      endif
      !remobilization reduces growth    
      DO NE=1,NumPlantChemElms
        Root1stNetGrowthElms(NE)=Root1stNetGrowthElms(NE)-Root1stKill(NE)
      ENDDO 
    endif
  endif
 
  !  write(433,*)'exitremob',Root1stNetGrowthElms(ielmc),Root1stKill(ielmc),RCO2MDef1stStructPayoff_Oltd,litrflxt(ielmc),&
  !    'net',Root1stKill(ielmc)-RCO2MDef1stStructPayoff_Oltd-litrflxt(ielmc)
  !
  !     CONSUMPTION OF NON-STRUCTURAL C,N,P BY PRIMARY ROOTS
  !
  !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
  !     Rmaint1st_CO2=root maintenance respiration
  !     RNonstCO2_Oltd=respiration from non-structural C limited by O2
  !     RootMycoNonst4GrowC_Oltd=total non-structural C used in growth and respn ltd by O2
  !     RCO2Nonst4Nassim_Oltd=respiration for N assimilation unltd,ltd by O2
  !     Frac2Senes1=fraction of primary root C to be remobilized
  !     Root2ndStrutRemob(:)=remobilization of C,N,P from senescing root
  !     RootMycoNonst4Grow_Oltd(ielmn),RootMycoNonst4Grow_Oltd(ielmp)=nonstructural N,P ltd by O2 used in growth
  ! add CO2 respiraiton from remobilization  
  RCO2T1st_Oltd      = RCO2T1st_Oltd+RCO2MDef1stStructPayoff_Oltd
  RCO2T1st_OUltd      = RCO2T1st_OUltd+RCO2MDef1stStructPayoff_OUltd
  Root1stcylc(ielmc) = Root1stcylc(ielmc)-RCO2MDef1stStructPayoff_Oltd
  call PrintInfo('end '//subname)
  end associate
  end subroutine RemobilizePrimeActRoots

!----------------------------------------------------------------------------------------------------
  subroutine PrimeRootsExtension(yearIJ,L,Lnext,N,NR,NZ,RootNetGrowthElms,Root1stExtenzPP)
  !
  !Description
  !Grow primary roots at the tip in layer L or downward to layer Lnext
  !Currently, all primary root axes are assumed to only grow downward.
  !In the real world, they may grow with certain angles.
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: L,Lnext  !current and next soil layer to support root growth
  integer, intent(in) :: N    !root/myco indicator
  integer, intent(in) :: NR   !root axis 
  integer, intent(in) :: NZ   !pft number
  real(r8), intent(in) :: RootNetGrowthElms(NumPlantChemElms)  !whole population biomass growth for primary roots
  real(r8), intent(in) :: Root1stExtenzPP                      !per plant primary root axial extension
  real(r8) :: FGROL,FGROZ                                      !fraction of root extension in current and next lower soil layer
  real(r8) :: FGROLA
  integer :: NE
  REAL(R8) :: XFRE(NumPlantChemElms)
  character(len=*), parameter :: subname='PrimeRootsExtension'

! begin_execution
  associate(                                                                      &
    ZERO4Groth_pft                 => plt_biom%ZERO4Groth_pft                    ,& !input  :threshold zero for plang growth calculation, [-]
    CumSoilThickness_vr            => plt_site%CumSoilThickness_vr               ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    DLYR3                          => plt_site%DLYR3                             ,& !input  :vertical thickness of soil layer, [m]
    MaxNumRootLays                 => plt_site%MaxNumRootLays                    ,& !input  :maximum root layer number,[-]
    NGTopRootLayer_pft             => plt_morph%NGTopRootLayer_pft               ,& !input  :soil layer at planting depth, [-]
    SeedDepth_pft                  => plt_morph%SeedDepth_pft                    ,& !input  :seeding depth, [m]
    RootVH2O_pvr                   => plt_morph%RootVH2O_pvr                     ,& !input  :root layer volume water, [m2 d-2]
    Root1stMaxRadius1_pft          => plt_morph%Root1stMaxRadius1_pft            ,& !input  :maximum radius of primary roots, [m]
    Root1stMaxRadius_pft           => plt_morph%Root1stMaxRadius_pft             ,& !input  :maximum radius of primary roots, [m]
    PlantPopuLive_pft              => plt_site%PlantPopuLive_pft                 ,& !input  :plant population, [d-2]    
    Num1stAxesPerStructRootX_pft   => plt_morph%Num1stAxesPerStructRootX_pft     ,& !input  :number of primary root axes per NR, [d-2]
    RootCRRadius0_rpvr             => plt_morph%RootCRRadius0_rpvr               ,& !inoput: initial radius of roots that may undergo secondary growth, [m]
    Root1stRadius_rpvr             => plt_morph%Root1stRadius_rpvr               ,& !inoput: root layer radius for each primary axes,  [m]
    RootAge_rpvr                   => plt_morph%RootAge_rpvr                     ,& !inoput :root age,[h]
    RootMyco1stStrutElms_rpvr      => plt_biom%RootMyco1stStrutElms_rpvr         ,& !inoput :root layer element primary axes, [g d-2]
    Root1stActStructElms_rpvr      => plt_biom%Root1stActStructElms_rpvr         ,& !inoput :root layer active zone element in primary axes, [g d-2]
    RootMyco1stElm_raxs            => plt_biom%RootMyco1stElm_raxs               ,& !inoput :root C biomass of primary axes, [g d-2]
    PopuRootMycoC_pvr              => plt_biom% PopuRootMycoC_pvr                ,& !inoput :root layer C, [gC d-2]
    RootMycoNonstElms_rpvr         => plt_biom%RootMycoNonstElms_rpvr            ,& !inoput :root layer nonstructural element, [g d-2]
    PSIRootTurg_vr                 => plt_ew%PSIRootTurg_vr                      ,& !inoput :root turgor water potential, [Mpa]
    PSIRootOSMO_vr                 => plt_ew%PSIRootOSMO_vr                      ,& !inoput :root osmotic water potential, [Mpa]
    PSIRoot_pvr                    => plt_ew%PSIRoot_pvr                         ,& !inoput :root total water potential, [Mpa]
    Root1stRadius_pvr              => plt_morph%Root1stRadius_pvr                ,& !inoput :root layer radius primary axes, [m]
    Root1stLenPP_rpvr              => plt_morph%Root1stLenPP_rpvr                ,& !inoput :mean primary root axis length in soil layer per plant, [m d-2]
    Root1stDepz_raxes              => plt_morph%Root1stDepz_raxes                ,& !inoput :root layer depth, [m]
    NRoot1stTipLay_raxes           => plt_morph%NRoot1stTipLay_raxes              & !output :soil layer number for deepest root axes, [-]
  )
  call PrintInfo('beg '//subname)
  !  PRIMARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
  !
  !  Root1stExtenzPP=primary root length extension, downward. 
  !  RootMycoNonst4Grow_Oltd=primary root C growth ltd by O2
  !  PP=PFT population
  !  FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
  !  RTDP1=primary root depth from soil surface
  !  SeedDepth_pft=seeding depth
  !  Frac2Senes1=fraction of primary root C to be remobilized
  !  Root1stLenPP_rpvr=primary root length
  !  RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
  !  DLYR=soil layer thickness
  !  only the non-woody part is contributing to the elongation
  ! elongation due to positive growth of non-woody part

  !
  !     ALLOCATE PRIMARY ROOT GROWTH TO CURRENT
  !     AND NEXT SOIL LAYER WHEN PRIMARY ROOTS EXTEND ACROSS LOWER
  !     BOUNDARY OF CURRENT LAYER
  !
  !
  !     FGROL,FGROZ=fraction of Root1stExtenzPP in current,next lower soil layer  
  ! Question, 03/29/2024, Jinyun Tang: needs double check the following calculation 
  ! if FGROL < 1.0, then the extension is all in current layer, meaning FGROZ=0.0
  !  
  IF(Root1stExtenzPP.GT.ZERO4Groth_pft(NZ) .AND. L.LT.MaxNumRootLays)THEN
    !root grow donward, and has not reached the maximum rooting layer
    FGROL=AZMAX1(AMIN1(1.0_r8,(CumSoilThickness_vr(L)-Root1stDepz_raxes(NR,NZ))/Root1stExtenzPP))        
    FGROZ=AZMAX1(1.0_r8-FGROL)
  ELSE
    FGROL = 1.0_r8
    FGROZ = 0._r8
  ENDIF

  FGROLA=AMIN1(1._r8,Root1stExtenzPP*FGROL/DLYR3(L))

  !
  !     UPDATE STATE VARIABLES FOR PRIMARY ROOT LENGTH, GROWTH
  !     AND AXIS NUMBER
  !

  !     Root1stLenPP_rpvr=primary root length per plant
  !
  Root1stDepz_raxes(NR,NZ) = Root1stDepz_raxes(NR,NZ)+Root1stExtenzPP
  if(Root1stExtenzPP.LT.0._r8) &
    Root1stDepz_raxes(NR,NZ)=MAX(SeedDepth_pft(NZ),Root1stDepz_raxes(NR,NZ))
  
  DO NE=1,NumPlantChemElms
    RootMyco1stElm_raxs(NE,NR,NZ)         = RootMyco1stElm_raxs(NE,NR,NZ)+RootNetGrowthElms(NE)
    Root1stActStructElms_rpvr(NE,L,NR,NZ) = Root1stActStructElms_rpvr(NE,L,NR,NZ)+RootNetGrowthElms(NE)*FGROL
    RootMyco1stStrutElms_rpvr(NE,L,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)+RootNetGrowthElms(NE)*FGROL
  ENDDO
 
  RootCRRadius0_rpvr(L,NR,NZ) = AMAX1(Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ))
  Root1stRadius_rpvr(L,NR,NZ) = RootCRRadius0_rpvr(L,NR,NZ)
  Root1stLenPP_rpvr(L,NR,NZ)=Root1stLenPP_rpvr(L,NR,NZ)+Root1stExtenzPP*FGROL
  !write(9960,*)yearIJ%I*1000+yearIJ%J/24.,L,NR,Root1stLenPP_rpvr(L,NR,NZ),Root1stExtenzPP*FGROL,DLYR3(L)
  !
  !     TRANSFER STRUCTURAL, NONSTRUCTURAL C,N,P INTO NEXT SOIL LAYER
  !     WHEN PRIMARY ROOT EXTENDS ACROSS LOWER BOUNDARY
  !     OF CURRENT SOIL LAYER
  !
  !     PSIRoot_pvr,PSIRootTurg_vr,PSIRootOSMO_vr=root total,turgor,osmotic water potential
  !
  !root grows into next soil layer
  IF(FGROZ.GT.0.0_r8)THEN
    RootAge_rpvr(Lnext,NR,NZ)=1._r8
    DO NE=1,NumPlantChemElms
      Root1stActStructElms_rpvr(NE,Lnext,NR,NZ) = Root1stActStructElms_rpvr(NE,Lnext,NR,NZ)+RootNetGrowthElms(NE)*FGROZ
      RootMyco1stStrutElms_rpvr(NE,Lnext,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,Lnext,NR,NZ)+RootNetGrowthElms(NE)*FGROZ
    ENDDO

    PopuRootMycoC_pvr(N,Lnext,NZ)   = PopuRootMycoC_pvr(N,Lnext,NZ)+Root1stActStructElms_rpvr(ielmc,Lnext,NR,NZ)
    Root1stLenPP_rpvr(Lnext,NR,NZ)  = Root1stLenPP_rpvr(Lnext,NR,NZ)+Root1stExtenzPP*FGROZ
    Root1stRadius_pvr(N,Lnext,NZ)   = Root1stRadius_pvr(N,L,NZ)
    RootCRRadius0_rpvr(Lnext,NR,NZ) = RootCRRadius0_rpvr(L,NR,NZ)

    !divide and extend into next layer
    DO NE=1,NumPlantChemElms
      XFRE(NE)                              = FGROZ*RootMycoNonstElms_rpvr(NE,N,L,NZ)
      RootMycoNonstElms_rpvr(NE,N,L,NZ)     = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ)-XFRE(NE))
      RootMycoNonstElms_rpvr(NE,N,Lnext,NZ) = AZMAX1(RootMycoNonstElms_rpvr(NE,N,Lnext,NZ)+XFRE(NE))
    ENDDO
    PSIRoot_pvr(N,Lnext,NZ)      = PSIRoot_pvr(N,L,NZ)
    PSIRootOSMO_vr(N,Lnext,NZ)   = PSIRootOSMO_vr(N,L,NZ)
    PSIRootTurg_vr(N,Lnext,NZ)   = PSIRootTurg_vr(N,L,NZ)
    NRoot1stTipLay_raxes(NR,NZ)  = MIN(MAX(NGTopRootLayer_pft(NZ),Lnext),MaxNumRootLays)
  ENDIF
  !Retreat geometry is determined by depth, not a negative increment in one layer.
  !The caller then reconciles pools using L, the old tip layer.
  if(Root1stExtenzPP.LT.0._r8)call UpdateRetreatedRootGeometry(L,NR,NZ)
  call PrintInfo('end '//subname)
  end associate
  end subroutine PrimeRootsExtension

!----------------------------------------------------------------------------------------------------
  subroutine PrimeRootsWithdraw(L,NR,NZ,N,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr)
  implicit none
  integer, intent(in) :: L,NR,NZ,N

  real(r8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  integer :: LL,NN,NE,NTG,M
  real(r8) :: XFRD,FracRootSinkL
  real(r8) :: XFRE
  logical :: RootDepzChk
  real(r8), parameter :: scalp=0.99999_r8
  character(len=*), parameter :: subname='PrimeRootsWithdraw'

! begin_execution
  associate(                                                                      &
    ZERO4Groth_pft                 => plt_biom%ZERO4Groth_pft                    ,& !input  :threshold zero for plang growth calculation, [-]
    CumSoilThickness_vr            => plt_site%CumSoilThickness_vr               ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    ZEROS2                         => plt_site%ZEROS2                            ,& !input  :threshold zero for numerical stability,[-]
    DLYR3                          => plt_site%DLYR3                             ,& !input  :vertical thickness of soil layer, [m]
    VLSoilPoreMicP_vr              => plt_soilchem%VLSoilPoreMicP_vr             ,& !input  :volume of soil layer, [m3 d-2]
    LitrfallElms_pvr               => plt_bgcr%LitrfallElms_pvr                  ,& !input  :plant LitrFall element, [g d-2 h-1]
    Root1stLigStructElms_rpvr      => plt_biom%Root1stLigStructElms_rpvr         ,& !inoput :root layer lignified zone element in primary axes, [g d-2]
    Root1stDepz_raxes              => plt_morph%Root1stDepz_raxes                ,& !input  :root layer depth, [m]
    NGTopRootLayer_pft             => plt_morph%NGTopRootLayer_pft               ,& !input  :soil layer at planting depth, [-]
    iPlantNfixType_pft             => plt_morph%iPlantNfixType_pft               ,& !input  :N2 fixation type,[-]
    Myco_pft                       => plt_morph%Myco_pft                         ,& !input  :mycorrhizal type (no or yes),[-]
    SeedDepth_pft                  => plt_morph%SeedDepth_pft                    ,& !input  :seeding depth, [m]
    PlantPopuLive_pft              => plt_site%PlantPopuLive_pft                 ,& !input  :plant population, [d-2]
    Num1stAxesPerStructRootXPOP_pft   => plt_morph%Num1stAxesPerStructRootXPOP_pft     ,& !output :population primary root axes number on one structrual axis, [d-2]        
    Root1stAxesTipDepz2Surf_pft    => plt_morph%Root1stAxesTipDepz2Surf_pft      ,& !input  :plant primary depth relative to column surface, [m]
    Num1stAxesPerStructRootX_pft   => plt_morph%Num1stAxesPerStructRootX_pft     ,& !input  :primary root axes number, [d-2]
    icwood                         => pltpar%icwood                              ,& !input  :group id of coarse woody litter
    k_woody_comp                   => pltpar%k_woody_comp                        ,& !input  :woody litter complex id
    RootAge_rpvr                   => plt_morph%RootAge_rpvr                     ,& !inoput :root age,[h]
    RootMyco1stStrutElms_rpvr      => plt_biom%RootMyco1stStrutElms_rpvr         ,& !inoput :root layer element primary axes, [g d-2]
    RootMyco1stElm_raxs           => plt_biom%RootMyco1stElm_raxs               ,& !inout :axis primary structural C/N/P
    RootMediumStructElms_rpvr     => plt_biom%RootMediumStructElms_rpvr         ,& !inout :medium-root C/N/P
    RootMediumLength_rpvr         => plt_morph%RootMediumLength_rpvr            ,& !inout :medium-root length
    RootMediumXNum_rpvr           => plt_morph%RootMediumXNum_rpvr              ,& !inout :medium-root axes
    Root1stActStructElms_rpvr      => plt_biom%Root1stActStructElms_rpvr         ,& !inoput :root layer active zone element in primary axes, [g d-2]
    RootMyco2ndStrutElms_rpvr      => plt_biom%RootMyco2ndStrutElms_rpvr         ,& !inoput :root layer element secondary axes, [g d-2]
    RootMycoNonstElms_rpvr         => plt_biom%RootMycoNonstElms_rpvr            ,& !inoput :root layer nonstructural element, [g d-2]
    PopuRootMycoC_pvr              => plt_biom% PopuRootMycoC_pvr                ,& !inoput :root layer C, [gC d-2]
    RootNodulStrutElms_rpvr        => plt_biom%RootNodulStrutElms_rpvr           ,& !inoput :root layer nodule element, [g d-2]
    PlantElmAllocMat4Litr          => plt_soilchem%PlantElmAllocMat4Litr         ,& !input  :litter kinetic fraction, [-]
    RootNodulNonstElms_rpvr        => plt_biom%RootNodulNonstElms_rpvr           ,& !inoput :root layer nonstructural element, [g d-2]
    RootGasLossDisturb_pft         => plt_bgcr%RootGasLossDisturb_pft            ,& !inoput :gaseous flux fron root disturbance, [g d-2 h-1]
    trcg_rootml_pvr                => plt_rbgc%trcg_rootml_pvr                   ,& !inoput :root gas content, [g d-2]
    trcs_rootml_pvr                => plt_rbgc%trcs_rootml_pvr                   ,& !inoput :root aqueous content, [g d-2]
    Root1stXNumL_pvr               => plt_morph%Root1stXNumL_pvr                 ,& !inoput :root layer number primary axes, [d-2]
    Root2ndXNumL_rpvr              => plt_morph%Root2ndXNumL_rpvr                ,& !inoput :root layer number axes, [d-2]
    Root1stLenPP_rpvr              => plt_morph%Root1stLenPP_rpvr                ,& !inoput :mean primary root axis length in soil layer, [m d-2]
    Root2ndLen_rpvr                => plt_morph%Root2ndLen_rpvr                  ,& !inoput :population secondary root length, [m d-2]
    Root2ndXNum_rpvr               => plt_morph%Root2ndXNum_rpvr                 ,& !inoput :root layer number secondary axes, [d-2]
    NRoot1stTipLay_raxes           => plt_morph%NRoot1stTipLay_raxes              & !output :maximum soil layer number for root axes, [-]
  )
  !     TRANSFER PRIMARY ROOT C,N,P TO NEXT SOIL LAYER ABOVE THE
  !     CURRENT SOIL LAYER WHEN NEGATIVE PRIMARY ROOT GROWTH FORCES
  !     WITHDRAWAL FROM THE CURRENT SOIL LAYER AND ALL SECONDARY ROOTS
  !     IN THE CURRENT SOIL LAYER HAVE BEEN LOST
  !
  !     NRoot1stTipLay_raxes=deepest root layer
  !     VLSoilPoreMicP_vr=soil layer volume excluding macropore, rocks
  !     Root1stAxesTipDepz2Surf_pft(NR,NZ)=primary root depth from soil surface
  !     CumSoilThickness_vr=depth from soil surface to layer bottom
  !     SeedDepth_pft=seeding depth
  !     FracRootSinkL=fraction of primary+secondary root sink strength in axis
  !     WTRTD=root C mass
  !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root

  call PrintInfo('beg '//subname)

  D5115: DO LL=L,NGTopRootLayer_pft(NZ)+1,-1
    !
    !test root presence in layer LL-1, or root axes above seeding depth
    RootDepzChk=Root1stDepz_raxes(NR,NZ).LE.CumSoilThickness_vr(LL-1)
    !  
    !EXIT if LL-1 is not active soil or has root, no withdraw
    !
    IF(VLSoilPoreMicP_vr(LL-1).LE.ZEROS2 .OR. .not. RootDepzChk)CYCLE
    !
    IF(RootSinkC_vr(N,LL).GT.ZERO4Groth_pft(NZ))THEN
      FracRootSinkL=(Root1stSink_pvr(LL,NR)+Root2ndSink_pvr(N,LL,NR))/RootSinkC_vr(N,LL)
    ELSE
      FracRootSinkL=1.0_r8
    ENDIF

    FracRootSinkL=MAX(0._r8,MIN(1._r8,FracRootSinkL))

    DO NE=1,NumPlantChemElms
      Root1stActStructElms_rpvr(NE,LL-1,NR,NZ) = Root1stActStructElms_rpvr(NE,LL-1,NR,NZ)+Root1stActStructElms_rpvr(NE,LL,NR,NZ)
      RootMyco1stStrutElms_rpvr(NE,LL-1,NR,NZ) = Root1stActStructElms_rpvr(NE,LL-1,NR,NZ) &
        +Root1stLigStructElms_rpvr(NE,LL-1,NR,NZ)
      Root1stActStructElms_rpvr(NE,LL,NR,NZ)   = 0._r8
      RootMyco1stStrutElms_rpvr(NE,LL,NR,NZ)   = 0._r8

      !add lignified coarse root to litter
      if(Root1stLigStructElms_rpvr(NE,LL,NR,NZ).GT.0._r8)then
        DO M=1,jsken
          LitrfallElms_pvr(NE,M,k_woody_comp,LL,NZ) = LitrfallElms_pvr(NE,M,k_woody_comp,LL,NZ)+PlantElmAllocMat4Litr(NE,icwood,M,NZ)*Root1stLigStructElms_rpvr(NE,LL,NR,NZ)
        ENDDO
        RootMyco1stElm_raxs(NE,NR,NZ)=RootMyco1stElm_raxs(NE,NR,NZ)-Root1stLigStructElms_rpvr(NE,LL,NR,NZ)
        Root1stLigStructElms_rpvr(NE,LL,NR,NZ)=0._r8
      endif
      RootMediumStructElms_rpvr(NE,LL-1,NR,NZ)=RootMediumStructElms_rpvr(NE,LL-1,NR,NZ) &
        +RootMediumStructElms_rpvr(NE,LL,NR,NZ)
      RootMediumStructElms_rpvr(NE,LL,NR,NZ)=0._r8
    ENDDO
    RootMediumLength_rpvr(LL-1,NR,NZ)=RootMediumLength_rpvr(LL-1,NR,NZ)+RootMediumLength_rpvr(LL,NR,NZ)
    RootMediumLength_rpvr(LL,NR,NZ)=0._r8
    RootMediumXNum_rpvr(LL-1,NR,NZ)=RootMediumXNum_rpvr(LL-1,NR,NZ)+RootMediumXNum_rpvr(LL,NR,NZ)
    RootMediumXNum_rpvr(LL,NR,NZ)=0._r8

    D5110: DO NN=1,Myco_pft(NZ)
      
      DO NE=1,NumPlantChemElms
        RootMyco2ndStrutElms_rpvr(NE,NN,LL-1,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,NN,LL-1,NR,NZ)+RootMyco2ndStrutElms_rpvr(NE,NN,LL,NR,NZ)

        RootMyco2ndStrutElms_rpvr(NE,NN,LL,NR,NZ) = 0._r8
        XFRE = FracRootSinkL*RootMycoNonstElms_rpvr(NE,NN,LL,NZ)
        call ExchFluxLimiter(RootMycoNonstElms_rpvr(NE,NN,LL,NZ),RootMycoNonstElms_rpvr(NE,NN,LL-1,NZ),XFRE)

        RootMycoNonstElms_rpvr(NE,NN,LL,NZ)   = RootMycoNonstElms_rpvr(NE,NN,LL,NZ)-XFRE
        RootMycoNonstElms_rpvr(NE,NN,LL-1,NZ) = RootMycoNonstElms_rpvr(NE,NN,LL-1,NZ)+XFRE
      ENDDO

      ! Move geometry with the surviving biomass for roots and mycorrhizae.
      Root2ndLen_rpvr(NN,LL-1,NR,NZ) = Root2ndLen_rpvr(NN,LL-1,NR,NZ)+Root2ndLen_rpvr(NN,LL,NR,NZ)
      Root2ndLen_rpvr(NN,LL,NR,NZ) = 0._r8
      Root2ndXNum_rpvr(NN,LL-1,NR,NZ) = Root2ndXNum_rpvr(NN,LL-1,NR,NZ)+Root2ndXNum_rpvr(NN,LL,NR,NZ)
      Root2ndXNum_rpvr(NN,LL,NR,NZ) = 0._r8
      Root2ndXNumL_rpvr(NN,LL-1,NZ) = SUM(Root2ndXNum_rpvr(NN,LL-1,1:plt_morph%NumStructuralRootAxes_pft(NZ),NZ))
      Root2ndXNumL_rpvr(NN,LL,NZ) = SUM(Root2ndXNum_rpvr(NN,LL,1:plt_morph%NumStructuralRootAxes_pft(NZ),NZ))

      XFRD = FracRootSinkL*PopuRootMycoC_pvr(NN,LL,NZ)
      XFRD = AMAX1(AMIN1(PopuRootMycoC_pvr(NN,LL,NZ),XFRD),-PopuRootMycoC_pvr(NN,LL-1,NZ))*scalp
      PopuRootMycoC_pvr(NN,LL,NZ)   = PopuRootMycoC_pvr(NN,LL,NZ)-XFRD
      PopuRootMycoC_pvr(NN,LL-1,NZ) = PopuRootMycoC_pvr(NN,LL-1,NZ)+XFRD
      !
      !     WITHDRAW GASES IN PRIMARY ROOTS
      !
      !     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2, O2, CH4, N2O, NH3, H2
      !     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2, O2, CH4, N2O, NH3, H2
      !     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2, O2, CH4, N2O, NH3, H2
      !     FracRootSinkL=fraction of primary root sink strength in axis
      !
      DO NTG=idg_beg,idg_NH3
        RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-FracRootSinkL &
          *(trcg_rootml_pvr(NTG,NN,LL,NZ)+trcs_rootml_pvr(NTG,NN,LL,NZ))
        trcg_rootml_pvr(NTG,NN,LL,NZ) = (1.0_r8-FracRootSinkL)*trcg_rootml_pvr(NTG,NN,LL,NZ)
        trcs_rootml_pvr(NTG,NN,LL,NZ) = (1.0_r8-FracRootSinkL)*trcs_rootml_pvr(NTG,NN,LL,NZ)
      ENDDO
    ENDDO D5110
    !
    !     RESET ROOT NUMBER AND PRIMARY ROOT LENGTH
    !
    !     Root2ndXNum_rpvr,Root2ndXNumL_rpvr=number of secondary root axes
    !     Root1stXNumL_pvr=number of primary root axes
    !     Root1stLenPP_rpvr=primary root length
    !     CumSoilThickness_vr=depth from soil surface to layer bottom
    !     SeedDepth_pft=seeding depth
    !
    Root1stXNumL_pvr(LL,NZ)     = MAX(0._r8,Root1stXNumL_pvr(LL,NZ)-Num1stAxesPerStructRootXPOP_pft(NZ))

    !
    !     WITHDRAW C,N,P FROM ROOT NODULES IN LEGUMES
    !
    !     iPlantNfixType_pft=N2 fixation: 1,2,3=rapid to slow root symbiosis
    !     FracRootSinkL=fraction of primary root sink strength in axis
    !     WTNDL,WTNDLN,WTNDLP=root bacterial C,N,P mass
    !     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in root bacteria
    !
    IF(is_root_N2fix(iPlantNfixType_pft(NZ)))THEN
      DO NE=1,NumPlantChemElms
        XFRE = FracRootSinkL*RootNodulStrutElms_rpvr(NE,LL,NZ)
        XFRE = AMAX1(AMIN1(RootNodulStrutElms_rpvr(NE,LL,NZ),XFRE),-RootNodulStrutElms_rpvr(NE,LL-1,NZ))*scalp
        RootNodulStrutElms_rpvr(NE,LL,NZ) = RootNodulStrutElms_rpvr(NE,LL,NZ)-XFRE
        RootNodulStrutElms_rpvr(NE,LL-1,NZ) = RootNodulStrutElms_rpvr(NE,LL-1,NZ)+XFRE

        XFRE = FracRootSinkL*RootNodulNonstElms_rpvr(NE,LL,NZ)
        XFRE = AMAX1(AMIN1(RootNodulNonstElms_rpvr(NE,LL,NZ),XFRE),-RootNodulNonstElms_rpvr(NE,LL-1,NZ))*scalp
        RootNodulNonstElms_rpvr(NE,LL,NZ)   = RootNodulNonstElms_rpvr(NE,LL,NZ)-XFRE
        RootNodulNonstElms_rpvr(NE,LL-1,NZ) = RootNodulNonstElms_rpvr(NE,LL-1,NZ)+XFRE
      ENDDO
    ENDIF
    RootAge_rpvr(LL-1,NR,NZ)    = 0._r8
  ENDDO D5115
  call PrintInfo('end '//subname)
  end associate
  end subroutine PrimeRootsWithdraw
!----------------------------------------------------------------------------------------------------

  function morphogen_signal(Abase, dist_scaled)result(A)
  !
  !Description:
  !compute morphogen signal at scaled distance dist_scaled
  implicit none
  real(r8), intent(in) :: Abase        !base morephogen signal, [0-1]
  real(r8), intent(in) :: dist_scaled  !scaled distance, [-]

  real(r8) :: a

  if(dist_scaled<6._r8)then
    A=Abase+(1._r8-Abase)*exp(-dist_scaled)
  else
  A=Abase
  endif
  end function morphogen_signal
!----------------------------------------------------------------------------------------------------
  subroutine CytoKininDynamics(I,J,NZ,Root2ndSink_pvr)
  !
  !Steady-state model of cytokinin concentration in coarse roots
  !
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)  
  real(r8) :: QH2Oroots_rpvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)   !water flux along primary root axes
  integer  :: L,N,NR,L1

  character(len=*), parameter :: subname='CytoKininDynamics'
  real(r8), parameter :: phi_min =0.4_r8  
  real(r8), parameter :: phi_max =0.7_r8  !asymptotic limit of the conductive area fraction, [m2/m2]
  real(r8), parameter :: R_cri = 0.05     !Critical radius where the woody radius is considered 95% mature, [m]
  real(r8) :: LumenFraction, SingleMRootLumenArea
  real(r8) :: MediumRootCount(JZ1), MediumMass, MediumMassLeft, MediumInMass, ExportMass
  real(r8) :: FineVolume, FineMass, FineMassLeft, FineExport, FlowVolume
  real(r8), parameter :: dtCyto = 1._r8 !hourly cytokinin update, [h]
  real(r8) :: Cytokinin1stConc_copy(JZ1)
  real(r8) :: d_effective_hourly(JZ1)
  real(r8) :: d_water_hourly
  integer  :: NU_loc
  real(r8) :: sap_area,sap_area_root,sapAreaLeaf
  real(r8) :: Xe(JZ1),lumenVolC,lumenVolM
  real(r8) :: Root2ndSinkT,QH2OMediumRoots,fMR,FineRootWeight
  real(r8) :: Root2ndCumSinkProf(JZ1)
  real(r8) :: fMR_numer(pltpar%jroots,JZ1),fMR_denorm(pltpar%jroots,JZ1)
  real(r8), parameter :: dRootExp=0.55_r8
  real(r8) :: tlumenArea(JZ1)
  associate(                                                                      &    
    RPlantRootH2OUptk_pvr          => plt_ew%RPlantRootH2OUptk_pvr               ,& !input  :root water uptake, [m3 d-2 h-1]
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft            ,& !input  :plant growth type (vascular, non-vascular),[-]
    NU                             => plt_site%NU                                ,& !input  :current soil surface layer number, [-]
    NK                             => plt_site%NK                                ,& !input  :total number of root soil layers
    Myco_pft                       => plt_morph%Myco_pft                         ,& !input  :mycorrhizal type (no or yes),[-]
    NumStructuralRootAxes_pft      => plt_morph%NumStructuralRootAxes_pft        ,& !input  :number of structural root axes,[-]
    TKS_vr                         => plt_ew%TKS_vr                              ,& !input  :mean annual soil temperature, [K]
    DLYR3                          => plt_site%DLYR3                             ,& !input  :vertical thickness of soil layer, [m]
    AREA3                          => plt_site%AREA3                             ,& !input  :grid area size, [m2]
    PlantPopuLive_pft              => plt_site%PlantPopuLive_pft                 ,& !input  :plant population, [d-2]
    iPlant2ndGrothPattern_pft      => plt_pheno%iPlant2ndGrothPattern_pft        ,& !input  :plant expression of secondary growth, [-]
    Radius95pctMature_pft          => plt_morph%Radius95pctMature_pft            ,& !input  :Critical radius where the woody radius is considered 95% mature, [m]
    Num1stAxesPerStructRootXPOP_pft   => plt_morph%Num1stAxesPerStructRootXPOP_pft     ,& !input  :population primary root axes number, [d-2]
    SapFlowVlinear_pvr             => plt_ew%SapFlowVlinear_pvr                  ,& !output :Sap flow mean linear velocity, [m h-1]
    CumSoilThickness_vr            => plt_site%CumSoilThickness_vr               ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    Root1stRadius_rpvr             => plt_morph%Root1stRadius_rpvr               ,& !input  :root layer radius for each primary axes, [m]
    NGTopRootLayer_pft             => plt_morph%NGTopRootLayer_pft               ,& !input  :soil layer at planting depth, [-]
    MaxSoilLays4Root_pft           => plt_morph%MaxSoilLays4Root_pft             ,& !input  :maximum soil layer number for all root axes,[-]
    Root2ndProdCytok_rpvr          => plt_rbgc%Root2ndProdCytok_rpvr             ,& !input  :cytokinin production rate due to fine root/myco elongation, [gC CK h-1]
    RootMRProdCytok_rpvr           => plt_rbgc%RootMRProdCytok_rpvr              ,& !input  :cytokinin production rate due to medium roots metabolism, [gC CK h-1]
    RootMediumVH2O_rpvr            => plt_morph%RootMediumVH2O_rpvr              ,& !output :population medium-root lumen volume for axis NR, [m3 H2O]
    Root2ndVH2O_rpvr               => plt_morph%Root2ndVH2O_rpvr                 ,& !input  :water-occupied 2nd root volume, [m3 H2O m-3]
    Root1stVH2O_rpvr               => plt_morph%Root1stVH2O_rpvr                 ,& !input  :water-occupied 1st root volume, [m3 H2O m-3]
    xylemPhi_min_pft               => plt_morph%xylemPhi_min_pft                 ,& !input  :the fraction found in the youngest xylem that as lumen for tree, [m2/m2]
    xylemPhi_max_pft               => plt_morph%xylemPhi_max_pft                 ,& !input  :asymptotic limit fraction of the xyxlem area as lumen for tree, [m2/m2]
    xylemPhi_mean_pft              => plt_morph%xylemPhi_mean_pft                ,& !input  :the mean are fraction found in the root xylem as lumen for non-tree roots, [m2/m2]
    RootSinkScalar_pvr             => plt_morph%RootSinkScalar_pvr               ,& !input  :root sink scalar to account for the missing of intermediate size roots, [0-1]
    RootMediumStructElms_rpvr      => plt_biom%RootMediumStructElms_rpvr         ,& !input :medium-root structural elements, [g d-2]
    RootMediumLength_rpvr          => plt_morph%RootMediumLength_rpvr            ,& !input  :total medium root length in layer for for axis NR, [m d-2]      
    RootMediumLength_pvr           => plt_morph%RootMediumLength_pvr             ,& !output :within layer mean length of individual medium roots, [m]
    SapFlowVLinear_rpvr            => plt_ew%SapFlowVLinear_rpvr                 ,& !output :linear sap flow for primary root axis normalized by lumen area, [m h-1]
    Cytokinin2ndConc_rpvr          => plt_rbgc%Cytokinin2ndConc_rpvr             ,& !output :cytokinin concentration in fine roots, [gC m-3 H2O]
    Cytokinin1stConc_rpvr          => plt_rbgc%Cytokinin1stConc_rpvr             ,& !output :cytokinin concentration in corase roots, [gC m-3 H2O]
    CytokininMRConc_rpvr           => plt_rbgc%CytokininMRConc_rpvr              ,& !output :cytokinin concentration in medium size roots, [gC m-3 H2O]
    CRootLumenArea_rpvr            => plt_morph%CRootLumenArea_rpvr              ,& !output :coarse roots lumen area for root axes, [m2]    
    CRootLumenArea_pvr             => plt_morph%CRootLumenArea_pvr               ,& !output :roots lumen area for coarse root, [m2]    
    MRootLumenArea_rpvr            => plt_morph%MRootLumenArea_rpvr              ,& !output :population medium-root lumen area for axis NR, [m2]
    RootMediumXNum_rpvr            => plt_morph%RootMediumXNum_rpvr              ,& !inoput :Number of medium size root axes in layer for structrual axes, [d-2]    
    RootFineFrac2Med_rpvr           => plt_morph%RootFineFrac2Med_rpvr           ,& !fine-axis-count-weighted medium-root routing by category, [-]
    MRootLumenArea_pvr             => plt_morph%MRootLumenArea_pvr               ,& !output :Lumen area for medium size root, [m2]
    Root2ndXNum_rpvr               => plt_morph%Root2ndXNum_rpvr                 ,& !input  :root layer number secondary axes, [d-2]
    RootMediumRadius_rpvr          => plt_morph%RootMediumRadius_rpvr            ,& !input  :root layer radius for medium size axes, [m d-2]          
    Root2ndXNumL_rpvr              => plt_morph%Root2ndXNumL_rpvr                 & !input  :within soil layer number of whole population 2nd root axes, [d-2]
  )
  RootFineFrac2Med_rpvr(:,:,NZ) = 0._r8
  
  call PrintInfo('beg '//subname)
  
  !fine roots and mycorrhizae are assumed to functioning in parallel
  MRootLumenArea_rpvr(:,:,NZ) = 0._r8
  RootMediumVH2O_rpvr(:,:,NZ) = 0._r8
  MRootLumenArea_pvr(:,NZ) = 0._r8
  RootMediumLength_pvr(:,NZ) = 0._r8
  QH2Oroots_rpvr = 0._r8
  MediumRootCount = 0._r8
  if(NumStructuralRootAxes_pft(NZ).EQ.0)return
  !1. compute the water flux and linear sapflow for each NR, assuming no compressibility
  !2. compute fine root cytokinin net production 
  !3. compute the xylem cross section area

  NU_loc                   = NU
  SapFlowVlinear_pvr(:,NZ) = 0._r8
  tlumenArea(:)            = 0._r8
  CRootLumenArea_pvr(:,NZ)  = 0._r8
  if(is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))then          
    !sum up the total sink strength
    Root2ndSinkT = sum(Root2ndSink_pvr(ipltroot,NU:MaxSoilLays4Root_pft(NZ),:))
    !get the sap area demanded by canopy
    sapAreaLeaf  = CalcPotentialSapArea(I,J,NZ)
  endif
  DNR100: DO NR=1,NumStructuralRootAxes_pft(NZ)

    Root1stVH2O_rpvr(MIN(MaxSoilLays4Root_pft(NZ)+1,NK),NR,NZ)    = 0._r8
    SapFlowVLinear_rpvr(MIN(MaxSoilLays4Root_pft(NZ)+1,NK),NR,NZ) = 0._r8
    Root2ndCumSinkProf(:) = 0._R8
    DO L=MaxSoilLays4Root_pft(NZ),NU_loc,-1
      RootSinkScalar_pvr(L,NR,NZ) = 1._r8
      L1                    = MIN(L+1,NK)
      Root2ndCumSinkProf(L) = Root2ndCumSinkProf(L1)+Root2ndSink_pvr(ipltroot,L,NR)
      !decay and make a copy
      Cytokinin1stConc_rpvr(L,NR,NZ)=Cytokinin1stConc_rpvr(L,NR,NZ)*exp(-kDCytoC)
      
      if(isclose(Root1stRadius_rpvr(L,NR,NZ),0._r8))then
        SapFlowVLinear_rpvr(L,NR,NZ)=SapFlowVLinear_rpvr(L1,NR,NZ)
        CRootLumenArea_rpvr(L,NR,NZ) = CRootLumenArea_rpvr(L1,NR,NZ)
        cycle
      endif
      
      if(is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))then        
        LumenFraction = xylemPhi_min_pft(NZ)+(xylemPhi_max_pft(NZ)-xylemPhi_min_pft(NZ))*(1._r8-sfexp(-3._r8*Root1stRadius_rpvr(L,NR,NZ)/Radius95pctMature_pft(NZ)))
        sap_area_root = GetCoarseRootXylemSecArea(dmax,Root1stRadius_rpvr(L,NR,NZ))*Num1stAxesPerStructRootXPOP_pft(NZ)
        ! Retain geometric sapwood area when no fine-root sink remains.
        ! Normalize canopy-supported area only when total sink strength is positive.
        sap_area = sap_area_root
        if(Root1stRadius_rpvr(L,NR,NZ).GE.dmax .and. Root2ndSinkT.GT.0._r8)then
          sap_area = AMIN1(sapAreaLeaf*Root2ndCumSinkProf(L)/Root2ndSinkT, sap_area_root)
        endif
        !Both candidate sapwood areas represent the population; retain a per-plant lumen-area floor.
        CRootLumenArea_rpvr(L,NR,NZ) = AMAX1(sap_area*LumenFraction,1.E-9_R8*PlantPopuLive_pft(NZ)) &
          /RootSinkScalar_pvr(L,NR,NZ)
        ! Match the living-tissue filter used for the gas-transport geometry.
        if(RootMediumStructElms_rpvr(ielmc,L,NR,NZ).GT.0._r8 .and. &
           RootMediumXNum_rpvr(L,NR,NZ).GT.0._r8 .and. &
           RootMediumLength_rpvr(L,NR,NZ).GT.0._r8)then
          !medium size roots
          LumenFraction = xylemPhi_min_pft(NZ)+(xylemPhi_max_pft(NZ)-xylemPhi_min_pft(NZ))*(1._r8-sfexp(-3._r8*RootMediumRadius_rpvr(L,NR,NZ)/Radius95pctMature_pft(NZ)))
          sap_area_root = GetCoarseRootXylemSecArea(dmax,RootMediumRadius_rpvr(L,NR,NZ))
          SingleMRootLumenArea = sap_area_root*LumenFraction
          !Count and total length already include the population; apply each only once.
          MRootLumenArea_rpvr(L,NR,NZ) = SingleMRootLumenArea*RootMediumXNum_rpvr(L,NR,NZ)
          RootMediumVH2O_rpvr(L,NR,NZ) = SingleMRootLumenArea*RootMediumLength_rpvr(L,NR,NZ)
          MRootLumenArea_pvr(L,NZ)     = MRootLumenArea_pvr(L,NZ)+MRootLumenArea_rpvr(L,NR,NZ)
          RootMediumLength_pvr(L,NZ)   = RootMediumLength_pvr(L,NZ)+RootMediumLength_rpvr(L,NR,NZ)
          MediumRootCount(L)          = MediumRootCount(L)+RootMediumXNum_rpvr(L,NR,NZ)
        else
          MRootLumenArea_rpvr(L,NR,NZ) = 0._r8
          RootMediumVH2O_rpvr(L,NR,NZ) = 0._r8
        endif
      else        
        CRootLumenArea_rpvr(L,NR,NZ) = PICON*Root1stRadius_rpvr(L,NR,NZ)**2*Num1stAxesPerStructRootXPOP_pft(NZ)*xylemPhi_mean_pft(NZ)
        MRootLumenArea_rpvr(L,NR,NZ) = 0._r8
      endif     

      CRootLumenArea_pvr(L,NZ)  = CRootLumenArea_pvr(L,NZ)+CRootLumenArea_rpvr(L,NR,NZ)
      Root1stVH2O_rpvr(L,NR,NZ) = Root1stVH2O_rpvr(L1,NR,NZ)

      DO N=1,Myco_pft(NZ)  
        !roots actively taking up soil water 
        if(Root2ndXNumL_rpvr(N,L,NZ).GT.0._r8 .and. RPlantRootH2OUptk_pvr(N,L,NZ) .LT.0._r8)then
          QH2Oroots_rpvr(N,L,NR)=-RPlantRootH2OUptk_pvr(N,L,NZ)*Root2ndXNum_rpvr(N,L,NR,NZ)/Root2ndXNumL_rpvr(N,L,NZ)
        else
          QH2Oroots_rpvr(N,L,NR)=0._r8
        endif
        Root1stVH2O_rpvr(L,NR,NZ) = Root1stVH2O_rpvr(L,NR,NZ)+QH2Oroots_rpvr(N,L,NR)
  
        !production
        if(Root2ndVH2O_rpvr(N,L,NR,NZ).GT.0._r8)then
          Cytokinin2ndConc_rpvr(N,L,NR,NZ)=Cytokinin2ndConc_rpvr(N,L,NR,NZ)+Root2ndProdCytok_rpvr(N,L,NR,NZ)/Root2ndVH2O_rpvr(N,L,NR,NZ)
          !linear decay 
          Cytokinin2ndConc_rpvr(N,L,NR,NZ)=AZMAX1(Cytokinin2ndConc_rpvr(N,L,NR,NZ)*sfexp(-kDCytof(N)))
        else
          Cytokinin2ndConc_rpvr(N,L,NR,NZ)=0._r8
        endif
      ENDDO
      SapFlowVLinear_rpvr(L,NR,NZ) = Root1stVH2O_rpvr(L,NR,NZ)/CRootLumenArea_rpvr(L,NR,NZ)
      tlumenArea(L)                = tlumenArea(L)+CRootLumenArea_rpvr(L,NR,NZ)
      SapFlowVlinear_pvr(L,NZ)     = SapFlowVlinear_pvr(L,NZ)+Root1stVH2O_rpvr(L,NR,NZ)
    ENDDO  
  ENDDO DNR100

  DO L=NU_loc,MaxSoilLays4Root_pft(NZ)
    if(MediumRootCount(L).GT.0._r8)then
      RootMediumLength_pvr(L,NZ)=RootMediumLength_pvr(L,NZ)/MediumRootCount(L)
    endif
  ENDDO
  
  DO L=MaxSoilLays4Root_pft(NZ),NU_loc,-1
    IF(.not.isclose(tlumenArea(L),0._r8))then
      SapFlowVlinear_pvr(L,NZ)=AZERO(SapFlowVlinear_pvr(L,NZ)/tlumenArea(L))
    else
      L1=MIN(L+1,NK)
      SapFlowVlinear_pvr(L,NZ)=AZERO(SapFlowVlinear_pvr(L1,NZ))  
    endif    
  ENDDO
  
  !3. export, because of the no compressibility assumption, the initial cytokinin concentration is due to 
  !advection from lower layers
  fMR_denorm=0._r8;fMR_numer=0._r8
  DNR200: DO NR=1,NumStructuralRootAxes_pft(NZ)    
    !make a copy
    DO L=MaxSoilLays4Root_pft(NZ),NU_loc,-1      
      L1=MIN(L+1,NK)            
      lumenVolC=CRootLumenArea_rpvr(L,NR,NZ)*DLYR3(L)
      lumenVolM=RootMediumVH2O_rpvr(L,NR,NZ)
      fMR=0._r8
      if(lcoarseroot .and. lumenVolM.GT.0._r8)then
        !The local routing fraction is geometric, including when water flow is zero.
        !Approximate primary length by layer thickness times the population axis count,
        !matching the population lengths used for fine-root branching.
        fMR = RootMediumLength_rpvr(L,NR,NZ) &
          /(RootMediumLength_rpvr(L,NR,NZ)+DLYR3(L)*Num1stAxesPerStructRootXPOP_pft(NZ))
      endif
      ! Match the fine-axis weighting used for water uptake and layer conductance.
      ! Groups without medium roots still contribute their fine axes to the denominator.
      DO N=1,Myco_pft(NZ)
        FineRootWeight = MAX(0._r8,Root2ndXNum_rpvr(N,L,NR,NZ))
        fMR_numer(N,L) = fMR_numer(N,L)+FineRootWeight*fMR
        fMR_denorm(N,L) = fMR_denorm(N,L)+FineRootWeight
      ENDDO
      QH2OMediumRoots=0._r8
      MediumInMass=0._r8
      if(lumenVolC.GT.0._r8)then
        DO N=1,Myco_pft(NZ)
          QH2OMediumRoots = QH2OMediumRoots+QH2Oroots_rpvr(N,L,NR)*fMR
          FineVolume = Root2ndVH2O_rpvr(N,L,NR,NZ)
          FlowVolume = AZMAX1(QH2Oroots_rpvr(N,L,NR))*dtCyto
          FineExport = 0._r8
          if(FineVolume.GT.0._r8 .and. FlowVolume.GT.0._r8)then
            !Production and decay were applied in DNR100. Debit each donor once.
            FineMass = AZMAX1(Cytokinin2ndConc_rpvr(N,L,NR,NZ))*FineVolume
            FineMassLeft = FineMass*(FineVolume/(FineVolume+FlowVolume))
            FineExport = FineMass-FineMassLeft
            Cytokinin2ndConc_rpvr(N,L,NR,NZ) = FineMassLeft/FineVolume
          endif
          !Split the same exported mass; without medium roots fMR is zero.
          MediumInMass = MediumInMass+fMR*FineExport
          Cytokinin1stConc_rpvr(L,NR,NZ) = Cytokinin1stConc_rpvr(L,NR,NZ) &
            +(1._r8-fMR)*FineExport/lumenVolC
        ENDDO
      endif
      if(lumenVolM.GT.0._r8)then
        !Population mass after sources and decay; local metabolism continues at zero flow.
        MediumMass = AZMAX1(CytokininMRConc_rpvr(L,NR,NZ)*lumenVolM+MediumInMass &
          +RootMRProdCytok_rpvr(L,NR,NZ)*dtCyto)*sfexp(-kDCytof(1)*dtCyto)

        if(QH2OMediumRoots.GT.0._r8 .and. lumenVolC.GT.0._r8)then
          !Implicit well-mixed outflow: export cannot exceed the donor inventory.
          MediumMassLeft = AZMAX1(MediumMass*(lumenVolM/(lumenVolM+QH2OMediumRoots*dtCyto)))
          ExportMass = AZMAX1(MediumMass-MediumMassLeft)
          MediumMass = MediumMassLeft
          Cytokinin1stConc_rpvr(L,NR,NZ) = Cytokinin1stConc_rpvr(L,NR,NZ)+ExportMass/lumenVolC
        endif
        CytokininMRConc_rpvr(L,NR,NZ)=MediumMass/lumenVolM
      else
        CytokininMRConc_rpvr(L,NR,NZ)=0._r8
      endif
      Cytokinin1stConc_rpvr(L,NR,NZ) = AZMAX1(Cytokinin1stConc_rpvr(L,NR,NZ))
      CRootLumenArea_rpvr(L,NR,NZ) = CRootLumenArea_rpvr(L,NR,NZ) +1.e-12_r8      
    ENDDO

    Cytokinin1stConc_copy(:) = Cytokinin1stConc_rpvr(:,NR,NZ)    
    !counter number of grids    
    L1=MaxSoilLays4Root_pft(NZ)-NU_loc+1
    !reverse the grid ordinates
    DO L=MaxSoilLays4Root_pft(NZ),NU_loc,-1
      Xe(L)=CumSoilThickness_vr(MaxSoilLays4Root_pft(NZ))-CumSoilThickness_vr(L-1)
    ENDDO

    call advect_remap_mass_loss(L1, dt=1._r8, xr=Xe(MaxSoilLays4Root_pft(NZ):NU_loc:-1), &
      c     = Cytokinin1stConc_copy(MaxSoilLays4Root_pft(NZ):NU_loc:-1),     &
      Areas = CRootLumenArea_rpvr(MaxSoilLays4Root_pft(NZ):NU_loc:-1,NR,NZ), &
      ur    = SapFlowVLinear_rpvr(MaxSoilLays4Root_pft(NZ):NU_loc:-1,NR,NZ), &
      c_new = Cytokinin1stConc_rpvr(MaxSoilLays4Root_pft(NZ):NU_loc:-1,NR,NZ))

    !solver area scaled diffusion of cytokinin in 1st roots
    Cytokinin1stConc_copy(:) = Cytokinin1stConc_rpvr(:,NR,NZ)    

    !compute diffusivity
    DO L=NU_loc,MaxSoilLays4Root_pft(NZ)
      d_water_hourly               = (4.3e-5_r8 * (1.03_r8**(TKS_vr(L) - 298.15_r8))) / 24.0_r8
      d_effective_hourly(L)        = d_water_hourly / 1.2_r8  ! 1.2 is the fixed xylem tortuosity      
    ENDDO

    call solve_root_diffusion_step(L1, dt=1._r8, c_init=Cytokinin1stConc_copy(NU_loc:MaxSoilLays4Root_pft(NZ)), &
      lumen_areas=CRootLumenArea_rpvr(NU_loc:MaxSoilLays4Root_pft(NZ),NR,NZ), &
      layer_thicknesses=DLYR3(NU_loc:MaxSoilLays4Root_pft(NZ)), &
      d_effective=d_effective_hourly(NU_loc:MaxSoilLays4Root_pft(NZ)), &
      c_next=Cytokinin1stConc_rpvr(NU_loc:MaxSoilLays4Root_pft(NZ),NR,NZ))
  ENDDO DNR200
  
  DO L=NU,MaxSoilLays4Root_pft(NZ)
    DO N=1,Myco_pft(NZ)
      if(fMR_denorm(N,L).GT.0._r8)then
        RootFineFrac2Med_rpvr(N,L,NZ)=fMR_numer(N,L)/fMR_denorm(N,L)
      endif
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine CytoKininDynamics

!----------------------------------------------------------------------------------------------------
  function GetCoarseRootXylemSecArea(dmax,R)result(area)
  !
  !compute coarse root xylem cross section area 
  implicit none
  real(r8), intent(in) :: dmax  !Sapwood Depth, or the radius when heartwood start to form, species dependent [m]
  real(r8), intent(in) :: R     !coarse root radius [m]

  real(r8) :: area
  real(r8) :: d,hb
  real(r8), parameter :: Rcri=1.e-3_r8       !critical radius for very young roots [m]

  if(R<Rcri)then     !for soybean and very young seedling
    hb   = 0.3_r8*R
    area = PiCON*(R-hb)**2
  elseif(R<dmax)then
    hb   = 0.1_r8*R
    area = PiCON*(R-hb)**2
  else    
    d    = dmax*tanh(R/dmax)
    area = PICON*(2._r8*R*d-d**2)
  endif
  end function GetCoarseRootXylemSecArea
!----------------------------------------------------------------------------------------------------
  function CalcPotentialSapArea(I,J,NZ)result(sapAreaLeaf)
  !
  !compute leaf supported total xylem area
  !using the basal stalk radius
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  real(r8), parameter :: huber_cef=3._r8    !humber coefficient
  real(r8) :: sapAreaGeom
  real(r8) :: sapAreaLeaf
  real(r8) :: DTransportTube
  real(r8) :: sapAreaLeafMax
  character(len=*), parameter :: subname='CalcPotentialSapArea'

  associate(                                                              &
    CanopyLeafAreaMAX_pft      => plt_morph%CanopyLeafAreaMAX_pft        ,& !input  :running maximum leaf area, [m2 d-2]
    lreset_laimax_pft          => plt_morph%lreset_laimax_pft            ,& !input  :toggle to reset max leaf area, [-]
    PlantPopuLive_pft          => plt_site%PlantPopuLive_pft             ,& !input  :live population in grid cell, [d-2]
    StalkAveRadius_pft         => plt_morph%StalkAveRadius_pft           ,& !input  :main stalk radius,[m]    
    CanopyLeafArea_pft         => plt_morph%CanopyLeafArea_pft            & !input  :plant canopy leaf area, [m2 d-2]
  )
  call PrintInfo('beg '//subname)
  DTransportTube = AMIN1(ZSTX,FSTK*StalkAveRadius_pft(NZ))
  !Match the population basis of CanopyLeafAreaMAX_pft before applying the nonlinear area limit.
  sapAreaGeom    = PICON*(2._r8*StalkAveRadius_pft(NZ)*DTransportTube-DTransportTube**2)*PlantPopuLive_pft(NZ)

  !obtain the running maximum of leaf area
  if(lreset_laimax_pft(NZ))then
    CanopyLeafAreaMax_pft(NZ) = MAX(sapAreaGeom/huber_cef,1.e-6_r8)
  else
    CanopyLeafAreaMax_pft(NZ) = AMAX1(CanopyLeafAreaMax_pft(NZ),CanopyLeafArea_pft(NZ))
  endif
  !  
  sapAreaLeafMax = huber_cef*CanopyLeafAreaMax_pft(NZ)
  if(sapAreaLeafMax.GT.1.e-5_r8)then
    sapAreaLeaf    = sapAreaLeafMax*(1._r8-sfexp(-sapAreaGeom/sapAreaLeafMax))
  else
    sapAreaLeaf    = sapAreaLeafMax
  endif
  call PrintInfo('end '//subname)
  end associate

  end function CalcPotentialSapArea
!----------------------------------------------------------------------------------------------------

  subroutine GetMediumRootPholTranspResit(yearIJ,NZ,SingleMedRootPholeResist_vr)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), intent(out) :: SingleMedRootPholeResist_vr(JZ1,pltpar%MaxNumRootAxes) !single medium root resistance for phloem transport
  integer :: L,NR
  real(r8) :: AreaTranspt,DTransptTube

  character(len=*), parameter :: subname='GetMediumRootPholTranspResit'


  call PrintInfo('beg '//subname)
  associate(                                                             &
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft               ,& !input  :threshold zero for plang growth calculation, [-]  
    NRoot1stTipLay_raxes       => plt_morph%NRoot1stTipLay_raxes        ,& !input  :maximum soil layer number for root axes, [-]        
    Root1stMaxRadius1_pft      => plt_morph%Root1stMaxRadius1_pft       ,& !input  :root radius primary axes, [m]            
    NumStructuralRootAxes_pft  => plt_morph%NumStructuralRootAxes_pft   ,& !input  :number of structural root axes,[-]        
    iPlant2ndGrothPattern_pft  => plt_pheno%iPlant2ndGrothPattern_pft   ,& !input  :plant expression of secondary growth, [-]    
    iPlantRootProfile_pft      => plt_pheno%iPlantRootProfile_pft       ,& !input  :plant growth type (vascular, non-vascular),[-]    
    RootMediumMeanLength_rpvr  => plt_morph%RootMediumMeanLength_rpvr   ,& !input  :mean medium-root length per structural axis, [m]
    Root1stLenPP_rpvr          => plt_morph%Root1stLenPP_rpvr           ,& !input  :length per primary-root axis in layer, [m]
    RootMediumRadius_rpvr      => plt_morph%RootMediumRadius_rpvr        & !input  :root layer radius for medium size axes, [m d-2]      
  )
  call plt_morph%RefreshMediumRootMeanLength(NZ)
  SingleMedRootPholeResist_vr=0._r8
  if(.not.lcoarseroot .or. .not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)))return

  DO NR=1,NumStructuralRootAxes_pft(NZ)    
    DO L=1,NRoot1stTipLay_raxes(NR,NZ)       
      !Also evaluate a nascent axis for establishment; this creates no living count.
      if(Root1stLenPP_rpvr(L,NR,NZ).GT.0._r8)then
        DTransptTube                   = AMIN1(ZSTX,AMAX1(FSTK*RootMediumRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(ipltroot,NZ)))
        AreaTranspt                    = 2._r8*AMAX1(RootMediumRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(ipltroot,NZ))*DTransptTube-DTransptTube**2
        SingleMedRootPholeResist_vr(L,NR) = AMAX1(0.5_r8*RootMediumMeanLength_rpvr(L,NR,NZ),Root2ndTipLen4uptk)/AreaTranspt
      endif
    enddo
  enddo
  call PrintInfo('end '//subname)
  end associate
  end subroutine GetMediumRootPholTranspResit
!----------------------------------------------------------------------------------------------------

  subroutine GetPrimaryRootPholTranspResist(yearIJ,NZ,Root1stLocDepz_vr,RootEffDepz_vr,R1stPholemResist_vr)
  !
  !compute the 
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), intent(in) :: Root1stLocDepz_vr(pltpar%MaxNumRootAxes,JZ1)  
  real(r8), intent(in) :: RootEffDepz_vr(pltpar%MaxNumRootAxes,JZ1)  
  real(r8), intent(out) :: R1stPholemResist_vr(JZ1,pltpar%MaxNumRootAxes)       ![m-1] resistance scalar for phloem transport in structrual roots

  real(r8) :: AreaTranspt(JZ1)  
  
  real(r8) :: DTransptTube
  integer :: NR,L
  character(len=*), parameter :: subname='GetPrimaryRootPholTranspResist'

  associate(                                                             &
    Root1stDepz_raxes          => plt_morph%Root1stDepz_raxes           ,& !input  :root layer depth, [m]  
    NGTopRootLayer_pft         => plt_morph%NGTopRootLayer_pft          ,& !input  :soil layer at planting depth, [-]        
    HypocotHeight_pft          => plt_morph%HypocotHeight_pft           ,& !input  :cotyledon height, [m]    
    NU                         => plt_site%NU                           ,& !input  :current soil surface layer number, [-]    
    MaxNumRootLays             => plt_site%MaxNumRootLays               ,& !input  :maximum root layer number,[-]    
    Root1stMaxRadius1_pft      => plt_morph%Root1stMaxRadius1_pft       ,& !input  :root radius primary axes, [m]        
    CumSoilThickness_vr        => plt_site%CumSoilThickness_vr          ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]    
    NumStructuralRootAxes_pft  => plt_morph%NumStructuralRootAxes_pft   ,& !input  :number of structural root axes,[-]    
    SeedDepth_pft              => plt_morph%SeedDepth_pft               ,& !input  :seeding depth, [m]        
    NRoot1stTipLay_raxes       => plt_morph%NRoot1stTipLay_raxes        ,& !input  :maximum soil layer number for root axes, [-]        
    CanopyHeight4WatUptake_pft => plt_morph%CanopyHeight4WatUptake_pft  ,& !input  :canopy height, [m]    
    StalkAveRadius_pft         => plt_morph%StalkAveRadius_pft          ,& !input  :main stalk radius,[m]    
    fRootTube_rpvr             => plt_morph%fRootTube_rpvr              ,& !output :fraction of root cross section for transport,[-]        
    DLYR3                      => plt_site%DLYR3                        ,& !input  :vertical thickness of soil layer, [m]
    Root1stRadius_rpvr         => plt_morph%Root1stRadius_rpvr           & !input  :root layer radius for each primary axes, [m]
  )

  call PrintInfo('beg '//subname)

  DO NR=1,NumStructuralRootAxes_pft(NZ)    
    DO L=1,NRoot1stTipLay_raxes(NR,NZ)
      if(RootEffDepz_vr(NR,L).LE.SeedDepth_pft(NZ) .OR. L.LT.NGTopRootLayer_pft(NZ) .OR. L.EQ.1)then                
        Root1stRadius_rpvr(L,NR,NZ) = AZMAX1(StalkAveRadius_pft(NZ),Root1stMaxRadius1_pft(ipltroot,NZ))
        DTransptTube                = AMIN1(ZSTX,AMAX1(FSTK*Root1stRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(ipltroot,NZ)))
        AreaTranspt(L)              = 2._r8*AMAX1(Root1stRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(ipltroot,NZ))*DTransptTube-DTransptTube**2
        R1stPholemResist_vr(L,NR)   = (CanopyHeight4WatUptake_pft(NZ)+RootEffDepz_vr(NR,L))/AreaTranspt(L)        
      else
        !R^2-(R-d)^2=(2R-d)*d
        DTransptTube           = AMIN1(ZSTX,AMAX1(FSTK*Root1stRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(ipltroot,NZ)))
        AreaTranspt(L)         = 2._r8*AMAX1(Root1stRadius_rpvr(L,NR,NZ),Root1stMaxRadius1_pft(ipltroot,NZ))*DTransptTube-DTransptTube**2
        R1stPholemResist_vr(L,NR) = R1stPholemResist_vr(L-1,NR)+0.5_r8*(Root1stLocDepz_vr(NR,L-1)/AreaTranspt(L-1)+Root1stLocDepz_vr(NR,L)/AreaTranspt(L))
      endif  

      if(Root1stRadius_rpvr(L,NR,NZ).GT.0._r8)then
        fRootTube_rpvr(L,NR,NZ) = 1._r8-(1._r8-AMIN1(DTransptTube/Root1stRadius_rpvr(L,NR,NZ),FSTK))**2
      else
        fRootTube_rpvr(L,NR,NZ) = 1._r8/FRTX
      endif
    ENDDO
  ENDDO
  
  call PrintInfo('end '//subname)

  end associate
  end subroutine GetPrimaryRootPholTranspResist

!----------------------------------------------------------------------------------------------------

  subroutine GetInLayerRootDepth(NZ,Root1stLocDepz_vr,RootEffDepz_vr)
  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(out) :: Root1stLocDepz_vr(MaxNumRootAxes,JZ1)
  real(r8), intent(out) :: RootEffDepz_vr(MaxNumRootAxes,JZ1)

  integer :: NR,L

  associate(                                                             &
    Root1stDepz_raxes          => plt_morph%Root1stDepz_raxes           ,& !input  :root layer depth, [m]
    CumSoilThickness_vr        => plt_site%CumSoilThickness_vr          ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    DLYR3                      => plt_site%DLYR3                        ,& !input  :vertical thickness of soil layer, [m]
    MaxNumRootLays             => plt_site%MaxNumRootLays               ,& !input  :lowest layer included in root mass accounting
    HypocotHeight_pft          => plt_morph%HypocotHeight_pft           ,& !input  :cotyledon height, [m]
    MaxSoilLays4Root_pft       => plt_morph%MaxSoilLays4Root_pft        ,& !input  :maximum soil layer number for all root axes,[-]  
    NGTopRootLayer_pft         => plt_morph%NGTopRootLayer_pft               ,& !input  :layer at planting depth
    SeedDepth_pft              => plt_morph%SeedDepth_pft               ,& !input  :seeding depth, [m]
    NRoot1stTipLay_raxes       => plt_morph%NRoot1stTipLay_raxes        ,& !input  :maximum soil layer number for root axes, [-]
    NumStructuralRootAxes_pft  => plt_morph%NumStructuralRootAxes_pft    & !input  :number of structural root axes,[-]
  )

  Root1stLocDepz_vr = 0._r8
  RootEffDepz_vr    = 0._r8
  ! Root depth may extend beyond the bottom of the permitted rooting domain.
  ! Keep its tip in the last root layer; NK includes an extra soil layer whose
  ! reserves are excluded from SumRootBiome and must not supply the shoot.
  MaxSoilLays4Root_pft(NZ) = MIN(MaxSoilLays4Root_pft(NZ),MaxNumRootLays)
  DO NR  = 1, NumStructuralRootAxes_pft(NZ)
    DO L=1,MaxNumRootLays
      IF(L.GT.NGTopRootLayer_pft(NZ) .AND. Root1stDepz_raxes(NR,NZ).LE.CumSoilThickness_vr(L-1))CYCLE
      IF(Root1stDepz_raxes(NR,NZ).LE.CumSoilThickness_vr(L) .OR. L.EQ.MaxNumRootLays)then
        MaxSoilLays4Root_pft(NZ)    = MAX(MaxSoilLays4Root_pft(NZ),L)
        NRoot1stTipLay_raxes(NR,NZ) = L
      endif
      Root1stLocDepz_vr(NR,L) = AZMAX1(Root1stDepz_raxes(NR,NZ)-RTDPX-CumSoilThickness_vr(L-1)) !distance from root tipbase to layer L-1 bottom
      Root1stLocDepz_vr(NR,L) = AZMAX1(AMIN1(DLYR3(L),Root1stLocDepz_vr(NR,L))- AZMAX1(SeedDepth_pft(NZ)-HypocotHeight_pft(NZ)-CumSoilThickness_vr(L-1)))
      RootEffDepz_vr(NR,L)    = AMAX1(SeedDepth_pft(NZ),CumSoilThickness_vr(L-1))+0.5_r8*Root1stLocDepz_vr(NR,L)
    ENDDO
  ENDDO
  end associate
  end subroutine GetInLayerRootDepth

!----------------------------------------------------------------------------------------------------
  subroutine UpdateRetreatedRootGeometry(OldTipLayer,NR,NZ)
  !Rebuild primary geometry from actual depth in the column coordinate system.
  implicit none
  integer, intent(in) :: OldTipLayer,NR,NZ
  integer :: LL,TipLayer
  real(r8) :: Depth,SeedDepth

  associate( &
    Root1stDepz_raxes             => plt_morph%Root1stDepz_raxes             ,& !input :primary-root axis depth, [m]
    SeedDepth_pft                 => plt_morph%SeedDepth_pft                 ,& !input :planting depth, [m]
    NGTopRootLayer_pft            => plt_morph%NGTopRootLayer_pft            ,& !input :layer at planting depth
    CumSoilThickness_vr           => plt_site%CumSoilThickness_vr           ,& !input :depth to layer bottom, [m]
    DLYR3                         => plt_site%DLYR3                         ,& !input :layer thickness, [m]
    SoilBulkDensity_vr            => plt_soilchem%SoilBulkDensity_vr            ,& !input :soil bulk density, [Mg m-3]
    ZERO                          => plt_site%ZERO                          ,& !input :numerical zero threshold
    Root1stLenPP_rpvr             => plt_morph%Root1stLenPP_rpvr             ,& !inout :primary-root length per axis in layer, [m]
    Root1stLenLoc_rpvr            => plt_morph%Root1stLenLoc_rpvr            ,& !inout :local primary-root depth in layer, [m]
    NRoot1stTipLay_raxes          => plt_morph%NRoot1stTipLay_raxes          ,& !inout :primary-root tip layer
    Root1stAxesTipDepz2Surf_pft   => plt_morph%Root1stAxesTipDepz2Surf_pft    & !inout :primary-root tip depth relative to column surface, [m]
  )

  Depth=Root1stDepz_raxes(NR,NZ)
  SeedDepth=SeedDepth_pft(NZ)
  TipLayer=NGTopRootLayer_pft(NZ)
  DO LL=NGTopRootLayer_pft(NZ),OldTipLayer
    Root1stLenPP_rpvr(LL,NR,NZ)=MAX(0._r8,MIN(Depth,CumSoilThickness_vr(LL)) &
      -MAX(SeedDepth,CumSoilThickness_vr(LL-1)))
    Root1stLenLoc_rpvr(LL,NR,NZ)=MAX(0._r8, &
      MIN(DLYR3(LL),Depth-CumSoilThickness_vr(LL-1)))
    if(Depth.GT.CumSoilThickness_vr(LL-1))TipLayer=LL
  ENDDO
  NRoot1stTipLay_raxes(NR,NZ)=TipLayer
  Root1stAxesTipDepz2Surf_pft(NR,NZ)=Depth
  if(SoilBulkDensity_vr(TipLayer).GT.ZERO) &
    Root1stAxesTipDepz2Surf_pft(NR,NZ)=Depth-CumSoilThickness_vr(0)
  end associate
  end subroutine UpdateRetreatedRootGeometry
!----------------------------------------------------------------------------------------------------
  pure function func_cyto(cyto_conc,enh_cyto_pft)result(fctyok)
  implicit none
  real(r8), intent(in) :: cyto_conc
  real(r8), intent(in) :: enh_cyto_pft
  real(r8) :: fctyok
  real(r8), parameter :: Km_cyto = 1.2_r8     !~20nM, after upscaling the cytokinin production rate by 1.e3
  real(r8), parameter :: Ki_cyto = 12._r8     !~200 nM, after upscaling the cytokinin production rate by 1.e3

  fctyok     = enh_cyto_pft*cyto_conc/(Km_cyto+cyto_conc*(1._r8+cyto_conc/Ki_cyto))
  fctyok     = AMIN1(fctyok,1._r8)     
  end function func_cyto
!----------------------------------------------------------------------------------------------------
  subroutine SummarizeRootSink(yearIJ,NZ,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,Root1stSinkTip,&
    RootMSink_pvr,RootSinkC,fUSE4MR_pvr)
  !
  !DESCRIPTION
  !summarize root sink strength for growth allocation and shoot-root coupling
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8),INTENT(OUT) :: RootSinkC_vr(pltpar%jroots,JZ1)       !conductance for nonstrucal CNP allocation
  real(r8),intent(out) :: Root1stSink_pvr(JZ1,pltpar%MaxNumRootAxes)
  real(r8),intent(out) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8),intent(out) :: Root1stSinkTip(pltpar%MaxNumRootAxes)
  real(r8),intent(out) :: RootMSink_pvr(JZ1,pltpar%MaxNumRootAxes)  
  real(r8),INTENT(OUT) :: RootSinkC(pltpar%jroots)
  real(r8),intent(out) :: fUSE4MR_pvr(JZ1,pltpar%MaxNumRootAxes) !fine-root activation of medium-root branching, [0-1]
  real(r8) :: R1stPholemResist_vr(JZ1,MaxNumBranches)    
  integer  :: N,L,K,NR,NE,ntu,LTip
  REAL(R8) :: Root1stLocDepz_vr(MaxNumRootAxes,JZ1)  !local rooting depth in layer L
  real(r8) :: RecoRootMycoC4Nup,CUPRO,CUPRC
  real(r8) :: DistCanopyPrimeRootTip       !distance between primary root tip and canopy top, which characterizes source-sink distance
  real(r8) :: RootEffDepz,PrimaryRootCd4NonstTP
  real(r8) :: FineRootCd4NonstTP,rscal,Root1stTipRadius
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: RCO2flx,dUptakeN,dUptakeP,DTransptTube,AreaTranspt,dist2Tip
  real(r8) :: DistRootEffDepz,dist_scaled,BaseSink
  real(r8) :: RootEffDepz_vr(MaxNumRootAxes,JZ1)    
  real(r8) :: SingleMedRootPholeResist_vr(JZ1,pltpar%MaxNumRootAxes)
  REAL(R8) :: CRootSinkProf(JZ1,MaxNumRootAxes)  !fraction profile of fine roots in soil
  real(r8) :: dSapFlow,cytoc_scal
  real(r8), parameter :: Km_cyto  = 1.2_r8     !~20nM, after upscaling the cytokinin production rate by 1.e3 
  real(r8), parameter :: Ki_cyto = 12._r8      !~200 nM, after upscaling the cytokinin production rate by 1.e3  
  real(r8) :: Km_app, Ki_app,CMFR,CFR,UsedConnCd,ConnCd  
  real(r8) :: fctyok,MediumRootCd4NonstTP
  real(r8) :: PotentialMediumAxes,MediumAxesForGrowth
  real(r8), parameter :: MinMediumEstablishmentFrac=1.e-4_r8 !seed for growth allocation only, [-]

  real(r8), parameter :: CMassCost4NiUptk=0.86_r8  !=12/14.
  logical :: checkCoarseRootLay
  character(len=*), parameter :: subname='SummarizeRootSink'

  associate(                                                                      &
    ZERO4Groth_pft                 => plt_biom%ZERO4Groth_pft                    ,& !input  :threshold zero for plang growth calculation, [-]
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft            ,& !input  :plant growth type (vascular, non-vascular),[-]
    VLSoilPoreMicP_vr              => plt_soilchem%VLSoilPoreMicP_vr             ,& !input  :volume of soil layer, [m3 d-2]
    NGTopRootLayer_pft             => plt_morph%NGTopRootLayer_pft               ,& !input  :soil layer at planting depth, [-]
    ZEROS2                         => plt_site%ZEROS2                            ,& !input  :threshold zero for numerical stability,[-]
    NU                             => plt_site%NU                                ,& !input  :current soil surface layer number, [-]
    CumSoilThickness_vr            => plt_site%CumSoilThickness_vr               ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    ZERO                           => plt_site%ZERO                              ,& !input  :threshold zero for numerical stability, [-]
    Root1stMaxRadius1_pft          => plt_morph%Root1stMaxRadius1_pft            ,& !input  :root radius primary axes, [m]
    DLYR3                          => plt_site%DLYR3                             ,& !input  :vertical thickness of soil layer, [m]
    RootOUlmNutUptake_pvr          => plt_rbgc%RootOUlmNutUptake_pvr             ,& !input  :root uptake of NH4 band unconstrained by O2, [g d-2 h-1]
    RootCUlmNutUptake_pvr          => plt_rbgc%RootCUlmNutUptake_pvr             ,& !input  :root uptake of NH4 band unconstrained by root nonstructural C, [g d-2 h-1]
    CanopyHeight4WatUptake_pft     => plt_morph%CanopyHeight4WatUptake_pft       ,& !input  :canopy height, [m]
    RootMatureAge_pft              => plt_morph%RootMatureAge_pft                ,& !input  : Root maturation age, [h]
    RootAge_rpvr                   => plt_morph%RootAge_rpvr                     ,& !input  : root age,[h]
    NRoot1stTipLay_raxes           => plt_morph%NRoot1stTipLay_raxes             ,& !input  :maximum soil layer number for root axes, [-]
    PSIRoot_pvr                    => plt_ew%PSIRoot_pvr                         ,& !input  :root total water potential, [Mpa]
    Myco_pft                       => plt_morph%Myco_pft                         ,& !input  :mycorrhizal type (no or yes),[-]
    Root1stRadius_rpvr             => plt_morph%Root1stRadius_rpvr               ,& !input  :root layer radius for each primary axes, [m]
    Cytokinin1stConc_rpvr          => plt_rbgc%Cytokinin1stConc_rpvr             ,& !output :cytokinin concentration in corase roots, [gC m-3 H2O]
    CRootLumenArea_rpvr            => plt_morph%CRootLumenArea_rpvr              ,& !output :coarse roots lumen area for root axes, [m2]
    iPlantTurnoverPattern_pft      => plt_pheno%iPlantTurnoverPattern_pft        ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    Root1stDepz_raxes              => plt_morph%Root1stDepz_raxes                ,& !input  :root layer depth, [m]
    HypocotHeight_pft              => plt_morph%HypocotHeight_pft                ,& !input  :cotyledon height, [m]
    iPlant2ndGrothPattern_pft      => plt_pheno%iPlant2ndGrothPattern_pft        ,& !input  :plant expression of secondary growth, [-]
    Root1stMaxRadius_pft           => plt_morph%Root1stMaxRadius_pft             ,& !input  :maximum radius of primary roots, [m]
    Root2ndRadius_rpvr             => plt_morph%Root2ndRadius_rpvr               ,& !input  :root layer radius secondary axes, [m]
    Root2ndXNum_rpvr               => plt_morph%Root2ndXNum_rpvr                 ,& !input  :root layer number secondary axes, [d-2]
    GroSrcRootStress_pvr           => plt_rbgc%GroSrcRootStress_pvr              ,& !input  :root growth stress due to nutrient and water, [-]
    Root2ndEffLen4uptk_rpvr        => plt_morph%Root2ndEffLen4uptk_rpvr          ,& !input  :Layer effective root length four resource uptake, [m]
    PlantPopuLive_pft              => plt_site%PlantPopuLive_pft                 ,& !input  :plant population, [d-2]
    DistRootEffDepz_pvr            => plt_morph%DistRootEffDepz_pvr              ,& !output :Effective shoot-root transport depth, [m]
    SapFlowVLinear_rpvr            => plt_ew%SapFlowVLinear_rpvr                 ,& !input  :linear sap flow for primary root axis, [m h-1]
    RootMyco2ndStrutElms_rpvr      => plt_biom%RootMyco2ndStrutElms_rpvr         ,& !inoput :root layer element secondary axes, [g d-2]
    SeedDepth_pft                  => plt_morph%SeedDepth_pft                    ,& !input  :seeding depth, [m]
    enh_cyto_pft                   => plt_morph%enh_cyto_pft                     ,& !input  :cytokinin sensitivity of coarse root thickening
    Cytokinin2ndConc_rpvr          => plt_rbgc%Cytokinin2ndConc_rpvr             ,& !input  :cytokinin concentration in fine roots, [gC m-3 H2O]
    MaxSoilLays4Root_pft           => plt_morph%MaxSoilLays4Root_pft             ,& !input  :maximum soil layer number for all root axes,[-]
    NumStructuralRootAxes_pft      => plt_morph%NumStructuralRootAxes_pft        ,& !input  :number of structural root axes,[-]
    RootMediumXNum_rpvr            => plt_morph%RootMediumXNum_rpvr              ,& !input  :number of medium root axes in soil layer, [# d-2]
    RootMediumStructElms_rpvr      => plt_biom%RootMediumStructElms_rpvr         ,& !input :population medium-root structural elements, [g d-2]
    Root1stLenPP_rpvr              => plt_morph%Root1stLenPP_rpvr                ,& !input :length per primary-root axis in layer, [m]
    MediumRootBranchFreq_pft       => plt_morph%MediumRootBranchFreq_pft         ,& !input :potential medium-root branching frequency, [m-1]
    Root2ndLen_rpvr                => plt_morph%Root2ndLen_rpvr                  ,& !input :population fine-root/mycorrhizal length, [m d-2]
    RootHPO4Uptake_pft             => plt_rbgc%RootHPO4Uptake_pft                ,& !inoput :total root uptake of HPO4, [g d-2 h-1]
    RootH2PO4Uptake_pft            => plt_rbgc%RootH2PO4Uptake_pft               ,& !inoput :total root uptake of PO4, [g d-2 h-1]
    RootNH4Uptake_pft              => plt_rbgc%RootNH4Uptake_pft                 ,& !inoput :total root uptake of NH4, [g d-2 h-1]
    RootNO3Uptake_pft              => plt_rbgc%RootNO3Uptake_pft                 ,& !inoput :total root uptake of NO3, [g d-2 h-1]
    fRootTube_rpvr                 => plt_morph%fRootTube_rpvr                   ,& !output :fraction of root for transport,[-]
    Num1stAxesPerStructRootX_pft   => plt_morph%Num1stAxesPerStructRootX_pft     ,& !input :primary root axes number, [d-2]
    Num1stAxesPerStructRootXPOP_pft=> plt_morph%Num1stAxesPerStructRootXPOP_pft  ,& !output :population primary root axes number on one structrual axis, [d-2]        
    RootMyco1stSinkC_rpvr          => plt_rbgc%RootMyco1stSinkC_rpvr             ,& !input  :primary root C sink, [gC d-2 h-1]
    RootMyco2ndSinkC_rpvr          => plt_rbgc%RootMyco2ndSinkC_rpvr             ,& !input  :fine root/myco carbon sink, [gC d-2 h-1]
    Root2ndSinkWeight_pvr          => plt_morph%Root2ndSinkWeight_pvr            ,& !output :Secondary root nonst element sink profile, [d-2]
    RootMediumMeanLength_rpvr      => plt_morph%RootMediumMeanLength_rpvr        ,& !input  :mean medium-root length per structural axis, [m]
    RootMediumLength_rpvr          => plt_morph%RootMediumLength_rpvr            ,& !input  :root layer length for medium size axes, [m d-2]      
    Root1stSinkWeight_pvr          => plt_morph%Root1stSinkWeight_pvr            ,& !output :primary root nonst element sink profile, [d-2]
    RootMSinkWeight_pvr            => plt_morph%RootMSinkWeight_pvr              ,& !output :medium size roots nonst element sink profile, [d-2]
    Root1stTipSinkWeight_pft       => plt_morph%Root1stTipSinkWeight_pft         ,& !output :primary root tip nonst element sink, [d-2]
    fctyok_scalar_rpvr             => plt_morph%fctyok_scalar_rpvr               ,& !output :cytokinin scalar for corase root sink, [-]
    RootNutUptakeN_pft             => plt_rbgc%RootNutUptakeN_pft                ,& !inoput :total N uptake by plant roots, [gN d-h2 h-1]
    flag2ndGrowth_pvr              => plt_morph%flag2ndGrowth_pvr                ,& !inoput  :flag for secondary growth of primary roots, [-]
    RootNutUptakeP_pft             => plt_rbgc%RootNutUptakeP_pft                ,& !inoput :total P uptake by plant roots, [gP d-h2 h-1]
    RootMycoNonstElms_rpvr         => plt_biom%RootMycoNonstElms_rpvr            ,& !inoput :root layer nonstructural element, [g d-2]
    RootNutUptake_pvr              => plt_rbgc%RootNutUptake_pvr                 ,& !inoput :root uptake of Nutrient band, [g d-2 h-1]
    RootCO2EmisPot_pvr             => plt_rbgc%RootCO2EmisPot_pvr                ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootRespPotent_pvr             => plt_rbgc%RootRespPotent_pvr                ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]
    RootCO2Autor_pvr               => plt_rbgc%RootCO2Autor_pvr                   & !inoput :root respiration constrained by O2, [g d-2 h-1]
  )

!     FOR ROOTS (N=1) AND MYCORRHIZAE (N=2) IN EACH SOIL LAYER
  call PrintInfo('beg '//subname)

  RootSinkC_vr               = 0._R8
  Root1stSink_pvr            = 0._r8
  Root2ndSink_pvr            = 0._r8
  RootMSink_pvr              = 0._r8  
  Root1stSinkTip             = 0._r8  
  RootSinkC                  = 0._r8
  RCO2flx                    = 0._r8
  fRootTube_rpvr             = 0._r8
  fctyok_scalar_rpvr(:,:,NZ) = 0._r8
  fUSE4MR_pvr                = 0._r8

  call GetInLayerRootDepth(NZ,Root1stLocDepz_vr,RootEffDepz_vr)

  call GetPrimaryRootPholTranspResist(yearIJ,NZ,Root1stLocDepz_vr,RootEffDepz_vr,R1stPholemResist_vr)

  call GetMediumRootPholTranspResit(yearIJ,NZ,SingleMedRootPholeResist_vr)

!  call SumRootBiome(yearIJ,NZ,mass_inital)
  flag2ndGrowth_pvr(:,:,NZ)=.false.
  !
  D4995: DO N=1,Myco_pft(NZ)
    Ltip=NU
    D4990: DO L=MaxSoilLays4Root_pft(NZ),NU,-1
      !
      !     RESPIRATION FROM NUTRIENT UPTAKE CALCULATED IN 'UPTAKE':
      !     ACTUAL, O2-UNLIMITED AND C-UNLIMITED
      !
      !     VLSoilPoreMicP_vr=soil layer volume excluding macropore, rocks
      !     RecoRootMycoC4Nup=C respiration for nutrient uptake
      !     CUPRO,CUPRC=RecoRootMycoC4Nup unlimited by O2,root nonstructural C
      !     why is 0.86? it refers to C cost for N asimilation
      
      IF(VLSoilPoreMicP_vr(L).GT.ZEROS2)THEN
        !accumulate C cost for nitrogen (NO3,NO3B,NH4,NH4B) and phosphorus (H1PO4,H2PO4,H1PO4B,H2PO4B) uptake

        RecoRootMycoC4Nup=0._r8;CUPRO=0._r8;CUPRC=0._r8
        do ntu=ids_nutb_beg+1,ids_nuts_end
          if(ntu==ids_NO2B .or. ntu==ids_NO2)cycle
          RecoRootMycoC4Nup = RecoRootMycoC4Nup+RootNutUptake_pvr(ntu,N,L,NZ)
          CUPRO             = CUPRO+RootOUlmNutUptake_pvr(ntu,N,L,NZ)
          CUPRC             = CUPRC+RootCUlmNutUptake_pvr(ntu,N,L,NZ)
        enddo  
        !
        !0.86 gC per gN or gP uptake 
        RecoRootMycoC4Nup = CMassCost4NiUptk*RecoRootMycoC4Nup
        CUPRO             = CMassCost4NiUptk*CUPRO
        CUPRC             = CMassCost4NiUptk*CUPRC
        !
        !
        !     RootCO2Autor_pvr=total root respiration
        !     RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
        !     RecoRootMycoC4Nup=C respiration for nutrient uptake
        !     CUPRO,CUPRC=RecoRootMycoC4Nup unlimited by O2,root nonstructural C
        !     CPOOLR=non-structural C mass in root
        !
        RootRespPotent_pvr(N,L,NZ) = RootRespPotent_pvr(N,L,NZ)+CUPRO
        RootCO2EmisPot_pvr(N,L,NZ) = RootCO2EmisPot_pvr(N,L,NZ)+CUPRC

        !correct small magnitude values to avoid numerical downflow
        RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=AZERO(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))

        if(RootMycoNonstElms_rpvr(ielmc,N,L,NZ).LT.RecoRootMycoC4Nup .and. RecoRootMycoC4Nup.GT.0._r8)then
          rscal=RootMycoNonstElms_rpvr(ielmc,N,L,NZ)/RecoRootMycoC4Nup*0.99999_r8          
          
          do ntu=ids_nutb_beg+1,ids_nuts_end
            !skip NO2(-) and NO2(-) band 
            if(ntu==ids_NO2B .or. ntu==ids_NO2)cycle
            RootNutUptake_pvr(ntu,N,L,NZ) = RootNutUptake_pvr(ntu,N,L,NZ)*rscal
          enddo
          RecoRootMycoC4Nup=RecoRootMycoC4Nup*rscal
        endif

        !root/myco non-structrual C used for root N/P uptake
        RootCO2Autor_pvr(N,L,NZ)             = RootCO2Autor_pvr(N,L,NZ)-RecoRootMycoC4Nup
        RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-RecoRootMycoC4Nup
        RCO2flx                              = RCO2flx-RecoRootMycoC4Nup
        !
        !     EXUDATION AND UPTAKE OF C, N AND P TO/FROM SOIL AND ROOT
        !     OR MYCORRHIZAL NON-STRUCTURAL C,N,P POOLS
        !
        !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
        !     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
        !     RootNH4Uptake_pvr,RootNH4BUptake_pvr,RUPN03,RootNO3BUptake_pvr=uptake from non-band,band of NH4,NO3
        !     RootH2PO4Uptake_pvr,RootH2PO4BUptake_pvr,RootHPO4Uptake_pvr,RootNutUptake_pvr=uptake from non-band,band of H2PO4,HPO4
        !
        dUptakeN= RootNutUptake_pvr(ids_NH4,N,L,NZ) +RootNutUptake_pvr(ids_NH4B,N,L,NZ) &
           +RootNutUptake_pvr(ids_NO3,N,L,NZ)+RootNutUptake_pvr(ids_NO3B,N,L,NZ)
        dUptakeP=RootNutUptake_pvr(ids_H2PO4,N,L,NZ)+RootNutUptake_pvr(ids_H2PO4B,N,L,NZ) &
           +RootNutUptake_pvr(ids_H1PO4,N,L,NZ)+RootNutUptake_pvr(ids_H1PO4B,N,L,NZ)

        RootNH4Uptake_pft(NZ)   = RootNH4Uptake_pft(NZ)  +(RootNutUptake_pvr(ids_NH4,N,L,NZ)+RootNutUptake_pvr(ids_NH4B,N,L,NZ))
        RootNO3Uptake_pft(NZ)   = RootNO3Uptake_pft(NZ)  +(RootNutUptake_pvr(ids_NO3,N,L,NZ)+RootNutUptake_pvr(ids_NO3B,N,L,NZ))
        RootH2PO4Uptake_pft(NZ) = RootH2PO4Uptake_pft(NZ) +(RootNutUptake_pvr(ids_H2PO4,N,L,NZ)+RootNutUptake_pvr(ids_H2PO4B,N,L,NZ))
        RootHPO4Uptake_pft(NZ)  = RootHPO4Uptake_pft(NZ)  +(RootNutUptake_pvr(ids_H1PO4,N,L,NZ)+RootNutUptake_pvr(ids_H1PO4B,N,L,NZ))

        RootMycoNonstElms_rpvr(ielmn,N,L,NZ)=RootMycoNonstElms_rpvr(ielmn,N,L,NZ)+dUptakeN
        RootMycoNonstElms_rpvr(ielmp,N,L,NZ)=RootMycoNonstElms_rpvr(ielmp,N,L,NZ)+dUptakeP

        RootNutUptakeN_pft(NZ)=RootNutUptakeN_pft(NZ)+dUptakeN
        RootNutUptakeP_pft(NZ)=RootNutUptakeP_pft(NZ)+dUptakeP
        
        IF(N.EQ.ipltroot)DistRootEffDepz_pvr(L,NZ)=0._r8
        D4985: DO NR=1,NumStructuralRootAxes_pft(NZ)
          !
          !do not go beyond lowest soil layer
          if(L.GT.NRoot1stTipLay_raxes(NR,NZ))cycle
          !
          !     PRIMARY ROOT SINK STRENGTH FROM ROOT RADIUS AND ROOT DEPTH
          !
          !          
          IF(N.EQ.ipltroot)THEN
            DistRootEffDepz           = RootEffDepz_vr(NR,L)+CanopyHeight4WatUptake_pft(NZ)
            DistRootEffDepz_pvr(L,NZ) = DistRootEffDepz_pvr(L,NZ)+DistRootEffDepz
            !
            !     SECONDARY ROOT SINK STRENGTH FROM ROOT RADIUS, ROOT AXIS NUMBER,
            !     AND ROOT LENGTH IN SERIES WITH PRIMARY ROOT SINK STRENGTH
            !
            Root1stTipRadius=AMAX1(Root1stMaxRadius1_pft(N,NZ),(1.0_r8+PSIRoot_pvr(N,L,NZ)/EMODR)*Root1stMaxRadius_pft(N,NZ))
            !

            !for coarse root dynamics, it must be woody plant, meets age and tip distance criteria, and 
            !lignification is on (assuming 2nd growth is modeled together with lignification)
            checkCoarseRootLay= is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ)) .and. &
                RootAge_rpvr(L,NR,NZ).GT.RootMatureAge_pft(NZ) .and. RootEffDepz_vr(NR,L).LT.Root1stDepz_raxes(NR,NZ) .and. &
                Root1stDepz_raxes(NR,NZ).GT.15._r8*Root1stTipRadius .and. lcoarseroot

            if(checkCoarseRootLay)then
              cytoc_scal = Cytokinin1stConc_rpvr(L,NR,NZ)     
              fctyok     = func_cyto(cytoc_scal,enh_cyto_pft(NZ))
            else
              fctyok=0._r8
            endif
            fctyok_scalar_rpvr(L,NR,NZ) = fctyok

            !primary root conductance for nonstructural biomass transport
            PrimaryRootCd4NonstTP = Num1stAxesPerStructRootXPOP_pft(NZ)/R1stPholemResist_vr(L,NR)            

            IF(DistRootEffDepz.GT.ZERO)THEN
              !fine root conductance for nonstructural biomass transport
              FineRootCd4NonstTP = safe_adb(Root2ndXNum_rpvr(N,L,NR,NZ)*Root2ndRadius_rpvr(N,L,NZ)**2,Root2ndEffLen4uptk_rpvr(N,L,NZ))
              !     
              MediumRootCd4NonstTP = 0._r8
              !Potential axes seed establishment without entering living transport geometry.
              PotentialMediumAxes = AZMAX1(MediumRootBranchFreq_pft(NZ)*Root1stLenPP_rpvr(L,NR,NZ) &
                *Num1stAxesPerStructRootXPOP_pft(NZ))
              MediumAxesForGrowth = 0._r8
              if(RootMediumStructElms_rpvr(ielmc,L,NR,NZ).GT.0._r8 .and. &
                 RootMediumLength_rpvr(L,NR,NZ).GT.0._r8) &
                MediumAxesForGrowth=AZMAX1(RootMediumXNum_rpvr(L,NR,NZ))
              if(MediumAxesForGrowth.LE.0._r8) &
                MediumAxesForGrowth=MinMediumEstablishmentFrac*PotentialMediumAxes
              IF(MediumAxesForGrowth.GT.0._r8 .and. FineRootCd4NonstTP.GT.ZERO4Groth_pft(NZ) .and. &
                 RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ).GT.ZERO4Groth_pft(NZ) .and. &
                 Root2ndLen_rpvr(N,L,NR,NZ).GT.0._r8 .and. PrimaryRootCd4NonstTP.GT.0._r8)then
                if(SingleMedRootPholeResist_vr(L,NR).GT.0._r8)then
                  MediumRootCd4NonstTP = MediumAxesForGrowth/SingleMedRootPholeResist_vr(L,NR)
                  ConnCd               = 1._r8/(1._r8/PrimaryRootCd4NonstTP+1._r8/MediumRootCd4NonstTP)
                  UsedConnCd           = 1._r8 / (1._r8/ConnCd + 1._r8/FineRootCd4NonstTP)
                  RootMSink_pvr(L,NR)  = f_mroot(NZ,L,NR,RootEffDepz_vr(NR,L))*UsedConnCd
                  !No activation floor: absent fine roots leave both growth quantities zero.
                  fUSE4MR_pvr(L,NR)    = MIN(1._r8,AZMAX1(UsedConnCd/ConnCd))
                endif
              endif

              IF(FineRootCd4NonstTP.GT.ZERO4Groth_pft(NZ))THEN
                ! Both lengths refer to an individual root; the mean is refreshed above.
                IF(MediumRootCd4NonstTP.GT.0._r8 .and. &
                   RootMediumStructElms_rpvr(ielmc,L,NR,NZ).GT.0._r8 .and. &
                   RootMediumLength_rpvr(L,NR,NZ).GT.0._r8 .and. &
                   RootMediumMeanLength_rpvr(L,NR,NZ).GT.Root2ndTipLen4uptk)then
                  !Use the same population lengths as fine-root branching and cytokinin routing:
                  !medium-root length and layer thickness times the population primary-axis count.
                  CMFR = 1._r8/(1._r8/MediumRootCd4NonstTP &
                    +(RootMediumLength_rpvr(L,NR,NZ)+DLYR3(L)*Num1stAxesPerStructRootXPOP_pft(NZ)) &
                    /(FineRootCd4NonstTP*RootMediumLength_rpvr(L,NR,NZ)))
                  CFR  = FineRootCd4NonstTP*DLYR3(L)*Num1stAxesPerStructRootXPOP_pft(NZ) &
                    /(RootMediumLength_rpvr(L,NR,NZ)+DLYR3(L)*Num1stAxesPerStructRootXPOP_pft(NZ))
                  Root2ndSink_pvr(N,L,NR)=(1._r8-fctyok)/(1._r8/PrimaryRootCd4NonstTP+1._r8/(CMFR+CFR))
                else                  
                  Root2ndSink_pvr(N,L,NR)=(1._r8-fctyok)*PrimaryRootCd4NonstTP*FineRootCd4NonstTP/(PrimaryRootCd4NonstTP+FineRootCd4NonstTP)  
                endif  
              ELSE
                Root2ndSink_pvr(N,L,NR)=0._r8
              ENDIF              
            ELSE
              Root2ndSink_pvr(N,L,NR) = 0._r8
              RootMSink_pvr(L,NR)     = 0._r8
            ENDIF
            !            
            Ltip=NRoot1stTipLay_raxes(NR,NZ)

            IF(Root1stDepz_raxes(NR,NZ).GT.CumSoilThickness_vr(L-1))THEN                                            
              If(checkCoarseRootLay)THEN
                !tip and corase root co-exist in layer
                flag2ndGrowth_pvr(L,NR,NZ) = .true.
                Root1stSink_pvr(L,NR)      = PrimaryRootCd4NonstTP*fctyok
                RootSinkC_vr(N,L)          = RootSinkC_vr(N,L)+Root1stSink_pvr(L,NR)+RootMSink_pvr(L,NR)
              ELSE
                Root1stSink_pvr(L,NR) = 0._r8                
              ENDIF
              !
              if(Root1stDepz_raxes(NR,NZ).LE.CumSoilThickness_vr(L))THEN
                NRoot1stTipLay_raxes(NR,NZ) = L
                Ltip                        = L
                flag2ndGrowth_pvr(L,NR,NZ) = checkCoarseRootLay .and. Root1stSink_pvr(L,NR).GT.0._r8
                Root1stSinkTip(NR)         = RTSK(iPlantRootProfile_pft(NZ))*Num1stAxesPerStructRootXPOP_pft(NZ)/R1stPholemResist_vr(L,NR)
                RootSinkC_vr(N,L)          = RootSinkC_vr(N,L)+Root1stSinkTip(NR)
              endif
            ENDIF            
          ELSE
            if(Root2ndSink_pvr(ipltroot,L,NR).GT.0._r8)then
              !mycorrhizae
              Root2ndSink_pvr(N,L,NR)=safe_adb(Root2ndXNum_rpvr(N,L,NR,NZ)*Root2ndRadius_rpvr(N,L,NZ)**2,Root2ndEffLen4uptk_rpvr(N,L,NZ))
            else
              Root2ndSink_pvr(N,L,NR)=0._r8
            endif
          ENDIF                   
          RootSinkC_vr(N,L) = RootSinkC_vr(N,L)+Root2ndSink_pvr(N,L,NR)
        ENDDO D4985       
        !get the mean root-shoot distance
        IF(N.EQ.ipltroot .and. DistRootEffDepz_pvr(L,NZ).GT.0._r8)DistRootEffDepz_pvr(L,NZ)=DistRootEffDepz_pvr(L,NZ)/NumStructuralRootAxes_pft(NZ) 
        RootSinkC(N)      = RootSinkC(N)+RootSinkC_vr(N,L)
      ENDIF
    ENDDO D4990
    
    Root2ndSinkWeight_pvr(:,N,NZ)=0._r8
    if(RootSinkC(N).GT.0._r8)then
      DO L=NU,MaxSoilLays4Root_pft(NZ)
        DO NR=1,NumStructuralRootAxes_pft(NZ)
          Root2ndSinkWeight_pvr(L,N,NZ)=Root2ndSinkWeight_pvr(L,N,NZ)+Root2ndSink_pvr(N,L,NR)
        ENDDO
        Root2ndSinkWeight_pvr(L,N,NZ)=Root2ndSinkWeight_pvr(L,N,NZ)/RootSinkC(N)        
      ENDDO
    endif
  ENDDO D4995

  Root1stSinkWeight_pvr(:,NZ)  = 0._r8
  Root1stTipSinkWeight_pft(NZ) = 0._r8
  RootMSinkWeight_pvr(:,NZ)    = 0._R8
  if(RootSinkC(ipltroot).GT.0._r8)then    
    DO L=NU,MaxSoilLays4Root_pft(NZ)
      DO NR=1,NumStructuralRootAxes_pft(NZ)
        Ltip=NRoot1stTipLay_raxes(NR,NZ)  
        Root1stSinkWeight_pvr(L,NZ)=Root1stSinkWeight_pvr(L,NZ)+Root1stSink_pvr(L,NR)
        if(L.eq.Ltip)then
          Root1stTipSinkWeight_pft(NZ) = Root1stTipSinkWeight_pft(NZ)+Root1stSinkTip(NR)
        endif
        RootMSinkWeight_pvr(L,NZ) =RootMSinkWeight_pvr(L,NZ)+ AZMAX1(RootMSink_pvr(L,NR))
      ENDDO
      RootMSinkWeight_pvr(L,NZ) =RootMSinkWeight_pvr(L,NZ)/RootSinkC(ipltroot)
      Root1stSinkWeight_pvr(L,NZ)=Root1stSinkWeight_pvr(L,NZ)/RootSinkC(ipltroot)
    ENDDO
    if(Root1stTipSinkWeight_pft(NZ).GT.0._r8) &
      Root1stTipSinkWeight_pft(NZ)=Root1stTipSinkWeight_pft(NZ)/RootSinkC(ipltroot)
  endif   
  
!  call SumRootBiome(yearIJ,NZ,mass_finale)
  !write(1114,*)yearIJ%I*1000+yearIJ%J/24.,1,2,0.0,plt_biom%RootMycoNonstElms_rpvr(ielmc,1,2,NZ),'af'//subname
  call PrintInfo('end '//subname)
  end associate
  end subroutine SummarizeRootSink
!----------------------------------------------------------------------------------------------------
  subroutine CalcCum2ndRootProfile(yearIJ,NZ,CRootSinkProf)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), intent(out) :: CRootSinkProf(JZ1,MaxNumRootAxes) !normalized coarse root sink profile based on pipe-theory

  real(r8) :: cumfrootmassC(JZ1) !depth-cumulative fine root C mass
  real(r8) :: SumFRootmassC
  INTEGER :: NR,L1,L,N

  character(len=*), parameter :: subname='CalcCum2ndRootProfile'
  associate(                                                              &
    Myco_pft                   => plt_morph%Myco_pft                     ,& !input  :mycorrhizal type (no or yes),[-]    
    NU                         => plt_site%NU                            ,& !input  :current soil surface layer number, [-]
    NK                         => plt_site%NK                            ,& !input  :current soil root bottom layer number, [-]        
    RootMyco2ndStrutElms_rpvr  => plt_biom%RootMyco2ndStrutElms_rpvr     ,& !input  :root layer element secondary axes, [g d-2]
    Root1stActStructElms_rpvr  => plt_biom%Root1stActStructElms_rpvr     ,& !input  :root layer active zone element in primary axes, [g d-2]
    MaxSoilLays4Root_pft       => plt_morph%MaxSoilLays4Root_pft         ,& !input  :maximum soil layer number for all root axes,[-]
    NRoot1stTipLay_raxes       => plt_morph%NRoot1stTipLay_raxes         ,& !inoput :soil layer number for deepest root axes, [-]    
    NumStructuralRootAxes_pft       => plt_morph%NumStructuralRootAxes_pft          & !input  :number of structural root axes,[-]
  )
  CALL PrintInfo('beg '//subname)

  CRootSinkProf(:,:) = 0._r8
  DO NR=1,NumStructuralRootAxes_pft(NZ)
    cumfrootmassC        = 0._r8
    SumFRootmassC        = 0._r8
    !sum up fine roots/myco
    DO L = NRoot1stTipLay_raxes(NR,NZ),NU,-1
      L1 = MIN(NK,L+1)      
      cumfrootmassC(L) = cumfrootmassC(L1)
      DO N=1,Myco_pft(NZ)
        cumfrootmassC(L) = cumfrootmassC(L)+ RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)
      ENDDO
      SumFRootmassC = SumFRootmassC+cumfrootmassC(L)
    ENDDO

    !root tip is not lignified, add to avoid division by zero
    DO L = NRoot1stTipLay_raxes(NR,NZ),NU,-1      
      cumfrootmassC(L)=cumfrootmassC(L)+Root1stActStructElms_rpvr(ielmc,NRoot1stTipLay_raxes(NR,NZ),NR,NZ)
    enddo

    if(SumFRootmassC.GT.0._r8)then
      DO L = NRoot1stTipLay_raxes(NR,NZ),NU,-1
        CRootSinkProf(L,NR) = cumfrootmassC(L)/cumfrootmassC(NRoot1stTipLay_raxes(NR,NZ))
      ENDDO
    endif
  ENDDO
  CALL PrintInfo('end '//subname)
  end associate
  end subroutine CalcCum2ndRootProfile
!----------------------------------------------------------------------------------------------------
  subroutine RootCheck(I,J,NZ,info)
  !
  !root axes are sorted by root order and then also by vertical layers
  !check primary root, with and without vertical distribution
  implicit none
  integer , intent(in) :: I,J,NZ
  character(len=*), intent(in) :: info
  integer :: NR,L,N,NE
  real(r8) :: dRootMyco1stElm_raxs(NumPlantChemElms)

  associate(                                                          &
    NumStructuralRootAxes_pft      => plt_morph%NumStructuralRootAxes_pft      ,& !input  :number of structural root axes,[-]
    Myco_pft                  => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    NL                        => plt_site%NL                         ,& !input  :lowest soil layer number,[-]
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr  ,& !input  :root layer element primary axes, [g d-2]
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs         & !input  :root C primary axes, [g d-2]
  )
  DO NR=1,NumStructuralRootAxes_pft(NZ)
    
    dRootMyco1stElm_raxs(:)=0._r8
    DO L=NU,NL
      DO NE=1,NumPlantChemElms
        dRootMyco1stElm_raxs(NE)=dRootMyco1stElm_raxs(NE)+RootMyco1stStrutElms_rpvr(NE,L,NR,NZ)  
      ENDDO
    ENDDO
!    write(111,*)trim(info),' rootcheck:',I+J/24.,NZ,N,NR,&
!        dRootMyco1stElm_raxs(:),RootMyco1stElm_raxs(:,NR,NZ),dRootMyco1stElm_raxs(:)-RootMyco1stElm_raxs(:,NR,NZ)
    dRootMyco1stElm_raxs(:)=dRootMyco1stElm_raxs(:)-RootMyco1stElm_raxs(:,NR,NZ)   
    if(any(abs(dRootMyco1stElm_raxs(:))>1.e-8_r8))stop
    
  ENDDO
  end associate
  end subroutine RootCheck
!----------------------------------------------------------------------------------------------------

  subroutine UpdateRootElongationZoneAge(I,J,NR,NZ,ElongateRate)
  !
  !Description
  !
  implicit none
  integer, intent(in) :: I,J,NR,NZ
  real(r8), intent(in)    :: ElongateRate  !Root length Growth rate (>0) or Retreat rate (<0)
  real(r8), parameter :: DLength=0.05_r8  !segment length [m]  
  real(r8), parameter :: dt=1._r8
  character(len=*), parameter :: subname='UpdateRootElongationZoneAge'
  real(r8) :: growth, retreat, current_len, current_age
  real(r8) :: new_total_len, new_mean_age, excess_len
  integer  :: k, curr_idx
  
  associate(                                                     &
    RootSegAges_raxes      => plt_morph%RootSegAges_raxes      , & !inoput: age of different active root segments
    NActiveRootSegs_raxes  => plt_morph%NActiveRootSegs_raxes  , & !inoput: number of active root segments
    RootSegBaseDepth_raxes => plt_morph%RootSegBaseDepth_raxes , & !inoput: base depth of different root axes
    NMaxRootSegs           => pltpar%NMaxRootSegs              , & !input:  Maximum number of root segments for age tracking
    RootSeglengths_raxes   => plt_morph%RootSeglengths_raxes   , & !inoput: root length in each segment
    IndRootSegBase_raxes   => plt_morph%IndRootSegBase_raxes   , & !inoput: Index for base segment under tracking
    IndRootSegTip_raxes    => plt_morph%IndRootSegTip_raxes      & !inoput: Index for tip segment under tracking
  )

  call PrintInfo('beg '//subname)
  ! ---------------------------------------------------------
  ! 1. AGING STEP
  ! ---------------------------------------------------------
  if (NActiveRootSegs_raxes(NR,NZ) > 0) then
    curr_idx = IndRootSegBase_raxes(NR,NZ)
    do k = 1, NActiveRootSegs_raxes(NR,NZ)
      RootSegAges_raxes(curr_idx,NR,NZ) = RootSegAges_raxes(curr_idx,NR,NZ) + dt
      
      ! Circular Step (1-based logic):
      ! ((Current - 1 + 1) % NMaxRootSegs) + 1
      curr_idx = mod(curr_idx, NMaxRootSegs) + 1
    end do
  end if

  ! ---------------------------------------------------------
  ! 2. HANDLE GROWTH (Positive Rate)
  ! ---------------------------------------------------------
  if (ElongateRate > 0.0_r8) then
    growth = ElongateRate * dt
    
    ! Initialize if empty
    if (NActiveRootSegs_raxes(NR,NZ) == 0) then
        NActiveRootSegs_raxes(NR,NZ)  = 1
        IndRootSegBase_raxes(NR,NZ)   = 1
        IndRootSegTip_raxes(NR,NZ)    = 1
        RootSeglengths_raxes(1,NR,NZ) = 0.0_r8
        RootSegAges_raxes(1,NR,NZ)    = 0.0_r8
    end if
      
    ! Get current tip data
    current_len = RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ)
    current_age = RootSegAges_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ)
    
    ! Weighted Mean Age Logic
    new_total_len = current_len + growth
    new_mean_age  = (current_age * current_len) / new_total_len
    
    ! Update Tip
    RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) = new_total_len
    RootSegAges_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ)    = new_mean_age
      
    ! Check for Overflow (Segment > DLength)
    if (RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) > DLength) then
      excess_len = RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) - DLength
      
      ! Cap the current full segment
      RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) = DLength
      
      ! Move Tip Forward (1-based circular)
      IndRootSegTip_raxes(NR,NZ) = mod(IndRootSegTip_raxes(NR,NZ), NMaxRootSegs) + 1
      
      ! Check for Buffer Full
      if (NActiveRootSegs_raxes(NR,NZ) < NMaxRootSegs) then
        NActiveRootSegs_raxes(NR,NZ) = NActiveRootSegs_raxes(NR,NZ) + 1
      else
        ! Buffer Full: Tip overwrites Base
        ! Base must move forward to escape
        ! Update RootSegBaseDepth_raxes(NR) because we lose the geometry of the old base
        RootSegBaseDepth_raxes(NR,NZ) = RootSegBaseDepth_raxes(NR,NZ) + RootSeglengths_raxes(IndRootSegBase_raxes(NR,NZ),NR,NZ)
        
        IndRootSegBase_raxes(NR,NZ) = mod(IndRootSegBase_raxes(NR,NZ), NMaxRootSegs) + 1
      end if
      
      ! Initialize the new tip
      RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) = excess_len
      RootSegAges_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ)    = 0.0_r8
    end if
      
  ! ---------------------------------------------------------
  ! 3. HANDLE RETREAT (Negative Rate)
  ! ---------------------------------------------------------
  else if (ElongateRate < 0.0_r8) then
    retreat = abs(ElongateRate) * dt
    
    do while (retreat > 0.0_r8 .and. NActiveRootSegs_raxes(NR,NZ) > 0)
      current_len = RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ)
      
      if (retreat >= current_len) then
        ! Consume entire segment
        retreat = retreat - current_len
        RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) = 0.0_r8
        RootSegAges_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ)    = 0.0_r8
        
        if (NActiveRootSegs_raxes(NR,NZ) > 1) then
          ! Move Tip Backward (1-based circular)
          ! Logic: ((Current - 1 - 1 + N) % N) + 1
          ! Simplified: ((Current - 2 + N) % N) + 1
          IndRootSegTip_raxes(NR,NZ) = mod(IndRootSegTip_raxes(NR,NZ) - 2 + NMaxRootSegs, NMaxRootSegs) + 1
          NActiveRootSegs_raxes(NR,NZ) = NActiveRootSegs_raxes(NR,NZ) - 1
        else
          ! Root is gone
          NActiveRootSegs_raxes(NR,NZ) = 0
          retreat                   = 0.0_r8
        end if
      else
        ! Partial consumption
        RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) = RootSeglengths_raxes(IndRootSegTip_raxes(NR,NZ),NR,NZ) - retreat
        retreat = 0.0_r8
      end if
    end do
  end if
  call PrintInfo('end '//subname)
  end associate      
  end subroutine UpdateRootElongationZoneAge
!----------------------------------------------------------------------------------------------------

  subroutine FindMatureRootSegs(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  character(len=*), parameter :: subname='FindMatureRootSegs'
  real(r8), parameter :: DLength=0.05_r8  !segment length [m]    
  real(r8), parameter :: RootLminMat=0.05_r8 !minimum length between growing tip and zone of secondary growth [m]
  INTEGER :: k,NR,KD
  logical :: foundMatureSeg

  associate(                                                      &
    NumStructuralRootAxes_pft    => plt_morph%NumStructuralRootAxes_pft   , & !input: number of structural root axes,[-]  
    RootSegAges_raxes       => plt_morph%RootSegAges_raxes      , & !input: age of different active root segments
    NActiveRootSegs_raxes   => plt_morph%NActiveRootSegs_raxes  , & !inoput: number of active root segments
    RootSegBaseDepth_raxes  => plt_morph%RootSegBaseDepth_raxes , & !inoput: base depth of different root axes, [m]
    Root1stDepz_raxes       => plt_morph%Root1stDepz_raxes      , & !input : root layer depth, [m]    
    NMaxRootSegs            => pltpar%NMaxRootSegs              , & !input : Maximum number of root segments for age tracking
    RootMatureAge_pft       => plt_morph%RootMatureAge_pft      , & !input : Root maturation age, [h]
    SeedDepth_pft           => plt_morph%SeedDepth_pft          , & !input  :seeding depth, [m]
    RootSeglengths_raxes    => plt_morph%RootSeglengths_raxes   , & !inoput: root length in each segment, [m]
    IndRootSegBase_raxes    => plt_morph%IndRootSegBase_raxes   , & !inoput: Index for base segment under tracking
    IndRootSegTip_raxes     => plt_morph%IndRootSegTip_raxes      & !inoput: Index for tip segment under tracking
  )
  call PrintInfo('beg '//subname)
  DO NR=1,NumStructuralRootAxes_pft(NZ)
    foundMatureSeg=.false.

    !the tracking array is of shape [1,2,...->tip,...,base->,...,max]
    IF(IndRootSegBase_raxes(NR,NZ) > IndRootSegTip_raxes(NR,NZ))then
      DO K=IndRootSegTip_raxes(NR,NZ),1,-1
        if(RootSegAges_raxes(K,NR,NZ)>=RootMatureAge_pft(NZ))then
          foundMatureSeg=.true.
          KD=K-1
          RootSegBaseDepth_raxes(NR,NZ) = RootSegBaseDepth_raxes(NR,NZ)+DLength*KD
          RootSegBaseDepth_raxes(NR,NZ) = min(RootSegBaseDepth_raxes(NR,NZ)+RootSeglengths_raxes(K,NR,NZ),Root1stDepz_raxes(NR,NZ))
          IF(Root1stDepz_raxes(NR,NZ).LT.RootLminMat)RootSegBaseDepth_raxes(NR,NZ) = SeedDepth_pft(NZ)            

          RootSeglengths_raxes(1:K,NR,NZ) = 0._r8
          RootSegAges_raxes(1:K,NR,NZ)    = 0._r8
          IndRootSegBase_raxes(NR,NZ)     = K+1
          NActiveRootSegs_raxes(NR,NZ)    = NActiveRootSegs_raxes(NR,NZ)-KD-1
          exit
        ENDIF
      ENDDO
      if(foundMatureSeg)then
        if(IndRootSegBase_raxes(NR,NZ)>IndRootSegTip_raxes(NR,NZ))then
          !the elongation zone is dead.
          NActiveRootSegs_raxes(NR,NZ)=0
          RootSegBaseDepth_raxes(NR,NZ)=min(RootSegBaseDepth_raxes(NR,NZ),Root1stDepz_raxes(NR,NZ))
          IF(Root1stDepz_raxes(NR,NZ).LT.RootLminMat)RootSegBaseDepth_raxes(NR,NZ) = SeedDepth_pft(NZ)            
        endif
        cycle      
      endif  

      DO K=NMaxRootSegs,IndRootSegBase_raxes(NR,NZ),-1
        if(RootSegAges_raxes(K,NR,NZ)>=RootMatureAge_pft(NZ))then
          foundMatureSeg                = .true.
          Kd=(K-IndRootSegBase_raxes(NR,NZ))
          RootSegBaseDepth_raxes(NR,NZ) = RootSegBaseDepth_raxes(NR,NZ)+DLength*Kd
          RootSegBaseDepth_raxes(NR,NZ) = min(RootSegBaseDepth_raxes(NR,NZ)+RootSeglengths_raxes(K,NR,NZ),Root1stDepz_raxes(NR,NZ))
          IF(Root1stDepz_raxes(NR,NZ).LT.RootLminMat)RootSegBaseDepth_raxes(NR,NZ) = SeedDepth_pft(NZ)            

          RootSegAges_raxes(IndRootSegBase_raxes(NR,NZ):K,NR,NZ)    = 0._r8
          RootSeglengths_raxes(IndRootSegBase_raxes(NR,NZ):K,NR,NZ) = 0._r8

          IndRootSegBase_raxes(NR,NZ)  = K+1
          NActiveRootSegs_raxes(NR,NZ) = NActiveRootSegs_raxes(NR,NZ)-KD-1

          exit  
        endif
      ENDDO
    else
      DO K=IndRootSegTip_raxes(NR,NZ),IndRootSegBase_raxes(NR,NZ),-1
        if(RootSegAges_raxes(K,NR,NZ)>=RootMatureAge_pft(NZ))then
          foundMatureSeg=.true.
          KD=K-IndRootSegBase_raxes(NR,NZ)
          RootSegBaseDepth_raxes(NR,NZ) = RootSegBaseDepth_raxes(NR,NZ)+DLength*KD          
          RootSegBaseDepth_raxes(NR,NZ) = RootSegBaseDepth_raxes(NR,NZ)+RootSeglengths_raxes(K,NR,NZ)
          IF(Root1stDepz_raxes(NR,NZ).LT.RootLminMat)RootSegBaseDepth_raxes(NR,NZ) = SeedDepth_pft(NZ)            

          RootSegAges_raxes(IndRootSegBase_raxes(NR,NZ):K,NR,NZ)    = 0._r8
          RootSeglengths_raxes(IndRootSegBase_raxes(NR,NZ):K,NR,NZ) = 0._r8

          NActiveRootSegs_raxes(NR,NZ) = NActiveRootSegs_raxes(NR,NZ)-KD-1
          IndRootSegBase_raxes(NR,NZ)  = K+1
          exit
        ENDIF
      ENDDO
    endif
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine FindMatureRootSegs  
!----------------------------------------------------------------------------------------------------
  subroutine ComputeSoilStress(I,J,Pstress_vr)
  implicit none
  INTEGER, intent(in) :: I,J
  real(r8), intent(out) :: Pstress_vr(JZ1)     !Effective Vertical Stress
  integer :: L
  character(len=*), parameter :: subname='ComputeSoilStress'
  associate(                                                     &
    NU                         => plt_site%NU                  , & !input : current soil surface layer number, [-]
    SoilWeightStress_vr        => plt_site%SoilWeightStress_vr , & !input : soil weight stress on root thickening, [MPa]
    SoilSuctStress_vr          => plt_site%SoilSuctStress_vr   , & !input : soil suction stress on root thickening, [MPa]
    rSat_vr                    => plt_site%rSat_vr               & !input : relative soil saturation, [-]            
  )
  call PrintInfo('beg '//subname)

  Pstress_vr(1)=SoilWeightStress_vr(1)*0.5_r8+rSat_vr(1)*abs(SoilSuctStress_vr(1))

  DO L=2,JZ1
    Pstress_vr(L)=Pstress_vr(L-1)+(SoilWeightStress_vr(L-1)+SoilWeightStress_vr(L))*0.5_r8+rSat_vr(L)*abs(SoilSuctStress_vr(L))
  enddo
  call PrintInfo('end '//subname)

  end associate
  end subroutine ComputeSoilStress
!----------------------------------------------------------------------------------------------------

  pure function CRootVolInflator(CRootRadius_scaled,gamma)result(ans)
  !
  !(the coarse + intermediate root active volume) /(coarse root active volume)
  implicit none
  real(r8), intent(in) :: CRootRadius_scaled !
  real(r8), intent(in) :: gamma              !scaling factor [0.5, 1)
  real(r8) :: ans,x,denorm1,denorm2
  real(r8) :: d,s
  real(r8), parameter :: d0=2._r8

  if(CRootRadius_scaled.LE.1._r8)then
    ans=1._r8
  else  
    x       = CRootRadius_scaled
    d       = x-d0
    s       = d*d/(d*d+d0*d0)
    denorm1 = 2._r8*(x**(gamma-1._r8)-1._r8)/(gamma-1._r8)
    denorm2 = (x**(gamma-2._r8)-1._r8)/(gamma-2._r8)
    ans     = 1._r8+s/(gamma*(denorm1-denorm2))
  endif  
  end function CRootVolInflator
!----------------------------------------------------------------------------------------------------
  real(r8) function f_mroot(NZ,L,NR,RootEffDepz)
  implicit none
  INTEGER, intent(in) :: NZ,NR,L
  real(r8), intent(in) :: RootEffDepz

  REAL(R8), parameter :: z_full_ref = 0.40_r8 ![m]
  real(r8), parameter :: x_full_ref = 0.20_r8 ![m] behind tip
  real(r8), parameter :: day2hour = 24._r8 ![h]
  real(r8) :: RTSK_val
  real(r8) :: dist_behind_tip
  real(r8) :: z_full,z_on
  real(r8) :: x_full,x_on
  real(r8) :: age_on,age_full

  associate(                                                         &
    iPlantRootProfile_pft => plt_pheno%iPlantRootProfile_pft        ,& !input :plant growth type (vascular, non-vascular),[-]
    Root1stDepz_raxes     => plt_morph%Root1stDepz_raxes            ,& !input :root layer depth, [m]
    RootMatureAge_pft     => plt_morph%RootMatureAge_pft            ,& !input :Root maturation age, [h]
    RootAge_rpvr          => plt_morph%RootAge_rpvr                  & !input :root age,[h]
  )

  RTSK_val        = RTSK(iPlantRootProfile_pft(NZ))
  dist_behind_tip = Root1stDepz_raxes(NR,NZ) - RootEffDepz

  z_full = z_full_ref * RTSK_val**0.5_r8
  z_on   = 0.3_r8 * z_full

  x_full = x_full_ref * RTSK_val**0.25_r8
  x_on   = 0.3_r8 * x_full

  age_on   = AZMAX1(RootMatureAge_pft(NZ),      30._r8 * day2hour) * RTSK_val**0.25_r8
  age_full = AZMAX1(3._r8*RootMatureAge_pft(NZ), 90._r8 * day2hour) * RTSK_val**0.25_r8

  f_mroot = smoothstep(z_on, z_full, Root1stDepz_raxes(NR,NZ)) * & !rooting depth
    smoothstep(x_on, x_full, dist_behind_tip) * &                  !behind root tip
    smoothstep(age_on, age_full, RootAge_rpvr(L,NR,NZ))            !rooting age

  end associate
  end function f_mroot

  ![tail]
end module RootMod
