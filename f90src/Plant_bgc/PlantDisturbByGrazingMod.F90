module PlantDisturbByGrazingMod
  use data_kind_mod, only : r8 => DAT_KIND_R8,yearIJ_type
  use minimathmod, only : isclose,AZMAX1
  use ElmIDMod
  use EcosimConst
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantBiomassAPIData, only : plt_biom
  use PlantDisturbanceAPIData, only : plt_distb
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantBGCPars
implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
! end_include_section
  public :: AbvgBiomRemovalByGrazing
  public :: RemoveStandDeadByGrazing
  public :: ApplyBiomRemovalByGrazing
  public :: CutBranchNonstalByGrazing
  public :: GrazingPlant

! GY,GZ=partitioning of grazed material to removal,respiration
!
  real(r8),PARAMETER :: GY=1._r8
  real(r8),parameter :: GZ=1._r8-GY  !respiraiton ratio of grazers
contains
![header]
!------------------------------------------------------------------------------------------

  subroutine AbvgBiomRemovalByGrazing(yearIJ,NZ,TotalElmnt2Litr,TotalElmntRemoval)
  
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), intent(in) :: TotalElmntRemoval(NumPlantChemElms)  
  real(r8), intent(in) :: TotalElmnt2Litr(NumPlantChemElms)  
  character(len=*), parameter :: subname='AbvgBiomRemovalByGrazing'
  integer :: NE
  real(r8) :: dRespC

  associate(                                                      &
    Eco_AutoR_CumYr_col     => plt_bgcr%Eco_AutoR_CumYr_col,      &
    CanopyRespC_CumYr_pft   => plt_bgcr%CanopyRespC_CumYr_pft,    &
    ECO_ER_col              => plt_bgcr%ECO_ER_col,               &
    GrossResp_pft           => plt_bgcr%GrossResp_pft,            &
    EcoHavstElmnt_CumYr_col => plt_distb%EcoHavstElmnt_CumYr_col, &
    EcoHavstElmnt_CumYr_pft => plt_distb%EcoHavstElmnt_CumYr_pft  &
  )
  call PrintInfo('beg '//subname)
  !     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
  !     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem from all PFT
  !     GY=growth yield of grazers
  !     TotalElmntRemoval(ielmc),TotalElmntRemoval(ielmn),TotalElmntRemoval(ielmp)=total C,N,P removed
  !     TotalElmnt2Litr(ielmc),TotalElmnt2Litr(ielmn),TotalElmnt2Litr(ielmp)=total C,N,P to litter
  !     ECO_ER_col=ecosystem respiration
  !     Eco_AutoR_CumYr_col=total autotrophic respiration

  EcoHavstElmnt_CumYr_pft(ielmc,NZ) = EcoHavstElmnt_CumYr_pft(ielmc,NZ)+GY*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
  EcoHavstElmnt_CumYr_col(ielmc)    = EcoHavstElmnt_CumYr_col(ielmc)+GY*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))

  DO NE=2,NumPlantChemElms
    EcoHavstElmnt_CumYr_pft(NE,NZ) = EcoHavstElmnt_CumYr_pft(NE,NZ)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
    EcoHavstElmnt_CumYr_col(NE)    = EcoHavstElmnt_CumYr_col(NE)+TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
  ENDDO
  !GZ : fraction of removal into respiration by herbivory/grazing
  dRespC                    = GZ*(TotalElmntRemoval(ielmc)-TotalElmnt2Litr(ielmc))
  GrossResp_pft(NZ)         = GrossResp_pft(NZ)-dRespC
  CanopyRespC_CumYr_pft(NZ) = CanopyRespC_CumYr_pft(NZ)-dRespC
!     Eco_NBP_CumYr_col=Eco_NBP_CumYr_col+GY*(TotalElmnt2Litr(ielmc)-TotalElmntRemoval(ielmc))
!     CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-dRespC
  ECO_ER_col          = ECO_ER_col-dRespC
  Eco_AutoR_CumYr_col = Eco_AutoR_CumYr_col-dRespC
  call PrintInfo('end '//subname)
  end associate
  end subroutine AbvgBiomRemovalByGrazing

!------------------------------------------------------------------------------------------

  subroutine RemoveStandDeadByGrazing(yearIJ,NZ,FracStdeadLeft,FHVSH)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer,  intent(in) :: NZ
  real(r8), intent(out) :: FracStdeadLeft  !fraction of remaining C after harvest
  real(r8), intent(out) :: FHVSH  !
  real(r8) :: HarvestedStdeadC              !harvested standing dead C, [gC d-2]

  associate(                                                      &
    CanopyCutProxy_pft      => plt_distb%CanopyCutProxy_pft,      &
    NU                      => plt_site%NU,                       &
    AREA3                   => plt_site%AREA3,                    &
    StandDeadStrutElms_pft  => plt_biom%StandDeadStrutElms_pft,   &
    FracBiomHarvsted        => plt_distb%FracBiomHarvsted,        &
    THIN_pft                => plt_distb%THIN_pft,                &
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft            &
  )

  IF(StandDeadStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
    HarvestedStdeadC = CanopyCutProxy_pft(NZ)*THIN_pft(NZ)*0.45_r8/24.0_r8*AREA3(NU)*FracBiomHarvsted(iHarvst_pft,iplthvst_stdead,NZ)
    FracStdeadLeft   = AZMAX1(1._r8-HarvestedStdeadC/StandDeadStrutElms_pft(ielmc,NZ))
    FHVSH  = FracStdeadLeft
  ELSE
    FracStdeadLeft = 1.0_r8
    FHVSH          = 1.0_r8
  ENDIF
  end associate  
  end subroutine RemoveStandDeadByGrazing
!------------------------------------------------------------------------------------------

  subroutine ApplyBiomRemovalByGrazing(I,J,NZ,EHVST21,EHVST22,EHVST23,EHVST24,&
    CanopyNonstElmRemoval,LeafElmntRemoval,FineNonleafElmntRemoval,WoodyElmntRemoval,StandeadElmntRemoval,&
    CanopyNonstElm2Litr,LeafElmnt2Litr,FineNonleafElmnt2Litr,WoodyElmnt2Litr,StandeadElmnt2Litr)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  real(r8), intent(in) :: EHVST21
  real(r8), intent(in) :: EHVST22
  real(r8), intent(in) :: EHVST23
  real(r8), intent(in) :: EHVST24
  real(r8), intent(in) :: CanopyNonstElmRemoval(NumPlantChemElms)
  real(r8), intent(in) :: LeafElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: FineNonleafElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: WoodyElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: StandeadElmntRemoval(NumPlantChemElms)
  real(r8), intent(out) :: CanopyNonstElm2Litr(NumPlantChemElms)  
  real(r8), intent(out) :: LeafElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: FineNonleafElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: WoodyElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: StandeadElmnt2Litr(NumPlantChemElms)
  real(r8) :: EHVST21h,EHVST22h,EHVST23h,EHVST24h
  integer :: NE
  associate(                                        &
    NU               => plt_site%NU,                &
    AREA3            => plt_site%AREA3,             &
    FERT             => plt_distb%FERT,             &
    FracBiomHarvsted => plt_distb%FracBiomHarvsted, &
    IYTYP            => plt_distb%IYTYP             &
  )
  EHVST21h = 1._r8-FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ)*0.5_r8
  EHVST22h = 1._r8-FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ)*0.5_r8
  EHVST23h = 1._r8-FracBiomHarvsted(iHarvst_col,iplthvst_stalk,NZ)*0.5_r8
  EHVST24h = 1._r8-FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ)*0.5_r8

  CanopyNonstElm2Litr(ielmc)   = CanopyNonstElmRemoval(ielmc)*EHVST21
  LeafElmnt2Litr(ielmc)        = LeafElmntRemoval(ielmc)*EHVST21
  FineNonleafElmnt2Litr(ielmc) = FineNonleafElmntRemoval(ielmc)*EHVST22
  WoodyElmnt2Litr(ielmc)       = WoodyElmntRemoval(ielmc)*EHVST23
  StandeadElmnt2Litr(ielmc)    = StandeadElmntRemoval(ielmc)*EHVST24

  DO NE=2,NumPlantChemElms
    CanopyNonstElm2Litr(NE)   = CanopyNonstElmRemoval(NE)*EHVST21h
    LeafElmnt2Litr(NE)        = LeafElmntRemoval(NE)*EHVST21h
    FineNonleafElmnt2Litr(NE) = FineNonleafElmntRemoval(NE)*EHVST22h
    WoodyElmnt2Litr(NE)       = WoodyElmntRemoval(NE)*EHVST23h
    StandeadElmnt2Litr(NE)    = StandeadElmntRemoval(NE)*EHVST24h
  ENDDO
  !
  !     ADD MANURE FROM GRAZING TO NEXT DAY FERTILIZER
  !
  !     FERT=fertilizer type from fertilizer input file
  !     IYTYP=fertilizer release type from fertilizer input file
  !
  FERT(ifert_plant_manuC)=FERT(ifert_plant_manuC)+(CanopyNonstElm2Litr(ielmc)+LeafElmnt2Litr(ielmc)+&
    FineNonleafElmnt2Litr(ielmc)+WoodyElmnt2Litr(ielmc)+StandeadElmnt2Litr(ielmc))/AREA3(NU)
  FERT(ifert_plant_manuN)=FERT(ifert_plant_manuN)+(CanopyNonstElm2Litr(ielmn)+LeafElmnt2Litr(ielmn)+&
    FineNonleafElmnt2Litr(ielmn)+WoodyElmnt2Litr(ielmn)+StandeadElmnt2Litr(ielmn))/AREA3(NU)*0.5_r8    
  FERT(ifert_N_urea)=FERT(ifert_N_urea)+(CanopyNonstElm2Litr(ielmn)+LeafElmnt2Litr(ielmn)+&
    FineNonleafElmnt2Litr(ielmn)+WoodyElmnt2Litr(ielmn)+StandeadElmnt2Litr(ielmn))/AREA3(NU)*0.5_r8
  FERT(ifert_plant_manuP)=FERT(ifert_plant_manuP)+(CanopyNonstElm2Litr(ielmp)+LeafElmnt2Litr(ielmp)+&
    FineNonleafElmnt2Litr(ielmp)+WoodyElmnt2Litr(ielmp)+StandeadElmnt2Litr(ielmp))/AREA3(NU)
  IYTYP=imanure_grazing
  end associate
  end subroutine ApplyBiomRemovalByGrazing


!--------------------------------------------------------------------------------

  subroutine GrazingPlant(I,J,NZ,HarvestedLeafC,HarvestedShethC,HarvestedEarC,HarvestedGrainC,&
    GrazedCanopyNonstC,HarvestedStalkC,HarvestedStalkRsrvC,HarvestedPetoleC,GrazedCanopyNoduleC,LeafLayerC_brch)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(out):: HarvestedLeafC
  real(r8), intent(out):: HarvestedShethC
  real(r8), intent(out):: HarvestedEarC
  real(r8), intent(out):: HarvestedGrainC
  real(r8), intent(out):: GrazedCanopyNonstC
  real(r8), intent(out):: HarvestedStalkC
  real(r8), intent(out):: HarvestedStalkRsrvC
  real(r8), intent(out) :: HarvestedPetoleC
  real(r8), intent(out) :: GrazedCanopyNoduleC  
  REAL(R8), intent(out) :: LeafLayerC_brch(NumCanopyLayers1,JP1,JP1)  
  real(r8) :: totShootC,TotPhytomassRemoved
  real(r8) :: GrazingDemndLeafC,GrazingMeetLeafC,GrazedLeafNonstC,GrazedLeafNoduleC,GrazingUnMetLeafC,GrazedPhytoLeafC_pft
  real(r8) :: GrazedDmndStalkC  
  real(r8) :: TotStalkC  
  real(r8) :: GrazingMeetPetoleC,GrazedPetoleNonstC,GrazedPetoleNoduleC,GrazingDmndEarC,GrazedEarC,GrazingDmdGrainC,GrazedGrainC  
  real(r8) :: GrazingDmndPetole,GrazingDmndHuskC,GrazedHuskC,GrazingDmndStalkRsrvC
  real(r8) :: GrazedStalkRsrvC,GrazedPhytoStalkC_pft,GrazingDmndStalkC
  real(r8) :: CCPOLX  !grazing preference based on nonstructural C
  real(r8) :: CCPLNX  !grazing preference based on canopy nodule amount
  integer :: K,NB,L

  associate(                                                             &
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft,          &
    LeafLayerElms_node          => plt_biom%LeafLayerElms_node,          &
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft,              &
    CanopyCutProxy_pft         => plt_distb%CanopyCutProxy_pft,        &
    iHarvstType_pft             => plt_distb%iHarvstType_pft,            &
    LeafStrutElms_pft           => plt_biom%LeafStrutElms_pft,           &
    HuskStrutElms_pft           => plt_biom%HuskStrutElms_pft,           &
    ShootElms_pft               => plt_biom%ShootElms_pft,               &
    AREA3                       => plt_site%AREA3,                       &
    NU                          => plt_site%NU,                          &
    FracBiomHarvsted            => plt_distb%FracBiomHarvsted,           &
    THIN_pft                    => plt_distb%THIN_pft,                   &
    CanopyNoduleNonstCConc_pft  => plt_biom%CanopyNoduleNonstCConc_pft,  &
    GrainStrutElms_pft          => plt_biom%GrainStrutElms_pft,          &
    StalkStrutElms_pft          => plt_biom%StalkStrutElms_pft,          &
    StalkRsrvElms_pft           => plt_biom%StalkRsrvElms_pft,           &
    EarStrutElms_pft            => plt_biom%EarStrutElms_pft,            &
    PetolShethStrutElms_pft         => plt_biom%PetolShethStrutElms_pft,         &
    CanopyNonstElmConc_pft      => plt_biom%CanopyNonstElmConc_pft,      &
    fTCanopyGroth_pft           => plt_pheno%fTCanopyGroth_pft,          &
    AvgCanopyBiomC2Graze_pft    => plt_biom%AvgCanopyBiomC2Graze_pft     &
  )
  !
  !     AvgCanopyBiomC2Graze_pft=section mean shoot C density times local cell area
  !     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
  !          iHarvstType_pft=3:reduction of clumping factor
  !          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
  !     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
  !          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
  !     TotPhytomassRemoved=total phytomass grazed, removed
  !     fTCanopyGroth_pft=temperature function for canopy growth
  !     CCPOLP=nonstructural C concentration in canopy
  !     CanopyNoduleNonstCConc_pft=nonstructural C concentration in canopy nodules
  ! 24. has unit of hour, 0.45 is biomass C fraction
  ! Both shoot C and its section mean are scaled to this cell area.
  IF(AvgCanopyBiomC2Graze_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
    TotPhytomassRemoved=CanopyCutProxy_pft(NZ)*THIN_pft(NZ)*0.45_r8/24.0_r8 &
      *AREA3(NU)*ShootElms_pft(ielmc,NZ)/AvgCanopyBiomC2Graze_pft(NZ)
  ELSE
    TotPhytomassRemoved=0._r8
  ENDIF

  !herbivory, means fauna/insect cutting?
  IF(iHarvstType_pft(NZ).EQ.iharvtyp_herbivo)THEN
    TotPhytomassRemoved=TotPhytomassRemoved*fTCanopyGroth_pft(NZ)
  ENDIF
  CCPOLX = CanopyNonstElmConc_pft(ielmc,NZ)/(1.0_r8+CanopyNonstElmConc_pft(ielmc,NZ))
  CCPLNX = CanopyNoduleNonstCConc_pft(NZ)/(1.0_r8+CanopyNoduleNonstCConc_pft(NZ))
  !
  !     LEAF,BACTERIA GRAZED,REMOVED
  !
  !     FracBiomHarvsted(iHarvst_pft,1:4)=fraction of leaf,non-foliar,woody, standing dead removed from PFT
  !     FracBiomHarvsted(iHarvst_col,1:4=fraction of leaf,non-foliar,woody, standing dead removed from ecosyst
  !     WHVSL*,WHVSC*,WHVSN=leaf,nonstructural,bacteria removed
  !     WTLF=PFT leaf C mass
  !     GrazingUnMetLeafC=grazing requirement unmet by leaf
  ! grazing from leaf, 
  GrazingDemndLeafC = TotPhytomassRemoved*FracBiomHarvsted(iHarvst_pft,iplthvst_leaf,NZ)
  GrazingMeetLeafC  = AMIN1(LeafStrutElms_pft(ielmc,NZ),GrazingDemndLeafC)
  HarvestedLeafC       = GrazingMeetLeafC*(1._r8-CCPOLX)
  
  GrazedLeafNonstC     = GrazingMeetLeafC*CCPOLX
  GrazedLeafNoduleC    = GrazingMeetLeafC*CCPLNX
  GrazingUnMetLeafC    = AZMAX1(GrazingDemndLeafC-GrazingMeetLeafC)
  GrazedPhytoLeafC_pft = TotPhytomassRemoved*FracBiomHarvsted(iHarvst_pft,iplthvst_finenonleaf,NZ)
  !
  !     OTHER NON-FOLIAR GRAZED,REMOVED
  !
  !     WTSHE,WTHSK,WTEAR,WTGR=PFT PetolSheth,husk,ear,grain C mass
  !     WHVSH*,WHVHS*,WHVEA*,WHVGR*,WHVSC*=
  !            PetolSheth,husk,ear,grain,nonstructural C removed
  !     GrazingUnMetLeafC=grazing requirement unmet by non-foliar removal
  !
  totShootC=PetolShethStrutElms_pft(ielmc,NZ)+HuskStrutElms_pft(ielmc,NZ)+EarStrutElms_pft(ielmc,NZ)+GrainStrutElms_pft(ielmc,NZ)
  
  IF(totShootC.GT.ZERO4Groth_pft(NZ))THEN
    GrazingDmndPetole   = GrazedPhytoLeafC_pft*PetolShethStrutElms_pft(ielmc,NZ)/totShootC+GrazingUnMetLeafC
    GrazingMeetPetoleC  = AMIN1(PetolShethStrutElms_pft(ielmc,NZ),GrazingDmndPetole)
    HarvestedPetoleC    = GrazingMeetPetoleC*(1._r8-CCPOLX)
    GrazedPetoleNonstC  = GrazingMeetPetoleC*CCPOLX
    GrazedPetoleNoduleC = GrazingMeetPetoleC*CCPLNX

    GrazingUnMetLeafC = AZMAX1(GrazingDmndPetole-GrazingMeetPetoleC)
    GrazingDmndHuskC  = GrazedPhytoLeafC_pft*HuskStrutElms_pft(ielmc,NZ)/totShootC+GrazingUnMetLeafC
    GrazedHuskC       = AMIN1(HuskStrutElms_pft(ielmc,NZ),GrazingDmndHuskC)
    HarvestedShethC   = GrazedHuskC

    GrazingUnMetLeafC = AZMAX1(GrazingDmndHuskC-GrazedHuskC)
    GrazingDmndEarC   = GrazedPhytoLeafC_pft*EarStrutElms_pft(ielmc,NZ)/totShootC+GrazingUnMetLeafC
    GrazedEarC        = AMIN1(EarStrutElms_pft(ielmc,NZ),GrazingDmndEarC)
    HarvestedEarC     = GrazedEarC

    GrazingUnMetLeafC = AZMAX1(GrazingDmndEarC-GrazedEarC)
    GrazingDmdGrainC  = GrazedPhytoLeafC_pft*GrainStrutElms_pft(ielmc,NZ)/totShootC+GrazingUnMetLeafC
    GrazedGrainC      = AMIN1(GrainStrutElms_pft(ielmc,NZ),GrazingDmdGrainC)
    HarvestedGrainC   = GrazedGrainC
    GrazingUnMetLeafC = AZMAX1(GrazingDmdGrainC-GrazedGrainC)
  ELSE
    HarvestedPetoleC    = 0._r8
    GrazedPetoleNonstC  = 0._r8
    GrazedPetoleNoduleC = 0._r8
    HarvestedShethC     = 0._r8
    HarvestedEarC       = 0._r8
    HarvestedGrainC     = 0._r8
    GrazingUnMetLeafC   = GrazingUnMetLeafC+GrazedPhytoLeafC_pft
  ENDIF

  GrazedCanopyNonstC    = GrazedLeafNonstC+GrazedPetoleNonstC
  GrazedCanopyNoduleC   = GrazedLeafNoduleC+GrazedPetoleNoduleC
  GrazedPhytoStalkC_pft = TotPhytomassRemoved*FracBiomHarvsted(iHarvst_pft,iplthvst_stalk,NZ)
  !
  !     STALK GRAZED, REMOVED
  !
  !     WTSTK,WTRSV=stalk,reserve C mass
  !     WHVST*,WHVRV*=stalk,reserve C removed
  !     GrazingUnMetLeafC=grazing requirement unmet by stalk,reserve
  !
  TotStalkC=StalkStrutElms_pft(ielmc,NZ)+StalkRsrvElms_pft(ielmc,NZ)

  IF(TotStalkC.GT.GrazedPhytoStalkC_pft+GrazingUnMetLeafC)THEN
    GrazingDmndStalkC = GrazedPhytoStalkC_pft*StalkStrutElms_pft(ielmc,NZ)/TotStalkC+GrazingUnMetLeafC
    GrazedDmndStalkC  = AMIN1(StalkStrutElms_pft(ielmc,NZ),GrazingDmndStalkC)
    HarvestedStalkC      = GrazedDmndStalkC

    GrazingUnMetLeafC     = AZMAX1(GrazingDmndStalkC-GrazedDmndStalkC)
    GrazingDmndStalkRsrvC = GrazedPhytoStalkC_pft*StalkRsrvElms_pft(ielmc,NZ)/TotStalkC+GrazingUnMetLeafC
    GrazedStalkRsrvC      = AMIN1(StalkRsrvElms_pft(ielmc,NZ),GrazingDmndStalkRsrvC)
    HarvestedStalkRsrvC      = GrazedStalkRsrvC
    GrazingUnMetLeafC     = AZMAX1(GrazingDmndStalkRsrvC-GrazedStalkRsrvC)
  ELSE
    HarvestedStalkC      = 0._r8
    HarvestedStalkRsrvC  = 0._r8
    GrazingUnMetLeafC = AZMAX1(GrazedPhytoStalkC_pft)
    !
    !     ALLOCATE UNMET DEMAND FOR GRAZING TO LEAF,PetolSheth,HUSK
    !     EAR,GRAIN
    !
    !     WHVSL*,WHVSC*,WHVSN=leaf,nonstructural,bacteria removed
    !     WHVSH*,WHVHS,WHVEA,WHVGR,WHVSC=
    !            PetolSheth,husk,ear,grain,nonstructural C removed
    !
    IF(GrazingUnMetLeafC.GT.0.0_r8)THEN
      GrazingMeetLeafC      = AMIN1(LeafStrutElms_pft(ielmc,NZ)-HarvestedLeafC-GrazedLeafNonstC,GrazingUnMetLeafC)
      HarvestedLeafC = HarvestedLeafC+GrazingMeetLeafC*(1._r8-CCPOLX)

      GrazedLeafNonstC = GrazedLeafNonstC+GrazingMeetLeafC*CCPOLX
      GrazedLeafNoduleC = GrazedLeafNoduleC+GrazingMeetLeafC*CCPLNX

      GrazingUnMetLeafC = AZMAX1(GrazingUnMetLeafC-GrazingMeetLeafC)

      IF(totShootC.GT.ZERO4Groth_pft(NZ))THEN
        GrazingDmndPetole = GrazingUnMetLeafC*PetolShethStrutElms_pft(ielmc,NZ)/totShootC
        GrazingMeetPetoleC     = AMIN1(PetolShethStrutElms_pft(ielmc,NZ),GrazingDmndPetole)
        HarvestedPetoleC            = HarvestedPetoleC+GrazingMeetPetoleC*(1._r8-CCPOLX)
        GrazedPetoleNonstC            = GrazedPetoleNonstC+GrazingMeetPetoleC*CCPOLX
        GrazedPetoleNoduleC            = GrazedPetoleNoduleC+GrazingMeetPetoleC*CCPLNX

        GrazingUnMetLeafC = AZMAX1(GrazingUnMetLeafC-GrazingMeetPetoleC)
        GrazingDmndHuskC  = GrazingUnMetLeafC*HuskStrutElms_pft(ielmc,NZ)/totShootC
        GrazedHuskC       = AMIN1(HuskStrutElms_pft(ielmc,NZ),GrazingDmndHuskC)
        HarvestedShethC      = HarvestedShethC+GrazedHuskC

        GrazingUnMetLeafC = AZMAX1(GrazingUnMetLeafC-GrazedHuskC)
        GrazingDmndEarC   = GrazingUnMetLeafC*EarStrutElms_pft(ielmc,NZ)/totShootC
        GrazedEarC        = AMIN1(EarStrutElms_pft(ielmc,NZ),GrazingDmndEarC)
        HarvestedEarC        = HarvestedEarC+GrazedEarC

        GrazingUnMetLeafC = AZMAX1(GrazingDmndEarC-GrazedEarC)
        GrazingDmdGrainC  = GrazingUnMetLeafC*GrainStrutElms_pft(ielmc,NZ)/totShootC
        GrazedGrainC      = AMIN1(GrainStrutElms_pft(ielmc,NZ),GrazingDmdGrainC)
        HarvestedGrainC      = HarvestedGrainC+GrazedGrainC
        
        GrazingUnMetLeafC = AZMAX1(GrazingDmdGrainC-GrazedGrainC)
      ENDIF
    ENDIF
  ENDIF
  !     ALL HARVEST REMOVALS
  !
  !     LeafLayerC_brch=branch leaf C mass in canopy layer
  !
  D9860: DO NB=1,NumOfBranches_pft(NZ)
    DO  L=1,NumCanopyLayers1
      DO  K=0,MaxNodesPerBranch1
        LeafLayerC_brch(L,NB,NZ)=0._r8
      enddo
    enddo
  ENDDO D9860

  D9870: DO NB=1,NumOfBranches_pft(NZ)
    DO  L=1,NumCanopyLayers1
      DO  K=0,MaxNodesPerBranch1
        LeafLayerC_brch(L,NB,NZ)=LeafLayerC_brch(L,NB,NZ)+LeafLayerElms_node(ielmc,L,K,NB,NZ)
      enddo
    enddo
  ENDDO D9870
  end associate
  end subroutine GrazingPlant

!--------------------------------------------------------------------------------

  subroutine CutBranchNonstalByGrazing(I,J,NB,NZ,GrazedCanopyNonstC,CanopyNonstElmCopy_brch,&
    GrazedCanopyNoduleC,CanopyNodulNonstElmCopy_brch,CanopyNonstElmAfhvst_brch,CanopyNodulNonstElmAfhvst,&
    CanopyNodulStrutElmAfhvst,FracLeafSheathLeft)
  implicit none
  integer, intent(in)   :: I,J
  integer, intent(in)   :: NB,NZ
  real(r8), intent(in)  :: GrazedCanopyNonstC
  real(r8), intent(in)  :: CanopyNonstElmCopy_brch(NumPlantChemElms)
  real(r8), intent(in)  :: GrazedCanopyNoduleC
  real(r8), intent(in)  :: FracLeafSheathLeft !fraction of branch leaf/sheath C remaining after grazing
  real(r8), intent(in)  :: CanopyNodulNonstElmCopy_brch(NumPlantChemElms)
  REAL(R8), intent(out) :: CanopyNonstElmAfhvst_brch(NumPlantChemElms)        !canopy nonstructural element after harvest
  REAL(R8), intent(out) :: CanopyNodulNonstElmAfhvst(NumPlantChemElms)   !canopy nodule nonstrucal element after harvest
  real(r8), intent(out) :: CanopyNodulStrutElmAfhvst(NumPlantChemElms)  !canopy nodule structural element after harvest
  real(r8) :: HarvestedCanopyNoduleC
  real(r8) :: HarvestedCanopyNonstC_brch
  real(r8) :: FracSheath_brch,fracHavested  
  integer  :: NE
  associate(                                                         &
    CanopyLeafSheathC_pft     => plt_biom%CanopyLeafSheathC_pft,     &
    CanopyNodulStrutElms_brch => plt_biom%CanopyNodulStrutElms_brch, &
    CanopyLeafSheathC_brch    => plt_biom%CanopyLeafSheathC_brch,    &
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch,      &
    CanopyNodulNonstElms_brch => plt_biom%CanopyNodulNonstElms_brch, &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,            &
    ZERO4LeafVar_pft          => plt_biom%ZERO4LeafVar_pft           &
  )
  IF(CanopyLeafSheathC_pft(NZ).GT.ZERO4LeafVar_pft(NZ))THEN
    FracSheath_brch=AZMAX1(CanopyLeafSheathC_brch(NB,NZ))/CanopyLeafSheathC_pft(NZ)
    !
    !there is nonstructural carbon pool to remove
    IF(CanopyNonstElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      HarvestedCanopyNonstC_brch  = AZMAX1(GrazedCanopyNonstC)*FracSheath_brch
      CanopyNonstElmAfhvst_brch(ielmc) = AZMAX1(CanopyNonstElmCopy_brch(ielmc)-HarvestedCanopyNonstC_brch)
      fracHavested=AMIN1(HarvestedCanopyNonstC_brch/CanopyNonstElmCopy_brch(ielmc),1._R8)
      DO NE=2,NumPlantChemElms
        CanopyNonstElmAfhvst_brch(NE)=CanopyNonstElmCopy_brch(NE)*(1._r8-fracHavested)
      ENDDO
    ELSE
      !Depleted C does not imply that the remaining N/P reserves were grazed.
      CanopyNonstElmAfhvst_brch(:)=CanopyNonstElmCopy_brch(:)*FracLeafSheathLeft
    ENDIF

    IF(CanopyNodulNonstElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      HarvestedCanopyNoduleC=AZMAX1(GrazedCanopyNoduleC)*FracSheath_brch
      CanopyNodulNonstElmAfhvst(ielmc)=AZMAX1(CanopyNodulNonstElmCopy_brch(ielmc)-HarvestedCanopyNoduleC)
      fracHavested=AMIN1(HarvestedCanopyNoduleC/CanopyNodulNonstElmCopy_brch(ielmc),1._R8)
      DO NE=2,NumPlantChemElms
        CanopyNodulNonstElmAfhvst(NE)=CanopyNodulNonstElmCopy_brch(NE)*(1._r8-fracHavested)
      ENDDO
      
      DO NE=1,NumPlantChemElms
        CanopyNodulStrutElmAfhvst(NE)=CanopyNodulStrutElms_brch(NE,NB,NZ)*(1._r8-fracHavested)
      ENDDO
    ELSE
      !Use host tissue removal when nodule C cannot define a removal fraction.
      CanopyNodulNonstElmAfhvst(:)=CanopyNodulNonstElmCopy_brch(:)*FracLeafSheathLeft
      CanopyNodulStrutElmAfhvst(:)=CanopyNodulStrutElms_brch(:,NB,NZ)*FracLeafSheathLeft
    ENDIF
  ELSE
    !A missing canopy allocation denominator must not manufacture removal.
    CanopyNonstElmAfhvst_brch(:)=CanopyNonstElmCopy_brch(:)*FracLeafSheathLeft
    CanopyNodulNonstElmAfhvst(:)=CanopyNodulNonstElmCopy_brch(:)*FracLeafSheathLeft
    CanopyNodulStrutElmAfhvst(:)=CanopyNodulStrutElms_brch(:,NB,NZ)*FracLeafSheathLeft
  ENDIF
  end associate
  end subroutine CutBranchNonstalByGrazing
  ![tail]
end module PlantDisturbByGrazingMod