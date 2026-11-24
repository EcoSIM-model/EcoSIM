module PlantDisturbByFireMod
  use data_kind_mod, only : r8 => DAT_KIND_R8,yearIJ_type
  use ElmIDMod
  use PlantAPICommonData
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantSoilChemistryAPIData, only : plt_soilchem
  use PlantAllometryAPIData, only : plt_allom
  use PlantBiomassAPIData, only : plt_biom
  use PlantDisturbanceAPIData, only : plt_distb
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantRootBGCAPIData, only : plt_rbgc
  use PlantBGCPars
  use EcoSimConst
  use PlantMathFuncMod
  implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  public :: StageRootRemovalByFire
  public :: AbvGrndLiterFallByFire
  public :: AbvgBiomRemovalByFire
  public :: ShootBiomRemovalByFire
  public :: InitPlantFireMod
  public :: RootRemovalLbyFire
  real(r8) :: EFIRE(2,5:5) !emission efficiency of N and P along with C loss due to fire

contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine InitPlantFireMod
  implicit none
  EFIRE=reshape(real((/0.917,0.167/),r8),shape(EFIRE))
  end subroutine InitPlantFireMod


!----------------------------------------------------------------------------------------------------
  subroutine StageRootRemovalByFire(yearIJ,NZ,L,FFIRE,DCORP,FracLeftThin)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  integer, intent(in) :: L
  real(r8), intent(in):: DCORP  !disturbance rate
  real(r8),intent(out) :: FFIRE(NumPlantChemElms) !fraction of C/N/P removed by fire combustion
  real(r8), intent(out):: FracLeftThin            !fraction left after fire
  character(len=*), parameter :: subname='StageRootRemovalByFire'

  associate(                                                 &
    iSoilDisturbType_col => plt_distb%iSoilDisturbType_col  ,& !input  :soil disturbance type, [-]
    CSoilOrgM_vr         => plt_soilchem%CSoilOrgM_vr       ,& !input  :soil organic C content, [gC kg soil-1]
    FracBiomHarvsted     => plt_distb%FracBiomHarvsted      ,& !input  :harvest efficiency, [-]
    iHarvstType_pft      => plt_distb%iHarvstType_pft       ,& !input  :type of harvest,[-]
    THETW_vr             => plt_soilchem%THETW_vr            & !input  :volumetric water content, [m3 m-3]
  )
  call PrintInfo('beg '//subname)
  
  !soil too wet, or insufficient soil organic matter or not  or not fire disturbance 
  IF(THETW_vr(L).GT.VolMaxSoilMoist4Fire .OR. CSoilOrgM_vr(ielmc,L).LE.FORGC .OR. iSoilDisturbType_col.NE.itill_fire)THEN
    FracLeftThin              = 1.0_r8
    FFIRE(1:NumPlantChemElms) = 0._r8
  ELSE
    FracLeftThin=1.0_r8-DCORP*FracBiomHarvsted(iHarvst_pft,iplthvst_stalk,NZ) &
      *AMIN1(1.0_r8,(CSoilOrgM_vr(ielmc,L)-FORGC)/(orgcden-FORGC))

    FFIRE(ielmc) = FracBiomHarvsted(iHarvst_col,iplthvst_stalk,NZ)
    FFIRE(ielmn) = FFIRE(ielmc)*EFIRE(1,iHarvstType_pft(NZ)) !N
    FFIRE(ielmp) = FFIRE(ielmc)*EFIRE(2,iHarvstType_pft(NZ)) !P
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine StageRootRemovalByFire

!----------------------------------------------------------------------------------------------------
  subroutine RemoveRootByFire(yearIJ,NZ,HarvestedBiomass,FFIRE)
  !
  !upon combustion:
  !plant C is turned into CO2, or CH4
  !plant N is turned into NH3, or N2O (which is set to zero)
  !plant P is turned in PO4.
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), INTENT(IN) :: HarvestedBiomass(NumPlantChemElms)
  real(r8), intent(in) :: FFIRE(NumPlantChemElms)

  character(len=*), parameter :: subname='RemoveRootByFire'
  real(r8) :: dFireE(NumPlantChemElms),dCO2,dCH4,dNH3
  integer :: NE

  associate(                                                     &
    CO2NetFix_pft        => plt_bgcr%CO2NetFix_pft,              & !inoput :canopy net CO2 exchange,      [gC d-2 h-1]
    PO4byFire_CumYr_pft  => plt_distb%PO4byFire_CumYr_pft,       & !inoput :plant PO4 emission from fire, [g d-2 ]
    N2ObyFire_CumYr_pft  => plt_distb%N2ObyFire_CumYr_pft,       & !inoput :plant N2O emission from fire, [g d-2 ]
    NH3byFire_CumYr_pft  => plt_distb%NH3byFire_CumYr_pft,       & !inoput :plant NH3 emission from fire, [g d-2 ]
    O2ByFire_CumYr_pft   => plt_distb%O2ByFire_CumYr_pft,        & !inoput :plant O2 uptake from fire,    [g d-2 ]
    CH4ByFire_CumYr_pft  => plt_distb%CH4ByFire_CumYr_pft,       & !inoput :plant CH4 emission from fire, [g d-2 ]
    CO2ByFire_CumYr_pft  => plt_distb%CO2ByFire_CumYr_pft,       & !inoput :plant CO2 emission from fire, [g d-2 ]
    PlantElmDistLoss_pft => plt_distb%PlantElmDistLoss_pft      ,& !inoput :plant loss to disturbance,    [g d-2 h-1]       
    RootLost2Fire_pft    => plt_distb%RootLost2Fire_pft         ,& !inoput :plant root biomass lost by fire, [g d-2 h-1]
    FireLossE_pft        => plt_distb%FireLossE_pft      ,       & !inoput :plant element lost by fire, [g d-2 h-1]
    Eco_NBP_CumYr_col    => plt_bgcr%Eco_NBP_CumYr_col           & !inoput :total NBP, >0 into ecosystem,  [g d-2]
  )
  call PrintInfo('beg '//subname)
  if(FFIRE(ielmc).LE.0._r8)return
  DO NE=1,NumPlantChemElms
    dFireE(NE)    = FFIRE(NE)*HarvestedBiomass(NE)
  ENDDO

  dCO2                    = (1._r8-FrcAsCH4byFire)*dFireE(ielmc)
  dCH4                    = FrcAsCH4byFire*dFireE(ielmc)
  dNH3                    = dFireE(ielmn)
  CO2ByFire_CumYr_pft(NZ) = CO2ByFire_CumYr_pft(NZ)-dCO2
  CH4ByFire_CumYr_pft(NZ) = CH4ByFire_CumYr_pft(NZ)-dCH4
  O2ByFire_CumYr_pft(NZ)  = O2ByFire_CumYr_pft(NZ)-dCO2*2.667_r8
  NH3byFire_CumYr_pft(NZ) = NH3byFire_CumYr_pft(NZ)-dNH3
  N2ObyFire_CumYr_pft(NZ) = N2ObyFire_CumYr_pft(NZ)-0.0_r8
  PO4byFire_CumYr_pft(NZ) = PO4byFire_CumYr_pft(NZ)-dFireE(ielmp)
  CO2NetFix_pft(NZ)       = CO2NetFix_pft(NZ)-dCO2
  !
  FireLossE_pft(:,NZ)        = FireLossE_pft(:,NZ)+dFireE(:)
  PlantElmDistLoss_pft(:,NZ) = PlantElmDistLoss_pft(:,NZ)+dFireE(:)
  RootLost2Fire_pft(:,NZ)    = RootLost2Fire_pft(:,NZ) +dFireE(:)
  !
  !should CO2 be removed here?
  Eco_NBP_CumYr_col       = Eco_NBP_CumYr_col-dCH4
  call PrintInfo('end '//subname)
  end associate
  end subroutine RemoveRootByFire

!----------------------------------------------------------------------------------------------------
  subroutine AbvGrndLiterFallByFire(yearIJ,NZ,CanopyNonstElm2Litr,StandeadElmntOffEcosystem, &
    FineNonleafElmOffEcosystem,LeafElmnt2Litr,LeafElmntOffEcosystem,NonstructElmntOffEcosystem,&
    WoodyElmntOffEcosystem,WoodyElmnt2Litr,StandeadElmnt2Litr,PetolShethElmntHarv2Litr,&
    FineNonleafElmnt2Litr,LeafElmntHarv2Litr,StandeadElmntHarv2Litr,WoodyElmntHarv2Litr)
  use PlantBalMod,        only: SumPlantBiome    
  implicit none 
  type(yearIJ_type),intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), intent(in) :: LeafElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: NonstructElmntOffEcosystem(NumPlantChemElms)  
  real(r8), intent(in) :: CanopyNonstElm2Litr(NumPlantChemElms)  
  real(r8), intent(in) :: StandeadElmntOffEcosystem(NumPlantChemElms)  
  real(r8), intent(in) :: FineNonleafElmOffEcosystem(NumPlantChemElms)  
  real(r8), intent(in) :: LeafElmntOffEcosystem(NumPlantChemElms)  
  real(r8), intent(in) :: WoodyElmntOffEcosystem(NumPlantChemElms)  
  real(r8), intent(in) :: WoodyElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: StandeadElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: PetolShethElmntHarv2Litr(NumPlantChemElms)
  real(r8), intent(in) :: FineNonleafElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: LeafElmntHarv2Litr(NumPlantChemElms)
  real(r8), intent(in) :: StandeadElmntHarv2Litr(NumPlantChemElms)
  real(r8), intent(in) :: WoodyElmntHarv2Litr(NumPlantChemElms)
  character(len=*), parameter :: subname='AbvGrndLiterFallByFire'

  integer  :: M,NE
  real(r8) :: dDeadE
  real(r8) :: FineLiterElms(NumPlantChemElms)
  real(r8) :: Stddead2LitrElms(NumPlantChemElms)
  real(r8) :: Wood2LitrElms(NumPlantChemElms)
  real(r8) :: tvegE(NumPlantChemElms)
  associate(                                                             &
    PlantElmAllocMat4Litr      => plt_soilchem%PlantElmAllocMat4Litr    ,& !input  :litter kinetic fraction, [-]
    icwood                     => pltpar%icwood                         ,& !input  :group id of coarse woody litter
    istalk                     => pltpar%istalk                         ,& !input  :group id of plant stalk litter group
    ilignin                    => pltpar%ilignin                        ,& !input  :kinetic id of litter component as lignin
    inonstruct                 => pltpar%inonstruct                     ,& !input  :group id of plant nonstructural litter
    k_fine_comp                => pltpar%k_fine_comp                    ,& !input  :fine litter complex id
    k_woody_comp               => pltpar%k_woody_comp                   ,& !input  :woody litter complex id
    ifoliar                    => pltpar%ifoliar                        ,& !input  :group id of plant foliar litter
    inonfoliar                 => pltpar%inonfoliar                     ,& !input  :group id of plant non-foliar litter group
    FracWoodStalkElmAlloc2Litr => plt_allom%FracWoodStalkElmAlloc2Litr  ,& !input  :woody element allocation,[-]
    iPlant2ndGrothPattern_pft  => plt_pheno%iPlant2ndGrothPattern_pft   ,& !input  :plant expression of secondary growth, [-]                
    iPlantTurnoverPattern_pft  => plt_pheno%iPlantTurnoverPattern_pft   ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    iPlantRootProfile_pft      => plt_pheno%iPlantRootProfile_pft       ,& !input  :plant growth type (vascular, non-vascular),[-]
    LitrfallElms_pvr           => plt_bgcr%LitrfallElms_pvr             ,& !inoput :plant LitrFall element, [g d-2 h-1]
    StandDeadCompKElms_pft     => plt_biom%StandDeadCompKElms_pft        & !inoput :standing dead element fraction, [g d-2]
  )

  call PrintInfo('beg '//subname)
  FineLiterElms    = FineNonleafElmnt2Litr(:)+PetolShethElmntHarv2Litr(:)
  Stddead2LitrElms = StandeadElmnt2Litr(:)+StandeadElmntHarv2Litr(:)
  Wood2LitrElms    = WoodyElmnt2Litr(:)+WoodyElmntHarv2Litr(:)

  D6485: DO M=1,jsken
    LitrfallElms_pvr(ielmc,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(ielmc,M,k_fine_comp,0,NZ) &
      +PlantElmAllocMat4Litr(ielmc,inonstruct,M,NZ)*CanopyNonstElm2Litr(ielmc) &
      +PlantElmAllocMat4Litr(ielmc,ifoliar,M,NZ)   *(LeafElmnt2Litr(ielmc)+LeafElmntHarv2Litr(ielmc)) &
      +PlantElmAllocMat4Litr(ielmc,inonfoliar,M,NZ)*FineLiterElms(ielmc)
    
    DO NE=2,NumPlantChemElms
      LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) &
        +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*NonstructElmntOffEcosystem(NE) &
        +PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)   *LeafElmntOffEcosystem(NE) &
        +PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*FineNonleafElmOffEcosystem(NE)

      !pyrogenic material 
      LitrfallElms_pvr(NE,ilignin,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,ilignin,k_fine_comp,0,NZ) &
        +PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*(CanopyNonstElm2Litr(NE)-NonstructElmntOffEcosystem(NE)) &
        +PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)   *(LeafElmnt2Litr(NE)+LeafElmntHarv2Litr(NE)-LeafElmntOffEcosystem(NE)) &
        +PlantElmAllocMat4Litr(NE,inonfoliar,M,NZ)*(FineLiterElms(NE) -FineNonleafElmOffEcosystem(NE))
    ENDDO

    !all above ground is subject to fire
    IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR.  &      !all aboveground turnedover
      (.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ),iPlant2ndGrothPattern_pft(NZ))))THEN !herbaceous
      
      LitrfallElms_pvr(ielmc,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(ielmc,M,k_fine_comp,0,NZ)+&
        PlantElmAllocMat4Litr(ielmc,istalk,M,NZ)*(Wood2LitrElms(ielmc)+Stddead2LitrElms(ielmc))

      DO NE=2,NumPlantChemElms
        LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)+&
          PlantElmAllocMat4Litr(NE,istalk,M,NZ)*(WoodyElmntOffEcosystem(NE)+StandeadElmntOffEcosystem(NE))

        !pyrogenic material 
        LitrfallElms_pvr(NE,ilignin,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,ilignin,k_fine_comp,0,NZ) &
          +PlantElmAllocMat4Litr(NE,istalk,M,NZ)*(Wood2LitrElms(NE)-WoodyElmntOffEcosystem(NE) &
          +Stddead2LitrElms(NE)-StandeadElmntOffEcosystem(NE))
      ENDDO
      !not grass, has woody component 
    ELSE

      !add wood to standing dead
      StandDeadCompKElms_pft(ielmc,M,NZ)=StandDeadCompKElms_pft(ielmc,M,NZ)+ &
        PlantElmAllocMat4Litr(ielmc,icwood,M,NZ)*Wood2LitrElms(ielmc)

      DO NE=2,NumPlantChemElms  
        dDeadE                          = PlantElmAllocMat4Litr(NE,icwood,M,NZ)*WoodyElmntOffEcosystem(NE)
        StandDeadCompKElms_pft(NE,M,NZ) = StandDeadCompKElms_pft(NE,M,NZ)+dDeadE
      ENDDO
        
      LitrfallElms_pvr(ielmc,M,k_woody_comp,0,NZ)=LitrfallElms_pvr(ielmc,M,k_woody_comp,0,NZ) &
        +PlantElmAllocMat4Litr(ielmc,istalk,M,NZ)*FracWoodStalkElmAlloc2Litr(ielmc,k_woody_comp)&
        *Stddead2LitrElms(ielmc) 
      
      DO NE=2,NumPlantChemElms  
        LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,0,NZ) &
          +PlantElmAllocMat4Litr(NE,istalk,M,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp) &
          *StandeadElmntOffEcosystem(NE)           
        !
        !pyrogenic material 
        LitrfallElms_pvr(NE,ilignin,k_woody_comp,0,NZ)=LitrfallElms_pvr(NE,ilignin,k_woody_comp,0,NZ) &
          +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_woody_comp) &
          *(Wood2LitrElms(NE)-WoodyElmntOffEcosystem(NE)+Stddead2LitrElms(NE)-StandeadElmntOffEcosystem(NE))
          
      ENDDO    

      LitrfallElms_pvr(ielmc,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(ielmc,M,k_fine_comp,0,NZ) &
        +PlantElmAllocMat4Litr(ielmc,istalk,M,NZ)*FracWoodStalkElmAlloc2Litr(ielmc,k_fine_comp) &
        *Stddead2LitrElms(ielmc)

      DO NE=2,NumPlantChemElms    
        LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) &
          +PlantElmAllocMat4Litr(NE,istalk,M,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp) &
          *StandeadElmntOffEcosystem(NE)

        !pyrogenic material         
        LitrfallElms_pvr(NE,ilignin,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,ilignin,k_fine_comp,0,NZ) &
          +PlantElmAllocMat4Litr(NE,icwood,M,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_fine_comp) &
          *(Wood2LitrElms(NE)-WoodyElmntOffEcosystem(NE)+Stddead2LitrElms(NE)-StandeadElmntOffEcosystem(NE)) 
          
      ENDDO  

    ENDIF
  ENDDO D6485

  call PrintInfo('end '//subname)
  end associate
  end subroutine AbvGrndLiterFallByFire

!----------------------------------------------------------------------------------------------------
  subroutine AbvgBiomRemovalByFire(yearIJ,NZ,TotalElmnt2Litr,TotalElmntRemoval)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), intent(in) :: TotalElmnt2Litr(NumPlantChemElms)
  real(r8), intent(in) :: TotalElmntRemoval(NumPlantChemElms)

  character(len=*), parameter :: subname='AbvgBiomRemovalByFire'
  real(r8) :: dFireE(NumPlantChemElms)
  integer :: NE

  associate(                                                      &
    Eco_NBP_CumYr_col        => plt_bgcr%Eco_NBP_CumYr_col       ,& !inoput :total NBP, [g d-2]
    CO2NetFix_pft            => plt_bgcr%CO2NetFix_pft           ,& !inoput :canopy net CO2 exchange, [gC d-2 h-1]
    FireLossE_pft            => plt_distb%FireLossE_pft          ,& !inoput :plant element lost by fire, [g d-2 h-1]    
    NH3byFire_CumYr_pft      => plt_distb%NH3byFire_CumYr_pft    ,& !inoput :plant NH3 emission from fire, [g d-2 ]
    PO4byFire_CumYr_pft      => plt_distb%PO4byFire_CumYr_pft    ,& !inoput :plant PO4 emission from fire, [g d-2 ]
    CH4ByFire_CumYr_pft      => plt_distb%CH4ByFire_CumYr_pft    ,& !inoput :plant CH4 emission from fire, [g d-2 ]
    O2ByFire_CumYr_pft       => plt_distb%O2ByFire_CumYr_pft     ,& !inoput :plant O2 uptake from fire, [g d-2 ]
    PlantElmDistLoss_pft     => plt_distb%PlantElmDistLoss_pft   ,& !inoput :plant loss to disturbance,    [g d-2 h-1]           
    N2ObyFire_CumYr_pft      => plt_distb%N2ObyFire_CumYr_pft    ,& !inoput :plant N2O emission from fire, [g d-2 ]
    CO2ByFire_CumYr_pft      => plt_distb%CO2ByFire_CumYr_pft     & !inoput :plant CO2 emission from fire, [g d-2 ]
  )
  call PrintInfo('beg '//subname)
  DO NE=1,NumPlantChemElms
    dFireE(NE)=TotalElmntRemoval(NE)-TotalElmnt2Litr(NE)
  ENDDO

  CO2ByFire_CumYr_pft(NZ) = CO2ByFire_CumYr_pft(NZ)-(1._r8-FrcAsCH4byFire)*dFireE(ielmc)
  CH4ByFire_CumYr_pft(NZ) = CH4ByFire_CumYr_pft(NZ)-FrcAsCH4byFire*dFireE(ielmc)
  !
  O2ByFire_CumYr_pft(NZ)  = O2ByFire_CumYr_pft(NZ)-(1._r8-FrcAsCH4byFire)*dFireE(ielmc)*2.667_r8
  !
  NH3byFire_CumYr_pft(NZ) = NH3byFire_CumYr_pft(NZ)-dFireE(ielmn)
  N2ObyFire_CumYr_pft(NZ) = N2ObyFire_CumYr_pft(NZ)-0.0_r8
  PO4byFire_CumYr_pft(NZ) = PO4byFire_CumYr_pft(NZ)-dFireE(ielmp)
  !
  CO2NetFix_pft(NZ)       = CO2NetFix_pft(NZ)-(1._r8-FrcAsCH4byFire)*dFireE(ielmc)
  Eco_NBP_CumYr_col       = Eco_NBP_CumYr_col-FrcAsCH4byFire*dFireE(ielmc)

  FireLossE_pft(:,NZ)        = FireLossE_pft(:,NZ)+dFireE(:)
  PlantElmDistLoss_pft(:,NZ) = PlantElmDistLoss_pft(:,NZ)+dFireE(:)
  call PrintInfo('end '//subname)
  end associate
  end subroutine AbvgBiomRemovalByFire

!----------------------------------------------------------------------------------------------------
  subroutine ShootBiomRemovalByFire(I,J,NZ,EHVST21,EHVST22, EHVST23, EHVST24,&
    StandeadElmntRemoval,CanopyNonstElmRemoval,LeafElmntRemoval,WoodyElmntRemoval,&
    FineNonleafElmntRemoval,CanopyNonstElm2Litr,NonstructElmntOffEcosystem,&
    LeafElmntOffEcosystem,FineNonleafElmOffEcosystem,WoodyElmntOffEcosystem,&
    StandeadElmntOffEcosystem,LeafElmnt2Litr,FineNonleafElmnt2Litr,&
    WoodyElmnt2Litr,StandeadElmnt2Litr)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ
  real(r8), intent(in) :: EHVST21      !fraction of fire-removed leaf left as litter
  real(r8), intent(in) :: EHVST22      !fraction of fire-removed fine non-leaf left as litter
  real(r8), intent(in) :: EHVST23      !fraction of fire-removed stalk left as litter
  real(r8), intent(in) :: EHVST24      !fraction of fire-removed standing dead left as litter
  real(r8), intent(in) :: StandeadElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: CanopyNonstElmRemoval(NumPlantChemElms)
  real(r8), intent(in) :: LeafElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: WoodyElmntRemoval(NumPlantChemElms)
  real(r8), intent(in) :: FineNonleafElmntRemoval(NumPlantChemElms)  
  real(r8), intent(out) :: CanopyNonstElm2Litr(NumPlantChemElms)
  real(r8), intent(out) :: NonstructElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: LeafElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: FineNonleafElmOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: WoodyElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: StandeadElmntOffEcosystem(NumPlantChemElms)
  real(r8), intent(out) :: LeafElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: FineNonleafElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: WoodyElmnt2Litr(NumPlantChemElms)
  real(r8), intent(out) :: StandeadElmnt2Litr(NumPlantChemElms)
  character(len=*), parameter :: subname='ShootBiomRemovalByFire'


  associate(                                         &
    iHarvstType_pft  => plt_distb%iHarvstType_pft   ,& !input  :type of harvest,[-]
    FracBiomHarvsted => plt_distb%FracBiomHarvsted   & !input  :harvest efficiency, [-]
  )
  call PrintInfo('beg '//subname)
  NonstructElmntOffEcosystem=0._r8;LeafElmntOffEcosystem=0._r8
  FineNonleafElmOffEcosystem=0._r8;WoodyElmntOffEcosystem=0._r8
  StandeadElmntOffEcosystem=0._r8

  !leaf
  CanopyNonstElm2Litr(ielmc) = CanopyNonstElmRemoval(ielmc)*EHVST21
  CanopyNonstElm2Litr(ielmn) = CanopyNonstElmRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ))
  CanopyNonstElm2Litr(ielmp) = CanopyNonstElmRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ))

  NonstructElmntOffEcosystem(ielmn) = CanopyNonstElmRemoval(ielmn)*EHVST21
  NonstructElmntOffEcosystem(ielmp) = CanopyNonstElmRemoval(ielmp)*EHVST21

  LeafElmnt2Litr(ielmc)        = LeafElmntRemoval(ielmc)*EHVST21
  LeafElmnt2Litr(ielmn)        = LeafElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ))
  LeafElmnt2Litr(ielmp)        = LeafElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ))

  LeafElmntOffEcosystem(ielmn) = LeafElmntRemoval(ielmn)*EHVST21
  LeafElmntOffEcosystem(ielmp) = LeafElmntRemoval(ielmp)*EHVST21
  
  !non-leaf
  FineNonleafElmnt2Litr(ielmc)=FineNonleafElmntRemoval(ielmc)*EHVST22
  FineNonleafElmnt2Litr(ielmn)=FineNonleafElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ))
  FineNonleafElmnt2Litr(ielmp)=FineNonleafElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ))

  FineNonleafElmOffEcosystem(ielmn)=FineNonleafElmntRemoval(ielmn)*EHVST22
  FineNonleafElmOffEcosystem(ielmp)=FineNonleafElmntRemoval(ielmp)*EHVST22

  !stalk
  WoodyElmnt2Litr(ielmc)        = WoodyElmntRemoval(ielmc)*EHVST23
  WoodyElmnt2Litr(ielmn)        = WoodyElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_stalk,NZ))
  WoodyElmnt2Litr(ielmp)        = WoodyElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_stalk,NZ))

  WoodyElmntOffEcosystem(ielmn) = WoodyElmntRemoval(ielmn)*EHVST23
  WoodyElmntOffEcosystem(ielmp) = WoodyElmntRemoval(ielmp)*EHVST23

  !standing dead
  StandeadElmnt2Litr(ielmc) = StandeadElmntRemoval(ielmc)*EHVST24
  StandeadElmnt2Litr(ielmn) = StandeadElmntRemoval(ielmn)*(1._r8-EFIRE(1,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ))
  StandeadElmnt2Litr(ielmp) = StandeadElmntRemoval(ielmp)*(1._r8-EFIRE(2,iHarvstType_pft(NZ))*FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ))

  StandeadElmntOffEcosystem(ielmn) = StandeadElmntRemoval(ielmn)*EHVST24
  StandeadElmntOffEcosystem(ielmp) = StandeadElmntRemoval(ielmp)*EHVST24

  call PrintInfo('end '//subname)
  end associate
  end subroutine ShootBiomRemovalByFire
!----------------------------------------------------------------------------------------------------
  subroutine RootRemovalLbyFire(yearIJ,N,L,NZ,FracLeftThin,XHVST1)

  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer , intent(in) :: N,L,NZ
  real(r8), intent(out) :: FracLeftThin  !FRACTION of biomass remains alive after fire removal
  real(r8), intent(out):: XHVST1         !fraction of biomass removed by disturbance
  character(len=*), parameter :: subname='RootRemovalLbyFire'
  real(r8) :: HarvestedBiomass(NumPlantChemElms)
  real(r8) :: FFIRE(NumPlantChemElms)
  real(r8) :: XFFIRE(NumPlantChemElms)
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
    RootMediumStructElms_rpvr => plt_biom%RootMediumStructElms_rpvr ,& !input :medium root structural C/N/P, [g d-2]
    Root1stActStructElms_rpvr => plt_biom%Root1stActStructElms_rpvr  ,& !inoput :Root layer primary axes Active zone structrual element, [g d-2]    
    Root1stLigStructElms_rpvr => plt_biom%Root1stLigStructElms_rpvr  ,& !inoput :root layer lignified zone element in primary axes, [g d-2]        
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
  IF(iHarvstType_pft(NZ).NE.iharvtyp_fire)THEN
    FracLeftThin              = 1.0_r8-THIN_pft(NZ)
    FFIRE(1:NumPlantChemElms) = 0._r8
  ELSE
    !fire
    call StageRootRemovalByFire(yearIJ,NZ,L,FFIRE,DCORP,FracLeftThin)
  ENDIF
  !
  XHVST1 = 1._r8-FracLeftThin

  !fraction to litter after fire
  XFFIRE = 1._r8-FFIRE(:) 

  D3385: DO M=1,jsken
    !nonstructural root biomass
    DO NE=1,NumPlantChemElms
      HarvestedBiomass(NE)=XHVST1*PlantElmAllocMat4Litr(NE,inonstruct,M,NZ)*AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ))
      LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+XFFIRE(NE)*HarvestedBiomass(NE)
    ENDDO

    call RemoveRootByFire(yearIJ,NZ,HarvestedBiomass,FFIRE)
    !
    D3960: DO NR=1,NumStructuralRootAxes_pft(NZ)
      if(N.eq.ipltroot)THEN
        !Active primary and medium roots share the woody/fine litter partition.
        DO NE=1,NumPlantChemElms
          HarvestedBiomass(NE)=XHVST1*PlantElmAllocMat4Litr(NE,icwood,M,NZ)*(AZMAX1(Root1stActStructElms_rpvr(NE,L,NR,NZ)) &
            +AZMAX1(RootMediumStructElms_rpvr(NE,L,NR,NZ))) &
            *FracRootElmAllocm(NE,k_woody_comp)
          LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)+XFFIRE(NE)*HarvestedBiomass(NE)      
        ENDDO
        call RemoveRootByFire(yearIJ,NZ,HarvestedBiomass,FFIRE)        

        !woody roots
        DO NE=1,NumPlantChemElms
          HarvestedBiomass(NE)=XHVST1*PlantElmAllocMat4Litr(NE,iroot,M,NZ)*AZMAX1(Root1stLigStructElms_rpvr(NE,L,NR,NZ))
          LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)+XFFIRE(NE)*HarvestedBiomass(NE)
        ENDDO
        call RemoveRootByFire(yearIJ,NZ,HarvestedBiomass,FFIRE)      

        DO NE=1,NumPlantChemElms
          HarvestedBiomass(NE)=XHVST1*PlantElmAllocMat4Litr(NE,iroot,M,NZ)*(AZMAX1(Root1stActStructElms_rpvr(NE,L,NR,NZ)) &
            +AZMAX1(RootMediumStructElms_rpvr(NE,L,NR,NZ))) &
            *FracRootElmAllocm(NE,k_fine_comp)
          LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+XFFIRE(NE)*HarvestedBiomass(NE)

        ENDDO
        call RemoveRootByFire(yearIJ,NZ,HarvestedBiomass,FFIRE)            
      ENDIF

      DO NE=1,NumPlantChemElms
        HarvestedBiomass(NE)=XHVST1*PlantElmAllocMat4Litr(NE,icwood,M,NZ)*AZMAX1(RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)) &
          *FracRootElmAllocm(NE,k_fine_comp)
        LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+XFFIRE(NE)*HarvestedBiomass(NE)     
      ENDDO

      call RemoveRootByFire(yearIJ,NZ,HarvestedBiomass,FFIRE)
      !fine roots
      DO NE=1,NumPlantChemElms
        HarvestedBiomass(NE)=XHVST1*PlantElmAllocMat4Litr(NE,iroot,M,NZ)*AZMAX1(RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)) &
          *FracRootElmAllocm(NE,k_woody_comp)
        LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_woody_comp,L,NZ)+XFFIRE(NE)*HarvestedBiomass(NE)
      ENDDO

      CALL RemoveRootByFire(yearIJ,NZ,HarvestedBiomass,FFIRE)

    enddo D3960
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
  end subroutine RootRemovalLbyFire
  ![tail]
  end module PlantDisturbByFireMod