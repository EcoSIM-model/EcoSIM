module StartqMod
  use data_kind_mod,    only : r8 => DAT_KIND_R8
  use minimathmod,      only : AZMAX1, isclose
  use DebugToolMod,     only : PrintInfo
  use abortutils,       only : endrun
  use UnitMod,          only : units
  use EcoSiMParDataMod, only : pltpar
  use PlantBGCPars,     only : BlkDActCoarseRoots,BlkDLigCoarseRoots,Rax_ref
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use EcosimConst
  use TracerIDMod
  use EcoSIMCtrlDataType
  use PlantDataRateType
  use ClimForcDataType
  use PlantTraitDataType
  use PlantMgmtDataType
  use CanopyDataType
  use CanopyRadDataType
  use RootDataType
  use EcoSIMHistMod
  use GridDataType
  use EcoSIMConfig
  use PlantMathFuncMod
  implicit none

  private
  character(len=*),private, parameter :: mod_filename =&
   __FILE__
  public :: startq
  public :: InitDimensionsandUptake
  public :: InitPlantPhenoMorphoBio
  contains

  SUBROUTINE startq(NHWQ,NHEQ,NVNQ,NVSQ,NZ1Q,NZ2Q)
!
!     THIS SUBROUTINE INITIALIZES ALL PLANT VARIABLES
!
  implicit none
  integer, intent(in) :: NHWQ,NHEQ,NVNQ,NVSQ
  integer, intent(in) :: NZ1Q,NZ2Q
  character(len=*), parameter :: subname='startq'
  integer :: NY,NX,K,L,M,NZ,NZ2X
!     begin_execution
!
!     INITIALIZE SHOOT GROWTH VARIABLES
!
!     IsPlantActive_pft=PFT flag:0=not active,1=active
!     iYearPlanting_pft,iDayPlanting_pft,iYearPlantHarvest_pft,iDayPlantHarvest_pft=year,day of planting,arvesting
!     PPI,PPX=initial,current population (m-2)
!     CF,ClumpFactorInit_pft=current,initial clumping factor
!     H2OCuticleResist_pft=cuticular resistance to water (h m-1)
!     CO2CuticleResist_pft=cuticular resistance to CO2 (s m-1)
!     CNWS,rProteinC2LeafP_pft=protein:N,protein:P ratios
!    RootProteinCMax_pft=maximum root protein concentration (g g-1)
!     O2I=intercellular O2 concentration in C3,C4 PFT (umol mol-1)
!
  call PrintInfo('beg '//subname)
  D9995: DO NX=NHWQ,NHEQ
    D9990: DO NY=NVNQ,NVSQ
      NZ2X=MIN(NZ2Q,NP_col(NY,NX))
      D9985: DO NZ=NZ1Q,NZ2X
        FireReSet_pft(NZ,NY,NX)=ifalse
        IF(IsPlantActive_pft(NZ,NY,NX).EQ.iFalse)THEN

          call InitShootGrowth(NZ,NY,NX)

          call PlantLitterFractions(NZ,NY,NX)

          call PFTThermalAcclimation(NZ,NY,NX)

          call InitDimensionsandUptake(NZ,NY,NX)

          call InitPlantPhenoMorphoBio(NZ,NY,NX)

          call InitMassBalance(NZ,NY,NX)

          call InitRootMychorMorphoBio(NZ,NY,NX)

          call InitSeedMorphoBio(NZ,NY,NX)

          call InitPlantHeatandWater(NZ,NY,NX)

          if(IGO==0 .and. is_cold_run())then
            call InitStandingDead(NZ,NY,NX)
          endif
          
        ENDIF
        ZERO4Groth_pft(NZ,NY,NX)   = ZERO*PlantPopuLive_pft(NZ,NY,NX)
        ZERO4Uptk_pft(NZ,NY,NX)    = ZERO*PlantPopuLive_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        ZERO4LeafVar_pft(NZ,NY,NX) = AMIN1(ZERO*PlantPopuLive_pft(NZ,NY,NX)*1.0E+06_r8,1.e-8_r8)
      ENDDO D9985
!
!     FILL OUT UNUSED ARRAYS
!
      D9986: DO NZ=NP0_col(NY,NX)+1,JP
        SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX) = 0._r8
        LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)     = 0._r8
        StandDeadStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)         = 0._r8
        StandDeadSurfArea_pft(NZ,NY,NX)                             = 0._r8
        CanopyHeightDead_pft(NZ,NY,NX)                              = 0._R8

        D6401: DO L=1,NL_col(NY,NX)
          DO  K=1,pltpar%NumOfPlantLitrCmplxs
            DO  M=1,jskenc
              LitrfallElms_pvr(1:NumPlantChemElms,M,K,L,NZ,NY,NX)=0._r8
            enddo
          enddo
        ENDDO D6401
      ENDDO D9986
    ENDDO D9990
  ENDDO D9995
  call PrintInfo('end '//subname)
  END subroutine startq
!------------------------------------------------------------------------------------------

  SUBROUTINE InitStandingDead(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  character(len=*), parameter :: subname='InitStandingDead'
  real(r8) :: WTSTDX
  integer :: M, NE

  associate(                                   &
    icwood          => pltpar%icwood           &
  )
  CALL PrintInfo('beg '//subname)

  StandDeadStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)         = 0._r8    
  WTSTDX                                                      = StandingDeadInitC_pft(NZ,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
  D155: DO M=1,jskenc
    StandDeadCompKElms_pft(ielmc,M,NZ,NY,NX)=WTSTDX*PlantElmAllocMat4Litr(ielmc,icwood,M,NZ,NY,NX)
    StandDeadCompKElms_pft(ielmn,M,NZ,NY,NX)=WTSTDX*rNCStalk_pft(NZ,NY,NX)*PlantElmAllocMat4Litr(ielmn,icwood,M,NZ,NY,NX)
    StandDeadCompKElms_pft(ielmp,M,NZ,NY,NX)=WTSTDX*rPCStalk_pft(NZ,NY,NX)*PlantElmAllocMat4Litr(ielmp,icwood,M,NZ,NY,NX)
    
    DO NE=1,NumPlantChemElms
      StandDeadStrutElms_pft(NE,NZ,NY,NX)=StandDeadStrutElms_pft(NE,NZ,NY,NX)+StandDeadCompKElms_pft(NE,M,NZ,NY,NX)
    ENDDO  
  ENDDO D155
  StandDeadCompKElms_pft(:,jskenp1,NZ,NY,NX)             = 0._r8
  CanopyHeightDead_pft(NZ,NY,NX)                         = 0._R8
  StandDeadSurfArea_pft(NZ,NY,NX)                        = 0._r8
  CanopySurfAreaProfDead_pft(1:NumCanopyLayers,NZ,NY,NX) = 0._r8

  ENGYX_pft(NZ,NY,NX)             = 0._r8
  DeltaTKC_pft(NZ,NY,NX)          = 0._r8
  TdegCCanopy_pft(NZ,NY,NX)       = ATCA_col(NY,NX)
  TKC_pft(NZ,NY,NX)               = units%Celcius2Kelvin(TdegCCanopy_pft(NZ,NY,NX))
  TCGroth_pft(NZ,NY,NX)           = TdegCCanopy_pft(NZ,NY,NX)
  TKGroth_pft(NZ,NY,NX)           = units%Celcius2Kelvin(TCGroth_pft(NZ,NY,NX))
  fTCanopyGroth_pft(NZ,NY,NX)     = 1.0_r8
  PSICanopy_pft(NZ,NY,NX)         = -1.0E-03_r8
  PSICanopyOsmo_pft(NZ,NY,NX)     = OrganOsmoPsi0pt_pft(NZ,NY,NX)+PSICanopy_pft(NZ,NY,NX)
  PSICanopyTurg_pft(NZ,NY,NX)     = AZMAX1(PSICanopy_pft(NZ,NY,NX)-PSICanopyOsmo_pft(NZ,NY,NX))
  Transpiration_pft(NZ,NY,NX)     = 0._r8
  FracPARads2Canopy_pft(NZ,NY,NX)     = 0._r8
  FracPARads2LiveCanopy_pft(NZ,NY,NX) = 0._r8

  call PrintInfo('end '//subname)
  end associate
  end subroutine InitStandingDead
!------------------------------------------------------------------------------------------

  subroutine InitShootGrowth(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX

  iYearPlanting_pft(NZ,NY,NX)     = iPlantingYear_pft(NZ,NY,NX) !planting year
  iDayPlanting_pft(NZ,NY,NX)      = iPlantingDay_pft(NZ,NY,NX)  !planting day
  iYearPlantHarvest_pft(NZ,NY,NX) = iHarvestYear_pft(NZ,NY,NX)
  iDayPlantHarvest_pft(NZ,NY,NX)  = iHarvestDay_pft(NZ,NY,NX)
  PPI_pft(NZ,NY,NX)               = PPatSeeding_pft(NZ,NY,NX)
  PPX_pft(NZ,NY,NX)               = PPI_pft(NZ,NY,NX)
  ClumpFactor_pft(NZ,NY,NX)       = ClumpFactorInit_pft(NZ,NY,NX)       !clumping factor
  
  H2OCuticleResist_pft(NZ,NY,NX) = CuticleResist_pft(NZ,NY,NX)/3600.0_r8
  CO2CuticleResist_pft(NZ,NY,NX) = CuticleResist_pft(NZ,NY,NX)*1.56_r8
!  rProteinC2RootN_pft(NZ,NY,NX)  = 2.15_r8
!  rProteinC2LeafN_pft(NZ,NY,NX)  = 2.6_r8
!  rProteinC2LeafP_pft(NZ,NY,NX)  = 25.0_r8
  RootProteinCMax_pft(NZ,NY,NX) = rNCRoot_pft(NZ,NY,NX)*rProteinC2RootN_pft(NZ,NY,NX)    
  IF(iPlantPhotosynsType_pft(NZ,NY,NX).EQ.ic3_photo)THEN
    O2I_pft(NZ,NY,NX)=2.10E+05_r8
  ELSE
    O2I_pft(NZ,NY,NX)=3.96E+05_r8
  ENDIF
  end subroutine InitShootGrowth
!------------------------------------------------------------------------------------------

  subroutine PlantLitterFractions(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer :: N,M
  real(r8) :: CNOPC(jskenc),CPOPC(jskenc)
  REAL(R8) :: CNOPCT,CPOPCT,litter_sum
  character(len=256) :: msg

  associate(                                   &
    iprotein        => pltpar%iprotein,        &
    icarbhyro       => pltpar%icarbhyro,       &
    icellulos       => pltpar%icellulos,       &
    ilignin         => pltpar%ilignin,         &
    NumLitterGroups => pltpar%NumLitterGroups, &
    inonstruct      => pltpar%inonstruct,      &
    ifoliar         => pltpar%ifoliar,         &
    inonfoliar      => pltpar%inonfoliar,      &
    istalk          => pltpar%istalk,          &
    iroot           => pltpar%iroot,           &
    icwood          => pltpar%icwood           &
  )
!
!     FRACTIONS OF PLANT LITTER ALLOCATED TO KINETIC COMPONENTS
!     PROTEIN(*,1),CH2O(*,2),CELLULOSE(*,3),LIGNIN(*,4) IN SOIL LITTER
!
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!
!     NONSTRUCTURAL
!
  PlantElmAllocMat4Litr(ielmc,inonstruct,iprotein,NZ,NY,NX)  = 0.0_r8
  PlantElmAllocMat4Litr(ielmc,inonstruct,icarbhyro,NZ,NY,NX) = 0.67_r8
  PlantElmAllocMat4Litr(ielmc,inonstruct,icellulos,NZ,NY,NX) = 0.33_r8
  PlantElmAllocMat4Litr(ielmc,inonstruct,ilignin,NZ,NY,NX)   = 0.0_r8
  !
  !     NON-VASCULAR (E.G. MOSSES)
  !
  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ,NY,NX)))THEN
    PlantElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ,NY,NX)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ,NY,NX) = 0.25_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ,NY,NX) = 0.30_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ,NY,NX)   = 0.38_r8

    PlantElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ,NY,NX)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ,NY,NX) = 0.25_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ,NY,NX) = 0.30_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ,NY,NX)   = 0.38_r8
  !
  !     LEGUMES
  !
  ELSEIF(iPlantNfixType_pft(NZ,NY,NX).NE.0)THEN
    PlantElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ,NY,NX)  = 0.16_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ,NY,NX) = 0.38_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ,NY,NX) = 0.34_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ,NY,NX)   = 0.12_r8

    PlantElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ,NY,NX)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ,NY,NX) = 0.41_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ,NY,NX) = 0.37_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ,NY,NX)   = 0.15_r8
  !
  !     ANNUALS, GRASSES, SHRUBS
  !
  ELSEIF(iPlantTurnoverPattern_pft(NZ,NY,NX).EQ.0 .OR. &
    (.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ,NY,NX),iPlant2ndGrothPattern_pft(NZ,NY,NX))))THEN
    PlantElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ,NY,NX)  = 0.08_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ,NY,NX) = 0.41_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ,NY,NX) = 0.36_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ,NY,NX)   = 0.15_r8

    PlantElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ,NY,NX)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ,NY,NX) = 0.41_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ,NY,NX) = 0.36_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ,NY,NX)   = 0.16_r8
  !
  !     DECIDUOUS TREES
  !
  ELSEIF(iPlantTurnoverPattern_pft(NZ,NY,NX).EQ.1 .OR. iPlantTurnoverPattern_pft(NZ,NY,NX).GE.3)THEN
    PlantElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ,NY,NX)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ,NY,NX) = 0.34_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ,NY,NX) = 0.36_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ,NY,NX)   = 0.23_r8

    PlantElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ,NY,NX)  = 0.0_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ,NY,NX) = 0.045_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ,NY,NX) = 0.660_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ,NY,NX)   = 0.295_r8
  !
  !     CONIFEROUS TREES
  !
  ELSE
    PlantElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ,NY,NX)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ,NY,NX) = 0.25_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ,NY,NX) = 0.38_r8
    PlantElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ,NY,NX)   = 0.30_r8

    PlantElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ,NY,NX)  = 0.0_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ,NY,NX) = 0.045_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ,NY,NX) = 0.660_r8
    PlantElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ,NY,NX)   = 0.295_r8
  ENDIF
  !
  !     FRACTIONS OF WOODY LITTER ALLOCATED TO
  !     PROTEIN, CH2O, CELLULOSE, LIGNIN
  !
  !     NON-VASCULAR
  !
  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ,NY,NX)))THEN
    PlantElmAllocMat4Litr(ielmc,istalk,iprotein,NZ,NY,NX)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,istalk,icarbhyro,NZ,NY,NX) = 0.25_r8
    PlantElmAllocMat4Litr(ielmc,istalk,icellulos,NZ,NY,NX) = 0.30_r8
    PlantElmAllocMat4Litr(ielmc,istalk,ilignin,NZ,NY,NX)   = 0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ,NY,NX).EQ.0 .OR. &
    (.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ,NY,NX),iPlant2ndGrothPattern_pft(NZ,NY,NX))))THEN
    PlantElmAllocMat4Litr(ielmc,istalk,iprotein,NZ,NY,NX)  = 0.03_r8
    PlantElmAllocMat4Litr(ielmc,istalk,icarbhyro,NZ,NY,NX) = 0.25_r8
    PlantElmAllocMat4Litr(ielmc,istalk,icellulos,NZ,NY,NX) = 0.57_r8
    PlantElmAllocMat4Litr(ielmc,istalk,ilignin,NZ,NY,NX)   = 0.15_r8
!
!     DECIDUOUS AND CONIFEROUS TREES
!
  ELSE
    PlantElmAllocMat4Litr(ielmc,istalk,iprotein,NZ,NY,NX)  = 0.0_r8
    PlantElmAllocMat4Litr(ielmc,istalk,icarbhyro,NZ,NY,NX) = 0.045_r8
    PlantElmAllocMat4Litr(ielmc,istalk,icellulos,NZ,NY,NX) = 0.660_r8
    PlantElmAllocMat4Litr(ielmc,istalk,ilignin,NZ,NY,NX)   = 0.295_r8
  ENDIF
!
!     FRACTIONS OF FINE ROOT LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN PC&E 25:601-608
!
!     NON-VASCULAR
!
  IF(is_root_bryophyte(iPlantRootProfile_pft(NZ,NY,NX)))THEN
    PlantElmAllocMat4Litr(ielmc,iroot,iprotein,NZ,NY,NX)  = 0.07_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icarbhyro,NZ,NY,NX) = 0.25_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icellulos,NZ,NY,NX) = 0.30_r8
    PlantElmAllocMat4Litr(ielmc,iroot,ilignin,NZ,NY,NX)   = 0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ,NY,NX).EQ.0 .OR. &
    (.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ,NY,NX),iPlant2ndGrothPattern_pft(NZ,NY,NX))))THEN
    PlantElmAllocMat4Litr(ielmc,iroot,iprotein,NZ,NY,NX)  = 0.057_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icarbhyro,NZ,NY,NX) = 0.263_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icellulos,NZ,NY,NX) = 0.542_r8
    PlantElmAllocMat4Litr(ielmc,iroot,ilignin,NZ,NY,NX)   = 0.138_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ,NY,NX).EQ.1 .OR. iPlantTurnoverPattern_pft(NZ,NY,NX).GE.3)THEN
    PlantElmAllocMat4Litr(ielmc,iroot,iprotein,NZ,NY,NX)  = 0.059_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icarbhyro,NZ,NY,NX) = 0.308_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icellulos,NZ,NY,NX) = 0.464_r8
    PlantElmAllocMat4Litr(ielmc,iroot,ilignin,NZ,NY,NX)   = 0.169_r8
!
!     CONIFEROUS TREES
!
  ELSE
    PlantElmAllocMat4Litr(ielmc,iroot,iprotein,NZ,NY,NX)  = 0.059_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icarbhyro,NZ,NY,NX) = 0.308_r8
    PlantElmAllocMat4Litr(ielmc,iroot,icellulos,NZ,NY,NX) = 0.464_r8
    PlantElmAllocMat4Litr(ielmc,iroot,ilignin,NZ,NY,NX)   = 0.169_r8
  ENDIF
!
!     COARSE WOODY LITTER FROM BOLES AND ROOTS
!
  PlantElmAllocMat4Litr(ielmc,icwood,iprotein,NZ,NY,NX)  = 0.00_r8
  PlantElmAllocMat4Litr(ielmc,icwood,icarbhyro,NZ,NY,NX) = 0.045_r8
  PlantElmAllocMat4Litr(ielmc,icwood,icellulos,NZ,NY,NX) = 0.660_r8
  PlantElmAllocMat4Litr(ielmc,icwood,ilignin,NZ,NY,NX)   = 0.295_r8
!
!     Validate that carbon litter fractions across kinetic components sum to one.
!
  D95: DO N=0,NumLitterGroups
    litter_sum = SUM(PlantElmAllocMat4Litr(ielmc,N,1:jskenc,NZ,NY,NX))
    IF(.not.isclose(litter_sum,1.0_r8))THEN
      write(msg,'(A,I0,A,I0,A,I0,A,I0,A,ES16.8)') &
        'Error: PlantElmAllocMat4Litr C fractions for PFT ',NZ, &
        ', NY ',NY,', NX ',NX,', litter group ',N,' sum to ',litter_sum
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
    CNOPCT = 0.0_r8
    CPOPCT = 0.0_r8
    D100: DO M  = 1, jskenc
      CNOPCT=CNOPCT+PlantElmAllocMat4Litr(ielmc,N,M,NZ,NY,NX)*CNOPC(M)
      CPOPCT=CPOPCT+PlantElmAllocMat4Litr(ielmc,N,M,NZ,NY,NX)*CPOPC(M)
    ENDDO D100
    D105: DO M=1,jskenc
      PlantElmAllocMat4Litr(ielmn,N,M,NZ,NY,NX)=PlantElmAllocMat4Litr(ielmc,N,M,NZ,NY,NX)*CNOPC(M)/CNOPCT
      PlantElmAllocMat4Litr(ielmp,N,M,NZ,NY,NX)=PlantElmAllocMat4Litr(ielmc,N,M,NZ,NY,NX)*CPOPC(M)/CPOPCT
    ENDDO D105
  ENDDO D110
  !
  !     CONCURRENT NODE GROWTH
  !
  !     FracGroth2Node_pft=scales node number for perennial vegetation (e.g. trees)
  !     NumCogrowthNode_pft=number of concurrently growing nodes
  !
  ! deciduous (0) or shallow root (not tree)
  IF(iPlantTurnoverPattern_pft(NZ,NY,NX).EQ.0 .OR. &
    (.not.is_plant_woody_vascular(iPlantRootProfile_pft(NZ,NY,NX),iPlant2ndGrothPattern_pft(NZ,NY,NX))))THEN
    FracGroth2Node_pft(NZ,NY,NX)=1.0_r8

    IF(MatureGroup_pft(NZ,NY,NX).LE.6)THEN
      NumCogrowthNode_pft(NZ,NY,NX) = 2
    ELSEIF(MatureGroup_pft(NZ,NY,NX).LE.10)THEN
      NumCoGrowthNode_pft(NZ,NY,NX)=3
    ELSEIF(MatureGroup_pft(NZ,NY,NX).LE.15)THEN
      NumCogrowthNode_pft(NZ,NY,NX)=4
    ELSE
      NumCoGrowthNode_pft(NZ,NY,NX)=5
    ENDIF
    !non-deciduous  
  ELSE
    FracGroth2Node_pft(NZ,NY,NX) =AMAX1(1.0_r8,0.04_r8/RateRefLeafAppearance_pft(NZ,NY,NX))
    NumCoGrowthNode_pft(NZ,NY,NX)=24
  ENDIF
  end associate
  end subroutine PlantLitterFractions
!------------------------------------------------------------------------------------------

  subroutine PFTThermalAcclimation(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
  real(r8), parameter :: TCZD = 5.0_r8        !basal value for threshold temperature for spring leafout/dehardening	oC
  real(r8), parameter :: TCXD = 12.0_r8       !basal value for threshold temperature for autumn leafoff/hardening	oC

!
!     PFT THERMAL ACCLIMATION
!
!     ZTYP,PlantInitThermoAdaptZone_pft=dynamic,initial thermal adaptation zone from PFT file
!     TempOffset_pft=shift in Arrhenius curve for thermal adaptation (oC)
!     TCZ,TC4LeafOff_pft=threshold temperature for leafout,leafoff
!     HTC=high temperature threshold for grain number loss (oC)
!     SeedTempSens_pft=sensitivity to HTC (seeds oC-1 above HTC)
!
  rPlantThermoAdaptZone_pft(NZ,NY,NX) = PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  TempOffset_pft(NZ,NY,NX)            = 2.667_r8*(2.5_r8-rPlantThermoAdaptZone_pft(NZ,NY,NX))
  TC4LeafOut_pft(NZ,NY,NX)            = TCZD-TempOffset_pft(NZ,NY,NX)
  TC4LeafOff_pft(NZ,NY,NX)            = AMIN1(15.0,TCXD-TempOffset_pft(NZ,NY,NX))
  IF(iPlantPhotosynsType_pft(NZ,NY,NX).EQ.3)THEN
    IF(DATAP(NZ,NY,NX)(1:4).EQ.'soyb')THEN
      HighTempLimitSeed_pft(NZ,NY,NX) = 30.0_r8+3.0_r8*rPlantThermoAdaptZone_pft(NZ,NY,NX)
      SeedTempSens_pft(NZ,NY,NX)      = 0.002_r8
    ELSE
      HighTempLimitSeed_pft(NZ,NY,NX) = 27.0_r8+3.0_r8*rPlantThermoAdaptZone_pft(NZ,NY,NX)
      SeedTempSens_pft(NZ,NY,NX)      = 0.002_r8
    ENDIF
  ELSE
    HighTempLimitSeed_pft(NZ,NY,NX) = 27.0_r8+3.0_r8*rPlantThermoAdaptZone_pft(NZ,NY,NX)
    SeedTempSens_pft(NZ,NY,NX)      = 0.005_r8
  ENDIF
  end subroutine PFTThermalAcclimation
!------------------------------------------------------------------------------------------

  subroutine InitDimensionsandUptake(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
  INTEGER :: L,N,NR
!
!     SEED CHARACTERISTICS
!
  call calc_seed_geometry(SeedCMass_pft(NZ,NY,NX),SeedWidth2LenRatio_pft(NZ,NY,NX),SeedVolumeMean_pft(NZ,NY,NX),&
    SeedMeanLen_pft(NZ,NY,NX),SeedAreaMean_pft(NZ,NY,NX))
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
  !     RSRR,RSRA=radial,axial root resistivity (m2 MPa-1 h-1)
  !
  SeedDepth_pft(NZ,NY,NX)=PlantinDepz_pft(NZ,NY,NX)+SoilSurfDepZ_col(NY,NX)
  DO L=NU_col(NY,NX),NL_col(NY,NX)
    IF(isclose(SeedDepth_pft(NZ,NY,NX),CumSoilThickness_vr(L,NY,NX)))then
      SeedDepth_pft(NZ,NY,NX)=AMAX1(CumSoilThickness_vr(L,NY,NX)-ppmc,ppmc)
      exit
    endif
  ENDDO
  ! Preserve the revised planting depth relative to the soil surface.
  PlantinDepz_pft(NZ,NY,NX)=SeedDepth_pft(NZ,NY,NX)-SoilSurfDepZ_col(NY,NX)

  D9795: DO L=NU_col(NY,NX),NL_col(NY,NX)
    !find the seeding layer, i.e. the first root layer
    IF(SeedDepth_pft(NZ,NY,NX).GE.CumSoilThickness_vr(L-1,NY,NX) &
      .AND. SeedDepth_pft(NZ,NY,NX).LT.CumSoilThickness_vr(L,NY,NX))THEN
      NGTopRootLayer_pft(NZ,NY,NX)   = L
      NMaxRootBotLayer_pft(NZ,NY,NX) = L
      D9790: DO NR=1,pltpar%MaxNumRootAxes
        NRoot1stTipLay_raxes(NR,NZ,NY,NX)=L
      ENDDO D9790
    ENDIF
  ENDDO D9795  
  CNRTS_pft(NZ,NY,NX)                      = rNCRoot_pft(NZ,NY,NX)*RootBiomGrosYld_pft(NZ,NY,NX)
  CPRTS_pft(NZ,NY,NX)                      = rPCRootr_pft(NZ,NY,NX)*RootBiomGrosYld_pft(NZ,NY,NX)
  Root1stMaxRadius_pft(imycorrhz,NZ,NY,NX) = 8.0E-06_r8     !transport hyphae
  Root2ndMaxRadius_pft(imycorrhz,NZ,NY,NX) = 2.0E-06_r8     !absorbing hyphae
  RootPorosity_pft(imycorrhz,NZ,NY,NX)     = RootPorosity_pft(1,NZ,NY,NX)
  VmaxNH4Root_pft(imycorrhz,NZ,NY,NX)      = VmaxNH4Root_pft(1,NZ,NY,NX)
  KmNH4Root_pft(imycorrhz,NZ,NY,NX)        = KmNH4Root_pft(1,NZ,NY,NX)
  CMinNH4Root_pft(imycorrhz,NZ,NY,NX)      = CMinNH4Root_pft(1,NZ,NY,NX)
  VmaxNO3Root_pft(imycorrhz,NZ,NY,NX)      = VmaxNO3Root_pft(1,NZ,NY,NX)
  KmNO3Root_pft(imycorrhz,NZ,NY,NX)        = KmNO3Root_pft(1,NZ,NY,NX)
  CminNO3Root_pft(imycorrhz,NZ,NY,NX)      = CminNO3Root_pft(1,NZ,NY,NX)
  VmaxPO4Root_pft(imycorrhz,NZ,NY,NX)      = VmaxPO4Root_pft(1,NZ,NY,NX)
  KmPO4Root_pft(imycorrhz,NZ,NY,NX)        = KmPO4Root_pft(1,NZ,NY,NX)
  CMinPO4Root_pft(imycorrhz,NZ,NY,NX)      = CMinPO4Root_pft(1,NZ,NY,NX)
  RootRadialResist_pft(imycorrhz,NZ,NY,NX) = 1.0E+04
  Root2ndAxialResist_pft(imycorrhz,NZ,NY,NX)  = 1.0E+12 
  !
  !     RootPoreTortu4Gas_pft=tortuosity for gas transport
  !     RootRaidus_rpft=path length for radial diffusion within root (m)
  !     FineRootVolPerMassC_pft=volume:C ratio (m3 gC-1)
  !     Root1stSpecLen_pft,Root2ndSpecLen_pft=specific primary,secondary root length (m g-1)
  !     Root1stXSecArea_pft,Root2ndXSecArea_pft=specific primary,secondary root area (m2 g-1)
  !
  CRootActVolPerMassC_pft(NZ,NY,NX) = 1.e-6_r8/(BlkDActCoarseRoots*(1._r8-RootPorosity_pft(ipltroot,NZ,NY,NX)))
  CRootLigVolPerMassC_pft(NZ,NY,NX) = 1.e-6_r8/(BlkDLigCoarseRoots*(1._r8-RootPorosity_pft(ipltroot,NZ,NY,NX)))
  D500: DO N=1,2
    RootPoreTortu4Gas_pft(N,NZ,NY,NX)   = RootPorosity_pft(N,NZ,NY,NX)**1.33_r8
    RootRaidus_rpft(N,NZ,NY,NX)         = LOG(1.0_r8/SQRT(AMAX1(0.01_r8,RootPorosity_pft(N,NZ,NY,NX))))
    FineRootVolPerMassC_pft(N,NZ,NY,NX) = 1._r8/(50.e3_r8*(1.0_r8-RootPorosity_pft(N,NZ,NY,NX)))         ![50 kgC m-3 root biomass]
    Root1stSpecLen_pft(N,NZ,NY,NX)      = FineRootVolPerMassC_pft(N,NZ,NY,NX)/(PICON*Root1stMaxRadius_pft(N,NZ,NY,NX)**2)
    Root2ndSpecLen_pft(N,NZ,NY,NX)      = FineRootVolPerMassC_pft(N,NZ,NY,NX)/(PICON*Root2ndMaxRadius_pft(N,NZ,NY,NX)**2)
    Root1stMaxRadius1_pft(N,NZ,NY,NX)   = Root1stMaxRadius_pft(N,NZ,NY,NX)
!    2*SQRT(0.25*(1.0-RootPorosity_pft(N,NZ,NY,NX)))
    Root2ndMaxRadius1_pft(N,NZ,NY,NX)=Root2ndMaxRadius_pft(N,NZ,NY,NX)
!    2*SQRT(0.25*(1.0-RootPorosity_pft(N,NZ,NY,NX)))
    Root1stXSecArea_pft(N,NZ,NY,NX) = PICON*Root1stMaxRadius1_pft(N,NZ,NY,NX)**2
    Root2ndXSecArea_pft(N,NZ,NY,NX) = PICON*Root2ndMaxRadius1_pft(N,NZ,NY,NX)**2
  ENDDO D500 
  end subroutine InitDimensionsandUptake
!------------------------------------------------------------------------------------------

  subroutine InitPlantPhenoMorphoBio(NZ,NY,NX)
  !
  use GridConsts, only : MaxNumBranches,NumGrowthStages,MaxNumRootAxes
  implicit none
  integer, intent(in) :: NZ, NY, NX
  character(len=*), parameter :: subname='InitPlantPhenoMorphoBio'
  integer :: K,L,M,N,NB
  real(r8) :: NumVessels
  call PrintInfo('beg '//subname)
  !
  !     INITIALIZE PLANT PHENOLOGY
  !
  !     PP=population (grid cell-1)
  !
  CanopyLeafAreaMAX_pft(NZ,NY,NX)           = 0._R8
  PlantPopuLive_pft(NZ,NY,NX)               = PPX_pft(NZ,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
  PlantPopuDead_pft(NZ,NY,NX)               = PlantPopuLive_pft(NZ,NY,NX)
  doInitPlant_pft(NZ,NY,NX)                 = ifalse
  isPlantShootAlive_pft(NZ,NY,NX)           = iTrue
  isPlantRootAlive_pft(NZ,NY,NX)            = iTrue
  BranchNumber_pft(NZ,NY,NX)                = 0
  NumOfBranches_pft(NZ,NY,NX)               = 0
  HypocotHeight_pft(NZ,NY,NX)               = 0._r8
  CanopyHeightLive_pft(NZ,NY,NX)            = 0._r8
  Cytokinin2ndConc_rpvr(:,:,:,NZ,NY,NX)     = 0._r8
  Days4FalseBreak_pft(NZ,NY,NX)             = 0
  RootSingleVesselRstAxial_pft(NZ,NY,NX)    = AlphaVesselAxialResist_pft(NZ,NY,NX)*Rax_ref*(1.e-6_r8/RootVesselRadius_pft(NZ,NY,NX))**4
  RootSingleVesselArea_pft(NZ,NY,NX)        = PICON*RootVesselRadius_pft(NZ,NY,NX)**2
  NumVessels                                = 0.2_r8*PICON*Root2ndMaxRadius_pft(ipltroot,NZ,NY,NX)**2/RootSingleVesselArea_pft(NZ,NY,NX)
  Root2ndAxialResist_pft(ipltroot,NZ,NY,NX) = RootSingleVesselRstAxial_pft(NZ,NY,NX)/NumVessels
    
  CRootLumenArea_pvr(:,NZ,NY,NX)            = 0.2_r8*PICON*Root1stMaxRadius_pft(ipltroot,NZ,NY,NX)**2
  D10: DO NB=1,MaxNumBranches
    doInitLeafOut_brch(NB,NZ,NY,NX)                = iTrue
    EnablePlantLeafOut_brch(NB,NZ,NY,NX)           = iTrue
    doPlantLeaveOff_brch(NB,NZ,NY,NX)              = iTrue
    Prep4Literfall_brch(NB,NZ,NY,NX)               = ifalse
    Hours4LiterfalAftMature_brch(NB,NZ,NY,NX)      = 0
    MatureGroup_brch(NB,NZ,NY,NX)                  = MatureGroup_pft(NZ,NY,NX)
    ShootNodeNum_brch(NB,NZ,NY,NX)                 = ShootNodeNumAtPlanting_pft(NZ,NY,NX)
    ShootNodeNumAtInitFloral_brch(NB,NZ,NY,NX)      = ShootNodeNum_brch(NB,NZ,NY,NX)
    ShootNodeNumAtAnthesis_brch(NB,NZ,NY,NX)         = 0._r8
    NumOfLeaves_brch(NB,NZ,NY,NX)                  = 0._r8
    LeafNumberAtFloralInit_brch(NB,NZ,NY,NX)       = 0._r8
    KLeafNumber_brch(NB,NZ,NY,NX)                  = 1
    KMinNumLeaf4GroAlloc_brch(NB,NZ,NY,NX)         = 1
    KHiestGroLeafNode_brch(NB,NZ,NY,NX)            = 1
    KLowestGroLeafNode_brch(NB,NZ,NY,NX)           = 0
    NodeNumNormByMatgrp_brch(NB,NZ,NY,NX)          = 0._r8
    ReprodNodeNumNormByMatrgrp_brch(NB,NZ,NY,NX)   = 0._r8
    TotalNodeNumNormByMatgrp_brch(NB,NZ,NY,NX)     = 0._r8
    TotReproNodeNumNormByMatrgrp_brch(NB,NZ,NY,NX) = 0._r8
    dReproNodeNumNormByMatG_brch(NB,NZ,NY,NX)      = 0._r8
    Hours4LenthenPhotoPeriod_brch(NB,NZ,NY,NX)     = 0._r8
    Hours4ShortenPhotoPeriod_brch(NB,NZ,NY,NX)     = 0._r8
    Hours4Leafout_brch(NB,NZ,NY,NX)                = Hours4LenthenPhotoPeriod_brch(NB,NZ,NY,NX)
    Hours4LeafOff_brch(NB,NZ,NY,NX)                = Hours4ShortenPhotoPeriod_brch(NB,NZ,NY,NX)
    Hours2LeafOut_brch(NB,NZ,NY,NX)                = 0._r8
    RubiscoActivity_brch(NB,NZ,NY,NX)              = 1.0_r8
    GrainFillDowreg_brch(NB,NZ,NY,NX)             = 1.0_r8
    HourFailGrainFill_brch(NB,NZ,NY,NX)            = 0
    HoursDoingRemob_brch(NB,NZ,NY,NX)              = 0
    BranchNumerID_brch(NB,NZ,NY,NX)                 = 0
    isPlantBranchAlive_brch(NB,NZ,NY,NX)            = iFalse
    D15: DO M=1,NumGrowthStages
      iPlantCalendar_brch(M,NB,NZ,NY,NX)=0
    ENDDO D15
  ENDDO D10
!
!     INITIALIZE PLANT MORPHOLOGY AND BIOMASS
!
  HoursTooLowPsiCan_pft(NZ,NY,NX)                                           = 0._r8
  ChillHours_pft(NZ,NY,NX)                                                  = 0._r8
  CanopyNonstElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)        = 0._r8
  CanopyNodulNonstElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)   = 0._r8
  ShootElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)              = 0._r8
  PetolShethStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)    = 0._r8
  StalkStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)         = 0._r8
  LeafStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)          = 0._r8
  StalkRsrvElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)          = 0._r8
  HuskStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)          = 0._r8
  GrainStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)         = 0._r8
  EarStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)           = 0._r8
  CanopyNodulStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)   = 0._r8
  LeafElmntRemobFlx_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)      = 0._r8
  LeafSenescInitialElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)      = -1._r8
  PetolShethChemElmRemobFlx_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX) = 0._r8
  PetolSenescInitialElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX) = -1._r8
  SenecStalkStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)    = 0._r8
  
  D25: DO NB=1,MaxNumBranches
    SapwoodBiomassC_brch(NB,NZ,NY,NX)    = 0._r8
    CanopyLeafSheathC_brch(NB,NZ,NY,NX)  = 0._r8
    PotentialSeedSites_brch(NB,NZ,NY,NX) = 0._r8
    SetNumberSeeds_brch(NB,NZ,NY,NX)     = 0._r8
    SingleGrainMeanBiomC_brch(NB,NZ,NY,NX) = 0._r8
    LeafAreaLive_brch(NB,NZ,NY,NX)       = 0._r8
    NH3Dep2Can_brch(NB,NZ,NY,NX)         = 0._r8
    LeafAreaDying_brch(NB,NZ,NY,NX)      = 0._r8
    CanPBranchHeight(NB,NZ,NY,NX)        = 0._r8
    
    D5: DO L=1,NumCanopyLayers
      CanopyStalkSurfArea_lbrch(L,NB,NZ,NY,NX)=0._r8
      DO N=1,NumLeafInclinationClasses
        StemAreaZsec_brch(N,L,NB,NZ,NY,NX)=0._r8
      enddo
    ENDDO D5

    DO K=0,MaxNodesPerBranch
      LeafArea_node(K,NB,NZ,NY,NX)                               = 0._r8
      StalkNodeHeight_brch(K,NB,NZ,NY,NX)                        = 0._r8
      StalkNodeVertLength_brch(K,NB,NZ,NY,NX)                    = 0._r8
      PetoleLength_node(K,NB,NZ,NY,NX)                           = 0._r8
      LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)       = 0._r8
      PetolShethElmntNode_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX) = 0._r8
      StructInternodeElms_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX) = 0._r8
      LeafProteinC_node(K,NB,NZ,NY,NX)                           = 0._r8
      PetoleProteinC_node(K,NB,NZ,NY,NX)                         = 0._r8

      D55: DO L=1,NumCanopyLayers
        CanopyLeafArea_lnode(L,K,NB,NZ,NY,NX)=0._r8
        LeafLayerElms_node(1:NumPlantChemElms,L,K,NB,NZ,NY,NX)=0._r8
      ENDDO D55

      IF(K.NE.0)THEN
        CPOOL3_node(K,NB,NZ,NY,NX)                = 0._r8
        CMassCO2BundleSheath_node(K,NB,NZ,NY,NX)  = 0._r8
        CMassHCO3BundleSheath_node(K,NB,NZ,NY,NX) = 0._r8
        CPOOL4_node(K,NB,NZ,NY,NX)                = 0._r8
        D45: DO L=1,NumCanopyLayers
          DO N=1,NumLeafInclinationClasses
            LeafAreaZsec_brch(N,L,K,NB,NZ,NY,NX)=0._r8
          enddo
        ENDDO D45
      ENDIF
    enddo
  ENDDO D25

  D35: DO L=1,NumCanopyLayers
    CanopyLeafAreaZ_pft(L,NZ,NY,NX) = 0._r8
    CanopyLeafCLyr_pft(L,NZ,NY,NX)  = 0._r8
    CanopyStemSurfAreaZ_pft(L,NZ,NY,NX) = 0._r8
  ENDDO D35
  CanopyNonstElms_pft(1:NumPlantChemElms,NZ,NY,NX)    = 0._r8
  CanopyNonstElmConc_pft(1:NumPlantChemElms,NZ,NY,NX) = 0._r8

  ShootElms_pft(1:NumPlantChemElms,NZ,NY,NX)     = 0._r8
  LeafStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)      = 0._r8
  PetolShethStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)    = 0._r8
  StalkStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)     = 0._r8
  StalkRsrvElms_pft(1:NumPlantChemElms,NZ,NY,NX)   = 0._r8
  HuskStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)   = 0._r8
  EarStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)    = 0._r8
  GrainStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)  = 0._r8
  RootStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)   = 0._r8
  ShootNoduleElms_pft(1:NumPlantChemElms,NZ,NY,NX) = 0._r8
  RootNoduleElms_pft(1:NumPlantChemElms,NZ,NY,NX)  = 0._r8

  CanopyNoduleNonstCConc_pft(NZ,NY,NX) = 0._r8
  CanopySapwoodC_pft(NZ,NY,NX)       = 0._r8
  CanopyLeafSheathC_pft(NZ,NY,NX)     = 0._r8
  CanopyLeafArea_pft(NZ,NY,NX)       = 0._r8
  RootBiomCPerPlant_pft(NZ,NY,NX)    = 0._r8
  CanopyStemSurfArea_pft(NZ,NY,NX)       = 0._r8
  call PrintInfo('end '//subname)
  end subroutine InitPlantPhenoMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitMassBalance(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
  integer :: M,NE
  character(len=*), parameter :: subname='InitMassBalance'
  associate(                 &
    icwood => pltpar%icwood  &
  )

  call PrintInfo('beg '//subname)
  !
  !     INITIALIZE MASS BALANCE CHECKS
  !
  IF(is_cold_run().AND.is_first_year)THEN
    GrossCO2Fix_pft(NZ,NY,NX)                                   = 0._r8
    SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX) = 0._r8
    GrossResp_pft(NZ,NY,NX)                                     = 0._r8
    CanopyRespC_CumYr_pft(NZ,NY,NX)                             = 0._r8
    PlantExudElm_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)         = 0._r8
    RootUptk_N_CumYr_pft(NZ,NY,NX)                              = 0._r8
    RootUptk_P_CumYr_pft(NZ,NY,NX)                              = 0._r8
    LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)     = 0._r8
    PlantN2Fix_CumYr_pft(NZ,NY,NX)                              = 0._r8
    NH3Dep2Can_pft(NZ,NY,NX)                                    = 0._r8
    NH3Emis_CumYr_pft(NZ,NY,NX)                                 = 0._r8
    CO2ByFire_CumYr_pft(NZ,NY,NX)                               = 0._r8
    CH4ByFire_CumYr_pft(NZ,NY,NX)                               = 0._r8
    O2ByFire_CumYr_pft(NZ,NY,NX)                                = 0._r8
    NH3byFire_CumYr_pft(NZ,NY,NX)                               = 0._r8
    N2ObyFire_CumYr_pft(NZ,NY,NX)                               = 0._r8
    PO4byFire_CumYr_pft(NZ,NY,NX)                               = 0._r8
    EcoHavstElmntCum_pft(1:NumPlantChemElms,NZ,NY,NX)           = 0._r8
    EcoHavstElmnt_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)        = 0._r8
    NetCumElmntFlx2Plant_pft(1:NumPlantChemElms,NZ,NY,NX)       = 0._r8
    ETCanopy_CumYr_pft(NZ,NY,NX)                                = 0._r8
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine InitMassBalance
!------------------------------------------------------------------------------------------

  subroutine InitPlantHeatandWater(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
  real(r8) :: FDM
!
!     INITIALIZE PLANT HEAT AND WATER STATUS
!
!     VHeatCapCanopy_pft=canopy heat capacity (MJ m-3 K-1)
!     TdegCCanopy_pft,TKC=canopy temperature for growth (oC,K)
!     TCGroth_pft,TKGroth_pft=canopy temperature for phenology (oC,K)
!     PSICanopy_pft,PSICanopyOsmo_pft,PSICanopyTurg_pft=canopy total,osmotic,turgor water potl(MPa)
!
  FDM                           = get_FDM(PSICanopy_pft(NZ,NY,NX))
  CanopyBiomWater_pft(NZ,NY,NX) = ppmc*CanopyLeafSheathC_pft(NZ,NY,NX)/FDM
  VHeatCapCanopy_pft(NZ,NY,NX)  = cpw*(ShootElms_pft(ielmc,NZ,NY,NX)*StemSpecVolume_pft(NZ,NY,NX)+CanopyBiomWater_pft(NZ,NY,NX))

  end subroutine InitPlantHeatandWater
!------------------------------------------------------------------------------------------

  subroutine InitRootMychorMorphoBio(NZ,NY,NX)
  use GridConsts, only : MaxNumRootAxes
  implicit none
  integer, intent(in) :: NZ, NY, NX
  character(len=*), parameter :: subname='InitRootMychorMorphoBio'
  integer :: K,L,M,N,NR
  REAL(R8) :: CCO2A
  REAL(R8) :: CCO2P
  REAL(R8) :: COXYA
  REAL(R8) :: COXYP
  call PrintInfo('beg '//subname)
  !
  !     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) MORPHOLOGY AND BIOMASS
  !
  !     PSIRoot_pvr,PSIRootOSMO_vr,PSIRootTurg_vr=root,myco total,osmotic,turgor water potl(MPa)
  !     CO2A,CO2P=root,myco gaseous,aqueous CO2 content (g)
  !     OXYA,OXYP=root,myco gaseous,aqueous O2 content (g)
  !
  NumStructuralRootAxes_pft(NZ,NY,NX) = 0
  RootNH4Uptake_pft(NZ,NY,NX)    = 0._r8
  RootNO3Uptake_pft(NZ,NY,NX)    = 0._r8
  RootH2PO4Uptake_pft(NZ,NY,NX)  = 0._r8
  RootHPO4Uptake_pft(NZ,NY,NX)   = 0._r8
  RootN2Fix_pft(NZ,NY,NX)        = 0._r8

  DO NR=1,MaxNumRootAxes
    RootSegBaseDepth_raxes(NR,NZ,NY,NX)                 = SeedDepth_pft(NZ,NY,NX)
    Root1stDepz_raxes(NR,NZ,NY,NX)                        = SeedDepth_pft(NZ,NY,NX)
    RootMyco1stElm_raxs(1:NumPlantChemElms,NR,NZ,NY,NX) = 0._r8
  ENDDO  

  D20: DO L=1,NL_col(NY,NX)
    D40: DO N=1,pltpar%jroots
      RPlantRootH2OUptk_pvr(N,L,NZ,NY,NX)                       = 0._r8
      PSIRoot_pvr(N,L,NZ,NY,NX)                                  = -0.01_r8
      PSIRootOSMO_vr(N,L,NZ,NY,NX)                               = OrganOsmoPsi0pt_pft(NZ,NY,NX)+PSIRoot_pvr(N,L,NZ,NY,NX)
      PSIRootTurg_vr(N,L,NZ,NY,NX)                               = AZMAX1(PSIRoot_pvr(N,L,NZ,NY,NX)-PSIRootOSMO_vr(N,L,NZ,NY,NX))
      RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ,NY,NX)    = 0._r8
      RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ,NY,NX) = 0._r8
      RootProteinConc_rpvr(N,L,NZ,NY,NX)                         = RootProteinCMax_pft(NZ,NY,NX)
      RootMycoActiveBiomC_pvr(N,L,NZ,NY,NX)                      = 0._r8
      PopuRootMycoC_pvr(N,L,NZ,NY,NX)                            = 0._r8
      RootProteinC_pvr(N,L,NZ,NY,NX)             = 0._r8
      Root2ndXNumL_rpvr(N,L,NZ,NY,NX)              = 0._r8
      RootTotLenPerPlant_pvr(N,L,NZ,NY,NX)          = 0._r8
      RootLenDensPerPlant_pvr(N,L,NZ,NY,NX)      = 0._r8
      RootPoreVol_pvr(N,L,NZ,NY,NX)              = 0._r8
      RootVH2O_pvr(N,L,NZ,NY,NX)                 = 0._r8
      Root1stRadius_pvr(N,L,NZ,NY,NX)            = Root1stMaxRadius_pft(N,NZ,NY,NX)
      Root2ndRadius_rpvr(N,L,NZ,NY,NX)            = Root2ndMaxRadius_pft(N,NZ,NY,NX)
      RootSAreaPerPlant_pvr(N,L,NZ,NY,NX)         = 0._r8
      Root2ndEffLen4uptk_rpvr(N,L,NZ,NY,NX)            = 1.0E-03
      RootNutUptake_pvr(ids_NH4,N,L,NZ,NY,NX)    = 0._r8
      RootNutUptake_pvr(ids_NO3,N,L,NZ,NY,NX)    = 0._r8
      RootNutUptake_pvr(ids_H2PO4,N,L,NZ,NY,NX)  = 0._r8
      RootNutUptake_pvr(ids_H1PO4,N,L,NZ,NY,NX)  = 0._r8
      RootNutUptake_pvr(ids_NH4B,N,L,NZ,NY,NX)   = 0._r8
      RootNutUptake_pvr(ids_NO3B,N,L,NZ,NY,NX)   = 0._r8
      RootNutUptake_pvr(ids_H2PO4B,N,L,NZ,NY,NX) = 0._r8
      RootNutUptake_pvr(ids_H1PO4B,N,L,NZ,NY,NX) = 0._r8
      RootO2Dmnd4Resp_pvr(N,L,NZ,NY,NX)          = 0._r8
      RootNH4DmndSoil_pvr(N,L,NZ,NY,NX)          = 0._r8
      RootNH4DmndBand_pvr(N,L,NZ,NY,NX)          = 0._r8
      RootNO3DmndSoil_pvr(N,L,NZ,NY,NX)          = 0._r8
      RootNO3DmndBand_pvr(N,L,NZ,NY,NX)          = 0._r8
      RootH2PO4DmndSoil_pvr(N,L,NZ,NY,NX)        = 0._r8
      RootH1PO4DmndSoil_pvr(N,L,NZ,NY,NX)        = 0._r8
      RootH2PO4DmndBand_pvr(N,L,NZ,NY,NX)        = 0._r8
      RootH1PO4DmndBand_pvr(N,L,NZ,NY,NX)        = 0._r8
      CCO2A                                           = CCO2EI_gperm3_col(NY,NX)
      CCO2P                                           = 0.030_r8*EXP(-2.621_r8-0.0317_r8*ATCA_col(NY,NX))*CO2EI_col(NY,NX)
      trcg_rootml_pvr(idg_CO2,N,L,NZ,NY,NX)           = CCO2A*RootPoreVol_pvr(N,L,NZ,NY,NX)
      trcs_rootml_pvr(idg_CO2,N,L,NZ,NY,NX)           = CCO2P*RootVH2O_pvr(N,L,NZ,NY,NX)
      trcg_air2root_flx_pvr(idg_CO2,N,L,NZ,NY,NX)     = 0._r8
      trcg_Root_gas2aqu_flx_vr(idg_CO2,N,L,NZ,NY,NX)  = 0._r8
      RootUptkSoiSol_pvr(idg_CO2,N,L,NZ,NY,NX)         = 0._r8
      RCO2Emis2Root_rpvr(N,L,NZ,NY,NX)                   = 0._r8
      COXYA                                           = AtmGasCgperm3_col(idg_O2,NY,NX)
      COXYP                                           = 0.032_r8*EXP(-6.175_r8-0.0211_r8*ATCA_col(NY,NX))*OXYE_col(NY,NX)
      trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ,NY,NX) = 0._r8
      trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ,NY,NX) = 0._r8
      trcg_rootml_pvr(idg_O2,N,L,NZ,NY,NX)            = COXYA*RootPoreVol_pvr(N,L,NZ,NY,NX)
      trcs_rootml_pvr(idg_O2,N,L,NZ,NY,NX)            = COXYP*RootVH2O_pvr(N,L,NZ,NY,NX)
      RAutoRootO2Limter_rpvr(N,L,NZ,NY,NX)            = 1.0_r8

      D30: DO NR=1,MaxNumRootAxes
        Root2ndXNum_rpvr(N,L,NR,NZ,NY,NX)                             = 0._r8
        Root2ndLen_rpvr(N,L,NR,NZ,NY,NX)                              = 0._r8
        RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ,NY,NX) = 0._r8
      ENDDO D30
    ENDDO D40

    Root1stXNumL_pvr(L,NZ,NY,NX)             = 0._r8      
    DO NR=1,MaxNumRootAxes
      RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,L,NR,NZ,NY,NX) = 0._r8
      Root1stActStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ,NY,NX) = 0._r8
      Root1stLigStructElms_rpvr(1:NumPlantChemElms,L,NR,NZ,NY,NX) = 0._r8
      Root1stLenPP_rpvr(L,NR,NZ,NY,NX)                              = 0._r8
      RootAge_rpvr(L,NR,NZ,NY,NX)                                 = 0._r8
    ENDDO

    D6400: DO K=1,pltpar%NumOfPlantLitrCmplxs
      DO  M=1,jskenc
        LitrfallElms_pvr(1:NumPlantChemElms,M,K,L,NZ,NY,NX)=0._r8
      enddo
    ENDDO D6400

    RootNodulNonstElms_rpvr(1:NumPlantChemElms,L,NZ,NY,NX) = 0._r8
    RootNodulStrutElms_rpvr(1:NumPlantChemElms,L,NZ,NY,NX) = 0._r8
    RootN2Fix_pvr(L,NZ,NY,NX)                              = 0._r8

  ENDDO D20

  RootNutUptake_pvr(ids_NH4,1:2,NL_col(NY,NX)+1:JZ,NZ,NY,NX)    = 0._r8
  RootNutUptake_pvr(ids_NH4B,1:2,NL_col(NY,NX)+1:JZ,NZ,NY,NX)   = 0._r8
  RootNutUptake_pvr(ids_H2PO4,1:2,NL_col(NY,NX)+1:JZ,NZ,NY,NX)  = 0._r8
  RootNutUptake_pvr(ids_H2PO4B,1:2,NL_col(NY,NX)+1:JZ,NZ,NY,NX) = 0._r8
  RootLenDensPerPlant_pvr(1:2,NL_col(NY,NX)+1:JZ,NZ,NY,NX)      = 0._r8
  call PrintInfo('end '//subname)
  end subroutine InitRootMychorMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitSeedMorphoBio(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
  character(len=*), parameter :: subname='InitSeedMorphoBio'
  REAL(R8) :: FDM
  INTEGER :: NR,L,NB,NE

  call PrintInfo('beg '//subname)
!
!     INITIALIZE SEED MORPHOLOGY AND BIOMASS
!
!
  SeedPlantedElm_pft(ielmc,NZ,NY,NX)    = SeedCMass_pft(NZ,NY,NX)*PlantPopuLive_pft(NZ,NY,NX)

  LeafStrutElms_brch(ielmn,1,NZ,NY,NX)  = rNCGrain_pft(NZ,NY,NX)*LeafStrutElms_brch(ielmc,1,NZ,NY,NX)
  LeafStrutElms_brch(ielmp,1,NZ,NY,NX)  = rPCGrain_pft(NZ,NY,NX)*LeafStrutElms_brch(ielmc,1,NZ,NY,NX)
  CanopyLeafSheathC_brch(1,NZ,NY,NX)    = LeafStrutElms_brch(ielmc,1,NZ,NY,NX)+PetolShethStrutElms_brch(ielmc,1,NZ,NY,NX)
  CanopyLeafSheathC_pft(NZ,NY,NX)        = CanopyLeafSheathC_pft(NZ,NY,NX)+CanopyLeafSheathC_brch(1,NZ,NY,NX)
    
  WatHeldOnCanopy_pft(NZ,NY,NX)          = 0._r8
  CanopyNonstElms_brch(ielmn,1,NZ,NY,NX) = rNCGrain_pft(NZ,NY,NX)*CanopyNonstElms_brch(ielmc,1,NZ,NY,NX)
  CanopyNonstElms_brch(ielmp,1,NZ,NY,NX) = rPCGrain_pft(NZ,NY,NX)*CanopyNonstElms_brch(ielmc,1,NZ,NY,NX)
  
  RootMyco1stStrutElms_rpvr(ielmn,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)= &
    rNCGrain_pft(NZ,NY,NX)*RootMyco1stStrutElms_rpvr(ielmc,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)
  RootMyco1stStrutElms_rpvr(ielmp,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)= &
    rPCGrain_pft(NZ,NY,NX)*RootMyco1stStrutElms_rpvr(ielmc,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)
  Root1stActStructElms_rpvr(:,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)=RootMyco1stStrutElms_rpvr(:,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)   
  RootMyco1stElm_raxs(ielmn,1,NZ,NY,NX)=rNCGrain_pft(NZ,NY,NX)*RootMyco1stElm_raxs(ielmc,1,NZ,NY,NX)
  RootMyco1stElm_raxs(ielmp,1,NZ,NY,NX)=rPCGrain_pft(NZ,NY,NX)*RootMyco1stElm_raxs(ielmc,1,NZ,NY,NX)
  RootMycoActiveBiomC_pvr(ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)= &
    RootMyco1stStrutElms_rpvr(ielmc,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)
  PopuRootMycoC_pvr(ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)= &
    RootMyco1stStrutElms_rpvr(ielmc,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)
  RootProteinC_pvr(ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)= &
    RootMycoActiveBiomC_pvr(ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)*RootProteinCMax_pft(NZ,NY,NX)
  RootMycoNonstElms_rpvr(ielmn,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)= rNCGrain_pft(NZ,NY,NX)&
    *RootMycoNonstElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)
  RootMycoNonstElms_rpvr(ielmp,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)=rPCGrain_pft(NZ,NY,NX) &
    *RootMycoNonstElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)

  !shoots
  NB=1
  DO NE=1,NumPlantChemElms
    ShootElms_brch(NE,NB,NZ,NY,NX)=LeafStrutElms_brch(NE,NB,NZ,NY,NX) &
      +PetolShethStrutElms_brch(NE,NB,NZ,NY,NX)+StalkStrutElms_brch(NE,NB,NZ,NY,NX) &
      +StalkRsrvElms_brch(NE,NB,NZ,NY,NX)+HuskStrutElms_brch(NE,NB,NZ,NY,NX)&
      +EarStrutElms_brch(NE,NB,NZ,NY,NX)+GrainStrutElms_brch(NE,NB,NZ,NY,NX) &
      +CanopyNonstElms_brch(NE,NB,NZ,NY,NX)              
  ENDDO
  
  RootElms_pft(1:NumPlantChemElms,NZ,NY,NX)=0._r8

  call PrintInfo('end '//subname)
  end subroutine InitSeedMorphoBio

  end module StartqMod
