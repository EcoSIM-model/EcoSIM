module UptakesMod
  ! Description:
  ! Solve plant water uptake, canopy energy balance and coupled gas exchange.
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use abortutils, only : endrun, iulog
  use data_kind_mod , only : r8 => DAT_KIND_R8,yearIJ_type
  use data_const_mod, only: GravAcceleration=>DAT_CONST_G
  use StomatesMod,    only: StomatalDynamics, PhotosynsDiag
  use EcoSIMCtrlMod,  only: etimer, ldo_sp_mode, ats_cpl_mode
  use UnitMod,        only: units
  use PlantBalMod,    only: SumPlantRootGas
  use MiniFuncMod,    only: gOC_to_m3_OM
  use PlantBGCPars
  use minimathmod
  use DebugToolMod
  use EcosimConst
  use EcoSIMSolverPar
  use UptakePars
  use NutUptakeMod
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantPhotosynthesisAPIData, only : plt_photo
  use PlantRadiationAPIData, only : plt_rad
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantSoilChemistryAPIData, only : plt_soilchem
  use PlantBiomassAPIData, only : plt_biom
  use PlantEnergyWaterAPIData, only : plt_ew
  use PlantRootBGCAPIData, only : plt_rbgc
  use PlantMathFuncMod
  use ElmIDMod, only : itrue
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public :: RootUptakes
  public :: InitUptake
  public :: SolvePassiveCanopyEnergy,PassiveCanopyResult
  public :: CanopyStemThermalCarbon,CanopyExchangeAreaFraction
  public :: CanopyVaporResistance


  type :: PassiveCanopyResult
    ! Description:
    ! Accepted one-hour passive surface state and integrated energy/water fluxes.
    real(r8) :: temperature,water,snow,capacity,resistance ! K, m3, m3 SWE, MJ/K, h/m
    real(r8) :: net_radiation,sensible,latent,vapor_enthalpy ! MJ per PFT/cell/hour
    real(r8) :: water_flux,snow_flux ! m3 per PFT/cell/hour, positive into canopy
  end type PassiveCanopyResult

  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine InitUptake

  implicit none

  call InitUptakePars

  end subroutine InitUptake


!----------------------------------------------------------------------------------------------------
  subroutine RootUptakes(yearIJ,LiveExchange_pft,PassiveExchange_pft)
  ! Description:
  ! Compute canopy atmospheric exchange and root energy/C/N/P transfers. Live
  ! exchange follows biological/hydraulic eligibility; standing dead exchanges
  ! heat and intercepted water independently, using the same exposure geometry.
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  logical, intent(out) :: PassiveExchange_pft(JP1) ! PFTs producing passive exchanges before management [-]
  logical, intent(out) :: LiveExchange_pft(JP1) ! PFTs producing live exchanges this hour, before management [-]


  character(len=*), parameter :: subname='RootUptakes'
  integer :: NN,N,NZ,K,L
  real(r8) :: FracGridPFTCovered
  real(r8) :: CdHeatCanopyAir
  real(r8) :: CanopyMassC4H2OStorage                 !Canopy C mass, g/m2
  real(r8) :: TotalSoilPSIMPa_vr(JZ1)     !total soil matric pressure, removing elevation adjustment [MPa]
  real(r8) :: RadialMeanLen_rvr(pltpar%jroots,JZ1)
  real(r8) :: FineRootRadius_rvr(pltpar%jroots,JZ1)
  real(r8) :: SoilResist4H2O_rvr(pltpar%jroots,JZ1)
  real(r8) :: SoilRootResist4H2O_pvr(pltpar%jroots,JZ1)
  real(r8) :: AllRootC_vr(JZ1)
  real(r8) :: FracPRoot4Uptake_pvr(pltpar%jroots,JZ1,JP1)
  real(r8) :: FracMinRoot4Uptake_rpvr(pltpar%jroots,JZ1,JP1)
  real(r8) :: FracSoilLBy1stRoots_pvr(JZ1,JP1)
  real(r8) :: RootEffLen4Absorption_pvr(pltpar%jroots,JZ1)
  real(r8) :: AirMicPore4Fill_vr(JZ1)  !
  real(r8) :: WatAvail4Uptake_vr(JZ1)
  real(r8) :: TKCX,CdH2ORootxSoil,PrecpHeatbyCanopy
  real(r8) :: VHeatCapCanopyPrev_pft  !canopy heat capacity [MJ/cell/K]
  real(r8) :: DIFF
  real(r8) :: PSIGravCanopyHeight               !gravitation potential at effective canopy height [MPa]
  real(r8) :: cumPRootH2OUptake,HeatEvapSens
  real(r8) :: CumPlantHeatLoss2Soil,RawH2OAtm2CanopySinkZ_pft,RawHeatAtm2CanopySinkZ_pft
  real(r8) :: FDMP
  real(r8) :: DeadHeatCapacity,LiveAreaFraction,DeadAreaColumn,PFTArea,PriorTemperature
  logical  :: HydroActivePlant
  logical  :: HydraulicPathAvailable
  logical  :: SoiLayerHasRoot_rvr(pltpar%jroots,JZ1)
!     begin_execution
  call PrintInfo('beg '//subname)
  associate(                                                        &
    HeatSensVapAir2Canopy_col => plt_ew%HeatSensVapAir2Canopy_col, &
    HeatPrecAir2Canopy_col => plt_ew%HeatPrecAir2Canopy_col, &
    AREA3                     => plt_site%AREA3                    ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    canopy_exposed_col        => plt_ew%canopy_exposed_col         ,& !input  :hourly exposed-canopy diagnostic, [-]
    CanopyBiomWater_pft       => plt_ew%CanopyBiomWater_pft        ,& !input  :canopy water content, [m3 d-2]
    CanopyStemSurfArea_pft    => plt_morph%CanopyStemSurfArea_pft  ,& !input  :plant stem area, [m2 d-2]
    RawH2OCanopy2Atm_pft      => plt_photo%RawH2OCanopy2Atm_pft    ,& !input  :canopy boundary layer resistance, [h m-1]
    RawHeatCanopy2Atm_pft     => plt_photo%RawHeatCanopy2Atm_pft   ,& !input  :canopy sensible-heat boundary layer resistance, [h m-1]
    CanopyLeafArea_col        => plt_morph%CanopyLeafArea_col      ,& !input  :grid canopy leaf area, [m2 d-2]
    CanopyLeafSheathC_pft     => plt_biom%CanopyLeafSheathC_pft    ,& !input  :canopy leaf + sheath C, [g d-2]
    CanopySapwoodC_pft        => plt_biom%CanopySapwoodC_pft       ,& !input  :canopy active stalk C, [g d-2]
    CumSoilThickness_vr       => plt_site%CumSoilThickness_vr      ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    FracPARads2Canopy_pft     => plt_rad%FracPARads2Canopy_pft     ,& !input  :fraction of incoming PAR absorbed by canopy, [-]
    HeatXAir2PCan_pft         => plt_ew%HeatXAir2PCan_pft          ,& !input  :canopy sensible heat flux, [MJ d-2 h-1]
    IsPlantActive_pft         => plt_pheno%IsPlantActive_pft       ,& !input  :flag for living pft, [-]
    LeafStalkAreaAll_col      => plt_morph%LeafStalkAreaAll_col    ,& !input  :stalk area of combined, each PFT canopy,[m^2 d-2]
    LeafStalkAreaAct_pft      => plt_morph%LeafStalkAreaAct_pft    ,& !input  :plant leaf+stem/stalk area, [m2 d-2]
    StandDeadSurfAreaAct_pft  => plt_morph%StandDeadSurfAreaAct_pft,& !input  :exposed standing dead area, [m2]
    StandDeadSurfArea_pft     => plt_morph%StandDeadSurfArea_pft, & !input: total dead area, including buried stems [m2]
    StandDeadStrutElms_pft    => plt_biom%StandDeadStrutElms_pft   ,& !input  :aboveground standing dead C/N/P, [g]
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft       ,& !input  :number of main branch,[-]
    NP                        => plt_site%NP                       ,& !input  :current number of plant species,[-]
    NU                        => plt_site%NU                       ,& !input  :current soil surface layer number, [-]
    PlantPopu_col             => plt_site%PlantPopu_col            ,& !input  :total plant population, [plants d-2]
    PlantPopuLive_pft         => plt_site%PlantPopuLive_pft        ,& !input  :plant population, [d-2]
    RainIntcptByCanopy_pft    => plt_ew%RainIntcptByCanopy_pft     ,& !input  :water flux onto canopy, [m3 d-2 h-1]
    SnowIntcptByCanopy_pft    => plt_ew%SnowIntcptByCanopy_pft     ,& !input  :snow flux onto canopy, [m3 d-2 h-1]
    Root1stDepz_raxes         => plt_morph%Root1stDepz_raxes       ,& !input  :root layer depth, [m]
    SeedDepth_pft             => plt_morph%SeedDepth_pft           ,& !input  :seeding depth, [m]
    TKC_pft                   => plt_ew%TKC_pft                    ,& !input  :canopy temperature, [K]
    StemSpecVolume_pft        => plt_morph%StemSpecVolume_pft      ,& !input  :stalk specific volume, [m3 gC-1]
    TairK                     => plt_ew%TairK                      ,& !input  :air temperature, [K]
    WatHeldOnCanopy_pft       => plt_ew%WatHeldOnCanopy_pft        ,& !input  :canopy surface water content, [m3 d-2]
    SnowOnCanopy_pft          => plt_ew%SnowOnCanopy_pft           ,& !input  :snow water equivalent held on canopy, [m3 d-2]
    ZERO4LeafVar_pft          => plt_biom%ZERO4LeafVar_pft         ,& !input  :threshold zero for leaf calculation, [-]
    ZEROS                     => plt_site%ZEROS                    ,& !input  :threshold zero for numerical stability,[-]
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch     ,& !input  :plant growth stage, [-]
    RawIsoTAtm2CanopySinkZ_col=> plt_ew%RawIsoTAtm2CanopySinkZ_col ,& !input  :isothermal aerodynamic resistance between canopy top and wind ref height in atmosphere, [h m-1]
    RawIsoTCanopy2Atm_pft     => plt_ew%RawIsoTCanopy2Atm_pft      ,& !input  :!isothermal aerodynamic resistance from canopy to free atmosphere (where wind speed is measured)
    Air_Heat_Latent_store_col => plt_ew%Air_Heat_Latent_store_col  ,& !inoput :total latent heat flux x boundary layer resistance, [MJ m-1]
    Air_Heat_Sens_store_col   => plt_ew%Air_Heat_Sens_store_col    ,& !inoput :total sensible heat flux x boundary layer resistance, [MJ m-1]
    PSICanopy_pft             => plt_ew%PSICanopy_pft              ,& !inoput :canopy total water potential, [Mpa]
    CanopyEvapTransLHeat_pft  => plt_ew%CanopyEvapTransLHeat_pft   ,& !output  :canopy latent heat flux, [MJ d-2 h-1]
    Transpiration_pft         => plt_ew%Transpiration_pft          ,& !output :canopy transpiration, [m2 d-2 h-1]
    TdegCCanopy_pft           => plt_ew%TdegCCanopy_pft            ,& !output :canopy temperature, [oC]
    VapXAir2Canopy_pft        => plt_ew%VapXAir2Canopy_pft          & !output :canopy evaporation, [m2 d-2 h-1]
  )

  HeatSensVapAir2Canopy_col=0._r8
  HeatPrecAir2Canopy_col=0._r8
  LiveExchange_pft=.false.
  PassiveExchange_pft=.false.
  call ZeroNutrientUptake

  call PrepNutrientUptake(yearIJ,TotalSoilPSIMPa_vr,AllRootC_vr,AirMicPore4Fill_vr,WatAvail4Uptake_vr)
  !
  !     IF PLANT SPECIES EXISTS
  if(NP.GT.0)THEN
    DeadAreaColumn   = sum(StandDeadSurfAreaAct_pft(1:NP))

    DO NZ=1,NP
      DeadHeatCapacity = cpo*gOC_to_m3_OM(CanopyStemThermalCarbon(StandDeadStrutElms_pft(ielmc,NZ), &
        StandDeadSurfArea_pft(NZ),StemSpecVolume_pft(NZ)))
      PFTArea          = LeafStalkAreaAct_pft(NZ)+StandDeadSurfAreaAct_pft(NZ)
      LiveAreaFraction = 1._r8

      if(PFTArea>ZEROS)LiveAreaFraction=LeafStalkAreaAct_pft(NZ)/PFTArea

      FracGridPFTCovered=0._r8

      if(LeafStalkAreaAll_col>ZEROS) &
        FracGridPFTCovered=CanopyExchangeAreaFraction(PFTArea,LeafStalkAreaAll_col, &
          CanopyLeafArea_col,DeadAreaColumn,AREA3(NU))

      IF(IsPlantActive_pft(NZ).EQ.iTrue .AND. PlantPopuLive_pft(NZ).GT.0.0_r8)THEN
        LiveExchange_pft(NZ)=.true.

        call PhotosynsDiag(yearIJ%I,yearIJ%J,NZ)

        call UpdateCanopyProperty(yearIJ%I,yearIJ%J,NZ)

        !     STOMATE=solve for minimum canopy stomatal resistance

        CALL StomatalDynamics(yearIJ,NZ)
        !
        call UpdateRootProperty(yearIJ,NZ,RadialMeanLen_rvr,FineRootRadius_rvr,AllRootC_vr,FracPRoot4Uptake_pvr,&
          FracMinRoot4Uptake_rpvr,FracSoilLBy1stRoots_pvr,RootEffLen4Absorption_pvr)
        !
        !     CALCULATE CANOPY WATER STATUS FROM CONVERGENCE SOLUTION FOR
        !     TRANSPIRATION - ROOT WATER UPTAKE = CHANGE IN CANOPY WATER CONTENT
        !     [SpecStalkVolume]=[m3/gC], assuming mean sapwood thicknees is 2.5 mm.
        CanopyMassC4H2OStorage = CanopyLeafSheathC_pft(NZ)+CanopyStemThermalCarbon(CanopySapwoodC_pft(NZ), &
          CanopyStemSurfArea_pft(NZ),StemSpecVolume_pft(NZ))

        ! A negative reservoir cannot be repaired by changing pressure. Reject it
        ! before exchanging water; silently clipping it would create water.
        if (.not.ieee_is_finite(CanopyBiomWater_pft(NZ)) .or. CanopyBiomWater_pft(NZ)<0._r8 &
            .or. .not.ieee_is_finite(CanopyMassC4H2OStorage)) then
          write(iulog,*)subname,': invalid canopy hydraulic storage; year/day/hour/PFT, water, storage C: ', &
            yearIJ%year,yearIJ%I,yearIJ%J,NZ,CanopyBiomWater_pft(NZ),CanopyMassC4H2OStorage
          call endrun('Invalid canopy hydraulic storage; use a valid checkpoint or an accounted water correction',__LINE__)
        endif

        CanopyMassC4H2OStorage=MAX(0._r8,CanopyMassC4H2OStorage)
        VHeatCapCanopyPrev_pft = cpo*gOC_to_m3_OM(CanopyMassC4H2OStorage)+DeadHeatCapacity &
          +cpw*(WatHeldOnCanopy_pft(NZ)+CanopyBiomWater_pft(NZ)) &
          +cps*SnowOnCanopy_pft(NZ)

        if(ldo_sp_mode)then
          ! Prescribed PFTs do not initialize the full-model branch calendar.
          HydroActivePlant=LeafStalkAreaAct_pft(NZ)>ZERO4LeafVar_pft(NZ)
        else
          HydroActivePlant=(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0)               &  !plant emerged
            .AND.(LeafStalkAreaAct_pft(NZ).GT.ZERO4LeafVar_pft(NZ).AND.FracPARads2Canopy_pft(NZ).GT.0.0_r8)  &  !active canopy
            .AND.(Root1stDepz_raxes(1,NZ).GT.SeedDepth_pft(NZ)+CumSoilThickness_vr(0))                          !active root
        endif
        !
        ! LeafStalkAreaAct is snow-exposed per-PFT area prepared at the same hour
        ! as the column flag and precipitation interception. Keep this geometry
        ! fixed through ground substeps, even if snow depth changes during them.
        HydroActivePlant=HydroActivePlant .and. canopy_exposed_col

        IF(HydroActivePlant)THEN

          call CalcPlantHydroResistance(yearIJ,NZ,RadialMeanLen_rvr,FineRootRadius_rvr,RootEffLen4Absorption_pvr,&
            SoilResist4H2O_rvr,SoilRootResist4H2O_pvr,CdH2ORootxSoil,PSIGravCanopyHeight,SoiLayerHasRoot_rvr)

          ! Second-stage eligibility: biological activity alone does not imply
          ! a connected hydraulic reservoir. Zero water is allowed to refill;
          ! small positive conductance and storage C retain their physical scale.
          if (.not.ieee_is_finite(CdH2ORootxSoil) .or. CdH2ORootxSoil<0._r8) then
            write(iulog,*)subname,': invalid root conductance; day/hour/PFT, conductance: ', &
              yearIJ%I,yearIJ%J,NZ,CdH2ORootxSoil
            call endrun('Invalid soil-canopy hydraulic conductance',__LINE__)
          endif

          HydraulicPathAvailable=CanopyMassC4H2OStorage>0._r8 .and. &
            CdH2ORootxSoil>0._r8 .and. any(SoiLayerHasRoot_rvr)
          if (.not.HydraulicPathAvailable) SoiLayerHasRoot_rvr=.false.
          !
          !       INITIALIZE CANOPY WATER POTENTIAL, OTHER VARIABLES USED IN ENERGY
          !       BALANCE THAT DON'T NEED TO BE RECALCULATED DURING CONVERGENCE
          !
          if (HydraulicPathAvailable) then
            PSICanopy_pft(NZ) = AMIN1(-ppmc,0.667_r8*PSICanopy_pft(NZ))
          else
            !canopy evaporation/condensation can still occur with standing dead
            PSICanopy_pft(NZ) = TotalSoilPSIMPa_vr(NU)
          endif
          PrecpHeatbyCanopy      = (RainIntcptByCanopy_pft(NZ)*cpw+SnowIntcptByCanopy_pft(NZ)*cps)*TairK
          !
          IF(LeafStalkAreaAll_col.GT.ZEROS)THEN
            !the grid has significant canopy (leaf+steam) area
            FracGridPFTCovered=CanopyExchangeAreaFraction(PFTArea,LeafStalkAreaAll_col, &
              CanopyLeafArea_col,DeadAreaColumn,AREA3(NU))
          ELSEIF(PlantPopu_col.GT.ZEROS)THEN
            !total population is > 0
            FracGridPFTCovered=PlantPopuLive_pft(NZ)/PlantPopu_col
          ELSE
            FracGridPFTCovered=1.0_r8/NP
          ENDIF

          TKCX   = TKC_pft(NZ)
          !
          !     CONVERGENCE SOLUTION
          !
          NN=CanopyEnergyH2OIter_func(yearIJ,NZ,FracGridPFTCovered,CanopyMassC4H2OStorage,&
            TotalSoilPSIMPa_vr,CdHeatCanopyAir,DIFF,cumPRootH2OUptake,CumPlantHeatLoss2Soil,&
            HeatEvapSens,FDMP,SoilRootResist4H2O_pvr,FracPRoot4Uptake_pvr,AirMicPore4Fill_vr,&
            WatAvail4Uptake_vr,TKCX,CdH2ORootxSoil,VHeatCapCanopyPrev_pft,PrecpHeatbyCanopy,&
            PSIGravCanopyHeight,SoiLayerHasRoot_rvr,HydraulicPathAvailable,LiveAreaFraction)
          !
          !     IF CONVERGENCE NOT ACHIEVED (RARE), SET DEFAULT
          !     TEMPERATURES, ENERGY FLUXES, WATER POTENTIALS, RESISTANCES
          !
          ! Reject stalled solves and exchanges that would overdraw the reservoir.
          if (.not.ieee_is_finite(DIFF)) DIFF=1._r8
          if (.not.ieee_is_finite(Transpiration_pft(NZ)-cumPRootH2OUptake)) DIFF=1._r8
          if (CanopyBiomWater_pft(NZ)+Transpiration_pft(NZ)-cumPRootH2OUptake<0._r8) DIFF=MAX(1._r8,DIFF)
          IF(NN.GE.MaxIterNum .or. DIFF.GE.DIFFY)THEN
            ! Reject unaccepted root exchange, preserving the intentional lag
            ! of earlier transfers. The wet canopy still needs a conservative
            ! thermal solution, rather than zero fluxes and an arbitrary reset.
            DIFF=MAX(DIFF,1._r8)
            call HandlingDivergence(yearIJ,NZ,TotalSoilPSIMPa_vr,DIFF,FDMP,cumPRootH2OUptake,HeatEvapSens, &
              CumPlantHeatLoss2Soil)
            TKC_pft(NZ)=TKCX
            call PassiveCanopyExchange(NZ,FracGridPFTCovered,CanopyMassC4H2OStorage)
            PassiveExchange_pft(NZ)=.true.
          else
            call UpdatePlantWaterVars(NZ,CdHeatCanopyAir,CanopyMassC4H2OStorage,TotalSoilPSIMPa_vr,SoilResist4H2O_rvr,SoilRootResist4H2O_pvr,&
              TKCX,VHeatCapCanopyPrev_pft,PrecpHeatbyCanopy,cumPRootH2OUptake,CumPlantHeatLoss2Soil,HeatEvapSens,SoiLayerHasRoot_rvr,DeadHeatCapacity)
            HeatSensVapAir2Canopy_col=HeatSensVapAir2Canopy_col+HeatEvapSens
            HeatPrecAir2Canopy_col=HeatPrecAir2Canopy_col+PrecpHeatbyCanopy
          endif
          !
          !     DEFAULT VALUES IF PLANT SPECIES DOES NOT EXIST
          !
        ELSE
          PriorTemperature=TKC_pft(NZ)
          call HandleBareSoil(NZ,TotalSoilPSIMPa_vr,FDMP,VHeatCapCanopyPrev_pft)
          if(StandDeadStrutElms_pft(ielmc,NZ)>0._r8)then
            TKC_pft(NZ)=PriorTemperature
            call PassiveCanopyExchange(NZ,FracGridPFTCovered,CanopyMassC4H2OStorage)
            PassiveExchange_pft(NZ)=.true.
          endif
        ENDIF
        !
        RawH2OAtm2CanopySinkZ_pft  = RawH2OCanopy2Atm_pft(NZ)*RawIsoTAtm2CanopySinkZ_col/RawIsoTCanopy2Atm_pft(NZ)
        RawHeatAtm2CanopySinkZ_pft = RawHeatCanopy2Atm_pft(NZ)*RawIsoTAtm2CanopySinkZ_col/RawIsoTCanopy2Atm_pft(NZ)
        Air_Heat_Latent_store_col  = Air_Heat_Latent_store_col+CanopyEvapTransLHeat_pft(NZ)*RawH2OAtm2CanopySinkZ_pft
        Air_Heat_Sens_store_col    = Air_Heat_Sens_store_col+HeatXAir2PCan_pft(NZ)*RawHeatAtm2CanopySinkZ_pft
        TdegCCanopy_pft(NZ)        = units%Kelvin2Celcius(TKC_pft(NZ))

        if(.not.ldo_sp_mode) then
          call SetCanopyGrowthFuncs(yearIJ,NZ)

          call PlantNutientO2Uptake(yearIJ,NZ,FDMP,RadialMeanLen_rvr,FineRootRadius_rvr,FracPRoot4Uptake_pvr,&
            FracMinRoot4Uptake_rpvr,FracSoilLBy1stRoots_pvr,RootEffLen4Absorption_pvr)
        endif
      ELSE
        ! This PFT skips all root uptake, so discard its previous O2 demand.
        plt_rbgc%RootO2Dmnd4Resp_pvr(:,:,NZ) = 0._r8
        if(StandDeadStrutElms_pft(ielmc,NZ)>0._r8)then
          call UpdateCanopyProperty(yearIJ%I,yearIJ%J,NZ)
          call PassiveCanopyExchange(NZ,FracGridPFTCovered,0._r8)
          PassiveExchange_pft(NZ)=.true.

          RawH2OAtm2CanopySinkZ_pft  = RawH2OCanopy2Atm_pft(NZ)*RawIsoTAtm2CanopySinkZ_col/RawIsoTCanopy2Atm_pft(NZ)
          RawHeatAtm2CanopySinkZ_pft = RawHeatCanopy2Atm_pft(NZ)*RawIsoTAtm2CanopySinkZ_col/RawIsoTCanopy2Atm_pft(NZ)
          Air_Heat_Latent_store_col  = Air_Heat_Latent_store_col+CanopyEvapTransLHeat_pft(NZ)*RawH2OAtm2CanopySinkZ_pft
          Air_Heat_Sens_store_col    = Air_Heat_Sens_store_col+HeatXAir2PCan_pft(NZ)*RawHeatAtm2CanopySinkZ_pft
          TdegCCanopy_pft(NZ)        = units%Kelvin2Celcius(TKC_pft(NZ))
        endif
      ENDIF
    ENDDO
  ENDIF

  end associate
  call PrintInfo('end '//subname)
  END subroutine RootUptakes

!------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------
  function CanopyExchangeAreaFraction(PFTArea,TotalArea,LeafArea,DeadArea,GroundArea) result(Fraction)
    ! Description:
    ! Partition effective projected cover among exposed PFTs. The inherited 0.5
    ! area-to-cover coefficient uses live leaf plus dead area continuously.
    real(r8), intent(in) :: PFTArea,TotalArea ! PFT and column exposed live+dead surface area [m2/cell]
    real(r8), intent(in) :: LeafArea,DeadArea ! column leaf area and exposed dead stem area [m2/cell]
    real(r8), intent(in) :: GroundArea ! column ground area [m2/cell]
    real(r8) :: Fraction ! PFT effective exchange area / ground area [-]
    character(len=*), parameter :: subname='CanopyExchangeAreaFraction'
    call PrintInfo('beg '//subname)
    Fraction=0._r8
    if(TotalArea>0._r8)Fraction=PFTArea/TotalArea*MIN(1._r8,0.5_r8*(LeafArea+DeadArea)/GroundArea)
    call PrintInfo('end '//subname)
  end function CanopyExchangeAreaFraction

!------------------------------------------------------------------------------------------
  function CanopyVaporResistance(NeutralResistance,Rib,AirTemperature,CanopyTemperature) result(Resistance)
    ! Description:
    ! Apply the common Richardson correction and physical minimum to canopy
    ! vapor transport. Iteration damping belongs to the caller, not this law.
    real(r8), intent(in) :: NeutralResistance ! neutral canopy-to-air resistance [h/m]
    real(r8), intent(in) :: Rib ! Richardson-number prefactor [1/K]
    real(r8), intent(in) :: AirTemperature,CanopyTemperature ! air and canopy temperatures [K]
    real(r8) :: Resistance ! stability-corrected vapor resistance [h/m]
    real(r8) :: Ri
    character(len=*), parameter :: subname='CanopyVaporResistance'
    call PrintInfo('beg '//subname)
    Ri=RichardsonNumber(Rib,AirTemperature,CanopyTemperature)
    Resistance=MAX(RawMinH2OCanopy2Atm_pft,NeutralResistance/(1._r8-10._r8*Ri))
    call PrintInfo('end '//subname)
  end function CanopyVaporResistance
!------------------------------------------------------------------------------------------
  subroutine SolvePassiveCanopyEnergy(T0,Ta,qa,Cdry,Wtissue,W0,S0,Rain,Snowfall,Area,SW,LW,F,Ra,Rib,Result)
    ! Description:
    ! Solve a one-hour passive canopy energy balance with accepted, donor-limited
    ! liquid/snow vapor fluxes and their sensible enthalpy. No root or stomatal
    ! state is accessed. All extensive quantities are per PFT and grid cell.
    real(r8), intent(in) :: T0,Ta                    ! previous canopy and reference-air temperatures [K]
    real(r8), intent(in) :: qa                       ! reference vapor water-equivalent concentration [m3/m3]
    real(r8), intent(in) :: Cdry                     ! live plus dead dry-material heat capacity [MJ/K]
    real(r8), intent(in) :: Wtissue,W0,S0            ! tissue liquid, intercepted liquid and snow SWE [m3]
    real(r8), intent(in) :: Rain,Snowfall            ! intercepted precipitation during this hour [m3]
    real(r8), intent(in) :: Area                     ! effective projected exchange area [m2]
    real(r8), intent(in) :: SW,LW,F                  ! absorbed SW, incident LW [MJ], one-side emission factor [MJ/K4]
    real(r8), intent(in) :: Ra,Rib                   ! neutral vapor resistance [h/m], Richardson prefactor [1/K]
    type(PassiveCanopyResult), intent(out) :: Result ! accepted state and integrated fluxes
    character(len=*), parameter :: subname='SolvePassiveCanopyEnergy'
    type(PassiveCanopyResult) :: Trial
    real(r8) :: Cold,Eold,Qprec,Tlo,Thi,Tmid,Flo,Fhi,Residual,Tolerance
    integer :: iteration
    integer, parameter :: itmax=100

    call PrintInfo('beg '//subname)
    if(any(.not.ieee_is_finite([T0,Ta,qa,Cdry,Wtissue,W0,S0,Rain,Snowfall,Area,SW,LW,F,Ra,Rib])) &
      .or. MIN(Cdry,Wtissue,W0,S0,Rain,Snowfall,Area,F)<0._r8 .or. Ra<=0._r8)then
      write(iulog,*)subname,': invalid passive canopy inputs'
      call endrun('Invalid passive canopy state',__LINE__)
    endif

    Cold      = Cdry+cpw*(Wtissue+W0)+cps*S0
    Eold      = Cold*T0
    Qprec     = (cpw*Rain+cps*Snowfall)*Ta
    Tolerance = 1.e-12_r8*MAX(1._r8,ABS(Eold),ABS(SW)+ABS(LW)+ABS(Qprec))

    if(Area==0._r8 .and. SW==0._r8 .and. LW==0._r8 .and. F==0._r8 &
      .and. Rain==0._r8 .and. Snowfall==0._r8)then
      call EvaluatePassiveFluxes(T0,Result,Residual)
      call PrintInfo('end '//subname)
      return
    endif

    Tlo=100._r8;Thi=400._r8
    call EvaluatePassiveFluxes(Tlo,Trial,Flo)
    call EvaluatePassiveFluxes(Thi,Trial,Fhi)
    if(Flo>0._r8 .or. Fhi<0._r8)then
      write(iulog,*)subname,': passive temperature not bracketed',T0,Ta,Flo,Fhi
      call endrun('Passive canopy energy balance has no physical temperature bracket',__LINE__)
    endif

    do iteration=1,itmax
      Tmid=0.5_r8*(Tlo+Thi)
      call EvaluatePassiveFluxes(Tmid,Result,Residual)
      if(ABS(Residual)<=Tolerance)exit
      if(Residual>0._r8)then
        Thi=Tmid
      else
        Tlo=Tmid
      endif
    enddo

    if(ABS(Residual)>Tolerance)then
      write(iulog,*)subname,': unconverged energy residual [MJ]',Residual
      call endrun('Passive canopy energy balance did not converge',__LINE__)
    endif
    call PrintInfo('end '//subname)
  contains
    subroutine EvaluatePassiveFluxes(T,State,Error)
      ! Description:
      ! Evaluate the nonlinear boundary and mass/storage terms at one trial temperature.
      real(r8), intent(in) :: T ! trial temperature [K]
      type(PassiveCanopyResult), intent(out) :: State ! trial accepted donor-limited fluxes
      real(r8), intent(out) :: Error ! stored-energy change minus accepted boundary heat [MJ]
      character(len=*), parameter :: subname='EvaluatePassiveFluxes'
      real(r8) :: Cd,Potential,DonorTemperature

      call PrintInfo('beg '//subname)
      associate(                        &
        Tout   => State%temperature,    &
        Wout   => State%water,          &
        Sout   => State%snow,           &
        Cout   => State%capacity,       &
        Rout   => State%resistance,     &
        RnOut  => State%net_radiation,  &
        Hout   => State%sensible,       &
        LEout  => State%latent,         &
        AdvOut => State%vapor_enthalpy, &
        dWout  => State%water_flux,     &
        dSout  => State%snow_flux       &
      )
      Rout      = CanopyVaporResistance(Ra,Rib,Ta,T)
      Cd        = Area/(Rout+ResistLeafSurf)
      ! As in the live canopy path, condensation remains liquid below freezing.
      ! Use liquid latent heat/enthalpy consistently; frost/freezing is not represented.
      Potential = Cd*(qa-vapsat(T))
      dWout     = MAX(Potential,-(W0+Rain))
      dSout     = MAX(MIN(0._r8,Potential-dWout),-(S0+Snowfall))
      Wout      = W0+Rain+dWout
      Sout      = S0+Snowfall+dSout
      Tout      = T
      Cout      = Cdry+cpw*(Wtissue+Wout)+cps*Sout
      Hout      = 1.25e-3_r8*Area/(Rout/BndlResistH2OToHeatRatio)*(Ta-T)
      LEout     = EvapLHTC*dWout+SublmHTC*dSout

      DonorTemperature = T
      if(dWout>0._r8)DonorTemperature=Ta
      AdvOut = cpw*DonorTemperature*dWout+cps*T*dSout
      RnOut  = SW+LW-2._r8*F*T**4
      Error  = Cout*T-Eold-(RnOut+Hout+LEout+AdvOut+Qprec)
      end associate
      call PrintInfo('end '//subname)
    end subroutine EvaluatePassiveFluxes
  end subroutine SolvePassiveCanopyEnergy
!------------------------------------------------------------------------------------------
  subroutine PassiveCanopyExchange(NZ,Cover,LiveDryCarbon)
    ! Description:
    ! Advance the passive thermal/wet-surface state at the existing PFT temperature.
    ! Include dead material once, suppress root/transpiration transfers, and retain
    ! material/temperature through burial. No biological activity flags are changed.
    integer, intent(in) :: NZ             ! named PFT slot [-]
    real(r8), intent(in) :: Cover         ! effective projected fraction of the grid cell [-]
    real(r8), intent(in) :: LiveDryCarbon ! thermally represented live dry carbon [gC/cell]
    character(len=*), parameter :: subname='PassiveCanopyExchange'
    type(PassiveCanopyResult) :: Result
    real(r8) :: DryCapacity,OldCapacity,Told,Area,SW,LW,F,Qprec
    call PrintInfo('beg '//subname)

    associate(                                                 &
      HeatSensVapAir2Canopy_col => plt_ew%HeatSensVapAir2Canopy_col, &
      HeatPrecAir2Canopy_col => plt_ew%HeatPrecAir2Canopy_col, &
      Dead               => plt_biom%StandDeadStrutElms_pft,   &
      DeadArea           => plt_morph%StandDeadSurfArea_pft, &
      StemVolume         => plt_morph%StemSpecVolume_pft, &
      Ta                 => plt_ew%TairK,                      &
      qa                 => plt_ew%VPA,                        &
      Exposed            => plt_ew%canopy_exposed_col,         &
      GridArea           => plt_site%AREA3,                    &
      NU                 => plt_site%NU,                       &
      Tissue             => plt_ew%CanopyBiomWater_pft,        &
      Water              => plt_ew%WatHeldOnCanopy_pft,        &
      Snow               => plt_ew%SnowOnCanopy_pft,           &
      Rain               => plt_ew%RainIntcptByCanopy_pft,     &
      Snowfall           => plt_ew%SnowIntcptByCanopy_pft,     &
      Tc                 => plt_ew%TKC_pft,                    &
      Ttrial             => plt_ew%TKCanopy_pft,               &
      Offset             => plt_ew%DeltaTKC_pft,               &
      Capacity           => plt_ew%VHeatCapCanopy_pft,         &
      Ra                 => plt_ew%RawIsoTCanopy2Atm_pft,      &
      Rib                => plt_ew%RIB,                        &
      Rwater             => plt_photo%RawH2OCanopy2Atm_pft,    &
      Rheat              => plt_photo%RawHeatCanopy2Atm_pft,   &
      Fraction           => plt_rad%FracPARads2Canopy_pft,     &
      AbsSW              => plt_rad%RadSWCanopyAbsorption_pft, &
      SkyLW              => plt_rad%LWRadSky_col,              &
      GroundLW           => plt_rad%LWRadGrnd_col,             &
      Emissivity         => plt_ew%EMS_Modify_Scalar_col,      &
      Rn                 => plt_rad%RadNet2Canopy_pft,         &
      LWout              => plt_rad%LWRadCanopy_pft,           &
      H                  => plt_ew%HeatXAir2PCan_pft,          &
      LE                 => plt_ew%CanopyEvapTransLHeat_pft,   &
      Storage            => plt_ew%HeatStorCanopy_pft,         &
      Evap               => plt_ew%VapXAir2Canopy_pft,         &
      LiquidFlux         => plt_ew%VapXAir2CanopyLiq_pft,      &
      SnowFlux           => plt_ew%SnoSub2AirCanopy_pft,       &
      Dew                => plt_ew%QdewCanopy_pft,             &
      Transp             => plt_ew%Transpiration_pft,          &
      RootFlux           => plt_ew%RPlantRootH2OUptk_pvr,      &
      RootStress         => plt_ew%RootH2OUptkStress_pvr,      &
      RootCd             => plt_ew%CdH2ORootxSoil_pft,         &
      AcceptedT          => Result%temperature,                &
      AcceptedWater      => Result%water,                      &
      AcceptedSnow       => Result%snow,                       &
      AcceptedCapacity   => Result%capacity,                   &
      AcceptedResistance => Result%resistance,                 &
      AcceptedRn         => Result%net_radiation,              &
      AcceptedH          => Result%sensible,                   &
      AcceptedLE         => Result%latent,                     &
      AcceptedAdv        => Result%vapor_enthalpy,             &
      AcceptedLiquidFlux => Result%water_flux,                 &
      AcceptedSnowFlux   => Result%snow_flux                   &
    )
      Told=Tc(NZ)
      ! A cold-start/legacy dead-only slot may have no previously initialized temperature.
      if(Told==spval .or. Told<=0._r8)Told=Ta

      DryCapacity = cpo*gOC_to_m3_OM(LiveDryCarbon+ &
        CanopyStemThermalCarbon(Dead(ielmc,NZ),DeadArea(NZ),StemVolume(NZ)))
      OldCapacity = DryCapacity+cpw*(Tissue(NZ)+Water(NZ))+cps*Snow(NZ)
      Area        = 0._r8
      SW          = 0._r8
      LW          = 0._r8
      F           = 0._r8
      if(Exposed .and. Cover>0._r8)then
        Area = Cover*GridArea(NU)
        SW   = AbsSW(NZ)
        LW   = (SkyLW+GroundLW)*Fraction(NZ)
        F    = EMMC*stefboltz_const*Fraction(NZ)*GridArea(NU)*Emissivity
      endif

      call SolvePassiveCanopyEnergy(Told,Ta,qa,DryCapacity,Tissue(NZ),Water(NZ),Snow(NZ),Rain(NZ),Snowfall(NZ), &
        Area,SW,LW,F,Ra(NZ),Rib,Result)

      Qprec              = (cpw*Rain(NZ)+cps*Snowfall(NZ))*Ta
      HeatSensVapAir2Canopy_col=HeatSensVapAir2Canopy_col+AcceptedAdv
      HeatPrecAir2Canopy_col=HeatPrecAir2Canopy_col+Qprec
      Water(NZ)          = AcceptedWater
      Snow(NZ)           = AcceptedSnow
      Tc(NZ)             = AcceptedT
      Ttrial(NZ)         = Tc(NZ)
      Offset(NZ)         = Tc(NZ)-Ta
      Capacity(NZ)       = AcceptedCapacity
      Rwater(NZ)         = AcceptedResistance
      Rheat(NZ)          = Rwater(NZ)/BndlResistH2OToHeatRatio
      Rn(NZ)             = AcceptedRn
      LWout(NZ)          = F*Tc(NZ)**4
      H(NZ)              = AcceptedH
      LE(NZ)             = AcceptedLE
      Storage(NZ)        = OldCapacity*Told-Capacity(NZ)*Tc(NZ)+AcceptedAdv+Qprec
      LiquidFlux(NZ)     = AcceptedLiquidFlux
      SnowFlux(NZ)       = AcceptedSnowFlux
      Evap(NZ)           = LiquidFlux(NZ)+SnowFlux(NZ)
      Dew(NZ)            = MAX(0._r8,LiquidFlux(NZ))
      Transp(NZ)         = 0._r8
      RootCd(NZ)         = 0._r8
      RootFlux(:,:,NZ)   = 0._r8
      RootStress(:,:,NZ) = 0._r8
    end associate
    call PrintInfo('end '//subname)
  end subroutine PassiveCanopyExchange

!----------------------------------------------------------------------------------------------------
  subroutine PrepNutrientUptake(yearIJ,TotalSoilPSIMPa_vr,AllRootC_vr,AirMicPore4Fill_vr,WatAvail4Uptake_vr)
!
!     prepare for uptake calculation
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  real(r8), intent(out) :: TotalSoilPSIMPa_vr(JZ1)   !total soil matric pressure after removing elevation adjustment, [MPa]
  real(r8), intent(out) :: AllRootC_vr(JZ1)          !total root C profile, [gC d-2]
  real(r8), intent(out) :: AirMicPore4Fill_vr(JZ1)
  real(r8), intent(out) :: WatAvail4Uptake_vr(JZ1)

  character(len=*), parameter :: subname='PrepNutrientUptake'
  integer :: NZ, L, N

  associate(                                                          &
    ALT                        => plt_site%ALT                       ,& !input  :altitude of the grid cell, [m]
    CanopyLeafArea_pft         => plt_morph%CanopyLeafArea_pft       ,& !input  :plant canopy leaf area, [m2 d-2]
    ElvAdjstedSoilH2OPSIMPa_vr => plt_ew%ElvAdjstedSoilH2OPSIMPa_vr  ,& !input  :soil micropore total water potential, [MPa]
    MaxNumRootLays             => plt_site%MaxNumRootLays            ,& !input  :maximum root layer number,[-]
    Myco_pft                   => plt_morph%Myco_pft                 ,& !input  :mycorrhizal type (no or yes),[-]
    NK                         => plt_site%NK                        ,& !input  :current hydrologically active layer, [-]
    NP                         => plt_site%NP                        ,& !input  :current number of plant species,[-]
    NP0                        => plt_site%NP0                       ,& !input  :intitial number of plant species,[-]
    NU                         => plt_site%NU                        ,& !input  :current soil surface layer number, [-]
    PopuRootMycoC_pvr          => plt_biom%PopuRootMycoC_pvr         ,& !input  :root layer C, [gC d-2]
    SoilBulkDensity_vr         => plt_soilchem%SoilBulkDensity_vr    ,& !input  :soil bulk density, [Mg m-3]
    SoilWatAirDry_vr           => plt_soilchem%SoilWatAirDry_vr      ,& !input  :air-dry water content, [m3 m-3]
    VLMicP_vr                  => plt_soilchem%VLMicP_vr             ,& !input  :total volume in micropores, [m3 d-2]
    VLSoilMicP_vr              => plt_soilchem%VLSoilMicP_vr         ,& !input  :total micropore volume in layer, [m3 d-2]
    VLWatMicPM_vr              => plt_site%VLWatMicPM_vr             ,& !input  :soil micropore water content, [m3 d-2]
    VLWatMicP_vr               => plt_soilchem%VLWatMicP_vr          ,& !input  :soil micropore water content, [m3 d-2]
    VLiceMicP_vr               => plt_soilchem%VLiceMicP_vr          ,& !input  :soil micropore ice content, [m3 d-2]
    ZERO                       => plt_site%ZERO                      ,& !input  :threshold zero for numerical stability, [-]
    LWRadCanopy_pft            => plt_rad%LWRadCanopy_pft            ,& !output :canopy longwave radiation, [MJ d-2 h-1]
    RadNet2Canopy_pft          => plt_rad%RadNet2Canopy_pft           & !output :canopy net radiation, [MJ d-2 h-1]
  )
  !
  !     RESET TOTAL UPTAKE ARRAYS
  !
  !
  call PrintInfo('beg '//subname)
  !
  ! NPH is the last iteration from solving for soil heat-moisture hydrothermal dynamics
  D9000: DO L=NU,NK
    !remove elevation dependence
    TotalSoilPSIMPa_vr(L)=ElvAdjstedSoilH2OPSIMPa_vr(L)-mGravAccelerat*ALT
    IF(SoilBulkDensity_vr(L).GT.ZERO)THEN
      WatAvail4Uptake_vr(L) = VLWatMicPM_vr(NPH,L)-SoilWatAirDry_vr(L)*VLSoilMicP_vr(L)     !maximum amount of water for uptake
      AirMicPore4Fill_vr(L) = AZMAX1(VLMicP_vr(L)-VLWatMicP_vr(L)-VLiceMicP_vr(L))   !air volume in soil
    ELSE
      WatAvail4Uptake_vr(L) = VLWatMicPM_vr(NPH,L)
      AirMicPore4Fill_vr(L) = 0.0_r8
    ENDIF
  ENDDO D9000

  DO L=NU,NK
    AllRootC_vr(L)=0.0_r8
    D9005: DO NZ=1,NP
      DO  N=1,Myco_pft(NZ)
      !turn off to include dead plants as well
!     IF(IsPlantActive_pft(NZ).EQ.iTrue .AND. PlantPopuLive_pft(NZ).GT.0.0)THEN
        AllRootC_vr(L)=AllRootC_vr(L)+AZMAX1(PopuRootMycoC_pvr(N,L,NZ))
!     ENDIF
      enddo
    ENDDO D9005
  ENDDO

  call PrintInfo('end '//subname)
  end associate
  end subroutine PrepNutrientUptake

!----------------------------------------------------------------------------------------------------
  subroutine UpdateCanopyProperty(I,J,NZ)
  ! Description:
  ! Update canopy resistance and trial temperature using exposed live/dead geometry.
  implicit none
  integer, intent(in) :: I,J,NZ

  character(len=*), parameter :: subname='UpdateCanopyProperty'
  real(r8) :: ALFZ       !attenuation coefficient
  real(r8) :: TFRADP     !total fraction of PAR on canopy
  real(r8) :: RawCanopy2SinkZ_pft(JP1) !isothermal aerodynamic resistance between current pft top to the sink height
  real(r8) :: ExchangeHeight(JP1),ExchangeArea(JP1)
  integer :: NB,K,L,N,NZZ

  call PrintInfo('beg '//subname)
  associate(                                                               &
    RawIsoTAtm2CanopySinkZ_col    => plt_ew%RawIsoTAtm2CanopySinkZ_col    ,& !input  :isothermal aerodynamic resistance between canopy top and wind ref height in atmosphere, [h m-1]
    RawIsoTSurf2CanopyHScal_col   => plt_ew%RawIsoTSurf2CanopyHScal_col   ,& !input  :scalar for isothermal aerodynamic resistance between zero-sink height and ground surface, [h m-1]
    CanopyHeight_col              => plt_morph%CanopyHeight_col           ,& !input  :canopy height , [m]
    CanopyHeightLive_pft          => plt_morph%CanopyHeightLive_pft       ,& !input  :canopy height, [m]
    DeltaTKC_pft                  => plt_ew%DeltaTKC_pft                  ,& !input  :change in canopy temperature, [K]
    FracPARads2Canopy_pft         => plt_rad%FracPARads2Canopy_pft        ,& !input  :fraction of incoming PAR absorbed by canopy, [-]
    KoppenClimZone                => plt_site%KoppenClimZone              ,& !input  :Koppen climate zone for the grid,[-]
    LeafAreaZsec_brch             => plt_morph%LeafAreaZsec_brch          ,& !input  :leaf surface area, [m2 d-2]
    LeafArea_node                 => plt_morph%LeafArea_node              ,& !input  :leaf area, [m2 d-2]
    LeafProteinC_node             => plt_biom%LeafProteinC_node           ,& !input  :layer leaf protein N, [g d-2]
    LeafStalkAreaAct_pft          => plt_morph%LeafStalkAreaAct_pft       ,& !input  :plant leaf+stem/stalk area, [m2 d-2]
    DeadArea                     => plt_morph%StandDeadSurfAreaAct_pft, &
    DeadHeight                   => plt_morph%CanopyHeightDead_pft, &
    NP                            => plt_site%NP                          ,& !input  :current number of plant species,[-]
    NumOfBranches_pft             => plt_morph%NumOfBranches_pft          ,& !input  :number of branches,[-]
    RoughnessLength               => plt_ew%RoughnessLength               ,& !input  :canopy surface roughness height, [m]
    TairK                         => plt_ew%TairK                         ,& !input  :air temperature, [K]
    ZERO                          => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]
    ZERO4Groth_pft                => plt_biom%ZERO4Groth_pft              ,& !input  :threshold zero for plang growth calculation, [-]
    ZeroPlaneDisplacem_col        => plt_ew%ZeroPlaneDisplacem_col        ,& !input  :zero plane displacement height, [m]
    KMinNumLeaf4GroAlloc_brch     => plt_morph%KMinNumLeaf4GroAlloc_brch  ,& !output :NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION,[-]
    RawIsoTCanopy2Atm_pft         => plt_ew%RawIsoTCanopy2Atm_pft         ,& !output :!isothermal aerodynamic resistance from canopy to free atmosphere (where wind speed is measured) 
    TKCanopy_pft                  => plt_ew%TKCanopy_pft                   & !output :canopy temperature, [K]
  )
!
!     APPLY CLUMPING FACTOR TO LEAF SURFACE AREA DEFINED BY
!     INCLINATION N, LAYER L, NODE K, BRANCH NB, SPECIES NZ,
!     N-S POSITION NY, E-W POSITION NX(AZIMUTH M ASSUMED UNIFORM)
!

!  D500: DO NB=1,NumOfBranches_pft(NZ)
!    D550: DO K=1,MaxNodesPerBranch1
!
!     NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
!
!     LeafArea_node=leaf area
!     LeafProteinC_node=leaf protein content
!
!      IF(LeafArea_node(K,NB,NZ).GT.ZERO4Groth_pft(NZ) .AND. LeafProteinC_node(K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
!        KMinNumLeaf4GroAlloc_brch(NB,NZ)=K
!      ENDIF

!    ENDDO D550
!  ENDDO D500
  !
  !     DIFFUSIVE RESISTANCE OF OTHER TALLER CANOPIES TO HEAT AND VAPOR
  !     TRANSFER OF CURRENT CANOPY ADDED TO BOUNDARY LAYER RESISTANCE
  !     OF TALLEST CANOPY CALCULATED IN 'HOUR1'
  !
  !     KoppenClimZone=Koppen climate zone
  !     ZC,ZT,RoughnessLength=PFT canopy,grid biome,surface roughness height
  !     FracPARads2Canopy_pft=fraction of radiation received by each PFT canopy
  !
  ExchangeHeight(1:NP)=CanopyHeightLive_pft(1:NP)
  ExchangeArea(1:NP)=LeafStalkAreaAct_pft(1:NP)+DeadArea(1:NP)
  do NZZ=1,NP
    if(DeadArea(NZZ)>0._r8 .and. ExchangeArea(NZZ)>0._r8) &
      ExchangeHeight(NZZ)=(LeafStalkAreaAct_pft(NZZ)*CanopyHeightLive_pft(NZZ) &
        +DeadArea(NZZ)*DeadHeight(NZZ))/ExchangeArea(NZZ)
  enddo
  IF(ExchangeArea(NZ).GT.0.0_r8)THEN
    IF(KoppenClimZone.GE.0)THEN
      TFRADP=0.0_r8
      D700: DO NZZ=1,NP
        IF(ExchangeHeight(NZZ).GT.ExchangeHeight(NZ)+RoughnessLength)THEN
          TFRADP=TFRADP+FracPARads2Canopy_pft(NZZ)
        ENDIF
      ENDDO D700
      !attenuation coefficient
      ALFZ=2.0_r8*TFRADP
      !
      !RACX:=minimum boundary layer resistances of canopy [h/m]
      !
      IF(RawIsoTSurf2CanopyHScal_col.GT.ZERO .AND. CanopyHeight_col.GT.ZERO .AND. ALFZ.GT.ZERO)THEN
        !eq.(25) from Choudhury and Monteith (1988)
        !resistance between current pft canopy height and sink height (=zero plane displacement + roughness length)
        RawCanopy2SinkZ_pft(NZ)=AMIN1(RACX,AZMAX1(CanopyHeight_col*EXP(ALFZ) &
          /(ALFZ/RawIsoTSurf2CanopyHScal_col)*(EXP(-ALFZ*ExchangeHeight(NZ)/CanopyHeight_col) &
          -EXP(-ALFZ*(ZeroPlaneDisplacem_col+RoughnessLength)/CanopyHeight_col))))
      ELSE
        RawCanopy2SinkZ_pft(NZ)=0.0_r8
      ENDIF
    ELSE
      RawCanopy2SinkZ_pft(NZ)=0.0_r8
    ENDIF
  ELSE
    RawCanopy2SinkZ_pft(NZ)=RACX
  ENDIF
  !isothermal aerodynamic resistance from canopy to free atmosphere (where wind speed is measured) 
  RawIsoTCanopy2Atm_pft(NZ)=RawIsoTAtm2CanopySinkZ_col+RawCanopy2SinkZ_pft(NZ)
  !
  !     INITIALIZE CANOPY TEMPERATURE WITH CURRENT AIR TEMPERATURE AND
  !     LAST HOUR'S CANOPY-AIR TEMPERATURE DIFFERENCE, AND CALL A
  !     subroutine TO CALCULATE MINIMUM CANOPY STOMATAL RESISTANCE
  !     FOR SUBSEQUENT USE IN ENERGY EXCHANGE CALCULATIONS
  !
  ! TKCanopy_pft = initial estimate of canopy temperature TKC
  ! TairK        = current air temperature
  ! DeltaTKC_pft = TKC-TairK from previous hour
  !
  TKCanopy_pft(NZ)=TairK+DeltaTKC_pft(NZ)

  end associate
  call PrintInfo('end '//subname)
  end subroutine UpdateCanopyProperty

!----------------------------------------------------------------------------------------------------
  subroutine UpdateRootProperty(yearIJ,NZ,RadialMeanLen_rvr,FineRootRadius_rvr,AllRootC_vr,&
    FracPRoot4Uptake_pvr,FracMinRoot4Uptake_rpvr,FracSoilLBy1stRoots_pvr,RootEffLen4Absorption_pvr)
  !
  ! Description:
  ! Update root characterization for pft NZ
  !
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer , intent(in) :: NZ
  real(r8), intent(in) :: AllRootC_vr(JZ1)
  real(r8), intent(out) :: RadialMeanLen_rvr(pltpar%jroots,JZ1)
  real(r8), intent(out) :: FineRootRadius_rvr(pltpar%jroots,JZ1)
  real(r8), intent(out) :: FracPRoot4Uptake_pvr(pltpar%jroots,JZ1,JP1)
  real(r8), intent(out) :: FracMinRoot4Uptake_rpvr(pltpar%jroots,JZ1,JP1)
  real(r8), intent(out) :: FracSoilLBy1stRoots_pvr(JZ1,JP1)                   !fraction of soil layer occupied by primary roots
  real(r8), intent(out) :: RootEffLen4Absorption_pvr(pltpar%jroots,JZ1)       !ratio of lateral root area to root radius, [m]

  character(len=*), parameter :: subname='UpdateRootProperty'
  real(r8) :: RootDepZ         !primary root depth
  real(r8) :: RootLocDepz      !Primary Root occupied thickness in current layer
  integer :: N,L,NR
  
  associate(                                                          &
    CumSoilThickness_vr        => plt_site%CumSoilThickness_vr       ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]
    DLYR3                      => plt_site%DLYR3                     ,& !input  :vertical thickness of soil layer, [m]
    FracSoiAsMicP_vr           => plt_site%FracSoiAsMicP_vr          ,& !input  :micropore fraction, [-]
    HypocotHeight_pft          => plt_morph%HypocotHeight_pft        ,& !input  :cotyledon height, [m]
    MaxSoilLays4Root_pft       => plt_morph%MaxSoilLays4Root_pft     ,& !input  :maximum soil layer number for all root axes,[-]
    Myco_pft                   => plt_morph%Myco_pft                 ,& !input  :mycorrhizal type (no or yes),[-]
    NU                         => plt_site%NU                        ,& !input  :current soil surface layer number, [-]
    NumStructuralRootAxes_pft  => plt_morph%NumStructuralRootAxes_pft,& !input  :number of structural root axes,[-]
    Root2ndAxialResist_pft     => plt_morph%Root2ndAxialResist_pft   ,& !input  :root axial resistivity, [MPa h m-4]    
    RootRadialResist_pft       => plt_morph%RootRadialResist_pft     ,& !input  :root radial resistivity, [MPa h m-1]        
    PlantPopuLive_pft          => plt_site%PlantPopuLive_pft         ,& !input  :plant population, [d-2]
    PopuRootMycoC_pvr          => plt_biom% PopuRootMycoC_pvr        ,& !input  :root layer C, [gC d-2]
    Root1stDepz_raxes          => plt_morph%Root1stDepz_raxes        ,& !input  :root layer depth, [m]
    Root2ndRadius_rpvr         => plt_morph%Root2ndRadius_rpvr       ,& !input  :root layer diameter secondary axes, [m]    
    Root2ndMaxRadius1_pft      => plt_morph%Root2ndMaxRadius1_pft    ,& !input  :root diameter secondary axes, [m]
    Root1stRadius_pvr          => plt_morph%Root1stRadius_pvr        ,& !input  :root layer diameter primary axes, [m]    
    Root1stMaxRadius_pft       => plt_morph%Root1stMaxRadius_pft      ,& !output :maximum radius of primary roots, [m]
    Root2ndMaxRadius_pft       => plt_morph%Root2ndMaxRadius_pft     ,& !input  :maximum radius of secondary roots, [m]
    RootLenDensPerPlant_pvr    => plt_morph%RootLenDensPerPlant_pvr  ,& !input  :root layer length density, [m m-3]
    RootAbsorbLenPerPlant_pvr  => plt_morph%RootAbsorbLenPerPlant_pvr,& !input  :total absorptive root length per plant in layer, [m p-1]
    RootPorosity_pft           => plt_morph%RootPorosity_pft         ,& !input  :root porosity, [m3 m-3]
    RootVH2O_pvr               => plt_morph%RootVH2O_pvr             ,& !input  :root layer volume water, [m2 d-2]
    SeedDepth_pft              => plt_morph%SeedDepth_pft            ,& !input  :seeding depth, [m]
    ZERO                       => plt_site%ZERO                      ,& !input  :threshold zero for numerical stability, [-]
    ZEROS                      => plt_site%ZEROS                      & !input  :threshold zero for numerical stability,[-]
  )
  call PrintInfo('beg '//subname)
  if(ldo_sp_mode)then
    !no mycorrhizae considered for prescribed phenology mode
    N=ipltroot
    DO  L=NU,MaxSoilLays4Root_pft(NZ)
      FracSoilLBy1stRoots_pvr(L,NZ)=1.0_r8
    ENDDO
    if(ats_cpl_mode)then
        !the following number should be revised
        !Root2ndAxialResist_pft(N,NZ) = 500.0_r8
        !RootRadialResist_pft(N,NZ) = 400000.0_r8
        DO N=1,Myco_pft(NZ)
            !test vals
            !Root2ndMaxRadius_pft(N,NZ) = 0.0002_r8
            DO L=NU,MaxSoilLays4Root_pft(NZ)
                Root2ndradius_rpvr(N,L,NZ) = Root2ndMaxRadius_pft(N,NZ)
                Root1stRadius_pvr(N,L,NZ) = Root1stMaxRadius_pft(N,NZ)
            ENDDO
        ENDDO
    endif
  else
    !obtain plant rooting depth
    RootDepZ=0.0_r8
    D2005: DO NR=1,NumStructuralRootAxes_pft(NZ)
      RootDepZ=AMAX1(RootDepZ,Root1stDepz_raxes(NR,NZ))
    ENDDO D2005

    DO  L=NU,MaxSoilLays4Root_pft(NZ)
      IF(L.EQ.NU)THEN
        FracSoilLBy1stRoots_pvr(L,NZ)=1.0_r8
      ELSE
        IF(DLYR3(L).GT.ZERO)THEN
          RootLocDepz = AZMAX1(RootDepZ-CumSoilThickness_vr(L-1)) !root covered depth
          RootLocDepz = AZMAX1(AMIN1(DLYR3(L),RootLocDepz)-AZMAX1(SeedDepth_pft(NZ)-CumSoilThickness_vr(L-1)-HypocotHeight_pft(NZ)))          
          FracSoilLBy1stRoots_pvr(L,NZ) = RootLocDepz/DLYR3(L)
        ELSE
          FracSoilLBy1stRoots_pvr(L,NZ)=0.0_r8
        ENDIF
      ENDIF
    ENDDO
  endif
  
  D2000: DO N=1,Myco_pft(NZ)

    DO  L=NU,MaxSoilLays4Root_pft(NZ)      
      IF(AllRootC_vr(L).GT.ZEROS)THEN
        FracPRoot4Uptake_pvr(N,L,NZ)=PopuRootMycoC_pvr(N,L,NZ)/AllRootC_vr(L)
      ELSE
        FracPRoot4Uptake_pvr(N,L,NZ)=1.0_r8
      ENDIF
      FracMinRoot4Uptake_rpvr(N,L,NZ)=FMN*AZMAX1(FracPRoot4Uptake_pvr(N,L,NZ))
      
      IF(RootAbsorbLenPerPlant_pvr(N,L,NZ).GT.ZERO .AND. FracSoilLBy1stRoots_pvr(L,NZ).GT.ZERO .and. .not.ldo_sp_mode)THEN        
        FineRootRadius_rvr(N,L)=AMAX1(Root2ndMaxRadius1_pft(N,NZ),SQRT(RootVH2O_pvr(N,L,NZ) &
          /((1.0_r8-RootPorosity_pft(N,NZ))*PICON*RootAbsorbLenPerPlant_pvr(N,L,NZ)*PlantPopuLive_pft(NZ))))
        
        !cylinder, pi*r^2*L=vol
        RadialMeanLen_rvr(N,L)=AMAX1(1.001_r8*FineRootRadius_rvr(N,L),&
          SQRT(FracSoilLBy1stRoots_pvr(L,NZ)*FracSoiAsMicP_vr(L)/(PICON*RootAbsorbLenPerPlant_pvr(N,L,NZ))))

        RootEffLen4Absorption_pvr(N,L)=TwoPiCON*RootAbsorbLenPerPlant_pvr(N,L,NZ)/FracSoilLBy1stRoots_pvr(L,NZ)

      ELSE
        FineRootRadius_rvr(N,L)        = Root2ndMaxRadius_pft(N,NZ)
        RadialMeanLen_rvr(N,L)         = 1.001_r8*FineRootRadius_rvr(N,L)
        RootEffLen4Absorption_pvr(N,L) = TwoPiCON*RootAbsorbLenPerPlant_pvr(N,L,NZ)
      ENDIF
      
    enddo
  ENDDO D2000
  call PrintInfo('end '//subname)
  end associate
  end subroutine UpdateRootProperty

!----------------------------------------------------------------------------------------------------
  subroutine HandlingDivergence(yearIJ,NZ,TotalSoilPSIMPa_vr,DIFF,FDMP,cumPRootH2OUptake,HeatEvapSens, &
    CumPlantHeatLoss2Soil)
  ! Description:
  ! Reject unaccepted root exchange and prepare the conservative passive canopy
  ! fallback. Report routine fallback use only when verbose logging is enabled.
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer  , intent(in) :: NZ
  REAL(R8) , INTENT(IN) :: TotalSoilPSIMPa_vr(JZ1)
  real(r8) , intent(in) :: DIFF   !divergence check
  real(r8), intent(inout):: FDMP
  real(r8), intent(inout) :: cumPRootH2OUptake,HeatEvapSens
  real(r8), intent(inout) :: CumPlantHeatLoss2Soil
  character(len=*), parameter :: subname='HandlingDivergence'
  real(r8) :: APSILT
  real(r8) :: CCPOLT
  real(r8) :: FTHRM
  real(r8) :: OSWT,Stomata_Stress

  integer :: N,L
! begin_execution
  call PrintInfo('beg '//subname)
  associate(                                                               &
    NY                          => plt_site%NY                            ,& ! input: column row index [-]
    NX                          => plt_site%NX                            ,& ! input: column index [-]
    AREA3                       => plt_site%AREA3                         ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    OrganOsmoPsi0pt_pft         => plt_ew%OrganOsmoPsi0pt_pft             ,& !input  :Organ osmotic potential when canopy water potential = 0 MPa, [MPa]
    CanopyNonstOMConc_pft       => plt_biom%CanopyNonstOMConc_pft         ,& !input  :canopy nonstructural C+N+P concentration, [g gC-1]
    FracPARads2Canopy_pft       => plt_rad%FracPARads2Canopy_pft          ,& !input  :fraction of incoming PAR absorbed by canopy, [-]
    MaxSoilLays4Root_pft        => plt_morph%MaxSoilLays4Root_pft         ,& !input  :maximum soil layer number for all root axes,[-]
    EMS_Modify_Scalar_col       => plt_ew%EMS_Modify_Scalar_col           ,& !input  :canopy longwave radiation emissivity scalar      
    Myco_pft                    => plt_morph%Myco_pft                     ,& !input  :mycorrhizal type (no or yes),[-]
    NGTopRootLayer_pft          => plt_morph%NGTopRootLayer_pft           ,& !input  :soil layer at planting depth, [-]
    NU                          => plt_site%NU                            ,& !input  :current soil surface layer number, [-]
    PSICanopyOsmo_pft           => plt_ew%PSICanopyOsmo_pft               ,& !input  :canopy osmotic water potential, [Mpa]
    PSICanopyTurg_pft           => plt_ew%PSICanopyTurg_pft               ,& !input  :plant canopy turgor water potential, [MPa]
    PSIRootOSMO_vr              => plt_ew%PSIRootOSMO_vr                  ,& !input  :root osmotic water potential, [Mpa]
    PSIRootTurg_vr              => plt_ew%PSIRootTurg_vr                  ,& !input  :root turgor water potential, [Mpa]
    RawIsoTCanopy2Atm_pft       => plt_ew%RawIsoTCanopy2Atm_pft           ,& !input  :isothermal aerodynamic resistance from canopy to free atmosphere (where wind speed is measured) 
    RootNonstOMConc_pvr         => plt_biom%RootNonstOMConc_pvr           ,& !input  :root/mycorrhizal nonstructural C+N+P concentration, [g gC-1]
    ShootElms_pft               => plt_biom%ShootElms_pft                 ,& !input  :canopy shoot structural chemical element mass, [g d-2]
    TKS_vr                      => plt_ew%TKS_vr                          ,& !input  :mean annual soil temperature, [K]
    TairK                       => plt_ew%TairK                           ,& !input  :air temperature, [K]
    iYearCurrent                => plt_site%iYearCurrent                  ,& !input  :current year,[-]
    DeltaTKC_pft                => plt_ew%DeltaTKC_pft                    ,& !inoput :change in canopy temperature, [K]
    PSICanopy_pft               => plt_ew%PSICanopy_pft                   ,& !output :canopy total water potential, [Mpa]
    PSIRoot_pvr                 => plt_ew%PSIRoot_pvr                     ,& !output :root total water potential, [Mpa]
    TKC_pft                     => plt_ew%TKC_pft                         ,& !output :canopy temperature, [K]
    RPlantRootH2OUptk_pvr       => plt_ew%RPlantRootH2OUptk_pvr           ,& !output :whole population root water uptake, [m3 d-2 h-1]
    RootH2OUptkStress_pvr       => plt_ew%RootH2OUptkStress_pvr           ,& !output :pontential root water uptake rate, [m3 d-2 h-1]
    CanPStomaResistH2O_pft      => plt_photo%CanPStomaResistH2O_pft       ,& !output :canopy stomatal resistance, [h m-1]
    RawH2OCanopy2Atm_pft        => plt_photo%RawH2OCanopy2Atm_pft         ,& !output :canopy boundary layer resistance, [h m-1]
    RawHeatCanopy2Atm_pft        => plt_photo%RawHeatCanopy2Atm_pft         ,& !output :canopy sensible-heat boundary layer resistance, [h m-1]
    LWRadCanopy_pft             => plt_rad%LWRadCanopy_pft                ,& !output :canopy longwave radiation, [MJ d-2 h-1]
    VHeatCapCanopy_pft          => plt_ew%VHeatCapCanopy_pft               & !output :canopy heat capacity, [MJ d-2 K-1]
  )
    if(do_debug())write(iulog,'(a,6i6)') &
      'Using passive canopy fallback; year/day/hour/row/column/PFT: ',iYearCurrent,yearIJ%I,yearIJ%J,NY,NX,NZ

    IF(DIFF.GE.DIFFY)THEN
      cumPRootH2OUptake=0._r8;HeatEvapSens=0._r8
      CumPlantHeatLoss2Soil=0._r8
      plt_rad%RadNet2Canopy_pft(NZ)       = 0.0_r8
      plt_ew%CanopyEvapTransLHeat_pft(NZ) = 0.0_r8
      plt_ew%HeatXAir2PCan_pft(NZ)        = 0.0_r8
      plt_ew%HeatStorCanopy_pft(NZ)       = 0.0_r8
      plt_ew%VapXAir2Canopy_pft(NZ)       = 0.0_r8
      plt_ew%Transpiration_pft(NZ)        = 0.0_r8
      plt_ew%SnoSub2AirCanopy_pft(NZ)     = 0.0_r8
      plt_ew%VapXAir2CanopyLiq_pft(NZ)    = 0.0_R8
      plt_ew%QdewCanopy_pft(NZ)           = 0.0_R8

      TKC_pft(NZ)         = TairK+DeltaTKC_pft(NZ)
      FTHRM               = EMMC*stefboltz_const*FracPARads2Canopy_pft(NZ)*AREA3(NU)*EMS_Modify_Scalar_col
      LWRadCanopy_pft(NZ) = FTHRM*TKC_pft(NZ)**4
      PSICanopy_pft(NZ)   = TotalSoilPSIMPa_vr(NGTopRootLayer_pft(NZ))

      CCPOLT  = CanopyNonstOMConc_pft(NZ)

      CALL update_osmo_turg_pressure(PSICanopy_pft(NZ),CCPOLT,OrganOsmoPsi0pt_pft(NZ),TKC_pft(NZ), &
        PSICanopyOsmo_pft(NZ),PSICanopyTurg_pft(NZ),FDMP)

      RawH2OCanopy2Atm_pft(NZ)      = RawIsoTCanopy2Atm_pft(NZ)
      RawHeatCanopy2Atm_pft(NZ) = RawH2OCanopy2Atm_pft(NZ)/BndlResistH2OToHeatRatio
      VHeatCapCanopy_pft(NZ)     = cpw*(ShootElms_pft(ielmc,NZ)*10.0E-06_r8)
      DeltaTKC_pft(NZ)           = 0.0_r8
      CanPStomaResistH2O_pft(NZ) = CalcStomataResist4H2O(NZ)
        
      D4290: DO N=1,Myco_pft(NZ)
        DO  L=NU,MaxSoilLays4Root_pft(NZ)
          PSIRoot_pvr(N,L,NZ) = TotalSoilPSIMPa_vr(L)
          if(ldo_sp_mode)then
            CCPOLT=0.1_r8
          else
            CCPOLT   = RootNonstOMConc_pvr(N,L,NZ)
          endif
          CALL update_osmo_turg_pressure(PSIRoot_pvr(N,L,NZ),CCPOLT,OrganOsmoPsi0pt_pft(NZ),TKS_vr(L),&
            PSIRootOSMO_vr(N,L,NZ),PSIRootTurg_vr(N,L,NZ))

          RPlantRootH2OUptk_pvr(N,L,NZ)=0.0_r8
          RootH2OUptkStress_pvr(N,L,NZ) =0._r8
      enddo
      ENDDO D4290
    ENDIF

  end associate
  call PrintInfo('end '//subname)
  end subroutine HandlingDivergence

!----------------------------------------------------------------------------------------------------
  function CanopyEnergyH2OIter_func(yearIJ,NZ,FracGridPFTCovered,CanopyMassC4H2OStorage,TotalSoilPSIMPa_vr,  &
    CdHeatCanopyAir,DIFF,cumPRootH2OUptake,CumPlantHeatLoss2Soil,HeatEvapSens,FDMP,                       &
    SoilRootResist4H2O_pvr,FracPRoot4Uptake_pvr,AirMicPore4Fill_vr,WatAvail4Uptake_vr,TKCX,CdH2ORootxSoil,&
    VHeatCapCanopyPrev_pft,PrecpHeatbyCanopy,PSIGravCanopyHeight,SoiLayerHasRoot_rvr,HydraulicPathAvailable,LiveAreaFraction) result(NN)
  ! Description:
  !The relationship between plant root water uptake and canopy water storage is approximated using the water tower concept. That is
  !Root at each soil layer interacts directly with the canopy water reservoir
  !Negative flux represents water taken up from soil
  !Positive flux represents water release into soil
  !Two processses together support hydraulic redistribution and transpiration.
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer  , intent(in) :: NZ
  real(r8) , intent(in) :: FracGridPFTCovered
  real(r8) , intent(in) :: CanopyMassC4H2OStorage
  real(r8) , intent(in) :: TotalSoilPSIMPa_vr(JZ1)
  real(r8) , intent(in) :: SoilRootResist4H2O_pvr(pltpar%jroots,JZ1)
  real(r8) , intent(in) :: FracPRoot4Uptake_pvr(pltpar%jroots,JZ1,JP1)
  real(r8) , intent(in) :: AirMicPore4Fill_vr(JZ1),WatAvail4Uptake_vr(JZ1)
  real(r8) , intent(in) :: TKCX
  real(r8) , intent(in) :: CdH2ORootxSoil            !total conductance, [m H2O d-2 h-1 MPa-1]
  real(r8) , intent(in) :: VHeatCapCanopyPrev_pft    !canopy heat capacity at previous time step, MJ/K
  real(r8) , intent(in) :: PrecpHeatbyCanopy         !heat added to canopy by precipitation [MJ]
  real(r8) , intent(in) :: PSIGravCanopyHeight       !Graviational water potential at effective canopy height [MPa]
  logical  , intent(in) :: SoiLayerHasRoot_rvr(pltpar%jroots,JZ1)
  logical  , intent(in) :: HydraulicPathAvailable
  real(r8) , intent(out):: CdHeatCanopyAir           !area-scaled conductance for canopy-air sensible heat flux
  real(r8) , intent(out):: DIFF
  real(r8) , intent(out):: cumPRootH2OUptake         !total root water uptake from soil  (<0 addition to canopy)
  real(r8) , intent(out):: CumPlantHeatLoss2Soil     !total heat associated to root water uptake from soil  (<0 addition to canopy)
  real(r8) , intent(out):: HeatEvapSens   !sensible heat due to evaporation/condensation  [MJ]
  real(r8) , intent(out):: FDMP
  real(r8) :: APSILT
  real(r8) :: CCPOLT
  real(r8) :: DTHS1
  real(r8) :: DIFFZ,DIFFU,DPSI
  real(r8) :: EX,PTransPre
  real(r8) :: FTHRM
  real(r8) :: LWRad2Canopy   !long wave radiation (from sky + ground) on canopy [MJ]
  real(r8) :: HeatAdd2Can
  real(r8) :: EffPFTGridArea4H2O,EffPFTGridArea4Heat,PSICanPPre
  real(r8) :: CdVapCanopyAir
  real(r8) :: PSILC               !canopy height adjusted canopy water potential, [MPa]
  real(r8) :: RSSUX,RSSU,RA1,RSSZ
  real(r8) :: TKC1,TKCY
  real(r8) :: HeatResidual ! accepted canopy energy-budget residual [MJ/cell]
  real(r8), parameter :: CanopyHeatTolerance=1.e-9_r8 ! accepted hourly energy closure [MJ/cell]
  real(r8) :: cumPRootH2OUptakePre
  real(r8) :: SymplasmicWaterPrev,SymplasmicWater  !previous, current biomass-bounded water
  real(r8) :: VPLeaf !Vapor pressure in an average leaf (or canopy)
  real(r8) :: scal
  real(r8) :: XC,Stomata_Stress, dCanopyAvailWater
  real(r8) :: psi_H2O_flx
  real(r8) :: CanopyAvailWater      !available canopy water for canopy ET and hydraulic redistribution
  real(r8) :: CanopyAvailSnow       !available canopy snow for sublimation
  real(r8) :: CumHeatPlant2Soil,CumHeatSoil2Plant
  real(r8) :: CumWaterPlant2Soil,CumWaterSoil2Plant
  real(r8) :: dRadNetdT
  integer :: I, J

  real(r8), intent(in) :: LiveAreaFraction ! exposed live / total live+dead surface area [-]
  character(len=*), parameter :: subname='CanopyEnergyH2OIter_func'
  integer  :: IC
  logical :: LIterationExit
  real(r8) :: DPSI_old,DTKUpdateRate,VHCPYZ,DIFFW
  real(r8), parameter :: VHCPYM=0.838E-04_r8
  real(r8), parameter :: DTKUpdateRate1=0.5_r8
  real(r8), parameter :: DTKUpdateRate2=0.1_r8
  character(len=64) :: fmt
!     return variables
  integer :: NN
!     local variables
  integer :: N,L
!     begin_execution

  call PrintInfo('beg '//subname)
  associate(                                                               &
    AREA3                       => plt_site%AREA3                         ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    OrganOsmoPsi0pt_pft         => plt_ew%OrganOsmoPsi0pt_pft             ,& !input  :Organ osmotic potential when canopy water potential = 0 MPa, [MPa]
    CanopyBiomWater_pft         => plt_ew%CanopyBiomWater_pft             ,& !input  :canopy water content, [m3 d-2]
    CanopyNonstOMConc_pft       => plt_biom%CanopyNonstOMConc_pft         ,& !input  :canopy nonstructural C+N+P concentration, [g gC-1]
    FracPARads2Canopy_pft       => plt_rad%FracPARads2Canopy_pft          ,& !input  :fraction of incoming PAR absorbed by canopy, [-]
    LWRadGrnd_col               => plt_rad%LWRadGrnd_col                  ,& !input  :longwave radiation emitted by ground surface, [MJ m-2 h-1]
    LWRadSky_col                => plt_rad%LWRadSky_col                   ,& !input  :sky longwave radiation , [MJ d-2 h-1]
    MaxSoilLays4Root_pft        => plt_morph%MaxSoilLays4Root_pft         ,& !input  :maximum soil layer number for all root axes,[-]
    EMS_Modify_Scalar_col       => plt_ew%EMS_Modify_Scalar_col           ,& !input  :canopy longwave radiation emissivity scalar      
    Myco_pft                    => plt_morph%Myco_pft                     ,& !input  :mycorrhizal type (no or yes),[-]
    NU                          => plt_site%NU                            ,& !input  :current soil surface layer number, [-]
    PSICanopyOsmo_pft           => plt_ew%PSICanopyOsmo_pft               ,& !input  :canopy osmotic water potential, [Mpa]
    PSICanopyTurg_pft           => plt_ew%PSICanopyTurg_pft               ,& !input  :plant canopy turgor water potential, [MPa]
    RainIntcptByCanopy_pft      => plt_ew%RainIntcptByCanopy_pft          ,& !input  :water flux into canopy, [m3 d-2 h-1]
    SnowIntcptByCanopy_pft      => plt_ew%SnowIntcptByCanopy_pft          ,& !input  :snow flux onto canopy, [m3 d-2 h-1]
    RIB                         => plt_ew%RIB                             ,& !input  :Richardson number for calculating boundary layer resistance, [-]
    RadSWCanopyAbsorption_pft   => plt_rad%RadSWCanopyAbsorption_pft      ,& !input  :canopy absorbed shortwave radiation, [MJ d-2 h-1]
    RawIsoTCanopy2Atm_pft       => plt_ew%RawIsoTCanopy2Atm_pft           ,& !input  :canopy roughness height, [m]
    TKS_vr                      => plt_ew%TKS_vr                          ,& !input  :mean annual soil temperature, [K]
    TairK                       => plt_ew%TairK                           ,& !input  :air temperature, [K]
    VPA                         => plt_ew%VPA                             ,& !input  :vapor concentration, [m3 m-3]
    WatHeldOnCanopy_pft         => plt_ew%WatHeldOnCanopy_pft             ,& !input  :canopy surface water content, [m3 d-2]
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft                ,& !input  :threshold zero for plang growth calculation, [-]
    ZERO4LeafVar_pft            => plt_biom%ZERO4LeafVar_pft              ,& !input  :threshold zero for leaf calculation, [-]
    RPlantRootH2OUptk_pvr       => plt_ew%RPlantRootH2OUptk_pvr           ,& !inoput :root water uptake, [m3 d-2 h-1]
    PSICanopy_pft               => plt_ew%PSICanopy_pft                   ,& !inoput :canopy total water potential, [Mpa]
    TKCanopy_pft                => plt_ew%TKCanopy_pft                    ,& !inoput :canopy/vegetation temperature, [K]
    CanPStomaResistH2O_pft      => plt_photo%CanPStomaResistH2O_pft       ,& !output :canopy stomatal resistance, [h m-1]
    RawH2OCanopy2Atm_pft        => plt_photo%RawH2OCanopy2Atm_pft         ,& !output :canopy boundary layer resistance, [h m-1]
    RawHeatCanopy2Atm_pft       => plt_photo%RawHeatCanopy2Atm_pft        ,& !output :canopy sensible-heat boundary layer resistance, [h m-1]
    CanopyEvapTransLHeat_pft    => plt_ew%CanopyEvapTransLHeat_pft        ,& !output :canopy latent heat flux, [MJ d-2 h-1]
    LWRadCanopy_pft             => plt_rad%LWRadCanopy_pft                ,& !output :canopy longwave radiation, [MJ d-2 h-1]
    RadNet2Canopy_pft           => plt_rad%RadNet2Canopy_pft              ,& !output :canopy net radiation, [MJ d-2 h-1]
    TKC_pft                     => plt_ew%TKC_pft                         ,& !output :canopy temperature, [K]
    SnowOnCanopy_pft            => plt_ew%SnowOnCanopy_pft                ,& !input  :snow water equivalent held on canopy, [m3 d-2]        
    Transpiration_pft           => plt_ew%Transpiration_pft               ,& !output :canopy transpiration, [m2 d-2 h-1]
    VHeatCapCanopy_pft          => plt_ew%VHeatCapCanopy_pft              ,& !output :canopy heat capacity, [MJ d-2 K-1]
    VapXAir2CanopyLiq_pft       => plt_ew%VapXAir2CanopyLiq_pft           ,& !output :canopy evaporation, [m2 d-2 h-1]
    VapXAir2Canopy_pft          => plt_ew%VapXAir2Canopy_pft              ,& !output :canopy evaporation+sublimation, [m2 d-2 h-1]    
    SnoSub2AirCanopy_pft        => plt_ew%SnoSub2AirCanopy_pft            ,& !output :canopy snow sublimation,[m3 d-2 h-1]        
    DeltaTKC_pft                => plt_ew%DeltaTKC_pft                    ,& !output :change in canopy temperature, [K]
    QdewCanopy_pft              => plt_ew%QdewCanopy_pft                  ,& !output :dew fall on to canopy, [m3 H2O d-2 h-1]
    RootH2OUptkStress_pvr       => plt_ew%RootH2OUptkStress_pvr            & !output :root water uptake stress indicated by rate, [m3 d-2 h-1]
  )
  I=yearIJ%I;J=yearIJ%J
  !CCPOLT: total nonstructural canopy C,N,P concentration
  !FTHRM:coefficient for LW emitted by canopy
  !LWRad2Canopy:long-wave absorbed by canopy
  FTHRM        = EMMC*stefboltz_const*FracPARads2Canopy_pft(NZ)*AREA3(NU)*EMS_Modify_Scalar_col
  LWRad2Canopy = (LWRadSky_col+LWRadGrnd_col)*FracPARads2Canopy_pft(NZ)

  CCPOLT = CanopyNonstOMConc_pft(NZ)
  !      
  CanopyAvailWater    = CanopyBiomWater_pft(NZ)+WatHeldOnCanopy_pft(NZ)+RainIntcptByCanopy_pft(NZ)
  CanopyAvailSnow     = SnowOnCanopy_pft(NZ)+SnowIntcptByCanopy_pft(NZ)

  EffPFTGridArea4H2O  = FracGridPFTCovered*AREA3(NU)                    !area coverd by pft m2
  EffPFTGridArea4Heat = FracGridPFTCovered*AREA3(NU)*1.25E-03_r8        !1.25e-3 is heat capacity of air [MJ/(m3 K)], assuming air volume is not affected by plant mass
  RA1                 = RawIsoTCanopy2Atm_pft(NZ)                   !canopy isothermal boundary later resistance
  !
  IC                    = 0
  XC                    = 0.5_r8   !learning/updating rate for canopy temperature, not used now
  LIterationExit        = .false.
  PSICanPPre            = 0.0_r8
  PTransPre             = 0.0_r8
  cumPRootH2OUptakePre  = 0.0_r8
  SymplasmicWaterPrev   = 0.0_r8
  DPSI_old              = 0._r8
  cumPRootH2OUptake     = 0.0_r8
  CumPlantHeatLoss2Soil = 0._r8
  if (.not.HydraulicPathAvailable) then
    RPlantRootH2OUptk_pvr(:,:,NZ)=0._r8
    RootH2OUptkStress_pvr(:,:,NZ)=0._r8
  endif
  DIFF=1._r8
  VHCPYZ                = VHCPYM*AREA3(NU)
  D4000: DO NN=1,MaxIterNum
    !
    !     NET RADIATION FROM ABSORBED SW AND NET LW
    !
    !     LWRadCanopy_pft=LW emitted by canopy
    !     DTHS1=net LW absorbed by canopy
    !
    TKC1                  = TKCanopy_pft(NZ)
    LWRadCanopy_pft(NZ)   = FTHRM*TKC1**4              !long wave radiation
    DTHS1                 = LWRad2Canopy-LWRadCanopy_pft(NZ)*2.0_r8        !net long wave radiation to canopy
    RadNet2Canopy_pft(NZ) = RadSWCanopyAbsorption_pft(NZ)+DTHS1  !total radiation to canopy
    dRadNetdT             = -8._r8*FTHRM*TKC1**3
    !
    !     BOUNDARY LAYER RESISTANCE FROM RICHARDSON NUMBER
    !
    !     RI=Ricardson's number
    !     RA=canopy boundary layer resistance
    !     EvapConductCanopy,CdHeatCanopyAir=canopy latent,sensible heat conductance
    !
    ! Relax toward the same resistance used by the passive scalar solve.
    RawH2OCanopy2Atm_pft(NZ)  = MAX(RawMinH2OCanopy2Atm_pft,0.9_r8*RA1, &
      MIN(1.1_r8*RA1,CanopyVaporResistance(RawIsoTCanopy2Atm_pft(NZ),RIB,TairK,TKC1)))
    RA1                       = RawH2OCanopy2Atm_pft(NZ)
    RawHeatCanopy2Atm_pft(NZ) = RawH2OCanopy2Atm_pft(NZ)/BndlResistH2OToHeatRatio
    CdVapCanopyAir            = EffPFTGridArea4H2O/RawH2OCanopy2Atm_pft(NZ)
    CdHeatCanopyAir           = EffPFTGridArea4Heat/RawHeatCanopy2Atm_pft(NZ) !canopy-air heat conductance for sensible heat
    !
    !     CANOPY WATER AND OSMOTIC POTENTIALS
    !
    !     OrganOsmoPsi0pt_pft=osmotic potential at PSICanopy_pft=0 from PFT file
    !     PSICanopyOsmo_pft,PSICanopyTurg_pft=canopy osmotic,turgor water potential
    !
    call update_osmo_turg_pressure(PSICanopy_pft(NZ),CCPOLT,OrganOsmoPsi0pt_pft(NZ),TKC1 &
       ,PSICanopyOsmo_pft(NZ),PSICanopyTurg_pft(NZ),FDMP)
    !
    !     CANOPY STOMATAL RESISTANCE
    !
    CanPStomaResistH2O_pft(NZ) = CalcStomataResist4H2O(NZ)
    !
    !     CANOPY VAPOR PRESSURE AND EVAPORATION OF INTERCEPTED WATER
    !     OR TRANSPIRATION OF UPTAKEN WATER
    !
    !     VPLeaf,VPA=vapor pressure inside canopy, in atmosphere
    !     TKC1=canopy/vegetation temperature
    !     EX=canopy-atmosphere water flux
    !     RA,ResistLeafSurf=canopy boundary layer,surface resistance
    !     VapXAir2Canopy_pft=water flux to,from air to canopy surfaces, [ton H2O/(h*m2)]
    !     WatHeldOnCanopy_pft=water volume on canopy surfaces
    !     EP=water flux to,from inside canopy
    !     CanopyEvapTransLHeat_pft=canopy latent heat flux
    !     HeatEvapSens=convective heat flux from CanopyEvapTransLHeat_pft, <0, into air
    !     VAP=latent heat of evaporation
    !     EvapConductCanopy=aerodynamic conductance, [m/h]
    !
    VPLeaf = vapsat(TKC1)*EXP(18.0_r8*AMAX1(PSICanopy_pft(NZ),-5000._r8)/(RGASC*TKC1))  !ton H2O/m3
    !
    !area-scaled air to canopy water vap gradient, [m/h]*[ton/m3] = [ton H2O/(h*m2)] < 0, transpire into air
    EX                     = CdVapCanopyAir*(VPA-VPLeaf)
    if(LiveAreaFraction<1._r8) &
      EX=CdVapCanopyAir*(VPA-(LiveAreaFraction*VPLeaf+(1._r8-LiveAreaFraction)*vapsat(TKC1)))
    VapXAir2Canopy_pft(NZ) = EX*RawH2OCanopy2Atm_pft(NZ)/(RawH2OCanopy2Atm_pft(NZ)+ResistLeafSurf)
    !
    !Dew condensation to canopy >0 to canopy
    QdewCanopy_pft(NZ)     = 0._r8
    IF(EX.GT.0.0_r8)THEN
      QdewCanopy_pft(NZ)        = VapXAir2Canopy_pft(NZ)
      VapXAir2CanopyLiq_pft(NZ) = VapXAir2Canopy_pft(NZ)
      SnoSub2AirCanopy_pft(NZ)  = 0._r8
      EX                        = 0.0_r8
      HeatEvapSens              = VapXAir2Canopy_pft(NZ)*cpw*TairK               !enthalpy of condensed water add to canopy, MJ/(h)
      
    ELSEIF(EX.LE.0.0_r8)THEN 
      !canopy loses water, and there is canopy-held water  
      if(WatHeldOnCanopy_pft(NZ).GT.0.0_r8 .OR. SnowOnCanopy_pft(NZ).GT.0._r8)THEN
        !evaporation, and there is water stored in canopy
        !VapXAir2Canopy_pft <0._r8, off canopy as evaporation, cannot be more than WatHeldOnCanopy_pft(NZ)
        !
        !evaporation
        VapXAir2CanopyLiq_pft(NZ) = AMAX1(VapXAir2Canopy_pft(NZ),-WatHeldOnCanopy_pft(NZ))
        VapXAir2Canopy_pft(NZ)    = VapXAir2Canopy_pft(NZ)-VapXAir2CanopyLiq_pft(NZ)
        HeatEvapSens              = VapXAir2CanopyLiq_pft(NZ)*cpw*TKC1             !enthalpy of evaporated water leaving canopy

        !sublimation 
        if(VapXAir2Canopy_pft(NZ).LE.0._r8)then
          SnoSub2AirCanopy_pft(NZ) = AMAX1(VapXAir2Canopy_pft(NZ),-SnowOnCanopy_pft(NZ))
          HeatEvapSens = HeatEvapSens+SnoSub2AirCanopy_pft(NZ)*cps*TKC1
        endif
        !
        VapXAir2Canopy_pft(NZ)=VapXAir2CanopyLiq_pft(NZ)+SnoSub2AirCanopy_pft(NZ)
        EX                    = EX-VapXAir2Canopy_pft(NZ)                        !demand for transpiration
      ELSE
        VapXAir2Canopy_pft(NZ)    = 0._r8
        VapXAir2CanopyLiq_pft(NZ) = 0.0_r8
        SnoSub2AirCanopy_pft(NZ)  = 0._r8
        HeatEvapSens              = 0._r8
      ENDIF
    ENDIF
    !
    !Transpiration_pft<0 means canopy lose water through transpiration
    !CanopyEvapTransLHeat_pft:latent heat flux, negative means into atmosphere
    !
    if(LiveAreaFraction<1._r8) &
      EX=MIN(0._r8,LiveAreaFraction*(CdVapCanopyAir*(VPA-VPLeaf)-VapXAir2Canopy_pft(NZ)))
    Transpiration_pft(NZ)        = EX*RawH2OCanopy2Atm_pft(NZ)/(RawH2OCanopy2Atm_pft(NZ)+CanPStomaResistH2O_pft(NZ))
    if (.not.HydraulicPathAvailable) Transpiration_pft(NZ)=0._r8
    CanopyEvapTransLHeat_pft(NZ) = (Transpiration_pft(NZ)+VapXAir2CanopyLiq_pft(NZ))*EvapLHTC+SnoSub2AirCanopy_pft(NZ)*SublmHTC
    HeatEvapSens                 = HeatEvapSens + Transpiration_pft(NZ)*cpw*TKC1
    !
    !     SENSIBLE + STORAGE HEAT FROM RN, LE AND CONVECTIVE HEAT FLUXES
    !
    !     HeatSenAddStore=initial estimate of sensible+storage heat flux
    !     PrecpHeatbyCanopy=convective heat flux from precip to canopy
    !
    HeatAdd2Can=RadNet2Canopy_pft(NZ)+CanopyEvapTransLHeat_pft(NZ)+HeatEvapSens+PrecpHeatbyCanopy-CumPlantHeatLoss2Soil
    !
    !     SOLVE FOR CANOPY TEMPERATURE CAUSED BY SENSIBLE + STORAGE HEAT
    !
    !     VHeatCapCanopy_pft=canopy heat capacity
    !     TKCY=equilibrium canopy temperature for HeatSenAddStore
    !
    !   VHeatCapCanopyPrev_pft= canopy heat capacity, including biomass, held and bounded water
    !   assuming transpiration does not change water binded by the canopy, i.e. CanopyBiomWater_pft(NZ)
    !   canopy bound water is changing
    !
    VHeatCapCanopy_pft(NZ)=VHeatCapCanopyPrev_pft+cpw*(RainIntcptByCanopy_pft(NZ)+VapXAir2CanopyLiq_pft(NZ) &
       +Transpiration_pft(NZ)-cumPRootH2OUptake)+cps*(SnowIntcptByCanopy_pft(NZ)+SnoSub2AirCanopy_pft(NZ))
    !   
    ! new canopy temperature 
    ! Linearize radiation about the current trial temperature TKC1; TKCX
    ! belongs only to the retained beginning inventory.
    TKCY=(TKCX*VHeatCapCanopyPrev_pft-dRadNetdT*TKC1+HeatAdd2Can+TairK*CdHeatCanopyAir) &
      /(VHeatCapCanopy_pft(NZ)+CdHeatCanopyAir-dRadNetdT)

    ! Limit each nonlinear step, not the final departure from air. A fixed
    ! air-temperature bound can exclude the heat-balanced solution.
    TKCY=AMIN1(TKC1+5.0_r8,AMAX1(TKC1-5.0_r8,TKCY))
    !
    !     RESET CANOPY TEMPERATURE FOR NEXT ITERATION
    !
    !     XC,IC=magnitude,direction of change in canopy temp for next cycle
    !
    IF((IC.EQ.0 .AND. TKCY.GT.TKC1)        & !canopy temperature shifts from cooling or no change to warming
      .OR. (IC.EQ.1 .AND. TKCY.LT.TKC1))THEN !canopy temperature shifts from warming to cooling
      !reducing update rate when temperature update change direction
      XC=0.5_r8*XC
    ENDIF
    !
    !obtain the updating rate
    if(VHeatCapCanopy_pft(NZ).GT.VHCPYZ)then
      !fast update if heat capacity is beyond threshold
      DTKUpdateRate=XC*DTKUpdateRate1
    else
      !slow update when canopy heat capacity is small.
      DTKUpdateRate=XC*DTKUpdateRate2    
    endif
    !
    TKCanopy_pft(NZ)=TKC1+DTKUpdateRate*(TKCY-TKC1) !maxium change 5 K
    !
    IF(TKCY.GT.TKC1)THEN
      !warming
      IC=1
    ELSE
      !cooling or no change
      IC=0
    ENDIF
    ! Match the capacity and accepted fluxes committed by UpdatePlantWaterVars.
    VHeatCapCanopy_pft(NZ)=VHeatCapCanopyPrev_pft+cpw*(RainIntcptByCanopy_pft(NZ)+VapXAir2CanopyLiq_pft(NZ) &
      +Transpiration_pft(NZ)-cumPRootH2OUptake)+cps*(SnowIntcptByCanopy_pft(NZ)+SnoSub2AirCanopy_pft(NZ))
    HeatResidual=TKCX*VHeatCapCanopyPrev_pft-TKCanopy_pft(NZ)*VHeatCapCanopy_pft(NZ) &
      +RadNet2Canopy_pft(NZ)+CanopyEvapTransLHeat_pft(NZ)+HeatEvapSens+PrecpHeatbyCanopy &
      -CumPlantHeatLoss2Soil+CdHeatCanopyAir*(TairK-TKCanopy_pft(NZ))
    ! A disconnected canopy still exchanges heat and intercepted water with
    ! the atmosphere. Its biomass reservoir and root fluxes remain unchanged.
    if (.not.HydraulicPathAvailable) then
      if (ABS(HeatResidual)<CanopyHeatTolerance) then
        DIFF=0._r8
        exit D4000
      endif
      cycle D4000
    endif
    !
    !     IF CONVERGENCE CRITERION IS MET OR ON EVERY TENTH ITERATION,
    !     PROCEED TO WATER BALANCE
    !
    !
    IF(ABS(TKCY-TKC1).LT.0.05_r8 .OR. (NN/10)*10.EQ.NN)THEN
      PSILC   = PSICanopy_pft(NZ)+PSIGravCanopyHeight
      !
      !     ROOT WATER UPTAKE FROM SOIL-CANOPY WATER POTENTIALS,
      !     SOIL + ROOT HYDRAULIC RESISTANCES
      !
      !     SoiLayerHasRoot_rvr=rooted layer flag
      !     RPlantRootH2OUptk_pvr=root water uptake from soil layer > 0, [m3 H2O/h] = [MPa] /[], => []=[m3 H2O/(MPa h)]=[m3 H2O/(kg m-1 h-2 h)]
      !     =[kg H2O/(kg m-1 h-1)]=[h m] over grid, =[h/m] specific area
      !     WatAvail4Uptake_vr,AirMicPore4Fill_vr=water volume available for uptake,air volume
      !     FracPRoot4Uptake_pvr=PFT fraction of biome root mass
      !     PSILC=height corrected canopy water potential 
      !     SoilRootResist4H2O_pvr=total soil+root resistance
      !     cumPRootH2OUptake=total root water uptake from soil  (<0 addition to canopy)
      !     CumPlantHeatLoss2Soil (<0.), add heat to canopy 
      !
      CumHeatPlant2Soil  = 0._r8;CumHeatSoil2Plant  = 0._r8
      CumWaterPlant2Soil = 0._r8;CumWaterSoil2Plant = 0._r8
      D4200: DO N=1,Myco_pft(NZ)
        D4201: DO L=NU,MaxSoilLays4Root_pft(NZ)
          !only update layers with roots
          IF(SoiLayerHasRoot_rvr(N,L))THEN
            !psi_H2O_flx<0 plant doing active H2O uptake, [mH2O d-2]            
            psi_H2O_flx=(PSILC-TotalSoilPSIMPa_vr(L))/SoilRootResist4H2O_pvr(N,L)

            RootH2OUptkStress_pvr(N,L,NZ) = AZMAX1(psi_H2O_flx)
            
            !if psi_H2O_flx >0, AMIN1 ensures it cannot fill more than available air-space
            !if psi_H2O_flx <0, AMAX1 ensures it cannot take more than available water

            RPlantRootH2OUptk_pvr(N,L,NZ) = AMAX1(AZMIN1(-WatAvail4Uptake_vr(L)*FracPRoot4Uptake_pvr(N,L,NZ)), &
              AMIN1(psi_H2O_flx, AirMicPore4Fill_vr(L)*FracPRoot4Uptake_pvr(N,L,NZ)))
            
            !plant roots/myco lose water to soil > 0, release water slowly back to soil
            IF(RPlantRootH2OUptk_pvr(N,L,NZ).GT.0.0_r8)THEN              

              !plant moves heat from canopy to soil
              CumHeatPlant2Soil  = CumHeatPlant2Soil+RPlantRootH2OUptk_pvr(N,L,NZ)*cpw*TKC1
              CumWaterPlant2Soil = CumWaterPlant2Soil+RPlantRootH2OUptk_pvr(N,L,NZ)
            ELSE
              !plant/myco gains water from soil < 0
              CumHeatSoil2Plant  = CumHeatSoil2Plant+RPlantRootH2OUptk_pvr(N,L,NZ)*cpw*TKS_vr(L)
              CumWaterSoil2Plant = CumWaterSoil2Plant+RPlantRootH2OUptk_pvr(N,L,NZ)
            ENDIF
          ELSE
            RootH2OUptkStress_pvr(N,L,NZ)=0._r8
            RPlantRootH2OUptk_pvr(N,L,NZ)=0.0_r8
          ENDIF
        enddo D4201
      ENDDO D4200

      !check for water availability
      dCanopyAvailWater=CanopyAvailWater-CumWaterSoil2Plant+Transpiration_pft(NZ)+VapXAir2CanopyLiq_pft(NZ) 
      dCanopyAvailWater=dCanopyAvailWater+CanopyAvailSnow+SnoSub2AirCanopy_pft(NZ) 
      !
      !plant does not have enough water to support release to soil
      if(dCanopyAvailWater.GE.CumWaterPlant2Soil)then
        !sufficient water to release into soil
        cumPRootH2OUptake     = CumWaterPlant2Soil+CumWaterSoil2Plant
        CumPlantHeatLoss2Soil = CumHeatPlant2Soil+CumHeatSoil2Plant
      elseif(CumWaterPlant2Soil.GT.1.e-12_r8)then
        !insufficient water, but still release to soil, downscale water release to soil
        scal=AZMAX1(dCanopyAvailWater)/CumWaterPlant2Soil*0.9999_r8
        D4202: DO N=1,Myco_pft(NZ)
          D4203: DO L=NU,MaxSoilLays4Root_pft(NZ)        
            if(RPlantRootH2OUptk_pvr(N,L,NZ)>0._r8)then
              RPlantRootH2OUptk_pvr(N,L,NZ)=RPlantRootH2OUptk_pvr(N,L,NZ)*scal
            endif
          ENDDO D4203
        ENDDO D4202
        cumPRootH2OUptake     = scal*CumWaterPlant2Soil+CumWaterSoil2Plant
        CumPlantHeatLoss2Soil = scal*CumHeatPlant2Soil+CumHeatSoil2Plant        
      else
        if(ldo_sp_mode)then
          !prescribed phenology: no release to the soil, so the plant receives
          !what it takes from the soil and the per-layer uptake
          !(RPlantRootH2OUptk_pvr, the soil sink) and the plant agree
          cumPRootH2OUptake     = CumWaterSoil2Plant
          CumPlantHeatLoss2Soil = CumHeatSoil2Plant
        else
          !full plant model unchanged (same inconsistency, not yet reviewed)
          cumPRootH2OUptake     = 0._r8
          CumPlantHeatLoss2Soil = 0._r8
        endif
      endif
      !
      !     TEST TRANSPIRATION - ROOT WATER UPTAKE VS. CHANGE IN CANOPY
      !     WATER STORAGE
      !
      !     SymplasmicWater  = total water can be held in biomass
      !     CanopyBiomWater_pft= canopy water content held in biomass
      !     DIFFZ= extra canopy water binding capacity
      !     DIFFU= change to canopy binded water
      !     DIFFU-DIFFZ=residual of canopy binded water,
      !     >0 gain binded water (increase leaf P, more positive), <0 lose binded water (more uptake, decrease leaf water P)
      !     DIFF=normalized difference between DIFFZ and DIFFU
      !     5.0E-03=acceptance criterion for DIFF
      !     RSSZ=change in canopy water potl vs change in canopy water cnt
      !     RSSU=change in canopy water potl vs change in transpiration
      !
      SymplasmicWater = 1.e-6_r8*CanopyMassC4H2OStorage/FDMP
      !
      DIFFZ = SymplasmicWater-CanopyBiomWater_pft(NZ)    !biomass water deficit, due to change in canopy water potential and C biomass
      DIFFU = Transpiration_pft(NZ)-cumPRootH2OUptake    !root water uptake excess, [m H2O h-1 d-2]
      DIFFW = DIFFU-DIFFZ
      !
      !ideally, the difference between DIFFZ and DIFFU should be as small as possible '
      !
      IF(.not.isclose(cumPRootH2OUptake,0.0_r8))THEN
        DIFF=ABS((DIFFW)/cumPRootH2OUptake)
      ELSEIF(SymplasmicWater.GT.0._r8)THEN
        DIFF=ABS(DIFFW)/SymplasmicWater
      ELSE
        DIFF=0._r8  
      ENDIF

      ! Water convergence also requires closure of the accepted heat budget.
    ! Match the capacity and accepted fluxes committed by UpdatePlantWaterVars.
    VHeatCapCanopy_pft(NZ)=VHeatCapCanopyPrev_pft+cpw*(RainIntcptByCanopy_pft(NZ)+VapXAir2CanopyLiq_pft(NZ) &
      +Transpiration_pft(NZ)-cumPRootH2OUptake)+cps*(SnowIntcptByCanopy_pft(NZ)+SnoSub2AirCanopy_pft(NZ))
    HeatResidual=TKCX*VHeatCapCanopyPrev_pft-TKCanopy_pft(NZ)*VHeatCapCanopy_pft(NZ) &
      +RadNet2Canopy_pft(NZ)+CanopyEvapTransLHeat_pft(NZ)+HeatEvapSens+PrecpHeatbyCanopy &
      -CumPlantHeatLoss2Soil+CdHeatCanopyAir*(TairK-TKCanopy_pft(NZ))
      IF(DIFF.LT.DIFFY)THEN
        IF(LIterationExit.and.ABS(HeatResidual)<CanopyHeatTolerance)EXIT
        if(.not.LIterationExit)then
          LIterationExit=.true.
          CALL StomatalDynamics(yearIJ,NZ,RawH2OCanopy2Atm_pft(NZ))
        endif
        CYCLE
      ENDIF
      !
      !resistance based on biome-associated water
      IF(ABS(SymplasmicWater-SymplasmicWaterPrev).GT.ZERO4Groth_pft(NZ))THEN
        RSSZ=ABS((PSICanopy_pft(NZ)-PSICanPPre)/(SymplasmicWater-SymplasmicWaterPrev))
      ELSEIF(CdH2ORootxSoil.GT.ZERO4Groth_pft(NZ))THEN
        RSSZ=1.0_r8/CdH2ORootxSoil   !resistance
      ELSE
        RSSZ=ZERO4LeafVar_pft(NZ)
      ENDIF

      !resistance based on transpiration flux
      IF(ABS(Transpiration_pft(NZ)-PTransPre).GT.ZERO4Groth_pft(NZ))THEN
        RSSUX=ABS((PSICanopy_pft(NZ)-PSICanPPre)/(Transpiration_pft(NZ)-PTransPre))
        IF(CdH2ORootxSoil.GT.ZERO4Groth_pft(NZ))THEN
          RSSU=AMIN1(1.0_r8/CdH2ORootxSoil,RSSUX)  !resistance for uptake
        ELSE
          RSSU=RSSUX
        ENDIF
        !
        !resistance based on root water uptake
      ELSEIF(ABS(cumPRootH2OUptake-cumPRootH2OUptakePre).GT.ZERO4Groth_pft(NZ))THEN
        RSSUX=ABS((PSICanopy_pft(NZ)-PSICanPPre)/(cumPRootH2OUptake-cumPRootH2OUptakePre))
        IF(CdH2ORootxSoil.GT.ZERO4Groth_pft(NZ))THEN
          RSSU=AMIN1(1.0_r8/CdH2ORootxSoil,RSSUX)  !resistance
        ELSE
          RSSU=RSSUX
        ENDIF
      ELSEIF(CdH2ORootxSoil.GT.ZERO4Groth_pft(NZ))THEN
        RSSU=1.0_r8/CdH2ORootxSoil
      ELSE
        RSSU=ZERO4LeafVar_pft(NZ)
      ENDIF
      !
      !     CHANGE IN CANOPY WATER POTENTIAL REQUIRED TO BRING AGREEMENT
      !     BETWEEN TRANSPIRATION - ROOT WATER UPTAKE AND CHANGE IN CANOPY
      !     WATER STORAGE
      !
      !     DPSI=change in PSICanopy_pft for next convergence cycle, MPa = [Cd] * m H2O d-2, [Cd]= MPa d2/m
      !     1.0E-03=acceptance criterion for DPSI
      !
      DPSI=MAX(-ABS(PSICanopy_pft(NZ)), &
        MIN(MIN(RSSZ,RSSU)*DIFFW,ABS(PSICanopy_pft(NZ))))
      !
      !     IF CONVERGENCE CRITERION IS MET THEN FINISH,
      !     OTHERWISE START NEXT ITERATION WITH CANOPY WATER POTENTIAL
      !     TRANSPIRATION, UPTAKE AND WATER CONTENT FROM CURRENT ITERATION
      !
      IF((NN.GE.30 .AND. ABS(DPSI).LT.1.0E-03_r8 .and. DIFF.LT.DIFFY) .OR. NN.GE.MaxIterNum)then
        IF(LIterationExit)EXIT
        LIterationExit=.true.
        CALL StomatalDynamics(yearIJ,NZ,RawH2OCanopy2Atm_pft(NZ))
      ELSE
        !prepare for next iteration
        PSICanPPre           = PSICanopy_pft(NZ)
        PTransPre            = Transpiration_pft(NZ)
        cumPRootH2OUptakePre = cumPRootH2OUptake
        SymplasmicWaterPrev    = SymplasmicWater

        PSICanopy_pft(NZ) = AMIN1(PSICanopy_pft(NZ)*0.5_r8,PSICanopy_pft(NZ)+0.5_r8*DPSI)
        DPSI_old          = DPSI
        XC                = 0.50_r8!
      ENDIF
    ENDIF
  ENDDO D4000

  call update_osmo_turg_pressure(PSICanopy_pft(NZ),CCPOLT,OrganOsmoPsi0pt_pft(NZ),TKC1 &
      ,PSICanopyOsmo_pft(NZ),PSICanopyTurg_pft(NZ),FDMP)  

  !
  !     FINAL CANOPY TEMPERATURE, DIFFERENCE WITH AIR TEMPERATURE
  !
  !     TKC=final estimate of canopy temperature TKCanopy_pft
  !     TairK=current air temperature
  !     DeltaTKC_pft=TKC-TairK for next hour
  !
  TKC_pft(NZ)      = TKCanopy_pft(NZ)
  DeltaTKC_pft(NZ) = TKC_pft(NZ)-TairK
  
  end associate
  call PrintInfo('end '//subname)
  end function CanopyEnergyH2OIter_func

!----------------------------------------------------------------------------------------------------
  subroutine CalcPlantHydroResistance(yearIJ,NZ,RadialMeanLen_rvr,FineRootRadius_rvr,RootEffLen4Absorption_pvr,&
    SoilResist4H2O_rvr,SoilRootResist4H2O_pvr,CdH2ORootxSoil,PSIGravCanopyHeight,SoiLayerHasRoot_rvr)

  implicit none
  type(yearIJ_type),intent(in) :: yearIJ
  integer, intent(in)   :: NZ
  real(r8), intent(in)  :: RadialMeanLen_rvr(pltpar%jroots,JZ1)   !Mean radial path length for water and nutrient uptake [m]
  real(r8), intent(in)  :: FineRootRadius_rvr(pltpar%jroots,JZ1)
  real(r8), intent(in)  :: RootEffLen4Absorption_pvr(pltpar%jroots,JZ1)
  real(r8), intent(out) :: SoilResist4H2O_rvr(pltpar%jroots,JZ1)         !soil resistance for water uptake
  real(r8), intent(out) :: SoilRootResist4H2O_pvr(pltpar%jroots,JZ1)   !added soil and root resistance, [MPa h d2 m-3]
  real(r8), intent(out) :: CdH2ORootxSoil                              !total root-soil conductance for plant H2O uptake
  real(r8), intent(out) :: PSIGravCanopyHeight    !gravimetric water potential at CanopyHeight4WatUptake_pft, [MPa]
  logical , intent(out) :: SoiLayerHasRoot_rvr(pltpar%jroots,JZ1)

  character(len=*), parameter :: subname='CalcPlantHydroResistance in '//trim(mod_filename)
  real(r8) :: RootRadialResist_rvr(pltpar%jroots,JZ1)      !radial root resistance for water uptake, [MPa h m-3]    
  real(r8) :: Root2ndAxialResist_rvr(pltpar%jroots,JZ1)    !secondary root axial resistance for water uptake,[MPa h m-3]      
  real(r8) :: Root1stAxialResist_rvr(pltpar%jroots,JZ1)    !primary root axial resistance for water uptake,[MPa h m-3]      
  real(r8) :: RootMediumAxialResist_rvr(pltpar%jroots,JZ1) !medium size root axial resistance for water uptake, [MPa h m-3]
  real(r8) :: FRADW 
  real(r8) :: FRAD1   !total transport vessels in primary roots
  real(r8) :: FRAD2
  real(r8) :: FRADM   !medium roots xylem/tracheid vessels
  real(r8) :: DTransptTube, AreaTranspt
  real(r8) :: RSSL,Root2ndSurfArea,condC,condM
  real(r8) :: StalkAxialResist     ![kg H2O m-5 h-1] =[m-2 h-1 H2O]
  integer :: N, L
  associate(                                                                  &
    DLYR3                       => plt_site%DLYR3                            ,& !input  :vertical thickness of soil layer, [m]  
    CanopyHeightLive_pft        => plt_morph%CanopyHeightLive_pft            ,& !input  :canopy height, [m]
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft           ,& !input  :plant growth type (vascular, non-vascular),[-]    
    CumSoilThickMidL_vr         => plt_site%CumSoilThickMidL_vr              ,& !input  :depth to middle of soil layer from surface of grid cell, [m]
    HYCDMicP4RootUptake_vr      => plt_soilchem%HYCDMicP4RootUptake_vr       ,& !input  :soil micropore hydraulic conductivity for root water uptake, [m MPa-1 h-1]
    MaxSoilLays4Root_pft        => plt_morph%MaxSoilLays4Root_pft            ,& !input  :maximum soil layer number for all root axes,[-]
    Root1stMaxRadius1_pft       => plt_morph%Root1stMaxRadius1_pft           ,& !input  :root diameter primary axes, [m]            
    Myco_pft                    => plt_morph%Myco_pft                        ,& !input  :mycorrhizal type (no or yes),[-]
    NU                          => plt_site%NU                               ,& !input  :current soil surface layer number, [-]
    PSICanopy_pft               => plt_ew%PSICanopy_pft                      ,& !input  :canopy total water potential, [Mpa]
    RootAge_rpvr                => plt_morph%RootAge_rpvr                    ,& !inoput :root age,[h]
    RootSingleVesselRstaxial_pft=> plt_morph%RootSingleVesselRstaxial_pft    ,& !input  :axial resistance for a single 1-m water transport vessel, [MPa h m-4]
    RootMatureAge_pft           => plt_morph%RootMatureAge_pft               ,& !input : Root maturation age, [h]
    PlantPopuLive_pft           => plt_site%PlantPopuLive_pft                ,& !input  :plant population, [d-2]
    Root1stRadius_pvr           => plt_morph%Root1stRadius_pvr               ,& !input  :root layer diameter primary axes, [m]
    Root1stXNumL_pvr            => plt_morph%Root1stXNumL_pvr                ,& !input  :root layer number primary axes, [d-2]
    Root2ndMaxRadius_pft        => plt_morph%Root2ndMaxRadius_pft            ,& !input  :maximum radius of secondary roots, [m]
    Root2ndEffLen4uptk_rpvr     => plt_morph%Root2ndEffLen4uptk_rpvr         ,& !input  :Layer effective root length four resource uptake, [m]
    Root2ndRadius_rpvr          => plt_morph%Root2ndRadius_rpvr              ,& !input  :root layer diameter secondary axes, [m]
    Root2ndXNumL_rpvr           => plt_morph%Root2ndXNumL_rpvr               ,& !input  :root layer number axes, [d-2]
    StalkAxialResist_pft        => plt_morph%StalkAxialResist_pft            ,& !input  :Stalk axial resistivity, [MPa h m-4]
    Root2ndAxialResist_pft      => plt_morph%Root2ndAxialResist_pft          ,& !input  :Root axial resistivity, [MPa h m-4]
    RootRadialResist_pft        => plt_morph%RootRadialResist_pft            ,& !input  :root radial resistivity, [MPa h m-1]    
    RootSingleVesselArea_pft    => plt_morph%RootSingleVesselArea_pft        ,& !input  :Single root vessel cross section area, [m2]
    RootLenDensPerPlant_pvr     => plt_morph%RootLenDensPerPlant_pvr         ,& !input  :root layer length density, [m m-3]
    RootAbsorbLenPerPlant_pvr   => plt_morph%RootAbsorbLenPerPlant_pvr       ,& !input  :total absorptive root length per plant in layer, [m p-1]
    THETW_vr                    => plt_soilchem%THETW_vr                     ,& !input  :volumetric water content, [m3 m-3]
    VLMicP_vr                   => plt_soilchem%VLMicP_vr                    ,& !input  :total volume in micropores, [m3 d-2]
    VLSoilPoreMicP_vr           => plt_soilchem%VLSoilPoreMicP_vr            ,& !input  :volume of soil layer, [m3 d-2]
    CRootLumenArea_pvr          => plt_morph%CRootLumenArea_pvr              ,& !input  :coarse roots lumen area, [m2]       
    MRootLumenArea_pvr          => plt_morph%MRootLumenArea_pvr              ,& !input  :medium roots lumen area, [m2]
    VLWatMicPM_vr               => plt_site%VLWatMicPM_vr                    ,& !input  :soil micropore water content, [m3 d-2]
    ZERO                        => plt_site%ZERO                             ,& !input  :threshold zero for numerical stability, [-]
    RootFineFrac2Med_rpvr       => plt_morph%RootFineFrac2Med_rpvr           ,& !fine-axis-count-weighted medium-root routing by category, [-]
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft                   ,& !input  :threshold zero for plang growth calculation, [-]
    ZEROS2                      => plt_site%ZEROS2                           ,& !input  :threshold zero for numerical stability,[-]
    CdH2ORootxSoil_pft          => plt_ew%CdH2ORootxSoil_pft                 ,& !output :total root and soil conductance for plant root water uptake, [mH2O h-1 d-2 MPa-1]
    RootMediumLength_pvr        => plt_morph%RootMediumLength_pvr            ,& !input  :root layer mean length for medium size axes, [m d-2]              
    RootResist4H2O_pvr          => plt_ew%RootResist4H2O_pvr                 ,& !output :total root (axial+radial) resistance for water uptake,[MPa h m-3]
    RootRadialKond2H2O_pvr      => plt_ew%RootRadialKond2H2O_pvr             ,& !output :radial root conductance for water uptake, [m5 H2O h-1 MPa-1]
    RootAxialKond2H2O_pvr       => plt_ew%RootAxialKond2H2O_pvr              ,& !output :axial root conductance for water uptake, [m3 H2O h-1 MPa-1]
    CanopyHeight4WatUptake_pft  => plt_morph%CanopyHeight4WatUptake_pft       & !output :canopy height, [m]
  )

  !     GRAVIMETRIC WATER POTENTIAL FROM CANOPY HEIGHT
  !
  !     CanopyHeight4WatUptake_pft=canopy height for water uptake
  !     FRADW=conducting elements of stalk relative to those of primary root, using the water poential to 
  !      mimic the stalk area expansion.
  !     PSICanopy_pft=canopy total water potential
  !     EMODW=wood modulus of elasticity (MPa)
  call PrintInfo('beg '//subname)

  CdH2ORootxSoil_pft(NZ)         = 0.0_r8
  CanopyHeight4WatUptake_pft(NZ) = 0.80_r8*CanopyHeightLive_pft(NZ)
  PSIGravCanopyHeight            = mGravAccelerat*CanopyHeight4WatUptake_pft(NZ)
  !
  !Hagen-Poiseuille law of stalk xylem/tracheid vessels
  FRADW  = 1.0E+04_r8*(AMAX1(0.5_r8,1.0_r8+PSICanopy_pft(NZ)/EMODW))**4
  !
  !     SOIL AND ROOT HYDRAULIC RESISTANCES TO ROOT WATER UPTAKE
  !
  !      VLSoilPoreMicP_vr,VLWatMicPM,THETW=soil,water volume,content
  !     RootLenDensPerPlant_pvr,RootAbsorbLenPerPlant_pvr=root length density,root length per plant
  !     HYCDMicP4RootUptake_vr=soil hydraulic conductivity for root uptake
  !     SoiLayerHasRoot_rvr:1=rooted,0=not rooted
  !     N:1=root,2=mycorrhizae
  !
  D3880: DO N=1,Myco_pft(NZ)
    DO  L=NU,MaxSoilLays4Root_pft(NZ)

      SoiLayerHasRoot_rvr(N,L)=VLSoilPoreMicP_vr(L).GT.ZEROS2         &
        .AND. VLWatMicPM_vr(NPH,L).GT.ZEROS2                          &
        .AND. RootAbsorbLenPerPlant_pvr(N,L,NZ).GT.ZERO               &
        .AND. HYCDMicP4RootUptake_vr(L).GT.ZERO                       &
        .AND. Root1stXNumL_pvr(L,NZ).GT.ZERO4Groth_pft(NZ)            &
        .AND. Root2ndXNumL_rpvr(N,L,NZ).GT.ZERO4Groth_pft(NZ)         &
        .AND. THETW_vr(L).GT.ZERO

      if(SoiLayerHasRoot_rvr(N,L))THEN

        !
        !     SOIL HYDRAULIC RESISTANCE FROM RADIAL UPTAKE GEOMETRY
        !     AND SOIL HYDRAULIC CONDUCTIVITY
        !
        !     SoiH2OResist=soil hydraulic resistance
        !     PP=plant population
        !     RadialMeanLen_rvr=mean radial path length of water and nutrient uptake
        !     FineRootRadius_rvr,RootEffLen4Absorption_pvr=root radius,surface/radius area
        
        RSSL                    = LOG(RadialMeanLen_rvr(N,L)/FineRootRadius_rvr(N,L))/(RootEffLen4Absorption_pvr(N,L)*PlantPopuLive_pft(NZ))
        SoilResist4H2O_rvr(N,L) = RSSL/HYCDMicP4RootUptake_vr(L)

        !
        !     RADIAL ROOT RESISTANCE FROM ROOT AREA AND RADIAL RESISTIVITY
        !     ENTERED IN 'READQ'
        !
        !     RootRadialResist_rvr=radial resistance
        !     RootRadialResist_pft=radial resistivity from PFT file
        !     VLMicP,VLWatMicPM=soil micropore,water volume
        !eq. (31) in Grant, 1998
        Root2ndSurfArea           = TwoPiCON*Root2ndRadius_rpvr(N,L,NZ)*RootAbsorbLenPerPlant_pvr(N,L,NZ)*PlantPopuLive_pft(NZ)
        RootRadialResist_rvr(N,L) = RootRadialResist_pft(N,NZ)*VLMicP_vr(L)/(Root2ndSurfArea*VLWatMicPM_vr(NPH,L))     ![MPa h m-1]/[m2]=[MPa h m-3]
        !
        !     ROOT AXIAL RESISTANCE FROM RADII AND LENGTHS OF PRIMARY AND
        !     SECONDARY ROOTS AND FROM AXIAL RESISTIVITY ENTERED IN 'READQ'
        !
        !     FRAD1,FRAD2=primary,secondary root radius relative to maximum, poisuelle flow conductance scalar
        !     secondary radius from PFT file Root2ndMaxRadius_pft at which Root2ndAxialResist_pft is defined
        !     Root1stAxialResist_rvr,Root2ndAxialResist_rvr=axial resistance of primary,secondary roots,[MPa h d2 m-3] =[MPa h m-1]=[]*m
        !     Root2ndEffLen4uptk_rpvr=Layer effective root length four resource uptake, [m]
        ! apply the Poiseuille relationship (Aguirrezabal et al., 1993, Grant, 1998)

        FRAD1 = AZMAX1(CRootLumenArea_pvr(L,NZ)/RootSingleVesselArea_pft(NZ),1._r8)
        FRADM = AZMAX1(MRootLumenArea_pvr(L,NZ)/RootSingleVesselArea_pft(NZ),1._R8)
        if(N.eq.ipltroot)then
          FRAD2                       = (Root2ndRadius_rpvr(N,L,NZ)/Root2ndMaxRadius_pft(N,NZ))**2
        else
          FRAD2                       = (Root2ndRadius_rpvr(N,L,NZ)/Root2ndMaxRadius_pft(N,NZ))**4
        endif
        StalkAxialResist            = AZMAX1(StalkAxialResist_pft(NZ)*CanopyHeight4WatUptake_pft(NZ)/(FRADW*Root1stXNumL_pvr(L,NZ)))
        Root1stAxialResist_rvr(N,L) = StalkAxialResist+RootSingleVesselRstaxial_pft(NZ)*CumSoilThickMidL_vr(L)/FRAD1

        if(FRADM.GT.0._r8)then
          RootMediumAxialResist_rvr(N,L)=0.5_r8*RootSingleVesselRstaxial_pft(NZ)*RootMediumLength_pvr(L,NZ)/FRADM
        ELSE
          RootMediumAxialResist_rvr(N,L)=0._r8
        endif
        Root2ndAxialResist_rvr(N,L) = Root2ndAxialResist_pft(N,NZ)*Root2ndEffLen4uptk_rpvr(N,L,NZ)/(FRAD2*Root2ndXNumL_rpvr(N,L,NZ))


        !
        !     TOTAL ROOT RESISTANCE = SOIL + RADIAL + AXIAL
        !
        !     RootResist=root radial+axial resistance
        !     SoilRootResist4H2O_pvr=total soil+root resistance, [MPa h d2 m-3]
        !     CdH2ORootxSoil=total soil+root conductance for all layers
        ! assuming all roots work in parallel
        RootRadialKond2H2O_pvr(N,L,NZ) = 1._r8/RootRadialResist_rvr(N,L)
        if(RootFineFrac2Med_rpvr(N,L,NZ).GT.0._R8)THEN
          condM=1._r8/(Root2ndAxialResist_rvr(N,L)/RootFineFrac2Med_rpvr(N,L,NZ)+RootMediumAxialResist_rvr(N,L))
          condC=(1._r8-RootFineFrac2Med_rpvr(N,L,NZ))/Root2ndAxialResist_rvr(N,L)
          RootAxialKond2H2O_pvr(N,L,NZ)  = DLYR3(L)**2/(Root1stAxialResist_rvr(N,L)+1._r8/(condC+condM))     !plant size-scaled axial root conductance to H2O        
          RootResist4H2O_pvr(N,L,NZ)     = RootRadialResist_rvr(N,L)+Root1stAxialResist_rvr(N,L)+1._r8/(condC+condM)
        else
          RootAxialKond2H2O_pvr(N,L,NZ)  = DLYR3(L)**2/(Root1stAxialResist_rvr(N,L)+Root2ndAxialResist_rvr(N,L))     !plant size-scaled axial root conductance to H2O        
          RootResist4H2O_pvr(N,L,NZ)     = RootRadialResist_rvr(N,L)+Root1stAxialResist_rvr(N,L)+Root2ndAxialResist_rvr(N,L)
        endif
        SoilRootResist4H2O_pvr(N,L)    = SoilResist4H2O_rvr(N,L)+RootResist4H2O_pvr(N,L,NZ)
        CdH2ORootxSoil_pft(NZ)         = CdH2ORootxSoil_pft(NZ)+1.0_r8/SoilRootResist4H2O_pvr(N,L)   !assume different root layers work in parallel
      ENDIF
    enddo
  ENDDO D3880
  ! Conductance is PFT-size dependent; the generic AZMAX1 cutoff (1e-12)
  ! can erase a valid small plant's connection while layer uptake remains on.
  CdH2ORootxSoil = CdH2ORootxSoil_pft(NZ)
  call PrintInfo('end '//subname)
  end associate
  end subroutine CalcPlantHydroResistance

!----------------------------------------------------------------------------------------------------
  subroutine HandleBareSoil(NZ,TotalSoilPSIMPa_vr,FDMP,CanopyHeatCapacity)
  ! Description:
  ! Clear atmospheric/root-water fluxes when exchange is inactive. Retain the
  ! canopy temperature and actual biomass/water capacity in both plant modes;
  ! resetting either without an accepted heat transfer changes the inventory.
  use data_const_mod, only : spval => DAT_CONST_SPVAL
  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: TotalSoilPSIMPa_vr(JZ1)
  real(r8), intent(out):: FDMP
  real(r8), intent(in) :: CanopyHeatCapacity ! current biomass/water heat capacity [MJ K-1]
  character(len=*), parameter :: subname='HandleBareSoil'
  integer :: N,L
  real(r8) :: APSILT
  real(r8) :: CCPOLT
  real(r8) :: FTHRM
  real(r8) :: OSWT

! begin_execution
  call PrintInfo('beg '//subname)
  associate(                                                                      &
    AREA3                          => plt_site%AREA3                             ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    OrganOsmoPsi0pt_pft            => plt_ew%OrganOsmoPsi0pt_pft                 ,& !input  :Organ osmotic potential when canopy water potential = 0 MPa, [MPa]
    CanopyHeightLive_pft           => plt_morph%CanopyHeightLive_pft             ,& !input  :canopy height, [m]
    CanopyNonstOMConc_pft          => plt_biom%CanopyNonstOMConc_pft             ,& !input  :canopy nonstructural C+N+P concentration, [g gC-1]
    FracPARads2Canopy_pft          => plt_rad%FracPARads2Canopy_pft              ,& !input  :fraction of incoming PAR absorbed by canopy, [-]
    MaxSoilLays4Root_pft           => plt_morph%MaxSoilLays4Root_pft             ,& !input  :maximum soil layer number for all root axes,[-]
    Myco_pft                       => plt_morph%Myco_pft                         ,& !input  :mycorrhizal type (no or yes),[-]
    NGTopRootLayer_pft             => plt_morph%NGTopRootLayer_pft               ,& !input  :soil layer at planting depth, [-]
    NU                             => plt_site%NU                                ,& !input  :current soil surface layer number, [-]
    PSICanopyOsmo_pft              => plt_ew%PSICanopyOsmo_pft                   ,& !input  :canopy osmotic water potential, [Mpa]
    PSIRootOSMO_vr                 => plt_ew%PSIRootOSMO_vr                      ,& !input  :root osmotic water potential, [Mpa]
    PSIRootTurg_vr                 => plt_ew%PSIRootTurg_vr                      ,& !input  :root turgor water potential, [Mpa]
    EMS_Modify_Scalar_col          => plt_ew%EMS_Modify_Scalar_col               ,& !input  :canopy longwave radiation emissivity scalar
    RawIsoTCanopy2Atm_pft          => plt_ew%RawIsoTCanopy2Atm_pft               ,& !input  :canopy roughness height, [m]
    RootNonstOMConc_pvr            => plt_biom%RootNonstOMConc_pvr               ,& !input  :root/mycorrhizal nonstructural C+N+P concentration, [g gC-1]
    ShootElms_pft                  => plt_biom%ShootElms_pft                     ,& !input  :canopy shoot structural chemical element mass, [g d-2]
    SnowDepth                      => plt_ew%SnowDepth                           ,& !input  :snowpack depth, [m]
    MaxNumRootLays                 => plt_site%MaxNumRootLays                    ,& !input  :maximum root layer number,[-]
    TKS_vr                         => plt_ew%TKS_vr                              ,& !input  :mean annual soil temperature, [K]
    TKSnow                         => plt_ew%TKSnow                              ,& !input  :snow temperature, [K]
    TairK                          => plt_ew%TairK                               ,& !input  :air temperature, [K]
    ZERO                           => plt_site%ZERO                              ,& !input  :threshold zero for numerical stability, [-]
    PSICanopyTurg_pft              => plt_ew%PSICanopyTurg_pft                   ,& !input  :plant canopy turgor water potential, [MPa]
    PSICanopy_pft                  => plt_ew%PSICanopy_pft                       ,& !output :canopy total water potential, [Mpa]
    PSIRoot_pvr                    => plt_ew%PSIRoot_pvr                         ,& !output :root total water potential, [Mpa]
    TKC_pft                        => plt_ew%TKC_pft                             ,& !output :canopy temperature, [K]
    TKCanopy_pft                   => plt_ew%TKCanopy_pft                        ,& !output :canopy iteration temperature, [K]
    RPlantRootH2OUptk_pvr          => plt_ew%RPlantRootH2OUptk_pvr               ,& !output :whole population root water uptake, [m3 d-2 h-1]
    RootH2OUptkStress_pvr          => plt_ew%RootH2OUptkStress_pvr               ,& !output :pontential root water uptake rate, [m3 d-2 h-1]
    CanPStomaResistH2O_pft         => plt_photo%CanPStomaResistH2O_pft           ,& !output :canopy stomatal resistance, [h m-1]
    RawH2OCanopy2Atm_pft           => plt_photo%RawH2OCanopy2Atm_pft             ,& !output :canopy boundary layer resistance, [h m-1]
    RawHeatCanopy2Atm_pft          => plt_photo%RawHeatCanopy2Atm_pft            ,& !output :canopy sensible-heat boundary layer resistance, [h m-1]
    DeltaTKC_pft                   => plt_ew%DeltaTKC_pft                        ,& !output :change in canopy temperature, [K]
    CanopyEvapTransLHeat_pft       => plt_ew%CanopyEvapTransLHeat_pft            ,& !output :canopy latent heat flux, [MJ d-2 h-1]
    HeatStorCanopy_pft             => plt_ew%HeatStorCanopy_pft                  ,& !output :canopy storage heat flux, [MJ d-2 h-1]
    HeatXAir2PCan_pft              => plt_ew%HeatXAir2PCan_pft                   ,& !output :canopy sensible heat flux, [MJ d-2 h-1]
    LWRadCanopy_pft                => plt_rad%LWRadCanopy_pft                    ,& !output :canopy longwave radiation, [MJ d-2 h-1]
    QdewCanopy_pft                 => plt_ew%QdewCanopy_pft                      ,& !output :dew fall on to canopy, [m3 H2O d-2 h-1]
    RadNet2Canopy_pft              => plt_rad%RadNet2Canopy_pft                  ,& !output :canopy net radiation, [MJ d-2 h-1]
    Transpiration_pft              => plt_ew%Transpiration_pft                   ,& !output :canopy transpiration, [m2 d-2 h-1]
    VHeatCapCanopy_pft             => plt_ew%VHeatCapCanopy_pft                  ,& !output :canopy heat capacity, [MJ d-2 K-1]
    SnoSub2AirCanopy_pft           => plt_ew%SnoSub2AirCanopy_pft                ,& !output :canopy snow sublimation,[m3 d-2 h-1]
    VapXAir2CanopyLiq_pft          => plt_ew%VapXAir2CanopyLiq_pft               ,& !output :canopy evaporation, [m2 d-2 h-1]
    VapXAir2Canopy_pft             => plt_ew%VapXAir2Canopy_pft                   & !output :canopy evaporation+sublimation, [m2 d-2 h-1]
  )
  RadNet2Canopy_pft(NZ)        = 0.0_r8
  CanopyEvapTransLHeat_pft(NZ) = 0.0_r8
  HeatXAir2PCan_pft(NZ)        = 0.0_r8
  HeatStorCanopy_pft(NZ)       = 0.0_r8
  VapXAir2Canopy_pft(NZ)       = 0.0_r8
  SnoSub2AirCanopy_pft(NZ)     = 0.0_r8
  VapXAir2CanopyLiq_pft(NZ)    = 0.0_r8
  Transpiration_pft(NZ)        = 0.0_r8
  QdewCanopy_pft(NZ)           = 0.0_r8
  
  TKCanopy_pft(NZ)=TKC_pft(NZ)

  FTHRM               = EMMC*stefboltz_const*FracPARads2Canopy_pft(NZ)*AREA3(NU)*EMS_Modify_Scalar_col
  LWRadCanopy_pft(NZ) = FTHRM*TKC_pft(NZ)**4
  PSICanopy_pft(NZ)   = TotalSoilPSIMPa_vr(NGTopRootLayer_pft(NZ))
  CCPOLT              = CanopyNonstOMConc_pft(NZ)

  call update_osmo_turg_pressure(PSICanopy_pft(NZ),CCPOLT,OrganOsmoPsi0pt_pft(NZ),TKC_pft(NZ)&
    ,PSICanopyOsmo_pft(NZ),PSICanopyTurg_pft(NZ),FDMP)
  !
  CanPStomaResistH2O_pft(NZ) = CalcStomataResist4H2O(NZ)
  !
  RawH2OCanopy2Atm_pft(NZ)  = RawIsoTCanopy2Atm_pft(NZ)
  RawHeatCanopy2Atm_pft(NZ) = RawH2OCanopy2Atm_pft(NZ)/BndlResistH2OToHeatRatio
  VHeatCapCanopy_pft(NZ)    = CanopyHeatCapacity
  DeltaTKC_pft(NZ)          = TKC_pft(NZ)-TairK

  DO N=1,Myco_pft(NZ)
    DO  L=NU,MaxNumRootLays
      PSIRoot_pvr(N,L,NZ)              = TotalSoilPSIMPa_vr(L)
      if(ldo_sp_mode)then
        CCPOLT =0.4_r8
      else
        CCPOLT = RootNonstOMConc_pvr(N,L,NZ)
      endif
      RPlantRootH2OUptk_pvr(N,L,NZ) = 0.0_r8
      RootH2OUptkStress_pvr(N,L,NZ)   = 0._r8
      call update_osmo_turg_pressure(PSIRoot_pvr(N,L,NZ),CCPOLT,OrganOsmoPsi0pt_pft(NZ),TKS_vr(L),&
        PSIRootOSMO_vr(N,L,NZ),PSIRootTurg_vr(N,L,NZ))
    enddo
  ENDDO
  
  end associate
  call PrintInfo('end '//subname)
  end subroutine HandleBareSoil

!----------------------------------------------------------------------------------------------------
  subroutine UpdatePlantWaterVars(NZ,CdHeatCanopyAir,CanopyMassC4H2OStorage,TotalSoilPSIMPa_vr,SoilResist4H2O_rvr,SoilRootResist4H2O_pvr,&
    TKCX,VHeatCapCanopyPrev_pft,PrecpHeatbyCanopy,cumPRootH2OUptake,CumPlantHeatLoss2Soil,HeatEvapSens,SoiLayerHasRoot_rvr,DeadHeatCapacity)
  ! Description:
  !Update canopy heat and water states after the numerical iterations.
  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: CdHeatCanopyAir
  real(r8), intent(in) :: CanopyMassC4H2OStorage         ![gC/m2]
  real(r8), intent(in) :: TotalSoilPSIMPa_vr(JZ1)    !Elevation adjusted soil water potential [MPa]
  real(r8), intent(in) :: SoilResist4H2O_rvr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: SoilRootResist4H2O_pvr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: TKCX,VHeatCapCanopyPrev_pft
  real(r8), intent(in) :: PrecpHeatbyCanopy
  real(r8), intent(in) :: cumPRootH2OUptake
  real(r8), intent(in) :: HeatEvapSens                     !sensible heat associated with evaporation [MJ/h]
  real(r8), intent(in) :: CumPlantHeatLoss2Soil
  logical , intent(in) :: SoiLayerHasRoot_rvr(pltpar%jroots,JZ1)  !indicator of root prescence

  real(r8), intent(in) :: DeadHeatCapacity ! standing-dead dry-material capacity [MJ/K]
  character(len=*), parameter :: subname='UpdatePlantWaterVars'
  real(r8) :: CCPOLT
  real(r8) :: OSWT
  integer :: N,L
  call PrintInfo('beg '//subname)
  associate(                                                          &
    OrganOsmoPsi0pt_pft       => plt_ew%OrganOsmoPsi0pt_pft          ,& !input  :Organ osmotic potential when canopy water potential = 0 MPa, [MPa]
    MaxSoilLays4Root_pft      => plt_morph%MaxSoilLays4Root_pft      ,& !input  :maximum soil layer number for all root axes,[-]
    Myco_pft                  => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    PSICanopy_pft             => plt_ew%PSICanopy_pft                ,& !input  :canopy total water potential, [Mpa]
    PSIRootOSMO_vr            => plt_ew%PSIRootOSMO_vr               ,& !input  :root osmotic water potential, [Mpa]
    PSIRootTurg_vr            => plt_ew%PSIRootTurg_vr               ,& !input  :root turgor water potential, [Mpa]
    RainIntcptByCanopy_pft    => plt_ew%RainIntcptByCanopy_pft       ,& !input  :water flux into canopy, [m3 d-2 h-1]
    SnowIntcptByCanopy_pft    => plt_ew%SnowIntcptByCanopy_pft       ,& !input  :Snow intercepted by canopy, [m3 d-2 h-1]
    RootNonstOMConc_pvr       => plt_biom%RootNonstOMConc_pvr        ,& !input  :root/mycorrhizal nonstructural C+N+P concentration, [g gC-1]
    TKCanopy_pft              => plt_ew%TKCanopy_pft                 ,& !input  :canopy temperature, [K]
    RootResist4H2O_pvr        => plt_ew%RootResist4H2O_pvr           ,& !input  :total root (axial+radial) resistance for water uptake
    TKS_vr                    => plt_ew%TKS_vr                       ,& !input  :mean annual soil temperature, [K]
    TairK                     => plt_ew%TairK                        ,& !input  :air temperature, [K]
    Transpiration_pft         => plt_ew%Transpiration_pft            ,& !input  :canopy transpiration, [m3 d-2 h-1]
    VapXAir2Canopy_pft        => plt_ew%VapXAir2Canopy_pft           ,& !input  :canopy evaporation+sublimation, [m3 d-2 h-1]
    VapXAir2CanopyLiq_pft     => plt_ew%VapXAir2CanopyLiq_pft        ,& !input  :canopy evaporation, [m3 d-2 h-1]
    SnoSub2AirCanopy_pft      => plt_ew%SnoSub2AirCanopy_pft         ,& !input  :canopy snow sublimation,[m3 d-2 h-1]
    CanopyBiomWater_pft       => plt_ew%CanopyBiomWater_pft          ,& !inoput :canopy water content, [m3 d-2]
    WatHeldOnCanopy_pft       => plt_ew%WatHeldOnCanopy_pft          ,& !inoput :canopy surface water content, [m3 d-2]
    SnowOnCanopy_pft          => plt_ew%SnowOnCanopy_pft             ,& !input  :snow water equivalent held on canopy, [m3 d-2]    
    PSIRoot_pvr               => plt_ew%PSIRoot_pvr                  ,& !output :root total water potential, [Mpa]
    VHeatCapCanopy_pft        => plt_ew%VHeatCapCanopy_pft           ,& !output :canopy heat capacity, [MJ d-2 K-1]
    HeatStorCanopy_pft        => plt_ew%HeatStorCanopy_pft           ,& !output :canopy storage heat flux, [MJ d-2 h-1]
    HeatXAir2PCan_pft         => plt_ew%HeatXAir2PCan_pft             & !output :canopy sensible heat flux, [MJ d-2 h-1]
  )
  !
  !     CANOPY SURFACE WATER STORAGE, SENSIBLE AND STORAGE HEAT FLUXES
  !     (NOT EXPLICITLY CALCULATED IN CONVERGENCE SOLUTION)
  !
  !     VOLWP,WatHeldOnCanopy_pft=water volume in canopy,on canopy surfaces
  !     HeatXAir2PCan_pft=canopy sensible,storage heat fluxes
  !     HeatStorCanopy_pft=Radiation-latent heat-sens, i.e. (>0) residual heat from air to canopy
  !     VHeatCapCanopyPrev_pft,VHeatCapCanopy_pft=previous,current canopy heat capacity
  !     CdHeatCanopyAir=canopy sensible heat conductance
  !     HeatEvapSens=convective heat flux from latent heat flux (>0)
  !     PrecpHeatbyCanopy=convective heat flux from precip to canopy
  !     cumPRootH2OUptake < 0., add water to canopy 
  !     Transpiration_pft < 0., lost from canopy, [m3 d-2 h-1]
  CanopyBiomWater_pft(NZ) = CanopyBiomWater_pft(NZ)+Transpiration_pft(NZ)-cumPRootH2OUptake
  WatHeldOnCanopy_pft(NZ) = WatHeldOnCanopy_pft(NZ)+RainIntcptByCanopy_pft(NZ)+VapXAir2CanopyLiq_pft(NZ)
  SnowOnCanopy_pft(NZ)    = SnowOnCanopy_pft(NZ)+SnowIntcptByCanopy_pft(NZ)+SnoSub2AirCanopy_pft(NZ)
  VHeatCapCanopy_pft(NZ)  = cpo*gOC_to_m3_OM(CanopyMassC4H2OStorage)+DeadHeatCapacity+cpw*(WatHeldOnCanopy_pft(NZ)+CanopyBiomWater_pft(NZ)) &
    +cps*SnowOnCanopy_pft(NZ)
  HeatXAir2PCan_pft(NZ)   = CdHeatCanopyAir*(TairK-TKCanopy_pft(NZ))
  HeatStorCanopy_pft(NZ)  = TKCX*VHeatCapCanopyPrev_pft-TKCanopy_pft(NZ)*VHeatCapCanopy_pft(NZ)+HeatEvapSens+PrecpHeatbyCanopy
  !
  !     ROOT TOTAL, OSMOTIC AND TURGOR WATER POTENTIALS
  !
  !     PSIRoot_pvr,PSICanopy_pft=root,canopy total water potential
  !     TotalSoilPSIMPa_vr=total soil water potential PSIST adjusted for surf elevn
  !     SoilResist4H2O_rvr,SoilRootResist4H2O_pvr,RootResist=soil,soil+root,root radial+axial resistance
  !     PSIRootOSMO_vr,PSIRootTurg_vr=root osmotic,turgor water potential
  !     OrganOsmoPsi0pt_pft=osmotic potential at PSIRoot_pvr=0 from PFT file

  !compute root pressure assuming zero water storage capacity in root, cf. Eq. (36) Grant (1998), Ecological modelling.
  D4505: DO N=1,Myco_pft(NZ)
    D4510: DO L=NU,MaxSoilLays4Root_pft(NZ)
      IF(SoiLayerHasRoot_rvr(N,L))THEN
        !assuming steady-state of root water status, so that (psi_c-psi_root)/r_root=(psi_root-psi_soil)/r_soil
        !
        PSIRoot_pvr(N,L,NZ)=AZMIN1((TotalSoilPSIMPa_vr(L)*RootResist4H2O_pvr(N,L,NZ) &
          +PSICanopy_pft(NZ)*SoilResist4H2O_rvr(N,L))/SoilRootResist4H2O_pvr(N,L))
      ELSE
        PSIRoot_pvr(N,L,NZ)=TotalSoilPSIMPa_vr(L)
      ENDIF
      !obtain total reserve
      if(ldo_sp_mode)then
        CCPOLT=0.4_r8
      else
        CCPOLT=RootNonstOMConc_pvr(N,L,NZ)
      endif

      CALL update_osmo_turg_pressure(PSIRoot_pvr(N,L,NZ),CCPOLT,OrganOsmoPsi0pt_pft(NZ),TKS_vr(L),&
        PSIRootOSMO_vr(N,L,NZ),PSIRootTurg_vr(N,L,NZ))

    ENDDO D4510
  ENDDO D4505
  end associate
  call PrintInfo('end '//subname)
  end subroutine UpdatePlantWaterVars

!----------------------------------------------------------------------------------------------------
  subroutine SetCanopyGrowthFuncs(yearIJ,NZ)

  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ

  character(len=*), parameter :: subname='SetCanopyGrowthFuncs'
  real(r8) :: ACTV,RTK,STK,TKGO,TKSO
  integer :: L
  associate(                                               &
    MainBranchNum_pft   => plt_morph%MainBranchNum_pft    ,& !input  :number of main branch,[-]
    MaxNumRootLays      => plt_site%MaxNumRootLays        ,& !input  :maximum root layer number,[-]
    NU                  => plt_site%NU                    ,& !input  :current soil surface layer number, [-]
    PSICanopy_pft       => plt_ew%PSICanopy_pft           ,& !input  :canopy total water potential, [Mpa]
    TCChill4Seed_pft    => plt_pheno%TCChill4Seed_pft     ,& !input  :temperature below which seed set is adversely affected, [oC]
    TKC_pft             => plt_ew%TKC_pft                 ,& !input  :canopy temperature, [K]
    TKS_vr              => plt_ew%TKS_vr                  ,& !input  :mean annual soil temperature, [K]
    TdegCCanopy_pft     => plt_ew%TdegCCanopy_pft         ,& !input  :canopy temperature, [oC]
    TempOffset_pft      => plt_pheno%TempOffset_pft       ,& !input  :adjustment of Arhhenius curves for plant thermal acclimation, [oC]
    iPlantCalendar_brch => plt_pheno%iPlantCalendar_brch  ,& !input  :plant growth stage, [-]
    ChillHours_pft      => plt_photo%ChillHours_pft       ,& !inoput :chilling effect on CO2 fixation, [-]
    PSICanPDailyMin_pft => plt_ew%PSICanPDailyMin_pft     ,& !inoput :minimum daily canopy water potential, [MPa]
    TKGroth_pft         => plt_pheno%TKGroth_pft          ,& !output :canopy growth temperature, [K]
    TCGroth_pft         => plt_pheno%TCGroth_pft          ,& !output :canopy growth temperature, [oC]
    fTCanopyGroth_pft   => plt_pheno%fTCanopyGroth_pft    ,& !output :canopy temperature growth function, [-]
    fTgrowRootP_vr      => plt_pheno%fTgrowRootP_vr        & !output :root layer temperature growth functiom, [-]
  )
  !
  !     SET CANOPY GROWTH TEMPERATURE FROM SOIL SURFACE
  !     OR CANOPY TEMPERATURE DEPENDING ON GROWTH STAGE
  !
  call PrintInfo('beg '//subname)
  IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).EQ.0)THEN
    !before seed emergence
    TKGroth_pft(NZ)=TKS_vr(NU)
  ELSE
    !after seed emergence
    TKGroth_pft(NZ)=TKC_pft(NZ)
  ENDIF
  TCGroth_pft(NZ)=real_truncate(units%Kelvin2Celcius(TKGroth_pft(NZ)),1.e-3_r8)
  !
  !     ARRHENIUS FUNCTION FOR CANOPY AND ROOT GROWTH WITH OFFSET
  !     FOR ZONE OF THERMAL ADAPTATION ENTERED IN 'READQ'
  !
  !     TKGroth_pft,TKGO=canopy temperature,canopy temp used in Arrhenius eqn
  !     TKS_vr,TKSO=soil temperature,soil temp used in Arrhenius eqn
  !     TempOffset_pft=shift in Arrhenius curve for thermal adaptation
  !     fTCanopyGroth_pft,fTgrowRootP_vr=temperature function for canopy,root growth (25 oC =1)
  !     RGAS,710.0=gas constant,enthalpy
  !     62500,197500,222500=energy of activn,high,low temp inactivn(KJ mol-1)
  !     PSICanPDailyMin=minimum daily canopy water potential
  !
  TKGO                  = real_truncate(TKGroth_pft(NZ)+TempOffset_pft(NZ),1.e-3_r8)
  fTCanopyGroth_pft(NZ) = calc_canopy_grow_tempf(TKGO)

  D100: DO L=NU,MaxNumRootLays
    TKSO                 = real_truncate(TKS_vr(L)+TempOffset_pft(NZ),1.e-3_r8)
    fTgrowRootP_vr(L,NZ) = calc_root_grow_tempf(TKSO)
  ENDDO D100

  PSICanPDailyMin_pft(NZ)=AMIN1(PSICanPDailyMin_pft(NZ),PSICanopy_pft(NZ))
  !
  !     DIURNAL CHILLING
  !
  !     TCChill4Seed_pft=chilling temperature from PFT file
  !     CHILL=accumulated chilling hours used to limit CO2 fixn in stomate.f
  !

  IF(TdegCCanopy_pft(NZ).LT.TCChill4Seed_pft(NZ))THEN
    ChillHours_pft(NZ)=AMIN1(24.0_r8,ChillHours_pft(NZ)+1.0_r8)
  ELSE
    ChillHours_pft(NZ)=AZMAX1(ChillHours_pft(NZ)-1.0_r8)
  ENDIF

  call PrintInfo('end '//subname)
  end associate
  end subroutine SetCanopyGrowthFuncs
  ![tail]
end module UptakesMod
