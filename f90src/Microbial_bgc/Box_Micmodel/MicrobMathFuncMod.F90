module MicrobMathFuncMod
  use EcosimConst
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use EcoSIMSolverPar,      only: dts_gas, NPH, NPT
  use MicFLuxTypeMod,       only: micfluxtype
  use MicForcTypeMod,       only: micforctype
  use MicStateTraitTypeMod, only: micsttype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: PrintInfo
  use ElmIDMod,             only: ielmc, ielmn, iLbiom_kinetic, iLbiom_struct
  use MicrobeDiagTypes,     only: Cumlate_Flux_Diag_type, Microbe_Diag_type, &
                                  Microbe_Flux_type, Microbe_State_type
  use minimathmod,          only: AZMAX1, safe_adb, fixEXConsumpFlux, &
                                  real_truncate
  use NitroPars,            only: BIOS, FMN, ORAD, RMOM
  implicit none

  character(len=*), parameter, private :: mod_filename=&
  __FILE__

  real(r8), parameter, private :: ZERO=1.0E-15_r8

  contains
!------------------------------------------------------------------------


  subroutine LimitReserveNutrientTransfers(Reserve,Supply,GrowthTransfer,MineralTransfer)
  !Growth and mineralization share one reserve donor. Call immediately before
  !updating biomass so transfers also respect any intervening priming exchange.
  !Supply includes organic uptake, fixation, recycling and underlying-soil fluxes.
  !Positive mineral uptake has already passed donor/competition limits: retain it.
  !Reduce release proportionally across compartments; never invent extra uptake.
  implicit none
  real(r8), intent(in) :: Reserve,Supply
  real(r8), intent(inout) :: GrowthTransfer(:),MineralTransfer(:)
  real(r8) :: Available,Growth,Release,Scale

  Available=MAX(0._r8,Reserve+Supply+SUM(MAX(0._r8,MineralTransfer)))
  Growth=SUM(GrowthTransfer)
  IF(Growth.GT.Available)THEN
    GrowthTransfer=GrowthTransfer*(Available/Growth)
  ENDIF
  Available=MAX(0._r8,Available-SUM(GrowthTransfer))
  Release=SUM(MAX(0._r8,-MineralTransfer))
  IF(Release.GT.Available)THEN
    Scale=Available/Release
    WHERE(MineralTransfer.LT.0._r8)MineralTransfer=MineralTransfer*Scale
  ENDIF
  end subroutine LimitReserveNutrientTransfers
!------------------------------------------------------------------------

  subroutine LitterSoilNutrientUptake(Nutrient,Demand,Capacity,KM,CMin, &
    NonbandShare,PrevBandDemand,BiomShare,micfor,Potential,Uptake)
  !
  !Description:  
  !Paired uptake from the soil beneath litter. Positive flux immobilizes nutrient;
  !negative demand releases it. No compartment can borrow the other's allocation.
  !
  implicit none
  integer, intent(in) :: Nutrient !1 NH4, 2 NO3, 3 H2PO4, 4 HPO4
  real(r8), intent(in) :: Demand,Capacity,KM,CMin,NonbandShare,PrevBandDemand,BiomShare
  type(micforctype), intent(in) :: micfor
  real(r8), intent(out) :: Potential(2),Uptake(2) !Nonband, fertilizer band
  real(r8) :: Frac(2),Share(2),ExcessConc(2),Available(2),FracSum

  Potential = 0._r8
  Uptake    = 0._r8
  Frac      = MAX(0._r8,micfor%LitrSoilNutFrac(:,Nutrient))
  FracSum   = SUM(Frac)
  IF(FracSum.LE.0._r8)RETURN
  Frac=Frac/FracSum
  IF(Demand.LE.0._r8)THEN
    Uptake=Demand*Frac
    RETURN
  ENDIF

  IF(micfor%VOLWU.LE.micfor%ZEROS)RETURN
  Share(1)=MIN(1._r8,MAX(0._r8,NonbandShare))
  IF(micfor%LitrSoilBandDemandPrev(Nutrient).GT.micfor%ZEROS)THEN
    Share(2)=MIN(1._r8,MAX(FMN,PrevBandDemand/micfor%LitrSoilBandDemandPrev(Nutrient)))
  ELSE
    Share(2)=MIN(1._r8,MAX(FMN,BiomShare))
  ENDIF

  ExcessConc = MAX(0._r8,micfor%LitrSoilNutConc(:,Nutrient)-CMin)
  Potential  = MIN(Demand,MAX(0._r8,Capacity))*Frac*ExcessConc/(ExcessConc+KM)
  Available  = MAX(0._r8,micfor%LitrSoilNutPool(:,Nutrient)-CMin*micfor%VOLWU*Frac)
  Uptake     = MIN(Potential,Share*Available)
  end subroutine LitterSoilNutrientUptake
!------------------------------------------------------------------------

  subroutine ReserveDenitrifMaintenance(MaintDemand,PrimaryResp,DenitResp, &
    AerobicRespFrac,DenitRespFrac,FracDenitResp4Maint,MaintDeficit)
  !
  !Description:
  !Maintenance demand retains its existing aerobic respiratory-C basis.
  !The yield relation (1/eta-1)=energy yield/growth energy cost converts
  !denitrification respiration to that basis using EO2X and ENOX, which
  !share EOMC. This leaves the primary pathway's maintenance convention intact.
  implicit none
  real(r8), intent(in) :: MaintDemand,PrimaryResp,DenitResp
  real(r8), intent(in) :: AerobicRespFrac,DenitRespFrac
  real(r8), intent(out) :: FracDenitResp4Maint !Fraction of denitrification respiration reserved for maintenance
  real(r8), optional, intent(out) :: MaintDeficit !Remaining demand in primary respiratory-C units
  real(r8) :: RemainingDemand,DenitMaintYield,DenitMaintSupply

  RemainingDemand     = AZMAX1(MaintDemand-PrimaryResp)
  FracDenitResp4Maint = 0._r8

  IF(RemainingDemand.GT.0._r8 .AND. DenitResp.GT.0._r8)THEN
    IF(AerobicRespFrac.GT.0._r8 .AND. AerobicRespFrac.LT.1._r8 &
      .AND. DenitRespFrac.GT.0._r8 .AND. DenitRespFrac.LT.1._r8)THEN
      DenitMaintYield  = (1._r8/DenitRespFrac-1._r8)/(1._r8/AerobicRespFrac-1._r8)
      DenitMaintSupply = DenitResp*DenitMaintYield
      IF(DenitMaintSupply.GE.RemainingDemand)THEN
        FracDenitResp4Maint = RemainingDemand/DenitMaintSupply
        RemainingDemand     = 0._r8
      ELSE
        FracDenitResp4Maint = 1._r8
        RemainingDemand     = AZMAX1(RemainingDemand-DenitMaintSupply)
      ENDIF
    ENDIF
  ENDIF
  IF(PRESENT(MaintDeficit))MaintDeficit=RemainingDemand
  end subroutine ReserveDenitrifMaintenance
!------------------------------------------------------------------------

  subroutine MicrobPhysTempFun(TKSO, TSensGrowth, TSensMaintR)
  !
  !the physiological temperature dependence of microbes
  implicit none
  real(r8), intent(in) :: TKSO
  real(r8), intent(out):: TSensGrowth !temperature sensitivity for growth respiration
  real(r8), intent(out):: TSensMaintR !temperature sensitivity for maintenance respiration
  real(r8) :: RTK,STK
  real(r8) :: ACTV,ACTVM

  RTK=RGASC*TKSO
  STK=710.0_r8*TKSO
  ACTV=1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  TSensGrowth=EXP(25.229_r8-62500._r8/RTK)/ACTV

  ACTVM=1+EXP((195000._r8-STK)/RTK)+EXP((STK-232500._r8)/RTK)
  TSensMaintR=EXP(25.214_r8-62500._r8/RTK)/ACTVM

  end subroutine MicrobPhysTempFun
!------------------------------------------------------------------------


  pure function TranspBasedsubstrateUptake(S_conc,diffusc, KM, V_max, zeros)result(uptake)
  !
  !transport based substrate uptake
  !assuming balance between MM-based uptake and supply by diffusion
  !q=V_max*s/(KM+s)=D(S_ext-S)
  !solve for q.
  implicit none
  real(r8), intent(in) :: S_conc   !substrate concentration
  real(r8), intent(in) :: diffusc !diffusion coefficient
  real(r8), intent(in) :: KM      !half saturation parameter
  real(r8), intent(in) :: V_max   !maximum uptake rate
  real(r8), optional, intent(in) :: zeros !threshold for active uptake
  real(r8) :: uptake
  real(r8) :: X,B,C
  real(r8) :: zero1

  if(present(zeros))then
    zero1=zeros
  else
    zero1=zero
  endif

  !obtain uptake flux
  X=diffusc*S_conc
  IF(X.GT.ZERO1)THEN
    B=-V_max-diffusc*KM-X
    C=X*V_max
    uptake=(-B-SQRT(B*B-4.0_r8*C))/2.0_r8
  ELSE
    uptake=0.0_r8
  ENDIF
  end function TranspBasedsubstrateUptake

!------------------------------------------------------------------------------------------

  subroutine AerobicHeterO2Uptake(I,J,NGL,N,K,FOXYX,OXKX,micfor,micstt,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: OXKX
  real(r8), intent(in) :: FOXYX  !preallocated O2 for uptake
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type),intent(inout) :: nmics
  type(micfluxtype), intent(inout) :: micflx
  integer  :: M,MX
  real(r8) :: COXYS1,DIFOX
  real(r8) :: B,C,O2AquaDiffusvity1
  real(r8) :: OXYG1,OXYS1
  real(r8) :: RUPMAX
  real(r8) :: ROXYFX       !gas dissolution flux
  real(r8) :: ROXYLX       !aqueous flux transport flux
  real(r8) :: RRADO,RMPOX
  real(r8) :: ROXDFQ  !gas dissolution
  real(r8) :: THETW1,VOLWOX
  real(r8) :: VOLPOX
  real(r8) :: X
  real(r8) :: VOLWPM,VOLOXM
  ! begin_execution
  associate(                                                 &
    OxyLimterHeter         => nmics%OxyLimterHeter,          & !Actual/potential O2 uptake ratio; 1 means no O2 restriction [-]
    OMActHeter             => nmics%OMActHeter,              & !Active microbial C biomass by heterotrophic guild and complex K
    RO2UptkHeter           => nmicf%RO2UptkHeter,            & !Realized total O2 uptake by heterotrophic guild and complex
    RespGrossHeter         => nmicf%RespGrossHeter,          & !Gross respiration C equivalent from the primary heterotrophic pathway; not referenced here
    RO2DmndHeter           => nmicf%RO2DmndHeter,            & !Total guild O2 demand before O2 limitation
    O2_irrig_conc          => micfor%O2_irrig_conc,          & !Dissolved O2 concentration in irrigation water
    O2_rain_conc           => micfor%O2_rain_conc,           & !Dissolved O2 concentration in rainwater
    COXYE                  => micfor%COXYE,                  & !Atmospheric gas-phase O2 concentration
    RO2GasXchangePrev      => micfor%RO2GasXchangePrev,      & !Previous-hour gaseous O2 exchange; negated when applied as a supply
    RO2AquaXchangePrev     => micfor%RO2AquaXchangePrev,     & !Previous-hour aqueous O2 exchange; negated when applied as a supply
    Irrig2LitRSurf_col     => micfor%Irrig2LitRSurf_col,     & !Irrigation water input to surface litter, carrying dissolved O2
    Rain2LitRSurf          => micfor%Rain2LitRSurf,          & !Rainwater input to surface litter, carrying dissolved O2
    litrm                  => micfor%litrm,                  & !True for the surface litter layer
    O2AquaDiffusvity       => micfor%O2AquaDiffusvity,       & !Aqueous O2 diffusivity before transport-substep scaling
    VLSoilPoreMicP         => micfor%VLSoilPoreMicP,         & !Layer micropore volume used in water and aerobic-uptake calculations
    VLSoilMicP             => micfor%VLSoilMicP,             & !Bulk volume associated with the layer micropore domain
    ZERO                   => micfor%ZERO,                   & !Small dimensionless or concentration threshold used by the routine
    ZEROS                  => micfor%ZEROS,                  & !Small mass or flux threshold used by the routine
    VLsoiAirPM             => micfor%VLsoiAirPM,             & !Soil air volume at each outer transport substep M
    VLWatMicP              => micfor%VLWatMicP,              & !Layer micropore water volume used for nutrient donor thresholds; not referenced here
    VLWatMicPM             => micfor%VLWatMicPM,             & !Micropore water volume at each outer transport substep M
    THETPM                 => micfor%THETPM,                 & !Air-filled soil pore fraction at each outer transport substep M [-]
    DiffusivitySolutEff    => micfor%DiffusivitySolutEff,    & !Gas-water exchange coefficient at each transport substep M
    FILM                   => micfor%FILM,                   & !Water-film thickness for microbial O2 diffusion at transport substep M
    TortMicPM              => micfor%TortMicPM,              & !Aqueous diffusion tortuosity factor at transport substep M [-]
    OXYG                   => micstt%OXYG,                   & !Gas-phase O2 donor pool
    OXYS                   => micstt%OXYS,                   & !Dissolved O2 donor pool
    COXYS                  => micstt%COXYS,                  & !Dissolved O2 concentration; not referenced here
    O2GSolubility          => micstt%O2GSolubility,          & !Equilibrium aqueous-to-gas O2 concentration ratio [-]
    COXYG                  => micstt%COXYG,                  & !Soil gas-phase O2 concentration
    RNO2DmndReduxSoilHeter => micflx%RNO2DmndReduxSoilHeter, & !Potential nonband NO2-N reduction demand used for competition; not referenced here
    RNO2DmndReduxBandHeter => micflx%RNO2DmndReduxBandHeter, & !Potential fertilizer-band NO2-N reduction demand used for competition; not referenced here
    REcoUptkSoilO2M        => micflx%REcoUptkSoilO2M         & !Accumulated microbial O2 uptake in each outer transport substep M
  )

  IF(RO2DmndHeter(NGL,K).GT.ZEROS .AND. FOXYX.GT.ZERO)THEN
    IF(.not.litrm .OR. VLSoilPoreMicP.GT.ZEROS)THEN
      !
      !write(*,*)'MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC'
      !     POPULATION
      !
      RUPMAX            = RO2DmndHeter(NGL,K)*dts_gas
      ROXYFX            = -RO2GasXchangePrev*dts_gas*FOXYX
      O2AquaDiffusvity1 = O2AquaDiffusvity*dts_gas

      IF(.not.litrm)THEN
        OXYG1  = OXYG*FOXYX
        ROXYLX = -RO2AquaXchangePrev*dts_gas*FOXYX
      ELSE
        !litter layer
        OXYG1  = COXYG*VLsoiAirPM(1)*FOXYX
        ROXYLX = -(RO2AquaXchangePrev+Rain2LitRSurf*O2_rain_conc+Irrig2LitRSurf_col*O2_irrig_conc)*dts_gas*FOXYX
      ENDIF
      OXYS1=OXYS*FOXYX

      if(OXYS1 <= 0._r8 .and. ROXYLX > 0._r8)ROXYLX=0._r8
      IF(OXYG1 <=0._r8 .and. ROXYFX>0._R8)ROXYFX=0._r8
!
      !write(*,*)'O2 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP'
!     TO MAINTAIN AQUEOUS O2 CONCENTRATION DURING REDUCTION
!
      D420: DO M=1,NPH
        !
        !     ACTUAL REDUCTION OF AQUEOUS BY AEROBES CALCULATED
        !     FROM MASS FLOW PLUS DIFFUSION = ACTIVE UPTAKE
        !     COUPLED WITH DISSOLUTION OF GASEOUS O2 DURING REDUCTION
        !     OF AQUEOUS O2 FROM DISSOLUTION RATE CONSTANT 'DiffusivitySolutEff'
        !     CALCULATED IN 'WATSUB'
        !
        !     VLWatMicPM,VLsoiAirPM,VLSoilPoreMicP=water, air and total volumes
        !     ORAD=microbial radius,FILM=water film thickness
        !     DIFOX=aqueous O2 diffusion, TortMicPM=tortuosity
        !     BIOS=microbial number, OMA=active biomass
        !     O2GSolubility=O2 solubility, OXKX=Km for O2 uptake
        !     OXYS,COXYS=aqueous O2 amount, concentration
        !     OXYG,COXYG=gaseous O2 amount, concentration
        !     RMPOX,REcoUptkSoilO2M=O2 uptake
        !
        THETW1 = AZMAX1(safe_adb(VLWatMicPM(M),VLSoilMicP))
        RRADO  = ORAD*(FILM(M)+ORAD)/FILM(M)
        DIFOX  = TortMicPM(M)*O2AquaDiffusvity1*12.57_r8*BIOS*OMActHeter(NGL,K)*RRADO
        VOLWOX = VLWatMicPM(M)*O2GSolubility
        VOLPOX = VLsoiAirPM(M)
        VOLWPM = VOLWOX+VOLPOX
        VOLOXM = VLWatMicPM(M)*FOXYX
        D425: DO MX=1,NPT
          call fixEXConsumpFlux(OXYG1,ROXYFX)
          call fixEXConsumpFlux(OXYS1,ROXYLX)

          COXYS1 = AMIN1(COXYE*O2GSolubility,safe_adb(OXYS1,VOLOXM))

          !obtain O2 uptake flux
          if(OXYS1<=ZEROS)then
            RMPOX=0.0_r8
          else
            RMPOX=TranspBasedsubstrateUptake(COXYS1,DIFOX, OXKX, RUPMAX, ZEROS)
          endif

          !Credit only O2 present in the aqueous donor during this substep.
          !O2 supplied by dissolution below is available in the next substep.
          RMPOX = MIN(MAX(0._r8,RMPOX),MAX(0._r8,OXYS1))
          OXYS1 = OXYS1-RMPOX

          !apply dissolution-volatilization
          IF(THETPM(M).GT.AirFillPore_Min.AND.VOLPOX.GT.ZEROS)THEN
            ROXDFQ=DiffusivitySolutEff(M)*(AMAX1(ZEROS,OXYG1)*VOLWOX-OXYS1*VOLPOX)/VOLWPM
            ROXDFQ=AMAX1(AMIN1(ROXDFQ,OXYG1),-OXYS1)
          ELSE
            ROXDFQ=0.0_r8
          ENDIF
          OXYG1 = OXYG1-ROXDFQ
          OXYS1 = OXYS1+ROXDFQ
          !accumulate upatke
          RO2UptkHeter(NGL,K) = RO2UptkHeter(NGL,K)+RMPOX
          REcoUptkSoilO2M(M)     = REcoUptkSoilO2M(M)+RMPOX
        ENDDO D425

      ENDDO D420
      !
      !     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (OxyLimterHeter)
      !
      !     OxyLimterHeter=ratio of O2-limited to O2-unlimited uptake
      !     RVMX4,RVNHB,RNO2DmndReduxSoilHeter,RNO2DmndReduxBandHeter=NH3,NO2 oxidation in non-band, band
      !
      OxyLimterHeter(NGL,K)=AMIN1(1.0_r8,AZMAX1(RO2UptkHeter(NGL,K)/RO2DmndHeter(NGL,K)))

    ELSE
      RO2UptkHeter(NGL,K)   = RO2DmndHeter(NGL,K)
      OxyLimterHeter(NGL,K) = 1.0_r8
    ENDIF
  ELSE
    RO2UptkHeter(NGL,K)   = 0.0_r8
    OxyLimterHeter(NGL,K) = 1.0_r8
  ENDIF
  !
  !     RespGrossHeter,RGOMP=O2-limited, O2-unlimited respiration
  !     RCO2ProdHeter,RAcetateProdHeter,RCH4ProdHeter,RH2ProdHeter=CO2,acetate,CH4,H2 production from RespGrossHeter
  !     RO2Uptk4RespHeter=O2-limited O2 uptake
  !     RSMetaOxidSoilAutor,RSMetaOxidBandAutor=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation
  !

  end associate
  end subroutine AerobicHeterO2Uptake

!------------------------------------------------------------------------------------------
  subroutine StageFuncGuild(N,NGL,K,TotActMicrobiom,FOQC,FOQA,micfor,naqfdiag,nmicdiag,nmics)
  implicit none
  integer, intent(in) :: N   !functional group id
  integer, intent(in) :: NGL !functional guild id
  integer, intent(in) :: K   !complex id
  real(r8),intent(in) :: TotActMicrobiom         !total active microbial biomass
  type(micforctype), intent(in) :: micfor
  real(r8),intent(out) :: FOQC                   !fraction of DOC acetate demand over all microbial demand, soil/band
  real(r8),intent(out) :: FOQA                   !fraction of acetate demand over all microbial demand, soil/band
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_Diag_type),intent(inout) :: nmicdiag
  type(Microbe_State_type), intent(inout):: nmics
  real(r8) :: WatStressMicb
  associate(                                                 &
    PSISoilMatricP         => micfor%PSISoilMatricP,         & !Soil matric water potential controlling microbial water stress
    ZEROS                  => micfor%ZEROS,                  & !Small mass or flux threshold used by the routine
    RDOMEcoDmndPrev        => micfor%RDOMEcoDmndPrev,        & !Previous-hour ecosystem DOC demand in each complex; competition denominator
    RAcetateEcoDmndPrev    => micfor%RAcetateEcoDmndPrev,    & !Previous-hour ecosystem acetate demand in each complex; competition denominator
    RAcetateUptkHeterPrev  => micfor%RAcetateUptkHeterPrev,  & !Previous-hour guild acetate uptake/demand used for competition
    mid_Aerob_Fungi        => micpar%mid_Aerob_Fungi,        & !Functional-group identifier for aerobic fungi
    mid_Facult_DenitBacter => micpar%mid_Facult_DenitBacter, & !Functional-group identifier for facultative denitrifying bacteria
    GrowthEnvScalHeter     => nmics%GrowthEnvScalHeter,      & !Temperature and water-potential multiplier on heterotrophic growth [-]
    FracHeterBiomOfActK    => nmics%FracHeterBiomOfActK,     & !Guild fraction of active heterotrophic biomass in complex K [-]
    RDOCUptkHeterPrev      => micfor%RDOCUptkHeterPrev,      & !Previous-hour guild DOC uptake/demand used for competition
    FracOMActHeter         => nmics%FracOMActHeter,          & !Guild/complex fraction of total active microbial C in the layer [-]
    OMActHeter             => nmics%OMActHeter,              & !Active microbial C biomass by heterotrophic guild and complex K
    TempMaintRHeter        => nmics%TempMaintRHeter,         & !Guild temperature multiplier on heterotrophic maintenance [-]
    TOMEK                  => nmicdiag%TOMEK          ,      & !Total active heterotrophic C/N/P in each substrate complex K
    FracNO2ReduxHeter      => nmics%FracNO2ReduxHeter,       & !Guild share of active nitrite-consuming biomass for competition [-]
    WSensGroHeter          => nmics%WSensGroHeter         ,  & !Guild soil-water-potential multiplier on heterotrophic growth [-]
    TSensGroHeter          => nmics%TSensGroHeter         ,  & !Guild temperature multiplier on heterotrophic growth [-]
    TSensMaintR            => nmicdiag%TSensMaintR,          & !Layer temperature response for microbial maintenance [-]
    TotBiomNO2Consumers    => nmicdiag%TotBiomNO2Consumers,  & !Active C summed over staged nitrite-competing microbial groups
    TSensGrowth            => nmicdiag%TSensGrowth           & !Layer temperature response for microbial growth [-]
  )

  ! WatStressMicb=water potential (PSISoilMatricP_vr) effect on microbial respiration
  ! OXKX=Km for O2 uptake
  ! OXKM=Km for heterotrophic O2 uptake set in starts.f
  !  GrowthEnvScalHeter=combined temp and water stress effect on growth respiration
  !  TempMaintRHeter=temperature effect on maintenance respiration
  !
  !different guilds can have different temperature and moisture sensitivity
  IF(N.EQ.mid_Aerob_Fungi)THEN
    WatStressMicb=EXP(0.1_r8*AMAX1(PSISoilMatricP,-500._r8))
  ELSE
    WatStressMicb=EXP(0.2_r8*AMAX1(PSISoilMatricP,-500._r8))
  ENDIF

  WSensGroHeter(NGL,K)=WatStressMicb
  TSensGroHeter(NGL,K)=TSensGrowth

  GrowthEnvScalHeter(NGL,K) = WSensGroHeter(NGL,K)*TSensGroHeter(NGL,K)
  TempMaintRHeter(NGL,K)    = TSensMaintR

  ! FracOMActHeter,FracNO2ReduxHeter=fraction of total active biomass C,N in each N and K

  IF(TotActMicrobiom.GT.ZEROS)THEN
    FracOMActHeter(NGL,K)=OMActHeter(NGL,K)/TotActMicrobiom
  ELSE
    FracOMActHeter(NGL,K)=1.0_r8
  ENDIF

  IF(TotBiomNO2Consumers.GT.ZEROS .and. N.EQ.mid_Facult_DenitBacter)THEN
    FracNO2ReduxHeter(NGL,K)=OMActHeter(NGL,K)/TotBiomNO2Consumers
  ELSE
    FracNO2ReduxHeter(NGL,K)=1.0_r8
  ENDIF

  IF(TOMEK(ielmc,K).GT.ZEROS)THEN
    FracHeterBiomOfActK(NGL,K)=OMActHeter(NGL,K)/TOMEK(ielmc,K)
  ELSE
    FracHeterBiomOfActK(NGL,K)=1.0_r8
  ENDIF
  !
  IF(RDOMEcoDmndPrev(K).GT.ZEROS)THEN
    FOQC=AMAX1(FMN,RDOCUptkHeterPrev(NGL,K)/RDOMEcoDmndPrev(K))
  ELSE
    FOQC=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
  ENDIF

  naqfdiag%TFOQC=naqfdiag%TFOQC+FOQC

  IF(RAcetateEcoDmndPrev(K).GT.ZEROS)THEN
    FOQA=AMAX1(FMN,RAcetateUptkHeterPrev(NGL,K)/RAcetateEcoDmndPrev(K))
  ELSE
    FOQA=AMAX1(FMN,FracHeterBiomOfActK(NGL,K))
  ENDIF
  naqfdiag%TFOQA  = naqfdiag%TFOQA+FOQA

  end associate
  end subroutine StageFuncGuild

!------------------------------------------------------------------------------------------

  subroutine StageAutotroph(NGL,N,TOMEAutoKC,micfor,nmics,nmicdiag)

  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(in) :: TOMEAutoKC
  type(micforctype), intent(in) :: micfor
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  real(r8) :: WatStressMicb

  associate(                                             &
    PSISoilMatricP      => micfor%PSISoilMatricP,        & !Soil matric water potential controlling microbial water stress
    OMActAutor          => nmics%OMActAutor,             & !Active microbial C biomass by autotrophic guild
    FracOMActAutor      => nmics%FracOMActAutor,         & !Guild fraction of total active microbial C in the layer [-]
    FracNO2XupAutor     => nmics%FracNO2XupAutor,        & !Staged autotrophic nitrite-competition weight [-]
    FracAutorBiomOfActK => nmics%FracAutorBiomOfActK,    & !Guild fraction of active biomass in the autotrophic assemblage [-]
    ZEROS               => micfor%ZEROS,                 & !Small mass or flux threshold used by the routine
    TSensMaintRAutor    => nmics%TSensMaintRAutor,       & !Guild temperature multiplier on autotrophic maintenance [-]
    GrowthEnvScalAutor  => nmics%GrowthEnvScalAutor,     & !Temperature and water-potential multiplier on autotrophic growth [-]
    TotBiomNO2Consumers => nmicdiag%TotBiomNO2Consumers, & !Active C summed over staged nitrite-competing microbial groups
    TotActMicrobiom     => nmicdiag%TotActMicrobiom,     & !Layer total active microbial C across heterotrophs and autotrophs
    WSensGroAutor       => nmics%WSensGroAutor,          & !Guild soil-water-potential multiplier on autotrophic growth [-]
    TSensGroAutor       => nmics%TSensGroAutor,          & !Guild temperature multiplier on autotrophic growth [-]
    TSensGrowth         => nmicdiag%TSensGrowth,         & !Layer temperature response for microbial growth [-]
    TSensMaintR         => nmicdiag%TSensMaintR          & !Layer temperature response for microbial maintenance [-]
  )
  !replace with trait specific parameterization
  WatStressMicb           = EXP(0.2_r8*PSISoilMatricP)
  WSensGroAutor(NGL)      = WatStressMicb
  TSensGroAutor(NGL)      = TSensGrowth

  GrowthEnvScalAutor(NGL) = AZMAX1(TSensGroAutor(NGL)*WSensGroAutor(NGL))
  TSensMaintRAutor(NGL)   = AZMAX1(TSensMaintR)

! FracOMActHeter,FOMN=fraction of total active biomass C,N in each N and K

  IF(TotActMicrobiom.GT.ZEROS)THEN
    FracOMActAutor(NGL)=OMActAutor(NGL)/TotActMicrobiom
  ELSE
    FracOMActAutor(NGL)=1.0_r8
  ENDIF

  !Use the same consumer set as StageBGCEnvironCondition's NO2 biomass total.
  IF(TotBiomNO2Consumers.GT.ZEROS .AND. (N.EQ.micpar%mid_AutoAmmoniaOxidBacter &
    .OR. N.EQ.micpar%mid_AutoNitriteOxidBacter .OR. N.EQ.micpar%mid_AutoAMONC10))THEN
    FracNO2XupAutor(NGL)=OMActAutor(NGL)/TotBiomNO2Consumers
  ELSE
    FracNO2XupAutor(NGL)=1.0_r8
  ENDIF

  IF(TOMEAutoKC.GT.ZEROS)THEN
    FracAutorBiomOfActK(NGL)=OMActAutor(NGL)/TOMEAutoKC
  ELSE
    FracAutorBiomOfActK(NGL)=1.0_r8
  ENDIF
  end associate
  end subroutine StageAutotroph
!------------------------------------------------------------------------------------------
  subroutine CalcRespMaintAutor(I,J,NGL,RMOMK,micfor,micstt,micflx,nmicf,nmics)
  implicit none
  integer, intent(in) :: I,J,NGL
  real(r8), intent(in) :: RMOMK(2)  !effect of low microbial C concentration on maintenance respiration
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  character(len=*),parameter :: subname='CalcRespMaintAutor'
  REAL(R8) :: FPH,RMOMX
  integer :: MID1

  associate(                                             &
    OMActAutor           => nmics%OMActAutor,            & !Active microbial C biomass by autotrophic guild; not referenced here
    OMN2Autor            => nmics%OMN2Autor,             & !Active structural autotrophic N by guild
    TSensMaintRAutor     => nmics%TSensMaintRAutor,      & !Guild temperature multiplier on autotrophic maintenance [-]
    RMaintDmndAutor      => nmicf%RMaintDmndAutor,       & !Maintenance-C demand by live biomass compartment and autotrophic guild
    RMaintRespAutor      => micflx%RMaintRespAutor,      & !Total hourly autotrophic guild maintenance-C demand
    pH                   => micfor%pH,                   & !Layer pH used in microbial rate and nutrient-speciation responses
    mBiomeAutor          => micstt%mBiomeAutor,          & !C/N/P pools indexed by element and flattened guild/biomass compartment
    JGniA                => micpar%JGniA,                & !First guild index for each autotrophic functional group; not referenced here
    JGnfA                => micpar%JGnfA                 & !Last guild index for each autotrophic functional group; not referenced here
  )
  call PrintInfo('beg '//subname)

  MID1                   = micpar%get_micb_id(iLbiom_kinetic,NGL)
  FPH                    = 1.0_r8+AZMAX1(0.25_r8*(6.5_r8-PH))
  RMOMX                  = RMOM*TSensMaintRAutor(NGL)*FPH
  RMaintDmndAutor(iLbiom_kinetic,NGL) = mBiomeAutor(ielmn,MID1)*RMOMX*RMOMK(iLbiom_kinetic)
  RMaintDmndAutor(iLbiom_struct,NGL)  = OMN2Autor(NGL)*RMOMX*RMOMK(iLbiom_struct)
  !
  !     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
  !
  !     RMaintRespAutor=total maintenance respiration, as measured by C needs to be respired
  !     RGrowthRespAutor=growth respiration
  !     RMaintDefcitcitAutor=senescence respiration
  !
  RMaintRespAutor(NGL)      = RMaintDmndAutor(iLbiom_kinetic,NGL)+RMaintDmndAutor(iLbiom_struct,NGL)
  call PrintInfo('end '//subname)
  end associate
  end subroutine CalcRespMaintAutor
!------------------------------------------------------------------------------------------

  subroutine CalcRespMaintHeter(NGL,K,RMOMK,micfor,micstt,micflx,nmicf,nmics)
  implicit none
  integer, intent(in) :: NGL,K
  real(r8), intent(in) :: RMOMK(2)  !effect of low microbial C concentration on maintenance respiration
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  character(len=*),parameter :: subname='CalcRespMaintHeter'
  REAL(R8) :: FPH,RMOMX
  integer :: MID1

  associate(                                             &
    OMN2                 => nmics%OMN2,                  & !Active structural heterotrophic N by guild and complex K
    TempMaintRHeter      => nmics%TempMaintRHeter,       & !Guild temperature multiplier on heterotrophic maintenance [-]
    RMaintDmndHeter      => nmicf%RMaintDmndHeter,       & !Maintenance-C demand by live compartment, heterotrophic guild and complex
    RMaintRespHeter      => nmicf%RMaintRespHeter,       & !Total hourly heterotrophic maintenance-C demand by guild and complex
    pH                   => micfor%pH,                   & !Layer pH used in microbial rate and nutrient-speciation responses
    mBiomeHeter          => micstt%mBiomeHeter           & !C/N/P pools indexed by element, flattened guild/compartment and complex K
  )
  call PrintInfo('beg '//subname)

  !     RMOMK=effect of low microbial C concentration on mntc respn

  FPH                                  = 1.0_r8+AZMAX1(0.25_r8*(6.5_r8-PH))
  RMOMX                                = RMOM*TempMaintRHeter(NGL,K)*FPH
  MID1                                 = micpar%get_micb_id(iLbiom_kinetic,NGL)
  RMaintDmndHeter(iLbiom_kinetic,NGL,K) = mBiomeHeter(ielmn,MID1,K)*RMOMX*RMOMK(iLbiom_kinetic)
  RMaintDmndHeter(iLbiom_struct,NGL,K)  = OMN2(NGL,K)*RMOMX*RMOMK(iLbiom_struct)
  !
  !     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
  !
  !     RMaintRespHeter=total maintenance respiration
  !     RGrowthRespHeter=growth respiration
  !     RMaintDefcitcitHeter=senescence respiration
  !
  RMaintRespHeter(NGL,K)      = RMaintDmndHeter(iLbiom_kinetic,NGL,K)+RMaintDmndHeter(iLbiom_struct,NGL,K)
  call PrintInfo('end '//subname)
  end associate
  end subroutine CalcRespMaintHeter

end module MicrobMathFuncMod
