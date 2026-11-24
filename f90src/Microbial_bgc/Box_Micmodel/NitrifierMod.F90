module NitrifierMod
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: PrintInfo
  use minimathmod,          only: AZMAX1
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod,    only: CalcRespMaintAutor, StageAutotroph

  implicit none

  private
  public :: AmmoniaOxiDenitCatabolism
  public :: AmmoniaOxidizerCatabolism
  public :: NitriteOxidizerCatabolism

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains
!------------------------------------------------------------------------------------------

  subroutine AmmoniaOxiDenitCatabolism(I,J,N,VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics, micflx)
  !
  !2NO2(-) + NH3 -> 1.5N2O + 2OH(-) + 0.5H2O, molar based
  !the energy used is used to assimilate CO2 for biomass
  !nitrate-ammonifying bacteria

  implicit none
  integer, intent(in) :: I,J,N
  real(r8), intent(in) :: VOLWZ
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx
  real(r8) :: FNO2S,FNO2B
  REAL(R8) :: FNO2,FNB2
  real(r8) :: FVMXDX
  real(r8) :: ROXYD,RDNOT
  real(r8) :: VMXDXS
  real(r8) :: VMXDXB
  real(r8) :: VMXDXT
  real(r8) :: VMXD4
  real(r8) :: VMXD4S
  real(r8) :: VMXD4B
  real(r8) :: ZNO2SX,ZNO2BX,XCO2
  real(r8) :: FNH4,FNB4,NH4AvailSoil,NH4AvailBand
  real(r8) :: NH4RemainingSoil,NH4RemainingBand
  integer  :: NGL
!     begin_execution
  associate(                                                &
    FracNO2XupAutor        => nmics%FracNO2XupAutor,        & !Staged autotrophic nitrite-competition weight [-]
    FracOMActAutor        => nmics%FracOMActAutor,        & !Guild fraction of active microbial C; fallback NH4 competition weight [-]
    RNH4EcoDmndSoilPrev    => micfor%RNH4EcoDmndSoilPrev,  & !Previous-hour ecosystem nonband NH4-N demand
    RNH4EcoDmndBandPrev    => micfor%RNH4EcoDmndBandPrev,  & !Previous-hour ecosystem fertilizer-band NH4-N demand
    VLNH4                 => micfor%VLNH4,               & !Nonband NH4 fraction used in fallback competition [-]
    VLNHB                 => micfor%VLNHB,               & !Fertilizer-band NH4 fraction used in fallback competition [-]
    ZNH4S                 => micstt%ZNH4S,               & !Nonband NH4-N donor pool before microbial uptake
    ZNH4B                 => micstt%ZNH4B,               & !Fertilizer-band NH4-N donor pool before microbial uptake
    RNH3OxidAutorPrev      => micflx%RNH3OxidAutorPrev,    & !Previous-hour nonband NH4-N oxidation demand used for competition
    RNH3OxidAutorBandPrev  => micflx%RNH3OxidAutorBandPrev,& !Previous-hour fertilizer-band NH4-N oxidation demand used for competition
    RNH3OxidAutor          => micflx%RNH3OxidAutor,        & !Potential nonband NH4-N demand from aerobic oxidation plus denitrification
    RNH3OxidAutorBand      => micflx%RNH3OxidAutorBand,    & !Potential fertilizer-band NH4-N demand from aerobic oxidation plus denitrification
    RO2Dmnd4GrossRespAutor => nmicf%RO2Dmnd4GrossRespAutor, & !Potential O2 demand supporting autotrophic gross respiration
    OMActAutor             => nmics%OMActAutor,             & !Active microbial C biomass by autotrophic guild
    RO2Uptk4RespAutor      => nmicf%RO2Uptk4RespAutor,      & !Realized O2 uptake attributed to autotrophic gross respiration
    RNO3UptkAutor          => nmicf%RNO3UptkAutor,          & !Autotrophic nitrate-N reduction bookkeeping; zero for nitrifiers in current code
    RNOxReduxAutorSoil     => nmicf%RNOxReduxAutorSoil,     & !Nonband NO2-N reduction by NH3 oxidation; contributes to N2O production
    RNOxReduxAutorBand     => nmicf%RNOxReduxAutorBand,     & !Fertilizer-band NO2-N reduction by NH3 oxidation; contributes to N2O production
    RNOxReduxRespAutorLim  => nmicf%RNOxReduxRespAutorLim,  & !C-equivalent respiration supported by autotrophic nitrite reduction
    RSMetaOxidSoilAutor    => nmicf%RSMetaOxidSoilAutor,    & !Nonband NH3-N oxidation including the added denitrification contribution
    RSMetaOxidBandAutor    => nmicf%RSMetaOxidBandAutor,    & !Fertilizer-band NH3-N oxidation including the added denitrification contribution
    RTotNH3OxidSoilAutor   => nmicf%RTotNH3OxidSoilAutor,   & !Ammonia-N oxidation to nitrite summed over nonband nitrifier guilds
    RTotNH3OxidBandAutor   => nmicf%RTotNH3OxidBandAutor,   & !Ammonia-N oxidation to nitrite summed over fertilizer-band nitrifier guilds
    RNO2EcoUptkSoilPrev    => micfor%RNO2EcoUptkSoilPrev,   & !Previous-hour ecosystem nonband NO2-N uptake; competition denominator
    VLNO3                  => micfor%VLNO3,                 & !Nonband fraction for nitrate/nitrite pools and uptake capacity [-]
    VLNOB                  => micfor%VLNOB,                 & !Fertilizer-band fraction for nitrate/nitrite pools and uptake capacity [-]
    RNO2EcoUptkBandPrev    => micfor%RNO2EcoUptkBandPrev,   & !Previous-hour ecosystem fertilizer-band NO2-N uptake; competition denominator
    ZEROS                  => micfor%ZEROS,                 & !Small mass or flux threshold used by the routine
    ZEROS2                 => micfor%ZEROS2,                & !Small water-volume threshold used by the routine
    CCO2S                  => micstt%CCO2S,                 & !Dissolved CO2-C concentration for substrate saturation
    CNO2B                  => micstt%CNO2B,                 & !Dissolved NO2-N concentration in fertilizer-band soil
    CNO2S                  => micstt%CNO2S,                 & !Dissolved NO2-N concentration in nonband soil
    ZNO2B                  => micstt%ZNO2B,                 & !NO2-N pool in fertilizer-band soil
    ZNO2S                  => micstt%ZNO2S,                 & !NO2-N pool in nonband soil
    JGniA                  => micpar%JGniA,                 & !First guild index for each autotrophic functional group
    JGnfA                  => micpar%JGnfA,                 & !Last guild index for each autotrophic functional group
    RNO2XupAutorPrev       => micflx%RNO2XupAutorPrev,      & !Previous-hour autotrophic nonband NO2-N redox uptake; reaction depends on functional group; used for competition
    RNO2XupAutorBandPrev   => micflx%RNO2XupAutorBandPrev,  & !Previous-hour autotrophic fertilizer-band NO2-N redox uptake; reaction depends on functional group; used for competition
    RNO2XupAutor           => micflx%RNO2XupAutor,          & !Potential nonband NO2-N demand for nitrifier denitrification
    RNO2XupAutorBand       => micflx%RNO2XupAutorBand       & !Potential fertilizer-band NO2-N demand for nitrifier denitrification
  )
  !
  !     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
  !     POPULATIONS
  !
  !     FNO2,FNB2=fraction of total biological demand for NO2
  !
  !     CCO2S=aqueous CO2 concentration
  !
  RTotNH3OxidSoilAutor = SUM(RSMetaOxidSoilAutor(JGniA(N):JGnfA(N)))
  RTotNH3OxidBandAutor = SUM(RSMetaOxidBandAutor(JGniA(N):JGnfA(N)))
  !Aerobic oxidation has already spent part of each compartment's donor pool.
  NH4RemainingSoil = MAX(0._r8,ZNH4S-RTotNH3OxidSoilAutor)
  NH4RemainingBand = MAX(0._r8,ZNH4B-RTotNH3OxidBandAutor)
  XCO2                 = CCO2S/(CCO2S+CCKM)
  DO NGL=JGniA(N),JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8 .or. RO2Dmnd4GrossRespAutor(NGL).LE.0.0_r8)cycle

    IF(RNO2EcoUptkSoilPrev.GT.ZEROS)THEN
      FNO2=AMAX1(FMN,RNO2XupAutorPrev(NGL)/RNO2EcoUptkSoilPrev)
    ELSE
      FNO2=AMAX1(FMN,FracNO2XupAutor(NGL)*VLNO3)
    ENDIF
    IF(RNO2EcoUptkBandPrev.GT.ZEROS)THEN
      FNB2=AMAX1(FMN,RNO2XupAutorBandPrev(NGL)/RNO2EcoUptkBandPrev)
    ELSE
      FNB2=AMAX1(FMN,FracNO2XupAutor(NGL)*VLNOB)
    ENDIF
    naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
    naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2

    !Reuse the NH4 allocation used by aerobic ammonia oxidation; denitrification
    !can spend only the remainder of that allocation, in the same compartment.
    IF(RNH4EcoDmndSoilPrev.GT.ZEROS)THEN
      FNH4=AMAX1(FMN,RNH3OxidAutorPrev(NGL)/RNH4EcoDmndSoilPrev)
    ELSE
      FNH4=AMAX1(FMN,VLNH4*FracOMActAutor(NGL))
    ENDIF
    IF(RNH4EcoDmndBandPrev.GT.ZEROS)THEN
      FNB4=AMAX1(FMN,RNH3OxidAutorBandPrev(NGL)/RNH4EcoDmndBandPrev)
    ELSE
      FNB4=AMAX1(FMN,VLNHB*FracOMActAutor(NGL))
    ENDIF
    NH4AvailSoil=MAX(0._r8,MIN(NH4RemainingSoil,FNH4*ZNH4S-RSMetaOxidSoilAutor(NGL)))
    NH4AvailBand=MAX(0._r8,MIN(NH4RemainingBand,FNB4*ZNH4B-RSMetaOxidBandAutor(NGL)))
    !
    !     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
    !     ACTIVE NITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2 AND CO2
    !     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
    !     NOT ACCEPTED BY O2
    !
    !     ROXYD=O2 demand RO2Dmnd4RespHeter not met by O2 uptake RO2Uptk4RespHeter
    !     VMXD4=demand for NO2-N reduction
    !     VMXDXS,VMXDXB=maximum NO2 reduction in non-band, band
    !     FNO2S,FNO2B=fractions of total NO2 in non-band, band
    !     CNO2S,CNO2B=NO2 concentrations in non-band, band
    !     Z2KM=Km for NO2 uptake
    !     FVMXDX=nonlinear effect of product inhibition for NOx reduction
    !     VMKI=product inhibition for NOx reduction
    !     VMXD4S,VMXD4B=substrate-unlimited NO2 reduction in non-band,band
    !     RNO2ReduxHeterSoil,RNO2ReduxHeterBand=substrate-limited NO2 reduction in non-band,band
    !     RNOxReduxRespDenitUlm,RNOxReduxRespDenitLim=total substrate-unltd,-ltd respn from NO2 reduction
    !     ECNO=efficiency CO2 conversion to biomass
    !     ECHZ=growth respiration efficiency
    !     RSMetaOxidSoilAutor,RSMetaOxidBandAutor=total O2-limited (1)NH4,(2)NO2,(3)CH4 oxidation
    !one O2 accepts 4e, one NO2(-) accepts 2e

    FNO2S  = VLNO3
    FNO2B  = VLNOB
    ROXYD  = AZMAX1(RO2Dmnd4GrossRespAutor(NGL)-RO2Uptk4RespAutor(NGL))

    !why is 0.875? 2NO2(-) + NH3 -> 1.5N2O+ 2OH(-)+0.5H2O, 0.875=14/16, N(+3)->N(+), O->O(2-),
    VMXD4  = 0.875_r8*ROXYD*XCO2
    VMXDXS = FNO2S*VMXD4*CNO2S/(CNO2S+Z2KM)
    VMXDXB = FNO2B*VMXD4*CNO2B/(CNO2B+Z2KM)
    VMXDXT = VMXDXS+VMXDXB

    IF(VOLWZ.GT.ZEROS2)THEN
      FVMXDX=1.0_r8/(1.0_r8+VMXDXT/(VMKI*VOLWZ))
    ELSE
      FVMXDX=0.0_r8
    ENDIF
    VMXD4S = VMXDXS*FVMXDX
    VMXD4B = VMXDXB*FVMXDX

    !Guild share of existing NO2 plus current-hour aerobic NH3 oxidation.
    !The production totals are fixed before the loop to avoid counting
    !NH3 consumed by denitrification below as new NO2 production.
    ZNO2SX                  = AZMAX1(ZNO2S+RTotNH3OxidSoilAutor)*FNO2
    ZNO2BX                  = AZMAX1(ZNO2B+RTotNH3OxidBandAutor)*FNB2
    !Two units of NO2-N require one unit of NH4-N. Limit the reaction before
    !deriving respiration, biomass supply, and N2O production from its rate.
    RNOxReduxAutorSoil(NGL) = AZMAX1(AMIN1(VMXD4S,ZNO2SX,2._r8*NH4AvailSoil))
    RNOxReduxAutorBand(NGL) = AZMAX1(AMIN1(VMXD4B,ZNO2BX,2._r8*NH4AvailBand))
    NH4RemainingSoil=MAX(0._r8,NH4RemainingSoil-0.5_r8*RNOxReduxAutorSoil(NGL))
    NH4RemainingBand=MAX(0._r8,NH4RemainingBand-0.5_r8*RNOxReduxAutorBand(NGL))

    !total NO2 reduced
    RDNOT                      = RNOxReduxAutorSoil(NGL)+RNOxReduxAutorBand(NGL)

    !C-biomass yield from the catabolic energy
    !ENOX: respiraiton coefficient = 1/(1+G_c/G_a), where G_c is Gibbs free energy of the catabolic reaction
    !and G_a is the Gibbs free energy of anabolic reaction, 1.5CO2+NH3(aq)+OH(-)->1.5CH2O+NO2(-)+0.5H2O.
    !ECNO: efficiency of CO2 conversion into biomass through nitrate reduction
    !ECNO refers to the fraction of NH3 used for making CH2O.
    !NO2(-) reduction is used to synthesize CH2O, and
    !2NO2(-) + CH2O -> N2O(aq) + CO2 + 2OH(-)
    RNOxReduxRespAutorLim(NGL) = RDNOT*ECNO*ENOX
    RNO3UptkAutor(NGL)         = 0.0_r8              !currently no NO3 reduction by nitrifiers
    RNO2XupAutor(NGL)          = VMXD4S
    RNO2XupAutorBand(NGL)      = VMXD4B
    !Include both pathways in next hour's NH4 competition demand. Keep these
    !potential demands separate from the donor-limited realized oxidation above.
    RNH3OxidAutor(NGL)         = RNH3OxidAutor(NGL)+0.5_r8*VMXD4S
    RNH3OxidAutorBand(NGL)     = RNH3OxidAutorBand(NGL)+0.5_r8*VMXD4B

    !NH4 oxidation by NO2(-), NH3+2NO2(-) -> 1.5N2O+2OH(-)+0.5H2O
    !NH4 -> N2O, 2NO2-> N2O
    RSMetaOxidSoilAutor(NGL)=RSMetaOxidSoilAutor(NGL)+RNOxReduxAutorSoil(NGL)/2._r8
    RSMetaOxidBandAutor(NGL)=RSMetaOxidBandAutor(NGL)+RNOxReduxAutorBand(NGL)/2._r8

  ENDDO
  end associate
  end subroutine AmmoniaOxiDenitCatabolism

!------------------------------------------------------------------------------------------
  subroutine AmmoniaOxidizerCatabolism(I,J,N,RMOMK,TOMEAutoKC,VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)
  !
  !Description:
  ! autotrophic NH3 oxidizer
  !NH3 + 1.5O2 -> NO2(-) + H+ + H2O
  !it first converts CO2 into CH2O (anabolic reaction) via
  !NH3 + 1.5CO2+ OH(-) -> 1.5CH2O + NO2(-) + 0.5H2O
  !CH2O is then respired to produce gross respiration
  implicit none
  integer,  intent(in)  :: I,J,N
  real(r8), intent(in) :: RMOMK(2)
  REAL(R8), INTENT(IN)  :: VOLWZ
  real(r8), intent(in)  :: TOMEAutoKC
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type),INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='AmmoniaOxidizerCatabolism'
  real(r8)  :: RGOMP         !O2-unlimited/potential respiration [gC h-1], as a measure of gross respiraiton/energy for maintenance+growth
  real(r8)  :: RVOXP         !potential NH3 oxidation, [gN h-1 d-2]
  real(r8)  :: RVOXPA        !potential oxidation of non-band soil NH3 [gN h-1 d-2]
  real(r8)  :: RVOXPB        !potential oxidation of band soil NH3 [gN h-1 d-2]
  real(r8) :: FNH4S,FNHBS
  real(r8) :: FNH4,FNB4
  real(r8) :: FCN4S,FCN4B
  real(r8) :: RNNH4,RNNHB
  real(r8) :: VMXX,VMX4S
  real(r8) :: VMX4B
  real(r8) :: ZNFN4S,ZNFN4B
  real(r8) :: VMAX,XCO2
  integer  :: NGL

!     begin_execution
  associate(                                                &
    GrowthEnvScalAutor     => nmics%GrowthEnvScalAutor,     & !Temperature and water-potential multiplier on autotrophic growth [-]
    FBiomNutStoiScalAutor  => nmics%FBiomNutStoiScalAutor,  & !Combined N/P stoichiometric multiplier on guild metabolic capacity [-]
    FSBSTAutor             => nmicdiag%FSBSTAutor,          & !Guild primary-substrate saturation diagnostic [-]
    FracOMActAutor         => nmics%FracOMActAutor,         & !Guild fraction of total active microbial C in the layer [-]
    OMActAutor             => nmics%OMActAutor,             & !Active microbial C biomass by autotrophic guild
    TSensGroAutor          => nmics%TSensGroAutor,          & !Guild temperature multiplier on autotrophic growth [-]
    RO2Dmnd4GrossRespAutor => nmicf%RO2Dmnd4GrossRespAutor, & !Potential O2 demand supporting autotrophic gross respiration
    ECHZAutor              => nmicf%ECHZAutor,              & !Guild respiration fraction used to convert growth respiration to C uptake [-]
    RSMetaOxidSoilAutor    => nmicf%RSMetaOxidSoilAutor,    & !Nonband potential aerobic NH3-N oxidation; O2-limited later
    RSMetaOxidBandAutor    => nmicf%RSMetaOxidBandAutor,    & !Fertilizer-band potential aerobic NH3-N oxidation; O2-limited later
    RespGrossAutor         => nmicf%RespGrossAutor,         & !Gross respiration C equivalent by autotrophic guild
    VLNH4                  => micfor%VLNH4,                 & !Nonband fraction for ammonium/ammonia pools and uptake capacity [-]
    VLNHB                  => micfor%VLNHB,                 & !Fertilizer-band fraction for ammonium/ammonia pools and uptake capacity [-]
    ZEROS                  => micfor%ZEROS,                 & !Small mass or flux threshold used by the routine
    ZEROS2                 => micfor%ZEROS2,                & !Small water-volume threshold used by the routine
    RNH4EcoDmndSoilPrev    => micfor%RNH4EcoDmndSoilPrev,   & !Previous-hour ecosystem NH4-N demand in nonband soil; competition denominator
    RNH4EcoDmndBandPrev    => micfor%RNH4EcoDmndBandPrev,   & !Previous-hour ecosystem NH4-N demand in fertilizer-band soil; competition denominator
    ZNFN0                  => micstt%ZNFN0,                 & !Nitrification-inhibition reference established with fertilizer application
    CCO2S                  => micstt%CCO2S,                 & !Dissolved CO2-C concentration for substrate saturation
    ZNFNI                  => micstt%ZNFNI,                 & !Temperature-decaying nitrification-inhibition state used with ZNFN0
    CNH3S                  => micstt%CNH3S,                 & !Dissolved NH3-N concentration in nonband soil
    CNH3B                  => micstt%CNH3B,                 & !Dissolved NH3-N concentration in fertilizer-band soil
    CNH4S                  => micstt%CNH4S,                 & !Dissolved NH4-N concentration in nonband soil
    CNH4B                  => micstt%CNH4B,                 & !Dissolved NH4-N concentration in fertilizer-band soil
    ZNH4S                  => micstt%ZNH4S,                 & !NH4-N pool in nonband soil
    ZNH4B                  => micstt%ZNH4B,                 & !NH4-N pool in fertilizer-band soil
    JGniA                  => micpar%JGniA,                 & !First guild index for each autotrophic functional group
    JGnfA                  => micpar%JGnfA,                 & !Last guild index for each autotrophic functional group
    RNH3OxidAutorPrev      => micflx%RNH3OxidAutorPrev,     & !Previous-hour nonband ammonia-N oxidation by nitrifier guilds; used for competition
    RNH3OxidAutorBandPrev  => micflx%RNH3OxidAutorBandPrev, & !Previous-hour fertilizer-band ammonia-N oxidation by nitrifier guilds; used for competition
    RNH3OxidAutor          => micflx%RNH3OxidAutor,         & !Nonband ammonia-N oxidation by nitrifier guilds
    RNH3OxidAutorBand      => micflx%RNH3OxidAutorBand,     & !Fertilizer-band ammonia-N oxidation by nitrifier guilds
    RO2MetaDmndAutor       => micflx%RO2MetaDmndAutor       & !Total autotrophic O2 demand from respiration and substrate oxidation
  )
!
!     FACTOR TO REGULATE COMPETITION FOR NH4 AMONG DIFFERENT
!     MICROBIAL AND ROOT POPULATIONS FNH4
!
!     FNH4,FNB4=frac of total biol demand for NH4 in non-band, band
!
!     CCO2S=aqueous CO2 concentration
!
  call PrintInfo('beg '//subname)
  XCO2    = CCO2S/(CCO2S+CCKM)
  DO NGL  = JGniA(N), JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle
    call StageAutotroph(NGL,N,TOMEAutoKC,micfor,nmics,nmicdiag)

    call CalcRespMaintAutor(I,J,NGL,RMOMK,micfor,micstt,micflx,nmicf,nmics)

    FNH4S=VLNH4
    FNHBS=VLNHB
    IF(RNH4EcoDmndSoilPrev.GT.ZEROS)THEN
      FNH4=AMAX1(FMN,RNH3OxidAutorPrev(NGL)/RNH4EcoDmndSoilPrev)
    ELSE
      FNH4=AMAX1(FMN,VLNH4*FracOMActAutor(NGL))
    ENDIF
    IF(RNH4EcoDmndBandPrev.GT.ZEROS)THEN
      FNB4=AMAX1(FMN,RNH3OxidAutorBandPrev(NGL)/RNH4EcoDmndBandPrev)
    ELSE
      FNB4=AMAX1(FMN,VLNHB*FracOMActAutor(NGL))
    ENDIF
!    naqfdiag%TFNH4X = naqfdiag%TFNH4X+FNH4
!    naqfdiag%TFNH4B = naqfdiag%TFNH4B+FNB4
    !
    !     NITRIFICATION INHIBITION
    !
    !     ZNFN0=inhibition when fertilizer added
    !     ZNFNI=reduction in inhibition since fertilizer added
    !     CNH4S,CNH4B=NH4 concentrations in non-band, band
    !     TSensGrowth=temperature effect
    !     RNFNI=rate constant for inhibition decline
    !     ZHKI=inhibition from high CNH4
    !     ZNFN4S,ZNFN4B=inhibition in non-band, band
    !
    IF(ZNFN0.GT.ZEROS)THEN
      ZNFNI=ZNFNI*(1.0_r8-RNFNI*TSensGroAutor(NGL))
      ZNFN4S=ZNFN0-ZNFNI/(1.0_r8+CNH4S/ZHKI)
      ZNFN4B=ZNFN0-ZNFNI/(1.0_r8+CNH4B/ZHKI)
    ELSE
      ZNFN4S=1.0_r8
      ZNFN4B=1.0_r8
    ENDIF
    !
    !     NH3 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
    !     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
    !     NH3 CONCENTRATIONS IN BAND AND NON-BAND SOIL ZONES
    !
    !     ECHZ=growth respiration efficiency
    !     VMXX=potential NH3 oxidation, VMXNH3Oxi=specific oxidation
    !     TFNG=temperature+water limitation, FBiomNutStoiScalAutorr=N,P limitation
    !     XCO2=aqueous CO2 limitation, OMA=active biomass
    !     VMAX= non-substrate limited NH3 oxidation
    !     VHKI=nonlinear increase in VMAX with VMXNH3Oxi
    !     FNH4S,FNHBS=fractions of NH4 in non-band, band
    !     CNH4S,CNH4B=NH4 concentration in non-band, band
    !     ZHKM=Km for NH4 uptake
    !     FNH4,FNB4=fractions of total NH4 demand in non-band, band
    !     ZNH4S,ZNH4B=NH4 amount in non-band, band
    !     RNNH4,RNNHB=NH3 oxidation in non-band, band
    !     RGOMP=O2-unlimited respiration
    !     ECNH=efficiency CO2 conversion to biomass
    !     RVMX4,RVMXB=nitrifier demand for NH4 in non-band, band
    !
    ECHZAutor(NGL) = EO2X
    VMXX           = VMXNH3Oxi*GrowthEnvScalAutor(NGL)*FBiomNutStoiScalAutor(NGL)*XCO2*OMActAutor(NGL)
    IF(VOLWZ.GT.ZEROS2)THEN
      VMAX=VMXX/(1.0_r8+VMXX/(VHKI*VOLWZ))
    ELSE
      VMAX=0.0_r8
    ENDIF

    FCN4S                  = FNH4S*CNH3S/(CNH3S+ZHKM)
    FCN4B                  = FNHBS*CNH3B/(CNH3B+ZHKM)
    FSBSTAutor(NGL)        = FCN4S+FCN4B
    !non-band soil NH3 uptake
    VMX4S                  = VMAX*FCN4S
    !banded soil NH3 uptake
    VMX4B                  = VMAX*FCN4B
    !substrate-limited non-band soil NH3 uptake
    RNNH4                  = AZMAX1(AMIN1(VMX4S,FNH4*ZNH4S))*ZNFN4S
    !substrate-limited banded soil NH3 uptake
    RNNHB                  = AZMAX1(AMIN1(VMX4B,FNB4*ZNH4B))*ZNFN4B
    !total NH3 uptake
    RVOXP                  = RNNH4+RNNHB
    RVOXPA                 = RNNH4
    RVOXPB                 = RNNHB

    !NH3+CO2-> CH2O+O2->CO2, a fraction ECNH of CO2 taken up is converted into CH2O, of which ECHZAutor(NGL) is respired
    !ECNH energy transfer efficiency, BELSER et al. (1984) reported the ratio is about 0.09
    !ECNH means the fraction of NH3 used to generate CH2O
    !RVOXP*ECNH represents CH2O produced driven by energy from RVOXP
    !RGOMP represents the potential respiraiton from burning the CH2O generated above.
    RGOMP                  = RVOXP*ECNH*ECHZAutor(NGL)  !relevant CO2 needs to be taken up
    RNH3OxidAutor(NGL)     = VMX4S
    RNH3OxidAutorBand(NGL) = VMX4B
    !
    !     O2 DEMAND FROM NH3 OXIDATION
    !
    !     RO2Dmnd4RespHeter=O2 demand from respiration by nitrifiers
    !     ROXYP,RO2Dmnd4RespHeter=O2 demand from respiration + NH3 oxidation
    ! C+O2 -> CO2,  respiration for growth and maintenance, 2.667=32./12.
    ! NH3+1.5O2-> NO2(-)+H2O+H(+), 1.5*32/14.=3.249, energy for CO2 reduction into biomass
    !
    RO2Dmnd4GrossRespAutor(NGL) = 2.667_r8*RGOMP
    RO2MetaDmndAutor(NGL)       = RO2Dmnd4GrossRespAutor(NGL)+3.429_r8*RVOXP
    RespGrossAutor(NGL)         = RGOMP     !this is CO2 production before O2 limitation
    RSMetaOxidSoilAutor(NGL)    = RVOXPA    !NH3 oxidation in non-band soil
    RSMetaOxidBandAutor(NGL)    = RVOXPB    !NH3 oxidation in band soil
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine AmmoniaOxidizerCatabolism
!------------------------------------------------------------------------------------------
  subroutine NitriteOxidizerCatabolism(I,J,N,RMOMK,TOMEAutoKC,micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)
  !
  !nitrite oxidation
  !NO2(-) + 0.5O2 -> NO3(-),
  implicit none
  integer,  intent(in) :: I,J,N
  real(r8), intent(in) :: TOMEAutoKC
  real(r8), intent(in) :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='NitriteOxidizerCatabolism'
  real(r8) :: FNH4S,FNHBS
  real(r8) :: RGOMP,RVOXP
  real(r8) :: RVOXPA,RVOXPB
  REAL(R8) :: fno2,FNB2
  real(r8) :: FCN2S,FCN2B
  real(r8) :: RNNO2,RNNOB
  real(r8) :: VMX2S,VMX2B
  real(r8) :: VMAX
  REAL(r8) :: XCO2
  integer  :: NGL

!     begin_execution
  associate(                                                &
    GrowthEnvScalAutor     => nmics%GrowthEnvScalAutor,     & !Temperature and water-potential multiplier on autotrophic growth [-]
    FBiomNutStoiScalAutor  => nmics%FBiomNutStoiScalAutor,  & !Combined N/P stoichiometric multiplier on guild metabolic capacity [-]
    FSBSTAutor             => nmicdiag%FSBSTAutor,          & !Guild substrate-response factor; larger values mean less limitation [-]
    FracNO2XupAutor        => nmics%FracNO2XupAutor,        & !Staged autotrophic nitrite-competition weight [-]
    OMActAutor             => nmics%OMActAutor,             & !Active microbial C biomass by autotrophic guild
    RO2Dmnd4GrossRespAutor => nmicf%RO2Dmnd4GrossRespAutor, & !Potential O2 demand supporting autotrophic gross respiration
    RespGrossAutor         => nmicf%RespGrossAutor,         & !Gross respiration C equivalent by autotrophic guild
    ECHZAutor              => nmicf%ECHZAutor,              & !Guild respiration fraction used to convert growth respiration to C uptake [-]
    RSMetaOxidBandAutor    => nmicf%RSMetaOxidBandAutor,    & !Fertilizer-band NO2-N oxidation to NO3-N
    RSMetaOxidSoilAutor    => nmicf%RSMetaOxidSoilAutor,    & !Nonband NO2-N oxidation to NO3-N
    VLNH4                  => micfor%VLNH4,                 & !Nonband fraction for ammonium/ammonia pools and uptake capacity [-]; not referenced here
    VLNHB                  => micfor%VLNHB,                 & !Fertilizer-band fraction for ammonium/ammonia pools and uptake capacity [-]; not referenced here
    VLNO3                  => micfor%VLNO3,                 & !Nonband fraction for nitrate/nitrite pools and uptake capacity [-]
    VLNOB                  => micfor%VLNOB,                 & !Fertilizer-band fraction for nitrate/nitrite pools and uptake capacity [-]
    ZEROS                  => micfor%ZEROS,                 & !Small mass or flux threshold used by the routine
    RNO2EcoUptkSoilPrev    => micfor%RNO2EcoUptkSoilPrev,   & !Previous-hour ecosystem nonband NO2-N uptake; competition denominator
    RNO2EcoUptkBandPrev    => micfor%RNO2EcoUptkBandPrev,   & !Previous-hour ecosystem fertilizer-band NO2-N uptake; competition denominator
    JGniA                  => micpar%JGniA,                 & !First guild index for each autotrophic functional group
    JGnfA                  => micpar%JGnfA,                 & !Last guild index for each autotrophic functional group
    CCO2S                  => micstt%CCO2S,                 & !Dissolved CO2-C concentration for substrate saturation
    CNO2S                  => micstt%CNO2S,                 & !Dissolved NO2-N concentration in nonband soil
    CNO2B                  => micstt%CNO2B,                 & !Dissolved NO2-N concentration in fertilizer-band soil
    ZNO2S                  => micstt%ZNO2S,                 & !NO2-N pool in nonband soil
    ZNO2B                  => micstt%ZNO2B,                 & !NO2-N pool in fertilizer-band soil
    RNO2XupAutorPrev       => micflx%RNO2XupAutorPrev,      & !Previous-hour autotrophic nonband NO2-N redox uptake; reaction depends on functional group; used for competition
    RNO2XupAutorBandPrev   => micflx%RNO2XupAutorBandPrev,  & !Previous-hour autotrophic fertilizer-band NO2-N redox uptake; reaction depends on functional group; used for competition
    RNO2XupAutor           => micflx%RNO2XupAutor,          & !Nonband NO2-N uptake for oxidation to NO3-N
    RNO2XupAutorBand       => micflx%RNO2XupAutorBand,      & !Fertilizer-band NO2-N uptake for oxidation to NO3-N
    RO2MetaDmndAutor       => micflx%RO2MetaDmndAutor       & !Total autotrophic O2 demand from respiration and substrate oxidation
  )
!     FACTOR TO REGULATE COMPETITION FOR NO2 AMONG DIFFERENT
!     MICROBIAL POPULATIONS
!
!     FNO2=fraction of total biological demand for NO2 in non-band, band
!
!     CCO2S=aqueous CO2 concentration
!
  call PrintInfo('beg '//subname)
  XCO2=CCO2S/(CCO2S+CCKM)
  DO NGL=JGniA(N),JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle

    call StageAutotroph(NGL,N,TOMEAutoKC,micfor,nmics,nmicdiag)

    call CalcRespMaintAutor(I,J,NGL,RMOMK,micfor,micstt,micflx,nmicf,nmics)

    IF(RNO2EcoUptkSoilPrev.GT.ZEROS)THEN
      FNO2=AMAX1(FMN,RNO2XupAutorPrev(NGL)/RNO2EcoUptkSoilPrev)
    ELSE
      FNO2=AMAX1(FMN,FracNO2XupAutor(NGL)*VLNO3)
    ENDIF

    IF(RNO2EcoUptkBandPrev.GT.ZEROS)THEN
      FNB2=AMAX1(FMN,RNO2XupAutorBandPrev(NGL)/RNO2EcoUptkBandPrev)
    ELSE
      FNB2=AMAX1(FMN,FracNO2XupAutor(NGL)*VLNOB)
    ENDIF

    naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
    naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2
    !
    !     NO2 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
    !     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
    !     NO2 CONCENTRATIONS
    !
    !     ECHZ=growth respiration efficiency
    !     VMAX= non-substrate limited NH3 oxidation
    !     VMXN=specific oxidation
    !     TFNG=temperature+water limitation, FBiomNutStoiScalAutorr=N,P limitation
    !     XCO2=aqueous CO2 limitation, OMA=active biomass
    !     OMA=active biomass
    !     CNO2S,CNO2B=NO2 concentration in non-band, band
    !     ZNKM=Km for NO2 uptake
    !     FNO2,FNB2=fractions of total NO2 demand in non-band, band
    !     ZNO2S,ZNO2B=NO2 amount in non-band, band
    !     RNNO2,RNNOB=NO2 oxidation in non-band, band
    !     RGOMP=O2-unlimited respiration
    !     ECNO=efficiency CO2 conversion to biomass
    !     RNO2DmndReduxSoilHeter_vr,RNO2DmndReduxBandHeter_vr=nitrifier demand for NO2 in non-band, band
    !

    VMAX=GrowthEnvScalAutor(NGL)*FBiomNutStoiScalAutor(NGL)*XCO2*OMActAutor(NGL)*VMXNO2Oxi
    ECHZAutor(NGL)  = EO2X
    FCN2S           = VLNO3*CNO2S/(CNO2S+ZNKM)
    FCN2B           = VLNOB*CNO2B/(CNO2B+ZNKM)
    FSBSTAutor(NGL) = FCN2S+FCN2B
    VMX2S           = VMAX*FCN2S
    VMX2B           = VMAX*FCN2B
    RNNO2           = AZMAX1(AMIN1(VMX2S,FNO2*ZNO2S))
    RNNOB           = AZMAX1(AMIN1(VMX2B,FNB2*ZNO2B))
    RVOXP           = RNNO2+RNNOB   !total NO2(-) to be oxidized, including those used for creating CH2O from CO2?
    RVOXPA          = RNNO2
    RVOXPB          = RNNOB

    !CO2 uptake due to NO2 uptake, and a fraction of ECNO becomes CH2O. RVOXP*ECNO is carbon fixed using the energy from NO2(-) oxidation by O2.
    !2NO2(-)+CO2 + H2O -> CH2O + 2NO3(-)
    !ECNO is about 0.1, BELSER et al. (1984) reported the ratio is about 0.02
    !CH2O + 2.667 O2 -> CO2 + H2O, it is assumed first all NO2 taken up and react with O2 to generate energy used to produce CH2O,
    !and CH2O is oxidized by O2 to produce gross respiraiton, in this case, ECNO = NO2(-) used for CH2O /(NO2(-) used for energy + NO2(-)used for CH2O)
    RGOMP                  = RVOXP*ECNO*ECHZAutor(NGL)     !relevant CO2 needs be taken up
    RNO2XupAutor(NGL)      = VMX2S
    RNO2XupAutorBand(NGL)  = VMX2B
    !
    !     O2 DEMAND FROM NO2 OXIDATION
    !
    !     RO2Dmnd4RespHeter=O2 demand from respiration by nitrifiers
    !     ROXYP,RO2Dmnd4RespHeter=O2 demand from respiration + NO2 oxidation
    !from mole-based NO2(-) + 0.5O2 -> NO3(-), have O/N=16/14=1.143
    RO2Dmnd4GrossRespAutor(NGL) = 2.667_r8*RGOMP
    RO2MetaDmndAutor(NGL)       = RO2Dmnd4GrossRespAutor(NGL)+1.143_r8*RVOXP
    RespGrossAutor(NGL)         = RGOMP     !CO2 production before O2 limitation
    RSMetaOxidSoilAutor(NGL)    = RVOXPA    !NO2(-)+0.5O2 -> NO3(-)
    RSMetaOxidBandAutor(NGL)    = RVOXPB    !NO2(-)+0.5O2 -> NO3(-)
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine NitriteOxidizerCatabolism


end module NitrifierMod
