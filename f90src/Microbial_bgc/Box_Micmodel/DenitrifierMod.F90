module DenitrifierMod
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use minimathmod,          only: AZMAX1, safe_adb
  use TracerIDMod
  use EcosimConst,          only: RGASC
  use NitroPars
  use MicrobeDiagTypes

  implicit none

  private
  public :: HeteroDenitrificCatabolism

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains

!------------------------------------------------------------------------------------------

  subroutine HeteroDenitrificCatabolism(N,K,VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)

  !Description
  !FACULTATIVE denitrifcation
  !(CH2O)6  + 6O2 -> 6CO2 +6 H2O
  !(CH2O)6 + 12NO3(-) -> 6CO2 + 12NO2(-) + 6H2O, 12*14/(6*32) = 7/8=0.875
  !(CH2O)6 + 12NO2(-) -> 6CO2 + 6N2O + 12OH(-), 12*14/(6*32)=7/8=0.875
  !(CH2O)6 + 12N2O    -> 6CO2 + 12N2 + 6H2O,  24*14/(6*32)=7/4 = 1.75
  !Denitrifiers do not use acetate (which is not right)
  !the reduction of NO2 into NO is not considered
  !Ref: The microbial nitrogen-cycling network, Kuypers et al., 2018
  !maintenance respiration is computed in aerobicBacteria module. 
  implicit none
  integer, intent(in) :: N,K
  real(r8), intent(in) :: VOLWZ            !volume of water to support biogeochemistry
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_State_type),intent(inout) :: ncplxs
  type(micfluxtype), intent(inout) :: micflx
  character(len=*), parameter :: subname='HeteroDenitrificCatabolism'
  integer :: NGL
  real(r8) :: FNO3S,FNO3B
  real(r8) :: FNO2S,FNO2B
  REAL(R8) :: FNO2,FNB2
  real(r8) :: FNO3,FNB3
  real(r8) :: FVMXDX
  REAL(R8) :: FN2O
  real(r8) :: OQCZ3
  real(r8) :: OQCD3
  real(r8) :: OQCD3S
  real(r8) :: OQCD3B
  real(r8) :: OQCZ2
  real(r8) :: OQCD2
  real(r8) :: OQCD2S
  real(r8) :: OQCD2B
  real(r8) :: OQCZ1
  real(r8) :: OQCD1
  real(r8) :: ROXYD
  real(r8) :: RNO3UptkSoil,RNO3UptkBand,RDNOX
  real(r8) :: RDNOT
  real(r8) :: RGOM3X
  real(r8) :: RNO2UptkSoil,RNO2UptkBand,RDN2X,RDN2T,RGOM2X,RDN2OX,RGOM1X
  real(r8) :: RNOxDOCReduxRespDenitLim1,RNOxDOCReduxRespDenitLim2,RNOxDOCReduxRespDenitLim3
  real(r8) :: RNOxAcetReduxRespDenitLim3,RNOxAcetReduxRespDenitLim2,RNOxAcetReduxRespDenitLim1
  real(r8) :: VMXD3
  real(r8) :: VMXDXS
  real(r8) :: VMXDXB
  real(r8) :: VMXDXT
  real(r8) :: VMXD3S,VMXD3B,VMXD2,VMXD2S,VMXD2B,VMXD1
  real(r8) :: VMXD1S
  real(r8) :: ZNO3SX,ZNO3BX
  real(r8) :: ZNO2SX,ZNO2BX
  real(r8) :: Z2OSX,FODC,OQAZ1,OQAZ2,OQAZ3
  real(r8), parameter :: eQNO3toOxy=12._r8/28._r8   !2NO3(-)+CH2O-> 2NO2(-) + CO2 + H2O, NO3(-) as e- acceptor by denitrifiers
  real(r8), parameter :: eQNO2toOxy=12._r8/28._r8   !2NO2(-)+CH2O-> N2O + CO2 + 2OH(-),  NO2(-) as e- acceptor by denitrifiers
  real(r8), parameter :: eQN2OtoOxy=12._r8/56._r8   !2N2O + CH2O -> 2N2+CO2+H2O, N2O as e- acceptor by denitrifiers

! begin_execution
  associate(                                                         &
    OxyLimterHeter             => nmics%OxyLimterHeter,              & !Actual/potential O2 uptake ratio; 1 means no O2 restriction [-]
    FracOMActHeter             => nmics%FracOMActHeter,              & !Guild/complex fraction of total active microbial C in the layer [-]
    FracNO2ReduxHeter          => nmics%FracNO2ReduxHeter,           & !Guild/complex share of active nitrite-consuming biomass [-]
    RO2Dmnd4RespHeter          => nmicf%RO2Dmnd4RespHeter,           & !Potential O2 demand supporting heterotrophic gross respiration
    RO2Uptk4RespHeter          => nmicf%RO2Uptk4RespHeter,           & !Realized O2 uptake attributed to heterotrophic gross respiration
    RNO3ReduxHeterSoil         => nmicf%RNO3ReduxHeterSoil,          & !Realized nonband NO3-N reduction to NO2-N by heterotrophic guilds
    RNO3ReduxHeterBand         => nmicf%RNO3ReduxHeterBand,          & !Realized fertilizer-band NO3-N reduction to NO2-N by heterotrophic guilds
    RNO2ReduxHeterSoil         => nmicf%RNO2ReduxHeterSoil,          & !Realized nonband NO2-N reduction to N2O-N by heterotrophic guilds
    RNO2ReduxHeterBand         => nmicf%RNO2ReduxHeterBand,          & !Realized fertilizer-band NO2-N reduction to N2O-N by heterotrophic guilds
    RN2OReduxHeter             => nmicf%RN2OReduxHeter,              & !Realized N2O-N reduction to N2 by guild and complex
    RNOxDOCReduxRespDenitLim   => nmicf%RNOxDOCReduxRespDenitLim,    & !DOC-C respiration supported by realized denitrification
    RNOxAcetReduxRespDenitLim  => nmicf%RNOxAcetReduxRespDenitLim,   & !Acetate-C respiration supported by realized denitrification
    RNOxReduxRespDenitUlm      => nmicf%RNOxReduxRespDenitUlm,       & !Potential denitrification-supported respiration before acceptor limitation
    FOQC                       => nmicf%FOQC,                        & !Guild share of DOC demand used to allocate the donor pool [-]
    FOQA                       => nmicf%FOQA,                        & !Guild share of acetate demand used to allocate the donor pool [-]
    RGOCP                      => nmicf%RGOCP,                       & !DOC-supported potential respiration before O2 limitation
    RGOAP                      => nmicf%RGOAP,                       & !Acetate-supported potential respiration before O2 limitation
    BulkSOMC                   => ncplxs%BulkSOMC,                   & !Colonized solid, residue and sorbed C summed for each complex K; not referenced here
    RNO2EcoUptkSoilPrev        => micfor%RNO2EcoUptkSoilPrev,        & !Previous-hour ecosystem nonband NO2-N uptake; competition denominator
    RN2OEcoUptkSoilPrev        => micfor%RN2OEcoUptkSoilPrev,        & !Previous-hour ecosystem N2O-N uptake; competition denominator
    RNO3EcoDmndSoilPrev        => micfor%RNO3EcoDmndSoilPrev,        & !Previous-hour ecosystem NO3-N demand in nonband soil; competition denominator
    VLNO3                      => micfor%VLNO3,                      & !Nonband fraction for nitrate/nitrite pools and uptake capacity [-]
    ZERO                       => micfor%ZERO,                       & !Small dimensionless or concentration threshold used by the routine
    ZEROS                      => micfor%ZEROS,                      & !Small mass or flux threshold used by the routine
    ZEROS2                     => micfor%ZEROS2,                     & !Small water-volume threshold used by the routine
    RNO2EcoUptkBandPrev        => micfor%RNO2EcoUptkBandPrev,        & !Previous-hour ecosystem fertilizer-band NO2-N uptake; competition denominator
    RNO3EcoDmndBandPrev        => micfor%RNO3EcoDmndBandPrev,        & !Previous-hour ecosystem NO3-N demand in fertilizer-band soil; competition denominator
    VLNOB                      => micfor%VLNOB,                      & !Fertilizer-band fraction for nitrate/nitrite pools and uptake capacity [-]
    CNO3B                      => micstt%CNO3B,                      & !Dissolved NO3-N concentration in fertilizer-band soil
    CNO3S                      => micstt%CNO3S,                      & !Dissolved NO3-N concentration in nonband soil
    CZ2OS                      => micstt%CZ2OS,                      & !Dissolved N2O-N concentration for denitrification
    Z2OS                       => micstt%Z2OS,                       & !Dissolved N2O-N donor pool
    ZNO2B                      => micstt%ZNO2B,                      & !NO2-N pool in fertilizer-band soil
    ZNO2S                      => micstt%ZNO2S,                      & !NO2-N pool in nonband soil
    ZNO3B                      => micstt%ZNO3B,                      & !NO3-N pool in fertilizer-band soil
    ZNO3S                      => micstt%ZNO3S,                      & !NO3-N pool in nonband soil
    CNO2B                      => micstt%CNO2B,                      & !Dissolved NO2-N concentration in fertilizer-band soil
    CNO2S                      => micstt%CNO2S,                      & !Dissolved NO2-N concentration in nonband soil
    FracBulkSOMC               => micstt%FracBulkSOMC,               & !Complex share of colonized solid, residue and sorbed C [-]
    CH2GS                      => micstt%CH2GS,                      & !Dissolved H2 concentration used in energy-yield and saturation calculations; not referenced here
    DOM                        => micstt%DOM,                        & !Dissolved organic pools by species (DOC, DON, DOP, acetate) and complex K
    RNO3ReduxDmndSoilHeterPrev => micflx%RNO3ReduxDmndSoilHeterPrev, & !Previous-hour potential nonband NO3-N reduction demand used for competition; used for competition
    RNO3ReduxDmndBandHeterPrev => micflx%RNO3ReduxDmndBandHeterPrev, & !Previous-hour potential fertilizer-band NO3-N reduction demand used for competition; used for competition
    RNO2DmndReduxSoilHeterPrev => micflx%RNO2DmndReduxSoilHeterPrev, & !Previous-hour potential nonband NO2-N reduction demand used for competition; used for competition
    RNO2DmndReduxBandHeterPrev => micflx%RNO2DmndReduxBandHeterPrev, & !Previous-hour potential fertilizer-band NO2-N reduction demand used for competition; used for competition
    RN2ODmndReduxHeterPrev     => micflx%RN2ODmndReduxHeterPrev,     & !Previous-hour guild N2O-N reduction demand used for competition
    RNO3ReduxDmndSoilHeter     => micflx%RNO3ReduxDmndSoilHeter,     & !Potential nonband NO3-N reduction demand used for competition
    RNO2DmndReduxSoilHeter     => micflx%RNO2DmndReduxSoilHeter,     & !Potential nonband NO2-N reduction demand used for competition
    RN2ODmndReduxHeter         => micflx%RN2ODmndReduxHeter,         & !Guild N2O-N reduction demand for denitrification before donor limitation
    RNO2DmndReduxBandHeter     => micflx%RNO2DmndReduxBandHeter,     & !Potential fertilizer-band NO2-N reduction demand used for competition
    RNO3ReduxDmndBandHeter     => micflx%RNO3ReduxDmndBandHeter      & !Potential fertilizer-band NO3-N reduction demand used for competition
  )

  !
  ! FACTOR TO CONSTRAIN NO3 UPAKE AMONG COMPETING MICROBIAL
  ! AND ROOT POPULATIONS
  !
  ! FNO3,FNB3=fraction of total biological demand for NO3
  !

  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    IF(RO2Dmnd4RespHeter(NGL,K).LE.0.0_r8)cycle

    FNO3S = VLNO3
    FNO3B = VLNOB
    IF(RNO3EcoDmndSoilPrev.GT.ZEROS)THEN
      FNO3=AMAX1(FMN,RNO3ReduxDmndSoilHeterPrev(NGL,K)/RNO3EcoDmndSoilPrev)
    ELSE
      FNO3=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNO3)
    ENDIF

    IF(RNO3EcoDmndBandPrev.GT.ZEROS)THEN
      FNB3=AMAX1(FMN,RNO3ReduxDmndBandHeterPrev(NGL,K)/RNO3EcoDmndBandPrev)
    ELSE
      FNB3=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNOB)
    ENDIF
    naqfdiag%TFNO3X=naqfdiag%TFNO3X+FNO3
    naqfdiag%TFNO3B=naqfdiag%TFNO3B+FNB3
    !
    !     NO3 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
    !     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO3
    !     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
    !     NOT ACCEPTED BY O2 IN BAND AND NON-BAND SOIL ZONES
    !
    !     ROXYD=O2 demand RO2Dmnd4RespHeter not met by O2 uptake RO2Uptk4RespHeter
    !     VMXD3=demand for NO3-N reduction
    !     VMXDXS,VMXDXB=maximum NO3 reduction in non-band, band
    !     FNO3S,FNO3B=fractions of total NO3 in non-band, band
    !     CNO3S,CNO3B=NO3 concentrations in non-band, band
    !     Z3KM,Z2KM=Km for NO3, NO2 uptake
    !     FVMXDX=nonlinear effect of product inhibition for NOx reduction
    !     VMKI=product inhibition for NOx reduction
    !     VMXD3S,VMXD3B=substrate-unlimited NO3 reduction in non-band,band
    !     OQCD3S,OQCD3B=DOC limitation to NO3 reduction in non-band, band
    !     RNO3ReduxHeterSoil,RNO3ReduxHeterBand=substrate-limited NO3 reduction in non-band,band
    !     RGOM3X,RNOxDOCReduxRespDenitLim3=substrate-unltd,-ltd respn from NO3 reduction
    !     RNO3ReduxDmndSoilHeter,RNO3ReduxDmndBandHeter=demand for NO3 reduction in non-band,band
    !non-band soil
    !oxygen deficit
    ROXYD = AZMAX1(RO2Dmnd4RespHeter(NGL,K)-RO2Uptk4RespHeter(NGL,K))

    !NO3 demand
    VMXD3 = 0.875_r8*ROXYD
    IF(CNO3S.GT.ZERO)THEN
      VMXDXS=FNO3S*VMXD3*CNO3S/(CNO3S+Z3KM)/(1.0_r8+(CNO2S*Z3KM)/(CNO3S*Z2KM))
    ELSE
      VMXDXS=0.0_r8
    ENDIF
    !band-soil
    IF(CNO3B.GT.ZERO)THEN
      VMXDXB=FNO3B*VMXD3*CNO3B/(CNO3B+Z3KM)/(1.0_r8+(CNO2B*Z3KM)/(CNO3B*Z2KM))
    ELSE
      VMXDXB=0.0_r8
    ENDIF

    VMXDXT=VMXDXS+VMXDXB
    IF(VOLWZ.GT.ZEROS2 .AND. FracBulkSOMC(K).GT.ZERO)THEN
      FVMXDX=1.0_r8/(1.0_r8+VMXDXT/(VMKI*VOLWZ*FracBulkSOMC(K)))
    ELSE
      FVMXDX=0.0_r8
    ENDIF

    VMXD3S                        = VMXDXS*FVMXDX
    VMXD3B                        = VMXDXB*FVMXDX
    RNO3ReduxDmndSoilHeter(NGL,K) = VMXD3S
    RNO3ReduxDmndBandHeter(NGL,K) = VMXD3B

    !DOC to be oxidized by NO3(-)
    OQCZ3  = AZMAX1(DOM(idom_doc,K)*FOQC(NGL,K)-RGOCP(NGL,K)*OxyLimterHeter(NGL,K))  !maximum DOC available for oxidation by NO3(-)
    OQAZ3  = AZMAX1(DOM(idom_acetate,K)*FOQA(NGL,K)-RGOAP(NGL,K)*OxyLimterHeter(NGL,K))
    OQCD3  = (OQCZ3+OQAZ3)/eQNO3toOxy     !NO3-N demand for DOC oxidation
    OQCD3S = OQCD3*FNO3S                  !NO3-N-soil demand for DOC oxidation
    OQCD3B = OQCD3*FNO3B                  !NO3-N-band demand for DOC oxidation
    FODC   = safe_adb(OQCZ3,OQCZ3+OQAZ3)
    ZNO3SX                        = ZNO3S*FNO3
    ZNO3BX                        = ZNO3B*FNB3
    RNO3UptkSoil                  = AZMAX1(AMIN1(ZNO3SX,VMXD3S))     !substrate-limited uptake in soil
    RNO3UptkBand                  = AZMAX1(AMIN1(ZNO3BX,VMXD3B))     !substrate-limited uptake in band
    !apply oxidation-demand limitation
    RNO3ReduxHeterSoil(NGL,K)     = AZMAX1(AMIN1(RNO3UptkSoil,OQCD3S))   !NO3-N-soil demand for DOC oxidation, NO3-NO2
    RNO3ReduxHeterBand(NGL,K)     = AZMAX1(AMIN1(RNO3UptkBand,OQCD3B))   !NO3-N-band demand for DOC oxidation
    RDNOX                         = RNO3UptkSoil+RNO3UptkBand
    RDNOT                         = RNO3ReduxHeterSoil(NGL,K)+RNO3ReduxHeterBand(NGL,K)
    RGOM3X                        = eQNO3toOxy*RDNOX                     !NO3-N supported potential DOC oxidation
    RNOxDOCReduxRespDenitLim3     = eQNO3toOxy*RDNOT*FODC                !NO3-N supported actual DOC oxidation
    RNOxAcetReduxRespDenitLim3    = eQNO3toOxy*RDNOT*(1._r8-FODC)        !NO3-N supported actual acetate oxidation
    !
    !     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
    !     POPULATIONS
    !
    !     FNO2,FNB2=fraction of total biological demand for NO2
    !
    FNO2S=VLNO3
    FNO2B=VLNOB
    IF(RNO2EcoUptkSoilPrev.GT.ZEROS)THEN
      FNO2=AMAX1(FMN,RNO2DmndReduxSoilHeterPrev(NGL,K)/RNO2EcoUptkSoilPrev)
    ELSE
      FNO2=AMAX1(FMN,FracNO2ReduxHeter(NGL,K)*VLNO3)
    ENDIF

    IF(RNO2EcoUptkBandPrev.GT.ZEROS)THEN
      FNB2=AMAX1(FMN,RNO2DmndReduxBandHeterPrev(NGL,K)/RNO2EcoUptkBandPrev)
    ELSE
      FNB2=AMAX1(FMN,FracNO2ReduxHeter(NGL,K)*VLNOB)
    ENDIF

    naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
    naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2
    !
    !     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
    !     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2
    !     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
    !     NOT ACCEPTED BY O2 AND NO3 IN BAND AND NON-BAND SOIL ZONES
    !
    !     VMXD2=demand for NO2-N reduction
    !     VMXDXS,VMXDXB=maximum NO2 reduction in non-band, band
    !     FNO2S,FNO2B=fractions of total NO2 in non-band, band
    !     CNO2S,CNO2B=NO2 concentrations in non-band, band
    !     Z2KM,Z1KM=Km for NO2, N2O uptake
    !     FVMXDX=nonlinear effect of product inhibition for NOx reduction
    !     VMKI=product inhibition for NOx reduction
    !     VMXD2S,VMXD2B=substrate-unlimited NO2 reduction in non-band,band
    !     OQCD2S,OQCD2B=DOC limitation to NO2 reduction in non-band, band
    !     RNO2ReduxHeterSoil,RNO2ReduxHeterBand=substrate-limited NO2 reduction in non-band,band
    !     RGOM2X,RNOxDOCReduxRespDenitLim2=substrate-unltd,-ltd respn from NO2 reduction
    !NO2(-) demand for C oxidation
    VMXD2=VMXD3-RDNOT
    IF(CNO2S.GT.ZERO)THEN
      VMXDXS=FNO2S*VMXD2*CNO2S/(CNO2S+Z2KM)/(1.0_r8+(CZ2OS*Z2KM)/(CNO2S*Z1KM))
    ELSE
      VMXDXS=0.0_r8
    ENDIF
    IF(CNO2B.GT.ZERO)THEN
      VMXDXB=FNO2B*VMXD2*CNO2B/(CNO2B+Z2KM)/(1.0_r8+(CZ2OS*Z2KM)/(CNO2B*Z1KM))
    ELSE
      VMXDXB=0.0_r8
    ENDIF

    VMXDXT=VMXDXS+VMXDXB
    IF(VOLWZ.GT.ZEROS2.AND.FracBulkSOMC(K).GT.ZERO)THEN
      FVMXDX=1.0_r8/(1.0_r8+VMXDXT/(VMKI*VOLWZ*FracBulkSOMC(K)))
    ELSE
      FVMXDX=0.0_r8
    ENDIF
    VMXD2S                        = VMXDXS*FVMXDX
    VMXD2B                        = VMXDXB*FVMXDX
    OQCZ2                         = AZMAX1(OQCZ3-RNOxDOCReduxRespDenitLim3)  !maximum available DOC for oxidation by NO2(-)
    OQAZ2                         = AZMAX1(OQAZ3-RNOxAcetReduxRespDenitLim3)
    FODC                          = safe_adb(OQCZ2,OQCZ2+OQAZ2)
    OQCD2                         = (OQCZ2+OQAZ2)/eQNO2toOxy
    OQCD2S                        = OQCD2*FNO3S
    OQCD2B                        = OQCD2*FNO3B
    ZNO2SX                        = (ZNO2S+RNO3ReduxHeterSoil(NGL,K))*FNO2
    ZNO2BX                        = (ZNO2B+RNO3ReduxHeterBand(NGL,K))*FNB2
    RNO2UptkSoil                  = AZMAX1(AMIN1(ZNO2SX,VMXD2S))
    RNO2UptkBand                  = AZMAX1(AMIN1(ZNO2BX,VMXD2B))
    RNO2ReduxHeterSoil(NGL,K)     = AZMAX1(AMIN1(VMXD2S,OQCD2S,ZNO2SX))
    RNO2ReduxHeterBand(NGL,K)     = AZMAX1(AMIN1(VMXD2B,OQCD2B,ZNO2BX))
    RDN2X                         = RNO2UptkSoil+RNO2UptkBand
    RDN2T                         = RNO2ReduxHeterSoil(NGL,K)+RNO2ReduxHeterBand(NGL,K)
    RGOM2X                        = eQNO2toOxy*RDN2X
    RNOxDOCReduxRespDenitLim2     = eQNO2toOxy*RDN2T*FODC
    RNOxAcetReduxRespDenitLim2    = eQNO2toOxy*RDN2T*(1._r8-FODC)
    RNO2DmndReduxSoilHeter(NGL,K) = VMXD2S
    RNO2DmndReduxBandHeter(NGL,K) = VMXD2B
    !
    !     FACTOR TO CONSTRAIN N2O UPAKE AMONG COMPETING MICROBIAL
    !     AND ROOT POPULATIONS
    !
    !     FN2O=fraction of total biological demand for N2O
    !
    IF(RN2OEcoUptkSoilPrev.GT.ZEROS)THEN
      FN2O=AMAX1(FMN,RN2ODmndReduxHeterPrev(NGL,K)/RN2OEcoUptkSoilPrev)
    ELSE
      FN2O=AMAX1(FMN,FracOMActHeter(NGL,K))
    ENDIF
    naqfdiag%TFN2OX=naqfdiag%TFN2OX+FN2O
    !
    !     N2O REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
    !     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS N2O
    !     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
    !     NOT ACCEPTED BY O2, NO3 AND NO2 IN BAND AND NON-BAND SOIL ZONES
    !
    !     VMXD1=demand for N2O-N reduction
    !     VMXDXS=maximum N2O reduction
    !     CZ2OS=N2O concentrations
    !     Z1KM=Km for N2O uptake
    !     FVMXDX=nonlinear effect of product inhibition for NOx reduction
    !     VMKI=product inhibition for NOx reduction
    !     VMXD1S=substrate-unlimited N2O reduction
    !     OQCD1=DOC limitation to N2O reduction
    !     RDN2O=substrate-limited N2O reduction
    !     RGOM1X,RNOxDOCReduxRespDenitLim1=substrate-unltd,-ltd  respn from N2O reduction
    !     RNOxReduxRespDenitUlm,RNOxDOCReduxRespDenitLim=total substrate-unltd,-ltd respn from NOx reduction
    !     RN2ODmndReduxHeter=demand for N2O reduction
    !
    VMXD1  = (VMXD2-RDN2T)*2.0_r8
    VMXDXS = VMXD1*CZ2OS/(CZ2OS+Z1KM)
    IF(VOLWZ.GT.ZEROS2 .AND. FracBulkSOMC(K).GT.ZERO)THEN
      FVMXDX=1.0_r8/(1.0_r8+VMXDXS/(VMKI*VOLWZ*FracBulkSOMC(K)))
    ELSE
      FVMXDX=0.0_r8
    ENDIF
    
    VMXD1S                           = VMXDXS*FVMXDX
    OQCZ1                            = AZMAX1(OQCZ2-RNOxDOCReduxRespDenitLim2)   !maximum available DOC for oxidation by N2O
    OQAZ1                            = AZMAX1(OQAZ2-RNOxAcetReduxRespDenitLim2)   !maximum available DOC for oxidation by N2O
    FODC                             = safe_adb(OQCZ1,OQCZ1+OQAZ1)
    OQCD1                            = (OQCZ1+OQAZ1)/eQN2OtoOxy
    Z2OSX                            = (Z2OS+RDN2T)*FN2O
    RDN2OX                           = AZMAX1(AMIN1(Z2OSX,VMXD1S))
    RN2OReduxHeter(NGL,K)            = AZMAX1(AMIN1(VMXD1S,OQCD1,Z2OSX))
    RGOM1X                           = eQN2OtoOxy*RDN2OX
    RNOxDOCReduxRespDenitLim1        = eQN2OtoOxy*RN2OReduxHeter(NGL,K)*FODC
    RNOxAcetReduxRespDenitLim1       = eQN2OtoOxy*RN2OReduxHeter(NGL,K)*(1._r8-FODC)
    RNOxReduxRespDenitUlm(NGL,K)     = RGOM3X+RGOM2X+RGOM1X
    RNOxDOCReduxRespDenitLim(NGL,K)  = RNOxDOCReduxRespDenitLim3+RNOxDOCReduxRespDenitLim2+RNOxDOCReduxRespDenitLim1
    RNOxAcetReduxRespDenitLim(NGL,K) = RNOxAcetReduxRespDenitLim3+RNOxAcetReduxRespDenitLim2+RNOxAcetReduxRespDenitLim1
    RN2ODmndReduxHeter(NGL,K)        = VMXD1S
  ENDDO
  end associate
  end subroutine HeteroDenitrificCatabolism
end module DenitrifierMod
