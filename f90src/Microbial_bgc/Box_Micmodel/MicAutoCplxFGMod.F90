module MicAutoCPLXMod
! USES:
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use minimathmod,          only: safe_adb, AZMAX1, fixEXConsumpFlux
  use MicForcTypeMod,       only: micforctype
  use MicFluxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: PrintInfo
  use MicrobeDiagTypes
  use ElmIDMod
  use TracerIDMod
  use EcoSIMSolverPar
  use EcoSimConst
  use NitroPars
  use MicrobMathFuncMod
  use MethanogenMod,        only: H2MethanogensCatabolism
  use MethanotrophMod,      only: AeroMethanotrophCatabolism, AMONC10Catabolism, &
                                  AMOANME2dCatabolism
  use NitrifierMod,         only: AmmoniaOxiDenitCatabolism, &
                                  AmmoniaOxidizerCatabolism, NitriteOxidizerCatabolism
  implicit none

  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: ActiveAutotrophs
  public :: AutotrophAnabolicUpdate
  contains

!------------------------------------------------------------------------------------------
  subroutine ActiveAutotrophs(I,J,N,SPOMK, RMOMK, &
    micfor,micstt,micflx,naqfdiag,nmicf,nmics,ncplxf,ncplxs,nmicdiag)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: N
  real(r8), intent(in) :: SPOMK(2)
  real(r8), intent(in)  :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout) :: ncplxs
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='ActiveAutotrophs'
  integer  :: M
  real(r8) :: COMC
  real(r8) :: FOQC,FOQA
  real(r8) :: RGOCP
  real(r8) :: RGOMP
  real(r8) :: RVOXP
  real(r8) :: RVOXPA    !oxidation in soil, for NH3, NO2 and CH4
  real(r8) :: RVOXPB    !oxidation in band, specifically for NH3 and NO2
  real(r8) :: RGrowthRespAutor
  real(r8) :: RMaintDefcitcitAutor
  real(r8) :: RMaintRespAutor

! begin_execution
  associate(                                                         &
    TOMEAutoK                  => ncplxs%TOMEAutoK,                  & !Total active autotrophic C/N/P summed over guilds
    ZNH4T                      => nmicdiag%ZNH4T,                    & !NH4-N pool in band plus nonband soil
    ZNO3T                      => nmicdiag%ZNO3T,                    & !NO3-N pool in band plus nonband soil
    ZNO2T                      => nmicdiag%ZNO2T,                    & !NO2-N pool in band plus nonband soil
    H2P4T                      => nmicdiag%H2P4T,                    & !H2PO4-P pool in band plus nonband soil
    H1P4T                      => nmicdiag%H1P4T,                    & !HPO4-P pool in band plus nonband soil
    VOLWZ                      => nmicdiag%VOLWZ,                    & !Effective water volume supporting microbial activity and rate constraints
    mid_AutoAmmoniaOxidBacter  => micpar%mid_AutoAmmoniaOxidBacter,  & !Functional-group identifier for ammonia oxidizers
    mid_AutoNitriteOxidBacter  => micpar%mid_AutoNitriteOxidBacter,  & !Functional-group identifier for nitrite oxidizers
    mid_AutoH2GenoCH4GenArchea => micpar%mid_AutoH2GenoCH4GenArchea, & !Functional-group identifier for hydrogenotrophic methanogens
    mid_AutoAeroCH4OxiBacter   => micpar%mid_AutoAeroCH4OxiBacter,   & !Functional-group identifier for aerobic methane oxidizers
    mid_AutoAMOANME2D          => micpar%mid_AutoAMOANME2D       ,   & !Functional-group identifier for nitrate-dependent ANME-2d methanotrophs
    mid_AutoAMONC10            => micpar%mid_AutoAMONC10           , & !Functional-group identifier for nitrite-dependent NC10 methanotrophs
    ZEROS                      => micfor%ZEROS,                      & !Small mass or flux threshold used by the routine
    SoilMicPMassLayer          => micfor%SoilMicPMassLayer,          & !Soil mass associated with the current layer micropore domain; not referenced here
    litrm                      => micfor%litrm,                      & !True for the surface litter layer
    VLSoilPoreMicP             => micfor%VLSoilPoreMicP              & !Layer micropore volume used in water and aerobic-uptake calculations
  )
  !
  !
  !     RESPIRATION RATES BY AUTOTROPHS 'RGOMP' FROM SPECIFIC
  !     OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
  !     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL
  !     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
  !     MICROBIAL COMPETITION FACTOR. N=(1) NH4 OXIDIZERS (2) NO2
  !     OXIDIZERS,(3) CH4 OXIDIZERS, (5) H2TROPHIC METHANOGENS
  !
  !
  call PrintInfo('beg '//subname)
  if (N.eq.mid_AutoAmmoniaOxidBacter)then
    ! NH3 OXIDIZERS
    call AmmoniaOxidizerCatabolism(I,J,N,RMOMK,TOMEAutoK(ielmc),VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)

  elseif (N.eq.mid_AutoNitriteOxidBacter)then
    ! NO2 OXIDIZERS
    call NitriteOxidizerCatabolism(I,J,N,RMOMK,TOMEAutoK(ielmc),micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)

  elseif (N.eq.mid_AutoH2GenoCH4GenArchea)then
    ! H2TROPHIC METHANOGENS
    call H2MethanogensCatabolism(I,J,N,RMOMK,TOMEAutoK(ielmc),micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)

  elseif (N.eq.mid_AutoAeroCH4OxiBacter)then
    ! METHANOTROPHS
    call AeroMethanotrophCatabolism(I,J,N,RMOMK,TOMEAutoK(ielmc),micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)
  elseif (N.eq.mid_AutoAMONC10)then
    call AMONC10Catabolism(I,J,N,RMOMK,TOMEAutoK(ielmc),micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)    
  elseif (N.eq.mid_AutoAMOANME2D)then
    call AMOANME2dCatabolism(I,J,N,RMOMK,TOMEAutoK(ielmc),VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)
    
  ENDIF

  IF(micpar%is_aerobic_autor(N))then
    call AerobicAutorO2Uptake(I,J,N,micfor,micstt,nmicf,nmics,micflx,naqfdiag)
  ENDif
  !
  !
  ! RO2UptkHeter, ROXYP=O2-limited, O2-unlimited rates of O2 uptake
  ! RUPMX=O2-unlimited rate of O2 uptake
  ! FOXYX=fraction of O2 uptake by N,K relative to total
  ! dts_gas=1/(NPH*NPT)
  ! ROXYF,ROXYL=net O2 gaseous, aqueous fluxes from previous hour
  ! O2AquaDiffusvity=aqueous O2 diffusivity
  ! OXYG,OXYS=gaseous, aqueous O2 amounts
  ! Rain2LitRSurf,Irrig2LitRSurf_col=surface water flux from precipitation, irrigation
  ! O2_rain_conc,O2_irrig_conc=O2 concentration in Rain2LitRSurf,Irrig2LitRSurf_col
  !
  !
  !     AUTOTROPHIC DENITRIFICATION by ammonia-oxidizer
  !
  IF(N.EQ.mid_AutoAmmoniaOxidBacter .AND. (.not.litrm .OR. VLSoilPoreMicP.GT.ZEROS))THEN  
    call AmmoniaOxiDenitCatabolism(I,J,N,VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,micflx)
  ENDIF
  !
  !     BIOMASS DECOMPOSITION AND MINERALIZATION
  !
  ! FACTORS CONSTRAINING DOC, ACETATE, O2, NH4, NO3, PO4 UPTAKE
  ! AMONG COMPETING MICROBIAL AND ROOT POPULATIONS IN SOIL LAYERS
  !

  call BiomNutMinerMobilAutor(I,J,N,ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,micflx,nmicf,nmics,naqfdiag)
  !
  call GatherAutotrophRespiration(I,J,N,micfor,micflx,nmicf,nmics)
  !
  call GatherAutotrophAnabolicFlux(I,J,N,micflx,spomk,rmomk,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)
  call PrintInfo('end '//subname)
  end associate
  end subroutine ActiveAutotrophs


!------------------------------------------------------------------------------------------

  subroutine SubstrateCompetAuto(NGL,N,FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    micfor,naqfdiag,nmicf,nmics,micflx)
  !
  !Description:
  !Substrate competition for autotrophs
  !  
  implicit none
  integer, intent(in) :: NGL,N
  real(r8), intent(out):: FNH4X
  real(r8),intent(out) :: FNB3X,FNB4X,FNO3X
  real(r8),intent(out) :: FPO4X,FPOBX,FP14X,FP1BX
  type(micforctype), intent(in) :: micfor
  type(Cumlate_Flux_Diag_type),INTENT(INOUT)::  naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx
! begin_execution
  associate(                                                   &
    FracOMActAutor          => nmics%FracOMActAutor,           & !Guild fraction of total active microbial C in the layer [-]
    AttenfNH4Autor          => micflx%AttenfNH4Autor,          & !Litter-microbial share of NH4-N uptake from underlying soil [-]
    AttenfNO3Autor          => micflx%AttenfNO3Autor,          & !Litter-microbial share of NO3-N uptake from underlying soil [-]
    AttenfH1PO4Autor        => micflx%AttenfH1PO4Autor,        & !Litter-microbial share of HPO4-P uptake from underlying soil [-]
    AttenfH2PO4Autor        => micflx%AttenfH2PO4Autor,        & !Litter-microbial share of H2PO4-P uptake from underlying soil [-]
    RNH4EcoDmndSoilPrev     => micfor%RNH4EcoDmndSoilPrev,     & !Previous-hour ecosystem NH4-N demand in nonband soil; competition denominator
    RNH4EcoDmndBandPrev     => micfor%RNH4EcoDmndBandPrev,     & !Previous-hour ecosystem NH4-N demand in fertilizer-band soil; competition denominator
    RNO3EcoDmndSoilPrev     => micfor%RNO3EcoDmndSoilPrev,     & !Previous-hour ecosystem NO3-N demand in nonband soil; competition denominator
    RNO3EcoDmndBandPrev     => micfor%RNO3EcoDmndBandPrev,     & !Previous-hour ecosystem NO3-N demand in fertilizer-band soil; competition denominator
    RH2PO4EcoDmndSoilPrev   => micfor%RH2PO4EcoDmndSoilPrev,   & !Previous-hour ecosystem H2PO4-P demand in nonband soil; competition denominator
    RH2PO4EcoDmndBandPrev   => micfor%RH2PO4EcoDmndBandPrev,   & !Previous-hour ecosystem H2PO4-P demand in fertilizer-band soil; competition denominator
    RH1PO4EcoDmndSoilPrev   => micfor%RH1PO4EcoDmndSoilPrev,   & !Previous-hour ecosystem HPO4-P demand in nonband soil; competition denominator
    RH1PO4EcoDmndBandPrev   => micfor%RH1PO4EcoDmndBandPrev,   & !Previous-hour ecosystem HPO4-P demand in fertilizer-band soil; competition denominator
    RNH4EcoDmndLitrPrev     => micfor%RNH4EcoDmndLitrPrev,     & !Previous-hour ecosystem NH4-N demand in underlying soil accessed by litter microbes; competition denominator
    RNO3EcoDmndLitrPrev     => micfor%RNO3EcoDmndLitrPrev,     & !Previous-hour ecosystem NO3-N demand in underlying soil accessed by litter microbes; competition denominator
    RH1PO4EcoDmndLitrPrev   => micfor%RH1PO4EcoDmndLitrPrev,   & !Previous-hour ecosystem HPO4-P demand in underlying soil accessed by litter microbes; competition denominator
    RH2PO4EcoDmndLitrPrev   => micfor%RH2PO4EcoDmndLitrPrev,   & !Previous-hour ecosystem H2PO4-P demand in underlying soil accessed by litter microbes; competition denominator
    RDOMEcoDmndPrev         => micfor%RDOMEcoDmndPrev,         & !Previous-hour ecosystem DOC demand in each complex; competition denominator; not referenced here
    RAcetateEcoDmndPrev     => micfor%RAcetateEcoDmndPrev,     & !Previous-hour ecosystem acetate demand in each complex; competition denominator; not referenced here
    SoilMicPMassLayer0      => micfor%SoilMicPMassLayer0,      & !Surface-litter soil-mass reference used in litter/soil exchange conditions; not referenced here
    Lsurf                   => micfor%Lsurf,                   & !True for the surface soil layer beneath litter; not referenced here
    litrm                   => micfor%litrm,                   & !True for the surface litter layer
    ZEROS                   => micfor%ZEROS,                   & !Small mass or flux threshold used by the routine
    VLNO3                   => micfor%VLNO3,                   & !Nonband fraction for nitrate/nitrite pools and uptake capacity [-]
    VLNOB                   => micfor%VLNOB,                   & !Fertilizer-band fraction for nitrate/nitrite pools and uptake capacity [-]
    VLPO4                   => micfor%VLPO4,                   & !Nonband fraction for phosphate pools and uptake capacity [-]
    VLPOB                   => micfor%VLPOB,                   & !Fertilizer-band fraction for phosphate pools and uptake capacity [-]
    VLNHB                   => micfor%VLNHB,                   & !Fertilizer-band fraction for ammonium/ammonia pools and uptake capacity [-]
    VLNH4                   => micfor%VLNH4,                   & !Nonband fraction for ammonium/ammonia pools and uptake capacity [-]
    RH1PO4UptkLitrAutorPrev => micflx%RH1PO4UptkLitrAutorPrev, & !Previous-hour potential HPO4-P uptake by autotrophic guilds from underlying soil accessed by litter microbes; used for competition
    RH2PO4UptkLitrAutorPrev => micflx%RH2PO4UptkLitrAutorPrev, & !Previous-hour potential H2PO4-P uptake by autotrophic guilds from underlying soil accessed by litter microbes; used for competition
    RNO3UptkLitrAutorPrev   => micflx%RNO3UptkLitrAutorPrev,   & !Previous-hour potential NO3-N uptake by autotrophic guilds from underlying soil accessed by litter microbes; used for competition
    RNH4UptkLitrAutorPrev   => micflx%RNH4UptkLitrAutorPrev,   & !Previous-hour potential NH4-N uptake by autotrophic guilds from underlying soil accessed by litter microbes; used for competition
    RH1PO4UptkBandAutorPrev => micflx%RH1PO4UptkBandAutorPrev, & !Previous-hour potential HPO4-P uptake by autotrophic guilds from fertilizer-band soil; used for competition
    RH1PO4UptkSoilAutorPrev => micflx%RH1PO4UptkSoilAutorPrev, & !Previous-hour potential HPO4-P uptake by autotrophic guilds from nonband soil; used for competition
    RH2PO4UptkBandAutorPrev => micflx%RH2PO4UptkBandAutorPrev, & !Previous-hour potential H2PO4-P uptake by autotrophic guilds from fertilizer-band soil; used for competition
    RH2PO4UptkSoilAutorPrev => micflx%RH2PO4UptkSoilAutorPrev, & !Previous-hour potential H2PO4-P uptake by autotrophic guilds from nonband soil; used for competition
    RNO3UptkBandAutorPrev   => micflx%RNO3UptkBandAutorPrev,   & !Previous-hour potential NO3-N uptake by autotrophic guilds from fertilizer-band soil; used for competition
    RNO3UptkSoilAutorPrev   => micflx%RNO3UptkSoilAutorPrev,   & !Previous-hour potential NO3-N uptake by autotrophic guilds from nonband soil; used for competition
    RNH4UptkBandAutorPrev   => micflx%RNH4UptkBandAutorPrev,   & !Previous-hour potential NH4-N uptake by autotrophic guilds from fertilizer-band soil; used for competition
    RNH4UptkSoilAutorPrev   => micflx%RNH4UptkSoilAutorPrev    & !Previous-hour potential NH4-N uptake by autotrophic guilds from nonband soil; used for competition
  )
! F*=fraction of substrate uptake relative to total uptake from
! previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
! NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
! POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
! oxidation, OQA=acetate oxidation
!
  
  IF(RNH4EcoDmndSoilPrev.GT.ZEROS)THEN
    FNH4X=AMAX1(FMN,RNH4UptkSoilAutorPrev(NGL)/RNH4EcoDmndSoilPrev)
  ELSE
    FNH4X=AMAX1(FMN,FracOMActAutor(NGL)*VLNH4)
  ENDIF
  IF(RNH4EcoDmndBandPrev.GT.ZEROS)THEN
    FNB4X=AMAX1(FMN,RNH4UptkBandAutorPrev(NGL)/RNH4EcoDmndBandPrev)
  ELSE
    FNB4X=AMAX1(FMN,FracOMActAutor(NGL)*VLNHB)
  ENDIF
  IF(RNO3EcoDmndSoilPrev.GT.ZEROS)THEN
    FNO3X=AMAX1(FMN,RNO3UptkSoilAutorPrev(NGL)/RNO3EcoDmndSoilPrev)
  ELSE
    FNO3X=AMAX1(FMN,FracOMActAutor(NGL)*VLNO3)
  ENDIF
  IF(RNO3EcoDmndBandPrev.GT.ZEROS)THEN
    FNB3X=AMAX1(FMN,RNO3UptkBandAutorPrev(NGL)/RNO3EcoDmndBandPrev)
  ELSE
    FNB3X=AMAX1(FMN,FracOMActAutor(NGL)*VLNOB)
  ENDIF
  IF(RH2PO4EcoDmndSoilPrev.GT.ZEROS)THEN
    FPO4X=AMAX1(FMN,RH2PO4UptkSoilAutorPrev(NGL)/RH2PO4EcoDmndSoilPrev)
  ELSE
    FPO4X=AMAX1(FMN,FracOMActAutor(NGL)*VLPO4)
  ENDIF
  IF(RH2PO4EcoDmndBandPrev.GT.ZEROS)THEN
    FPOBX=AMAX1(FMN,RH2PO4UptkBandAutorPrev(NGL)/RH2PO4EcoDmndBandPrev)
  ELSE
    FPOBX=AMAX1(FMN,FracOMActAutor(NGL)*VLPOB)
  ENDIF
  IF(RH1PO4EcoDmndSoilPrev.GT.ZEROS)THEN
    FP14X=AMAX1(FMN,RH1PO4UptkSoilAutorPrev(NGL)/RH1PO4EcoDmndSoilPrev)
  ELSE
    FP14X=AMAX1(FMN,FracOMActAutor(NGL)*VLPO4)
  ENDIF
  IF(RH1PO4EcoDmndBandPrev.GT.ZEROS)THEN
    FP1BX=AMAX1(FMN,RH1PO4UptkBandAutorPrev(NGL)/RH1PO4EcoDmndBandPrev)
  ELSE
    FP1BX=AMAX1(FMN,FracOMActAutor(NGL)*VLPOB)
  ENDIF

!  naqfdiag%TFNH4X = naqfdiag%TFNH4X+FNH4X
!  naqfdiag%TFNO3X = naqfdiag%TFNO3X+FNO3X
!  naqfdiag%TFPO4X = naqfdiag%TFPO4X+FPO4X
!  naqfdiag%TFP14X = naqfdiag%TFP14X+FP14X
!  naqfdiag%TFNH4B = naqfdiag%TFNH4B+FNB4X
!  naqfdiag%TFNO3B = naqfdiag%TFNO3B+FNB3X
!  naqfdiag%TFPO4B = naqfdiag%TFPO4B+FPOBX
!  naqfdiag%TFP14B = naqfdiag%TFP14B+FP1BX
  !
  ! FACTORS CONSTRAINING NH4, NO3, PO4 UPTAKE AMONG COMPETING
  ! MICROBIAL POPULATIONS IN SURFACE RESIDUE
  ! F*=fraction of substrate uptake relative to total uptake from
  ! previous hour in surface litter, labels as for soil layers above
  !
  !litter layer
  !All litter complexes and autotrophs tap the same underlying-soil pool.
  !Use their shared layer biomass denominator when previous demand is absent.
  IF(litrm)THEN
    IF(RNH4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfNH4Autor(NGL)=AMAX1(FMN,RNH4UptkLitrAutorPrev(NGL)/RNH4EcoDmndLitrPrev)
    ELSE
      AttenfNH4Autor(NGL)=AMAX1(FMN,FracOMActAutor(NGL))
    ENDIF
    IF(RNO3EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfNO3Autor(NGL)=AMAX1(FMN,RNO3UptkLitrAutorPrev(NGL)/RNO3EcoDmndLitrPrev)
    ELSE
      AttenfNO3Autor(NGL)=AMAX1(FMN,FracOMActAutor(NGL))
    ENDIF
    IF(RH2PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfH2PO4Autor(NGL)=AMAX1(FMN,RH2PO4UptkLitrAutorPrev(NGL)/RH2PO4EcoDmndLitrPrev)
    ELSE
      AttenfH2PO4Autor(NGL)=AMAX1(FMN,FracOMActAutor(NGL))
    ENDIF
    IF(RH1PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfH1PO4Autor(NGL)=AMAX1(FMN,RH1PO4UptkLitrAutorPrev(NGL)/RH1PO4EcoDmndLitrPrev)
    ELSE
      AttenfH1PO4Autor(NGL)=AMAX1(FMN,FracOMActAutor(NGL))
    ENDIF
  ENDIF
  !top soil layer
  !diagnostics off
!  IF(Lsurf.AND.SoilMicPMassLayer0.GT.ZEROS)THEN
!    naqfdiag%TFNH4X=naqfdiag%TFNH4X+micfor%AttenfNH4AutorR(NGL)
!    naqfdiag%TFNO3X=naqfdiag%TFNO3X+micfor%AttenfNO3AutorR(NGL)
!    naqfdiag%TFPO4X=naqfdiag%TFPO4X+micfor%AttenfH2PO4AutorR(NGL)
!    naqfdiag%TFP14X=naqfdiag%TFP14X+micfor%AttenfH1PO4AutorR(NGL)
!  ENDIF
  end associate
  end subroutine SubstrateCompetAuto
!------------------------------------------------------------------------------------------

  subroutine GatherAutotrophAnabolicFlux(I,J,N,micflx,spomk,rmomk,micfor,micstt,nmicf,nmics,ncplxf,ncplxs)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: N
  real(r8), intent(in) :: spomk(2)
  real(r8), intent(in) :: RMOMK(2)
  type(MicForcType), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(in) :: micflx  
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout) :: ncplxs

  character(len=*), parameter :: subname='GatherAutotrophAnabolicFlux'
  integer :: M,K,MID3,MID,MID1,NE,idom,NGL
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: CCC,CGOMX,CGOMD
  real(r8) :: FracDenitResp4Maint,MaintDenitResp
  real(r8) :: CXC,RCCE(NumPlantChemElms)
  real(r8) :: CGOXC
  real(r8) :: C3C,CNC,CPC
  real(r8) :: CGOMZ
  real(r8) :: SPOMX
  real(r8) :: FRM
  real(r8) :: AvailableBiomass(NumPlantChemElms)
!     begin_execution
  associate(                                                                    &
    rCNBiomeActAutor                 => nmics%rCNBiomeActAutor,                 & !Active autotrophic biomass nutrient:C ratios (N:C and P:C)
    GrowthEnvScalAutor               => nmics%GrowthEnvScalAutor,               & !Temperature and water-potential multiplier on autotrophic growth [-]
    OMActAutor                       => nmics%OMActAutor,                       & !Active microbial C biomass by autotrophic guild
    DOMuptk4GrothAutor               => nmicf%DOMuptk4GrothAutor,               & !Guild elemental uptake; C source is CO2 or CH4 according to functional group
    NonstX2stBiomAutor               => nmicf%NonstX2stBiomAutor,               & !C/N/P transfer from guild reserves into kinetic and structural biomass
    RespGrossAutor                   => nmicf%RespGrossAutor,                   & !Gross respiration C equivalent by autotrophic guild
    RNOxReduxRespAutorLim            => nmicf%RNOxReduxRespAutorLim,            & !C-equivalent respiration supported by autotrophic nitrite reduction
    RMaintDmndAutor                  => nmicf%RMaintDmndAutor,                  & !Maintenance-C demand by live biomass compartment and autotrophic guild
    RkillLitfalOMAutor               => nmicf%RkillLitfalOMAutor,               & !Unrecycled ordinary-mortality C/N/P from autotrophic biomass
    RkillLitrfal2HumOMAutor          => nmicf%RkillLitrfal2HumOMAutor,          & !Ordinary-mortality C/N/P routed to humified material from autotrophic biomass
    RkillLitrfal2ResduOMAutor        => nmicf%RkillLitrfal2ResduOMAutor,        & !Ordinary-mortality C/N/P routed to microbial residue from autotrophic biomass
    RMaintDefcitLitrfalOMAutor       => nmicf%RMaintDefcitLitrfalOMAutor,       & !Unrecycled starvation-derived C/N/P from autotrophic biomass
    RMaintDefLitrfal2HumOMAutor      => nmicf%RMaintDefLitrfal2HumOMAutor,      & !Starvation-derived C/N/P routed to humified material from autotrophic biomass
    RMaintDefLitrfal2ResduOMAutor    => nmicf%RMaintDefLitrfal2ResduOMAutor,    & !Starvation-derived C/N/P routed to microbial residue from autotrophic biomass
    RKillOMAutor                     => nmicf%RKillOMAutor,                     & !Ordinary mortality C/N/P withdrawal from autotrophic biomass
    RkillRecycOMAutor                => nmicf%RkillRecycOMAutor,                & !Ordinary-mortality C/N/P recycled to reserves from autotrophic biomass
    RMaintDefcitKillOMAutor          => nmicf%RMaintDefcitKillOMAutor,          & !Maintenance-starvation C/N/P withdrawal from autotrophic biomass
    RMaintDefcitRecycOMAutor         => nmicf%RMaintDefcitRecycOMAutor,         & !Starvation recycling: C respired, N/P returned to reserves from autotrophic biomass
    Resp4NFixAutor                   => nmicf%Resp4NFixAutor,                   & !Autotrophic respiration-C cost of N2 fixation (currently set to zero)
    ECHZAutor                        => nmicf%ECHZAutor,                        & !Guild respiration fraction used to convert growth respiration to C uptake [-]
    rNCOMCAutor                      => micpar%rNCOMCAutor,                     & !Target autotrophic N:C ratios by compartment and guild
    rPCOMCAutor                      => micpar%rPCOMCAutor,                     & !Target autotrophic P:C ratios by compartment and guild
    JGniA                            => micpar%JGniA,                           & !First guild index for each autotrophic functional group
    JGnfA                            => micpar%JGnfA,                           & !Last guild index for each autotrophic functional group
    FL                               => micpar%FL,                              & !Target fractions of active biomass in kinetic and structural compartments [-]
    ZEROS                            => micfor%ZEROS,                           & !Small mass or flux threshold used by the routine
    ZERO                             => micfor%ZERO,                            & !Small dimensionless or concentration threshold used by the routine
    RGrowthRespAutor                 => micflx%RGrowthRespAutor,                & !Autotrophic gross respiration remaining after maintenance
    RMaintDefcitcitAutor             => micflx%RMaintDefcitcitAutor,            & !Guild maintenance-C deficit after available gross respiration
    RMaintRespAutor                  => micflx%RMaintRespAutor,                 & !Total hourly autotrophic guild maintenance-C demand
    mBiomeAutor                      => micstt%mBiomeAutor,                     & !C/N/P pools indexed by element and flattened guild/biomass compartment
    EHUM                             => micstt%EHUM                             & !Fraction of microbial litterfall routed to humified organic matter [-]
  )
  call PrintInfo('beg '//subname)
  !     DOC, DON, DOP AND ACETATE UPTAKE DRIVEN BY GROWTH RESPIRATION
  !     FROM O2, NOX AND C REDUCTION
  !
  !     CGOMX=DOC+acetate uptake from aerobic growth respiration
  !     CGOMD=DOC+acetate uptake from denitrifier growth respiration
  !     RMaintRespAutor=maintenance respiration
  !     RespGrossHeter=total respiration
  !     RNOxReduxRespDenitLim=respiration for denitrifcation
  !     Resp4NFixHeter=respiration for N2 fixation
  !     ECHZ,ENOX=growth respiration efficiencies for CO2, O2, and NOx reduction
  !     CGOMC,CGOQC,CGOAC=total DOC+acetate, DOC, acetate uptake heterotrophs
  !     CGOMC=total CO2,CH4 uptake (autotrophs)
  !     CGOMN,CGOMP=DON, DOP uptake
  !     FGOCP,FGOAP=DOC,acetate/(DOC+acetate)
  !     OQN,OPQ=DON,DOP
  !     FOMK=faction of OMA in total OMA
  !     FCN,FCP=limitation from N,P
  !
  DO NGL=JGniA(N),JGnfA(N)  
    IF(OMActAutor(NGL).LE.0.0_r8)cycle        

    !potential growth respiraiton-respiraiton for N2-fixation 
    DOMuptk4GrothAutor(idom_beg:idom_end,NGL)=0._r8
    CGOMX = AMIN1(RMaintRespAutor(NGL),RespGrossAutor(NGL))+Resp4NFixAutor(NGL)+(RGrowthRespAutor(NGL)-Resp4NFixAutor(NGL))/ECHZAutor(NGL)
    call ReserveDenitrifMaintenance(RMaintRespAutor(NGL),RespGrossAutor(NGL), &
      RNOxReduxRespAutorLim(NGL),EO2X,ENOX,FracDenitResp4Maint)
    !CO2 supplies maintenance respiration plus uptake for the remaining growth.
    MaintDenitResp=RNOxReduxRespAutorLim(NGL)*FracDenitResp4Maint
    CGOMD=MaintDenitResp+(RNOxReduxRespAutorLim(NGL)-MaintDenitResp)/ENOX

    !C entering the biomass/respiration pathway, supplied by CO2 or CH4.
    !For aerobic methanotrophs, this is CH4 uptake for maintenance and growth.
    !For H2 methanogens, this is CO2-derived C entering the biomass/respiration
    !pathway; it excludes direct CO2-to-CH4 conversion represented by RVOXP.
    DOMuptk4GrothAutor(ielmc,NGL)=CGOMX+CGOMD

    !
    !     TRANSFER UPTAKEN C,N,P FROM STORAGE TO ACTIVE BIOMASS
    !
    !     OMC,OMN,OMP=nonstructural C,N,P
    !     CCC,CNC,CPC=C:N:P ratios used to calculate C,N,P recycling
    !     rNCOMC,rPCOMC=maximum microbial N:C, P:C ratios
    !     RCCC,RCCN,RCCP=C,N,P recycling fractions
    !     RCCZ,RCCY=min, max C recycling fractions
    !     RCCX,RCCQ=max N,P recycling fractions
    !
    MID1=micpar%get_micb_id(iLbiom_kinetic,NGL);MID3=micpar%get_micb_id(iLbiom_reserve,NGL)
    IF(mBiomeAutor(ielmc,MID3).GT.ZEROS.AND.mBiomeAutor(ielmc,MID1).GT.ZEROS)THEN
      CCC=AZMAX1(AMIN1(1.0_r8 &
        ,mBiomeAutor(ielmn,MID3)/(mBiomeAutor(ielmn,MID3)+mBiomeAutor(ielmc,MID3)*rNCOMCAutor(iLbiom_reserve,NGL)) &
        ,mBiomeAutor(ielmp,MID3)/(mBiomeAutor(ielmp,MID3)+mBiomeAutor(ielmc,MID3)*rPCOMCAutor(iLbiom_reserve,NGL))))
      CXC = mBiomeAutor(ielmc,MID3)/mBiomeAutor(ielmc,MID1)
      C3C = 1.0_r8/(1.0_r8+CXC/CKC)
      CNC = AZMAX1(AMIN1(1.0_r8 &
        ,mBiomeAutor(ielmc,MID3)/(mBiomeAutor(ielmc,MID3)+mBiomeAutor(ielmn,MID3)/rNCOMCAutor(iLbiom_reserve,NGL))))
      CPC=AZMAX1(AMIN1(1.0_r8 &
        ,mBiomeAutor(ielmc,MID3)/(mBiomeAutor(ielmc,MID3)+mBiomeAutor(ielmp,MID3)/rPCOMCAutor(iLbiom_reserve,NGL))))
      RCCC = RCCZ+AMAX1(CCC,C3C)*RCCY
      RCCN = CNC*RCCX
      RCCP = CPC*RCCQ
    ELSE
      RCCC = RCCZ
      RCCN = 0.0_r8
      RCCP = 0.0_r8
    ENDIF
    RCCE(ielmc)=RCCC
    RCCE(ielmn)=(RCCN+(1.0_r8-RCCN)*RCCC)
    RCCE(ielmp)=(RCCP+(1.0_r8-RCCP)*RCCC)    
    !
    !     MICROBIAL ASSIMILATION OF NONSTRUCTURAL C,N,P
    !
    !     CGOMZ=transfer from nonstructural to structural microbial C
    !     TFNG=temperature+water stress function
    !     OMGR=rate constant for transferring nonstructural to structural C
    !     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
    !     FL=partitioning between labile and resistant microbial components
    !     OMC,OMN,OMP=nonstructural microbial C,N,P
    !
    MID3  = micpar%get_micb_id(iLbiom_reserve,NGL)
    CGOMZ = GrowthEnvScalAutor(NGL)*OMGR*AZMAX1(mBiomeAutor(ielmc,MID3))

    DO M = 1, 2
      NonstX2stBiomAutor(ielmc,M,NGL)=FL(M)*CGOMZ
      IF(mBiomeAutor(ielmc,MID3).GT.ZEROS)THEN
        NonstX2stBiomAutor(ielmn,M,NGL)=AMIN1(FL(M)*AZMAX1(mBiomeAutor(ielmn,MID3)) &
          ,NonstX2stBiomAutor(ielmc,M,NGL)*mBiomeAutor(ielmn,MID3)/mBiomeAutor(ielmc,MID3))
        NonstX2stBiomAutor(ielmp,M,NGL)=AMIN1(FL(M)*AZMAX1(mBiomeAutor(ielmp,MID3)) &
          ,NonstX2stBiomAutor(ielmc,M,NGL)*mBiomeAutor(ielmp,MID3)/mBiomeAutor(ielmc,MID3))
      ELSE
        NonstX2stBiomAutor(ielmn,M,NGL)=0.0_r8
        NonstX2stBiomAutor(ielmp,M,NGL)=0.0_r8
      ENDIF

    !
    !     MICROBIAL DECOMPOSITION FROM BIOMASS, SPECIFIC DECOMPOSITION
    !     RATE, TEMPERATURE
    !
    !     SPOMX=rate constant for microbial decomposition
    !     SPOMC=basal decomposition rate
    !     SPOMK=effect of low microbial C concentration on microbial decay
    !     RXOMC,RXOMN,RXOMP=microbial C,N,P decomposition
    !     RDOMC,RDOMN,RDOMP=microbial C,N,P LitrFall
    !     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling
    !
      MID   = micpar%get_micb_id(M,NGL)
      SPOMX = SQRT(GrowthEnvScalAutor(NGL))*SPOMC(M)*SPOMK(M)
      
      DO NE=1,NumPlantChemElms
        RKillOMAutor(NE,M,NGL)=AZMAX1(AMIN1(mBiomeAutor(NE,MID),mBiomeAutor(NE,MID)*SPOMX))
            
        RkillRecycOMAutor(NE,M,NGL)=RKillOMAutor(NE,M,NGL)*RCCE(NE)

        RkillLitfalOMAutor(NE,M,NGL)=AZMAX1(RKillOMAutor(NE,M,NGL)-RkillRecycOMAutor(NE,M,NGL))
    !
    !     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
    !     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
    !
    !     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P LitrFall to humus
    !     EHUM=humus transfer fraction from hour1.f
    !     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P LitrFall to residue
    !
        RkillLitrfal2HumOMAutor(NE,M,NGL)=RkillLitfalOMAutor(NE,M,NGL)*EHUM
    !
    !     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
    !
        RkillLitrfal2ResduOMAutor(NE,M,NGL)=RkillLitfalOMAutor(NE,M,NGL)-RkillLitrfal2HumOMAutor(NE,M,NGL)
      ENDDO
      
    ENDDO
    
    !
    !     MICROBIAL DECOMPOSITION WHEN MAINTENANCE RESPIRATION
    !     EXCEEDS UPTAKE
    !
    !     OMC,OMN,OMP=microbial C,N,P
    !     RMaintRespAutor=total maintenance respiration
    !     RMaintDefcitcitAutor=senescence respiration
    !     RCCC=C recycling fraction
    !     RXMMC,RXMMN,RXMMP=microbial C,N,P loss from senescence
    !     RMaintDmndHeter=maintenance respiration
    !     CNOMA,CPOMA=N:C,P:C ratios of active biomass
    !     RDMMC,RDMMN,RDMMP=microbial C,N,P LitrFall from senescence
    !     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
    !
    IF(RMaintDefcitcitAutor(NGL).GT.ZEROS.AND.RMaintRespAutor(NGL).GT.ZEROS.AND.RCCC.GT.ZERO)THEN
      FRM=RMaintDefcitcitAutor(NGL)/RMaintRespAutor(NGL)
      DO  M=1,2
        !Cap C/N/P withdrawal by this compartment's own donor pools.
        MID=micpar%get_micb_id(M,NGL)
        !Ordinary mortality and starvation share this compartment's donor pools.
        AvailableBiomass=MAX(0._r8,mBiomeAutor(1:NumPlantChemElms,MID)-RKillOMAutor(1:NumPlantChemElms,M,NGL))
        RMaintDefcitKillOMAutor(ielmc,M,NGL)=AMIN1(AvailableBiomass(ielmc),AZMAX1(FRM*RMaintDmndAutor(M,NGL)/RCCC))
        RMaintDefcitKillOMAutor(ielmn,M,NGL)=AMIN1(AvailableBiomass(ielmn),AZMAX1(RMaintDefcitKillOMAutor(ielmc,M,NGL)*rCNBiomeActAutor(ielmn,NGL)))
        RMaintDefcitKillOMAutor(ielmp,M,NGL)=AMIN1(AvailableBiomass(ielmp),AZMAX1(RMaintDefcitKillOMAutor(ielmc,M,NGL)*rCNBiomeActAutor(ielmp,NGL)))
        DO NE=1,NumPlantChemElms
          RMaintDefcitRecycOMAutor(NE,M,NGL)   = RMaintDefcitKillOMAutor(NE,M,NGL)*RCCE(NE)
          RMaintDefcitLitrfalOMAutor(NE,M,NGL) = AZMAX1(RMaintDefcitKillOMAutor(NE,M,NGL)-RMaintDefcitRecycOMAutor(NE,M,NGL))
          !
          !     HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
          !     PRODUCTS
          !
          !     RHMMC,RHMMN,RHMMC=transfer of senesence LitrFall C,N,P to humus
          !     EHUM=humus transfer fraction
          !     RCMMC,RCMMN,RCMMC=transfer of senesence LitrFall C,N,P to residue
          !
          RMaintDefLitrfal2HumOMAutor(NE,M,NGL)   = RMaintDefcitLitrfalOMAutor(NE,M,NGL)*EHUM
          RMaintDefLitrfal2ResduOMAutor(NE,M,NGL) = RMaintDefcitLitrfalOMAutor(NE,M,NGL)-RMaintDefLitrfal2HumOMAutor(NE,M,NGL)
        ENDDO
      ENDDO
    ELSE
      DO  M=1,2
        DO NE=1,NumPlantChemElms
          RMaintDefcitKillOMAutor(NE,M,NGL)          = 0.0_r8
          RMaintDefcitLitrfalOMAutor(NE,M,NGL)       = 0.0_r8
          RMaintDefcitRecycOMAutor(NE,M,NGL)         = 0.0_r8
          RMaintDefLitrfal2HumOMAutor(NE,M,NGL)   = 0.0_r8
          RMaintDefLitrfal2ResduOMAutor(NE,M,NGL) = 0.0_r8
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine GatherAutotrophAnabolicFlux
!------------------------------------------------------------------------------------------

  subroutine AerobicAutorO2Uptake(I,J,N,micfor,micstt,nmicf,nmics,micflx,naqfdiag)
  implicit none
  integer, intent(in) :: I,J,N     !functional group id
  type(MicForcType), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag  
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type),intent(inout) :: nmics
  type(micfluxtype), intent(inout) :: micflx
  real(r8) :: FOXYX,OXKX
  integer  :: M,MX,NGL
  real(r8) :: COXYS1,DIFOX
  real(r8) :: B,C,O2AquaDiffusvity1
  real(r8) :: OXYG1,OXYS1
  real(r8) :: RUPMX,RO2DmndX
  real(r8) :: ROXYLX
  real(r8) :: RRADO,RMPOX,ROXDFQ
  real(r8) :: THETW1,VOLWOX
  real(r8) :: VOLPOX
  real(r8) :: X,VOLOXM
  real(r8) :: VOLWPM

  ! begin_execution
  associate(                                                       &
    fLimO2Autor               => nmics%fLimO2Autor,                & !Actual/potential O2 uptake ratio; 1 means no O2 restriction [-]
    OMActAutor                => nmics%OMActAutor,                 & !Active microbial C biomass by autotrophic guild
    FracOMActAutor            => nmics%FracOMActAutor,             & !Guild fraction of total active microbial C in the layer [-]
    RO2UptkAutor              => nmicf%RO2UptkAutor,               & !Realized total O2 uptake by autotrophic guild
    RespGrossAutor            => nmicf%RespGrossAutor,             & !Gross respiration C equivalent by autotrophic guild
    RO2Dmnd4GrossRespAutor    => nmicf%RO2Dmnd4GrossRespAutor,     & !Potential O2 demand supporting autotrophic gross respiration
    RO2Uptk4RespAutor         => nmicf%RO2Uptk4RespAutor,          & !Realized O2 uptake attributed to autotrophic gross respiration
    RCO2ProdAutor             => nmicf%RCO2ProdAutor,              & !CO2-C production by autotrophic guild
    RCH4ProdAutor             => nmicf%RCH4ProdAutor,              & !CH4-C production by autotrophic guild
    RSMetaOxidSoilAutor       => nmicf%RSMetaOxidSoilAutor,        & !Nonband catabolic substrate oxidation; substrate depends on functional group
    RSMetaOxidBandAutor       => nmicf%RSMetaOxidBandAutor,        & !Fertilizer-band catabolic substrate oxidation by autotrophic guild
    mid_AutoAmmoniaOxidBacter => micpar%mid_AutoAmmoniaOxidBacter, & !Functional-group identifier for ammonia oxidizers
    mid_AutoNitriteOxidBacter => micpar%mid_AutoNitriteOxidBacter, & !Functional-group identifier for nitrite oxidizers
    RO2GasXchangePrev         => micfor%RO2GasXchangePrev,         & !Previous-hour gaseous O2 exchange; negated when applied as a supply
    RO2MetaDmndAutorPrev      => micflx%RO2MetaDmndAutorPrev,      & !Previous-hour autotrophic guild O2 demand used for competition
    RO2MetaDmndAutor          => micflx%RO2MetaDmndAutor,          & !Total autotrophic O2 demand from respiration and substrate oxidation
    COXYE                     => micfor%COXYE,                     & !Atmospheric gas-phase O2 concentration
    RO2EcoDmndPrev            => micfor%RO2EcoDmndPrev,            & !Previous-hour ecosystem O2 demand; competition denominator
    O2_rain_conc              => micfor%O2_rain_conc,              & !Dissolved O2 concentration in rainwater
    O2_irrig_conc             => micfor%O2_irrig_conc,             & !Dissolved O2 concentration in irrigation water
    Irrig2LitRSurf_col        => micfor%Irrig2LitRSurf_col,        & !Irrigation water input to surface litter, carrying dissolved O2
    Rain2LitRSurf             => micfor%Rain2LitRSurf,             & !Rainwater input to surface litter, carrying dissolved O2
    litrm                     => micfor%litrm,                     & !True for the surface litter layer
    O2AquaDiffusvity          => micfor%O2AquaDiffusvity,          & !Aqueous O2 diffusivity before transport-substep scaling
    RO2AquaXchangePrev        => micfor%RO2AquaXchangePrev,        & !Previous-hour aqueous O2 exchange; negated when applied as a supply
    VLSoilPoreMicP            => micfor%VLSoilPoreMicP,            & !Layer micropore volume used in water and aerobic-uptake calculations
    VLSoilMicP                => micfor%VLSoilMicP,                & !Bulk volume associated with the layer micropore domain
    VLsoiAirPM                => micfor%VLsoiAirPM,                & !Soil air volume at each outer transport substep M
    VLWatMicPM                => micfor%VLWatMicPM,                & !Micropore water volume at each outer transport substep M
    FILM                      => micfor%FILM,                      & !Water-film thickness for microbial O2 diffusion at transport substep M
    THETPM                    => micfor%THETPM,                    & !Air-filled soil pore fraction at each outer transport substep M [-]
    TortMicPM                 => micfor%TortMicPM,                 & !Aqueous diffusion tortuosity factor at transport substep M [-]
    ZERO                      => micfor%ZERO,                      & !Small dimensionless or concentration threshold used by the routine
    ZEROS                     => micfor%ZEROS,                     & !Small mass or flux threshold used by the routine
    DiffusivitySolutEff       => micfor%DiffusivitySolutEff,       & !Gas-water exchange coefficient at each transport substep M
    O2GSolubility             => micstt%O2GSolubility,             & !Equilibrium aqueous-to-gas O2 concentration ratio [-]
    JGniA                     => micpar%JGniA,                     & !First guild index for each autotrophic functional group
    JGnfA                     => micpar%JGnfA,                     & !Last guild index for each autotrophic functional group
    OXYG                      => micstt%OXYG,                      & !Gas-phase O2 donor pool
    OXYS                      => micstt%OXYS,                      & !Dissolved O2 donor pool
    COXYG                     => micstt%COXYG,                     & !Soil gas-phase O2 concentration
    REcoUptkSoilO2M           => micflx%REcoUptkSoilO2M,           & !Accumulated microbial O2 uptake in each outer transport substep M
    RNH3OxidAutor             => micflx%RNH3OxidAutor,             & !Nonband ammonia-N oxidation by nitrifier guilds
    RNH3OxidAutorBand         => micflx%RNH3OxidAutorBand,         & !Fertilizer-band ammonia-N oxidation by nitrifier guilds
    RNO2XupAutor              => micflx%RNO2XupAutor,              & !Autotrophic nonband NO2-N redox uptake; reaction depends on functional group
    RNO2XupAutorBand          => micflx%RNO2XupAutorBand           & !Autotrophic fertilizer-band NO2-N redox uptake; reaction depends on functional group
  )

  DO NGL=JGniA(N),JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle    
    IF(RO2EcoDmndPrev.GT.ZEROS)THEN
      FOXYX=AMAX1(FMN,RO2MetaDmndAutorPrev(NGL)/RO2EcoDmndPrev)
    ELSE
      FOXYX=AMAX1(FMN,FracOMActAutor(NGL))
    ENDIF
    naqfdiag%TFOXYX   = naqfdiag%TFOXYX+FOXYX
    OXKX              = OXKA
    RO2UptkAutor(NGL) = 0._r8

    IF(RO2MetaDmndAutor(NGL).GT.ZEROS .AND. FOXYX.GT.ZERO)THEN
      IF(.not.litrm .OR. VLSoilPoreMicP.GT.ZEROS)THEN
        !
        !write(*,*)'MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC'
        !     POPULATION
        !
        RUPMX             = RO2MetaDmndAutor(NGL)*dts_gas    !
        RO2DmndX          = -RO2GasXchangePrev*dts_gas*FOXYX    !O2 demand
        O2AquaDiffusvity1 = O2AquaDiffusvity*dts_gas
        IF(.not.litrm)THEN
          OXYG1  = OXYG*FOXYX
          ROXYLX = -RO2AquaXchangePrev*dts_gas*FOXYX
        ELSE
          OXYG1  = COXYG*VLsoiAirPM(1)*FOXYX
          ROXYLX = -(RO2AquaXchangePrev+Rain2LitRSurf*O2_rain_conc &
            +Irrig2LitRSurf_col*O2_irrig_conc)*dts_gas*FOXYX
        ENDIF
        OXYS1=OXYS*FOXYX
        !Aqueous transport removal depends on the dissolved O2 donor pool.
        if(OXYS1<=0._r8 .and. ROXYLX>0._r8)ROXYLX=0._r8
        !
            !write(*,*)'O2 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP'
        !     TO MAINTAIN AQUEOUS O2 CONCENTRATION DURING REDUCTION
        !
        DO  M=1,NPH
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
          DIFOX  = TortMicPM(M)*O2AquaDiffusvity1*12.57_r8*BIOS*OMActAutor(NGL)*RRADO
          VOLWOX = VLWatMicPM(M)*O2GSolubility
          VOLPOX = VLsoiAirPM(M)
          VOLWPM = VOLWOX+VOLPOX
          VOLOXM = VLWatMicPM(M)*FOXYX
          !oxygen uptake in the layer
          DO  MX=1,NPT
            call fixEXConsumpFlux(OXYG1,RO2DmndX)
            call fixEXConsumpFlux(OXYS1,ROXYLX)
            COXYS1 = AMIN1(COXYE*O2GSolubility,safe_adb(OXYS1,VOLOXM))

            !solve for uptake flux
            IF(OXYS1<=ZEROS)THEN
              RMPOX=0.0_r8
            else
              RMPOX=TranspBasedsubstrateUptake(COXYS1,DIFOX, OXKX, RUPMX, ZEROS)
            ENDIF

            !Credit only O2 present in the aqueous donor during this substep.
            !O2 supplied by dissolution below is available in the next substep.
            RMPOX = MIN(MAX(0._r8,RMPOX),MAX(0._r8,OXYS1))
            OXYS1 = OXYS1-RMPOX

            !apply volatilization-dissolution
            IF(THETPM(M).GT.AirFillPore_Min.AND.VOLPOX.GT.ZEROS)THEN
              ROXDFQ=DiffusivitySolutEff(M)*(AMAX1(ZEROS,OXYG1)*VOLWOX-OXYS1*VOLPOX)/VOLWPM
              ROXDFQ=AMAX1(AMIN1(OXYG1,ROXDFQ),-OXYS1)
            ELSE
              ROXDFQ=0.0_r8
            ENDIF
            OXYG1 = OXYG1-ROXDFQ
            OXYS1 = OXYS1+ROXDFQ
            !accumulate uptake flux
            RO2UptkAutor(NGL)  = RO2UptkAutor(NGL)+RMPOX
            REcoUptkSoilO2M(M) = REcoUptkSoilO2M(M)+RMPOX
          ENDDO
        ENDDO
        !
        !     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (OxyLimterHeter)
        !
        !     OxyLimterHeter=ratio of O2-limited to O2-unlimited uptake
        !     RVMX4,RVNHB,RNO2DmndReduxSoilHeter_vr,RNO2DmndReduxBandHeter_vr=NH3,NO2 oxidation in non-band, band
        !
        fLimO2Autor(NGL)=AMIN1(1.0,AZMAX1(RO2UptkAutor(NGL)/RO2MetaDmndAutor(NGL)))
        IF(N.EQ.mid_AutoAmmoniaOxidBacter)THEN
          RNH3OxidAutor(NGL)     = RNH3OxidAutor(NGL)*fLimO2Autor(NGL)
          RNH3OxidAutorBand(NGL) = RNH3OxidAutorBand(NGL)*fLimO2Autor(NGL)
        ELSEIF(N.EQ.mid_AutoNitriteOxidBacter)THEN
          RNO2XupAutor(NGL)     = RNO2XupAutor(NGL)*fLimO2Autor(NGL)
          RNO2XupAutorBand(NGL) = RNO2XupAutorBand(NGL)*fLimO2Autor(NGL)
        ENDIF
      ELSE
        RO2UptkAutor(NGL) = RO2MetaDmndAutor(NGL)
        fLimO2Autor(NGL)  = 1.0_r8
      ENDIF
    ELSE
      fLimO2Autor(NGL)  = 1.0_r8
    ENDIF
    !
    !     RespGrossHeter,RGOMP=O2-limited, O2-unlimited respiration
    !     RCO2X,RAcetateProdHeter,RCH4ProdHeter,RH2ProdHeter=CO2,acetate,CH4,H2 production from RespGrossHeter
    !     RO2Uptk4RespHeter=O2-limited O2 uptake
    !     RSMetaOxidSoilAutor,RSMetaOxidBandAutor=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation

    !NH3 oxidizer assimilate CO2, CH4 oxidizer produces CO2, nitrite oxidizer assimilates CO2
    RespGrossAutor(NGL)      = RespGrossAutor(NGL)*fLimO2Autor(NGL)
    RCO2ProdAutor(NGL)       = RespGrossAutor(NGL)
    RCH4ProdAutor(NGL)       = 0.0_r8
    RO2Uptk4RespAutor(NGL)   = RO2Dmnd4GrossRespAutor(NGL)*fLimO2Autor(NGL)
    RSMetaOxidSoilAutor(NGL) = RSMetaOxidSoilAutor(NGL)*fLimO2Autor(NGL)
    RSMetaOxidBandAutor(NGL) = RSMetaOxidBandAutor(NGL)*fLimO2Autor(NGL)
  ENDDO
  end associate
  end subroutine AerobicAutorO2Uptake



!------------------------------------------------------------------------------------------
  subroutine BiomNutMinerMobilAutor(I,J,N,ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,micflx,nmicf,&
    nmics,naqfdiag)
  !
  !Description:
  !Do biomass nutrient mobilization and immobilization.   
  implicit none
  integer, intent(in) :: I,J,N
  real(r8), intent(in) :: ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type),intent(inout) :: nmics
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag  
  real(r8)  :: FNH4X
  real(r8)  :: FNB3X,FNB4X,FNO3X
  real(r8)  :: FPO4X,FPOBX,FP14X,FP1BX  
  real(r8) :: CNH4X,CNH4Y,CNO3X,CNO3Y
  real(r8) :: CH2PX,CH2PY
  real(r8) :: CH1PX,CH1PY
  real(r8) :: FNH4S,FNHBS
  real(r8) :: FNO3S,FNO3B
  real(r8) :: FH1PS,FH1PB
  real(r8) :: FH2PS,FH2PB
  real(r8) :: H2POM,H2PBM
  real(r8) :: H1POM,H1PBM
  real(r8) :: H1P4M,H2P4M
  real(r8) :: RNetNH4MinPotent
  real(r8) :: RINHX,RNetNO3Dmnd,RINOX,RNetH2PO4MinPotent,RIPOX,RNetH1PO4Dmnd
  real(r8) :: RIP1X,RNetNH4MinPotentLitr,RNetNO3DmndLitr,RNetH2PO4MinPotentLitr,RNetH1PO4DmndLitr
  real(r8) :: ZNH4M,ZNHBM
  real(r8) :: ZNO3M
  real(r8) :: ZNOBM
  integer :: MID3,NGL

  real(r8) :: LitrPotential(2),LitrUptake(2) !Nonband and band soil fluxes
!     begin_execution
  associate(                                             &
   GrowthEnvScalAutor    => nmics%GrowthEnvScalAutor,    & !Temperature and water-potential multiplier on autotrophic growth [-]
   OMActAutor            => nmics%OMActAutor,            & !Active microbial C biomass by autotrophic guild
   AttenfNH4Autor        => micflx%AttenfNH4Autor,       & !Litter-microbial share of NH4-N uptake from underlying soil [-]
   AttenfNO3Autor        => micflx%AttenfNO3Autor,       & !Litter-microbial share of NO3-N uptake from underlying soil [-]
   AttenfH2PO4Autor      => micflx%AttenfH2PO4Autor,     & !Litter-microbial share of H2PO4-P uptake from underlying soil [-]
   AttenfH1PO4Autor      => micflx%AttenfH1PO4Autor,     & !Litter-microbial share of HPO4-P uptake from underlying soil [-]
   RNH4TransfSoilAutor   => nmicf%RNH4TransfSoilAutor,   & !Net NH4-N transfer from nonband soil to microbes; positive immobilization
   RNO3TransfSoilAutor   => nmicf%RNO3TransfSoilAutor,   & !Net NO3-N transfer from nonband soil to microbes; positive immobilization
   RH2PO4TransfSoilAutor => nmicf%RH2PO4TransfSoilAutor, & !Net H2PO4-P transfer from nonband soil to microbes; positive immobilization
   RNH4TransfBandAutor   => nmicf%RNH4TransfBandAutor,   & !Net NH4-N transfer from fertilizer-band soil to microbes; positive immobilization
   RNO3TransfBandAutor   => nmicf%RNO3TransfBandAutor,   & !Net NO3-N transfer from fertilizer-band soil to microbes; positive immobilization
   RH2PO4TransfBandAutor => nmicf%RH2PO4TransfBandAutor, & !Net H2PO4-P transfer from fertilizer-band soil to microbes; positive immobilization
   RNH4TransfLitrAutor   => nmicf%RNH4TransfLitrAutor,   & !Net NH4-N transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
   RNO3TransfLitrAutor   => nmicf%RNO3TransfLitrAutor,   & !Net NO3-N transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
   RH2PO4TransfLitrAutor => nmicf%RH2PO4TransfLitrAutor, & !Net H2PO4-P transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
   RH1PO4TransfSoilAutor => nmicf%RH1PO4TransfSoilAutor, & !Net HPO4-P transfer from nonband soil to microbes; positive immobilization
   RH1PO4TransfBandAutor => nmicf%RH1PO4TransfBandAutor, & !Net HPO4-P transfer from fertilizer-band soil to microbes; positive immobilization
   RH1PO4TransfLitrAutor => nmicf%RH1PO4TransfLitrAutor, & !Net HPO4-P transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
   rNCOMCAutor           => micpar%rNCOMCAutor,          & !Target autotrophic N:C ratios by compartment and guild
   rPCOMCAutor           => micpar%rPCOMCAutor,          & !Target autotrophic P:C ratios by compartment and guild
   VLNH4                 => micfor%VLNH4,                & !Nonband fraction for ammonium/ammonia pools and uptake capacity [-]
   VLNHB                 => micfor%VLNHB,                & !Fertilizer-band fraction for ammonium/ammonia pools and uptake capacity [-]
   VLWatMicP             => micfor%VLWatMicP,            & !Layer micropore water volume used for nutrient donor thresholds
   VOLWU                 => micfor%VOLWU,                & !Water volume in the soil beneath surface litter
   VLNO3                 => micfor%VLNO3,                & !Nonband fraction for nitrate/nitrite pools and uptake capacity [-]
   VLNOB                 => micfor%VLNOB,                & !Fertilizer-band fraction for nitrate/nitrite pools and uptake capacity [-]
   VLPO4                 => micfor%VLPO4,                & !Nonband fraction for phosphate pools and uptake capacity [-]
   VLPOB                 => micfor%VLPOB,                & !Fertilizer-band fraction for phosphate pools and uptake capacity [-]
   litrm                 => micfor%litrm,                & !True for the surface litter layer
   mBiomeAutor           => micstt%mBiomeAutor,          & !C/N/P pools indexed by element and flattened guild/biomass compartment
   ! Underlying-soil pools and concentrations for supplemental litter uptake.
   ZNH4TU                => micstt%ZNH4TU, &               !NH4-N pool in underlying soil, band plus nonband
   ZNO3TU                => micstt%ZNO3TU, &               !NO3-N pool in underlying soil, band plus nonband
   H2P4TU                => micstt%H2P4TU, &               !H2PO4-P pool in underlying soil, band plus nonband
   H1P4TU                => micstt%H1P4TU, &               !HPO4-P pool in underlying soil, band plus nonband
   CNH4SU                => micstt%CNH4SU, &               !Dissolved NH4-N concentration in underlying nonband soil
   CNH4BU                => micstt%CNH4BU, &               !Dissolved NH4-N concentration in underlying fertilizer-band soil
   CNO3SU                => micstt%CNO3SU, &               !Dissolved NO3-N concentration in underlying nonband soil
   CNO3BU                => micstt%CNO3BU, &               !Dissolved NO3-N concentration in underlying fertilizer-band soil
   CH2P4U                => micstt%CH2P4U, &               !Dissolved H2PO4-P concentration in underlying nonband soil
   CH2P4BU               => micstt%CH2P4BU, &              !Dissolved H2PO4-P concentration in underlying fertilizer-band soil
   CH1P4U                => micstt%CH1P4U, &               !Dissolved HPO4-P concentration in underlying nonband soil
   CH1P4BU               => micstt%CH1P4BU, &              !Dissolved HPO4-P concentration in underlying fertilizer-band soil
   ZNH4S                 => micstt%ZNH4S,                & !NH4-N pool in nonband soil
   ZNH4B                 => micstt%ZNH4B,                & !NH4-N pool in fertilizer-band soil
   ZNO3S                 => micstt%ZNO3S,                & !NO3-N pool in nonband soil
   ZNO3B                 => micstt%ZNO3B,                & !NO3-N pool in fertilizer-band soil
   CNO3S                 => micstt%CNO3S,                & !Dissolved NO3-N concentration in nonband soil
   CNO3B                 => micstt%CNO3B,                & !Dissolved NO3-N concentration in fertilizer-band soil
   CNH4S                 => micstt%CNH4S,                & !Dissolved NH4-N concentration in nonband soil
   CNH4B                 => micstt%CNH4B,                & !Dissolved NH4-N concentration in fertilizer-band soil
   CH2P4                 => micstt%CH2P4,                & !Dissolved H2PO4-P concentration in nonband soil
   CH2P4B                => micstt%CH2P4B,               & !Dissolved H2PO4-P concentration in fertilizer-band soil
   H2PO4                 => micstt%H2PO4,                & !H2PO4-P pool in nonband soil
   H2POB                 => micstt%H2POB,                & !H2PO4-P pool in fertilizer-band soil
   CH1P4                 => micstt%CH1P4,                & !Dissolved HPO4-P concentration in nonband soil
   CH1P4B                => micstt%CH1P4B,               & !Dissolved HPO4-P concentration in fertilizer-band soil
   H1PO4                 => micstt%H1PO4,                & !HPO4-P pool in nonband soil
   H1POB                 => micstt%H1POB,                & !HPO4-P pool in fertilizer-band soil
   JGniA                 => micpar%JGniA,                & !First guild index for each autotrophic functional group
   JGnfA                 => micpar%JGnfA,                & !Last guild index for each autotrophic functional group
   RNH4UptkSoilAutor     => micflx%RNH4UptkSoilAutor,    & !Potential NH4-N uptake by autotrophic guilds from nonband soil
   RNH4UptkBandAutor     => micflx%RNH4UptkBandAutor,    & !Potential NH4-N uptake by autotrophic guilds from fertilizer-band soil
   RNO3UptkSoilAutor     => micflx%RNO3UptkSoilAutor,    & !Potential NO3-N uptake by autotrophic guilds from nonband soil
   RNO3UptkBandAutor     => micflx%RNO3UptkBandAutor,    & !Potential NO3-N uptake by autotrophic guilds from fertilizer-band soil
   NetNH4Mineralize      => micflx%NetNH4Mineralize,     & !Net mineral N exchange (NH4 plus NO3); positive immobilization, negative release
   RH2PO4UptkSoilAutor   => micflx%RH2PO4UptkSoilAutor,  & !Potential H2PO4-P uptake by autotrophic guilds from nonband soil
   RH2PO4UptkBandAutor   => micflx%RH2PO4UptkBandAutor,  & !Potential H2PO4-P uptake by autotrophic guilds from fertilizer-band soil
   NetPO4Mineralize      => micflx%NetPO4Mineralize,     & !Net phosphate exchange; positive immobilization, negative mineralization
   RH1PO4UptkSoilAutor   => micflx%RH1PO4UptkSoilAutor,  & !Potential HPO4-P uptake by autotrophic guilds from nonband soil
   RH1PO4UptkBandAutor   => micflx%RH1PO4UptkBandAutor,  & !Potential HPO4-P uptake by autotrophic guilds from fertilizer-band soil
   RH1PO4UptkLitrAutor   => micflx%RH1PO4UptkLitrAutor,  & !Potential HPO4-P uptake by autotrophic guilds from underlying soil accessed by litter microbes
   RNH4UptkLitrAutor     => micflx%RNH4UptkLitrAutor,    & !Potential NH4-N uptake by autotrophic guilds from underlying soil accessed by litter microbes
   RNO3UptkLitrAutor     => micflx%RNO3UptkLitrAutor,    & !Potential NO3-N uptake by autotrophic guilds from underlying soil accessed by litter microbes
   RH2PO4UptkLitrAutor   => micflx%RH2PO4UptkLitrAutor   & !Potential H2PO4-P uptake by autotrophic guilds from underlying soil accessed by litter microbes
  )
!     MINERALIZATION-IMMOBILIZATION OF NH4 IN SOIL FROM MICROBIAL
!     C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RNetNH4MinPotent=NH4 mineralization (-ve) or immobilization (+ve) demand
!     OMC,OMN=microbial nonstructural C,N
!     rNCOMC=maximum microbial N:C ratio
!     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
!     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
!     minimum NH4 concentration and Km for NH4 uptake
!     RINHX=microbially limited NH4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     VLWatMicP=water content
!     ZNH4M,ZNHBM=NH4 not available for uptake in non-band, band
!     FNH4X,FNB4X=fractions of biological NH4 demand in non-band, band
!     RNH4imobilSoilHeter,RNH4imobilBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
!     NetNH4Mineralize=total NH4 net mineraln (-ve) or immobiln (+ve)
!
  DO NGL=JGniA(N),JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle      
    call SubstrateCompetAuto(NGL,N,FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
        micfor,naqfdiag,nmicf,nmics,micflx)

    FNH4S=VLNH4
    FNHBS=VLNHB
    MID3=micpar%get_micb_id(iLbiom_reserve,NGL)
    RNetNH4MinPotent=mBiomeAutor(ielmc,MID3)*rNCOMCAutor(iLbiom_reserve,NGL)-mBiomeAutor(ielmn,MID3)
    IF(RNetNH4MinPotent.GT.0.0_r8)THEN
      CNH4X                    = AZMAX1(CNH4S-Z4MN)
      CNH4Y                    = AZMAX1(CNH4B-Z4MN)
      RINHX                    = AMIN1(RNetNH4MinPotent,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*Z4MX)
      RNH4UptkSoilAutor(NGL)   = FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
      RNH4UptkBandAutor(NGL)   = FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
      ZNH4M                    = Z4MN*VLWatMicP*FNH4S
      ZNHBM                    = Z4MN*VLWatMicP*FNHBS
      RNH4TransfSoilAutor(NGL) = AMIN1(FNH4X*AZMAX1((ZNH4S-ZNH4M)),RNH4UptkSoilAutor(NGL))
      RNH4TransfBandAutor(NGL) = AMIN1(FNB4X*AZMAX1((ZNH4B-ZNHBM)),RNH4UptkBandAutor(NGL))
    ELSE
      RNH4UptkSoilAutor(NGL)   = 0.0_r8
      RNH4UptkBandAutor(NGL)   = 0.0_r8
      RNH4TransfSoilAutor(NGL) = RNetNH4MinPotent*FNH4S
      RNH4TransfBandAutor(NGL) = RNetNH4MinPotent*FNHBS
    ENDIF
    NetNH4Mineralize=NetNH4Mineralize+(RNH4TransfSoilAutor(NGL)+RNH4TransfBandAutor(NGL))
!
!     MINERALIZATION-IMMOBILIZATION OF NO3 IN SOIL FROM MICROBIAL
!     C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RNetNO3Dmnd=NO3 immobilization (+ve) demand
!     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
!     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
!     min NO3 concentration and Km for NO3 uptake
!     RINOX=microbially limited NO3 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNO3S,FNO3B=fractions of NO3 in non-band, band
!     VLWatMicP=water content
!     ZNO3M,ZNOBM=NO3 not available for uptake in non-band, band
!     FNO3X,FNB3X=fractions of biological NO3 demand in non-band, band
!     RNO3imobilSoilHeter,RNO3imobilBandHeter=substrate-limited NO3 immobiln in non-band, band
!     NetNH4Mineralize=total net NH4+NO3 mineraln (-ve) or immobiln (+ve)
!
    FNO3S=VLNO3
    FNO3B=VLNOB
    RNetNO3Dmnd=AZMAX1(RNetNH4MinPotent-RNH4TransfSoilAutor(NGL)-RNH4TransfBandAutor(NGL))
    IF(RNetNO3Dmnd.GT.0.0_r8)THEN
      CNO3X                    = AZMAX1(CNO3S-ZOMN)
      CNO3Y                    = AZMAX1(CNO3B-ZOMN)
      RINOX                    = AMIN1(RNetNO3Dmnd,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*ZOMX)
      RNO3UptkSoilAutor(NGL)   = FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
      RNO3UptkBandAutor(NGL)   = FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
      ZNO3M                    = ZOMN*VLWatMicP*FNO3S
      ZNOBM                    = ZOMN*VLWatMicP*FNO3B
      RNO3TransfSoilAutor(NGL) = AMIN1(FNO3X*AZMAX1((ZNO3S-ZNO3M)),RNO3UptkSoilAutor(NGL))
      RNO3TransfBandAutor(NGL) = AMIN1(FNB3X*AZMAX1((ZNO3B-ZNOBM)),RNO3UptkBandAutor(NGL))
    ELSE
      RNO3UptkSoilAutor(NGL)   = 0.0_r8
      RNO3UptkBandAutor(NGL)   = 0.0_r8
      RNO3TransfSoilAutor(NGL) = RNetNO3Dmnd*FNO3S
      RNO3TransfBandAutor(NGL) = RNetNO3Dmnd*FNO3B
    ENDIF
    NetNH4Mineralize=NetNH4Mineralize+(RNO3TransfSoilAutor(NGL)+RNO3TransfBandAutor(NGL))
!
!     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SOIL FROM MICROBIAL
!     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RNetH2PO4MinPotent=H2PO4 mineralization (-ve) or immobilization (+ve) demand
!     OMC,OMP=microbial nonstructural C,P
!     rPCOMC=maximum microbial P:C ratio
!     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
!     min H2PO4 concentration and Km for H2PO4 uptake
!     RIPOX=microbially limited H2PO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
!     RH2PO4DmndSoilHeter_vr,RH2PO4DmndBandHeter_vr=substrate-unlimited H2PO4 mineraln-immobiln
!     H2POM,H2PBM=H2PO4 not available for uptake in non-band, band
!     VOLW=water content
!     FPO4X,FPOBX=fractions of biol H2PO4 demand in non-band, band
!     RH2PO4imobilSoilHeter,RH2PO4imobilBandHeter=substrate-limited H2PO4 mineraln-immobn in non-band, band
!     NetPO4Mineralize=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
    FH2PS=VLPO4
    FH2PB=VLPOB
    MID3=micpar%get_micb_id(iLbiom_reserve,NGL)
    RNetH2PO4MinPotent=(mBiomeAutor(ielmc,MID3)*rPCOMCAutor(iLbiom_reserve,NGL)-mBiomeAutor(ielmp,MID3))
    IF(RNetH2PO4MinPotent.GT.0.0)THEN
      CH2PX=AZMAX1(CH2P4-HPMN)
      CH2PY=AZMAX1(CH2P4B-HPMN)
      RIPOX=AMIN1(RNetH2PO4MinPotent,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*HPMX)
      RH2PO4UptkSoilAutor(NGL)=FH2PS*RIPOX*CH2PX/(CH2PX+HPKU)
      RH2PO4UptkBandAutor(NGL)=FH2PB*RIPOX*CH2PY/(CH2PY+HPKU)
      H2POM=HPMN*VLWatMicP*FH2PS
      H2PBM=HPMN*VLWatMicP*FH2PB
      RH2PO4TransfSoilAutor(NGL)=AMIN1(FPO4X*AZMAX1((H2PO4-H2POM)),RH2PO4UptkSoilAutor(NGL))
      RH2PO4TransfBandAutor(NGL)=AMIN1(FPOBX*AZMAX1((H2POB-H2PBM)),RH2PO4UptkBandAutor(NGL))
    ELSE
      RH2PO4UptkSoilAutor(NGL)=0.0_r8
      RH2PO4UptkBandAutor(NGL)=0.0_r8
      RH2PO4TransfSoilAutor(NGL)=RNetH2PO4MinPotent*FH2PS
      RH2PO4TransfBandAutor(NGL)=RNetH2PO4MinPotent*FH2PB
    ENDIF
    NetPO4Mineralize=NetPO4Mineralize+(RH2PO4TransfSoilAutor(NGL)+RH2PO4TransfBandAutor(NGL))
!
!     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SOIL FROM MICROBIAL
!     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
!
!     RNetH1PO4Dmnd=HPO4 mineralization (-ve) or immobilization (+ve) demand
!     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
!     min HPO4 concentration and Km for HPO4 uptake
!     RIP1X=microbially limited HPO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH1PS,FH1PB=fractions of HPO4 in non-band, band
!     RH1PO4DmndSoilHeter_vr,RH1PO4DmndBandHeter_vr=substrate-unlimited HPO4 mineraln-immobiln
!     H1POM,H1PBM=HPO4 not available for uptake in non-band, band
!     VOLW=water content
!     FP14X,FP1BX=fractions of biol HPO4 demand in non-band, band
!     RH1PO4imobilSoilHeter,RH1PO4imobilBandHeter=substrate-limited HPO4 mineraln-immobn in non-band, band
!     NetPO4Mineralize=total H2PO4+HPO4 net mineraln (-ve) or immobiln (+ve)
!
    FH1PS=VLPO4
    FH1PB=VLPOB
    RNetH1PO4Dmnd=0.1_r8*AZMAX1(RNetH2PO4MinPotent-RH2PO4TransfSoilAutor(NGL)-RH2PO4TransfBandAutor(NGL))
    IF(RNetH1PO4Dmnd.GT.0.0_r8)THEN
      CH1PX=AZMAX1(CH1P4-HPMN)
      CH1PY=AZMAX1(CH1P4B-HPMN)
      RIP1X=AMIN1(RNetH1PO4Dmnd,BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*HPMX)
      RH1PO4UptkSoilAutor(NGL)=FH1PS*RIP1X*CH1PX/(CH1PX+HPKU)
      RH1PO4UptkBandAutor(NGL)=FH1PB*RIP1X*CH1PY/(CH1PY+HPKU)
      H1POM=HPMN*VLWatMicP*FH1PS
      H1PBM=HPMN*VLWatMicP*FH1PB
      RH1PO4TransfSoilAutor(NGL)=AMIN1(FP14X*AZMAX1((H1PO4-H1POM)),RH1PO4UptkSoilAutor(NGL))
      RH1PO4TransfBandAutor(NGL)=AMIN1(FP1BX*AZMAX1((H1POB-H1PBM)),RH1PO4UptkBandAutor(NGL))
    ELSE
      RH1PO4UptkSoilAutor(NGL)=0.0_r8
      RH1PO4UptkBandAutor(NGL)=0.0_r8
      RH1PO4TransfSoilAutor(NGL)=RNetH1PO4Dmnd*FH1PS
      RH1PO4TransfBandAutor(NGL)=RNetH1PO4Dmnd*FH1PB
    ENDIF
    NetPO4Mineralize=NetPO4Mineralize+(RH1PO4TransfSoilAutor(NGL)+RH1PO4TransfBandAutor(NGL))
!
!     MINERALIZATION-IMMOBILIZATION OF NH4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RNetNH4MinPotentLitr=NH4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
!     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
!     minimum NH4 concentration and Km for NH4 uptake
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNH4S,FNHBS=fractions of NH4 in non-band, band
!     RNH4DmndLitrHeter=substrate-unlimited NH4 mineraln-immobiln
!     VOLW=water content
!     ZNH4M=NH4 not available for uptake
!     AttenfNH4Heter=fractions of biological NH4 demand
!     RNH4imobilLitrHeter=substrate-limited NH4 mineraln-immobiln
!     NetNH4Mineralize=total NH4 net mineraln (-ve) or immobiln (+ve)
!
    ! These transfers are charged to the underlying soil, not the litter pool.
    IF(litrm)THEN
      RNetNH4MinPotentLitr=RNetNH4MinPotent-RNH4TransfSoilAutor(NGL)-RNO3TransfSoilAutor(NGL) &
        -RNH4TransfBandAutor(NGL)-RNO3TransfBandAutor(NGL)

      call LitterSoilNutrientUptake(1,RNetNH4MinPotentLitr, &
        BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*Z4MX, &
        Z4KU,Z4MN,AttenfNH4Autor(NGL), &
        micflx%RNH4UptkLitrBandAutorPrev(NGL),nmics%FracOMActAutor(NGL),micfor,LitrPotential,LitrUptake)
      RNH4UptkLitrAutor(NGL)=LitrPotential(1)
      micflx%RNH4UptkLitrBandAutor(NGL)=LitrPotential(2)
      RNH4TransfLitrAutor(NGL)=SUM(LitrUptake)
      micflx%tRNH4MicrbImobilSoil=micflx%tRNH4MicrbImobilSoil+LitrUptake(1)
      micflx%tRNH4MicrbImobilBand=micflx%tRNH4MicrbImobilBand+LitrUptake(2)
      NetNH4Mineralize=NetNH4Mineralize+RNH4TransfLitrAutor(NGL)
!
!     MINERALIZATION-IMMOBILIZATION OF NO3 IN SURFACE RESIDUE FROM
!     MICROBIAL C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RNetNO3DmndLitr=NH4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CNO3SU,CNO3BU=aqueous NO3 concentrations in non-band, band
!     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
!     minimum NO3 concentration and Km for NO3 uptake
!     RNO3DmndLitrHeter_col=microbially limited NO3 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FNO3S,FNO3B=fractions of NO3 in non-band, band
!     RNO3imobilLitrHeter=substrate-unlimited NO3 immobiln
!     VOLWU=water content
!     ZNO3M=NO3 not available for uptake
!     AttenfNO3Heter=fraction of biological NO3 demand
!     RNO3imobilLitrHeter=substrate-limited NO3 immobiln
!     NetNH4Mineralize=total NH4+NO3 net mineraln (-ve) or immobiln (+ve)
!
      RNetNO3DmndLitr=AZMAX1(RNetNH4MinPotentLitr-RNH4TransfLitrAutor(NGL))

      call LitterSoilNutrientUptake(2,RNetNO3DmndLitr, &
        BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*ZOMX, &
        ZOKU,ZOMN,AttenfNO3Autor(NGL), &
        micflx%RNO3UptkLitrBandAutorPrev(NGL),nmics%FracOMActAutor(NGL),micfor,LitrPotential,LitrUptake)
      RNO3UptkLitrAutor(NGL)=LitrPotential(1)
      micflx%RNO3UptkLitrBandAutor(NGL)=LitrPotential(2)
      RNO3TransfLitrAutor(NGL)=SUM(LitrUptake)
      micflx%tRNO3MicrbImobilSoil=micflx%tRNO3MicrbImobilSoil+LitrUptake(1)
      micflx%tRNO3MicrbImobilBand=micflx%tRNO3MicrbImobilBand+LitrUptake(2)
      NetNH4Mineralize=NetNH4Mineralize+RNO3TransfLitrAutor(NGL)
!
!     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SURFACE RESIDUE FROM
!     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
!     ZONES OF SOIL SURFACE
!
!     RNetH2PO4MinPotentLitr=H2PO4 mineralization (-ve) or immobilization (+ve) demand
!     NU=surface layer number
!     CH2P4U,CH2P4BU=aqueous H2PO4 concentrations in non-band, band
!     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
!     minimum H2PO4 concentration and Km for H2PO4 uptake
!     RH2PO4DmndLitrHeter=microbially limited H2PO4 demand
!     BIOA=microbial surface area, OMA=active biomass
!     TFNG=temp+water stress
!     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
!     RH2PO4DmndLitrHeter=substrate-unlimited H2PO4 mineraln-immobiln
!     VOLWU=water content
!     H2P4M=H2PO4 not available for uptake
!     AttenfH2PO4Heter=fractions of biological H2PO4 demand
!     RH2PO4imobilLitrHeter=substrate-limited H2PO4 mineraln-immobiln
!     NetPO4Mineralize=total H2PO4 net mineraln (-ve) or immobiln (+ve)
!
      !Subtract all P already exchanged with litter before tapping topsoil.
      RNetH2PO4MinPotentLitr=RNetH2PO4MinPotent-RH2PO4TransfSoilAutor(NGL) &
        -RH2PO4TransfBandAutor(NGL)-RH1PO4TransfSoilAutor(NGL)-RH1PO4TransfBandAutor(NGL)

      call LitterSoilNutrientUptake(3,RNetH2PO4MinPotentLitr, &
        BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*HPMX, &
        HPKU,HPMN,AttenfH2PO4Autor(NGL), &
        micflx%RH2PO4UptkLitrBandAutorPrev(NGL),nmics%FracOMActAutor(NGL),micfor,LitrPotential,LitrUptake)
      RH2PO4UptkLitrAutor(NGL)=LitrPotential(1)
      micflx%RH2PO4UptkLitrBandAutor(NGL)=LitrPotential(2)
      RH2PO4TransfLitrAutor(NGL)=SUM(LitrUptake)
      micflx%tRH2PO4MicrbImobilSoil=micflx%tRH2PO4MicrbImobilSoil+LitrUptake(1)
      micflx%tRH2PO4MicrbImobilBand=micflx%tRH2PO4MicrbImobilBand+LitrUptake(2)
      NetPO4Mineralize=NetPO4Mineralize+RH2PO4TransfLitrAutor(NGL)
      !
      !     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SURFACE RESIDUE FROM
      !     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
      !     ZONES OF SOIL SURFACE
      !
      !     RNetH1PO4DmndLitr=HPO4 mineralization (-ve) or immobilization (+ve) demand
      !     NU=surface layer number
      !     CH1P4U,CH1P4BU=aqueous HPO4 concentrations in non-band, band
      !     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
      !     minimum HPO4 concentration and Km for HPO4 uptake
      !     RH1PO4DmndLitrHeter_col=microbially limited HPO4 demand
      !     BIOA=microbial surface area, OMA=active biomass
      !     TFNG=temp+water stress
      !     FH1PS,FH1PB=fractions of HPO4 in non-band, band
      !     RH1PO4DmndLitrHeter_col=substrate-unlimited HPO4 mineraln-immobiln
      !     VOLWU=water content
      !     H1P4M=HPO4 not available for uptake
      !     AttenfH1PO4Heter=fraction of biological HPO4 demand
      !     RH1PO4imobilLitrHeter=substrate-limited HPO4 minereraln-immobiln
      !     NetPO4Mineralize=total HPO4 net mineraln (-ve) or immobiln (+ve)
      !
      FH1PS = VLPO4
      FH1PB = VLPOB
      RNetH1PO4DmndLitr=0.1_r8*AZMAX1(RNetH2PO4MinPotentLitr-RH2PO4TransfLitrAutor(NGL))

      call LitterSoilNutrientUptake(4,RNetH1PO4DmndLitr, &
        BIOA*OMActAutor(NGL)*GrowthEnvScalAutor(NGL)*HPMX, &
        HPKU,HPMN,AttenfH1PO4Autor(NGL), &
        micflx%RH1PO4UptkLitrBandAutorPrev(NGL),nmics%FracOMActAutor(NGL),micfor,LitrPotential,LitrUptake)
      RH1PO4UptkLitrAutor(NGL)=LitrPotential(1)
      micflx%RH1PO4UptkLitrBandAutor(NGL)=LitrPotential(2)
      RH1PO4TransfLitrAutor(NGL)=SUM(LitrUptake)
      micflx%tRH1PO4MicrbImobilSoil=micflx%tRH1PO4MicrbImobilSoil+LitrUptake(1)
      micflx%tRH1PO4MicrbImobilBand=micflx%tRH1PO4MicrbImobilBand+LitrUptake(2)
      NetPO4Mineralize=NetPO4Mineralize+RH1PO4TransfLitrAutor(NGL)
    ENDIF
  ENDDO
  end associate
  end subroutine BiomNutMinerMobilAutor
!------------------------------------------------------------------------------------------

  subroutine GatherAutotrophRespiration(I,J,N,micfor,micflx,nmicf,nmics)
  implicit none
  integer, intent(in) :: I,J,N
  type(micforctype), intent(in) :: micfor

  type(micfluxtype), intent(inout) :: micflx  
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  character(len=*), parameter :: subname='GatherAutotrophRespiration'
  real(r8) :: RGN2P,FracDenitResp4Maint
  integer  :: NGL
!     begin_execution
  associate(                                             &
    OMActAutor           => nmics%OMActAutor,            & !Active microbial C biomass by autotrophic guild
    RNOxReduxRespAutorLim => nmicf%RNOxReduxRespAutorLim, & !C-equivalent respiration supported by nitrifier denitrification
    RespGrossAutor       => nmicf%RespGrossAutor,        & !Gross respiration C equivalent by autotrophic guild
    Resp4NFixAutor       => nmicf%Resp4NFixAutor,        & !Autotrophic respiration-C cost of N2 fixation (currently set to zero)
    RN2FixAutor          => nmicf%RN2FixAutor,           & !Autotrophic guild N2 fixation flux (set to zero in current respiration gathering)
    RGrowthRespAutor     => micflx%RGrowthRespAutor,     & !Autotrophic gross respiration remaining after maintenance
    RMaintDefcitcitAutor => micflx%RMaintDefcitcitAutor, & !Guild maintenance-C deficit after available gross respiration
    JGniA                => micpar%JGniA,                & !First guild index for each autotrophic functional group
    JGnfA                => micpar%JGnfA,                & !Last guild index for each autotrophic functional group
    RMaintRespAutor      => micflx%RMaintRespAutor       & !Total hourly autotrophic guild maintenance-C demand
  )
!     pH EFFECT ON MAINTENANCE RESPIRATION
!
!     FPH=pH effect on maintenance respiration
!     RMOM=specific maintenance respiration rate
!     TempMaintRHeter=temperature effect on maintenance respiration
!     OMN=microbial N biomass
!
  call PrintInfo('beg '//subname)
  DO NGL=JGniA(N),JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle      
    RGrowthRespAutor(NGL)     = AZMAX1(RespGrossAutor(NGL)-RMaintRespAutor(NGL))
    !Primary growth respiration remains separate from the denitrification
    !pathway; both can supply maintenance on the same reference energy basis.
    call ReserveDenitrifMaintenance(RMaintRespAutor(NGL),RespGrossAutor(NGL), &
      RNOxReduxRespAutorLim(NGL),EO2X,ENOX,FracDenitResp4Maint,RMaintDefcitcitAutor(NGL))
    !
    !     N2 FIXATION: N=(6) AEROBIC, (7) ANAEROBIC
    !     FROM GROWTH RESPIRATION, FIXATION ENERGY REQUIREMENT,
    !     MICROBIAL N REQUIREMENT IN LABILE (1) AND
    !     RESISTANT (2) FRACTIONS
    !
    !     RGN2P=respiration to meet N2 fixation demand
    !     OMC,OMN=microbial nonstructural C,N
    !     rNCOMC=maximum microbial N:C ratio
    !     EN2F=N2 fixation yield per unit nonstructural C
    !     RGrowthRespAutor=growth respiration
    !     Resp4NFixHeter=respiration for N2 fixation
    !     CZ2GS=aqueous N2 concentration
    !     ZFKM=Km for N2 uptake
    !     OMGR*OMC(3,NGL,N,K)=nonstructural C limitation to Resp4NFixHeter
    !     RN2FixHeter=N2 fixation rate
    !
    RN2FixAutor(NGL)    = 0.0_r8
    Resp4NFixAutor(NGL) = 0.0_r8
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine GatherAutotrophRespiration
!------------------------------------------------------------------------------------------

  subroutine AutotrophAnabolicUpdate(micfor,micstt,nmicf,nmicdiag,micflx)

  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_Diag_type), intent(inout) :: nmicdiag    
  type(micfluxtype), intent(inout) :: micflx
  character(len=*), parameter :: subname='AutotrophAnabolicUpdate'

  real(r8) :: CGROMC   !C for microbial biomass growth
  integer :: N,M,NGL,MID,MID3,NE
  real(r8) :: ReserveSupply,MineralTransfer(4),PreviousMineralTransfer
  associate(                                                                &
    DOMuptk4GrothAutor             => nmicf%DOMuptk4GrothAutor,             & !Guild elemental uptake; C source is CO2 or CH4 according to functional group
    NonstX2stBiomAutor             => nmicf%NonstX2stBiomAutor,             & !C/N/P transfer from guild reserves into kinetic and structural biomass
    Resp4NFixAutor                 => nmicf%Resp4NFixAutor,                 & !Autotrophic respiration-C cost of N2 fixation (currently set to zero)
    RespGrossAutor                 => nmicf%RespGrossAutor,                 & !Gross respiration C equivalent by autotrophic guild
    RNOxReduxRespAutorLim          => nmicf%RNOxReduxRespAutorLim,          & !C-equivalent respiration supported by autotrophic nitrite reduction
    RNO3TransfSoilAutor            => nmicf%RNO3TransfSoilAutor,            & !Net NO3-N transfer from nonband soil to microbes; positive immobilization
    RCO2ProdAutor                  => nmicf%RCO2ProdAutor,                  & !CO2-C production by autotrophic guild
    ROMProdCO2Autor                => nmicf%ROMProdCO2Autor,                & !CO2-C release from autotrophic biomass consumed to meet maintenance deficits
    RGrowthCAutor                  => nmicf%RGrowthCAutor,                  & !Net substrate-derived C credited to autotrophic guild reserves
    RCO2XumpAutor                  => nmicf%RCO2XumpAutor,                  & !CO2-C uptake for guild metabolism and biomass, including methanogenic CH4 production
    RH2PO4TransfSoilAutor          => nmicf%RH2PO4TransfSoilAutor,          & !Net H2PO4-P transfer from nonband soil to microbes; positive immobilization
    RNH4TransfBandAutor            => nmicf%RNH4TransfBandAutor,            & !Net NH4-N transfer from fertilizer-band soil to microbes; positive immobilization
    RNO3TransfBandAutor            => nmicf%RNO3TransfBandAutor,            & !Net NO3-N transfer from fertilizer-band soil to microbes; positive immobilization
    RH2PO4TransfBandAutor          => nmicf%RH2PO4TransfBandAutor,          & !Net H2PO4-P transfer from fertilizer-band soil to microbes; positive immobilization
    RkillLitrfal2HumOMAutor        => nmicf%RkillLitrfal2HumOMAutor,        & !Ordinary-mortality C/N/P routed to humified material from autotrophic biomass
    RMaintDefLitrfal2HumOMAutor    => nmicf%RMaintDefLitrfal2HumOMAutor,    & !Starvation-derived C/N/P routed to humified material from autotrophic biomass
    RN2FixAutor                    => nmicf%RN2FixAutor,                    & !Autotrophic guild N2 fixation flux (set to zero in current respiration gathering)
    RKillOMAutor                   => nmicf%RKillOMAutor,                   & !Ordinary mortality C/N/P withdrawal from autotrophic biomass
    RkillRecycOMAutor              => nmicf%RkillRecycOMAutor,              & !Ordinary-mortality C/N/P recycled to reserves from autotrophic biomass
    RMaintDefcitKillOMAutor        => nmicf%RMaintDefcitKillOMAutor,        & !Maintenance-starvation C/N/P withdrawal from autotrophic biomass
    RMaintDefcitRecycOMAutor       => nmicf%RMaintDefcitRecycOMAutor,       & !Starvation recycling: C respired, N/P returned to reserves from autotrophic biomass
    RNH4TransfLitrAutor            => nmicf%RNH4TransfLitrAutor,            & !Net NH4-N transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RNO3TransfLitrAutor            => nmicf%RNO3TransfLitrAutor,            & !Net NO3-N transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RH2PO4TransfLitrAutor          => nmicf%RH2PO4TransfLitrAutor,          & !Net H2PO4-P transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RH1PO4TransfSoilAutor          => nmicf%RH1PO4TransfSoilAutor,          & !Net HPO4-P transfer from nonband soil to microbes; positive immobilization
    RH1PO4TransfBandAutor          => nmicf%RH1PO4TransfBandAutor,          & !Net HPO4-P transfer from fertilizer-band soil to microbes; positive immobilization
    RH1PO4TransfLitrAutor          => nmicf%RH1PO4TransfLitrAutor,          & !Net HPO4-P transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RNH4TransfSoilAutor            => nmicf%RNH4TransfSoilAutor,            & !Net NH4-N transfer from nonband soil to microbes; positive immobilization
    litrm                          => micfor%litrm,                         & !True for the surface litter layer
    ElmAllocmatMicrblitr2POM       => micfor%ElmAllocmatMicrblitr2POM,      & !Partition of humified microbial litter into receiving solid components
    ElmAllocmatMicrblitr2POMU      => micfor%ElmAllocmatMicrblitr2POMU,     & !Underlying-soil partition of humified litter from surface microbes
    SolidOM                        => micstt%SolidOM,                       & !Solid C/N/P pools by substrate component and complex K
    mBiomeAutor                    => micstt%mBiomeAutor,                   & !C/N/P pools indexed by element and flattened guild/biomass compartment
    SOMHumProtein                  => micstt%SOMHumProtein,                 & !C/N/P transferred to the first humus component in underlying soil
    SOMHumCarbohyd                 => micstt%SOMHumCarbohyd,                & !C/N/P transferred to the second humus component in underlying soil
    mid_AutoAmmoniaOxidBacter      => micpar%mid_AutoAmmoniaOxidBacter,     & !Functional-group identifier for ammonia oxidizers
    mid_AutoAMONC10                => micpar%mid_AutoAMONC10          ,     & !Functional-group identifier for nitrite-dependent NC10 methanotrophs
    mid_AutoAMOANME2D              => micpar%mid_AutoAMOANME2D            , & !Functional-group identifier for nitrate-dependent ANME-2d methanotrophs
    mid_AutoNitriteOxidBacter      => micpar%mid_AutoNitriteOxidBacter,     & !Functional-group identifier for nitrite oxidizers
    mid_AutoH2GenoCH4GenArchea     => micpar%mid_AutoH2GenoCH4GenArchea ,   & !Functional-group identifier for hydrogenotrophic methanogens
    JGniA                          => micpar%JGniA,                         & !First guild index for each autotrophic functional group
    JGnfA                          => micpar%JGnfA,                         & !Last guild index for each autotrophic functional group
    NumMicbAFunGrupsPerCmplx       => micpar%NumMicbAFunGrupsPerCmplx,      & !Number of autotrophic functional groups
    icarbhyro                      => micpar%icarbhyro,                     & !Carbohydrate/second solid-component index
    iprotein                       => micpar%iprotein,                      & !Protein/first solid-component index
    k_humus                        => micpar%k_humus,                       & !Humus complex receiving humified microbial C/N/P
    is_activeMicrbFungrpAutor      => micpar%is_activeMicrbFungrpAutor      & !Activation flags for autotrophic functional groups
  )
  call PrintInfo('beg '//subname)
  DO  N=1,NumMicbAFunGrupsPerCmplx
    IF(is_activeMicrbFungrpAutor(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        !Reconcile reserve withdrawals before crediting structural biomass.
        MID3=micpar%get_micb_id(iLbiom_reserve,NGL)
        DO NE=ielmn,ielmp
          ReserveSupply=DOMuptk4GrothAutor(NE,NGL) &
            +SUM(RkillRecycOMAutor(NE,1:2,NGL))+SUM(RMaintDefcitRecycOMAutor(NE,1:2,NGL))
          IF(NE.EQ.ielmn)THEN
            ReserveSupply=ReserveSupply+RN2FixAutor(NGL)
            IF(litrm)ReserveSupply=ReserveSupply+RNH4TransfLitrAutor(NGL)+RNO3TransfLitrAutor(NGL)
            MineralTransfer=[RNH4TransfSoilAutor(NGL),RNH4TransfBandAutor(NGL),RNO3TransfSoilAutor(NGL),RNO3TransfBandAutor(NGL)]
          ELSE
            IF(litrm)ReserveSupply=ReserveSupply+RH2PO4TransfLitrAutor(NGL)+RH1PO4TransfLitrAutor(NGL)
            MineralTransfer=[RH2PO4TransfSoilAutor(NGL),RH2PO4TransfBandAutor(NGL),RH1PO4TransfSoilAutor(NGL),RH1PO4TransfBandAutor(NGL)]
          ENDIF
          PreviousMineralTransfer=SUM(MineralTransfer)
          call LimitReserveNutrientTransfers(mBiomeAutor(NE,MID3),ReserveSupply, &
            NonstX2stBiomAutor(NE,1:2,NGL),MineralTransfer)
          IF(NE.EQ.ielmn)THEN
            RNH4TransfSoilAutor(NGL)=MineralTransfer(1)
            RNH4TransfBandAutor(NGL)=MineralTransfer(2)
            RNO3TransfSoilAutor(NGL)=MineralTransfer(3)
            RNO3TransfBandAutor(NGL)=MineralTransfer(4)
            micflx%NetNH4Mineralize=micflx%NetNH4Mineralize+SUM(MineralTransfer)-PreviousMineralTransfer
          ELSE
            RH2PO4TransfSoilAutor(NGL)=MineralTransfer(1)
            RH2PO4TransfBandAutor(NGL)=MineralTransfer(2)
            RH1PO4TransfSoilAutor(NGL)=MineralTransfer(3)
            RH1PO4TransfBandAutor(NGL)=MineralTransfer(4)
            micflx%NetPO4Mineralize=micflx%NetPO4Mineralize+SUM(MineralTransfer)-PreviousMineralTransfer
          ENDIF
        ENDDO

        DO  M=1,2
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            mBiomeAutor(NE,MID)=mBiomeAutor(NE,MID)+NonstX2stBiomAutor(NE,M,NGL)-RKillOMAutor(NE,M,NGL)-RMaintDefcitKillOMAutor(NE,M,NGL)
          ENDDO

!     HUMIFICATION PRODUCTS
!
!     ElmAllocmatMicrblitr2POM=fractions allocated to humic vs fulvic humus
!     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P LitrFall to humus
!     RHMMC,RHMMN,RHMMC=transfer of senesence LitrFall C,N,P to humus
!
          IF(.not.litrm)THEN
            DO NE=1,NumPlantChemElms
              SolidOM(NE,iprotein,k_humus)=SolidOM(NE,iprotein,k_humus)+ElmAllocmatMicrblitr2POM(1) &
                *(RkillLitrfal2HumOMAutor(NE,M,NGL)+RMaintDefLitrfal2HumOMAutor(NE,M,NGL))
              SolidOM(NE,icarbhyro,k_humus)=SolidOM(NE,icarbhyro,k_humus)+ElmAllocmatMicrblitr2POM(2)&
                *(RkillLitrfal2HumOMAutor(NE,M,NGL)+RMaintDefLitrfal2HumOMAutor(NE,M,NGL))
            ENDDO
          ELSE
            DO NE=1,NumPlantChemElms
              SOMHumProtein(NE)=SOMHumProtein(NE)+ElmAllocmatMicrblitr2POMU(1) &
                *(RkillLitrfal2HumOMAutor(NE,M,NGL)+RMaintDefLitrfal2HumOMAutor(NE,M,NGL))
              SOMHumCarbohyd(NE)=SOMHumCarbohyd(NE)+ElmAllocmatMicrblitr2POMU(2) &
                *(RkillLitrfal2HumOMAutor(NE,M,NGL)+RMaintDefLitrfal2HumOMAutor(NE,M,NGL))
            ENDDO
          ENDIF
        ENDDO
        !
        !     INPUTS TO NONSTRUCTURAL POOLS
        !
        !     CGOMC=total DOC+acetate uptake
        !     RespGrossHeter=total respiration
        !     RNOxReduxRespDenitLim=respiration for denitrifcation
        !     Resp4NFixHeter=respiration for N2 fixation
        !     RCO2X=total CO2 emission
        !     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
        !     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling
        !     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
        !     CGOMN,CGOMP=DON, DOP uptake
        !
        CGROMC             = DOMuptk4GrothAutor(ielmc,NGL)-RespGrossAutor(NGL)-RNOxReduxRespAutorLim(NGL)-Resp4NFixAutor(NGL)
        RGrowthCAutor(NGL) = CGROMC
        RCO2ProdAutor(NGL) = RCO2ProdAutor(NGL)+Resp4NFixAutor(NGL)
        MID3               = micpar%get_micb_id(iLbiom_reserve,NGL)

        if(N.eq.mid_AutoAMONC10)then
          !environmental CO2 is assimilated for C biomass
          RCO2XumpAutor(NGL)= CGROMC
        elseif(N.EQ.mid_AutoAMOANME2D)then
          !recyle some CO2
          RCO2ProdAutor(NGL)=RCO2ProdAutor(NGL)-CGROMC
        elseif(N.eq.mid_AutoH2GenoCH4GenArchea)then
          !CO2 supplies carbon for both methane production and biomass growth.
          RCO2XumpAutor(NGL)=nmicf%RCH4ProdAutor(NGL)+CGROMC
          !Additional H2 consumed in biomass synthesis: CO2 + 2H2 -> CH2O + H2O.
          nmicdiag%RH2UptkAutor=nmicdiag%RH2UptkAutor+0.333_r8*CGROMC
        elseif(N.eq.mid_AutoNitriteOxidBacter .or. N.eq.mid_AutoAmmoniaOxidBacter)then
          !For NH3/NO2 oxidizers, CO2 supports both respiration and biomass growth.
          RCO2XumpAutor(NGL)=DOMuptk4GrothAutor(ielmc,NGL)
        endif

        DO M=1,2
          DO NE=1,NumPlantChemElms
            mBiomeAutor(NE,MID3)=mBiomeAutor(NE,MID3)-NonstX2stBiomAutor(NE,M,NGL)+RkillRecycOMAutor(NE,M,NGL)
          ENDDO

          !TODO: Resolve starvation respiration with pathway-specific redox budgets.
          !This CO2 release currently has no matching electron-acceptor debit.
          !Aerobic respiration needs available O2 and consistent demand/uptake;
          !anaerobic groups need their own substrates and reduced products,
          !not an O2 charge. Preserve unrespired C in organic pools and avoid
          !double-counting maintenance already supported by catabolism.
          !See the corresponding starvation update in MicBGCFGMod.
          !C is respired as CO2 while N and P are recycled.
          mBiomeAutor(ielmn,MID3) = mBiomeAutor(ielmn,MID3)+RMaintDefcitRecycOMAutor(ielmn,M,NGL)
          mBiomeAutor(ielmp,MID3) = mBiomeAutor(ielmp,MID3)+RMaintDefcitRecycOMAutor(ielmp,M,NGL)
          RCO2ProdAutor(NGL)      = RCO2ProdAutor(NGL)+RMaintDefcitRecycOMAutor(ielmc,M,NGL)
          ROMProdCO2Autor(NGL)    =  ROMProdCO2Autor(NGL)+RMaintDefcitRecycOMAutor(ielmc,M,NGL)
        ENDDO
        
        mBiomeAutor(ielmc,MID3)=mBiomeAutor(ielmc,MID3)+CGROMC
        mBiomeAutor(ielmn,MID3)=mBiomeAutor(ielmn,MID3)+DOMuptk4GrothAutor(ielmn,NGL) &
          +RNH4TransfSoilAutor(NGL)+RNH4TransfBandAutor(NGL)+RNO3TransfSoilAutor(NGL) &
          +RNO3TransfBandAutor(NGL)+RN2FixAutor(NGL)
        
        mBiomeAutor(ielmp,MID3)=mBiomeAutor(ielmp,MID3)+DOMuptk4GrothAutor(ielmp,NGL) &
          +RH2PO4TransfSoilAutor(NGL)+RH2PO4TransfBandAutor(NGL)+RH1PO4TransfSoilAutor(NGL)&
          +RH1PO4TransfBandAutor(NGL)
        IF(litrm)THEN
          mBiomeAutor(ielmn,MID3)=mBiomeAutor(ielmn,MID3)+RNH4TransfLitrAutor(NGL)+RNO3TransfLitrAutor(NGL)
          mBiomeAutor(ielmp,MID3)=mBiomeAutor(ielmp,MID3)+RH2PO4TransfLitrAutor(NGL)+RH1PO4TransfLitrAutor(NGL)
        ENDIF
      enddo
    ENDIF
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine AutotrophAnabolicUpdate

end module MicAutoCPLXMod
