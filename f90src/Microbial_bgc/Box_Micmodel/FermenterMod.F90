module FermenterMod
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use minimathmod
  use TracerIDMod
  use EcosimConst
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod,    only: CalcRespMaintHeter, StageFuncGuild

  implicit none

  private
  public :: AcetogFermentCatabolism

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains

!------------------------------------------------------------------------------------------

  subroutine AcetogFermentCatabolism(N,K,RMOMK,micfor,micstt,naqfdiag,ncplxs,nmicf,nmics,micflx,nmicdiag)
  !
  !Description:
  !Fermentation and acetogenic N2 fixers
  !(CH2O)6 +2H2O-> 2CO2 + 2(CH2O)2 + 4H2, mole based
  !(CH2O)6 -> 2CO2 + 2/3 (CH2O)2 + 8/(72)H2, mass based
  !fermenters only take up DOC/glucose
  !it can be fermenters or anaerobic N2 fixers
  implicit none
  integer, intent(in) :: N,K
  real(r8), intent(in) :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(OMCplx_State_type), intent(inout) :: ncplxs
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout):: nmics
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  real(r8) :: WatStressMicb
  real(r8) :: RGOMP   !potential respiration for metabolism
  integer  :: NGL
  real(r8) :: GH2X,GH2F
  real(r8) :: GOAX,GOAF
  real(r8) :: GHAX
  REAL(R8) :: oxyi
  real(r8) :: RGOFX,RGOFY,RGOFZ

  real(r8), parameter :: GlucoseC=72._r8  !glucose has 72 gC/mol

  ! begin_execution
  associate(                                            &
    FBiomStoiScalarHeter => nmics%FBiomStoiScalarHeter, & !Combined N/P stoichiometric multiplier on guild metabolic capacity [-]
    OMActHeter           => nmics%OMActHeter,           & !Active microbial C biomass by heterotrophic guild and complex K
    FSBSTHeter           => nmicdiag%FSBSTHeter,        & !Guild substrate-response factor; larger values mean less limitation [-]
    GrowthEnvScalHeter   => nmics%GrowthEnvScalHeter,   & !Temperature and water-potential multiplier on heterotrophic growth [-]
    RO2Dmnd4RespHeter    => nmicf%RO2Dmnd4RespHeter,    & !Potential O2 demand supporting heterotrophic gross respiration; zero for this anaerobic pathway
    RO2DmndHeter         => nmicf%RO2DmndHeter,         & !Total guild O2 demand before O2 limitation; zero for this anaerobic pathway
    ECHZHeter            => nmicf%ECHZHeter,            & !Guild respiration fraction used to convert growth respiration to C uptake [-]
    FOQC                 => nmicf%FOQC,                 & !Guild share of DOC demand used to allocate the donor pool [-]
    FOQA                 => nmicf%FOQA,                 & !Guild share of acetate demand used to allocate the donor pool [-]
    RCH4ProdHeter        => nmicf%RCH4ProdHeter,        & !CH4-C production by heterotrophic guild and complex
    RO2Uptk4RespHeter    => nmicf%RO2Uptk4RespHeter,    & !Realized O2 uptake attributed to heterotrophic gross respiration; zero for this anaerobic pathway
    RCO2ProdHeter        => nmicf%RCO2ProdHeter,        & !CO2-C production by heterotrophic guild and complex
    RespGrossHeter       => nmicf%RespGrossHeter,       & !Gross respiration C equivalent from the primary heterotrophic pathway
    RH2ProdHeter         => nmicf%RH2ProdHeter,         & !H2 production by heterotrophic guild and complex
    ROQC4HeterMicrobAct  => nmicf%ROQC4HeterMicrobAct,  & !Guild activity proxy for substrate hydrolysis, with DOC concentration unconstrained
    RAcetateProdHeter    => nmicf%RAcetateProdHeter,    & !Acetate-C production by heterotrophic guild and complex
    TotActMicrobiom      => nmicdiag%TotActMicrobiom,   & !Layer total active microbial C across heterotrophs and autotrophs
    FGOCP                => nmicf%FGOCP,                & !DOC-supported fraction of total primary guild respiration [-]
    FGOAP                => nmicf%FGOAP,                & !Acetate-supported fraction of total primary guild respiration [-]
    TKS                  => micfor%TKS,                 & !Layer absolute temperature [K]
    PSISoilMatricP       => micfor%PSISoilMatricP,      & !Soil matric water potential controlling microbial water stress; not referenced here
    ZERO                 => micfor%ZERO,                & !Small dimensionless or concentration threshold used by the routine
    DOM                  => micstt%DOM,                 & !Dissolved organic pools by species (DOC, DON, DOP, acetate) and complex K
    CH2GS                => micstt%CH2GS,               & !Dissolved H2 concentration used in energy-yield and saturation calculations
    COXYS                => micstt%COXYS,               & !Dissolved O2 concentration
    RO2DmndHetert        => micflx%RO2DmndHetert,       & !Guild O2 demand retained for substrate-competition accounting; zero for this anaerobic pathway
    RDOCUptkHeter        => micflx%RDOCUptkHeter,       & !Potential DOC uptake used in guild substrate-competition accounting
    RAcetateUptkHeter    => micflx%RAcetateUptkHeter,   & !Potential acetate uptake used in guild substrate-competition accounting
    mid_fermentor        => micpar%mid_fermentor,       & !Functional-group identifier for fermenters
    CDOM                 => ncplxs%CDOM                 & !Dissolved substrate concentrations by DOM species and complex K
  )

  !
  !     ENERGY YIELD FROM FERMENTATION DEPENDS ON H2 AND
  !     ACETATE CONCENTRATION
  !
  !     GH2F=energy yield of acetotrophic methanogenesis per g C
  !     GHAX=H2 effect on energy yield of fermentation
  !     GOAX=acetate effect on energy yield of fermentation
  !     ECHZHeter=growth respiration efficiency of fermentation

  !loop over all guilds of a given functional group
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    IF(OMActHeter(NGL,K).LE.0.0_r8)cycle
    !prepare trait parameters
    call StageFuncGuild(N,NGL,K,TotActMicrobiom,FOQC(NGL,K),FOQA(NGL,K),micfor,naqfdiag,nmicdiag,nmics)

    call CalcRespMaintHeter(NGL,K,RMOMK,micfor,micstt,micflx,nmicf,nmics)

    GH2X = RGASC*1.E-3_r8*TKS*LOG((AMAX1(1.0E-05_r8,CH2GS)/H2KI)**4)
    GH2F = GH2X/GlucoseC    !
    GOAX = RGASC*1.E-3_r8*TKS*LOG((AMAX1(ZERO,CDOM(idom_acetate,K))/OAKI)**2)
    GOAF = GOAX/GlucoseC
    GHAX = GH2F+GOAF
    IF(N.EQ.mid_fermentor)THEN
      ECHZHeter(NGL,K)=AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GCHX-GHAX))/EOMF)))
    ELSE
        !dizotrophs, i.e. N2 fixers
      ECHZHeter(NGL,K)=AMAX1(ENFX,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GCHX-GHAX))/EOMN)))
    ENDIF
    !
    !     RESPIRATION RATES BY HETEROTROPHIC ANAEROBES 'RGOMP' FROM
    !     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
    !     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL
    !     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
    !     MICROBIAL COMPETITION FACTOR
    !
    !     OXYI=O2 inhibition of fermentation
    !     FBiomStoiScalarHeter=N,P limitation on respiration
    !     VMXF=maximum respiration rate by fermenters
    !     WatStressMicb=water stress effect on respiration
    !     OMA=active fermenter biomass
    !     TSensGrowth=temp stress effect, FOQC=OQC limitation
    !     RFOMP=O2-unlimited respiration of DOC
    !     ROQC4HeterMicrobAct=microbial respiration used to represent microbial activity
    !
    OXYI  = 1.0_r8-1.0_r8/(1.0_r8+EXP(1.0_r8*AMAX1(-COXYS+2.5_r8,-50._r8)))
    FSBSTHeter(NGL,K) = CDOM(idom_doc,K)/(CDOM(idom_doc,K)+OQKM)*OXYI
    RGOFY             = AZMAX1(FBiomStoiScalarHeter(NGL,K)*OMActHeter(NGL,K))*VMXF*GrowthEnvScalHeter(NGL,K)
    RGOFZ             = RGOFY*FSBSTHeter(NGL,K)
    RGOFX             = AZMAX1(DOM(idom_doc,K)*FOQC(NGL,K)*ECHZHeter(NGL,K))

    !potential respiration to expense
    RGOMP                      = AMIN1(RGOFX,RGOFZ)
    FGOCP(NGL,K)               = 1.0_r8
    FGOAP(NGL,K)               = 0.0_r8
    RO2Dmnd4RespHeter(NGL,K)   = 0.0_r8   !demand no oxygen
    RO2DmndHeter(NGL,K)        = 0.0_r8
    RO2DmndHetert(NGL,K)       = 0.0_r8
    RDOCUptkHeter(NGL,K)       = RGOFZ    !potential DOC (unlimited) uptake flux
    RAcetateUptkHeter(NGL,K)   = 0.0_r8
    ROQC4HeterMicrobAct(NGL,K) = RGOFY*OXYI    !DOC/temperature-unlimited fermentation rate for microbial activity calculation
    naqfdiag%tCResp4H2Prod     = naqfdiag%tCResp4H2Prod+RGOMP

    !fermentation  (CH2O)6 -> 2CO2 + 2(CH2O)2
    RespGrossHeter(NGL,K)    = RGOMP
    RAcetateProdHeter(NGL,K)   = 0.667_r8*RespGrossHeter(NGL,K)
    RCO2ProdHeter(NGL,K)     = AZMAX1(RespGrossHeter(NGL,K)-RAcetateProdHeter(NGL,K))
    RCH4ProdHeter(NGL,K)     = 0.0_r8
    RO2Uptk4RespHeter(NGL,K) = RO2Dmnd4RespHeter(NGL,K)
    RH2ProdHeter(NGL,K)      = 0.111_r8*RespGrossHeter(NGL,K)

  ENDDO
!
  end associate
  end subroutine AcetogFermentCatabolism
end module FermenterMod
