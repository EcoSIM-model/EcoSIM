module AerobicBacteriaMod
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: PrintInfo
  use minimathmod
  use ElmIDMod
  use TracerIDMod
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod,    only: AerobicHeterO2Uptake, CalcRespMaintHeter, &
                                  StageFuncGuild

  implicit none

  private
  public :: AerobicHeteroBactCatabolism

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains

!------------------------------------------------------------------------------------------

  subroutine AerobicHeteroBactCatabolism(I,J,N,K,RMOMK,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
  !
  !Description
  !catabolism of aerobic heterotrophs
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: N,K
  real(r8), intent(in) :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout) :: nmics
  type(OMCplx_State_type),intent(inout):: ncplxs
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='AerobicHeteroBactCatabolism'
  real(r8) :: WatStressMicb  !moisture sensivity of microbial activity
  integer :: NGL
  real(r8) :: RGOMP  !total DOC/acetate C uptake for potential respiraiton
  real(r8) :: EO2Q   !respiraiton efficiency, i.e. fraction of 1 gC (uptake) used for respiraiton
  real(r8) :: OXKX,FOXYX
  real(r8) :: FSBSTC,FSBSTA
  real(r8) :: RGOCY,RGOCZ,RGOAZ
  real(r8) :: RGOCX,RGOAX

!     begin_execution
  associate(                                                               &
    OMActHeter             => nmics%OMActHeter,                            & !Active microbial C biomass by heterotrophic guild and complex K
    FBiomStoiScalarHeter   => nmics%FBiomStoiScalarHeter,                  & !Combined N/P stoichiometric multiplier on guild metabolic capacity [-]
    GrowthEnvScalHeter     => nmics%GrowthEnvScalHeter,                    & !Temperature and water-potential multiplier on heterotrophic growth [-]
    WSensGroHeter          => nmics%WSensGroHeter,                         & !Guild soil-water-potential multiplier on heterotrophic growth [-]; not referenced here
    TSensGroHeter          => nmics%TSensGroHeter,                         & !Guild temperature multiplier on heterotrophic growth [-]; not referenced here
    TempMaintRHeter        => nmics%TempMaintRHeter,                       & !Guild temperature multiplier on heterotrophic maintenance [-]; not referenced here
    FracOMActHeter         => nmics%FracOMActHeter,                        & !Guild/complex fraction of total active microbial C in the layer [-]
    FracNO2ReduxHeter      => nmics%FracNO2ReduxHeter,                     & !Guild share of active nitrite-consuming biomass for competition [-]; not referenced here
    FracHeterBiomOfActK    => nmics%FracHeterBiomOfActK,                   & !Guild fraction of active heterotrophic biomass in complex K [-]; not referenced here
    FSBSTHeter             => nmicdiag%FSBSTHeter,                         & !Guild substrate-response factor; larger values mean less limitation [-]
    RO2Dmnd4RespHeter      => nmicf%RO2Dmnd4RespHeter,                     & !Potential O2 demand supporting heterotrophic gross respiration
    OxyLimterHeter         => nmics%OxyLimterHeter,                        & !Actual/potential O2 uptake ratio; 1 means no O2 restriction [-]
    RO2DmndHeter           => nmicf%RO2DmndHeter,                          & !Total guild O2 demand before O2 limitation
    ROQC4HeterMicrobAct    => nmicf%ROQC4HeterMicrobAct,                   & !Guild activity proxy for substrate hydrolysis, with DOC concentration unconstrained
    ECHZHeter              => nmicf%ECHZHeter,                             & !Guild respiration fraction used to convert growth respiration to C uptake [-]
    RCO2ProdHeter          => nmicf%RCO2ProdHeter,                         & !CO2-C production by heterotrophic guild and complex
    RCH4ProdHeter          => nmicf%RCH4ProdHeter,                         & !CH4-C production by heterotrophic guild and complex
    RO2Uptk4RespHeter      => nmicf%RO2Uptk4RespHeter,                     & !Realized O2 uptake attributed to heterotrophic gross respiration
    RH2ProdHeter           => nmicf%RH2ProdHeter,                          & !H2 production by heterotrophic guild and complex
    RGOCP                  => nmicf%RGOCP,                                 & !DOC-supported potential respiration before O2 limitation
    RGOAP                  => nmicf%RGOAP,                                 & !Acetate-supported potential respiration before O2 limitation
    FOQC                   => nmicf%FOQC,                                  & !Guild share of DOC demand used to allocate the donor pool [-]
    FOQA                   => nmicf%FOQA,                                  & !Guild share of acetate demand used to allocate the donor pool [-]
    FGOCP                  => nmicf%FGOCP,                                 & !DOC-supported fraction of total primary guild respiration [-]
    RespGrossHeter         => nmicf%RespGrossHeter,                        & !Gross respiration C equivalent from the primary heterotrophic pathway
    RAcetateProdHeter      => nmicf%RAcetateProdHeter,                     & !Acetate-C production by heterotrophic guild and complex
    FGOAP                  => nmicf%FGOAP,                                 & !Acetate-supported fraction of total primary guild respiration [-]
    ZEROS                  => micfor%ZEROS,                                & !Small mass or flux threshold used by the routine
    RO2EcoDmndPrev         => micfor%RO2EcoDmndPrev,                       & !Previous-hour ecosystem O2 demand; competition denominator
    RDOMEcoDmndPrev        => micfor%RDOMEcoDmndPrev,                      & !Previous-hour ecosystem DOC demand in each complex; competition denominator; not referenced here
    RAcetateEcoDmndPrev    => micfor%RAcetateEcoDmndPrev,                  & !Previous-hour ecosystem acetate demand in each complex; competition denominator; not referenced here
    RDOCUptkHeterPrev      => micfor%RDOCUptkHeterPrev,                    & !Previous-hour guild DOC uptake/demand used for competition; not referenced here
    RAcetateUptkHeterPrev  => micfor%RAcetateUptkHeterPrev,                & !Previous-hour guild acetate uptake/demand used for competition; not referenced here
    PSISoilMatricP         => micfor%PSISoilMatricP,                       & !Soil matric water potential controlling microbial water stress; not referenced here
    DOM                    => micstt%DOM,                                  & !Dissolved organic pools by species (DOC, DON, DOP, acetate) and complex K
    mid_HeterAerobBacter   => micpar%mid_HeterAerobBacter,                 & !Functional-group identifier for aerobic heterotrophic bacteria
    mid_Facult_DenitBacter => micpar%mid_Facult_DenitBacter,               & !Functional-group identifier for facultative denitrifying bacteria
    mid_HeterAerobN2Fixer  => micpar%mid_HeterAerobN2Fixer,                & !Functional-group identifier for aerobic heterotrophic N2 fixers
    TOMEK                  => nmicdiag%TOMEK,                              & !Total active heterotrophic C/N/P in each substrate complex K; not referenced here
    TSensGrowth            => nmicdiag%TSensGrowth,                        & !Layer temperature response for microbial growth [-]; not referenced here
    TSensMaintR            => nmicdiag%TSensMaintR,                        & !Layer temperature response for microbial maintenance [-]; not referenced here
    TotBiomNO2Consumers    => nmicdiag%TotBiomNO2Consumers,                & !Active C summed over staged nitrite-competing microbial groups; not referenced here
    RO2DmndHetert          => micflx%RO2DmndHetert,                        & !Guild O2 demand retained for substrate-competition accounting
    RDOCUptkHeter          => micflx%RDOCUptkHeter,                        & !Potential DOC uptake used in guild substrate-competition accounting
    RAcetateUptkHeter      => micflx%RAcetateUptkHeter,                    & !Potential acetate uptake used in guild substrate-competition accounting
    RO2DmndHetertPrev      => micflx%RO2DmndHetertPrev,                    & !Previous-hour guild O2 demand used for competition
    tRespGrossHeterUlm     => naqfdiag%tRespGrossHeterUlm,                 & !Total gross heterotrophic respiration before O2 limitation
    TotActMicrobiom        => nmicdiag%TotActMicrobiom,                    & !Layer total active microbial C across heterotrophs and autotrophs
    tRGOXP                 => micflx%tRGOXP,                               & !Accumulated substrate-pool-limited potential respiration capacity
    tRGOZP                 => micflx%tRGOZP,                               & !Accumulated kinetic potential respiration capacity before donor/O2 limits
    FOCA                   => ncplxs%FOCA,                                 & !DOC fraction of DOC plus acetate in each substrate complex [-]
    FOAA                   => ncplxs%FOAA,                                 & !Acetate-C fraction of DOC plus acetate in each substrate complex [-]
    CDOM                   => ncplxs%CDOM                                  & !Dissolved substrate concentrations by DOM species and complex K
  )
  call PrintInfo('beg '//subname)
  !     ENERGY YIELDS OF O2 REDOX REACTIONS
  !     E* = growth respiration efficiency calculated in PARAMETERS
  !
  !loop over all guilds of a given functional group
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    IF(OMActHeter(NGL,K).LE..0_r8)cycle

    !prepare trait parameters
    call StageFuncGuild(N,NGL,K,TotActMicrobiom,FOQC(NGL,K),FOQA(NGL,K),micfor,naqfdiag,nmicdiag,nmics)

    call CalcRespMaintHeter(NGL,K,RMOMK,micfor,micstt,micflx,nmicf,nmics)

    OXKX  = OXKM
    IF(RO2EcoDmndPrev.GT.ZEROS)THEN
      FOXYX=AMAX1(FMN,RO2DmndHetertPrev(NGL,K)/RO2EcoDmndPrev)
    ELSE
      FOXYX=AMAX1(FMN,FracOMActHeter(NGL,K))
    ENDIF
    naqfdiag%TFOXYX = naqfdiag%TFOXYX+FOXYX

    ! N=OBLIGATE AEROBIC bacteria,
    IF(N.EQ.mid_HeterAerobBacter)THEN
      EO2Q=EO2X
    ! FACULTATIVE ANAEROBES,
    ELSEIF(N.EQ.mid_Facult_DenitBacter)THEN
      EO2Q=EO2D
    !aerobic N2 FIXERS
    ELSEIF(N.EQ.mid_HeterAerobN2Fixer)THEN
      EO2Q=ENFX
    ENDIF
    !
    ! O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
    ! 'RGO*Z'FROM SPECIFIC RESPIRATION RATE, ACTIVE BIOMASS, DOC OR
    ! ACETATE CONCENTRATION,MICROBIAL C:N:P FACTOR, AND TEMPERATURE
    ! FOLLOWED BY POTENTIAL RESPIRATION RATES 'RGO*P' WITH UNLIMITED
    ! SUBSTRATE USED FOR MICROBIAL COMPETITION FACTOR

    ! COQC,COQA=DOC,DOA concentration, FOCA,FOAA=DOC,DOA vs DOC+DOA
    ! FBiomStoiScalarHeter=N,P limitation,VMXO=specific respiration rate
    ! WatStressMicb=water stress effect, OMA=active biomass
    ! TSensGrowth=temp stress effect,FOQC,FOQA=OQC,OQA limitation
    ! RGOMP=O2-unlimited respiration of DOC+DOA
    ! RGOCP,RGOAP,RGOMP=O2-unlimited respiration of DOC, DOA, DOC+DOA
    !
    FSBSTC            = CDOM(idom_doc,K)/(CDOM(idom_doc,K)+OQKM)
    FSBSTA            = CDOM(idom_acetate,K)/(CDOM(idom_acetate,K)+OQKA)
    FSBSTHeter(NGL,K) = FOCA(K)*FSBSTC+FOAA(K)*FSBSTA

    RGOCY  = AZMAX1(FBiomStoiScalarHeter(NGL,K)*OMActHeter(NGL,K))*VMXO*GrowthEnvScalHeter(NGL,K)
    RGOCZ  = RGOCY*FSBSTC*FOCA(K) !MM uptake of DOC
    RGOAZ  = RGOCY*FSBSTA*FOAA(K) !MM uptake of acetate

    !obtain kinetically unlimited DOM/acetate uptake
    RGOCX = AZMAX1(DOM(idom_doc,K)*FOQC(NGL,K)*EO2Q)             !DOC respiration
    RGOAX = AZMAX1(DOM(idom_acetate,K)*FOQA(NGL,K)*EO2A)         !acetate respiraiton

    !obtain the final uptake
    RGOCP(NGL,K) = AMIN1(RGOCX,RGOCZ)      !DOC respiraiton
    RGOAP(NGL,K) = AMIN1(RGOAX,RGOAZ)      !acetate respiraiton
    RGOMP        = RGOCP(NGL,K)+RGOAP(NGL,K)      !total C respiration before O2 limitation

    tRGOXP = tRGOXP+RGOCX+RGOAX
    tRGOZP = tRGOZP+RGOCZ+RGOAZ            !potential C oxidation without C and O2 limitation
    IF(RGOMP.GT.ZEROS)THEN
      FGOCP(NGL,K) = RGOCP(NGL,K)/RGOMP
      FGOAP(NGL,K) = RGOAP(NGL,K)/RGOMP
    ELSE
      FGOCP(NGL,K) = 1.0_r8
      FGOAP(NGL,K) = 0.0_r8
    ENDIF
    !
    ! ENERGY YIELD AND O2 DEMAND FROM DOC AND ACETATE OXIDATION
    ! BY HETEROTROPHIC AEROBES

    ! ECHZHeter=growth respiration yield, averaged over acetate and DOC/glucose
    ! RO2Dmnd4RespHeter,RO2DmndHeter,RO2DmndHetert=O2 demand from DOC,DOA oxidation
    ! RDOCUptkHeter,RAcetateUptkHeter=DOC,DOA demand from DOC,DOA oxidation
    ! ROQC4HeterMicrobAct=microbial respiration used to represent microbial activity
    ! CH2O+O2 -> CO2 + H2O, (32/12.=2.667)
    !Respiration-weighted conversion: total uptake is the sum of each
    !donor's respiration divided by its own respiration fraction.
    ECHZHeter(NGL,K)         = 1._r8/(FGOCP(NGL,K)/EO2Q+FGOAP(NGL,K)/EO2A)
    RO2Dmnd4RespHeter(NGL,K) = 2.667_r8*RGOMP                 !O2 demand
    RO2DmndHeter(NGL,K)      = RO2Dmnd4RespHeter(NGL,K)

    !make a copy for flux limiter
    RO2DmndHetert(NGL,K)       = RO2DmndHeter(NGL,K)
    RDOCUptkHeter(NGL,K)       = RGOCZ               !potential DOC (O2-unlimited) uptake flux
    RAcetateUptkHeter(NGL,K)   = RGOAZ               !potential acetate (O2-unlimited) uptake flux

    call AerobicHeterO2Uptake(I,J,NGL,N,K,FOXYX,OXKX,micfor,micstt,nmicf,nmics,micflx)

    ROQC4HeterMicrobAct(NGL,K) = RGOCY*OxyLimterHeter(NGL,K)  !C demand for oxidation
    RespGrossHeter(NGL,K)      = RGOMP*OxyLimterHeter(NGL,K)  !actual respiration O2-limited
    RCO2ProdHeter(NGL,K)       = RespGrossHeter(NGL,K)
    RAcetateProdHeter(NGL,K)     = 0.0_r8
    RCH4ProdHeter(NGL,K)       = 0.0_r8
    RO2Uptk4RespHeter(NGL,K)   = RO2Dmnd4RespHeter(NGL,K)*OxyLimterHeter(NGL,K)
    RH2ProdHeter(NGL,K)        = 0.0_r8
    tRespGrossHeterUlm         = tRespGrossHeterUlm+RGOMP
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine AerobicHeteroBactCatabolism

end module AerobicBacteriaMod
