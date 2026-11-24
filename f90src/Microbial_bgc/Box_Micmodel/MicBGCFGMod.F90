module MicBGCMod
!!
! DESCRIPTION:
! codes to do soil biological transfOMBioResduations
!
! USES:
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use abortutils,           only: endrun,   destroy
  use EcoSIMCtrlMod,        only: etimer
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: DebugPrint,PrintInfo
  use minimathmod            
  use ElmIDMod
  use TracerIDMod
  use MicAutoCPLXMod
  use EcosimConst
  use EcoSIMSolverPar
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod
  use AerobicBacteriaMod,   only: AerobicHeteroBactCatabolism
  use CyanoBacteriaMod,     only: CyanoBacteriaCatabolism
  use DenitrifierMod,       only: HeteroDenitrificCatabolism
  use FermenterMod,         only: AcetogFermentCatabolism
  use FungiMod,             only: AerobicFungiCatabolism
  use MethanogenMod,        only: AcetoMethanogenCatabolism
  implicit none

  private

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  integer :: jcplx,NumMicbHFunGrupsPerCmplx,NumMicbAFunGrupsPerCmplx,jsken,ndbiomcp,nlbiomcp
  integer, pointer :: JGniA(:)
  integer, pointer :: JGnfA(:)
  integer, pointer :: JGniH(:)
  integer, pointer :: JGnfH(:)
!
  public :: initNitro1Layer, SoilBGCOneLayer

  contains

!------------------------------------------------------------------------------------------

  subroutine initNitro1Layer
!
! DESCRIPTION:
! initialize single layer microibal bgc model
  implicit none

  jcplx =micpar%jcplx
  NumMicbHFunGrupsPerCmplx  =micpar%NumMicbHFunGrupsPerCmplx
  NumMicbAFunGrupsPerCmplx  =micpar%NumMicbAFunGrupsPerCmplx  
  jsken =micpar%jsken
  ndbiomcp = micpar%ndbiomcp
  nlbiomcp = micpar%nlbiomcp

  JGniA => micpar%JGniA
  JGnfA => micpar%JGnfA
  JGniH => micpar%JGniH
  JGnfH => micpar%JGnfH

  call initNitroPars

  end subroutine initNitro1Layer

!------------------------------------------------------------------------------------------

  subroutine SoilBGCOneLayer(I,J,micfor,micstt,micflx,naqfdiag,nmicdiag)
  implicit none
  integer, intent(in) :: I,J
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type), intent(out) :: naqfdiag
  type(Microbe_Diag_type),intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='SoilBGCOneLayer'
!local variables
  integer :: LL,K,KL,NGL
  integer :: M,N

  type(Microbe_State_type) :: nmics
  type(Microbe_Flux_type) :: nmicf
  type(OMCplx_Flux_type) :: ncplxf
  type(OMCplx_State_type) :: ncplxs
  real(r8) :: totOMbeg(1:NumPlantChemElms)
  real(r8) :: totOMend(1:NumPlantChemElms)
  real(r8) :: domsb(1:NumPlantChemElms),somb(1:NumPlantChemElms)
  real(r8) :: sorbomb(1:NumPlantChemElms),biomheterb(1:NumPlantChemElms)
  real(r8) :: biomautob(1:NumPlantChemElms),biomresb(1:NumPlantChemElms)
  real(r8) :: domse(1:NumPlantChemElms),some(1:NumPlantChemElms)
  real(r8) :: sorbome(1:NumPlantChemElms),biomhetere(1:NumPlantChemElms)
  real(r8) :: biomautoe(1:NumPlantChemElms),biomrese(1:NumPlantChemElms)

! begin_execution
  call PrintInfo('beg '//subname)
  call nmicf%Init(jcplx)
  call nmics%Init(jcplx)
  call ncplxf%Init()
  call ncplxs%Init()
  call naqfdiag%ZeroOut()
  call micflx%ZeroOut()
  call nmicdiag%ZeroOut()

  call StageBGCEnvironCondition(I,J,micfor,KL,micstt,naqfdiag,nmicdiag,nmics,ncplxs)

!  call SumOneLayer('b',KL,micstt,micfor%litrM,totOMbeg,domsb,somb,sorbomb,biomheterb,biomautob,biomresb)
!
  ! write(*,*)'ActiveMicrobeCatabolism'
  call ActiveMicrobeCatabolism(I,J,KL,micfor,micstt,micflx,nmicdiag, &
    naqfdiag,nmicf,nmics,ncplxf,ncplxs)
  !
  !write(*,*)'ChemoDenitrification'
  call ChemoDenitrification(micfor,micstt,nmicdiag,naqfdiag,micflx)

        !
        !write(*,*)'PRIMING of DOC,DON,DOP BETWEEN LITTER AND NON-LITTER C'
  call OMTransferForPriming(KL,micfor,micstt,nmicf,nmics,ncplxf,ncplxs,nmicdiag)
        !
        !     TRANSFER ALL PRIMING AMONG ALL K
        !
        !     TMicHeterActivity=total respiration of DOC+DOA in soil layer
        !     ROQC4HeterMicActCmpK=total respiration of DOC+DOA in substrate complex
        !     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
        !     OMC,OMN,OMP=microbial C,N,P
        !
!
    !write(*,*)'DECOMPOSITION OF ORGANIC SUBSTRATES'
  call SolidOMDecomposition(I,J,KL,micfor,micstt,naqfdiag,nmicdiag,ncplxf,ncplxs,micflx)

  call RedistDecompProduct(micfor,KL,nmicdiag,nmicf,ncplxf,ncplxs,micstt)

  call RDOMSorption(KL,micfor,micstt,nmicf,ncplxf,ncplxs)

  call AutotrophAnabolicUpdate(micfor,micstt,nmicf,nmicdiag,micflx)

  call HeterotrophAnabolicUpdate(I,J,micfor,micstt,nmicf,micflx)
  
  call MicrobialLitterColonization(I,J,KL,micfor,micstt,ncplxf,ncplxs,nmicdiag)
  
  !     AGGREGATE ALL TRANSFOMBioResduATIONS CALCULATED ABOVE FOR EACH N,K
  !
  call AggregateTransfOMBioResdue(KL,micfor,micstt,nmicdiag,naqfdiag,nmicf,ncplxf,micflx)

!  call SumOneLayer('e',KL,micstt,micfor%litrM,totOMend,domse,some,sorbome,biomhetere,biomautoe,biomrese)

  micstt%TotActMicrobiom=nmicdiag%TotActMicrobiom
  !obtain total gross respiration by heterotrophs
  naqfdiag%tRespGrossHeter =sum(nmicf%RespGrossHeter)
  call nmics%destroy()
  call nmicf%destroy()
  call ncplxf%destroy()
  call ncplxs%destroy()
  call PrintInfo('end '//subname)
  end subroutine SoilBGCOneLayer
!------------------------------------------------------------------------------------------
  subroutine SumOneLayer(tag,KL,micstt,litrM,toms,doms,som,sorbom,biomheter,biomauto,biomres)
  !
  !Description:
  !summarize all organic matter states in the specified layer
  implicit none
  character(len=*), intent(in) :: tag
  integer, intent(in) :: KL
  type(micsttype), intent(in) :: micstt  
  logical, intent(in) :: litrM
  real(r8), intent(out) :: toms(1:NumPlantChemElms)
  real(r8), intent(out) :: doms(1:NumPlantChemElms)
  real(r8), intent(out) :: sorbom(1:NumPlantChemElms)
  real(r8), intent(out) :: som(1:NumPlantChemElms)
  real(r8), intent(out) :: biomheter(1:NumPlantChemElms)
  real(r8), intent(out) :: biomauto(1:NumPlantChemElms)
  real(r8), intent(out) :: biomres(1:NumPlantChemElms)
  integer :: K,M,NBM,NE,N,NGL,NLB,MID

!     begin_execution
  associate(                                      &
    SOMPomProtein     => micstt%SOMPomProtein,    & !C/N/P transferred to the first POM component in underlying soil
    SOMHumProtein     => micstt%SOMHumProtein,    & !C/N/P transferred to the first humus component in underlying soil
    SOMHumCarbohyd    => micstt%SOMHumCarbohyd,   & !C/N/P transferred to the second humus component in underlying soil
    SolidOM           => micstt%SolidOM,          & !Solid C/N/P pools by substrate component and complex K
    OMBioResdu        => micstt%OMBioResdu,       & !Microbial residue C/N/P pools by residue compartment and complex K
    SorbedOM          => micstt%SorbedOM,         & !Sorbed organic pools by DOM species and complex K
    DOM               => micstt%DOM,              & !Dissolved organic pools by species (DOC, DON, DOP, acetate) and complex K
    mBiomeHeter       => micstt%mBiomeHeter,      & !C/N/P pools indexed by element, flattened guild/compartment and complex K
    mBiomeAutor       => micstt%mBiomeAutor,      & !C/N/P pools indexed by element and flattened guild/biomass compartment
    NumLiveAutoBioms  => micpar%NumLiveAutoBioms, & !Number of flattened live autotrophic biomass entries (guilds times compartments); not referenced here
    NumLiveHeterBioms => micpar%NumLiveHeterBioms & !Number of flattened live heterotrophic biomass entries per complex; not referenced here
  )

  toms=0._r8
  som=0._r8;biomheter=0._r8;biomauto=0._r8;biomres=0._r8;doms=0._r8;sorbom=0._r8
  if(litrM)then
    DO NE=1,NumPlantChemElms
      toms(NE)=toms(NE)+SOMPomProtein(NE)+SOMHumProtein(NE)+SOMHumCarbohyd(NE)
      som(NE)=som(NE)+SOMPomProtein(NE)+SOMHumProtein(NE)+SOMHumCarbohyd(NE)
    ENDDO
  endif
  
  DO K=1,KL
   !add solid organic matter
    DO M=1,jsken
      DO NE=1,NumPlantChemElms
        TOMS(NE)=TOMS(NE)+SolidOM(NE,M,K)
        som(NE)=som(NE)+SolidOM(NE,M,K)
      ENDDO
    ENDDO

    !add dead microbial biomass
    DO NBM=1,ndbiomcp
      DO NE=1,NumPlantChemElms
        TOMS(NE)    = TOMS(NE)+OMBioResdu(NE,NBM,K)
        biomres(NE) = biomres(NE)+OMBioResdu(NE,NBM,K)
      ENDDO
    ENDDO

    !add sorbed om and DOM
    DO NE=1,NumPlantChemElms
      TOMS(NE)   = TOMS(NE)+SorbedOM(NE,K)
      TOMS(NE)   = TOMS(NE)+DOM(NE,K)
      sorbom(NE) = sorbom(NE)+SorbedOM(NE,K)
      doms(NE)   = doms(NE)+DOM(NE,K)
    ENDDO

    !add live heterotrophic biomass     
    DO N=1,NumMicbHFunGrupsPerCmplx
      if(.not.micpar%is_activeMicrbFungrpHeter(N))cycle
      DO NGL=JGniH(n),JGnfH(n)
        DO nlb=1,micpar%nlbiomcp
          MID=micpar%get_micb_id(nlb,NGL)                
          DO NE=1,NumPlantChemElms
            TOMS(NE)=TOMS(NE)+mBiomeHeter(NE,mid,K)
            biomheter(NE)=biomheter(NE)+mBiomeHeter(NE,mid,K)
          ENDDO
        ENDDO  
      ENDDO
    ENDDO 
    sorbom(ielmc) = sorbom(ielmc)+SorbedOM(idom_acetate,K)
    doms(ielmc)   = doms(ielmc)+DOM(idom_acetate,K)
    TOMS(ielmc)   = TOMS(ielmc)+DOM(idom_acetate,K)+SorbedOM(idom_acetate,K)
  ENDDO

  !add live autotrophic biomass
  DO N=1,NumMicbAFunGrupsPerCmplx
    if(.not.micpar%is_activeMicrbFungrpAutor(N))cycle
    DO NGL=JGniA(N),JGnfA(N)
      DO nlb=1,micpar%nlbiomcp
        MID=micpar%get_micb_id(nlb,NGL)    
        DO NE=1,NumPlantChemElms
          TOMS(NE)=TOMS(NE)+mBiomeAutor(NE,mid)      
          biomauto(NE)=biomauto(NE)+mBiomeAutor(NE,mid)      
        ENDDO 
      ENDDO    
    ENDDO
  ENDDO

  end associate

  end subroutine SumOneLayer
!------------------------------------------------------------------------------------------

  subroutine StageBGCEnvironCondition(I,J,micfor,KL,micstt,naqfdiag,nmicdiag,nmics,ncplxs)
  implicit none
  integer, intent(in) :: I,J
  type(micforctype), intent(in) :: micfor
  integer, intent(out) :: KL                     !total number of complexes to do calculation
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(OMCplx_State_type),intent(inout):: ncplxs
  character(len=*), parameter :: subname='StageBGCEnvironCondition'
  real(r8) :: TBulkOMC  
  integer  :: K,NE
  integer  :: M,N,NGL,MID1,MID2
  real(r8) :: ORGCL
  real(r8) :: TKSO
  real(r8) :: TSolidOMC,TSorbedOMC
  real(r8) :: theta0=0.5_r8            !tuning parameter for litter layer moist volume 
!     begin_execution
  associate(                                                       &
    rCNBiomeActHeter          => nmics%rCNBiomeActHeter,           & !Active heterotrophic biomass nutrient:C ratios (N:C and P:C)
    OMActHeter                => nmics%OMActHeter,                 & !Active microbial C biomass by heterotrophic guild and complex K
    OMC2                      => nmics%OMC2,                       & !Active structural heterotrophic C by guild and complex K
    OMN2                      => nmics%OMN2,                       & !Active structural heterotrophic N by guild and complex K
    FOM2                      => nmics%FOM2,                       & !Active fraction of structural heterotrophic biomass by guild and complex [-]
    FCN                       => nmics%FCN,                        & !N-status multiplier on heterotrophic activity by guild and complex [-]
    FCP                       => nmics%FCP,                        & !P-status multiplier on heterotrophic activity by guild and complex [-]
    FBiomStoiScalarHeter      => nmics%FBiomStoiScalarHeter,       & !Combined N/P stoichiometric multiplier on guild metabolic capacity [-]
    BulkSOMC                  => ncplxs%BulkSOMC,                  & !Colonized solid, residue and sorbed C summed for each complex K
    TOMEAutoK                 => ncplxs%TOMEAutoK,                 & !Total active autotrophic C/N/P summed over guilds
    TOMEK                     => nmicdiag%TOMEK,                   & !Total active heterotrophic C/N/P in each substrate complex K
    FOCA                      => ncplxs%FOCA,                      & !DOC fraction of DOC plus acetate in each substrate complex [-]
    FOAA                      => ncplxs%FOAA,                      & !Acetate-C fraction of DOC plus acetate in each substrate complex [-]
    rCNDOM                    => ncplxs%rCNDOM,                    & !DON:DOC ratio in each substrate complex [-]
    rCPDOM                    => ncplxs%rCPDOM,                    & !DOP:DOC ratio in each substrate complex [-]
    CDOM                      => ncplxs%CDOM,                      & !Dissolved substrate concentrations by DOM species and complex K
    OMBioResduK               => ncplxs%OMBioResduK,               & !Total microbial residue C in each complex K
    SolidOMCK                 => ncplxs%SolidOMCK,                 & !Total solid organic C in each complex K
    SolidOMActK               => ncplxs%SolidOMActK,               & !Total microbially colonized solid C in each complex K
    tMaxNActMicrbK            => ncplxs%tMaxNActMicrbK,            & !Maximum N content of active heterotrophic biomass in each complex
    tMaxPActMicrbK            => ncplxs%tMaxPActMicrbK,            & !Maximum P content of active heterotrophic biomass in each complex
    TOMBioResdu               => nmicdiag%TOMBioResdu,             & !Total microbial residue C summed over substrate complexes
    TotActMicrobiom           => nmicdiag%TotActMicrobiom,         & !Layer total active microbial C across heterotrophs and autotrophs
    TotBiomNO2Consumers       => nmicdiag%TotBiomNO2Consumers,     & !Active C summed over staged nitrite-competing microbial groups
    XCO2                      => nmicdiag%XCO2,                    & !Dissolved CO2 saturation factor [-]
    TSensGrowth               => nmicdiag%TSensGrowth,             & !Layer temperature response for microbial growth [-]
    TSensMaintR               => nmicdiag%TSensMaintR,             & !Layer temperature response for microbial maintenance [-]
    ThetaLitr                 => nmicdiag%ThetaLitr,               & !Effective water content for surface-litter microbial activity
    ThetaZ                    => nmicdiag%ThetaZ,                  & !Effective soil water content above the unavailable-water threshold
    VOLWZ                     => nmicdiag%VOLWZ,                   & !Effective water volume supporting microbial activity and rate constraints
    ZNH4T                     => nmicdiag%ZNH4T,                   & !NH4-N pool in band plus nonband soil
    ZNO3T                     => nmicdiag%ZNO3T,                   & !NO3-N pool in band plus nonband soil
    ZNO2T                     => nmicdiag%ZNO2T,                   & !NO2-N pool in band plus nonband soil
    H2P4T                     => nmicdiag%H2P4T,                   & !H2PO4-P pool in band plus nonband soil
    H1P4T                     => nmicdiag%H1P4T,                   & !HPO4-P pool in band plus nonband soil
    rCNBiomeActAutor          => nmics%rCNBiomeActAutor,           & !Active autotrophic biomass nutrient:C ratios (N:C and P:C)
    OMActAutor                => nmics%OMActAutor,                 & !Active microbial C biomass by autotrophic guild
    OMC2Autor                 => nmics%OMC2Autor,                  & !Active structural autotrophic C by guild
    OMN2Autor                 => nmics%OMN2Autor,                  & !Active structural autotrophic N by guild
    FOM2Autor                 => nmics%FOM2Autor,                  & !Active fraction of structural autotrophic biomass by guild [-]
    FCNAutor                  => nmics%FCNAutor,                   & !N-status multiplier on autotrophic guild activity [-]
    FCPAutor                  => nmics%FCPAutor,                   & !P-status multiplier on autotrophic guild activity [-]
    FBiomNutStoiScalAutor     => nmics%FBiomNutStoiScalAutor,      & !Combined N/P stoichiometric multiplier on guild metabolic capacity [-]
    rNCOMCAutor               => micpar%rNCOMCAutor,               & !Target autotrophic N:C ratios by compartment and guild
    rPCOMCAutor               => micpar%rPCOMCAutor,               & !Target autotrophic P:C ratios by compartment and guild
    rNCOMC                    => micpar%rNCOMC,                    & !Target heterotrophic N:C ratios by compartment, guild and complex
    rPCOMC                    => micpar%rPCOMC,                    & !Target heterotrophic P:C ratios by compartment, guild and complex
    FL                        => micpar%FL,                        & !Target fractions of active biomass in kinetic and structural compartments [-]
    k_humus                   => micpar%k_humus,                   & !Humus complex index
    k_POM                     => micpar%k_POM,                     & !Particulate-organic-matter complex index
    is_activeMicrbFungrpAutor => micpar%is_activeMicrbFungrpAutor, & !Activation flags for autotrophic functional groups
    mid_Facult_DenitBacter    => micpar%mid_Facult_DenitBacter,    & !Functional-group identifier for facultative denitrifying bacteria
    mid_AutoAmmoniaOxidBacter => micpar%mid_AutoAmmoniaOxidBacter, & !Functional-group identifier for ammonia oxidizers
    litrm                     => micfor%litrm,                     & !True for the surface litter layer
    VLSoilPoreMicP            => micfor%VLSoilPoreMicP,            & !Layer micropore volume used in water and aerobic-uptake calculations
    VWatLitRHoldCapcity       => micfor%VWatLitRHoldCapcity,       & !Surface-litter water-holding capacity
    VLWatMicP                 => micfor%VLWatMicP,                 & !Layer micropore water volume used for nutrient donor thresholds; not referenced here
    VOLW0                     => micfor%VOLW0,                     & !Surface-litter volumetric liquid water content
    THETY                     => micfor%THETY,                     & !Water-content threshold subtracted from available soil water
    VLitR                     => micfor%VLitR,                     & !Surface litter volume; not referenced here
    VLSoilMicP                => micfor%VLSoilMicP,                & !Bulk volume associated with the layer micropore domain; not referenced here
    POROS                     => micfor%POROS,                     & !Soil porosity used in the effective-water calculation [-]
    ZEROS                     => micfor%ZEROS,                     & !Small mass or flux threshold used by the routine
    FieldCapacity             => micfor%FieldCapacity,             & !Soil field-capacity water content used to constrain effective water volume
    THETW                     => micfor%THETW,                     & !Soil liquid water content used in the effective-water calculation
    TKS                       => micfor%TKS,                       & !Layer absolute temperature [K]
    TempOffset                => micfor%TempOffset,                & !Temperature offset for acclimation of microbial responses
    VLWatMicPM                => micfor%VLWatMicPM,                & !Micropore water volume at each outer transport substep M
    ZEROS2                    => micfor%ZEROS2,                    & !Small water-volume threshold used by the routine
    ZERO                      => micfor%ZERO,                      & !Small dimensionless or concentration threshold used by the routine
    CCO2S                     => micstt%CCO2S,                     & !Dissolved CO2-C concentration for substrate saturation
    tOMActC                   => micstt%tOMActC,                   & !Total active heterotrophic C summed over substrate complexes
    SolidOM                   => micstt%SolidOM,                   & !Solid C/N/P pools by substrate component and complex K
    SolidOMAct                => micstt%SolidOMAct,                & !Microbially colonized solid C by substrate component and complex K
    TSolidOMActC              => micstt%TSolidOMActC,              & !Layer total colonized solid organic C
    TSolidOMC                 => micstt%TSolidOMC,                 & !Layer total solid organic C
    OMBioResdu                => micstt%OMBioResdu,                & !Microbial residue C/N/P pools by residue compartment and complex K
    SorbedOM                  => micstt%SorbedOM,                  & !Sorbed organic pools by DOM species and complex K
    mBiomeHeter               => micstt%mBiomeHeter,               & !C/N/P pools indexed by element, flattened guild/compartment and complex K
    DOM                       => micstt%DOM,                       & !Dissolved organic pools by species (DOC, DON, DOP, acetate) and complex K
    H1PO4                     => micstt%H1PO4,                     & !HPO4-P pool in nonband soil
    H1POB                     => micstt%H1POB,                     & !HPO4-P pool in fertilizer-band soil
    H2PO4                     => micstt%H2PO4,                     & !H2PO4-P pool in nonband soil
    H2POB                     => micstt%H2POB,                     & !H2PO4-P pool in fertilizer-band soil
    ZNH4B                     => micstt%ZNH4B,                     & !NH4-N pool in fertilizer-band soil
    ZNH4S                     => micstt%ZNH4S,                     & !NH4-N pool in nonband soil
    ZNO2B                     => micstt%ZNO2B,                     & !NO2-N pool in fertilizer-band soil
    ZNO2S                     => micstt%ZNO2S,                     & !NO2-N pool in nonband soil
    ZNO3B                     => micstt%ZNO3B,                     & !NO3-N pool in fertilizer-band soil
    ZNO3S                     => micstt%ZNO3S,                     & !NO3-N pool in nonband soil
    mBiomeAutor               => micstt%mBiomeAutor,               & !C/N/P pools indexed by element and flattened guild/biomass compartment
    FracBulkSOMC              => micstt%FracBulkSOMC               & !Complex share of colonized solid, residue and sorbed C [-]
  )
  call PrintInfo('beg '//subname)
! get KL, the number of mic-om complexes

!
!     TEMPERATURE FUNCTIONS FOR GROWTH AND MAINTENANCE
!     WITH TempOffset FOR THERMAL ADAPTATION
  IF(litrm)THEN
    ! surface litter layer
    KL=micpar%NumOfLitrCmplxs
    IF(VWatLitRHoldCapcity.GT.ZEROS2)THEN
      ThetaLitr = AMIN1(VOLW0/VWatLitRHoldCapcity,1._r8)
      ThetaZ    = AZMAX1(AMIN1(ThetaLitr,FieldCapacity)-THETY)
!      VOLWZ     = ThetaZ**2/(theta0+ThetaZ)*VWatLitRHoldCapcity    !effective water volume, the definition seems not very accurate
!      VOLWZ     = ThetaZ**0.4_r8*VWatLitRHoldCapcity 
      VOLWZ     = ThetaZ*VWatLitRHoldCapcity 
    ELSE
      VOLWZ=0.0_r8
    ENDIF
  ELSE
!     non-surface layer
    KL     = micpar%jcplx
    ThetaZ = AZMAX1((AMIN1(AMAX1(0.5_r8*POROS,FieldCapacity),THETW)-THETY))
!    VOLWZ  = ThetaZ**2/(theta0+ThetaZ)*VLSoilPoreMicP  !effective water volume
!    VOLWZ  = ThetaZ**0.4_r8*VLSoilPoreMicP  !effective water volume
    VOLWZ  = ThetaZ*VLSoilPoreMicP  !effective water volume
  ENDIF
  VOLWZ=real_truncate(VOLWZ,1.e-5_r8)

!     TKS=soil temperature
!     TempOffset=adjustment for acclimation based on MAT in starts.f
!     8.313,710.0=gas constant,enthalpy
!     62500=activation energy
!     197500,195000 low temp inactivation for growth,maintenance
!     222500,232500 high temp inactivation for growth,maintenance
!     TSensGrowth,TSensMaintR=temperature function for growth,maintenance respiration
! the offset could be micobial guild/group specific
  TKSO=TKS+TempOffset

  call MicrobPhysTempFun(TKSO, TSensGrowth, TSensMaintR)

!
!     TOTAL MINERAL NH4, NO3 AND PO4
!
!     allocate NH4, NO3, HPO4, H2PO4 to non-band and band fractions
!
  ZNH4T=AZMAX1(ZNH4S)+AZMAX1(ZNH4B)
  ZNO3T=AZMAX1(ZNO3S)+AZMAX1(ZNO3B)
  H1P4T=AZMAX1(H1PO4)+AZMAX1(H1POB)
  H2P4T=AZMAX1(H2PO4)+AZMAX1(H2POB)
  ZNO2T=AZMAX1(ZNO2S)+AZMAX1(ZNO2B)
!
!     CCO2S=aqueous CO2 concentration
!
  XCO2=CCO2S/(CCO2S+CCKM)
!
!     TOTAL SUBSTRATE
!
!     TSolidOMC=total SOC, TSolidOMActC=total colonized SOC
!     TOMBioResdu=total microbial residue, TSorbedOMC=total adsorbed C
!     in each K:
!     SolidOMCK=total SOC n each K, SolidOMActK=total colonized SOC
!     OMBioResduK=total microbial residue, OHCT=total adsorbed C
!
  TSolidOMC    = 0.0_r8
  TSolidOMActC = 0.0_r8
  TOMBioResdu  = 0.0_r8
  TSorbedOMC   = 0.0_r8
!
!     TOTAL SOLID SUBSTRATE
!
  DO  K=1,KL
    SolidOMCK(K)   =0.0_r8
    SolidOMActK(K) =0.0_r8
    DO M=1,jsken
      SolidOMCK(K)   = SolidOMCK(K)+SolidOM(ielmc,M,K)
      SolidOMActK(K) = SolidOMActK(K)+SolidOMAct(M,K)
    enddo
    TSolidOMC    = TSolidOMC+SolidOMCK(K)
    TSolidOMActC = TSolidOMActC+SolidOMActK(K)
  enddo
!
!     TOTAL BIORESIDUE
!
  DO  K=1,KL
    OMBioResduK(K)=0.0_r8
    DO  M=1,ndbiomcp
      OMBioResduK(K)=OMBioResduK(K)+OMBioResdu(ielmc,M,K)
    ENDDO
    TOMBioResdu=TOMBioResdu+OMBioResduK(K)
!
!     TOTAL ADSORBED AND DISSOLVED SUBSTRATE
!
!     BulkSOMC=total SOC
!
    TSorbedOMC=TSorbedOMC+SorbedOM(ielmc,K)+SorbedOM(idom_acetate,K)
  enddo

  D860: DO K=1,KL
    BulkSOMC(K)=SolidOMActK(K)+OMBioResduK(K)+SorbedOM(ielmc,K)+SorbedOM(idom_acetate,K)
  ENDDO D860
  
  TBulkOMC=TSolidOMActC+TOMBioResdu+TSorbedOMC
!
!     C:N AND C:P RATIOS OF TOTAL BIOMASS
!     CNOMA,CPOMA=N,P contents of active biomass OMA
!     FCN,FCP=effects of N,P limitations on biomass activity
!
  TotActMicrobiom     = 0.0_r8
  TotBiomNO2Consumers = 0.0_r8
  D890: DO K = 1, jcplx
    IF(.not.litrm .OR. (K.NE.k_POM .AND. K.NE.k_humus))THEN
      ! the omb complexes, three biomass components, labile, recalcitrant and reserve
      D895: DO N=1,NumMicbHFunGrupsPerCmplx
        DO NGL=JGniH(n),JGnfH(n)
          MID1=micpar%get_micb_id(iLbiom_kinetic,NGL)
          
          IF(mBiomeHeter(ielmc,MID1,K).GT.ZEROS)THEN
            rCNBiomeActHeter(ielmn,NGL,K)=AZMAX1(mBiomeHeter(ielmn,MID1,K)/mBiomeHeter(ielmc,MID1,K))
            rCNBiomeActHeter(ielmp,NGL,K)=AZMAX1(mBiomeHeter(ielmp,MID1,K)/mBiomeHeter(ielmc,MID1,K))
          ELSE
            rCNBiomeActHeter(ielmn,NGL,K)=rNCOMC(iLbiom_kinetic,NGL,K)
            rCNBiomeActHeter(ielmp,NGL,K)=rPCOMC(iLbiom_kinetic,NGL,K)
          ENDIF
          OMActHeter(NGL,K)           = AZMAX1(mBiomeHeter(ielmc,MID1,K)/FL(iLbiom_kinetic))
          FCN(NGL,K)                  = AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(rCNBiomeActHeter(ielmn,NGL,K)/rNCOMC(iLbiom_kinetic,NGL,K))))
          FCP(NGL,K)                  = AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(rCNBiomeActHeter(ielmp,NGL,K)/rPCOMC(iLbiom_kinetic,NGL,K))))
          FBiomStoiScalarHeter(NGL,K) = AMIN1(FCN(NGL,K),FCP(NGL,K))

!       TOTAL BIOMASS
!       OMC2=active biomass in recalcitrant fraction
!
          TotActMicrobiom=TotActMicrobiom+OMActHeter(NGL,K)
          IF(N.EQ.mid_Facult_DenitBacter .AND. micpar%is_activeMicrbFungrpHeter(N))THEN
            TotBiomNO2Consumers=TotBiomNO2Consumers+OMActHeter(NGL,K)
          ENDIF
          MID2=micpar%get_micb_id(iLbiom_struct,NGL)
          OMC2(NGL,K)=AZMAX1(AMIN1(OMActHeter(NGL,K)*FL(iLbiom_struct),mBiomeHeter(ielmc,MID2,K)))
          IF(mBiomeHeter(ielmc,MID2,K).GT.ZEROS)THEN
            FOM2(NGL,K)=AZMAX1(OMC2(NGL,K)/mBiomeHeter(ielmc,MID2,K))
            OMN2(NGL,K)=AZMAX1(FOM2(NGL,K)*mBiomeHeter(ielmn,MID2,K))
          ELSE
            FOM2(NGL,K)=0.0_r8
            OMN2(NGL,K)=0.0_r8
          ENDIF
        ENDDO
      ENDDO D895
    ENDIF
  ENDDO D890

! the autotrohpic complex
  DO N=1,NumMicbAFunGrupsPerCmplx
    IF(is_activeMicrbFungrpAutor(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        MID1=micpar%get_micb_id(iLbiom_kinetic,NGL)
        IF(mBiomeAutor(ielmc,MID1).GT.ZEROS)THEN
          rCNBiomeActAutor(ielmn,NGL)=AZMAX1(mBiomeAutor(ielmn,MID1)/mBiomeAutor(ielmc,MID1))
          rCNBiomeActAutor(ielmp,NGL)=AZMAX1(mBiomeAutor(ielmp,MID1)/mBiomeAutor(ielmc,MID1))
        ELSE
          rCNBiomeActAutor(ielmn,NGL)=rNCOMCAutor(iLbiom_kinetic,NGL)
          rCNBiomeActAutor(ielmp,NGL)=rPCOMCAutor(iLbiom_kinetic,NGL)
        ENDIF
        OMActAutor(NGL)           = AZMAX1(mBiomeAutor(ielmc,MID1)/FL(iLbiom_kinetic))
        FCNAutor(NGL)             = AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(rCNBiomeActAutor(ielmn,NGL)/rNCOMCAutor(iLbiom_kinetic,NGL))))
        FCPAutor(NGL)             = AMIN1(1.0_r8,AMAX1(0.50_r8,SQRT(rCNBiomeActAutor(ielmp,NGL)/rPCOMCAutor(iLbiom_kinetic,NGL))))
        FBiomNutStoiScalAutor(NGL) = AMIN1(FCNAutor(NGL),FCPAutor(NGL))
      !
      !       TOTAL BIOMASS
      !       OMC2=active biomass in recalcitrant fraction
      !
        TotActMicrobiom=TotActMicrobiom+OMActAutor(NGL)

        !All active NO2 consumers share the same layer-wide fallback budget.
        IF(N.EQ.mid_AutoAmmoniaOxidBacter .OR. N.EQ.micpar%mid_AutoNitriteOxidBacter &
          .OR. N.EQ.micpar%mid_AutoAMONC10)THEN
          TotBiomNO2Consumers=TotBiomNO2Consumers+OMActAutor(NGL)
        ENDIF

        MID2           = micpar%get_micb_id(iLbiom_struct,NGL)
        OMC2Autor(NGL) = AZMAX1(AMIN1(OMActAutor(NGL)*FL(iLbiom_struct),mBiomeAutor(ielmc,MID2)))
        IF(mBiomeAutor(ielmc,MID2).GT.ZEROS)THEN
          FOM2Autor(NGL)=AZMAX1(OMC2Autor(NGL)/mBiomeAutor(ielmc,MID2))
          OMN2Autor(NGL)=AZMAX1(FOM2Autor(NGL)*mBiomeAutor(ielmn,MID2))
        ELSE
          FOM2Autor(NGL)=0.0_r8
          OMN2Autor(NGL)=0.0_r8
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  tOMActC=0._r8
  D690: DO K=1,KL
    TOMEK(:,K)        = 0.0_r8
    tMaxNActMicrbK(K) = 0.0_r8
    tMaxPActMicrbK(K) = 0.0_r8
    D685: DO N=1,NumMicbHFunGrupsPerCmplx
      DO NGL=JGniH(N),JGnfH(N)
        if(OMActHeter(NGL,K)>ZEROS)THEN
          TOMEK(ielmc,K)    = TOMEK(ielmc,K)+OMActHeter(NGL,K)
          DO NE=2,NumPlantChemElms
            TOMEK(NE,K)    = TOMEK(NE,K)+OMActHeter(NGL,K)*rCNBiomeActHeter(NE,NGL,K)
          ENDDO
          tMaxNActMicrbK(K) = tMaxNActMicrbK(K)+OMActHeter(NGL,K)*rNCOMC(iLbiom_kinetic,NGL,K)   !maximum total N in active micb
          tMaxPActMicrbK(K) = tMaxPActMicrbK(K)+OMActHeter(NGL,K)*rPCOMC(iLbiom_kinetic,NGL,K)   !maximum total P in active micb
        ENDIF
      ENDDO
    ENDDO D685
    tOMActC=tOMActC+TOMEK(ielmc,K)
  ENDDO D690
  
  TOMEAutoK(:)      = 0._r8
  DO N=1,NumMicbAFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      if(OMActAutor(NGL)>ZEROS)then
        TOMEAutoK(ielmc)  = TOMEAutoK(ielmc)+OMActAutor(NGL)
        TOMEAutoK(ielmn)  = TOMEAutoK(ielmn)+OMActAutor(NGL)*rCNBiomeActAutor(ielmn,NGL)
        TOMEAutoK(ielmp)  = TOMEAutoK(ielmp)+OMActAutor(NGL)*rCNBiomeActAutor(ielmp,NGL)
      endif
    ENDDO
  ENDDO
!
!     FracBulkSOMC=fraction of total SOC in each substrate complex K
!
  D790: DO K=1,KL
    IF(TBulkOMC.GT.ZEROS)THEN
      FracBulkSOMC(K)=BulkSOMC(K)/TBulkOMC
    ELSE
      FracBulkSOMC(K)=1.0_r8
    ENDIF
    !
    !     DOC CONCENTRATIONS
    !
    !     COQC,COQA=aqueous DOC,acetate concentrations
    !     VLWatMicPM=soil water content, FracBulkSOMC=fraction of total SOC
    !     occupied by each substrate complex K
    !
    IF(VLWatMicPM(NPH).GT.ZEROS2)THEN
      IF(FracBulkSOMC(K).GT.ZERO)THEN
        CDOM(idom_doc,K)     = AZMAX1(DOM(idom_doc,K)/(VLWatMicPM(NPH)*FracBulkSOMC(K)))
        CDOM(idom_acetate,K) = AZMAX1(DOM(idom_acetate,K)/(VLWatMicPM(NPH)*FracBulkSOMC(K)))
      ELSE
        CDOM(idom_doc,K)     = AZMAX1(DOM(idom_doc,K)/VLWatMicPM(NPH))
        CDOM(idom_acetate,K) = AZMAX1(DOM(idom_acetate,K)/VLWatMicPM(NPH))
      ENDIF
    ELSE
      CDOM(idom_doc,K)     = 0.0_r8
      CDOM(idom_acetate,K) = 0.0_r8
    ENDIF
!
!     rCNDOM,rCPDOM=DON:DOC,DOP:DOC,FOCA,FOAA=DOC,DOA:(DOC+DOA)
!
    IF(DOM(idom_doc,K).GT.ZEROS)THEN
      rCNDOM(K)=AZMAX1(DOM(idom_don,K)/DOM(idom_doc,K))
      rCPDOM(K)=AZMAX1(DOM(idom_dop,K)/DOM(idom_doc,K))
    ELSE
      rCNDOM(K)=0.0_r8
      rCPDOM(K)=0.0_r8
    ENDIF
    IF(DOM(idom_doc,K).GT.ZEROS.AND.DOM(idom_acetate,K).GT.ZEROS)THEN
      FOCA(K)=DOM(idom_doc,K)/(DOM(idom_doc,K)+DOM(idom_acetate,K))
      FOAA(K)=1.0_r8-FOCA(K)
    ELSEIF(DOM(idom_doc,K).GT.ZEROS)THEN
      FOCA(K)=1.0_r8
      FOAA(K)=0.0_r8
    ELSE
      FOCA(K)=0.0_r8
      FOAA(K)=1.0_r8
    ENDIF
  ENDDO D790
!
  call PrintInfo('end '//subname)
  end associate
  end subroutine StageBGCEnvironCondition

!------------------------------------------------------------------------------------------
  subroutine GetMicrobDensFactorHeter(N,K,micfor, micstt, ORGCL, SPOMK, RMOMK)

  implicit none  
  integer, intent(in) :: N,K
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt    
  real(r8), intent(in) :: ORGCL
  real(r8), intent(out) :: SPOMK(2)
  real(r8), intent(out)  :: RMOMK(2)    

  integer :: NGL,M,MID  
  real(r8) :: TOMCNK(2)   !total microbial C in labile (M=1) and recalcitrant (M=2) biomass components 
  real(r8) :: COMC

  associate(                             &
    ZEROS       => micfor%ZEROS,         & !Small mass or flux threshold used by the routine
    mBiomeHeter => micstt%mBiomeHeter    & !C/N/P pools indexed by element, flattened guild/compartment and complex K
  )

  TOMCNK(:)=0.0_r8
  DO NGL=JGniH(N),JGnfH(N)
    DO M=1,2
      MID=micpar%get_micb_id(M,NGL)          
      TOMCNK(M)=TOMCNK(M)+mBiomeHeter(ielmc,MID,K)
    ENDDO
  ENDDO  

!     ADJUST MCROBIAL GROWTH AND DECOMPOSITION RATES FOR BIOMASS
!
!     COMC=microbial C concentration relative to substrate
!     SPOMK=effect of microbial C concentration on microbial decay
!     RMOMK=effect of microbial C concentration on maintenance respn
!

  IF(ORGCL.GT.ZEROS)THEN
    DO M=1,2
      COMC     = TOMCNK(M)/ORGCL
      SPOMK(M) = COMC/(COMC+COMKI)
      RMOMK(M) = COMC/(COMC+COMKM)
    ENDDO
  ELSE
    DO M=1,2
      SPOMK(M)=1.0_r8
      RMOMK(M)=1.0_r8
    ENDDO
  ENDIF
  end associate
  end subroutine GetMicrobDensFactorHeter
!------------------------------------------------------------------------------------------
  subroutine GetMicrobDensFactorAutor(N, micfor, micstt, ORGCL, SPOMK, RMOMK)

  implicit none  
  integer, intent(in) :: N
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt  
  real(r8), intent(in) :: ORGCL
  real(r8), intent(out) :: SPOMK(2)
  real(r8), intent(out)  :: RMOMK(2)    

  integer :: NGL,M,MID  
  real(r8) :: TOMCNK(2),COMC

  associate(                             &
    ZEROS       => micfor%ZEROS,         & !Small mass or flux threshold used by the routine
    mBiomeAutor => micstt%mBiomeAutor    & !C/N/P pools indexed by element and flattened guild/biomass compartment
  )

  TOMCNK(:)=0.0_r8
  DO NGL=JGniA(N),JGnfA(N)
    DO M=1,2
      MID=micpar%get_micb_id(M,NGL)          
      TOMCNK(M)=TOMCNK(M)+mBiomeAutor(ielmc,MID)
    ENDDO
  ENDDO  

!     ADJUST MCROBIAL GROWTH AND DECOMPOSITION RATES FOR BIOMASS
!
!     COMC=microbial C concentration relative to substrate
!     SPOMK=effect of microbial C concentration on microbial decay
!     RMOMK=effect of microbial C concentration on maintenance respn
!

  IF(ORGCL.GT.ZEROS)THEN
    DO M=1,2
      COMC=TOMCNK(M)/ORGCL
      SPOMK(M)=COMC/(COMC+COMKI)
      RMOMK(M)=COMC/(COMC+COMKM)
    ENDDO
  ELSE
    DO M=1,2
      SPOMK(M)=1.0_r8
      RMOMK(M)=1.0_r8
    ENDDO
  ENDIF
  end associate
  end subroutine GetMicrobDensFactorAutor
!------------------------------------------------------------------------------------------

  subroutine ActiveMicrobeCatabolism(I,J,KL,micfor,micstt,micflx,nmicdiag,naqfdiag, &
    nmicf, nmics,ncplxf,ncplxs)
  !
  !  Description:
  !
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: KL                !number of complexes involved in calculation
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout):: nmics
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout):: ncplxs
  integer :: K,M,N,NGL,MID
  real(r8) :: SPOMK(2)
  real(r8) :: RMOMK(2)    
  character(len=*), parameter :: subname='ActiveMicrobeCatabolism'
  REAL(R8) :: OXKX
  real(r8) :: ORGCL
! begin_execution
  associate(                                                     &
  TempMaintRHeter           => nmics%TempMaintRHeter,            & !Guild temperature multiplier on heterotrophic maintenance [-]; not referenced here
  GrowthEnvScalHeter        => nmics%GrowthEnvScalHeter,         & !Temperature and water-potential multiplier on heterotrophic growth [-]; not referenced here
  OMActHeter                => nmics%OMActHeter,                 & !Active microbial C biomass by heterotrophic guild and complex K; not referenced here
  TDOMUptkHeter             => ncplxf%TDOMUptkHeter,             & !Total heterotrophic DOM uptake by species and complex K
  k_humus                   => micpar%k_humus,                   & !Humus complex index; not referenced here
  k_POM                     => micpar%k_POM,                     & !Particulate-organic-matter complex index; not referenced here
  is_activeMicrbFungrpAutor => micpar%is_activeMicrbFungrpAutor, & !Activation flags for autotrophic functional groups
  is_activeMicrbFungrpHeter => micpar%is_activeMicrbFungrpHeter, & !Activation flags for heterotrophic functional groups
  litrm                     => micfor%litrm,                     & !True for the surface litter layer; not referenced here
  H1PO4                     => micstt%H1PO4,                     & !HPO4-P pool in nonband soil; not referenced here
  H1POB                     => micstt%H1POB,                     & !HPO4-P pool in fertilizer-band soil; not referenced here
  H2PO4                     => micstt%H2PO4,                     & !H2PO4-P pool in nonband soil; not referenced here
  H2POB                     => micstt%H2POB,                     & !H2PO4-P pool in fertilizer-band soil; not referenced here
  ORGC                      => micfor%ORGC,                      & !Layer organic-C reference used in microbial density calculations
  SoilMicPMassLayer         => micfor%SoilMicPMassLayer,         & !Soil mass associated with the current layer micropore domain
  mBiomeAutor               => micstt%mBiomeAutor,               & !C/N/P pools indexed by element and flattened guild/biomass compartment; not referenced here
  mBiomeHeter               => micstt%mBiomeHeter                & !C/N/P pools indexed by element, flattened guild/compartment and complex K; not referenced here
  )
  call PrintInfo('beg '//subname)
  !Heterotrophic microbes
  TDOMUptkHeter(:,:) = 0.0_r8
  ORGCL              = AMIN1(1.0E+05_r8*SoilMicPMassLayer,ORGC)

  !heterotrophs
  D760: DO K=1,KL

    DO  N=1,NumMicbHFunGrupsPerCmplx
      if(.not.is_activeMicrbFungrpHeter(N))cycle
      call GetMicrobDensFactorHeter(N,K,micfor, micstt, ORGCL,SPOMK,RMOMK)
      
      call ActiveHeterotrophsK(I,J,N,K,SPOMK,RMOMK,&
        micfor,micstt,naqfdiag,nmicf,nmics,ncplxf,ncplxs,micflx,nmicdiag)
    ENDDO
  ENDDO D760

! Autotrophic microbes
  
  DO  N=1,NumMicbAFunGrupsPerCmplx
    IF(.not.is_activeMicrbFungrpAutor(N))cycle

    call GetMicrobDensFactorAutor(N,micfor, micstt, ORGCL,SPOMK,RMOMK)

    call ActiveAutotrophs(I,J,N,SPOMK, RMOMK, &
      micfor,micstt,micflx,naqfdiag,nmicf,nmics,ncplxf,ncplxs,nmicdiag)

  ENDDO

  !summarize microbial activity as a proxy for hydrolysis
  DO  K=1,KL
    DO  N=1,NumMicbHFunGrupsPerCmplx
      DO NGL=JGniH(N),JGnfH(N)
        nmicdiag%ROQC4HeterMicActCmpK(K)=nmicdiag%ROQC4HeterMicActCmpK(K)+nmicf%ROQC4HeterMicrobAct(NGL,K)
      enddo
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine ActiveMicrobeCatabolism
!------------------------------------------------------------------------------------------
  subroutine ActiveHeterotrophsK(I,J,N,K,SPOMK,RMOMK,&
     micfor,micstt,naqfdiag,nmicf,nmics,ncplxf,ncplxs,micflx,nmicdiag)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: K,N
  real(r8), intent(in) :: SPOMK(2)        
  real(r8), intent(in) :: RMOMK(2)      
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout) :: ncplxs
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Diag_type), intent(inout) :: nmicdiag  
  character(len=*), parameter :: subname='ActiveHeterotrophsK'
  real(r8) :: FNH4X
  real(r8) :: FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX
  real(r8) :: RMaintDefcitcitHeter
  integer  :: NGL

! begin_execution
  associate(                                             &
    ZNH4T               => nmicdiag%ZNH4T,               & !NH4-N pool in band plus nonband soil
    ZNO3T               => nmicdiag%ZNO3T,               & !NO3-N pool in band plus nonband soil
    ZNO2T               => nmicdiag%ZNO2T,               & !NO2-N pool in band plus nonband soil
    H2P4T               => nmicdiag%H2P4T,               & !H2PO4-P pool in band plus nonband soil
    H1P4T               => nmicdiag%H1P4T,               & !HPO4-P pool in band plus nonband soil
    VOLWZ               => nmicdiag%VOLWZ,               & !Effective water volume supporting microbial activity and rate constraints
    TotBiomNO2Consumers => nmicdiag%TotBiomNO2Consumers, & !Active C summed over staged nitrite-competing microbial groups; not referenced here
    OMActHeter          => nmics%OMActHeter,             & !Active microbial C biomass by heterotrophic guild and complex K
    litrm               => micfor%litrm,                 & !True for the surface litter layer
    ZEROS               => micfor%ZEROS,                 & !Small mass or flux threshold used by the routine
    VLSoilPoreMicP      => micfor%VLSoilPoreMicP         & !Layer micropore volume used in water and aerobic-uptake calculations
  )
  call PrintInfo('beg '//subname)
  !
  ! HETEROTROPHIC BIOMASS RESPIRATION

  IF(micpar%is_aerobic_hetr(N))THEN
    !   RESPIRATION BY HETEROTROPHIC AEROBES:
    !   N=(1)OBLIGATE AEROBES,(2)FACULTATIVE ANAEROBES,(3)FUNGI, (6)N2 FIXERS

    IF(N.EQ.micpar%mid_Aerob_Fungi)THEN
      !(3)FUNGI
      call AerobicFungiCatabolism(I,J,N,K,RMOMK,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
    ELSEIF(N.EQ.micpar%mid_HeterMixtCynoBacter)THEN
      call CyanoBacteriaCatabolism(I,J,N,K,RMOMK,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
    else  
      !heterotrophic aerboic bacteria, including facultative denitrifiers
      call AerobicHeteroBactCatabolism(I,J,N,K,RMOMK,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
    endif

    IF(N.EQ.micpar%mid_Facult_DenitBacter .AND. (.not.litrm .OR. VLSoilPoreMicP.GT.ZEROS))THEN
      !non-litter layer denitrifcation
      call HeteroDenitrificCatabolism(N,K,VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx)
    ENDIF

  ELSEIF(micpar%is_anaerobic_hetr(N))THEN
    !     RESPIRATION BY HETEROTROPHIC ANAEROBES:
    !     N=(4)ACETOGENIC FERMENTERS (7) ACETOGENIC N2 FIXERS
    call AcetogFermentCatabolism(N,K,RMOMK,micfor,micstt,naqfdiag,ncplxs,nmicf,nmics,micflx,nmicdiag)

  ELSEIF(N.EQ.micpar%mid_HeterAcetoCH4GenArchea)THEN
    !     ENERGY YIELD FROM ACETOTROPHIC METHANOGENESIS
    call AcetoMethanogenCatabolism(N,K,RMOMK,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
  ENDIF

  !
  DO NGL=JGniH(N),JGnfH(N)        
    IF(OMActHeter(NGL,K).LE.0.0_r8)cycle
    !     BIOMASS DECOMPOSITION AND MINERALIZATION
    !
    call SubstrateAttenf4Compet(NGL,N,K,FNH4X, FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
      micfor,naqfdiag,nmicf,nmics,micflx)

    call BiomassMineralization(NGL,N,K,FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX, &
      ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt, nmicf,nmics,micflx)

    call GatherHeterotrophRespiration(I,J,NGL,N,K,RMaintDefcitcitHeter, &
      micfor,micstt,nmicf,nmics,micflx)

    call GatherHetertrophAnabolicFlux(I,J,NGL,N,K,&
      RMaintDefcitcitHeter,SPOMK,micfor,micstt,nmicf, nmics,ncplxf,ncplxs,micflx)
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine ActiveHeterotrophsK

!------------------------------------------------------------------------------------------

  subroutine ChemoDenitrification(micfor,micstt,nmicdiag,naqfdiag,micflx)
  implicit none
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  type(Cumlate_Flux_Diag_type), INTENT(INOUT):: naqfdiag
  type(micfluxtype), intent(inout) :: micflx
  real(r8) :: H_1p_conc,CHNO2,CHNOB
  real(r8) :: FNO3S,FNO3B
  REAL(R8) :: FNO2,FNB2
  real(r8) :: VMXC4S,VMXC4B
  character(len=*), parameter :: subname='ChemoDenitrification'
!     begin_execution
  associate(                                               &
    TSensGrowth           => nmicdiag%TSensGrowth,         & !Layer temperature response for microbial growth [-]
    RNO2ReduxSoilChemo    => naqfdiag%RNO2ReduxSoilChemo,  & !Realized chemical nonband NO2-N consumption
    RNO2ReduxBandChemo    => naqfdiag%RNO2ReduxBandChemo,  & !Realized chemical fertilizer-band NO2-N consumption
    RN2OProdSoilChemo     => naqfdiag%RN2OProdSoilChemo,   & !N2O-N production by chemical nitrite decomposition in nonband soil
    RN2OProdBandChemo     => naqfdiag%RN2OProdBandChemo,   & !N2O-N production by chemical nitrite decomposition in fertilizer-band soil
    RNO3ProdSoilChemo     => naqfdiag%RNO3ProdSoilChemo,   & !NO3-N production by chemical nitrite decomposition in nonband soil
    RNO3ProdBandChemo     => naqfdiag%RNO3ProdBandChemo,   & !NO3-N production by chemical nitrite decomposition in fertilizer-band soil
    RNO2ReduxChemo        => naqfdiag%RNO2ReduxChemo,      & !DON-N production from chemical nitrite reduction, summed over band and nonband
    RNO2EcoUptkSoilPrev   => micfor%RNO2EcoUptkSoilPrev,   & !Previous-hour ecosystem nonband NO2-N uptake; competition denominator
    RNO2EcoUptkBandPrev   => micfor%RNO2EcoUptkBandPrev,   & !Previous-hour ecosystem fertilizer-band NO2-N uptake; competition denominator
    ph                    => micfor%pH,                    & !Layer pH used in microbial rate and nutrient-speciation responses
    VLWatMicPM            => micfor%VLWatMicPM,            & !Micropore water volume at each outer transport substep M
    ZEROS                 => micfor%ZEROS,                 & !Small mass or flux threshold used by the routine
    ZERO                  => micfor%ZERO,                  & !Small dimensionless or concentration threshold used by the routine
    VLNOB                 => micfor%VLNOB,                 & !Fertilizer-band fraction for nitrate/nitrite pools and uptake capacity [-]
    VLNO3                 => micfor%VLNO3,                 & !Nonband fraction for nitrate/nitrite pools and uptake capacity [-]
    CNO2S                 => micstt%CNO2S,                 & !Dissolved NO2-N concentration in nonband soil
    CNO2B                 => micstt%CNO2B,                 & !Dissolved NO2-N concentration in fertilizer-band soil
    ZNO2B                 => micstt%ZNO2B,                 & !NO2-N pool in fertilizer-band soil
    ZNO2S                 => micstt%ZNO2S,                 & !NO2-N pool in nonband soil
    RNO2DmndSoilChemoPrev => micflx%RNO2DmndSoilChemoPrev, & !Previous-hour potential chemical nonband NO2-N consumption used for competition; used for competition
    RNO2DmndBandChemoPrev => micflx%RNO2DmndBandChemoPrev, & !Previous-hour potential chemical fertilizer-band NO2-N consumption used for competition; used for competition
    RNO2DmndSoilChemo     => micflx%RNO2DmndSoilChemo,     & !Potential chemical nonband NO2-N consumption used for competition
    RNO2DmndBandChemo     => micflx%RNO2DmndBandChemo      & !Potential chemical fertilizer-band NO2-N consumption used for competition
  )
  call PrintInfo('beg '//subname)
!
!     FNO2,FNB2=fraction of total NO2 demand in non-band,band
!     VMXC4S,VMXC4B=substrate-unlimited NO2 reduction in non-band,band
!     CHNO2,CHNOB=nitrous acid concentration in non-band,band
!     VLWatMicPM=soil water content
!     FNO3S,FNO3B=fractions of NO2 in non-band,band
!     TSensGrowth=temperature stress function
!     RNO2ReduxSoilChemo,RNO2ReduxBandChemo=substrate-limited nitrous acid reduction in non-band,band
!     RN2OProdSoilChemo,RN2OProdBandChemo=N2O production from nitrous acid reduction in non-band,band
!     RNO3ProdSoilChemo,RNO3ProdBandChemo=NO3 production from nitrous acid reduction in non-band,band
!     RNO2ReduxChemo=DON production from nitrous acid reduction
!     RNO2DmndSoilChemo,RNO2DmndBandChemo=demand for NO2 reduction in non-band,band
!     nitrous acid concn CHNO2
  H_1p_conc = AMAX1(ZERO,10.0_r8**(-(PH-3.0_r8)))
  CHNO2     = CNO2S*H_1p_conc/0.5_r8
  CHNOB     = CNO2B*H_1p_conc/0.5_r8

  IF(RNO2EcoUptkSoilPrev.GT.ZEROS)THEN
    FNO2=AMAX1(FMN,RNO2DmndSoilChemoPrev/RNO2EcoUptkSoilPrev)
  ELSE
    FNO2=FMN*VLNO3
  ENDIF
  IF(RNO2EcoUptkBandPrev.GT.ZEROS)THEN
    FNB2=AMAX1(FMN,RNO2DmndBandChemoPrev/RNO2EcoUptkBandPrev)
  ELSE
    FNB2=FMN*VLNOB
  ENDIF
  naqfdiag%TFNO2X    = naqfdiag%TFNO2X+FNO2
  naqfdiag%TFNO2B    = naqfdiag%TFNO2B+FNB2
  FNO3S              = VLNO3
  FNO3B              = VLNOB
  VMXC4S             = 7.5E-02_r8*CHNO2*VLWatMicPM(NPH)*FNO3S*TSensGrowth
  VMXC4B             = 7.5E-02_r8*CHNOB*VLWatMicPM(NPH)*FNO3B*TSensGrowth
  RNO2ReduxSoilChemo = AZMAX1(AMIN1(ZNO2S*FNO2,VMXC4S))
  RNO2ReduxBandChemo = AZMAX1(AMIN1(ZNO2B*FNB2,VMXC4B))
  RN2OProdSoilChemo  = 0.10_r8*RNO2ReduxSoilChemo
  RN2OProdBandChemo  = 0.10_r8*RNO2ReduxBandChemo
  RNO3ProdSoilChemo  = 0.80_r8*RNO2ReduxSoilChemo
  RNO3ProdBandChemo  = 0.80_r8*RNO2ReduxBandChemo
  RNO2ReduxChemo     = 0.10_r8*(RNO2ReduxSoilChemo+RNO2ReduxBandChemo)
  RNO2DmndSoilChemo  = VMXC4S
  RNO2DmndBandChemo  = VMXC4B
  call PrintInfo('end '//subname)
  end associate
  end subroutine ChemoDenitrification
!------------------------------------------------------------------------------------------

  subroutine OMTransferForPriming(KL,micfor,micstt,nmicf,nmics,ncplxf,ncplxs,nmicdiag)

  implicit none
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout) :: nmics
  type(OMCplx_Flux_type), intent(inout):: ncplxf
  type(OMCplx_State_type),intent(inout):: ncplxs
  type(Microbe_Diag_type), intent(inout) :: nmicdiag  
  integer  :: K,M,N,KK,NGL,MID,idom,NE
  real(r8) :: OSRT
  real(r8) :: XFRK,XFROM(idom_beg:idom_end)
  real(r8) :: XFME
  character(len=*), parameter :: subname='OMTransferForPriming'

!     begin_execution
  associate(                                               &
    GrowthEnvScalHeter   => nmics%GrowthEnvScalHeter,      & !Temperature and water-potential multiplier on heterotrophic growth [-]
    XferBiomeHeterK      => nmicf%XferBiomeHeterK,         & !Net C/N/P transfer into each guild/complex from inter-complex priming
    ROQC4HeterMicActCmpK => nmicdiag%ROQC4HeterMicActCmpK, & !Complex-total heterotrophic activity proxy driving organic-matter hydrolysis
    XferRespHeterK       => ncplxf%XferRespHeterK,         & !Net redistribution of heterotrophic activity proxy among complexes
    XferDOMK             => ncplxf%XferDOMK,               & !Net transfer into each DOM species/complex from inter-complex priming
    BulkSOMC             => ncplxs%BulkSOMC,               & !Colonized solid, residue and sorbed C summed for each complex K
    TMicHeterActivity    => micstt%TMicHeterActivity,      & !Layer-total heterotrophic activity proxy used in hydrolysis and priming
    DOM                  => micstt%DOM,                    & !Dissolved organic pools by species (DOC, DON, DOP, acetate) and complex K
    mBiomeHeter          => micstt%mBiomeHeter,            & !C/N/P pools indexed by element, flattened guild/compartment and complex K
    ZEROS                => micfor%ZEROS,                  & !Small mass or flux threshold used by the routine
    TScal4Difsvity       => micfor%TScal4Difsvity          & !Temperature multiplier on inter-complex priming transfers [-]
  )
  call PrintInfo('beg '//subname)
!
!     BulkSOMC=total SOC in each K
!     XFRK,XFRC,XFRN,XFRP,XFRA=transfer of respiration,DOC,DON,DOP,acetate
!     between each K and KK, FPRIM=priming transfer rate constant
!     TScal4Difsvity=temperature effect on priming transfers
!     ROQC4HeterMicActCmpK,OQC,OQN,OQP=respiration,DOC,DON,DOP
!     XferRespHeterK,XOQCZ,XOQNZ,XOQPZ,XOQAZ=total XFRK,XFRC,XFRN,XFRP,XFRA for all K
!
  D795: DO K=1,KL
    IF(K.LE.KL-1)THEN
      D800: DO KK=K+1,KL
        OSRT=BulkSOMC(K)+BulkSOMC(KK)
        IF(BulkSOMC(K).GT.ZEROS.AND.BulkSOMC(KK).GT.ZEROS)THEN
          XFRK=FPRIM*TScal4Difsvity*(ROQC4HeterMicActCmpK(K)*BulkSOMC(KK)-ROQC4HeterMicActCmpK(KK)*BulkSOMC(K))/OSRT
          DO idom=idom_beg,idom_end
            XFROM(idom)=AMAX1(AMIN1(FPRIM*TScal4Difsvity*(DOM(idom,K)*BulkSOMC(KK)-DOM(idom,KK)*BulkSOMC(K))/OSRT,DOM(idom,K)),-DOM(idom,KK))
          ENDDO
          IF(ROQC4HeterMicActCmpK(K)+XferRespHeterK(K)-XFRK.GT.0.0_r8 .AND. ROQC4HeterMicActCmpK(KK)+XferRespHeterK(KK)+XFRK.GT.0.0_r8)THEN
            XferRespHeterK(K)  = XferRespHeterK(K)-XFRK
            XferRespHeterK(KK) = XferRespHeterK(KK)+XFRK
          ENDIF
          DO iDOM=idom_beg,idom_end
            IF(DOM(idom,K)+XferDOMK(idom,K)-XFROM(idom).GT.0.0_r8 &
              .AND.DOM(idom,KK)+XferDOMK(idom,KK)+XFROM(idom).GT.0.0_r8)THEN
              XferDOMK(idom,K)  = XferDOMK(idom,K)-XFROM(idom)
              XferDOMK(idom,KK) = XferDOMK(idom,KK)+XFROM(idom)
            ENDIF
          ENDDO
!
!     PRIMING of MICROBIAL C,N,P BETWEEN LITTER AND NON-LITTER C
!
!     XFMC,XFMN,XFMP=transfer of microbial C,N,P
!     between each K and KK, FPRIMM=priming transfer rate constant
!     GrowthEnvScalHeter=temperature+water effect
!     OMC,OMN,OMP=microbial C,N,P
!     BulkSOMC=total SOC in each K
!     XOMCZ,XOMNZ,XOMPZ=total microbial C,N,P transfer for all K
!
          D850: DO N=1,NumMicbHFunGrupsPerCmplx
            DO  M=1,nlbiomcp
              DO NGL=JGniH(N),JGnfH(N)
                MID=micpar%get_micb_id(M,NGL)
                DO NE=1,NumPlantChemElms
                  XFME=FPRIMM*GrowthEnvScalHeter(NGL,K)*(mBiomeHeter(NE,MID,K)*BulkSOMC(KK) &
                    -mBiomeHeter(NE,MID,KK)*BulkSOMC(K))/OSRT
                  IF(mBiomeHeter(NE,MID,K)+XferBiomeHeterK(NE,M,NGL,K)-XFME.GT.0.0_r8 &
                    .AND.mBiomeHeter(NE,MID,KK)+XferBiomeHeterK(NE,M,NGL,KK)+XFME.GT.0.0_r8)THEN
                    XferBiomeHeterK(NE,M,NGL,K)  = XferBiomeHeterK(NE,M,NGL,K)-XFME
                    XferBiomeHeterK(NE,M,NGL,KK) = XferBiomeHeterK(NE,M,NGL,KK)+XFME
                  ENDIF
                ENDDO
              enddo
            enddo
          ENDDO D850
        ENDIF
      ENDDO D800
    ENDIF
  ENDDO D795
!
!     TRANSFER ALL PRIMING AMONG ALL K
!
!     TMicHeterActivity=total respiration of DOC+DOA in soil layer
!     ROQC4HeterMicActCmpK=total respiration of DOC+DOA in substrate complex
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     OMC,OMN,OMP=microbial C,N,P
!
  TMicHeterActivity=0.0_r8
  D840: DO K=1,KL
    ROQC4HeterMicActCmpK(K) = ROQC4HeterMicActCmpK(K)+XferRespHeterK(K)
    TMicHeterActivity       = TMicHeterActivity+ROQC4HeterMicActCmpK(K)
        
    DO idom=idom_beg,idom_end
      DOM(idom,K)=DOM(idom,K)+XferDOMK(idom,K)
    ENDDO
    DO  N=1,NumMicbHFunGrupsPerCmplx
      DO  M=1,nlbiomcp
        do NGL=JGniH(N),JGnfH(N)
          MID=micpar%get_micb_id(M,NGL)        
          DO NE=1,NumPlantChemElms
            mBiomeHeter(NE,MID,K)=AZERO(mBiomeHeter(NE,MID,K)+XferBiomeHeterK(NE,M,NGL,K))
          ENDDO
        enddo
      enddo
    enddo
  ENDDO D840
  call PrintInfo('end '//subname)
  end associate
  end subroutine OMTransferForPriming
!------------------------------------------------------------------------------------------

  subroutine RDOMSorption(KL,micfor,micstt,nmicf,ncplxf,ncplxs)
  !
  ! Description:
  !DOC ADSORPTION - DESORPTION  
  implicit none
  integer, intent(in) :: KL                 !from 1 to 5
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type),intent(inout) :: ncplxs
  real(r8) :: AECX
  real(r8) :: OQEX(idom_beg:idom_end)
  real(r8) :: OHEX(idom_beg:idom_end)
  integer :: idom,K
  real(r8) :: VLSoilPoreMicPX
  real(r8) :: VLSoilPoreMicPW,VOLCX,VOLCW,VOLAX,VOLAW
!     begin_execution
  associate(                                            &
    RDOMSorp             => ncplxf%RDOMSorp,            & !DOM transfer from solution to sorbed pools; positive adsorption, negative desorption
    BulkSOMC             => ncplxs%BulkSOMC,            & !Colonized solid, residue and sorbed C summed for each complex K; not referenced here
    FOCA                 => ncplxs%FOCA,                & !DOC fraction of DOC plus acetate in each substrate complex [-]
    FOAA                 => ncplxs%FOAA,                & !Acetate-C fraction of DOC plus acetate in each substrate complex [-]
    SoilMicPMassLayer    => micfor%SoilMicPMassLayer,   & !Soil mass associated with the current layer micropore domain
    ZERO                 => micfor%ZERO,                & !Small dimensionless or concentration threshold used by the routine
    ZEROS2               => micfor%ZEROS2,              & !Small water-volume threshold used by the routine
    ZEROS                => micfor%ZEROS,               & !Small mass or flux threshold used by the routine
    litrm                => micfor%litrm,               & !True for the surface litter layer
    VLWatMicPM           => micfor%VLWatMicPM,          & !Micropore water volume at each outer transport substep M
    FracBulkSOMC         => micstt%FracBulkSOMC,        & !Complex share of colonized solid, residue and sorbed C [-]
    DOM                  => micstt%DOM,                 & !Dissolved organic pools by species (DOC, DON, DOP, acetate) and complex K
    SorbedOM             => micstt%SorbedOM,            & !Sorbed organic pools by DOM species and complex K
    AEC                  => micfor%AEC                  & !Anion exchange capacity controlling DOM sorption
  )
  !TSORP=0._r8
!     VLWatMicPM=soil water content, FracBulkSOMC=fraction of total SOC
!     AEC,AECX=anion exchange capacity
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
!     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!     TSORP,HSORP=sorption rate constant and coefficient for OHC
!     FOCA,FOAA=fractions of DOC and acetate vs. DOC+acetate
!     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of OQC,acetate,DON,DOP
!
  DO K=1,KL
    IF(VLWatMicPM(NPH).GT.ZEROS2.AND.FracBulkSOMC(K).GT.ZERO)THEN
      IF(litrm)THEN
        AECX=0.5E+03_r8
      ELSE
        AECX=AEC
      ENDIF
      DO idom=idom_beg,idom_end
        !RedistDecompProduct already debited microbial uptake via SubstrateDribbling.
        !Use the remaining dissolved pool for DOC, DON, DOP, and acetate sorption.
        OQEX(idom) = AMAX1(ZEROS,DOM(idom,K))                       !free DOM
        OHEX(idom) = AMAX1(ZEROS,SorbedOM(idom,K))                   !adsorbed DOM
      ENDDO

      VLSoilPoreMicPX = SoilMicPMassLayer*AECX*HSORP*FracBulkSOMC(K)
      VLSoilPoreMicPW = VLWatMicPM(NPH)*FracBulkSOMC(K)
      IF(FOCA(K).GT.ZERO)THEN
        VOLCX                = FOCA(K)*VLSoilPoreMicPX
        VOLCW                = FOCA(K)*VLSoilPoreMicPW
        RDOMSorp(idom_doc,K) = TSORP*(OQEX(idom_doc)*VOLCX-OHEX(idom_doc)*VOLCW)/(VOLCX+VOLCW)
      ELSE
        RDOMSorp(idom_doc,K)=TSORP*(OQEX(idom_doc)*VLSoilPoreMicPX-OHEX(idom_doc)*VLSoilPoreMicPW)/(VLSoilPoreMicPX+VLSoilPoreMicPW)
      ENDIF

      IF(FOAA(K).GT.ZERO)THEN
        VOLAX                    = FOAA(K)*VLSoilPoreMicPX
        VOLAW                    = FOAA(K)*VLSoilPoreMicPW
        RDOMSorp(idom_acetate,K) = TSORP*(OQEX(idom_acetate)*VOLAX-OHEX(idom_acetate)*VOLAW)/(VOLAX+VOLAW) 
      ELSE
        RDOMSorp(idom_acetate,K)=TSORP*(OQEX(idom_acetate)*VLSoilPoreMicPX-OHEX(idom_acetate)*VLSoilPoreMicPW)/(VLSoilPoreMicPX+VLSoilPoreMicPW)
      ENDIF
      RDOMSorp(idom_don,K)=TSORP*(OQEX(idom_don)*VLSoilPoreMicPX-OHEX(idom_don)*VLSoilPoreMicPW)/(VLSoilPoreMicPX+VLSoilPoreMicPW)
      RDOMSorp(idom_dop,K)=TSORP*(OQEX(idom_dop)*VLSoilPoreMicPX-OHEX(idom_dop)*VLSoilPoreMicPW)/(VLSoilPoreMicPX+VLSoilPoreMicPW)
    ELSE
      DO idom=idom_beg,idom_end
        RDOMSorp(idom,K)=0.0_r8
      enddo
    ENDIF
    DO idom=idom_beg,idom_end
      DOM(idom,K)      = DOM(idom,K)-RDOMSorp(idom,K)
      SorbedOM(idom,K) = SorbedOM(idom,K)+RDOMSorp(idom,K)
    ENDDO
  ENDDO
  end associate
  end subroutine RDOMSorption
!------------------------------------------------------------------------------------------

  subroutine SolidOMDecomposition(I,J,KL,micfor,micstt,naqfdiag,nmicdiag,ncplxf,ncplxs,micflx)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), INTENT(INOUT) :: naqfdiag
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type),intent(inout) :: ncplxs
  type(micfluxtype), intent(inout) :: micflx  
  integer  :: M,NE,idom,K
  real(r8) :: CNOMX,CPOMX
  real(r8) :: COQCK,COSC
  real(r8) :: CPR,CNR
  real(r8) :: DCKD
  real(r8) :: DFNS
  real(r8) :: OQCI
  real(r8) :: RHOSCM,dHyd
  real(r8) :: FCNK(1:jcplx),FCPK(1:jcplx)

  character(len=*), parameter :: subname='SolidOMDecomposition'

!     begin_execution
  associate(                                               &
    tRHydlySOM           => micflx%tRHydlySOM,             & !Solid-substrate hydrolysis summed by element over components and complexes
    RHydlySOCK           => micflx%RHydlySOCK,             & !Total solid-organic-C hydrolysis in each complex K
    tRHydlyBioReSOM      => micflx%tRHydlyBioReSOM,        & !Microbial residue hydrolysis summed by element over compartments and complexes
    tRHydlySoprtOM       => micflx%tRHydlySoprtOM,         & !Sorbed-substrate hydrolysis summed by element over complexes
    RHydlysSolidOM       => ncplxf%RHydlysSolidOM,         & !Solid-substrate C/N/P hydrolysis by component and complex K
    RHumifySolidOM       => ncplxf%RHumifySolidOM,         & !Hydrolyzed solid C/N/P transferred to particulate organic matter
    RDecmpProdDOM        => ncplxf%RDecmpProdDOM,          & !Solid-substrate hydrolysis products entering DOM after humification
    RHydlysBioResduOM    => ncplxf%RHydlysBioResduOM,      & !Microbial residue C/N/P hydrolysis by compartment and complex K
    RHydlysSorptOM       => ncplxf%RHydlysSorptOM,         & !Sorbed-substrate hydrolysis by DOM species and complex K
    ROQC4HeterMicActCmpK => nmicdiag%ROQC4HeterMicActCmpK, & !Complex-total heterotrophic activity proxy driving organic-matter hydrolysis
    BulkSOMC             => ncplxs%BulkSOMC,               & !Colonized solid, residue and sorbed C summed for each complex K
    TOMEK                => nmicdiag%TOMEK,                & !Total active heterotrophic C/N/P in each substrate complex K
    tMaxNActMicrbK       => ncplxs%tMaxNActMicrbK,         & !Maximum N content of active heterotrophic biomass in each complex
    tMaxPActMicrbK       => ncplxs%tMaxPActMicrbK,         & !Maximum P content of active heterotrophic biomass in each complex
    VOLWZ                => nmicdiag%VOLWZ,                & !Effective water volume supporting microbial activity and rate constraints
    TSensGrowth          => nmicdiag%TSensGrowth,          & !Layer temperature response for microbial growth [-]; not referenced here
    RHydrolysisScalCmpK  => nmicdiag%RHydrolysisScalCmpK,  & !Microbial activity, density and product-inhibition scaling for hydrolysis
    iprotein             => micpar%iprotein,               & !Protein/first solid-component index
    icarbhyro            => micpar%icarbhyro,              & !Carbohydrate/second solid-component index
    icellulos            => micpar%icellulos,              & !Cellulose solid-component index
    ilignin              => micpar%ilignin,                & !Lignin solid-component index
    SPOSC                => micpar%SPOSC,                  & !Specific solid-substrate hydrolysis coefficients by component and complex K
    k_POM                => micpar%k_POM,                  & !Particulate-organic-matter complex index
    CNRH                 => micpar%CNRH,                   & !Target humified-material N:C ratios by receiving complex
    CPRH                 => micpar%CPRH,                   & !Target humified-material P:C ratios by receiving complex
    CDOM                 => ncplxs%CDOM,                   & !Dissolved substrate concentrations by DOM species and complex K
    EPOC                 => micstt%EPOC,                   & !Fraction controlling transfer of decomposed litter C to particulate organic matter
    CNOSC                => micstt%CNOSC,                  & !Solid substrate N:C ratios by component and complex K; not referenced here
    CPOSC                => micstt%CPOSC,                  & !Solid substrate P:C ratios by component and complex K; not referenced here
    SorbedOM             => micstt%SorbedOM,               & !Sorbed organic pools by DOM species and complex K
    SolidOM              => micstt%SolidOM,                & !Solid C/N/P pools by substrate component and complex K
    SolidOMAct           => micstt%SolidOMAct,             & !Microbially colonized solid C by substrate component and complex K
    OMBioResdu           => micstt%OMBioResdu,             & !Microbial residue C/N/P pools by residue compartment and complex K
    VLSoilMicP           => micfor%VLSoilMicP,             & !Bulk volume associated with the layer micropore domain
    ZEROS                => micfor%ZEROS,                  & !Small mass or flux threshold used by the routine
    ZEROS2               => micfor%ZEROS2,                 & !Small water-volume threshold used by the routine
    litrm                => micfor%litrm,                  & !True for the surface litter layer
    SoilMicPMassLayer    => micfor%SoilMicPMassLayer       & !Soil mass associated with the current layer micropore domain
  )
  call PrintInfo('beg '//subname)
!     FCPK=N,P limitation to microbial activity in each K
!     CNOMX,CPOMX=N:C,P:C ratios relative to set maximum values
!     COQCK=aqueous concentration of microbial activity
!     DCKD=Km for decomposition of SOC at current COQCK
!     DCKM0,DCKML=Km for decomposition of SOC at zero COQCK
!     DCKI=inhibition of decomposition by microbial concentration
!     BulkSOMC=total SOC
!     COSC=concentration of total SOC
!     SoilMicPMassLayer,VLSoilPoreMicP=mass, volume of soil layer
!     DFNS=effect of microbial concentration on decomposition
!     OQCI=DOC product inhibition for decomposition
!     OQKI=DOC product inhibition constant for decomposition
!

  DO K=1,KL

    IF(TOMEK(ielmc,K).GT.ZEROS)THEN
      CNOMX   = TOMEK(ielmn,K)/tMaxNActMicrbK(K)
      CPOMX   = TOMEK(ielmp,K)/tMaxPActMicrbK(K)
      FCNK(K) = AMIN1(1.0_r8,AMAX1(0.50_r8,CNOMX))
      FCPK(K) = AMIN1(1.0_r8,AMAX1(0.50_r8,CPOMX))
    ELSE
      FCNK(K) = 1.0_r8
      FCPK(K) = 1.0_r8
    ENDIF
  !
  !     AQUEOUS CONCENTRATION OF BIOMASS TO CACULATE INHIBITION
  !     CONSTANT FOR DECOMPOSITION
  !
    IF(VOLWZ.GT.ZEROS2)THEN
      COQCK=AMIN1(0.1E+06_r8,ROQC4HeterMicActCmpK(K)/VOLWZ)
    ELSE
      COQCK=0.1E+06_r8
    ENDIF
    IF(litrm)THEN
      DCKD=DCKM0*(1.0_r8+COQCK/DCKI)
    ELSE
      DCKD=DCKML*(1.0_r8+COQCK/DCKI)
    ENDIF
    
    IF(BulkSOMC(K).GT.ZEROS)THEN
      IF(SoilMicPMassLayer.GT.ZEROS)THEN
        COSC=BulkSOMC(K)/SoilMicPMassLayer
      ELSE
        COSC=BulkSOMC(K)/VLSoilMicP
      ENDIF
      DFNS = COSC/(COSC+DCKD)
      OQCI = 1.0_r8/(1.0_r8+CDOM(idom_doc,K)/OQKI)
  !
  !     C, N, P DECOMPOSITION RATE OF SOLID SUBSTRATES 'RDOS*' FROM
  !     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
  !     TEMPERATURE, SUBSTRATE C:N, C:P
  !
  !     CNS,CPS=N:C,P:C ratios of SOC
  !     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP
  !     OSA,OSN,OSP=active biomass C,N,P
  !     SPOSC=specific decomposition rate constant
  !     ROQC4HeterMicActCmpK=total respiration of DOC+DOA used to represent microbial activity
  !     DFNS=effect of microbial concentration on decomposition
  !     OQCI=DOC product inhibition for decomposition
  !     TSensGrowth=temperature stress effect
  !     BulkSOMC=total SOC
  !     FCNK,FCPK=N,P limitation to microbial activity in each K
  !
      RHydrolysisScalCmpK(K)=AZMAX1(ROQC4HeterMicActCmpK(K))*DFNS*OQCI/BulkSOMC(K)
      
      D785: DO M=1,jsken
        IF(SolidOM(ielmc,M,K).GT.ZEROS)THEN                  
          RHydlysSolidOM(ielmc,M,K) = SolidOMAct(M,K)*AMIN1(0.5_r8,SPOSC(M,K)*RHydrolysisScalCmpK(K))            
          dHyd=AMIN1(1._r8,RHydlysSolidOM(ielmc,M,K)/SolidOM(ielmc,M,K))      
          RHydlysSolidOM(ielmn,M,K) = AZERO(SolidOM(ielmn,M,K))*dHyd/FCNK(K)          
          RHydlysSolidOM(ielmp,M,K) = AZERO(SolidOM(ielmp,M,K))*dHyd/FCPK(K)
                
          DO NE=1,NumPlantChemElms
            RHydlysSolidOM(NE,M,K) = AMIN1(RHydlysSolidOM(NE,M,K),SolidOM(NE,M,K))
            tRHydlySOM(NE)         = tRHydlySOM(NE)+RHydlysSolidOM(NE,M,K)
          ENDDO
          RHydlySOCK(K)=RHydlySOCK(K)+RHydlysSolidOM(ielmc,M,K)
        ELSE

          DO NE    = 1, NumPlantChemElms
            RHydlysSolidOM(NE,M,K)=0.0_r8
          ENDDO
        ENDIF
      ENDDO D785
      call PrintInfo('end 785')      
  !
  !     HUMIFICATION OF DECOMPOSED RESIDUE LIGNIN WITH PROTEIN,
  !     CH2O AND CELLULOSE 'RHOS*' WITH REMAINDER 'RCOS*' TO DOC,DON,DOP
  !
  !     RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
  !     RDOSC,RDOSN,RDOSP=decomposition of SOC,SON,SOP
  !     CNRH,CPRH=N:C,P:C in POC
  !     EPOC=fraction of RDOSC allocated to POC from hour1.f
  !     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
  !   
      IF(K.LE.micpar%NumOfLitrCmplxs)THEN
        !litter complexes
        RHumifySolidOM(ielmc,ilignin,K)=AZMAX1(AMIN1(RHydlysSolidOM(ielmn,ilignin,K)/CNRH(k_POM) &
          ,RHydlysSolidOM(ielmp,ilignin,K)/CPRH(k_POM),EPOC*RHydlysSolidOM(ielmc,ilignin,K)))
        RHOSCM                           = 0.10_r8*RHumifySolidOM(ielmc,ilignin,K)
        RHumifySolidOM(ielmc,iprotein,K) = AZMAX1(AMIN1(RHydlysSolidOM(ielmc,iprotein,K) &
          ,RHydlysSolidOM(ielmn,iprotein,K)/CNRH(k_POM),RHydlysSolidOM(ielmp,iprotein,K)/CPRH(k_POM),RHOSCM))
        RHumifySolidOM(ielmc,icarbhyro,K)=AZMAX1(AMIN1(RHydlysSolidOM(ielmc,icarbhyro,K) &
          ,RHydlysSolidOM(ielmn,icarbhyro,K)/CNRH(k_POM),RHydlysSolidOM(ielmp,icarbhyro,K)/CPRH(k_POM),RHOSCM))
        RHumifySolidOM(ielmc,icellulos,K)=AZMAX1(AMIN1(RHydlysSolidOM(ielmc,icellulos,K) &
          ,RHydlysSolidOM(ielmn,icellulos,K)/CNRH(k_POM),RHydlysSolidOM(ielmp,icellulos,K)/CPRH(k_POM),RHOSCM-RHumifySolidOM(ielmc,icarbhyro,K)))

        D805: DO M=1,jsken
          RHumifySolidOM(ielmn,M,K) = AMIN1(RHydlysSolidOM(ielmn,M,K),RHumifySolidOM(ielmc,M,K)*CNRH(k_POM))
          RHumifySolidOM(ielmp,M,K) = AMIN1(RHydlysSolidOM(ielmp,M,K),RHumifySolidOM(ielmc,M,K)*CPRH(k_POM))
          DO NE=1,NumPlantChemElms
            RDecmpProdDOM(NE,M,K)=RHydlysSolidOM(NE,M,K)-RHumifySolidOM(NE,M,K)
          ENDDO
        ENDDO D805
      ELSE
        !non-litter complexes
        D810: DO M=1,jsken
          DO NE=1,NumPlantChemElms      
            RHumifySolidOM(NE,M,K) = 0.0_r8
            RDecmpProdDOM(NE,M,K)  = RHydlysSolidOM(NE,M,K)
          ENDDO
        ENDDO D810
      ENDIF
      call PrintInfo('end 805')      
    ELSE
      D780: DO M=1,jsken
        DO NE=1,NumPlantChemElms    
          RHydlysSolidOM(NE,M,K) = 0.0_r8
          RHumifySolidOM(NE,M,K) = 0.0_r8
          RDecmpProdDOM(NE,M,K)  = 0.0_r8
        ENDDO  
      ENDDO D780
    ENDIF
  !
  !     C, N, P DECOMPOSITION RATE OF BIORESIDUE 'RDOR*' FROM
  !     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
  !     TEMPERATURE, SUBSTRATE C:N, C:P
  !
  !     ORC,ORN,ORP=microbial residue C,N,P
  !     CNR,CPR=N:C,P:C ratios of microbial residue
  !     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
  !     SPORC=specific decomposition rate constant for microbial residue
  !     ROQC4HeterMicActCmpK=total respiration of DOC+DOA used to represent microbial activity
  !     DFNS=effect of microbial concentration on decomposition
  !     OQCI=DOC product inhibition for decomposition
  !     TSensGrowth=temperature stress effect
  !     BulkSOMC=total SOC
  !     FCNK,FCPK=N,P limitation to microbial activity in each K
  !
    IF(BulkSOMC(K).GT.ZEROS)THEN
      RHydrolysisScalCmpK(K)=AZMAX1(ROQC4HeterMicActCmpK(K))*DFNS*OQCI/BulkSOMC(K)
      D775: DO M=1,ndbiomcp
        IF(OMBioResdu(ielmc,M,K).GT.ZEROS)THEN
          dHyd=AMIN1(1._r8,SPORC(M)*RHydrolysisScalCmpK(K))
          RHydlysBioResduOM(ielmc,M,K) = AZERO(OMBioResdu(ielmc,M,K))*dHyd
          RHydlysBioResduOM(ielmn,M,K) = AZERO(OMBioResdu(ielmn,M,K))*dHyd/FCNK(K)
          RHydlysBioResduOM(ielmp,M,K) = AZERO(OMBioResdu(ielmp,M,K))*dHyd/FCPK(K)

          DO NE=1,NumPlantChemElms
            RHydlysBioResduOM(NE,M,K) = AMIN1(RHydlysBioResduOM(NE,M,K),OMBioResdu(NE,M,K))
            tRHydlyBioReSOM(NE)       = tRHydlyBioReSOM(NE)+RHydlysBioResduOM(NE,M,K)
          ENDDO
        ELSE
          DO NE=1,NumPlantChemElms      
            RHydlysBioResduOM(NE,M,K)=0.0_r8
          ENDDO  
        ENDIF
      ENDDO D775
      call PrintInfo('end 775')      

    ELSE
      D776: DO M=1,ndbiomcp
        DO NE=1,NumPlantChemElms    
          RHydlysBioResduOM(NE,M,K)=0.0_r8
        ENDDO  
      ENDDO D776
    ENDIF
  !
  !     C, N, P DECOMPOSITION RATE OF SORBED SUBSTRATES 'RDOH*' FROM
  !     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
  !     TEMPERATURE, SUBSTRATE C:N, C:P
  !
  !     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
  !     rCNSorbOM,rCPSorbOM=N:C,P:C ratios of adsorbed C,N,P
  !     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
  !     SPOHC=specific decomposition rate constant for adsorbed C
  !     ROQC4HeterMicActCmpK=total respiration of DOC+DOA used to represent microbial activity
  !     DFNS=effect of microbial concentration on decomposition
  !     OQCI=DOC product inhibition for decomposition
  !     TSensGrowth=temperature stress effect
  !     BulkSOMC=total SOC
  !     FCNK,FCPK=N,P limitation to microbial activity in each K
  !
    IF(BulkSOMC(K).GT.ZEROS)THEN
      RHydrolysisScalCmpK(K)=AZMAX1(ROQC4HeterMicActCmpK(K))*DFNS*OQCI/BulkSOMC(K)    
      IF(SorbedOM(ielmc,K).GT.ZEROS)THEN
        dHyd=AMIN1(1._r8, SPOHC*RHydrolysisScalCmpK(K))
        RHydlysSorptOM(ielmc,K)        = SorbedOM(ielmc,K)*dHyd
        RHydlysSorptOM(ielmn,K)        = AZERO(SorbedOM(ielmn,K))*dHyd/FCNK(K)
        RHydlysSorptOM(ielmp,K)        = AZERO(SorbedOM(ielmp,K))*dHyd/FCPK(K)

        DO NE=1,NumPlantChemElms
          RHydlysSorptOM(NE,K) = AMIN1(SorbedOM(NE,K),RHydlysSorptOM(NE,K))
          tRHydlySoprtOM(NE)   = tRHydlySoprtOM(NE)+RHydlysSorptOM(NE,K)
        ENDDO
      ELSE
        DO NE=1,NumPlantChemElms
          RHydlysSorptOM(NE,K)=0.0_r8
        ENDDO
      ENDIF
      !Acetate has its own donor and hydrolysis rate; it does not require
      !sorbed DOC. Retain its existing rate law without DOC product inhibition.
      RHydlysSorptOM(idom_acetate,K) = AZERO(SorbedOM(idom_acetate,K)) &
        *AZMAX1(AMIN1(1._r8,SPOHA*ROQC4HeterMicActCmpK(K)*DFNS/BulkSOMC(K)))
      !The elemental C diagnostic includes both DOC-C and acetate-C release.
      tRHydlySoprtOM(ielmc) = tRHydlySoprtOM(ielmc)+RHydlysSorptOM(idom_acetate,K)
    ELSE

      DO idom=idom_beg,idom_end          
        RHydlysSorptOM(idom,K)=0.0_r8
      ENDDO
    ENDIF
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine SolidOMDecomposition
!------------------------------------------------------------------------------------------

  subroutine RedistDecompProduct(micfor,KL,nmicdiag,nmicf,ncplxf,ncplxs,micstt)

  implicit none
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Diag_type), intent(in) :: nmicdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout):: ncplxs
  integer :: K,M,N,NGL,NE,idom,N1,N2
  integer :: ResidueFallbackComplex !Receiving complex when existing microbial residue is negligible
  real(r8) :: FORC(0:jcplx)
  real(r8) :: dflux,scal
!     begin_execution
  associate(                                                                    &
    k_POM                            => micpar%k_POM,                           & !Particulate-organic-matter complex index
    TDOMUptkHeter                    => ncplxf%TDOMUptkHeter,                   & !Total heterotrophic DOM uptake by species and complex K
    DOMUptk4GrothHeter               => nmicf%DOMUptk4GrothHeter,               & !Guild elemental uptake; C includes DOC, acetate and cyanobacterial CO2 fixation
    RMetabDOCUptkHeter               => nmicf%RMetabDOCUptkHeter,               & !Total DOC-C uptake for primary metabolism and denitrification
    RMetabAcetUptkHeter              => nmicf%RMetabAcetUptkHeter,              & !Total acetate-C uptake for primary metabolism and denitrification
    RkillLitrfal2ResduOMHeter        => nmicf%RkillLitrfal2ResduOMHeter,        & !Ordinary-mortality C/N/P routed to microbial residue from heterotrophic biomass; not referenced here
    RMaintDefLitrfal2ResduOMHeter    => nmicf%RMaintDefLitrfal2ResduOMHeter,    & !Starvation-derived C/N/P routed to microbial residue from heterotrophic biomass; not referenced here
    RCCMEAutor                       => nmicf%RCCMEAutor,                       & !Autotrophic mortality/starvation residue C/N/P routed to each receiving complex
    RAcetateProdHeter                => nmicf%RAcetateProdHeter,                & !Acetate-C production by heterotrophic guild and complex
    RCCMEHeter                       => nmicf%RCCMEHeter,                       & !Heterotrophic mortality/starvation residue C/N/P summed for each receiving complex
    RHydlysSolidOM                   => ncplxf%RHydlysSolidOM,                  & !Solid-substrate C/N/P hydrolysis by component and complex K
    RHumifySolidOM                   => ncplxf%RHumifySolidOM,                  & !Hydrolyzed solid C/N/P transferred to particulate organic matter
    RDecmpProdDOM                    => ncplxf%RDecmpProdDOM,                   & !Solid-substrate hydrolysis products entering DOM after humification
    RHydlysBioResduOM                => ncplxf%RHydlysBioResduOM,               & !Microbial residue C/N/P hydrolysis by compartment and complex K
    RHydlysSorptOM                   => ncplxf%RHydlysSorptOM,                  & !Sorbed-substrate hydrolysis by DOM species and complex K
    RDOMSorp                         => ncplxf%RDOMSorp,                        & !DOM transfer from solution to sorbed pools; positive adsorption, negative desorption; not referenced here
    OMBioResduK                      => ncplxs%OMBioResduK,                     & !Total microbial residue C in each complex K
    TOMBioResdu                      => nmicdiag%TOMBioResdu,                   & !Total microbial residue C summed over substrate complexes
    RMaintDefLitrfal2ResduOMAutor    => nmicf%RMaintDefLitrfal2ResduOMAutor,    & !Starvation-derived C/N/P routed to microbial residue from autotrophic biomass
    RkillLitrfal2ResduOMAutor        => nmicf%RkillLitrfal2ResduOMAutor,        & !Ordinary-mortality C/N/P routed to microbial residue from autotrophic biomass
    SolidOM                          => micstt%SolidOM,                         & !Solid C/N/P pools by substrate component and complex K
    iprotein                         => micpar%iprotein,                        & !Protein/first solid-component index
    SolidOMAct                       => micstt%SolidOMAct,                      & !Microbially colonized solid C by substrate component and complex K; not referenced here
    SOMPomProtein                    => micstt%SOMPomProtein,                   & !C/N/P transferred to the first POM component in underlying soil
    DOM                              => micstt%DOM,                             & !Dissolved organic pools by species (DOC, DON, DOP, acetate) and complex K
    DOM_MicP_drib                    => micstt%DOM_MicP_drib ,                  & !DOM pool correction recorded for downstream mass accounting
    OMBioResdu                       => micstt%OMBioResdu,                      & !Microbial residue C/N/P pools by residue compartment and complex K
    SorbedOM                         => micstt%SorbedOM,                        & !Sorbed organic pools by DOM species and complex K
    ZEROS                            => micfor%ZEROS,                           & !Small mass or flux threshold used by the routine
    NumHetetr1MicCmplx               => micpar%NumHetetr1MicCmplx ,             & !Number of heterotrophic guilds in one organic-matter complex
    is_activeMicrbFungrpHeter        => micpar%is_activeMicrbFungrpHeter,       & !Activation flags for heterotrophic functional groups
    Litrm                            => micfor%litrm                            & !True for the surface litter layer
  )
!
!     REDISTRIBUTE AUTOTROPHIC DECOMPOSITION PRODUCTS AMONG
!     HETEROTROPHIC SUBSTRATE-MICROBE complexES
!
!     FORC=fraction of total microbial residue
!     OMBioResduK=microbial residue
!     RCCMEheter,RCCMN,RCCMP=transfer of auto LitrFall C,N,P to each hetero K
!     RkillLitrfal2ResduOMHeter,RkillLitrfal2ResduOMHeter,RCOMP=transfer of microbial C,N,P LitrFall to residue
!     RMaintDefLitrfal2ResduOMHeter,RCMMN,RMaintDefLitrfal2ResduOMHeter=transfer of senesence LitrFall C,N,P to residue
!
  !Surface litter has no POM complex: seed its fine-litter residue pool.
  !Mineral soil retains POM as the fallback; both destinations lie in 1:KL.
  ResidueFallbackComplex=k_POM
  IF(litrm)ResidueFallbackComplex=micpar%k_fine_comp
  D1690: DO K=1,KL
    IF(TOMBioResdu.GT.ZEROS)THEN
      !partition dead autotrophs proportionally
      FORC(K)=OMBioResduK(K)/TOMBioResdu
    ELSE
      IF(K.EQ.ResidueFallbackComplex)THEN
        FORC(K)=1.0_r8
      ELSE
        FORC(K)=0.0_r8
      ENDIF
    ENDIF
    D1685: DO N=1,NumMicbAFunGrupsPerCmplx
      if(.not.micpar%is_activeMicrbFungrpAutor(N))cycle
      D1680: DO M=1,ndbiomcp
        DO NGL=JGniA(N),JGnfA(N)
          DO NE=1,NumPlantChemElms
            RCCMEAutor(NE,M,K)=RCCMEAutor(NE,M,K)+(RkillLitrfal2ResduOMAutor(NE,M,NGL)+RMaintDefLitrfal2ResduOMAutor(NE,M,NGL))*FORC(K)
          ENDDO
        ENDDO
      ENDDO D1680
    ENDDO D1685
  ENDDO D1690
!
!   REDISTRIBUTE C,N AND P TRANSFOMBioResduATIONS AMONG STATE
!   VARIABLES IN SUBSTRATE-MICROBE complexES
!

  D590: DO K=1,KL
    D580: DO M=1,jsken
      !
      !     SUBSTRATE DECOMPOSITION PRODUCTS
      !
      !     OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
      !     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP
      !     OQC,OQN,OQP,OQA=DOC,DON,DOP
      !     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
      !
      DO NE=1,NumPlantChemElms
        SolidOM(NE,M,K)=SolidOM(NE,M,K)-RHydlysSolidOM(NE,M,K)
      ENDDO

!     SolidOMAct(M,K)=SolidOMAct(M,K)-RHydlysSolidOM(ielmc,M,K)
      DO NE=1,NumPlantChemElms
        DOM(NE,K)=DOM(NE,K)+RDecmpProdDOM(NE,M,K)
      ENDDO
      !
      !     LIGNIFICATION PRODUCTS
      !
      !       RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
      !
      IF(.not.litrm)THEN
       ! add to POM carbonhydrate from humification
 !      SolidOMAct(1,k_POM)=SolidOMAct(1,k_POM)+RHumifySolidOM(ielmc,M,K)
        DO NE=1,NumPlantChemElms
          SolidOM(NE,iprotein,k_POM)=SolidOM(NE,iprotein,k_POM)+RHumifySolidOM(NE,M,K)
        ENDDO
      ELSE
        DO NE=1,NumPlantChemElms
          SOMPomProtein(NE)=SOMPomProtein(NE)+RHumifySolidOM(NE,M,K)
        ENDDO
      ENDIF
    ENDDO D580
    !
    !     MICROBIAL RESIDUE DECOMPOSITION PRODUCTS
    !
    !     ORC,ORN,ORP=microbial residue C,N,P
    !     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
    !     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
    !
    D575: DO M=1,ndbiomcp
      DO NE=1,NumPlantChemElms    
        OMBioResdu(NE,M,K) = OMBioResdu(NE,M,K)-RHydlysBioResduOM(NE,M,K)
        DOM(NE,K)          = DOM(NE,K)+RHydlysBioResduOM(NE,M,K)
      ENDDO
    ENDDO D575

    N1=1; N2=NumHetetr1MicCmplx
    DO  N=1,NumMicbHFunGrupsPerCmplx
      if(micpar%is_anaerobic_hetr(N))then
        DO NGL=JGniH(N),JGnfH(N)
          DOM(idom_acetate,K) = DOM(idom_acetate,K)+RAcetateProdHeter(NGL,K)
        ENDDO
      endif
    ENDDO

!    write(151,*)K,'b',DOM(idom_doc,K)+SorbedOM(idom_doc,K)+DOM(idom_acetate,K)+SorbedOM(idom_acetate,K),&
!      DOM(idom_doc,K),SorbedOM(idom_doc,K),DOM(idom_acetate,K),SorbedOM(idom_acetate,K)    
    !update adsorbed OM    
    DO idom=idom_beg,idom_end
      dflux            = RHydlysSorptOM(idom,K)
      DOM(idom,K)      = DOM(idom,K)+dflux
      SorbedOM(idom,K) = SorbedOM(idom,K)-dflux
    ENDDO
!    write(151,*)K,'e',DOM(idom_doc,K)+SorbedOM(idom_doc,K)+DOM(idom_acetate,K)+SorbedOM(idom_acetate,K), &
!      DOM(idom_doc,K),SorbedOM(idom_doc,K),DOM(idom_acetate,K),SorbedOM(idom_acetate,K)    
    !
    !     MICROBIAL UPTAKE OF DISSOLVED C, N, P
    !
    !     RMetabDOCUptkHeter,RMetabAcetUptkHeter,DOMUptk4GrothHeter,DOMUptk4GrothHeter=DOC,acetate,DON,DOP uptake
    !     RAcetateProdHeter=acetate production from fermentation
    
    call SubstrateDribbling(n1,n2,RMetabDOCUptkHeter(N1:N2,K),DOM_MicP_drib(idom_doc,K),DOM(idom_doc,K))

    call SubstrateDribbling(n1,n2,DOMUptk4GrothHeter(ielmn,N1:N2,K),DOM_MicP_drib(idom_don,K),DOM(idom_don,K))

    call SubstrateDribbling(n1,n2,DOMUptk4GrothHeter(ielmp,N1:N2,K),DOM_MicP_drib(idom_dop,K),DOM(idom_dop,K))

    call SubstrateDribbling(n1,n2,RMetabAcetUptkHeter(N1:N2,K),DOM_MicP_drib(idom_acetate,K),DOM(idom_acetate,K))
!   
!     MICROBIAL DECOMPOSITION PRODUCTS
!
    D570: DO N=1,NumMicbHFunGrupsPerCmplx
      if(.not.is_activeMicrbFungrpHeter(N))cycle
      DO NGL=JGniH(N),JGnfH(N)
        TDOMUptkHeter(idom_doc,K)     = TDOMUptkHeter(idom_doc,K)+RMetabDOCUptkHeter(NGL,K)
        TDOMUptkHeter(idom_don,K)     = TDOMUptkHeter(idom_don,K)+DOMUptk4GrothHeter(ielmn,NGL,K)
        TDOMUptkHeter(idom_dop,K)     = TDOMUptkHeter(idom_dop,K)+DOMUptk4GrothHeter(ielmp,NGL,K)
        TDOMUptkHeter(idom_acetate,K) = TDOMUptkHeter(idom_acetate,K)+RMetabAcetUptkHeter(NGL,K)        
      enddo
    ENDDO D570

    !These totals already include all contributing guilds; apply once per complex.
    D565: DO M=1,ndbiomcp
      DO NE=1,NumPlantChemElms
        OMBioResdu(NE,M,K)=OMBioResdu(NE,M,K)+RCCMEHeter(NE,M,K)+RCCMEAutor(NE,M,K)
      ENDDO
    ENDDO D565

  ENDDO D590
  end associate
  end subroutine RedistDecompProduct
!------------------------------------------------------------------------------------------

  subroutine HeterotrophAnabolicUpdate(I,J,micfor,micstt,nmicf,micflx)
  implicit none
  integer, intent(in) :: I,J
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx
  character(len=*), parameter :: subname='HeterotrophAnabolicUpdate'
  integer  :: K,M,N,NGL,MID3,MID,NE
  real(r8) ::CGROMC,dmassC
!     begin_execution
  real(r8) :: ReserveSupply,MineralTransfer(4),PreviousMineralTransfer
  associate(                                                                &
    DOMUptk4GrothHeter             => nmicf%DOMUptk4GrothHeter,             & !Guild elemental uptake; C includes DOC, acetate and cyanobacterial CO2 fixation
    NonstX2stBiomHeter             => nmicf%NonstX2stBiomHeter,             & !C/N/P transfer from reserves into kinetic/structural biomass by guild and complex
    RespGrossHeter                 => nmicf%RespGrossHeter,                 & !Gross respiration C equivalent from the primary heterotrophic pathway
    RNOxDOCReduxRespDenitLim       => nmicf%RNOxDOCReduxRespDenitLim,       & !DOC-C respiration supported by realized denitrification
    RNOxAcetReduxRespDenitLim      => nmicf%RNOxAcetReduxRespDenitLim,      & !Acetate-C respiration supported by realized denitrification
    RNO3imobilSoilHeter            => nmicf%RNO3imobilSoilHeter,            & !Net NO3-N transfer from nonband soil to microbes; positive immobilization
    RCO2ProdHeter                  => nmicf%RCO2ProdHeter,                  & !CO2-C production by heterotrophic guild and complex
    RH2PO4imobilSoilHeter          => nmicf%RH2PO4imobilSoilHeter,          & !Net H2PO4-P transfer from nonband soil to microbes; positive immobilization
    RNH4imobilBandHeter            => nmicf%RNH4imobilBandHeter,            & !Net NH4-N transfer from fertilizer-band soil to microbes; positive immobilization
    RNO3imobilBandHeter            => nmicf%RNO3imobilBandHeter,            & !Net NO3-N transfer from fertilizer-band soil to microbes; positive immobilization
    RH2PO4imobilBandHeter          => nmicf%RH2PO4imobilBandHeter,          & !Net H2PO4-P transfer from fertilizer-band soil to microbes; positive immobilization
    RkillLitrfal2HumOMHeter        => nmicf%RkillLitrfal2HumOMHeter,        & !Ordinary-mortality C/N/P routed to humified material from heterotrophic biomass
    RMaintDefLitrfal2HumOMHeter    => nmicf%RMaintDefLitrfal2HumOMHeter,    & !Starvation-derived C/N/P routed to humified material from heterotrophic biomass
    RN2FixHeter                    => nmicf%RN2FixHeter,                    & !Heterotrophic guild N2-N fixation into biomass
    RKillOMHeter                   => nmicf%RKillOMHeter,                   & !Ordinary mortality C/N/P withdrawal from heterotrophic biomass
    RkillRecycOMHeter              => nmicf%RkillRecycOMHeter,              & !Ordinary-mortality C/N/P recycled to reserves from heterotrophic biomass
    RMaintDefcitKillOMHeter        => nmicf%RMaintDefcitKillOMHeter,        & !Maintenance-starvation C/N/P withdrawal from heterotrophic biomass
    RMaintDefcitRecycOMHeter       => nmicf%RMaintDefcitRecycOMHeter,       & !Starvation recycling: C respired, N/P returned to reserves from heterotrophic biomass
    RNH4imobilLitrHeter            => nmicf%RNH4imobilLitrHeter,            & !Net NH4-N transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RNO3imobilLitrHeter            => nmicf%RNO3imobilLitrHeter,            & !Net NO3-N transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RH2PO4imobilLitrHeter          => nmicf%RH2PO4imobilLitrHeter,          & !Net H2PO4-P transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RH1PO4imobilSoilHeter          => nmicf%RH1PO4imobilSoilHeter,          & !Net HPO4-P transfer from nonband soil to microbes; positive immobilization
    RH1PO4imobilBandHeter          => nmicf%RH1PO4imobilBandHeter,          & !Net HPO4-P transfer from fertilizer-band soil to microbes; positive immobilization
    RH1PO4imobilLitrHeter          => nmicf%RH1PO4imobilLitrHeter,          & !Net HPO4-P transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RNH4imobilSoilHeter            => nmicf%RNH4imobilSoilHeter,            & !Net NH4-N transfer from nonband soil to microbes; positive immobilization
    k_POM                          => micpar%k_POM,                         & !Particulate-organic-matter complex index
    k_humus                        => micpar%k_humus,                       & !Humus complex index
    mBiomeHeter                    => micstt%mBiomeHeter,                   & !C/N/P pools indexed by element, flattened guild/compartment and complex K
    SolidOM                        => micstt%SolidOM,                       & !Solid C/N/P pools by substrate component and complex K
    SOMHumProtein                  => micstt%SOMHumProtein,                 & !C/N/P transferred to the first humus component in underlying soil
    SOMHumCarbohyd                 => micstt%SOMHumCarbohyd,                & !C/N/P transferred to the second humus component in underlying soil
    ElmAllocmatMicrblitr2POM       => micfor%ElmAllocmatMicrblitr2POM,      & !Partition of humified microbial litter into receiving solid components
    ElmAllocmatMicrblitr2POMU      => micfor%ElmAllocmatMicrblitr2POMU,     & !Underlying-soil partition of humified litter from surface microbes
    icarbhyro                      => micpar%icarbhyro,                     & !Carbohydrate/second solid-component index
    iprotein                       => micpar%iprotein,                      & !Protein/first solid-component index
    Litrm                          => micfor%litrm,                         & !True for the surface litter layer
    NetCAssimhr                    => micflx%NetCAssimhr,                   & !Accumulated net heterotrophic C assimilation after metabolic costs
    GrosAssimhr                    => micflx%GrosAssimhr,                   & !Accumulated gross heterotrophic substrate-C uptake
    NetNH4Mineralize               => micflx%NetNH4Mineralize,              & !Net mineral N exchange (NH4 plus NO3); positive immobilization, negative release
    NetPO4Mineralize               => micflx%NetPO4Mineralize               & !Net phosphate exchange; positive immobilization, negative mineralization
  )
!
!     OMC,OMN,OMP=microbial C,N,P
!     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
!     RXOMC,RXOMN,RKillOMHeter=microbial C,N,P decomposition
!     RMaintDefcitKillOMHeter,RXMMN,RXMMP=microbial C,N,P loss from senescence
!

  D550: DO K=1,jcplx
    IF(.not.litrm .OR. (K.NE.k_POM .AND. K.NE.k_humus))THEN
      DO  N=1,NumMicbHFunGrupsPerCmplx
        if(.not.micpar%is_activeMicrbFungrpHeter(N))cycle
        DO NGL=JGniH(N),JGnfH(N)
          !Reconcile reserve withdrawals before crediting structural biomass.
          MID3=micpar%get_micb_id(iLbiom_reserve,NGL)
          DO NE=ielmn,ielmp
            ReserveSupply=DOMuptk4GrothHeter(NE,NGL,K) &
              +SUM(RkillRecycOMHeter(NE,1:2,NGL,K))+SUM(RMaintDefcitRecycOMHeter(NE,1:2,NGL,K))
            IF(NE.EQ.ielmn)THEN
              ReserveSupply=ReserveSupply+RN2FixHeter(NGL,K)
              IF(litrm)ReserveSupply=ReserveSupply+RNH4imobilLitrHeter(NGL,K)+RNO3imobilLitrHeter(NGL,K)
              MineralTransfer=[RNH4imobilSoilHeter(NGL,K),RNH4imobilBandHeter(NGL,K),RNO3imobilSoilHeter(NGL,K),RNO3imobilBandHeter(NGL,K)]
            ELSE
              IF(litrm)ReserveSupply=ReserveSupply+RH2PO4imobilLitrHeter(NGL,K)+RH1PO4imobilLitrHeter(NGL,K)
              MineralTransfer=[RH2PO4imobilSoilHeter(NGL,K),RH2PO4imobilBandHeter(NGL,K),RH1PO4imobilSoilHeter(NGL,K),RH1PO4imobilBandHeter(NGL,K)]
            ENDIF
            PreviousMineralTransfer=SUM(MineralTransfer)
            call LimitReserveNutrientTransfers(mBiomeHeter(NE,MID3,K),ReserveSupply, &
              NonstX2stBiomHeter(NE,1:2,NGL,K),MineralTransfer)
            IF(NE.EQ.ielmn)THEN
              RNH4imobilSoilHeter(NGL,K)=MineralTransfer(1)
              RNH4imobilBandHeter(NGL,K)=MineralTransfer(2)
              RNO3imobilSoilHeter(NGL,K)=MineralTransfer(3)
              RNO3imobilBandHeter(NGL,K)=MineralTransfer(4)
              micflx%NetNH4Mineralize=micflx%NetNH4Mineralize+SUM(MineralTransfer)-PreviousMineralTransfer
            ELSE
              RH2PO4imobilSoilHeter(NGL,K)=MineralTransfer(1)
              RH2PO4imobilBandHeter(NGL,K)=MineralTransfer(2)
              RH1PO4imobilSoilHeter(NGL,K)=MineralTransfer(3)
              RH1PO4imobilBandHeter(NGL,K)=MineralTransfer(4)
              micflx%NetPO4Mineralize=micflx%NetPO4Mineralize+SUM(MineralTransfer)-PreviousMineralTransfer
            ENDIF
          ENDDO

          D540: DO M=1,2
            MID=micpar%get_micb_id(M,NGL)     
            DO NE=1,NumPlantChemElms     
              mBiomeHeter(NE,MID,K)=mBiomeHeter(NE,MID,K)+NonstX2stBiomHeter(NE,M,NGL,K) &
                -RKillOMHeter(NE,M,NGL,K)-RMaintDefcitKillOMHeter(NE,M,NGL,K)
            ENDDO  
!
!     HUMIFICATION PRODUCTS
!
!     ElmAllocmatMicrblitr2POM=fractions allocated to humic vs fulvic humus
!     RHOMC,RHOMN,RkillLitrfal2HumOMHeter=transfer of microbial C,N,P LitrFall to humus
!     RMaintDefLitrfal2HumOMHeter,RHMMN,RMaintDefLitrfal2HumOMHeter=transfer of senesence LitrFall C,N,P to humus
!
            IF(.not.litrm)THEN              
              DO NE=1,NumPlantChemElms
                !add as protein
                SolidOM(NE,iprotein,k_humus)=SolidOM(NE,iprotein,k_humus) &
                  +ElmAllocmatMicrblitr2POM(iprotein)*(RkillLitrfal2HumOMHeter(NE,M,NGL,K)&
                  +RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K))
                !add as carbon hydro
                SolidOM(NE,icarbhyro,k_humus)=SolidOM(NE,icarbhyro,k_humus) &
                  +ElmAllocmatMicrblitr2POM(icarbhyro)*(RkillLitrfal2HumOMHeter(NE,M,NGL,K)&
                  +RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K))
              ENDDO
            ELSE
              DO NE=1,NumPlantChemElms
                SOMHumProtein(NE)=SOMHumProtein(NE)+ElmAllocmatMicrblitr2POMU(iprotein) &
                  *(RkillLitrfal2HumOMHeter(NE,M,NGL,K)+RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K))
                SOMHumCarbohyd(NE)=SOMHumCarbohyd(NE)+ElmAllocmatMicrblitr2POMU(icarbhyro) &
                  *(RkillLitrfal2HumOMHeter(NE,M,NGL,K)+RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K))                  
              ENDDO
            ENDIF
          ENDDO D540
          !
          !     INPUTS TO NONSTRUCTURAL POOLS
          !
          !     DOMUptk4GrothHeter=total DOC+acetate uptake
          !     RespGrossHeter=total respiration, including produciton of CO2 and acetate
          !     RNOxDOCReduxRespDenitLim=DOC respiration for denitrifcation
          !     Resp4NFixHeter=respiration for N2 fixation
          !     RCO2ProdHeter=total CO2 emission
          !     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
          !     R3OMC,R3OMN,RkillRecycOMHeter=microbial C,N,P recycling
          !     RMaintDefcitRecycOMHeter,R3MMN,R3MMP=microbial C,N,P recycling from senescence
          !     DOMUptk4GrothHeter,DOMUptk4GrothHeter=DON, DOP uptake
          !     RNH4imobilSoilHeter,RNH4imobilBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
          !     RNO3imobilSoilHeter,RNO3imobilBandHeter=substrate-limited NO3 immobiln in non-band, band
          !     RH2PO4imobilSoilHeter,RH2PO4imobilBandHeter=substrate-limited H2PO4 mineraln-immobn in non-band, band
          !     RH1PO4imobilSoilHeter,RH1PO4imobilBandHeter=substrate-limited HPO4 mineraln-immobn in non-band, band
          !     RNH4imobilLitrHeter,RNO3imobilLitrHeter =substrate-limited NH4,NO3 mineraln-immobiln
          !     RH2PO4imobilLitrHeter,RH1PO4imobilLitrHeter=substrate-limited H2PO4,HPO4 mineraln-immobiln
          !
          !Both denitrification substrates contribute uptake and respiratory C loss.
          CGROMC = DOMUptk4GrothHeter(ielmc,NGL,K)-RespGrossHeter(NGL,K) &
            -RNOxDOCReduxRespDenitLim(NGL,K)-RNOxAcetReduxRespDenitLim(NGL,K)
          !
          !N2 fixation redirects part of existing growth respiration; its
          !cost already reduces growth-associated uptake in the anabolic flux.
          !RespGrossHeter and catabolic products include that respiration.
          !Cyanobacterial extra photosynthate respiration is also included
          !upstream with its matching O2 flux. Do not debit or emit it again.
          !
          NetCAssimhr = NetCAssimhr+CGROMC
          GrosAssimhr = GrosAssimhr+DOMUptk4GrothHeter(ielmc,NGL,K)
          MID3        = micpar%get_micb_id(iLbiom_reserve,NGL)

          D555: DO M = 1, 2
            DO NE=1,NumPlantChemElms
              mBiomeHeter(NE,MID3,K)=mBiomeHeter(NE,MID3,K)-NonstX2stBiomHeter(NE,M,NGL,K)+RkillRecycOMHeter(NE,M,NGL,K)
            ENDDO
            DO NE=2,NumPlantChemElms
              mBiomeHeter(NE,MID3,K)=mBiomeHeter(NE,MID3,K)+RMaintDefcitRecycOMHeter(NE,M,NGL,K)
            ENDDO
            !TODO: Resolve starvation respiration with pathway-specific redox budgets.
            !This CO2 release currently has no matching electron-acceptor debit.
            !Aerobic respiration needs available O2 and consistent demand/uptake;
            !anaerobic groups need their own substrates and reduced products,
            !not an O2 charge. Preserve unrespired C in organic pools and avoid
            !double-counting maintenance already supported by catabolism.
            !See the corresponding starvation update in MicAutoCplxFGMod.
            !respire
            RCO2ProdHeter(NGL,K)=RCO2ProdHeter(NGL,K)+RMaintDefcitRecycOMHeter(ielmc,M,NGL,K)
          ENDDO D555

          mBiomeHeter(ielmc,MID3,K) = AZERO(mBiomeHeter(ielmc,MID3,K)+CGROMC)
          mBiomeHeter(ielmn,MID3,K) = mBiomeHeter(ielmn,MID3,K)+DOMUptk4GrothHeter(ielmn,NGL,K) &
            +RNH4imobilSoilHeter(NGL,K)+RNH4imobilBandHeter(NGL,K)+RNO3imobilSoilHeter(NGL,K) &
            +RNO3imobilBandHeter(NGL,K)+RN2FixHeter(NGL,K)
          mBiomeHeter(ielmp,MID3,K)=mBiomeHeter(ielmp,MID3,K)+DOMUptk4GrothHeter(ielmp,NGL,K) &
            +RH2PO4imobilSoilHeter(NGL,K)+RH2PO4imobilBandHeter(NGL,K)+RH1PO4imobilSoilHeter(NGL,K) &
            +RH1PO4imobilBandHeter(NGL,K)
         
          IF(litrm)THEN
            mBiomeHeter(ielmn,MID3,K)=mBiomeHeter(ielmn,MID3,K)+RNH4imobilLitrHeter(NGL,K)+RNO3imobilLitrHeter(NGL,K)
            mBiomeHeter(ielmp,MID3,K)=mBiomeHeter(ielmp,MID3,K)+RH2PO4imobilLitrHeter(NGL,K)+RH1PO4imobilLitrHeter(NGL,K)
          ENDIF
        enddo
      ENDDO
    ENDIF
  ENDDO D550
  end associate
  end subroutine HeterotrophAnabolicUpdate
!------------------------------------------------------------------------------------------

  subroutine MicrobialLitterColonization(I,J,KL,micfor,micstt,ncplxf,ncplxs,nmicdiag)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type),intent(inout) :: ncplxs
  type(Microbe_Diag_type),intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='MicrobialLitterColonization'

  integer  :: K,M
  real(r8) :: DOSAK

!     begin_execution
  associate(                                               &
    ROQC4HeterMicActCmpK => nmicdiag%ROQC4HeterMicActCmpK, & !Complex-total heterotrophic activity proxy driving organic-matter hydrolysis
    SolidOMCK            => ncplxs%SolidOMCK,              & !Total solid organic C in each complex K
    SolidOMActK          => ncplxs%SolidOMActK,            & !Total microbially colonized solid C in each complex K
    ZEROS                => micfor%ZEROS,                  & !Small mass or flux threshold used by the routine
    SolidOMAct           => micstt%SolidOMAct,             & !Microbially colonized solid C by substrate component and complex K
    SolidOM              => micstt%SolidOM,                & !Solid C/N/P pools by substrate component and complex K
    DOSA                 => micpar%DOSA                    & !Rate coefficient for microbial colonization of solid substrate
  )

  call PrintInfo('beg '//subname)
!     SolidOMCK,SolidOMActK,OSCX=total,colonized,uncolonized SOC
!     OSA,OSC=colonized,total litter
!     DOSA=rate constant for litter colonization
!     ROQC4HeterMicActCmpK=total respiration of DOC+DOA used to represent microbial activity
!
  D475: DO K=1,KL
    SolidOMCK(K)   = 0.0_r8
    SolidOMActK(K) = 0.0_r8
    DO  M=1,jsken
      SolidOMCK(K)    = SolidOMCK(K)+SolidOM(ielmc,M,K)
      SolidOMActK(K) = SolidOMActK(K)+SolidOMAct(M,K)
    enddo
  ENDDO D475

  D480: DO K=1,KL
    IF(SolidOMCK(K).GT.ZEROS)THEN
      DOSAK=DOSA(K)*AZMAX1(ROQC4HeterMicActCmpK(K))

      D485: DO M=1,jsken
        SolidOMAct(M,K)=AMIN1(SolidOM(ielmc,M,K),SolidOMAct(M,K)+DOSAK*SolidOM(ielmc,M,K)/SolidOMCK(K))
      ENDDO D485
    ELSE
      D490: DO M=1,jsken
        SolidOMAct(M,K)=AMIN1(SolidOM(ielmc,M,K),SolidOMAct(M,K))
      ENDDO D490
    ENDIF
  ENDDO D480
  call PrintInfo('end '//subname)
  end associate
  end subroutine MicrobialLitterColonization
!------------------------------------------------------------------------------------------

  subroutine AggregateTransfOMBioResdue(KL,micfor,micstt,nmicdiag,naqfdiag,nmicf,ncplxf,micflx)
  implicit none
  integer,intent(in) :: KL
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_Diag_type),intent(in) :: nmicdiag
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(micfluxtype), intent(inout) :: micflx

  integer  :: K,M,N,NGL,NE,idom
  real(r8) :: NC10NO2ReduxSoil, NC10NO2ReduxBand
  real(r8) :: NH4toN2OSoil, NH4toN2OBand
  real(r8) :: O2RespUptake,O2RespDemand
  real(r8) :: DONChemoPerComplex
!     begin_execution
  NC10NO2ReduxSoil = 0._r8
  NC10NO2ReduxBand = 0._r8
  NH4toN2OSoil     = 0._r8
  NH4toN2OBand     = 0._r8
  associate(                                                       &
    DOMUptk4GrothHeter        => nmicf%DOMUptk4GrothHeter,         & !Guild elemental uptake; C includes DOC, acetate and cyanobacterial CO2 fixation
    RMetabDOCUptkHeter        => nmicf%RMetabDOCUptkHeter,         & !Total DOC-C uptake for primary metabolism and denitrification
    RMetabAcetUptkHeter       => nmicf%RMetabAcetUptkHeter,        & !Total acetate-C uptake for primary metabolism and denitrification
    RO2UptkHeter              => nmicf%RO2UptkHeter,               & !Realized total O2 uptake by heterotrophic guild and complex
    RO2Uptk4RespHeter         => nmicf%RO2Uptk4RespHeter,         & !Gross respiratory O2 consumption, including photosynthetic supply
    RO2Dmnd4RespHeter         => nmicf%RO2Dmnd4RespHeter,         & !Gross potential respiratory O2 demand
    RO2DmndHeter              => nmicf%RO2DmndHeter,               & !Total guild O2 demand before O2 limitation
    RNO3ReduxHeterSoil        => nmicf%RNO3ReduxHeterSoil,         & !Realized nonband NO3-N reduction to NO2-N by heterotrophic guilds
    RNO3ReduxHeterBand        => nmicf%RNO3ReduxHeterBand,         & !Realized fertilizer-band NO3-N reduction to NO2-N by heterotrophic guilds
    RNO2ReduxHeterSoil        => nmicf%RNO2ReduxHeterSoil,         & !Realized nonband NO2-N reduction to N2O-N by heterotrophic guilds
    RNO2ReduxHeterBand        => nmicf%RNO2ReduxHeterBand,         & !Realized fertilizer-band NO2-N reduction to N2O-N by heterotrophic guilds
    RN2OReduxHeter            => nmicf%RN2OReduxHeter,             & !Realized N2O-N reduction to N2 by guild and complex
    RNOxDOCReduxRespDenitLim  => nmicf%RNOxDOCReduxRespDenitLim,   & !DOC-C respiration supported by realized denitrification
    RNOxAcetReduxRespDenitLim => nmicf%RNOxAcetReduxRespDenitLim,  & !Acetate-C respiration supported by realized denitrification
    RNH4imobilSoilHeter       => nmicf%RNH4imobilSoilHeter,        & !Net NH4-N transfer from nonband soil to microbes; positive immobilization
    RNO3imobilSoilHeter       => nmicf%RNO3imobilSoilHeter,        & !Net NO3-N transfer from nonband soil to microbes; positive immobilization
    RH2PO4imobilSoilHeter     => nmicf%RH2PO4imobilSoilHeter,      & !Net H2PO4-P transfer from nonband soil to microbes; positive immobilization
    RNH4imobilBandHeter       => nmicf%RNH4imobilBandHeter,        & !Net NH4-N transfer from fertilizer-band soil to microbes; positive immobilization
    RNO3imobilBandHeter       => nmicf%RNO3imobilBandHeter,        & !Net NO3-N transfer from fertilizer-band soil to microbes; positive immobilization
    RH2PO4imobilBandHeter     => nmicf%RH2PO4imobilBandHeter,      & !Net H2PO4-P transfer from fertilizer-band soil to microbes; positive immobilization
    RH2ProdHeter              => nmicf%RH2ProdHeter,               & !H2 production by heterotrophic guild and complex
    RNH4imobilLitrHeter       => nmicf%RNH4imobilLitrHeter,        & !Net NH4-N transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RNO3imobilLitrHeter       => nmicf%RNO3imobilLitrHeter,        & !Net NO3-N transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RH2PO4imobilLitrHeter     => nmicf%RH2PO4imobilLitrHeter,      & !Net H2PO4-P transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RCO2ProdHeter             => nmicf%RCO2ProdHeter,              & !CO2-C production by heterotrophic guild and complex
    RAcetateProdHeter         => nmicf%RAcetateProdHeter,          & !Acetate-C production by heterotrophic guild and complex
    RCH4ProdHeter             => nmicf%RCH4ProdHeter,              & !CH4-C production by heterotrophic guild and complex
    TRDOM2DIE                 => micflx%TRDOM2DIE,                 & !Net conversion from organic to inorganic C/N/P for balance accounting
    RSMetaOxidSoilAutor       => nmicf%RSMetaOxidSoilAutor,        & !Nonband catabolic substrate oxidation; substrate depends on functional group
    RSMetaOxidBandAutor       => nmicf%RSMetaOxidBandAutor,        & !Fertilizer-band catabolic substrate oxidation by autotrophic guild
    RH1PO4imobilSoilHeter     => nmicf%RH1PO4imobilSoilHeter,      & !Net HPO4-P transfer from nonband soil to microbes; positive immobilization
    RH1PO4imobilBandHeter     => nmicf%RH1PO4imobilBandHeter,      & !Net HPO4-P transfer from fertilizer-band soil to microbes; positive immobilization
    RH1PO4imobilLitrHeter     => nmicf%RH1PO4imobilLitrHeter,      & !Net HPO4-P transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RCO2XumpAutor             => nmicf%RCO2XumpAutor,              & !CO2-C uptake for guild metabolism and biomass, including methanogenic CH4 production
    RN2FixHeter               => nmicf%RN2FixHeter,                & !Heterotrophic guild N2-N fixation into biomass
    RDecmpProdDOM             => ncplxf%RDecmpProdDOM,             & !Solid-substrate hydrolysis products entering DOM after humification
    RHydlysBioResduOM         => ncplxf%RHydlysBioResduOM,         & !Microbial residue C/N/P hydrolysis by compartment and complex K
    RHydlysSorptOM            => ncplxf%RHydlysSorptOM,            & !Sorbed-substrate hydrolysis by DOM species and complex K
    RDOMSorp                  => ncplxf%RDOMSorp,                  & !DOM transfer from solution to sorbed pools; positive adsorption, negative desorption
    TSensGrowth               => nmicdiag%TSensGrowth,             & !Layer temperature response for microbial growth [-]
    RH2UptkAutor              => nmicdiag%RH2UptkAutor,            & !Total H2 uptake by hydrogenotrophic methanogens, including biomass synthesis
    RNO2ReduxSoilChemo        => naqfdiag%RNO2ReduxSoilChemo,      & !Realized chemical nonband NO2-N consumption
    RNO2ReduxBandChemo        => naqfdiag%RNO2ReduxBandChemo,      & !Realized chemical fertilizer-band NO2-N consumption
    RN2OProdSoilChemo         => naqfdiag%RN2OProdSoilChemo,       & !N2O-N production by chemical nitrite decomposition in nonband soil
    RN2OProdBandChemo         => naqfdiag%RN2OProdBandChemo,       & !N2O-N production by chemical nitrite decomposition in fertilizer-band soil
    RNO3ProdSoilChemo         => naqfdiag%RNO3ProdSoilChemo,       & !NO3-N production by chemical nitrite decomposition in nonband soil
    RNO3ProdBandChemo         => naqfdiag%RNO3ProdBandChemo,       & !NO3-N production by chemical nitrite decomposition in fertilizer-band soil
    VOLWZ                     => nmicdiag%VOLWZ,                   & !Effective water volume supporting microbial activity and rate constraints
    DOMuptk4GrothAutor        => nmicf%DOMuptk4GrothAutor,         & !Guild elemental uptake; C source is CO2 or CH4 according to functional group
    RNOxReduxAutorBand        => nmicf%RNOxReduxAutorBand,         & !Autotrophic fertilizer-band oxidized-N reduction; reactant/product depend on functional group
    RNO3UptkAutor             => nmicf%RNO3UptkAutor,              & !Autotrophic nitrate-N reduction bookkeeping; zero for nitrifiers in current code
    RCH4ProdAutor             => nmicf%RCH4ProdAutor,              & !CH4-C production by autotrophic guild
    RGrowthCAutor             => nmicf%RGrowthCAutor,              & !Net substrate-derived C credited to autotrophic guild reserves
    RNOxReduxAutorSoil        => nmicf%RNOxReduxAutorSoil,         & !Autotrophic nonband oxidized-N reduction; reactant/product depend on functional group
    RO2UptkAutor              => nmicf%RO2UptkAutor,               & !Realized total O2 uptake by autotrophic guild
    RNOxReduxRespAutorLim     => nmicf%RNOxReduxRespAutorLim,      & !C-equivalent respiration supported by autotrophic nitrite reduction
    RCO2ProdAutor             => nmicf%RCO2ProdAutor,              & !CO2-C production by autotrophic guild
    RH1PO4TransfLitrAutor     => nmicf%RH1PO4TransfLitrAutor,      & !Net HPO4-P transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RH2PO4TransfLitrAutor     => nmicf%RH2PO4TransfLitrAutor,      & !Net H2PO4-P transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RNO3TransfLitrAutor       => nmicf%RNO3TransfLitrAutor,        & !Net NO3-N transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RNH4TransfLitrAutor       => nmicf%RNH4TransfLitrAutor,        & !Net NH4-N transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
    RN2FixAutor               => nmicf%RN2FixAutor,                & !Autotrophic guild N2 fixation flux (set to zero in current respiration gathering)
    RH1PO4TransfBandAutor     => nmicf%RH1PO4TransfBandAutor,      & !Net HPO4-P transfer from fertilizer-band soil to microbes; positive immobilization
    RH2PO4TransfBandAutor     => nmicf%RH2PO4TransfBandAutor,      & !Net H2PO4-P transfer from fertilizer-band soil to microbes; positive immobilization
    RNO3TransfBandAutor       => nmicf%RNO3TransfBandAutor,        & !Net NO3-N transfer from fertilizer-band soil to microbes; positive immobilization
    RNH4TransfBandAutor       => nmicf%RNH4TransfBandAutor,        & !Net NH4-N transfer from fertilizer-band soil to microbes; positive immobilization
    RH2PO4TransfSoilAutor     => nmicf%RH2PO4TransfSoilAutor,      & !Net H2PO4-P transfer from nonband soil to microbes; positive immobilization
    RNH4TransfSoilAutor       => nmicf%RNH4TransfSoilAutor,        & !Net NH4-N transfer from nonband soil to microbes; positive immobilization
    RNO3TransfSoilAutor       => nmicf%RNO3TransfSoilAutor,        & !Net NO3-N transfer from nonband soil to microbes; positive immobilization
    RH1PO4TransfSoilAutor     => nmicf%RH1PO4TransfSoilAutor,      & !Net HPO4-P transfer from nonband soil to microbes; positive immobilization
    TSens4MicbGrwoth          => micstt%TSens4MicbGrwoth,          & !Layer microbial growth-temperature response exported for coupling
    VWatMicrobAct             => micstt%VWatMicrobAct,             & !Effective microbial-activity water volume exported for coupling
    litrm                     => micfor%litrm,                     & !True for the surface litter layer
    Lsurf                     => micfor%Lsurf,                     & !True for the surface soil layer beneath litter
    k_POM                     => micpar%k_POM,                     & !Particulate-organic-matter complex index
    k_humus                   => micpar%k_humus,                   & !Humus complex index
    mid_AutoAmmoniaOxidBacter => micpar%mid_AutoAmmoniaOxidBacter, & !Functional-group identifier for ammonia oxidizers
    mid_AutoAeroCH4OxiBacter  => micpar%mid_AutoAeroCH4OxiBacter,  & !Functional-group identifier for aerobic methane oxidizers
    mid_AutoNitriteOxidBacter => micpar%mid_AutoNitriteOxidBacter, & !Functional-group identifier for nitrite oxidizers
    is_activeMicrbFungrpAutor => micpar%is_activeMicrbFungrpAutor, & !Activation flags for autotrophic functional groups
    mid_AutoAMONC10           => micpar%mid_AutoAMONC10,           & !Functional-group identifier for nitrite-dependent NC10 methanotrophs
    mid_AutoAMOANME2D         => micpar%mid_AutoAMOANME2D    ,     & !Functional-group identifier for nitrate-dependent ANME-2d methanotrophs
    RCH4MetaDmndAutor         => micflx%RCH4MetaDmndAutor,         & !Guild CH4 demand used in competition accounting
    RCH4UptkAutor             => micflx%RCH4UptkAutor,             & !Complete net microbial CH4-C uptake; negative values represent production
    RCO2NetUptkMicb           => micflx%RCO2NetUptkMicb,           & !Net microbial CO2-C uptake; negative values represent release
    RH2NetUptkMicb            => micflx%RH2NetUptkMicb,            & !Net microbial H2 uptake; negative values represent production
    RN2NetUptkMicb            => micflx%RN2NetUptkMicb,            & !Net N2 exchange from reduction pathways; negative values represent production
    RN2ONetUptkMicb           => micflx%RN2ONetUptkMicb,           & !Net microbial/chemical N2O-N uptake; negative values represent production
    RO2UptkMicb               => micflx%RO2UptkMicb,               & !Total microbial O2 uptake in the layer
    RH1PO4MicbReliz2Band      => micflx%RH1PO4MicbReliz2Band,      & !Net microbial HPO4-P release to fertilizer-band soil; positive supply, negative uptake
    RH1PO4MicbReliz2Soil      => micflx%RH1PO4MicbReliz2Soil,      & !Net microbial HPO4-P release to nonband soil; positive supply, negative uptake
    RH2PO4MicbReliz2Band      => micflx%RH2PO4MicbReliz2Band,      & !Net microbial H2PO4-P release to fertilizer-band soil; positive supply, negative uptake
    RH2PO4MicbReliz2Soil      => micflx%RH2PO4MicbReliz2Soil,      & !Net microbial H2PO4-P release to nonband soil; positive supply, negative uptake
    MicrbN2Fix                => micflx%MicrbN2Fix,                & !Total microbial N2-N fixation; positive means uptake
    RNH4MicbReliz2Band        => micflx%RNH4MicbReliz2Band,        & !Net microbial NH4-N release to fertilizer-band soil; positive supply, negative uptake
    RNH4MicbReliz2Soil        => micflx%RNH4MicbReliz2Soil,        & !Net microbial NH4-N release to nonband soil; positive supply, negative uptake
    RNO2MicbReliz2Band        => micflx%RNO2MicbReliz2Band,        & !Net microbial NO2-N release to fertilizer-band soil; positive supply, negative uptake
    RNO2MicbReliz2Soil        => micflx%RNO2MicbReliz2Soil,        & !Net microbial NO2-N release to nonband soil; positive supply, negative uptake
    RNO3MicbReliz2Band        => micflx%RNO3MicbReliz2Band,        & !Net microbial NO3-N release to fertilizer-band soil; positive supply, negative uptake
    RNO3MicbReliz2Soil        => micflx%RNO3MicbReliz2Soil,        & !Net microbial NO3-N release to nonband soil; positive supply, negative uptake
    REcoDOMProd               => micflx%REcoDOMProd                & !Net DOM production by species and complex after microbial uptake
  )
  
  D650: DO K=1,KL
    IF(.not.litrm .OR. (K.NE.k_POM .AND. K.NE.k_humus))THEN
      DO N=1,NumMicbHFunGrupsPerCmplx
        DO NGL=JGniH(N),JGnfH(N)
          naqfdiag%tRNH4MicrbImobilSoil   = naqfdiag%tRNH4MicrbImobilSoil+RNH4imobilSoilHeter(NGL,K)
          naqfdiag%tRNO3MicrbImobilSoil   = naqfdiag%tRNO3MicrbImobilSoil+RNO3imobilSoilHeter(NGL,K)
          naqfdiag%tRH2PO4MicrbImobilSoil = naqfdiag%tRH2PO4MicrbImobilSoil+RH2PO4imobilSoilHeter(NGL,K)
          naqfdiag%tRH1PO4MicrbImobilSoil = naqfdiag%tRH1PO4MicrbImobilSoil+RH1PO4imobilSoilHeter(NGL,K)  !> 0 uptake 
          naqfdiag%tRNH4MicrbImobilBand   = naqfdiag%tRNH4MicrbImobilBand+RNH4imobilBandHeter(NGL,K)
          naqfdiag%tRNO3MicrbImobilBand   = naqfdiag%tRNO3MicrbImobilBand+RNO3imobilBandHeter(NGL,K)
          naqfdiag%tRH2PO4MicrbImobilBand = naqfdiag%tRH2PO4MicrbImobilBand+RH2PO4imobilBandHeter(NGL,K)
          naqfdiag%tRH1PO4MicrbImobilBand = naqfdiag%tRH1PO4MicrbImobilBand+RH1PO4imobilBandHeter(NGL,K)
          naqfdiag%TFixN2                 = naqfdiag%TFixN2+RN2FixHeter(NGL,K)
      !
          naqfdiag%tRCO2MicrbProd     = naqfdiag%tRCO2MicrbProd+RCO2ProdHeter(NGL,K)
          naqfdiag%tRCH4MicrbProd     = naqfdiag%tRCH4MicrbProd+RCH4ProdHeter(NGL,K)
          naqfdiag%tRNOxMicrbRedux    = naqfdiag%tRNOxMicrbRedux+RNOxDOCReduxRespDenitLim(NGL,K)+RNOxAcetReduxRespDenitLim(NGL,K)
          naqfdiag%tRO2MicrbUptk      = naqfdiag%tRO2MicrbUptk+RO2UptkHeter(NGL,K)
          naqfdiag%TReduxNO3Soil      = naqfdiag%TReduxNO3Soil+RNO3ReduxHeterSoil(NGL,K)       !NO3->NO2
          naqfdiag%TReduxNO3Band      = naqfdiag%TReduxNO3Band+RNO3ReduxHeterBand(NGL,K)       !NO3->NO2
          naqfdiag%TDeniReduxNO2Soil  = naqfdiag%TDeniReduxNO2Soil+RNO2ReduxHeterSoil(NGL,K)   !NO2->N2O
          naqfdiag%TDeniReduxNO2Band  = naqfdiag%TDeniReduxNO2Band+RNO2ReduxHeterBand(NGL,K)   !NO2->N2O
          naqfdiag%TReduxNO2toN2OSoil = naqfdiag%TReduxNO2toN2OSoil+RNO2ReduxHeterSoil(NGL,K)
          naqfdiag%TReduxNO2toN2OBand = naqfdiag%TReduxNO2toN2OBand+RNO2ReduxHeterBand(NGL,K)
          naqfdiag%TReduxN2OtoN2      = naqfdiag%TReduxN2OtoN2+RN2OReduxHeter(NGL,K)               !N2O -> N2
          naqfdiag%TProdH2            = naqfdiag%TProdH2+RH2ProdHeter(NGL,K)
          O2RespUptake=RO2UptkHeter(NGL,K)
          O2RespDemand=RO2DmndHeter(NGL,K)
          IF(N.EQ.micpar%mid_HeterMixtCynoBacter)THEN
            !Net photosynthetic O2 release belongs in gas exchange, not in
            !the actual/potential respiratory O2 ratio exported to the API.
            O2RespUptake=RO2Uptk4RespHeter(NGL,K)
            O2RespDemand=RO2Dmnd4RespHeter(NGL,K)
          ENDIF
          naqfdiag%tRO2UptkHeterG = naqfdiag%tRO2UptkHeterG+O2RespUptake
          naqfdiag%tRO2DmndHeterG = naqfdiag%tRO2DmndHeterG+O2RespDemand
          nmicf%RO2UptkHeterG(NGL)=nmicf%RO2UptkHeterG(NGL)+O2RespUptake
          nmicf%RO2DmndHeterG(NGL)=nmicf%RO2DmndHeterG(NGL)+O2RespDemand
          TRDOM2DIE(ielmc)            = TRDOM2DIE(ielmc)+RCO2ProdHeter(NGL,K)+RCH4ProdHeter(NGL,K)
        ENDDO
      ENDDO
    ENDIF
  ENDDO D650

  DO  N=1,NumMicbAFunGrupsPerCmplx
    IF(is_activeMicrbFungrpAutor(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        naqfdiag%tRNH4MicrbImobilSoil   = naqfdiag%tRNH4MicrbImobilSoil+RNH4TransfSoilAutor(NGL)
        naqfdiag%tRNO3MicrbImobilSoil   = naqfdiag%tRNO3MicrbImobilSoil+RNO3TransfSoilAutor(NGL)
        naqfdiag%tRH2PO4MicrbImobilSoil = naqfdiag%tRH2PO4MicrbImobilSoil+RH2PO4TransfSoilAutor(NGL)
        naqfdiag%tRH1PO4MicrbImobilSoil = naqfdiag%tRH1PO4MicrbImobilSoil+RH1PO4TransfSoilAutor(NGL)
        naqfdiag%tRNH4MicrbImobilBand   = naqfdiag%tRNH4MicrbImobilBand+RNH4TransfBandAutor(NGL)
        naqfdiag%tRNO3MicrbImobilBand   = naqfdiag%tRNO3MicrbImobilBand+RNO3TransfBandAutor(NGL)
        naqfdiag%tRH2PO4MicrbImobilBand = naqfdiag%tRH2PO4MicrbImobilBand+RH2PO4TransfBandAutor(NGL)
        naqfdiag%tRH1PO4MicrbImobilBand = naqfdiag%tRH1PO4MicrbImobilBand+RH1PO4TransfBandAutor(NGL)
        naqfdiag%TFixN2                 = naqfdiag%TFixN2+RN2FixAutor(NGL)
        naqfdiag%tRCO2MicrbProd       = naqfdiag%tRCO2MicrbProd+RCO2ProdAutor(NGL)
        naqfdiag%tRCH4MicrbProd       = naqfdiag%tRCH4MicrbProd+RCH4ProdAutor(NGL)
        naqfdiag%tRNOxMicrbRedux      = naqfdiag%tRNOxMicrbRedux+RNOxReduxRespAutorLim(NGL)
        naqfdiag%tRO2MicrbUptk        = naqfdiag%tRO2MicrbUptk+RO2UptkAutor(NGL)
        naqfdiag%TReduxNO3Soil        = naqfdiag%TReduxNO3Soil+RNO3UptkAutor(NGL)
        
        IF(N.EQ.mid_AutoAmmoniaOxidBacter)THEN
          !Nitrite-N consumed in NO2 -> N2O.
          naqfdiag%TReduxNO2toN2OSoil = naqfdiag%TReduxNO2toN2OSoil &
            +RNOxReduxAutorSoil(NGL)
          naqfdiag%TReduxNO2toN2OBand = naqfdiag%TReduxNO2toN2OBand &
            +RNOxReduxAutorBand(NGL)

          !Additional ammonia-N incorporated into N2O.
          NH4toN2OSoil = NH4toN2OSoil+0.5_r8*RNOxReduxAutorSoil(NGL)
          NH4toN2OBand = NH4toN2OBand+0.5_r8*RNOxReduxAutorBand(NGL)

          !Total nitrifier N2O-N production, including ammonia-N.
          naqfdiag%TNitNO2Redux2N2OSoil = naqfdiag%TNitNO2Redux2N2OSoil &
            +1.5_r8*RNOxReduxAutorSoil(NGL)
          naqfdiag%TNitNO2Redux2N2OBand = naqfdiag%TNitNO2Redux2N2OBand &
            +1.5_r8*RNOxReduxAutorBand(NGL)

        ELSEIF(N.EQ.mid_AutoAMONC10)THEN
          !NC10: NO2-N consumed equals N2-N produced.
          NC10NO2ReduxSoil = NC10NO2ReduxSoil+RNOxReduxAutorSoil(NGL)
          NC10NO2ReduxBand = NC10NO2ReduxBand+RNOxReduxAutorBand(NGL)

        ELSEIF(N.EQ.mid_AutoAMOANME2D)THEN
          !ANME-2d: NO3-N consumed equals NO2-N produced.
          naqfdiag%TReduxNO3Soil = naqfdiag%TReduxNO3Soil+RNOxReduxAutorSoil(NGL)
          naqfdiag%TReduxNO3Band = naqfdiag%TReduxNO3Band+RNOxReduxAutorBand(NGL)
        ENDIF

      ENDDO
    ENDIF
  ENDDO

  IF(Lsurf)THEN
    naqfdiag%tRNH4MicrbImobilSoil   = naqfdiag%tRNH4MicrbImobilSoil+micfor%tRNH4MicrbImobilSoil
    naqfdiag%tRNH4MicrbImobilBand   = naqfdiag%tRNH4MicrbImobilBand+micfor%tRNH4MicrbImobilBand
    naqfdiag%tRNO3MicrbImobilSoil   = naqfdiag%tRNO3MicrbImobilSoil+micfor%tRNO3MicrbImobilSoil
    naqfdiag%tRNO3MicrbImobilBand   = naqfdiag%tRNO3MicrbImobilBand+micfor%tRNO3MicrbImobilBand
    naqfdiag%tRH2PO4MicrbImobilSoil = naqfdiag%tRH2PO4MicrbImobilSoil+micfor%tRH2PO4MicrbImobilSoil
    naqfdiag%tRH2PO4MicrbImobilBand = naqfdiag%tRH2PO4MicrbImobilBand+micfor%tRH2PO4MicrbImobilBand
    naqfdiag%tRH1PO4MicrbImobilSoil = naqfdiag%tRH1PO4MicrbImobilSoil+micfor%tRH1PO4MicrbImobilSoil
    naqfdiag%tRH1PO4MicrbImobilBand = naqfdiag%tRH1PO4MicrbImobilBand+micfor%tRH1PO4MicrbImobilBand
  ENDIF

  ! tRCO2GrothAutor=total CO2 uptake by autotrophs, ammonia oxidizer
  ! nitrite oxidizer, and hydrogenotrophic methanogens,
  ! all of which involves CO2 for both energy and C biomass.
  D645: DO N=1,NumMicbAFunGrupsPerCmplx
    IF(micpar%is_CO2_autotroph(N))THEN
      DO NGL=JGniA(N),JGnfA(N)
        naqfdiag%tRCO2GrothAutor=naqfdiag%tRCO2GrothAutor+RCO2XumpAutor(NGL)
      ENDDO
    ENDIF
  ENDDO D645

!
!     ALLOCATE AGGREGATED TRANSFOMBioResduATIONS INTO ARRAYS TO UPDATE
!     STATE VARIABLES IN 'REDIST'
!
!     RCO2NetUptkMicb=net CO2 uptake, < 0 producig CO2
!     tRCO2MicrbProd total CO2 emission by heterotrophs reducing O2
!     tRNOxMicrbRedux=total CO2 emission by denitrifiers reducing NOx
!     RSMetaOxidSoilAutor(3)=CH4 oxidation
!     RCH4UptkAutor=net CH4 uptake, >0, means uptake
!     DOMUptk4GrothHeter=total CH4 uptake by autotrophs
!     tRCH4MicrbProd=total CH4 emission
!     RH2NetUptkMicb=net H2 uptake
!     RH2UptkAutor,TProdH2=total H2 uptake, emission
!     RO2UptkMicb,tRO2MicrbUptk=total O2 uptake
!     RN2NetUptkMicb=total N2 production
!     TReduxN2OtoN2=total N2O reduction
!     RN2ONetUptkMicb=total N2O uptake
!     NO2(-) -> N2O
!     TReduxNO2toN2OSoil,TReduxNO2toN2OBand=total NO2 reduction in non-band,band
!     RN2OProdSoilChemo,RN2OProdBandChemo=nitrous acid reduction in non-band,band
!
  RCO2NetUptkMicb = naqfdiag%tRCO2GrothAutor-naqfdiag%tRCO2MicrbProd-naqfdiag%tRNOxMicrbRedux
  RCH4UptkAutor   = -naqfdiag%tRCH4MicrbProd

  DO  N=1,NumMicbAFunGrupsPerCmplx
    if(N.eq.mid_AutoAeroCH4OxiBacter)then
      DO NGL=JGniA(N),JGnfA(N)
        !use CH4 for both energy and biomass
        !CO2 is released from CH4 oxidation
        RCO2NetUptkMicb = RCO2NetUptkMicb-RSMetaOxidSoilAutor(NGL)
        !total CH4 uptake includes catabolic oxidation and biomass synthesis
        RCH4UptkAutor          = RCH4UptkAutor+RSMetaOxidSoilAutor(NGL)+DOMuptk4GrothAutor(ielmc,NGL)
        RCH4MetaDmndAutor(NGL) = RCH4MetaDmndAutor(NGL)+DOMuptk4GrothAutor(ielmc,NGL)
        TRDOM2DIE(ielmc)       = TRDOM2DIE(ielmc)-RGrowthCAutor(NGL)
      ENDDO
    elseif(N.eq.mid_AutoAMONC10 .or. N.eq.mid_AutoAMOANME2D)then
      DO NGL=JGniA(N),JGnfA(N)
        RCH4UptkAutor    = RCH4UptkAutor+RSMetaOxidSoilAutor(NGL)
        TRDOM2DIE(ielmc) = TRDOM2DIE(ielmc)-RGrowthCAutor(NGL)
      ENDDO
    ENDIF
  ENDDO

  !>0. microbial uptake
  RH2NetUptkMicb  = RH2UptkAutor-naqfdiag%TProdH2
  RO2UptkMicb     = naqfdiag%tRO2MicrbUptk
  RN2NetUptkMicb  = -naqfdiag%TReduxN2OtoN2-NC10NO2ReduxSoil-NC10NO2ReduxBand
  RN2ONetUptkMicb = -naqfdiag%TReduxNO2toN2OSoil-naqfdiag%TReduxNO2toN2OBand-RN2OProdSoilChemo &
    -RN2OProdBandChemo-NH4toN2OSoil-NH4toN2OBand+naqfdiag%TReduxN2OtoN2
!
  D655: DO K=1,jcplx
    D660: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        REcoDOMProd(NE,K)=REcoDOMProd(NE,K)+RDecmpProdDOM(NE,M,K)
      ENDDO
    ENDDO D660

    D665: DO M=1,ndbiomcp
      DO NE=1,NumPlantChemElms
        REcoDOMProd(NE,K)=REcoDOMProd(NE,K)+RHydlysBioResduOM(NE,M,K)
      ENDDO
    ENDDO D665
    DO NE=1,NumPlantChemElms
      REcoDOMProd(NE,K)=REcoDOMProd(NE,K)+RHydlysSorptOM(NE,K)
    ENDDO
    REcoDOMProd(idom_acetate,K)=REcoDOMProd(idom_acetate,K)+RHydlysSorptOM(idom_acetate,K)
    D670: DO N=1,NumMicbHFunGrupsPerCmplx
      DO NGL=JGniH(N),JGnfH(N)
        REcoDOMProd(idom_doc,K)     = REcoDOMProd(idom_doc,K)-RMetabDOCUptkHeter(NGL,K)
        REcoDOMProd(idom_don,K)     = REcoDOMProd(idom_don,K)-DOMUptk4GrothHeter(ielmn,NGL,K)
        REcoDOMProd(idom_dop,K)     = REcoDOMProd(idom_dop,K)-DOMUptk4GrothHeter(ielmp,NGL,K)
        REcoDOMProd(idom_acetate,K) = REcoDOMProd(idom_acetate,K)-RMetabAcetUptkHeter(NGL,K)+RAcetateProdHeter(NGL,K)
      ENDDO
    ENDDO D670
    DO idom=idom_beg,idom_end
      REcoDOMProd(idom,K)=REcoDOMProd(idom,K)-RDOMSorp(idom,K)
    ENDDO
  ENDDO D655

  !Chemical nitrite reduction yields 80% NO3-N, 10% N2O-N and 10% DON-N.
  !The layer-wide reaction has no complex-specific donor: distribute its
  !DON equally among active complexes (litter excludes POM and humus).
  !Credit both the updated pool exported by MicAPIRecv and its DOM flux.
  DONChemoPerComplex = naqfdiag%RNO2ReduxChemo/REAL(KL,r8)
  micstt%DOM(idom_don,1:KL) = micstt%DOM(idom_don,1:KL)+DONChemoPerComplex
  REcoDOMProd(idom_don,1:KL) = REcoDOMProd(idom_don,1:KL)+DONChemoPerComplex
!
!     RNH4MicbReliz2Soil,RNH4MicbReliz2Band=net change in NH4 in band,non-band
!     tRNH4MicrbImobilSoil,tRNH4MicrbImobilBand=total NH4 mineraln-immobn in non-band,band
!     RSMetaOxidSoilAutor(1),RSMetaOxidBandAutor(1)=total NH4 oxidation in non-band,band
!     RNO3MicbReliz2Soil,RNO3MicbReliz2Band=net change in NO3 in band,non-band
!     tRNO3MicrbImobilSoil,tRNO3MicrbImobilBand=total NO3 immobn in non-band,band
!     RSMetaOxidSoilAutor(2),RSMetaOxidBandAutor(2)=total NO2 oxidation in non-band,band
!     TReduxNO3Soil,TReduxNO3Band=total NO3 reduction in non-band,band
!     RNO3ProdSoilChemo,RNO3ProdBandChemo=NO3 production from nitrous acid reduction in non-band,band
!     RNO2MicbReliz2Soil,RNO2MicbReliz2Band=net change in NO3 in band,non-band
!     TReduxNO2toN2OSoil,TReduxNO2toN2OBand=total NO2 reduction in non-band,band
!     RNO2ReduxSoilChemo,RNO2ReduxBandChemo=substrate-limited nitrous acid reduction in non-band,band
!     RH2PO4MicbReliz2Soil,RH2PO4MicbReliz2Band=net change in H2PO4 in band,non-band
!     tRH2PO4MicrbImobilSoil,tRH2PO4MicrbImobilBand=total H2PO4 mineraln-immobn in non-band,band
!     RH1PO4MicbReliz2Soil,RH1PO4MicbReliz2Band=net change in HPO4 in band,non-band
!     tRH1PO4MicrbImobilSoil,tRH1PO4MicrbImobilBand=total HPO4 mineraln-immobn in non-band,band
!     MicrbN2Fix=total N2 fixation
!     XZHYS=total H+ production
!     TFixN2=total N2 fixation
!

  micflx%TRDOM2DIE(ielmn)=micflx%TRDOM2DIE(ielmn)+naqfdiag%tRNH4MicrbImobilSoil &
    +naqfdiag%tRNH4MicrbImobilBand+naqfdiag%tRNO3MicrbImobilSoil &
    +naqfdiag%tRNO3MicrbImobilBand+naqfdiag%TFixN2+naqfdiag%RNO2ReduxChemo
  micflx%TRDOM2DIE(ielmp)=micflx%TRDOM2DIE(ielmp)+naqfdiag%tRH1PO4MicrbImobilSoil &
    +naqfdiag%tRH2PO4MicrbImobilSoil+naqfdiag%tRH1PO4MicrbImobilBand &
    +naqfdiag%tRH2PO4MicrbImobilBand


  RNH4MicbReliz2Soil=-naqfdiag%tRNH4MicrbImobilSoil
  RNO3MicbReliz2Soil=-naqfdiag%tRNO3MicrbImobilSoil-naqfdiag%TReduxNO3Soil+RNO3ProdSoilChemo
  RNO2MicbReliz2Soil=+naqfdiag%TReduxNO3Soil-naqfdiag%TReduxNO2toN2OSoil &
    -NC10NO2ReduxSoil-RNO2ReduxSoilChemo
  RH2PO4MicbReliz2Soil=-naqfdiag%tRH2PO4MicrbImobilSoil
  RH1PO4MicbReliz2Soil=-naqfdiag%tRH1PO4MicrbImobilSoil     !< 0 uptake
  RNH4MicbReliz2Band=-naqfdiag%tRNH4MicrbImobilBand
  RNO3MicbReliz2Band=-naqfdiag%tRNO3MicrbImobilBand-naqfdiag%TReduxNO3Band+RNO3ProdBandChemo
  RNO2MicbReliz2Band=naqfdiag%TReduxNO3Band-naqfdiag%TReduxNO2toN2OBand &
    -NC10NO2ReduxBand-RNO2ReduxBandChemo

  !mid_AutoAmmoniaOxidBacter=1, mid_AutoNitriteOxidBacter=2, mid_AutoAeroCH4OxiBacter=3
  DO NGL=JGniA(mid_AutoAmmoniaOxidBacter),JGnfA(mid_AutoAmmoniaOxidBacter)
    RNH4MicbReliz2Soil=RNH4MicbReliz2Soil-RSMetaOxidSoilAutor(NGL)   !some NH3 -> NO2, some NH3-> N2O
    RNH4MicbReliz2Band=RNH4MicbReliz2Band-RSMetaOxidBandAutor(NGL)

    !on mole-basis 2NO2(-) + NH3 -> 1.5N2O + 2OH(-) + 0.5H2O, RNOxReduxAutorSoil(NGL)/2._r8 goes to N2O
    RNO2MicbReliz2Soil=RNO2MicbReliz2Soil+RSMetaOxidSoilAutor(NGL)-RNOxReduxAutorSoil(NGL)/2._r8
    RNO2MicbReliz2Band=RNO2MicbReliz2Band+RSMetaOxidBandAutor(NGL)-RNOxReduxAutorBand(NGL)/2._r8
  ENDDO

  DO NGL=JGniA(mid_AutoNitriteOxidBacter),JGnfA(mid_AutoNitriteOxidBacter)
    naqfdiag%tNO2OxiAuto=naqfdiag%tNO2OxiAuto+RSMetaOxidSoilAutor(NGL)+RSMetaOxidBandAutor(NGL)
    RNO3MicbReliz2Soil=RNO3MicbReliz2Soil+RSMetaOxidSoilAutor(NGL)
    RNO2MicbReliz2Soil=RNO2MicbReliz2Soil-RSMetaOxidSoilAutor(NGL)
    RNO3MicbReliz2Band=RNO3MicbReliz2Band+RSMetaOxidBandAutor(NGL)
    RNO2MicbReliz2Band=RNO2MicbReliz2Band-RSMetaOxidBandAutor(NGL)
  ENDDO

  RH2PO4MicbReliz2Band = -naqfdiag%tRH2PO4MicrbImobilBand
  RH1PO4MicbReliz2Band = -naqfdiag%tRH1PO4MicrbImobilBand
  MicrbN2Fix           = naqfdiag%TFixN2
  !how to aggregate the mean T sensitivity?
  TSens4MicbGrwoth     = TSensGrowth
  VWatMicrobAct        = VOLWZ

  DO NGL=JGNiA(mid_AutoAeroCH4OxiBacter),jGnfA(mid_AutoAeroCH4OxiBacter)
    naqfdiag%tCH4OxiAero=naqfdiag%tCH4OxiAero+RSMetaOxidSoilAutor(NGL)
  ENDDO

  DO NGL=JGNiA(mid_AutoAMONC10),jGnfA((mid_AutoAMONC10))
    naqfdiag%tCH4OxiANMO=naqfdiag%tCH4OxiANMO+RSMetaOxidSoilAutor(NGL)
  ENDDO

  DO NGL=JGNiA(mid_AutoAMOANME2D),jGnfA(mid_AutoAMOANME2D)
    naqfdiag%tCH4OxiANMO=naqfdiag%tCH4OxiANMO+RSMetaOxidSoilAutor(NGL)    
  ENDDO  
  naqfdiag%tRNH3Oxi=nmicf%RTotNH3OxidSoilAutor+nmicf%RTotNH3OxidBandAutor

  end associate
  end subroutine AggregateTransfOMBioResdue
!------------------------------------------------------------------------------------------

  subroutine SubstrateAttenf4Compet(NGL,N,K,FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,&
    micfor,naqfdiag,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(out):: FNH4X,FNB4X            !fraction of NH4 demand over all plant+microbial demand, soil/band 
  real(r8),intent(out) :: FNO3X,FNB3X            !fraction of NO3 demand over all plant+microbial demand, soil/band
  real(r8),intent(out) :: FPO4X,FPOBX            !fraction of H2PO4 demand over all plant+microbial demand, soil/band
  real(r8), intent(out):: FP14X,FP1BX            !fraction of H1PO4 demand over all plant+microbial demand, soil/band
 
  type(micforctype), intent(in) :: micfor
  type(Cumlate_Flux_Diag_type),INTENT(INOUT)::  naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx
! begin_execution
  associate(                                                   &
    FracOMActHeter          => nmics%FracOMActHeter,           & !Guild/complex fraction of total active microbial C in the layer [-]
    AttenfNH4Heter          => micflx%AttenfNH4Heter,          & !Litter-microbial share of NH4-N uptake from underlying soil [-]
    AttenfNO3Heter          => micflx%AttenfNO3Heter,          & !Litter-microbial share of NO3-N uptake from underlying soil [-]
    AttenfH1PO4Heter        => micflx%AttenfH1PO4Heter,        & !Litter-microbial share of HPO4-P uptake from underlying soil [-]
    AttenfH2PO4Heter        => micflx%AttenfH2PO4Heter,        & !Litter-microbial share of H2PO4-P uptake from underlying soil [-]
    RO2DmndHetertPrev       => micflx%RO2DmndHetertPrev,       & !Previous-hour guild O2 demand used for competition; not referenced here
    RNH4DmndSoilHeterPrev   => micflx%RNH4DmndSoilHeterPrev,   & !Previous-hour potential NH4-N uptake by heterotrophic guilds from nonband soil; used for competition
    RNH4DmndBandHeterPrev   => micflx%RNH4DmndBandHeterPrev,   & !Previous-hour potential NH4-N uptake by heterotrophic guilds from fertilizer-band soil; used for competition
    RNO3DmndSoilHeterPrev   => micflx%RNO3DmndSoilHeterPrev,   & !Previous-hour potential NO3-N uptake by heterotrophic guilds from nonband soil; used for competition
    RNO3DmndBandHeterPrev   => micflx%RNO3DmndBandHeterPrev,   & !Previous-hour potential NO3-N uptake by heterotrophic guilds from fertilizer-band soil; used for competition
    RH2PO4DmndSoilHeterPrev => micflx%RH2PO4DmndSoilHeterPrev, & !Previous-hour potential H2PO4-P uptake by heterotrophic guilds from nonband soil; used for competition
    RH2PO4DmndBandHeterPrev => micflx%RH2PO4DmndBandHeterPrev, & !Previous-hour potential H2PO4-P uptake by heterotrophic guilds from fertilizer-band soil; used for competition
    RH1PO4DmndSoilHeterPrev => micflx%RH1PO4DmndSoilHeterPrev, & !Previous-hour potential HPO4-P uptake by heterotrophic guilds from nonband soil; used for competition
    RH1PO4DmndBandHeterPrev => micflx%RH1PO4DmndBandHeterPrev, & !Previous-hour potential HPO4-P uptake by heterotrophic guilds from fertilizer-band soil; used for competition
    RNH4DmndLitrHeterPrev   => micflx%RNH4DmndLitrHeterPrev,   & !Previous-hour potential NH4-N uptake by heterotrophic guilds from underlying soil accessed by litter microbes; used for competition
    RNO3DmndLitrHeterPrev   => micflx%RNO3DmndLitrHeterPrev,   & !Previous-hour potential NO3-N uptake by heterotrophic guilds from underlying soil accessed by litter microbes; used for competition
    RH2PO4DmndLitrHeterPrev => micflx%RH2PO4DmndLitrHeterPrev, & !Previous-hour potential H2PO4-P uptake by heterotrophic guilds from underlying soil accessed by litter microbes; used for competition
    RH1PO4DmndLitrHeterPrev => micflx%RH1PO4DmndLitrHeterPrev, & !Previous-hour potential HPO4-P uptake by heterotrophic guilds from underlying soil accessed by litter microbes; used for competition
    litrm                   => micfor%litrm,                   & !True for the surface litter layer
    VLNH4                   => micfor%VLNH4,                   & !Nonband fraction for ammonium/ammonia pools and uptake capacity [-]
    VLNHB                   => micfor%VLNHB,                   & !Fertilizer-band fraction for ammonium/ammonia pools and uptake capacity [-]
    VLNOB                   => micfor%VLNOB,                   & !Fertilizer-band fraction for nitrate/nitrite pools and uptake capacity [-]
    VLNO3                   => micfor%VLNO3,                   & !Nonband fraction for nitrate/nitrite pools and uptake capacity [-]
    VLPOB                   => micfor%VLPOB,                   & !Fertilizer-band fraction for phosphate pools and uptake capacity [-]
    VLPO4                   => micfor%VLPO4,                   & !Nonband fraction for phosphate pools and uptake capacity [-]
    ZEROS                   => micfor%ZEROS,                   & !Small mass or flux threshold used by the routine
    RNH4EcoDmndSoilPrev     => micfor%RNH4EcoDmndSoilPrev,     & !Previous-hour ecosystem NH4-N demand in nonband soil; competition denominator
    RNH4EcoDmndBandPrev     => micfor%RNH4EcoDmndBandPrev,     & !Previous-hour ecosystem NH4-N demand in fertilizer-band soil; competition denominator
    RNO3EcoDmndSoilPrev     => micfor%RNO3EcoDmndSoilPrev,     & !Previous-hour ecosystem NO3-N demand in nonband soil; competition denominator
    RNH4EcoDmndLitrPrev     => micfor%RNH4EcoDmndLitrPrev,     & !Previous-hour ecosystem NH4-N demand in underlying soil accessed by litter microbes; competition denominator
    RNO3EcoDmndLitrPrev     => micfor%RNO3EcoDmndLitrPrev,     & !Previous-hour ecosystem NO3-N demand in underlying soil accessed by litter microbes; competition denominator
    RH1PO4EcoDmndLitrPrev   => micfor%RH1PO4EcoDmndLitrPrev,   & !Previous-hour ecosystem HPO4-P demand in underlying soil accessed by litter microbes; competition denominator
    RH2PO4EcoDmndLitrPrev   => micfor%RH2PO4EcoDmndLitrPrev,   & !Previous-hour ecosystem H2PO4-P demand in underlying soil accessed by litter microbes; competition denominator
    RNO3EcoDmndBandPrev     => micfor%RNO3EcoDmndBandPrev,     & !Previous-hour ecosystem NO3-N demand in fertilizer-band soil; competition denominator
    RH2PO4EcoDmndSoilPrev   => micfor%RH2PO4EcoDmndSoilPrev,   & !Previous-hour ecosystem H2PO4-P demand in nonband soil; competition denominator
    RH2PO4EcoDmndBandPrev   => micfor%RH2PO4EcoDmndBandPrev,   & !Previous-hour ecosystem H2PO4-P demand in fertilizer-band soil; competition denominator
    RH1PO4EcoDmndSoilPrev   => micfor%RH1PO4EcoDmndSoilPrev,   & !Previous-hour ecosystem HPO4-P demand in nonband soil; competition denominator
    RH1PO4EcoDmndBandPrev   => micfor%RH1PO4EcoDmndBandPrev,   & !Previous-hour ecosystem HPO4-P demand in fertilizer-band soil; competition denominator
    RDOMEcoDmndPrev         => micfor%RDOMEcoDmndPrev,         & !Previous-hour ecosystem DOC demand in each complex; competition denominator; not referenced here
    RAcetateEcoDmndPrev     => micfor%RAcetateEcoDmndPrev,     & !Previous-hour ecosystem acetate demand in each complex; competition denominator; not referenced here
    Lsurf                   => micfor%Lsurf,                   & !True for the surface soil layer beneath litter; not referenced here
    SoilMicPMassLayer0      => micfor%SoilMicPMassLayer0       & !Surface-litter soil-mass reference used in litter/soil exchange conditions; not referenced here
  )
! F*=fraction of substrate uptake relative to total uptake from
! previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
! NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
! POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
! oxidation, OQA=acetate oxidation
!
  
  IF(RNH4EcoDmndSoilPrev.GT.ZEROS)THEN
    FNH4X=AMAX1(FMN,RNH4DmndSoilHeterPrev(NGL,K)/RNH4EcoDmndSoilPrev)
  ELSE
    FNH4X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNH4)
  ENDIF
  IF(RNH4EcoDmndBandPrev.GT.ZEROS)THEN
    FNB4X=AMAX1(FMN,RNH4DmndBandHeterPrev(NGL,K)/RNH4EcoDmndBandPrev)
  ELSE
    FNB4X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNHB)
  ENDIF
  IF(RNO3EcoDmndSoilPrev.GT.ZEROS)THEN
    FNO3X=AMAX1(FMN,RNO3DmndSoilHeterPrev(NGL,K)/RNO3EcoDmndSoilPrev)
  ELSE
    FNO3X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNO3)
  ENDIF
  IF(RNO3EcoDmndBandPrev.GT.ZEROS)THEN
    FNB3X=AMAX1(FMN,RNO3DmndBandHeterPrev(NGL,K)/RNO3EcoDmndBandPrev)
  ELSE
    FNB3X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLNOB)
  ENDIF
  IF(RH2PO4EcoDmndSoilPrev.GT.ZEROS)THEN
    FPO4X=AMAX1(FMN,RH2PO4DmndSoilHeterPrev(NGL,K)/RH2PO4EcoDmndSoilPrev)
  ELSE
    FPO4X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPO4)
  ENDIF
  IF(RH2PO4EcoDmndBandPrev.GT.ZEROS)THEN
    FPOBX=AMAX1(FMN,RH2PO4DmndBandHeterPrev(NGL,K)/RH2PO4EcoDmndBandPrev)
  ELSE
    FPOBX=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPOB)
  ENDIF
  IF(RH1PO4EcoDmndSoilPrev.GT.ZEROS)THEN
    FP14X=AMAX1(FMN,RH1PO4DmndSoilHeterPrev(NGL,K)/RH1PO4EcoDmndSoilPrev)
  ELSE
    FP14X=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPO4)
  ENDIF
  IF(RH1PO4EcoDmndBandPrev.GT.ZEROS)THEN
    FP1BX=AMAX1(FMN,RH1PO4DmndBandHeterPrev(NGL,K)/RH1PO4EcoDmndBandPrev)
  ELSE
    FP1BX=AMAX1(FMN,FracOMActHeter(NGL,K)*VLPOB)
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
  !All litter complexes and autotrophs tap the same underlying-soil pool.
  !Use their shared layer biomass denominator when previous demand is absent.
  IF(litrm)THEN
    IF(RNH4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfNH4Heter(NGL,K)=AMAX1(FMN,RNH4DmndLitrHeterPrev(NGL,K)/RNH4EcoDmndLitrPrev)
    ELSE
      AttenfNH4Heter(NGL,K)=AMAX1(FMN,FracOMActHeter(NGL,K))
    ENDIF
    IF(RNO3EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfNO3Heter(NGL,K)=AMAX1(FMN,RNO3DmndLitrHeterPrev(NGL,K)/RNO3EcoDmndLitrPrev)
    ELSE
      AttenfNO3Heter(NGL,K)=AMAX1(FMN,FracOMActHeter(NGL,K))
    ENDIF
    IF(RH2PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfH2PO4Heter(NGL,K)=AMAX1(FMN,RH2PO4DmndLitrHeterPrev(NGL,K)/RH2PO4EcoDmndLitrPrev)
    ELSE
      AttenfH2PO4Heter(NGL,K)=AMAX1(FMN,FracOMActHeter(NGL,K))
    ENDIF
    IF(RH1PO4EcoDmndLitrPrev.GT.ZEROS)THEN
      AttenfH1PO4Heter(NGL,K)=AMAX1(FMN,RH1PO4DmndLitrHeterPrev(NGL,K)/RH1PO4EcoDmndLitrPrev)
    ELSE
      AttenfH1PO4Heter(NGL,K)=AMAX1(FMN,FracOMActHeter(NGL,K))
    ENDIF
  ENDIF

!diagnostics off
!  IF(Lsurf .AND. K.NE.micpar%k_POM .AND. K.NE.micpar%k_humus .AND. SoilMicPMassLayer0.GT.ZEROS)THEN
!    naqfdiag%TFNH4X=naqfdiag%TFNH4X+micfor%AttenfNH4HeterR(NGL,K)
!    naqfdiag%TFNO3X=naqfdiag%TFNO3X+micfor%AttenfNO3HeterR(NGL,K)
!    naqfdiag%TFPO4X=naqfdiag%TFPO4X+micfor%AttenfH2PO4HeterR(NGL,K)
!    naqfdiag%TFP14X=naqfdiag%TFP14X+micfor%AttenfH1PO4HeterR(NGL,K)
!  ENDIF
  end associate
  end subroutine SubstrateAttenf4Compet

!------------------------------------------------------------------------------------------

  subroutine BiomassMineralization(NGL,N,K,FNH4X,FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX, &
    ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T,micfor,micstt,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: FNH4X
  real(r8), intent(in) :: FNB3X,FNB4X,FNO3X
  real(r8), intent(in) :: FPO4X,FPOBX,FP14X,FP1BX
  real(r8), intent(in) :: ZNH4T,ZNO3T,ZNO2T,H2P4T,H1P4T
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type),intent(inout) :: nmics
  character(len=*), parameter :: subname='BiomassMineralization'

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
  real(r8) :: RINHP
  real(r8) :: RINHX,RINOP,RINOX,RIPOP,RIPOX,RIP1P
  real(r8) :: RIP1X,RINHPR,RINOPR,RIPOPR,RIP1PR
  real(r8) :: ZNH4M,ZNHBM
  real(r8) :: ZNO3M
  real(r8) :: ZNOBM
  integer :: MID3

  real(r8) :: LitrPotential(2),LitrUptake(2) !Nonband and band soil fluxes
!     begin_execution
  associate(                                                       &
   GrowthEnvScalHeter              => nmics%GrowthEnvScalHeter,    & !Temperature and water-potential multiplier on heterotrophic growth [-]
   OMActHeter                      => nmics%OMActHeter,            & !Active microbial C biomass by heterotrophic guild and complex K
   AttenfNH4Heter                  => micflx%AttenfNH4Heter,       & !Litter-microbial share of NH4-N uptake from underlying soil [-]
   AttenfNO3Heter                  => micflx%AttenfNO3Heter,       & !Litter-microbial share of NO3-N uptake from underlying soil [-]
   AttenfH2PO4Heter                => micflx%AttenfH2PO4Heter,     & !Litter-microbial share of H2PO4-P uptake from underlying soil [-]
   AttenfH1PO4Heter                => micflx%AttenfH1PO4Heter,     & !Litter-microbial share of HPO4-P uptake from underlying soil [-]
   RNH4imobilSoilHeter             => nmicf%RNH4imobilSoilHeter,   & !Net NH4-N transfer from nonband soil to microbes; positive immobilization
   RNO3imobilSoilHeter             => nmicf%RNO3imobilSoilHeter,   & !Net NO3-N transfer from nonband soil to microbes; positive immobilization
   RH2PO4imobilSoilHeter           => nmicf%RH2PO4imobilSoilHeter, & !Net H2PO4-P transfer from nonband soil to microbes; positive immobilization
   RNH4imobilBandHeter             => nmicf%RNH4imobilBandHeter,   & !Net NH4-N transfer from fertilizer-band soil to microbes; positive immobilization
   RNO3imobilBandHeter             => nmicf%RNO3imobilBandHeter,   & !Net NO3-N transfer from fertilizer-band soil to microbes; positive immobilization
   RH2PO4imobilBandHeter           => nmicf%RH2PO4imobilBandHeter, & !Net H2PO4-P transfer from fertilizer-band soil to microbes; positive immobilization
   RNH4imobilLitrHeter             => nmicf%RNH4imobilLitrHeter,   & !Net NH4-N transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
   RNO3imobilLitrHeter             => nmicf%RNO3imobilLitrHeter,   & !Net NO3-N transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
   RH2PO4imobilLitrHeter           => nmicf%RH2PO4imobilLitrHeter, & !Net H2PO4-P transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
   RH1PO4imobilSoilHeter           => nmicf%RH1PO4imobilSoilHeter, & !Net HPO4-P transfer from nonband soil to microbes; positive immobilization
   RH1PO4imobilBandHeter           => nmicf%RH1PO4imobilBandHeter, & !Net HPO4-P transfer from fertilizer-band soil to microbes; positive immobilization
   RH1PO4imobilLitrHeter           => nmicf%RH1PO4imobilLitrHeter, & !Net HPO4-P transfer from underlying soil accessed by litter microbes to microbes; positive immobilization
   rNCOMC                          => micpar%rNCOMC,               & !Target heterotrophic N:C ratios by compartment, guild and complex
   rPCOMC                          => micpar%rPCOMC,               & !Target heterotrophic P:C ratios by compartment, guild and complex
   litrm                           => micfor%litrm,                & !True for the surface litter layer
   VLNH4                           => micfor%VLNH4,                & !Nonband fraction for ammonium/ammonia pools and uptake capacity [-]
   VLNHB                           => micfor%VLNHB,                & !Fertilizer-band fraction for ammonium/ammonia pools and uptake capacity [-]
   VLNO3                           => micfor%VLNO3,                & !Nonband fraction for nitrate/nitrite pools and uptake capacity [-]
   VLNOB                           => micfor%VLNOB,                & !Fertilizer-band fraction for nitrate/nitrite pools and uptake capacity [-]
   VLPOB                           => micfor%VLPOB,                & !Fertilizer-band fraction for phosphate pools and uptake capacity [-]
   VLWatMicP                       => micfor%VLWatMicP,            & !Layer micropore water volume used for nutrient donor thresholds
   VOLWU                           => micfor%VOLWU,                & !Water volume in the soil beneath surface litter
   VLPO4                           => micfor%VLPO4,                & !Nonband fraction for phosphate pools and uptake capacity [-]
   ZNH4B                           => micstt%ZNH4B,                & !NH4-N pool in fertilizer-band soil
   ZNH4S                           => micstt%ZNH4S,                & !NH4-N pool in nonband soil
   ZNO3B                           => micstt%ZNO3B,                & !NO3-N pool in fertilizer-band soil
   ZNO3S                           => micstt%ZNO3S,                & !NO3-N pool in nonband soil
   ZNH4TU                          => micstt%ZNH4TU,               & !NH4-N pool in underlying soil, band plus nonband
   ZNO3TU                          => micstt%ZNO3TU,               & !NO3-N pool in underlying soil, band plus nonband
   H1P4TU                          => micstt%H1P4TU,               & !HPO4-P pool in underlying soil, band plus nonband
   H2P4TU                          => micstt%H2P4TU,               & !H2PO4-P pool in underlying soil, band plus nonband
   CNH4BU                          => micstt%CNH4BU,               & !Dissolved NH4-N concentration in underlying fertilizer-band soil
   CNH4SU                          => micstt%CNH4SU,               & !Dissolved NH4-N concentration in underlying nonband soil
   CH2P4U                          => micstt%CH2P4U,               & !Dissolved H2PO4-P concentration in underlying nonband soil
   CH2P4BU                         => micstt%CH2P4BU,              & !Dissolved H2PO4-P concentration in underlying fertilizer-band soil
   CH1P4U                          => micstt%CH1P4U,               & !Dissolved HPO4-P concentration in underlying nonband soil
   CH1P4BU                         => micstt%CH1P4BU,              & !Dissolved HPO4-P concentration in underlying fertilizer-band soil
   CH2P4                           => micstt%CH2P4,                & !Dissolved H2PO4-P concentration in nonband soil
   CH2P4B                          => micstt%CH2P4B,               & !Dissolved H2PO4-P concentration in fertilizer-band soil
   CNH4B                           => micstt%CNH4B,                & !Dissolved NH4-N concentration in fertilizer-band soil
   CNH4S                           => micstt%CNH4S,                & !Dissolved NH4-N concentration in nonband soil
   CH1P4                           => micstt%CH1P4,                & !Dissolved HPO4-P concentration in nonband soil
   CH1P4B                          => micstt%CH1P4B,               & !Dissolved HPO4-P concentration in fertilizer-band soil
   H1PO4                           => micstt%H1PO4,                & !HPO4-P pool in nonband soil
   H1POB                           => micstt%H1POB,                & !HPO4-P pool in fertilizer-band soil
   H2PO4                           => micstt%H2PO4,                & !H2PO4-P pool in nonband soil
   H2POB                           => micstt%H2POB,                & !H2PO4-P pool in fertilizer-band soil
   CNO3B                           => micstt%CNO3B,                & !Dissolved NO3-N concentration in fertilizer-band soil
   CNO3S                           => micstt%CNO3S,                & !Dissolved NO3-N concentration in nonband soil
   CNO3SU                          => micstt%CNO3SU,               & !Dissolved NO3-N concentration in underlying nonband soil
   CNO3BU                          => micstt%CNO3BU,               & !Dissolved NO3-N concentration in underlying fertilizer-band soil
   mBiomeHeter                     => micstt%mBiomeHeter,          & !C/N/P pools indexed by element, flattened guild/compartment and complex K
   RNO3DmndBandHeter               => micflx%RNO3DmndBandHeter,    & !Potential NO3-N uptake by heterotrophic guilds from fertilizer-band soil
   RNO3DmndSoilHeter               => micflx%RNO3DmndSoilHeter,    & !Potential NO3-N uptake by heterotrophic guilds from nonband soil
   RNH4DmndSoilHeter               => micflx%RNH4DmndSoilHeter,    & !Potential NH4-N uptake by heterotrophic guilds from nonband soil
   RNH4DmndBandHeter               => micflx%RNH4DmndBandHeter,    & !Potential NH4-N uptake by heterotrophic guilds from fertilizer-band soil
   RH2PO4DmndSoilHeter             => micflx%RH2PO4DmndSoilHeter,  & !Potential H2PO4-P uptake by heterotrophic guilds from nonband soil
   RH2PO4DmndBandHeter             => micflx%RH2PO4DmndBandHeter,  & !Potential H2PO4-P uptake by heterotrophic guilds from fertilizer-band soil
   RH1PO4DmndSoilHeter             => micflx%RH1PO4DmndSoilHeter,  & !Potential HPO4-P uptake by heterotrophic guilds from nonband soil
   RH1PO4DmndBandHeter             => micflx%RH1PO4DmndBandHeter,  & !Potential HPO4-P uptake by heterotrophic guilds from fertilizer-band soil
   RNH4DmndLitrHeter               => micflx%RNH4DmndLitrHeter,    & !Potential NH4-N uptake by heterotrophic guilds from underlying soil accessed by litter microbes
   RNO3DmndLitrHeter               => micflx%RNO3DmndLitrHeter,    & !Potential NO3-N uptake by heterotrophic guilds from underlying soil accessed by litter microbes
   RH2PO4DmndLitrHeter             => micflx%RH2PO4DmndLitrHeter,  & !Potential H2PO4-P uptake by heterotrophic guilds from underlying soil accessed by litter microbes
   RH1PO4DmndLitrHeter             => micflx%RH1PO4DmndLitrHeter,  & !Potential HPO4-P uptake by heterotrophic guilds from underlying soil accessed by litter microbes
   RNiDemand                       => micflx%RNiDemand          ,  & !Accumulated mineral-N requirement of heterotrophic biomass
   RPiDemand                       => micflx%RPiDemand          ,  & !Accumulated mineral-P requirement of heterotrophic biomass
   NetNH4Mineralize                => micflx%NetNH4Mineralize,     & !Net mineral N exchange (NH4 plus NO3); positive immobilization, negative release
   NetPO4Mineralize                => micflx%NetPO4Mineralize      & !Net phosphate exchange; positive immobilization, negative mineralization
  )
  call PrintInfo('beg '//subname)
  !     MINERALIZATION-IMMOBILIZATION OF NH4 IN SOIL FROM MICROBIAL
  !     C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
  !
  !     RINHP=NH4 mineralization (-ve) or immobilization (+ve) demand
  !     OMC,OMN=microbial nonstructural C,N
  !     rNCOMC=maximum microbial N:C ratio
  !     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
  !     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
  !     minimum NH4 concentration and Km for NH4 uptake
  !     RINHX=microbially limited NH4 demand
  !     BIOA=microbial surface area, OMA=active biomass
  !     GrowthEnvScalHeter=temp+water stress
  !     FNH4S,FNHBS=fractions of NH4 in non-band, band
  !     RNH4DmndSoilHeter,RNH4DmndBandHeter=substrate-unlimited NH4 mineraln-immobiln
  !     VOLW=water content
  !     ZNH4M,ZNHBM=NH4 not available for uptake in non-band, band
  !     FNH4X,FNB4X=fractions of biological NH4 demand in non-band, band
  !     RNH4imobilSoilHeter,RNH4imobilBandHeter=substrate-limited NH4 mineraln-immobiln in non-band, band
  !     NetNH4Mineralize=total NH4 net mineraln (-ve) or immobiln (+ve)
  ! update may be needed, May 17th, 2023, jyt.
  FNH4S     = VLNH4
  FNHBS     = VLNHB
  MID3      = micpar%get_micb_id(iLbiom_reserve,NGL)
  RINHP     = (mBiomeHeter(ielmc,MID3,K)*rNCOMC(iLbiom_reserve,NGL,K)-mBiomeHeter(ielmn,MID3,K))
  RNiDemand = RNiDemand+RINHP
  
  IF(RINHP.GT.0.0_r8)THEN
    !immobilization
    CNH4X                      = AZMAX1(CNH4S-Z4MN)
    CNH4Y                      = AZMAX1(CNH4B-Z4MN)
    RINHX                      = AMIN1(RINHP,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*Z4MX)
    RNH4DmndSoilHeter(NGL,K)   = FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
    RNH4DmndBandHeter(NGL,K)   = FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
    !Local donor thresholds use this layer; VOLWU is for underlying-soil access.
    ZNH4M                      = Z4MN*VLWatMicP*FNH4S
    ZNHBM                      = Z4MN*VLWatMicP*FNHBS
    RNH4imobilSoilHeter(NGL,K) = AMIN1(FNH4X*AZMAX1((ZNH4S-ZNH4M)),RNH4DmndSoilHeter(NGL,K))
    RNH4imobilBandHeter(NGL,K) = AMIN1(FNB4X*AZMAX1((ZNH4B-ZNHBM)),RNH4DmndBandHeter(NGL,K))
    !mineralization  (<0._r8)
  ELSE
    !mineralization
    RNH4DmndSoilHeter(NGL,K)   = 0.0_r8
    RNH4DmndBandHeter(NGL,K)   = 0.0_r8
    RNH4imobilSoilHeter(NGL,K) = RINHP*FNH4S
    RNH4imobilBandHeter(NGL,K) = RINHP*FNHBS
  ENDIF

  NetNH4Mineralize=NetNH4Mineralize+(RNH4imobilSoilHeter(NGL,K)+RNH4imobilBandHeter(NGL,K))
  !
  !     MINERALIZATION-IMMOBILIZATION OF NO3 IN SOIL FROM MICROBIAL
  !     C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
  !
  !     RINOP=NO3 immobilization (+ve) demand
  !     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
  !     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
  !     min NO3 concentration and Km for NO3 uptake
  !     RINOX=microbially limited NO3 demand
  !     BIOA=microbial surface area, OMA=active biomass
  !     GrowthEnvScalHeter=temp+water stress
  !     FNO3S,FNO3B=fractions of NO3 in non-band, band
  !     RNO3DmndSoilHeter,RNO3DmndBandHeter=substrate-unlimited NO3 immobiln
  !     VOLW=water content
  !     ZNO3M,ZNOBM=NO3 not available for uptake in non-band, band
  !     FNO3X,FNB3X=fractions of biological NO3 demand in non-band, band
  !     RNO3imobilSoilHeter,RNO3imobilBandHeter=substrate-limited NO3 immobiln in non-band, band
  !     NetNH4Mineralize=total net NH4+NO3 mineraln (-ve) or immobiln (+ve)
  !
  FNO3S = VLNO3
  FNO3B = VLNOB
  RINOP = AZMAX1(RINHP-RNH4imobilSoilHeter(NGL,K)-RNH4imobilBandHeter(NGL,K))
  !immobilization
  IF(RINOP.GT.0.0_r8)THEN
    CNO3X                      = AZMAX1(CNO3S-ZOMN)
    CNO3Y                      = AZMAX1(CNO3B-ZOMN)
    RINOX                      = AMIN1(RINOP,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*ZOMX)
    RNO3DmndSoilHeter(NGL,K)   = FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
    RNO3DmndBandHeter(NGL,K)   = FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
    ZNO3M                      = ZOMN*VLWatMicP*FNO3S
    ZNOBM                      = ZOMN*VLWatMicP*FNO3B
    RNO3imobilSoilHeter(NGL,K) = AMIN1(FNO3X*AZMAX1((ZNO3S-ZNO3M)),RNO3DmndSoilHeter(NGL,K))
    RNO3imobilBandHeter(NGL,K) = AMIN1(FNB3X*AZMAX1((ZNO3B-ZNOBM)),RNO3DmndBandHeter(NGL,K))
    !mineralization, 
  ELSE
    RNO3DmndSoilHeter(NGL,K)   = 0.0_r8
    RNO3DmndBandHeter(NGL,K)   = 0.0_r8
    RNO3imobilSoilHeter(NGL,K) = 0.0_r8
    RNO3imobilBandHeter(NGL,K) = 0.0_r8
  ENDIF
  NetNH4Mineralize=NetNH4Mineralize+(RNO3imobilSoilHeter(NGL,K)+RNO3imobilBandHeter(NGL,K))
  !
  !     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SOIL FROM MICROBIAL
  !     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
  !
  !     RIPOP=H2PO4 mineralization (-ve) or immobilization (+ve) demand
  !     OMC,OMP=microbial nonstructural C,P
  !     rPCOMC=maximum microbial P:C ratio
  !     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
  !     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
  !     min H2PO4 concentration and Km for H2PO4 uptake
  !     RIPOX=microbially limited H2PO4 demand
  !     BIOA=microbial surface area, OMA=active biomass
  !     GrowthEnvScalHeter=temp+water stress
  !     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
  !     RH2PO4DmndSoilHeter,RH2PO4DmndBandHeter=substrate-unlimited H2PO4 mineraln-immobiln
  !     H2POM,H2PBM=H2PO4 not available for uptake in non-band, band
  !     VOLW=water content
  !     FPO4X,FPOBX=fractions of biol H2PO4 demand in non-band, band
  !     RH2PO4imobilSoilHeter,RH2PO4imobilBandHeter=substrate-limited H2PO4 mineraln-immobn in non-band, band
  !     NetPO4Mineralize=total H2PO4 net mineraln (-ve) or immobiln (+ve)
  !
  FH2PS     = VLPO4
  FH2PB     = VLPOB
  MID3      = micpar%get_micb_id(iLbiom_reserve,NGL)
  RIPOP     = (mBiomeHeter(ielmc,MID3,K)*rPCOMC(iLbiom_reserve,NGL,K)-mBiomeHeter(ielmp,MID3,K))
  RPiDemand = RPiDemand+RIPOP
  !immobilization
  IF(RIPOP.GT.0.0_r8)THEN
    CH2PX                        = AZMAX1(CH2P4-HPMN)
    CH2PY                        = AZMAX1(CH2P4B-HPMN)
    RIPOX                        = AMIN1(RIPOP,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*HPMX)
    RH2PO4DmndSoilHeter(NGL,K)   = FH2PS*RIPOX*CH2PX/(CH2PX+HPKU)
    RH2PO4DmndBandHeter(NGL,K)   = FH2PB*RIPOX*CH2PY/(CH2PY+HPKU)
    H2POM                        = HPMN*VLWatMicP*FH2PS
    H2PBM                        = HPMN*VLWatMicP*FH2PB
    RH2PO4imobilSoilHeter(NGL,K) = AMIN1(FPO4X*AZMAX1((H2PO4-H2POM)),RH2PO4DmndSoilHeter(NGL,K))
    RH2PO4imobilBandHeter(NGL,K) = AMIN1(FPOBX*AZMAX1((H2POB-H2PBM)),RH2PO4DmndBandHeter(NGL,K))
  !mineralization  
  ELSE
    RH2PO4DmndSoilHeter(NGL,K)   = 0.0_r8
    RH2PO4DmndBandHeter(NGL,K)   = 0.0_r8
    RH2PO4imobilSoilHeter(NGL,K) = RIPOP*FH2PS
    RH2PO4imobilBandHeter(NGL,K) = RIPOP*FH2PB
  ENDIF
  NetPO4Mineralize=NetPO4Mineralize+(RH2PO4imobilSoilHeter(NGL,K)+RH2PO4imobilBandHeter(NGL,K))
  !
  !     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SOIL FROM MICROBIAL
  !     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
  !
  !     RIP1P=HPO4 mineralization (-ve) or immobilization (+ve) demand
  !     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
  !     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
  !     min HPO4 concentration and Km for HPO4 uptake
  !     RIP1X=microbially limited HPO4 demand
  !     BIOA=microbial surface area, OMA=active biomass
  !     GrowthEnvScalHeter=temp+water stress
  !     FH1PS,FH1PB=fractions of HPO4 in non-band, band
  !     RH1PO4DmndSoilHeter,RH1PO4DmndBandHeter=substrate-unlimited HPO4 mineraln-immobiln
  !     H1POM,H1PBM=HPO4 not available for uptake in non-band, band
  !     VOLW=water content
  !     FP14X,FP1BX=fractions of biol HPO4 demand in non-band, band
  !     RH1PO4imobilSoilHeter,RH1PO4imobilBandHeter=substrate-limited HPO4 mineraln-immobn in non-band, band uptake (> 0)
  !     NetPO4Mineralize=total H2PO4+HPO4 net mineraln (-ve) or immobiln (+ve)
  !
  FH1PS=VLPO4
  FH1PB=VLPOB
  ! why 0.1 here?  
  RIP1P=0.1_r8*AZMAX1(RIPOP-RH2PO4imobilSoilHeter(NGL,K)-RH2PO4imobilBandHeter(NGL,K))
  !immobilization
  IF(RIP1P.GT.0.0_r8)THEN
    CH1PX                        = AZMAX1(CH1P4-HPMN)
    CH1PY                        = AZMAX1(CH1P4B-HPMN)
    RIP1X                        = AMIN1(RIP1P,BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*HPMX)
    RH1PO4DmndSoilHeter(NGL,K)   = FH1PS*RIP1X*CH1PX/(CH1PX+HPKU)
    RH1PO4DmndBandHeter(NGL,K)   = FH1PB*RIP1X*CH1PY/(CH1PY+HPKU)
    H1POM                        = HPMN*VLWatMicP*FH1PS
    H1PBM                        = HPMN*VLWatMicP*FH1PB
    RH1PO4imobilSoilHeter(NGL,K) = AMIN1(FP14X*AZMAX1((H1PO4-H1POM)),RH1PO4DmndSoilHeter(NGL,K))
    RH1PO4imobilBandHeter(NGL,K) = AMIN1(FP1BX*AZMAX1((H1POB-H1PBM)),RH1PO4DmndBandHeter(NGL,K))
  !mineralization  
  ELSE
    RH1PO4DmndSoilHeter(NGL,K)   = 0.0_r8
    RH1PO4DmndBandHeter(NGL,K)   = 0.0_r8
    RH1PO4imobilSoilHeter(NGL,K) = 0.0_r8
    RH1PO4imobilBandHeter(NGL,K) = 0.0_r8
  ENDIF
  NetPO4Mineralize=NetPO4Mineralize+(RH1PO4imobilSoilHeter(NGL,K)+RH1PO4imobilBandHeter(NGL,K))
  !
  !     MINERALIZATION-IMMOBILIZATION OF NH4 IN SURFACE RESIDUE FROM
  !     MICROBIAL C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL
  !     ZONES OF SOIL SURFACE
  !
  !     RINHPR=NH4 mineralization (-ve) or immobilization (+ve) demand
  !     NU=surface layer number
  !     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
  !     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
  !     minimum NH4 concentration and Km for NH4 uptake
  !     BIOA=microbial surface area, OMA=active biomass
  !     GrowthEnvScalHeter=temp+water stress
  !     FNH4S,FNHBS=fractions of NH4 in non-band, band
  !     RNH4DmndLitrHeter=substrate-unlimited NH4 mineraln-immobiln
  !     VOLW=water content
  !     ZNH4M=NH4 not available for uptake
  !     AttenfNH4Heter=fractions of biological NH4 demand
  !     RNH4imobilLitrHeter=substrate-limited NH4 mineraln-immobiln
  !     NetNH4Mineralize=total NH4 net mineraln (-ve) or immobiln (+ve)
  !when there is not sufficient uptake from the litter layer, take it from soil
  IF(litrm)THEN
    RINHPR=RINHP-RNH4imobilSoilHeter(NGL,K)-RNO3imobilSoilHeter(NGL,K) &
      -RNH4imobilBandHeter(NGL,K)-RNO3imobilBandHeter(NGL,K)
    !immobilization by tap into the top soil layer

    call LitterSoilNutrientUptake(1,RINHPR, &
        BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*Z4MX, &
        Z4KU,Z4MN,AttenfNH4Heter(NGL,K), &
        micflx%RNH4DmndLitrBandHeterPrev(NGL,K),nmics%FracOMActHeter(NGL,K),micfor,LitrPotential,LitrUptake)
    RNH4DmndLitrHeter(NGL,K)=LitrPotential(1)
    micflx%RNH4DmndLitrBandHeter(NGL,K)=LitrPotential(2)
    RNH4imobilLitrHeter(NGL,K)=SUM(LitrUptake)
    micflx%tRNH4MicrbImobilSoil=micflx%tRNH4MicrbImobilSoil+LitrUptake(1)
    micflx%tRNH4MicrbImobilBand=micflx%tRNH4MicrbImobilBand+LitrUptake(2)
    NetNH4Mineralize=NetNH4Mineralize+RNH4imobilLitrHeter(NGL,K)
    !
    !     MINERALIZATION-IMMOBILIZATION OF NO3 IN SURFACE RESIDUE FROM
    !     MICROBIAL C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL
    !     ZONES OF SOIL SURFACE
    !
    !     RINOPR=NH4 mineralization (-ve) or immobilization (+ve) demand
    !     NU=surface layer number
    !     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
    !     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
    !     minimum NO3 concentration and Km for NO3 uptake
    !     RNO3DmndLitrHeter=microbially limited NO3 demand
    !     BIOA=microbial surface area, OMA=active biomass
    !     GrowthEnvScalHeter=temp+water stress
    !     FNO3S,FNO3B=fractions of NO3 in non-band, band
    !     RNO3imobilLitrHeter=substrate-unlimited NO3 immobiln
    !     VOLW=water content
    !     ZNO3M=NO3 not available for uptake
    !     AttenfNO3Heter=fraction of biological NO3 demand
    !     RNO3imobilLitrHeter=substrate-limited NO3 immobiln
    !     NetNH4Mineralize=total NH4+NO3 net mineraln (-ve) or immobiln (+ve)
    !
    RINOPR=AZMAX1(RINHPR-RNH4imobilLitrHeter(NGL,K))
    !immobilization by tapping into the top soil layer

    call LitterSoilNutrientUptake(2,RINOPR, &
        BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*ZOMX, &
        ZOKU,ZOMN,AttenfNO3Heter(NGL,K), &
        micflx%RNO3DmndLitrBandHeterPrev(NGL,K),nmics%FracOMActHeter(NGL,K),micfor,LitrPotential,LitrUptake)
    RNO3DmndLitrHeter(NGL,K)=LitrPotential(1)
    micflx%RNO3DmndLitrBandHeter(NGL,K)=LitrPotential(2)
    RNO3imobilLitrHeter(NGL,K)=SUM(LitrUptake)
    micflx%tRNO3MicrbImobilSoil=micflx%tRNO3MicrbImobilSoil+LitrUptake(1)
    micflx%tRNO3MicrbImobilBand=micflx%tRNO3MicrbImobilBand+LitrUptake(2)
    NetNH4Mineralize=NetNH4Mineralize+RNO3imobilLitrHeter(NGL,K)
    !
    !     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SURFACE RESIDUE FROM
    !     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
    !     ZONES OF SOIL SURFACE
    !
    !     RIPOPR=H2PO4 mineralization (-ve) or immobilization (+ve) demand
    !     NU=surface layer number
    !     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
    !     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
    !     minimum H2PO4 concentration and Km for H2PO4 uptake
    !     RH2PO4DmndLitrHeter=microbially limited H2PO4 demand
    !     BIOA=microbial surface area, OMA=active biomass
    !     GrowthEnvScalHeter=temp+water stress
    !     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
    !     RH2PO4DmndLitrHeter=substrate-unlimited H2PO4 mineraln-immobiln
    !     VOLW=water content
    !     H2P4M=H2PO4 not available for uptake
    !     AttenfH2PO4Heter=fractions of biological H2PO4 demand
    !     RH2PO4imobilLitrHeter=substrate-limited H2PO4 mineraln-immobiln
    !     NetPO4Mineralize=total H2PO4 net mineraln (-ve) or immobiln (+ve)
    !
    !Subtract all P already exchanged with litter before tapping topsoil.
    RIPOPR=RIPOP-RH2PO4imobilSoilHeter(NGL,K)-RH2PO4imobilBandHeter(NGL,K) &
      -RH1PO4imobilSoilHeter(NGL,K)-RH1PO4imobilBandHeter(NGL,K)
    !immobilization by tapping into top soil layer

    call LitterSoilNutrientUptake(3,RIPOPR, &
        BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*HPMX, &
        HPKU,HPMN,AttenfH2PO4Heter(NGL,K), &
        micflx%RH2PO4DmndLitrBandHeterPrev(NGL,K),nmics%FracOMActHeter(NGL,K),micfor,LitrPotential,LitrUptake)
    RH2PO4DmndLitrHeter(NGL,K)=LitrPotential(1)
    micflx%RH2PO4DmndLitrBandHeter(NGL,K)=LitrPotential(2)
    RH2PO4imobilLitrHeter(NGL,K)=SUM(LitrUptake)
    micflx%tRH2PO4MicrbImobilSoil=micflx%tRH2PO4MicrbImobilSoil+LitrUptake(1)
    micflx%tRH2PO4MicrbImobilBand=micflx%tRH2PO4MicrbImobilBand+LitrUptake(2)
    NetPO4Mineralize=NetPO4Mineralize+RH2PO4imobilLitrHeter(NGL,K)
    !
    !     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SURFACE RESIDUE FROM
    !     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
    !     ZONES OF SOIL SURFACE
    !
    !     RIP1PR=HPO4 mineralization (-ve) or immobilization (+ve) demand
    !     NU=surface layer number
    !     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
    !     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
    !     minimum HPO4 concentration and Km for HPO4 uptake
    !     RH1PO4DmndLitrHeter=microbially limited HPO4 demand
    !     BIOA=microbial surface area, OMA=active biomass
    !     GrowthEnvScalHeter=temp+water stress
    !     FH1PS,FH1PB=fractions of HPO4 in non-band, band
    !     RH1PO4DmndLitrHeter=substrate-unlimited HPO4 mineraln-immobiln
    !     VOLW=water content
    !     H1P4M=HPO4 not available for uptake
    !     AttenfH1PO4Heter=fraction of biological HPO4 demand
    !     RH1PO4imobilLitrHeter=substrate-limited HPO4 minereraln-immobiln
    !     NetPO4Mineralize=total HPO4 net mineraln (-ve) or immobiln (+ve)
    !
    FH1PS=VLPO4
    FH1PB=VLPOB
    RIP1PR=0.1_r8*AZMAX1(RIPOPR-RH2PO4imobilLitrHeter(NGL,K))
    !immobilization

    call LitterSoilNutrientUptake(4,RIP1PR, &
        BIOA*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)*HPMX, &
        HPKU,HPMN,AttenfH1PO4Heter(NGL,K), &
        micflx%RH1PO4DmndLitrBandHeterPrev(NGL,K),nmics%FracOMActHeter(NGL,K),micfor,LitrPotential,LitrUptake)
    RH1PO4DmndLitrHeter(NGL,K)=LitrPotential(1)
    micflx%RH1PO4DmndLitrBandHeter(NGL,K)=LitrPotential(2)
    RH1PO4imobilLitrHeter(NGL,K)=SUM(LitrUptake)
    micflx%tRH1PO4MicrbImobilSoil=micflx%tRH1PO4MicrbImobilSoil+LitrUptake(1)
    micflx%tRH1PO4MicrbImobilBand=micflx%tRH1PO4MicrbImobilBand+LitrUptake(2)
    NetPO4Mineralize=NetPO4Mineralize+RH1PO4imobilLitrHeter(NGL,K)
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine BiomassMineralization

!------------------------------------------------------------------------------------------

  subroutine GatherHeterotrophRespiration(I,J,NGL,N,K,RMaintDefcitcitHeter,micfor,micstt,nmicf,nmics,micflx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NGL,N,K
  real(r8), intent(out) :: RMaintDefcitcitHeter   !deficit for maintenance respiraiton, [gC d-2 h-1]
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(micfluxtype), intent(inout) :: micflx    
  character(len=*), parameter :: subname='GatherHeterotrophRespiration'
  integer :: MID3,MID1
  REAL(R8) :: FPH,RMOMX
  real(r8) :: RGN2P,DenitResp,FracDenitResp4Maint
  real(r8) :: PhotoGrowthC,OrganicGrowthResp,PhotoNFixResp,OrganicNFixResp,PhotoResp
!     begin_execution
  associate(                                                   &
    RCO2FixCyano           => nmicf%RCO2FixCyano,            & !Available photosynthetic C, including maintenance supply
    fPhotoR                => nmicf%fPhotoR,                 & !Photosynthate fraction of actual primary respiration
    RCO2ProdHeter          => nmicf%RCO2ProdHeter,           & !Net CO2 release, including the photosynthetic uptake credit
    RO2UptkHeter           => nmicf%RO2UptkHeter,            & !Signed net external O2 uptake, negative for release
    RO2Uptk4RespHeter      => nmicf%RO2Uptk4RespHeter,       & !Actual gross respiratory O2 consumption
    RO2Dmnd4RespHeter      => nmicf%RO2Dmnd4RespHeter,       & !Potential gross respiratory O2 consumption
    REcoUptkSoilO2M        => micflx%REcoUptkSoilO2M,       & !Signed external O2 exchange accumulated by transport substep
    TempMaintRHeter         => nmics%TempMaintRHeter,          & !Guild temperature multiplier on heterotrophic maintenance [-]; not referenced here
    OMN2                    => nmics%OMN2,                     & !Active structural heterotrophic N by guild and complex K; not referenced here
    Resp4NFixHeter          => nmicf%Resp4NFixHeter,           & !Respiration-C cost of heterotrophic N2 fixation
    RNOxDOCReduxRespDenitLim => nmicf%RNOxDOCReduxRespDenitLim, & !DOC-C respiration supported by denitrification
    RNOxAcetReduxRespDenitLim => nmicf%RNOxAcetReduxRespDenitLim, & !Acetate-C respiration supported by denitrification
    RespGrossHeter          => nmicf%RespGrossHeter,           & !Gross respiration C equivalent from the primary heterotrophic pathway
    RMaintRespHeter         => nmicf%RMaintRespHeter,          & !Total hourly heterotrophic maintenance-C demand by guild and complex
    RMaintDmndHeter         => nmicf%RMaintDmndHeter,          & !Maintenance-C demand by live compartment, heterotrophic guild and complex; not referenced here
    RNH4imobilSoilHeter     => nmicf%RNH4imobilSoilHeter,      & !Net NH4-N transfer from nonband soil to microbes; positive immobilization; not referenced here
    RNO3imobilSoilHeter     => nmicf%RNO3imobilSoilHeter,      & !Net NO3-N transfer from nonband soil to microbes; positive immobilization; not referenced here
    RGrowthRespHeter        => micflx%RGrowthRespHeter    ,    & !Post-maintenance capacity; cyanobacteria also include residual photosynthetic C
    RN2FixHeter             => nmicf%RN2FixHeter,              & !Heterotrophic guild N2-N fixation into biomass
    rNCOMC                  => micpar%rNCOMC,                  & !Target heterotrophic N:C ratios by compartment, guild and complex
    rPCOMC                  => micpar%rPCOMC,                  & !Target heterotrophic P:C ratios by compartment, guild and complex; not referenced here
    pH                      => micfor%pH,                      & !Layer pH used in microbial rate and nutrient-speciation responses; not referenced here
    ZEROS                   => micfor%ZEROS,                   & !Small mass or flux threshold used by the routine
    mid_HeterMixtCynoBacter => micpar%mid_HeterMixtCynoBacter, & !Functional-group identifier for mixotrophic cyanobacteria
    mid_HeterAerobN2Fixer   => micpar%mid_HeterAerobN2Fixer,   & !Functional-group identifier for aerobic heterotrophic N2 fixers
    mid_HeterAnaerobN2Fixer => micpar%mid_HeterAnaerobN2Fixer, & !Functional-group identifier for anaerobic heterotrophic N2 fixers
    CZ2GS                   => micstt%CZ2GS,                   & !Dissolved N2-N concentration for nitrogen fixation
    mBiomeHeter             => micstt%mBiomeHeter              & !C/N/P pools indexed by element, flattened guild/compartment and complex K
  )
  call PrintInfo('beg '//subname)
  !     pH EFFECT ON MAINTENANCE RESPIRATION
  !
  !     FPH=pH effect on maintenance respiration
  !     RMOM=specific maintenance respiration rate
  !     TempMaintRHeter=temperature effect on maintenance respiration
  !     OMN=microbial N biomass
  !
  !
  RGrowthRespHeter(NGL,K) = AZMAX1(RespGrossHeter(NGL,K)-RMaintRespHeter(NGL,K))
  IF(N.EQ.mid_HeterMixtCynoBacter)THEN
    !Photosynthetic growth is the residual after maintenance. It can also
    !support N2 fixation; retain that capacity separately from actual CO2 loss.
    PhotoResp=RespGrossHeter(NGL,K)*fPhotoR(NGL,K)
    PhotoGrowthC=AZMAX1(RCO2FixCyano(NGL,K)-PhotoResp)
    OrganicGrowthResp=RGrowthRespHeter(NGL,K)
    RGrowthRespHeter(NGL,K)=OrganicGrowthResp+PhotoGrowthC
  ENDIF
  !
  !Keep primary growth respiration separate; reserve denitrification energy
  !for any remaining maintenance before its growth uptake is calculated below.
  !
  DenitResp=RNOxDOCReduxRespDenitLim(NGL,K)+RNOxAcetReduxRespDenitLim(NGL,K)
  call ReserveDenitrifMaintenance(RMaintRespHeter(NGL,K),RespGrossHeter(NGL,K), &
    DenitResp,EO2X,ENOX,FracDenitResp4Maint,RMaintDefcitcitHeter)
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
  !     RGrowthRespHeter=growth respiration
  !     Resp4NFixHeter=respiration for N2 fixation
  !     CZ2GS=aqueous N2 concentration
  !     ZFKM=Km for N2 uptake
  !     OMGR*OMC(3,NGL,N,K)=nonstructural C limitation to Resp4NFixHeter
  !     RN2FixHeter=N2 fixation rate, [gN d-2 h-1]
  !
  IF(N.EQ.mid_HeterAerobN2Fixer .OR. N.EQ.mid_HeterAnaerobN2Fixer .or. N.eq.mid_HeterMixtCynoBacter)THEN
    MID3  = micpar%get_micb_id(iLbiom_reserve,NGL)
    RGN2P = AZMAX1(mBiomeHeter(ielmc,MID3,K)*rNCOMC(iLbiom_reserve,NGL,K)-mBiomeHeter(ielmn,MID3,K))/EN2F(N)
    IF(RGrowthRespHeter(NGL,K).GT.ZEROS)THEN
      Resp4NFixHeter(NGL,K)=AMIN1(RGrowthRespHeter(NGL,K)*RGN2P/(RGrowthRespHeter(NGL,K)+RGN2P) &
        *CZ2GS/(CZ2GS+ZFKM),OMGR*mBiomeHeter(ielmc,MID3,K))
      IF(N.EQ.mid_HeterMixtCynoBacter)THEN
        !Organic growth respiration can be redirected to N2 fixation without
        !respiring its C twice. Additional photosynthate respiration uses only
        !C and internally produced O2 left after primary metabolism.
        OrganicNFixResp=MIN(OrganicGrowthResp,Resp4NFixHeter(NGL,K))
        PhotoNFixResp=MIN(PhotoGrowthC,AZMAX1(-RO2UptkHeter(NGL,K))/2.667_r8, &
          AZMAX1(Resp4NFixHeter(NGL,K)-OrganicNFixResp))
        Resp4NFixHeter(NGL,K)=OrganicNFixResp+PhotoNFixResp
        RespGrossHeter(NGL,K)=RespGrossHeter(NGL,K)+PhotoNFixResp
        IF(RespGrossHeter(NGL,K).GT.0._r8) &
          fPhotoR(NGL,K)=(PhotoResp+PhotoNFixResp)/RespGrossHeter(NGL,K)
        RCO2ProdHeter(NGL,K)=RCO2ProdHeter(NGL,K)+PhotoNFixResp
        RO2UptkHeter(NGL,K)=RO2UptkHeter(NGL,K)+2.667_r8*PhotoNFixResp
        RO2Uptk4RespHeter(NGL,K)=RO2Uptk4RespHeter(NGL,K)+2.667_r8*PhotoNFixResp
        RO2Dmnd4RespHeter(NGL,K)=RO2Dmnd4RespHeter(NGL,K)+2.667_r8*PhotoNFixResp
        REcoUptkSoilO2M(1:NPH)=REcoUptkSoilO2M(1:NPH)+2.667_r8*PhotoNFixResp/REAL(NPH,r8)
      ENDIF
      RN2FixHeter(NGL,K)=Resp4NFixHeter(NGL,K)*EN2F(N)  
    ELSE
      Resp4NFixHeter(NGL,K) = 0.0_r8
      RN2FixHeter(NGL,K)    = 0._r8
    ENDIF    
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine GatherHeterotrophRespiration
!------------------------------------------------------------------------------------------

  subroutine GatherHetertrophAnabolicFlux(I,J,NGL,N,K,RMaintDefcitcitHeter,&
    SPOMK,micfor,micstt,nmicf,nmics,ncplxf,ncplxs,micflx)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NGL,N,K
  real(r8), intent(in) :: RMaintDefcitcitHeter  !maintenance deficit, [gC d-2 h-1]
  real(r8), intent(in) :: SPOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_Flux_type), intent(inout) :: ncplxf
  type(OMCplx_State_type), intent(inout) :: ncplxs
  type(micfluxtype), intent(inout) :: micflx  
  character(len=*), parameter :: subname='GatherHetertrophAnabolicFlux'

  integer :: M,MID1,MID3,MID,NE
  real(r8) :: C3C,CNC,CPC
  real(r8) :: CCC,CXC,RCCE(NumPlantChemElms)  

  real(r8) :: RCCC,RCCN,RCCP  
  real(r8) :: CGOMX,CGOMD,AGOMD
  real(r8) :: DenitResp,FracDenitResp4Maint,MaintDOC,MaintAcet
  real(r8) :: OrganicResp,OrganicMaint,OrganicGrowthResp,OrganicNFixResp
  real(r8) :: PrimaryDOCUptake,PrimaryAcetUptake
  real(r8) :: DOCRespFrac,AcetRespFrac
  real(r8) :: CGOXC              !DOC+acetate uptake flux, [gC d-2 h-1]
  real(r8) :: CGOMZ
  real(r8) :: SPOMX
  real(r8) :: FracMaintDeficit
  real(r8) :: AvailableBiomass(NumPlantChemElms)
!     begin_execution
  associate(                                                                    &
    rCNBiomeActHeter                 => nmics%rCNBiomeActHeter,                 & !Active heterotrophic biomass nutrient:C ratios (N:C and P:C)
    GrowthEnvScalHeter               => nmics%GrowthEnvScalHeter,               & !Temperature and water-potential multiplier on heterotrophic growth [-]
    FCN                              => nmics%FCN,                              & !N-status multiplier on heterotrophic activity by guild and complex [-]
    FCP                              => nmics%FCP,                              & !P-status multiplier on heterotrophic activity by guild and complex [-]
    FracHeterBiomOfActK              => nmics%FracHeterBiomOfActK,              & !Guild fraction of active heterotrophic biomass in complex K [-]
    DOMUptk4GrothHeter               => nmicf%DOMUptk4GrothHeter,               & !Guild elemental uptake; C includes DOC, acetate and cyanobacterial CO2 fixation
    RMetabDOCUptkHeter               => nmicf%RMetabDOCUptkHeter,               & !Total DOC-C uptake for primary metabolism and denitrification
    RMaintRespHeter                  => nmicf%RMaintRespHeter,                  & !Total hourly heterotrophic maintenance-C demand by guild and complex
    FGOCP                            => nmicf%FGOCP,                            & !DOC-supported fraction of total primary guild respiration [-]
    FGOAP                            => nmicf%FGOAP,                            & !Acetate-supported fraction of total primary guild respiration [-]
    RMetabAcetUptkHeter              => nmicf%RMetabAcetUptkHeter,              & !Total acetate-C uptake for primary metabolism and denitrification
    ECHZHeter                        => nmicf%ECHZHeter,                        & !Guild respiration fraction used to convert growth respiration to C uptake [-]
    NonstX2stBiomHeter               => nmicf%NonstX2stBiomHeter,               & !C/N/P transfer from reserves into kinetic/structural biomass by guild and complex
    RespGrossHeter                   => nmicf%RespGrossHeter,                   & !Gross respiration C equivalent from the primary heterotrophic pathway
    RCO2FixCyano                     => nmicf%RCO2FixCyano,                     & !Photosynthetic CO2-C fixation by cyanobacterial guild and complex
    RGrowthRespHeter                 => micflx%RGrowthRespHeter        ,        & !Post-maintenance capacity; cyanobacteria also include residual photosynthetic C
    RNOxDOCReduxRespDenitLim         => nmicf%RNOxDOCReduxRespDenitLim,         & !DOC-C respiration supported by realized denitrification
    RNOxAcetReduxRespDenitLim        => nmicf%RNOxAcetReduxRespDenitLim,        & !Acetate-C respiration supported by realized denitrification
    RMaintDmndHeter                  => nmicf%RMaintDmndHeter,                  & !Maintenance-C demand by live compartment, heterotrophic guild and complex
    fPhotoR                          => nmicf%fPhotoR,                          & !Photosynthate-supported fraction of cyanobacterial gross respiration [-]
    RkillLitfalOMHeter               => nmicf%RkillLitfalOMHeter,               & !Unrecycled ordinary-mortality C/N/P from heterotrophic biomass
    RkillLitrfal2HumOMHeter          => nmicf%RkillLitrfal2HumOMHeter,          & !Ordinary-mortality C/N/P routed to humified material from heterotrophic biomass
    RkillLitrfal2ResduOMHeter        => nmicf%RkillLitrfal2ResduOMHeter,        & !Ordinary-mortality C/N/P routed to microbial residue from heterotrophic biomass
    RCCMEHeter                       => nmicf%RCCMEHeter,                       & !Heterotrophic mortality/starvation residue C/N/P summed for each receiving complex
    RMaintDefLitrfal2HumOMHeter      => nmicf%RMaintDefLitrfal2HumOMHeter,      & !Starvation-derived C/N/P routed to humified material from heterotrophic biomass
    RMaintDefLitrfal2ResduOMHeter    => nmicf%RMaintDefLitrfal2ResduOMHeter,    & !Starvation-derived C/N/P routed to microbial residue from heterotrophic biomass
    RMaintDefcitLitrfalOMHeter       => nmicf%RMaintDefcitLitrfalOMHeter,       & !Unrecycled starvation-derived C/N/P from heterotrophic biomass
    RKillOMHeter                     => nmicf%RKillOMHeter,                     & !Ordinary mortality C/N/P withdrawal from heterotrophic biomass
    RkillRecycOMHeter                => nmicf%RkillRecycOMHeter,                & !Ordinary-mortality C/N/P recycled to reserves from heterotrophic biomass
    RMaintDefcitKillOMHeter          => nmicf%RMaintDefcitKillOMHeter,          & !Maintenance-starvation C/N/P withdrawal from heterotrophic biomass
    RMaintDefcitRecycOMHeter         => nmicf%RMaintDefcitRecycOMHeter,         & !Starvation recycling: C respired, N/P returned to reserves from heterotrophic biomass
    Resp4NFixHeter                   => nmicf%Resp4NFixHeter,                   & !Respiration-C cost of heterotrophic N2 fixation
    TDOMUptkHeter                    => ncplxf%TDOMUptkHeter,                   & !Total heterotrophic DOM uptake by species and complex K; not referenced here
    rCNDOM                           => ncplxs%rCNDOM,                          & !DON:DOC ratio in each substrate complex [-]
    rCPDOM                           => ncplxs%rCPDOM,                          & !DOP:DOC ratio in each substrate complex [-]
    rNCOMC                           => micpar%rNCOMC,                          & !Target heterotrophic N:C ratios by compartment, guild and complex
    rPCOMC                           => micpar%rPCOMC,                          & !Target heterotrophic P:C ratios by compartment, guild and complex
    FL                               => micpar%FL,                              & !Target fractions of active biomass in kinetic and structural compartments [-]
    EHUM                             => micstt%EHUM,                            & !Fraction of microbial litterfall routed to humified organic matter [-]
    mBiomeHeter                      => micstt%mBiomeHeter,                     & !C/N/P pools indexed by element, flattened guild/compartment and complex K
    DOM                              => micstt%DOM,                             & !Dissolved organic pools by species (DOC, DON, DOP, acetate) and complex K
    CDOMuptk1                        => micflx%CDOMuptk1,                       & !Total DOC plus acetate C uptake assigned to the primary growth pathway
    CDOMuptk2                        => micflx%CDOMuptk2,                       & !Accumulated DOC-C uptake supporting denitrification (excludes acetate uptake)
    tROMT                            => micflx%tROMT,                           & !Accumulated heterotrophic maintenance-C demand
    tGROMO                           => micflx%tGROMO,                          & !Accumulated realized gross heterotrophic respiration C
    ZEROS                            => micfor%ZEROS,                           & !Small mass or flux threshold used by the routine
    ZERO                             => micfor%ZERO                             & !Small dimensionless or concentration threshold used by the routine
  )
  call PrintInfo('beg '//subname)
  !     DOC, DON, DOP AND ACETATE UPTAKE DRIVEN BY GROWTH RESPIRATION
  !     FROM O2, NOX AND C REDUCTION
  !
  !     CGOMX=DOC+acetate uptake from aerobic growth respiration
  !     CGOMD=DOC  uptake from denitrifier growth respiration
  !     RMaintRespHeter=maintenance respiration
  !     RespGrossHeter=total respiration
  !     RNOxDOCReduxRespDenitLim=respiration generated by denitrifcation
  !     Resp4NFixHeter=respiration for N2 fixation
  !     ECHZHeter,ENOX=growth respiration efficiencies for O2, NOx reduction
  !     DOMUptk4GrothHeter,RMetabDOCUptkHeter,RMetabAcetUptkHeter=total DOC+acetate, DOC, acetate uptake(heterotrophs
  !     DOMUptk4GrothHeter=total CO2,CH4 uptake (autotrophs)
  !     DOMUptk4GrothHeter,DOMUptk4GrothHeter=DON, DOP uptake
  !     FracHeterBiomOfActK=faction of OMActHeterin total OMA
  !     rCNDOM,rCPDOM=DON/DOC, DOP/DOC
  !     FCN,FCP=limitation from N,P
  ! gross respiration equals to maintenance+respiraiton for N-fixation + growth respiraiton
  if(N.eq.micpar%mid_HeterMixtCynoBacter)THEN
    !Photosynthetic C is a separate input, used first for maintenance.
    !Only organic respiration is converted into DOC/acetate uptake. The
    !N2-fixation diversion consumes organic growth respiration first.
    OrganicResp=AZMAX1(RespGrossHeter(NGL,K)*(1._r8-fPhotoR(NGL,K)))
    OrganicMaint=MIN(OrganicResp,AZMAX1(RMaintRespHeter(NGL,K) &
      -RespGrossHeter(NGL,K)*fPhotoR(NGL,K)))
    OrganicGrowthResp=AZMAX1(OrganicResp-OrganicMaint)
    OrganicNFixResp=MIN(OrganicGrowthResp,Resp4NFixHeter(NGL,K))
  ELSE
    OrganicMaint=AMIN1(RMaintRespHeter(NGL,K),RespGrossHeter(NGL,K))
    OrganicGrowthResp=RGrowthRespHeter(NGL,K)
    OrganicNFixResp=Resp4NFixHeter(NGL,K)
  endif

  !Match each donor's respiration fraction to its catabolic allocation cap.
  !Fermenters and acetotrophic methanogens retain their single-donor yield.
  DOCRespFrac=ECHZHeter(NGL,K)
  AcetRespFrac=ECHZHeter(NGL,K)
  IF(N.EQ.micpar%mid_HeterAerobBacter)THEN
    DOCRespFrac=EO2X
    AcetRespFrac=EO2A
  ELSEIF(N.EQ.micpar%mid_Facult_DenitBacter)THEN
    DOCRespFrac=EO2D
    AcetRespFrac=EO2A
  ELSEIF(N.EQ.micpar%mid_Aerob_Fungi)THEN
    DOCRespFrac=EO2G
    AcetRespFrac=EO2A
  ELSEIF(N.EQ.micpar%mid_HeterAerobN2Fixer .OR. N.EQ.micpar%mid_HeterMixtCynoBacter)THEN
    DOCRespFrac=ENFX
    AcetRespFrac=EO2A
  ENDIF
  !Respiration shares are not uptake shares when donor yields differ.
  !Maintenance and N2-fixation allocations earn no growth-associated uptake.
  PrimaryDOCUptake=FGOCP(NGL,K)*(OrganicMaint+OrganicNFixResp &
    +(OrganicGrowthResp-OrganicNFixResp)/DOCRespFrac)
  PrimaryAcetUptake=FGOAP(NGL,K)*(OrganicMaint+OrganicNFixResp &
    +(OrganicGrowthResp-OrganicNFixResp)/AcetRespFrac)
  CGOMX=PrimaryDOCUptake+PrimaryAcetUptake

  DenitResp=RNOxDOCReduxRespDenitLim(NGL,K)+RNOxAcetReduxRespDenitLim(NGL,K)
  call ReserveDenitrifMaintenance(RMaintRespHeter(NGL,K),RespGrossHeter(NGL,K), &
    DenitResp,EO2X,ENOX,FracDenitResp4Maint)
  !
  !Reserve maintenance proportionally from DOC and acetate respiration. Only
  !the remaining respiration earns growth-associated substrate uptake.
  MaintDOC=RNOxDOCReduxRespDenitLim(NGL,K)*FracDenitResp4Maint
  MaintAcet = RNOxAcetReduxRespDenitLim(NGL,K)*FracDenitResp4Maint
  CGOMD     = MaintDOC+(RNOxDOCReduxRespDenitLim(NGL,K)-MaintDOC)/ENOX
  AGOMD     = MaintAcet+(RNOxAcetReduxRespDenitLim(NGL,K)-MaintAcet)/ENOX
  CDOMuptk1 = CDOMuptk1+CGOMX !DOC used for growth
  CDOMuptk2 = CDOMuptk2+CGOMD !DOC used for denitrifcation
  tROMT     = tROMT+RMaintRespHeter(NGL,K)
  tGROMO    = tGROMO+RespGrossHeter(NGL,K)

  !Total organic C uptake includes acetate used for denitrification.
  DOMUptk4GrothHeter(ielmc,NGL,K) = CGOMX+CGOMD+AGOMD
  !Fixation supplies maintenance as well as growth; never drop it when the
  !maintenance demand exhausts available respiratory capacity.
  if(N.eq.micpar%mid_HeterMixtCynoBacter)then
    DOMUptk4GrothHeter(ielmc,NGL,K) = DOMUptk4GrothHeter(ielmc,NGL,K)+RCO2FixCyano(NGL,K)
  endif
  RMetabDOCUptkHeter(NGL,K)  = PrimaryDOCUptake+CGOMD
  RMetabAcetUptkHeter(NGL,K) = PrimaryAcetUptake+AGOMD
  CGOXC                      = RMetabDOCUptkHeter(NGL,K)+RMetabAcetUptkHeter(NGL,K)

  !obtain organic nutrient uptake
  DOMUptk4GrothHeter(ielmn,NGL,K)=AZMAX1(AMIN1(DOM(idom_don,K)*FracHeterBiomOfActK(NGL,K),CGOXC*rCNDOM(K)/FCN(NGL,K)))
  DOMUptk4GrothHeter(ielmp,NGL,K)=AZMAX1(AMIN1(DOM(idom_dop,K)*FracHeterBiomOfActK(NGL,K),CGOXC*rCPDOM(K)/FCP(NGL,K)))

  !
  !     TRANSFER UPTAKEN C,N,P FROM STORAGE/nonstructural TO ACTIVE BIOMASS
  !
  !     OMC,OMN,OMP=nonstructural C,N,P
  !     CCC,CNC,CPC=C:N:P ratios used to calculate C,N,P recycling
  !     rNCOMC,rPCOMC=maximum microbial N:C, P:C ratios
  !     RCCC,RCCN,RCCP=C,N,P recycling fractions
  !     RCCZ,RCCY=min, max C recycling fractions
  !     RCCX,RCCQ=max N,P recycling fractions
  !
  !
  !     MICROBIAL ASSIMILATION OF NONSTRUCTURAL C,N,P
  !
  !     CGOMZ=transfer from nonstructural to structural microbial C
  !     GrowthEnvScalHeter=temperature+water stress function
  !     OMGR=rate constant for transferring nonstructural to structural C
  !     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
  !     FL=partitioning between labile and resistant microbial components
  !     OMC,OMN,OMP=nonstructural microbial C,N,P
  !
  MID3  = micpar%get_micb_id(iLbiom_reserve,NGL)
  CGOMZ = GrowthEnvScalHeter(NGL,K)*OMGR*AZMAX1(mBiomeHeter(ielmc,MID3,K))
  
  ! M=1:labile, 2, recalcitrant
  D745: DO M=1,2
    !reserve flux to kinetic and structural biomass
    NonstX2stBiomHeter(ielmc,M,NGL,K)=FL(M)*CGOMZ
    IF(mBiomeHeter(ielmc,MID3,K).GT.ZEROS)THEN
      Do NE=2,NumPlantChemElms
        NonstX2stBiomHeter(NE,M,NGL,K)=AZMAX1(mBiomeHeter(NE,MID3,K))*AMIN1(FL(M), &
          NonstX2stBiomHeter(ielmc,M,NGL,K)/mBiomeHeter(ielmc,MID3,K))
      ENDDO  
    ELSE
      NonstX2stBiomHeter(2:NumPlantChemElms,M,NGL,K)=0.0_r8
    ENDIF
    !
    !     MICROBIAL DECOMPOSITION FROM BIOMASS, SPECIFIC DECOMPOSITION
    !     RATE, TEMPERATURE
    !
    !     SPOMX=rate constant for microbial decomposition
    !     SPOMC=basal decomposition/mortality rate
    !     SPOMK=effect of low microbial C concentration on microbial decay
    !     RXOMC,RXOMN,RKillOMHeter=microbial C,N,P decomposition
    !     RkillLitfalOMHeter,RDOMN,RDOMP=microbial C,N,P LitrFall
    !     R3OMC,R3OMN,RkillRecycOMHeter=microbial C,N,P recycling
    !

    MID1=micpar%get_micb_id(iLbiom_kinetic,NGL)

    IF(mBiomeHeter(ielmc,MID3,K).GT.ZEROS .AND. mBiomeHeter(ielmc,MID1,K).GT.ZEROS)THEN
      CCC=AZMAX1(AMIN1(1.0_r8 &
        ,mBiomeHeter(ielmn,MID3,K)/(mBiomeHeter(ielmn,MID3,K)+mBiomeHeter(ielmc,MID3,K)*rNCOMC(iLbiom_reserve,NGL,K)) &
        ,mBiomeHeter(ielmp,MID3,K)/(mBiomeHeter(ielmp,MID3,K)+mBiomeHeter(ielmc,MID3,K)*rPCOMC(iLbiom_reserve,NGL,K))))
      CXC  = mBiomeHeter(ielmc,MID3,K)/mBiomeHeter(ielmc,MID1,K)
      C3C  = 1.0_r8/(1.0_r8+CXC/CKC)

      CNC  = AZMAX1(AMIN1(1.0_r8,mBiomeHeter(ielmc,MID3,K)/(mBiomeHeter(ielmc,MID3,K)+mBiomeHeter(ielmn,MID3,K)/rNCOMC(iLbiom_reserve,NGL,K))))
      CPC  = AZMAX1(AMIN1(1.0_r8,mBiomeHeter(ielmc,MID3,K)/(mBiomeHeter(ielmc,MID3,K)+mBiomeHeter(ielmp,MID3,K)/rPCOMC(iLbiom_reserve,NGL,K))))
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

    MID   = micpar%get_micb_id(M,NGL)
    SPOMX = SQRT(GrowthEnvScalHeter(NGL,K))*SPOMC(M)*SPOMK(M)
    DO NE=1,NumPlantChemElms
      RKillOMHeter(NE,M,NGL,K)=AZMAX1(AMIN1(mBiomeHeter(NE,MID,K),mBiomeHeter(NE,MID,K)*SPOMX))

      RkillRecycOMHeter(NE,M,NGL,K)= RKillOMHeter(NE,M,NGL,K)*RCCE(NE)

      RkillLitfalOMHeter(NE,M,NGL,K)=AZMAX1(RKillOMHeter(NE,M,NGL,K)-RkillRecycOMHeter(NE,M,NGL,K))
      !
      !     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
      !     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
      !
      !     RHOMC,RHOMN,RkillLitrfal2HumOMHeter=transfer of microbial C,N,P LitrFall to humus
      !     EHUM=humus transfer fraction from hour1.f
      !     RkillLitrfal2ResduOMHeter,RkillLitrfal2ResduOMHeter,RCOMP=transfer of microbial C,N,P LitrFall to residue
      !
      RkillLitrfal2HumOMHeter(NE,M,NGL,K)=RkillLitfalOMHeter(NE,M,NGL,K)*EHUM
      !
      !     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
      !
      RkillLitrfal2ResduOMHeter(NE,M,NGL,K)=RkillLitfalOMHeter(NE,M,NGL,K)-RkillLitrfal2HumOMHeter(NE,M,NGL,K)

    ENDDO
  ENDDO D745

  !=========================================================================================================
  !     MICROBIAL DECOMPOSITION/renomalization WHEN MAINTENANCE RESPIRATION
  !     EXCEEDS UPTAKE
  !
  !     OMC,OMN,OMP=microbial C,N,P
  !     RMaintRespHeter=total maintenance respiration
  !     RMaintDefcitcitHeter=senescence respiration
  !     RCCC=C recycling fraction
  !     RMaintDefcitKillOMHeter,RXMMN,RXMMP=microbial C,N,P loss from senescence
  !     RMaintDmndHeter=maintenance respiration
  !     CNOMA,CPOMA=N:C,P:C ratios of active biomass
  !     RDMMC,RMaintDefcitLitrfalOMHeter,RDMMP=microbial C,N,P LitrFall from senescence
  !     RMaintDefcitRecycOMHeter,R3MMN,R3MMP=microbial C,N,P recycling from senescence
  !

  IF(RMaintDefcitcitHeter.GT.ZEROS .AND. RMaintRespHeter(NGL,K).GT.ZEROS .AND. RCCC.GT.ZERO)THEN
    FracMaintDeficit=RMaintDefcitcitHeter/RMaintRespHeter(NGL,K)
    D730: DO M=1,2
      MID                                    = micpar%get_micb_id(M,NGL)
      !Ordinary mortality and starvation share this compartment's donor pools.
      AvailableBiomass=MAX(0._r8,mBiomeHeter(1:NumPlantChemElms,MID,K)-RKillOMHeter(1:NumPlantChemElms,M,NGL,K))
      RMaintDefcitKillOMHeter(ielmc,M,NGL,K) = AMIN1(AvailableBiomass(ielmc),AZMAX1(FracMaintDeficit*RMaintDmndHeter(M,NGL,K)/RCCC))
      RMaintDefcitKillOMHeter(ielmn,M,NGL,K) = AMIN1(AvailableBiomass(ielmn),AZMAX1(RMaintDefcitKillOMHeter(ielmc,M,NGL,K)*rCNBiomeActHeter(ielmn,NGL,K)))
      RMaintDefcitKillOMHeter(ielmp,M,NGL,K) = AMIN1(AvailableBiomass(ielmp),AZMAX1(RMaintDefcitKillOMHeter(ielmc,M,NGL,K)*rCNBiomeActHeter(ielmp,NGL,K)))
        
      DO NE=1,NumPlantChemElms
        RMaintDefcitRecycOMHeter(NE,M,NGL,K)   = RMaintDefcitKillOMHeter(NE,M,NGL,K)*RCCE(NE)
        RMaintDefcitLitrfalOMHeter(NE,M,NGL,K) = AZMAX1(RMaintDefcitKillOMHeter(NE,M,NGL,K)-RMaintDefcitRecycOMHeter(NE,M,NGL,K))
        !
        ! HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
        ! PRODUCTS
        !
        ! RMaintDefLitrfal2HumOMHeter,RHMMN,RMaintDefLitrfal2HumOMHeter=transfer of senesence LitrFall C,N,P to humus
        ! EHUM=humus transfer fraction
        ! RMaintDefLitrfal2ResduOMHeter,RCMMN,RMaintDefLitrfal2ResduOMHeter=transfer of senesence LitrFall C,N,P to residue
        !

        RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K)   = RMaintDefcitLitrfalOMHeter(NE,M,NGL,K)*EHUM
        RMaintDefLitrfal2ResduOMHeter(NE,M,NGL,K) = RMaintDefcitLitrfalOMHeter(NE,M,NGL,K)-RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K)

      ENDDO
    ENDDO D730
  ELSE
    D720: DO M=1,2
      DO NE=1,NumPlantChemElms
        RMaintDefcitKillOMHeter(NE,M,NGL,K)       = 0.0_r8
        RMaintDefcitLitrfalOMHeter(NE,M,NGL,K)    = 0.0_r8
        RMaintDefcitRecycOMHeter(NE,M,NGL,K)      = 0.0_r8
        RMaintDefLitrfal2HumOMHeter(NE,M,NGL,K)   = 0.0_r8
        RMaintDefLitrfal2ResduOMHeter(NE,M,NGL,K) = 0.0_r8
      ENDDO

    ENDDO D720
  ENDIF

  !Ordinary mortality contributes residue even without a maintenance deficit.
  !Accumulate each guild once, including any maintenance-deficit residue.
  DO M=1,2
    DO NE=1,NumPlantChemElms
      RCCMEHeter(NE,M,K)=RCCMEHeter(NE,M,K) &
        +RkillLitrfal2ResduOMHeter(NE,M,NGL,K) &
        +RMaintDefLitrfal2ResduOMHeter(NE,M,NGL,K)
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine GatherHetertrophAnabolicFlux

end module MicBGCMod
