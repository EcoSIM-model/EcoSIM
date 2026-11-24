module CyanoBacteriaMod
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
  use EcosimConst
  use EcoSIMSolverPar
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod,    only: AerobicHeterO2Uptake, CalcRespMaintHeter, &
                                  StageFuncGuild
  implicit none

  private
  public :: CyanoBacteriaCatabolism

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains
!------------------------------------------------------------------------------------------
  subroutine CyanoBacteriaCatabolism(I,J,N,K,RMOMK,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
  !

  implicit none
  integer, intent(in) :: I,J,N,K
  real(r8), intent(in) :: RMOMK(2)
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(OMCplx_State_type), intent(inout) :: ncplxs
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='CyanoBacteriaCatabolism'

  real(r8), parameter :: K_low_photo= 6.e-3_r8   ![h-1]
  integer :: NGL
  real(r8) :: EO2Q
  real(r8) :: fBIOM
  real(r8) :: fPAR,xCO2,r_photo,f_low_photo
  real(r8) :: RPhotoResp,PhotoO2Surplus,NetO2Release
  real(r8) :: RGOMP,RHeterRespUlm
  
  real(r8) :: FSBSTC,FSBSTA
  real(r8) :: RGOCY,RGOCZ,RGOAZ
  real(r8) :: RGOCX,RGOAX
  real(r8) :: OXKX,FOXYX
  real(r8) :: WatStressMicb
  real(r8) :: f_cyno_heter
  associate(                                                               &
    OMActHeter             => nmics%OMActHeter,                            & !Active microbial C biomass by heterotrophic guild and complex K
    FBiomStoiScalarHeter   => nmics%FBiomStoiScalarHeter,                  & !Combined N/P stoichiometric multiplier on guild metabolic capacity [-]
    GrowthEnvScalHeter     => nmics%GrowthEnvScalHeter,                    & !Temperature and water-potential multiplier on heterotrophic growth [-]
    WSensGroHeter          => nmics%WSensGroHeter,                         & !Guild soil-water-potential multiplier on heterotrophic growth [-]; not referenced here
    TSensGroHeter          => nmics%TSensGroHeter,                         & !Guild temperature multiplier on heterotrophic growth [-]; not referenced here
    TempMaintRHeter        => nmics%TempMaintRHeter,                       & !Guild temperature multiplier on heterotrophic maintenance [-]; not referenced here
    FracOMActHeter         => nmics%FracOMActHeter,                        & !Guild/complex fraction of total active microbial C in the layer [-]
    FracHeterBiomOfActK    => nmics%FracHeterBiomOfActK,                   & !Guild fraction of active heterotrophic biomass in complex K [-]; not referenced here
    FSBSTHeter             => nmicdiag%FSBSTHeter,                         & !Guild substrate-response factor; larger values mean less limitation [-]
    TotActMicrobiom        => nmicdiag%TotActMicrobiom,                    & !Layer total active microbial C across heterotrophs and autotrophs
    TOMEK                  => nmicdiag%TOMEK,                              & !Total active heterotrophic C/N/P in each substrate complex K; not referenced here
    TSensGrowth            => nmicdiag%TSensGrowth,                        & !Layer temperature response for microbial growth [-]; not referenced here
    TSensMaintR            => nmicdiag%TSensMaintR,                        & !Layer temperature response for microbial maintenance [-]; not referenced here
    RO2Dmnd4RespHeter      => nmicf%RO2Dmnd4RespHeter,                     & !Gross potential respiratory O2 demand (before the photosynthetic supply credit)
    OxyLimterHeter         => nmics%OxyLimterHeter,                        & !Actual/potential O2 uptake ratio; 1 means no O2 restriction [-]
    RO2DmndHeter           => nmicf%RO2DmndHeter,                          & !Total guild O2 demand before O2 limitation
    ROQC4HeterMicrobAct    => nmicf%ROQC4HeterMicrobAct,                   & !Guild activity proxy for substrate hydrolysis, with DOC concentration unconstrained
    ECHZHeter              => nmicf%ECHZHeter,                             & !Guild respiration fraction used to convert growth respiration to C uptake [-]
    RCO2ProdHeter          => nmicf%RCO2ProdHeter,                         & !CO2-C production by heterotrophic guild and complex
    RCH4ProdHeter          => nmicf%RCH4ProdHeter,                         & !CH4-C production by heterotrophic guild and complex
    RO2Uptk4RespHeter      => nmicf%RO2Uptk4RespHeter,                     & !Gross respiratory O2 consumption, including internally supplied O2
    RO2UptkHeter           => nmicf%RO2UptkHeter,                          & !Net external O2 uptake; negative means photosynthetic O2 release
    REcoUptkSoilO2M        => micflx%REcoUptkSoilO2M,                     & !Net external O2 exchange accumulated by transport substep
    RCO2FixCyano           => nmicf%RCO2FixCyano   ,                       & !Photosynthetic CO2-C fixation by cyanobacterial guild and complex
    RH2ProdHeter           => nmicf%RH2ProdHeter,                          & !H2 production by heterotrophic guild and complex
    RGOCP                  => nmicf%RGOCP,                                 & !DOC-supported potential respiration before O2 limitation
    RGOAP                  => nmicf%RGOAP,                                 & !Acetate-supported potential respiration before O2 limitation
    FOQC                   => nmicf%FOQC,                                  & !Guild share of DOC demand used to allocate the donor pool [-]
    FOQA                   => nmicf%FOQA,                                  & !Guild share of acetate demand used to allocate the donor pool [-]
    FGOCP                  => nmicf%FGOCP,                                 & !DOC fraction of organic-substrate respiration; excludes photosynthate [-]
    FGOAP                  => nmicf%FGOAP,                                 & !Acetate fraction of organic-substrate respiration; excludes photosynthate [-]
    fPhotoR                => nmicf%fPhotoR,                               & !Photosynthate-supported fraction of cyanobacterial gross respiration [-]
    RespGrossHeter         => nmicf%RespGrossHeter,                        & !Gross respiration C equivalent from the primary heterotrophic pathway
    RAcetateProdHeter      => nmicf%RAcetateProdHeter,                     & !Acetate-C production by heterotrophic guild and complex
    RO2DmndHetert          => micflx%RO2DmndHetert,                        & !Guild O2 demand retained for substrate-competition accounting
    RDOCUptkHeter          => micflx%RDOCUptkHeter,                        & !Potential DOC uptake used in guild substrate-competition accounting
    RAcetateUptkHeter      => micflx%RAcetateUptkHeter,                    & !Potential acetate uptake used in guild substrate-competition accounting
    RMaintRespHeter        => nmicf%RMaintRespHeter,                       & !Total hourly maintenance-C demand; photosynthate supplies it first
    tRGOXP                 => micflx%tRGOXP,                               & !Accumulated substrate-pool-limited potential respiration capacity
    tRGOZP                 => micflx%tRGOZP,                               & !Accumulated kinetic potential respiration capacity before donor/O2 limits
    RDOMEcoDmndPrev        => micfor%RDOMEcoDmndPrev,                      & !Previous-hour ecosystem DOC demand in each complex; competition denominator; not referenced here
    RAcetateEcoDmndPrev    => micfor%RAcetateEcoDmndPrev,                  & !Previous-hour ecosystem acetate demand in each complex; competition denominator; not referenced here
    RDOCUptkHeterPrev      => micfor%RDOCUptkHeterPrev,                    & !Previous-hour guild DOC uptake/demand used for competition; not referenced here
    RAcetateUptkHeterPrev  => micfor%RAcetateUptkHeterPrev,                & !Previous-hour guild acetate uptake/demand used for competition; not referenced here
    RO2EcoDmndPrev         => micfor%RO2EcoDmndPrev,                       & !Previous-hour ecosystem O2 demand; competition denominator
    RO2DmndHetertPrev      => micflx%RO2DmndHetertPrev,                    & !Previous-hour guild O2 demand used for competition
    PSISoilMatricP         => micfor%PSISoilMatricP,                       & !Soil matric water potential controlling microbial water stress; not referenced here
    PAR_rad                => micfor%PAR_rad,                              & !Photosynthetically active radiation for cyanobacterial photosynthesis
    ZEROS                  => micfor%ZEROS,                                & !Small mass or flux threshold used by the routine
    CCO2S                  => micstt%CCO2S,                                & !Dissolved CO2-C concentration for substrate saturation
    DOM                    => micstt%DOM,                                  & !Dissolved organic pools by species (DOC, DON, DOP, acetate) and complex K
    FOCA                   => ncplxs%FOCA,                                 & !DOC fraction of DOC plus acetate in each substrate complex [-]
    FOAA                   => ncplxs%FOAA,                                 & !Acetate-C fraction of DOC plus acetate in each substrate complex [-]
    CDOM                   => ncplxs%CDOM,                                 & !Dissolved substrate concentrations by DOM species and complex K
    tRespGrossHeterUlm     => naqfdiag%tRespGrossHeterUlm                  & !Total gross heterotrophic respiration before O2 limitation
  )
  call PrintInfo('beg '//subname)

  !reduction of heterotrophic rate compared to actual heterotrophs
  f_cyno_heter=0.1_r8
  if(PAR_rad.LE.ZEROS)then
    f_cyno_heter=0.01_r8
  endif

  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)

    IF(OMActHeter(NGL,K).LE.0.0_r8)cycle

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

    EO2Q=ENFX
    fBIOM=AZMAX1(FBiomStoiScalarHeter(NGL,K)*OMActHeter(NGL,K))

    !photosynthesis
    !6CO2 + 6H2O + light -> (CH2O)6 + 6O2

    if(PAR_rad.GT.ZEROS)then

      fPAR=AZMAX1(fPAR_func(PAR_rad))
      xCO2 = CCO2S/(CCO2S+CCKM)
      r_photo=fPAR*xCO2*GrowthEnvScalHeter(NGL,K)
      RCO2FixCyano(NGL,K) = r_photo * fBIOM

      !Fixation is the available photosynthetic C, not respiration itself.
      !Use it for maintenance first; the residual is available for growth
      !and N2 fixation below, without an extra growth-respiration charge.
      RPhotoResp=MIN(RMaintRespHeter(NGL,K),RCO2FixCyano(NGL,K))

      !(CH2O)6  + 6O2 -> 6CO2 +6 H2O
      f_low_photo = K_low_photo / (r_photo + K_low_photo)
    else
      f_low_photo         = 1._r8
      RCO2FixCyano(NGL,K) = 0._r8
      RPhotoResp          = 0._r8
    endif

    !grow on DOC and acetate
    FSBSTC            = CDOM(idom_doc,K)/(CDOM(idom_doc,K)+OQKM)
    FSBSTA            = CDOM(idom_acetate,K)/(CDOM(idom_acetate,K)+OQKA)
    FSBSTHeter(NGL,K) = FOCA(K)*FSBSTC+FOAA(K)*FSBSTA

    RGOCY  = fBIOM*VMXO*GrowthEnvScalHeter(NGL,K)*f_low_photo*f_cyno_heter
    RGOCZ  = RGOCY*FSBSTC*FOCA(K)  !MM uptake of DOC
    RGOAZ  = RGOCY*FSBSTA*FOAA(K)  !MM uptake of acetate

    !obtain kinetically unlimited DOM/acetate uptake
    RGOCX = AZMAX1(DOM(idom_doc,K)*FOQC(NGL,K)*EO2Q)      !DOC respiration
    RGOAX = AZMAX1(DOM(idom_acetate,K)*FOQA(NGL,K)*EO2A)  !acetate respiraiton

    !obtain the final uptake
    RGOCP(NGL,K) = AMIN1(RGOCX,RGOCZ)             !DOC respiraiton
    RGOAP(NGL,K) = AMIN1(RGOAX,RGOAZ)             !acetate respiraiton
    RHeterRespUlm = RGOCP(NGL,K)+RGOAP(NGL,K)
    RGOMP        = RPhotoResp+RHeterRespUlm      !total C respiration before O2 limitation

    tRGOXP = tRGOXP+RGOCX+RGOAX
    tRGOZP = tRGOZP+RGOCZ+RGOAZ            !potential C oxidation without C and O2 limitation
    IF(RHeterRespUlm.GT.0._r8)THEN
      FGOCP(NGL,K) = RGOCP(NGL,K)/RHeterRespUlm
      FGOAP(NGL,K) = RGOAP(NGL,K)/RHeterRespUlm
    ELSE
      FGOCP(NGL,K) = 1.0_r8
      FGOAP(NGL,K) = 0.0_r8      
    ENDIF

    !Photosynthesis supplies its own respiration first, then organic-C
    !respiration. Only the remaining positive demand competes for soil O2.
    PhotoO2Surplus=2.667_r8*AZMAX1(RCO2FixCyano(NGL,K)-RPhotoResp)
    RO2Dmnd4RespHeter(NGL,K) = 2.667_r8*RGOMP
    RO2DmndHeter(NGL,K)      = AZMAX1(2.667_r8*RHeterRespUlm-PhotoO2Surplus)
    !Respiration-weighted conversion to organic substrate uptake.
    ECHZHeter(NGL,K)         = 1._r8/(FGOCP(NGL,K)/EO2Q+FGOAP(NGL,K)/EO2A)

    !make a copy for flux limiter
    RO2DmndHetert(NGL,K)       = RO2DmndHeter(NGL,K)
    RDOCUptkHeter(NGL,K)       = RGOCZ               !potential DOC (O2-unlimited) uptake flux
    RAcetateUptkHeter(NGL,K)   = RGOAZ               !potential acetate (O2-unlimited) uptake flux

    call AerobicHeterO2Uptake(I,J,NGL,N,K,FOXYX,OXKX,micfor,micstt,nmicf,nmics,micflx)

    !The uptake solver limits only external demand. Include the internal
    !photosynthetic supply when limiting organic respiration, even in anoxia.
    IF(RHeterRespUlm.GT.0._r8)THEN
      OxyLimterHeter(NGL,K)=MIN(1._r8, &
        (RO2UptkHeter(NGL,K)+PhotoO2Surplus)/(2.667_r8*RHeterRespUlm))
    ENDIF
    !Export unused photosynthetic O2 through both hourly and substep budgets.
    !Keep RO2DmndHetert nonnegative for next hour's competition weights.
    NetO2Release=AZMAX1(PhotoO2Surplus-2.667_r8*RHeterRespUlm)
    RO2UptkHeter(NGL,K)=RO2UptkHeter(NGL,K)-NetO2Release
    REcoUptkSoilO2M(1:NPH)=REcoUptkSoilO2M(1:NPH)-NetO2Release/REAL(NPH,r8)

    ROQC4HeterMicrobAct(NGL,K) = RGOCY*OxyLimterHeter(NGL,K) !C demand for oxidation
    RespGrossHeter(NGL,K)      = RPhotoResp+RHeterRespUlm*OxyLimterHeter(NGL,K)
    
    IF(RespGrossHeter(NGL,K).GT.0._R8)THEN
      fPhotoR(NGL,K)   = RPhotoResp/RespGrossHeter(NGL,K) 
    ELSE
      fPhotoR(NGL,K)   = 0._R8
    ENDIF  

    RCO2ProdHeter(NGL,K)       = RespGrossHeter(NGL,K)-RCO2FixCyano(NGL,K)
    RAcetateProdHeter(NGL,K)   = 0.0_r8
    RCH4ProdHeter(NGL,K)       = 0.0_r8
    RO2Uptk4RespHeter(NGL,K)   = 2.667_r8*RespGrossHeter(NGL,K)
    RH2ProdHeter(NGL,K)        = 0.0_r8
    !The diagnostic is the maximum respiratory C supply, including fixed C
    !that may later support N2 fixation or remain as residual growth.
    tRespGrossHeterUlm         = tRespGrossHeterUlm+RCO2FixCyano(NGL,K)+RHeterRespUlm

  ENDDO

  call PrintInfo('end '//subname)
  end associate

  end subroutine CyanoBacteriaCatabolism
!------------------------------------------------------------------------------------------
  function fPAR_func(PAR_RAD)result(ans)
  implicit none
  real(r8), intent(in) :: PAR_RAD                ![umol photon m-2 s-1]
  real(r8), parameter :: alpha_cyno = 1.1e-4_r8  ! [h-1] / [umol photon m-2 s-1]
  real(r8), parameter :: beta_cyno  = 3.3e-4_r8  ! inverse [PAR]
  real(r8), parameter :: gamma_cyno = 1.5e-3_r8  ! inverse [PAR]

  real(r8) :: ans   !photosynthesis rate [h-1]


  ans = alpha_cyno * PAR_rad * (1._r8 - beta_cyno*PAR_rad) &
        / (1._r8 + gamma_cyno*PAR_rad)

  end function fPAR_func

  end module CyanoBacteriaMod
