module MethanogenMod
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: PrintInfo
  use minimathmod,          only: AZMAX1, real_truncate
  use TracerIDMod
  use EcosimConst,          only: RGASC
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod,    only: CalcRespMaintAutor, CalcRespMaintHeter, &
                                  StageAutotroph, StageFuncGuild

  implicit none

  private
  public :: AcetoMethanogenCatabolism, H2MethanogensCatabolism

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains

!------------------------------------------------------------------------------------------

  subroutine AcetoMethanogenCatabolism(N,K,RMOMK,micfor,micstt,naqfdiag,nmicf,nmics,ncplxs,micflx,nmicdiag)
  implicit none
  integer, intent(in) :: N,K
  real(r8), intent(in) :: RMOMK(2)

  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout):: nmics
  type(OMCplx_State_type),intent(in):: ncplxs
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  integer :: NGL
  reaL(r8) :: RGOMP         !substrate-limited potential respiration
  real(r8) :: GOMX,GOMM
  real(r8) :: RGOGY,RGOGZ
  real(r8) :: RGroMax   !kinetically unlimited acetate uptake
  real(r8) :: FNH4X
  real(r8) :: FNB3X,FNB4X,FNO3X,FPO4X,FPOBX,FP14X,FP1BX,FOQA
  real(r8)  :: WatStressMicb

! begin_execution
  associate(                                            &
    FBiomStoiScalarHeter => nmics%FBiomStoiScalarHeter, & !Combined N/P stoichiometric multiplier on guild metabolic capacity [-]
    OMActHeter           => nmics%OMActHeter,           & !Active microbial C biomass by heterotrophic guild and complex K
    GrowthEnvScalHeter   => nmics%GrowthEnvScalHeter,   & !Temperature and water-potential multiplier on heterotrophic growth [-]
    FSBSTHeter           => nmicdiag%FSBSTHeter,        & !Guild substrate-response factor; larger values mean less limitation [-]
    RO2DmndHeter         => nmicf%RO2DmndHeter,         & !Total guild O2 demand before O2 limitation; zero for this anaerobic pathway
    RO2Dmnd4RespHeter    => nmicf%RO2Dmnd4RespHeter,    & !Potential O2 demand supporting heterotrophic gross respiration; zero for this anaerobic pathway
    ROQC4HeterMicrobAct  => nmicf%ROQC4HeterMicrobAct,  & !Guild activity proxy for substrate hydrolysis, with DOC concentration unconstrained
    RAcetateProdHeter    => nmicf%RAcetateProdHeter,    & !Acetate-C production by heterotrophic guild and complex
    RO2Uptk4RespHeter    => nmicf%RO2Uptk4RespHeter,    & !Realized O2 uptake attributed to heterotrophic gross respiration; zero for this anaerobic pathway
    RCO2ProdHeter        => nmicf%RCO2ProdHeter,        & !CO2-C production by heterotrophic guild and complex
    RH2ProdHeter         => nmicf%RH2ProdHeter,         & !H2 production by heterotrophic guild and complex
    RCH4ProdHeter        => nmicf%RCH4ProdHeter,        & !CH4-C production by heterotrophic guild and complex
    FOQC                 => nmicf%FOQC,                 & !Guild share of DOC demand used to allocate the donor pool [-]
    ECHZHeter            => nmicf%ECHZHeter,            & !Guild respiration fraction used to convert growth respiration to C uptake [-]
    FGOCP                => nmicf%FGOCP,                & !DOC-supported fraction of total primary guild respiration [-]
    RespGrossHeter       => nmicf%RespGrossHeter,       & !Gross respiration C equivalent from the primary heterotrophic pathway
    FGOAP                => nmicf%FGOAP,                & !Acetate-supported fraction of total primary guild respiration [-]
    TotActMicrobiom      => nmicdiag%TotActMicrobiom,   & !Layer total active microbial C across heterotrophs and autotrophs
    CDOM                 => ncplxs%CDOM,                & !Dissolved substrate concentrations by DOM species and complex K
    DOM                  => micstt%DOM,                 & !Dissolved organic pools by species (DOC, DON, DOP, acetate) and complex K
    RO2DmndHetert        => micflx%RO2DmndHetert,       & !Guild O2 demand retained for substrate-competition accounting; zero for this anaerobic pathway
    RDOCUptkHeter        => micflx%RDOCUptkHeter,       & !Potential DOC uptake used in guild substrate-competition accounting
    RAcetateUptkHeter    => micflx%RAcetateUptkHeter,   & !Potential acetate uptake used in guild substrate-competition accounting
    PSISoilMatricP       => micfor%PSISoilMatricP,      & !Soil matric water potential controlling microbial water stress
    ZERO                 => micfor%ZERO,                & !Small dimensionless or concentration threshold used by the routine
    TKS                  => micfor%TKS                  & !Layer absolute temperature [K]
  )
  !
  !     GOMX=acetate effect on energy yield
  !     ECHZHeter=growth respiration efficiency of aceto. methanogenesis
  !
  !loop over all guilds of a given functional group
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    IF(OMActHeter(NGL,K).LE.0.0_r8)cycle
    WatStressMicb=real_truncate(EXP(0.2_r8*AMAX1(PSISoilMatricP,-500._r8)),1.e-3_r8)
    !prepare parameters
    call StageFuncGuild(N,NGL,K,TotActMicrobiom,FOQC(NGL,K),FOQA,micfor,naqfdiag,nmicdiag,nmics)

    call CalcRespMaintHeter(NGL,K,RMOMK,micfor,micstt,micflx,nmicf,nmics)
    !1.E-3 is to convert J into kJ
    GOMX = RGASC*1.E-3_r8*TKS*LOG((AMAX1(ZERO,CDOM(idom_acetate,K))/OAKI))
    GOMM = GOMX/24.0_r8
    ECHZHeter(NGL,K) = AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GC4X+GOMM))/EOMH)))
    !
    !     RESPIRATION RATES BY ACETOTROPHIC METHANOGENS 'RGOMP' FROM
    !     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
    !     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL C
    !     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
    !     MICROBIAL COMPETITION FACTOR
    !
    !     COQA=DOA concentration
    !     OQKAM=Km for acetate uptake,FBiomStoiScalarHeter=N,P limitation
    !     VMXCH4gAcet=specific respiration rate
    !     WatStressMicb=water stress effect, OMA=active biomass
    !     TSensGrowth=temp stress effect, FOQA= acetate limitation
    !     RGroMax=substrate-limited respiration of acetate
    !     RGroMax=competition-limited respiration of acetate
    !     OQA=acetate, FOQA=fraction of biological demand for acetate
    !     RGOMP=O2-unlimited respiration of acetate
    !     ROXY*=O2 demand, RDOCUptkHeter,ROQCA=DOC, acetate demand
    !     ROQC4HeterMicrobAct=microbial respiration used to represent microbial activity
    !
    FSBSTHeter(NGL,K)          = CDOM(idom_acetate,K)/(CDOM(idom_acetate,K)+OQKAM)
    RGOGY                      = FBiomStoiScalarHeter(NGL,K)*VMXCH4gAcet*OMActHeter(NGL,K)*GrowthEnvScalHeter(NGL,K)
    RGOGZ                      = RGOGY*FSBSTHeter(NGL,K)
    RGroMax                    = AZMAX1(DOM(idom_acetate,K)*FOQA*ECHZHeter(NGL,K))
    RGOMP                      = AMIN1(RGroMax,RGOGZ)
    FGOCP(NGL,K)               = 0.0_r8
    FGOAP(NGL,K)               = 1.0_r8
    RO2Dmnd4RespHeter(NGL,K)   = 0.0_r8
    RO2DmndHeter(NGL,K)        = 0.0_r8
    RO2DmndHetert(NGL,K)       = 0.0_r8
    RDOCUptkHeter(NGL,K)       = 0.0_r8
    RAcetateUptkHeter(NGL,K)   = RGOGZ
    ROQC4HeterMicrobAct(NGL,K) = 0.0_r8

    !given CH3COOH -> CH4+CO2, 0.5 is into CH4.
    naqfdiag%tCH4ProdAceto=naqfdiag%tCH4ProdAceto+0.5_r8*RGOMP

    ! CH3COOH -> CO2 + CH4
    RespGrossHeter(NGL,K)   = RGOMP
    RCO2ProdHeter(NGL,K)    = 0.50_r8*RespGrossHeter(NGL,K)
    RAcetateProdHeter(NGL,K)  = 0.0_r8
    RCH4ProdHeter(NGL,K)    = AZMAX1(RespGrossHeter(NGL,K)-RCO2ProdHeter(NGL,K))
    RO2Uptk4RespHeter(NGL,K)= RO2Dmnd4RespHeter(NGL,K)
    RH2ProdHeter(NGL,K)     = 0.0_r8
  ENDDO
  end associate
  end subroutine AcetoMethanogenCatabolism

!------------------------------------------------------------------------------------------

  subroutine H2MethanogensCatabolism(I,J,N,RMOMK,TOMEAutoKC,micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)
  !
  !Hydrogenotrophic CH4 production
  !CO2 + 4H2 -> CH4 + 2H2O
  !H2 is produced from fermentation
  !CO2+0.667H2-> CH4+ 3H2O
  !use CO2 for both energy generation and C biomass
  implicit none
  integer, intent(in) :: I,J, N
  real(r8), intent(in) :: RMOMK(2)
  real(r8), intent(in) :: TOMEAutoKC
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type), intent(inout) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='H2MethanogensCatabolism'
  real(r8) :: GH2X,GH2H
  real(r8) :: H2GSX       !Remaining shared H2 supply, including fermentation
  real(r8) :: RespPerCatabolicC,H2PerCatabolicC,H2PerGrowthRespC
  real(r8) :: RVOXMaxH2,GrowthC,H2Demand
  real(r8) :: VMAX
  REAL(R8) :: XCO2
  real(r8) :: RGOMP,RVOXP,GH2C
  real(r8) :: ECH2         !Efficiency of converting CO2 into biomass (CH2O) by hydrogenotrophic methanogen
  real(r8), parameter :: GCHA=38.9/12._r8 !Gibbs free energy for anabolic reaction, CO2(aq)+2H2(aq)->CH2O + H2O,[kJ (gC)-1]

  integer  :: NGL

  associate(                                                &
    GrowthEnvScalAutor     => nmics%GrowthEnvScalAutor,     &  !Temperature and water-potential multiplier on autotrophic growth [-]
    FBiomNutStoiScalAutor  => nmics%FBiomNutStoiScalAutor,  &  !Combined N/P stoichiometric multiplier on guild metabolic capacity [-]
    FSBSTAutor             => nmicdiag%FSBSTAutor,          &  !Dissolved H2 concentration saturation factor for methanogenesis [-]
    OMActAutor             => nmics%OMActAutor,             &  !Active microbial C biomass by autotrophic guild
    RO2Dmnd4GrossRespAutor => nmicf%RO2Dmnd4GrossRespAutor, &  !Potential O2 demand supporting autotrophic gross respiration; zero for this anaerobic pathway
    RO2Uptk4RespAutor      => nmicf%RO2Uptk4RespAutor,      &  !Realized O2 uptake attributed to autotrophic gross respiration; zero for this anaerobic pathway
    RespGrossAutor         => nmicf%RespGrossAutor,         &  !Gross respiration C equivalent by autotrophic guild
    RO2UptkAutor           => nmicf%RO2UptkAutor,           &  !Realized total O2 uptake by autotrophic guild; zero for this anaerobic pathway
    RCO2ProdAutor          => nmicf%RCO2ProdAutor,          &  !CO2-C production by autotrophic guild; not referenced here
    RCH4ProdAutor          => nmicf%RCH4ProdAutor,          &  !CH4-C production by autotrophic guild
    ECHZAutor              => nmicf%ECHZAutor,              &  !Guild respiration fraction used to convert growth respiration to C uptake [-]
    JGniA                  => micpar%JGniA,                 &  !First guild index for each autotrophic functional group; not referenced here
    JGnfA                  => micpar%JGnfA,                 &  !Last guild index for each autotrophic functional group; not referenced here
    TKS                    => micfor%TKS,                   &  !Layer absolute temperature [K]
    CH2GS                  => micstt%CH2GS,                 &  !Dissolved H2 concentration used in energy-yield and saturation calculations
    CCO2S                  => micstt%CCO2S,                 &  !Dissolved CO2-C concentration for substrate saturation
    H2GS                   => micstt%H2GS,                  &  !Dissolved H2 donor pool
    RH2UptkAutor           => nmicdiag%RH2UptkAutor,        &  !Accumulated H2 uptake for methane production; biomass H2 is added later
    RMaintRespAutor        => micflx%RMaintRespAutor,        & !Total hourly autotrophic guild maintenance-C demand
    RO2MetaDmndAutor       => micflx%RO2MetaDmndAutor       &  !Total autotrophic O2 demand from respiration and substrate oxidation; zero for this anaerobic pathway
  )
  !     begin_execution
  !
  !     CO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
  !     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND H2
  !
  !     GH2H=energy yield of hydrogenotrophic methanogenesis per g C
  !     ECHZ=growth respiration efficiency of hydrogen. methanogenesis
  !     VMAX=substrate-unlimited H2 oxidation rate
  !     H2GSX=aqueous H2 (H2GS) + total H2 from fermentation (tCResp4H2Prod)
  !     CH2GS=H2 concentration, H2KM=Km for H2 uptake
  !     RGOMP=H2 oxidation, ROXY*=O2 demand
  !
  !     and energy yield of hydrogenotrophic
  !     methanogenesis GH2X at ambient H2 concentration CH2GS
  !     CCO2S=aqueous CO2 concentration
  !
  call PrintInfo('beg '//subname)
  XCO2         = CCO2S/(CCO2S+CCKM)
  RH2UptkAutor = 0.0_r8
  !0.111 gH2 per gC fermented: C6H12O6 + 2H2O -> 2C2H4O2 + 4H2 + 2CO2.
  H2GSX = AZMAX1(H2GS+0.111_r8*naqfdiag%tCResp4H2Prod)
  !TODO: add an explicit H2 competition factor to allocate this supply among guilds.
  !Until then, guilds draw sequentially from the shared remaining H2GSX budget.
  DO NGL=micpar%JGniA(N),micpar%JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle
    call StageAutotroph(NGL,N,TOMEAutoKC,micfor,nmics,nmicdiag)

    call CalcRespMaintAutor(I,J,NGL,RMOMK,micfor,micstt,micflx,nmicf,nmics)

    !Use catabolic reaction: CO2(aq)+4H2(aq) -> CH4 + 2H2O, 8/12=0.667, 1.5=12/8,
    !to drive anabolic reaction: CO2(aq)+2H2(aq) -> CH2O + H2O,

    GH2X = RGASC*1.E-3_r8*TKS*LOG((AMAX1(1.0E-05_r8,CH2GS)/H2KI)**2)
    GH2C = RGASC*1.E-3_r8*TKS*LOG((AMAX1(1.0E-05_r8,CH2GS)/H2KI)**4)/12._r8
    GH2H = GH2X/12.0_r8
    ECH2 = AMIN1(1._r8,AZMAX1((GCOX+GH2C)/(GCHA+GH2H)))
    !biomass yield as measured based on C, using respiration CH2O +2H2 -> CH4 + H2O for energy
    ECHZAutor(NGL)  = AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+AZMAX1((GCOX+GH2H))/EOMH)))
    VMAX            = OMActAutor(NGL)*VMXCH4gH2*GrowthEnvScalAutor(NGL)*FBiomNutStoiScalAutor(NGL)*XCO2
    FSBSTAutor(NGL) = CH2GS/(CH2GS+H2KM)

    !For catabolic rate P=RVOXP, R=RespPerCatabolicC*P is respiration.
    !Methane C is P+R; growth C is MAX(0,R-maintenance)*(1/ECHZ-1).
    !Use the same H2:C factors as the methane and later biomass flux updates.
    RespPerCatabolicC = ECH2*ECHZAutor(NGL)
    H2PerCatabolicC   = 0.667_r8*(1._r8+RespPerCatabolicC)
    H2PerGrowthRespC  = 0.333_r8*(1._r8/ECHZAutor(NGL)-1._r8)

    !Solve the H2 budget on the maintenance-only or growth branch.
    RVOXMaxH2 = H2GSX/H2PerCatabolicC
    IF(RespPerCatabolicC*RVOXMaxH2.GT.RMaintRespAutor(NGL))THEN
      RVOXMaxH2 = (H2GSX+H2PerGrowthRespC*RMaintRespAutor(NGL)) &
        /(H2PerCatabolicC+H2PerGrowthRespC*RespPerCatabolicC)
    ENDIF
    RVOXP = AZMAX1(AMIN1(RVOXMaxH2,VMAX*FSBSTAutor(NGL)))
    RGOMP = RVOXP*RespPerCatabolicC
    GrowthC = AZMAX1(RGOMP-RMaintRespAutor(NGL))*(1._r8/ECHZAutor(NGL)-1._r8)
    H2Demand = 0.667_r8*(RVOXP+RGOMP)+0.333_r8*GrowthC
    !Reserve both methane and biomass demand before considering the next guild.
    !Biomass H2 is added to RH2UptkAutor later in AutotrophAnabolicUpdate.
    H2GSX = AZMAX1(H2GSX-H2Demand)

    RO2Dmnd4GrossRespAutor(NGL) = 0.0_r8
    RO2MetaDmndAutor(NGL)       = 0.0_r8

    !obtains CO2 uptake for energy generation
    RespGrossAutor(NGL)    = RGOMP
    RCH4ProdAutor(NGL)     = RVOXP+RGOMP
    naqfdiag%tCH4ProdH2    = naqfdiag%tCH4ProdH2+RCH4ProdAutor(NGL)
    RO2Uptk4RespAutor(NGL) = 0._r8
    RH2UptkAutor           = RH2UptkAutor+0.667_r8*RCH4ProdAutor(NGL)
    RO2UptkAutor(NGL)      = 0.0_r8
  ENDDO
!
  call PrintInfo('end '//subname)
  end associate
  end subroutine H2MethanogensCatabolism
end module MethanogenMod
