module NoduleBGCMod

  use data_kind_mod, only: r8 => DAT_KIND_R8
  use minimathmod,   only: safe_adb, AZMAX1
  use DebugToolMod,  only: PrintInfo
  use abortutils  ,  only: endrun
  use EcosimConst
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantSoilChemistryAPIData, only : plt_soilchem
  use PlantAllometryAPIData, only : plt_allom
  use PlantBiomassAPIData, only : plt_biom
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantRootBGCAPIData, only : plt_rbgc
  use PlantMathFuncMod
  use PlantBGCPars
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: CanopyNoduleBiochemistry
  public :: RootNodulBiochemistry
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine CanopyNoduleBiochemistry(I,J,NZ,NB,TFN5,WFNG)
  !
  !Nodule update for branch NB on pft NZ
  implicit none
  integer, intent(in) :: I,J,NZ,NB
  real(r8), intent(in) :: TFN5,WFNG
  character(len=*), parameter :: subname='CanopyNoduleBiochemistry'
  integer :: M
  real(r8) :: ZPOOLD
  real(r8) :: CCC,CNC,CPC
  REAL(R8) :: cpoolt
  real(r8) :: CPOOLD
  real(r8) :: NodulNonstElmConc(NumPlantChemElms)
  real(r8) :: NodulUseNonstElms_brch(NumPlantChemElms)
  real(r8) :: CCNDLB
  real(r8) :: FCNPF
  real(r8) :: FXRNX
  real(r8) :: NoduleBiomCGrowth
  real(r8) :: PPOOLD
  real(r8) :: RCO2T
  real(r8) :: RespNonst_Oltd,Rmaint,RXNDL
  real(r8) :: RGNDL,RSNDL
  real(r8) :: RGN2P,Rauto4Nfix_brch
  real(r8) :: RCanopyN2Fix
  real(r8) :: NoduleElmDecayLoss(NumPlantChemElms)
  real(r8) :: NoduleElmntDecay2Litr(NumPlantChemElms)
  real(r8) :: NodulElmDecayRecyc(NumPlantChemElms)
  real(r8) :: NoduleCResp
  real(r8) :: NodulELmLoss2Senes(NumPlantChemElms)
  real(r8) :: NodulELmSenes2Litr(NumPlantChemElms)
  real(r8) :: NodulELmSenes2Recyc(NumPlantChemElms)
  real(r8) :: SPNDLI
  real(r8) :: SPNDX
  real(r8) :: WTLSB1,WTNDB1,WTLSBT
  real(r8) :: XFRE(NumPlantChemElms),NoduleInfectE(NumPlantChemElms)
  REAL(R8) :: RCCC,RCCN,RCCP
  integer :: NE
!     begin_execution
  associate(                                                          &
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    ZERO                      => plt_site%ZERO                       ,& !input  :threshold zero for numerical stability, [-]
    AREA3                     => plt_site%AREA3                      ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    k_fine_comp               => pltpar%k_fine_comp                  ,& !input  :fine litter complex id
    PlantElmAllocMat4Litr     => plt_soilchem%PlantElmAllocMat4Litr  ,& !input  :litter kinetic fraction, [-]
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft        ,& !input  :N2 fixation type,[-]
    fTCanopyGroth_pft         => plt_pheno%fTCanopyGroth_pft         ,& !input  :canopy temperature growth function, [-]
    ifoliar                   => pltpar%ifoliar                      ,& !input  :group id of plant foliar litter
    NoduGrowthYield_pft       => plt_allom%NoduGrowthYield_pft       ,& !input  :nodule growth yield, [g g-1]
    rNCNodule_pft             => plt_allom%rNCNodule_pft             ,& !input  :nodule N:C ratio, [gN gC-1]
    rPCNoduler_pft            => plt_allom%rPCNoduler_pft            ,& !input  :nodule P:C ratio, [gP gC-1]
    CanopyLeafSheathC_brch    => plt_biom%CanopyLeafSheathC_brch     ,& !input  :plant branch leaf + sheath C, [g d-2]
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft             ,& !input  :threshold zero for plang growth calculation, [-]
    ZERO4LeafVar_pft          => plt_biom%ZERO4LeafVar_pft           ,& !input  :threshold zero for leaf calculation, [-]
    CanopyGrosRCO2_pft        => plt_bgcr%CanopyGrosRCO2_pft         ,& !inoput :canopy plant+nodule autotrophic respiraiton, [gC d-2]
    ECO_ER_col                => plt_bgcr%ECO_ER_col                 ,& !inoput :ecosystem respiration, [g d-2 h-1]
    CanopyRespC_CumYr_pft     => plt_bgcr%CanopyRespC_CumYr_pft      ,& !inoput :total autotrophic respiration, [gC d-2 ]
    Eco_AutoR_CumYr_col       => plt_bgcr%Eco_AutoR_CumYr_col        ,& !inoput :ecosystem autotrophic respiration, [g d-2 h-1]
    CO2NetFix_pft             => plt_bgcr%CO2NetFix_pft              ,& !inoput :canopy net CO2 exchange, [gC d-2 h-1]
    LitrfallElms_pvr          => plt_bgcr%LitrfallElms_pvr           ,& !inoput :plant LitrFall element, [g d-2 h-1]
    CanopyN2Fix_pft           => plt_rbgc%CanopyN2Fix_pft            ,& !output :total canopy N2 fixation, [g d-2 h-1]            
    NodulInfectElms_pft       => plt_bgcr%NodulInfectElms_pft        ,& !inoput :nodule infection chemical element mass, [g d-2]
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch       ,& !inoput :branch nonstructural element, [g d-2]
    CanopyNodulNonstElms_brch => plt_biom%CanopyNodulNonstElms_brch  ,& !inoput :branch nodule nonstructural element, [g d-2]
    RNFixCO2_pft              => plt_bgcr%RNFixCO2_pft               ,& !inoput :CO2 respired by nodules, [gC d-2]
    CanopyNodulStrutElms_brch => plt_biom%CanopyNodulStrutElms_brch   & !inoput :branch nodule structural element, [g d-2]
  )
  call PrintInfo('beg '//subname)
!     iPlantNfixType_pft=N2 fixation: 4,5,6=rapid to slow canopy symbiosis
!
  
  IF(is_canopy_N2fix(iPlantNfixType_pft(NZ)))THEN
    !
    !     INITIAL INFECTION
    !
    !     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
    !     NodulBiomCatInfection=initial bacterial mass at infection
    !     AREA=grid cell area
    !     CNND,rPCNoduler_pft=bacterial N:C,P:C ratio from PFT file
    !
    
    IF(CanopyNodulStrutElms_brch(ielmc,NB,NZ).LE.0.0_r8)THEN    
      NoduleInfectE(ielmc)=NodulBiomCatInfection*AREA3(NU)
      NoduleInfectE(ielmn)=NodulBiomCatInfection*AREA3(NU)*rNCNodule_pft(NZ)
      NoduleInfectE(ielmp)=NodulBiomCatInfection*AREA3(NU)*rPCNoduler_pft(NZ)

      DO NE=1,NumPlantChemElms
        CanopyNodulStrutElms_brch(NE,NB,NZ) = CanopyNodulStrutElms_brch(NE,NB,NZ)+NoduleInfectE(NE)
        NodulInfectElms_pft(NE,NZ)          = NodulInfectElms_pft(NE,NZ)+NoduleInfectE(NE)
      ENDDO
    ENDIF
    !
    !     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
    !     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
    !     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
    !     AND TEMPERATURE
    !
    !     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
    !     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
    !     NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmn),NodulNonstElmConc(ielmp)=nonstructural C,N,P concn in bacteria
    !     CNKI,CPKI=nonstructural N,P inhibition constant on growth
    !     FCNPF=N,P constraint to bacterial activity
    !
    IF(CanopyNodulStrutElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      DO NE=1,NumPlantChemElms
        NodulNonstElmConc(NE)=AZMAX1(CanopyNodulNonstElms_brch(NE,NB,NZ)/CanopyNodulStrutElms_brch(ielmc,NB,NZ))
      ENDDO
    ELSE
      NodulNonstElmConc(1:NumPlantChemElms)=1.0_r8
    ENDIF

    IF(NodulNonstElmConc(ielmc).GT.ZERO)THEN
      CCC=AZMAX1(AMIN1(1.0_r8,safe_adb(NodulNonstElmConc(ielmn),NodulNonstElmConc(ielmn)+NodulNonstElmConc(ielmc)*CNKI) &
        ,safe_adb(NodulNonstElmConc(ielmp),NodulNonstElmConc(ielmp)+NodulNonstElmConc(ielmc)*CPKI)))
      CNC=AZMAX1(AMIN1(1.0_r8,safe_adb(NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmc)+NodulNonstElmConc(ielmn)/CNKI)))
      CPC=AZMAX1(AMIN1(1.0_r8,safe_adb(NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmc)+NodulNonstElmConc(ielmp)/CPKI)))
    ELSE
      CCC=0._r8
      CNC=0._r8
      CPC=0._r8
    ENDIF
    IF(CanopyNodulStrutElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FCNPF=AMIN1(1.0_r8 &
        ,SQRT(CanopyNodulStrutElms_brch(ielmn,NB,NZ)/(CanopyNodulStrutElms_brch(ielmc,NB,NZ)*rNCNodule_pft(NZ))) &
        ,SQRT(CanopyNodulStrutElms_brch(ielmp,NB,NZ)/(CanopyNodulStrutElms_brch(ielmc,NB,NZ)*rPCNoduler_pft(NZ))))
    ELSE
      FCNPF=1.0_r8
    ENDIF
    !MM factor for nodule maintenance respiration
    SPNDLI=NodulNonstElmConc(ielmc)/(NodulNonstElmConc(ielmc)+SPNDLK)
    !
    !     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
    !     NON-STRUCTURAL C:N:P
    !
    !     RespNonst_Oltd=respiration from non-structural C
    !     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
    !     VMXO=specific respiration rate by bacterial N2 fixers
    !     WTNDB=bacterial C mass
    !     fTCanopyGroth_pft=temperature function for canopy growth
    !     FCNPF=N,P constraint to bacterial activity
    !     WFNG=growth function of canopy water potential
    !
    !Cap the temperature- and water-adjusted demand against available nonstructural C.
    RespNonst_Oltd=AZMAX1(AMIN1(CanopyNodulNonstElms_brch(ielmc,NB,NZ),&
      VMXO*CanopyNodulStrutElms_brch(ielmc,NB,NZ)*FCNPF*fTCanopyGroth_pft(NZ)*WFNG))*SPNDLI

    !     CPOOLNX=CanopyNodulNonstElms_brch(ielmc,NB,NZ)
    !     VMXOX=VMXO*CanopyNodulStrutElms_brch(ielmc,NB,NZ)*FCNPF*fTCanopyGroth_pft(NZ)*WFNG
    !
    !     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
    !     NODULE STRUCTURAL N
    !
    !     Rmaint=bacterial maintenance respiration
    !     RmSpecPlant=specific maintenance respiration rate (g C g-1 N h-1)
    !     TFN5=temperature function for canopy maintenance respiration
    !     WTNDBN=bacterial N mass
    !
    Rmaint=AZMAX1(RmSpecPlant*TFN5*CanopyNodulStrutElms_brch(ielmn,NB,NZ))
    !
    !     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
    !     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
    !
    !     RXNDL=difference between non-structural C respn and mntc respn
    !     RGNDL=growth respiration unlimited by N,P
    !     RSNDL=excess maintenance respiration
    !
    RXNDL=RespNonst_Oltd-Rmaint
    RGNDL=AZMAX1(RXNDL)
    RSNDL=AZMAX1(-RXNDL)
    !
    !     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
    !     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
    !     CONSTRAINED BY MICROBIAL N REQUIREMENT
    !
    !     RGN2P=respiration requirement to maintain bacterial N:C ratio
    !     WTNDB,WTNDBN=bacterial C,N mass
    !     rNCNodule_pft=bacterial N:C ratio from PFT file
    !     EN2F=N fixation yield from C oxidation (g N g-1 C)
    !     RGNDL=growth respiration unlimited by N,P
    !     Rauto4Nfix_brch=respiration for N2 fixation
    !     RCanopyN2Fix,CanopyN2Fix_pft=branch,total N2 fixation
    !
    RGN2P=AZMAX1(CanopyNodulStrutElms_brch(ielmc,NB,NZ)*rNCNodule_pft(NZ)-&
      CanopyNodulStrutElms_brch(ielmn,NB,NZ))/EN2F
      
    !N fixation needs to meet with two requirements, 1: nodule bacteria are not N-saturated (RGN2P>0), 2.there is extra respiration for fixation (RGNDL>0)
    IF(RGNDL.GT.ZERO4Groth_pft(NZ))THEN
      Rauto4Nfix_brch=RGNDL*RGN2P/(RGNDL+RGN2P)
    ELSE
      Rauto4Nfix_brch=0._r8
    ENDIF
    RCanopyN2Fix        = Rauto4Nfix_brch*EN2F
    CanopyN2Fix_pft(NZ) = CanopyN2Fix_pft(NZ)+RCanopyN2Fix
    !
    !     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION AND LEAKAGE
    !
    !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !     RCCZN,RCCYN=min,max fractions for bacteria C recycling
    !     RCCXN,RCCQN=max fractions for bacteria N,P recycling
    !     CanopyLeafSheathC_brch=leaf+PetolSheth mass
    !     CCNDLB=bacteria:leaf+PetolSheth ratio
    !     RDNDBX=effect of CCNDLB on bacteria decomposition rate
    !     SPNDX=specific bacterial decomposition rate at current CCNDLB
    !     SPNDL=specific decomposition rate by bacterial N2 fixers
    !     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
    !     NoduleElmDecayLoss(ielmc),NoduleElmDecayLoss(ielmn),NoduleElmDecayLoss(ielmp)=bacterial C,N,P loss from decomposition
    !     NoduleElmntDecay2Litr(ielmc),NoduleElmntDecay2Litr(ielmn),NoduleElmntDecay2Litr(ielmp)=bacterial C,N,P decomposition to LitrFall
    !     NodulElmDecayRecyc(ielmc),NodulElmDecayRecyc(ielmn),NodulElmDecayRecyc(ielmp)=bacterial C,N,P decomposition to recycling
    !
    RCCC  = RCCZN+CCC*RCCYN
    RCCN  = CNC*RCCXN
    RCCP  = CPC*RCCQN
    SPNDX = SPNDL*SQRT(fTCanopyGroth_pft(NZ)*WFNG)
    DO NE=1,NumPlantChemElms
      NoduleElmDecayLoss(NE)=SPNDX*CanopyNodulStrutElms_brch(NE,NB,NZ)
    ENDDO

    NoduleElmntDecay2Litr(ielmc)=NoduleElmDecayLoss(ielmc)*(1.0_r8-RCCC)
    NoduleElmntDecay2Litr(ielmn)=NoduleElmDecayLoss(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
    NoduleElmntDecay2Litr(ielmp)=NoduleElmDecayLoss(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)

    DO NE=1,NumPlantChemElms
      NodulElmDecayRecyc(NE)=NoduleElmDecayLoss(NE)-NoduleElmntDecay2Litr(NE)
    ENDDO
    !
    !     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
    !     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
    !     ENTERED IN 'READQ'
    !
    !     NodulUseNonstElms_brch(ielmc)=total non-structural C used in bacterial growth and growth respiration
    !     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
    !     Rmaint=bacterial maintenance respiration
    !     RespNonst_Oltd=respiration from non-structural C
    !     NodulElmDecayRecyc(ielmc)=bacterial C decomposition to recycling
    !     RGNDL=growth respiration ltd by O2
    !     Rauto4Nfix_brch=respiration for N2 fixation
    !     NoduleBiomCGrowth=bacterial growth
    !     NoduGrowthYield_pft=bacterial growth yield
    !     NoduleCResp=bacterial respiration for growth and N2 fixation
    !     NodulUseNonstElms_brch(ielmn),NodulUseNonstElms_brch(ielmp)=nonstructural N,P used in growth
    !     CNND,rPCNoduler_pft=bacterial N:C,P:C ratio from PFT file
    !     NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmn),NodulNonstElmConc(ielmp)=nonstructural C,N,P concn in bacteria
    !     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
    !
    NodulUseNonstElms_brch(ielmc)=AMIN1(CanopyNodulNonstElms_brch(ielmc,NB,NZ)-AMIN1(Rmaint,RespNonst_Oltd) &
      -Rauto4Nfix_brch+NodulElmDecayRecyc(ielmc),(RGNDL-Rauto4Nfix_brch)/(1.0_r8-NoduGrowthYield_pft(NZ)))
    NoduleBiomCGrowth=NodulUseNonstElms_brch(ielmc)*NoduGrowthYield_pft(NZ)
    NoduleCResp=Rauto4Nfix_brch+NodulUseNonstElms_brch(ielmc)*(1.0_r8-NoduGrowthYield_pft(NZ))
    NodulUseNonstElms_brch(ielmn)=AZMAX1(AMIN1(CanopyNodulNonstElms_brch(ielmn,NB,NZ),NoduleBiomCGrowth*rNCNodule_pft(NZ))) &
      *NodulNonstElmConc(ielmn)/(NodulNonstElmConc(ielmn)+CZKM)
    NodulUseNonstElms_brch(ielmp)=AZMAX1(AMIN1(CanopyNodulNonstElms_brch(ielmp,NB,NZ),NoduleBiomCGrowth*rPCNoduler_pft(NZ))) &
      *NodulNonstElmConc(ielmp)/(NodulNonstElmConc(ielmp)+CPKM)
    !
    !     NODULE SENESCENCE
    !
    !     RSNDL=excess maintenance respiration
    !     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
    !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !     NodulELmLoss2Senes(ielmc),NodulELmLoss2Senes(ielmc),NodulELmLoss2Senes(ielmp)=bacterial C,N,P loss from senescence
    !     NodulELmSenes2Litr(ielmc),NodulELmSenes2Litr(ielmc),NodulELmSenes2Litr(ielmp)=bacterial C,N,P senescence to LitrFall
    !     NodulELmSenes2Recyc(ielmc),NodulELmSenes2Recyc(ielmc),NodulELmSenes2Recyc(ielmp)=bacterial C,N,P senescence to recycling
    !
    IF(RSNDL.GT.0.0_r8 .AND. CanopyNodulStrutElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ) .AND. RCCC.GT.ZERO)THEN
      NodulELmLoss2Senes(ielmc)=RSNDL/RCCC
      NodulELmLoss2Senes(ielmn)=NodulELmLoss2Senes(ielmc)*CanopyNodulStrutElms_brch(ielmn,NB,NZ)&
        /CanopyNodulStrutElms_brch(ielmc,NB,NZ)
      NodulELmLoss2Senes(ielmp)=NodulELmLoss2Senes(ielmc)*CanopyNodulStrutElms_brch(ielmp,NB,NZ)&
        /CanopyNodulStrutElms_brch(ielmc,NB,NZ)

      NodulELmSenes2Litr(ielmc)=NodulELmLoss2Senes(ielmc)*(1.0_r8-RCCC)
      NodulELmSenes2Litr(ielmn)=NodulELmLoss2Senes(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
      NodulELmSenes2Litr(ielmp)=NodulELmLoss2Senes(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
      DO NE=1,NumPlantChemElms
        NodulELmSenes2Recyc(NE)=NodulELmLoss2Senes(NE)-NodulELmSenes2Litr(NE)
      ENDDO
    ELSE
      NodulELmLoss2Senes(1:NumPlantChemElms)  = 0._r8
      NodulELmSenes2Litr(1:NumPlantChemElms)  = 0._r8
      NodulELmSenes2Recyc(1:NumPlantChemElms) = 0._r8
    ENDIF
    !
    !     TOTAL NODULE RESPIRATION
    !
    !     RCO2T=total C respiration
    !     Rmaint=bacterial maintenance respiration
    !     RespNonst_Oltd=respiration from non-structural C
    !     NoduleCResp=bacterial respiration for growth and N2 fixation
    !     NodulELmSenes2Recyc(ielmc)=bacterial C senescence to recycling
    !     CanopyGrosRCO2_pft,CanopyRespC_CumYr_pft=total,above-ground PFT respiration
    !     CO2NetFix_pft=PFT net CO2 fixation
    !     ECO_ER_col=ecosystem respiration
    !     Eco_AutoR_CumYr_col=total autotrophic respiration
    !
    RCO2T                     = AMIN1(Rmaint,RespNonst_Oltd)+NoduleCResp+NodulELmSenes2Recyc(ielmc)
    RNFixCO2_pft(NZ)          = RNFixCO2_pft(NZ)+RCO2T
    CanopyGrosRCO2_pft(NZ)    = CanopyGrosRCO2_pft(NZ)-RCO2T
    CanopyRespC_CumYr_pft(NZ) = CanopyRespC_CumYr_pft(NZ)-RCO2T
    CO2NetFix_pft(NZ)         = CO2NetFix_pft(NZ)-RCO2T
    ECO_ER_col                = ECO_ER_col-RCO2T
    Eco_AutoR_CumYr_col       = Eco_AutoR_CumYr_col-RCO2T
    !
    !     NODULE LitrFall CAUSED BY REMOBILIZATION
    !
    !     CSNC,ZSNC,PSNC=C,N,P LitrFall from decomposition and senescence
    !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
    !     NoduleElmntDecay2Litr(ielmc),NoduleElmntDecay2Litr(ielmn),NoduleElmntDecay2Litr(ielmp)=bacterial C,N,P decomposition to LitrFall
    !     NodulELmSenes2Litr(ielmc),NodulELmSenes2Litr(ielmc),NodulELmSenes2Litr(ielmp)=bacterial C,N,P senescence to LitrFall
    !
    D6470: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,0,NZ) &
          +PlantElmAllocMat4Litr(NE,ifoliar,M,NZ)*(NoduleElmntDecay2Litr(NE)+NodulELmSenes2Litr(NE))
      ENDDO
    ENDDO D6470
    !
    !     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
    !
    !     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
    !     Rmaint=bacterial maintenance respiration
    !     RespNonst_Oltd=respiration from non-structural C
    !     Rauto4Nfix_brch=respiration for N2 fixation
    !     NodulUseNonstElms_brch(ielmc)=total non-structural C used in bacterial growth and growth respiration
    !     NodulElmDecayRecyc(ielmc),NodulElmDecayRecyc(ielmn),NodulElmDecayRecyc(ielmp)=bacterial C,N,P decomposition to recycling
    !     NodulELmSenes2Recyc(ielmc),NodulELmSenes2Recyc(ielmc),NodulELmSenes2Recyc(ielmp)=bacterial C,N,P senescence to recycling
    !     NodulUseNonstElms_brch(ielmn),NodulUseNonstElms_brch(ielmp)=nonstructural N,P used in growth
    !     RCanopyN2Fix=branch N2 fixation
    !
    CanopyNodulNonstElms_brch(ielmc,NB,NZ)=CanopyNodulNonstElms_brch(ielmc,NB,NZ)-AMIN1(Rmaint,RespNonst_Oltd)-Rauto4Nfix_brch&
      -NodulUseNonstElms_brch(ielmc)+NodulElmDecayRecyc(ielmc)
    CanopyNodulNonstElms_brch(ielmn,NB,NZ)=CanopyNodulNonstElms_brch(ielmn,NB,NZ) &
      -NodulUseNonstElms_brch(ielmn)+NodulElmDecayRecyc(ielmn)+NodulELmSenes2Recyc(ielmn)+RCanopyN2Fix
    CanopyNodulNonstElms_brch(ielmp,NB,NZ)=CanopyNodulNonstElms_brch(ielmp,NB,NZ) &
      -NodulUseNonstElms_brch(ielmp)+NodulElmDecayRecyc(ielmp)+NodulELmSenes2Recyc(ielmp)
!
!     UPDATE STATE VARIABLES FOR NODULE C, N, P
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     NoduleBiomCGrowth=bacterial growth
!     NoduleElmDecayLoss(ielmc),NoduleElmDecayLoss(ielmn),NoduleElmDecayLoss(ielmp)=bacterial C,N,P loss from decomposition
!     NodulELmLoss2Senes(ielmc),NodulELmLoss2Senes(ielmc),NodulELmLoss2Senes(ielmp)=bacterial C,N,P loss from senescence
!     NodulUseNonstElms_brch(ielmn),NodulUseNonstElms_brch(ielmp)=nonstructural N,P used in growth
!
    CanopyNodulStrutElms_brch(ielmc,NB,NZ)=CanopyNodulStrutElms_brch(ielmc,NB,NZ) &
      +NoduleBiomCGrowth-NoduleElmDecayLoss(ielmc)-NodulELmLoss2Senes(ielmc)
    DO NE=2,NumPlantChemElms  
      CanopyNodulStrutElms_brch(NE,NB,NZ)=CanopyNodulStrutElms_brch(NE,NB,NZ) &
        +NodulUseNonstElms_brch(NE)-NoduleElmDecayLoss(NE)-NodulELmLoss2Senes(NE)
    ENDDO  
    !
    !
    !     TRANSFER NON-STRUCTURAL C,N,P BETWEEN BRANCH AND NODULES
    !     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
    !
    !     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
    !     CanopyLeafSheathC_brch=leaf+PetolSheth C mass
    !     WTNDB=bacterial C mass
    !     NodulBiomCatInfection=initial bacterial mass at infection
    !     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
    !     CCNGB=parameter to calculate nonstructural C,N,P exchange
    !     CCNDLB=bacteria:leaf+PetolSheth ratio
    !     XFRE(ielmc),XFRE(ielmn),XFRE(ielmp)=nonstructural C,N,P transfer
    !     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
    !
    IF(CanopyNonstElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ) .AND. CanopyLeafSheathC_brch(NB,NZ).GT.ZERO4LeafVar_pft(NZ))THEN
      CCNDLB = CanopyNodulStrutElms_brch(ielmc,NB,NZ)/CanopyLeafSheathC_brch(NB,NZ)
      WTLSB1 = CanopyLeafSheathC_brch(NB,NZ)
      WTNDB1 = AMIN1(CanopyLeafSheathC_brch(NB,NZ),AMAX1(NodulBiomCatInfection*AREA3(NU),CanopyNodulStrutElms_brch(ielmc,NB,NZ)))
      WTLSBT = WTLSB1+WTNDB1
      IF(WTLSBT.GT.ZERO4Groth_pft(NZ))THEN
        FXRNX=FXRN(iPlantNfixType_pft(NZ))/(1.0_r8+CCNDLB/CCNGB)
!    2/(1.0+CCNDLB/(CCNGB*FXRN(iPlantNfixType_pft(NZ))))
        CPOOLD=(CanopyNonstElms_brch(ielmc,NB,NZ)*WTNDB1-CanopyNodulNonstElms_brch(ielmc,NB,NZ)*WTLSB1)/WTLSBT
        XFRE(ielmc)=FXRNX*CPOOLD
        CanopyNonstElms_brch(ielmc,NB,NZ)=CanopyNonstElms_brch(ielmc,NB,NZ)-XFRE(ielmc)
        CanopyNodulNonstElms_brch(ielmc,NB,NZ)=CanopyNodulNonstElms_brch(ielmc,NB,NZ)+XFRE(ielmc)

        CPOOLT=CanopyNonstElms_brch(ielmc,NB,NZ)+CanopyNodulNonstElms_brch(ielmc,NB,NZ)
        IF(CPOOLT.GT.ZERO4Groth_pft(NZ))THEN
          ZPOOLD=(CanopyNonstElms_brch(ielmn,NB,NZ)*CanopyNodulNonstElms_brch(ielmc,NB,NZ) &
            -CanopyNodulNonstElms_brch(ielmn,NB,NZ)*CanopyNonstElms_brch(ielmc,NB,NZ))/CPOOLT
          XFRE(ielmn)=FXRNX*ZPOOLD
          PPOOLD=(CanopyNonstElms_brch(ielmp,NB,NZ)*CanopyNodulNonstElms_brch(ielmc,NB,NZ) &
            -CanopyNodulNonstElms_brch(ielmp,NB,NZ)*CanopyNonstElms_brch(ielmc,NB,NZ))/CPOOLT
          XFRE(ielmp)=FXRNX*PPOOLD
          DO NE=2,NumPlantChemElms
            CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)-XFRE(NE)
            CanopyNodulNonstElms_brch(NE,NB,NZ)=CanopyNodulNonstElms_brch(NE,NB,NZ)+XFRE(NE)
          ENDDO
        ENDIF

      ENDIF
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine CanopyNoduleBiochemistry

!----------------------------------------------------------------------------------------------------
  subroutine RootNodulBiochemistry(I,J,NZ,TFN6_vr)
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: TFN6_vr(JZ1)
  integer :: L,M
  real(r8) :: ZPOOLD
  real(r8) :: CCC,CNC,CPC
  real(r8) :: CPOOLT
  real(r8) :: CPOOLD
  real(r8) :: NodulNonstElmConc(NumPlantChemElms)
  real(r8) :: NodulNonstElm4Groth_rvr(NumPlantChemElms),CPOOLNX
  real(r8) :: CCNDLR
  real(r8) :: FCNPF
  real(r8) :: FXRNX
  real(r8) :: NoduleBiomCGrowth
  real(r8) :: PPOOLD
  real(r8) :: RCO2T
  real(r8) :: RCO2TM
  real(r8) :: RespNonst_Oltd,Rmaint,RXNDL
  real(r8) :: RGNDL,RSNDL
  real(r8) :: RGN2P,Rauto4Nfix_rvr
  real(r8) :: NoduleElmDecayLoss(NumPlantChemElms)
  real(r8) :: NoduleElmntDecay2Litr(NumPlantChemElms)
  real(r8) :: NodulElmDecayRecyc(NumPlantChemElms)
  real(r8) :: NoduleCResp
  real(r8) :: NodulELmLoss2Senes(NumPlantChemElms)
  real(r8) :: NodulELmSenes2Litr(NumPlantChemElms)
  real(r8) :: NodulELmSenes2Recyc(NumPlantChemElms)
  real(r8) :: RespNonst_OUltd,RXNDLM,RGNDLM
  real(r8) :: RSNDLM
  real(r8) :: NoduleRemobilizableC !structural C available for respiration after background decay
  real(r8) :: NoduleMaintStructC_OUltd,NoduleMaintStructC_Oltd !potential/actual structural C respired for maintenance
  real(r8) :: SPNDLI
  real(r8) :: SPNDX
  real(r8) :: WTRTD1,WTNDL1,WTRTDT
  real(r8) :: XFRE(NumPlantChemElms)
  real(r8) :: NoduleInfectE(NumPlantChemElms)
  real(r8) :: RootNoduleElms_beg(NumPlantChemElms)
  real(r8) :: RootNoduleElms_end(NumPlantChemElms)
  real(r8) :: RCCC,RCCN,RCCP,resp,err
  real(r8) :: NoduleLitrfall(NumPlantChemElms)
  real(r8) :: dELoss
  integer  :: NE
  character(len=*), parameter :: subname='RootNodulBiochemistry'
!     begin_execution
  associate(                                                       &
    NU                      => plt_site%NU                        ,& !input  :current soil surface layer number, [-]
    AREA3                   => plt_site%AREA3                     ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    ZERO                    => plt_site%ZERO                      ,& !input  :threshold zero for numerical stability, [-]
    fTgrowRootP_vr          => plt_pheno%fTgrowRootP_vr           ,& !input  :root layer temperature growth functiom, [-]
    MaxNumRootLays          => plt_site%MaxNumRootLays            ,& !input  :maximum root layer number,[-]
    NoduGrowthYield_pft     => plt_allom%NoduGrowthYield_pft      ,& !input  :nodule growth yield, [g g-1]
    rNCNodule_pft           => plt_allom%rNCNodule_pft            ,& !input  :nodule N:C ratio, [gN gC-1]
    rPCNoduler_pft          => plt_allom%rPCNoduler_pft           ,& !input  :nodule P:C ratio, [gP gC-1]
    k_fine_comp             => pltpar%k_fine_comp                 ,& !input  :fine litter complex id
    iroot                   => pltpar%iroot                       ,& !input  :group id of plant root litter
    fRootGrowPSISense_pvr   => plt_pheno%fRootGrowPSISense_pvr    ,& !input  :water stress to plant root growth, [-]
    RAutoRootO2Limter_rpvr  => plt_rbgc%RAutoRootO2Limter_rpvr    ,& !input  :O2 constraint to root respiration (0-1), [-]
    PopuRootMycoC_pvr       => plt_biom% PopuRootMycoC_pvr        ,& !input  :root layer C, [gC d-2]
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft            ,& !input  :threshold zero for plang growth calculation, [-]
    ZERO4LeafVar_pft        => plt_biom%ZERO4LeafVar_pft          ,& !input  :threshold zero for leaf calculation, [-]
    PlantElmAllocMat4Litr   => plt_soilchem%PlantElmAllocMat4Litr ,& !input  :litter kinetic fraction, [-]
    iPlantNfixType_pft      => plt_morph%iPlantNfixType_pft       ,& !input  :N2 fixation type,[-]
    NMaxRootBotLayer_pft    => plt_morph%NMaxRootBotLayer_pft     ,& !input  :maximum soil layer number for all root axes, [-]
    RootRespPotent_pvr      => plt_rbgc%RootRespPotent_pvr        ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]
    RootCO2EmisPot_pvr      => plt_rbgc%RootCO2EmisPot_pvr        ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootCO2Autor_pvr        => plt_rbgc%RootCO2Autor_pvr          ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    NodulInfectElms_pft     => plt_bgcr%NodulInfectElms_pft       ,& !inoput :nodule infection chemical element mass, [g d-2]
    LitrfallElms_pvr        => plt_bgcr%LitrfallElms_pvr          ,& !inoput :plant LitrFall element, [g d-2 h-1]
    RootN2Fix_pft           => plt_rbgc%RootN2Fix_pft             ,& !inoput :total root N2 fixation, [g d-2 h-1]
    RootNodulStrutElms_rpvr => plt_biom%RootNodulStrutElms_rpvr   ,& !inoput :layered root nodule element, [g d-2]
    RootNodulNonstElms_rpvr => plt_biom%RootNodulNonstElms_rpvr   ,& !inoput :layered root nonstructural element, [g d-2]
    RootMycoNonstElms_rpvr  => plt_biom%RootMycoNonstElms_rpvr    ,& !inoput :layered root nonstructural element, [g d-2]
    RCO2Nodule_pvr          => plt_bgcr%RCO2Nodule_pvr            ,& !inoput :layered root nodule respiration, [gC d-2 h-1]
    RNFixCO2_pft            => plt_bgcr%RNFixCO2_pft              ,& !inoput :CO2 respired by nodules, [gC d-2]    
    RootN2Fix_pvr           => plt_bgcr%RootN2Fix_pvr              & !output :root N2 fixation, [gN d-2 h-1]
  )

  call PrintInfo('beg '//subname)
  !     iPlantNfixType_pft=N2 fixation: 1,2,3=rapid to slow root symbiosis
  !     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
  !     NodulBiomCatInfection=initial bacterial mass at infection
  !     AREA=grid cell area
  !     CNND,rPCNoduler_pft=bacterial N:C,P:C ratio from PFT file
  !
!  call SumRootNoduleBiome(I,J,NZ,RootNoduleElms_beg)
  NoduleLitrfall(:)=0._r8
  IF(is_root_N2fix(iPlantNfixType_pft(NZ)))THEN
    D5400: DO L=NU,MaxNumRootLays
      if(L.GT.NMaxRootBotLayer_pft(NZ))then
        DO NE=1,NumPlantChemElms
          if(RootNodulStrutElms_rpvr(NE,L,NZ).GT.0._r8 .or. RootNodulNonstElms_rpvr(NE,L,NZ).GT.0._r8)then      
            D6371: DO M=1,jsken
              dELoss=PlantElmAllocMat4Litr(NE,iroot,M,NZ)*AZMAX1(RootNodulStrutElms_rpvr(NE,L,NZ)+RootNodulNonstElms_rpvr(NE,L,NZ))
              LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+dELoss
              NoduleLitrfall(NE)=  NoduleLitrfall(NE)+dELoss
            ENDDO D6371
            RootNodulStrutElms_rpvr(NE,L,NZ)=0._r8
            RootNodulNonstElms_rpvr(NE,L,NZ)=0._r8
          endif  
        ENDDO
        cycle
      endif
      IF(PopuRootMycoC_pvr(ipltroot,L,NZ).GT.ZERO4LeafVar_pft(NZ))THEN
        !
        !     INITIAL INFECTION
        !
        IF(RootNodulStrutElms_rpvr(ielmc,L,NZ).LE.0.0_r8)THEN
          NoduleInfectE(ielmc)=NodulBiomCatInfection*AREA3(NU)
          NoduleInfectE(ielmn)=NodulBiomCatInfection*AREA3(NU)*rNCNodule_pft(NZ)
          NoduleInfectE(ielmp)=NodulBiomCatInfection*AREA3(NU)*rPCNoduler_pft(NZ)

          DO NE=1,NumPlantChemElms
            RootNodulStrutElms_rpvr(NE,L,NZ) = RootNodulStrutElms_rpvr(NE,L,NZ)+NoduleInfectE(NE)
            NodulInfectElms_pft(NE,NZ)       = NodulInfectElms_pft(NE,NZ)+NoduleInfectE(NE)
          ENDDO                      
        ENDIF
        !
        !     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
        !     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
        !     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
        !     AND TEMPERATURE
        !
        !     NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmn),NodulNonstElmConc(ielmp)=nonstructural C,N,P concn in bacteria
        !     CNKI,CPKI=nonstructural N,P inhibition constant on growth
        !     FCNPF=N,P constraint to bacterial activity
        !
        IF(RootNodulStrutElms_rpvr(ielmc,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
          DO NE=1,NumPlantChemElms
            NodulNonstElmConc(NE)=AZMAX1(RootNodulNonstElms_rpvr(NE,L,NZ)/RootNodulStrutElms_rpvr(ielmc,L,NZ))
          ENDDO
        ELSE
          NodulNonstElmConc(1:NumPlantChemElms)=1.0_r8
        ENDIF

        IF(NodulNonstElmConc(ielmc).GT.ZERO)THEN
          CCC=AZMAX1(AMIN1(1.0_r8 &
            ,safe_adb(NodulNonstElmConc(ielmn),NodulNonstElmConc(ielmn)+NodulNonstElmConc(ielmc)*CNKI) &
            ,safe_adb(NodulNonstElmConc(ielmp),NodulNonstElmConc(ielmp)+NodulNonstElmConc(ielmc)*CPKI)))
          CNC=AZMAX1(AMIN1(1.0_r8,safe_adb(NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmc)+NodulNonstElmConc(ielmn)/CNKI)))
          CPC=AZMAX1(AMIN1(1.0_r8,safe_adb(NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmc)+NodulNonstElmConc(ielmp)/CPKI)))
        ELSE
          CCC=0._r8
          CNC=0._r8
          CPC=0._r8
        ENDIF

        IF(RootNodulStrutElms_rpvr(ielmc,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
          FCNPF=AMIN1(1.0_r8 &
            ,SQRT(RootNodulStrutElms_rpvr(ielmn,L,NZ)/(RootNodulStrutElms_rpvr(ielmc,L,NZ)*rNCNodule_pft(NZ))) &
            ,SQRT(RootNodulStrutElms_rpvr(ielmp,L,NZ)/(RootNodulStrutElms_rpvr(ielmc,L,NZ)*rPCNoduler_pft(NZ))))
        ELSE
          FCNPF=1.0_r8
        ENDIF
        SPNDLI=NodulNonstElmConc(ielmc)/(NodulNonstElmConc(ielmc)+SPNDLK)
        !
        !     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
        !     NON-STRUCTURAL C:N:P
        !
        !     RespNonst_OUltd=respiration from non-structural C unltd by O2
        !     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
        !     VMXO=specific respiration rate by bacterial N2 fixers
        !     WTNDL=bacterial C mass
        !     fTgrowRootP_vr=temperature function for root growth
        !     FCNPF=N,P constraint to bacterial activity
        !     fRootGrowPSISense=growth function of root water potential
!
        CPOOLNX         = RootNodulNonstElms_rpvr(ielmc,L,NZ)
        !Cap the adjusted demand before applying the physiological oxygen factor below.
        RespNonst_OUltd = AZMAX1(AMIN1(RootNodulNonstElms_rpvr(ielmc,L,NZ) &
          ,VMXO*RootNodulStrutElms_rpvr(ielmc,L,NZ)*FCNPF*fTgrowRootP_vr(L,NZ) &
          *fRootGrowPSISense_pvr(ipltroot,L,NZ)))*SPNDLI
        !
        !     O2-LIMITED NODULE RESPIRATION FROM 'WFR' IN 'UPTAKE'
        !
        !     RespNonst_Oltd=respiration from non-structural C ltd by O2
        !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes, assuming nodule is associated with roots
        !
        RespNonst_Oltd=RespNonst_OUltd*RAutoRootO2Limter_rpvr(ipltroot,L,NZ)
        !
        !     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
        !     NODULE STRUCTURAL N
        !
        !     Rmaint=bacterial maintenance respiration
        !     RmSpecPlant=specific maintenance respiration rate (g C g-1 N h-1)
        !     TFN6_vr=temperature function for root maintenance respiration
        !     WTNDLN=bacterial N mass
        !
        Rmaint=AZMAX1(RmSpecPlant*TFN6_vr(L)*RootNodulStrutElms_rpvr(ielmn,L,NZ))
        !
        !     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
        !     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
        !
        !     RXNDLM,RXNDL=difference between non-structural C respn and mntc respn unltd,ltd by O2
        !     RGNDLM,RGNDL=growth respiration unlimited by N,P and unltd,ltd by O2
        !     RSNDLM,RSNDL=excess maintenance respiration unltd,ltd by O2
        !
        RXNDLM = RespNonst_OUltd-Rmaint
        RXNDL  = RespNonst_Oltd-Rmaint
        RGNDLM = AZMAX1(RXNDLM)
        RGNDL  = AZMAX1(RXNDL)     !respiration available for growth and N fixation
        RSNDLM = AZMAX1(-RXNDLM)
        RSNDL  = AZMAX1(-RXNDL)
        !
        !     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
        !     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
        !     CONSTRAINED BY MICROBIAL N REQUIREMENT
        !
        !     RGN2P=respiration requirement to maintain bacterial N:C ratio
        !     WTNDL,WTNDLN=bacterial C,N mass
        !     rNCNodule_pft=bacterial N:C ratio from PFT file
        !     EN2F=N fixation yield from C oxidation (g N g-1 C)
        !     RGNDL=growth respiration unlimited by N,P
        !     Rauto4Nfix_rvr=respiration for N2 fixation
        !     RootN2Fix_pvr,RootN2Fix_pft=layer,total root N2 fixation
        !
        RGN2P=AZMAX1(RootNodulStrutElms_rpvr(ielmc,L,NZ)*rNCNodule_pft(NZ)-RootNodulStrutElms_rpvr(ielmn,L,NZ))/EN2F
        IF(RGNDL.GT.ZERO4Groth_pft(NZ))THEN
          Rauto4Nfix_rvr=RGNDL*RGN2P/(RGNDL+RGN2P)
        ELSE
          Rauto4Nfix_rvr=0._r8
        ENDIF
        RootN2Fix_pvr(L,NZ) = Rauto4Nfix_rvr*EN2F
        RootN2Fix_pft(NZ)   = RootN2Fix_pft(NZ)+RootN2Fix_pvr(L,NZ)
        
        !
        !     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION
        !
        !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
        !     RCCZN,RCCYN=min,max fractions for bacteria C recycling
        !     RCCXN,RCCQN=max fractions for bacteria N,P recycling
        !     WTRTD=root C mass
        !     CCNDLR=bacteria:root ratio
        !     RDNDLX=effect of CCNDLR on bacteria decomposition rate
        !     CCNKR=Km for bacterial vs root mass in decomposition
        !     SPNDX=specific bacterial decomposition rate at current CCNDLR
        !     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
        !     NoduleElmDecayLoss(ielmc),NoduleElmDecayLoss(ielmn),NoduleElmDecayLoss(ielmp)=bacterial C,N,P loss from decomposition
        !
        RCCC  = RCCZN+CCC*RCCYN
        RCCN  = CNC*RCCXN
        RCCP  = CPC*RCCQN
        SPNDX = SPNDL*SQRT(fTgrowRootP_vr(L,NZ)*fRootGrowPSISense_pvr(ipltroot,L,NZ))
        DO NE=1,NumPlantChemElms
          NoduleElmDecayLoss(NE)=SPNDX*RootNodulStrutElms_rpvr(NE,L,NZ)
        ENDDO

        NoduleElmntDecay2Litr(ielmc)=NoduleElmDecayLoss(ielmc)*(1.0_r8-RCCC)
        NoduleElmntDecay2Litr(ielmn)=NoduleElmDecayLoss(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
        NoduleElmntDecay2Litr(ielmp)=NoduleElmDecayLoss(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
        DO NE=1,NumPlantChemElms
          NodulElmDecayRecyc(NE)=NoduleElmDecayLoss(NE)-NoduleElmntDecay2Litr(NE)
        ENDDO
        !
        !     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
        !     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
        !     ENTERED IN 'READQ'
        !
        !     NodulNonstElm4Groth_rvr(ielmc)=total non-structural C used in bacterial growth and growth respiration
        !     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
        !     Rmaint=bacterial maintenance respiration
        !     RespNonst_Oltd=respiration from non-structural C
        !     NodulElmDecayRecyc(ielmc)=bacterial C decomposition to recycling
        !     RGNDL=growth respiration ltd by O2
        !     Rauto4Nfix_brch=respiration for N2 fixation
        !     NoduleBiomCGrowth=bacterial growth
        !     NoduGrowthYield_pft=bacterial growth yield
        !     NoduleCResp=bacterial respiration for growth and N2 fixation
        !     NodulNonstElm4Groth_rvr(ielmn),NodulNonstElm4Groth_rvr(ielmp)=nonstructural N,P used in growth
        !     CNND,rPCNoduler_pft=bacterial N:C,P:C ratio from PFT file
        !     NodulNonstElmConc(ielmc),NodulNonstElmConc(ielmn),NodulNonstElmConc(ielmp)=nonstructural C,N,P concn in bacteria
        !     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
        !
        NodulNonstElm4Groth_rvr(ielmc)=AMIN1(RootNodulNonstElms_rpvr(ielmc,L,NZ)-AMIN1(Rmaint,RespNonst_Oltd) &
          -Rauto4Nfix_rvr+NodulElmDecayRecyc(ielmc),(RGNDL-Rauto4Nfix_rvr)/(1.0_r8-NoduGrowthYield_pft(NZ)))
        NoduleBiomCGrowth=NodulNonstElm4Groth_rvr(ielmc)*NoduGrowthYield_pft(NZ)
        NoduleCResp=Rauto4Nfix_rvr+NodulNonstElm4Groth_rvr(ielmc)*(1.0_r8-NoduGrowthYield_pft(NZ))
        NodulNonstElm4Groth_rvr(ielmn)=AZMAX1(AMIN1(RootNodulNonstElms_rpvr(ielmn,L,NZ),NoduleBiomCGrowth*rNCNodule_pft(NZ)))&
          *NodulNonstElmConc(ielmn)/(NodulNonstElmConc(ielmn)+CZKM)
        NodulNonstElm4Groth_rvr(ielmp)=AZMAX1(AMIN1(RootNodulNonstElms_rpvr(ielmp,L,NZ),NoduleBiomCGrowth*rPCNoduler_pft(NZ)))&
          *NodulNonstElmConc(ielmp)/(NodulNonstElmConc(ielmp)+CPKM)
        !
        !     NODULE SENESCENCE
        !
        !     RSNDL=excess maintenance respiration
        !     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
        !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
        !     NodulELmLoss2Senes(:)=bacterial C,N,P loss from senescence
        !     NodulELmSenes2Litr(:)=bacterial C,N,P senescence to LitrFall
        !     NodulELmSenes2Recyc(:)=bacterial C,N,P senescence to recycling
        !
        NoduleMaintStructC_OUltd=0._r8
        NoduleMaintStructC_Oltd=0._r8
        IF(RootNodulStrutElms_rpvr(ielmc,L,NZ).GT.ZERO4Groth_pft(NZ) .AND. RCCC.GT.ZERO)THEN
          !Background decay and maintenance must not consume the same structural C.
          NoduleRemobilizableC     = AZMAX1(RootNodulStrutElms_rpvr(ielmc,L,NZ)-NoduleElmDecayLoss(ielmc))*RCCC
          NoduleMaintStructC_OUltd = AMIN1(RSNDLM,NoduleRemobilizableC)

          !Apply the physiological oxygen factor to structural respiration as for other roots.
          NoduleMaintStructC_Oltd=AMIN1(RSNDL,NoduleRemobilizableC)*RAutoRootO2Limter_rpvr(ipltroot,L,NZ)

          !Derive all actual C/N/P losses, litter, and recycling from this one payment.
          NodulELmLoss2Senes(ielmc) = NoduleMaintStructC_Oltd/RCCC
          NodulELmLoss2Senes(ielmn) = NodulELmLoss2Senes(ielmc)*RootNodulStrutElms_rpvr(ielmn,L,NZ)/RootNodulStrutElms_rpvr(ielmc,L,NZ)
          NodulELmLoss2Senes(ielmp) = NodulELmLoss2Senes(ielmc)*RootNodulStrutElms_rpvr(ielmp,L,NZ)/RootNodulStrutElms_rpvr(ielmc,L,NZ)

          NodulELmSenes2Litr(ielmc)=NodulELmLoss2Senes(ielmc)*(1.0_r8-RCCC)
          NodulELmSenes2Litr(ielmn)=NodulELmLoss2Senes(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
          NodulELmSenes2Litr(ielmp)=NodulELmLoss2Senes(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
          DO NE=1,NumPlantChemElms
            NodulELmSenes2Recyc(NE)=NodulELmLoss2Senes(NE)-NodulELmSenes2Litr(NE)
          ENDDO
        ELSE
          NodulELmLoss2Senes(1:NumPlantChemElms)=0._r8
          NodulELmSenes2Litr(1:NumPlantChemElms)=0._r8
          NodulELmSenes2Recyc(1:NumPlantChemElms)=0._r8
        ENDIF
        !
        !     TOTAL NODULE RESPIRATION is added to root respiration
        !
        !Potential respiration uses its own structural payment, independent of oxygen stress.
        RCO2TM                            = AMIN1(Rmaint,RespNonst_OUltd)+RGNDLM+NoduleMaintStructC_OUltd
        RCO2T                             = AMIN1(Rmaint,RespNonst_Oltd)+NoduleCResp+NodulELmSenes2Recyc(ielmc)
        RootRespPotent_pvr(ipltroot,L,NZ) = RootRespPotent_pvr(ipltroot,L,NZ)+RCO2TM
        RootCO2EmisPot_pvr(ipltroot,L,NZ) = RootCO2EmisPot_pvr(ipltroot,L,NZ)+RCO2T
        RootCO2Autor_pvr(ipltroot,L,NZ)   = RootCO2Autor_pvr(ipltroot,L,NZ)-RCO2T
        RCO2Nodule_pvr(L,NZ)              = RCO2Nodule_pvr(L,NZ)+RCO2T
        RNFixCO2_pft(NZ)                  = RNFixCO2_pft(NZ)+RCO2T
        !
        !     NODULE LitrFall CAUSED BY REMOBILIZATION
        !
        !     CSNC,ZSNC,PSNC=C,N,P LitrFall from decomposition and senescence
        !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
        !     NoduleElmntDecay2Litr(:)=bacterial C,N,P decomposition to LitrFall
        !     NodulELmSenes2Litr(:)=bacterial C,N,P senescence to LitrFall
        !
        D6370: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            dELoss=PlantElmAllocMat4Litr(NE,iroot,M,NZ)*AZMAX1(NoduleElmntDecay2Litr(NE)+NodulELmSenes2Litr(NE))
            LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)=LitrfallElms_pvr(NE,M,k_fine_comp,L,NZ)+dELoss
            NoduleLitrfall(NE)=  NoduleLitrfall(NE)+dELoss
          ENDDO
        ENDDO D6370
        !
        !     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
        !
        !     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
        !     Rmaint=bacterial maintenance respiration
        !     RespNonst_Oltd=respiration from non-structural C
        !     Rauto4Nfix_rvr=respiration for N2 fixation
        !     NodulNonstElm4Groth_rvr(ielmc)=total non-structural C used in bacterial growth and growth respiration
        !     NodulElmDecayRecyc(:)=bacterial C,N,P decomposition to recycling
        !     NodulELmSenes2Recyc(:)=bacterial C,N,P senescence to recycling
        !     NodulNonstElm4Groth_rvr(:)=nonstructural N,P used in growth
        !     RootN2Fix_pvr=root N2 fixation
        !
        Do NE=1,NumPlantChemElms
          RootNodulNonstElms_rpvr(NE,L,NZ)=RootNodulNonstElms_rpvr(NE,L,NZ)-NodulNonstElm4Groth_rvr(NE)+NodulElmDecayRecyc(NE)
        ENDDO
        RootNodulNonstElms_rpvr(ielmc,L,NZ)=RootNodulNonstElms_rpvr(ielmc,L,NZ)-AMIN1(Rmaint,RespNonst_Oltd)-Rauto4Nfix_rvr
        RootNodulNonstElms_rpvr(ielmn,L,NZ)=RootNodulNonstElms_rpvr(ielmn,L,NZ)+NodulELmSenes2Recyc(ielmn)+RootN2Fix_pvr(L,NZ)
        RootNodulNonstElms_rpvr(ielmp,L,NZ)=RootNodulNonstElms_rpvr(ielmp,L,NZ)+NodulELmSenes2Recyc(ielmp)
        !
        !     UPDATE STATE VARIABLES FOR NODULE C, N, P
        !
        !     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
        !     NoduleBiomCGrowth=bacterial growth
        !     NoduleElmDecayLoss(:)=bacterial C,N,P loss from decomposition
        !     NodulELmLoss2Senes(:)=bacterial C,N,P loss from senescence
        !     NodulNonstElm4Groth_rvr(ielmn),NodulNonstElm4Groth_rvr(ielmp)=nonstructural N,P used in growth
        !
        DO NE=1,NumPlantChemElms  
          RootNodulStrutElms_rpvr(NE,L,NZ)=RootNodulStrutElms_rpvr(NE,L,NZ)-NoduleElmDecayLoss(NE)-NodulELmLoss2Senes(NE)
        ENDDO  

        RootNodulStrutElms_rpvr(ielmc,L,NZ)=RootNodulStrutElms_rpvr(ielmc,L,NZ)+NoduleBiomCGrowth
        DO NE=2,NumPlantChemElms  
          RootNodulStrutElms_rpvr(NE,L,NZ)=RootNodulStrutElms_rpvr(NE,L,NZ)+NodulNonstElm4Groth_rvr(NE)
        ENDDO  
        !
        !     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND NODULES
        !     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
        !
        !     CPOOLR,ZPOOLR,PPOOLR=root non-structural C,N,P mass
        !     WTRTD=root C mass
        !     WTNDL=bacterial C mass
        !     NodulBiomCatInfection=initial bacterial mass at infection
        !     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
        !     CCNGR=parameter to calculate nonstructural C,N,P exchange
        !     CCNDLR=bacteria:root ratio
        !     XFRE(ielmc),XFRE(ielmn),XFRE(ielmp)=nonstructural C,N,P transfer
        !     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
        !
        IF(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ).GT.ZERO4Groth_pft(NZ) &
          .AND. PopuRootMycoC_pvr(ipltroot,L,NZ).GT.ZERO4LeafVar_pft(NZ))THEN
          CCNDLR = RootNodulStrutElms_rpvr(ielmc,L,NZ)/PopuRootMycoC_pvr(ipltroot,L,NZ)
          WTRTD1 = PopuRootMycoC_pvr(ipltroot,L,NZ)
          WTNDL1 = AMIN1(PopuRootMycoC_pvr(ipltroot,L,NZ),AMAX1(NodulBiomCatInfection*AREA3(NU),RootNodulStrutElms_rpvr(ielmc,L,NZ)))
          WTRTDT = WTRTD1+WTNDL1
          IF(WTRTDT.GT.ZERO4Groth_pft(NZ))THEN
            FXRNX=FXRN(iPlantNfixType_pft(NZ))/(1.0_r8+CCNDLR/CCNGR)  ! /(1.0+CCNDLR/(CCNGR*FXRN(iPlantNfixType_pft(NZ))))
            CPOOLD=(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)*WTNDL1-RootNodulNonstElms_rpvr(ielmc,L,NZ)*WTRTD1)/WTRTDT
            XFRE(ielmc)=FXRNX*CPOOLD
            RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)-XFRE(ielmc)
            RootNodulNonstElms_rpvr(ielmc,L,NZ)         = RootNodulNonstElms_rpvr(ielmc,L,NZ)+XFRE(ielmc)
            CPOOLT                                      = RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)+RootNodulNonstElms_rpvr(ielmc,L,NZ)

            IF(CPOOLT.GT.ZERO4Groth_pft(NZ))THEN
              ZPOOLD=(RootMycoNonstElms_rpvr(ielmn,ipltroot,L,NZ)*RootNodulNonstElms_rpvr(ielmc,L,NZ) &
                -RootNodulNonstElms_rpvr(ielmn,L,NZ)* RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRE(ielmn)=FXRNX*ZPOOLD
              PPOOLD=(RootMycoNonstElms_rpvr(ielmp,ipltroot,L,NZ)*RootNodulNonstElms_rpvr(ielmc,L,NZ) &
                -RootNodulNonstElms_rpvr(ielmp,L,NZ)* RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRE(ielmp)=FXRNX*PPOOLD
              DO NE=2,NumPlantChemElms
                RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)-XFRE(NE)
                RootNodulNonstElms_rpvr(NE,L,NZ)         = RootNodulNonstElms_rpvr(NE,L,NZ)+XFRE(NE)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO D5400
  ENDIF

  call PrintInfo('end '//subname)
  return
  call SumRootNoduleBiome(I,J,NZ,RootNoduleElms_end)

  end associate
  end subroutine RootNodulBiochemistry
!----------------------------------------------------------------------------------------------------
  subroutine SumRootNoduleBiome(I,J,NZ,RootNoduleElms_pft)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(out):: RootNoduleElms_pft(NumPlantChemElms)
  integer :: NE

  associate(                                                       &
    NU                      => plt_site%NU                        ,& !input  :current soil surface layer number, [-]  
    NMaxRootBotLayer_pft    => plt_morph%NMaxRootBotLayer_pft     ,& !input  :maximum soil layer number for all root axes, [-]
    RootNodulStrutElms_rpvr => plt_biom%RootNodulStrutElms_rpvr   ,& !input  :layered root nodule element, [g d-2]    
    RootMycoNonstElms_rpvr  => plt_biom%RootMycoNonstElms_rpvr    ,& !input  :layered root nonstructural element, [g d-2]    
    RootNodulNonstElms_rpvr => plt_biom%RootNodulNonstElms_rpvr    & !inoput :layered root nonstructural element, [g d-2]
  ) 

  DO NE=1, NumPlantChemElms
    RootNoduleElms_pft(NE)  = sum(RootNodulStrutElms_rpvr(NE,NU:NMaxRootBotLayer_pft(NZ),NZ))+sum(RootNodulNonstElms_rpvr(NE,NU:NMaxRootBotLayer_pft(NZ),NZ)) &
      + sum(RootMycoNonstElms_rpvr(NE,ipltroot,NU:NMaxRootBotLayer_pft(NZ),NZ))
  ENDDO    

  end associate
  end subroutine SumRootNoduleBiome

  ![tail]
end module NoduleBGCMod
