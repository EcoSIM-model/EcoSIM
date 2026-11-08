  module stomatesMod
  use data_kind_mod, only: r8 => DAT_KIND_R8, yearIJ_type
  use DebugToolMod,  only: PrintInfo
  use EcosimConst
  use minimathmod
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantPhotosynthesisAPIData, only : plt_photo
  use PlantRadiationAPIData, only : plt_rad
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantBiomassAPIData, only : plt_biom
  use PlantEnergyWaterAPIData, only : plt_ew
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantBGCPars
  use EcoSIMCtrlMod , only : etimer,lverb 
  implicit none

  private
  character(len=*), parameter :: mod_filename = &
  __FILE__


  real(r8), parameter :: Hours2KillAnuals(0:5)=real((/336.0,672.0,672.0,672.0,672.0,672.0/),r8)  !number of hours with no grain fill to terminate annuals

  public :: StomatalDynamics
  public :: PhotosynsDiag
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine PhotosynsDiag(I,J,NZ)  
  !
  !Description:
  !Diagnose photosynthesis related variables
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='PhotosynsDiag'
  real(r8) :: BundlSheathRubiscoC, BundlSheathChlC,MesophyllChlC, MesophyllPEPC, MesophyllRubiscoC
  real(r8) :: tLeafArea,tLeafC
  integer :: K,NB

  associate(                                                                   &
    ZERO                            => plt_site%ZERO                          ,& !input  :threshold zero for numerical stability, [-]  
    NumOfBranches_pft               => plt_morph%NumOfBranches_pft            ,& !input  :number of branches,[-]
    VmaxSpecRubCarboxyRef_pft       => plt_photo%VmaxSpecRubCarboxyRef_pft    ,& !input  :rubisco carboxylase activity at 25 oC, [umol g-1 h-1 m-2]    
    VmaxPEPCarboxyRef_pft           => plt_photo%VmaxPEPCarboxyRef_pft        ,& !input  :PEP carboxylase activity at 25 oC [umol g-1 h-1]    
    ProteinCperm2LeafArea_node      => plt_photo%ProteinCperm2LeafArea_node   ,& !output :Protein C per m2 of leaf aera, [gC (leaf area m-2)]        
    iPlantPhotosynsType_pft         => plt_photo%iPlantPhotosynsType_pft      ,& !input  :plant photosynthetic type (C3 or C4),[-]      
    VmaxRubOxyRef_pft               => plt_photo%VmaxRubOxyRef_pft            ,& !input  :rubisco oxygenase activity at 25 oC, [umol g-1 h-1]    
    iPlantPhenolType_pft            => plt_pheno%iPlantPhenolType_pft         ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]  
    Hours4Leafout_brch              => plt_pheno%Hours4Leafout_brch           ,& !input  :heat requirement for spring leafout/dehardening, [h]    
    Hours4LeafOff_brch              => plt_pheno%Hours4LeafOff_brch           ,& !input  :cold requirement for autumn leafoff/hardening, [h]    
    HourReq4LeafOff_brch            => plt_pheno%HourReq4LeafOff_brch         ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]    
    HourReq4LeafOut_brch            => plt_pheno%HourReq4LeafOut_brch         ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    ZERO4Groth_pft                  => plt_biom%ZERO4Groth_pft                ,& !input  :threshold zero for plang growth calculation, [-]    
    LeafElmntNode_brch              => plt_biom%LeafElmntNode_brch            ,& !input  :branch leaf element, [g d-2]    
    LeafPEP2Protein_pft             => plt_photo%LeafPEP2Protein_pft          ,& !input  :leaf PEP carboxylase content, [gC gC-1]    
    fMesophyllChlProtein_pft        => plt_photo%fMesophyllChlProtein_pft     ,& !input  :fraction of Chl-bound protein in mesophyll cell, [gC gC-1]    
    LeafProtein2Chl_pft             => plt_photo%LeafProtein2Chl_pft          ,& !input  :fraction of leaf protein that is chlorophyll-binded, [gC gC-1]    
    LeafRubisco2Protein_pft         => plt_photo%LeafRubisco2Protein_pft      ,& !input  :leaf rubisco content, [gC gC-1]    
    SpecLeafChlAct_pft              => plt_photo%SpecLeafChlAct_pft           ,& !input  :cholorophyll activity at 25 oC, [umol g-1 h-1]    
    LeafArea_node                   => plt_morph%LeafArea_node                ,& !input  :node leaf area, [m2 d-2]    
    LeafProteinC_node               => plt_biom%LeafProteinC_node             ,& !input  :node leaf protein N, [g d-2]    
    VcMaxPEPCarboxyRef_brch         => plt_photo%VcMaxPEPCarboxyRef_brch      ,& !output :reference maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]          
    VoMaxRubiscoRef_brch            => plt_photo%VoMaxRubiscoRef_brch         ,& !output :maximum rubisco oxygenation rate at reference temperature, [umol g-1 h-1]
    VcMaxRubiscoRef_brch            => plt_photo%VcMaxRubiscoRef_brch         ,& !output :maximum rubisco carboxylation rate at reference temperature, [umol g-1 h-1]
    CanopyVcMaxRubisco25C_pft       => plt_photo%CanopyVcMaxRubisco25C_pft    ,& !output :Canopy VcMax for rubisco carboxylation, [umol s-1 m-2]
    CanopyVoMaxRubisco25C_pft       => plt_photo%CanopyVoMaxRubisco25C_pft    ,& !output :Canopy VoMax for rubisco oxygenation, [umol s-1 m-2]
    CanopyVcMaxPEP25C_pft           => plt_photo%CanopyVcMaxPEP25C_pft        ,& !output :Canopy VcMax in PEP C4 fixation, [umol s-1 m-2]
    ElectronTransptJmax25C_pft      => plt_photo%ElectronTransptJmax25C_pft   ,& !output :Canopy Jmax at reference temperature, [umol e- s-1 m-2]
    LeafProteinC_brch               => plt_biom%LeafProteinC_brch             ,& !output :Protein C for the branches, [gC protein d-2]
    LeafProteinCperm2LA_pft         => plt_biom%LeafProteinCperm2LA_pft       ,& !output :Protein C for the plant, [gC protein m-2 leaf area]
    LeafC3ChlC_brch                 => plt_biom%LeafC3ChlC_brch               ,& !output :Bundle sheath C4/mesophyll C3 chlorophyll C for the branches, [gC chlorophyll d-2]     
    LeafC4ChlC_brch                 => plt_biom%LeafC4ChlC_brch               ,& !output :Mesophyll chlorophyll C for the branches, [gC chlorophyll d-2]     
    LeafRubiscoC_brch               => plt_biom%LeafRubiscoC_brch             ,& !output :Bundle sheath C4/mesophyll C3 Rubisco C for the branches, [gC Rubisco d-2]
    LeafPEPC_brch                   => plt_biom%LeafPEPC_brch                 ,& !output :PEP C for the branches, [gC PEP d-2]
    LeafC3ChlCperm2LA_pft           => plt_biom%LeafC3ChlCperm2LA_pft         ,& !output :Bundle sheath C4/mesophyll C3 chlorophyll C for the branches, [gC chlorophyll d-2]     
    LeafC4ChlCperm2LA_pft           => plt_biom%LeafC4ChlCperm2LA_pft         ,& !output :Mesophyll chlorophyll C for the branches, [gC chlorophyll d-2]     
    LeafRubiscoCperm2LA_pft         => plt_biom%LeafRubiscoCperm2LA_pft       ,& !output :Bundle sheath C4/mesophyll C3 Rubisco C for the branches, [gC Rubisco d-2]
    LeafPEPCperm2LA_pft             => plt_biom%LeafPEPCperm2LA_pft           ,& !output :PEP C for the branches, [gC PEP d-2]
    SpecificLeafArea_pft            => plt_biom%SpecificLeafArea_pft          ,& !ouptut :specifc leaf area per g C of leaf mass, [m2 leaf area (gC leaf C)-1]
    ElectronTransptJmaxRef_brch     => plt_photo%ElectronTransptJmaxRef_brch   & !output :Branch Maximum electron transport rate at reference temperature, [umol e- s-1]    
  )
  call PrintInfo('beg '//subname)
  tLeafArea=0._r8;tLeafC=0._r8
  LeafProteinCperm2LA_pft(NZ)    = 0._r8
  CanopyVcMaxRubisco25C_pft(NZ)  = 0._r8
  CanopyVoMaxRubisco25C_pft(NZ)  = 0._r8
  CanopyVcMaxPEP25C_pft(NZ)      = 0._r8
  ElectronTransptJmax25C_pft(NZ) = 0._r8
  LeafRubiscoCperm2LA_pft(NZ)    = 0._r8
  LeafPEPCperm2LA_pft(NZ)        = 0._r8
  LeafC3ChlCperm2LA_pft(NZ)      = 0._r8
  LeafC4ChlCperm2LA_pft(NZ)      = 0._r8

  DO NB=1,NumOfBranches_pft(NZ)
!
!     FEEDBACK ON CO2 FIXATION
!
!     iPlantPhenolType_pft=phenology type from PFT file
!     Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
!
    VcMaxRubiscoRef_brch(NB,NZ)        = 0._r8
    VoMaxRubiscoRef_brch(NB,NZ)        = 0._r8
    VcMaxPEPCarboxyRef_brch(NB,NZ)     = 0._r8
    ElectronTransptJmaxRef_brch(NB,NZ) = 0._r8
    LeafProteinC_brch(NB,NZ)           = 0._r8
    LeafRubiscoC_brch(NB,NZ)           = 0._r8
    LeafPEPC_brch(NB,NZ)               = 0._r8
    LeafC3ChlC_brch(NB,NZ)             = 0._r8
    LeafC4ChlC_brch(NB,NZ)             = 0._r8

    IF(iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen &
      .OR. Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ) &
      .OR. Hours4LeafOff_brch(NB,NZ).LT.HourReq4LeafOff_brch(NB,NZ))THEN
      
      DO K=1,MaxNodesPerBranch1          
        IF(LeafArea_node(K,NB,NZ).GT.ZERO4Groth_pft(NZ) .AND. LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
          LeafProteinC_brch(NB,NZ)            = LeafProteinC_brch(NB,NZ)+LeafProteinC_node(K,NB,NZ)
          tLeafArea                           = tLeafArea+LeafArea_node(K,NB,NZ)
          tLeafC                              = tLeafC+LeafElmntNode_brch(ielmc,K,NB,NZ)
          ProteinCperm2LeafArea_node(K,NB,NZ) = LeafProteinC_node(K,NB,NZ)/LeafArea_node(K,NB,NZ)

          IF(ProteinCperm2LeafArea_node(K,NB,NZ).GT.ZERO)THEN
            IF(iPlantPhotosynsType_pft(NZ).EQ.ic4_photo)THEN        
              BundlSheathRubiscoC = LeafProteinC_node(K,NB,NZ)*LeafRubisco2Protein_pft(NZ)
              BundlSheathChlC     = LeafProteinC_node(K,NB,NZ)*LeafProtein2Chl_pft(NZ)*(1._r8-fMesophyllChlProtein_pft(NZ))/5.5_r8        !1 gC chl for 5.5 gC protein
              MesophyllPEPC       = LeafProteinC_node(K,NB,NZ)*LeafPEP2Protein_pft(NZ)
              MesophyllChlC       = LeafProteinC_node(K,NB,NZ)*LeafProtein2Chl_pft(NZ)*fMesophyllChlProtein_pft(NZ)/3.7_r8                !1 gC chl for 3.7 gC protein

              VcMaxPEPCarboxyRef_brch(NB,NZ) = VcMaxPEPCarboxyRef_brch(NB,NZ) +VmaxPEPCarboxyRef_pft(NZ)*MesophyllPEPC
              VcMaxRubiscoRef_brch(NB,NZ)    = VcMaxRubiscoRef_brch(NB,NZ)+VmaxSpecRubCarboxyRef_pft(NZ)*BundlSheathRubiscoC
              VoMaxRubiscoRef_brch(NB,NZ)    = VoMaxRubiscoRef_brch(NB,NZ)+VmaxRubOxyRef_pft(NZ)*BundlSheathRubiscoC

              LeafPEPC_brch(NB,NZ)     = LeafPEPC_brch(NB,NZ)+MesophyllPEPC
              
              LeafC4ChlC_brch(NB,NZ)   = LeafC4ChlC_brch(NB,NZ)+MesophyllChlC              
              LeafRubiscoC_brch(NB,NZ) = LeafRubiscoC_brch(NB,NZ)+BundlSheathRubiscoC
              
              LeafC3ChlC_brch(NB,NZ)   = LeafC3ChlC_brch(NB,NZ)+BundlSheathChlC

            ELSE IF(iPlantPhotosynsType_pft(NZ).EQ.ic3_photo)then    
              MesophyllChlC     = LeafProteinC_node(K,NB,NZ)*LeafProtein2Chl_pft(NZ)/3.5_r8        ! 1 gC chl for 3.5 gC protein
              MesophyllRubiscoC = LeafProteinC_node(K,NB,NZ)*LeafRubisco2Protein_pft(NZ)
              
              VcMaxRubiscoRef_brch(NB,NZ) = VcMaxRubiscoRef_brch(NB,NZ)+VmaxSpecRubCarboxyRef_pft(NZ)*MesophyllRubiscoC
              VoMaxRubiscoRef_brch(NB,NZ) = VoMaxRubiscoRef_brch(NB,NZ)+VmaxRubOxyRef_pft(NZ)*MesophyllRubiscoC

              LeafRubiscoC_brch(NB,NZ) = LeafRubiscoC_brch(NB,NZ)+MesophyllRubiscoC
              
              LeafC3ChlC_brch(NB,NZ)   = LeafC3ChlC_brch(NB,NZ)+MesophyllChlC
    
            ENDIF
            ElectronTransptJmaxRef_brch(NB,NZ) = ElectronTransptJmaxRef_brch(NB,NZ)+SpecLeafChlAct_pft(NZ)*MesophyllChlC
          ENDIF
        ENDIF  
      ENDDO

      CanopyVcMaxRubisco25C_pft(NZ)  = CanopyVcMaxRubisco25C_pft(NZ)+VcMaxRubiscoRef_brch(NB,NZ)
      CanopyVoMaxRubisco25C_pft(NZ)  = CanopyVoMaxRubisco25C_pft(NZ)+VoMaxRubiscoRef_brch(NB,NZ)
      CanopyVcMaxPEP25C_pft(NZ)      = CanopyVcMaxPEP25C_pft(NZ) +VcMaxPEPCarboxyRef_brch(NB,NZ)
      ElectronTransptJmax25C_pft(NZ) = ElectronTransptJmax25C_pft(NZ)+ElectronTransptJmaxRef_brch(NB,NZ)

      LeafProteinCperm2LA_pft(NZ) = LeafProteinCperm2LA_pft(NZ)+LeafProteinC_brch(NB,NZ)
      LeafRubiscoCperm2LA_pft(NZ) = LeafRubiscoCperm2LA_pft(NZ)+LeafRubiscoC_brch(NB,NZ)
      LeafPEPCperm2LA_pft(NZ)     = LeafPEPCperm2LA_pft(NZ)+LeafPEPC_brch(NB,NZ)
      LeafC3ChlCperm2LA_pft(NZ)   = LeafC3ChlCperm2LA_pft(NZ)+LeafC3ChlC_brch(NB,NZ)
      LeafC4ChlCperm2LA_pft(NZ)   = LeafC4ChlCperm2LA_pft(NZ)+LeafC4ChlC_brch(NB,NZ)
    endif
  ENDDO

  if(tLeafArea>5.e-3_r8)then
    LeafProteinCperm2LA_pft(NZ) = LeafProteinCperm2LA_pft(NZ)/tLeafArea
    LeafRubiscoCperm2LA_pft(NZ) = LeafRubiscoCperm2LA_pft(NZ)/tLeafArea
    LeafPEPCperm2LA_pft(NZ)     = LeafPEPCperm2LA_pft(NZ)/tLeafArea
    LeafC3ChlCperm2LA_pft(NZ)   = LeafC3ChlCperm2LA_pft(NZ)/tLeafArea
    LeafC4ChlCperm2LA_pft(NZ)   = LeafC4ChlCperm2LA_pft(NZ)/tLeafArea

    CanopyVcMaxRubisco25C_pft(NZ)  = CanopyVcMaxRubisco25C_pft(NZ)/tLeafArea
    CanopyVoMaxRubisco25C_pft(NZ)  = CanopyVoMaxRubisco25C_pft(NZ)/tLeafArea
    CanopyVcMaxPEP25C_pft(NZ)      = CanopyVcMaxPEP25C_pft(NZ) /tLeafArea
    ElectronTransptJmax25C_pft(NZ) = ElectronTransptJmax25C_pft(NZ)/tLeafArea
    SpecificLeafArea_pft(NZ)       = tLeafArea/tLeafC
  else 
    LeafProteinCperm2LA_pft(NZ)    = 0._r8
    CanopyVcMaxRubisco25C_pft(NZ)  = 0._r8
    CanopyVoMaxRubisco25C_pft(NZ)  = 0._r8
    CanopyVcMaxPEP25C_pft(NZ)      = 0._r8
    ElectronTransptJmax25C_pft(NZ) = 0._r8
    LeafRubiscoCperm2LA_pft(NZ)    = 0._r8
    LeafPEPCperm2LA_pft(NZ)        = 0._r8
    LeafC3ChlCperm2LA_pft(NZ)      = 0._r8
    LeafC4ChlCperm2LA_pft(NZ)      = 0._r8
    SpecificLeafArea_pft(NZ)       = 0._r8
  endif
  call PrintInfo('end '//subname)
  end associate
  end subroutine PhotosynsDiag
!----------------------------------------------------------------------------------------------------
  subroutine StomatalDynamics(yearIJ,NZ,RawH2OCanopy2Atm_pft)
!
!     THIS subroutine CALCULATES CANOPY STOMATAL RESISTANCE AT MAXIMUM
!     CANOPY TURGOR FOR USE IN ENERGY BALANCE EQUATIONS IN 'UPTAKE'
!
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), optional, intent(in) :: RawH2OCanopy2Atm_pft

  character(len=*), parameter :: subname='StomatalDynamics'
  integer :: K,L,M,NB,N,I,J
  REAL(R8):: CanopyCO2BndlResist_pft
  real(r8), parameter :: RawMinH2OCanopy2Atm = 5.56E-03_r8 ![h/m]
  real(r8):: RI
!     begin_execution
  associate(                                                             &
    RIB                         => plt_ew%RIB                           ,& !input  :Richardson number for calculating boundary layer resistance, [-]
    RawIsoTCanopy2Atm_pft       => plt_ew%RawIsoTCanopy2Atm_pft         ,& !input  :canopy isothermal boundary later resistance, [h m-1]
    TairK                       => plt_ew%TairK                         ,& !input  :air temperature, [K]
    TKCanopy_pft                => plt_ew%TKCanopy_pft                  ,& !input  :canopy temperature, [K]
    CO2E                        => plt_site%CO2E                        ,& !input  :atmospheric CO2 concentration, [umol mol-1]
    CanopyLeafArea_pft          => plt_morph%CanopyLeafArea_pft         ,& !input  :plant canopy leaf area, [m2 d-2]
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft              ,& !input  :threshold zero for plang growth calculation, [-]
    NetCO2Flx2Canopy_col        => plt_bgcr%NetCO2Flx2Canopy_col        ,& !input  :total net canopy CO2 exchange, [g d-2 h-1]
    SineSunInclinationAngle_col => plt_rad%SineSunInclinationAngle_col  ,& !input  :sine of solar angle, [-]
    CanopyCi2CaRatio_pft        => plt_photo%CanopyCi2CaRatio_pft       ,& !input  :Ci:Ca ratio, [-]
    H2OCuticleResist_pft        => plt_photo%H2OCuticleResist_pft       ,& !input  :maximum stomatal resistance to vapor, [s h-1]
    CanopyGasCO2_pft            => plt_photo%CanopyGasCO2_pft           ,& !inoput :canopy gaesous CO2 concentration, [umol mol-1]
    AirConc_pft                 => plt_photo%AirConc_pft                ,& !output :total gas concentration, [mol m-3]
    CanopyMinStomaResistH2O_pft => plt_photo%CanopyMinStomaResistH2O_pft,& !output :canopy minimum stomatal resistance, [s m-1]
    LeafIntracellularCO2_pft    => plt_photo%LeafIntracellularCO2_pft    & !output :leaf gaseous CO2 concentration, [umol m-3]
  )
  call PrintInfo('beg '//subname)
  I=yearIJ%I;J=yearIJ%J
!
!     CANOPY TEMPERATURE + OFFSET FOR THERMAL ADAPTATION FROM 'READQ'
!
!     CANOPY BOUNDARY LAYER RESISTANCE
!
!     RI=Richardson's number
!     RIB=canopy isothermal Richardson's number
!     TairK,TKCanopy_pft=air,canopy temperature
!     RawIsoTCanopy2Atm_pft=canopy isothermal boundary later resistance
!     CanopyCO2BndlResist_pft=canopy boundary layer resistance to CO2, h/m
!     AirConc_pft=number of moles of air per m3
!
  if(present(RawH2OCanopy2Atm_pft))then
    CanopyCO2BndlResist_pft  = BndlResistCO2ToH2ORatio*RawH2OCanopy2Atm_pft
  else
    RI                       = RichardsonNumber(RIB,TairK,TKCanopy_pft(NZ))
    CanopyCO2BndlResist_pft  = BndlResistCO2ToH2ORatio &
      *AMAX1(RawMinH2OCanopy2Atm,RawIsoTCanopy2Atm_pft(NZ)/(1.0_r8-10.0_r8*RI))
  endif

  AirConc_pft(NZ)          = GetMolAirPerm3(TKCanopy_pft(NZ))    !assuming pressure is one atmosphere

  ! For prescribed phenolgoy, the canopy CO2 concentration should be computed differently
  !
  !     CANOPY CO2 CONCENTRATION FROM CO2 INFLUXES AND EFFLUXES 
  !
  !     CanopyGasCO2_pft,CO2E=CO2 concentrations in canopy air,atmosphere, umol mol-1 (ppmv)
  !     NetCO2Flx2Canopy_col=net CO2 flux in canopy air from soil,plants, g d-2 h-1, set to zero for prescribed phenolgoy
  !     assuming steady state, canopy CO2 concentration is computed with mass balance. 
  !     how 8.33E+04 is determined. 
  CanopyGasCO2_pft(NZ) = CO2E-8.33E+04_r8*NetCO2Flx2Canopy_col*CanopyCO2BndlResist_pft/AirConc_pft(NZ)
  CanopyGasCO2_pft(NZ) = AMIN1(CO2E+200.0_r8,AZMAX1(CO2E-200.0_r8,CanopyGasCO2_pft(NZ)))
  !
  !     MESOPHYLL CO2 CONCENTRATION FROM CI:CA RATIO ENTERED IN 'READQ' 
  !
  !     LeafIntracellularCO2_pft=intercellular CO2 concentration
  !     CanopyCi2CaRatio_pft=intercellular:atmospheric CO2 concn ratio from PFT file, parameter
  !     SineSunInclinationAngle_col=sine of solar angle
  !     CanopyLeafArea_pft=PFT leaf area 
  !
  LeafIntracellularCO2_pft(NZ)=CanopyCi2CaRatio_pft(NZ)*CanopyGasCO2_pft(NZ)

  IF(SineSunInclinationAngle_col.GT.0.0_r8 .AND. CanopyLeafArea_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
!
    call PhotoActivePFT(I,J,NZ)
  ELSE
!
    CanopyMinStomaResistH2O_pft(NZ)=H2OCuticleResist_pft(NZ)
  ENDIF
  call PrintInfo('end '//subname) 
  RETURN
  end associate
  END subroutine StomatalDynamics

!----------------------------------------------------------------------------------------------------
  subroutine C3FixCO2(LP,I,J,K,N,M,L,NB,NZ,PAR_zsec,Tau_rad,CH2O)
  implicit none
  integer, intent(in) :: LP,I,J,K,N,M,L,NB,NZ
  real(r8),intent(in) :: PAR_zsec,Tau_rad
  real(r8), intent(inout) :: CH2O
  real(r8) :: ETLF,EGRO,PARX,PARJ
  real(r8) :: VL
  character(len=*), parameter :: subname='C3FixCO2'
!     begin_execution
  associate(                                                                   &
    CO2lmtRubiscoCarboxyRate_node => plt_photo%CO2lmtRubiscoCarboxyRate_node  ,& !input  :carboxylation rate, [umol m-2 s-1]
    RubiscoCarboxyEff_node        => plt_photo%RubiscoCarboxyEff_node         ,& !input  :carboxylation efficiency, [umol umol-1]
    LigthSatCarboxyRate_node      => plt_photo%LigthSatCarboxyRate_node       ,& !input  :maximum light carboxylation rate under saturating CO2, [umol m-2 s-1]
    LeafEffArea_zsec           => plt_photo%LeafEffArea_zsec            ,& !input  :leaf irradiated surface area, [m2 d-2]
    RubiscoActivity_brch          => plt_photo%RubiscoActivity_brch            & !input  :branch down-regulation of CO2 fixation, [-]
  )
!
  call PrintInfo('beg '//subname)
!
!     QNTM=quantum efficiency
!     PAR=direct PAR flux
!     LigthSatCarboxyRate_node=light saturated e- transport rate
!     ETLF=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO=light-limited rubisco carboxylation rate
!
  PARX = QNTM*PAR_zsec
  PARJ = PARX+LigthSatCarboxyRate_node(K,NB,NZ)
  ETLF = (PARJ-SQRT(PARJ*PARJ-CURV4*PARX*LigthSatCarboxyRate_node(K,NB,NZ)))/CURV2
  EGRO = ETLF*RubiscoCarboxyEff_node(K,NB,NZ)
!
!     C3 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=rubisco carboxylation rate limited by light,CO2,N,P
!     CO2lmtRubiscoCarboxyRate_node=rubisco carboxylation rate limited by CO2
!     EGRO=light-limited rubisco carboxylation rate
!     RubiscoActivity_brch=N,P feedback inhibition on C3 CO2 fixation
!     CH2O=total rubisco carboxylation rate
!     LeafEffArea_zsec=unself-shaded leaf surface area
!     TAU_DirectSunLit=fraction of direct radiation transmitted from layer above
!
  VL   = AMIN1(CO2lmtRubiscoCarboxyRate_node(K,NB,NZ),EGRO)*RubiscoActivity_brch(NB,NZ)
  CH2O = CH2O+VL*LeafEffArea_zsec(N,L,K,NB,NZ)*TAU_Rad

  call PrintInfo('end '//subname)
  end associate
  end subroutine C3FixCO2

!----------------------------------------------------------------------------------------------------
  subroutine C3PhotosynsCanopyLayerL(I,J,L,K,NB,NZ,CH2O)
  implicit none
  integer, intent(in) :: I,J,L,K,NB,NZ
  real(r8), intent(inout) :: CH2O
  integer :: N,M,LP
  real(r8) :: PAR_zsec,Tau_rad
!     begin_execution
  associate(                                                &
    ZERO4Groth_pft       => plt_biom%ZERO4Groth_pft        ,& !input  :threshold zero for plang growth calculation, [-]
    LeafEffArea_zsec  => plt_photo%LeafEffArea_zsec  ,& !input  :leaf irradiated surface area, [m2 d-2]
    RadTotPARAbsorption_zsec       => plt_rad%RadTotPARAbsorption_zsec         ,& !input  :sunlit incoming PAR, [umol m-2 s-1]
    RadDifPARAbsorption_zsec       => plt_rad%RadDifPARAbsorption_zsec         ,& !input  :shade incoming PAR, [umol m-2 s-1]
    LeafAreaZsec_brch    => plt_morph%LeafAreaZsec_brch    ,& !input  :leaf surface area, [m2 d-2]    
    TAU_DirectSunLit  => plt_rad%TAU_DirectSunLit    ,& !input  :fraction of radiation intercepted by canopy layer, [-]
    TAU_DirectSunSha          => plt_rad%TAU_DirectSunSha             & !input  :fraction of radiation transmitted by canopy layer, [-]
  )
!     FOR EACH INCLINATION AND AZIMUTH CLASS
!
  DO N=1,NumLeafInclinationClasses1

    DO M=1,NumOfSkyAzimuthSects1
      IF(LeafEffArea_zsec(N,L,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
!
        DO LP=1,2
          IF(LP==1)THEN
            ! SUNLIT LEAVES
            PAR_zsec = RadTotPARAbsorption_zsec(N,M,L,NZ)
            Tau_rad  = TAU_DirectSunLit(L+1)
          else
            ! shade          
            PAR_zsec = RadDifPARAbsorption_zsec(N,M,L,NZ)
            Tau_rad  = TAU_DirectSunSha(L+1)
          ENDIF
    !
          if(PAR_zsec>0._r8)call C3FixCO2(LP,I,J,K,N,M,L,NB,NZ,PAR_zsec,Tau_rad,CH2O)     
     
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  end associate
  end subroutine C3PhotosynsCanopyLayerL

!----------------------------------------------------------------------------------------------------
  subroutine C3Photosynthesis(I,J,K,NB,NZ,CH2O,TFN_Carboxy,TFN_Oxygen,TFN_eTranspt,Km4RubOxy,ProteinCLeafAreaDensity)
  implicit none
  integer, intent(in) :: I,J,K,NB,NZ
  real(r8), intent(inout) :: CH2O
  real(r8), intent(in) :: TFN_Carboxy     !temperature sensitivity of carboyxlase
  real(r8), intent(in) :: TFN_Oxygen      !temperature sensitivity of oxygenase
  real(r8), intent(in) :: TFN_eTranspt    !temperature sensitivity of electron transport
  real(r8), intent(in) :: ProteinCLeafAreaDensity
  real(r8), intent(in) :: Km4RubOxy
  integer :: L
  real(r8) :: MesophyllChlDensity,MesophyllRubiscoSurfDensity
  real(r8) :: VOGRO,VoMaxRubiscoRef_node,VcMaxRubiscoRef_node,ElectronTransptJmaxRef_node
  character(len=*), parameter :: subname='C3Photosynthesis'
!     begin_execution
  associate(                                                                   &
    CanopyLeafArea_lnode          => plt_morph%CanopyLeafArea_lnode           ,& !input  :layer/node/branch leaf area, [m2 d-2]
    ZERO4Groth_pft                => plt_biom%ZERO4Groth_pft                  ,& !input  :threshold zero for plang growth calculation, [-]
    O2L_pft                       => plt_photo%O2L_pft                        ,& !input  :leaf aqueous O2 concentration, [uM]
    aquCO2Intraleaf_pft           => plt_photo%aquCO2Intraleaf_pft            ,& !input  :leaf aqueous CO2 concentration, [uM]
    LeafProtein2Chl_pft           => plt_photo%LeafProtein2Chl_pft            ,& !input  :fraction of leaf protein that is chlorophyll-binded, [gC gC-1]    
    Km4LeafaqCO2_pft              => plt_photo%Km4LeafaqCO2_pft               ,& !input  :leaf aqueous CO2 Km no O2, [uM]
    Km4RubiscoCarboxy_pft         => plt_photo%Km4RubiscoCarboxy_pft          ,& !input  :leaf aqueous CO2 Km ambient O2, [uM]
    SpecLeafChlAct_pft            => plt_photo%SpecLeafChlAct_pft             ,& !input  :cholorophyll activity at 25 oC, [umol g-1 h-1]
    VmaxSpecRubCarboxyRef_pft     => plt_photo%VmaxSpecRubCarboxyRef_pft      ,& !input  :rubisco carboxylase activity at 25 oC, [umol g-1 h-1 m-2]
    VmaxRubOxyRef_pft             => plt_photo%VmaxRubOxyRef_pft              ,& !input  :rubisco oxygenase activity at 25 oC, [umol g-1 h-1]
    LeafProteinC_node             => plt_biom%LeafProteinC_node               ,& !input  :node leaf protein C, [g d-2]    
    LeafArea_node                 => plt_morph%LeafArea_node                  ,& !input  :node leaf area, [m2 d-2]    
    LeafRubisco2Protein_pft       => plt_photo%LeafRubisco2Protein_pft        ,& !input  :leaf rubisco content, [gC gC-1]
    TAU_DirectSunLit           => plt_rad%TAU_DirectSunLit              ,& !input  :fraction of radiation intercepted by canopy layer, [-]
    TAU_DirectSunSha                   => plt_rad%TAU_DirectSunSha                      ,& !input  :fraction of radiation transmitted by canopy layer, [-]
    ElectronTransptJmaxRef_brch   => plt_photo%ElectronTransptJmaxRef_brch    ,& !inoput :Branch Maximum electron transport rate at reference temperature, [umol e- s-1]
    Vmax4RubiscoCarboxy_node      => plt_photo%Vmax4RubiscoCarboxy_node       ,& !output :maximum dark carboxylation rate under saturating CO2, [umol m-2 s-1]
    LigthSatCarboxyRate_node      => plt_photo%LigthSatCarboxyRate_node       ,& !output :maximum light carboxylation rate under saturating CO2, [umol m-2 s-1]
    CO2lmtRubiscoCarboxyRate_node => plt_photo%CO2lmtRubiscoCarboxyRate_node  ,& !output :carboxylation rate, [umol m-2 s-1]
    RubiscoCarboxyEff_node        => plt_photo%RubiscoCarboxyEff_node         ,& !output :carboxylation efficiency, [umol umol-1]
    CO2CompenPoint_node           => plt_photo%CO2CompenPoint_node             & !output :leaf node CO2 compensation point, [uM]
  )
  call PrintInfo('beg '//subname)
!
!     SURFICIAL DENSITY OF RUBISCO AND ITS LeafProtein2Chl_pftOROPHYLL
!
!     MesophyllRubiscoSurfDensity=surficial density of rubisco in mesophyll, [gC rubisco (gC protein)-1]*[gC protein m-2 leaf]
!     MesophyllChlDensity=surficial density of chlorophyll in esophyll,[gC chol (gC protein)-1]*[gC protein m-2 leaf]
!     LeafProtein2Chl_pft=mass ratio of mesophyll cholorophyll to leaf protein
!     ProteinCLeafAreaDensity=leaf protein surficial density, [g protein (m2 LA)-1]
!
  MesophyllRubiscoSurfDensity = LeafRubisco2Protein_pft(NZ)*ProteinCLeafAreaDensity
  MesophyllChlDensity         = LeafProtein2Chl_pft(NZ)*ProteinCLeafAreaDensity/3.5_r8          !1gC chl for 3.5 gC protein


!
!     CO2-LIMITED C3 CARBOXYLATION RATES
!
!     Vmax4RubiscoCarboxy_node=rubisco carboxylation rate unlimited by CO2
!     VmaxSpecRubCarboxyRef_pft=specific rubisco carboxylation activity from PFT file
!     TFN_Carboxy=temperature function for carboxylation
!     MesophyllRubiscoSurfDensity=surficial density of rubisco in mesophyll
!     VOGRO=rubisco oxygenation rate
!     TFN_Oxygen=temperature function for oxygenation
!     CO2CompenPoint_node=C3 CO2 compensation point (uM)
!     aquCO2Intraleaf_pft,O2L=intercellular CO2,O2 concentrations (uM)
!     Km4LeafaqCO2_pft,Km4RubiscoCarboxy_pft=Km for rubisco carboxylation without,with O2
!     Km4RubOxy=Km for rubisco oxygenation
!     CO2lmtRubiscoCarboxyRate_node=rubisco carboxylation rate limited by CO2
!
  VcMaxRubiscoRef_node                   = VmaxSpecRubCarboxyRef_pft(NZ)*MesophyllRubiscoSurfDensity
  VoMaxRubiscoRef_node                   = VmaxRubOxyRef_pft(NZ)*MesophyllRubiscoSurfDensity
  Vmax4RubiscoCarboxy_node(K,NB,NZ)      = VcMaxRubiscoRef_node*TFN_Carboxy
  VOGRO                                  = VoMaxRubiscoRef_node*TFN_Oxygen
  CO2CompenPoint_node(K,NB,NZ)           = 0.5_r8*O2L_pft(NZ)*VOGRO*Km4LeafaqCO2_pft(NZ)/(Vmax4RubiscoCarboxy_node(K,NB,NZ)*Km4RubOxy)
  ! Eq. (1) of Grant, 1989
  CO2lmtRubiscoCarboxyRate_node(K,NB,NZ) = AZMAX1(Vmax4RubiscoCarboxy_node(K,NB,NZ)&
    *(aquCO2Intraleaf_pft(NZ)-CO2CompenPoint_node(K,NB,NZ))/(aquCO2Intraleaf_pft(NZ)+Km4RubiscoCarboxy_pft(NZ)))
!
!     C3 ELECTRON TRANSFER RATES
!
!     LigthSatCarboxyRate_node=light-limited rubisco carboxylation rate, 
!     SpecLeafChlAct_pft=specific chlorophyll activity from PFT file
!     TFN_eTranspt=temperature function for e- transport
!     MesophyllChlDensity=surficial density of chlorophyll in mesophyll
!     RubiscoCarboxyEff_node=rubisco caboxylation efficiency
!     aquCO2Intraleaf_pft=intercellular CO2 concentrations (uM)
!     CO2CompenPoint_node=C3 CO2 compensation point (uM)
!     ELEC3=e- requirement for CO2 fixn by rubisco
!
  ElectronTransptJmaxRef_node          = SpecLeafChlAct_pft(NZ)*MesophyllChlDensity
  LigthSatCarboxyRate_node(K,NB,NZ)    = ElectronTransptJmaxRef_node*TFN_eTranspt
  RubiscoCarboxyEff_node(K,NB,NZ)      = AZMAX1((aquCO2Intraleaf_pft(NZ)-CO2CompenPoint_node(K,NB,NZ)) &
    /(ELEC3*aquCO2IntraLeaf_pft(NZ)+10.5_r8*CO2CompenPoint_node(K,NB,NZ)))

!
!     FOR EACH CANOPY LAYER
!
!     CanopyLeafArea_lnode=leaf area
!     LeafEffArea_zsec=unself-shaded leaf surface area
!

  DO L=NumCanopyLayers1,1,-1    
    IF(CanopyLeafArea_lnode(L,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      call C3PhotosynsCanopyLayerL(I,J,L,K,NB,NZ,CH2O)
    ENDIF
  ENDDO
  call PrintInfo('end '//subname)
  end associate

  end subroutine C3Photosynthesis

!----------------------------------------------------------------------------------------------------
  subroutine C4Photosynthesis(I,J,K,NB,NZ,CH2O,TFN_Carboxy,TFN_Oxygen,TFN_eTranspt,Km4RubOxy,ProteinCLeafAreaDensity)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: K       !leaf node id
  integer, intent(in) :: NB      !branch id
  integer, intent(in) :: NZ      !pft id
  real(r8), intent(inout) :: CH2O
  real(r8), intent(in) :: TFN_Carboxy        !temperature sensitivity of carboxylation
  real(r8), intent(in) :: TFN_Oxygen         !temperature sensitivity of oxygenation
  real(r8), intent(in) :: TFN_eTranspt       !temperature sensitivity of electron transport 
  real(r8), intent(in) :: Km4RubOxy
  real(r8), intent(in) :: ProteinCLeafAreaDensity
  integer :: L
  real(r8) :: CC4M
  real(r8) :: CCBS,MesophyllChlDensity
  real(r8) :: BundlSheathChlDensity
  real(r8) :: MesophyllPEPSurfDensity,BundlSheathRubiscoSurfDensity
  real(r8) :: ElectronTransptJmaxRef_node
  real(r8) :: VcMaxPEPCarboxyRef_node,VoMaxRubiscoRef_node,VcMaxRubiscoRef_node
  real(r8) :: VOGRO   !vmax 4 oxygenation in mesophyll
  character(len=*), parameter :: subname='C4Photosynthesis'
!     begin_execution
  associate(                                                                       &
    LeafElmntNode_brch              => plt_biom%LeafElmntNode_brch                ,& !input  :leaf element, [g d-2]
    ZERO4Groth_pft                  => plt_biom%ZERO4Groth_pft                    ,& !input  :threshold zero for plang growth calculation, [-]
    CanopyLeafArea_lnode            => plt_morph%CanopyLeafArea_lnode             ,& !input  :layer/node/branch leaf area, [m2 d-2]
    GrainFillDowreg_brch            => plt_photo%GrainFillDowreg_brch             ,& !input  :down-regulation of C4 photosynthesis, [-]
    fMesophyllChlProtein_pft        => plt_photo%fMesophyllChlProtein_pft         ,& !input  :fraction of Chl-bound protein in mesophyll cell, [gC gC-1]    
    O2L_pft                         => plt_photo%O2L_pft                          ,& !input  :leaf aqueous O2 concentration, [uM]
    aquCO2Intraleaf_pft             => plt_photo%aquCO2Intraleaf_pft              ,& !input  :leaf aqueous CO2 concentration, [uM]
    SpecLeafChlAct_pft              => plt_photo%SpecLeafChlAct_pft               ,& !input  :cholorophyll activity at 25 oC, [umol g-1 h-1]
    VmaxPEPCarboxyRef_pft           => plt_photo%VmaxPEPCarboxyRef_pft            ,& !input  :PEP carboxylase activity at 25 oC [umol g-1 h-1]
    LeafRubisco2Protein_pft         => plt_photo%LeafRubisco2Protein_pft          ,& !input  :leaf rubisco content, [gC gC-1]
    Km4PEPCarboxy_pft               => plt_photo%Km4PEPCarboxy_pft                ,& !input  :Km for PEP carboxylase activity, [uM]
    CMassCO2BundleSheath_node       => plt_photo%CMassCO2BundleSheath_node        ,& !input  :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CPOOL4_node                     => plt_photo%CPOOL4_node                      ,& !input  :leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
    Km4LeafaqCO2_pft                => plt_photo%Km4LeafaqCO2_pft                 ,& !input  :leaf aqueous CO2 Km no O2, [uM]
    Km4RubiscoCarboxy_pft           => plt_photo%Km4RubiscoCarboxy_pft            ,& !input  :leaf aqueous CO2 Km ambient O2, [uM]
    LeafProtein2Chl_pft             => plt_photo%LeafProtein2Chl_pft              ,& !input  :fraction of leaf protein that is chlorophyll-binded, [gC gC-1]    
    VmaxSpecRubCarboxyRef_pft       => plt_photo%VmaxSpecRubCarboxyRef_pft        ,& !input  :rubisco carboxylase activity at 25 oC, [umol g-1 h-1]
    VmaxRubOxyRef_pft               => plt_photo%VmaxRubOxyRef_pft                ,& !input  :rubisco oxygenase activity at 25 oC, [umol g-1 h-1]
    LeafPEP2Protein_pft             => plt_photo%LeafPEP2Protein_pft              ,& !input  :leaf PEP carboxylase content, [gC gC-1]
    LeafArea_node                   => plt_morph%LeafArea_node                    ,& !input  :leaf area on node, [m2 d-2]        
    NutrientCtrlonC4Carboxy_node    => plt_photo%NutrientCtrlonC4Carboxy_node     ,& !inoput :down-regulation of C4 photosynthesis, [-]
    ElectronTransptJmaxRef_brch     => plt_photo%ElectronTransptJmaxRef_brch      ,& !inoput :Branch Maximum electron transport rate at reference temperature, [umol e- s-1]
    RubiscoCarboxyEff_node          => plt_photo%RubiscoCarboxyEff_node           ,& !output :carboxylation efficiency, [umol umol-1]
    LigthSatCarboxyRate_node        => plt_photo%LigthSatCarboxyRate_node         ,& !output :maximum light carboxylation rate under saturating CO2, [umol m-2 s-1]
    CO2lmtRubiscoCarboxyRate_node   => plt_photo%CO2lmtRubiscoCarboxyRate_node    ,& !output :carboxylation rate, [umol m-2 s-1]
    Vmax4RubiscoCarboxy_node        => plt_photo%Vmax4RubiscoCarboxy_node         ,& !output :maximum dark carboxylation rate under saturating CO2, [umol m-2 s-1]
    CO2CompenPoint_node             => plt_photo%CO2CompenPoint_node              ,& !output :CO2 compensation point, [uM]
    Vmax4PEPCarboxy_node            => plt_photo%Vmax4PEPCarboxy_node             ,& !output :maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
    CO2lmtPEPCarboxyRate_node       => plt_photo%CO2lmtPEPCarboxyRate_node        ,& !output :C4 carboxylation rate, [umol m-2 s-1]
    C4CarboxyEff_node               => plt_photo%C4CarboxyEff_node                ,& !output :C4 carboxylation efficiency, [umol umol-1]
    LigthSatC4CarboxyRate_node      => plt_photo%LigthSatC4CarboxyRate_node        & !output :maximum light C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  )
  call PrintInfo('beg '//subname)
!
!     FEEDBACK ON C4 CARBOXYLATION FROM C4 NON-STRUCTURAL C
!
!     CC4M,CCBS=C4 nonstruct C concn in mesophyll,bundle sheath (uM)
!     CPOOL4_node,CMassCO2BundleSheath_node=C4 nonstructural C mass in mesophyll,bundle sheath
!     WGLF=leaf C mass
!     FWCBundlSheath,FWCMesophyll=leaf water content in bundle sheath, mesophyll
!     NutrientCtrlonC4Carboxy_node=N,P feedback inhibition on C4 CO2 fixation
!
  CC4M = AZMAX1(0.021E+09_r8*CPOOL4_node(K,NB,NZ)/(LeafElmntNode_brch(ielmc,K,NB,NZ)*FWCMesophyll))
  CCBS = AZMAX1(0.083E+09_r8*CMassCO2BundleSheath_node(K,NB,NZ)/(LeafElmntNode_brch(ielmc,K,NB,NZ)*FWCBundlSheath))
  
  NutrientCtrlonC4Carboxy_node(K,NB,NZ) = 1.0_r8/(1.0_r8+CC4M/C4KI_pepcarboxy)
  NutrientCtrlonC4Carboxy_node(K,NB,NZ) = NutrientCtrlonC4Carboxy_node(K,NB,NZ)*GrainFillDowreg_brch(NB,NZ)
!
!     SURFICIAL DENSITY OF LeafPEP2Protein_pftAND ITS LeafProtein2Chl_pftOROPHYLL
!
!     MesophyllPEPSurfDensity=surficial density of PEP carboxylase in mesophyll
!     MesophyllChlDensity=surficial density of chlorophyll in mesophyll
!     LeafPEP2Protein_pft=C mass ratio of PEP enzyme to total leaf protein
!     fMesophyllChlProtein_pft=fraction of leaf protein in mesophyll chlorophyll
!     ProteinCLeafAreaDensity=leaf protein surficial density, gC/m2 leaf
!
  MesophyllPEPSurfDensity = ProteinCLeafAreaDensity*LeafPEP2Protein_pft(NZ)
  MesophyllChlDensity     = ProteinCLeafAreaDensity*LeafProtein2Chl_pft(NZ)*fMesophyllChlProtein_pft(NZ)/3.7_r8   !1gC chl for 3.7 gC protein in mesophyll 

!
!     CO2-LIMITED C4 CARBOXYLATION RATES
!
!     Vmax4PEPCarboxy_node,CO2lmtPEPCarboxyRate_node=PEP carboxylation rate unlimited,limited by CO2
!     VmaxPEPCarboxyRef_pft=specific PEP carboxylase activity from PFT file
!     TFN_Carboxy=temperature function for carboxylation
!     MesophyllPEPSurfDensity=surficial density of PEP carboxylase in mesophyll
!     aquCO2Intraleaf_pft=intercellular CO2 concentrations (uM)
!     COMP4=C4 CO2 compensation point (uM)
!     Km4PEPCarboxy_pft=Km for VmaxPEPCarboxyRef_pft from PFT file (uM)
!
  VcMaxPEPCarboxyRef_node            = VmaxPEPCarboxyRef_pft(NZ)*MesophyllPEPSurfDensity
  Vmax4PEPCarboxy_node(K,NB,NZ)      = VcMaxPEPCarboxyRef_node*TFN_Carboxy
  CO2lmtPEPCarboxyRate_node(K,NB,NZ) = AZMAX1(Vmax4PEPCarboxy_node(K,NB,NZ) &
    *(aquCO2Intraleaf_pft(NZ)-COMP4)/(aquCO2Intraleaf_pft(NZ)+Km4PEPCarboxy_pft(NZ)))

!
!     C4 ELECTRON TRANSFER RATES
!
!     LigthSatC4CarboxyRate_node=light saturated e- transport rate
!     SpecLeafChlAct_pft=specific chlorophyll activity from PFT file
!     TFN_eTranspt=temperature function for e- transport
!     MesophyllChlDensity=surficial density of chlorophyll in mesophyll
!     C4CarboxyEff_node=PEP caboxylation efficiency
!     aquCO2Intraleaf_pft=intercellular CO2 concentrations (uM)
!     COMP4=C4 CO2 compensation point (uM)
!     ELEC4=e- requirement for CO2 fixn by PEP carboxylase
!
  ElectronTransptJmaxRef_node         = SpecLeafChlAct_pft(NZ)*MesophyllChlDensity  
  LigthSatC4CarboxyRate_node(K,NB,NZ) = ElectronTransptJmaxRef_node*TFN_eTranspt
  C4CarboxyEff_node(K,NB,NZ)          = AZMAX1((aquCO2Intraleaf_pft(NZ)-COMP4)/(ELEC4*aquCO2Intraleaf_pft(NZ)+10.5_r8*COMP4))


!
!     FOR EACH CANOPY LAYER
!
!     CanopyLeafArea_lnode=leaf area
!     LeafEffArea_zsec=unself-shaded leaf surface area
!
  DO L=NumCanopyLayers1,1,-1
    IF(CanopyLeafArea_lnode(L,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      call C4PhotosynsCanopyLayerL(I,J,L,K,NB,NZ,CH2O)
    ENDIF
  ENDDO
!
!     VARIABLES FOR C3 PHOTOSYNTHESIS DRIVEN BY C4
!
!     MesophyllRubiscoSurfDensity=surficial density of rubisco in bundle sheath, 
!     BundlSheathChlDensity=surficial density of chlorophyll in bundle sheath
!     LeafRubisco2Protein_pft=fraction of leaf protein in rubisco
!     LeafProtein2Chl_pft=fraction of leaf protein in bundle sheath chlorophyll
!     ProteinCLeafAreaDensity=leaf protein surficial density
!
  BundlSheathRubiscoSurfDensity = ProteinCLeafAreaDensity*LeafRubisco2Protein_pft(NZ)
  BundlSheathChlDensity         = ProteinCLeafAreaDensity*LeafProtein2Chl_pft(NZ)*(1._r8-fMesophyllChlProtein_pft(NZ))/5.5_r8  !1gC chl for 5.5 gC protein

!
!     CO2-LIMITED C3 CARBOXYLATION RATES
!
!     Vmax4RubiscoCarboxy_node=rubisco carboxylation rate unlimited by CO2
!     VmaxSpecRubCarboxyRef_pft=specific rubisco carboxylation activity from PFT file
!     TFN_Carboxy=temperature function for carboxylation
!     MesophyllRubiscoSurfDensity=surficial density of rubisco in bundle sheath
!     VOGRO=rubisco oxygenation rate
!     TFN_Oxygen=temperature function for oxygenation
!     CO2CompenPoint_node=C3 CO2 compensation point (uM)
!     aquCO2Intraleaf_pft,O2L=intercellular CO2,O2 concentrations (uM)
!     Km4LeafaqCO2_pft,Km4RubiscoCarboxy_pft=Km for rubisco carboxylation without,with O2
!     Km4RubOxy=Km for rubisco oxygenation
!     CO2lmtRubiscoCarboxyRate_node=rubisco carboxylation rate limited by CO2
!     CCBS=C4 nonstruct C concn in bundle sheath (uM)
!

  VcMaxRubiscoRef_node                   = VmaxSpecRubCarboxyRef_pft(NZ)*BundlSheathRubiscoSurfDensity
  Vmax4RubiscoCarboxy_node(K,NB,NZ)      = VcMaxRubiscoRef_node*TFN_Carboxy
  VoMaxRubiscoRef_node                   = VmaxRubOxyRef_pft(NZ)*BundlSheathRubiscoSurfDensity
  VOGRO                                  = VoMaxRubiscoRef_node*TFN_Oxygen
  CO2CompenPoint_node(K,NB,NZ)           = 0.5_r8*O2L_pft(NZ)*VOGRO*Km4LeafaqCO2_pft(NZ)/(Vmax4RubiscoCarboxy_node(K,NB,NZ)*Km4RubOxy)
  CO2lmtRubiscoCarboxyRate_node(K,NB,NZ) = AZMAX1(Vmax4RubiscoCarboxy_node(K,NB,NZ)* &
    (CCBS-CO2CompenPoint_node(K,NB,NZ))/(CCBS+Km4RubiscoCarboxy_pft(NZ)))
!
!     C3 ELECTRON TRANSFER RATES
!
!     LigthSatCarboxyRate_node=light-limited rubisco carboxylation rate
!     SpecLeafChlAct_pft=specific chlorophyll activity from PFT file
!     TFN_eTranspt=temperature function for e- transport
!     BundlSheathChlDensity=surficial density of chlorophyll in bundle sheath
!     RubiscoCarboxyEff_node=rubisco caboxylation efficiency
!     aquCO2Intraleaf_pft=intercellular CO2 concentrations (uM)
!     CO2CompenPoint_node=C3 CO2 compensation point (uM)
!     ELEC3=e- requirement for CO2 fixn by rubisco
!
  LigthSatCarboxyRate_node(K,NB,NZ) = SpecLeafChlAct_pft(NZ)*TFN_eTranspt*BundlSheathChlDensity
  RubiscoCarboxyEff_node(K,NB,NZ)   = AZMAX1((CCBS-CO2CompenPoint_node(K,NB,NZ))/(ELEC3*CCBS+10.5_r8*CO2CompenPoint_node(K,NB,NZ)))
  call PrintInfo('end '//subname)
  end associate
  end subroutine C4Photosynthesis

!----------------------------------------------------------------------------------------------------
  subroutine C4PhotosynsCanopyLayerL(I,J,L,K,NB,NZ,CH2O)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: L      !canopy layer id
  integer, intent(in) :: K      !leaf node id
  integer, intent(in) :: NB     !branch id
  integer, intent(in) :: NZ     !pft id
  real(r8), intent(inout) :: CH2O
  integer :: M,N,LP
  real(r8) :: PAR_zsec   !photosynthetically active radiation
  real(r8) :: TAU_rad    !PAR transmisivity
!     begin_execution
  associate(                                                &
    LeafEffArea_zsec  => plt_photo%LeafEffArea_zsec  ,& !input  :leaf irradiated surface area, [m2 d-2]
    ZERO4Groth_pft       => plt_biom%ZERO4Groth_pft        ,& !input  :threshold zero for plang growth calculation, [-]
    RadTotPARAbsorption_zsec       => plt_rad%RadTotPARAbsorption_zsec         ,& !input  :direct incoming PAR, [umol m-2 s-1]
    RadDifPARAbsorption_zsec       => plt_rad%RadDifPARAbsorption_zsec         ,& !input  :diffuse incoming PAR, [umol m-2 s-1]
    TAU_DirectSunSha          => plt_rad%TAU_DirectSunSha            ,& !input  :fraction of radiation transmitted by canopy layer, [-]
    TAU_DirectSunLit  => plt_rad%TAU_DirectSunLit     & !input  :fraction of radiation intercepted by canopy layer, [-]
  )
!
!     FOR EACH INCLINATION AND AZIMUTH CLASS
!
  DO  N=1,NumLeafInclinationClasses1
    DO  M=1,NumOfSkyAzimuthSects1
      IF(LeafEffArea_zsec(N,L,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
        DO LP=1,2
          ! SUNLIT LEAVES
          IF(LP==1)THEN
            PAR_zsec = RadTotPARAbsorption_zsec(N,M,L,NZ)
            TAU_rad  = TAU_DirectSunLit(L+1)
          ELSE
            ! SHADED LEAVES
            PAR_zsec = RadDifPARAbsorption_zsec(N,M,L,NZ)
            TAU_rad  = TAU_DirectSunSha(L+1)
          ENDIF
          if(PAR_zsec>0._r8)call C4FixCO2(I,J,K,N,M,L,NB,NZ,PAR_zsec,TAU_rad,CH2O)
        ENDDO  
      ENDIF
    ENDDO
  ENDDO
  end associate
  end subroutine C4PhotosynsCanopyLayerL

!----------------------------------------------------------------------------------------------------
  subroutine C4FixCO2(I,J,K,N,M,L,NB,NZ,PAR_zsec,TAU_rad,CH2O)

  implicit none
  integer, intent(in) :: I,J,K,N,M,L,NB,NZ
  real(r8),INTENT(IN) :: PAR_zsec,TAU_rad  
  real(r8), intent(inout) :: CH2O  
  real(r8) :: ETLF4,EGRO4,PARX,PARJ
  real(r8) :: VL
!     begin_execution
  associate(                                                                 &
    LeafEffArea_zsec          => plt_photo%LeafEffArea_zsec           ,& !input  :leaf irradiated surface area, [m2 d-2]
    NutrientCtrlonC4Carboxy_node => plt_photo%NutrientCtrlonC4Carboxy_node  ,& !input  :down-regulation of C4 photosynthesis, [-]
    C4CarboxyEff_node            => plt_photo%C4CarboxyEff_node             ,& !input  :C4 carboxylation efficiency, [umol umol-1]
    CO2lmtPEPCarboxyRate_node    => plt_photo%CO2lmtPEPCarboxyRate_node     ,& !input  :C4 carboxylation rate, [umol m-2 s-1]
    LigthSatC4CarboxyRate_node   => plt_photo%LigthSatC4CarboxyRate_node     & !input  :maximum light C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
  )
!
!     LIGHT-LIMITED CARBOXYLATION RATES
!
!     QNTM=quantum efficiency
!     RadDifPARAbsorption_zsec=diffuse PAR flux
!     LigthSatC4CarboxyRate_node=light saturated e- transport rate
!     ETLF4=light-limited e- transport rate
!     CURV=shape parameter for e- transport response to PAR
!     EGRO4=light-limited PEP carboxylation rate
!
  PARX  = QNTM*PAR_zsec
  PARJ  = PARX+LigthSatC4CarboxyRate_node(K,NB,NZ)
  ETLF4 = (PARJ-SQRT(PARJ*PARJ-CURV4*PARX*LigthSatC4CarboxyRate_node(K,NB,NZ)))/CURV2
  EGRO4 = ETLF4*C4CarboxyEff_node(K,NB,NZ)
!
!     C4 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
!
!     VL=PEP carboxylation rate limited by light,CO2,N,P
!     CO2lmtPEPCarboxyRate_node=PEP carboxylation rate limited by CO2
!     EGRO4=light-limited PEP carboxylation rate
!     CH2O=total PEP carboxylation rate
!     NutrientCtrlonC4Carboxy_node=N,P feedback inhibition on C4 CO2 fixation
!     LeafEffArea_zsec=unself-shaded leaf surface area
!     TAU_DirectSunSha=fraction of diffuse radiation transmitted from layer above
!  
  VL   = AMIN1(CO2lmtPEPCarboxyRate_node(K,NB,NZ),EGRO4)*NutrientCtrlonC4Carboxy_node(K,NB,NZ)
  CH2O = CH2O+VL*LeafEffArea_zsec(N,L,K,NB,NZ)*TAU_Rad
  end associate
  end subroutine C4FixCO2

!----------------------------------------------------------------------------------------------------
  subroutine LiveBranchPhotosynthesis(I,J,NB,NZ,CH2O,TFN_Carboxy,TFN_Oxygen,TFN_eTranspt,Km4RubOxy)
  implicit none
  integer, intent(in):: I,J,NB,NZ
  real(r8), intent(in) :: TFN_Carboxy
  real(r8), intent(in) :: TFN_Oxygen
  real(r8), intent(in) :: TFN_eTranspt
  real(r8), intent(in) :: Km4RubOxy  
  real(r8),intent(inout) :: CH2O

  character(len=*), parameter :: subname='LiveBranchPhotosynthesis'
  integer :: K  !leaf node id
  real(r8) :: ProteinCLeafAreaDensity   !protein area density, [gC protein m-2 LA]
!     begin_execution
  associate(                                                              &
    iPlantPhotosynsType_pft    => plt_photo%iPlantPhotosynsType_pft    ,& !input  :plant photosynthetic type (C3 or C4),[-]
    ZERO                        => plt_site%ZERO                         ,& !input  :threshold zero for numerical stability, [-]
    LeafElmntNode_brch          => plt_biom%LeafElmntNode_brch           ,& !input  :leaf element, [g d-2]
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft               ,& !input  :threshold zero for plang growth calculation, [-]
    LeafProteinC_node           => plt_biom%LeafProteinC_node            ,& !input  :layer leaf protein C, [g d-2]
    LeafArea_node               => plt_morph%LeafArea_node               ,& !input  :leaf area, [m2 d-2]
    Vmax4PEPCarboxy_node        => plt_photo%Vmax4PEPCarboxy_node        ,& !output :maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
    Vmax4RubiscoCarboxy_node    => plt_photo%Vmax4RubiscoCarboxy_node    ,& !output :maximum dark carboxylation rate under saturating CO2, [umol m-2 s-1]
    LeafProteinC_brch           => plt_biom%LeafProteinC_brch             & !output :Protein C for the branches, [gC protein d-2]
  )
  
  call PrintInfo('beg '//subname)
  DO K=1,MaxNodesPerBranch1
    IF(LeafArea_node(K,NB,NZ).GT.ZERO4Groth_pft(NZ) .AND. LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      ProteinCLeafAreaDensity   = LeafProteinC_node(K,NB,NZ)/LeafArea_node(K,NB,NZ)      
    ELSE
      ProteinCLeafAreaDensity=0.0_r8
    ENDIF
    
    IF(ProteinCLeafAreaDensity.GT.ZERO)THEN
      !
      !     iPlantPhotosynsType_pft=photosynthesis type:3=C3,4=C4 from PFT file
      !
      IF(iPlantPhotosynsType_pft(NZ).EQ.ic4_photo)THEN        
        call C4Photosynthesis(I,J,K,NB,NZ,CH2O,TFN_Carboxy,TFN_Oxygen,TFN_eTranspt,Km4RubOxy,ProteinCLeafAreaDensity)
      ELSE IF(iPlantPhotosynsType_pft(NZ).EQ.ic3_photo)then        

        call C3Photosynthesis(I,J,K,NB,NZ,CH2O,TFN_Carboxy,TFN_Oxygen,TFN_eTranspt,Km4RubOxy,ProteinCLeafAreaDensity)
      ENDIF
     !
    ELSE
      Vmax4PEPCarboxy_node(K,NB,NZ)     = 0.0_r8
      Vmax4RubiscoCarboxy_node(K,NB,NZ) = 0.0_r8
    ENDIF
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine LiveBranchPhotosynthesis

!----------------------------------------------------------------------------------------------------
  subroutine PhenoActiveBranch(I,J,NB,NZ,CH2O,TFN_Carboxy,TFN_Oxygen,TFN_eTranspt,Km4RubOxy)
  implicit none
  integer , intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: TFN_Carboxy
  real(r8), intent(in) :: TFN_Oxygen
  real(r8), intent(in) :: TFN_eTranspt
  real(r8), intent(in) :: Km4RubOxy
  real(r8), intent(inout) :: CH2O  

  character(len=*), parameter :: subname='PhenoActiveBranch'
  integer :: NE
  real(r8) :: CNS,CPS
!     begin_execution
  associate(                                                           &
    LeafPetoNonstElmConc_brch => plt_biom%LeafPetoNonstElmConc_brch   ,& !input  :branch nonstructural C concentration, [g d-2]
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft       ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft  ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    HourFailGrainFill_brch    => plt_pheno%HourFailGrainFill_brch     ,& !input  :flag to detect physiological maturity from grain fill, [-]
    Hours2LeafOut_brch        => plt_pheno%Hours2LeafOut_brch         ,& !input  :counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
    iPlantPhenolPattern_pft   => plt_pheno%iPlantPhenolPattern_pft    ,& !input  :plant growth habit: annual or perennial,[-]
    isPlantBranchAlive_brch    => plt_pheno%isPlantBranchAlive_brch     ,& !input  :flag to detect branch death, [-]
    ZERO                      => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]
    RubiscoActivity_brch      => plt_photo%RubiscoActivity_brch       ,& !inoput :branch down-regulation of CO2 fixation, [-]
    GrainFillDowreg_brch      => plt_photo%GrainFillDowreg_brch        & !output :grain fill down-regulation of annual plants, [-]
  )

  call PrintInfo('beg '//subname)
!
!     FEEDBACK ON C3 CARBOXYLATION FROM NON-STRUCTURAL C:N:P
!
!     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
!     RubiscoActivity_brch=N,P feedback inhibition on C3 CO2 fixation
!     CNKI,CPKI=nonstructural N,P inhibition constant on rubisco
!
  IF(LeafPetoNonstElmConc_brch(ielmc,NB,NZ).GT.ZERO)THEN
    CNS=LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmn,NB,NZ)&
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CNKI_rubisco)
    CPS=LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmp,NB,NZ)&
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CPKI_rubisco)
    RubiscoActivity_brch(NB,NZ)=AMIN1(CNS,CPS)
  ELSE
    RubiscoActivity_brch(NB,NZ)=1.0_r8
  ENDIF
!
!     CHILLING
!
!     CHILL=accumulated chilling hours used to limit CO2 fixn
!
!     RubiscoActivity_brch(NB,NZ)=RubiscoActivity_brch(NB,NZ)/(1.0_r8+0.25_r8*ChillHours_pft(NZ))
!
!     DEHARDENING OF EVERGREENS IN SPRING
!
!     ATRP=hours above threshold temperature for dehardening since leafout
!     Hours4ConiferSpringDeharden=hours to full dehardening of conifers in spring
! deciduous
  IF(iPlantPhenolType_pft(NZ).NE.iphenotyp_evgreen .AND. iPlantTurnoverPattern_pft(NZ).GE.2)THEN
    !conifer modification
    RubiscoActivity_brch(NB,NZ)=RubiscoActivity_brch(NB,NZ)*AZMAX1(AMIN1(1.0_r8,Hours2LeafOut_brch(NB,NZ)/(0.9_r8*Hours4ConiferSpringDeharden)))
  ENDIF
!
!     TERMINATION OF ANNUALS
!
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     HourFailGrainFill_brch=number of hours with no grain fill after start of grain fill
!     Hours2KillAnuals=number of hours with no grain fill to terminate annuals
!
  IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. HourFailGrainFill_brch(NB,NZ).GT.0.0_r8)THEN
    GrainFillDowreg_brch(NB,NZ)=AZMAX1(1.0_r8-HourFailGrainFill_brch(NB,NZ)/Hours2KillAnuals(iPlantPhenolType_pft(NZ)))
  ELSE
    GrainFillDowreg_brch(NB,NZ)=1.0_r8
  ENDIF
  RubiscoActivity_brch(NB,NZ)=RubiscoActivity_brch(NB,NZ)*GrainFillDowreg_brch(NB,NZ)
!
!     FOR EACH NODE
!
!     LeafArea_node,WGLF,LeafProteinC_node=leaf area,C mass,protein mass
!     ProteinCLeafAreaDensity=leaf protein C surficial density (gC protein m-2 leaf)
!

  IF(isPlantBranchAlive_brch(NB,NZ).EQ.iTrue)THEN
    call LiveBranchPhotosynthesis(I,J,NB,NZ,CH2O,TFN_Carboxy,TFN_Oxygen,TFN_eTranspt,Km4RubOxy)
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine PhenoActiveBranch

!----------------------------------------------------------------------------------------------------
  subroutine PrepPhotosynthesis(I,J,NZ,CH2O,TFN_Carboxy,TFN_Oxygen,TFN_eTranspt,Km4RubOxy)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(out) :: CH2O                !carboxylation rate
  real(r8), intent(out) :: TFN_Carboxy
  real(r8), intent(out) :: TFN_Oxygen
  real(r8), intent(out) :: TFN_eTranspt
  real(r8), intent(out) :: Km4RubOxy
  real(r8) :: ACTV,RTK
  real(r8) :: STK,TCCZ
  real(r8) :: TKCO
!     begin_execution
  associate(                                                           &
    TKCanopy_pft              => plt_ew%TKCanopy_pft                  ,& !input  :canopy temperature, [K]
    TempOffset_pft            => plt_pheno%TempOffset_pft             ,& !input  :adjustment of Arhhenius curves for plant thermal acclimation, [oC]
    XKO2_pft                  => plt_photo%XKO2_pft                   ,& !input  :Km for rubisco oxygenase activity, [uM]
    CanopyGasCO2_pft          => plt_photo%CanopyGasCO2_pft           ,& !input  :canopy gaesous CO2 concentration, [umol mol-1]
    AirConc_pft               => plt_photo%AirConc_pft                ,& !input  :total gas concentration, [mol m-3]
    LeafIntracellularCO2_pft  => plt_photo%LeafIntracellularCO2_pft   ,& !input  :leaf gaseous CO2 concentration, [umol m-3]
    O2I_pft                   => plt_photo%O2I_pft                    ,& !input  :leaf gaseous O2 concentration, [umol m-3]
    XKCO2_pft                 => plt_photo%XKCO2_pft                  ,& !input  :Km for rubisco carboxylase activity, [uM]
    CO2Solubility_pft         => plt_photo%CO2Solubility_pft          ,& !output :leaf CO2 solubility, [uM /umol mol-1]
    O2L_pft                   => plt_photo%O2L_pft                    ,& !output :leaf aqueous O2 concentration, [uM]
    DiffCO2Atmos2Intracel_pft => plt_photo%DiffCO2Atmos2Intracel_pft  ,& !output :gaesous CO2 concentration difference across stomates, [umol m-3]
    aquCO2Intraleaf_pft       => plt_photo%aquCO2Intraleaf_pft        ,& !output :leaf aqueous CO2 concentration, [uM]
    LeafO2Solubility_pft      => plt_photo%LeafO2Solubility_pft       ,& !output :leaf O2 solubility, [uM /umol mol-1]
    Km4RubiscoCarboxy_pft     => plt_photo%Km4RubiscoCarboxy_pft      ,& !output :leaf aqueous CO2 Km ambient O2, [uM]
    TFN_Carboxy_pft           => plt_photo%TFN_Carboxy_pft            ,& !output :temperature dependence of carboxylation, [-]
    TFN_Oxygen_pft            => plt_photo%TFN_Oxygen_pft             ,& !output :temperature dependence of oxygenation, [-]
    TFN_eTranspt_pft          => plt_photo%TFN_eTranspt_pft           ,& !output :temperature dependence of electron transport, [-]
    Km4LeafaqCO2_pft          => plt_photo%Km4LeafaqCO2_pft            & !output :leaf aqueous CO2 Km no O2, [uM]
  )
!
!     CO2 AND O2 AQUEOUS SOLUBILITY
!
!     TCCZ=canopy temperature
!     SCO2,LeafO2Solubility_pft=solubility of CO2,O2 (uM/(umol mol-1))
!     aquCO2Intraleaf_pft,O2L=intercellular CO2,O2 concentrations (uM)
!     DiffCO2Atmos2Intracel_pft=difference between atmosph and intercellular CO2 concn (umol m-3)
!
  TCCZ                          = TKCanopy_pft(NZ)-TC2K
  CO2Solubility_pft(NZ)         = EXP(-2.621_r8-0.0317_r8*TCCZ)
  LeafO2Solubility_pft(NZ)      = EXP(-6.175_r8-0.0211_r8*TCCZ)
  aquCO2Intraleaf_pft(NZ)       = LeafIntracellularCO2_pft(NZ)*CO2Solubility_pft(NZ)
  O2L_pft(NZ)                   = O2I_pft(NZ)*LeafO2Solubility_pft(NZ)
  DiffCO2Atmos2Intracel_pft(NZ) = AirConc_pft(NZ)*(CanopyGasCO2_pft(NZ)-LeafIntracellularCO2_pft(NZ))
!
!     ARRHENIUS FUNCTIONS FOR CARBOXYLATION AND OXYGENATION
!
!     TKCanopy_pft,TKCO=canopy temperature,canopy temp used in Arrhenius eqn
!     TempOffset_pft=shift in Arrhenius curve for thermal adaptation
!     TFN_Carboxy,TFN_Oxygen,TFN_eTranspt=temperature function for carboxylation,
!     oxygenation,e- transport (25 oC =1)
!     8.313,710.0=gas constant,enthalpy
!     197500,222500=energy of high,low temp inactivn(KJ mol-1)
!     65000,60000,43000=activation energy for carboxylation,
!     oxygenation,e- transport
!
  CH2O                 = 0.0_r8
  TKCO                 = TKCanopy_pft(NZ)+TempOffset_pft(NZ)
  RTK                  = RGASC*TKCO
  STK                  = 710.0_r8*TKCO
  ACTV                 = 1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  TFN_Carboxy          = EXP(26.237_r8-65000._r8/RTK)/ACTV
  TFN_Oxygen           = EXP(24.220_r8-60000._r8/RTK)/ACTV
  TFN_eTranspt         = EXP(17.362_r8-43000._r8/RTK)/ACTV
  TFN_Carboxy_pft(NZ)  = TFN_Carboxy
  TFN_Oxygen_pft(NZ)   = TFN_Oxygen
  TFN_eTranspt_pft(NZ) = TFN_eTranspt

!
!     M-M CONSTANT FOR CARBOXYLATION FROM 'READQ' ADJUSTED FOR TEMPERATURE
!
!     Km4LeafaqCO2_pft,Km4RubiscoCarboxy_pft=Km for rubisco carboxylation without,with O2
!     Km4RubOxy=Km for rubisco oxygenation
!
  Km4LeafaqCO2_pft(NZ)      = XKCO2_pft(NZ)*EXP(16.136_r8-40000._r8/RTK)
  Km4RubOxy                 = XKO2_pft(NZ)*EXP(8.067_r8-20000._r8/RTK)
  !Eq.(1) of Grant, 1989
  Km4RubiscoCarboxy_pft(NZ) = Km4LeafaqCO2_pft(NZ)*(1.0_r8+O2L_pft(NZ)/Km4RubOxy)

  end associate
  end subroutine PrepPhotosynthesis

!----------------------------------------------------------------------------------------------------
  subroutine PhotoActivePFT(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NZ

  integer :: NB,K
  real(r8) :: CH2O
  real(r8) :: RSX              !minimum canopy stomatal resistance to CO2
  real(r8) :: TFN_Carboxy,TFN_Oxygen,TFN_eTranspt,Km4RubOxy
  real(r8), parameter :: secsperhour=3600.0_r8
  character(len=*), parameter :: subname='PhotoActivePFT'

!     begin_execution
  associate(                                                               &
    HourReq4LeafOff_brch        => plt_pheno%HourReq4LeafOff_brch         ,& !input  :number of hours below set temperature required for autumn leafoff/hardening, [-]
    Hours4LeafOff_brch          => plt_pheno%Hours4LeafOff_brch           ,& !input  :cold requirement for autumn leafoff/hardening, [h]
    HourReq4LeafOut_brch        => plt_pheno%HourReq4LeafOut_brch         ,& !input  :hours above threshold temperature required for spring leafout/dehardening, [-]
    Hours4Leafout_brch          => plt_pheno%Hours4Leafout_brch           ,& !input  :heat requirement for spring leafout/dehardening, [h]
    iPlantPhenolType_pft        => plt_pheno%iPlantPhenolType_pft         ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft                ,& !input  :threshold zero for plang growth calculation, [-]
    NU                          => plt_site%NU                            ,& !input  :current soil surface layer number, [-]
    AREA3                       => plt_site%AREA3                         ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    DiffCO2Atmos2Intracel_pft   => plt_photo%DiffCO2Atmos2Intracel_pft    ,& !input  :gaesous CO2 concentration difference across stomates, [umol m-3]
    H2OCuticleResist_pft        => plt_photo%H2OCuticleResist_pft         ,& !input  :maximum stomatal resistance to vapor, [s h-1]
    FracPARads2LiveCanopy_pft   => plt_rad%FracPARads2LiveCanopy_pft      ,& !input  :fraction of incoming PAR absorbed by live canopy, [-]
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft            ,& !input  :number of branches,[-]
    Vmax4RubiscoCarboxy_node    => plt_photo%Vmax4RubiscoCarboxy_node     ,& !output :maximum dark carboxylation rate under saturating CO2, [umol m-2 s-1]
    RubiscoActivity_brch        => plt_photo%RubiscoActivity_brch         ,& !output :branch down-regulation of CO2 fixation, [-]
    GrainFillDowreg_brch       => plt_photo%GrainFillDowreg_brch        ,& !output :down-regulation of C4 photosynthesis, [-]
    Vmax4PEPCarboxy_node        => plt_photo%Vmax4PEPCarboxy_node         ,& !output :maximum dark C4 carboxylation rate under saturating CO2, [umol m-2 s-1]
    CanopyMinStomaResistH2O_pft => plt_photo%CanopyMinStomaResistH2O_pft   & !output :canopy minimum stomatal resistance, [s m-1]
  )
  call PrintInfo('beg '//subname)
  
  call PrepPhotosynthesis(I,J,NZ,CH2O,TFN_Carboxy,TFN_Oxygen,TFN_eTranspt,Km4RubOxy)
!
!     FOR EACH BRANCH
!
  DO NB=1,NumOfBranches_pft(NZ)
!
!     FEEDBACK ON CO2 FIXATION
!
!     iPlantPhenolType_pft=phenology type from PFT file
!     Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
!

    IF(iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen &
      .OR. Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ) &
      .OR. Hours4LeafOff_brch(NB,NZ).LT.HourReq4LeafOff_brch(NB,NZ))THEN
      
      !there are photosynthetically active leaves       
      call PhenoActiveBranch(I,J,NB,NZ,CH2O,TFN_Carboxy,TFN_Oxygen,TFN_eTranspt,Km4RubOxy)

    ELSE
      RubiscoActivity_brch(NB,NZ)  = 0.0_r8
      GrainFillDowreg_brch(NB,NZ) = 1.0_r8
      DO K=1,MaxNodesPerBranch1
        Vmax4PEPCarboxy_node(K,NB,NZ)     = 0.0_r8
        Vmax4RubiscoCarboxy_node(K,NB,NZ) = 0.0_r8
      ENDDO
    ENDIF
  ENDDO

!
!     MINIMUM CANOPY STOMATAL RESISTANCE FROM CO2 CONCENTRATION
!     DIFFERENCE DIVIDED BY TOTAL CO2 FIXATION
!
!     RSX,CanopyMinStomaResistH2O_pft=minimum canopy stomatal resistance to CO2,H2O (h m-1)
!     CH2O=total PEP(C4) or rubisco(C3) carboxylation rate
!     FracPARads2LiveCanopy_pft=fraction of radiation received by each live PFT canopy
!     DiffCO2Atmos2Intracel_pft=difference between atmosph and intercellular CO2 concn (umol m-3)
!     AREA=area of grid cell
!     RSMY=minimum stomatal resistance for CO2 uptake (h m-1)
! hourly time step
! determine minimum canopy stomatal resistance to CO2 based on CO2 gradient
  IF(CH2O.GT.ZERO4Groth_pft(NZ))THEN
    RSX=FracPARads2LiveCanopy_pft(NZ)*DiffCO2Atmos2Intracel_pft(NZ)*AREA3(NU)/(CH2O*secsperhour)
  ELSE
    RSX=H2OCuticleResist_pft(NZ)*1.56_r8
  ENDIF

  CanopyMinStomaResistH2O_pft(NZ)=AMIN1(H2OCuticleResist_pft(NZ),AMAX1(RSMY_stomaCO2,RSX/1.56_r8))

  call PrintInfo('end '//subname)
  end associate
  end subroutine PhotoActivePFT
  ![tail]
end module stomatesMod
