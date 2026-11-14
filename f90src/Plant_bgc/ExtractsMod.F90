module ExtractsMod
!!
! Description:
!     THIS SUBROUTINE AGGREGATES ALL SOIL-PLANT C,N,P EXCHANGES
!     FROM 'UPTAKE' AMD 'GROSUB' AND SENDS RESULTS TO 'REDIST'
!
  use data_kind_mod, only: r8 => DAT_KIND_R8,yearIJ_type
  use minimathmod,   only: AZMAX1
  use PlantBalMod,   only: SumPlantRootGas
  use EcoSIMCtrlMod, only: lverb, ldo_sp_mode
  use DebugToolMod,  only: PrintInfo
  use PlantSiteAPIData, only : plt_site
  use PlantRadiationAPIData, only : plt_rad
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantBiomassAPIData, only : plt_biom
  use PlantEnergyWaterAPIData, only : plt_ew
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantRootBGCAPIData, only : plt_rbgc
  use EcosimConst
  use PlantBGCPars
  use PlantAPICommonData
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  public :: extracts
  public :: ExtractRootWaterHeat,ExtractCanopyEnergyWater,ExtractCanopyInventory
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  SUBROUTINE extracts(YEARIJ,LiveExchange_pft,PassiveExchange_pft,HydraulicsApplied)
  ! Description:
  ! Aggregate accepted hourly exchanges even when management removed their living population.
  implicit none

  TYPE(yearIJ_type), intent(in) :: yearIJ
  logical, intent(in) :: PassiveExchange_pft(JP1) ! accepted passive exchange, recorded before management [-]
  logical, intent(in) :: LiveExchange_pft(JP1) ! PFT live exchanges accepted before management [-]
  logical, optional, intent(in) :: HydraulicsApplied ! water/energy fluxes already aggregated before soil transport [-]
  integer  :: I, J
  character(len=*), parameter :: subname='extracts'
  integer :: NZ,L
  logical :: IncludeHydraulics

  I=yearIJ%I;J=yearIJ%J
  call PrintInfo('beg '//subname)
  IncludeHydraulics=.true.
  if(present(HydraulicsApplied))IncludeHydraulics=.not.HydraulicsApplied
  if(.not.IncludeHydraulics)then
    ! Rebuild post-growth inventories, retaining the early hourly energy change.
    plt_ew%CanopyBiomWater_col=0._r8
    plt_ew%WatHeldOnCanopy_col=0._r8
    plt_ew%SnowOnCanopy_col=0._r8
    plt_ew%CanopyHeatStor_col=0._r8
  endif

  call ExtractLitterfall()

  call SumPlantRootGas(I,J)

  DO NZ=1,plt_site%NP
    IF(LiveExchange_pft(NZ))THEN
      if(IncludeHydraulics)call ExtractRootWaterHeat(NZ)
      if(.not.ldo_sp_mode)then
        call CalcTotalLeafArea(NZ)

        call TotalGasandSoluteUptake(I,J,NZ)
      endif
      call ExtractCanopyFluxes(I,J,NZ,IncludeHydraulics)
    elseif(PassiveExchange_pft(NZ))then
      ! Thermal/wet-surface exchange does not require living roots or population.
      if(IncludeHydraulics)call ExtractCanopyEnergyWater(NZ)
    ENDIF
  ENDDO

  ! Inventories survive inactive/retired slots independently of hourly flux eligibility.
  do NZ=1,plt_site%NP0
    call ExtractCanopyInventory(NZ)
  enddo

  call PrintInfo('end '//subname)  
  RETURN
  END subroutine extracts


!----------------------------------------------------------------------------------------------------
  subroutine ExtractLitterfall()
  ! Description:
  ! Aggregate litter, standing-dead stocks and cumulative element-balance diagnostics for all retained slots.
  implicit none
  character(len=*), parameter :: subname='ExtractLitterfall'
  integer :: NZ,L,K,M
  integer :: NE
  call PrintInfo('beg '//subname)
  associate(                                                          &
    LitrfallElms_pft          => plt_bgcr%LitrfallElms_pft           ,& !input  :plant element LitrFall, [g d-2 h-1]
    LitrfallElms_pvr          => plt_bgcr%LitrfallElms_pvr           ,& !input  :plant LitrFall element, [g d-2 h-1]
    MaxSoilLays4Root_pft      => plt_morph%MaxSoilLays4Root_pft      ,& !input  :maximum soil layer number for all root axes,[-]
    NP0                       => plt_site%NP0                        ,& !input  :intitial number of plant species,[-]
    StandDeadCompKElms_pft    => plt_biom%StandDeadCompKElms_pft     ,& !input: current standing-dead kinetic pools [g/cell]
    StandDeadStrutElms_pft    => plt_biom%StandDeadStrutElms_pft     ,& !input  :standing dead element, [g d-2]
    PlantElmBalCum_pft        => plt_site%PlantElmBalCum_pft         ,& !input: cumulative PFT C/N/P balance residual [g]
    PlantElemntStoreLandscape => plt_site%PlantElemntStoreLandscape  ,& !inout: sum of retained-PFT residuals [g]
    LitrFallStrutElms_col     => plt_bgcr%LitrFallStrutElms_col      ,& !inoput :total LitrFall structural element mass, [g d-2 h-1]
    LitrfalStrutElms_vr       => plt_bgcr%LitrfalStrutElms_vr        ,& !inoput :total LitrFall element, [g d-2 h-1]
    StandingDeadStrutElms_col => plt_biom%StandingDeadStrutElms_col  ,& !inoput :total standing dead biomass chemical element, [g d-2]
    CanopyLeafAareZ_col       => plt_morph%CanopyLeafAareZ_col       ,& !output :total leaf area, [m2 d-2]
    CanopyLeafArea_col        => plt_morph%CanopyLeafArea_col        ,& !output :grid canopy leaf area, [m2 d-2]
    CanopyStemAareZ_col       => plt_morph%CanopyStemAareZ_col       ,& !output :total stem area, [m2 d-2]
    StemArea_col              => plt_morph%StemArea_col              ,& !output :grid canopy stem area, [m2 d-2]
    tCanLeafC_clyr            => plt_biom%tCanLeafC_clyr              & !output :total leaf carbon mass in canopy layers, [gC d-2]
  )

  DO NZ=1,NP0
    ! Growth/management changes the component pools. Refresh their derived sum
    ! before inventory aggregation; hourly flux eligibility is already recorded.
    if(.not.ldo_sp_mode)then
      do NE=1,NumPlantChemElms
        StandDeadStrutElms_pft(NE,NZ)=sum(StandDeadCompKElms_pft(NE,:,NZ))
      enddo
    endif
    !
    !   TOTAL LitrFall OF ALL PLANT SPECIES
    !
    !   ZCSNC,ZZSNC,ZPSNC=total C,N,P LitrFall
    !   HCSNC,HZSNC,HPSNC=hourly PFT C,N,P LitrFall from grosub.f
    !   StandingDeadStrutElms_col=total standing dead C,N,P mass
    !   WTSTG=PFT standing dead C,N,P mass
    !   LitrfallElms_pvr,=cumulative PFT C,N,P LitrFall from grosub.f
    !   LitrfalStrutElms_vr,=cumulative total C,N,P LitrFall
    !
    DO NE=1,NumPlantChemElms
      LitrFallStrutElms_col(NE)     = LitrFallStrutElms_col(NE)+LitrfallElms_pft(NE,NZ)
      StandingDeadStrutElms_col(NE) = StandingDeadStrutElms_col(NE)+StandDeadStrutElms_pft(NE,NZ)
      if(.not.ldo_sp_mode)PlantElemntStoreLandscape(NE)=PlantElemntStoreLandscape(NE)+PlantElmBalCum_pft(NE,NZ)
    ENDDO

    DO  L=0,MaxSoilLays4Root_pft(NZ)
      DO K=1,pltpar%NumOfPlantLitrCmplxs
        DO NE=1,NumPlantChemElms
          DO  M=1,pltpar%jsken
            LitrfalStrutElms_vr(NE,M,K,L)=LitrfalStrutElms_vr(NE,M,K,L)+AZMAX1(LitrfallElms_pvr(NE,M,K,L,NZ))
          enddo
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CanopyLeafArea_col = 0._r8
  StemArea_col       = 0._r8

  DO L  = 1, NumCanopyLayers1
    CanopyLeafAareZ_col(L) = 0._r8
    tCanLeafC_clyr(L)      = 0._r8
    CanopyStemAareZ_col(L) = 0._r8
  ENDDO
  end associate
  call PrintInfo('end '//subname)
  end subroutine ExtractLitterfall

!----------------------------------------------------------------------------------------------------
  subroutine CalcTotalLeafArea(NZ)
!
!     TOTAL LEAF AREA OF ALL PLANT SPECIES
!
!     CanopyLeafAareZ_col,CanopyStemAareZ_col=total leaf,stalk area of combined canopy layer
!     CanopyLeafAreaZ_pft,CanopyStemSurfAreaZ_pft=PFT leaf,stalk area in canopy layer
!     tCanLeafC_clyr=total leaf C of combined canopy layer
!     CanopyLeafCLyr_pft=PFT leaf C in canopy layer
!
  implicit none
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='CalcTotalLeafArea'
  integer :: L
  associate(                                                       &
    CanopyLeafAreaZ_pft     => plt_morph%CanopyLeafAreaZ_pft      ,& !input  :canopy layer leaf area, [m2 d-2]
    CanopyLeafCLyr_pft      => plt_biom%CanopyLeafCLyr_pft        ,& !input  :canopy layer leaf C, [g d-2]
    CanopyStemSurfAreaZ_pft => plt_morph%CanopyStemSurfAreaZ_pft  ,& !input  :plant canopy layer stem area, [m2 d-2]
    CanopyLeafAareZ_col     => plt_morph%CanopyLeafAareZ_col      ,& !inoput :total leaf area, [m2 d-2]
    CanopyStemAareZ_col     => plt_morph%CanopyStemAareZ_col      ,& !inoput :total stem area, [m2 d-2]
    tCanLeafC_clyr          => plt_biom%tCanLeafC_clyr             & !inoput :total leaf carbon mass in canopy layers, [gC d-2]
  )
  call PrintInfo('beg '//subname)
  DO L=1,NumCanopyLayers1
    CanopyLeafAareZ_col(L) = CanopyLeafAareZ_col(L)+CanopyLeafAreaZ_pft(L,NZ)
    tCanLeafC_clyr(L)      = tCanLeafC_clyr(L)+CanopyLeafCLyr_pft(L,NZ)
    CanopyStemAareZ_col(L) = CanopyStemAareZ_col(L)+CanopyStemSurfAreaZ_pft(L,NZ)
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine CalcTotalLeafArea

!------------------------------------------------------------------------------------------
  subroutine ExtractRootWaterHeat(NZ)
  ! Description:
  ! Aggregate accepted root water and donor enthalpy independently of BGC.
  ! Standalone plant modes apply them before current soil transport; external
  ! soil coupling controls when it consumes the accepted transfers.
  integer, intent(in) :: NZ ! PFT whose live exchange was accepted [-]
  character(len=*), parameter :: subname='ExtractRootWaterHeat'
  integer :: N,L

  call PrintInfo('beg '//subname)
  associate(                                               &
    Flux              => plt_ew%RPlantRootH2OUptk_pvr    , & ! input: canopy-to-soil water [m3/cell/hour; root type, layer, PFT]
    Water             => plt_ew%TWaterPlantRoot2Soil_vr  , & ! inout: queued water [m3/cell/hour; layer]
    Heat              => plt_ew%THeatLossRoot2Soil_vr    , & ! inout: queued donor enthalpy [MJ/cell/hour; layer]
    CanopyTemperature => plt_ew%TKCanopy_pft             , & ! input: accepted canopy temperature [K; PFT]
    SoilTemperature   => plt_ew%TKS_vr                   , & ! input: soil donor temperature [K; layer]
    RootTypes         => plt_morph%Myco_pft              , & ! input: number of root/mycorrhizal types [-; PFT]
    TopLayer          => plt_site%NU                     , & ! input: current upper soil-layer index [-]
    RootBottom        => plt_site%MaxNumRootLays         , & ! input: deepest rooted layer index [-]
    SoilBottom        => plt_site%NK                       & ! input: deepest available soil-layer index [-]
  )

  ! Match the full-model aggregation's layer range and addition order.
  do L=TopLayer,MIN(RootBottom+1,SoilBottom)
    do N=1,RootTypes(NZ)
      Water(L)=Water(L)+Flux(N,L,NZ)
      if(Flux(N,L,NZ)>0._r8)then
        Heat(L)=Heat(L)+Flux(N,L,NZ)*cpw*CanopyTemperature(NZ)
      else
        Heat(L)=Heat(L)+Flux(N,L,NZ)*cpw*SoilTemperature(L)
      endif
    enddo
  enddo
  end associate
  call PrintInfo('end '//subname)
  end subroutine ExtractRootWaterHeat
!------------------------------------------------------------------------------------------
  subroutine TotalGasandSoluteUptake(I,J,NZ)
  ! Description:
  ! Aggregate root density, gases and nutrient exchange. Root water and heat
  ! are queued separately so prescribed phenology can use the same transfer path.

  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), parameter :: subname='TotalGasandSoluteUptake'
  integer :: N,L,K,idg,NE,ids
 
  call PrintInfo('beg '//subname)
  associate(                                                          &
    AREA3                     => plt_site%AREA3                      ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    Myco_pft                  => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    MaxNumRootLays            => plt_site%MaxNumRootLays             ,& !input  :maximum root layer number,[-]
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    PlantPopuLive_pft         => plt_site%PlantPopuLive_pft          ,& !input  :plant population, [d-2]
    RCO2Emis2Root_rpvr        => plt_rbgc%RCO2Emis2Root_rpvr         ,& !input  :aqueous CO2 flux from roots to root water, [g d-2 h-1]
    RootH1PO4DmndBand_pvr     => plt_rbgc%RootH1PO4DmndBand_pvr      ,& !input  :HPO4 demand in band by each root population, [g d-2 h-1]
    RootH1PO4DmndSoil_pvr     => plt_rbgc%RootH1PO4DmndSoil_pvr      ,& !input  :HPO4 demand in non-band by each root population, [g d-2 h-1]
    RootH2PO4DmndBand_pvr     => plt_rbgc%RootH2PO4DmndBand_pvr      ,& !input  :root uptake of H2PO4 band, [g d-2 h-1]
    RootH2PO4DmndSoil_pvr     => plt_rbgc%RootH2PO4DmndSoil_pvr      ,& !input  :root uptake of H2PO4 non-band, [g d-2 h-1]
    RootLenDensPerPlant_pvr   => plt_morph%RootLenDensPerPlant_pvr   ,& !input  :root layer length density, [m m-3]
    Soil2RootMycoExudE_pvr    => plt_rbgc%Soil2RootMycoExudE_pvr     ,& !input  :root uptake (+ve) - exudation (-ve) of DOE, [g d-2 h-1]
    RootNH4DmndBand_pvr       => plt_rbgc%RootNH4DmndBand_pvr        ,& !input  :root uptake of NO3 band unconstrained by NO3, [g d-2 h-1]
    RootNH4DmndSoil_pvr       => plt_rbgc%RootNH4DmndSoil_pvr        ,& !input  :root uptake of NH4 non-band unconstrained by NH4, [g d-2 h-1]
    RootNO3DmndBand_pvr       => plt_rbgc%RootNO3DmndBand_pvr        ,& !input  :root uptake of NO3 non-band unconstrained by NO3, [g d-2 h-1]
    RootNO3DmndSoil_pvr       => plt_rbgc%RootNO3DmndSoil_pvr        ,& !input  :root uptake of NH4 band unconstrained by NH4, [g d-2 h-1]
    RootNutUptake_pvr         => plt_rbgc%RootNutUptake_pvr          ,& !input  :root uptake of Nutrient band, [g d-2 h-1]
    RootO2Dmnd4Resp_pvr       => plt_rbgc%RootO2Dmnd4Resp_pvr        ,& !input  :root O2 demand from respiration, [g d-2 h-1]
    RootO2Uptk_pvr            => plt_rbgc%RootO2Uptk_pvr             ,& !input  :aqueous O2 flux from roots to root water, [g d-2 h-1]
    RootO2_TotSink_pvr        => plt_bgcr%RootO2_TotSink_pvr         ,& !input  :root O2 sink for autotrophic respiraiton, [gC d-2 h-1]
    RootUptkSoiSol_pvr        => plt_rbgc%RootUptkSoiSol_pvr         ,& !input  :aqueous CO2 flux from roots to soil water, [g d-2 h-1]
    trcg_air2root_flx_pvr     => plt_rbgc%trcg_air2root_flx_pvr      ,& !input  :gaseous tracer flux through roots, [g d-2 h-1]
    RootCO2Ar2RootX_rpvr      => plt_rbgc%RootCO2Ar2RootX_rpvr       ,& !input  :root/myco respiration released to root/myco, [gC d-2 h-1]    
    RootCO2Ar2RootX_pvr       => plt_rbgc%RootCO2Ar2RootX_pvr        ,& !input  :root respiration released to root, [gC d-2 h-1]        
    REcoH1PO4DmndBand_vr      => plt_bgcr%REcoH1PO4DmndBand_vr       ,& !inoput :HPO4 demand in band by all microbial, root, myco populations, [gP d-2 h-1]
    REcoH1PO4DmndSoil_vr      => plt_bgcr%REcoH1PO4DmndSoil_vr       ,& !inoput :HPO4 demand in non-band by all microbial, root, myco populations, [gP d-2 h-1]
    REcoH2PO4DmndBand_vr      => plt_bgcr%REcoH2PO4DmndBand_vr       ,& !inoput :total root + microbial PO4 uptake band, [gP d-2 h-1]
    REcoH2PO4DmndSoil_vr      => plt_bgcr%REcoH2PO4DmndSoil_vr       ,& !inoput :total root + microbial PO4 uptake non-band, [gP d-2 h-1]
    REcoNH4DmndBand_vr        => plt_bgcr%REcoNH4DmndBand_vr         ,& !inoput :total root + microbial NH4 uptake band, [gN d-2 h-1]
    REcoNH4DmndSoil_vr        => plt_bgcr%REcoNH4DmndSoil_vr         ,& !inoput :total root + microbial NH4 uptake non-band, [gN d-2 h-1]
    REcoNO3DmndBand_vr        => plt_bgcr%REcoNO3DmndBand_vr         ,& !inoput :total root + microbial NO3 uptake band, [gN d-2 h-1]
    REcoNO3DmndSoil_vr        => plt_bgcr%REcoNO3DmndSoil_vr         ,& !inoput :total root + microbial NO3 uptake non-band, [gN d-2 h-1]
    REcoO2DmndResp_vr         => plt_bgcr%REcoO2DmndResp_vr          ,& !inoput :total root + microbial O2 uptake, [g d-2 h-1]
    RUptkRootO2_vr            => plt_bgcr%RUptkRootO2_vr             ,& !inoput :total root internal O2 flux, [g d-2 h-1]
    RootCO2Emis2Root_vr       => plt_bgcr%RootCO2Emis2Root_vr        ,& !inoput :total root CO2 flux, [gC d-2 h-1]
    RootCO2Emis2Root_pvr      => plt_bgcr%RootCO2Emis2Root_pvr       ,& !inoput :root CO2 flux, [gC d-2 h-1]    
    RootO2_TotSink_vr         => plt_bgcr%RootO2_TotSink_vr          ,& !inoput :all root O2 sink for autotrophic respiraiton, [gC d-2 h-1]
    tRootMycoExud2Soil_vr     => plt_bgcr%tRootMycoExud2Soil_vr      ,& !inoput :total root element exchange, [g d-2 h-1]
    totRootLenDens_vr         => plt_morph%totRootLenDens_vr         ,& !inoput :total root length density, [m m-3]
    trcg_air2root_flx_vr      => plt_rbgc%trcg_air2root_flx_vr       ,& !inoput :total internal root gas flux, [gC d-2 h-1]    
    trcs_Soil2plant_uptake_pvr=> plt_rbgc%trcs_Soil2plant_uptake_pvr ,& !inoput :plant root-soil solute flux non-band, [g d-2 h-1]    
    trcs_Soil2plant_uptake_vr => plt_rbgc%trcs_Soil2plant_uptake_vr   & !inoput :total root-soil solute flux non-band, [g d-2 h-1]
  )
  DO L=NU,MIN(MaxNumRootLays+1,plt_site%NK)
    DO N=1,Myco_pft(NZ)  
      !
      !     TOTAL ROOT DENSITY
      !
      !     totRootLenDens_vr=total root length density
      !     RootLenDensPerPlant_pvr=PFT root length density per plant
      !     RPlantRootH2OUptk_pvr=total water uptake
      !     RPlantRootH2OUptk_pvr=PFT root water uptake
      !     THeatLossRoot2Soil_vr=total convective heat in root water uptake, 
      !     TKS=soil temperature
      !     PP=PFT population, this is dynamic, and can goes to zero
      !
      IF(N.EQ.ipltroot)THEN
        totRootLenDens_vr(L)=totRootLenDens_vr(L)+RootLenDensPerPlant_pvr(N,L,NZ)*PlantPopuLive_pft(NZ)/AREA3(L)
      ENDIF
      !
      !     TOTAL ROOT BOUNDARY GAS FLUXES
      !
      DO idg=idg_beg,idg_NH3
        trcg_air2root_flx_vr(idg,L)=trcg_air2root_flx_vr(idg,L)+trcg_air2root_flx_pvr(idg,N,L,NZ)
      ENDDO
      RootO2_TotSink_vr(L)       = RootO2_TotSink_vr(L) + RootO2_TotSink_pvr(N,L,NZ)  !total O2 uptake by roots to support root autotrophic respiraiton
      RootCO2Emis2Root_vr(L)     = RootCO2Emis2Root_vr(L)+RCO2Emis2Root_rpvr(N,L,NZ)
      RootCO2Emis2Root_pvr(L,NZ) = RootCO2Emis2Root_pvr(L,NZ)+RCO2Emis2Root_rpvr(N,L,NZ)
      RUptkRootO2_vr(L)          = RUptkRootO2_vr(L)+RootO2Uptk_pvr(N,L,NZ)

      !(>0) uptake from soil into roots
      DO idg=idg_beg,idg_end
        trcs_Soil2plant_uptake_vr(idg,L)     = trcs_Soil2plant_uptake_vr(idg,L)+RootUptkSoiSol_pvr(idg,N,L,NZ)
        trcs_Soil2plant_uptake_pvr(idg,L,NZ) = trcs_Soil2plant_uptake_pvr(idg,L,NZ)+RootUptkSoiSol_pvr(idg,N,L,NZ)
      ENDDO
      RootCO2Ar2RootX_pvr(L,NZ) = RootCO2Ar2RootX_pvr(L,NZ)+RootCO2Ar2RootX_rpvr(N,L,NZ)
 
      !Nutrients are sorted according to band|soil
      do ids=ids_NH4B,ids_nuts_end
        trcs_Soil2plant_uptake_vr(ids,L)=trcs_Soil2plant_uptake_vr(ids,L)+RootNutUptake_pvr(ids,N,L,NZ)
        trcs_Soil2plant_uptake_pvr(ids,L,NZ)=trcs_Soil2plant_uptake_pvr(ids,L,NZ)+RootNutUptake_pvr(ids,N,L,NZ)
      ENDDO
      !
      !     TOTAL ROOT C,N,P EXUDATION
      !
      !     tRootMycoExud2Soil_vr=total nonstructl C,N,P exchange
      !     RDFOMC,RDFOMN,RDFOMP=PFT nonstructl C,N,P exchange
      !
      DO K=1,jcplx
        DO NE=1,NumPlantChemElms
          tRootMycoExud2Soil_vr(NE,K,L)=tRootMycoExud2Soil_vr(NE,K,L)-Soil2RootMycoExudE_pvr(NE,N,K,L,NZ)
        ENDDO
      ENDDO
      !
      !     TOTAL ROOT O2, NH4, NO3, PO4 UPTAKE CONTRIBUTES TO
      !     TOTAL ROOT + MICROBIAL UPTAKE USED TO CALCULATE
      !     COMPETITION CONSTRAINTS
      !
      !
      REcoO2DmndResp_vr(L)    = REcoO2DmndResp_vr(L)+RootO2Dmnd4Resp_pvr(N,L,NZ)
      REcoNH4DmndSoil_vr(L)   = REcoNH4DmndSoil_vr(L)+RootNH4DmndSoil_pvr(N,L,NZ)
      REcoNO3DmndSoil_vr(L)   = REcoNO3DmndSoil_vr(L)+RootNO3DmndSoil_pvr(N,L,NZ)
      REcoH2PO4DmndSoil_vr(L) = REcoH2PO4DmndSoil_vr(L)+RootH2PO4DmndSoil_pvr(N,L,NZ)
      REcoH1PO4DmndSoil_vr(L) = REcoH1PO4DmndSoil_vr(L)+RootH1PO4DmndSoil_pvr(N,L,NZ)
      REcoNH4DmndBand_vr(L)   = REcoNH4DmndBand_vr(L)+RootNH4DmndBand_pvr(N,L,NZ)
      REcoNO3DmndBand_vr(L)   = REcoNO3DmndBand_vr(L)+RootNO3DmndBand_pvr(N,L,NZ)
      REcoH2PO4DmndBand_vr(L) = REcoH2PO4DmndBand_vr(L)+RootH2PO4DmndBand_pvr(N,L,NZ)
      REcoH1PO4DmndBand_vr(L) = REcoH1PO4DmndBand_vr(L)+RootH1PO4DmndBand_pvr(N,L,NZ)
    ENDDO
  ENDDO
  
  end associate
  call PrintInfo('end '//subname)
  end subroutine TotalGasandSoluteUptake

!----------------------------------------------------------------------------------------------------
  subroutine ExtractCanopyFluxes(I,J,NZ,IncludeEnergy)
  ! Description:
  ! Accumulate live-plant carbon/nutrient outputs and the shared canopy heat/water ledger.
  implicit none
  integer, intent(in) :: I,J,NZ
  logical, intent(in) :: IncludeEnergy ! aggregate physical fluxes only when not already exported [-]
  character(len=*), parameter :: subname='ExtractCanopyFluxes'
  integer :: NE,idg

  call PrintInfo('beg '//subname)
  associate(                                                       &
    CO2NetFix_pft            => plt_bgcr%CO2NetFix_pft,            &
    CanopyLeafArea_pft       => plt_morph%CanopyLeafArea_pft,      &
    CanopyStemSurfArea_pft   => plt_morph%CanopyStemSurfArea_pft,  &
    NH3Dep2Can_pft           => plt_bgcr%NH3Dep2Can_pft,           &
    PlantRootSoilElmNetX_pft => plt_rbgc%PlantRootSoilElmNetX_pft, &
    RootGasLossDisturb_pft   => plt_bgcr%RootGasLossDisturb_pft,   &
    CanopyLeafArea_col       => plt_morph%CanopyLeafArea_col,      &
    Canopy_NEE_col           => plt_bgcr%Canopy_NEE_col,           &
    LitrFallStrutElms_col    => plt_bgcr%LitrFallStrutElms_col,    &
    NH3Emis_CumYr_pft        => plt_bgcr%NH3Emis_CumYr_pft,        &
    StemArea_col             => plt_morph%StemArea_col,            &
    TRootGasLossDisturb_col  => plt_rbgc%TRootGasLossDisturb_col   &
    )
    if(IncludeEnergy)call ExtractCanopyEnergyWater(NZ)
    if(.not.ldo_sp_mode)then
      Canopy_NEE_col=Canopy_NEE_col+CO2NetFix_pft(NZ)
      CanopyLeafArea_col=CanopyLeafArea_col+CanopyLeafArea_pft(NZ)
      StemArea_col=StemArea_col+CanopyStemSurfArea_pft(NZ)
      DO NE=1,NumPlantChemElms
        LitrFallStrutElms_col(NE)=LitrFallStrutElms_col(NE)-PlantRootSoilElmNetX_pft(NE,NZ)
      ENDDO
      DO idg=idg_beg,idg_NH3
        TRootGasLossDisturb_col(idg)=TRootGasLossDisturb_col(idg)+RootGasLossDisturb_pft(idg,NZ)
      ENDDO
      NH3Emis_CumYr_pft(NZ)=NH3Emis_CumYr_pft(NZ)+NH3Dep2Can_pft(NZ)
    endif
  end associate
  call PrintInfo('end '//subname)
  end subroutine ExtractCanopyFluxes

  !------------------------------------------------------------------------------------------
  subroutine ExtractCanopyEnergyWater(NZ)
    ! Description:
    ! Accumulate accepted canopy heat/water fluxes independently of living-plant
    ! carbon and nutrient accounting. Inventories are aggregated separately.
    integer, intent(in) :: NZ ! PFT slot [-]
    character(len=*), parameter :: subname='ExtractCanopyEnergyWater'
    call PrintInfo('beg '//subname)
    associate( &
          Eco_NetRad_col => plt_rad%Eco_NetRad_col, &
          RadNet2Canopy_pft => plt_rad%RadNet2Canopy_pft, &
          Eco_Heat_Latent_col => plt_ew%Eco_Heat_Latent_col, &
          CanopyEvapTransLHeat_pft => plt_ew%CanopyEvapTransLHeat_pft, &
          Eco_Heat_Sens_col => plt_ew%Eco_Heat_Sens_col, &
          HeatXAir2PCan_pft => plt_ew%HeatXAir2PCan_pft, &
          Eco_Heat_GrndSurf_col => plt_ew%Eco_Heat_GrndSurf_col, &
          HeatStorCanopy_pft => plt_ew%HeatStorCanopy_pft, &
          ETCanopy_CumYr_pft => plt_ew%ETCanopy_CumYr_pft, &
          Transpiration_pft => plt_ew%Transpiration_pft, &
          VapXAir2Canopy_pft => plt_ew%VapXAir2Canopy_pft, &
          QVegET_col => plt_ew%QVegET_col, &
          VapXAir2Canopy_col => plt_ew%VapXAir2Canopy_col, &
          LWRadCanG => plt_ew%LWRadCanG, &
          LWRadCanopy_pft => plt_rad%LWRadCanopy_pft)
      Eco_NetRad_col         = Eco_NetRad_col+RadNet2Canopy_pft(NZ)
      Eco_Heat_Latent_col    = Eco_Heat_Latent_col+CanopyEvapTransLHeat_pft(NZ)
      Eco_Heat_Sens_col      = Eco_Heat_Sens_col+HeatXAir2PCan_pft(NZ)
      Eco_Heat_GrndSurf_col  = Eco_Heat_GrndSurf_col+HeatStorCanopy_pft(NZ)
      ETCanopy_CumYr_pft(NZ) = ETCanopy_CumYr_pft(NZ)+Transpiration_pft(NZ)+VapXAir2Canopy_pft(NZ)
      QVegET_col             = QVegET_col+Transpiration_pft(NZ)+VapXAir2Canopy_pft(NZ)
      VapXAir2Canopy_col     = VapXAir2Canopy_col+VapXAir2Canopy_pft(NZ)
      LWRadCanG              = LWRadCanG+LWRadCanopy_pft(NZ)

    end associate
    call PrintInfo('end '//subname)
  end subroutine ExtractCanopyEnergyWater

!------------------------------------------------------------------------------------------
  subroutine ExtractCanopyInventory(NZ)
  ! Description:
  ! Retain water, snow and sensible energy in every named slot, including PFTs
  ! removed after exchange. No old exchange rates are accumulated here.
  integer, intent(in) :: NZ ! retained PFT slot [-]
  character(len=*), parameter :: subname='ExtractCanopyInventory'
  real(r8) :: Energy,WetCapacity
  call PrintInfo('beg '//subname)
  associate( &
    Population     => plt_site%PlantPopuLive_pft,      &
    Dead           => plt_biom%StandDeadStrutElms_pft, &
    HeatRemoved    => plt_ew%CanopyHeatLoss2Dist_col,  &
    Tissue         => plt_ew%CanopyBiomWater_pft,      &
    Water          => plt_ew%WatHeldOnCanopy_pft,      &
    Snow           => plt_ew%SnowOnCanopy_pft,         &
    TissueTotal    => plt_ew%CanopyBiomWater_col,      &
    WaterTotal     => plt_ew%WatHeldOnCanopy_col,      &
    SnowTotal      => plt_ew%SnowOnCanopy_col,         &
    Capacity       => plt_ew%VHeatCapCanopy_pft,       &
    Temperature    => plt_ew%TKC_pft,                  &
    EnergyTotal    => plt_ew%CanopyHeatStor_col,       &
    EnergyChange   => plt_ew%HeatFlx2Canopy_col,       &
    PreviousEnergy => plt_ew%ENGYX_pft                 &
  )
  ! Standing dead can be removed after CutPlant rebuilt capacity. Export only
  ! the residual dry heat; retained liquid and snow keep their own capacity.
  if(Population(NZ)<=0._r8.and.Dead(ielmc,NZ)<=0._r8)then
    WetCapacity=cpw*(Tissue(NZ)+Water(NZ))+cps*Snow(NZ)
    if(Capacity(NZ)>WetCapacity)then
      HeatRemoved  = HeatRemoved+(Capacity(NZ)-WetCapacity)*Temperature(NZ)
      Capacity(NZ) = WetCapacity
    endif
  endif

  TissueTotal = TissueTotal+Tissue(NZ)
  WaterTotal  = WaterTotal+Water(NZ)
  SnowTotal   = SnowTotal+Snow(NZ)
  ! Unused cold-start slots can have a missing temperature and zero capacity.
  Energy=0._r8
  if(Capacity(NZ)>0._r8)Energy=Capacity(NZ)*Temperature(NZ)

  EnergyTotal        = EnergyTotal+Energy
  EnergyChange       = EnergyChange+Energy-PreviousEnergy(NZ)
  PreviousEnergy(NZ) = Energy
  end associate
  call PrintInfo('end '//subname)
  end subroutine ExtractCanopyInventory
![tail]
end module ExtractsMod
