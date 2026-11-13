module PlantAPI4Uptake
  ! Description:
  ! Transfer column and PFT state to/from the uptake workspace, preserving
  ! prescribed plant existence separately from snow-exposed exchange area.
  use EcoSIMCtrlMod, only: ldo_sp_mode
!
! interface to integrate the plant model
! for prescribed phenology
  use data_kind_mod,              only: r8 => DAT_KIND_R8
  use EcoSiMParDataMod,           only: micpar, pltpar
  use SoilPhysDataType,           only: SurfAlbedo_col
  use MiniMathMod,                only: AZMAX1
  use EcoSIMCtrlMod,              only: ats_cpl_mode
  use PlantSiteAPIData,           only: plt_site
  use PlantPhotosynthesisAPIData, only: plt_photo
  use PlantRadiationAPIData,      only: plt_rad
  use PlantMorphologyAPIData,     only: plt_morph
  use PlantPhenologyAPIData,      only: plt_pheno
  use PlantSoilChemistryAPIData,  only: plt_soilchem
  use PlantBiomassAPIData,        only: plt_biom
  use PlantEnergyWaterAPIData,    only: plt_ew
  use NumericalAuxMod
  use EcoSIMSolverPar
  use EcoSIMHistMod
  use SnowDataType
  use TracerIDMod
  use SurfLitterDataType
  use LandSurfDataType
  use SoilPropertyDataType
  use ChemTranspDataType
  use EcoSimSumDataType
  use SoilHeatDataType
  use SOMDataType
  use ClimForcDataType
  use EcoSIMCtrlDataType
  use GridDataType
  use RootDataType
  use SoilWaterDataType
  use CanopyDataType
  use PlantDataRateType
  use PlantTraitDataType
  use CanopyRadDataType
  use FlagDataType
  use EcosimBGCFluxType
  use FertilizerDataType
  use SoilBGCDataType
  use PlantMgmtDataType
  use PlantAPICommonData
  use PlantTraitsTransferMod, only: SendPlantTraits
  use PrescribedLeafBiochemMod, only: PreparePrescribedLeafBiochemistry
  use PlantBGCRatesAPIData, only: plt_bgcr
implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: PlantUptakeAPISend
  public :: PlantUPtakeAPIRecv

  contains
!------------------------------------------------------------------------------------------

  subroutine PlantUptakeAPISend(I,J,NY,NX)
  ! Description:
  ! Load the current column's plant, soil and atmospheric inputs into the uptake
  ! workspace, including its hourly canopy-exposure diagnostic.
  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer :: M,L,NZ,N,NB,K

  ! Current root transfers are accumulated after hydrology and applied next
  ! hour. Reset the SP workspace per column, including its unused layer slots.
  plt_ew%TWaterPlantRoot2Soil_vr=0._r8
  plt_ew%THeatLossRoot2Soil_vr=0._r8
  !sent variables also modified
  plt_site%PlantPopu_col           = PlantPopu_col(NY,NX)              !total plant population in the column
  plt_site%ZERO                    = ZERO                              !numerical threshold
  plt_site%ZERO2                   = ZERO2                             !numerical threshold
  plt_morph%LeafStalkAreaAll_col   = LeafStalkAreaAll_col(NY,NX)          !leaf+stalk area,  set as phenology input
  plt_ew%BulkFactor4Snow_col       = BulkFactor4Snow_col(NY,NX)
  plt_morph%CanopyLeafArea_col     = CanopyLeafArea_col(NY,NX)         !canopy leaf area, set as phenolgoy input
  plt_site%NL                      = NL_col(NY,NX)                     !lower node number of the total vertical number of soil layers used for water uptake
  plt_site%NP0                     = NP0_col(NY,NX)                    !total number of plants in the column
  plt_site%MaxNumRootLays          = NL_col(NY,NX)                     !assume roots exist in all layers
  plt_site%NP                      = NP_col(NY,NX)                     !total number of plants in the column, set to NP0_col for simplicity
  plt_site%NU                      = NU_col(NY,NX)                     !top soil layer indices
  plt_site%NK                      = NK_col(NY,NX)                     !total layers with roots,  can be set to NL_col
  plt_site%ALT                     = SurfElevation_col(NY,NX)                    !column altitude [m]
  plt_ew%TKSnow                    = TKSnow_snvr(1,NY,NX)              !surface snow temperature
  plt_ew%SnowDepth                 = SnowDepth_col(NY,NX)              !total snow depth
  plt_ew%canopy_exposed_col        = canopy_exposed_col(NY,NX)         !hourly exposed-canopy flag for this column
  plt_ew%RawIsoTAtm2CanopySinkZ_col   = RawIsoTAtm2CanopySinkZ_col(NY,NX)    !aboveground canopy resistance,  computed in the radiation code
  plt_ew%TairK                     = TairK_col(NY,NX)                  !air temperature
  plt_ew%RoughnessLength               = RoughnessLength_col(NY,NX)            !computed in the radiation coce
  plt_morph%CanopyHeight_col       = CanopyHeight_col(NY,NX)           !canopy height as prescribed input
  plt_rad%LWRadGrnd_col            = LWRadGrnd_col(NY,NX)              !computed in surface energy model
  plt_rad%LWRadSky_col             = LWRadSky_col(NY,NX)               !computed in surface energy model
  plt_ew%VPA                       = VPA_col(NY,NX)                    !vapor pressure as climate input
  plt_ew%RIB                       = RIB_col(NY,NX)                    !richardson stability number,  computed in surface radiation model
  plt_site%ZEROS2                  = ZEROS2(NY,NX)                     !numerical threshold
  plt_ew%Air_Heat_Sens_store_col   = Air_Heat_Sens_store_col(NY,NX)    !initialized to zero,  but is needed iterately to update heat source for surface energy balance calculation?
  plt_ew%Air_Heat_Latent_store_col = Air_Heat_Latent_store_col(NY,NX)  !initialized to zero,  but is needed iteratively for surface energy calculation
  plt_ew%EMS_modify_scalar_col     = EMS_modify_scalar_col(NY,NX)             !scaling factor for longwave radiation, should be 1 in the coupling, but needs to be in the API
  !plt_ew%ZERO4PlantDisplace_col    = ZERO4PlantDisplace_col(NY,NX)    
  plt_ew%ZeroPlaneDisplacem_col    = ZeroPlaneDisplacem_col(NY,NX)     !numerical threshold
  plt_rad%Eco_NetRad_col           = Eco_NetRad_col(NY,NX)             !net ecosystem radiation, updated through model iterations, reset to zero in hour1.F90
  plt_ew%Eco_Heat_Sens_col         = Eco_Heat_Sens_col(NY,NX)          !whole ecosystem sensible heat, updated iteratively, reset to zero in hour1.F90
  plt_ew%Eco_Heat_Latent_col       = Eco_Heat_Latent_col(NY,NX)        !whole ecosystem latent heat, updated iterately, reset to zero in hour1.F90
  plt_ew%Eco_Heat_GrndSurf_col     = Eco_Heat_GrndSurf_col(NY,NX)      !heat to ground surface, updated iterately, reset to zero in hour1.F90
  plt_ew%CanopyBiomWater_col       = CanopyBiomWater_col(NY,NX)              !canopy water content, updated iterately
  plt_ew%WatHeldOnCanopy_col       = WatHeldOnCanopy_col(NY,NX)        !water held on canopy, updated iterately
  plt_ew%SnowOnCanopy_col          = SnowOnCanopy_col(NY,NX)
  plt_ew%QVegET_col                = QVegET_col(NY,NX)                 !canopy evapotranspiration, reset to zero in hour1.F90
  plt_ew%VapXAir2Canopy_col        = VapXAir2Canopy_col(NY,NX)         !canopy evaporation, reset to zero in hour1.F90
  plt_ew%CanopyHeatStor_col        = CanopyHeatStor_col(NY,NX)         !canopy heat storage, reset to zero in hour1.F90
  plt_ew%HeatSensVapAir2Canopy_col        = HeatSensVapAir2Canopy_col(NY,NX)         !heat flux to canopy, reset to zero in hour1.F90
  plt_ew%HeatPrecAir2Canopy_col        = HeatPrecAir2Canopy_col(NY,NX)         !heat flux to canopy, reset to zero in hour1.F90
  plt_ew%HeatFlx2Canopy_col        = HeatFlx2Canopy_col(NY,NX)         !heat flux to canopy, reset to zero in hour1.F90
  plt_ew%LWRadCanG                 = LWRadCanG_col(NY,NX)              !longwave radiation onto groud, reset to zero in hour1.F90
  DO L=1,NK_col(NY,NX)
   plt_soilchem%VLSoilPoreMicP_vr(L)      = VLSoilPoreMicP_vr(L,NY,NX)             !total soil volume (associated with micropores) in layer L, which is soil volume when no macropore is considered
   plt_ew%ElvAdjstedSoilH2OPSIMPa_vr(L)   = ElvAdjstedSoilH2OPSIMPa_vr(L,NY,NX)    !elevation adjusted soil matric water potential, MPa
   plt_soilchem%SoilBulkDensity_vr(L)     = SoilBulkDensity_vr(L,NY,NX)            !soil bulk density
   plt_site%CumSoilThickness_vr(L)        = CumSoilThickness_vr(L,NY,NX)           !cumulative soil thickness starting from the soil surface (which is zero)
   plt_site%AREA3(L)                      = AREA_3D(3,L,NY,NX)                     !grid size
   plt_ew%TKS_vr(L)                       = TKS_vr(L,NY,NX)                        !soil temperature
   plt_soilchem%THETW_vr(L)               = THETW_vr(L,NY,NX)                      !relative saturation of soil micropores in layer L
   plt_soilchem%SoilWatAirDry_vr(L)       = SoilWatAirDry_vr(L,NY,NX)              !relative moisture content when considered as air-dry
   plt_soilchem%VLSoilMicP_vr(L)          = VLSoilMicP_vr(L,NY,NX)                 !volume of soil micropore,
   plt_soilchem%VLiceMicP_vr(L)           = VLiceMicP_vr(L,NY,NX)                  !volume of ice in soil layer
   plt_soilchem%VLWatMicP_vr(L)           = VLWatMicP_vr(L,NY,NX)                  !volume of water in soil layer
   plt_soilchem%VLMicP_vr(L)              = VLMicP_vr(L,NY,NX)                     !total micropore volume
   plt_site%FracSoiAsMicP_vr(L)           = FracSoiAsMicP_vr(L,NY,NX)              !fraction of soil as micropore
   plt_site%DLYR3(L)                      = DLYR_3D(3,L,NY,NX)                     !layer thickness
   plt_soilchem%HYCDMicP4RootUptake_vr(L) = HYCDMicP4RootUptake_vr(L,NY,NX)        !soil hydraulic conductivity for water uptake,
   plt_site%CumSoilThickMidL_vr(L)        = CumSoilThickMidL_vr(L,NY,NX)           !layer middle soil depth
   DO M=1,NPH
     plt_site%VLWatMicPM_vr(M,L)           = VLWatMicP_vr(L,NY,NX)                 !intermediate iteration soil moisture
   ENDDO

  ENDDO

  DO NZ=1,NP0_col(NY,NX)
    plt_ew%fSnowCanopy_pft(NZ)              =fSnowCanopy_pft(NZ,NY,NX)   
    Myco_pft(NZ,NY,NX)=1
    plt_ew%DeltaTKC_pft(NZ)                 = DeltaTKC_pft(NZ,NY,NX)
    
    if(ats_cpl_mode)then
      if(DATAPI(NZ,NY,NX).GT.ZERO4LeafVar_pft(NZ,NY,NX))then
          IsPlantActive_pft(NZ,NY,NX)=iTrue
        else
          IsPlantActive_pft(NZ,NY,NX)=iFalse
      endif
    elseif(ldo_sp_mode)then
      ! Burial disables atmospheric exchange, not the prescribed PFT itself.
      ! Keep it eligible for the zero-exchange fallback, which clears old fluxes
      ! while retaining tissue water and intercepted stores through snow cover.
      if(CanopyLeafArea_pft(NZ,NY,NX)+CanopyStemSurfArea_pft(NZ,NY,NX)>ZERO4LeafVar_pft(NZ,NY,NX))then
        IsPlantActive_pft(NZ,NY,NX)=iTrue
      else
        IsPlantActive_pft(NZ,NY,NX)=iFalse
      endif
    else
      if(LeafStalkAreaAct_pft(NZ,NY,NX).GT.ZERO4LeafVar_pft(NZ,NY,NX))then
        IsPlantActive_pft(NZ,NY,NX)=iTrue
      else
        IsPlantActive_pft(NZ,NY,NX)=iFalse
      endif
    endif
    plt_morph%StalkAxialResist_pft(NZ)         = StalkAxialResist_pft(NZ,NY,NX)
    plt_morph%RootSingleVesselRstaxial_pft(NZ) = RootSingleVesselRstaxial_pft(NZ,NY,NX)
    plt_morph%RootSingleVesselArea_pft(NZ)     = RootSingleVesselArea_pft(NZ,NY,NX)
    plt_ew%ETCanopy_CumYr_pft(NZ)              = ETCanopy_CumYr_pft(NZ,NY,NX)
    plt_ew%ENGYX_pft(NZ)                       = ENGYX_pft(NZ,NY,NX)
    plt_ew%VHeatCapCanopy_pft(NZ)                = VHeatCapCanopy_pft(NZ,NY,NX)
    plt_morph%NGTopRootLayer_pft(NZ)        = 1
    plt_morph%LeafStalkAreaAct_pft(NZ)      = LeafStalkAreaAct_pft(NZ,NY,NX)      !pft leaf+stalk area
    plt_morph%StandDeadSurfAreaAct_pft(NZ)  = StandDeadSurfAreaAct_pft(NZ,NY,NX)
    plt_morph%StandDeadSurfArea_pft(NZ)     = StandDeadSurfArea_pft(NZ,NY,NX)
    plt_biom%StandDeadStrutElms_pft(:,NZ)   = StandDeadStrutElms_pft(:,NZ,NY,NX)
    plt_morph%HypocotHeight_pft(NZ)         = HypocotHeight_pft(NZ,NY,NX)      !(approximately) plant height
    plt_morph%NumOfBranches_pft(NZ)         = NumOfBranches_pft(NZ,NY,NX)      !number of branches, can be set to one for simplicity
    plt_ew%RawIsoTCanopy2Atm_pft(NZ)        = RawIsoTCanopy2Atm_pft(NZ,NY,NX)    !canopy resistance, updated iteratively
    plt_morph%CanopyHeightLive_pft(NZ)      = CanopyHeightLive_pft(NZ,NY,NX)       !canopy height, set as prescribed input
    plt_morph%CanopyHeightDead_pft(NZ)      = CanopyHeightDead_pft(NZ,NY,NX)
    plt_morph%CanopyStemSurfArea_pft(NZ)    = CanopyStemSurfArea_pft(NZ,NY,NX)     !canopy stem area, set as prescribed input
    plt_morph%CanopyLeafArea_pft(NZ)        = CanopyLeafArea_pft(NZ,NY,NX)     !canopy leaf area, set as prescribed input
    plt_morph%MaxSoilLays4Root_pft(NZ)      = NK_col(NY,NX)                     !can be derived from root type, set to maximum for simplicity
    plt_morph%Myco_pft(NZ)                  = Myco_pft(NZ,NY,NX)                !no mycorrhizae
    plt_biom%CanopyLeafSheathC_pft(NZ)      = CanopyLeafSheathC_pft(NZ,NY,NX)    !need to convert from leaf area
    plt_biom%CanopySapwoodC_pft(NZ)         = CanopySapwoodC_pft(NZ,NY,NX)        !need to convert from stem area
    plt_biom%CanopyNonstElmConc_pft(:,NZ)=CanopyNonstElmConc_pft(:,NZ,NY,NX)
    ! SP prescribes the canopy total directly; no elemental partition is inferred.
    plt_biom%CanopyNonstOMConc_pft(NZ)=CanopyNonstOMConc_pft(NZ,NY,NX)
    ! SP prescribes the same RTNSC in every root layer; no elemental partition
    ! or nonstructural biomass pool is inferred from this osmotic input.
    if(ldo_sp_mode)then
      RootNonstOMConc_pvr(:,:,NZ,NY,NX)=PrescribedRootNonstOMConc_pft(NZ,NY,NX)
    else
      RootNonstOMConc_pvr(:,:,NZ,NY,NX)=sum(RootNonstructElmConc_pvr(:,:,:,NZ,NY,NX),dim=1)
    endif
    plt_biom%RootNonstOMConc_pvr(:,:,NZ)=0._r8
    plt_biom%RootMyco2ndStrutElms_pvr(:,:,:,NZ)=0._r8
    plt_biom%RootMyco2ndStrutElms_pvr(:,:,1:NL_col(NY,NX),NZ)= &
      RootMyco2ndStrutElms_pvr(:,:,1:NL_col(NY,NX),NZ,NY,NX)
    plt_site%PlantPopuLive_pft(NZ)          = PlantPopuLive_pft(NZ,NY,NX)
    plt_biom%ZERO4LeafVar_pft(NZ)           = ZERO4LeafVar_pft(NZ,NY,NX)
    plt_biom%ZERO4Groth_pft(NZ)             = ZERO4Groth_pft(NZ,NY,NX)
    plt_photo%RCS_pft(NZ)                   = RCS_pft(NZ,NY,NX)                    !read in pft parameter
    plt_photo%CuticleResist_pft(NZ)         = CuticleResist_pft(NZ,NY,NX)          !read in pft parameter
    plt_photo%H2OCuticleResist_pft(NZ)      = H2OCuticleResist_pft(NZ,NY,NX)       !set based on read in pft parameter
    DO N=1,Myco_pft(NZ,NY,NX)
      plt_morph%Root1stMaxRadius_pft(N,NZ) = Root1stMaxRadius_pft(N,NZ,NY,NX)
      plt_morph%Root2ndMaxRadius_pft(N,NZ) = Root2ndMaxRadius_pft(N,NZ,NY,NX)
      plt_morph%Root2ndAxialResist_pft(N,NZ)  = Root2ndAxialResist_pft(N,NZ,NY,NX)
      plt_morph%RootPorosity_pft(N,NZ)     = RootPorosity_pft(N,NZ,NY,NX)
      plt_morph%RootRadialResist_pft(N,NZ) = RootRadialResist_pft(N,NZ,NY,NX)
      plt_morph%Root1stMaxRadius1_pft(N,NZ) = Root1stMaxRadius1_pft(N,NZ,NY,NX)
      plt_morph%Root2ndMaxRadius1_pft(N,NZ) = Root2ndMaxRadius1_pft(N,NZ,NY,NX)
    ENDDO

    DO L=1,NK_col(NY,NX)
      plt_morph%Root1stXNumL_pvr(L,NZ)     = Root1stXNumL_pvr(L,NZ,NY,NX)
      plt_morph%CRootLumenArea_pvr(L,NZ)   = CRootLumenArea_pvr(L,NZ,NY,NX)
      plt_morph%RootFineFrac2Med_rpvr(:,L,NZ) = RootFineFrac2Med_rpvr(:,L,NZ,NY,NX)
      plt_morph%RootMediumLength_pvr(L,NZ) = RootMediumLength_pvr(L,NZ,NY,NX)
      plt_morph%MRootLumenArea_pvr(L,NZ)   = MRootLumenArea_pvr(L,NZ,NY,NX)
      DO N=1,Myco_pft(NZ,NY,NX)
        plt_morph%Root2ndXNumL_rpvr(N,L,NZ)   = Root2ndXNumL_rpvr(N,L,NZ,NY,NX)
        plt_biom%RootNonstructElmConc_pvr(:,N,L,NZ)=RootNonstructElmConc_pvr(:,N,L,NZ,NY,NX)
        plt_biom%RootNonstOMConc_pvr(N,L,NZ)=RootNonstOMConc_pvr(N,L,NZ,NY,NX)
        plt_morph%Root1stRadius_pvr(N,L,NZ) = Root1stRadius_pvr(N,L,NZ,NY,NX)
        plt_morph%Root2ndRadius_rpvr(N,L,NZ) = Root2ndRadius_rpvr(N,L,NZ,NY,NX)
        plt_morph%Root2ndEffLen4uptk_rpvr(N,L,NZ) = Root2ndEffLen4uptk_rpvr(N,L,NZ,NY,NX)
        plt_morph%RootLenDensPerPlant_pvr(N,L,NZ) = RootLenDensPerPlant_pvr(N,L,NZ,NY,NX)
        plt_morph%RootTotLenPerPlant_pvr(N,L,NZ)     = RootTotLenPerPlant_pvr(N,L,NZ,NY,NX)
        plt_morph%RootAbsorbLenPerPlant_pvr(N,L,NZ)=RootAbsorbLenPerPlant_pvr(N,L,NZ,NY,NX)
      ENDDO
    ENDDO
    !variables updated from surface energy calculation
    DO L=1,NumCanopyLayers
      plt_rad%RadSWCanopyLAbsroption_pft(L,NZ)   =RadSWCanopyLAbsroption_pft(L,NZ,NY,NX)
      DO  M=1,NumOfSkyAzimuthSects
        DO  N=1,NumLeafInclinationClasses
          plt_rad%RadTotPARAbsorption_zsec(N,M,L,NZ)    = RadTotPARAbsorption_zsec(N,M,L,NZ,NY,NX)
          plt_rad%RadDifPARAbsorption_zsec(N,M,L,NZ) = RadDifPARAbsorption_zsec(N,M,L,NZ,NY,NX)
        ENDDO
      ENDDO
      plt_rad%TAU_DirectSunLit(L) = TAU_DirectSunLit(L,NY,NX)
      plt_rad%TAU_DirectSunSha(L)         = TAU_DirectSunSha(L,NY,NX)
    ENDDO
    plt_rad%TAU_DirectSunLit(NumCanopyLayers+1) = TAU_DirectSunLit(NumCanopyLayers+1,NY,NX)
    plt_rad%TAU_DirectSunSha(NumCanopyLayers+1)         = TAU_DirectSunSha(NumCanopyLayers+1,NY,NX)

    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO K=1,MaxNodesPerBranch
        DO  L=1,NumCanopyLayers
          DO N=1,NumLeafInclinationClasses
             plt_morph%LeafAreaZsec_brch(N,L,K,NB,NZ)=LeafAreaZsec_brch(N,L,K,NB,NZ,NY,NX)
          ENDDO
        ENDDO
      ENDDO

      DO K=0,MaxNodesPerBranch
        plt_morph%LeafArea_node(K,NB,NZ)   = LeafArea_node(K,NB,NZ,NY,NX)
      ENDDO
    ENDDO
    !variables updated iteratively
    plt_ew%CanopyBiomWater_pft(NZ)         = CanopyBiomWater_pft(NZ,NY,NX)    !canopy biom water,    updated iteratively
    plt_ew%TKC_pft(NZ)                     = TKC_pft(NZ,NY,NX)                !canopy temeprature updated iteratively
    plt_ew%TKCanopy_pft(NZ)                = TKCanopy_pft(NZ,NY,NX)           !canopy temperature updated iteratively
    plt_photo%RawH2OCanopy2Atm_pft(NZ)        = RawH2OCanopy2Atm_pft(NZ,NY,NX)
    plt_photo%RawHeatCanopy2Atm_pft(NZ)        = RawHeatCanopy2Atm_pft(NZ,NY,NX)
    plt_ew%SnowOnCanopy_pft(NZ)            = SnowOnCanopy_pft(NZ,NY,NX)
    plt_ew%WatHeldOnCanopy_pft(NZ)         = WatHeldOnCanopy_pft(NZ,NY,NX)     !water held by canopy surface
    plt_rad%FracPARads2Canopy_pft(NZ)      = FracPARads2Canopy_pft(NZ,NY,NX)
    plt_rad%FracPARads2LiveCanopy_pft(NZ)  = FracPARads2LiveCanopy_pft(NZ,NY,NX)
    plt_ew%HeatXAir2PCan_pft(NZ)           = HeatXAir2PCan_pft(NZ,NY,NX)
    plt_rad%RadPARCanopyAbsorption_pft(NZ) = RadPARCanopyAbsorption_pft(NZ,NY,NX)      !computed from surface energy module
    plt_rad%RadSWCanopyAbsorption_pft(NZ)  = RadSWCanopyAbsorption_pft(NZ,NY,NX)       !computed from surface energy module
    plt_ew%RainIntcptByCanopy_pft(NZ)      = RainIntcptByCanopy_pft(NZ,NY,NX)  !computed from hour1, rainfall partition
    plt_ew%SnowIntcptByCanopy_pft(NZ)      = SnowIntcptByCanopy_pft(NZ,NY,NX)
    plt_pheno%IsPlantActive_pft(NZ)        = IsPlantActive_pft(NZ,NY,NX)       !lai >0,              active
    plt_ew%PSICanopy_pft(NZ)               = PSICanopy_pft(NZ,NY,NX)
    plt_pheno%TempOffset_pft(NZ)           = TempOffset_pft(NZ,NY,NX)             !set based on read in pft parameter
    plt_ew%OrganOsmoPsi0pt_pft(NZ)         = OrganOsmoPsi0pt_pft(NZ,NY,NX)          !read in pft parameter

  ENDDO
  ! HOUR1 clears the count. Prescribed phenology sets activity here instead
  ! of running PhenologyUpdate; keep restart filtering consistent with it.
  NumActivePlants_col(NY,NX) = sum(IsPlantActive_pft(1:NP0_col(NY,NX),NY,NX))
  plt_site%NumActivePlants = NumActivePlants_col(NY,NX)
  if(ldo_sp_mode)then
    ! SP uses atmospheric CO2 directly; it does not solve a canopy-air carbon
    ! inventory using the full plant model's respiratory and growth fluxes.
    plt_site%CO2E=CO2E_col(NY,NX)
    plt_bgcr%NetCO2Flx2Canopy_col=0._r8
    plt_rad%SineSunInclinationAngle_col=SineSunInclinationAngle_col(NY,NX)
    do NZ=1,NP0_col(NY,NX)
      call SendPlantTraits(NY,NX,NZ)
      plt_photo%O2I_pft(NZ)=O2I_pft(NZ,NY,NX)
      plt_photo%CO2CuticleResist_pft(NZ)=CO2CuticleResist_pft(NZ,NY,NX)
      plt_morph%ClumpFactorNow_pft(NZ)=ClumpFactorNow_pft(NZ,NY,NX)
      call PreparePrescribedLeafBiochemistry(I,J,NZ,PrescribedLeafProteinDensity_pft(NZ,NY,NX), &
        CanopyLeafAreaZ_pft(:,NZ,NY,NX))
    enddo
  endif
  end subroutine PlantUptakeAPISend
!------------------------------------------------------------------------------------------

  subroutine PlantUPtakeAPIRecv(I,J,NY,NX)
  ! Description:
  ! Transfer prescribed-canopy fluxes and persistent exchange state to native arrays.

  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer :: N,L,NZ

  Air_Heat_Latent_store_col(NY,NX) = plt_ew%Air_Heat_Latent_store_col
  Air_Heat_Sens_store_col(NY,NX)   = plt_ew%Air_Heat_Sens_store_col
  Eco_NetRad_col(NY,NX)            = plt_rad%Eco_NetRad_col
  Eco_Heat_Latent_col(NY,NX)       = plt_ew%Eco_Heat_Latent_col
  Eco_Heat_Sens_col(NY,NX)         = plt_ew%Eco_Heat_Sens_col
  Eco_Heat_GrndSurf_col(NY,NX)     = plt_ew%Eco_Heat_GrndSurf_col
  CanopyBiomWater_col(NY,NX)       = plt_ew%CanopyBiomWater_col
  WatHeldOnCanopy_col(NY,NX)       = plt_ew%WatHeldOnCanopy_col
  SnowOnCanopy_col(NY,NX)          = plt_ew%SnowOnCanopy_col
  QVegET_col(NY,NX)                = plt_ew%QVegET_col
  VapXAir2Canopy_col(NY,NX)        = plt_ew%VapXAir2Canopy_col
  CanopyHeatStor_col(NY,NX)        = plt_ew%CanopyHeatStor_col
  HeatSensVapAir2Canopy_col(NY,NX)        = plt_ew%HeatSensVapAir2Canopy_col
  HeatPrecAir2Canopy_col(NY,NX)        = plt_ew%HeatPrecAir2Canopy_col
  HeatFlx2Canopy_col(NY,NX)        = plt_ew%HeatFlx2Canopy_col
  LWRadCanG_col(NY,NX)             = plt_ew%LWRadCanG
  ! Preserve the accepted root-transfer queue for the next hydrology step.
  do L=1,NL_col(NY,NX)
    TWaterPlantRoot2Soil_vr(L,NY,NX)=plt_ew%TWaterPlantRoot2Soil_vr(L)
    THeatLossRoot2Soil_vr(L,NY,NX)=plt_ew%THeatLossRoot2Soil_vr(L)
  enddo
  DO NZ=1,NP0_col(NY,NX)
    CdH2ORootxSoil_pft(NZ,NY,NX)     = plt_ew%CdH2ORootxSoil_pft(NZ)
    ENGYX_pft(NZ,NY,NX)              = plt_ew%ENGYX_pft(NZ)
    ETCanopy_CumYr_pft(NZ,NY,NX)     = plt_ew%ETCanopy_CumYr_pft(NZ)
    DeltaTKC_pft(NZ,NY,NX)           = plt_ew%DeltaTKC_pft(NZ)
    CanopyEvapTransLHeat_pft(NZ,NY,NX)     = plt_ew%CanopyEvapTransLHeat_pft(NZ)
    HeatXAir2PCan_pft(NZ,NY,NX)      = plt_ew%HeatXAir2PCan_pft(NZ)
    HeatStorCanopy_pft(NZ,NY,NX)    = plt_ew%HeatStorCanopy_pft(NZ)
    VHeatCapCanopy_pft(NZ,NY,NX)    = plt_ew%VHeatCapCanopy_pft(NZ)
    TKC_pft(NZ,NY,NX)                = plt_ew%TKC_pft(NZ)
    Transpiration_pft(NZ,NY,NX)      = plt_ew%Transpiration_pft(NZ)
    VapXAir2Canopy_pft(NZ,NY,NX)     = plt_ew%VapXAir2Canopy_pft(NZ)
    LWRadCanopy_pft(NZ,NY,NX)        = plt_rad%LWRadCanopy_pft(NZ)
    RadNet2Canopy_pft(NZ,NY,NX)      = plt_rad%RadNet2Canopy_pft(NZ)
    TKCanopy_pft(NZ,NY,NX)           = plt_ew%TKCanopy_pft(NZ)
    PSICanopy_pft(NZ,NY,NX)          = plt_ew%PSICanopy_pft(NZ)
    !plant tissue water is carried between steps (as PlantAPIRecv does for the
    !full plant model); ROOTUPTAKES updates it by transpiration - root uptake
    CanopyBiomWater_pft(NZ,NY,NX)    = plt_ew%CanopyBiomWater_pft(NZ)
    WatHeldOnCanopy_pft(NZ,NY,NX)    = plt_ew%WatHeldOnCanopy_pft(NZ)
    SnowOnCanopy_pft(NZ,NY,NX)       = plt_ew%SnowOnCanopy_pft(NZ)
    PSICanopyOsmo_pft(NZ,NY,NX)      = plt_ew%PSICanopyOsmo_pft(NZ)
    PSICanopyTurg_pft(NZ,NY,NX)      = plt_ew%PSICanopyTurg_pft(NZ)
    CanPStomaResistH2O_pft(NZ,NY,NX) = plt_photo%CanPStomaResistH2O_pft(NZ)
    if(ldo_sp_mode)then
      ! Retain the biochemical capacity and resistance diagnostics used by
      ! StomatalDynamics. SP does not run ComputeGPP or accumulate fixed carbon.
      CanopyGasCO2_pft(NZ,NY,NX)=plt_photo%CanopyGasCO2_pft(NZ)
      CanopyMinStomaResistH2O_pft(NZ,NY,NX)=plt_photo%CanopyMinStomaResistH2O_pft(NZ)
      LeafProteinCperm2LA_pft(NZ,NY,NX)=plt_biom%LeafProteinCperm2LA_pft(NZ)
      LeafRubiscoCperm2LA_pft(NZ,NY,NX)=plt_biom%LeafRubiscoCperm2LA_pft(NZ)
      LeafPEPCperm2LA_pft(NZ,NY,NX)=plt_biom%LeafPEPCperm2LA_pft(NZ)
      LeafC3ChlCperm2LA_pft(NZ,NY,NX)=plt_biom%LeafC3ChlCperm2LA_pft(NZ)
      LeafC4ChlCperm2LA_pft(NZ,NY,NX)=plt_biom%LeafC4ChlCperm2LA_pft(NZ)
      CanopyVcMaxRubisco25C_pft(NZ,NY,NX)=plt_photo%CanopyVcMaxRubisco25C_pft(NZ)
      CanopyVoMaxRubisco25C_pft(NZ,NY,NX)=plt_photo%CanopyVoMaxRubisco25C_pft(NZ)
      CanopyVcMaxPEP25C_pft(NZ,NY,NX)=plt_photo%CanopyVcMaxPEP25C_pft(NZ)
      ElectronTransptJmax25C_pft(NZ,NY,NX)=plt_photo%ElectronTransptJmax25C_pft(NZ)
      TFN_Carboxy_pft(NZ,NY,NX)=plt_photo%TFN_Carboxy_pft(NZ)
      TFN_Oxygen_pft(NZ,NY,NX)=plt_photo%TFN_Oxygen_pft(NZ)
      TFN_eTranspt_pft(NZ,NY,NX)=plt_photo%TFN_eTranspt_pft(NZ)
      LeafIntracellularCO2_pft(NZ,NY,NX)=plt_photo%LeafIntracellularCO2_pft(NZ)
      aquCO2Intraleaf_pft(NZ,NY,NX)=plt_photo%aquCO2Intraleaf_pft(NZ)
      SpecificLeafArea_pft(NZ,NY,NX)=plt_biom%SpecificLeafArea_pft(NZ)
      Km4RubiscoCarboxy_pft(NZ,NY,NX)=plt_photo%Km4RubiscoCarboxy_pft(NZ)
      O2L_pft(NZ,NY,NX)=plt_photo%O2L_pft(NZ)
    endif
    ! Preserve the solved canopy resistances for history, restart and the next API send.
    RawIsoTCanopy2Atm_pft(NZ,NY,NX) = plt_ew%RawIsoTCanopy2Atm_pft(NZ)
    RawH2OCanopy2Atm_pft(NZ,NY,NX)   = plt_photo%RawH2OCanopy2Atm_pft(NZ)
    RawHeatCanopy2Atm_pft(NZ,NY,NX)  = plt_photo%RawHeatCanopy2Atm_pft(NZ)
    TdegCCanopy_pft(NZ,NY,NX)        = plt_ew%TdegCCanopy_pft(NZ)
    DO L=1,NK_col(NY,NX)
      DO N=1,Myco_pft(NZ,NY,NX)
        PSIRootOSMO_vr(N,L,NZ,NY,NX)        = plt_ew%PSIRootOSMO_vr(N,L,NZ)
        PSIRootTurg_vr(N,L,NZ,NY,NX)        = plt_ew%PSIRootTurg_vr(N,L,NZ)
        RPlantRootH2OUptk_pvr(N,L,NZ,NY,NX) = plt_ew%RPlantRootH2OUptk_pvr(N,L,NZ)
      ENDDO
    ENDDO
  ENDDO

  end subroutine PlantUPtakeAPIRecv

end module PlantAPI4Uptake
