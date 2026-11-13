module ATSEcoSIMAdvanceMod
  !
  !Description

  use data_kind_mod, only : r8 => DAT_KIND_R8, yearIJ_type
  use StartqMod,     only : InitPlantPhenoMorphoBio, InitDimensionsandUptake
  use SoilWaterDataType
  use SoilHeatDataType
  use SharedDataMod
  use GridDataType
  use GridConsts
  use SOMDataType
  USE SoilPhysDataType
  use SnowDataType
  use LandSurfDataType
  use PlantTraitDataType
  use PlantMgmtDataType
  use CanopyDataType
  !use PlantAPIData, only: CO2E, CH4E, OXYE, Z2GE, Z2OE, ZNH3E, &
  !    H2GE
  use ClimForcDataType, only : LWRadSky_col, TairK_col, &
      VPA_col, WindSpeedAtm_col, RainH, VPK_col
  use SoilPropertyDataType
  use HydroThermData
  use SurfPhysData
  use SnowPhysData
  use CanopyHydroMod
  use EcoSIMSolverPar, only : NPH, dts_HeatWatTP
  use UnitMod    , only : units
  use EcoSIMCtrlDataType
  use MiniMathMod
  use ClimForcDataType
  use PrescribePhenolMod
  use RootDataType
  use FlagDataType
  use ATSUtilsMod
  use EcoSIMHistMod
  use PlantTraitsTransferMod
  !use PlantAPIData
  use EcoSIMCtrlMod,      only: ldo_sp_mode
  use EcosimBGCFluxType,  only: Eco_NetRad_col, Eco_Heat_Latent_col, &
    Eco_Heat_Sens_col, Eco_Heat_GrndSurf_col


implicit none
  character(len=*), private, parameter :: mod_filename=&
  __FILE__
  public :: RunEcoSIMSurfaceBalance
  contains

  subroutine RunEcoSIMSurfaceBalance(NYS)
  ! Description:
  ! Stage and advance the surface balance using coupled soil states and the
  ! shared snow, litter and atmospheric exchange preparation.
  !
  use EcoSimConst
  use GridMod           , only : SetMeshATS
  use SurfPhysMod       , only : RunSurfacePhysModelM, StageSurfacePhysModel, UpdateSurfaceAtM, &
      SetHourlyAccumulatorsATS
  use SnowBalanceMod    , only : SnowMassUpdate
  use StartsMod         , only : set_ecosim_solver
  use SnowBalanceMod    , only : SnowMassUpdate, SnowpackLayering
  use SurfaceRadiationMod
  !use SoilHydroParaMod  , only : SetColdRunSoilStates
  use PlantMod
  use WthrMod
  use PlantInfoMod
  !InitPlantMod is not used: its InitPlantPhenoMorphoBio(NZ) would hide the
  !column versions from StartqMod (module level) called in the PFT loop
  use MiniFuncMod, only : GetDayLength
  !use InitEcoSIM
  implicit none
  integer, intent(in) :: NYS  !Number of columns?

  integer :: NY,NX,NZ,L,NHW,NHE,NVN,NVS, I, J, M, heat_vec_size, NPH_Test, LS, npfts, NP
  real(r8) :: Wat_next
  real(r8) :: YSIN(NumOfSkyAzimuthSects),YCOS(NumOfSkyAzimuthSects),SkyAzimuthAngle(NumOfSkyAzimuthSects)
  real(r8) :: ResistVapThroughLitr(JY,JX)
  real(r8) :: RainEkReducedKsat(JY,JX)
  real(r8) :: HeatFluxAir2Soi(JY,JX)
  real(r8) :: TopLayWatVol_col(JY,JX)
  real(r8) :: Qinfl2MicP_col(JY,JX)
  real(r8) :: HeatInfl2Soil(JY,JX)
  real(r8) :: Qinfl2MacP_col(JY,JX)
  real(r8) :: dLWRaddTKsoi1(JY,JX)
  !real(r8) :: SnowDepth_col(JY,JX)
  real(r8) :: PrecAsRain(JY,JX)
  real(r8) :: PrecAsSnow(JY,JX)
  real(r8) :: VLTSoiPore
  real(r8) :: VPS(JY,JX)
  real(r8) :: EMM, CLD
  real(r8) :: SkyLWEst  !EcoSIM estimate of incoming longwave [MJ m-2 h-1]
  real(r8) :: dTexp, Tfinal
  real(r8), PARAMETER :: TSNOW=-0.25_r8  !oC, threshold temperature for snowfall
  real(r8), PARAMETER :: surf_e_max=0.1_r8
  real(r8), PARAMETER :: evap_sat_floor=0.3_r8  !top-cell liquid saturation below which soil evaporation stops [-]
  real(r8) :: Qinfl2MicPM(JY,JX)
  real(r8) :: Hinfl2SoilM(JY,JX)
  real(r8) :: VLWat_test(JZ,JY,JX)
  real(r8) :: THETF, flux_per_area
  type(yearIJ_type) :: yearIJ

  !associate( RPlantRootH2OUptk_pvr     => plt_ew%RPlantRootH2OUptk_pvr         & !inoput :root water uptake, [m3 d-2 h-1]
  !)

  !All the necessary sizes are taken from GridConsts
  !real(r8) :: LeafAreaZsecLive_lpft(NumLeafInclinationClasses,NumCanopyLayers,JP)
  !real(r8) :: StemAreaZsecLive_lpft(NumLeafInclinationClasses,NumCanopyLayers,JP)

  allocate(a_VEG_test(5))
  
  NHW=1;NHE=1;NVN=1;NVS=NYS
  !clock from ATS: day of year I (1-365) and hour J (1-24); EcoSIM's hour J
  !covers [J-1, J) h, here the hour starting at this call. The sun position
  !(PrepHourlyWeather) and the plant model use them; they were fixed at I=1,
  !J=1 (sun always below the horizon, no shortwave) before.
  I=current_day+1;J=current_hour+1
  NPH_Test=1
  NX=1
  npfts=1
  !ldo_sp_mode = .True.
  !LeafAreaZsecLive_lpft(:,:,:) = 0.2
  !StemAreaZsecLive_lpft(:,:,:) = 0.05

  !what Day/Month is it?
  call ComputeDatefromATS(current_day, current_year, current_month, day_of_month, total_days_in_month)
  write(*,*) "(ATSEcoSIMAdvance) month: ", current_month, " day: ", day_of_month, " of ", total_days_in_month
  yearIJ%year=current_year
  yearIJ%J=J
  yearIJ%I=I

  call SetMeshATS(NHW,NVN,NHE,NVS)
  !call InitModules()
  !load NK_col here?
  !NK_col(NX,1) = 14

  do NY=1, NYS
    !Hourly reset of the canopy column totals. extracts (ExtractCanopyFluxes)
    !rebuilds these every hour by adding the per-PFT values to them. Standalone
    !EcoSIM zeroes them every hour in SetHourlyDiagnostics (Hour1Mod.F90),
    !which the ATS coupling (prescribed phenology) does not call; without this
    !block they become running totals. This block replaces that missing reset,
    !at the same point of the hour (before the canopy radiation and the plant
    !model). SnowOnCanopy_col is read by the canopy radiation (effective LAI);
    !the others are diagnostics or outputs. The per-PFT state (e.g.
    !SnowOnCanopy_pft, WatHeldOnCanopy_pft) is not reset.
    QVegET_col(NY,NX)            = 0._r8
    CanopyBiomWater_col(NY,NX)   = 0._r8
    WatHeldOnCanopy_col(NY,NX)   = 0._r8
    SnowOnCanopy_col(NY,NX)      = 0._r8
    VapXAir2Canopy_col(NY,NX)    = 0._r8
    Eco_NetRad_col(NY,NX)        = 0._r8
    Eco_Heat_Latent_col(NY,NX)   = 0._r8
    Eco_Heat_Sens_col(NY,NX)     = 0._r8
    Eco_Heat_GrndSurf_col(NY,NX) = 0._r8
    CanopyHeatStor_col(NY,NX)    = 0._r8
    HeatFlx2Canopy_col(NY,NX)    = 0._r8
    !outputs summed over PFTs at the end of the advance; zeroed so they hold
    !this hour's value (a_EvapCan was previously never reset, a running total)
    a_Transpiration(NY) = 0._r8
    a_EvapCan(NY)       = 0._r8
    call SetHourlyAccumulatorsATS(NY,NX)
    NK_col(NY,NX) = a_NL(NY)
    Myco_pft(1,NY,NX) = 1
    NP0_col(NY,NX) = 1
    DO L=a_NU(NY),a_NL(NY)
      a_SSWS(L,NY) = 0.0_r8
    enddo
  enddo
  write(*,*) "(ATS-EcoSIM Advance) Day: ", current_day, " Year: ", current_year

  !Weather from ATS: PrepHourlyWeather takes the IWTHR=-999 branch and starts
  !from the ATS incident shortwave of each column, which CalcRadiation splits into
  !direct/diffuse with the sun position (and the canopy model then between canopy
  !and ground). IWTHR is otherwise unset in ATS mode, which selected the
  !hourly-file branch (no data, zero shortwave) and overwrote the ATS value.
  IWTHR = -999
  do NY=1,NYS
    RadSWATS_col(NY,NX) = swrad(NY)*0.0036_r8   !W m-2 to MJ m-2 h-1
  enddo

  call PrepHourlyWeather(I,J,NHW,NHE,NVN,NVS)
  do NY=1,NYS
    NU_col(NY,NX)               = a_NU(NY)
    NL_col(NY,NX)               = a_NL(NY)

    ASP_col(NY,NX)=a_ASP(NY)
    !TairKClimMean_col(NY,NX) = a_ATKA(NY)
    !CO2E_col(NY,NX)      = atm_co2
    !CH4E_col(NY,NX)      = atm_ch4
    !OXYE_col(NY,NX)      = atm_o2
    !Z2GE_col(NY,NX)      = atm_n2
    !Z2OE_col(NY,NX)      = atm_n2o
    !ZNH3E_col(NY,NX)     = atm_nh3
    !H2GE_col(NY,NX)      = atm_H2
    TairK_col(NY,NX)      = tairc(NY)
    !convert VPA from ATS units (Pa) to EcoSIM (MPa)
    !VPA(NY,NX) = vpair(NY)/1.0e6_r8

    !site from ATS (PK parameters "solar noon [h]", "latitude [degrees]")
    SolarNoonHour_col(NY,NX) = site_solar_noon
    ALAT_col(NY,NX)       = site_latitude
    DayLensCurr_col(NY,NX) = GetDayLength(ALAT_col(NY,NX),I)   !day length [h] (DayMod in standalone)
    !VPS(NY,NX)              = vapsat0(TairK_col(NY,NX))*EXP(-ALTI_col(NY,NX)/7272.0_r8)
    VPK_col(NY,NX)          = vpair(NY)/1.0e3 !vapor pressure in kPa
    !VPK_col(NY,NX)          = AMIN1(VPK_col(NY,NX),VPS(NY,NX))
    VPA_col(NY,NX)              = VPK_col(NY,NX)*2.173E-03_r8/TairK_col(NY,NX)
    !convert WindSpeedAtm_col from ATS units (m s^-1) to EcoSIM (m h^-1)
    WindSpeedAtm_col(NY,NX) = uwind(NY)*3600.0_r8
    !converting radiation units from ATS (W m^-2) to EcoSIM (MJ m^-2 h^-1)
    RadSWGrnd_col(NY,NX) = swrad(NY)*0.0036_r8
    
    !Set snow and initialize snow layers
    SnowAlbedo_col(NY,NX) = a_SALB(NY)
    SoilAlbedo_col(NY,NX) = 0.2_r8
    SnowDepth_col(NY,NX) = surf_snow_depth(NY)
    !call InitSnowLayers(NY,NX) 
    
    !EMM = 2.445 !There is a more elaborate calcuation of sky emissivity but I don't think we'll need that yet
    !EMM = 0.684
    !Incoming longwave: use the ATS value (sunrad, [W m-2]) when present, as
    !standalone EcoSIM does for longwave from the weather file. Otherwise
    !estimate it as in WthrMod; cloudiness is unknown here, so CLD = 0.2 (clear)
    CLD = 0.2_r8
    EMM = 0.625_r8*AMAX1(1.0_r8,(1.0E+03_r8*VPK_col(NY,NX)/TairK_col(NY,NX))**0.131_r8)
    EMM = EMM*(1.0_r8+0.242_r8*CLD**0.583_r8)
    SkyLWEst = EMM*stefboltz_const*TairK_col(NY,NX)**4._r8   ![MJ m-2 h-1]
    if(sunrad(NY) > 0._r8)then
      SkyLonwRad_col(NY,NX) = sunrad(NY)*0.0036_r8   !W m-2 to MJ m-2 h-1
    else
      SkyLonwRad_col(NY,NX) = SkyLWEst
    endif
    LWRadSky_col(NY,NX) = SkyLonwRad_col(NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    TCA_col(NY,NX) = units%Kelvin2Celcius(TairK_col(NY,NX))
    !SineSunInclinationAngle_col(NY,NX) = 1.0_r8

    DO L=NU_col(NY,NX),NL_col(NY,NX)
      CumDepz2LayBottom_vr(L,NY,NX) = a_CumDepz2LayBottom_vr(L,NY)
      !Convert Bulk Density from ATS (kg m^-3) to EcoSIM (Mg m^-3)
      !SoilBulkDensity_vr(L,NY,NX)  = a_BKDSI(L,NY)/1.0e3_r8
      SoilBulkDensity_vr = a_RDENS(L,NY)*(1.0_r8-a_PORO(L,NY))/1.0e3_r8
      CSoilOrgM_vr(ielmc,L,NY,NX)  = a_CORGC(L,NY)
      CSoilOrgM_vr(ielmn,L,NY,NX)  = a_CORGN(L,NY)
      CSoilOrgM_vr(ielmp,L,NY,NX)  = a_CORGP(L,NY)
      !Convert ATS units (mols) to EcoSIM units (mH2O)
      VLWatMicP1_vr(L,NY,NX)       = a_WC(L,NY)/(a_LDENS(L,NY)*AREA_3D(3,NU_col(NY,NX),NY,NX))

      VLiceMicP1_vr(L,NY,NX)       = 0.0
      TKSoil1_vr(L,NY,NX)          = a_TEMP(L,NY)
      TKS_vr(L,NY,NX)              = a_TEMP(L,NY)
      VHeatCapacity1_vr(L,NY,NX)   = heat_capacity
      SoilFracAsMicP_vr(L,NY,NX)   = 1.0 !This is percentage of Soil void that is micropores (100% for now)
      !Convert Matric Pressure from ATS [Pa] to EcoSIM [MPa]
      PSISM1_vr(L,NY,NX)           = -1.0*a_MATP(L,NY)/1.0e6
      ElvAdjstedSoilH2OPSIMPa_vr(L,NY,NX) = PSISM1_vr(L,NY,NX) !needed for root uptake
      POROS_vr(L,NY,NX)            = a_PORO(L,NY)
      !AREA3(L,NY,NX)              = a_AREA3(L,NY)
      VLTSoiPore = VLSoilMicP_vr(L,NY,NX)

      !Additional parameters for phenology
      !VLSoilPoreMicP_vr(L,NY,NX) = SoilFracAsMicP_vr(L,NY,NX) !total soil volume (associated with micropores) in layer L, which is soil volume when no macropore is considered
      !VLSoilPoreMicP_vr(L,NY,NX) =
      !VLSoilMicP_vr(L,NY,NX)=VLSoilPoreMicP_vr(L,NY,NX)
      VLWatMicP_vr(L,NY,NX) = VLWatMicP1_vr(L,NY,NX) !soil moisture
      !HYCDMicP4RootUptake_vr(L,NY,NX) = PSISM1_vr(L,NY,NX)!soil hydraulic conductivity for water uptake

      !Hydraulic conductivity is hard-coded for now this will be fixed pending a decision on how best to
      ! get it out of the ATS WRM model
      HYCDMicP4RootUptake_vr(L,NY,NX) = 0.000571
      THETW_vr(L,NY,NX) = a_LSAT(L,NY) !relative Saturation of soil micropores in layerL

      !Parameters needed for water balance after phenology is added:
      !Do I need to compute VGeom or is it already there
      FracSoiAsMicP_vr(L,NY,NX)  = POROS_vr(L,NY,NX)*SoilFracAsMicP_vr(L,NY,NX)
      VGeomLayer_vr(L,NY,NX)     = AREA_3D(3,L,NY,NX)*DLYR_3D(3,L,NY,NX)
      VLSoilPoreMicP_vr(L,NY,NX) = VGeomLayer_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
      VLSoilMicP_vr(L,NY,NX)=VLSoilPoreMicP_vr(L,NY,NX)
      VLMicP_vr(L,NY,NX)    = POROS_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
      VLSoilMicPMass_vr(L,NY,NX) = SoilBulkDensity_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)

      WiltPoint_vr(L,NY,NX) = 0.25_r8
      FieldCapacity_vr(L,NY,NX) = 0.5_r8
      LOGFldCapacity_vr(L,NY,NX) = LOG(FieldCapacity_vr(L,NY,NX))
      LOGWiltPoint_vr(L,NY,NX)   = LOG(WiltPoint_vr(L,NY,NX))
      LOGPOROS_vr(L,NY,NX)=LOG(POROS_vr(L,NY,NX))
      PSD_vr(L,NY,NX)               = LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX)
      FCD_vr(L,NY,NX)               = LOGFldCapacity_vr(L,NY,NX)-LOGWiltPoint_vr(L,NY,NX)
      SRP_vr(L,NY,NX)               = 1.00_r8
      LOGPSIAtSat(NY,NX)            = LOG(-PSIPS)
      LOGPSIFLD_col(NY,NX)          = LOG(-pressure_at_field_capacity)
      LOGPSIMXD_col(NY,NX)          = LOGPSIFLD_col(NY,NX)-LOGPSIAtSat(NY,NX)

      !Top layer water volume:
      !Only the top-cell water above evap_sat_floor can evaporate within the
      !hour, so the hourly sink sent to ATS cannot dry the cell into the steep,
      !near-residual end of its retention curve
      TopLayWatVol_col(NY,NX) = VLWatMicP1_vr(NU_col(NY,NX),NY,NX) &
        *AZMAX1(1._r8-evap_sat_floor/AMAX1(a_LSAT(NU_col(NY,NX),NY),1.e-6_r8))

      !Adapting from SoilHydroParaMod.F90 for the no organic matter in soil case
      !THETF=AMIN1(POROS_vr(L,NY,NX),EXP((LOGPSIAtSat(NY,NX)-LOG(0.033_r8)) &
      !  *(LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX))/LOGPSIMXD_col(NY,NX)+LOGPOROS_vr(L,NY,NX)))
      !HYCDMicP4RootUptake_vr(L,NY,NX)=1.54_r8*((POROS_vr(L,NY,NX)-THETF)/THETF)**2

      !Root2ndMaxRadius_pft(1,L) = 0.0002_r8
      !Root2ndMaxRadius_pft(2,L) = 5.0e-06_r8
      !call SetColdRunSoilStates(I,J,L,NY,NX)

      IF(VLTSoiPore.GT.ZEROS2(NY,NX))THEN
        !fraction as water
        FracSoiPAsWat_vr(L,NY,NX)=AZMAX1t(VLWatMicP1_vr(L,NY,NX)/VLTSoiPore)
        !fraction as ice
        FracSoiPAsIce_vr(L,NY,NX)=AZMAX1t(VLiceMicP1_vr(L,NY,NX)/VLTSoiPore)
        !fraction as air
        FracAirFilledSoilPore_vr(L,NY,NX)=AZMAX1t(VLairMicP1_vr(L,NY,NX)/VLTSoiPore)
      ELSE
        FracSoiPAsWat_vr(L,NY,NX)=POROS_vr(L,NY,NX)
        FracSoiPAsIce_vr(L,NY,NX)=0.0_r8
        FracAirFilledSoilPore_vr(L,NY,NX)=0.0_r8
      ENDIF
    ENDDO

    !The surface energy balance (longwave emission, sensible heat, soil vapor
    !pressure, snow-soil conduction) reads the top soil layer temperature. Use
    !the ATS surface temperature there instead of the top subsurface cell
    !temperature, which can be much deeper than the surface and does not see the
    !cooling or warming applied to the surface by surf_e_source.
    if(a_TSURF(NY) > 0._r8)then
      TKSoil1_vr(NU_col(NY,NX),NY,NX) = a_TSURF(NY)
      TKS_vr(NU_col(NY,NX),NY,NX)     = a_TSURF(NY)
    endif

    !Setting the snow, if passed as total precipitation do full temp calc, else
    !Just set snow and rain to their given values
    !5/16 - I don't think we need to account for area here at all
    !       I am switching out PrecAsRain and just directly substituting
    !       with RainFalPrec
    if(p_bool)then
      IF(TCA_col(NY,NX).GT.TSNOW)THEN
        RainFalPrec_col(NY,NX)=p_total(NY)*3600.0*AREA_3D(3,NU_col(NY,NX),NY,NX) !convert from m/s to m/hr
        SnoFalPrec_col(NY,NX)=0.0_r8
      ELSE
        RainFalPrec_col(NY,NX)=0.0_r8
        SnoFalPrec_col(NY,NX)=p_total(NY)*3600.0*AREA_3D(3,NU_col(NY,NX),NY,NX) !convert to m SWE/s to m SWE/hr
      ENDIF
    else
      RainFalPrec_col(NY,NX)=p_rain(NY)*3600.0*AREA_3D(3,NU_col(NY,NX),NY,NX) !convert from m/s to m/hr
      SnoFalPrec_col(NY,NX)=p_snow(NY)*3600.0*AREA_3D(3,NU_col(NY,NX),NY,NX) !convert from m SWE/s to mSWE/hr
    endif
    !Set Prec equal to variables for after Irrigation and canopy processing
    !since that is not included yet
    PrecRainAndIrrig_col = RainFalPrec_col
    RainPrecThrufall_col = RainFalPrec_col


    !RainFalPrec_col(NY,NX)=PrecAsRain(NY,NX)*AREA(3,NU_col(NY,NX),NY,NX)
    !SnoFalPrec_col(NY,NX)=PrecAsSnow(NY,NX)*AREA(3,NU_col(NY,NX),NY,NX)
    POROS_vr(0,NY,NX) = 1.0

    !Fill in column-wise values needed for prescribed phenology
    !CanopyHeight_col(NY,NX) = 17.0
    !CanopyHeightLive_pft(1,NY,NX) = 17.0
    !tlai_day_pft(1,NY,NX) = a_LAI(NY)
    !tsai_day_pft(1,NY,NX) = a_SAI(NY)
    !irootType_col(NY,NX) = a_PFT(NY,1)
    !if(ldo_sp_mode) call PlantCanopyRadsModel(I,J,NY,NX,0.0_r8)
    !Fill number of plants from npfts
    if(a_PFT(1,NY).EQ.0.0)then
      NP0_col(NY,NX) = 0
      NP_col(NY,NX) = 0
    else
      NP0_col(NY,NX) = num_pfts
      NP_col(NY,NX) = num_pfts
    endif
    
    !if(ldo_sp_mode) call startq(NHW,NHE,NVN,NVS,1,num_pfts)

    do i=1,5
      a_VEG_test(i) = 1
    enddo
    a_VEG_test(1) = a_VEG(1)
    
    !prescribed phenology. Every hour: the canopy structure prescribed by ATS
    !and the plant traits (parameters only). The plant water/energy state
    !(canopy water potential, canopy temperatures, canopy heat storage, plant
    !water) is initialized once, at the first advance, and then carried from
    !hour to hour as in standalone EcoSIM; it is no longer reset every hour.
    if (ldo_sp_mode) then
      DO NZ=1,num_pfts
        if (a_PFT(NZ,NY) .GT. 0.0) then
            tlai_day_pft(NZ,NY,NX) = a_LAI(NZ,NY)*AREA_3D(3, NU_col(NY,NX), NY, NX)
            tsai_day_pft(NZ,NY,NX) = a_SAI(NZ,NY)*AREA_3D(3, NU_col(NY,NX), NY, NX)
            SnowOnCanopy_pft(NZ,NY,NX) = a_CanSnow(NZ,NY)
            SnowOnCanopy_col(NY,NX)    = SnowOnCanopy_col(NY,NX) + a_CanSnow(NZ,NY)

            irootType_col(NY,NX) = irootType_pft(NZ,NY,NX)
            CanopyHeight_col(NY,NX) = CanopyHeightLive_pft(NZ,NY,NX)
            !Load the PFT array from ATS and fill the plant traits
            ! based on that mapping
            DATAPI(NZ,NY,NX) = a_PFT(NZ,NY)
            call ReadPlantProperties(0,NZ,NY,NX,.false.)

            call InitDimensionsandUptake(NZ,NY,NX)

            call InitPlantPhenoMorphoBio(NZ,NY,NX)

            call SendPlantTraits(NY,NX,NZ)

            !Set cuticle resistances scaled from pft trait file
            H2OCuticleResist_pft(NZ,NY,NX) = CuticleResist_pft(NZ,NY,NX)/3600.0_r8
            CO2CuticleResist_pft(NZ,NY,NX) = CuticleResist_pft(NZ,NY,NX)*1.56_r8
            ATCA_col(NY,NX)                = 6.0_r8

            if (a_PlantInit(NY,1) .LT. 0.5_r8) call InitPrescribedPlantState(NZ,NY,NX)
        else
            !Empty PFT slot (zeros are packed at the end of a_PFT). NP_col still
            !includes it, so mark it inactive (DATAPI=0) with no canopy or roots
            DATAPI(NZ,NY,NX)               = 0
            tlai_day_pft(NZ,NY,NX)         = 0._r8
            tsai_day_pft(NZ,NY,NX)         = 0._r8
            irootType_pft(NZ,NY,NX)        = 0
            MeanSecondaryAxisLength_pft(NZ,NY,NX) = 0._r8
            CanopyNonstOMConc_pft(NZ,NY,NX) = 0._r8
            CanopyHeightLive_pft(NZ,NY,NX) = 0._r8
        endif
      enddo
      a_PlantInit(NY,1) = 1.0_r8
    endif
  ENDDO
  
  !Run precribe phenology interpolation which sets necessary variables
  !To do the phenology
  if(ldo_sp_mode) call PrescribePhenologyInterp(I, NHW, NHE, NVN, NVS)

  !PlantCanopyRadsModel processes the input radiation
  ! into radiation to ground and to canopy
  do NY=1,NYS
    if(ldo_sp_mode)then  
      call PlantCanopyRadsModel(I,J,NY,NX,0.0_r8)
      
      !SnowPrecThrufall_col(NY,NX) = SnoFalPrec_col(NY,NX)
      !RainPrecThrufall_col(NY,NX) = RainFalPrec_col(NY,NX)
      !TKSnowThrufall_col(NY,NX) = TairK_col(NY,NX)
        
      call CanopyInterceptprecip(NY,NX)
     else
       !If not using phenology set the thrufall snow and rain
       !to the total
       SnowPrecThrufall_col(NY,NX) = SnoFalPrec_col(NY,NX)
       RainPrecThrufall_col(NY,NX) = RainFalPrec_col(NY,NX)
       TKSnowThrufall_col(NY,NX) = TairK_col(NY,NX)
     endif
  enddo

  PSIAtFldCapacity_col = pressure_at_field_capacity
  PSIAtWiltPoint_col = pressure_at_wilting_point
  ! Reset hourly accumulators before shared staging so any staged transfer
  ! remains in the ledger consumed by the surface iterations.
  do NY=1, NYS
    Qinfl2MicP_col(NY,NX)=0.0_r8
    HeatInfl2Soil(NY,NX)=0.0_r8
    Qinflx2Soil_col(NY,NX)=0.0_r8
    HeatFlx2Grnd_col(NY,NX)=0.0_r8
    a_EvapSnow(NY)=0.0_r8
    a_Sublim(NY)=0.0_r8
    
    HeatFlowSno2SoiByWat_col(NY,NX)=0.0_r8
    HeatConvFlxSno2Soi_col(NY,NX)=0.0_r8
    HeatCndFlxSno2Soi_col(NY,NX)=0.0_r8
    HeatConvFlxLitr2Soi_col(NY,NX)=0.0_r8
    HeatCndFlxLitr2Soi_col(NY,NX)=0.0_r8
  
    HeatSensAir2Grnd_col(NY,NX)=0.0_r8
    Radnet2Grnd_col(NY,NX)=0.0_r8
    LatentHeatEvapAir2Grnd_col(NY,NX)=0.0_r8
    HeatSensVapAir2Grnd_col(NY,NX)=0.0_r8
  enddo

  !This does the surface energy and water balance. In the case of plants
  ! it processes the precipitation and radiation already separated out in
  ! PlantCanopyRadsModel
  call StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistVapThroughLitr)

  VHeatCapacity1_vr(0,1,1) = 0.0

  !Runs the main subcycling for the surface energy and water balance in
  !addition to snow.
  !
  DO M=1,NPH

    !call RunSurfacePhysModelM(I,J,M,NHE,NHW,NVS,NVN,ResistVapThroughLitr,&
    !  KSatReductByRainKineticEnergy,TopLayWatVol,HeatFluxAir2Soi,Qinfl2MicPM,Hinfl2SoilM)

    !call RunSurfacePhysModelM(I,J,M,NHE,NHW,NVS,NVN,ResistVapThroughLitr, &
    !KSatReductByRainKineticEnergy,TopLayWatVol_col,HeatFluxAir2Soi, &
    !Qinfl2MicP_col,HeatInfl2Soil)

    call RunSurfacePhysModelM(I,J,M,NHE,NHW,NVS,NVN,ResistVapThroughLitr,RainEkReducedKsat,&
      TopLayWatVol_col,HeatFluxAir2Soi,dLWRaddTKsoi1,Qinfl2MicP_col,HeatInfl2Soil,Qinfl2MacP_col)
    
    !call subroutine RunSurfacePhysModelM(I,J,M,NHE,NHW,NVS,NVN,ResistVapThroughLitr,RainEkReducedKsat,&
    !  TopLayWatVol_col,HeatFluxAir2Soi,Qinfl2MicP,HeatInfl2Soil,Qinfl2MacP)

    !sum fluxes over iterations in each column
    do NY=1, NYS
      HeatFlx2Grnd_col(NY,NX) = HeatFlx2Grnd_col(NY,NX)+HeatInfl2Soil(NY,NX)
      Qinflx2Soil_col(NY,NX)  = Qinflx2Soil_col(NY,NX)+Qinfl2MicP_col(NY,NX)
      !EVAPW_col/EVAPS_col are reset at every substep (InitSnowAccumsM), so
      !add them over the substeps to get the hourly snow evaporation/sublimation
      a_EvapSnow(NY)          = a_EvapSnow(NY)+EVAPW_col(NY,NX)
      a_Sublim(NY)            = a_Sublim(NY)+EVAPS_col(NY,NX)
    enddo

    !also update state variables for iteration M
    call UpdateSurfaceAtM(I,J,M,NHW,NHE,NVN,NVS)

  ENDDO

  !Does final update of the snow masses after the surface model is run
  do NY=1,NYS
    call SnowMassUpdate(I,J,NY,NX,Qinfl2MicPM(NY,NX),Hinfl2SoilM(NY,NX))
  ENDDO

  !Does energy and water exchange in the plant roots and canopy
  ! Specifically runs the CanopyEnergyH2Oiter_func
  if(ldo_sp_mode) call PlantModel(yearIJ,NHW,NHE,NVN,NVS)

  DO NY=1,NYS
    !for every column send the top layer to the transfer var
    !Convert heat and water flux flux from the subcycled value
    !The ATS units for water flux is m/s for energy it is MW/m^2
    !As EcoSIM is an hourly model, the conversion is done in ATS
    !write(*,*) "dts_HeatWatTP: ", dts_HeatWatTP
    !surf_e_source(NY) = (HeatFlx2Grnd_col(NY,1)-HeatFlowSno2SoiByWat_col(NY,1)) / (column_area(NY))
    !surf_e_source(NY) = MIN(HeatFlx2Grnd_col(NY,1) / (column_area(NY)), surf_e_max)
    !EcoSIM advects heat with water as cpw*T*Q (T in K, zero at 0 K) while ATS
    !internal energy is zero at TFice (IEM reference temperature, must match the
    !ATS input). Shift the energy of the exported water Qinflx2Soil_col to the
    !ATS reference, otherwise each m3 of infiltration adds ~cpw*273 MJ of heat.
    flux_per_area = (HeatFlx2Grnd_col(NY,1) - cpw*TFice*Qinflx2Soil_col(NY,1)) / column_area(NY)
    !Uncapped export:
    !surf_e_source(NY) = flux_per_area
    !Cap the magnitude at surf_e_max [MJ m-2 h-1], keeping the sign
    surf_e_source(NY) = SIGN(MIN(ABS(flux_per_area), surf_e_max), flux_per_area)
    
    dTexp = surf_e_source(NY)/heat_capacity
    Tfinal = TKsoil1_vr(1,NY,NX) + dTexp 
    
    !per grid cell [m3 d-2 h-1] to per area [m h-1]; ATS converts h-1 to s-1
    surf_w_source(NY) = Qinflx2Soil_col(NY,1) / column_area(NY)
    surf_snow_depth(NY) = SnowDepth_col(NY,1)
    !Now update subsurface flux from roots
    LWRadCanGPrev_col(NY,NX) = LWRadCanG_col(NY,NX)
    TLEX_col(NY,NX)             = Air_Heat_Latent_store_col(NY,NX)
    TSHX_col(NY,NX)             = Air_Heat_Sens_store_col(NY,NX)
    a_ET(NY) = QVegET_col(NY,NX) !canopy evapotranspiration

    !Now that the heat carryover variables are set to the coupler versions
    !They need to be cleared for the next iteration
    LWRadCanG_col(NY,NX) = 0.0
    Air_Heat_Latent_store_col(NY,NX) = 0.0
    Air_Heat_Sens_store_col(NY,NX) = 0.0

    !Sum over PFTs for total transpiration and canopy evap
    !and canopy on snow
    do NZ=1,num_pfts
      a_Transpiration(NY) = a_Transpiration(NY) + Transpiration_pft(NZ,NY,NX)
      a_EvapCan(NY)  = a_EvapCan(NY) + VapXAir2Canopy_pft(NZ,NY,NX)
      a_CanSnow(NZ,NY) = SnowOnCanopy_pft(NZ,NY,NX)
    enddo

    a_EvapGrnd(NY) = TEvapXAir2Toplay_col(NY,NX) !bare ground evaporation
    a_EvapLitr(NY) = TEvapXAir2LitR_col(NY,NX) !litter evaporation
    !a_EvapSnow, a_Sublim: added over the substeps in the M loop above
    !DO LS=1,JS
    !  a_TS(LS,L) = TKSnow1_snvr(LS,NY,NX)
    !ENDDO
    DO L=NU_col(NY,NX),NL_col(NY,NX)
      DO NZ=1,num_pfts
        a_SSWS(L,NY) = a_SSWS(L,NY) + RPlantRootH2OUptk_pvr(1,L,NZ,NY,NX)
      ENDDO
    ENDDO
  ENDDO

  call destroy(a_VEG_test)
  
  !end associate
  end subroutine RunEcoSIMSurfaceBalance

!------------------------------------------------------------------------------------------
  subroutine InitPrescribedPlantState(NZ,NY,NX)
  !
  !Description
  !Prescribed phenology: initial plant water and energy state of PFT NZ,
  !set at the first advance (when ATS's air temperature is available) and
  !then carried between hours and checkpointed (see ATSStateRegistryMod), as
  !in standalone EcoSIM (StartqMod).
  implicit none
  integer, intent(in) :: NZ,NY,NX

  !canopy energy state: canopy at air temperature, no stored heat
  ENGYX_pft(NZ,NY,NX)         = 0._r8
  DeltaTKC_pft(NZ,NY,NX)      = 0._r8
  TKC_pft(NZ,NY,NX)           = TairK_col(NY,NX)
  TKCanopy_pft(NZ,NY,NX)      = TairK_col(NY,NX)
  TdegCCanopy_pft(NZ,NY,NX)   = units%Kelvin2Celcius(TKC_pft(NZ,NY,NX))
  TCGroth_pft(NZ,NY,NX)       = TdegCCanopy_pft(NZ,NY,NX)
  TKGroth_pft(NZ,NY,NX)       = TKC_pft(NZ,NY,NX)
  fTCanopyGroth_pft(NZ,NY,NX) = 1.0_r8

  !canopy water potential: standalone initial value (StartqMod.F90)
  PSICanopy_pft(NZ,NY,NX)     = -1.0E-03_r8
  PSICanopyOsmo_pft(NZ,NY,NX) = OrganOsmoPsi0pt_pft(NZ,NY,NX)+PSICanopy_pft(NZ,NY,NX)
  PSICanopyTurg_pft(NZ,NY,NX) = AZMAX1(PSICanopy_pft(NZ,NY,NX)-PSICanopyOsmo_pft(NZ,NY,NX))

  !plant tissue water: unset, PrescribePhenologyInterp starts it in equilibrium
  !with the leaf carbon
  CanopyBiomWater_pft(NZ,NY,NX) = 0._r8
  end subroutine InitPrescribedPlantState

end module ATSEcoSIMAdvanceMod
