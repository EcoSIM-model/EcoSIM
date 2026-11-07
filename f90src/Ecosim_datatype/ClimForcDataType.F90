module ClimForcDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod, only : warming_exp
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8) :: RMAX       !maximum hourly radiation,	[MJ m-2 h-1]
  real(r8) :: TAVG1      !parameter to calculate hourly  air temperature from daily value,	[oC]
  real(r8) :: TAVG2      !parameter to calculate hourly  air temperature from daily value,	[oC]
  real(r8) :: TAVG3      !parameter to calculate hourly  air temperature from daily value,	[oC]
  real(r8) :: AMP1       !parameter to calculate hourly  air temperature from daily value,	[oC]
  real(r8) :: AMP2       !parameter to calculate hourly  air temperature from daily value,	[oC]
  real(r8) :: AMP3       !parameter to calculate hourly  air temperature from daily value,	[oC]
  real(r8) :: VAVG1      !parameter to calculate hourly  vapor pressure from daily value,	[kPa]
  real(r8) :: VAVG2      !parameter to calculate hourly  vapor pressure from daily value,	[kPa]
  real(r8) :: VAVG3      !parameter to calculate hourly  vapor pressure from daily value,	[kPa]
  real(r8) :: VMP1       !parameter to calculate hourly  vapor pressure from daily value,	[kPa]
  real(r8) :: VMP2       !parameter to calculate hourly  vapor pressure from daily value,	[kPa]
  real(r8) :: VMP3       !parameter to calculate hourly  vapor pressure from daily value,	[kPa]
  real(r8) :: SAZI                              !solar azimuth of solar angle, [-]
  real(r8) :: SCOS                              !cosine of solar angle, [-]
  real(r8) :: DOY                               !day of year, [-]

  real(r8) :: atm_co2_mon(12)    !monthly atmospheric O2, [ppmv]
  real(r8) :: atm_ch4_mon(12)    !monthly atmospheric CH4, [ppmv]
  real(r8) :: atm_n2o_mon(12)    !monthly atmospheric N2O, [ppmv]
  real(r8) :: TMPX(366)                         !maximum daily air temperature, [oC]
  real(r8) :: TMPN(366)                         !minimum daily air temperature, [oC]
  real(r8) :: SRAD(366)                         !daily solar radiation, [MJ m-2 d-1]
  real(r8) :: RAIN(366)                         !daily precipitation, [mm d-1 ]
  real(r8) :: WIND(366)                         !daily wind travel, [m d-1]
  real(r8) :: DWPT(2,366)                       !daily dewpoint temperature, [oC]
  real(r8) :: DryDepoOMC(366)                   !dry deposition of organic C, [gC/day/m2]

  real(r8) :: TMP_hrly(24,366)                   !hourly air temperature, [oC]
  real(r8) :: SWRad_hrly(24,366)                !hourly solar radiation, [MJ m-2 h-1]
  real(r8) :: RAINH(24,366)                     !hourly precipitation, [mm h-1]
  real(r8) :: WINDH(24,366)                     !hourly wind speed, [m h-1]
  real(r8) :: DWPTH(24,366)                     !hourly dewpoint temperature, [oC]
  real(r8) :: RadLWClm(24,366)                  !longwave radiation (MJ m-2 h-1)
  real(r8) :: PBOT_hrly(24,366)                 !hourly surface atmospheric pressure, [kPa]
  real(r8) :: DRAD(12)                          !change factor for radiation, [-]
  real(r8) :: DTMPX(12)                         !change factor for maximum temperature, [-]
  real(r8) :: DTMPN(12)                         !change factor for minimum temperature, [-]
  real(r8) :: DHUM(12)                          !change factor for humidity, [-]
  real(r8) :: DPREC(12)                         !change factor for precipitation, [-]
  real(r8) :: DWIND(12)                         !change factor for wind speed, [-]
  real(r8) :: DCN4R(12)                         !change factor for NH4 in precipitation, [-]
  real(r8) :: DCNOR(12)                         !change factor for NO3 in precipitation, [-]
  
  real(r8),target,allocatable ::  srad_scalar_col(:,:)     !solar radiation scalar due to open top chamber
  real(r8),target,allocatable ::  EMS_Modify_Scalar_col(:,:)      !emissivity scalar due to open top chamber
  real(r8),target,allocatable ::  Eco_RadSW_col(:,:)       !shortwave radiation absorbed by the ecosystem [MJ/h]
  real(r8),target,allocatable ::  RadSWATS_col(:,:)        !incident shortwave radiation from ATS (IWTHR=-999), [MJ m-2 h-1]
  real(r8),target,allocatable ::  TKS_ref_vr(:,:,:,:)      !reference tempeature profile from control run to warming experiment [K]
  real(r8),target,allocatable ::  TDTPX(:,:,:)                       !accumulated change  for maximum temperature, [-]
  real(r8),target,allocatable ::  TDTPN(:,:,:)                       !accumulated change  for minimum temperature, [-]
  real(r8),target,allocatable ::  TDRAD(:,:,:)                       !accumulated change  for radiation, [-]
  real(r8),target,allocatable ::  TDHUM(:,:,:)                       !accumulated change  for humidity, [-]
  real(r8),target,allocatable ::  TDPRC(:,:,:)                       !accumulated change  for precipitation, [-]
  real(r8),target,allocatable ::  TDWND(:,:,:)                       !accumulated change  for wind speed, [-]
  real(r8),target,allocatable ::  TDCN4(:,:,:)                       !accumulated change  for NH4 in precipitation, [-]
  real(r8),target,allocatable ::  TDCNO(:,:,:)                       !accumulated change  for NO3 in precipitation, [-]
  real(r8),target,allocatable ::  TCA_col(:,:)                           !air temperature, [oC]
  real(r8),target,allocatable ::  TairK_col(:,:)                         !air temperature, [K]
  real(r8),target,allocatable ::  WindSpeedAtm_col(:,:)                  !wind speed, [m h-1]
  real(r8),target,allocatable ::  VPA_col(:,:)                           !atmospheric vapor concentration, [m3 m-3]
  real(r8),target,allocatable ::  VPK_col(:,:)                           !atmospheric vapor pressure, [kPa]
  real(r8),target,allocatable ::  PBOT_col(:,:)                          !atmospheric pressure [kPa]
  real(r8),target,allocatable ::  DayLensCurr_col(:,:)                          !daylength, [h]
  real(r8),target,allocatable ::  DayLenthPrev_col(:,:)                          !daylength of previous day, [h]
  real(r8),target,allocatable ::  DayLenthMax_col(:,:)                          !maximum daylength, [h]
  real(r8),target,allocatable ::  OMEGA2Ground(:,:,:)                      !sine of solar beam on leaf surface, [-]
  real(r8),target,allocatable ::  LWRadSky_col(:,:)                  !sky longwave radiation , [MJ/h]
  real(r8),target,allocatable ::  TRAD_col(:,:)                          !total daily solar radiation, [MJ d-1]
  real(r8),target,allocatable ::  HUDX_col(:,:)                          !daily maximum vapor pressure , [kPa]
  real(r8),target,allocatable ::  HUDN_col(:,:)                          !daily minimum vapor pressure , [kPa]
  real(r8),target,allocatable ::  TWIND_col(:,:)                         !total daily wind travel, [m d-1]
!  real(r8),target,allocatable ::  PrecDaily_col(:,:)                !total daily precipitation, [m d-1]
  real(r8),target,allocatable ::  SkyLonwRad_col(:,:)                !sky longwave radiation , [MJ m-2 h-1]
  real(r8),target,allocatable ::  TempOffset_col(:,:)                !TempOffset_col for calculating temperature in Arrhenius curves, [oC]
  real(r8),target,allocatable ::  PrecDirect2Grnd_col(:,:)                     !direct precipitation at ground surface used to calculate soil erosion, [m h-1]
  real(r8),target,allocatable ::  PrecIndirect2Grnd_col(:,:)                         !indirect precipitation at ground surface used to calculate soil erosion, [m h-1]
  real(r8),target,allocatable ::  CO2EI_col(:,:)                         !initial atmospheric CO2 concentration, [umol mol-1]
  real(r8),target,allocatable ::  CCO2EI_gperm3_col(:,:)                        !initial atmospheric CO2 concentration, [gC m-3]

  real(r8),target,allocatable ::  AtmGasCgperm3_col(:,:,:)               !atmospheric gas concentration, [g m-3]
  real(r8),target,allocatable ::  AtmGmms_col(:,:,:)                     !atmospheric gas concentration, [umol mol-1]
  real(r8),target,allocatable ::  OXYE_col(:,:)                          !atmospheric O2 concentration, [umol mol-1]
  real(r8),target,allocatable ::  Z2OE_col(:,:)                          !atmospheric N2O concentration, [umol mol-1]
  real(r8),target,allocatable ::  Z2GE_col(:,:)                      !atmospheric N2 concentration, [umol mol-1]
  real(r8),target,allocatable ::  ZNH3E_col(:,:)                     !atmospheric NH3 concentration, [umol mol-1]
  real(r8),target,allocatable ::  CH4E_col(:,:)                      !atmospheric CH4 concentration, [umol mol-1]
  real(r8),target,allocatable ::  H2GE_col(:,:)                      !atmospheric H2 concentration, [umol mol-1]
  real(r8),target,allocatable ::  CO2E_col(:,:)                      !atmospheric CO2 concentration, [umol mol-1]
  real(r8),target,allocatable ::  ARGE_col(:,:)                      !atmospheric AR concentration, [umol mol-1]
  real(r8),target,allocatable ::  SolarNoonHourYM_col(:,:)           !year mean time of solar noon, [h]
  real(r8),target,allocatable ::  SolarNoonHour_col(:,:)             !time of solar noon each day, [h]
  real(r8),target,allocatable ::  RadSWDirect_col(:,:)               !direct shortwave radiation, [W m-2]
  real(r8),target,allocatable ::  RadSWDiffus_col(:,:)               !diffuse shortwave radiation, [W m-2]
  real(r8),target,allocatable ::  RadDirectPAR_col(:,:)              !direct PAR, [umol m-2 s-1]
  real(r8),target,allocatable ::  RadPARDiffus_col(:,:)              !diffuse PAR, [umol m-2 s-1]
  real(r8),target,allocatable ::  SineSunInclinationAngle_col(:,:)              !sine of solar angle, [-]
  real(r8),target,allocatable ::  SineSunInclAnglNxtHour_col(:,:)        !sine of solar angle next hour, [-]
  real(r8),target,allocatable ::  TLEX_col(:,:)                          !total latent heat flux x boundary layer resistance, [MJ m-1]
  real(r8),target,allocatable ::  TSHX_col(:,:)                          !total sensible heat flux x boundary layer resistance, [MJ m-1]
  real(r8),target,allocatable ::  Air_Heat_Latent_store_col(:,:)         !total latent heat flux x boundary layer resistance, [MJ m-1]
  real(r8),target,allocatable ::  Air_Heat_Sens_store_col(:,:)           !total sensible heat flux x boundary layer resistance, [MJ m-1]
  real(r8),target,allocatable ::  dSensHeatdTKQ_col(:,:)                 !hourly resistance estimate of weighted sensible-heat derivative [MJ m-1 K-1]
  real(r8),target,allocatable ::  TKQHeatRef_col(:,:)                     !Tq used for last completed hour [K]; zero means no completed-hour reference
  real(r8),target,allocatable ::  dLatHeatdVPQ_col(:,:)                   !weighted latent-energy derivative w.r.t. vapor concentration [MJ m2 ton-1]
  real(r8),target,allocatable ::  VPQLatHeatRef_col(:,:)                  !last completed-hour vapor concentration [ton m-3]; -1 means no reference
  real(r8),target,allocatable ::  SoilHeatSrcDepth_col(:,:)              !depth of soil heat sink/source, [m]
  real(r8),target,allocatable ::  TKSD_col(:,:)                          !temperature of soil heat sink/source, [oC]
  real(r8),target,allocatable ::  ATCAI_col(:,:)                         !initial mean annual air temperature, [oC]
  real(r8),target,allocatable ::  RadSWSolarBeam_col(:,:)                !shortwave radiation in solar beam, [MJ m-2 h-1]
  real(r8),target,allocatable ::  RadPARSolarBeam_col(:,:)               !PAR radiation in solar beam, [umol m-2 s-1]
  real(r8),target,allocatable ::  ATCA_col(:,:)                          !mean annual air temperature, [oC]
  real(r8),target,allocatable ::  ATCS_col(:,:)                          !mean annual soil temperature, [oC]
  real(r8),target,allocatable ::  TairKClimMean_col(:,:)                 !mean annual air temperature, [K]
  real(r8),target,allocatable ::  ATKS_col(:,:)                          !mean annual soil temperature, [K]
  real(r8),target,allocatable ::  RainFalPrec_col(:,:)                   !rainfall, [m3 d-2 h-1]
  real(r8),target,allocatable ::  SnoFalPrec_col(:,:)                    !snowfall, [m3 d-2 h-1]
  real(r8),target,allocatable ::  Irrigation_col(:,:)                    !irrigation, [m3 d-2 h-1]
  real(r8),target,allocatable ::  PrecAtm_col(:,:)                       !rainfall + snowfall, [m3 d-2 h-1]
  real(r8),target,allocatable ::  PrecRainAndIrrig_col(:,:)              !rainfall + irrigation, [m3 d-2 h-1]
  real(r8),target,allocatable ::  EnergyImpact4Erosion_col(:,:)          !cumulative rainfall energy impact on soil surface,[MJ d-2]
  real(r8),target,allocatable ::  pH_rain_col(:,:)                       !precipitation pH, [-]
  real(r8),target,allocatable ::  CN4RI_col(:,:)                         !precipitation initial NH4 concentration, [mol N m-3]
  real(r8),target,allocatable ::  CNORI_col(:,:)                         !precipitation initial NO3 concentration, [mol N m-3]
  real(r8),target,allocatable ::  NH4_rain_mole_conc(:,:)                !precipitation  NH4 concentration, [mol m-3]
  real(r8),target,allocatable ::  NO3_rain_mole_conc(:,:)                !precipitation  NO3 concentration, [mol m-3]
  real(r8),target,allocatable ::  H2PO4_rain_mole_conc(:,:)              !precipitation  H2PO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  trcg_rain_mole_conc_col(:,:,:)         !precipitation volatile concentration, [mol m-3]
  real(r8),target,allocatable ::  HPO4_rain_mole_conc_col(:,:)           !precipitation  HPO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  GDD_col(:,:)                           !growing degree day with base temperature at oC
  real(r8),target,allocatable ::  PrecHeat_col(:,:)                      !precipitation heat to surface, [MJ/d2/h]
  real(r8),target,allocatable ::  RainLitr_col(:,:)                      !water from aboveground falling litter,[m3 d-2]
  real(r8),target, allocatable :: trcs_solcoef_col(:,:,:)                !parameter for computing RGasSinkScalar_vr,[-]  
  real(r8),target,allocatable :: sp_popdens_pft(:,:,:) !prescribed annual population density, [plants m-2]
  real(r8),target,allocatable ::  tlai_mon_pft(:,:,:,:)                  !monthly lai for different pft used in prescribed phenology, [m2/m2]
  real(r8),target,allocatable ::  tsai_mon_pft(:,:,:,:)                  !monthly sai for different pft used in prescribed phenology, [m2/m2]
  real(r8),target,allocatable ::  tlai_day_pft(:,:,:)                    !daily prescribed leaf area for each PFT, [m2 cell-1]
  real(r8),target,allocatable ::  tsai_day_pft(:,:,:)                    !daily prescribed stem area for each PFT, [m2 cell-1]
  real(r8),target,allocatable ::  height_top_mon_pft(:,:,:,:)            !monthly canopy top height used in prescribed phenology, [m]
  real(r8),target,allocatable ::  height_bot_mon_pft(:,:,:,:)            !monthly canopy bottom height used in prescribed phenology, [m]
  REAL(R8),target,allocatable ::  NWetDep_col(:,:)                       !atmospheric inorganic N wet deposition, [gN m-2 h-1]
  real(r8),target,allocatable :: f_aerosol_DOM_col(:,:)                  !fraction of bioaerosol C as dissolved organic matter [0-1]
  real(r8),target,allocatable :: f_aerosol_DeadMB_col(:,:)               !fraction of bioaerosol C as dead microbial biomass [0-1]
  real(r8),target,allocatable :: f_aerosol_LiveMB_col(:,:)               !fraction of bioaerosol C as live microbial biomass [0-1]
  real(r8),target,allocatable :: f_aerosol_DeadNMB_col(:,:)              !fraction of bioaerosol C as dead nonmicrobial biomass [0-1]
  real(r8),target,allocatable :: f_aerosol_LichB_col(:,:)                !fraction of bioaerosol C as lichen spore, [0-1]
  real(r8),target,allocatable :: f_aerosol_MossB_col(:,:)                !fraction of bioaerosol C as moss spore, [0-1]
  real(r8),target,allocatable :: fAeroScalar_col(:,:)                    !optional per-topounit bioaerosol carbon deposition scalar read from grid file, [-]
  real(r8),target,allocatable :: SeedCDeposition_pft(:,:,:)              !seed C deposition for each PFT, [gC m-2]

  contains
!----------------------------------------------------------------------

  subroutine InitClimForcData()
  use TracerIDMod
  implicit none

  if(len(trim(warming_exp))>10)then
    allocate(TKS_ref_vr(8784,JZ,JY,JX));TKS_ref_vr=0._r8
  endif
  allocate(f_aerosol_DOM_col(JY,JX));f_aerosol_DOM_col=0._r8
  allocate(f_aerosol_DeadMB_col(JY,JX));f_aerosol_DeadMB_col=0._r8
  allocate(f_aerosol_LiveMB_col(JY,JX));f_aerosol_LiveMB_col=0._r8
  allocate(f_aerosol_DeadNMB_col(JY,JX));f_aerosol_DeadNMB_col=0._r8
  allocate(f_aerosol_LichB_col(JY,JX));f_aerosol_LichB_col=0._r8
  allocate(f_aerosol_MossB_col(JY,JX));f_aerosol_MossB_col=0._r8
  allocate(fAeroScalar_col(JY,JX));fAeroScalar_col=0._r8
  allocate(SeedCDeposition_pft(JP,JY,JX));SeedCDeposition_pft=0._r8
  allocate(NWetDep_col(JY,JX)); NWetDep_col=0._r8
  allocate(sp_popdens_pft(JP,JY,JX));sp_popdens_pft=0._r8
  allocate(tlai_mon_pft(12,JP,JY,JX));tlai_mon_pft=0._r8
  allocate(tsai_mon_pft(12,JP,JY,JX));tsai_mon_pft=0._r8
  allocate(tlai_day_pft(JP,JY,JX));tlai_day_pft=0._r8
  allocate(tsai_day_pft(JP,JY,JX));tsai_day_pft=0._r8
  allocate(height_top_mon_pft(12,JP,JY,JX)); height_top_mon_pft=0._r8
  allocate(height_bot_mon_pft(12,JP,JY,JX)); height_bot_mon_pft=0._r8

  allocate(EMS_Modify_Scalar_col(JY,JX)); EMS_Modify_Scalar_col = 1._r8
  allocate(srad_scalar_col(JY,JX));srad_scalar_col=1._r8
  allocate(trcs_solcoef_col(idg_beg:idg_NH3,JY,JX)); trcs_solcoef_col=0._r8
  allocate(Eco_RadSW_col(JY,JX)); Eco_RadSW_col=0._r8
  allocate(RadSWATS_col(JY,JX));  RadSWATS_col=0._r8
  allocate(GDD_col(JY,JX)); GDD_col=0._r8
  allocate(TDTPX(12,JY,JX));    TDTPX=0._r8
  allocate(TDTPN(12,JY,JX));    TDTPN=0._r8
  allocate(TDRAD(12,JY,JX));    TDRAD=0._r8
  allocate(TDHUM(12,JY,JX));    TDHUM=0._r8
  allocate(TDPRC(12,JY,JX));    TDPRC=0._r8
  allocate(TDWND(12,JY,JX));    TDWND=0._r8
  allocate(TDCN4(12,JY,JX));    TDCN4=0._r8
  allocate(TDCNO(12,JY,JX));    TDCNO=0._r8
  allocate(PrecHeat_col(JY,JX)); PrecHeat_col=0._r8
  allocate(RainLitr_col(JY,JX)); RainLitr_col=0._r8
  allocate(TCA_col(JY,JX));         TCA_col=0._r8
  allocate(TairK_col(JY,JX));         TairK_col=0._r8
  allocate(WindSpeedAtm_col(JY,JX));          WindSpeedAtm_col=0._r8
  allocate(VPA_col(JY,JX));         VPA_col=0._r8
  allocate(VPK_col(JY,JX));         VPK_col=0._r8
  allocate(PBOT_col(JY,JX));        PBOT_col=1.01325E+02_r8
  allocate(DayLensCurr_col(JY,JX));        DayLensCurr_col=0._r8
  allocate(DayLenthPrev_col(JY,JX));        DayLenthPrev_col=0._r8
  allocate(DayLenthMax_col(JY,JX));        DayLenthMax_col=0._r8
  allocate(OMEGA2Ground(NumOfSkyAzimuthSects,JY,JX));  OMEGA2Ground=0._r8
  allocate(LWRadSky_col(JY,JX));         LWRadSky_col=0._r8
  allocate(TRAD_col(JY,JX));        TRAD_col=0._r8
  allocate(HUDX_col(JY,JX));        HUDX_col=0._r8
  allocate(HUDN_col(JY,JX));        HUDN_col=0._r8
  allocate(TWIND_col(JY,JX));       TWIND_col=0._r8
!  allocate(PrecDaily_col(JY,JX));        PrecDaily_col=0._r8
  allocate(SkyLonwRad_col(JY,JX));        SkyLonwRad_col=0._r8
  allocate(TempOffset_col(JY,JX));      TempOffset_col=0._r8
  allocate(PrecDirect2Grnd_col(JY,JX));       PrecDirect2Grnd_col=0._r8
  allocate(PrecIndirect2Grnd_col(JY,JX));       PrecIndirect2Grnd_col=0._r8
  allocate(CO2EI_col(JY,JX));       CO2EI_col=0._r8
  allocate(CCO2EI_gperm3_col(JY,JX));      CCO2EI_gperm3_col=0._r8

  allocate(AtmGasCgperm3_col(idg_beg:idg_end,JY,JX)); AtmGasCgperm3_col=0._r8
  allocate(AtmGmms_col(idg_beg:idg_end,JY,JX)); AtmGmms_col=0._r8

  allocate(OXYE_col(JY,JX));        OXYE_col=0._r8
  allocate(Z2OE_col(JY,JX));        Z2OE_col=0._r8
  allocate(Z2GE_col(JY,JX));        Z2GE_col=0._r8
  allocate(ZNH3E_col(JY,JX));       ZNH3E_col=0._r8
  allocate(CH4E_col(JY,JX));        CH4E_col=0._r8
  allocate(H2GE_col(JY,JX));        H2GE_col=0._r8
  allocate(CO2E_col(JY,JX));        CO2E_col=0._r8
  allocate(ARGE_col(JY,JX));        ARGE_col=0._r8
  allocate(SolarNoonHourYM_col(JY,JX));       SolarNoonHourYM_col=0._r8
  allocate(SolarNoonHour_col(JY,JX)); SolarNoonHour_col=0._r8
  allocate(RadSWDirect_col(JY,JX));        RadSWDirect_col=0._r8
  allocate(RadSWDiffus_col(JY,JX));        RadSWDiffus_col=0._r8
  allocate(RadDirectPAR_col(JY,JX));        RadDirectPAR_col=0._r8
  allocate(RadPARDiffus_col(JY,JX));        RadPARDiffus_col=0._r8
  allocate(SineSunInclinationAngle_col(JY,JX));        SineSunInclinationAngle_col=0._r8
  allocate(SineSunInclAnglNxtHour_col(JY,JX));       SineSunInclAnglNxtHour_col=0._r8
  allocate(TLEX_col(JY,JX));        TLEX_col=0._r8
  allocate(TSHX_col(JY,JX));        TSHX_col=0._r8
  allocate(Air_Heat_Latent_store_col(JY,JX));        Air_Heat_Latent_store_col=0._r8
  allocate(Air_Heat_Sens_store_col(JY,JX));        Air_Heat_Sens_store_col=0._r8
  allocate(dSensHeatdTKQ_col(JY,JX)); dSensHeatdTKQ_col=0._r8
  allocate(TKQHeatRef_col(JY,JX));    TKQHeatRef_col=0._r8
  allocate(dLatHeatdVPQ_col(JY,JX));  dLatHeatdVPQ_col=0._r8
  allocate(VPQLatHeatRef_col(JY,JX)); VPQLatHeatRef_col=-1._r8
  allocate(SoilHeatSrcDepth_col(JY,JX));      SoilHeatSrcDepth_col=0._r8
  allocate(TKSD_col(JY,JX));        TKSD_col=0._r8
  allocate(ATCAI_col(JY,JX));       ATCAI_col=0._r8
  allocate(RadSWSolarBeam_col(JY,JX));         RadSWSolarBeam_col=0._r8
  allocate(RadPARSolarBeam_col(JY,JX));         RadPARSolarBeam_col=0._r8
  allocate(ATCA_col(JY,JX));        ATCA_col=0._r8
  allocate(ATCS_col(JY,JX));        ATCS_col=0._r8
  allocate(TairKClimMean_col(JY,JX));        TairKClimMean_col=0._r8
  allocate(ATKS_col(JY,JX));        ATKS_col=0._r8
  allocate(RainFalPrec_col(JY,JX));       RainFalPrec_col=0._r8
  allocate(Irrigation_col(JY,JX));    Irrigation_col=0._r8
  allocate(SnoFalPrec_col(JY,JX));       SnoFalPrec_col=0._r8
  allocate(PrecAtm_col(JY,JX));       PrecAtm_col=0._r8
  allocate(PrecRainAndIrrig_col(JY,JX));       PrecRainAndIrrig_col=0._r8
  allocate(EnergyImpact4Erosion_col(JY,JX));       EnergyImpact4Erosion_col=0._r8
  allocate(pH_rain_col(JY,JX));         pH_rain_col=0._r8
  allocate(CN4RI_col(JY,JX));       CN4RI_col=0._r8
  allocate(CNORI_col(JY,JX));       CNORI_col=0._r8
  allocate(NH4_rain_mole_conc(JY,JX));        NH4_rain_mole_conc=0._r8
  allocate(NO3_rain_mole_conc(JY,JX));        NO3_rain_mole_conc=0._r8
  allocate(H2PO4_rain_mole_conc(JY,JX));        H2PO4_rain_mole_conc=0._r8
  allocate(trcg_rain_mole_conc_col(idg_beg:idg_NH3,JY,JX));        trcg_rain_mole_conc_col=0._r8
  allocate(HPO4_rain_mole_conc_col(JY,JX));       HPO4_rain_mole_conc_col=0._r8

  end subroutine InitClimForcData

!----------------------------------------------------------------------
  subroutine DestructClimForcData
  use abortutils, only : destroy
  implicit none

  call destroy(f_aerosol_DOM_col)
  call destroy(f_aerosol_DeadMB_col)
  call destroy(f_aerosol_LiveMB_col)
  call destroy(f_aerosol_DeadNMB_col)
  call destroy(f_aerosol_MossB_col)
  call destroy(f_aerosol_LichB_col)
  call destroy(fAeroScalar_col)
  call destroy(SeedCDeposition_pft)
  call destroy(NWetDep_col)
  call destroy(height_bot_mon_pft)  
  call destroy(height_top_mon_pft)
  call destroy(sp_popdens_pft)
  call destroy(tlai_mon_pft)
  call destroy(tsai_mon_pft)
  call destroy(tlai_day_pft)
  call destroy(tsai_day_pft)
  call destroy(PrecHeat_col)
  call destroy(RainLitr_col)
  call destroy(TDTPX)
  call destroy(TDTPN)
  call destroy(TDRAD)
  call destroy(TDHUM)
  call destroy(TDPRC)
  call destroy(TDWND)
  call destroy(TDCN4)
  call destroy(TDCNO)
  call destroy(TCA_col)
  call destroy(TairK_col)
  call destroy(WindSpeedAtm_col)
  call destroy(VPA_col)
  call destroy(VPK_col)
  call destroy(PBOT_col)
  call destroy(DayLensCurr_col)
  call destroy(DayLenthPrev_col)
  call destroy(DayLenthMax_col)
  call destroy(OMEGA2Ground)
  call destroy(LWRadSky_col)
  call destroy(TRAD_col)
  call destroy(HUDX_col)
  call destroy(HUDN_col)
  call destroy(TWIND_col)
  call destroy(ARGE_col)
!  call destroy(PrecDaily_col)
  call destroy(SkyLonwRad_col)
  call destroy(TempOffset_col)
  call destroy(PrecDirect2Grnd_col)
  call destroy(PrecIndirect2Grnd_col)
  call destroy(CO2EI_col)
  call destroy(CCO2EI_gperm3_col)
  call destroy(srad_scalar_col)
  call destroy(EMS_Modify_Scalar_col)
  call destroy(AtmGasCgperm3_col)
  call destroy(AtmGmms_col)
  call destroy(TKS_ref_vr)
  call destroy(Eco_RadSW_col)
  call destroy(RadSWATS_col)
  call destroy(OXYE_col)
  call destroy(Z2OE_col)
  call destroy(Z2GE_col)
  call destroy(ZNH3E_col)
  call destroy(CH4E_col)
  call destroy(H2GE_col)
  call destroy(trcs_solcoef_col)
  call destroy(SolarNoonHourYM_col)
  call destroy(SolarNoonHour_col)
  call destroy(CO2E_col)
  call destroy(RadSWDirect_col)
  call destroy(RadSWDiffus_col)
  call destroy(RadDirectPAR_col)
  call destroy(RadPARDiffus_col)
  call destroy(SineSunInclinationAngle_col)
  call destroy(SineSunInclAnglNxtHour_col)
  call destroy(TLEX_col)
  call destroy(TSHX_col)
  call destroy(Air_Heat_Latent_store_col)
  call destroy(Air_Heat_Sens_store_col)
  call destroy(dSensHeatdTKQ_col)
  call destroy(TKQHeatRef_col)
  call destroy(dLatHeatdVPQ_col)
  call destroy(VPQLatHeatRef_col)
  call destroy(SoilHeatSrcDepth_col)
  call destroy(TKSD_col)
  call destroy(ATCAI_col)
  call destroy(RadSWSolarBeam_col)
  call destroy(RadPARSolarBeam_col)
  call destroy(ATCA_col)
  call destroy(ATCS_col)
  call destroy(TairKClimMean_col)
  call destroy(ATKS_col)
  call destroy(RainFalPrec_col)
  call destroy(Irrigation_col)
  call destroy(SnoFalPrec_col)
  call destroy(PrecAtm_col)
  call destroy(PrecRainAndIrrig_col)
  call destroy(EnergyImpact4Erosion_col)
  call destroy(pH_rain_col)
  call destroy(CN4RI_col)
  call destroy(CNORI_col)
  call destroy(NH4_rain_mole_conc)
  call destroy(NO3_rain_mole_conc)
  call destroy(H2PO4_rain_mole_conc)
  call destroy(trcg_rain_mole_conc_col)
  call destroy(HPO4_rain_mole_conc_col)
  end subroutine DestructClimForcData
end module ClimForcDataType
