module WthrMod
  !
  ! Description:
  ! code to process the weather forcing
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use MiniMathMod,       only: safe_adb, vapsat0, isclose
  use MiniFuncMod,       only: get_sun_declin,f_microbe_TOffset
  use EcoSIMCtrlMod,     only: etimer, frectyp,fixClime, ats_cpl_mode
  use PlantMgmtDataType, only: NP_col
  use MiniMathMod,       only: AZMAX1
  use UnitMod,           only: units
  use EcosysWarmingMod
  use DebugToolMod  
  use EcosimConst
  use CanopyRadDataType
  use GridConsts
  use FlagDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use LandSurfDataType
  use FertilizerDataType
  use PlantTraitDataType
  use AqueChemDatatype
  use IrrigationDataType
  use GridDataType
  use SoilWaterDataType
  use EcoSIMConfig
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  real(r8) :: AMP,CLD,DTA,DHR,DTS,EMM,RADZ,VPX,XJ

  !
  !     CDIR,CDIF=fraction of solar SW,sky diffuse radiation in visible
  !     PDIR,PDIF=PAR:SW ratio (umol s-1/(MJ h-1))
  !     TSNOW=temperature below which precipitation is snow (oC)
  !
  real(r8), parameter :: SolConst=4.896_r8 !MJ/m2/hour, solar constant ~ 1360 W/m2
  real(r8), PARAMETER :: CDIR=0.42_r8   !
  real(r8), parameter :: CDIF=0.58_r8
  real(r8), parameter :: PDIR=1269.4_r8
  real(r8), parameter :: PDIF=1269.4_r8
  real(r8), PARAMETER :: TSNOW=-0.25_r8  !oC, threshold temperature for snowfall

  public :: PrepHourlyWeather
  contains

  SUBROUTINe PrepHourlyWeather(I,J,NHW,NHE,NVN,NVS)
  !
  !     Description:
  !
  !     REINITIALIZES WEATHER VARIABLES USED IN OTHER
  !     SUBROUTINES
  !
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  character(len=*), parameter :: subname='PrepHourlyWeather'
  integer :: ITYPE,NX,NY,N,NZ,mon
  real(r8) :: PrecAsRain_col(JY,JX)
  real(r8) :: PrecAsSnow_col(JY,JX)
  real(r8) :: PRECII_col(JY,JX)   !surface irrigation
  real(r8) :: PRECUI_col(JY,JX)   !subsurface irrigation
  real(r8) :: RADN_col(JY,JX)     ![MJ/m2/h]
  real(r8) :: VPS(JY,JX)

  !     execution begins here
  call PrintInfo('beg '//subname)
  XJ=J
  DOY=I-1+XJ/24
  !
  !     SWITCH OUT ECOSYS WEATHER HERE IF CESM WEATHER IS READ IN
  !
  !     IF
  !
  !     IWTHR=weather data type in first(1) or second(2) scene
  !     ITYPE=weather data type:1=daily,2=hourly
  !

  ! Restore neutral modifiers for every forcing path before applying this
  ! hour's eligible OTC treatment, including after its season or dates end.
  EMS_Modify_Scalar_col(NVN:NVS,NHW:NHE)=1._r8
  srad_scalar_col(NVN:NVS,NHW:NHE)=1._r8
  ITYPE=IWTHR

  !
  !     CALCULATE HOURLY TEMPERATURE, RADIATION, WINDSPEED, VAPOR PRESSURE
  !     AND PRECIPITATION FROM DAILY WEATHER ARRAYS LOADED IN 'READS'
  !
  IF(ITYPE.EQ.1)THEN
!
    call DailyWeather(I,J,NHW,NHE,NVN,NVS,RADN_col,PrecAsRain_col,PrecAsSnow_col,VPS)
    !     CALCULATE HOURLY TEMPERATURE, RADIATION, WINDSPEED, VAPOR PRESSURE
    !     AND PRECIPITATION FROM HOURLY WEATHER ARRAYS LOADED IN 'READS'
    !
  else if (IWTHR.EQ.-999)then
    !ATS coupling: per-column incident shortwave from ATS. Precipitation is set
    !by the ATS driver after this call, so it is zeroed here.
    RADN_col       = RadSWATS_col
    PrecAsRain_col = 0._r8
    PrecAsSnow_col = 0._r8
  ELSE
    call HourlyWeather(I,J,NHW,NHE,NVN,NVS,RADN_col,PrecAsRain_col,PrecAsSnow_col,VPS)
  ENDIF

  ! Current geometry controls day/night eligibility; radiation below must
  ! use the eligible treatment's transmission and atmospheric temperature.
  call UpdateSolarGeometry(I,J,NHW,NHE,NVN,NVS,RADN_col)

  if(check_warming_dates(iYearCurrent,I,J))then    
    call apply_OTC_warming(I,J,NHW,NHE,NVN,NVS)
    call apply_IR_warming(I,J,NHW,NHE,NVN,NVS) 
  endif  

  call CalcRadiation(I,J,NHW,NHE,NVN,NVS,RADN_col,PRECUI_col,PRECII_col)
!
!  IF(ICLM.EQ.1.OR.ICLM.EQ.2)THEN
!    !perturb climate data
!    call CorrectClimate(I,J,NHW,NHE,NVN,NVS,PRECUI_col,PrecAsRain_col,PRECII_col,PrecAsSnow_col,VPS)
!  ENDIF
!
  mon=etimer%get_curr_mon()

  call SummaryClimateForc(I,J,NHW,NHE,NVN,NVS,PRECUI_col,PrecAsRain_col,PRECII_col,PrecAsSnow_col)
  call PrintInfo('end '//subname)
  END subroutine PrepHourlyWeather
!------------------------------------------------------------------------------------------

  subroutine DailyWeather(I,J,NHW,NHE,NVN,NVS,RADN_col,PrecAsRain_col,PrecAsSnow_col,VPS)
!
!     Description:
!
  implicit none
  integer,  intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8), intent(out) :: RADN_col(JY,JX),PrecAsRain_col(JY,JX)
  real(r8), intent(out) :: PrecAsSnow_col(JY,JX),VPS(JY,JX)
  integer :: NY,NX
  !     begin_execution

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      !
      !     KoppenClimZone=Koppen climate zone:-2=phytotron
      !     RADN_col=hourky SW radiation
      !     RMAX=maximum hourly radiation
      !     SolarNoonHour_col=time of solar noon
      !     DayLensCurr_col=current daylength
!
      IF(KoppenClimZone_col(NY,NX).NE.-2)THEN
!       not phytotron
        IF(DayLensCurr_col(NY,NX).GT.ZERO)THEN
          RADN_col(NY,NX)=AZMAX1(RMAX*SIN((J-(SolarNoonHour_col(NY,NX) &
            -DayLensCurr_col(NY,NX)/2.0_r8))*PICON/DayLensCurr_col(NY,NX)))
        ELSE
          RADN_col(NY,NX)=0.0_r8
        ENDIF
      ELSE
!       phytotron
        RADN_col(NY,NX)=RMAX/24.0_r8
      ENDIF
      !
      !     TCA,TairK=air temperature (oC,K)
      !     TAVG*,AMP*=daily averages, amplitudes from day.f
      !
      IF(J.LT.(SolarNoonHour_col(NY,NX)-DayLensCurr_col(NY,NX)*0.5_r8))THEN
        TCA_col(NY,NX)=TAVG1+AMP1*SIN(((J+SolarNoonHour_col(NY,NX)-3.0)*PICON &
          /(SolarNoonHour_col(NY,NX)+9.0-DayLensCurr_col(NY,NX)/2.0))+PICON2h)
      ELSEIF(J.GT.SolarNoonHour_col(NY,NX)+3)THEN
        TCA_col(NY,NX)=TAVG3+AMP3*SIN(((J-SolarNoonHour_col(NY,NX)-3.0)*PICON &
          /(SolarNoonHour_col(NY,NX)+9.0-DayLensCurr_col(NY,NX)/2.0))+PICON2h)
      ELSE
        TCA_col(NY,NX)=TAVG2+AMP2*SIN(((J-(SolarNoonHour_col(NY,NX) &
          -DayLensCurr_col(NY,NX)/2.0_r8))*PICON/(3.0_r8+DayLensCurr_col(NY,NX)/2.0_r8))-PICON2h)
      ENDIF
      TairK_col(NY,NX)=units%Celcius2Kelvin(TCA_col(NY,NX))
      if(abs(TairK_col(NY,NX))>400._r8)then
        print*,'air temperature problematic',TairK_col(NY,NX),TCA_col(NY,NX)
      endif
      !
      !     VPK,VPS=ambient,saturated vapor pressure
      !     VAVG*,VAMP*=daily averages, amplitudes from day.f
      !      ALTI=altitude
      !
      IF(J.LT.(SolarNoonHour_col(NY,NX)-DayLensCurr_col(NY,NX)/2))THEN
        VPK_col(NY,NX)=VAVG1+VMP1*SIN(((J+SolarNoonHour_col(NY,NX)-3.0_r8)*PICON &
          /(SolarNoonHour_col(NY,NX)+9.0_r8-DayLensCurr_col(NY,NX)/2.0_r8))+PICON2h)
      ELSEIF(J.GT.SolarNoonHour_col(NY,NX)+3)THEN
        VPK_col(NY,NX)=VAVG3+VMP3*SIN(((J-SolarNoonHour_col(NY,NX)-3.0_r8)*PICON &
          /(SolarNoonHour_col(NY,NX)+9.0_r8-DayLensCurr_col(NY,NX)/2.0_r8))+PICON2h)
      ELSE
        VPK_col(NY,NX)=VAVG2+VMP2*SIN(((J-(SolarNoonHour_col(NY,NX) &
          -DayLensCurr_col(NY,NX)/2.0_r8))*PICON /(3.0_r8+DayLensCurr_col(NY,NX)/2.0_r8))-PICON2h)
      ENDIF
      !VPS(NY,NX)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/TairK_col(NY,NX))) &
      VPS(NY,NX)     = vapsat0(TairK_col(ny,nx))*EXP(-AltitudeLandscape_col(NY,NX)/7272.0_r8)
      VPK_col(NY,NX) = AMIN1(VPS(NY,NX),VPK_col(NY,NX))
      PBOT_col(NY,NX)=1.01325E+02_r8*exp(-SurfElevation_col(NY,NX)/hpresc)
!
      !     UA=wind speed
!
      WindSpeedAtm_col(NY,NX)=AMAX1(3600.0_r8,WIND(I))
!
!     TSNOW=temperature below which precipitation is snow (oC)
!     PrecAsRain_col,PrecAsSnow_col=rainfall,snowfall
!
      IF(J.GE.13.AND.J.LE.24)THEN
        IF(TCA_col(NY,NX).GT.TSNOW)THEN
          PrecAsRain_col(NY,NX)=RAIN(I)/12.0_r8   !rainfall
          PrecAsSnow_col(NY,NX)=0.0_r8   !snow fall
        ELSE
          PrecAsRain_col(NY,NX)=0.0_r8
          PrecAsSnow_col(NY,NX)=RAIN(I)/12.0_r8   !precip as snowfall
    !     IF(PrecAsSnow_col(NY,NX).LT.0.25E-03_r8)PrecAsSnow_col(NY,NX)=0.0_r8
        ENDIF
      ELSE
        PrecAsRain_col(NY,NX)=0.0_r8
        PrecAsSnow_col(NY,NX)=0.0_r8
      ENDIF
    ENDDO
  enddo
!
  end subroutine DailyWeather
!------------------------------------------------------------------------------------------

  subroutine HourlyWeather(I,J,NHW,NHE,NVN,NVS,RADN_col,PrecAsRain_col,PrecAsSnow_col,VPS)
!
!     Description:
!
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8), intent(out) :: RADN_col(JY,JX)
  real(r8), intent(out) :: PrecAsRain_col(JY,JX),PrecAsSnow_col(JY,JX)
  real(r8), intent(out) :: VPS(JY,JX)
  integer :: NY,NX
  !     begin_execution

  DO  NX=NHW,NHE
    DO NY=NVN,NVS
 !
      !     RADN_col=SW radiation at horizontal surface,  [MJ/m2/hour]
      !     TCA,TairK=air temperature (oC,K)
      !     VPK,VPS=ambient,saturated vapor pressure, KPa
      !     UA=wind speed,  m/hr
      !     TSNOW=temperature below which precipitation is snow (oC)
      !     PrecAsRain_col,PrecAsSnow_col=rainfall,snowfall, m H2O /m2 /hr
!
      RADN_col(NY,NX)  = SWRad_hrly(J,I)
      TCA_col(NY,NX)   = TMP_hrly(J,I)
      TairK_col(NY,NX) = units%Celcius2Kelvin(TCA_col(NY,NX))

      !elevation corrected saturated air vapor pressure, KPa
      VPS(NY,NX)              = vapsat0(TairK_col(ny,nx))*EXP(-AltitudeLandscape_col(NY,NX)/7272.0_r8)
      VPK_col(NY,NX)          = AMIN1(DWPTH(J,I),VPS(NY,NX))
      WindSpeedAtm_col(NY,NX) = AMAX1(3600.0_r8,WINDH(J,I))
      PBOT_col(NY,NX)         = PBOT_hrly(J,I)
      !snowfall is determined by air tempeature
      IF(TCA_col(NY,NX).GT.TSNOW)THEN
        PrecAsRain_col(NY,NX) = RAINH(J,I)
        PrecAsSnow_col(NY,NX) = 0.0_r8
      ELSE
        PrecAsRain_col(NY,NX) = 0.0_r8
        PrecAsSnow_col(NY,NX) = RAINH(J,I)
      ENDIF
    enddo
  enddo

  end subroutine HourlyWeather
!------------------------------------------------------------------------------------------

  subroutine UpdateSolarGeometry(I,J,NHW,NHE,NVN,NVS,RADN_col)
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8), intent(inout) :: RADN_col(JY,JX) !Incident horizontal shortwave; cleared when geometry excludes sunlight
  character(len=*), parameter :: subname='UpdateSolarGeometry'
  integer :: NY,NX
  real(r8) :: AZI,DEC,DECLIN
  real(r8), parameter :: PICON12=PICON/12._r8

  call PrintInfo('beg '//subname)
  DECLIN=get_sun_declin(I)
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      IF(KoppenClimZone_col(NY,NX).GE.-1)THEN
        AZI=SIN(ALAT_col(NY,NX)*RadianPerDegree)*SIN(DECLIN*RadianPerDegree)
        DEC=COS(ALAT_col(NY,NX)*RadianPerDegree)*COS(DECLIN*RadianPerDegree)
        !check eq.(11.1) in Campbell and Norman, 1998, p168.
        if(fixClime)then
          !always assume the light is from zenith when using fixed climate forcing.
          SineSunInclinationAngle_col(NY,NX)       = AZMAX1(AZI+DEC*COS(PICON12*(SolarNoonHour_col(NY,NX)-11.5_r8)))
          SineSunInclAnglNxtHour_col(NY,NX) = AZMAX1(AZI+DEC*COS(PICON12*(SolarNoonHour_col(NY,NX)-12.5_r8)))
        else
          SineSunInclinationAngle_col(NY,NX)       = AZMAX1(AZI+DEC*COS(PICON12*(SolarNoonHour_col(NY,NX)-(J-0.5_r8))))
          SineSunInclAnglNxtHour_col(NY,NX) = AZMAX1(AZI+DEC*COS(PICON12*(SolarNoonHour_col(NY,NX)-(J+0.5_r8))))
        endif

        !IF(SineSunInclinationAngle_col(NY,NX).GT.0.0_r8 .AND. SineSunInclinationAngle_col(NY,NX).LT.TWILGT)SineSunInclinationAngle_col(NY,NX)=TWILGT
        IF(RADN_col(NY,NX).LE.0.0_r8)SineSunInclinationAngle_col(NY,NX)=0.0_r8
        IF(SineSunInclinationAngle_col(NY,NX).LE.-TWILGT)RADN_col(NY,NX)=0.0_r8

      ELSE
        IF(RADN_col(NY,NX).LE.0.0_r8)THEN
          SineSunInclinationAngle_col(NY,NX)=0.0_r8
        ELSE
          SineSunInclinationAngle_col(NY,NX)=1.0_r8
        ENDIF
        SineSunInclAnglNxtHour_col(NY,NX)=1.0_r8
      ENDIF
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine UpdateSolarGeometry
!------------------------------------------------------------------------------------------
  subroutine CalcRadiation(I,J,NHW,NHE,NVN,NVS,RADN_col,PRECUI_col,PRECII_col)
!
  !     Description:
!
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8), intent(inout) :: RADN_col(JY,JX)   !Incident shortwave radiation at surface horizontal plane, [MJ/hr/m2]
  real(r8), intent(out) :: PRECUI_col(JY,JX)
  real(r8), intent(out) :: PRECII_col(JY,JX)
  integer :: NY,NX
  real(r8) :: RADX    !solar constant at horizontal surface
  !     begin_execution
  !     CALCULATE DIRECT, DIFFUSE AND LONGWAVE RADIATION FROM
  !     INCOMING RADIATION READ IN 'READS', SOLAR ANGLE, HUMIDITY,
  !     TEMPERATURE AND CLOUDINESS
  !

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
!
!     IF OUTDOORS
!
!     SineSunInclinationAngle_col,SineSunInclAnglNxtHour_col=sine solar angle of current,next hour
!     
!     IETYP: koppen climate zone
      IF(KoppenClimZone_col(NY,NX).GE.-1)THEN
        RADX            = SolConst*AZMAX1(SineSunInclinationAngle_col(NY,NX))
        RADN_col(NY,NX) = AMIN1(RADX,RADN_col(NY,NX))
        !
        !     DIRECT VS DIFFUSE RADIATION IN SOLAR OR SKY BEAMS
        !
        !     RADZ=diffuse radiation at horizontal surface
        !     RADS,RADY,RAPS,RadPARDiffus_col=direct,diffuse SW,PAR in solar beam
        !
        RADZ                    = AMIN1(RADN_col(NY,NX),0.5_r8*(RADX-RADN_col(NY,NX)))
        RadSWDirect_col(NY,NX)  = safe_adb(RADN_col(NY,NX)-RADZ,SineSunInclinationAngle_col(NY,NX))
        RadSWDirect_col(NY,NX)  = AMIN1(4.167_r8,RadSWDirect_col(NY,NX))*srad_scalar_col(NY,NX)
        RadSWDiffus_col(NY,NX)  = RADZ/TotSineSkyAngles_grd*srad_scalar_col(NY,NX)
        RadDirectPAR_col(NY,NX) = RadSWDirect_col(NY,NX)*CDIR*PDIR  !MJ/m2/hr
        RadPARDiffus_col(NY,NX) = RadSWDiffus_col(NY,NX)*CDIF*PDIF  !MJ/m2/hr
        !
        !     ATMOSPHERIC RADIATIVE PROPERTIES 
        !  Duarte et al., 2006, AFM 139:171, Assessing daytime downward longwave radiation estimates for clear and cloudy skies in Southern Brazil
        !     CLD=cloudiness factor for EMM
        !     EMM=sky emissivity
        !     VPK,TairK=vapor pressure,temperature
        !
        IF(RADX.GT.ZERO)THEN
          CLD=AMIN1(1.0_r8,AMAX1(0.2_r8,2.33_r8-3.33_r8*RADN_col(NY,NX)/RADX))
        ELSE
          CLD=0.2_r8
        ENDIF
        EMM=0.625_r8*AMAX1(1.0_r8,(1.0E+03_r8*VPK_col(NY,NX)/TairK_col(NY,NX))**0.131_r8)
        EMM=EMM*(1.0_r8+0.242_r8*CLD**0.583_r8)
        !
        !     IF PHYTOTRON
        !
      ELSE
        CLD=0.0_r8
        EMM=0.96
      ENDIF
!
      !     LONGWAVE RADIATION
      !
      !     RadLWClm=longwave radiation
      !     THSX=longwave radiation from weather file or calculated from
      !     atmospheric properties

      IF(RadLWClm(J,I).GT.0.0_r8)THEN
        ! Match legacy ecosys: observed incoming longwave replaces the
        ! sky-emissivity estimate. Adding both terms double counts longwave.
        SkyLonwRad_col(NY,NX)=RadLWClm(J,I)
      ELSE
        SkyLonwRad_col(NY,NX)=EMM*stefboltz_const*TairK_col(NY,NX)**4
      ENDIF

!
      !     INSERT CESM WEATHER HERE
      !
      !     ELSE
      !     RADS=DIRECT SW RADIATION (MJ M-2 H-1)
      !     RadSWDiffus_col=INDIRECT SW RADIATION (MJ M-2 H-1)
      !     RadDirectPAR_col=DIRECT PAR (UMOL M-2 S-1)
      !     RadPARDiffus_col=INDIRECT PAR (UMOL M-2 S-1)
      !     THSX=LW RADIATION (MJ M-2 H-1)
      !     TCA=AIR TEMPERATURE (C)
      !     TairK=AIR TEMPERATURE (K)
      !     VPK=VAPOR PRESSURE (KPA)
      !     UA=WINDSPEED (M H-1)
      !     PrecAsRain_col(NY,NX)=RAIN (M H2O H-1 m-2)
      !     PrecAsSnow_col(NY,NX)=SNOW (M H2O H-1 m-2)
      !     SineSunInclinationAngle_col=SOLAR ANGLE CURRENT HOUR (SINE)
      !     SineSunInclAnglNxtHour_col=SOLAR ANGLE NEXT HOUR (SINE)
      !     ENDIF
      !
      !     ADD IRRIGATION
      !
      !     PRECII_col,PRECUI=surface,subsurface irrigation
      !     RRIG=irrigation from soil management file in reads.f
      !
!      WDPTHD=WDPTH(I,NY,NX)+CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX)
!      IF(WDPTHD.LE.CumDepz2LayBottom_vr(NU_col(NY,NX),NY,NX))THEN
        PRECII_col(NY,NX) = RRIG(J,I,NY,NX)   !surface irrigation
        PRECUI_col(NY,NX) = 0.0_r8
!      ELSE
!        PRECII_col(NY,NX)=0.0_r8
!        PRECUI_col(NY,NX)=RRIG(J,I,NY,NX)   !subsurface irrigation
!      ENDIF
    ENDDO
  ENDDO
  end subroutine CalcRadiation
!------------------------------------------------------------------------------------------

  subroutine CorrectClimate(I,J,NHW,NHE,NVN,NVS,PRECUI_col,PrecAsRain_col,PRECII_col,PrecAsSnow_col,VPS)
  !
  !     DESCRIPTION:
  !     IMPLEMENT CLIMATE CHANGES READ IN 'READS' TO HOURLY TEMPERATURE,
  !     RADIATION, WINDSPEED,VAPOR PRESSURE, PRECIPITATION, IRRIGATION
  !     AND CO2
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8), intent(inout) :: PRECUI_col(JY,JX),PrecAsRain_col(JY,JX),PRECII_col(JY,JX)
  real(r8), intent(inout) :: PrecAsSnow_col(JY,JX),VPS(JY,JX)
  integer :: NY,NX,NZ,N
  !     begin_execution
  !
  !     ICLM=changes to weather data (0=none,1=step,2=transient)
  !     N=season or month
  !
  !
  !     SEASONAL CHANGES
  !
  IF(I.GT.334.OR.I.LE.59)THEN
    N=1
  ELSEIF(I.GT.59.AND.I.LE.151)THEN
    N=2
  ELSEIF(I.GT.151.AND.I.LE.243)THEN
    N=3
  ELSE
    N=4
  ENDIF
  !

  D9925: DO NX=NHW,NHE
    D9920: DO NY=NVN,NVS
!
      !     TEMPERATURE CHANGE
      !
      !     TDTPX,TDTPN=change in max,min temperature
      !     DTA,AMP,DHR=change in daily average,amplitude of air temperature
      !     DHR=diurnal effect on AMP
      !
      IF(.not.isclose(TDTPX(N,NY,NX),0.0_r8).OR.(.not.isclose(TDTPN(N,NY,NX),0.0_r8)))THEN
        DTA=0.5_r8*(TDTPX(N,NY,NX)+TDTPN(N,NY,NX))
        AMP=0.5_r8*(TDTPX(N,NY,NX)-TDTPN(N,NY,NX))
        DHR=SIN(0.2618_r8*(J-(SolarNoonHour_col(NY,NX)+3.0_r8))+PICON2h)
        TCA_col(NY,NX)=TCA_col(NY,NX)+DTA+AMP*DHR
        TairK_col(NY,NX)=units%Celcius2Kelvin(TCA_col(NY,NX))
!
        !     ACCLIMATION TO GRADUAL CLIMATE CHANGE
        !
        !     DTS=change in daily average soil temperature
        !     ATCA,ATCS=mean annual air,soil temperature
        !     TempOffset_col=shift in Arrhenius curve for MFT activity in nitro.f
        !     TempOffset_pft=shift in Arrhenius curve for PFT activity in uptake.f
        !     rPlantThermoAdaptZone_pft=PFT thermal adaptation zone
        !     HTC=high temperature threshold for grain number loss (oC)
        !     GROUPI,ShootNodeNumAtPlanting_pft=node number at floral initiation,planting (maturity group)
        !
        IF(ICLM.EQ.2 .AND. J.EQ.1)THEN
          DTS                   = 0.5_r8*DTA
          ATCA_col(NY,NX)       = ATCAI_col(NY,NX)+DTA
          ATCS_col(NY,NX)       = ATCAI_col(NY,NX)+DTS
          TempOffset_col(NY,NX) = f_microbe_TOffset(ATCS_col(NY,NX))
          DO NZ=1,NP_col(NY,NX)
            rPlantThermoAdaptZone_pft(NZ,NY,NX)=PlantInitThermoAdaptZone_pft(NZ,NY,NX)+0.30_r8/2.667_r8*DTA
            TempOffset_pft(NZ,NY,NX)=2.667*(2.5-rPlantThermoAdaptZone_pft(NZ,NY,NX))
            !     TC4LeafOut_pft(NZ,NY,NX)=TCZD-TempOffset_pft(NZ,NY,NX)
            !     TC4LeafOff_pft(NZ,NY,NX)=AMIN1(15.0,TC4LeafOut_pft(NZ,NY,NX)+TCXD)
            IF(iPlantPhotosynsType_pft(NZ,NY,NX).EQ.3)THEN
              HighTempLimitSeed_pft(NZ,NY,NX)=27.0+3.0*rPlantThermoAdaptZone_pft(NZ,NY,NX)
            ELSE
              HighTempLimitSeed_pft(NZ,NY,NX)=30.0+3.0*rPlantThermoAdaptZone_pft(NZ,NY,NX)
            ENDIF
            MatureGroup_pft(NZ,NY,NX)=GROUPX_pft(NZ,NY,NX)+0.30_r8*DTA
            IF(iPlantTurnoverPattern_pft(NZ,NY,NX).NE.0)THEN
              MatureGroup_pft(NZ,NY,NX)=MatureGroup_pft(NZ,NY,NX)/25.0_r8
            ENDIF
            MatureGroup_pft(NZ,NY,NX)=MatureGroup_pft(NZ,NY,NX)-ShootNodeNumAtPlanting_pft(NZ,NY,NX)

          ENDDO
        ENDIF
        !
        !     ADJUST VAPOR PRESSURE FOR TEMPERATURE CHANGE
        !
        IF(isclose(DHUM(N),1.0_r8))THEN
          VPX=VPS(NY,NX)
          !VPS(NY,NX)=0.61*EXP(5360.0*(3.661E-03-1.0/TairK_col(NY,NX))) &
          vps(ny,ny)=vapsat0(TairK_col(ny,nx))*EXP(-AltitudeLandscape_col(NY,NX)/7272.0)
          VPK_col(NY,NX)=VPK_col(NY,NX)*VPS(NY,NX)/VPX
        ENDIF
      ENDIF
!
      !     CHANGES IN VAPOR PRESSURE,RADIATION,WIND SPEED,PRECIPITATION
      !     IRRIGATION ,CO2,NH4,NO3
      !
      !     TDRAD,TDWND,TDHUM=change in radiation,windspeed,vapor pressure
      !     TDPRC,TDIRRI=change in precipitation,irrigation
      !     TDCN4,TDCNO=change in atm CO2,NH4,NO3 concn in precipitation
!
      RadSWDirect_col(NY,NX)  = RadSWDirect_col(NY,NX)*TDRAD(N,NY,NX)
      RadSWDiffus_col(NY,NX)  = RadSWDiffus_col(NY,NX)*TDRAD(N,NY,NX)
      RadDirectPAR_col(NY,NX) = RadDirectPAR_col(NY,NX)*TDRAD(N,NY,NX)
      RadPARDiffus_col(NY,NX) = RadPARDiffus_col(NY,NX)*TDRAD(N,NY,NX)
      WindSpeedAtm_col(NY,NX) = WindSpeedAtm_col(NY,NX)*TDWND(N,NY,NX)
      VPK_col(NY,NX)          = AMIN1(VPS(NY,NX),VPK_col(NY,NX)*TDHUM(N,NY,NX))
      PrecAsRain_col(NY,NX)   = PrecAsRain_col(NY,NX)*TDPRC(N,NY,NX)
      PrecAsSnow_col(NY,NX)   = PrecAsSnow_col(NY,NX)*TDPRC(N,NY,NX)
      PRECII_col(NY,NX)       = PRECII_col(NY,NX)*TDIRI(N,NY,NX)
      PRECUI_col(NY,NX)       = PRECUI_col(NY,NX)*TDIRI(N,NY,NX)
      NH4_rain_mole_conc(NY,NX)    = CN4RI_col(NY,NX)*TDCN4(N,NY,NX)
      NO3_rain_mole_conc(NY,NX)    = CNORI_col(NY,NX)*TDCNO(N,NY,NX)
    ENDDO D9920
  ENDDO D9925
  end subroutine CorrectClimate
!------------------------------------------------------------------------------------------

  subroutine SummaryClimateForc(I,J,NHW,NHE,NVN,NVS,PRECUI_col,PrecAsRain_col,PRECII_col,PrecAsSnow_col)
  !
  !     DESCRIPTION:
  !     DAILY WEATHER TOTALS, MAXIMA AND MINIMA FOR DAILY OUTPUT
  !     CHECK AGAINST INPUTS FROM WEATHER FILE
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), intent(in):: PRECUI_col(JY,JX),PrecAsRain_col(JY,JX)
  real(r8), intent(in) :: PRECII_col(JY,JX),PrecAsSnow_col(JY,JX)
  integer :: NY,NX
  !
  !     begin_execution

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      IF(SineSunInclinationAngle_col(NY,NX).GT.0._r8)then
        TRAD_col(NY,NX)= RadSWDirect_col(NY,NX)*SineSunInclinationAngle_col(NY,NX)+RadSWDiffus_col(NY,NX)*TotSineSkyAngles_grd
      endif
      HUDX_col(NY,NX)  = AMAX1(HUDX_col(NY,NX),VPK_col(NY,NX))          !maximum humidity, vapor pressure, [KPa]
      HUDN_col(NY,NX)  = AMIN1(HUDN_col(NY,NX),VPK_col(NY,NX))          !minimum humidity, vapor pressure, [KPa]
      TWIND_col(NY,NX) = TWIND_col(NY,NX)+WindSpeedAtm_col(NY,NX)       !wind speed, [m/hr]
      VPA_col(NY,NX)   = VPK_col(NY,NX)*2.173E-03_r8/TairK_col(NY,NX)    !atmospheric vapor concentration, [m3 m-3],       2.173E-03_r8 = 18g/mol/(8.3142)
!
      !     WATER AND HEAT INPUTS TO GRID CELLS
      !
      !     AREA=area of grid cell (m2)
      !
      !     PRECR,PRECW=rainfall,snowfall
      !     PRECI,PRECU=surface,subsurface irrigation
      !     PRECA,PrecAtm_col=rain+irrigation,rain+snow
      !     THS=sky LW radiation
      !
      if(.not.ats_cpl_mode)then
        RainFalPrec_col(NY,NX)      = PrecAsRain_col(NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
        SnoFalPrec_col(NY,NX)       = PrecAsSnow_col(NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
      endif

      IrrigSurface_col(NY,NX)     = PRECII_col(NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
      IrrigSubsurf_col(NY,NX)     = PRECUI_col(NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
      Irrigation_col(NY,NX)       = IrrigSurface_col(NY,NX)+IrrigSubsurf_col(NY,NX)
      PrecRainAndIrrig_col(NY,NX) = RainFalPrec_col(NY,NX)+IrrigSurface_col(NY,NX)
      PrecAtm_col(NY,NX)          = RainFalPrec_col(NY,NX)+SnoFalPrec_col(NY,NX)
      LWRadSky_col(NY,NX)         = SkyLonwRad_col(NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
      PrecipAtm2LandSurf_col(NY,NX)=RainFalPrec_col(NY,NX)+SnoFalPrec_col(NY,NX)+IrrigSurface_col(NY,NX)
    ENDDO
  ENDDO
  end subroutine SummaryClimateForc

end module WthrMod
