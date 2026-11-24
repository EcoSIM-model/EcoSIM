
  module DayMod
  use data_kind_mod,      only: r8 => DAT_KIND_R8
  use minimathmod,        only: isLeap,  AZMAX1
  use MiniFuncMod,        only: GetDayLength, calculate_equation_of_time
  use PrescribePhenolMod, only: PrescribePhenologyInterp
  use SurfLitterDataType, only: XTillCorp_col
  use abortutils,         only: iulog
  use DebugToolMod,       only: PrintInfo
  use DepositionMod,     only : ApplyBioAerosol
  use CanopyRadDataType
  use EcosimConst  
  use EcoSIMCtrlMod
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SoilBGCDataType
  use ClimForcDataType
  use FertilizerDataType
  use PlantTraitDataType
  use PlantDataRateType
  use CanopyDataType
  use RootDataType
  use PlantMgmtDataType
  use SOMDataType
  use EcosimBGCFluxType
  use EcoSIMHistMod
  use SoilPropertyDataType
  use IrrigationDataType
  use SedimentDataType
  use GridDataType
  use EcoSIMConfig

  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: day, CompleteDailySetup
  contains

  SUBROUTINE day(I,NHW,NHE,NVN,NVS,DeferStateChanges)
!
!     THIS SUBROUTINE REINITIALIZES DAILY VARIABLES USED IN OTHER
!     SUBROUTINES E.G. LAND MANAGEMENT
!
  implicit none

  integer, intent(in) :: I
  integer, intent(in) :: NHW,NHE,NVN,NVS
  logical, optional, intent(in) :: DeferStateChanges !true until the restart checkpoint is loaded
  logical :: defer_updates
  integer :: NY,NX
  real(r8) :: eot,leapday
  character(len=*), parameter :: subname='day'
!     execution begins here
!     begin_execution
!
!     WRITE DATE
!
!     CDATE=DDMMYYYY
!
  call PrintInfo('beg '//subname)
  defer_updates=.false.
  if(present(DeferStateChanges))defer_updates=DeferStateChanges
  if(isLeap(iYearCurrent))then
    leapday=1._r8
  else
    leapday=0._r8
  endif

  eot=calculate_equation_of_time(I,leapday)
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      SolarNoonHour_col(NY,NX) = SolarNoonHourYM_col(NY,NX)+ eot
    ENDDO
  ENDDO    

  if(defer_updates)then
    ! Reconstruct date-derived forcing while skipping actions on initial state.
    call UpdateDailyAccumulators(I,NHW,NHE,NVN,NVS,StartNewDay=.false.)
  else
    call CompleteDailySetup(I,1,NHW,NHE,NVN,NVS)
  endif

  call PrintInfo('end '//subname)

  END subroutine day

!------------------------------------------------------------------------------------------     
  subroutine CompleteDailySetup(I,J,NHW,NHE,NVN,NVS)
  ! Run after checkpoint loading on the restart day; DAY uses J=1 normally.
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS !J is the first simulated hour, 1:24
  character(len=*), parameter :: subname='CompleteDailySetup'
  integer :: NY,NX
  logical :: StartNewDay

  call PrintInfo('beg '//subname)
  StartNewDay=J==1
  if(StartNewDay)then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        call ApplyBioAerosol(I,12,NY,NX)
      ENDDO
    ENDDO
  endif
  call UpdateDailyAccumulators(I,NHW,NHE,NVN,NVS,StartNewDay)
  call TillageandIrrigationEvents(I,NHW,NHE,NVN,NVS,RunIrrigation=StartNewDay)
  if(ldo_sp_mode)call PrescribePhenologyInterp(I,NHW,NHE,NVN,NVS)
  call PrintInfo('end '//subname)
  end subroutine CompleteDailySetup
!------------------------------------------------------------------------------------------

  subroutine UpdateDailyAccumulators(I, NHW, NHE, NVN, NVS,StartNewDay)
!     WRITE DAILY MAX MIN ACCUMULATORS FOR WEATHER VARIABLES

  implicit none
  integer, intent(in) :: I, NHW, NHE, NVN, NVS
  logical, intent(in) :: StartNewDay !false for intraday reconstruction; preserve restored daily state

!  real(r8) :: AZI
!  REAL(R8) :: DEC
  integer :: NE,NX,NY,I2,I3,ITYPE,N
  character(len=*), parameter :: subname='UpdateDailyAccumulators'

  call PrintInfo('beg '//subname)
  D955: DO NX=NHW,NHE
    D950: DO NY=NVN,NVS

      ! Reset daily accumulators only when beginning a new day. Intraday
      ! reconstruction must preserve the partial sums restored from checkpoint.
      if(StartNewDay)then
        HUDX_col(NY,NX)=0._r8
        HUDN_col(NY,NX)=100.0_r8
        TWIND_col(NY,NX)=0._r8
        SeasonalNonstCDayAve_pft(:,NY,NX) = 0._r8
      endif
!      PrecDaily_col(NY,NX)=0._r8
!
!
!     CALCULATE DAYLENGTH FROM SOLAR ANGLE
!
!     DayLenthPrev,DLYN=daylength of previous,current day
!     ALAT=latitude
!
      if(StartNewDay)DayLenthPrev_col(NY,NX)=DayLensCurr_col(NY,NX)

      DayLensCurr_col(NY,NX)=GetDayLength(ALAT_col(NY,NX),I)
!
!     TIME STEP OF WEARHER DATA
!     ITYPE 1=daily,2=hourly
!
      ITYPE=IWTHR

!
!     PARAMETERS FOR CALCULATING HOURLY RADIATION, TEMPERATURE
!     AND VAPOR PRESSURE FROM DAILY VALUES
!
!     DLYN=daylength
!     SRAD=daily radiation from weather file
!     RMAX=maximum hourly radiation
!     I2,I,I3=previous,current,next day
!     TMPX,TMPN=maximum,minimum daily temperature from weather file
!     DWPT=daily vapor pressure from weather file
!     TAVG*,AMP*,VAVG*,VMP*=daily avgs, amps to calc hourly values in wthr.f
!     KoppenClimZone=Koppen climate zone

      IF(ITYPE.EQ.1)THEN
        IF(KoppenClimZone_col(NY,NX).GE.-1)THEN
          IF(DayLensCurr_col(NY,NX).GT.ZERO)THEN
            RMAX=SRAD(I)/(DayLensCurr_col(NY,NX)*0.658_r8)
          ELSE
            RMAX=0._r8
          ENDIF
        ELSE
          RMAX=SRAD(I)
        ENDIF
        I2=I-1
        I3=I+1
        IF(I.EQ.1)I2=LYRX
        IF(I.EQ.IBEGIN)I2=I
        IF(I.EQ.DazCurrYear)I3=I
        TAVG1=(TMPX(I2)+TMPN(I))/2._r8
        TAVG2=(TMPX(I)+TMPN(I))/2._r8
        TAVG3=(TMPX(I)+TMPN(I3))/2._r8
        AMP1=TAVG1-TMPN(I)
        AMP2=TAVG2-TMPN(I)
        AMP3=TAVG3-TMPN(I3)
        VAVG1=(DWPT(1,I2)+DWPT(2,I))/2._r8
        VAVG2=(DWPT(1,I)+DWPT(2,I))/2._r8
        VAVG3=(DWPT(1,I)+DWPT(2,I3))/2._r8
        VMP1=VAVG1-DWPT(2,I)
        VMP2=VAVG2-DWPT(2,I)
        VMP3=VAVG3-DWPT(2,I3)
      ENDIF
!
!     MODIFIERS TO TEMPERATURE, RADIATION, WIND, HUMIDITY, PRECIPITATION,
!     IRRIGATION AND CO2 INPUTS FROM CLIMATE CHANGES ENTERED IN OPTION
!     FILE IN 'READS'
!
!     ICLM=type of climate change 1=step,2=incremental
!     TDTPX,TDTPN=change in max,min temperature
!     TDRAD,TDWND,TDHUM=change in radiation,windspeed,vapor pressure
!     TDPRC,TDIRRI=change in precipitation,irrigation
!     TDCN4,TDCNO=change in atm CO2,NH4,NO3 concn in precipitation
!
      if(StartNewDay)then
      D600: DO N=1,12
!
!     STEP CHANGES
!
        IF(ICLM.EQ.1)THEN
          TDTPX(N,NY,NX)=DTMPX(N)
          TDTPN(N,NY,NX)=DTMPN(N)
          TDRAD(N,NY,NX)=DRAD(N)
          TDWND(N,NY,NX)=DWIND(N)
          TDHUM(N,NY,NX)=DHUM(N)
          TDPRC(N,NY,NX)=DPREC(N)
          TDIRI(N,NY,NX)=DIRRI(N)
          TDCN4(N,NY,NX)=DCN4R(N)
          TDCNO(N,NY,NX)=DCNOR(N)
!
!     INCRENENTAL CHANGES
!
        ELSEIF(ICLM.EQ.2)THEN
! DazCurrYear: number of days in current year
          TDTPX(N,NY,NX)=TDTPX(N,NY,NX)+DTMPX(N)/DazCurrYear
          TDTPN(N,NY,NX)=TDTPN(N,NY,NX)+DTMPN(N)/DazCurrYear
          TDRAD(N,NY,NX)=TDRAD(N,NY,NX)+(DRAD(N)-1.0_r8)/DazCurrYear
          TDWND(N,NY,NX)=TDWND(N,NY,NX)+(DWIND(N)-1.0_r8)/DazCurrYear
          TDHUM(N,NY,NX)=TDHUM(N,NY,NX)+(DHUM(N)-1.0_r8)/DazCurrYear
          TDPRC(N,NY,NX)=TDPRC(N,NY,NX)+(DPREC(N)-1.0_r8)/DazCurrYear
          TDIRI(N,NY,NX)=TDIRI(N,NY,NX)+(DIRRI(N)-1.0_r8)/DazCurrYear
          TDCN4(N,NY,NX)=TDCN4(N,NY,NX)+(DCN4R(N)-1.0_r8)/DazCurrYear
          TDCNO(N,NY,NX)=TDCNO(N,NY,NX)+(DCNOR(N)-1.0_r8)/DazCurrYear
        ENDIF
      ENDDO D600
      endif
    ENDDO D950
  ENDDO D955
  call PrintInfo('end '//subname)
  END subroutine UpdateDailyAccumulators
!-----------------------------------------------------------------------------------------

  subroutine TillageandIrrigationEvents(I, NHW, NHE, NVN, NVS,RunIrrigation)
  !
  implicit none

  integer, intent(in) :: I, NHW, NHE, NVN, NVS
  logical, optional, intent(in) :: RunIrrigation !false on intraday restart; use the restored schedule
  logical :: decide_irrigation
  integer :: NY,NX,J,L
  real(r8) :: TWP,TVW,TFZ
  real(r8) :: CORP,DIRRA1,DIRRA2,FW,FZ,RR
  character(len=*), parameter :: subname='TillageandIrrigationEvents'

  call PrintInfo('beg '//subname)
  decide_irrigation=.true.
  if(present(RunIrrigation))decide_irrigation=RunIrrigation
  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
!
!     ATTRIBUTE MIXING COEFFICIENTS AND SURFACE ROUGHNESS PARAMETERS
!     TO TILLAGE EVENTS FROM SOIL MANAGEMENT FILE IN 'READS'
!
!     ITILL=soil disturbance type, 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     CORP=soil mixing fraction used in redist.f
!
      CORP=0._r8
      IF(iSoilDisturbType_col(I,NY,NX).LE.10)THEN
        ! type-1 tillage
        CORP=AMIN1(1.0_r8,AZMAX1(iSoilDisturbType_col(I,NY,NX)/10.0_r8))
      ELSEIF(iSoilDisturbType_col(I,NY,NX).LE.20)THEN
        ! type-2 tillage
        CORP=AMIN1(1.0_r8,AZMAX1((iSoilDisturbType_col(I,NY,NX)-10.0_r8)/10.0_r8))
      ENDIF
      !fraction of material not mixed
      XTillCorp_col(NY,NX)=1.0_r8-CORP
!      WRITE(990,*)'TILL',I,iSoilDisturbType_col(I,NY,NX),CORP,XTillCorp_col(NY,NX)
!2227  FORMAT(A8,2I4,12E12.4)
!
!     AUTOMATIC IRRIGATION IF SELECTED
!
!     DATA1(6)=irrigation file name
!     IIRRA=start,finish dates(1,2),hours(3,4) of automated irrigation
!     DIRRX=depth to which water depletion and rewatering is calculated(1)
!     DIRRA=depth to,at which irrigation is applied(1,2)
!     POROS,FC,WP=water content at saturation,field capacity,wilting point
!     CIRRA_col= fraction of FC to which irrigation will raise SWC
!     FW=fraction of soil layer in irrigation zone
!     FZ=SWC at which irrigation is triggered
!     VLSoilPoreMicP_vr,VOLW,VOLI=total,water,ice volume
!     iIrrigOpt_col=flag for irrigation criterion,0=SWC,1=canopy water potential
!     FIRRA_col=depletion of SWC from CIRRA_col to WP(iIrrigOpt_col=0),or minimum canopy
!     water potential(iIrrigOpt_col=1), to trigger irrigation
!     RR=total irrigation requirement [m3 per cell]
!     RRIG=hourly irrigation depth [m h-1]; WthrMod multiplies by cell area
!
      IF(Lirri_auto_col(NY,NX) .and. decide_irrigation)THEN
      !automated irrigation
        IF(I.GE.IIRRA(1,NY,NX).AND.I.LE.IIRRA(2,NY,NX))THEN
          TFZ    = 0._r8
          TWP    = 0._r8
          TVW    = 0._r8
          DIRRA1 = DIRRA(1,NY,NX)+CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX)
          DIRRA2 = DIRRA(2,NY,NX)+CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX)

          D165: DO L=NU_col(NY,NX),NL_col(NY,NX)
            IF(CumDepz2LayBottom_vr(L-1,NY,NX).LT.DIRRA1)THEN
              FW  = AMIN1(1.0_r8,(DIRRA1-CumDepz2LayBottom_vr(L-1,NY,NX))/(CumDepz2LayBottom_vr(L,NY,NX)-CumDepz2LayBottom_vr(L-1,NY,NX)))
              FZ  = AMIN1(POROS_vr(L,NY,NX),WiltPoint_vr(L,NY,NX)+CIRRA_col(NY,NX)*(FieldCapacity_vr(L,NY,NX)-WiltPoint_vr(L,NY,NX)))
              TFZ = TFZ+FW*FZ*VLSoilPoreMicP_vr(L,NY,NX)
              TWP = TWP+FW*WiltPoint_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
              TVW = TVW+FW*(VLWatMicP_vr(L,NY,NX)+VLiceMicP_vr(L,NY,NX))
            ENDIF
          ENDDO D165

          IF((iIrrigOpt_col(NY,NX).EQ.iIrrig_swc .AND. TVW.LT.TWP+FIRRA_col(NY,NX)*(TFZ-TWP)) &
            .OR. (iIrrigOpt_col(NY,NX).EQ.iIrrig_cwp .AND. PSICanPDailyMin_pft(1,NY,NX).LT.FIRRA_col(NY,NX)))THEN
            RR=AZMAX1(TFZ-TVW)
            IF(RR.GT.0.0_r8)THEN
              D170: DO J=IIRRA(3,NY,NX),IIRRA(4,NY,NX)
                RRIG(J,I,NY,NX)=RR/(AREA_3D(3,NU_col(NY,NX),NY,NX) &
                  *(IIRRA(4,NY,NX)-IIRRA(3,NY,NX)+1))
              ENDDO D170
              WDPTH(I,NY,NX)=DIRRA(2,NY,NX)
              WRITE(iulog,2222)'auto',iYearCurrent,I,IIRRA(3,NY,NX),IIRRA(4,NY,NX) &
                ,iIrrigOpt_col(NY,NX),RR,TFZ,TVW,TWP,FIRRA_col(NY,NX),PSICanPDailyMin_pft(1,NY,NX) &
                ,CIRRA_col(NY,NX),DIRRA1,WDPTH(I,NY,NX)

2222  FORMAT(A8,5I6,40E12.4)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO D9990
  ENDDO D9995
  call PrintInfo('end '//subname)
  end subroutine TillageandIrrigationEvents
END module DayMod
