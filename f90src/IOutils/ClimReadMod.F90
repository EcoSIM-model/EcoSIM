module ClimReadMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils , only : endrun
  use ncdio_pio
  use fileUtil   , only : open_safe
  use minimathmod, only : isLeap,AZMAX1,vapsat0
  use ClimForcDataType
  use EcoSIMHistMod
  use FlagDataType
  use EcoSimConst
  USE LandSurfDataType
  use EcoSIMCtrlDataType
  use EcoSIMCtrlMod
  use EcoSIMConfig
  use UnitMod, only : units
  use DebugToolMod, only : PrintInfo
implicit none
  private

  character(len=*), parameter :: mod_filename = &
  __FILE__
  CHARACTER(len=1) :: IVAR(20),VAR(50),TYP(50),CTYPE
  real(r8) :: DATK(50),OUT(50)
  real(r8) :: datav(40)
  integer,save :: IYRD,IFLGY,iYearCorrect
  integer :: IFLG3
  DATA IYRD,IFLGY,iYearCorrect/0,0,0/

  type, public :: atm_forc_type
    real(r8) :: Z0G
    real(r8) :: ZNOONG
    real(r8) :: PHRG
    real(r8) :: CN4RIG
    real(r8) :: CNORIG
    real(r8) :: CPORG
    real(r8) :: CALRG
    real(r8) :: CFERG
    real(r8) :: CCARG
    real(r8) :: CMGRG
    real(r8) :: CNARG
    real(r8) :: CKARG
    real(r8) :: CSORG
    real(r8) :: CCLRG
  end type atm_forc_type
  public :: ReadClim
  public :: ReadClimNC    !read climate data
  public :: ReadAtmosComposition
  public :: get_clm_years
  public :: read_soil_warming_Tref
  contains

!------------------------------------------------------------------------------------------

  subroutine interp3hourweather(I,J)
  !
  !interpolate 3-hourly weather to hourly
  !I is day, J is hour
  implicit none
  integer, intent(in) :: I,J

  integer :: II,JJ
  JJ=J-3
  II=I
  IF(JJ.EQ.0)THEN
    JJ=24
    II=I-1
  ENDIF
!
! INFILL 3-HOURLY WEATHER DATA
!
  IF(II.LT.1)THEN
    TMP_hrly(J-2,I)=TMP_hrly(J,I)
    TMP_hrly(J-1,I)=TMP_hrly(J,I)
    SWRad_hrly(J-2,I)=SWRad_hrly(J,I)
    SWRad_hrly(J-1,I)=SWRad_hrly(J,I)
    WINDH(J-2,I)=WINDH(J,I)
    WINDH(J-1,I)=WINDH(J,I)
    DWPTH(J-2,I)=DWPTH(J,I)
    DWPTH(J-1,I)=DWPTH(J,I)
    RAINH(J,I)=RAINH(J,I)/3.0_r8
    RAINH(J-2,I)=RAINH(J,I)
    RAINH(J-1,I)=RAINH(J,I)
    RadLWClm(J-2,I)=RadLWClm(J,I)
    RadLWClm(J-1,I)=RadLWClm(J,I)
  ELSE
    TMP_hrly(J-2,I)=0.667_r8*TMP_hrly(JJ,II)+0.333_r8*TMP_hrly(J,I)
    TMP_hrly(J-1,I)=0.333_r8*TMP_hrly(JJ,II)+0.667_r8*TMP_hrly(J,I)
    SWRad_hrly(J-2,I)=0.667_r8*SWRad_hrly(JJ,II)+0.333_r8*SWRad_hrly(J,I)
    SWRad_hrly(J-1,I)=0.333_r8*SWRad_hrly(JJ,II)+0.667_r8*SWRad_hrly(J,I)
    WINDH(J-2,I)=0.667_r8*WINDH(JJ,II)+0.333_r8*WINDH(J,I)
    WINDH(J-1,I)=0.333_r8*WINDH(JJ,II)+0.667_r8*WINDH(J,I)
    DWPTH(J-2,I)=0.667_r8*DWPTH(JJ,II)+0.333_r8*DWPTH(J,I)
    DWPTH(J-1,I)=0.333_r8*DWPTH(JJ,II)+0.667_r8*DWPTH(J,I)
    RAINH(J,I)=RAINH(J,I)/3.0_r8
    RAINH(J-2,I)=RAINH(J,I)
    RAINH(J-1,I)=RAINH(J,I)
    RadLWClm(J-2,I)=0.667_r8*RadLWClm(JJ,II)+0.333_r8*RadLWClm(J,I)
    RadLWClm(J-1,I)=0.333_r8*RadLWClm(JJ,II)+0.667_r8*RadLWClm(J,I)
  ENDIF
  end subroutine interp3hourweather

!------------------------------------------------------------------------------------------

  subroutine readhourweather(idat,DAT,IYEAR,I,J,IH,go60,NN,NI,TTYPE)

!     read hourly weather data, from ascii file
!
!     DERIVE DAY I AND HOUR J FROM TIME VARIABLES IVAR
!
  implicit none
  integer, intent(in) :: idat(:)
  real(r8),intent(inout) :: dat(:)
  integer, intent(inout) :: I      !julian day
  integer, intent(inout) :: J      !hour
  integer, intent(inout) :: IH
  integer, intent(in) :: NN,NI,iyear
  character(len=1),intent(in) :: TTYPE
  logical, intent(out) :: go60

  real(r8) :: tempK
  integer :: LPY
  integer :: K
  integer :: M
  integer :: N

  go60=.false.

  D190: DO K=1,NI
    IF(IVAR(K).EQ.'M')THEN
    !  month
      M=IDAT(K)
    ELSEIF(IVAR(K).EQ.'D')THEN
     ! julian day
      N=IDAT(K)
    ELSEIF(IVAR(K).EQ.'H')THEN
     ! hour
      J=IDAT(K)
    ENDIF
    IF(IVAR(K).EQ.'Y')THEN
     ! year information
      IFLGY=1
    ENDIF
  ENDDO D190

  IF(IFLGY.EQ.1.AND.iYearCorrect.LT.iYearCurrent)then
    GO60=.true.
    return
  endif
  IF(CTYPE.EQ.'J')THEN
    I=N
  ELSE
    LPY=0
    IF(isLeap(IYEAR).and.M.GT.2)LPY=1
    IF(M.EQ.1)THEN
      I=N
    ELSE
      I=30*(M-1)+ICOR(M-1)+N+LPY
    ENDIF
  ENDIF

  IF(J.GT.24.AND.(J/100)*100.NE.J)THEN
  !time information given as long hour
    D80: DO K=1,NN
      DATK(K)=DATK(K)+DAT(K)
    ENDDO D80
    IH=IH+1
    go60=.true.
    return
  ENDIF

  IF(J.GT.24)J=INT(J/100)
  IF(J.EQ.0)THEN
    J=24
    I=I-1
    IF(I.LT.1)then
      GO60=.true.
      return
    endif
  ENDIF
!
!     DERIVE START DATE FROM TIME VARIABLES
!
  IF(IFLG3.EQ.0)THEN
    IBEGIN=N
    ISTART=MAX(ISTART,IBEGIN)
    IFLG3=1
  ENDIF

  RadLWClm(J,I)=0.0_r8
!
!     CONVERT HOURLY WEATHER VARIABLES TO MODEL UNITS
!     AND ENTER INTO MODEL ARRAYS
!
!     TMPH=temperature (oC)
!     SWRad_hrly=solar radiation (MJ m-2 h-1)
!     WINDH=windspeed (m h-1)
!     DWPTH=vapor pressure (kPa)
!     RAINH=precipitation (m h-1)
!     RadLWClm=longwave radiation (MJ m-2 h-1)
!
  D95: DO K=1,NN
!
!       TEMPERATURE
!
    IF(VAR(K).EQ.'T')THEN
      IF(TYP(K).EQ.'F')THEN
! Fahrenheit to celcius
        TMP_hrly(J,I)=((DAT(K)+DATK(K))/IH-32.0_r8)*0.556_r8
      ELSEIF(TYP(K).EQ.'K')THEN
! Temperature given as kelvin
        TMP_hrly(J,I)=units%Kelvin2Celcius((DAT(K)+DATK(K))/IH)
      ELSE
        TMP_hrly(J,I)=(DAT(K)+DATK(K))/IH
      ENDIF
!
!         SOLAR RADIATION
!
    ELSEIF(VAR(K).EQ.'R')THEN
      !input is W m-2, to MJ m-2 per hour (3600 seconds per hour * 1.e-6)
      IF(TYP(K).EQ.'W')THEN
        SWRad_hrly(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.0036_r8)
      !input is 1.e4 J m-2/hour        
      ELSEIF(TYP(K).EQ.'J')THEN
        SWRad_hrly(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.01_r8)
      ! input is kJ m-2 per hour  
      ELSEIF(TYP(K).EQ.'K')THEN
        SWRad_hrly(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.001_r8)
      !input is Roentgen, which is 0.457 milligray, be careful for this  
      ELSEIF(TYP(K).EQ.'P')THEN
        SWRad_hrly(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.0036_r8*0.457_r8)
      !input is   
!     ELSEIF(TYP(K).EQ.'M')THEN
!     SWRad_hrly(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*3.6*0.457)
      !input is MJ/hour
      ELSE
        SWRad_hrly(J,I)=AZMAX1((DAT(K)+DATK(K))/IH)
      ENDIF
!
!     WIND SPEED
!
    ELSEIF(VAR(K).EQ.'W')THEN
      ! input is m s-1, into m per hour
      IF(TYP(K).EQ.'S')THEN
        WINDH(J,I)=(DAT(K)+DATK(K))/IH*3600.0_r8
      ! input is km per hour into m per hour
      ELSEIF(TYP(K).EQ.'H')THEN
        WINDH(J,I)=(DAT(K)+DATK(K))/IH*1000.0_r8
      ! input is mile per hour into m per hour, it should be 1609.34 though        
      ELSEIF(TYP(K).EQ.'M')THEN
        WINDH(J,I)=(DAT(K)+DATK(K))/IH*1609.34_r8
      !input is m/hour  
      ELSE
        WINDH(J,I)=(DAT(K)+DATK(K))/IH
      ENDIF
    ELSEIF(VAR(K).EQ.'H')THEN
!
!       VAPOR PRESSURE
!
      ! input is celcius degree of dew point
      IF(TYP(K).EQ.'D')THEN
        tempK=units%Celcius2Kelvin((DAT(K)+DATK(K))/IH)
        DWPTH(J,I)=vapsat0(tempK)
      ! input is Fahrenheit dew point
      ELSEIF(TYP(K).EQ.'F')THEN       
        DAT(K)=(DAT(K)-32.0_r8)*0.556_r8
        tempK=units%Celcius2Kelvin((DAT(K)+DATK(K))/IH)
        DWPTH(J,I)=vapsat0(tempK)
      ! input is relative humidity, [0, 1], return value kPa
      ELSEIF(TYP(K).EQ.'H')THEN
        tempK=units%Celcius2Kelvin(TMP_hrly(J,I))
        DWPTH(J,I)=vapsat0(tempK)*AZMAX1(AMIN1(1.0_r8,(DAT(K)+DATK(K))/IH))
      ! input is relative humidity [0, 100],return value kPa        
      ELSEIF(TYP(K).EQ.'R')THEN     
        tempK=units%Celcius2Kelvin(TMP_hrly(J,I))
        DWPTH(J,I)=vapsat0(tempK)*AZMAX1(AMIN1(1._r8,0.01_r8*(DAT(K)+DATK(K))/IH))
      ! input is molar mixing ratio? [0,100]        
      ELSEIF(TYP(K).EQ.'S')THEN
        DWPTH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH)*0.0289_r8/18.0_r8*101.325_r8 &
          *EXP(-ALTIG/7272.0_r8)*288.15_r8/units%Celcius2Kelvin(TMP_hrly(J,I))
      !input is molar mixing ratio [0, 1] ALTIG is doing altitude correction,          
      ELSEIF(TYP(K).EQ.'G')THEN
        DWPTH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH)*28.9_r8/18.0_r8*101.325_r8 &
          *EXP(-ALTIG/7272.0_r8)*288.15_r8/units%Celcius2Kelvin(TMP_hrly(J,I))
      !input is hPa          
      ELSEIF(TYP(K).EQ.'M')THEN
        DWPTH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.1_r8)
      !input is kPa  
      ELSE
        DWPTH(J,I)=AZMAX1((DAT(K)+DATK(K))/IH)
      ENDIF
!
!       PRECIPITATION, [m h-1]
!
  ELSEIF(VAR(K).EQ.'P')THEN
    !input is mm/m2/hour
    IF(TYP(K).EQ.'M')THEN
    !   [mm m-2] (mm to m, 1.e-3)
      RAINH(J,I)=AZMAX1(DAT(K)+DATK(K))/1.0E+03_r8
    !input is cm/m2/hour   
    ELSEIF(TYP(K).EQ.'C')THEN
      !   [cm m-2]
      RAINH(J,I)=AZMAX1(DAT(K)+DATK(K))/1.0E+02_r8
    !input is inch/m2/hour  
    ELSEIF(TYP(K).EQ.'I')THEN
    !  [inch m-2] -> [m m-2]
      RAINH(J,I)=AZMAX1(DAT(K)+DATK(K))*0.0254_r8
    !input is m/m2/second  
    ELSEIF(TYP(K).EQ.'S')THEN
     !  [second m-2]
      IF(TTYPE.EQ.'H')THEN
        ! given in hour
        RAINH(J,I)=AZMAX1(DAT(K)+DATK(K))*3.6_r8
      ELSE
        ! given in half hour, 
        RAINH(J,I)=AZMAX1(DAT(K)+DATK(K))*1.8_r8
      ENDIF
    !input is m/  
    ELSE
      RAINH(J,I)=AZMAX1(DAT(K)+DATK(K))
    ENDIF
!
!       LONGWAVE RADIATION (OPTIONAL)
!
  ELSEIF(VAR(K).EQ.'L')THEN
    IF(TYP(K).EQ.'W')THEN
      ! watss m-2, into MJ per hour, *3600 seconds * 1.0-6 = 0.0036
      RadLWClm(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.0036_r8)
    ELSEIF(TYP(K).EQ.'J')THEN
      RadLWClm(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.01_r8)
    ELSEIF(TYP(K).EQ.'K')THEN
        ! kJ m-2, into MJ per hour, * 0.001
        RadLWClm(J,I)=AZMAX1((DAT(K)+DATK(K))/IH*0.001_r8)
      ELSE
        RadLWClm(J,I)=AZMAX1((DAT(K)+DATK(K))/IH)
      ENDIF
    ENDIF
    DATK(K)=0.0_r8
  ENDDO D95
  end subroutine readhourweather
!------------------------------------------------------------------------------------------

  subroutine readdayweather(idat,dat,IYEAR,I,go110,go60,NTX,NFX,NN,NI,IX)
!     read daily weather data
!
!     DERIVE DAY I FROM TIME VARIABLES IVAR
!
!     IWTHR=weather data type:1=daily,2=hourly for first(L=1) or second(L=2) scene
!
  implicit none
  integer, intent(in) :: idat(:)
  real(r8),intent(inout) :: dat(:)
  integer, intent(in) :: NTX,NFX,NN,NI,IYEAR
  integer, intent(out) :: I
  logical, intent(out) :: go110, go60
  INTEGER, INTENT(OUT) :: IX
  integer :: LPY
  integer :: K,M,N

  go110=.false.
  go60=.false.

  D160: DO K=1,NI
    IF(IVAR(K).EQ.'M')THEN
      M=IDAT(K)
    ELSEIF(IVAR(K).EQ.'D')THEN
      N=IDAT(K)
    ENDIF
    IF(IVAR(K).EQ.'Y')THEN
      IFLGY=1
! compute the shifted/correced year
      iYearCorrect=IDAT(K)+(NTX-1)*NFX
      IF(isLeap(IDAT(K)))then
        IYRD=366
      else
        IYRD=365
      endif
    ENDIF
  ENDDO D160

  IF(IFLGY.EQ.1.AND.iYearCorrect.LT.iYearCurrent)THEN
! if the year read in is not equal to year required
    GO60=.TRUE.
    RETURN
  ENDIF

! CTYPE for calendar format
  IF(CTYPE.EQ.'J')THEN
    I=N
  ELSE
    LPY=0
    IF(isLeap(IYEAR).and.M.GT.2)LPY=1
    IF(M.EQ.1)THEN
      I=N
    ELSE
      I=30*(M-1)+ICOR(M-1)+N+LPY
    ENDIF
  ENDIF
!
! DERIVE START DATE FROM TIME VARIABLES
!
  IF(IFLG3.EQ.0)THEN
    IBEGIN=I
    ISTART=MAX(ISTART,IBEGIN)
    IFLG3=1
  ENDIF
!
! CONVERT DAILY WEATHER VARIABLES TO MODEL UNITS
! AND ENTER INTO MODEL ARRAYS
!
! TMPX,TMPN=maximum,minimum temperature (OC)
! SRAD=solar radiation (MJ m-2 d-1)
! WIND=windspeed (m h-1)
! DWPT=vapor pressure (kPa)
! RAIN=precipitation (mm d-1)
!
  D65: DO K=1,NN
!
!   MAX,MIN REMPERATURE
!   maximum temperature
    IF(VAR(K).EQ.'M')THEN
      !input is Fahrenheit
      IF(TYP(K).EQ.'F')THEN
        TMPX(I)=(DAT(K)-32.0_r8)*0.556_r8
      !input is kelvin  
      ELSEIF(TYP(K).EQ.'K')THEN
        TMPX(I)=DAT(K)-273.15_r8
      !input is celcius  
      ELSE
        TMPX(I)=DAT(K)
      ENDIF
    !minimum temperature  
    ELSEIF(VAR(K).EQ.'N')THEN
      !input is Fahrenheit
      IF(TYP(K).EQ.'F')THEN
        TMPN(I)=(DAT(K)-32.0_r8)*0.556_r8
      !input is kelvin
      ELSEIF(TYP(K).EQ.'K')THEN
        TMPN(I)=DAT(K)-273.16_r8
      !input is celcius  
      ELSE
        TMPN(I)=DAT(K)
      ENDIF
!
!     SOLAR RADIATION
!
    ELSEIF(VAR(K).EQ.'R')THEN
      !input is 1.e4 Calories/day
      IF(TYP(K).EQ.'L')THEN
        SRAD(I)=AZMAX1(DAT(K)/23.87_r8)
      !input is 1.e4 Joule/day  
      ELSEIF(TYP(K).EQ.'J')THEN
        SRAD(I)=AZMAX1(DAT(K)*0.01_r8)
      !input is MJ/day 
      ELSE
        SRAD(I)=AZMAX1(DAT(K))
      ENDIF
!
!     WIND SPEED
!
    ELSEIF(VAR(K).EQ.'W')THEN
      !input is m/s
      IF(TYP(K).EQ.'S')THEN
        WIND(I)=ABS(DAT(K))*3600.0_r8
      !input is km/hour  
      ELSEIF(TYP(K).EQ.'H')THEN
        WIND(I)=ABS(DAT(K))*1000.0_r8
      !input is km/day  
      ELSEIF(TYP(K).EQ.'D')THEN
        WIND(I)=ABS(DAT(K))*1000.0_r8/24.0_r8
      !input is mile/hour  
      ELSEIF(TYP(K).EQ.'M')THEN
        WIND(I)=ABS(DAT(K))*1609.34_r8
      !input is m/day  
      ELSE
        WIND(I)=ABS(DAT(K))
      ENDIF
!
!     VAPOR PRESSURE (dew point)
!
    ELSEIF(VAR(K).EQ.'H')THEN
      !input is degree
      IF(TYP(K).EQ.'D')THEN
        DWPT(1,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/units%Celcius2Kelvin(DAT(K))))
        DWPT(2,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/units%Celcius2Kelvin(DAT(K))))
      !input is Fahrenheit  
      ELSEIF(TYP(K).EQ.'F')THEN
        DAT(K)=(DAT(K)-32.0_r8)*0.556_r8
        DWPT(1,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/units%Celcius2Kelvin(DAT(K))))
        DWPT(2,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/units%Celcius2Kelvin(DAT(K))))
      !input is relative humidity reported in absolute value
      ELSEIF(TYP(K).EQ.'H')THEN
        DAT(K)=AZMAX1(AMIN1(1.0_r8,DAT(K)))
        DWPT(1,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/units%Celcius2Kelvin((TMPN(I)+TMPX(I))/2._r8)))*DAT(K)
        DWPT(2,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/units%Celcius2Kelvin(TMPN(I))))
      !input is relative humidity  
      ELSEIF(TYP(K).EQ.'R')THEN
        DAT(K)=AZMAX1(AMIN1(100.0_r8,DAT(K)))
        DWPT(1,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/units%Celcius2Kelvin((TMPN(I)+TMPX(I))/2._r8)))*DAT(K)*0.01_r8
        DWPT(2,I)=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/units%Celcius2Kelvin(TMPN(I))))
     !input is specific humidity kg/kg   
      ELSEIF(TYP(K).EQ.'S')THEN
        DWPT(1,I)=AZMAX1(DAT(K))*0.0289_r8/18.0_r8*101.325_r8 &
          *EXP(-ALTIG/7272.0_r8)*288.15_r8/units%Celcius2Kelvin((TMPN(I)+TMPX(I))/2._r8)
        DWPT(2,I)=AZMAX1(DAT(K))*0.0289_r8/18.0_r8*101.325_r8 &
          *EXP(-ALTIG/7272.0_r8)*288.15_r8/units%Celcius2Kelvin(TMPN(I))
      !input is specific humidity g/g    
      ELSEIF(TYP(K).EQ.'G')THEN
        DWPT(1,I)=AZMAX1(DAT(K))*28.9_r8/18.0_r8*101.325_r8 &
          *EXP(-ALTIG/7272.0_r8)*288.15_r8/units%Celcius2Kelvin((TMPN(I)+TMPX(I))/2._r8)
        DWPT(2,I)=AZMAX1(DAT(K))*28.9_r8/18.0_r8*101.325_r8 &
          *EXP(-ALTIG/7272.0_r8)*288.15_r8/units%Celcius2Kelvin(TMPN(I))
      !input is hPa, -> kPa    
      ELSEIF(TYP(K).EQ.'M')THEN
        DWPT(1,I)=AZMAX1(DAT(K)*0.1_r8)
        DWPT(2,I)=AZMAX1(DAT(K)*0.1_r8)
      !input is kPa  
      ELSE
        DWPT(1,I)=AZMAX1(DAT(K))
        DWPT(2,I)=AZMAX1(DAT(K))
      ENDIF
!
!     PRECIPITATION
!
    ELSEIF(VAR(K).EQ.'P')THEN
      !input is mm /day
      IF(TYP(K).EQ.'M')THEN
        RAIN(I)=AZMAX1(DAT(K))/1.0E+03_r8
      !input is cm/day  
      ELSEIF(TYP(K).EQ.'C')THEN
        RAIN(I)=AZMAX1(DAT(K))/1.0E+02_r8
      !input is inch/day  
      ELSEIF(TYP(K).EQ.'I')THEN
        RAIN(I)=AZMAX1(DAT(K))*0.0254_r8
      !input is m/day  
      ELSE
        RAIN(I)=AZMAX1(DAT(K))
      ENDIF
    ENDIF
  ENDDO D65
  IX=I
  IF(IFLGY.EQ.1.AND.I.EQ.IYRD)THEN
    GO110=.true.
    return
  ENDIF
  GO60=.TRUE.
  end subroutine readdayweather

!------------------------------------------------------------------------------------------

  function getclimttype(ttype)result(ans)
!  return climate data time step
  implicit none
  character(len=1), intent(in):: ttype
  character(len=16) :: ans

  select case (ttype)
  case ('D')
    ans='daily'
  case ('H')
    ans='hourly'
  case ('3')
    ans='3-hourly'
  end select
  end function getclimttype

!------------------------------------------------------------------------------------------

  function getclimctype(ctype)result(ans)
! return climate calendar type
  implicit none
  character(len=1), intent(in) :: ctype
  character(len=16) :: ans

  select case(ctype)
  case ('J')
    ans='Julian day'
  case default
    ans='not julian day'
  end select
  end function getclimctype

!------------------------------------------------------------------------------------------

  subroutine ReadClim(iyear,clmfile,NTX,NFX,I,IX,TTYPE,atmf)
  !
  !read climate/weather data, 
  implicit none
  character(len=*), intent(in) :: clmfile
  integer, intent(in) :: iyear
  integer, intent(in) :: NTX,NFX
  integer, intent(out) :: I,IX
  CHARACTER(len=1),intent(out) :: TTYPE
  type(atm_forc_type), intent(out) :: atmf
  real(r8) :: DAT(50)
  integer :: IDAT(20)

  integer :: K, KK, IH, NI, NN, J
  integer :: LL,iostat
  LOGICAL :: GO110,GO60
  character(len=100):: line

! OPEN WEATHER file(3,
  call OPEN_safe(3,PREFIX,clmfile,'OLD',mod_filename,__LINE__)
  IFLG3=0
  READ(3,'(2A1,2I2,50A1)')TTYPE,CTYPE,NI,NN,(IVAR(K),K=1,NI),(VAR(K),K=1,NN)
  READ(3,'(50A1)')(TYP(K),K=1,NN)
  read(3,*)(datav(kk),kk=1,3)
  atmf%Z0G=datav(1)
  iFlagRaiseZ0GbyVeg=int(datav(2))
  atmf%ZNOONG=datav(3)

!  fourth line in the weather file
  read(3, '(A)') line
  READ(line,*,iostat=iostat)atmf%PHRG,atmf%CN4RIG,atmf%CNORIG,atmf%CPORG,atmf%CALRG,&
    atmf%CFERG,atmf%CCARG,atmf%CMGRG,atmf%CNARG,atmf%CKARG,atmf%CSORG,atmf%CCLRG
  if(iostat/=0)then
    print*,line
    call endrun("Error reading rainfall chem composition data in "//mod_filename, __LINE__)
  endif
  if(lverb)then
    write(*,'(40A)')('-',ll=1,40)
    write(*,*)'read weather file head from ',clmfile
    write(*,*)'time step format: TTYPE ',TTYPE,trim(getclimttype(TTYPE))
    write(*,*)'calendar format: CTYPE ',CTYPE,trim(getclimctype(CTYPE))
    write(*,*)'number of time variables: NI ',NI
    write(*,*)'time var type: IVAR ',(IVAR(K),K=1,NI)
    write(*,*)'number of weather data variables: NN ',NN
    write(*,*)'weather var type: VAR ',(VAR(K),K=1,NN)
    write(*,*)'weather var units: TYP ',(TYP(K),K=1,NN)
    write(*,*)'windspeed measurement height: Z0G',atmf%Z0G
    write(*,*)'flag for raising Z0G with vegn: IFLGW ',iFlagRaiseZ0GbyVeg
    write(*,*)'time of solar noon: ZNOONG ',atmf%ZNOONG
    write(*,*)'pH in precipitation: PHRG ',atmf%PHRG
    write(*,*)'NH4 conc in precip: CN4RIG ',atmf%CN4RIG
    write(*,*)'NO3 conc in precip: CNORIG ',atmf%CNORIG
    write(*,*)'H2PO4 conc in precip: CPORG ',atmf%CPORG
    write(*,*)'Al conc in precip: CALRG ',atmf%CALRG
    write(*,*)'Fe conc in precip: CFERG ',atmf%CFERG
    write(*,*)'Ca conc in precip: CCARG ',atmf%CCARG
    write(*,*)'Mg conc in precip: CMGRG ',atmf%CMGRG
    write(*,*)'Na conc in precip: CNARG ',atmf%CNARG
    write(*,*)'K conc in precip: CKARG ',atmf%CKARG
    write(*,*)'SO4 conc in precip: CSORG ',atmf%CSORG
    write(*,*)'Cl conc in precip: CCLRG ',atmf%CCLRG
    write(*,*)'weather data are in the format time,weather variable'
  endif
  D55: DO K=1,NN
    DATK(K)=0.0_r8
  ENDDO D55
  IH=1

! the file reading loop
  DO while(.TRUE.)
    !read one line data
    read(3,*,END=111)(datav(k),k=1,NI),(DAT(K),K=1,NN)

!   time information
    do k = 1, ni
      idat(k)=int(datav(k))
    enddo
!
    IF(TTYPE.EQ.'D')THEN
!   READ DAILY WEATHER DATA AND CONVERT TO MODEL UNITS
      if(lverb)write(*,*)'read daily weather file'
      IWTHR=1
      call readdayweather(idat,dat,iyear,I,GO110,GO60,NTX,NFX,NN,NI,IX)
      IF(GO60)cycle  !year mismatch, read the next line
      IF(GO110)EXIT
!!
    ELSE
!     READ HOURLY WEATHER DATA AND CONVERT TO MODEL UNITS
      IWTHR=2
      if(lverb)write(*,*)'read hourly weather file'
      call readhourweather(idat,dat,iyear,I,J,IH,go60,NN,NI,TTYPE)
      IF(GO60)cycle  !year mismatch,read the next year
      IH=1
      IX=I
      IF(TTYPE.EQ.'3')THEN
!       weather data is every 3 hrs, do interpolation
        call interp3hourweather(I,J)
      ENDIF

      IF(IFLGY.EQ.1.AND.I.EQ.IYRD.AND.J.EQ.24)THEN
        EXIT
      ENDIF
      cycle
    ENDIF
  enddo

111  CLOSE(3)
  IFLGY=0
  end subroutine ReadClim


!----------------------------------------------------------------------

  subroutine ReadClimNC(yearc,yeari,L,atmf)
  !
  !DESCRIPTION
  !read climate data from netcdf files
  use EcoSIMCtrlMod, only : clm_hour_file_in
  use EcoSIMCtrlMod, only : clm_day_file_in
  use abortutils, only : endrun
  implicit none
  integer, intent(in) :: yearc  ! current year of simulation
  integer, intent(in) :: yeari  ! current year of forcing data
  integer, intent(in) :: L
  type(atm_forc_type), intent(out) :: atmf

  type(file_desc_t) :: clm_nfid
  integer :: nyears,ngrid,LYR
  type(Var_desc_t) :: vardesc
  logical :: readvar
  real(r8), allocatable :: fdatam(:,:,:)
  real(r8), allocatable :: fdatam1(:,:)
  integer , allocatable :: idatav(:)
  real(r8), allocatable ::fdatav(:)
  integer :: year,J,II,JJ,I,irec
  integer, parameter :: iHourly=2
  integer, parameter :: iDaily=1
  
  IWTHR=get_forc_step_type(yeari)
  
  LYR=0

  if(IWTHR==iDaily)then
   !open file
    call ncd_pio_openfile(clm_nfid, clm_day_file_in, ncd_nowrite)

    nyears=get_dim_len(clm_nfid, 'year')

    do irec=1,nyears
      call ncd_getvar(clm_nfid,'year',irec,year)
      if(year==yeari)exit
    enddo
    ngrid=get_dim_len(clm_nfid,'ngrid')
    if(ngrid/=1)then
      call endrun("ngrid > 1 in "//mod_filename,__LINE__)
    endif
    allocate(idatav(ngrid))
    allocate(fdatav(ngrid))
    allocate(fdatam1(ngrid,366))
    I=365
    if(isleap(yeari))I=366
    !there is only grid
    call ncd_getvar(clm_nfid,'TMPX',irec,fdatam1); call reshape1(TMPX,fdatam1)
    call ncd_getvar(clm_nfid,'TMPN',irec,fdatam1); call reshape1(TMPN,fdatam1)
    call ncd_getvar(clm_nfid,'WIND',irec,fdatam1); call reshape1(WIND,fdatam1)    
    call ncd_getvar(clm_nfid,'RAIN',irec,fdatam1); call reshape1(RAIN,fdatam1)
    call ncd_getvar(clm_nfid,'SRAD',irec,fdatam1); call reshape1(SRAD,fdatam1)
    call ncd_getvar(clm_nfid,'DWPT',irec,fdatam1); DWPT(1,:)=fdatam1(1,:);DWPT(2,:)=DWPT(1,:)

    call ncd_getvar(clm_nfid,'Z0G',irec,fdatav); atmf%Z0G=fdatav(1)
    call ncd_getvar(clm_nfid,'ZNOONG',irec,fdatav); atmf%ZNOONG=fdatav(1)
    call ncd_getvar(clm_nfid,'PHRG',irec,fdatav); atmf%PHRG=fdatav(1)
    call ncd_getvar(clm_nfid,'CN4RIG',irec,fdatav); atmf%CN4RIG=fdatav(1)
    call ncd_getvar(clm_nfid,'CNORIG',irec,fdatav); atmf%CNORIG=fdatav(1)
    call ncd_getvar(clm_nfid,'CPORG',irec,fdatav); atmf%CPORG=fdatav(1)
    call ncd_getvar(clm_nfid,'CALRG',irec,fdatav); atmf%CALRG=fdatav(1)
    call ncd_getvar(clm_nfid,'CFERG',irec,fdatav); atmf%CFERG=fdatav(1)
    call ncd_getvar(clm_nfid,'CCARG',irec,fdatav); atmf%CCARG=fdatav(1)
    call ncd_getvar(clm_nfid,'CMGRG',irec,fdatav); atmf%CMGRG=fdatav(1)
    call ncd_getvar(clm_nfid,'CNARG',irec,fdatav); atmf%CNARG=fdatav(1)
    call ncd_getvar(clm_nfid,'CKARG',irec,fdatav); atmf%CKARG=fdatav(1)
    call ncd_getvar(clm_nfid,'CSORG',irec,fdatav); atmf%CSORG=fdatav(1)
    call ncd_getvar(clm_nfid,'CCLRG',irec,fdatav); atmf%CCLRG=fdatav(1)
    call ncd_getvar(clm_nfid,'IFLGW',irec,idatav); iFlagRaiseZ0GbyVeg=idatav(1)
    
    call check_var(clm_nfid, 'DRYDEPOC', vardesc, readvar)
    if(readvar)then
      !dry deposition is provided, now (year,day,ngrid)
      call ncd_getvar(clm_nfid,'DRYDEPOC',irec,fdatam1)
      DryDepoOMC(:) = fdatam1(:,1)  ! extract first (only) grid cell
    else
      DryDepoOMC=0._r8
    endif
    !fill in day 366
    if (.not. isLeap(yeari))then
      I=365
      TMPX(I+1)=TMPX(I)
      TMPN(I+1)=TMPN(I)
      WIND(I+1)=WIND(I)
      RAIN(I+1)=RAIN(I)
      SRAD(I+1)=SRAD(I)
      DWPT(:,I+1)=DWPT(:,I)
    endif

    deallocate(fdatam1)
    deallocate(fdatav)
    deallocate(idatav)

  elseif(IWTHR==iHourly)then
    !hourly data
    call ncd_pio_openfile(clm_nfid, clm_hour_file_in, ncd_nowrite)

    nyears=get_dim_len(clm_nfid, 'year')
    do irec=1,nyears
      call ncd_getvar(clm_nfid,'year',irec,year)
      if(year==yeari)exit
    enddo

    ! while the climate data may be for multiple grid, only the first grid of the data
    ! read is assigned to current climate arrays
    ngrid=get_dim_len(clm_nfid,'ngrid')
    if(ngrid/=1)then
      call endrun("ngrid > 1 in "//mod_filename,__LINE__)
    endif
    allocate(fdatam(ngrid,24,366))
    allocate(idatav(ngrid))
    allocate(fdatav(ngrid))
    allocate(fdatam1(ngrid,366))    
    I=365
    if(isleap(yeari))I=366

    call ncd_getvar(clm_nfid,'TMPH',irec,fdatam); call reshape2(TMP_hrly,fdatam)
    call ncd_getvar(clm_nfid,'WINDH',irec,fdatam); call reshape2(WINDH,fdatam)
    call ncd_getvar(clm_nfid,'DWPTH',irec,fdatam); call reshape2(DWPTH,fdatam)
    call ncd_getvar(clm_nfid,'RAINH',irec,fdatam); call reshape2(RAINH,fdatam)
    call ncd_getvar(clm_nfid,'SRADH',irec,fdatam); call reshape2(SWRad_hrly,fdatam)

    call check_var(clm_nfid, 'PATM', vardesc, readvar,print_err=.false.)
    if(readvar)then
      call ncd_getvar(clm_nfid,'PATM',irec,fdatam); call reshape2(PBOT_hrly,fdatam)
    else
      PBOT_hrly=1.01325E+02_r8
    endif

    call check_var(clm_nfid, 'DRYDEPOC', vardesc, readvar)
    if(readvar)then
      !dry deposition is provided, now (year,day,ngrid)
      call ncd_getvar(clm_nfid,'DRYDEPOC',irec,fdatam1)
      DryDepoOMC(:) = fdatam1(1,:)  ! extract first (only) grid cell
    else
      DryDepoOMC=0._r8
    endif

    call ncd_getvar(clm_nfid,'Z0G',irec,fdatav); atmf%Z0G=fdatav(1)
    call ncd_getvar(clm_nfid,'ZNOONG',irec,fdatav); atmf%ZNOONG=fdatav(1)
    call ncd_getvar(clm_nfid,'PHRG',irec,fdatav); atmf%PHRG=fdatav(1)
    call ncd_getvar(clm_nfid,'CN4RIG',irec,fdatav); atmf%CN4RIG=fdatav(1)
    call ncd_getvar(clm_nfid,'CNORIG',irec,fdatav); atmf%CNORIG=fdatav(1)
    call ncd_getvar(clm_nfid,'CPORG',irec,fdatav); atmf%CPORG=fdatav(1)
    call ncd_getvar(clm_nfid,'CALRG',irec,fdatav); atmf%CALRG=fdatav(1)
    call ncd_getvar(clm_nfid,'CFERG',irec,fdatav); atmf%CFERG=fdatav(1)
    call ncd_getvar(clm_nfid,'CCARG',irec,fdatav); atmf%CCARG=fdatav(1)
    call ncd_getvar(clm_nfid,'CMGRG',irec,fdatav); atmf%CMGRG=fdatav(1)
    call ncd_getvar(clm_nfid,'CNARG',irec,fdatav); atmf%CNARG=fdatav(1)
    call ncd_getvar(clm_nfid,'CKARG',irec,fdatav); atmf%CKARG=fdatav(1)
    call ncd_getvar(clm_nfid,'CSORG',irec,fdatav); atmf%CSORG=fdatav(1)
    call ncd_getvar(clm_nfid,'CCLRG',irec,fdatav); atmf%CCLRG=fdatav(1)
    call ncd_getvar(clm_nfid,'IFLGW',irec,idatav); iFlagRaiseZ0GbyVeg=idatav(1)

    call check_var(clm_nfid, 'XRADH', vardesc, readvar,print_err=.false.)
    if(readvar)then
      call ncd_getvar(clm_nfid,'XRADH',irec,fdatam); call reshape2(RadLWClm,fdatam)
    else
      RadLWClm=0._r8
    endif

    if (.not. isLeap(yeari))then
      I=365
      DO J=1,24
        TMP_hrly(J,I+1) = TMP_hrly(J,I)
        WINDH(J,I+1)= WINDH(J,I)
        DWPTH(J,I+1)= DWPTH(J,I)
        RAINH(J,I+1)= RAINH(J,I)
        SWRad_hrly(J,I+1)= SWRad_hrly(J,I)
        RadLWClm(J,I+1)= RadLWClm(J,I)
        PBOT_hrly(J,I+1)=PBOT_hrly(J,I)
      ENDDO
    endif

    if(isleap(yearc))LYR=1
    DO II=1,366
      DO JJ=1,24
        SWRad_hrly(JJ,II)=SWRad_hrly(JJ,II)*3600.e-6_r8  !convert from W/m2 to MJ/m^2/hour
        RadLWClm(JJ,II)=RadLWClm(JJ,II)*3600.e-6_r8       !convert from W/m2 to MJ/m^2/hour
        WINDH(JJ,II)=WINDH(JJ,II)*3600._r8        !convert from m/s to m/hour
        RAINH(JJ,II)=RAINH(JJ,II)*1.e-3_r8        !convert into m hr^-1
      enddo
    ENDDO

    if(lverb)then
    DO II=1,365+LYR
      DO JJ=1,24
        print*,TMP_hrly(JJ,II),SWRad_hrly(JJ,II),WINDH(JJ,II),RAINH(JJ,II),DWPTH(JJ,II)
      ENDDO
    ENDDO
    endif
    deallocate(fdatam)
    deallocate(fdatav)
    deallocate(idatav)
    deallocate(fdatam1)    
  endif
  call ncd_pio_closefile(clm_nfid)
  end subroutine ReadClimNC
!----------------------------------------------------------------------

  subroutine reshape2(arr2d,arr3d)
  implicit none
  real(r8), dimension(:,:), intent(out) :: arr2d
  real(r8), dimension(:,:,:),intent(in) :: arr3d

  integer :: ii,jj
  DO ii=1,366
    DO jj=1,24
      arr2d(jj,ii)=arr3d(1,jj,ii)
    ENDDO
  ENDDO
  end subroutine reshape2
!----------------------------------------------------------------------

  subroutine reshape1(arr1d,arr2d)
  implicit none
  real(r8), dimension(:), intent(OUT) :: arr1d
  real(r8), dimension(:,:),intent(IN) :: arr2d

  integer :: ii
  DO ii=1,366
    arr1d(ii)=arr2d(1,ii)
  ENDDO
  end subroutine reshape1


!------------------------------------------------------------------------------------------

  subroutine ReadAtmosComposition(yeari,NHW,NHE,NVN,NVS)
  !
  !DESCRIPTION
  !read in atmospheric concentrations for CO2, CH4, and N2O
  !for year yeari
  implicit none
  integer, intent(in) :: yeari
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8) :: atm_co2   !ppm
  real(r8) :: atm_ch4   !ppb
  real(r8) :: atm_n2o   !ppb
  type(file_desc_t) :: atm_ghg_nfid
  integer :: iyear,ntime
  real(r8) :: year0,yearN
  character(len=512) :: msg
  INTEGER :: NY,NX,K
  character(len=*), parameter :: subname='ReadAtmosComposition'

  call PrintInfo('beg '//subname)

  call ncd_pio_openfile(atm_ghg_nfid, atm_ghg_in, ncd_nowrite)

  iyear=0
  if(atm_co2_fix<=0._r8 .or. atm_ch4_fix<=0._r8 .or. atm_n2o_fix<=0._r8)then
    !time stamps are calendar_year+month/12: Jan of Y is Y+1/12, Dec of Y is Y+1
    ntime=get_dim_len(atm_ghg_nfid,'time')
    call ncd_getvar(atm_ghg_nfid,'year',1,year0)
    call ncd_getvar(atm_ghg_nfid,'year',ntime,yearN)
    if(nint((year0-int(year0))*12._r8)/=1)then
      call endrun(subname//': first record of '//trim(atm_ghg_in)// &
        ' is not January (expect year stamp = calendar_year+1/12)',__LINE__)
    endif
    iyear=(yeari-int(year0))*12
    if(iyear<0 .or. iyear+12>ntime)then
      write(msg,'(A,I0,A,I0,A,I0,A)')subname//': requested year ',yeari, &
        ' is outside the full-year coverage ',int(year0),'-',nint(yearN)-1,' of'
      call endrun(trim(msg)//' '//trim(atm_ghg_in),__LINE__)
    endif
  endif

  if(atm_co2_fix>0._r8)then
    atm_co2_mon=atm_co2_fix
  else
    DO k=1,12
      call ncd_getvar(atm_ghg_nfid,'CO2',iyear+k,atm_co2_mon(k))
    enddo
  endif  
  if(atm_ch4_fix>0._r8)then
    atm_ch4_mon=atm_ch4_fix
  else
    do K=1,12
      call ncd_getvar(atm_ghg_nfid,'CH4',iyear+k,atm_ch4_mon(k))
    enddo
  endif
  if(atm_n2o_fix>0._r8)then
    atm_n2o_mon=atm_n2o_fix
  else
    do k=1,12
      call ncd_getvar(atm_ghg_nfid,'N2O',iyear+k,atm_n2o_mon(k))
    enddo
  endif

  call ncd_pio_closefile(atm_ghg_nfid)
  
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      CO2EI_col(NY,NX)=atm_co2_mon(1)
      CO2E_col(NY,NX) =CO2EI_col(NY,NX)
      CH4E_col(NY,NX) =atm_ch4_mon(1)*1.e-3_r8  !ppb to ppm
      Z2OE_col(NY,NX) =atm_n2o_mon(1)*1.e-3_r8  !ppb to ppm
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine ReadAtmosComposition


!------------------------------------------------------------------------------------------
  subroutine get_clm_years()
  !Description  
  !obtain the forcing year
  use EcoSIMCtrlMod, only : clm_hour_file_in      !file for hourly climate forcing
  use EcoSIMCtrlMod, only : clm_day_file_in       !file for daily climate forcing
  use EcoSIMCtrlMod, only : yearf1,yearf2, nyeardal1             
  use fileUtil, only : file_exists  
  implicit none
  type(file_desc_t) :: clm_nfid

  yearf1=-999;nyeardal1=0
  if(len_trim(clm_day_file_in)/=0)then
    if(.not.file_exists(clm_day_file_in))then
      call endrun("Do not find clm_day_file_in file "//trim(clm_day_file_in)//" in "//mod_filename, __LINE__)
    endif

    call ncd_pio_openfile(clm_nfid, clm_day_file_in, ncd_nowrite)
    nyeardal1=get_dim_len(clm_nfid, 'year')
    
    call ncd_getvar(clm_nfid,'year',1,yearf1)        
    call ncd_pio_closefile(clm_nfid)
  endif

  yearf2=-999
  if(len_trim(clm_hour_file_in)/=0)then
    if(.not.file_exists(clm_hour_file_in))then
      call endrun("Do not find clm_hour_file_in file "//trim(clm_hour_file_in)//" in "//mod_filename, __LINE__)
    endif

    call ncd_pio_openfile(clm_nfid, clm_hour_file_in, ncd_nowrite)
    
    call ncd_getvar(clm_nfid,'year',1,yearf2)        

    call ncd_pio_closefile(clm_nfid)
  endif
  end subroutine get_clm_years
!-----------------------------------------------------------------------
  function get_forc_step_type(yeari)result(iclmtype)
  !
  !  yearf1<=year<yearf2, daily climate forcing
  !  yearf2<year, hourly climate forcing
  implicit none
  integer, intent(in) :: yeari  ! current year of forcing data
  integer :: iclmtype

  iclmtype=2
  if(nyeardal1>0)then
    !there is daily data provided
    if(yeari<yearf1)then
      call endrun("Do not find requested daily climate data in "//mod_filename, __LINE__)
    elseif(yeari<yearf2)then
      iclmtype=1
    else
      iclmtype=2
    endif  
    return
  endif

  if(yeari<yearf2)then
    call endrun("Do not find requested hourly climate data in "//mod_filename, __LINE__)
  endif
  return
  end function get_forc_step_type
!-----------------------------------------------------------------------

  subroutine read_soil_warming_Tref(year,NHW,NHE,NVN,NVS)
  !
  !Description
  !read soil warming file
  use abortutils,      only: destroy
  use EcosysWarmingMod, only: get_warming_fname, warming_hours_needed
  use GridMod,         only: get_col
  implicit none
  integer, intent(in) :: year
  integer, intent(in) :: NHW,NHE,NVN,NVS

  character(len=256) :: fname
  type(file_desc_t) :: ncid
  type(Var_desc_t) :: vardesc   ! Output variable descriptor
  logical           :: readvar   ! If variable exists or not
  integer           :: nsteps,nlevs,ncol
  integer           :: kk, icol, ny,nx, ll
  real(r8), allocatable    :: data(:,:)  
  
  if(.not.get_warming_fname(year,fname))return
  
  call ncd_pio_openfile (ncid, trim(fname), ncd_nowrite)
  call check_var(ncid, 'TEMP_vr', vardesc, readvar, print_err=.true.)
  nsteps = get_dim_len(ncid,'time')
  nlevs  = get_dim_len(ncid,'levsoi')
  ncol   = get_dim_len(ncid,'column')
  allocate(data(ncol,nlevs))

  !hours missing from the file would leave zeros or the previous year's values in TKS_ref_vr
  if(nsteps<warming_hours_needed(year))then
    write(*,*)'read_soil_warming_Tref: ',trim(fname),' has ',nsteps,' hourly records, needs ', &
      warming_hours_needed(year),' for year ',year
    call endrun('read_soil_warming_Tref: too few hourly records in '//trim(fname),__LINE__)
  endif

  !write(*,*)nsteps,nlevs,ncol
  do kk = 1,nsteps
    call ncd_getvar(ncid,'TEMP_vr', kk, data)

    DO NX=NHW,NHE
      DO NY=NVN,NVS
        icol=get_col(NY,NX)
        TKS_ref_vr(kk,1:nlevs,NY,NX)=data(icol,1:nlevs)+TC2K
      ENDDO
    ENDDO  
!    if(mod(kk,24*30)==0)then
!      print*,kk/24,TKS_ref_vr(kk,1:3,1,1),data(1,1:3)
!    endif
  enddo
  call ncd_pio_closefile(ncid)

  call destroy(data)

  end subroutine read_soil_warming_Tref

end module ClimReadMod
