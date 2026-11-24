module MiniFuncMod
  ! Description:
  ! Supply shared physical response functions and safe numerical combinations.

  use data_kind_mod     , only : r8 => DAT_KIND_R8
  use EcosimConst
  use MiniMathMod
implicit none

  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  
  contains
!------------------------------------------------------------------------------------------

  pure function FilmThickness(PSISM, is_top_layer)result(ans)
  ! Description:
  ! Evaluate the empirical water-film thickness. Limit its suction argument at
  ! the former saturated-soil value PSIPS, retaining a finite saturated film
  ! when a retention scheme reaches zero pressure. Soil pressure is unchanged.
  implicit none
  real(r8), intent(in) :: PSISM ! matric potential [MPa]
  logical, optional, intent(in) :: is_top_layer ! select the surface-film fit [-]
  real(r8) :: ans ! water-film thickness [m]
  real(r8) :: psi_film
  logical :: is_top_layer_loc
  is_top_layer_loc=.false.
  if(PRESENT(is_top_layer))is_top_layer_loc=is_top_layer
  psi_film=AMIN1(PSISM,PSIPS)
  if(is_top_layer_loc)then
    ans=AMAX1(ppmc,EXP(-13.650_r8-0.857_r8*LOG(-psi_film)))
  else
    ans=AMAX1(ppmc,EXP(-13.833_r8-0.857_r8*LOG(-psi_film)))
  endif
  end function FilmThickness
!------------------------------------------------------------------------------------------
  pure function Harmonicmean_safe(a,b)result(hm)
  implicit none
  real(r8), intent(in) :: a,b
  real(r8) :: hm
  real(r8) :: tiny_val
  tiny_val = TINY(1.0_r8)

  IF (a <= tiny_val .OR. b <= tiny_val) THEN
    hm = 0.0_r8
  ELSE  
    hm = 2.0_r8 * MIN(a,b) / (1.0_r8 + MIN(a,b) / MAX(a,b))  
  endif
  end function Harmonicmean_safe
!------------------------------------------------------------------------------------------


  pure function TEFAQUDIF(TK)result(ans)
!
! temperautre effect on aqueous diffusivity
  implicit none
  real(r8), intent(in) :: TK !temperature, [Kelvin]
  real(r8) :: ans

! In general, TK should never be less than zero, however, the
! following avoids segmentatio fault
  if(TK>0._r8)then
    ans =(TK/298.15_r8)**6._r8
  else
    ans=1._r8
  endif
  end function TEFAQUDIF
!------------------------------------------------------------------------------------------


  pure function TEFGASDIF(TK)result(ans)
!
! temperautre effect on gaseous diffusivity
  implicit none
  real(r8), intent(in) :: TK !temperature, [Kelvin]
  real(r8) :: ans

! in general, TK should never be less than zero, however, the
! following avoids segmentatio fault
  if(TK>0._r8)then
    ans =(TK/298.15_r8)**1.75_r8
  else
    ans=1._r8
  endif

  end function TEFGASDIF
!------------------------------------------------------------------------------------------

  pure function TortMicporeW(THETWT)result(ans)
!
! aqueous tortuosity in micropore
  implicit none
  real(r8), intent(in) :: THETWT     !relative saturation of micropore
  real(r8) :: ans
  ans =0.7_r8*THETWT**2._r8
  end function TortMicporeW
!------------------------------------------------------------------------------------------
  pure function TortMacporeW(THETWH)result(ans)
  implicit none
  real(r8), intent(in) :: THETWH   !realative saturation of macropore

  real(r8) :: ans

  ans=AMIN1(1.0_r8,2.8_r8*THETWH**3_r8)

  end function TortMacporeW

!------------------------------------------------------------------------------------------
  pure function fDiffusivitySolutEff(tempscalar,THETWA,Z3SR,is_litter)result(ans)
!
! Description
! compute coefficient for air-water gas transfer
!
  implicit none
  real(r8), intent(in) :: tempscalar           !rate scalar including temperature effect
  real(r8), intent(in) :: THETWA               !relative saturation
  real(r8), intent(in) :: Z3SR                 !saturation threshold
  logical, optional, intent(in) :: is_litter   !indicator if it is litter layer
  real(r8) :: ans
  real(r8), parameter :: Z1S=0.010_r8
  real(r8), parameter :: Z2SW=12.0_r8
  real(r8), parameter :: Z2SD=12.0_r8
  real(r8), parameter :: Z3SX=0.50_r8
  real(r8), parameter :: Z1R=0.01_r8
  real(r8), parameter :: Z2RW=12.0_r8
  real(r8), parameter :: Z2RD=12.0_r8
  real(r8), parameter :: Z3R=0.50_r8
  real(r8) :: Z3S
  logical :: is_litter_loc
! Z1S,Z2SW,Z2SD,Z3SX=parameters for air-water gas transfers in soil
! Z1R,Z2RW,Z2RD,Z3RX=parameters for air-water gas transfers in litter

  is_litter_loc=.false.
  if(present(is_litter))is_litter_loc=is_litter
  if(is_litter_loc)then
    IF(THETWA.GT.Z3R)THEN
      ans=AMAX1(0.0_r8,tempscalar/((Z1R**(-1._r8))*EXP(Z2RW*(THETWA-Z3R))))
    ELSE
      ans=AMIN1(1.0_r8,tempscalar/((Z1R**(-1._r8))*EXP(Z2RD*(THETWA-Z3R))))
    ENDIF
  else
    Z3S=AMAX1(Z3SX,Z3SR)
    IF(THETWA.GT.Z3S)THEN
      ans=AMAX1(0.0_r8,tempscalar/((Z1S**(-1._r8))*EXP(Z2SW*(THETWA-Z3S))))
    ELSE
      ans=AMIN1(1.0_r8,tempscalar/((Z1S**(-1._r8))*EXP(Z2SD*(THETWA-Z3S))))
    ENDIF
  endif
  end function fDiffusivitySolutEff
!------------------------------------------------------------------------------------------
  pure function f_microbe_TOffset(atcs)result(ans)
  !
  !computing the temeprature offset for microbial activity
  implicit none
  real(r8), intent(in) :: atcs   !mean annual soil temperature, [oC]
  real(r8) :: ans

  ans=0.333_r8*(12.5_r8-AMAX1(0.0_r8,AMIN1(25.0_r8,ATCS)))

  end function f_microbe_TOffset

!------------------------------------------------------------------------------------------
  function GetDayLength(ALAT,I,DECLIN)result(DYL)
! Description:
! CALCULATE MAXIMUM DAYLENTH FOR PLANT PHENOLOGY
  implicit none
  real(r8), intent(in) :: ALAT   !latitude in degree
  integer, intent(in) :: I      !day of year
  real(r8), optional, intent(out) :: DECLIN
  real(r8) :: DYL
  real(r8) :: AZI,DEC,DECLIN1
  !the following calculation makes Sep 22 and Mar 23 as the 
  !days that has solar declination close to zero.
  !0.9863=360./365

  DECLIN1=get_sun_declin(I)
  AZI=SIN(ALAT*RadianPerDegree)*SIN(DECLIN1*RadianPerDegree)
  DEC=COS(ALAT*RadianPerDegree)*COS(DECLIN1*RadianPerDegree)
  IF(AZI/DEC.GE.1.0_r8-TWILGT)THEN
    DYL=24.0_r8
  ELSEIF(AZI/DEC.LE.-1.0_r8+TWILGT)THEN
    DYL=0.0_r8
  ELSE
    DYL=12.0_r8*(1.0_r8+2.0_r8/PICON*ASIN(TWILGT+AZI/DEC))
  ENDIF
  if(present(DECLIN))DECLIN=DECLIN1
  end function GetDayLength
!------------------------------------------------------------------------------------------

  function get_sun_declin(I)result(DECLIN)
  implicit none
  integer, intent(in) :: I  !day of year
  real(r8)  :: XI     
  real(r8) :: decday
  real(r8) :: DECLIN

  XI=I
  IF(I.EQ.366)XI=365.5_r8
  DECDAY=XI+100._r8
  DECLIN=SIN((DECDAY*0.9863_r8)*RadianPerDegree)*(-23.47_r8)

  end function get_sun_declin

! -------------------------------------------------------------------------
! Function: calculate_equation_of_time
! 
! Approximation of the Equation of Time (EOT) in minutes.
! Formula: B = (360/365) * (d - 81) in degrees
!          EOT = 9.87 * sin(2B) - 7.53 * cos(B) - 1.5 * sin(B)
! -------------------------------------------------------------------------
  pure FUNCTION calculate_equation_of_time(day_of_year,leapday) RESULT(eot)
  implicit none

  INTEGER, INTENT(IN) :: day_of_year
  real(r8), intent(in) :: leapday

  REAL(r8) :: eot  ![h]
  REAL(r8) :: b, b_rad

  ! B calculation (convert to radians for trig functions)
  ! Ensure floating point division is used
  b = (360.0_r8 / (365.0_r8+leapday)) * REAL(day_of_year - 81, r8)
  
  ! Convert degrees to radians
  b_rad = b * (PICON / 180.0_r8)

  ! EOT calculation
  eot = 9.87_r8 * SIN(2.0_r8 * b_rad) - &
        7.53_r8 * COS(b_rad) - &
        1.50_r8 * SIN(b_rad)
  eot = eot/60._r8 
  END FUNCTION calculate_equation_of_time

!------------------------------------------------------------------------------------------

  pure function gOC_to_m3_OM(gram_OC)result(ans)
  !
  !Description
  !convert gram organic C into m3 OM
  implicit none
  real(r8), intent(in) :: gram_OC    ![gC]
  real(r8) :: ans

  !1.e-6*[gC]/([gC gOM-1]*[gOM cm-3])-> [cm3 OM]*1.e-6 -> [m3 OM]
  ans = 1.e-6_r8*AZMAX1(gram_OC)/(OMCMassFrac*DensitySolidOM)

  end function gOC_to_m3_OM  
end module MiniFuncMod
