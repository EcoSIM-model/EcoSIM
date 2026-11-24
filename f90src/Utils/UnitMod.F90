module UnitMod
!
!DESCRIPTION
! code for unit conversion
  use data_kind_mod , only : r8 => DAT_KIND_R8
  use data_const_mod, only : Tref => DAT_CONST_TKFRZ,spval => DAT_CONST_SPVAL
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  type, public :: unit_type    
    real(r8) :: ppmv    
    real(r8) :: ppbv
    real(r8) :: pptv
    real(r8) :: gram2kg
    real(r8) :: gram2Tg
    real(r8) :: gram2Pg
    real(r8) :: gram2ton
    real(r8) :: gram2Mg
    real(r8) :: cm3tom3
    real(r8) :: Pa2kPa    
    real(r8) :: Pa2MPa
    real(r8) :: Pa2Atmos
    real(r8) :: msq2hectare
    real(r8) :: acre2hectare
    real(r8) :: day2seconds
    real(r8) :: hour2seconds
    real(r8) :: molpm3tomolar    !mole per m3 to mol/L
    real(r8) :: km2mile
    real(r8) :: calorie2joule
  contains
    procedure, public :: Initailize => Init
    procedure, public :: Celcius2Kelvin
    procedure, public :: Kelvin2Celcius
    procedure, public :: Fahrenheit2Celcius
    procedure, public :: Celcius2Fahrenheit
    procedure, public :: get_SecondsPerDay
  end type unit_type

  type(unit_type), public, target :: units
contains

  subroutine Init(this)
  implicit none
  class(unit_type) :: this

  this%ppmv        = 1.e-6_r8
  this%ppbv        = 1.e-9_r8
  this%pptv        = 1.e-12_r8
  this%gram2kg     =1.e-3_r8
  this%Pa2kPa      =1.e-3_r8
  this%Pa2MPa      =1.e-3_r8
  this%gram2Mg     =1.e-6_r8
  this%gram2Pg     =1.e-15_r8
  this%gram2Tg     =1.e-12_r8
  this%gram2ton    =1.e-6_r8
  this%Pa2Atmos    =9.8692e-6_r8
  this%msq2hectare =1.e-4_r8
  this%acre2hectare=0.404686_r8
  this%day2seconds =86400._r8
  this%hour2seconds=3600._r8
  this%molpm3tomolar=1.e-3_r8
  this%km2mile     =0.62137119223733_r8
  this%calorie2joule=4.184_r8
  this%cm3tom3     =1.e-6_r8
  end subroutine Init

!------------------------------------------------------------------------
  elemental function Celcius2Kelvin(this,TC)result(TK)
  implicit none
  class(unit_type), intent(in) :: this
  real(r8), intent(in) :: TC
  real(r8) :: TK

  if(TC/=spval)then
    TK=TC+Tref
  else
    TK=spval
  endif

  end function Celcius2Kelvin
!------------------------------------------------------------------------
  elemental function get_SecondsPerDay(this)result(ans)
  implicit none
  class(unit_type), intent(in) :: this
  real(r8) :: ans

  ans = this%day2seconds

  end function get_SecondsPerDay
!------------------------------------------------------------------------
  elemental function Kelvin2Celcius(this,TK)result(TC)
  implicit none
  class(unit_type), intent(in) :: this
  real(r8), intent(in) :: TK
  real(r8) :: TC

  if(TK/=spval)then
    TC=TK-Tref
  else
    TC=spval
  endif
  
  end function Kelvin2Celcius
!------------------------------------------------------------------------

  elemental function Fahrenheit2Celcius(this,F)result(C)
  implicit none
  class(unit_type), intent(in) :: this
  real(r8), intent(in) :: f
  real(r8) :: C
  real(r8), parameter :: offset=32._r8

  C = 5._r8/9._r8*(F - offset)

  end  function Fahrenheit2Celcius
!------------------------------------------------------------------------

  elemental function Celcius2Fahrenheit(this,C)result(F)
  implicit none
  class(unit_type), intent(in) :: this
  real(r8), intent(in) :: C
  real(r8) :: f
  real(r8), parameter :: offset=32._r8

  F = 9._r8/5._r8*C+offset

  end  function Celcius2Fahrenheit
end module UnitMod  
