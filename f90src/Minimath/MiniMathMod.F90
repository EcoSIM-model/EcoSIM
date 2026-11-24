module minimathmod
!!
! Description:
! Some small subroutines/function to do safe math.

  use data_kind_mod, only: r8 => DAT_KIND_R8
  use abortutils,    only: iulog
  use EcoSimConst

  implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: safe_adb
  public :: p_adb
  public :: isclose         !test if two values a and b are close in magnitude
  public :: vapsat, vapsat0
  public :: isLeap,iisleap
  
  public :: AZMAX1,AZMIN1,AZMAX1t,AZMAX1d,AZMIN1d
  public :: GetMolAirPerm3
  public :: fSiLU
  public :: fixnegmass
  public :: fixEXConsumpFlux
  public :: yearday,isletter
  public :: dssign
  public :: symmetric_flux_limiter
  public :: flux_mass_limiter
  public :: AZERO,AZERO1  
  public :: SubstrateLimit
  public :: real_truncate
  public :: pMod
  public :: isAinsideBC,isABetweenBC,isALeftinBC,isARightinBC
  public :: sfexp
  public :: Viscosity_H2O
  public :: VapMass2KPa
  public :: attn_func
  public :: SubstrateDribbling  
  interface SubstrateDribbling
    module procedure SubstrateDribbling_vec
    module procedure SubstrateDribbling_scal
  end interface SubstrateDribbling
  interface AZMAX1
    module procedure AZMAX1_s
    module procedure AZMAX1_d
    module procedure AZMAX1_v
  end interface AZMAX1

  interface AZMIN1
    module procedure AZMIN1_s
    module procedure AZMIN1_d
  end interface AZMIN1
  
  public :: get_flux_scalar
  public :: addone
  public :: RichardsonNumber
  real(r8), parameter :: tiny_val=1.e-12_r8

  contains

!------------------------------------------------------------------------------------------

   pure function safe_adb(a,b)result(ans)
   !!
   ! Description:
   ! damp division by zero to zero
   implicit none
   real(r8), intent(in) :: a,b
   real(r8) :: ans

   if(abs(b)<tiny_val)then
     !assume complete dissipation
     ans=0._R8
   else
     ans = AZERO(a/b)
   endif

   return
   end function safe_adb
!------------------------------------------------------------------------------------------

   pure function p_adb(a,b)result(ans)
   !!
   ! Description:
   ! ans=max(0.,a/b)
   implicit none
   real(r8), intent(in) :: a,b
   real(r8) :: ans

   ans=AMAX1(0._r8,a/b)
   return
   end function p_adb
!------------------------------------------------------------------------------------------

  function dssign(snow)result(ans)
  implicit none
  real(r8), intent(in) :: snow
  real(r8) :: ans

  if(snow<=1.e-8_r8)then
    ans=0._r8
  else
    ans=1._r8
  endif
  end function dssign
!------------------------------------------------------------------------------------------
  pure function isAinsideBC(a,b,c)result(ans)
  implicit none
  real(r8), intent(in) :: a,b,C
  logical :: ans

  ans = a>b .and. a < c
  end function isAinsideBC
!------------------------------------------------------------------------------------------
  pure function isALeftinBC(a,b,c)result(ans)
  implicit none
  real(r8), intent(in) :: a,b,C
  logical :: ans

  ans = a>=b .and. a < c
  end function isALeftinBC
!------------------------------------------------------------------------------------------
  pure function isARightinBC(a,b,c)result(ans)
  implicit none
  real(r8), intent(in) :: a,b,C
  logical :: ans

  ans = a>b .and. a <= c
  end function isARightinBC
!------------------------------------------------------------------------------------------
  pure function isABetweenBC(a,b,c)result(ans)
  implicit none
  real(r8), intent(in) :: a,b,C
  logical :: ans

  ans = a>=b .and. a <= c
  end function isABetweenBC

!------------------------------------------------------------------------------------------

  pure function vapsat(tempK)result(ans)
  !
  ! Description
  ! compute saturated vapor pressure, based on temperature tempK (in K)
  ! ep=\mu R T= (m/M)/V0 RT, V0=1m3, m=ep*V0*M/(RT), [kPa]*[m3]*[g/mol]/[Pa m3]~kg 
  implicit none
  real(r8), intent(in) :: tempK

  real(r8) :: ans  !ton, i.e. (ans*10^3=kg/m3) in terms vapor concentration, 2.173~18/8.314
  
  ans=2.173E-03_r8/tempK*0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/tempK))
  end function vapsat

!------------------------------------------------------------------------------------------

  pure function vapsat0(tempK)result(ans)
  !
  ! Description
  ! compute saturated vapor pressure, based on temperature tempK (in K)
  implicit none
  real(r8), intent(in) :: tempK

  real(r8) :: ans  !(kPa)
  ans=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/tempK))
  end function vapsat0

!------------------------------------------------------------------------------------------
  pure function attn_func(tau)result(attn)
  implicit none
  real(r8), intent(in) :: tau
  real(r8) :: attn

  if (tau < 1.e-6_r8) then
    attn = 1._r8 - 0.5_r8*tau + tau*tau/6._r8
  else
    attn = (1._r8-sfexp(-tau))/tau
  endif  
  end function attn_func

!------------------------------------------------------------------------------------------

  pure function isLeap(year)result(ans)
!
! Description
! Determine if it is a leap year

  implicit none
  integer, intent(in) :: year
  logical :: ans

  ans =(mod(year,400)== 0) .or. (mod(year,4)==0 .and. mod(year,100)/=0)
  end function isLeap
!------------------------------------------------------------------------------------------

  pure function iisLeap(year)result(ans)
!
! Description
! Determine if it is a leap year

  implicit none
  integer, intent(in) :: year
  integer :: ans

  if (isLeap(year))then
    ans=1
  else
    ans=0
  endif
  end function iisLeap

!------------------------------------------------------------------------------------------

  pure function AZMAX1t(val)result(ans)
  implicit none
  real(r8), intent(in) :: val

  real(r8) :: ans

  ans=AMAX1(val,tiny_val)

  end function AZMAX1t
!------------------------------------------------------------------------------------------

  pure function AZMAX1d(val,tiny_val2)result(ans)
  !
  implicit none
  real(r8), intent(in) :: val
  real(r8), intent(in) :: tiny_val2 !positive tiny threshold value

  real(r8) :: ans

  if(val>tiny_val2)then
    ans=val 
  else
    ans=0._r8
  endif
  end function AZMAX1d

!--------------------------------------------------------------------------------
  pure function pMOD(a,b)result(c)
  !
  ! compute c=MOD(A,B)
  ! if C==0 and A/=0 then C=B
  implicit none
  integer, intent(in) :: a,B
  integer :: C
  
  if(a<0)then
    c=mod(a+b,b)
  else
    c=mod(a,b)
  endif

  if(c==0 .and. a/=0)c=b
  end function pMOD
!------------------------------------------------------------------------------------------

  pure function AZMIN1d(val,tiny_val2)result(ans)
  implicit none
  real(r8), intent(in) :: val
  real(r8), intent(in) :: tiny_val2 !positive tiny threshold value

  real(r8) :: ans

  if(val<-tiny_val2)then
    ans=val 
  else
    ans=0._r8
  endif
  end function AZMIN1d

!------------------------------------------------------------------------------------------
  pure function AZERO(val,tiny1)result(ans)
  implicit none
  real(r8), intent(in) :: val
  real(r8), optional, intent(in) :: tiny1
  real(r8) :: ans
  real(r8), parameter :: tiny_val1=1.e-10_r8
  real(r8) :: tiny  
  if(present(tiny1))then
    tiny=tiny1
  else
    tiny=tiny_val1
  endif  

  if(abs(val)>=tiny)then
    ans=val
  else  
    ans=0._r8
  endif

  end function AZERO

!------------------------------------------------------------------------------------------
  pure function AZERO1(val,tiny_val2)result(ans)
  implicit none
  real(r8), intent(in) :: val
  real(r8), intent(in) :: tiny_val2
  real(r8) :: ans

  if(abs(val)>=tiny_val2)then
    ans=val
  else  
    ans=0._r8
  endif

  end function AZERO1
!------------------------------------------------------------------------------------------

  pure function AZMAX1_s(val)result(ans)
  implicit none
  real(r8), intent(in) :: val

  real(r8) :: ans

  if(val>=tiny_val)then
    ans=val
  else  
    ans=0._r8
  endif

  end function AZMAX1_s  
!------------------------------------------------------------------------------------------

  pure function AZMAX1_d(val1,val2)result(ans)
  implicit none
  real(r8), intent(in) :: val1,val2

  real(r8) :: ans

  ans=AMAX1(0.0_r8,val1,val2)

  end function AZMAX1_d
!------------------------------------------------------------------------------------------

  pure function AZMAX1_v(val1)result(ans)
  implicit none
  real(r8), dimension(:), intent(in) :: val1
  
  real(r8) :: ans(size(val1))
  integer :: jj, sz
  sz = size(val1)
  DO jj=1,sz
    ans(jj)=AMAX1(0.0_r8,val1(jj))
  enddo

  end function AZMAX1_v


!------------------------------------------------------------------------------------------

  pure function AZMIN1_s(val)result(ans)
  implicit none
  real(r8), intent(in) :: val

  real(r8) :: ans

  ans=AMIN1(0.0_r8,val)

  end function AZMIN1_s


!------------------------------------------------------------------------------------------

  pure function AZMIN1_d(val1,val2)result(ans)
  implicit none
  real(r8), intent(in) :: val1,val2

  real(r8) :: ans

  ans=AMIN1(0.0_r8,val1,val2)

  end function AZMIN1_d

! ----------------------------------------------------------------------

  function addone(itemp)result(ans)
!
!  DESCRIPTION
! increase itemp by one
  implicit none
  integer, intent(inout) :: itemp

  integer :: ans

  itemp=itemp+1
  ans=itemp
  end function addone


! ----------------------------------------------------------------------
  pure function isclose(a,b)result(ans)
  !DESCRIPTION
  !determine if a is close to b in magnitude by relative magnitude tiny_val

  implicit none
  real(r8), intent(in) :: a,b
  real(r8) :: c,ac,bc
  logical :: ans
  
  c=max(abs(a),abs(b))  
  if (c<tiny_val) then
    ans=.True.
    return     
  endif

  if(abs(a+b)<1.e-20_r8)then
    ans=.false.
  else  
    ac=a/c;bc=b/c    
    ans=abs((ac-bc)/(ac+bc))<tiny_val
  endif
  end function isclose

! ----------------------------------------------------------------------

  pure function RichardsonNumber(RIB,TK1,TK2)result(ans)
  implicit none
  real(r8), intent(in) :: RIB  !isothermal RI
  real(r8), intent(in) :: TK1  !target temperature
  real(r8), intent(in) :: TK2  !reference temperature

  real(r8) :: ans !>0 stable, < 0 unstable

  ans = AMAX1(-0.3_r8,AMIN1(0.075_r8,RIB*(TK1-TK2)))
  
  end function RichardsonNumber
! ----------------------------------------------------------------------
  pure function GetMolAirPerm3(TKair,Patm_Pa)result(ans)
  implicit none
  real(r8), intent(in) :: TKair
  real(r8), optional, intent(in) :: Patm_Pa  !atmospheric pressure in Pascal

  real(r8) :: ans

  if(present(Patm_Pa))then
    ans = Patm_Pa/(TKair*RGASC)
  else
    ans= 1.01325E5_r8/(TKair*RGASC)
  endif

  end function GetMolAirPerm3
! ----------------------------------------------------------------------
  function fSiLU(x,b)result(ans)
  !Sigmoid linear unit function
  implicit none
  real(r8), intent(in) :: x
  real(r8),optional, intent(in) :: b
  real(r8) :: b_loc
  
  real(r8) :: ans
  
  if(present(b))then
    b_loc=b
  else
    b_loc=1._r8
  endif
  
  ans=x/(1._r8+exp(-b_loc*x))
  
  end function fSiLU
! ----------------------------------------------------------------------
  function fixnegmass(val,refcon)result(ans)
  implicit none
  real(r8), intent(in) :: val
  real(r8), optional, intent(in) :: refcon
  real(r8) :: ans

  ans=val

  if(present(refcon))then
    if(val<0._r8 .and. val>AMIN1(-refcon*1.e-3_r8,-1.e-5))ans=0._r8
  else
    if(val<0._r8 .and. val >-1.e-5_r8)ans=0._r8
  endif

  end function fixnegmass
! ----------------------------------------------------------------------

  subroutine fixEXConsumpFlux(mass,consum_flux,dsgn)
  implicit none
  real(r8), intent(inout) :: mass
  real(r8), intent(inout) :: consum_flux
  integer, optional, intent(in) :: dsgn  !sign of consumption flux
  integer :: dsgnl

  dsgnl=1
  if(present(dsgn))dsgnl=dsgn


  !cut off too small mass
  if(abs(mass)<1.e-12_r8)mass=0._r8    
  !return for zero flux
  if(isclose(consum_flux,0._r8))return

  
  if(dsgnl>0)then  
    !udpate as mass=mass-flux
    if(mass<consum_flux)then
      consum_flux = mass
      mass        = 0._r8
    else
      mass=mass-consum_flux  
    endif
    
  else  
    !update as mass=mass+flux
    if(mass<-consum_flux)then  
      consum_flux = -mass
      mass        = 0._r8
    else
      mass=mass+consum_flux
    endif
  endif
  end subroutine fixEXConsumpFlux

! ----------------------------------------------------------------------

  pure function yearday(year,month,day)result(doy)
  implicit none
  integer, intent(in) :: year
  integer, intent(in) :: month
  integer, intent(in) :: day

  integer, parameter :: daz(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  integer :: doy
  integer :: jj

  doy=0
  do jj = 1, month-1
    doy=doy+daz(jj)
  enddo
  if(month>2)then
    doy=doy+iisleap(year);
  endif
  doy=doy+day
  end function yearday
! ----------------------------------------------------------------------

  pure function isletter(c)result(ans)
  implicit none
  character(len=1), intent(in) :: c
  logical :: ans

  ans=(c>='a' .and. c<='z') .or. (c>='A' .and. c<='Z')

  end function isletter

! ----------------------------------------------------------------------

  pure elemental function symmetric_flux_limiter(flux,reference_flux)result(ans)
  ! Bound flux to [-abs(reference_flux),abs(reference_flux)]. Both arguments
  ! have the same units. The reference sign does not change the bounds;
  ! the flux retains its direction. Donor mass availability is handled separately.
  implicit none
  real(r8), intent(in) :: flux,reference_flux
  real(r8) :: ans

  ans=MIN(ABS(reference_flux),MAX(-ABS(reference_flux),flux))
  end function symmetric_flux_limiter
!------------------------------------------------------------------------------------------

  function flux_mass_limiter(flux,massa,massb)result(ans)
  !
  !limit flux by massa and massb 
  !assuming
  !massa=massa+flux
  !massb=massb-flux
  implicit none
  real(r8), intent(in) :: flux
  real(r8), intent(in) :: massa
  real(r8), intent(in) :: massb
  
  real(r8) :: ans
  if(flux>0._r8)then
    ans=AMIN1(flux,massb)-tiny_val
  else 
    ans=-AMIN1(-flux,massa)+tiny_val
  endif
  end function flux_mass_limiter

! ----------------------------------------------------------------------
  subroutine get_flux_scalar(x0,flux,x1,pscal)
  implicit none
  real(r8), intent(in) :: x0,flux
  real(r8), intent(out):: x1
  real(r8), intent(inout):: pscal
  real(r8) :: p1

  x1=x0
  if(isclose(flux,0._r8))return
  if(pscal>1.e-6_r8)x1=x1+flux  
  if(x1<0._r8)then
    p1=safe_adb(-x0,flux)
    pscal=AMIN1(0.999_r8,pscal,p1) 
    x1=0._r8
  endif

  end subroutine get_flux_scalar
! ----------------------------------------------------------------------

  subroutine SubstrateLimit(n1,n2,demand_flux,y,scal)
  implicit none
  integer, intent(in) :: n1,n2
  real(r8), intent(inout) :: demand_flux(n1:n2)  !consumption/demand flux
  real(r8), intent(inout) :: Y  
  real(r8), optional, intent(out)   :: scal
  real(r8) :: tDemand,scal1
  integer :: N

  tDemand = sum(demand_flux)
  scal1    = 1._r8
  if(y-tDemand>=tiny_val .or. tDemand<=tiny_val)then
    y=AZMAX1(y-tDemand)
  else
    scal1=safe_adb(y,tDemand)
    DO n=n1,n2
      demand_flux(n)=demand_flux(n)*scal1
    enddo
    y=0._r8
  endif
  if(present(scal))scal=scal1
  end subroutine SubstrateLimit

!------------------------------------------------------------------------
  subroutine SubstrateDribbling_scal(demand_flux,dribbling_flx,y)

  implicit none
  real(r8), intent(in) :: demand_flux        !consumption/demand flux
  real(r8), intent(inout) :: dribbling_flx
  real(r8), intent(inout) :: Y
  real(r8) :: tDemand

  tDemand = demand_flux + dribbling_flx

  if(y>tDemand)then
    y=y-tDemand
    dribbling_flx=0._r8
  else
    dribbling_flx=tDemand-y
    y=0._r8
  endif
  
  end subroutine SubstrateDribbling_scal
!------------------------------------------------------------------------
  subroutine SubstrateDribbling_vec(n1,n2,demand_flux,dribbling_flx,y)
  
  !
  !computing the dribbling flux to meet the non-negative requirement of 
  !variable y
  implicit none
  integer, intent(in) :: n1,n2
  real(r8), intent(in) :: demand_flux(n1:n2)  !consumption/demand flux
  real(r8), intent(inout) :: dribbling_flx    !>0
  real(r8), intent(inout) :: Y
  real(r8) :: tDemand

  tDemand = sum(demand_flux) + dribbling_flx

  if(y>tDemand)then
    y=y-tDemand
    dribbling_flx=0._r8
  else
    dribbling_flx=tDemand-y
    y=0._r8
  endif
  
  end subroutine SubstrateDribbling_vec
!------------------------------------------------------------------------

  pure function real_truncate(val,precision)result(ans)
  !
  !truncate val to precision
  implicit none
  real(r8), intent(in) :: val
  real(r8), intent(in) :: precision
  real(r8) :: ans

  ans=aint(val/precision)*precision

  end function real_truncate

!------------------------------------------------------------------------
  pure function Viscosity_H2O(TEMPC)result(viscwl)
  implicit none
  real(r8), intent(in) :: TEMPC  !temperautre [oC]
  real(r8), parameter :: VISCW=1.0E-06_r8    !water viscosity at oC, [Mg m-1 s]
  
  real(r8) :: viscwl
  
  VISCWL    = VISCW*EXP(0.533_r8-0.0267_r8*TEMPC)

  end function Viscosity_H2O
!------------------------------------------------------------------------

  pure function sfexp(a)result(ans)
  !
  !do filtered exponential function calculation
  implicit none
  real(r8), intent(in) :: a
  real(r8) :: ans

  if(a>35._r8)then
    ans=1.e15_r8 !return a large number
  elseif(a<-35._r8)then
    ans=1.e-15_r8
  else
    ans=exp(a)
  endif
  end function sfexp
!------------------------------------------------------------------------
  pure function VapMass2KPa(VP,TA,fwd)result(ans)
  implicit none
  real(r8), intent(in) :: VP ! pressure, trune kPa, false  m3/m3, 
  real(r8), intent(in) :: TA ! temperature K
  logical, intent(in) :: fwd !true, does kPa -> m3/m3; false, m3/m3 -> kPa
  real(r8) :: ans !true m3/m3, false kPa

  if(fwd)then
    ans =VP*2.173E-03_r8/TA
  else
    ans=VP*TA/2.173E-03_r8  
  endif

  end function VapMass2KPa
end module minimathmod
