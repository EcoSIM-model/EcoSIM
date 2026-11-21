module BETResistanceMod
  ! Description:
  ! Cache a native-to-BET cubic in log suction for evaporation resistance only.
  ! Material anchors use a reference isotherm; transient temperature changes only
  ! the BET pressure scale and four Hermite coefficients. Water state is external.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use EcoSimConst, only: RGASC
  use DebugToolMod, only: PrintInfo
  use abortutils, only: endrun,iulog
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  real(r8), parameter, public :: BETResistanceHumidity=.3_r8 ! maximum BET joining RH [-]
  real(r8), parameter, public :: BETReferenceTemperature=293.15_r8 ! reference isotherm [K]
  real(r8), parameter, public :: SoilBETCutoffFCFraction=.6_r8 ! default soil wet join, theta/field_capacity [-]
  real(r8), parameter, public :: LitterBETCutoffFCFraction=.5_r8 ! default litter wet join, theta/field_capacity [-]
  type, public :: BETResistanceParameters
    ! Description:
    ! Derived, reconstructible material coefficients; no prognostic water state.
    logical :: ready=.false.                 ! coefficients initialized [-]
    real(r8) :: theta_wet=0._r8              ! native/cubic joining content [m3 m-3]
    real(r8) :: potential_wet=0._r8          ! native pressure at wet join [MPa]
    real(r8) :: theta_bet=0._r8              ! cubic/BET joining content [m3 m-3]
    real(r8) :: monolayer=0._r8              ! BET monolayer on material bulk-volume basis [m3 m-3]
    real(r8) :: adsorption_c=0._r8           ! adjustable BET C, distinct from density ratio B [-]
    real(r8) :: coefficient(0:3)=0._r8       ! reference cubic in normalized content [ln(MPa)]
  end type BETResistanceParameters
  public :: MakeBETResistanceParameters,BETResistanceState,BETReferencePressure
contains
!------------------------------------------------------------------------------------------
  function BETReferencePressure(potential_wet) result(potential)
    ! Description:
    ! Select RH=0.3 at reference temperature unless the wet cutoff is already too
    ! dry. Then use three times its suction to preserve pressure ordering even at
    ! half the reference temperature. The saturation cutoff is never moved.
    real(r8), intent(in) :: potential_wet ! native wet-end pressure [MPa]
    real(r8) :: potential ! reference BET joining pressure [MPa]
    character(len=*), parameter :: subname='BETReferencePressure'

    call PrintInfo('beg '//subname)
    potential=RGASC*BETReferenceTemperature/18._r8*LOG(BETResistanceHumidity)
    if(potential>=2._r8*potential_wet)potential=3._r8*potential_wet
    call PrintInfo('end '//subname)
  end function BETReferencePressure
!------------------------------------------------------------------------------------------
  function MakeBETResistanceParameters(theta_wet,potential_wet,derivative_wet,theta_bet_candidate,adsorption_c) result(parameters)
    ! Description:
    ! Anchor the cubic at the prescribed native wet content and a BET isotherm.
    ! Start from the native inverse at the reference BET pressure. If necessary,
    ! raise its dry content to make the Hermite log-suction curve monotone for
    ! T>=Tref/2. This is an explicit shape constraint, not an adsorption-data fit.
    ! C and monolayer are material constants; no matching search runs per step.
    real(r8), intent(in) :: theta_wet           ! prescribed native wet join [m3 m-3]
    real(r8), intent(in) :: potential_wet       ! native wet joining pressure [MPa]
    real(r8), intent(in) :: derivative_wet      ! native dpsi/dtheta at wet join [MPa per volume fraction]
    real(r8), intent(in) :: theta_bet_candidate ! native inverse at reference BET joining pressure [m3 m-3]
    real(r8), intent(in) :: adsorption_c        ! configured BET adsorption parameter C [-]
    type(BETResistanceParameters) :: parameters ! cached reference material coefficients
    real(r8) :: thermal,pressure_bet,log_x,x,factor,match,delta,ratio,wet_elasticity
    real(r8) :: theta_bet,monolayer,width,ylo,yhi,dlo,dhi
    character(len=*), parameter :: subname='MakeBETResistanceParameters'

    call PrintInfo('beg '//subname)
    if(.not.all(ieee_is_finite([theta_wet,potential_wet,derivative_wet,theta_bet_candidate,adsorption_c]))) &
      call endrun('Nonfinite cubic BET resistance anchors')
    if(theta_wet<=0._r8.or.potential_wet>=0._r8.or.derivative_wet<0._r8.or.adsorption_c<=1._r8) &
      call endrun('Invalid cubic BET resistance anchors or C (must exceed one)')
    thermal        = RGASC*BETReferenceTemperature/18._r8
    pressure_bet   = BETReferencePressure(potential_wet)
    log_x          = pressure_bet/thermal
    x              = EXP(log_x)
    factor         = (1._r8-x)*(1._r8+(adsorption_c-1._r8)*x)
    match          = factor/(log_x*(1._r8+(adsorption_c-1._r8)*x*x))
    ylo            = LOG(-pressure_bet);yhi = LOG(-potential_wet)
    delta          = ylo-yhi-LOG(2._r8)
    wet_elasticity = -derivative_wet*theta_wet/potential_wet
    ! Normalized Hermite endpoint slopes alpha<=1, beta<=2 are sufficient for
    ! monotonicity. Design at Tref/2; higher temperature increases the secant.
    ratio=(-match)/((-match)+delta)
    if(wet_elasticity>0._r8)ratio=MAX(ratio,1._r8-2._r8*delta/wet_elasticity)
    theta_bet=MAX(theta_wet*ratio,MIN(.95_r8*theta_wet,MAX(0._r8,theta_bet_candidate)))
    if(theta_bet<=0._r8.or.theta_bet>=theta_wet)call endrun('Unresolved cubic BET moisture interval')
    ! In logarithmic form this also works when the chosen joining RH underflows.
    monolayer=EXP(LOG(theta_bet)+LOG(factor)-LOG(adsorption_c)-log_x)
    if(.not.ieee_is_finite(monolayer).or.monolayer<=0._r8)call endrun('Invalid cubic BET monolayer')
    width                    = theta_wet-theta_bet
    dlo                      = match/theta_bet;dhi = derivative_wet/potential_wet
    parameters%theta_wet     = theta_wet
    parameters%potential_wet = potential_wet
    parameters%theta_bet     = theta_bet
    parameters%monolayer     = monolayer
    parameters%adsorption_c  = adsorption_c
    parameters%coefficient   = [ylo,width*dlo,3._r8*(yhi-ylo)-width*(2._r8*dlo+dhi), &
      2._r8*(ylo-yhi)+width*(dlo+dhi)]
    parameters%ready=.true.
    call PrintInfo('end '//subname)
  end function MakeBETResistanceParameters
!------------------------------------------------------------------------------------------
  subroutine BETResistanceState(parameters,theta,temperature,potential,log_derivative,log_theta_over_humidity)
    ! Description:
    ! Evaluate cached cubic/BET retention and stable logarithmic transport factors.
    ! Cubic temperature updates need one log(T/Tref), no root search. The native
    ! wet value and slope remain fixed; the BET endpoint follows current T.
    ! Verify monotonicity analytically using extrema of the derivative quadratic.
    ! Positive contents only; the caller owns the exact-zero liquid limit.
    type(BETResistanceParameters), intent(in) :: parameters ! cached material isotherm
    real(r8), intent(in) :: theta ! positive liquid content [m3 m-3]
    real(r8), intent(in) :: temperature ! current material temperature [K]
    real(r8), intent(out) :: potential ! resistance-only pressure [MPa]
    real(r8), intent(out) :: log_derivative ! ln(dpsi/dtheta), native MPa/volume-fraction units
    real(r8), intent(out) :: log_theta_over_humidity ! ln(theta/RH), finite at positive hyperdry contents
    real(r8) :: y,b,x,log_x,log_y,factor,derivative_factor,log_capacity,thermal
    real(r8) :: c(0:3),shift,t,width,dydt,vertex,largest
    character(len=*), parameter :: subname='BETResistanceState'

    call PrintInfo('beg '//subname)
    associate( &
      material => parameters, & ! input: cached material coefficients [-]
      theta_wet => parameters%theta_wet, & ! input: prescribed wet join [m3 m-3]
      theta_bet => parameters%theta_bet, & ! input: derived BET join [m3 m-3]
      monolayer => parameters%monolayer, & ! input: reference monolayer [m3 m-3]
      adsorption_c => parameters%adsorption_c) ! input: configured adsorption C [-]
    if(.not.material%ready)call endrun('Uninitialized cubic BET material')
    if(.not.all(ieee_is_finite([theta,temperature])).or.theta<=0._r8.or.theta>theta_wet.or.temperature<=0._r8) &
      call endrun('Cubic BET content or temperature outside domain')
    thermal=RGASC*temperature/18._r8
    if(theta>theta_bet)then
      shift=LOG(temperature/BETReferenceTemperature)
      c=material%coefficient+[shift,0._r8,-3._r8*shift,2._r8*shift]
      largest=MAX(c(1),c(1)+2._r8*c(2)+3._r8*c(3))
      if(c(3)<0._r8)then
        vertex=-c(2)/(3._r8*c(3))
        if(vertex>0._r8.and.vertex<1._r8)largest=MAX(largest,c(1)+vertex*(2._r8*c(2)+3._r8*c(3)*vertex))
      endif
      if(largest>1.e-12_r8*MAX(1._r8,MAXVAL(ABS(c))))call endrun('Nonmonotone temperature-updated cubic BET curve')
      width=theta_wet-theta_bet
      t=(theta-theta_bet)/width
      y=c(0)+t*(c(1)+t*(c(2)+t*c(3)))
      dydt=c(1)+t*(2._r8*c(2)+3._r8*c(3)*t)
      potential=-EXP(y)
      log_derivative=-huge(1._r8)
      if(dydt<0._r8)log_derivative=y+LOG(-dydt)-LOG(width)
      log_theta_over_humidity=LOG(theta)-potential/thermal
    else
      log_capacity=LOG(monolayer)+LOG(adsorption_c)
      log_y=LOG(theta)-LOG(monolayer)
      if(log_y<.5_r8*LOG(tiny(1._r8)))then
        x=0._r8
        log_x=LOG(theta)-log_capacity
      else
        y=EXP(log_y)
        b=adsorption_c-(adsorption_c-2._r8)*y
        x=2._r8*y/(b+SQRT(b*b+4._r8*(adsorption_c-1._r8)*y*y))
        log_x=LOG(x)
      endif
      factor=(1._r8-x)*(1._r8+(adsorption_c-1._r8)*x)
      derivative_factor=1._r8
      if(x>=SQRT(tiny(1._r8)))derivative_factor=1._r8+(adsorption_c-1._r8)*x*x
      potential=thermal*log_x
      log_theta_over_humidity=log_capacity-LOG(factor)
      log_derivative=LOG(thermal)+2._r8*LOG(factor)-log_capacity-log_x-LOG(derivative_factor)
    endif
    end associate
    call PrintInfo('end '//subname)
  end subroutine BETResistanceState
end module BETResistanceMod
