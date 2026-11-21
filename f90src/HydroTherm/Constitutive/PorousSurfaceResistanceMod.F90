module PorousSurfaceResistanceMod
  ! Description:
  ! Combine parallel vapor diffusion and capillary liquid transport for the
  ! Tang-Riley surface boundary. Material-specific retention stays with callers.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use DebugToolMod, only: PrintInfo
  use abortutils, only: endrun,iulog
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: PorousSurfaceResistance
contains
!------------------------------------------------------------------------------------------
  function PorousSurfaceResistance(theta,thickness,hydraulic_conductivity,slope, &
    bulk_vapor_diffusivity,bunsen_coefficient,log_liquid_prefactor) result(resistance)
    ! Description:
    ! Evaluate r=(d/2)/(D_v_bulk+B*theta*K*dpsi/dtheta) in native hours and MPa.
    ! Gas and liquid conductances act in parallel; air fraction/tortuosity must
    ! already be included once in D_v_bulk. No atmospheric film or flux is added.
    ! Logarithmic sums preserve very dry limits without multiplying large B by K.
    ! An optional ln(B*theta*dpsi/dtheta) preserves BET cancellations before K.
    real(r8), intent(in) :: theta                  ! liquid content per bulk transport volume [m3 m-3]
    real(r8), intent(in) :: thickness              ! full porous control-volume thickness [m]
    real(r8), intent(in) :: hydraulic_conductivity ! pressure-gradient K [m2 MPa-1 h-1]
    real(r8), intent(in) :: slope                  ! dpsi/dtheta of this material [MPa per (m3 m-3)]
    real(r8), intent(in) :: bulk_vapor_diffusivity ! includes pore air/tortuosity [m2 h-1]
    real(r8), intent(in) :: bunsen_coefficient     ! equal-unit liquid/vapor density ratio [-]
    real(r8), optional, intent(in) :: log_liquid_prefactor ! ln(B*theta*dpsi/dtheta), native MPa; excludes K
    real(r8) :: resistance                         ! internal resistance, HUGE if blocked [h m-1]
    real(r8) :: log_liquid,log_vapor,log_transport,log_resistance
    logical :: liquid_transport
    character(len=*), parameter :: subname='PorousSurfaceResistance'

    call PrintInfo('beg '//subname)
    if(.not.all(ieee_is_finite([theta,thickness,hydraulic_conductivity,slope, &
      bulk_vapor_diffusivity,bunsen_coefficient])))then
      write(iulog,*)subname,': nonfinite input'
      call endrun('Invalid porous surface resistance input')
    endif
    if(MIN(theta,thickness,hydraulic_conductivity,slope,bulk_vapor_diffusivity,bunsen_coefficient)<0._r8)then
      write(iulog,*)subname,': negative input'
      call endrun('Invalid porous surface resistance input')
    endif
    resistance=0._r8
    if(thickness>0._r8)then
      liquid_transport=slope>0._r8.and.theta>0._r8.and.hydraulic_conductivity>0._r8.and.bunsen_coefficient>0._r8
      if(present(log_liquid_prefactor))then
        if(.not.ieee_is_finite(log_liquid_prefactor))call endrun('Nonfinite porous resistance log prefactor')
        liquid_transport=theta>0._r8.and.hydraulic_conductivity>0._r8.and.bunsen_coefficient>0._r8
        if(liquid_transport)log_liquid=log_liquid_prefactor+LOG(hydraulic_conductivity)
      elseif(liquid_transport)then
        log_liquid=LOG(bunsen_coefficient)+LOG(theta)+LOG(hydraulic_conductivity)+LOG(slope)
      endif
      if(bulk_vapor_diffusivity>0._r8)then
        log_vapor=LOG(bulk_vapor_diffusivity)
        log_transport=log_vapor
        if(liquid_transport)log_transport=MAX(log_liquid,log_vapor) &
          +LOG(1._r8+EXP(-ABS(log_liquid-log_vapor)))
      elseif(liquid_transport)then
        log_transport=log_liquid
      endif

      resistance=huge(1._r8)
      if(bulk_vapor_diffusivity>0._r8.or.liquid_transport)then
        log_resistance=LOG(thickness)-LOG(2._r8)-log_transport
        if(log_resistance<LOG(tiny(1._r8)))then
          resistance=0._r8
        elseif(log_resistance<LOG(huge(1._r8)))then
          resistance=EXP(log_resistance)
        endif
      endif
    endif
    call PrintInfo('end '//subname)
  end function PorousSurfaceResistance
end module PorousSurfaceResistanceMod
