module VanGenuchtenSoilHydraulicMod
  !
  ! Description:
  ! Implement van Genuchten retention paired with Mualem conductivity, m=(n-1)/n.
  ! Alpha is stored in MPa-1; theta_s is the current micropore-domain porosity.
  ! Reference: van Genuchten (1980), doi:10.2136/sssaj1980.03615995004400050002x.
  !
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use abortutils, only: endrun,iulog
  use DebugToolMod, only: PrintInfo
  use SoilPhysDataType, only: VGThetaR_vr,VGAlpha_vr,VGN_vr,VGL_vr
  use SoilHydraulicModelMod, only: SoilHydraulicModel,SoilHydraulicParameters,NumSoilMoistureClasses
  implicit none
  private
  type, extends(SoilHydraulicModel), public :: VanGenuchtenSoilHydraulicModel
    ! Description:
    ! Use spatial VG parameters stored with soil physical state; surface domains
    ! continue to use their separate implementations in the model adapters.
  contains
    procedure :: matric_potential      => VGMatricPotential
    procedure :: water_content         => VGWaterContent
    procedure :: conductivity_table    => VGConductivityTable
    procedure :: relative_conductivity => VGRelativeConductivity
  end type VanGenuchtenSoilHydraulicModel
contains
!------------------------------------------------------------------------------------------
  subroutine GetVGCoefficients(parameters,theta_r,alpha,n,ell,m)
    !
    ! Description:
    ! Read and validate one ordinary-soil layer's coefficients before logarithms.
    !
    type(SoilHydraulicParameters), intent(in) :: parameters ! layer indices and current porosity [-, m3 m-3]
    real(r8), intent(out) :: theta_r                        ! residual liquid content [m3 m-3 micropore soil]
    real(r8), intent(out) :: alpha                          ! inverse pressure scale [MPa-1]
    real(r8), intent(out) :: n,ell,m                        ! shape, connectivity, and constrained m exponents [-]
    character(len=*), parameter :: subname='GetVGCoefficients'

    call PrintInfo('beg '//subname)

    if(.not.allocated(VGThetaR_vr) .or. .not.allocated(VGAlpha_vr).or. &
       .not.allocated(VGN_vr) .or. .not.allocated(VGL_vr))then
       call endrun('Van Genuchten parameters are not initialized')
    endif
    associate(                   &
      L   => parameters%layer,   &
      NY  => parameters%row,     &
      NX  => parameters%column,  &
      phi => parameters%porosity &
    )
    if(L<1 .or. L>ubound(VGThetaR_vr,1) .or. NY<1 .or. NY>size(VGThetaR_vr,2) .or. &
        NX<1.or.NX>size(VGThetaR_vr,3)) call endrun('Invalid van Genuchten parameter indices')
    theta_r = VGThetaR_vr(L,NY,NX)
    alpha   = VGAlpha_vr(L,NY,NX)
    n       = VGN_vr(L,NY,NX)
    ell     = VGL_vr(L,NY,NX)
    if(.not.all(ieee_is_finite([theta_r,alpha,n,ell,phi])))then
      write(iulog,*) 'Nonfinite VG parameter at L,NY,NX=',L,NY,NX
      call endrun('Invalid van Genuchten parameters')
    endif
    if(theta_r<0._r8.or.theta_r>=phi.or.phi<=0._r8.or.phi>1._r8.or.alpha<=0._r8.or.n<=1._r8)then
      write(iulog,*) 'VG theta_r,alpha,n,l,porosity at L,NY,NX=',L,NY,NX,theta_r,alpha,n,ell,phi
      call endrun('Invalid van Genuchten parameters')
    endif
    m=(n-1._r8)/n
    if(m<=0._r8.or.m>=1._r8)call endrun('Van Genuchten n is not numerically resolvable')

    ! This bound admits fitted negative l while retaining a vanishing, monotone dry K.
    if(ell<=-2._r8/m)call endrun('Van Genuchten l must exceed -2/m')
    end associate
    call PrintInfo('end '//subname)
  end subroutine GetVGCoefficients
!------------------------------------------------------------------------------------------
  function VGMatricPotential(this,parameters,theta) result(potential)
    ! Description:
    ! Evaluate negative matric pressure from effective saturation. The model's dry
    ! pressure bound is applied to the diagnostic, without changing water inventory.
    class(VanGenuchtenSoilHydraulicModel), intent(in) :: this ! selected constitutive implementation [-]
    type(SoilHydraulicParameters), intent(in) :: parameters   ! one-layer context
    real(r8), intent(inout) :: theta                          ! liquid content; limited to porosity when saturated [m3 m-3]
    real(r8) :: potential                                     ! matric potential [MPa]
    real(r8) :: theta_r,alpha,n,ell,m,log_se,log_pressure
    character(len=*), parameter :: subname='VGMatricPotential'

    call PrintInfo('beg '//subname)
    associate(                              &
      phi       => parameters%porosity,     &
      dry_limit => parameters%dry_potential &
      )
    call GetVGCoefficients(parameters,theta_r,alpha,n,ell,m)

    if(theta>=phi)then
      theta=phi;potential=0._r8
    elseif(theta<=theta_r)then
      potential=dry_limit
    else
      if(theta>theta_r+0.5_r8*(phi-theta_r))then
        log_se=VGLog1p(-(phi-theta)/(phi-theta_r))
      else
        log_se=LOG((theta-theta_r)/(phi-theta_r))
      endif
      log_pressure=VGLogExpm1(-log_se/m)/n-LOG(alpha)
      if(log_pressure>=LOG(-dry_limit))then
        potential=dry_limit
      else
        potential=-EXP(log_pressure)
      endif
    endif
    end associate
    call PrintInfo('end '//subname)
  end function VGMatricPotential
!------------------------------------------------------------------------------------------
  function VGWaterContent(this,parameters,potential) result(theta)
    ! Description:
    ! Evaluate theta_r+(porosity-theta_r)*(1+(alpha*abs(psi))**n)**(-m),
    ! using log-domain arithmetic and the existing model dry-pressure bound.
    class(VanGenuchtenSoilHydraulicModel), intent(in) :: this ! selected constitutive implementation [-]
    type(SoilHydraulicParameters), intent(in) :: parameters   ! one-layer context
    real(r8), intent(in) :: potential                         ! matric potential [MPa]
    real(r8) :: theta                                         ! liquid content [m3 m-3 micropore soil]
    real(r8) :: theta_r,alpha,n,ell,m,x,log_se
    character(len=*), parameter :: subname='VGWaterContent'

    call PrintInfo('beg '//subname)
    call GetVGCoefficients(parameters,theta_r,alpha,n,ell,m)
    associate(phi => parameters%porosity,dry_limit => parameters%dry_potential)
      if(potential>=0._r8)then
        theta=phi
      else
        x=n*(LOG(alpha)+LOG(-AMAX1(potential,dry_limit)))
        if(x>0._r8)then
          log_se=-m*(x+VGLog1p(EXP(-x)))
        else
          log_se=-m*VGLog1p(EXP(x))
        endif
        theta=theta_r+(phi-theta_r)*EXP(log_se)
      endif
    end associate
    call PrintInfo('end '//subname)
  end function VGWaterContent
!------------------------------------------------------------------------------------------
  function VGRelativeConductivity(this,parameters,theta) result(relative_k)
    ! Description:
    ! Evaluate analytic Mualem K/Ksat=Se**l*[1-(1-Se**(1/m))**m]**2.
    ! Log-domain/asymptotic evaluation avoids cancellation near residual water.
    class(VanGenuchtenSoilHydraulicModel), intent(in) :: this ! selected constitutive implementation [-]
    type(SoilHydraulicParameters), intent(in) :: parameters   ! one-layer context
    real(r8), intent(in) :: theta                             ! liquid content [m3 m-3 micropore soil]
    real(r8) :: relative_k                                    ! analytic relative conductivity [-]
    real(r8) :: theta_r,alpha,n,ell,m,log_se,log_q,log_bracket,log_k
    logical :: underflows
    character(len=*), parameter :: subname='VGRelativeConductivity'

    call PrintInfo('beg '//subname)

    associate(phi => parameters%porosity)
    call GetVGCoefficients(parameters,theta_r,alpha,n,ell,m)
    if(theta>=phi)then
      relative_k=1._r8
    elseif(theta<=theta_r)then
      relative_k=0._r8
    else
      if(theta>theta_r+0.5_r8*(phi-theta_r))then
        log_se=VGLog1p(-(phi-theta)/(phi-theta_r))
      else
        log_se=LOG((theta-theta_r)/(phi-theta_r))
      endif
      if(log_se>=0._r8)then
        relative_k=1._r8
      else
        underflows=.false.
        if(ell>0._r8)underflows=-log_se>(-LOG(tiny(1._r8)))/ell
        if(underflows)then
          relative_k=0._r8
        else
          log_q=log_se/m
          if(log_q < -36._r8)then
            log_k=(ell+2._r8/m)*log_se+2._r8*LOG(m)
          else
            log_bracket=VGLogOneMinusExp(m*VGLogOneMinusExp(log_q))
            log_k=ell*log_se+2._r8*log_bracket
          endif
          if(log_k<LOG(tiny(1._r8)))then
            relative_k=0._r8
          else
            relative_k=EXP(AMIN1(0._r8,log_k))
          endif
        endif
      endif
    endif
    end associate
    call PrintInfo('end '//subname)
  end function VGRelativeConductivity
!------------------------------------------------------------------------------------------
  subroutine VGConductivityTable(this,parameters,ksat_vertical,ksat_horizontal,water_classes,pressure_classes,conductivity)
    ! Description:
    ! Sample the paired VG/Mualem functions on the existing total-saturation grid.
    ! Runtime indexing remains based on theta/porosity, not effective saturation.
    class(VanGenuchtenSoilHydraulicModel), intent(in) :: this         ! selected constitutive implementation [-]
    type(SoilHydraulicParameters), intent(in) :: parameters           ! one-layer context
    real(r8), intent(in) :: ksat_vertical,ksat_horizontal             ! saturated coefficients [m2 MPa-1 h-1]
    real(r8), intent(out) :: water_classes(NumSoilMoistureClasses)    ! liquid samples [m3 m-3]
    real(r8), intent(out) :: pressure_classes(NumSoilMoistureClasses) ! corresponding potentials [MPa]
    real(r8), intent(out) :: conductivity(3,NumSoilMoistureClasses)   ! horizontal 1:2, vertical 3 [m2 MPa-1 h-1]
    real(r8) :: relative_k
    integer :: k
    character(len=*), parameter :: subname='VGConductivityTable'

    call PrintInfo('beg '//subname)
    associate(phi => parameters%porosity)
      do k=1,NumSoilMoistureClasses
        water_classes(k)    = phi*(1._r8-(k-1._r8)/NumSoilMoistureClasses)
        pressure_classes(k) = this%matric_potential(parameters,water_classes(k))
        relative_k          = this%relative_conductivity(parameters,water_classes(k))
        conductivity(1:2,k) = ksat_horizontal*relative_k
        conductivity(3,k)   = ksat_vertical*relative_k
      enddo
    end associate
    call PrintInfo('end '//subname)
  end subroutine VGConductivityTable
!------------------------------------------------------------------------------------------
  function VGLog1p(x) result(value)
    ! Description:
    ! Evaluate log(1+x) accurately for small x, with x greater than -1.
    real(r8), intent(in) :: x  ! dimensionless logarithm increment [-]
    real(r8) :: value          ! log(1+x) [-]
    character(len=*), parameter :: subname='VGLog1p'

  call PrintInfo('beg '//subname)
  if(ABS(x)<1.e-4_r8)then
    value=x*(1._r8-x/2._r8+x*x/3._r8-x**3/4._r8+x**4/5._r8)
  else
    value=LOG(1._r8+x)
  endif
  call PrintInfo('end '//subname)
  end function VGLog1p
!------------------------------------------------------------------------------------------
  function VGExpm1(x) result(value)
    ! Description:
    ! Evaluate exp(x)-1 without losing small increments.
    real(r8), intent(in) :: x ! dimensionless exponent [-]
    real(r8) :: value         ! exp(x)-1 [-]
    character(len=*), parameter :: subname='VGExpm1'

  call PrintInfo('beg '//subname)
  if(ABS(x)<1.e-4_r8)then
    value=x*(1._r8+x/2._r8+x*x/6._r8+x**3/24._r8+x**4/120._r8)
  else
    value=EXP(x)-1._r8
  endif
  call PrintInfo('end '//subname)
  end function VGExpm1
!------------------------------------------------------------------------------------------
  function VGLogOneMinusExp(x) result(value)
    ! Description:
    ! Evaluate log(1-exp(x)) for negative x without cancellation near zero.
    real(r8), intent(in) :: x ! negative dimensionless exponent [-]
    real(r8) :: value ! log(1-exp(x)) [-]
    character(len=*), parameter :: subname='VGLogOneMinusExp'
    call PrintInfo('beg '//subname)
    if(x < -LOG(2._r8))then
      value=VGLog1p(-EXP(x))
    else
      value=LOG(-VGExpm1(x))
    endif
    call PrintInfo('end '//subname)
  end function VGLogOneMinusExp
!------------------------------------------------------------------------------------------
  function VGLogExpm1(x) result(value)
    ! Description:
    ! Evaluate log(exp(x)-1) for positive x without overflowing exp(x).
    real(r8), intent(in) :: x ! positive dimensionless exponent [-]
    real(r8) :: value         ! log(exp(x)-1) [-]
    character(len=*), parameter :: subname='VGLogExpm1'

  call PrintInfo('beg '//subname)
  if(x>LOG(2._r8))then
    value=x+VGLog1p(-EXP(-x))
  else
    value=LOG(VGExpm1(x))
  endif
  call PrintInfo('end '//subname)
end function VGLogExpm1
end module VanGenuchtenSoilHydraulicMod
