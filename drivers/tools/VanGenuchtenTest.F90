program VanGenuchtenTest
  ! Description:
  ! Check independent closed-form values, inverse/bounds, spatial coefficients,
  ! the fixed lookup grid, and stable dry/wet evaluation of the VG/Mualem pair.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use abortutils, only: iulog
  use DebugToolMod, only: PrintInfo
  use MiniFuncMod, only: FilmThickness
  use EcoSimConst, only: PSIPS
  use EcoSIMCtrlMod, only: lverb,soil_hydraulic_scheme,VanGenuchtenSchemeName
  use SoilPhysDataType, only: VGThetaR_vr,VGAlpha_vr,VGN_vr,VGL_vr
  use SoilHydraulicModelMod, only: SoilHydraulicModel,SoilHydraulicParameters,NumSoilMoistureClasses
  use SoilHydraulicSchemeMod, only: InitializeSoilHydraulicScheme,FinalizeSoilHydraulicScheme, &
    GetSoilHydraulicModel,GetSoilHydraulicSchemeName
  implicit none
  class(SoilHydraulicModel), pointer :: model
  type(SoilHydraulicParameters) :: parameters,other_layer
  real(r8) :: theta,psi,previous_k,k,back,water(100),pressure(100),conductivity(3,100)
  real(r8), parameter :: ns(6)=[1.01_r8,1.1_r8,1.3_r8,2._r8,4._r8,20._r8]
  real(r8), parameter :: alphas(5)=[1.e-6_r8,1._r8,100._r8,1.e4_r8,1.e200_r8]
  real(r8), parameter :: ses(12)=[1.e-12_r8,1.e-8_r8,1.e-4_r8,.01_r8,.1_r8,.25_r8, &
    .5_r8,.75_r8,.9_r8,.99_r8,1._r8-1.e-8_r8,1._r8-1.e-12_r8]
  character(len=32) :: mode
  integer :: i,j,a

  lverb=.false.
  call get_command_argument(1,mode)
  soil_hydraulic_scheme='van_genuchten_mualem'
  call InitializeSoilHydraulicScheme()
  call require('canonical VG selector',GetSoilHydraulicSchemeName()==VanGenuchtenSchemeName)
  model=>GetSoilHydraulicModel()
  parameters=MakeParameters(1)
  if(mode=='missing_parameters')then
    psi=model%water_content(parameters,-.01_r8)
    error stop 'missing VG parameters were accepted'
  endif
  allocate(VGThetaR_vr(0:2,1,1),VGAlpha_vr(0:2,1,1),VGN_vr(0:2,1,1),VGL_vr(0:2,1,1))
  VGThetaR_vr=.05_r8;VGAlpha_vr=100._r8;VGN_vr=2._r8;VGL_vr=.5_r8
  select case(trim(mode))
  case('invalid_n')
    VGN_vr(1,1,1)=1._r8
  case('invalid_alpha')
    VGAlpha_vr(1,1,1)=0._r8
  case('invalid_residual')
    VGThetaR_vr(1,1,1)=.45_r8
  case('invalid_l')
    VGL_vr(1,1,1)=-4._r8
  end select
  if(index(mode,'invalid_')==1)then
    psi=model%water_content(parameters,-.01_r8)
    error stop 'invalid VG parameters were accepted'
  endif

  theta=.25_r8
  psi=model%matric_potential(parameters,theta)
  call check_close('n=2 pressure reference',psi,-.017320508075688772935_r8)
  call check_close('n=2 Mualem reference',model%relative_conductivity(parameters,theta),.012691995684869118603_r8)
  call check_close('FC reference',model%water_content(parameters,-.033_r8),.16600295811314830881_r8)
  call check_close('WP reference',model%water_content(parameters,-1.5_r8),.052666607409382642893_r8)
  call check_close('dry bound',model%dry_content(parameters),.0500015999999999872_r8)
  call check_close('Ksat reference uses inverse',model%ksat_reference_content(parameters),.16600295811314830881_r8)
  call check_close('equilibrium capacity cap',model%equilibrium_content(parameters,psi,.20_r8),.20_r8)
  theta=1._r8
  psi=model%matric_potential(parameters,theta)
  call require('saturation adjustment',theta==.45_r8.and.psi==0._r8)
  call check_close('saturated surface film',FilmThickness(psi,.true.),FilmThickness(PSIPS,.true.))
  call check_close('saturated subsurface film',FilmThickness(psi),FilmThickness(PSIPS))
  call check_close('unsaturated film unchanged',FilmThickness(-.01_r8),EXP(-13.833_r8-0.857_r8*LOG(.01_r8)))
  theta=.05_r8
  call check_close('residual pressure bound',model%matric_potential(parameters,theta),-2500._r8)
  call require('residual conductivity',model%relative_conductivity(parameters,theta)==0._r8)
  call check_close('nonnegative pressure saturation',model%water_content(parameters,.1_r8),.45_r8)
  call check_close('lower pressure cap',model%water_content(parameters,-1.e30_r8),model%dry_content(parameters))

  call model%conductivity_table(parameters,1._r8,2._r8,water,pressure,conductivity)
  call check_close('class 51 conductivity',conductivity(3,51),.0067181768942885614293_r8)
  call require('dry classes have zero conductivity',all(conductivity(:,90:100)==0._r8))
  do i=1,NumSoilMoistureClasses
    call check_close('total saturation class grid',water(i),.45_r8*(1._r8-(i-1._r8)/100._r8))
    call check_close('horizontal scaling',conductivity(1,i),2._r8*conductivity(3,i))
  enddo

  other_layer=MakeParameters(2)
  VGAlpha_vr(2,1,1)=200._r8
  theta=.25_r8
  call check_close('spatial alpha',model%matric_potential(other_layer,theta),-.008660254037844386468_r8)

  ! Cover varied shape and scale, including pressure clipping and conductivity underflow.
  do j=1,size(ns)
    VGN_vr(1,1,1)=ns(j)
    do a=1,size(alphas)
      VGAlpha_vr(1,1,1)=alphas(a)
      previous_k=0._r8
      do i=1,size(ses)
        theta=.05_r8+.4_r8*ses(i)
        psi=model%matric_potential(parameters,theta)
        back=model%water_content(parameters,psi)
        call check_close('bounded forward/inverse',back,MAX(theta,model%dry_content(parameters)))
        k=model%relative_conductivity(parameters,theta)
        call require('finite monotone bounded K',ieee_is_finite(k).and.k>=previous_k.and.k<=1._r8)
        previous_k=k
        if(mode=='dump')write(*,'(a,*(g0,1x))') 'SAMPLE ',ns(j),alphas(a),.5_r8,theta,psi,back,k
      enddo
    enddo
  enddo
  VGN_vr(1,1,1)=2._r8;VGAlpha_vr(1,1,1)=100._r8;VGL_vr(1,1,1)=-1._r8
  previous_k=0._r8
  do i=1,1000
    theta=.05_r8+.4_r8*i/1000._r8
    k=model%relative_conductivity(parameters,theta)
    call require('negative connectivity remains monotone',ieee_is_finite(k).and.k>=previous_k.and.k<=1._r8)
    previous_k=k
  enddo
  call FinalizeSoilHydraulicScheme()
  write(*,*) 'VanGenuchtenTest passed'
contains
!------------------------------------------------------------------------------------------
  function MakeParameters(layer) result(parameters)
    ! Description:
    ! Construct a scalar context with a known porosity and model pressure bound.
    integer, intent(in) :: layer ! ordinary-soil layer index [-]
    type(SoilHydraulicParameters) :: parameters ! context; unused legacy caches remain valid
    character(len=*), parameter :: subname='MakeParameters'
    call PrintInfo('beg '//subname)
    parameters=SoilHydraulicParameters(porosity=.45_r8,field_capacity=.30_r8,log_porosity=LOG(.45_r8), &
      log_field_capacity=LOG(.30_r8),log_wet_content_ratio=LOG(1.5_r8),log_dry_content_ratio=LOG(2._r8), &
      log_saturation_potential=LOG(.0005_r8),log_field_potential=LOG(.033_r8), &
      log_wet_potential_ratio=LOG(66._r8),log_dry_potential_ratio=LOG(1.5_r8/.033_r8), &
      curvature=1._r8,saturation_potential=0._r8,dry_potential=-2500._r8, &
      saturation_tolerance=1.e-6_r8,layer=layer,row=1,column=1)
    call PrintInfo('end '//subname)
  end function MakeParameters
!------------------------------------------------------------------------------------------
  subroutine check_close(label,actual,expected)
    ! Description:
    ! Compare references and inverse checks with a fixed double-precision allowance.
    character(len=*), intent(in) :: label ! assertion description [-]
    real(r8), intent(in) :: actual,expected ! evaluated and reference values [same units]
    character(len=*), parameter :: subname='check_close'
    call PrintInfo('beg '//subname)
    if(.not.ieee_is_finite(actual).or.ABS(actual-expected)>3.e-12_r8*MAX(ABS(expected),tiny(1._r8)))then
      write(iulog,*) trim(label),': actual=',actual,' expected=',expected
      error stop 1
    endif
    call PrintInfo('end '//subname)
  end subroutine check_close
!------------------------------------------------------------------------------------------
  subroutine require(label,condition)
    ! Description:
    ! Check a physical bound or qualitative property independently of implementation.
    character(len=*), intent(in) :: label ! assertion description [-]
    logical, intent(in) :: condition ! required property [-]
    character(len=*), parameter :: subname='require'
    call PrintInfo('beg '//subname)
    if(.not.condition)then
      write(iulog,*) 'Failed: ',trim(label)
      error stop 1
    endif
    call PrintInfo('end '//subname)
  end subroutine require
end program VanGenuchtenTest
