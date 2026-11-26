program MoistureFactorTest
  ! Description:
  ! Check column parameter lifetime, wet-tail limits and effective K selection.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf,ieee_is_finite
  use GridConsts, only: JX,JY,JZ
  use EcoSIMCtrlMod, only: lverb
  use DebugToolMod, only: PrintInfo
  use SoilPropertyDataType, only: DMeff_col,DMeff_initial,InitSoilProperty,DestructSoilProperty
  use SoilPhysParaMod, only: moisture_factor,getMoistK
  implicit none
  character(len=*), parameter :: mod_filename = &
  __FILE__
  character(len=32) :: mode
  real(r8) :: saturation,f,theta_eff
  integer :: k

  lverb=.false.
  JX=2;JY=2;JZ=2
  call get_command_argument(1,mode)
  select case(trim(mode))
  case('negative')
    DMeff_initial=-.1_r8
  case('nan')
    DMeff_initial=ieee_value(0._r8,ieee_quiet_nan)
  case('infinite')
    DMeff_initial=ieee_value(0._r8,ieee_positive_inf)
  end select
  call InitSoilProperty()
  if(len_trim(mode)>0)error stop 'invalid DMeff was accepted'
  call require('default parameter is zero',all(DMeff_col==0._r8))
  do k=0,1000
    saturation=real(k,r8)/1000._r8
    call require('zero DMeff exactly retains moisture',moisture_factor(saturation,1._r8,1,1)==1._r8)
  enddo
  DMeff_col(1,1)=.5_r8
  DMeff_col(2,2)=.9_r8
  call require('unmodified second column',moisture_factor(.8_r8,1._r8,1,2)==1._r8)
  call require('independent column amplitudes',moisture_factor(.8_r8,1._r8,2,2)< &
    moisture_factor(.8_r8,1._r8,1,1))
  do k=0,1000
    saturation=real(k,r8)/1000._r8
    f=moisture_factor(saturation,1._r8,1,1)
    call require('finite bounded factor',ieee_is_finite(f).and.f>=0._r8.and.f<=1._r8)
    if(saturation<=.6_r8)call require('dry tail unchanged',f==1._r8)
    if(saturation>.6_r8.and.saturation<1._r8)call require('wet tail reduced',f<1._r8)
    call require('effective content within actual water',saturation*f<=saturation)
  enddo
  call require('exact wet cutoff unchanged',moisture_factor(.6_r8,1._r8,1,1)==1._r8)
  call require('supplied cutoff remains abrupt',moisture_factor(NEAREST(.6_r8,1._r8),1._r8,1,1)<.9_r8)
  call require('saturation unchanged',moisture_factor(1._r8,1._r8,1,1)==1._r8)
  call require('oversaturation bounded',moisture_factor(2._r8,1._r8,1,1)==1._r8)
  call require('empty or invalid pore support is neutral',moisture_factor(0._r8,0._r8,1,1)==1._r8)
  theta_eff=.8_r8*moisture_factor(.8_r8,1._r8,1,1)
  call require('effective moisture selects a drier K class',getMoistK(theta_eff,1._r8)>getMoistK(.8_r8,1._r8))
  DMeff_col(1,1)=8._r8
  call require('large test amplitudes cannot make water negative',moisture_factor(.8_r8,1._r8,1,1)==0._r8)
  call DestructSoilProperty()
  call require('column storage released',.not.allocated(DMeff_col))
  DMeff_initial=.5_r8
  call InitSoilProperty()
  call require('uniform test input initializes every column',all(DMeff_col==.5_r8))
  call DestructSoilProperty()
  write(*,*)'MoistureFactorTest passed'
contains
!------------------------------------------------------------------------------------------
  subroutine require(label,condition)
    ! Description:
    ! Fail with the name of the violated moisture-scaling invariant.
    character(len=*), intent(in) :: label ! invariant description [-]
    logical, intent(in) :: condition      ! expected truth value [-]
    character(len=*), parameter :: subname='require'
    call PrintInfo('beg '//subname)
    if(.not.condition)then
      write(*,*)trim(label)
      error stop 'MoistureFactorTest failed'
    endif
    call PrintInfo('end '//subname)
  end subroutine require
end program MoistureFactorTest
