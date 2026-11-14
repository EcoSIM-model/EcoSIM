program TillageWeightsTest
  ! Description:
  ! Verify conservative tillage transfers with independently rounded layer
  ! depths/thicknesses, partial tillage and an empty tillage interval.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use abortutils, only: iulog
  use DebugToolMod, only: PrintInfo
  use GridConsts, only: JZ
  use GridDataType, only: NU_col,NL_col,DLYR_3D,CumSoilThickness_vr
  use EcoSIMCtrlDataType, only: ZERO
  use TillageMixMod, only: DeriveTillageProf,MixExtensive1D,MixExtensive2D,MixExtensive3D,MixExtensive4D
  implicit none
  real(r8) :: FI(4),TI(4),Energy(0:4),Water(0:4),Elements(3,0:4)
  real(r8) :: Pools3(2,3,0:4),Pools4(2,2,3,0:4),BeforeElements(3),Before3(2,3),Before4(2,2,3)
  real(r8) :: BeforeEnergy,BeforeWater,Depth
  integer :: LL,L,NE,it

  JZ=4;ZERO=1.e-16_r8
  allocate(NU_col(1,1),NL_col(1,1),DLYR_3D(3,0:4,1,1),CumSoilThickness_vr(0:4,1,1))
  NU_col=1;NL_col=4;DLYR_3D=0._r8;CumSoilThickness_vr=0._r8
  ! Representative single-precision input rounding leaves a small mismatch
  ! between independently supplied cumulative depths and layer thicknesses.
  DLYR_3D(3,1:4,1,1)=[0.009999999776482582_r8,0.04000000096857548_r8, &
    0.05000000074505806_r8,0.04999999701976776_r8]
  CumSoilThickness_vr(1:4,1,1)=[0.009999999776482582_r8,0.04999999701976776_r8, &
    0.09999999403953552_r8,0.1499999910593033_r8]

  do it=1,2
    Depth=0.1_r8
    if(it==2)Depth=0.075_r8
    call DeriveTillageProf(1,1,Depth,LL,FI,TI)
    call check_close('recipient shares sum to one',sum(FI),1._r8)
    if(any(FI<0._r8).or.any(TI<0._r8).or.any(TI>1._r8))error stop 'invalid tillage fractions'
    if(it==2)then
      call check_close('partly tilled layer',TI(3), &
        (Depth-(CumSoilThickness_vr(3,1,1)-DLYR_3D(3,3,1,1)))/DLYR_3D(3,3,1,1))
      call check_close('untouched lower layer',TI(4),0._r8)
    endif
    Energy=[0.38956796391073312_r8,20._r8,30._r8,40._r8,50._r8]
    Water=[1.e-4_r8,0.01_r8,0.02_r8,0.03_r8,0.04_r8]
    do L=0,4
      do NE=1,3
        Elements(NE,L)=real(NE*(L+1),r8)*0.1_r8
      enddo
      Pools3(:,:,L)=real(L+1,r8)*0.03_r8
      Pools4(:,:,:,L)=real(L+1,r8)*0.07_r8
    enddo
    BeforeEnergy=sum(Energy);BeforeWater=sum(Water)
    BeforeElements=sum(Elements,dim=2);Before3=sum(Pools3,dim=3);Before4=sum(Pools4,dim=4)
    call MixExtensive1D(0.7_r8,TI,FI,1,LL,Energy,0.7_r8)
    call MixExtensive1D(0.7_r8,TI,FI,1,LL,Water,0.7_r8)
    call MixExtensive2D(0.7_r8,TI,FI,1,LL,Elements,0.7_r8)
    call MixExtensive3D(0.7_r8,TI,FI,1,LL,Pools3,0.7_r8)
    call MixExtensive4D(0.7_r8,TI,FI,1,LL,Pools4,0.7_r8)
    call check_close('heat retained by internal tillage',sum(Energy),BeforeEnergy)
    call check_close('water retained by internal tillage',sum(Water),BeforeWater)
    do NE=1,3
      call check_close('C/N/P retained by internal tillage',sum(Elements(NE,:)),BeforeElements(NE))
    enddo
    call check_close('three-dimensional pools retained',maxval(abs(sum(Pools3,dim=3)-Before3)),0._r8)
    call check_close('four-dimensional pools retained',maxval(abs(sum(Pools4,dim=4)-Before4)),0._r8)
  enddo

  call DeriveTillageProf(1,1,0._r8,LL,FI,TI)
  if(LL>=1.or.any(FI/=0._r8).or.any(TI/=0._r8))error stop 'nonempty zero-depth tillage profile'
  write(iulog,*)'PASS: tillage retains heat, water and C/N/P with rounded geometry'
contains
!------------------------------------------------------------------------------------------
  subroutine check_close(Label,Actual,Expected)
    ! Description:
    ! Compare conservative pool totals and normalized geometry weights.
    character(len=*), intent(in) :: Label ! invariant being checked [-]
    real(r8), intent(in) :: Actual ! computed total [units identified by Label]
    real(r8), intent(in) :: Expected ! reference total [same units as Actual]
    character(len=*), parameter :: subname='check_close'
    call PrintInfo('beg '//subname)
    if(.not.ieee_is_finite(Actual).or.abs(Actual-Expected)>2.e-13_r8*max(abs(Expected),1._r8))then
      write(iulog,*)Label,Actual,Expected
      error stop 'nonconservative tillage redistribution'
    endif
    call PrintInfo('end '//subname)
  end subroutine check_close
end program TillageWeightsTest
