program SurfaceVaporTest
  ! Description:
  ! Verify physical-volume source partition, phase donor limits, frost/dew,
  ! common ice-equilibrium vapor and fixed-store latent-heat derivatives.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use EcoSimConst, only: DENSICE,EvapLHTC,SublmHTC,TFice,LtHeatIceMelt,RGASC
  use SurfaceVaporExchangeMod, only: SurfaceVaporDensity,PartitionSurfaceVapor
  use MiniMathMod, only: vapsat
  use DebugToolMod, only: PrintInfo
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real(r8) :: Liquid,Ice,Derivative,LiquidPlus,IcePlus,Unused,Step,Expected
  integer :: k
  real(r8),parameter :: Demands(3)=[-0.1_r8,-3.8_r8,-10._r8]

  call PartitionSurfaceVapor(-0.1_r8,1._r8,3._r8,260._r8,1._r8,Liquid,Ice,Derivative)
  call close('liquid physical-volume weight',Liquid,-0.025_r8)
  call close('ice physical-volume weight, not mass weight',Ice,-0.075_r8)
  call close('unsaturated donor demand',Liquid+Ice,-0.1_r8)
  call close('mixed latent derivative',Derivative,0.25_r8*EvapLHTC+0.75_r8*SublmHTC)
  call PartitionSurfaceVapor(-10._r8,1._r8,3._r8,260._r8,1._r8,Liquid,Ice,Derivative)
  call close('liquid exhausted',Liquid,-1._r8)
  call close('ice exhausted on water-equivalent basis',Ice,-3._r8*DENSICE)
  call close('exhausted latent derivative',Derivative,0._r8)
  call PartitionSurfaceVapor(-10._r8,0._r8,3._r8,260._r8,0.5_r8,Liquid,Ice,Derivative)
  call close('ice-only has no liquid debit',Liquid,0._r8)
  call close('fractional ice donor safety',Ice,-1.5_r8*DENSICE)
  call PartitionSurfaceVapor(-0.1_r8,1._r8,0._r8,280._r8,1._r8,Liquid,Ice,Derivative)
  call close('liquid-only loss',Liquid,-0.1_r8)
  call close('liquid-only no ice',Ice,0._r8)
  call PartitionSurfaceVapor(-0.1_r8,0._r8,0._r8,260._r8,1._r8,Liquid,Ice,Derivative)
  call close('empty donor',abs(Liquid)+abs(Ice)+Derivative,0._r8)
  call PartitionSurfaceVapor(0.1_r8,0._r8,0._r8,260._r8,1._r8,Liquid,Ice,Derivative,0.1_r8)
  call close('frost on initially empty surface',Ice,0.1_r8)
  call close('frost latent derivative',Derivative,SublmHTC)
  call close('frost not liquid',Liquid,0._r8)
  call PartitionSurfaceVapor(0.1_r8,0._r8,0._r8,260._r8,1._r8,Liquid,Ice,Derivative,-0.05_r8)
  call close('cold adsorption below ice saturation stays liquid',Liquid,0.1_r8)
  call close('no frost below ice saturation',Ice,0._r8)
  call close('adsorption latent derivative',Derivative,EvapLHTC)
  call PartitionSurfaceVapor(0.1_r8,0._r8,0._r8,260._r8,1._r8,Liquid,Ice,Derivative,0.04_r8)
  call close('frost limited by ice supersaturation',Ice,0.04_r8)
  call close('remaining cold adsorption',Liquid,0.06_r8)
  call close('combined frost and adsorption derivative',Derivative,SublmHTC)
  call PartitionSurfaceVapor(0.1_r8,0._r8,0._r8,280._r8,1._r8,Liquid,Ice,Derivative)
  call close('dew on warm surface',Liquid,0.1_r8)
  call close('dew not ice',Ice,0._r8)

  ! Compare the heat derivative to finite differences, including a branch with
  ! ice capped but liquid still responding. Keep the stencil off limiter joins.
  Step=1.e-6_r8
  do k=1,size(Demands)
    call PartitionSurfaceVapor(Demands(k),1._r8,3._r8,260._r8,1._r8,Liquid,Ice,Derivative)
    call PartitionSurfaceVapor(Demands(k)+Step,1._r8,3._r8,260._r8,1._r8,LiquidPlus,IcePlus,Unused)
    Expected=((LiquidPlus-Liquid)*EvapLHTC+(IcePlus-Ice)*SublmHTC)/Step
    if(abs(Expected-Derivative)>1.e-5_r8)error stop 'latent derivative disagrees with accepted phase flux'
  enddo
  call close('liquid equilibrium retains Kelvin correction',SurfaceVaporDensity(280._r8,-1._r8,0._r8), &
    vapsat(280._r8)*EXP(-18._r8/(RGASC*280._r8)))
  Expected=vapsat(260._r8)*EXP(18._r8*LtHeatIceMelt*(260._r8-TFice)/(TFice*RGASC*260._r8))
  call close('ice sets common chemical potential',SurfaceVaporDensity(260._r8,-100._r8,1._r8),Expected)
  call close('ice amount does not scale equilibrium RH', &
    SurfaceVaporDensity(260._r8,-100._r8,1._r8),SurfaceVaporDensity(260._r8,-1._r8,0.01_r8))
  call close('continuous at melting temperature',SurfaceVaporDensity(TFice,0._r8,1._r8),vapsat(TFice))
  write(*,*)'PASS: surface phase-volume partition, vapor equilibrium and latent accounting'
contains
!------------------------------------------------------------------------------------------
  subroutine close(Label,Actual,Reference)
    ! Description:
    ! Check independent scalar mass and heat references with roundoff tolerance.
    character(len=*),intent(in) :: Label ! invariant being tested [-]
    real(r8),intent(in) :: Actual ! computed value [test-dependent native units]
    real(r8),intent(in) :: Reference ! independently calculated reference [same units]
    character(len=*),parameter :: subname='close'
    call PrintInfo('beg '//subname)
    if(.not.ieee_is_finite(Actual).or.abs(Actual-Reference)>1.e-12_r8*MAX(1.e-20_r8,abs(Reference)))then
      write(*,*)Label,Actual,Reference
      error stop 'surface vapor check failed'
    endif
    call PrintInfo('end '//subname)
  end subroutine close
end program SurfaceVaporTest
