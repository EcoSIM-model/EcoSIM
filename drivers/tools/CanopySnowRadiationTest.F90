program CanopySnowRadiationTest
  ! Description:
  ! Check snow-added area, burial/re-exposure, grid-size invariance, PFT radiation
  ! partition and the resistance law shared by live and passive canopy solves.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use PlantMathFuncMod, only: CanopySnowOptics,CanopyRadiationFractions
  use UptakesMod, only: InitUptake,CanopyVaporResistance,SolvePassiveCanopyEnergy,PassiveCanopyResult
  use DebugToolMod, only: PrintInfo
  use abortutils, only: iulog
  implicit none
  real(r8) :: Fraction,Area,ScaledFraction,ScaledArea,Ground,ScaledGround
  real(r8) :: Live(2),Dead(2),Snow(2),Canopy(2),LiveFraction(2),ScaledCanopy(2),ScaledLive(2)
  real(r8) :: Resistance,Expected
  integer :: SnowType(2)
  type(PassiveCanopyResult) :: Result

  call InitUptake()
  call CanopySnowOptics(4,654.8029119215743_r8,.6068384808278249_r8,Fraction,Area)
  call close_value('snow adds area rather than multiplying by cell-total SWE',Area,715.4867600043568_r8)
  call CanopySnowOptics(4,65480.29119215743_r8,60.68384808278249_r8,ScaledFraction,ScaledArea)
  call close_value('snow cover is invariant to cell subdivision',ScaledFraction,Fraction)
  call close_value('radiative area is extensive',ScaledArea,100._r8*Area)
  call CanopySnowOptics(4,0._r8,.6068384808278249_r8,Fraction,Area)
  if(Fraction/=0._r8.or.Area/=0._r8)error stop 'buried canopy intercepted radiation'
  call CanopySnowOptics(4,654.8029119215743_r8,.6068384808278249_r8,Fraction,Area)
  if(Fraction<=0._r8.or.Area<=654.8029119215743_r8)error stop 're-exposure used stale snow optics'
  call CanopySnowOptics(4,1.e-310_r8,0._r8,Fraction,Area)
  if(Fraction/=0._r8)error stop 'tiny snow-free canopy acquired snow cover'

  Live=[4._r8,8._r8];Dead=[2._r8,1._r8];Snow=[.01_r8,.02_r8];SnowType=[4,1]
  call CanopyRadiationFractions(Live,Dead,Snow,SnowType,100._r8,Canopy,LiveFraction,Ground)
  call close_value('column partition closes',SUM(Canopy)+Ground,1._r8)
  call close_value('different snow classes weight actual effective areas',Canopy(1)/Canopy(2),7._r8/9.4_r8)
  call close_value('live/dead partition within PFT 1',LiveFraction(1)/Canopy(1),4._r8/6._r8)
  call close_value('live/dead partition within PFT 2',LiveFraction(2)/Canopy(2),8._r8/9._r8)
  call CanopyRadiationFractions(100._r8*Live,100._r8*Dead,100._r8*Snow,SnowType,10000._r8, &
    ScaledCanopy,ScaledLive,ScaledGround)
  call close_value('column optical depth is invariant to grid size',ScaledGround,Ground)
  call close_value('PFT optical weight is invariant to grid size',ScaledCanopy(1),Canopy(1))
  call close_value('live optical weight is invariant to grid size',ScaledLive(2),LiveFraction(2))

  Snow=0._r8
  call CanopyRadiationFractions(Live,Dead,Snow,SnowType,100._r8,Canopy,LiveFraction,Ground)
  Expected=EXP(-.65_r8*15._r8/100._r8)
  call close_value('snow-free Beer law is preserved',Ground,Expected)
  call close_value('snow-free physical area partition',Canopy(1)/(1._r8-Ground),6._r8/15._r8)
  Live=0._r8;Dead=0._r8;Snow=.01_r8
  call CanopyRadiationFractions(Live,Dead,Snow,SnowType,100._r8,Canopy,LiveFraction,Ground)
  if(Ground/=1._r8.or.any(Canopy/=0._r8).or.any(LiveFraction/=0._r8)) &
    error stop 'buried snow pools changed the radiation partition'

  Resistance=CanopyVaporResistance(.05_r8,.12_r8,268.57_r8,300._r8)
  call close_value('unstable resistance reaches Richardson target',Resistance,.0125_r8)
  Resistance=CanopyVaporResistance(.05_r8,.12_r8,268.57_r8,250._r8)
  call close_value('stable resistance reaches Richardson target',Resistance,.2_r8)
  Resistance=CanopyVaporResistance(.001_r8,.12_r8,268.57_r8,300._r8)
  call close_value('physical resistance floor is preserved',Resistance,.00139_r8)
  call SolvePassiveCanopyEnergy(280._r8,275._r8,0._r8,1._r8,0._r8,0._r8,0._r8, &
    0._r8,0._r8,1._r8,1._r8,0._r8,0._r8,.05_r8,.12_r8,Result)
  call close_value('passive solve uses the common resistance law',Result%resistance, &
    CanopyVaporResistance(.05_r8,.12_r8,275._r8,Result%temperature))
  write(iulog,*)'PASS: canopy snow area, radiation partition, exposure and common resistance'
contains
!------------------------------------------------------------------------------------------
  subroutine close_value(Label,Actual,Expected)
    ! Description:
    ! Check a physical scaling or partition identity within roundoff.
    character(len=*), intent(in) :: Label ! physical invariant being checked [-]
    real(r8), intent(in) :: Actual,Expected ! calculated and reference values [same units]
    character(len=*), parameter :: subname='close_value'
    call PrintInfo('beg '//subname)
    if(ABS(Actual-Expected)>1.e-12_r8*MAX(1._r8,ABS(Expected)))then
      write(iulog,*)Label,Actual,Expected
      error stop 'canopy snow radiation invariant failed'
    endif
    call PrintInfo('end '//subname)
  end subroutine close_value
end program CanopySnowRadiationTest
