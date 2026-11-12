program StandingDeadCanopyTest
  ! Description:
  ! Check passive heat/water conservation, donor limits, mortality shell capacity and cover continuity.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use UptakesMod, only: SolvePassiveCanopyEnergy,PassiveCanopyResult,InitUptake
  use UptakesMod, only: CanopyStemThermalCarbon,CanopyExchangeAreaFraction
  use MiniMathMod, only: vapsat
  use MiniFuncMod, only: gOC_to_m3_OM
  use EcoSimConst, only: cpw,cps,cpo
  use DebugToolMod, only: PrintInfo
  use abortutils, only: iulog
  implicit none
  type(PassiveCanopyResult) :: Result,LiveStem
  real(r8) :: T0,Ta,qa,Cdry,Tissue,W0,S0,Rain,Snowfall
  real(r8) :: EffectiveLiveC,EffectiveDeadC,Cover0,CoverSmall,Partition
  call InitUptake()
  T0=290._r8;Ta=300._r8;qa=0._r8;Cdry=cpo*gOC_to_m3_OM(100._r8)
  Tissue=0._r8;W0=0._r8;S0=0._r8;Rain=0._r8;Snowfall=0._r8
  call check_case('dry standing dead',1._r8,0.4_r8,0._r8,0._r8)
  if(Result%water_flux/=0._r8 .or. Result%snow_flux/=0._r8)error stop 'dry canopy created evaporative water'
  if(Result%temperature<=T0 .or. Result%sensible==0._r8)error stop 'dead canopy did not exchange heat'
  W0=1.e-6_r8;S0=2.e-6_r8;Rain=3.e-6_r8;Snowfall=4.e-6_r8
  call check_case('liquid and snow donor limits',1._r8,1._r8,0._r8,0._r8)
  if(Result%water<0._r8 .or. Result%snow<0._r8)error stop 'negative intercepted water'
  if(Result%water_flux<-(W0+Rain) .or. Result%snow_flux<-(S0+Snowfall))error stop 'vapor overdraw'
  qa=2._r8*vapsat(Ta);W0=0._r8;S0=0._r8;Rain=0._r8;Snowfall=0._r8
  call check_case('dew onto initially dry material',1._r8,0._r8,0._r8,0._r8)
  if(Result%water_flux<=0._r8)error stop 'missing passive condensation'
  T0=260._r8;Ta=260._r8;qa=1.1_r8*vapsat(Ta)
  call check_case('supercooled liquid condensation',.01_r8,0._r8,0._r8,0._r8)
  if(Result%temperature>=273.15_r8.or.Result%water_flux<=0._r8.or.Result%snow_flux/=0._r8) &
    error stop 'inconsistent supercooled liquid approximation'
  T0=290._r8;Ta=300._r8
  W0=1.e-5_r8;S0=2.e-5_r8;Tissue=3.e-5_r8
  call check_case('buried material retains heat capacity',0._r8,0._r8,0._r8,0._r8)
  if(Result%temperature/=T0 .or. Result%water/=W0 .or. Result%snow/=S0)error stop 'burial changed isolated material'
  if(Result%sensible/=0._r8 .or. Result%latent/=0._r8)error stop 'buried canopy exchanged with air'

  ! A 500 gC sapwood pool and its 5000 gC whole dead stem both represent
  ! the same 0.3 m2 outer shell: mortality must not activate the interior.
  EffectiveLiveC=CanopyStemThermalCarbon(500._r8,.3_r8,2.e-6_r8)
  EffectiveDeadC=CanopyStemThermalCarbon(5000._r8,.3_r8,2.e-6_r8)
  if(ABS(EffectiveDeadC-375._r8)>1.e-10_r8)error stop 'dead capacity included internal wood'
  if(EffectiveLiveC/=EffectiveDeadC)error stop 'live/dead shell mass changed at mortality'
  Tissue=0._r8;W0=0._r8;S0=0._r8;qa=0._r8
  Cdry=cpo*gOC_to_m3_OM(EffectiveLiveC)
  call check_case('live stem shell',.15_r8,.02_r8,0._r8,0._r8)
  LiveStem=Result
  Cdry=cpo*gOC_to_m3_OM(EffectiveDeadC)
  call check_case('same stem after death',.15_r8,.02_r8,0._r8,0._r8)
  if(Result%temperature/=LiveStem%temperature .or. Result%sensible/=LiveStem%sensible) &
    error stop 'mortality changed stem thermal response without heat input'
  if(CanopyStemThermalCarbon(100._r8,.3_r8,2.e-6_r8)/=100._r8)error stop 'thin material lost thermal mass'

  ! Adding arbitrarily little dead area must not suddenly include live stems
  ! in the cover multiplier, for this PFT or for another PFT in the column.
  Cover0=CanopyExchangeAreaFraction(2._r8,2._r8,1._r8,0._r8,1._r8)
  CoverSmall=CanopyExchangeAreaFraction(2._r8,2._r8+1.e-10_r8,1._r8,1.e-10_r8,1._r8)
  if(ABS(Cover0-.5_r8)>1.e-14_r8 .or. ABS(CoverSmall-Cover0)>1.e-9_r8)error stop 'discontinuous canopy cover'
  Partition=CanopyExchangeAreaFraction(2._r8,4._r8,3._r8,1._r8,1._r8) &
    +CanopyExchangeAreaFraction(2._r8,4._r8,3._r8,1._r8,1._r8)
  if(Partition/=1._r8)error stop 'PFT effective areas exceed ground area'
  write(iulog,*)'PASS: passive canopy heat/water, dry/dew/snow donors and burial'
contains
!------------------------------------------------------------------------------------------
  subroutine check_case(label,Area,SW,LW,F)
    ! Description:
    ! Check the physical energy and water budgets of one accepted one-hour solve.
    character(len=*), intent(in) :: label ! scenario name [-]
    real(r8), intent(in) :: Area ! effective exchange area [m2]
    real(r8), intent(in) :: SW,LW,F ! incident energies [MJ] and one-side emission factor [MJ/K4]
    character(len=*), parameter :: subname='check_case'
    real(r8) :: Eold,Heat,Error
    call PrintInfo('beg '//subname)
    associate(T => Result%temperature,C => Result%capacity,Rn => Result%net_radiation,H => Result%sensible, &
      LE => Result%latent,Adv => Result%vapor_enthalpy,W => Result%water,S => Result%snow, &
      dW => Result%water_flux,dS => Result%snow_flux)
      call SolvePassiveCanopyEnergy(T0,Ta,qa,Cdry,Tissue,W0,S0,Rain,Snowfall,Area,SW,LW,F,.02_r8,.001_r8,Result)
      Eold=(Cdry+cpw*(Tissue+W0)+cps*S0)*T0
      Heat=Rn+H+LE+Adv+(cpw*Rain+cps*Snowfall)*Ta
      Error=C*T-Eold-Heat
      if(.not.ieee_is_finite(T) .or. ABS(Error)>2.e-11_r8*MAX(1._r8,ABS(Eold)))then
        write(iulog,*)label,' energy residual',Error
        error stop 'passive canopy energy failed'
      endif
      if(ABS(W-W0-Rain-dW)>1.e-18_r8 .or. ABS(S-S0-Snowfall-dS)>1.e-18_r8)error stop 'passive canopy water failed'
    end associate
    call PrintInfo('end '//subname)
  end subroutine check_case
end program StandingDeadCanopyTest
