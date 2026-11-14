program ThinLitterEnergyTest
  ! Description:
  ! Check accepted litter enthalpy before/after runoff and after litterfall,
  ! including small positive capacities, small cover and snow-hidden litter.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use abortutils, only: iulog
  use DebugToolMod, only: PrintInfo
  use GridConsts, only: JX,JY,JX0,JY0,JZ,JH,JV,JD
  use GridDataType, only: NUM_col
  use EcoSIMCtrlDataType, only: ZERO,ZEROS,ZEROS2
  use EcoSIMSolverPar, only: NPH,dts_HeatWatTP
  use EcoSimConst, only: cpo,cpw
  use MiniFuncMod, only: gOC_to_m3_OM
  use HydroThermData
  use SurfLitterDataType
  use SurfPhysData
  use SoilWaterDataType
  use SoilHeatDataType
  use SOMDataType, only: SoilOrgM_vr
  use SurfSoilDataType, only: HeatByRad2Surf_col,HeatSensAir2Surf_col,HeatNet2Surf_col, &
    VLWatHeldCapSurf_col,LWRadBySurf_col
  use ClimForcDataType, only: TairK_col
  use SnowPhysData, only: cumWatFlx2LitRByRunoff_col,cumHeatFlx2LitRByRunoff_col
  use SurfLitterPhysMod, only: UpdateLitRBe4RunoffM,UpdateLitRAftRunoffM,UpdateLitRPhys,SurfLitREnergyBalanceM
  use SurfLitterPhysMod, only: EquilibrateThinLitterM
  implicit none
  real(r8), parameter :: Covers(3)=[0._r8,5.e-4_r8,0.5_r8]
  real(r8), parameter :: Capacities(2)=[1.e-4_r8,1.e-8_r8]
  real(r8), parameter :: RainWater=2.e-8_r8,RunoffWater=1.e-8_r8
  real(r8) :: ExpectedHeat,LandscapeHeat
  real(r8) :: SurfaceFluxes(5) ! litter-soil and atmospheric transfers [MJ/cell or m3/cell]
  integer :: ic,iv

  JX=1;JY=1;JX0=1;JY0=1;JZ=1;JH=1;JV=1;JD=2
  NPH=2;dts_HeatWatTP=0.5_r8;ZERO=1.e-16_r8
  call InitHydroThermData()
  call InitSoilWater()
  call InitSoilHeatData()
  call InitSurfLitter(3)
  call InitSurfPhysData()
  allocate(NUM_col(1,1),ZEROS(1,1),ZEROS2(1,1))
  allocate(SoilOrgM_vr(3,0:1,1,1),TairK_col(1,1))
  allocate(VLWatHeldCapSurf_col(1,1),LWRadBySurf_col(1,1))
  allocate(HeatByRad2Surf_col(1,1),HeatSensAir2Surf_col(1,1),HeatNet2Surf_col(1,1))
  allocate(cumWatFlx2LitRByRunoff_col(1,1),cumHeatFlx2LitRByRunoff_col(1,1))
  NUM_col=1;ZEROS=1.e-16_r8;ZEROS2=1.e-16_r8

  do ic=1,size(Covers)
    do iv=1,size(Capacities)
      call reset_state(Covers(ic),Capacities(iv))
      ExpectedHeat=VHeatCapacity1_vr(0,1,1)*TKSoil1_vr(0,1,1)+cpw*RainWater*310._r8
      WatFLo2LitrM_col=RainWater
      HeatFLoByWat2LitRM_col=cpw*RainWater*310._r8
      call UpdateLitRBe4RunoffM(1,1,1,1,1)
      call check_close('accepted heat before runoff',VHeatCapacity1_vr(0,1,1)*TKSoil1_vr(0,1,1),ExpectedHeat)
      call check_close('accepted water before runoff',VLWatMicP1_vr(0,1,1),1.e-7_r8+RainWater)

      ExpectedHeat=ExpectedHeat+cpw*RunoffWater*300._r8
      cumWatFlx2LitRByRunoff_col=RunoffWater
      cumHeatFlx2LitRByRunoff_col=cpw*RunoffWater*300._r8
      call UpdateLitRAftRunoffM(1,1,1,1,1)
      call check_close('accepted runoff heat',VHeatCapacity_vr(0,1,1)*TKS_vr(0,1,1),ExpectedHeat)
      call check_close('accepted runoff water',VLWatMicP_vr(0,1,1),1.e-7_r8+RainWater+RunoffWater)

      ! The caller adds intercepted litterfall water before updating its heat.
      VLWatMicP_vr(0,1,1)=VLWatMicP_vr(0,1,1)+RainWater
      ExpectedHeat=ExpectedHeat+cpw*RainWater*310._r8
      LandscapeHeat=7._r8
      call UpdateLitRPhys(1,1,1,1,RainWater,cpw*RainWater*310._r8,LandscapeHeat)
      call check_close('litterfall heat with small or hidden footprint',VHeatCapacity_vr(0,1,1)*TKS_vr(0,1,1),ExpectedHeat)
      call check_close('no artificial equilibration heat',LandscapeHeat,7._r8)
    enddo
  enddo

  ! Entering the exposed-surface solver must not reset an inactive positive store.
  do iv=1,size(Capacities)
  call reset_state(0._r8,Capacities(iv))
  VHeatCapacity1_vr(1,1,1)=1._r8
  ExpectedHeat=VHeatCapacity1_vr(0,1,1)*TKSoil1_vr(0,1,1)
  call SurfLitREnergyBalanceM(1,1,1,1,1,0._r8,1.e-7_r8,0._r8,0._r8, &
    0._r8,0._r8,0._r8,0._r8,0._r8,SurfaceFluxes(1),SurfaceFluxes(2), &
    SurfaceFluxes(3),SurfaceFluxes(4),SurfaceFluxes(5))
  call check_close('inactive surface entry retains energy',VHeatCapacity1_vr(0,1,1)*TKSoil1_vr(0,1,1),ExpectedHeat)
  call check_close('inactive surface adds no atmospheric or water exchange',sum(abs(SurfaceFluxes)),0._r8)
  ! The shared transfer path also runs when snow skips the exposed solver.
  call EquilibrateThinLitterM(1,1)
  SurfaceFluxes(1)=HeatFlow2Soili_3D(3,1,1,1)
  call check_close('thin litter and soil approach common temperature', &
    (ExpectedHeat-SurfaceFluxes(1))/VHeatCapacity1_vr(0,1,1), &
    TKSoil1_vr(1,1,1)+SurfaceFluxes(1)/VHeatCapacity1_vr(1,1,1))
  call UpdateLitRBe4RunoffM(1,1,1,1,1)
  call check_close('thin contact conserves combined sensible heat', &
    VHeatCapacity1_vr(0,1,1)*TKSoil1_vr(0,1,1)+SurfaceFluxes(1),ExpectedHeat)
  enddo

  ! Litter freeze/thaw was accepted earlier, even when its exposed cover is zero.
  call reset_state(0._r8,1.e-4_r8)
  VLWatMicP_vr=VLWatMicP1_vr;VHeatCapacity_vr=VHeatCapacity1_vr;TKS_vr=TKSoil1_vr
  TLitrIceHeatFlxFrez_col=1.e-5_r8
  LandscapeHeat=7._r8
  call UpdateLitRPhys(1,1,1,1,0._r8,0._r8,LandscapeHeat)
  call check_close('hidden litter phase heat is tallied',THeatSoiThaw_col(1,1),1.e-5_r8)
  call check_close('accepted phase heat in landscape ledger',LandscapeHeat,7._r8+1.e-5_r8)
  write(iulog,*)'PASS: thin/hidden litter conserves accepted water and heat'
contains
!------------------------------------------------------------------------------------------
  subroutine reset_state(Cover,CapacityThreshold)
    ! Description:
    ! Initialize one thin wet litter store colder than the underlying soil.
    real(r8), intent(in) :: Cover ! exposed litter footprint fraction [-]
    real(r8), intent(in) :: CapacityThreshold ! significant-capacity threshold [MJ/cell/K]
    character(len=*), parameter :: subname='reset_state'
    call PrintInfo('beg '//subname)
    SoilOrgM_vr=0._r8;SoilOrgM_vr(1,0,1,1)=1.e-3_r8
    m3OM_col=gOC_to_m3_OM(1.e-3_r8)
    TKSoil1_vr=285._r8;TKSoil1_vr(0,1,1)=270._r8
    VLWatMicP1_vr=0._r8;VLWatMicP1_vr(0,1,1)=1.e-7_r8
    VLiceMicP1_vr=0._r8;VLairMicP1_vr=0._r8
    VHeatCapacity1_vr=cpw*VLWatMicP1_vr;VHeatCapacity1_vr(0,1,1)=VHeatCapacity1_vr(0,1,1)+cpo*m3OM_col(1,1)
    VLPoreLitR_col=1.e-5_r8;VLitR_col=1.e-5_r8;VLWatHeldCapSurf_col=1.e-5_r8;VWatLitRHoldCapcity_col=1.e-5_r8
    VHeatCapLitRMin_col=CapacityThreshold;FracSurfByLitR_col=Cover
    WatFLo2LitrM_col=0._r8;HeatFLoByWat2LitRM_col=0._r8;HeatFLoByWat2LitR_col=0._r8
    LitrIceFlxThaw_col=0._r8;LitrIceHeatFlxFrez_col=0._r8;TLitrIceHeatFlxFrez_col=0._r8
    cumWatFlx2LitRByRunoff_col=0._r8;cumHeatFlx2LitRByRunoff_col=0._r8
    CdLitRHSens_col=0._r8;LWEmscefLitR_col=0._r8
    HeatFlow2Soili_3D=0._r8
    HeatByRad2Surf_col=0._r8;HeatSensAir2Surf_col=0._r8;HeatNet2Surf_col=0._r8;LWRadBySurf_col=0._r8
    TKS_vr=TKSoil1_vr;VHeatCapacity_vr=VHeatCapacity1_vr
    VLWatMicP_vr=VLWatMicP1_vr;VLiceMicP_vr=VLiceMicP1_vr
    THeatSoiThaw_col=0._r8;TairK_col=293._r8
    call PrintInfo('end '//subname)
  end subroutine reset_state
!------------------------------------------------------------------------------------------
  subroutine check_close(Label,Actual,Expected)
    ! Description:
    ! Require finite, conservative inventories against independently mixed water inputs.
    character(len=*), intent(in) :: Label ! invariant being checked [-]
    real(r8), intent(in) :: Actual ! computed energy or water inventory [MJ/cell or m3/cell]
    real(r8), intent(in) :: Expected ! reference inventory [same units as Actual]
    character(len=*), parameter :: subname='check_close'
    call PrintInfo('beg '//subname)
    if(.not.ieee_is_finite(Actual).or.abs(Actual-Expected)>2.e-12_r8*max(abs(Expected),1.e-20_r8))then
      write(iulog,*)Label,Actual,Expected
      error stop 'thin litter lost accepted heat or water'
    endif
    call PrintInfo('end '//subname)
  end subroutine check_close
end program ThinLitterEnergyTest
