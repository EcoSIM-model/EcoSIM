program ThinSnowEnergyTest
  ! Check energy retention, redistribution and conservative thin-snow removal.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use abortutils, only: iulog
  use DebugToolMod, only: PrintInfo
  use GridConsts, only: JX,JY,JX0,JY0,JZ,JH,JV,JS
  use GridDataType, only: AREA_3D,NUM_col,NU_col
  use EcoSIMCtrlMod, only: lverb,fixWaterLevel,salt_model
  use EcoSIMCtrlDataType, only: ZERO,ZERO2,ZEROS,ZEROS2
  use EcoSIMSolverPar, only: NPH
  use EcoSimConst, only: cps,cpw,cpi,cpo,DENSICE,TFice,VLHeatCapSnoMin
  use HydroThermData, only: TKSnow0_snvr
  use ElmIDMod, only: ielmc
  use MiniFuncMod, only: gOC_to_m3_OM
  use SnowDataType
  use SnowPhysData
  use SnowBalanceMod, only: SnowMassUpdate,SnowpackLayering
  use SnowPhysMod, only: UpdateSnowAtM
  use ClimForcDataType, only: TairK_col,TCA_col
  use SoilHeatDataType, only: InitSoilHeatData,TKS_vr,VHeatCapacity_vr,XPhaseChangeHeatL_snvr
  use SoilWaterDataType, only: VLWatMicP_vr,VLWatMacP_vr,VLiceMicP_vr,VLiceMacP_vr,Qinflx2Soil_col
  use SoilPhysDataType, only: VHeatCapSolidSoil_vr
  use SoilPropertyDataType, only: SoilBulkDensity_vr
  use SOMDataType, only: SoilOrgM_vr
  use SurfLitterDataType, only: WatFLo2LitR_col
  use EcoSimSumDataType, only: HEATIN_lnd
  implicit none
  real(r8), parameter :: SnowAmounts(3)=[1.e-11_r8,0.999e-4_r8,1.001e-4_r8]
  real(r8) :: BeforeEnergy,BeforeWater,BeforeIce,Incoming,ReturnedWater,ReturnedHeat
  real(r8) :: LowerEnergy
  integer :: it,M,Receiver

  JX=1; JY=1; JX0=1; JY0=1; JZ=JS; JH=1; JV=1
  lverb=.false.; fixWaterLevel=.false.; salt_model=.false.; NPH=2
  ZERO=1.e-16_r8; ZERO2=1.e-8_r8
  call InitSnowData()
  call InitSnowPhysData()
  call InitSoilHeatData()
  allocate(TKSnow0_snvr(JS,1,1))
  allocate(ZEROS(1,1),ZEROS2(1,1),AREA_3D(3,0:JS,1,1),NUM_col(1,1),NU_col(1,1))
  allocate(TairK_col(1,1),TCA_col(1,1),SoilBulkDensity_vr(0:JS,1,1),SoilOrgM_vr(3,0:JS,1,1))
  allocate(VLWatMicP_vr(0:JS,1,1),VLWatMacP_vr(0:JS,1,1),VLiceMicP_vr(0:JS,1,1),VLiceMacP_vr(0:JS,1,1))
  allocate(VHeatCapSolidSoil_vr(0:JS,1,1),Qinflx2Soil_col(1,1),WatFLo2LitR_col(1,1))
  ZEROS=ZERO; ZEROS2=ZERO2; AREA_3D=1._r8; NUM_col=1; NU_col=1

  ! Positive stores on both sides of the former threshold retain energy.
  ! The smallest case also lies below the former geometric deletion cutoff.
  do it=1,size(SnowAmounts)
    call reset_state()
    call set_layer(1,SnowAmounts(it),260._r8)
    BeforeEnergy=sum(VLHeatCapSnow_snvr*TKSnow_snvr)
    call SnowMassUpdate(1,1,1,1)
    call check_close('hourly energy retention',sum(VLHeatCapSnow_snvr*TKSnow_snvr),BeforeEnergy)
    call check_close('hourly mass retention',VLDrySnoWE_snvr(1,1,1),SnowAmounts(it))
    call check_close('no prescribed external heat',HEATIN_lnd,12._r8)
  enddo

  ! Redistributed cold snow must carry its own energy into a thin receiver,
  ! including the last-substep branch, which updates a different capacity cache.
  do M=1,NPH
    call reset_state()
    call set_layer(1,5.e-5_r8,260._r8)
    VLDrySnoWE0_snvr=VLDrySnoWE_snvr; TKSnow0_snvr=TKSnow_snvr
    VLSnowHeatCapM_snvr(M,:,:,:)=VLHeatCapSnow_snvr
    Incoming=1.e-5_r8
    cumDrySnoByRedistM_col=Incoming
    cumHeatSnoByRedistM_col=cps*Incoming*250._r8
    BeforeEnergy=sum(VLHeatCapSnow_snvr*TKSnow_snvr)+cps*Incoming*250._r8
    call UpdateSnowAtM(1,1,M,1,1)
    call check_close('redistributed energy',sum(VLHeatCapSnow_snvr*TKSnow_snvr),BeforeEnergy)
    call check_close('redistributed mass',sum(VLDrySnoWE_snvr),5.e-5_r8+Incoming)
  enddo

  ! Remove thin frozen snow into soil, then litter. Conserve both sensible
  ! energy and latent inventory, while leaving an occupied lower layer intact.
  do Receiver=0,1
    call reset_state()
    call set_layer(1,5.e-5_r8,260._r8)
    call set_layer(2,2.e-4_r8,265._r8)
    TairK_col=280._r8
    if(Receiver==0)then
      SoilOrgM_vr(ielmc,0,1,1)=100._r8
      VHeatCapacity_vr(0,1,1)=cpo*gOC_to_m3_OM(100._r8)+cpw*VLWatMicP_vr(0,1,1)+cpi*VLiceMicP_vr(0,1,1)
    endif
    BeforeEnergy=sum(VLHeatCapSnow_snvr*TKSnow_snvr)+VHeatCapacity_vr(Receiver,1,1)*TKS_vr(Receiver,1,1)
    BeforeWater=sum(VLDrySnoWE_snvr)+sum(VLWatSnow_snvr)+DENSICE*sum(VLIceSnow_snvr) &
      +VLWatMicP_vr(Receiver,1,1)+DENSICE*VLiceMicP_vr(Receiver,1,1)
    BeforeIce=sum(VLDrySnoWE_snvr)+DENSICE*sum(VLIceSnow_snvr)+DENSICE*VLiceMicP_vr(Receiver,1,1)
    LowerEnergy=VLHeatCapSnow_snvr(2,1,1)*TKSnow_snvr(2,1,1)
    call SnowMassUpdate(1,1,1,1)
    call check_close('removal energy',sum(VLHeatCapSnow_snvr*TKSnow_snvr) &
      +VHeatCapacity_vr(Receiver,1,1)*TKS_vr(Receiver,1,1),BeforeEnergy)
    call check_close('removal water',sum(VLDrySnoWE_snvr)+sum(VLWatSnow_snvr)+DENSICE*sum(VLIceSnow_snvr) &
      +VLWatMicP_vr(Receiver,1,1)+DENSICE*VLiceMicP_vr(Receiver,1,1),BeforeWater)
    call check_close('removal frozen inventory',sum(VLDrySnoWE_snvr)+DENSICE*sum(VLIceSnow_snvr) &
      +DENSICE*VLiceMicP_vr(Receiver,1,1),BeforeIce)
    call check_close('lower layer energy retained',VLHeatCapSnow_snvr(2,1,1)*TKSnow_snvr(2,1,1),LowerEnergy)
    call check_close('no reset heat before removal',HEATIN_lnd,12._r8)
  enddo

  ! Layer merging conserves the temperatures' energy weights and water mass.
  call reset_state()
  call set_layer(1,1.5e-4_r8,260._r8)
  call set_layer(2,0.5e-4_r8,250._r8)
  BeforeEnergy=sum(VLHeatCapSnow_snvr*TKSnow_snvr)
  BeforeWater=sum(VLDrySnoWE_snvr)
  call SnowpackLayering(1,1,1,1)
  call check_close('merged energy',sum(VLHeatCapSnow_snvr*TKSnow_snvr),BeforeEnergy)
  call check_close('merged mass',sum(VLDrySnoWE_snvr),BeforeWater)
  call check_close('single occupied merged layer',real(nsnol_col(1,1),r8),1._r8)

  ! Optional coupling outputs are defined even when no transfer occurs.
  call reset_state()
  ReturnedWater=-1._r8; ReturnedHeat=-1._r8
  call SnowMassUpdate(1,1,1,1,QWatinfl2Mic=ReturnedWater)
  call SnowMassUpdate(1,1,1,1,QHeatInfl2Soil=ReturnedHeat)
  call check_close('no water transfer',ReturnedWater,0._r8)
  call check_close('no heat transfer',ReturnedHeat,0._r8)
  write(*,*)'PASS: thin snow energy, redistribution, removal, merging and optional outputs'
contains
!------------------------------------------------------------------------------------------
  subroutine reset_state()
  ! Set a cold one-column surface and zero all transfers used by these checks.
  character(len=*), parameter :: subname='reset_state'
  call PrintInfo('beg '//subname)
  VLDrySnoWE_snvr=0._r8; VLWatSnow_snvr=0._r8; VLIceSnow_snvr=0._r8
  VLHeatCapSnow_snvr=0._r8; SnoDens_snvr=0.1_r8; NewSnowDens_col=0.1_r8
  TKSnow_snvr=260._r8; TCSnow_snvr=TKSnow_snvr-TFice
  VLSnoDWIprev_snvr=0._r8; SnowThickL_snvr=0._r8; cumSnowDepz_col=0._r8
  VLSnoDWIMax_snvr=0.1_r8; VLHeatCapSnowMin_col=VLHeatCapSnoMin; nsnol_col=0
  SnoXfer2SnoLay_snvr=0._r8; WatXfer2SnoLay_snvr=0._r8; IceXfer2SnoLay_snvr=0._r8; HeatXfer2SnoLay_snvr=0._r8
  CumSno2SnowL_snvr=0._r8; CumWat2SnowL_snvr=0._r8; CumIce2SnowL_snvr=0._r8; CumHeat2SnowL_snvr=0._r8
  TDrysnoByRedist_col=0._r8; TIceSnoByRedist_col=0._r8; TWatSnoByRedist_col=0._r8; THeatSnoByRedist_col=0._r8
  cumDrySnoByRedistM_col=0._r8; cumWatSnoByRedistM_col=0._r8
  cumIceSnoByRedistM_col=0._r8; cumHeatSnoByRedistM_col=0._r8
  VLDrySnoWE0_snvr=0._r8; VLIceSnow0_snvr=0._r8; VLWatSnow0_snvr=0._r8
  THeatSnowThaw_col=0._r8; XPhaseChangeHeatL_snvr=0._r8; HEATIN_lnd=12._r8
  QSnowH2Oloss_col=0._r8; QSnowHeatLoss_col=0._r8; Qinflx2Soil_col=0._r8; WatFLo2LitR_col=0._r8
  TairK_col=240._r8; TCA_col=TairK_col-TFice
  SoilBulkDensity_vr=1._r8; SoilOrgM_vr=0._r8; VHeatCapSolidSoil_vr=1._r8
  VLWatMicP_vr=1.e-3_r8; VLWatMacP_vr=2.e-4_r8; VLiceMicP_vr=1.e-4_r8; VLiceMacP_vr=5.e-5_r8
  VHeatCapacity_vr=VHeatCapSolidSoil_vr+cpw*(VLWatMicP_vr+VLWatMacP_vr)+cpi*(VLiceMicP_vr+VLiceMacP_vr)
  VHeatCapacity_vr(0,:,:)=cpw*VLWatMicP_vr(0,:,:)+cpi*VLiceMicP_vr(0,:,:)
  TKS_vr=270._r8
  trcg_solsml_snvr=0._r8; trcn_solsml_snvr=0._r8
  call PrintInfo('end '//subname)
  end subroutine reset_state
!------------------------------------------------------------------------------------------
  subroutine set_layer(L,WaterWE,Temperature)
  ! Initialize a dry-snow layer and its energy/geometry caches.
  integer, intent(in) :: L             !snow layer index [1:JS]
  real(r8), intent(in) :: WaterWE      !dry-snow water-equivalent volume [m3/cell]
  real(r8), intent(in) :: Temperature  !snow temperature [K]
  character(len=*), parameter :: subname='set_layer'
  call PrintInfo('beg '//subname)
  VLDrySnoWE_snvr(L,1,1)=WaterWE
  VLHeatCapSnow_snvr(L,1,1)=cps*WaterWE
  TKSnow_snvr(L,1,1)=Temperature; TCSnow_snvr(L,1,1)=Temperature-TFice
  VLSnoDWIprev_snvr(L,1,1)=WaterWE/SnoDens_snvr(L,1,1)
  SnowThickL_snvr(L,1,1)=VLSnoDWIprev_snvr(L,1,1)
  nsnol_col(1,1)=nsnol_col(1,1)+1
  call PrintInfo('end '//subname)
  end subroutine set_layer
!------------------------------------------------------------------------------------------
  subroutine check_close(Label,Actual,Expected)
  ! Require finite values and relative closure of physical inventories.
  character(len=*), intent(in) :: Label !invariant being checked
  real(r8), intent(in) :: Actual       !computed inventory [units specified by Label]
  real(r8), intent(in) :: Expected     !reference inventory [same units as Actual]
  character(len=*), parameter :: subname='check_close'
  call PrintInfo('beg '//subname)
  if(.not.ieee_is_finite(Actual) .or. abs(Actual-Expected)>2.e-13_r8*max(abs(Expected),1.e-20_r8))then
    write(iulog,*)'FAIL: ',Label,Actual,Expected
    error stop 1
  endif
  call PrintInfo('end '//subname)
  end subroutine check_close
end program ThinSnowEnergyTest
