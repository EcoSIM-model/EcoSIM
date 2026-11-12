program CanopyWaterStoreTest
  ! Description:
  ! Check native cold-start passive energy, burial retention and removal drainage
  ! with liquid/snow donor enthalpy, including slots beyond the active NP count
  ! and inert wet stores with full throughfall in radiation-only mode.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use CanopyHydroMod, only: CanopyInterceptPrecip,SyncCanopyThermalInventory
  use CanopyDataType
  use PlantTraitDataType
  use SoilWaterDataType
  use ClimForcDataType
  use PlantMgmtDataType, only: NP_col,NP0_col
  use FlagDataType, only: IsPlantActive_pft
  use EcoSIMCtrlMod, only: ldo_radiation_test,ldo_sp_mode
  use EcoSIMCtrlDataType, only: ZEROS
  use EcoSIMSolverPar, only: NPH,dts_HeatWatTP
  use ElmIDMod, only: ifalse,ielmc
  use EcoSimConst, only: cpw,cps,cpo
  use MiniFuncMod, only: gOC_to_m3_OM
  use DebugToolMod, only: PrintInfo
  use abortutils, only: iulog
  implicit none
  real(r8) :: InitialEnergy,Capacity,W,S,T,ExpectedHeat

  allocate(NP_col(1,1),NP0_col(1,1),ZEROS(1,1),TairK_col(1,1),WindSpeedAtm_col(1,1))
  allocate(BulkFactor4Snow_col(1,1),fSnowCanopy_col(1,1),LeafStalkAreaAll_col(1,1))
  allocate(LeafStalkAreaAct_pft(2,1,1),StandDeadSurfAreaAct_pft(2,1,1))
  allocate(CanopyLeafArea_pft(2,1,1),CanopyStemSurfArea_pft(2,1,1),StandDeadSurfArea_pft(2,1,1))
  allocate(CanopyLeafSheathC_pft(2,1,1),CanopySapwoodC_pft(2,1,1))
  allocate(StandDeadStrutElms_pft(3,2,1,1),StemSpecVolume_pft(2,1,1),PlantPopuLive_pft(2,1,1))
  allocate(IsPlantActive_pft(2,1,1),iPlantRootProfile_pft(2,1,1),iPlantSnowIntercepType_pft(2,1,1))
  allocate(WatHeldOnCanopy_pft(2,1,1),SnowOnCanopy_pft(2,1,1),CanopyBiomWater_pft(2,1,1))
  allocate(VHeatCapCanopy_pft(2,1,1),ENGYX_pft(2,1,1),TKC_pft(2,1,1),TKCanopy_pft(2,1,1))
  allocate(RainIntcptByCanopy_pft(2,1,1),SnowIntcptByCanopy_pft(2,1,1),FracPARads2Canopy_pft(2,1,1))
  allocate(fSnowCanopy_pft(2,1,1),BulkFactor4Snow_pft(2,1,1))
  allocate(Rain2Canopy_col(1,1),RainIntceptByCanopy_col(1,1),SnowIntceptByCanopy_col(1,1))
  allocate(PrecRainAndIrrig_col(1,1),SnoFalPrec_col(1,1),RainPrecThrufall_col(1,1),SnowPrecThrufall_col(1,1))
  allocate(TKRainThrufall_col(1,1),TKSnowThrufall_col(1,1))
  allocate(RainCanopyUnload_col(1,1),SnowCanopyUnload_col(1,1),HeatCanopyUnload_col(1,1))

  NPH=30;dts_HeatWatTP=1._r8/30._r8
  NP_col=1;NP0_col=2;ZEROS=1.e-12_r8
  W=1.e-4_r8;S=2.e-4_r8;T=260._r8
  TairK_col=270._r8;WindSpeedAtm_col=0._r8
  LeafStalkAreaAll_col=0._r8;LeafStalkAreaAct_pft=0._r8;StandDeadSurfAreaAct_pft=0._r8
  CanopyLeafArea_pft=0._r8;CanopyStemSurfArea_pft=0._r8
  StandDeadSurfArea_pft=0._r8;StandDeadSurfArea_pft(1,1,1)=.3_r8
  StandDeadStrutElms_pft=0._r8;StandDeadStrutElms_pft(ielmc,1,1,1)=5000._r8
  CanopyLeafSheathC_pft=0._r8;CanopySapwoodC_pft=0._r8
  StemSpecVolume_pft=2.e-6_r8;PlantPopuLive_pft=0._r8;IsPlantActive_pft=ifalse
  iPlantRootProfile_pft=2;iPlantSnowIntercepType_pft=4
  WatHeldOnCanopy_pft=0._r8;WatHeldOnCanopy_pft(1,1,1)=W
  SnowOnCanopy_pft=0._r8;SnowOnCanopy_pft(1,1,1)=S
  CanopyBiomWater_pft=0._r8;TKC_pft=T;TKCanopy_pft=T
  VHeatCapCanopy_pft=0._r8;ENGYX_pft=0._r8
  FracPARads2Canopy_pft=0._r8;PrecRainAndIrrig_col=0._r8;SnoFalPrec_col=0._r8
  call reset_fluxes()
  call SyncCanopyThermalInventory(1,1,InitialEnergy)
  Capacity=cpo*gOC_to_m3_OM(375._r8)+cpw*W+cps*S
  call check_close('cold initial energy',InitialEnergy,Capacity*T)
  call check_close('cold previous energy',ENGYX_pft(1,1,1),Capacity*T)
  call CanopyInterceptPrecip(1,1)
  call check_close('buried liquid retained',WatHeldOnCanopy_pft(1,1,1),W)
  call check_close('buried snow retained',SnowOnCanopy_pft(1,1,1),S)
  call check_close('burial is not drainage',HeatCanopyUnload_col(1,1),0._r8)

  ! Existing energy must be replaced by a signed increment, including live-only SP.
  ldo_sp_mode=.true.
  StandDeadStrutElms_pft=0._r8
  PlantPopuLive_pft(1,1,1)=1._r8
  CanopyLeafSheathC_pft(1,1,1)=100._r8
  CanopySapwoodC_pft=0._r8
  ENGYX_pft(1,1,1)=0.2_r8
  ExpectedHeat=(cpo*gOC_to_m3_OM(100._r8)+cpw*W+cps*S)*T
  call SyncCanopyThermalInventory(1,1,InitialEnergy)
  call check_close('SP initial energy increment',InitialEnergy,ExpectedHeat-0.2_r8)
  call check_close('SP prescribed capacity energy',ENGYX_pft(1,1,1),ExpectedHeat)
  call SyncCanopyThermalInventory(1,1,InitialEnergy)
  call check_close('unchanged SP inventory has no source',InitialEnergy,0._r8)
  CanopyLeafSheathC_pft(1,1,1)=110._r8
  call SyncCanopyThermalInventory(1,1,InitialEnergy)
  call check_close('prescribed mass enthalpy gain',InitialEnergy,cpo*gOC_to_m3_OM(10._r8)*T)
  CanopyLeafSheathC_pft(1,1,1)=90._r8
  call SyncCanopyThermalInventory(1,1,InitialEnergy)
  call check_close('prescribed mass enthalpy loss',InitialEnergy,-cpo*gOC_to_m3_OM(20._r8)*T)
  ldo_sp_mode=.false.
  CanopyLeafSheathC_pft(1,1,1)=100._r8
  call SyncCanopyThermalInventory(1,1,InitialEnergy)
  call check_close('full-model live-only material enthalpy gain',InitialEnergy,cpo*gOC_to_m3_OM(10._r8)*T)
  call SyncCanopyThermalInventory(1,1,InitialEnergy)
  call check_close('unchanged full-model inventory has no source',InitialEnergy,0._r8)
  CanopyLeafSheathC_pft(1,1,1)=80._r8
  call SyncCanopyThermalInventory(1,1,InitialEnergy)
  call check_close('full-model material enthalpy loss',InitialEnergy,-cpo*gOC_to_m3_OM(20._r8)*T)
  PlantPopuLive_pft=0._r8

  ! Complete removal. Slot 1 is retained but is no longer inside NP.
  NP_col=0;StandDeadSurfArea_pft=0._r8;StandDeadStrutElms_pft=0._r8
  ! Inactive slots may retain old geometry; it cannot support surface water.
  CanopyLeafArea_pft(1,1,1)=1._r8;LeafStalkAreaAct_pft(1,1,1)=1._r8
  PrecRainAndIrrig_col=3.e-4_r8;SnoFalPrec_col=4.e-4_r8
  call reset_fluxes()
  call CanopyInterceptPrecip(1,1)
  ExpectedHeat=(cpw*W+cps*S)*T
  call check_close('removed liquid released',WatHeldOnCanopy_pft(1,1,1),0._r8)
  call check_close('removed snow released',SnowOnCanopy_pft(1,1,1),0._r8)
  call check_close('liquid mass to ground',RainPrecThrufall_col(1,1),4.e-4_r8)
  call check_close('snow mass to ground',SnowPrecThrufall_col(1,1),6.e-4_r8)
  call check_close('net liquid canopy transfer',RainIntceptByCanopy_col(1,1),-W)
  call check_close('net snow canopy transfer',SnowIntceptByCanopy_col(1,1),-S)
  call check_close('donor heat returned to ground',HeatCanopyUnload_col(1,1),ExpectedHeat)
  call check_close('recipient liquid enthalpy',cpw*TKRainThrufall_col(1,1)*RainPrecThrufall_col(1,1), &
    cpw*(W*T+3.e-4_r8*270._r8))
  call check_close('recipient snow enthalpy',cps*TKSnowThrufall_col(1,1)*SnowPrecThrufall_col(1,1), &
    cps*(S*T+4.e-4_r8*270._r8))
  call check_close('removed wet heat capacity',VHeatCapCanopy_pft(1,1,1),0._r8)
  call reset_fluxes()
  PrecRainAndIrrig_col=0._r8;SnoFalPrec_col=0._r8
  call CanopyInterceptPrecip(1,1)
  call check_close('no repeated release',RainPrecThrufall_col(1,1)+SnowPrecThrufall_col(1,1),0._r8)
  ! The production solver advances 0.99999 h of substeps. A release is a
  ! complete donor amount; its receiving rate must integrate to that amount.
  WatHeldOnCanopy_pft(1,1,1)=W;SnowOnCanopy_pft(1,1,1)=S
  VHeatCapCanopy_pft(1,1,1)=cpw*W+cps*S
  dts_HeatWatTP=.99999_r8/real(NPH,r8)
  call reset_fluxes()
  call CanopyInterceptPrecip(1,1)
  call check_close('integrated liquid release',RainPrecThrufall_col(1,1)*NPH*dts_HeatWatTP,W)
  call check_close('integrated snow release',SnowPrecThrufall_col(1,1)*NPH*dts_HeatWatTP,S)
  call check_close('integrated donor heat', &
    (cpw*TKRainThrufall_col(1,1)*RainPrecThrufall_col(1,1) &
    +cps*TKSnowThrufall_col(1,1)*SnowPrecThrufall_col(1,1))*NPH*dts_HeatWatTP,ExpectedHeat)
  ! Radiation-only mode has no plant solver to accept interception or unloading.
  ! Wet active and retired slots remain inert, even under rain, snow and wind.
  ldo_radiation_test=.true.
  NP_col=1;PlantPopuLive_pft(1,1,1)=1._r8
  LeafStalkAreaAll_col=1._r8;FracPARads2Canopy_pft=1._r8
  WindSpeedAtm_col=1.e5_r8
  WatHeldOnCanopy_pft=W;SnowOnCanopy_pft=S
  VHeatCapCanopy_pft=cpw*W+cps*S
  ExpectedHeat=sum(ENGYX_pft)
  RainIntcptByCanopy_pft=1._r8;SnowIntcptByCanopy_pft=1._r8
  PrecRainAndIrrig_col=3.e-4_r8;SnoFalPrec_col=4.e-4_r8
  call reset_fluxes()
  call CanopyInterceptPrecip(1,1)
  call check_close('radiation rain throughfall',RainPrecThrufall_col(1,1),3.e-4_r8)
  call check_close('radiation snow throughfall',SnowPrecThrufall_col(1,1),4.e-4_r8)
  call check_close('radiation liquid retained',sum(WatHeldOnCanopy_pft),2._r8*W)
  call check_close('radiation snow retained',sum(SnowOnCanopy_pft),2._r8*S)
  call check_close('radiation energy retained',sum(ENGYX_pft),ExpectedHeat)
  call check_close('radiation no rain interception',sum(RainIntcptByCanopy_pft),0._r8)
  call check_close('radiation no snow interception',sum(SnowIntcptByCanopy_pft),0._r8)
  call check_close('radiation no unloaded heat',HeatCanopyUnload_col(1,1),0._r8)
  call check_close('radiation rain temperature',TKRainThrufall_col(1,1),TairK_col(1,1))
  call check_close('radiation snow temperature',TKSnowThrufall_col(1,1),TairK_col(1,1))
  ldo_radiation_test=.false.
  write(iulog,*)'PASS: canopy inventory, conservative drainage and radiation-only throughfall'

contains
!------------------------------------------------------------------------------------------
  subroutine reset_fluxes()
    ! Description:
    ! Reproduce the hourly reset of canopy precipitation accumulators.
    character(len=*), parameter :: subname='reset_fluxes'
    call PrintInfo('beg '//subname)
    Rain2Canopy_col=0._r8;RainIntceptByCanopy_col=0._r8;SnowIntceptByCanopy_col=0._r8
    call PrintInfo('end '//subname)
  end subroutine reset_fluxes
!------------------------------------------------------------------------------------------
  subroutine check_close(label,actual,expected)
    ! Description:
    ! Compare independently derived donor/recipient masses and sensible energies.
    character(len=*), intent(in) :: label ! budget being checked [-]
    real(r8), intent(in) :: actual,expected ! values in the budget units
    character(len=*), parameter :: subname='check_close'
    call PrintInfo('beg '//subname)
    if(ABS(actual-expected)>1.e-12_r8)then
      write(iulog,*)label,actual,expected
      error stop 'canopy store budget failed'
    endif
    call PrintInfo('end '//subname)
  end subroutine check_close
end program CanopyWaterStoreTest
