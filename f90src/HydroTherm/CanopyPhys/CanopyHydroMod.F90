module CanopyHydroMod
  ! Description:
  ! Partition precipitation and canopy snow unloading using the exposed live
  ! and standing-dead surface area, preserving water routed to the ground.
  ! Radiation-only runs keep canopy stores inert and use full throughfall.
  ! Synchronize canopy thermal inventories before hourly exchange.
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use EcoSIMCtrlMod, only: ldo_radiation_test,ldo_sp_mode
  use abortutils,        only: endrun
  use MiniMathMod,       only: AZMAX1
  use PlantMgmtDataType, only: NP_col,NP0_col
  use ElmIDMod, only: ielmc
  use EcoSimConst, only: cpw,cps,cpo
  use EcoSIMSolverPar, only: NPH,dts_HeatWatTP
  use MiniFuncMod, only: gOC_to_m3_OM
  use PlantMathFuncMod, only: CanopyStemThermalCarbon,CanopySnowOptics
  use PlantBGCPars, only: CanopySnowCapacityPerArea,CanopySnowAreaCoeff
  USE EcoSIMCtrlDataType,ONLY: ZEROS
  use DebugToolMod  
  USE SoilWaterDataType
  use PlantTraitDataType
  use CanopyDataType
  use ClimForcDataType
implicit none

  private
  character(len=*), parameter :: mod_filename=&
  __FILE__
  public :: CanopyInterceptPrecip,SyncCanopyThermalInventory

  real(r8), parameter :: FoliarWatRetcap(0:3)=real((/5.0E-04,2.5E-04,2.5E-04,2.5E-04/),r8)                !bryophyte, grasses, deep roots, very deep roots
  contains
!------------------------------------------------------------------------------------------

  subroutine CanopyInterceptPrecip(NY,NX)
  ! Description:
  ! Route precipitation to ground in radiation-only mode; otherwise intercept
  ! it on exposed live/dead canopy and route excess to ground.
  implicit none
  integer, intent(in) :: NY,NX
  character(len=*), parameter :: subname='CanopyInterceptPrecip'
  integer :: NZ
  logical :: HasCanopyMaterial
  real(r8) :: CanopyWatHeldCap,CanopySnowHeldCap  !maximum rain/snow holding capacity by canopy (leaf+stem) [m3 H2O]~ [ton H2O]
  real(r8) :: rain2canopy_pft,snow2canopy_pft     !rain/snowfall onto canopy [m H2O/h]
  real(r8) :: SLArea                              !stem+leaf surface area [m2]
  real(r8) :: SnowWindUnload,SnowTempUnload,CanopySnowUnload_col,SnowUnload_pft
  real(r8) :: SnowRadiativeArea ! physical plus retained-snow interception area [m2/cell]
  real(r8) :: ENGYS,ENGYW,ReleasedWater,ReleasedSnow,ReleasedHeat,ReleaseInterval
!
!     CANOPY RETENTION OF PRECIPITATION
!
!     FoliarWatRetcap=foliar surface water retention capacity
!     CanopyLeafArea_pft,CanopyStemSurfArea_pft=leaf,stalk area of PFT
!     FLWC,TFLWC=water retention of PFT,combined canopy
!     PRECA=precipitation+irrigation
!     FracPARads2Canopy_pft=fraction of radiation received by each PFT canopy
!     VOLWC=canopy surface water retention
!
!     Snow and liquid interception have separate retention capacities.
  call PrintInfo('beg '//subname)
  ! These substeps deliberately total 0.99999 h in the legacy solver. Convert
  ! complete hourly donor releases to rates over the actual receiving interval.
  ReleaseInterval=real(NPH,r8)*dts_HeatWatTP
  if(ReleaseInterval<=0._r8)call endrun('Invalid canopy unloading interval',__LINE__)
  CanopySnowUnload_col       = 0._r8
  BulkFactor4Snow_col(NY,NX) = 0.0_r8
  fSnowCanopy_col(NY,NX)     = 0._r8
  ENGYS                      = 0._R8
  ENGYW                      = 0._r8
  RainCanopyUnload_col(NY,NX) = 0._r8
  SnowCanopyUnload_col(NY,NX) = 0._r8
  HeatCanopyUnload_col(NY,NX) = 0._r8
  DO NZ=1,NP0_col(NY,NX)

    ! This diagnostic sends radiation directly to ground and skips the plant
    ! energy/water solve. Route precipitation there too; retain old canopy stocks.
    if(ldo_radiation_test)then
      RainIntcptByCanopy_pft(NZ,NY,NX)=0._r8
      SnowIntcptByCanopy_pft(NZ,NY,NX)=0._r8
      cycle
    endif

    ! Total material, rather than exposed area, distinguishes removal from burial.
    ! Population also excludes residual live geometry left in an inactive slot.
    HasCanopyMaterial=StandDeadStrutElms_pft(ielmc,NZ,NY,NX)>0._r8.or. &
      (PlantPopuLive_pft(NZ,NY,NX)>0._r8.and. &
       CanopyLeafArea_pft(NZ,NY,NX)+CanopyStemSurfArea_pft(NZ,NY,NX)>ZEROS(NY,NX))
    if(.not.HasCanopyMaterial)then
      ReleasedWater=WatHeldOnCanopy_pft(NZ,NY,NX)
      ReleasedSnow=SnowOnCanopy_pft(NZ,NY,NX)
      ReleasedHeat=(cpw*ReleasedWater+cps*ReleasedSnow)*TKC_pft(NZ,NY,NX)
      RainCanopyUnload_col(NY,NX)=RainCanopyUnload_col(NY,NX)+ReleasedWater
      CanopySnowUnload_col=CanopySnowUnload_col+ReleasedSnow
      ENGYW=ENGYW+ReleasedWater*TKC_pft(NZ,NY,NX)
      ENGYS=ENGYS+ReleasedSnow*TKC_pft(NZ,NY,NX)
      HeatCanopyUnload_col(NY,NX)=HeatCanopyUnload_col(NY,NX)+ReleasedHeat
      WatHeldOnCanopy_pft(NZ,NY,NX)=0._r8
      SnowOnCanopy_pft(NZ,NY,NX)=0._r8
      VHeatCapCanopy_pft(NZ,NY,NX)=cpw*CanopyBiomWater_pft(NZ,NY,NX)
    endif

    SLArea=LeafStalkAreaAct_pft(NZ,NY,NX)+StandDeadSurfAreaAct_pft(NZ,NY,NX)
    if(HasCanopyMaterial.and.NZ<=NP_col(NY,NX).and.SLArea.GT.ZEROS(NY,NX))then
      CanopyWatHeldCap  = FoliarWatRetcap(iPlantRootProfile_pft(NZ,NY,NX))*SLArea
      CanopySnowHeldCap = CanopySnowCapacityPerArea(iPlantSnowIntercepType_pft(NZ,NY,NX))*SLArea

      snow2canopy_pft                  = AZMAX1(SnoFalPrec_col(NY,NX)*FracPARads2Canopy_pft(NZ,NY,NX))
      rain2canopy_pft                  = AZMAX1(PrecRainAndIrrig_col(NY,NX)*FracPARads2Canopy_pft(NZ,NY,NX))
      SnowIntcptByCanopy_pft(NZ,NY,NX) = AZMAX1(AMIN1(snow2canopy_pft,CanopySnowHeldCap-SnowOnCanopy_pft(NZ,NY,NX)))
      RainIntcptByCanopy_pft(NZ,NY,NX) = AZMAX1(AMIN1(rain2canopy_pft,CanopyWatHeldCap-WatHeldOnCanopy_pft(NZ,NY,NX)))
      Rain2Canopy_col(NY,NX)           = Rain2Canopy_col(NY,NX)+rain2canopy_pft
      RainIntceptByCanopy_col(NY,NX)   = RainIntceptByCanopy_col(NY,NX)+RainIntcptByCanopy_pft(NZ,NY,NX)
      SnowIntceptByCanopy_col(NY,NX)    = SnowIntceptByCanopy_col(NY,NX)+SnowIntcptByCanopy_pft(NZ,NY,NX)

      !wind-induced unloading
      SnowWindUnload= WindSpeedAtm_col(NY,NX)*SnowOnCanopy_pft(NZ,NY,NX)/1.56E5_r8
      !temperature-induced unloading
      SnowTempUnload = SnowOnCanopy_pft(NZ,NY,NX)*AZMAX1(TKCanopy_pft(NZ,NY,NX)-270._r8)/(1.87E5_r8*3600._r8)

      !update by unloading
      SnowUnload_pft                = AMIN1(SnowWindUnload+SnowTempUnload,SnowOnCanopy_pft(NZ,NY,NX))
      SnowOnCanopy_pft(NZ,NY,NX)    = AZMAX1(SnowOnCanopy_pft(NZ,NY,NX)-SnowUnload_pft)
      ENGYS                         = ENGYS+SnowUnload_pft*TKC_pft(NZ,NY,NX)
      HeatCanopyUnload_col(NY,NX) = HeatCanopyUnload_col(NY,NX)+cps*SnowUnload_pft*TKC_pft(NZ,NY,NX)
      
      CanopySnowUnload_col          = CanopySnowUnload_col+SnowUnload_pft
      call CanopySnowOptics(iPlantSnowIntercepType_pft(NZ,NY,NX),SLArea, &
        SnowOnCanopy_pft(NZ,NY,NX)+SnowIntcptByCanopy_pft(NZ,NY,NX), &
        fSnowCanopy_pft(NZ,NY,NX),SnowRadiativeArea)
      BulkFactor4Snow_pft(NZ,NY,NX) = CanopySnowAreaCoeff(iPlantSnowIntercepType_pft(NZ,NY,NX))
      BulkFactor4Snow_col(NY,NX)    = BulkFactor4Snow_col(NY,NX)+BulkFactor4Snow_pft(NZ,NY,NX)*SLArea
      fSnowCanopy_col(NY,NX)        = fSnowCanopy_col(NY,NX)+fSnowCanopy_pft(NZ,NY,NX)*SLArea
    else
      SnowIntcptByCanopy_pft(NZ,NY,NX) = 0._r8
      RainIntcptByCanopy_pft(NZ,NY,NX) = 0._r8
      fSnowCanopy_pft(NZ,NY,NX)        = 0._r8
      BulkFactor4Snow_pft(NZ,NY,NX)    = 0._R8
    endif
  ENDDO
  ReleasedWater=RainCanopyUnload_col(NY,NX)
  RainCanopyUnload_col(NY,NX)=ReleasedWater/ReleaseInterval
  SnowCanopyUnload_col(NY,NX)=CanopySnowUnload_col/ReleaseInterval
  RainPrecThrufall_col(NY,NX) = PrecRainAndIrrig_col(NY,NX)-RainIntceptByCanopy_col(NY,NX)+RainCanopyUnload_col(NY,NX)
  SnowPrecThrufall_col(NY,NX) = SnoFalPrec_col(NY,NX)-SnowIntceptByCanopy_col(NY,NX)+SnowCanopyUnload_col(NY,NX)
  IF(ENGYS.GT.ZEROS(NY,NX))THEN
    TKSnowThrufall_col(NY,NX) = (TairK_col(NY,NX)*(SnoFalPrec_col(NY,NX)-SnowIntceptByCanopy_col(NY,NX))+ENGYS/ReleaseInterval)/SnowPrecThrufall_col(NY,NX)
  else
    TKSnowThrufall_col(NY,NX) = TairK_col(NY,NX)
  endif
  TKRainThrufall_col(NY,NX)=TairK_col(NY,NX)
  if(RainPrecThrufall_col(NY,NX)>0._r8)then
    TKRainThrufall_col(NY,NX)=(TairK_col(NY,NX)*(PrecRainAndIrrig_col(NY,NX) &
      -RainIntceptByCanopy_col(NY,NX))+ENGYW/ReleaseInterval)/RainPrecThrufall_col(NY,NX)
  endif
  ! Net canopy transfers include unloading; ground throughfall already includes it.
  RainIntceptByCanopy_col(NY,NX)=RainIntceptByCanopy_col(NY,NX)-ReleasedWater
  SnowIntceptByCanopy_col(NY,NX)=SnowIntceptByCanopy_col(NY,NX)-CanopySnowUnload_col
  if(LeafStalkAreaAll_col(NY,NX).GT.ZEROS(NY,NX))THEN
    BulkFactor4Snow_col(NY,NX) = BulkFactor4Snow_col(NY,NX)/LeafStalkAreaAll_col(NY,NX)
    fSnowCanopy_col(NY,NX)     = fSnowCanopy_col(NY,NX) /LeafStalkAreaAll_col(NY,NX)
  ENDIF  
  call PrintInfo('end '//subname)
  end subroutine CanopyInterceptPrecip

!------------------------------------------------------------------------------------------
  subroutine SyncCanopyThermalInventory(NY,NX,EnergyAdjustment)
    ! Description:
    ! Rebuild capacity at the retained temperature from current live/dead material
    ! and water before interception. Return the signed stock change: the caller
    ! initializes cold-start inventories or books changed material enthalpy.
    integer, intent(in) :: NY,NX ! grid column indices [-]
    real(r8), intent(out) :: EnergyAdjustment ! signed canopy sensible inventory adjustment [MJ/cell]
    character(len=*), parameter :: subname='SyncCanopyThermalInventory'
    integer :: NZ
    real(r8) :: Carbon,Capacity,Temperature
    call PrintInfo('beg '//subname)
    EnergyAdjustment=0._r8
    do NZ=1,NP0_col(NY,NX)
      if(StandDeadStrutElms_pft(ielmc,NZ,NY,NX)<=0._r8 .and. &
        PlantPopuLive_pft(NZ,NY,NX)<=0._r8 .and. &
        .not.(ldo_sp_mode.and..not.ldo_radiation_test))cycle
      Carbon=CanopyStemThermalCarbon(StandDeadStrutElms_pft(ielmc,NZ,NY,NX), &
        StandDeadSurfArea_pft(NZ,NY,NX),StemSpecVolume_pft(NZ,NY,NX))
      if(PlantPopuLive_pft(NZ,NY,NX)>0._r8)then
        Carbon=Carbon+CanopyLeafSheathC_pft(NZ,NY,NX)+CanopyStemThermalCarbon( &
          CanopySapwoodC_pft(NZ,NY,NX),CanopyStemSurfArea_pft(NZ,NY,NX),StemSpecVolume_pft(NZ,NY,NX))
      endif
      Capacity=cpo*gOC_to_m3_OM(Carbon)+cpw*(CanopyBiomWater_pft(NZ,NY,NX) &
        +WatHeldOnCanopy_pft(NZ,NY,NX))+cps*SnowOnCanopy_pft(NZ,NY,NX)
      Temperature=TKC_pft(NZ,NY,NX)
      if(Temperature<=0._r8)Temperature=TairK_col(NY,NX)
      TKC_pft(NZ,NY,NX)=Temperature
      TKCanopy_pft(NZ,NY,NX)=Temperature
      VHeatCapCanopy_pft(NZ,NY,NX)=Capacity
      EnergyAdjustment=EnergyAdjustment+Capacity*Temperature-ENGYX_pft(NZ,NY,NX)
      ENGYX_pft(NZ,NY,NX)=Capacity*Temperature
    enddo
    call PrintInfo('end '//subname)
  end subroutine SyncCanopyThermalInventory

end module CanopyHydroMod
