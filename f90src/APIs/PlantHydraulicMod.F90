module PlantHydraulicMod
  ! Description:
  ! Advance plant water and energy exchange independently of plant
  ! growth, nutrient uptake and litter accounting. The soil driver consumes the
  ! accepted layer water and donor enthalpy before advancing soil transport.
  use data_kind_mod, only: r8 => DAT_KIND_R8,yearIJ_type
  use DebugToolMod, only: PrintInfo
  use PlantAPICommonData, only: JP1,JZ1,pltpar
  use PlantAPI4Uptake, only: PlantUptakeAPISend,PlantUptakeAPIRecv
  use PlantSiteAPIData, only: plt_site
  use EcoSIMCtrlMod, only: ldo_sp_mode
  use PlantMod, only: PreparePlantColumn
  use PlantDisturbMod, only: PrepLandscapeGrazing
  use PlantPhaseWorkMod
  use UptakesMod, only: RootWaterUptake
  use ExtractsMod, only: ExtractRootWaterHeat,ExtractCanopyEnergyWater,ExtractCanopyInventory
  implicit none
  private
  public :: PlantHydraulicModel
contains
!------------------------------------------------------------------------------------------
  subroutine PlantHydraulicModel(yearIJ,NHW,NHE,NVN,NVS)
    ! Description:
    ! Solve one hourly canopy/root exchange and aggregate only water/energy.
    ! Full mode retains prepared column state for its later chemical/BGC phase.
    type(yearIJ_type), intent(in) :: yearIJ ! model year, day of year and hour [-]
    integer, intent(in) :: NHW,NHE,NVN,NVS ! active west/east/north/south column bounds [-]
    character(len=*), parameter :: subname='PlantHydraulicModel'
    logical :: LiveExchange(JP1)
    logical :: PassiveExchange(JP1)
    real(r8) :: CanopyDryMatter_pft(JP1)                          ! accepted dry-matter fraction [-]
    real(r8) :: RadialMeanLen_pvr(pltpar%jroots,JZ1,JP1)          !
    real(r8) :: FineRootRadius_pvr(pltpar%jroots,JZ1,JP1)         ! [m]
    real(r8) :: FracPRoot4Uptake_pvr(pltpar%jroots,JZ1,JP1)       !
    real(r8) :: FracSoilLBy1stRoots_pvr(JZ1,JP1)                  ! [-]
    real(r8) :: RootEffLen4Absorption_rpvr(pltpar%jroots,JZ1,JP1) ! [m]
    integer :: NY,NX,NZ

    call PrintInfo('beg '//subname)
    associate(               &
      I => yearIJ%I        , & ! input: model day of year [-]
      J => yearIJ%J        , & ! input: model hour of day [-]
      NP => plt_site%NP    , & ! input: current PFT count in the column [-]
      NP0 => plt_site%NP0    & ! input: retained PFT slot count in the column [-]
    )
    if(.not.ldo_sp_mode)then
      call InitPlantPhaseWork(NHW,NHE,NVN,NVS)
      call PrepLandscapeGrazing(I,J,NHW,NHE,NVN,NVS)
    endif

    do NX=NHW,NHE
      do NY=NVN,NVS
        if(ldo_sp_mode)then
          call PlantUptakeAPISend(I,J,NY,NX)
        else
          call SwapPlantPhaseColumn(NY,NX)
          call PreparePlantColumn(yearIJ,NY,NX)
        endif
        call RootWaterUptake(yearIJ,LiveExchange,PassiveExchange,CanopyDryMatter_pft, &
          RadialMeanLen_pvr,FineRootRadius_pvr,FracPRoot4Uptake_pvr,FracSoilLBy1stRoots_pvr,RootEffLen4Absorption_rpvr)

        do NZ=1,NP
          if(LiveExchange(NZ))call ExtractRootWaterHeat(NZ)
          if(LiveExchange(NZ).or.PassiveExchange(NZ))call ExtractCanopyEnergyWater(NZ)
        enddo
        do NZ=1,NP0
          call ExtractCanopyInventory(NZ)
        enddo
        call PlantUptakeAPIRecv(I,J,NY,NX)
        if(.not.ldo_sp_mode)then
          WaterLive_pft(:,NY,NX)           = LiveExchange
          WaterPassive_pft(:,NY,NX)        = PassiveExchange
          WaterDryMatter_pft(:,NY,NX)      = CanopyDryMatter_pft
          WaterPath_pvr(:,:,:,NY,NX)       = RadialMeanLen_pvr
          WaterRadius_pvr(:,:,:,NY,NX)     = FineRootRadius_pvr
          WaterRootShare_pvr(:,:,:,NY,NX)  = FracPRoot4Uptake_pvr
          WaterRootDepth_pvr(:,:,NY,NX)    = FracSoilLBy1stRoots_pvr
          WaterRootLength_pvr(:,:,:,NY,NX) = RootEffLen4Absorption_rpvr
          call SavePlantHydraulicColumn(yearIJ,NY,NX)
        endif
      enddo
    enddo
    end associate
    call PrintInfo('end '//subname)
  end subroutine PlantHydraulicModel
end module PlantHydraulicMod
