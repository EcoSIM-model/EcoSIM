module PlantEnergyWaterAPIData
  ! Description:
  ! Owns the energywater API type, instance and lifecycle procedures.
  use DebugToolMod, only: PrintInfo
  use PlantAPICommonData
  implicit none
  private
  public :: plt_ew_init, plt_ew_destroy

  type, public :: plant_ew_type
  logical :: canopy_exposed_col=.false.      ! canopy area protrudes above snow; hourly diagnostic, [-]
  real(r8) :: SnowDepth                     !snowpack depth, [m]
  real(r8) :: VcumWatSnow_col               !water volume in snowpack, [m3 d-2]
  real(r8) :: VcumDrySnoWE_col              !snow volume in snowpack (water equivalent), [m3 d-2]
  real(r8) :: Air_Heat_Sens_store_col       !total sensible heat flux x boundary layer resistance, [MJ m-1]
  real(r8) :: VLHeatCapSnowMin_col          !minimum snowpack heat capacity, [MJ d-2 K-1]
  real(r8) :: VLHeatCapSurfSnow_col         !snowpack heat capacity, [MJ m-3 K-1]
  real(r8) :: CanopyHeatStor_col            !total canopy heat content, [MJ  d-2]
  real(r8) :: LWRadCanG                     !grid total canopy LW emission, [MJ d-2 h-1]
  real(r8) :: QVegET_col                    !total canopy evaporation + transpiration, [m3 d-2]
  real(r8) :: VapXAir2Canopy_col            !grid canopy evaporation, [m3 d-2]
  real(r8) :: HeatSensVapAir2Canopy_col=0._r8 ! accepted atmospheric transfer [MJ/cell/hour]
  real(r8) :: HeatPrecAir2Canopy_col=0._r8 ! accepted atmospheric transfer [MJ/cell/hour]
  real(r8) :: HeatFlx2Canopy_col            !total canopy heat flux, [MJ  d-2]
  real(r8) :: H2OLoss_CumYr_col             !total subsurface water flux, [m3 d-2]
  real(r8) :: VPA                           !vapor concentration, [m3 m-3]
  real(r8) :: EMS_Modify_Scalar_col                !canopy longwave radiation emissivity scalar
  real(r8) :: TairK                         !air temperature, [K]
  real(r8) :: CanopyBiomWater_col                 !total canopy water content stored with dry matter, [m3 d-2]
  real(r8) :: Eco_Heat_Latent_col           !ecosystem latent heat flux, [MJ d-2 h-1]
  real(r8) :: Air_Heat_Latent_store_col     !total latent heat flux x boundary layer resistance, [MJ m-1]
  real(r8) :: VcumIceSnow_col               !ice volume in snowpack, [m3 d-2]
  real(r8) :: TKSnow                        !snow temperature, [K]

  real(r8) :: WatHeldOnCanopy_col           !canopy surface water content, [m3 d-2]
  real(r8) :: fSnowCanopy_col               !fraction snow covered canopy, [-]
  real(r8) :: SnowOnCanopy_col              !canopy snow water content, [m3 d-2]
  real(r8) :: Eco_Heat_Sens_col             !ecosystem sensible heat flux, [MJ d-2 h-1]
  real(r8) :: RoughnessLength                   !canopy surface roughness height, [m]
  real(r8) :: BulkFactor4Snow_col           !grid bulking factor for canopy snow interception effect on radiation, [m2 (kg SWE)-1]
  real(r8) :: ZeroPlaneDisplacem_col        !zero plane displacement height, [m]
  real(r8) :: RawIsoTAtm2CanopySinkZ_col        !isothermal aerodynamic resistance between zero-sink height and wind ref height in atmosphere, [h m-1]
  real(r8) :: RawCanopyH2SinkZ_col           !isothermal aerodynamic resistance bewtween canopy height and zero sink height, [h m-1]
  real(r8) :: RawIsoTSurf2CanopyHScal_col   !scalar for isothermal aerodynamic resistance between zero-sink height and ground surface, [h m-1]
  real(r8) :: RIB                           !Richardson number for calculating boundary layer resistance, [-]
  real(r8) :: Eco_Heat_GrndSurf_col         !ecosystem storage heat flux, [MJ d-2 h-1]
  real(r8) :: CanopyHeatLoss2Dist_col           !canopy energy +/- due to disturbance, [MJ /d2]
  real(r8) :: QCanopyWatLoss2Dist_col           !canopy water +/- due to disturbance, [m3 H2O/d2]
  real(r8) :: TPlantRootH2OUptake_col       !total water uptake by roots, [m3 H2O/d2]
  real(r8), pointer :: CdH2ORootxSoil_pft(:)          => null()    !total root and soil conductance for plant root water uptake, [mH2O h-1 d-2 MPa-1]
  REAL(R8), pointer :: BulkFactor4Snow_pft(:)         => null()    !pft bulking factor for canopy snow interception effect on radiation, [m2 (ton SWE)-1]
  real(r8), pointer :: fSnowCanopy_pft(:)             => null()    !fraction of canopy is snow covered, [-]
  real(r8), pointer :: Transpiration_pft(:)           => null()    !canopy transpiration,                                         [m3 d-2 h-1]
  real(r8), pointer :: PSICanopyTurg_pft(:)           => null()    !plant canopy turgor water potential,                          [MPa]
  real(r8), pointer :: RainIntcptByCanopy_pft(:)      => null()    !water flux into canopy,                                       [m3 d-2 h-1]
  real(r8), pointer :: SnowIntcptByCanopy_pft(:)      => null()    !snow flux into canopy,                                        [m3 d-2 h-1]
  real(r8), pointer :: PSICanopy_pft(:)               => null()    !canopy total water potential,                                 [Mpa]
  real(r8), pointer :: VapXAir2Canopy_pft(:)          => null()    !canopy evaporation+sublimation,                                           [m3 d-2 h-1]
  real(r8), pointer :: VapXAir2CanopyLiq_pft(:)       => null()    !canopy evaporation, [m3 d-2 h-1]
  real(r8), pointer :: HeatStorCanopy_pft(:)          => null()    !canopy storage heat flux,                                     [MJ d-2 h-1]
  real(r8), pointer :: CanopyEvapTransLHeat_pft(:)    => null()    !canopy latent heat flux,                                      [MJ d-2 h-1]
  real(r8), pointer :: RawIsoTCanopy2Atm_pft(:)       => null()    !canopy isothermal boundary later resistance,                  [h m-1]
  real(r8), pointer :: TKS_vr(:)                      => null()    !mean annual soil temperature,                                 [K]
  real(r8), pointer :: PSICanPDailyMin_pft(:)         => null()    !minimum daily canopy water potential,                         [MPa]
  real(r8), pointer :: TdegCCanopy_pft(:)             => null()    !canopy temperature,                                           [oC]
  real(r8), pointer :: DeltaTKC_pft(:)                => null()    !change in canopy temperature,                                 [K]
  real(r8), pointer :: ENGYX_pft(:)                   => null()    !canopy heat storage from previous time step,                  [MJ d-2]
  real(r8), pointer :: TKC_pft(:)                     => null()    !canopy temperature,                                           [K]
  real(r8), pointer :: PSICanopyOsmo_pft(:)           => null()    !canopy osmotic water potential,                               [Mpa]
  real(r8), pointer :: OrganOsmoPsi0pt_pft(:)           => null()    !Organ osmotic potential when canopy water potential = 0 MPa, [MPa]
  real(r8), pointer :: HeatXAir2PCan_pft(:)           => null()    !canopy sensible heat flux,                                    [MJ d-2 h-1]
  real(r8), pointer :: TKCanopy_pft(:)                => null()    !canopy temperature,                                           [K]
  real(r8), pointer :: PSIRoot_pvr(:,:,:)             => null()    !root total water potential,                                   [Mpa]
  real(r8), pointer :: PSIRootOSMO_vr(:,:,:)          => null()    !root osmotic water potential,                                 [Mpa]
  real(r8), pointer :: PSIRootTurg_vr(:,:,:)          => null()    !root turgor water potential,                                  [Mpa]
  real(r8), pointer :: RPlantRootH2OUptk_pvr(:,:,:)   => null()    !root water uptake,                                            [m3 d-2 h-1]
  real(r8), pointer :: SapFlowVlinear_pvr(:,:)        => null()    !Sap flow mean linear velocity, [m h-1]
  real(r8), pointer :: SapFlowVLinear_rpvr(:,:,:)     => null()    !linear sap flow for primary root axis, [m h-1]
  real(r8), pointer :: RootH2OUptkStress_pvr(:,:,:)   => null()    !root water uptake stress indicated by rate,                   [m3 d-2 h-1]
  real(r8), pointer :: THeatLossRoot2Soil_vr(:)       => null()    !total root heat uptake,                                       [MJ d-2]
  real(r8), pointer :: TWaterPlantRoot2Soil_vr(:)     => null()    !total root water uptake,                                      [m3 d-2]
  real(r8), pointer :: SnowOnCanopy_pft(:)            => null()    !canopy held snow, [m3 d-2]
  real(r8), pointer :: SnoSub2AirCanopy_pft(:)        => null()    !canopy snow sublimation,[m3 d-2 h-1]
  real(r8), pointer :: WatHeldOnCanopy_pft(:)         => null()    !canopy surface water content,                                 [m3 d-2]
  real(r8), pointer :: CanopyBiomWater_pft(:)         => null()    !canopy water content,                                         [m3 d-2]
  real(r8), pointer :: VHeatCapCanopy_pft(:)          => null()    !canopy heat capacity,                                         [MJ d-2 K-1]
  real(r8), pointer :: ElvAdjstedSoilH2OPSIMPa_vr(:)  => null()    !soil micropore total water potential,                         [MPa]
  real(r8), pointer :: ETCanopy_CumYr_pft(:)          => null()    !total transpiration,                                          [m H2O d-2]
  real(r8), pointer :: QdewCanopy_pft(:)              => null()    !dew fall on to canopy,                                        [m3 H2O d-2 h-1]
  real(r8), pointer :: RootResist4H2O_pvr(:,:,:)      => null()    !total root (axial+radial) resistance for water uptake,        [MPa-1 h-1]
  real(r8), pointer :: RootRadialKond2H2O_pvr(:,:,:)  => null()    !radial root conductance for water uptake, [m3 H2O h-1 MPa-1]
  real(r8), pointer :: RootAxialKond2H2O_pvr(:,:,:)   => null()    !axial root conductance for water uptake, [m3 H2O h-1 MPa-1]

  contains
    procedure, public :: Init => plt_ew_init
    procedure, public :: Destroy=> plt_ew_destroy
  end type plant_ew_type


  type(plant_ew_type)       , public, target :: plt_ew        !plant energy and water type

contains

  subroutine  plt_ew_init(this)
  ! Description:
  ! Initialize the plant energy/water workspace and derived canopy exposure.
  implicit none
  class(plant_ew_type) :: this

  this%canopy_exposed_col=.false.
  allocate(this%RootResist4H2O_pvr(jroots,JZ1,JP1)); this%RootResist4H2O_pvr=spval
  allocate(this%RootRadialKond2H2O_pvr(jroots,JZ1,JP1));this%RootRadialKond2H2O_pvr=spval
  allocate(this%RootAxialKond2H2O_pvr(jroots,JZ1,JP1));this%RootAxialKond2H2O_pvr=spval
  allocate(this%QdewCanopy_pft(JP1)); this%QdewCanopy_pft=spval
  allocate(this%ETCanopy_CumYr_pft(JP1));this%ETCanopy_CumYr_pft=spval
  allocate(this%ElvAdjstedSoilH2OPSIMPa_vr(0:JZ1));this%ElvAdjstedSoilH2OPSIMPa_vr=spval
  allocate(this%THeatLossRoot2Soil_vr(0:JZ1));this%THeatLossRoot2Soil_vr=spval
  allocate(this%TKCanopy_pft(JP1));this%TKCanopy_pft=spval
  allocate(this%HeatXAir2PCan_pft(JP1));this%HeatXAir2PCan_pft=spval
  allocate(this%RainIntcptByCanopy_pft(JP1));this%RainIntcptByCanopy_pft=spval
  allocate(this%SnowIntcptByCanopy_pft(JP1));this%SnowIntcptByCanopy_pft=spval
  allocate(this%PSICanopyTurg_pft(JP1));this%PSICanopyTurg_pft=spval
  allocate(this%PSICanopy_pft(JP1));this%PSICanopy_pft=spval
  allocate(this%VapXAir2Canopy_pft(JP1));this%VapXAir2Canopy_pft=spval
  allocate(this%VapXAir2CanopyLiq_pft(JP1));this%VapXAir2CanopyLiq_pft=spval
  allocate(this%HeatStorCanopy_pft(JP1));this%HeatStorCanopy_pft=spval
  allocate(this%CanopyEvapTransLHeat_pft(JP1));this%CanopyEvapTransLHeat_pft=spval
  allocate(this%WatHeldOnCanopy_pft(JP1));this%WatHeldOnCanopy_pft=spval
  allocate(this%SnowOnCanopy_pft(JP1)); this%SnowOnCanopy_pft=spval
  allocate(this%SnoSub2AirCanopy_pft(JP1));this%SnoSub2AirCanopy_pft=spval
  allocate(this%VHeatCapCanopy_pft(JP1));this%VHeatCapCanopy_pft=spval
  allocate(this%CanopyBiomWater_pft(JP1));this%CanopyBiomWater_pft=spval
  allocate(this%PSIRoot_pvr(jroots,JZ1,JP1));this%PSIRoot_pvr=spval
  allocate(this%PSIRootOSMO_vr(jroots,JZ1,JP1));this%PSIRootOSMO_vr=spval
  allocate(this%PSIRootTurg_vr(jroots,JZ1,JP1));this%PSIRootTurg_vr=spval
  allocate(this%RPlantRootH2OUptk_pvr(jroots,JZ1,JP1));this%RPlantRootH2OUptk_pvr=spval
  allocate(this%SapFlowVlinear_pvr(JZ1,JP1)); this%SapFlowVlinear_pvr=0._r8
  allocate(this%SapFlowVLinear_rpvr(JZ1,MaxNumRootAxes,JP1)); this%SapFlowVLinear_rpvr=0._r8
  allocate(this%RootH2OUptkStress_pvr(jroots,JZ1,JP1)); this%RootH2OUptkStress_pvr=spval
  allocate(this%TWaterPlantRoot2Soil_vr(0:JZ1));this%TWaterPlantRoot2Soil_vr=spval
  allocate(this%Transpiration_pft(JP1));this%Transpiration_pft=spval
  allocate(this%CdH2ORootxSoil_pft(JP1)); this%CdH2ORootxSoil_pft=0._r8
  allocate(this%fSnowCanopy_pft(JP1)); this%fSnowCanopy_pft=spval
  allocate(this%BulkFactor4Snow_pft(JP1)); this%BulkFactor4Snow_pft=spval
  allocate(this%PSICanopyOsmo_pft(JP1));this%PSICanopyOsmo_pft=spval
  allocate(this%TKS_vr(0:JZ1));this%TKS_vr=spval
  allocate(this%OrganOsmoPsi0pt_pft(JP1));this%OrganOsmoPsi0pt_pft=spval
  allocate(this%RawIsoTCanopy2Atm_pft(JP1));this%RawIsoTCanopy2Atm_pft=spval
  allocate(this%DeltaTKC_pft(JP1));this%DeltaTKC_pft=spval
  allocate(this%TKC_pft(JP1));this%TKC_pft=spval
  allocate(this%ENGYX_pft(JP1));this%ENGYX_pft=spval
  allocate(this%TdegCCanopy_pft(JP1));this%TdegCCanopy_pft=spval
  allocate(this%PSICanPDailyMin_pft(JP1));this%PSICanPDailyMin_pft=spval

  end subroutine plt_ew_init

  subroutine plt_ew_destroy(this)
  ! Description:
  ! Release every array owned by this PlantEnergyWaterAPIData workspace.
  class(plant_ew_type), intent(inout) :: this ! owned API state; component units and dimensions as declared above
  character(len=*), parameter :: subname='plt_ew_destroy'

  call PrintInfo('beg '//subname)
  if(associated(this%CdH2ORootxSoil_pft))deallocate(this%CdH2ORootxSoil_pft);nullify(this%CdH2ORootxSoil_pft)
  if(associated(this%BulkFactor4Snow_pft))deallocate(this%BulkFactor4Snow_pft);nullify(this%BulkFactor4Snow_pft)
  if(associated(this%fSnowCanopy_pft))deallocate(this%fSnowCanopy_pft);nullify(this%fSnowCanopy_pft)
  if(associated(this%Transpiration_pft))deallocate(this%Transpiration_pft);nullify(this%Transpiration_pft)
  if(associated(this%PSICanopyTurg_pft))deallocate(this%PSICanopyTurg_pft);nullify(this%PSICanopyTurg_pft)
  if(associated(this%RainIntcptByCanopy_pft))deallocate(this%RainIntcptByCanopy_pft);nullify(this%RainIntcptByCanopy_pft)
  if(associated(this%SnowIntcptByCanopy_pft))deallocate(this%SnowIntcptByCanopy_pft);nullify(this%SnowIntcptByCanopy_pft)
  if(associated(this%PSICanopy_pft))deallocate(this%PSICanopy_pft);nullify(this%PSICanopy_pft)
  if(associated(this%VapXAir2Canopy_pft))deallocate(this%VapXAir2Canopy_pft);nullify(this%VapXAir2Canopy_pft)
  if(associated(this%VapXAir2CanopyLiq_pft))deallocate(this%VapXAir2CanopyLiq_pft);nullify(this%VapXAir2CanopyLiq_pft)
  if(associated(this%HeatStorCanopy_pft))deallocate(this%HeatStorCanopy_pft);nullify(this%HeatStorCanopy_pft)
  if(associated(this%CanopyEvapTransLHeat_pft))deallocate(this%CanopyEvapTransLHeat_pft);nullify(this%CanopyEvapTransLHeat_pft)
  if(associated(this%RawIsoTCanopy2Atm_pft))deallocate(this%RawIsoTCanopy2Atm_pft);nullify(this%RawIsoTCanopy2Atm_pft)
  if(associated(this%TKS_vr))deallocate(this%TKS_vr);nullify(this%TKS_vr)
  if(associated(this%PSICanPDailyMin_pft))deallocate(this%PSICanPDailyMin_pft);nullify(this%PSICanPDailyMin_pft)
  if(associated(this%TdegCCanopy_pft))deallocate(this%TdegCCanopy_pft);nullify(this%TdegCCanopy_pft)
  if(associated(this%DeltaTKC_pft))deallocate(this%DeltaTKC_pft);nullify(this%DeltaTKC_pft)
  if(associated(this%ENGYX_pft))deallocate(this%ENGYX_pft);nullify(this%ENGYX_pft)
  if(associated(this%TKC_pft))deallocate(this%TKC_pft);nullify(this%TKC_pft)
  if(associated(this%PSICanopyOsmo_pft))deallocate(this%PSICanopyOsmo_pft);nullify(this%PSICanopyOsmo_pft)
  if(associated(this%OrganOsmoPsi0pt_pft))deallocate(this%OrganOsmoPsi0pt_pft);nullify(this%OrganOsmoPsi0pt_pft)
  if(associated(this%HeatXAir2PCan_pft))deallocate(this%HeatXAir2PCan_pft);nullify(this%HeatXAir2PCan_pft)
  if(associated(this%TKCanopy_pft))deallocate(this%TKCanopy_pft);nullify(this%TKCanopy_pft)
  if(associated(this%PSIRoot_pvr))deallocate(this%PSIRoot_pvr);nullify(this%PSIRoot_pvr)
  if(associated(this%PSIRootOSMO_vr))deallocate(this%PSIRootOSMO_vr);nullify(this%PSIRootOSMO_vr)
  if(associated(this%PSIRootTurg_vr))deallocate(this%PSIRootTurg_vr);nullify(this%PSIRootTurg_vr)
  if(associated(this%RPlantRootH2OUptk_pvr))deallocate(this%RPlantRootH2OUptk_pvr);nullify(this%RPlantRootH2OUptk_pvr)
  if(associated(this%SapFlowVlinear_pvr))deallocate(this%SapFlowVlinear_pvr);nullify(this%SapFlowVlinear_pvr)
  if(associated(this%SapFlowVLinear_rpvr))deallocate(this%SapFlowVLinear_rpvr);nullify(this%SapFlowVLinear_rpvr)
  if(associated(this%RootH2OUptkStress_pvr))deallocate(this%RootH2OUptkStress_pvr);nullify(this%RootH2OUptkStress_pvr)
  if(associated(this%THeatLossRoot2Soil_vr))deallocate(this%THeatLossRoot2Soil_vr);nullify(this%THeatLossRoot2Soil_vr)
  if(associated(this%TWaterPlantRoot2Soil_vr))deallocate(this%TWaterPlantRoot2Soil_vr);nullify(this%TWaterPlantRoot2Soil_vr)
  if(associated(this%SnowOnCanopy_pft))deallocate(this%SnowOnCanopy_pft);nullify(this%SnowOnCanopy_pft)
  if(associated(this%SnoSub2AirCanopy_pft))deallocate(this%SnoSub2AirCanopy_pft);nullify(this%SnoSub2AirCanopy_pft)
  if(associated(this%WatHeldOnCanopy_pft))deallocate(this%WatHeldOnCanopy_pft);nullify(this%WatHeldOnCanopy_pft)
  if(associated(this%CanopyBiomWater_pft))deallocate(this%CanopyBiomWater_pft);nullify(this%CanopyBiomWater_pft)
  if(associated(this%VHeatCapCanopy_pft))deallocate(this%VHeatCapCanopy_pft);nullify(this%VHeatCapCanopy_pft)
  if(associated(this%ElvAdjstedSoilH2OPSIMPa_vr))deallocate(this%ElvAdjstedSoilH2OPSIMPa_vr);nullify(this%ElvAdjstedSoilH2OPSIMPa_vr)
  if(associated(this%ETCanopy_CumYr_pft))deallocate(this%ETCanopy_CumYr_pft);nullify(this%ETCanopy_CumYr_pft)
  if(associated(this%QdewCanopy_pft))deallocate(this%QdewCanopy_pft);nullify(this%QdewCanopy_pft)
  if(associated(this%RootResist4H2O_pvr))deallocate(this%RootResist4H2O_pvr);nullify(this%RootResist4H2O_pvr)
  if(associated(this%RootRadialKond2H2O_pvr))deallocate(this%RootRadialKond2H2O_pvr);nullify(this%RootRadialKond2H2O_pvr)
  if(associated(this%RootAxialKond2H2O_pvr))deallocate(this%RootAxialKond2H2O_pvr);nullify(this%RootAxialKond2H2O_pvr)
  call PrintInfo('end '//subname)
  end subroutine plt_ew_destroy
end module PlantEnergyWaterAPIData
