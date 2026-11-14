module EcoSIMDesctruct
  ! Release model state and transport workspaces at the end of a simulation.

  implicit none

  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: DestructEcoSIM
  contains

  subroutine DestructEcoSIM
  ! Pair model and workspace cleanup with the corresponding initialization paths.
  use EcoSimSumDataType   , only : DestructEcoSimSum
  use FlagDataType        , only : DestructFlagData
  use WatsubDataMod       , only : DestructWatSubData
  use ErosionMod          , only : DestructErosion
  use SoilHeatDatatype    , only : DestructSoilHeatData
  use SoilPhysDataType    , only : DestructSoilPhysData
  use SurfSoilDataType    , only : DestructSurfSoilData
  use GridDataType        , only : DestructGridData
  use RootDataType        , only : DestructRootData
  use SnowDataType        , only : DestructSnowData
  use PlantTraitDataType  , only : DestructPlantTraits
  use LandSurfDataType    , only : DestructLandSurfData
  use EcoSIMCtrlDataType  , only : DestructEcoSIMCtrlData
  use EcoSIMCtrlMod       , only : pft_nfid,salt_model, plant_model
  use ncdio_pio           , only : ncd_pio_closefile
  use EcosimBGCFluxType   , only : DestructEcosimBGCFluxData
  use EcoSIMHistMod       , only : DestructEcoSIMHistData
  use CanopyDataType      , only : DestructCanopyData
  use ClimForcDataType    , only : DestructClimForcData
  use AqueChemDatatype    , only : DestructAquaChem
  use MicrobialDataType   , only : DestructMicrobialData
  use ChemTranspDataType  , only : DestructChemTranspData
  use SOMDataType         , only : DestructSOMData
  use FertilizerDataType  , only : DestructFertilizerData
  use CanopyRadDataType   , only : DestructCanopyRad
  use SurfLitterDataType  , only : DestructSurfLitter
  use SoilPropertyDataType, only : DestructSoilProperty
  use PlantDataRateType   , only : DestructPlantRates
  use IrrigationDataType  , only : DestructIrrigation
  use SoilBGCDataType     , only : DestructSoilBGCData
  use SedimentDataType    , only : DestructSedimentData
  use SoilWaterDataType   , only : DestructSoilWater
  use PlantAPIData        , only : DestructPlantAPIData
  use PlantPhaseWorkMod   , only : DestroyPlantPhaseWork
  use PlantMgmtDataType   , only : DestructPlantMngmtData
  use InitSOMBGCMOD       , only : DestructSOMBGC
  use TranspNoSaltMod     , only : DestructTranspNoSalt
  use TranspSaltDataMod   , only : DestructTranspSaltData
  use SnowPhysData        , only : DestructSnowPhysData
  use HydroThermData      , only : DestructHydroThermData
  use BalanceCheckDataType, only : DestructBalanceCheckData
  use EcosysWarmingMod     , only : destructSoilWarming
  use TracerIDMod         , only : CleanUpTracerIDs
  use NumericalAuxMod     , only : DestructNumericAux
  use PlantTraitTableMod  , only : DestructPlantTraitTable
  use MicBGCAPI           , only : MicAPI_cleanup  
  use RadiationDataMod    , only : DestroyRadiationData
  use MLDataDiagType      , only : DestructMLDataType
  implicit none

  call MicAPI_cleanup

  call DestructMicrobialData

  call DestructChemTranspData

  call DestructSOMData

  call DestructFertilizerData

  call DestructCanopyRad

  call DestructSurfLitter

  call DestructIrrigation

  call DestructSedimentData

  call DestructSoilWater

  call DestroyPlantPhaseWork
  call DestructPlantAPIData

  call DestructPlantMngmtData

  call  DestructPlantRates

  call DestructAquaChem

  call DestructSOMBGC

  call DestructFlagData

  ! Non-salt transport is also active in salt-enabled standalone simulations.
  call DestructTranspNoSalt
  if(salt_model)call DestructTranspSaltData

  call DestructPlantTraitTable()

  call DestructEcoSIMCtrlData

  call DestructCanopyData

  call DestructEcoSIMHistData

  call DestructEcosimBGCFluxData

  call DestructLandSurfData

  call DestructGridData

  call DestructPlantTraits

  call DestructRootData

  call DestructSnowData

  call DestructSoilHeatData

  call DestructSoilPhysData

  call DestructSurfSoilData

  call DestructErosion

  call DestructSnowPhysData

  call DestructWatSubData

  call DestructHydroThermData

  call DestructEcoSimSum

  call DestructBalanceCheckData

  call destructSoilWarming()

  call DestructNumericAux()
  
  call DestructMLDataType
  
  call DestroyRadiationData
  call CleanUpTracerIDs
  end subroutine DestructEcoSIM

end module EcoSIMDesctruct
