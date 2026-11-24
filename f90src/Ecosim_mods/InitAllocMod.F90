module InitAllocMod
  ! Allocate model state and workspaces shared by standalone and coupled drivers.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod

implicit none
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  public :: InitAlloc
  contains
  subroutine InitAlloc()
  !
  !DESCRIPTION
  !allocate memeory for ecosim data

  use EcoSiMParDataMod    , only : micpar,pltpar
  use EcoSimSumDataType   , only : InitEcoSimSum
  use SnowDataType        , only : InitSnowData
  use FlagDataType        , only : InitFlagData
  use SurfSoilDataType    , only : InitSurfSoilData
  use SoilPhysDataType    , only : InitSoilPhysData
  use SoilHeatDatatype    , only : InitSoilHeatData
  use GridDataType        , only : InitGridData
  use RootDataType        , only : InitRootData
  use LandSurfDataType    , only : InitLandSurfData
  use EcosimBGCFluxType   , only : InitEcosimBGCFluxData
  use CanopyDataType      , only : InitCanopyData
  use EcoSIMCtrlDataType  , only : InitEcoSIMCtrlData
  use EcoSIMHistMod       , only : InitEcoSIMHistData
  use ClimForcDataType    , only : InitClimForcData
  use MicrobialDataType   , only : InitMicrobialData
  use SOMDataType         , only : InitSOMData
  use ChemTranspDataType  , only : InitChemTranspData
  use FertilizerDataType  , only : InitFertilizerData
  use CanopyRadDataType   , only : InitCanopyRad
  use InitSoluteMod       , only : InitSoluteProperty
  use TranspSaltDataMod   , only : InitTranspSaltData
  use AqueChemDatatype    , only : initaquachem
  use PlantDataRateType   , only : InitPlantRates
  use EcoSIMCtrlMod       , only : plant_model
  use PlantTraitDataType   , only : InitPlantTraits
  use SoilPropertyDataType, only : InitSoilProperty
  use SurfLitterDataType  , only : InitSurfLitter
  use SurfPhysData        , only : InitSurfPhysData
  use IrrigationDataType  , only : InitIrrigation
  use SoilBGCDataType     , only : InitSoilBGCData
  use SedimentDataType    , only : InitSedimentData
  use SoilWaterDataType   , only : InitSoilWater
  use PlantAPIData        , only : InitPlantAPIData
  use RadiationDataMod    , only : InitRadiationData
  use PlantMgmtDataType  , only : InitPlantMngmtData
  use InitSOMBGCMod       , only : InitSOMBGC
  use EcoSIMConfig        , only : jcplx1 => jcplxcm1
  use TracerIDMod         , only : InitTracerIDs
  use SnowPhysData        , only : InitSnowPhysData
  use HydroThermData      , only : InitHydroThermData
  use EcosysWarmingMod     , only : InitSoilWarming
  use NumericalAuxMod     , only : InitNumericAux
  use WatsubMod           , only : InitWatsub
  use BalanceCheckDataType, only : InitBalanceCheckData
  use PlantBGCPars          , only : InitPlantTraitTable
  use MLDataDiagType        , only : InitMLDataType
  use GridConsts  
  implicit none
  
! begin_execution

  call InitSOMBGC()

  call InitPlantMorphSize()

  if(plant_model)then    
    !and obtain plant NumGrowthStages and MaxNumRootAxes for output set up
    call InitPlantTraitTable(pltpar,NumGrowthStages,MaxNumRootAxes)
  endif  

  call InitGridData

  call InitTracerIDs(salt_model)

  call InitLandSurfData

  call InitEcoSIMCtrlData

  call InitCanopyData

  call InitCanopyRad

  call InitAquaChem

  call InitPlantMngmtData

  call InitPlantRates(micpar%NumOfPlantLitrCmplxs,pltpar%jroots)

  call InitSoilProperty

  call InitWatsub

  call InitSurfLitter(micpar%NumOfLitrCmplxs)

  call InitSedimentData

  call InitSoilWater

  call InitIrrigation

  call InitPlantAPIData()

  call InitMicrobialData

  call InitChemTranspData(salt_model)

  call InitSoilBGCData(pltpar%NumOfPlantLitrCmplxs)

  call InitSOMData(micpar%NumOfLitrCmplxs)

  call InitFertilizerData

  call InitPlantTraits(pltpar%NumOfPlantLitrCmplxs)

  call InitFlagData

  call InitEcoSimSum

  call InitRootData(2)

  call InitClimForcData

  call InitEcosimBGCFluxData

  call InitSnowData

  call InitEcoSIMHistData

  call InitSoilHeatData

  call InitSurfSoilData

  call InitHydroThermData

  call InitSnowPhysData

  call InitSoilPhysData

  call InitNumericAux()

  call InitRadiationData
  
  call InitMLDataType
  if(salt_model)then
    call InitSoluteProperty
    call InitTranspSaltData
  endif

  call InitSoilWarming
  
  call InitBalanceCheckData

  end subroutine InitAlloc
!------------------------------------------------------------------------------------------
  subroutine InitPlantMorphSize()
  use EcoSiMParDataMod, only : pltpar,micpar
  use GridConsts
  implicit none

  pltpar%JZ1                     = JZ
  pltpar%NumCanopyLayers1        = NumCanopyLayers
  pltpar%JP1                     = JP
  pltpar%NumOfLeafAzimuthSectors = NumOfLeafAzimuthSectors
  pltpar%NumOfSkyAzimuthSects1   = NumOfSkyAzimuthSects
  pltpar%NumLeafInclinationClasses1   = NumLeafInclinationClasses
  pltpar%MaxNodesPerBranch1      = MaxNodesPerBranch
  pltpar%iprotein                = micpar%iprotein
  pltpar%icarbhyro               = micpar%icarbhyro
  pltpar%icellulos               = micpar%icellulos
  pltpar%ilignin                 = micpar%ilignin
  pltpar%k_woody_comp            = micpar%k_woody_comp
  pltpar%k_fine_comp             = micpar%k_fine_comp
  
  !the following variable should be consistent with the soil bgc model
  pltpar%jcplx                = micpar%jcplx
  pltpar%NumOfPlantLitrCmplxs = micpar%NumOfPlantLitrCmplxs
  pltpar%jsken                = micpar%jsken
  pltpar%NumLitterGroups      = 5     !number of liter groups
  pltpar%MaxNumBranches       = 10    !number of branches
  pltpar%MaxNumRootAxes       = 10
  pltpar%NumGrowthStages      = 10
  pltpar%NumOfPlantMorphUnits = 7
  pltpar%NMaxRootSegs = 10    !
  NMaxRootSegs = pltpar%NMaxRootSegs
  MaxNumBranches       = pltpar%MaxNumBranches
  NumLitterGroups      = pltpar%NumLitterGroups
  MaxNumRootAxes       = pltpar%MaxNumRootAxes
  NumOfPlantMorphUnits = pltpar%NumOfPlantMorphUnits
  NumGrowthStages      = pltpar%NumGrowthStages
  end subroutine InitPlantMorphSize
end module InitAllocMod
