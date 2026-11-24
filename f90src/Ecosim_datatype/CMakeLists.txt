set(ECOSIM_DATATYPE_SOURCES
  GridDataType.F90
  MicrobialDataType.F90
  SOMDataType.F90
  SoilPhysDataType.F90
  FlagDataType.F90
  ChemTranspDataType.F90
  PlantTraitDataType.F90
  FertilizerDataType.F90
  CanopyRadDataType.F90
  EcoSIMCtrlDataType.F90
  SoilHeatDataType.F90
  SoilWaterDataType.F90
  SoilBGCDataType.F90
  PlantDataRateType.F90
  ClimForcDataType.F90
  LandSurfDataType.F90
  SnowDataType.F90
  PlantMgmtDataType.F90
  SurfLitterDataType.F90
  SurfSoilDataType.F90
  CanopyDataType.F90
  RootDataType.F90
  PlantTraitTableMod.F90
  NumericalAuxMod.F90
  EcoSimSumDataType.F90
  EcosimBGCFluxType.F90
  AqueChemDatatype.F90
  EcoSIMHistMod.F90
  SoilPropertyDataType.F90
  IrrigationDataType.F90
  SedimentDataType.F90
  BalanceCheckDataType.F90
  MLDataDiagType.F90
)

#include_directories(${CMAKE_BINARY_DIR}/f90src/Utils/)
#include_directories(${CMAKE_BINARY_DIR}/f90src/Modelconfig/)
#include_directories(${CMAKE_BINARY_DIR}/f90src/Mesh/)
include(add_ecosim_library)
add_ecosim_library(ecosim_datatype ${ECOSIM_DATATYPE_SOURCES})
#add_library(ecosim_datatype OBJECT ${ECOSIM_DATATYPE_SOURCES})

target_include_directories(ecosim_datatype PRIVATE 
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/)

target_include_directories(ecosim_datatype PUBLIC ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ecosim_datatype PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ecosim_datatype)

install(TARGETS ecosim_datatype DESTINATION ${CMAKE_BINARY_DIR}/lib)
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")  
  install(TARGETS ecosim_datatype DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
