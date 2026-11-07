set(TOOLS_SOURCES
	ClimTransformer.F90
  SoilManagementReader.F90
  PlantManagementReader.F90
  ClimReader.F90
  etimerTest.F90
  restartTest.F90
  NamelistTest.F90
  HFileTest.F90
  SoilWarmReadTest.F90
  GroundResistanceTest.F90
  ThinSnowEnergyTest.F90
  SaltTransportDataTest.F90
	EcoATSTest.F90
)
include(add_ecosim_library)
include(add_ecosim_executable)
add_ecosim_executable(ClimTransformer.x ClimTransformer.F90)
add_ecosim_executable(SoilManagementReader.x SoilManagementReader.F90)
add_ecosim_executable(PlantManagementReader.x PlantManagementReader.F90)
add_ecosim_executable(ClimReader.x ClimReader.F90)
add_ecosim_executable(etimerTest.x etimerTest.F90)
add_ecosim_executable(NamelistTest.x NamelistTest.F90)
add_ecosim_executable(HFileTest.x HFileTest.F90)
add_ecosim_executable(restartTest.x restartTest.F90)
#add_ecosim_executable(EcoATSTest.x EcoATSTest.F90)
add_ecosim_executable(SoilWarmReadTest.x SoilWarmReadTest.F90)

#target_include_directories(EcoATSTest.x PRIVATE
#  ${CMAKE_BINARY_DIR}/f90src/ATSUtils/
#  )

target_include_directories(ClimTransformer.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/)

target_include_directories(SoilManagementReader.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/)

target_include_directories(PlantManagementReader.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/)

target_include_directories(ClimReader.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/)

target_include_directories(etimerTest.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/)

target_include_directories(NamelistTest.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/)

target_include_directories(HFileTest.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/)

target_include_directories(restartTest.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/)

  target_include_directories(SoilWarmReadTest.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/Disturbances/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/)


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ClimTransformer.x DESTINATION bin)
  install(TARGETS SoilManagementReader.x DESTINATION bin)
  install(TARGETS PlantManagementReader.x DESTINATION bin)
  install(TARGETS ClimReader.x DESTINATION bin)
  install(TARGETS etimerTest.x DESTINATION bin)
  install(TARGETS NamelistTest.x DESTINATION bin)
  install(TARGETS HFileTest.x DESTINATION bin)
  install(TARGETS restartTest.x DESTINATION bin)
  #install(TARGETS EcoATSTest.x DESTINATION bin)
  install(TARGETS SoilWarmReadTest.x DESTINATION bin)

  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/)
endif()

add_ecosim_executable(GroundResistanceTest.x GroundResistanceTest.F90)
target_include_directories(GroundResistanceTest.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SurfPhys/)

add_ecosim_executable(ThinSnowEnergyTest.x ThinSnowEnergyTest.F90)
target_include_directories(ThinSnowEnergyTest.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys/)

add_ecosim_executable(SaltTransportDataTest.x SaltTransportDataTest.F90)
target_include_directories(SaltTransportDataTest.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_mods/)
