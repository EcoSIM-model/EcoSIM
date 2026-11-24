set(BALANCES_MOD_SOURCES
  RedistDataMod.F90
  LateralTranspMod.F90
  SoilLayerDynMod.F90
#  RunoffBalMod.F90
  ErosionBalMod.F90
	RedistMod.F90
  TillageMixMod.F90
)
include(add_ecosim_library)
add_ecosim_library(balances ${BALANCES_MOD_SOURCES})

target_include_directories(balances  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/   
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/  
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SurfPhys/
  ${CMAKE_BINARY_DIR}/f90src/Disturbances/  
  ${CMAKE_BINARY_DIR}/f90src/ModelDiags/    
  ${CMAKE_BINARY_DIR}/f90src/Plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Box_Micmodel  
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Layers_Micmodel)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};balances PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};balances)
install(TARGETS balances DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS balances DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
