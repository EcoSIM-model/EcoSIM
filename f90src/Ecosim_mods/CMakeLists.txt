set(ECOSIM_MOD_SOURCES
    # Allocate transport workspaces before either standalone or coupled stepping.
    ../Transport/Salt/TranspSaltDataMod.F90
    InitAllocMod.F90
	StarteMod.F90
	StartsMod.F90
	StartqMod.F90  
)
include(add_ecosim_library)
add_ecosim_library(ecosim_mods ${ECOSIM_MOD_SOURCES})

target_include_directories(ecosim_mods  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/   
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/ModelDiags/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SoilPhys
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SurfPhys
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData
  ${CMAKE_BINARY_DIR}/f90src/Disturbances
  ${CMAKE_BINARY_DIR}/f90src/APIData/
  ${CMAKE_BINARY_DIR}/f90src/APIs/
  ${CMAKE_BINARY_DIR}/f90src/Plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/Geochem/Box_chem
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Box_Micmodel
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Layers_Micmodel)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ecosim_mods PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ecosim_mods)
install(TARGETS ecosim_mods DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ecosim_mods DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
