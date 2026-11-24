set(APIS_SOURCES
  PlantColumnTransferMod.F90
  PlantStateTransferMod.F90
  PlantCanopyTransferMod.F90
  PlantRootTransferMod.F90
  PlantRootAxisTransferMod.F90
  PlantTraitsTransferMod.F90
	PlantAPI.F90
  PlantCanAPI.F90
  PlantMod.F90
  PlantPhaseWorkMod.F90
  PlantHydraulicMod.F90
  MicBGCAPI.F90
  PlantAPI4Uptake.F90  
  GeochemAPI.F90
  SurfPhysAPI.F90  
)

include(add_ecosim_library)
add_ecosim_library(APIs ${APIS_SOURCES})

target_include_directories(APIs PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/    
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/ModelDiags/
  ${CMAKE_BINARY_DIR}/f90src/APIData/
  ${CMAKE_BINARY_DIR}/f90src/Plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Box_Micmodel
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Layers_Micmodel
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/
  ${CMAKE_BINARY_DIR}/f90src/Disturbances/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_mods/
  ${CMAKE_BINARY_DIR}/f90src/Balances/
  ${CMAKE_BINARY_DIR}/f90src/Transport/Nonsalt
  ${CMAKE_BINARY_DIR}/f90src/Transport/Salt
  ${CMAKE_BINARY_DIR}/f90src/Geochem/Box_chem
  ${CMAKE_BINARY_DIR}/f90src/Geochem/Layers_chem)


set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};APIs PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};APIs)

install(TARGETS APIs DESTINATION ${CMAKE_BINARY_DIR}/lib)
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS APIs DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
