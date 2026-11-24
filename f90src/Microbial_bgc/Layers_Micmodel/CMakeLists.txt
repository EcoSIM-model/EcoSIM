set(LAYERS_SOURCES
  InitSOMBGCMod.F90
	SoilBGCNLayMod.F90
  SoilPARAttenuationMod.F90
)
include(add_ecosim_library)
add_ecosim_library(Layers_Micmodel ${LAYERS_SOURCES})

target_include_directories(Layers_Micmodel  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/     
  ${CMAKE_BINARY_DIR}/f90src/Mesh/  
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/ModelDiags/   
  ${CMAKE_BINARY_DIR}/f90src/IOutils/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Disturbances/
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Box_Micmodel)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};Layers_Micmodel PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};Layers_Micmodel)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS Layers_Micmodel DESTINATION ${CMAKE_BINARY_DIR}/lib)  
  install(TARGETS Layers_Micmodel DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
