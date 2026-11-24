set(SNOWPHYS_SOURCES
 SnowPhysData.F90
 SnowPhysMod.F90
 SnowBalanceMod.F90
 SnowTransportMod.F90
)
include(add_ecosim_library)
add_ecosim_library(SnowPhys ${SNOWPHYS_SOURCES})

target_include_directories(SnowPhys  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Layers_Micmodel    
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};SnowPhys PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};SnowPhys)
install(TARGETS SnowPhys DESTINATION ${CMAKE_BINARY_DIR}/lib)


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS SnowPhys DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
