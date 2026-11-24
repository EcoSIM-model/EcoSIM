set(SURFPHYS_SOURCES
  SurfPhysData.F90
  SurfLitterPhysMod.F90
  SurfPhysMod.F90
)
include(add_ecosim_library)
add_ecosim_library(SurfPhys ${SURFPHYS_SOURCES})

target_include_directories(SurfPhys  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys/
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Layers_Micmodel)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};SurfPhys PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};SurfPhys)
install(TARGETS SurfPhys DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS SurfPhys DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
