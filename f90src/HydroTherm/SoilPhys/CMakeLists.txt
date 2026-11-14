# PhysData owns this module now; discard an orphan from older configured builds.
file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/watsubdatamod.mod")

set(SOILPHYS_SOURCES
  SoilHydroParaMod.F90
  WatsubMod.F90
)

include(add_ecosim_library)
add_ecosim_library(SoilPhys ${SOILPHYS_SOURCES})

target_include_directories(SoilPhys  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/ModelDiags/
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/Disturbances/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Layers_Micmodel  
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SurfPhys/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};SoilPhys PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};SoilPhys)
install(TARGETS SoilPhys DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS SoilPhys DESTINATION lib)  
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
