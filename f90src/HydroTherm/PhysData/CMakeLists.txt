# Litter hydraulics now belongs to Constitutive; discard its old module artifact.
file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/litterhydraulicmod.mod")
# PhysPars is now compiled by modelpars; discard the previous module location.
file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/physpars.mod")

set(PHYSDATA_SOURCES
 WatsubDataMod.F90
 SurfaceAirExchangeMod.F90
 SurfaceVaporExchangeMod.F90
 AtmosCompositionMod.F90
 HydroThermData.F90
 SoilPhysParaMod.F90
 SoilFaceMappingMod.F90
)
include(add_ecosim_library)
add_ecosim_library(PhysData ${PHYSDATA_SOURCES})

target_include_directories(PhysData  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/  
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};PhysData PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};PhysData)

install(TARGETS PhysData DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS PhysData DESTINATION lib)  
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
