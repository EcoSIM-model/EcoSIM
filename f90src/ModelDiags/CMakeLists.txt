set(MODELDIAGS_MOD_SOURCES
  BalancesMod.F90
  HydrologyDiagMod.F90
  SoilDiagsMod.F90
  MicrobialDiagMod.F90
  MLDataDiagMod.F90
)
include(add_ecosim_library)
add_ecosim_library(ModelDiags ${MODELDIAGS_MOD_SOURCES})

target_include_directories(ModelDiags  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ModelDiags PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ModelDiags)
install(TARGETS ModelDiags DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ModelDiags DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
