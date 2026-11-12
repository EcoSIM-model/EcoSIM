set(CANOPYPHYS_SOURCES
 CanopyHydroMod.F90
)
include(add_ecosim_library)
add_ecosim_library(CanopyPhys ${CANOPYPHYS_SOURCES})

target_include_directories(CanopyPhys  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};CanopyPhys PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};CanopyPhys)
install(TARGETS CanopyPhys DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS CanopyPhys DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
