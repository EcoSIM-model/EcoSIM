set(DISTURB_MOD_SOURCES
	ErosionMod.F90
  FireMod.F90
  EcosysWarmingMod.F90
  PlantDisturbMod.F90
  SoilDisturbMod.F90
  FertilizerMod.F90
)
include(add_ecosim_library)
add_ecosim_library(disturb_mods ${DISTURB_MOD_SOURCES})

target_include_directories(disturb_mods  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  )

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};disturb_mods PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};disturb_mods)
install(TARGETS disturb_mods DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS disturb_mods DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
