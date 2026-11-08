set(SOIL_HYDRAULIC_SOURCES
  SoilHydraulicModelMod.F90
  EcoSIMSoilHydraulicMod.F90
  SoilHydraulicSchemeMod.F90
)
include(add_ecosim_library)
add_ecosim_library(SoilHydraulics ${SOIL_HYDRAULIC_SOURCES})
target_include_directories(SoilHydraulics PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};SoilHydraulics PARENT_SCOPE)
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS SoilHydraulics DESTINATION lib)
endif()
