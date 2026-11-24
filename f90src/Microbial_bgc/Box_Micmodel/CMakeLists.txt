set(BOX_MODEL_SOURCES
  MicrobMathFuncMod.F90
  AerobicBacteriaMod.F90
  CyanoBacteriaMod.F90
  FungiMod.F90
  FermenterMod.F90
  DenitrifierMod.F90
  MethanogenMod.F90
  MethanotrophMod.F90
  NitrifierMod.F90
  MicAutoCplxFGMod.F90
  MicBGCFGMod.F90
  MicFluxTypeMod.F90
  MicForcTypeMod.F90
  MicStateTraitTypeMod.F90
  MicrobeDiagTypes.F90
)

include(add_ecosim_library)
add_ecosim_library(Box_Micmodel ${BOX_MODEL_SOURCES})

target_include_directories(Box_Micmodel  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};Box_Micmodel PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};Box_Micmodel)
install(TARGETS Box_Micmodel DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS Box_Micmodel DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
