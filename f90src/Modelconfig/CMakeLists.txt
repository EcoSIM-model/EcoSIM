set(ECOSIM_CONFG_SOURCES
  EcoSIMSolverPar.F90
  EcoSIMCtrlMod.F90
  EcoSIMConfig.F90
  ElmIDMod.F90
  TracerIDMod.F90
)
include(add_ecosim_library)
add_ecosim_library(modelconfig ${ECOSIM_CONFG_SOURCES})
#add_library(modelconfig OBJECT ${ECOSIM_CONFG_SOURCES})
#install(TARGETS modelconfig OBJECTS DESTINATION lib)  

target_include_directories(modelconfig  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  PUBLIC
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};modelconfig PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};modelconfig)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS modelconfig DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
