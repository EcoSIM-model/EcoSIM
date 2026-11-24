set(DEBUGTOOLS_MOD_SOURCES
  DebugToolMod.F90
)
include(add_ecosim_library)
add_ecosim_library(DebugTools ${DEBUGTOOLS_MOD_SOURCES})

target_include_directories(DebugTools  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};DebugTools PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};DebugTools)
install(TARGETS DebugTools DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS DebugTools DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
