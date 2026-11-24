set(MINIMATH_SOURCES
  LinearAlgebraMod.F90
  MiniMathMod.F90
  MiniFuncMod.F90
)

include(add_ecosim_library)
add_ecosim_library(minimath ${MINIMATH_SOURCES})

target_include_directories(minimath  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/)



set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};minimath PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};minimath)
install(TARGETS minimath DESTINATION ${CMAKE_BINARY_DIR}/lib)  

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS minimath DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
