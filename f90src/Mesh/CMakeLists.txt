set(ECOSIM_MESH_SOURCES
  GridConsts.F90
  GridMod.F90
)

#include_directories(${CMAKE_BINARY_DIR}/f90src/Utils/)
#include_directories(${CMAKE_BINARY_DIR}/f90src/Minimath/)
#include_directories(${CMAKE_BINARY_DIR}/f90src/Modelconfig/)


include(add_ecosim_library)
add_ecosim_library(grid_mesh ${ECOSIM_MESH_SOURCES})
#add_library(mesh OBJECT ${ECOSIM_MESH_SOURCES})
#install(TARGETS mesh OBJECTS DESTINATION lib)

target_include_directories(grid_mesh PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};grid_mesh PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};grid_mesh)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS grid_mesh DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
