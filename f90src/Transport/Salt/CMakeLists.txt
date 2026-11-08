set(SALT_MOD_SOURCES
  IngridTranspMod.F90
	TranspSaltMod.F90
)
include(add_ecosim_library)
add_ecosim_library(Salt ${SALT_MOD_SOURCES})

target_include_directories(Salt  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/  
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_mods/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys/  
  ${CMAKE_BINARY_DIR}/f90src/APIs/
  ${CMAKE_BINARY_DIR}/f90src/Plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Box_model
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Layers_Micmodel)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};Salt PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};Salt)
install(TARGETS Salt DESTINATION ${CMAKE_BINARY_DIR}/lib)  

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS Salt DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
