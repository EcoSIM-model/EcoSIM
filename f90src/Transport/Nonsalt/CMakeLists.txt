set(NONSALT_MOD_SOURCES
  TranspNoSaltDataMod.F90
  InitNoSaltTransportMod.F90  
  TranspNoSaltFastMod.F90
  TranspNoSaltSlowMod.F90
	TranspNoSaltMod.F90
)
include(add_ecosim_library)
add_ecosim_library(Nonsalt ${NONSALT_MOD_SOURCES})

target_include_directories(Nonsalt  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/  
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/
  ${CMAKE_BINARY_DIR}/f90src/APIs/
  ${CMAKE_BINARY_DIR}/f90src/Plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Box_Micmodel
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Layers_Micmodel)


set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};Nonsalt PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};Nonsalt)
install(TARGETS Nonsalt DESTINATION ${CMAKE_BINARY_DIR}/lib)  

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS Nonsalt DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
