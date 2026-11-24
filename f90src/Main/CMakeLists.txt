set(ECOSIM_MAIN_SOURCES
  InitEcoSIM.F90
  EcoSIMDesctruct.F90
)
include(add_ecosim_library)
add_ecosim_library(ecosim_main ${ECOSIM_MAIN_SOURCES})

target_include_directories(ecosim_main  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/APIData/
  ${CMAKE_BINARY_DIR}/f90src/Plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/Disturbances/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SoilPhys
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Box_Micmodel
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Layers_Micmodel
  ${CMAKE_BINARY_DIR}/f90src/Modelforc/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_mods/
  ${CMAKE_BINARY_DIR}/f90src/Balances/
  ${CMAKE_BINARY_DIR}/f90src/Transport/Nonsalt
  ${CMAKE_BINARY_DIR}/f90src/Transport/Salt
  ${CMAKE_BINARY_DIR}/f90src/APIs/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ecosim_main PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ecosim_main)
install(TARGETS ecosim_main DESTINATION ${CMAKE_BINARY_DIR}/lib)  

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ecosim_main DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
