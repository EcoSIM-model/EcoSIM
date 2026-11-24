file(GLOB BOXGEOCHEM_SOURCES "*.F90")
#set(BOXGEOCHEM_SOURCES
#  SoluteChemDataType.F90
#  SaltChemEquilibriaMod.F90
#  ChemEquilibriaMod.F90  
#  InitSoluteMod.F90
#  SoluteParMod.F90
#)
include(add_ecosim_library)
add_ecosim_library(Box_chem ${BOXGEOCHEM_SOURCES})

target_include_directories(Box_chem  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/   
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/
  ${CMAKE_BINARY_DIR}/f90src/Plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/)


set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};Box_chem PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};Box_chem)

install(TARGETS Box_chem DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS Box_chem DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
