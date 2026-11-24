set(MODELFORC_MOD_SOURCES
  DepositionMod.F90
  YearMod.F90
	DayMod.F90
	Hour1Mod.F90
	WthrMod.F90
)
include(add_ecosim_library)
add_ecosim_library(modelforc ${MODELFORC_MOD_SOURCES})

target_include_directories(modelforc  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/
  ${CMAKE_BINARY_DIR}/f90src/Disturbances/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ModelDiags/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/
  ${CMAKE_BINARY_DIR}/f90src/ATSUtils/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/CanopyPhys
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SoilPhys
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData
  ${CMAKE_BINARY_DIR}/f90src/APIData/
  ${CMAKE_BINARY_DIR}/f90src/APIs/
  ${CMAKE_BINARY_DIR}/f90src/Plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/Geochem/Box_chem
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Box_Micmodel
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Layers_Micmodel
  )

target_link_libraries(modelforc APIData)


set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};modelforc PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};modelforc)
install(TARGETS modelforc DESTINATION ${CMAKE_BINARY_DIR}/lib)  

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS modelforc DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
