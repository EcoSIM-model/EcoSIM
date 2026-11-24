file(GLOB MODELPARS_SOURCES "*.F90")
#set(MODELPARS_SOURCES
#  NitroPars.F90
#  PlantBGCPars.F90
#  TracerPropMod.F90
#  MicBGCPars.F90
#  EcoSiMParDataMod.F90
#  ChemTracerParsMod.F90
#)
include(add_ecosim_library)
add_ecosim_library(modelpars ${MODELPARS_SOURCES})
#add_library(modelpars OBJECT ${MODELPARS_SOURCES})
#install(TARGETS modelpars OBJECTS DESTINATION lib)

target_include_directories(modelpars  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/  
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/)


set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};modelpars PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};modelpars)

message("install install install")
file(GLOB files_with_wildcard "${CMAKE_BINARY_DIR}/f90src/Modelpars/libmodelpars.*")
if (files_with_wildcard)
message("lib exists")
foreach(file ${files_with_wildcard})
  file(COPY "${file}" DESTINATION "${CMAKE_BINARY_DIR}/lib")
endforeach()

else()
message("no file CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
endif()
install(TARGETS modelpars DESTINATION lib)
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  message("install install install")
  install(TARGETS modelpars DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
