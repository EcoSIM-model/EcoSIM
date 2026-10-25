set(IOUTILS_SOURCES
    ReadManagementMod.F90
    ForcWriterMod.F90
  	ClimReadMod.F90
    PlantInfoMod.F90
    Read4SatePhenolMod.F90
    HistDataType.F90
    HistData/HistAllocateColumn.F90
    HistData/HistAllocatePlant.F90
    HistData/HistAllocateProfile.F90
    HistData/HistRegisterColumn.F90
    HistData/HistRegisterPlant.F90
    HistData/HistRegisterSoilBGC.F90
    HistData/HistRegisterSoilPhysRoot.F90
    HistData/HistUpdate.F90
    HistData/HistUpdateColumn.F90
    HistData/HistUpdateSoil.F90
    HistData/HistUpdatePlant.F90
    MicrobeInfoMod.F90
    HistFileMod.F90
    restUtilMod.F90
    RestartMemoryMod.F90
    RestartMod.F90
  	readimod.F90
	  readsmod.F90
    bhistMod.F90
)
include(add_ecosim_library)
add_ecosim_library(ioutils ${IOUTILS_SOURCES})
#add_library(ioutils OBJECT ${IOUTILS_SOURCES})
#install(TARGETS ioutils OBJECTS DESTINATION lib)

# HistFileMod.F90 is large enough that -O3 (CMAKE_Fortran_FLAGS_RELEASE) makes
# gfortran emit a branch/relocation that Apple's integrated assembler cannot
# encode on arm64 ("addend too big for relocation" / "fixup value out of
# range"), aborting the Release build. -O1 compiles cleanly and this file is
# I/O bookkeeping, not a hot loop, so the perf cost is negligible. The
# trailing flag wins for gfortran, so this appended -O1 overrides the
# target/directory -O3 for this one source file only.
set_source_files_properties(HistFileMod.F90 PROPERTIES COMPILE_OPTIONS "-O1")

target_include_directories(ioutils  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/   
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/APIData/
  ${CMAKE_BINARY_DIR}/f90src/ModelDiags/
  ${CMAKE_BINARY_DIR}/f90src/Plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Layers_Micmodel/  
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SoilPhys/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/
  ${CMAKE_BINARY_DIR}/f90src/Disturbances/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_mods/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ioutils PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ioutils)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ioutils DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
