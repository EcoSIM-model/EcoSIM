set(ECOSIM_UTIL_SOURCES
  data_kind_mod.F90
  data_const_mod.F90
  StringToolsMod.F90
  ncdio_pio.F90
  abortutils.F90
  UnitMod.F90
  TestMod.F90
#  clock.c
#  getfilename.c
  timings.F90
  fileUtil.F90
  EcoSimConst.F90
  ecosim_log_mod.F90
  ecosim_time_mod.F90
  ModelStatusType.F90

)

#configure_file(shr_assert.h shr_assert.h COPYONLY)
configure_file(dtypes.h dtypes.h COPYONLY)

include(add_ecosim_library)
add_ecosim_library(util ${ECOSIM_UTIL_SOURCES})
#add_library(util OBJECT ${ECOSIM_UTIL_SOURCES})    
#install(TARGETS util OBJECTS DESTINATION lib)

target_include_directories(util PUBLIC ${CMAKE_BINARY_DIR}/f90src/Utils)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};util PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};util)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS util DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
