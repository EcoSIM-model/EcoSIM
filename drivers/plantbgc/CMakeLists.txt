set(PLANT_DRIVER_SOURCES
  PlantTestMod.F90
)

include(add_ecosim_library)
add_ecosim_library(plant_driver ${PLANT_DRIVER_SOURCES})

target_include_directories(plant_driver PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/)


set(ECOSIM_LIBRARIES plant_driver;${ECOSIM_LIBRARIES} PARENT_SCOPE)
set(ECOSIM_LIBRARIES plant_driver;${ECOSIM_LIBRARIES})

include(add_ecosim_executable)
add_ecosim_executable(plant.x plantdriver.F90)

target_include_directories(plant.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/)


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS plant.x DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
