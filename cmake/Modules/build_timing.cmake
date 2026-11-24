# Time executed compile/link/archive commands for both scripted and direct builds.
option(ECOSIM_BUILD_TIMING "Record wall time and exit status of build commands" ON)
set(ECOSIM_BUILD_TIMING_LOG "${CMAKE_BINARY_DIR}/build-timings.log" CACHE FILEPATH
    "Append-only TSV build timing log")

if(ECOSIM_BUILD_TIMING)
  if(NOT UNIX OR NOT CMAKE_GENERATOR MATCHES "Makefiles|Ninja")
    message(WARNING "Per-command build timing requires Unix Makefiles or Ninja on Unix")
  else()
    find_program(ECOSIM_TIMING_PYTHON NAMES python3)
    if(NOT ECOSIM_TIMING_PYTHON)
      message(FATAL_ERROR "Build timing requires python3; alternatively set ECOSIM_BUILD_TIMING=OFF")
    endif()
    get_filename_component(_ecosim_timer "${CMAKE_CURRENT_LIST_DIR}/../build_timer.py" ABSOLUTE)
    foreach(_lang C CXX Fortran)
      # Retain existing launchers such as ccache inside the timer.
      set(_ecosim_compile_launcher
          "${ECOSIM_TIMING_PYTHON};${_ecosim_timer};run;--log;${ECOSIM_BUILD_TIMING_LOG};--phase;compile;--")
      if(CMAKE_${_lang}_COMPILER_LAUNCHER)
        list(APPEND _ecosim_compile_launcher ${CMAKE_${_lang}_COMPILER_LAUNCHER})
      endif()
      set(CMAKE_${_lang}_COMPILER_LAUNCHER "${_ecosim_compile_launcher}")
    endforeach()
    # RULE_LAUNCH_LINK also covers ar/ranlib and supports the project's CMake 3.5
    # minimum; the newer language-specific linker launchers do not cover archives.
    get_property(_ecosim_previous_link DIRECTORY PROPERTY RULE_LAUNCH_LINK)
    set_property(DIRECTORY PROPERTY RULE_LAUNCH_LINK
      "\"${ECOSIM_TIMING_PYTHON}\" \"${_ecosim_timer}\" run --log \"${ECOSIM_BUILD_TIMING_LOG}\" --phase link -- ${_ecosim_previous_link}")
    message(STATUS "Build timing log: ${ECOSIM_BUILD_TIMING_LOG}")
  endif()
endif()
