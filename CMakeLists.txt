# NOTE: this is the old CMake file for compiling with a UNIX Makefile
# Fortran submodule dependency scanning requires CMake 3.12.
cmake_minimum_required (VERSION 3.12)

# We declare the project here.
project (ecosim LANGUAGES NONE)

# Version numbers.
set (ECOSIM_MAJOR_VERSION 0)
set (ECOSIM_MINOR_VERSION 1)
set (ECOSIM_PATCH_VERSION 0)
set (ECOSIM_VERSION "${ECOSIM_MAJOR_VERSION}.${ECOSIM_MINOR_VERSION}.${ECOSIM_PATCH_VERSION}")

  # This tells CMake/Clang to target version 26.0 specifically
if (APPLE)
    # Execute a system command to get the current SDK version (e.g., 26.0 or 15.1)
    execute_process(
        COMMAND xcrun --show-sdk-version
        OUTPUT_VARIABLE CURRENT_SDK_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Only force the deployment target if it hasn't been set or is outdated
    if(NOT CMAKE_OSX_DEPLOYMENT_TARGET)
        set(CMAKE_OSX_DEPLOYMENT_TARGET "${CURRENT_SDK_VERSION}" CACHE STRING "OSX Deployment Target" FORCE)
    endif()

    message(STATUS "Setting Apple Deployment Target to SDK Version: ${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif()

set(TPL_INSTALL_PREFIX $ENV{AMANZI_TPLS_DIR})

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Options for building EcoSIM. These come from the xSDK compliance rules.
#ioption(USE_XSDK_DEFAULTS "Set to use xSDK defaults for options [ON]." ON)
#option(CMAKE_INSTALL_PREFIX "Sets installation prefix [/usr/local].")
#option(XSDK_ENABLE_DEBUG "Enables Debug mode builds [OFF]." OFF)
#option(BUILD_SHARED_LIBS "Builds shared libraries [OFF]." OFF)
#option(XSDK_WITH_NETCDF "Enables support for netcdf [OFF]." ON)
#option(TPL_NETCDF_LIBRARIES "List of absolute paths to netcdf link libraries [].")
#option(TPL_NETCDF_INCLUDE_DIRS "List of absolute paths to netcdf include directories [].")

# LANGUAGES NONE may leave the install prefix unset until languages are enabled
# below. Also replace /usr/local cached by older configurations, so repeated
# build_EcoSIM.sh runs install locally. Preserve other custom install prefixes.
if (NOT DEFINED CMAKE_INSTALL_PREFIX OR CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT
    OR CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/local" CACHE PATH "EcoSIM installation directory" FORCE)
endif()

## JDM: Not sure why we need this at all - is it used in running tests?
# set the number of cores to something reasonable so we don't over
# subscribe a laptop or exceed hpc login node restrictions. could
# switch to something more sophisticated like polymec if necessary.
set(NUMBER_OF_CORES 4)
math(EXPR NUM_BUILD_THREADS "${NUMBER_OF_CORES} + 1")

# This includes some custom macros - but does not call them
include(set_up_platform)
include(set_up_compilers)

# Make sure compilers are set. This must be done before enabling languages.
if (NOT CMAKE_C_COMPILER)
  if (NOT $ENV{CC} STREQUAL "")
    set(CMAKE_C_COMPILER $ENV{CC})
  else()
    set(CMAKE_C_COMPILER cc)
  endif()
endif()
if (NOT CMAKE_C_FLAGS)
  set(CMAKE_C_FLAGS $ENV{CFLAGS})
endif()
if (NOT CMAKE_CXX_COMPILER)
  if (NOT $ENV{CXX} STREQUAL "")
    set(CMAKE_CXX_COMPILER $ENV{CXX})
  else()
    set(CMAKE_CXX_COMPILER c++)
  endif()
endif()
if (NOT CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS $ENV{CXX_FLAGS})
endif()
if (NOT CMAKE_Fortran_COMPILER)
  if (NOT $ENV{FC} STREQUAL "")
    set(CMAKE_Fortran_COMPILER $ENV{FC})
  else()
    set(CMAKE_Fortran_COMPILER gfortran)
  endif()
endif()
if (NOT CMAKE_Fortran_FLAGS)
  set(CMAKE_Fortran_FLAGS $ENV{FCFLAGS})
endif()

# Calling the platform setup macro
set_up_platform()

# Enabling the languages for the project
enable_language(C)
enable_language(CXX)
enable_language(Fortran)

message("-- C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})")
message("-- CXX compiler is ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID})")
message("-- Fortran compiler is ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")

# Calling the custom macro to setup the compilers - not clear if this overrides flags we already set
set_up_compilers()
include(build_timing)

if (BGC)
  set (ECOSIM_BGC 1)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DECOSIM_BGC")
else()
  set (ECOSIM_BGC 0)
endif()

if (ECOSIM)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DECOSIM")
endif()

# JDM: Isn't this covered by the specialized macro we already called?
# Figure out the system type.
set(ECOSIM_HAVE_BOOL 1) # All reasonable C99 compilers have this now.
if (APPLE EQUAL 1)
  set(SYS_FLAGS "-DAPPLE=1")
  set(DYLIB_SUFFIX "dylib")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
else ()
  if (LINUX EQUAL 1)
    set(SYS_FLAGS "-DLINUX=1")
    set(DYLIB_SUFFIX "so")
  else()
    if (WIN32 EQUAL 1)
      set(ECOSIM_HAVE_BOOL 0) # MS doesn't have reasonable C compilers.
      set(SYS_FLAGS "-DWINDOWS=1")
      set(DYLIB_SUFFIX "dll")
    endif()
  endif ()
endif ()

# JDM: We are overriding things that are already provided, probably should wrap this?
# Here we make sure CMake-installed binaries use the correct runpath, and
# that the path is not stripped during installation.
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CMAKE_BINARY_DIR}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Basic libraries to be linked in.
# added the standard mathematical library for C and C++
#set(ECOSIM_LIBRARIES "")
#set(ECOSIM_TPLS m)

# JDM: We should provide these for Amanzi-ATS builds, but in general, yes - it's a pain but you should find them
if (${NEED_LAPACK})
  # NEED_LAPACK is set by set_up_platform()
  include(FindBLAS)
  include(FindLAPACK)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  if (${LAPACK_LIBRARY_DIR})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${LAPACK_LIBRARY_DIR}")
  endif()
  if (${BLAS_LIBRARY_DIR})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${BLAS_LIBRARY_DIR}")
  endif()
  set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};${LAPACK_LIBRARIES};${BLAS_LIBRARIES})
endif()

set(ECOSIM_NEED_PFUNIT FALSE)
set(ECOSIM_NEED_PETSC FALSE)
set(ECOSIM_INCLUDE_DIRS "" CACHE STRING "Include directory of EcoSIM" FORCE)

# JDM: Do we need this?
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Templates/ecosim.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/ecosim.cmake"
  @ONLY
)

# Other third-party libraries.
# adding in ATS tpl support:

if(ATS_ECOSIM)
  message(STATUS "AG: >>>>>> ATS_ECOSIM is set")
  #set(ATS_ECOSIM TRUE)
else()
  message(STATUS "AG: >>>>>> ATS_ECOSIM NOT set")
endif()

if (ATS_ECOSIM)
    message(STATUS "JDM: >>>>>> Using TPLs provided by Amanzi-ATS build")
    message(STATUS "JDM: >>>>>> TPL include and libraries definitions ... ")
    message(STATUS "JDM: >>>>>> HDF5")
    message(STATUS "JDM: >>>>>>   TPLS_HDF5_INCLUDE_DIRS = ${TPLS_HDF5_INCLUDE_DIRS}")
    message(STATUS "JDM: >>>>>>   TPLS_HDF5_LIBRARIES    = ${TPLS_HDF5_LIBRARIES}")
    message(STATUS "JDM: >>>>>> NETCDF")
    message(STATUS "JDM: >>>>>>   TPLS_NETCDF_INCLUDE_DIRS = ${TPLS_NETCDF_INCLUDE_DIRS}")
    message(STATUS "JDM: >>>>>>   TPLS_NETCDF_LIBRARIES    = ${TPLS_NETCDF_LIBRARIES}")
    message(STATUS "JDM: >>>>>> NETCDF_FORTRAN")
    message(STATUS "JDM: >>>>>>   TPLS_NETCDF_FORTRAN_INCLUDE_DIRS = ${TPLS_NETCDF_FORTRAN_INCLUDE_DIRS}")
    message(STATUS "JDM: >>>>>>   TPLS_NETCDF_FORTRAN_LIBRARIES    = ${TPLS_NETCDF_FORTRAN_LIBRARIES}")
    message(STATUS "JDM: >>>>>>>>>>>")
 
    set(ECOSIM_HAVE_NETCDF TRUE)
    #Needed to share with ATS (probably an easier way to do this)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${NETCDF_CFLAGS} -fPIC")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${NETCDF_FFLAGS} -fPIC -g -O0 -fbacktrace -fbounds-check")

    #What were these even for? are they neccesary?
    #include_directories("${TPL_INSTALL_PREFIX}/include")
    include_directories(${TPLS_NETCDF_FORTRAN_INCLUDE_DIRS} ${TPLS_NETCDF_INCLUDE_DIRS} ${TPLS_HDF5_INCLUDE_DIRS})
    #link_directories("/global2/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/lib")

    # Check for math library (libm)
    find_library(M_LIBRARY m)
    if (M_LIBRARY)
        set(ECOSIM_TPLS "${ECOSIM_TPLS} ${M_LIBRARY}")
    endif()

    # Check for zlib library
    find_library(ZLIB_LIBRARY z)
    if (ZLIB_LIBRARY)
        set(ECOSIM_TPLS "${ECOSIM_TPLS} ${ZLIB_LIBRARY}")
    endif()

    # Check for libzip library
    find_library(ZIP_LIBRARY zip)
    if (ZIP_LIBRARY)
        set(ECOSIM_TPLS "${ECOSIM_TPLS} ${ZIP_LIBRARY}")
    endif()

    # Check for libxml2 library
    find_library(XML2_LIBRARY xml2)
    if (XML2_LIBRARY)
        set(ECOSIM_TPLS "${ECOSIM_TPLS} ${XML2_LIBRARY}")
    endif()

    # Check for libcurl library
    find_library(CURL_LIBRARY curl)
    if (CURL_LIBRARY)
        set(ECOSIM_TPLS "${ECOSIM_TPLS} ${CURL_LIBRARY}")
    endif()

    set(ECOSIM_TPLS "${TPLS_NETCDF_FORTRAN_LIBRARIES}")
    
else()
    add_subdirectory(3rd-partylibs)
    
    #string(REPLACE "-ffpe-trap=invalid" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()

#string(REPLACE "-ffpe-trap=invalid" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
#message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
#message(STATUS "CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
#message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")

if (NOT ECOSIM_HAVE_NETCDF)
  list(APPEND ECOSIM_TPLS ${TPL_NETCDF_LIBRARIES})
  list(APPEND ECOSIM_INCLUDE_DIRS ${TPL_NETCDF_INCLUDE_DIRS})
endif()
# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.

message(STATUS "linking directories: ")
message(STATUS "linking include and lib in ${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/lib")

#list(APPEND ECOSIM_INCLUDE_DIRS "${PROJECT_BINARY_DIR}/include")
#include_directories(${ECOSIM_INCLUDE_DIRS})

# set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/local)

# JDM: I think this could some much earlier in the file.
enable_testing()

# Source code itself.
include_directories("${PROJECT_SOURCE_DIR}")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
    message(STATUS "CMAKE_BUILD_TYPE not specified, defaulting to Debug")
endif()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

add_subdirectory(f90src)
add_subdirectory(drivers)
