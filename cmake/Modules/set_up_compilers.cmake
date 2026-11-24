# This macro sets compiler specific flags
macro(set_up_compilers)
  #
  # General C compiler flags.
  #
  message("HOSTNAME ${HOSTNAME}")
  message("CMAKE_C_COMPILER_ID ${CMAKE_C_COMPILER_ID}")
  if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wpedantic -Wextra")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror-implicit-function-declaration")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function -Wno-cast-qual")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-but-set-variable -Wno-overlength-strings")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-discarded-qualifiers")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-conversion -Wno-maybe-uninitialized")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${NETCDF_CFLAGS} -ffpe-trap=invalid,zero,overflow,underflow")

    if (BUILD_SHARED_LIBS)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fPIC")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffpe-trap=invalid,zero,overflow,underflow")

    if (LINUX EQUAL 1)
      # Counter some of GCC's more recent stinginess on Linux.
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=200809L")
    endif()

  elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror-implicit-function-declaration ")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-builtin")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter ")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function")

  elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
    if (CMAKE_BUILD_TYPE MATCHES "Debug")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -g")
    else()
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")
    endif()
  endif()

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

  #
  # Fortran compiler flags.
  #
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_BUILD_TYPE MATCHES "Debug")
#     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W -Wall -std=gnu -pedantic -finit-local-zero -Wall -cpp -g -fbacktrace -fmax-stack-var-size=524288 ")
     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W -Wall -fbounds-check -std=gnu -pedantic -fcheck=all -cpp -g -fbacktrace -fmax-stack-var-size=524288 ")
    else()
#     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W -Wall -std=gnu -pedantic -finit-local-zero -cpp -O2 -fmax-stack-var-size=524288")
     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=gnu -fbounds-check -pedantic -finit-local-zero -cpp -O2 -fmax-stack-var-size=524288")

    endif()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
# --- Jing Tao ---
# UNDERFLOW removed from the trap list. Underflow is a BENIGN IEEE condition -- a
# result too small for the normal range decays gracefully into the subnormals -- not an
# error, and trapping it aborts arithmetic that is perfectly well defined. It is the
# proximate cause of 887 of the 1,129 EcoSIM_BioCON R3 crashes (79%), all in dying
# stands. Mechanism, reproduced in isolation with these exact flags: PlantMathFuncMod::
# ExchFluxLimiter caps XFRE at fromState*0.9999, so the caller in PlantBranchMod::
# WithdrawNutBranchLeaves computes fromState - XFRE = fromState*1e-4. Once a decaying
# pool falls below ~1e-304 that quotient lands among the subnormals, where reduced
# precision makes it INEXACT; tiny + inexact is exactly the IEEE underflow condition,
# so the trap fires and the run dies at a SUBTRACTION. Verified: the same program with
# 'underflow' removed returns 1e-309 and 1e-311 correctly and exits 0.
# invalid / zero / overflow are KEPT -- those are genuine errors and must still abort.
# Removing a trap cannot change any arithmetic RESULT, only whether execution
# continues, so a run that completed under the old flags is bit-for-bit unaffected.
# That is a claim to TEST, not assume -- see the inertness gate in
# use_cases/EcoSIM_BioCON/memory/phase_results/20260820a_*/RELAUNCH_FAILED_CASES.md
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-variable -ffpe-trap=invalid,zero,overflow")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-parameter -Wno-unused-function")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wuninitialized -Wno-unused-dummy-argument")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCPRINTEL")

    if (HOSTNAME MATCHES "cori")
      message("Using set up for cori")
      if (CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -debug -r8 -i4 -align dcommons -auto-scalar -fimf-arch-consistency=true -fpe-all=0")
      else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -mp1 -r8 -i4 -align dcommons -auto-scalar -fimf-arch-consistency=true")
      endif()
    elseif (HOSTNAME MATCHES "biomes")
      message("Using set up for biomes")
      if (CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -debug -r8 -i4 -align dcommons -auto-scalar -fimf-arch-consistency=true -fpe-all=0")
      else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -show -list -O2 -mp1 -r8 -i4 -align dcommons -auto-scalar -fimf-arch-consistency=true")
      endif()      
    elseif (HOSTNAME MATCHES "[scs]")
      message("Using set up for lawrencium")
      if (CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -debug -r8 -i4 -align dcommons -auto-scalar -fimf-arch-consistency=true -fpe-all=0")
      else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -mp1 -r8 -i4 -align dcommons -auto-scalar -fimf-arch-consistency=true")
      endif()
	  else()
      if (CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -debug -r8 -i4 -align dcommons -auto-scalar -fpe-all=0 -check all")
      else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O -show -check -mp1 -r8 -i4 -align dcommons -auto-scalar -fimf-arch-consistency=true")
      endif()
    endif()
  endif()

endmacro()
