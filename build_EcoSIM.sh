
# ############################################################################ #
#
# EcoSIM build script
#
# ############################################################################ #

# A couple of helper functions:

status_message()
{
  #Printing status messages in a nice format
  local GREEN='32m'
  if [ "${no_color}" -eq "${TRUE}" ]; then
    echo "[`date`] $1"
  else
    echo -e "[`date`]\033[$GREEN $1\033[m"
  fi
}

############# EDIT THESE! ##################
parallel_jobs=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)
debug=0
mpi=0
shared=0
verbose=1
sanitize=0
regression_test=0

precision="double"
prefix=""
systype=""

print_help() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  CC=<compiler>          Set C compiler"
    echo "  CXX=<compiler>         Set C++ compiler"
    echo "  FC=<compiler>          Set Fortran compiler"
    echo "  --regression_test      Enable regression test"
    echo "  --debug                Enable debug mode"
    echo "  -j, --jobs=<n>         Number of parallel build jobs (default: all cores)"    
    echo "  --mpi                  Enable MPI"
    echo "  --shared               Enable shared libraries"
    echo "  --verbose              Enable verbose output"
    echo "  --sanitize             Enable sanitizer"
    echo "  --clean                Cleans build directory"
    echo "  --help                 Display this help message"
    echo "Timing logs: build/<configuration>/build-timings.log and build-timings-summary.txt"
    exit 0
}


#Using command-line arguments
for arg in "$@"
do
  case $arg in
    --help|--h)
      print_help
      ;; 
    CC=*|--CC=*)
      CC="${arg#CC=}" 
      ;;
    CXX=*|--CXX=*)
      CXX="${arg#CXX=}"
      ;;
    FC=*|--FC=*)
      FC="${arg#FC=}"   
      ;;
    --regression_test)
      regression_test=1
      ;;
    --debug)
      debug=1
      ;;
    --mpi)
      mpi=1
      ;;
    --shared)
      shared=1
      ;;
    --verbose)
      verbose=1
      ;;                
    --sanitize)
      sanitize=1
      ;; 
    -j=*|--jobs=*)
      parallel_jobs="${arg#*=}"
      ;;      
    --clean)
      echo "Cleaning build directory..."
      rm -rf ./build
      rm -rf ./local
      exit 0
      ;;
  esac
done
#This is a little confusing, but we have to move into the build dir
#and then point cmake to the top-level CMakeLists file which will
#be a directory up from the build dir hence setting
#ecosim_source_dir to ../
ecosim_source_dir=$(pwd)

############## END EDIT ####################

cmake_binary=$(command -v cmake) || {
    echo "CMake was not found on PATH. Add its bin directory to PATH before running $0." >&2
    exit 127
}

cputype=$(uname -m | sed "s/\\ /_/g")
systype=$(uname -s)

BUILDDIR="${systype}-${cputype}"

if [ "$shared" -eq 1 ]; then
    BUILDDIR="${BUILDDIR}-shared"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DBUILD_SHARED_LIBS=ON"
fi

if [ "$precision" = "double" ]; then
    BUILDDIR="${BUILDDIR}-double"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DECOSIM_PRECISION=double"
else
    BUILDDIR="${BUILDDIR}-${precision}"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DECOSIM_PRECISION=${precision}"
fi

if [ "$debug" -eq 0 ]; then
    BUILDDIR="${BUILDDIR}-Release"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_BUILD_TYPE=Release"
else
    BUILDDIR="${BUILDDIR}-Debug"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_BUILD_TYPE=Debug"
fi

if [ -n "$prefix" ]; then
    echo "add prefix here ${prefix}"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_INSTALL_PREFIX:PATH=${prefix}"
fi

if [ "$systype" == "Darwin" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DAPPLE=1"
elif [ "$systype" == "Linux" ]; then
    #Do I want this?
    #CONFIG_FLAGS="${CONFIG_FLAGS} -DLINUX=1"
    echo "We're on linux"
fi

if [ "$sanitize" -eq 1 ]; then
    BUILDDIR="${BUILDDIR}-AddressSanitizer"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DADDRESS_SANITIZER=1"
fi

if [ -n "$ATS_ECOSIM" ]; then
    echo "Building ATS-EcoSIM"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DATS_ECOSIM=1" 
    BUILDDIR="./"
    #Having this automatically set to use mpi compilers
    if [ "$mpi" -eq 0 ]; then
    	mpi=1
    fi
else
  echo "Building standalone EcoSIM"    
fi

#if mpi is set use the mpi compilers, if not use default
if [ "$mpi" -eq 1 ]; then
  # Use MPI versions of the compilers
  if [ -n "$MPICC" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_C_COMPILER=${MPICC}"
  fi

  if [ -n "$MPICXX" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_CXX_COMPILER=${MPICXX}"
  fi

  if [ -n "$MPIF90" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_Fortran_COMPILER=${MPIF90}"
  fi
else
  # Use non-MPI versions of the compilers
  if [ -n "$CC" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_C_COMPILER=${CC}"
  fi

  if [ -n "$CXX" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_CXX_COMPILER=${CXX}"
  fi

  if [ -n "$FC" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_Fortran_COMPILER=${FC}"
  fi
fi

ecosim_build_dir="build/$BUILDDIR"

# Time top-level phases and group compiler/linker launcher records by invocation.
timing_python=$(command -v python3) || { echo "Build timing requires python3." >&2; exit 1; }
timing_script="${ecosim_source_dir}/cmake/build_timer.py"
export ECOSIM_BUILD_TIMING_LOG="${ecosim_source_dir}/${ecosim_build_dir}/build-timings.log"
export ECOSIM_BUILD_RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)-$$"
timing_summary="${ecosim_source_dir}/${ecosim_build_dir}/build-timings-summary.txt"

write_timing_summary() {
  local build_status=$?
  if [ -f "$ECOSIM_BUILD_TIMING_LOG" ]; then
    "$timing_python" "$timing_script" report --log "$ECOSIM_BUILD_TIMING_LOG" \
      --run-id "$ECOSIM_BUILD_RUN_ID" > "$timing_summary"
    # Keep per-command details in the file; print only the phase summary.
    sed '/^Compile\/link steps, slowest first:/,$d' "$timing_summary"
    echo "Timing summary: $timing_summary"
  fi
  return "$build_status"
}
trap write_timing_summary EXIT

timed_step() {
  local phase=$1
  shift
  "$timing_python" "$timing_script" run --phase "$phase" -- "$@"
}

cmd_configure="${cmake_binary} \
  ${CONFIG_FLAGS}
  ${ecosim_source_dir}"

# Check if the build exists
if [ ! -d "$ecosim_build_dir" ]; then
    echo "Directory $ecosim_build_dir does not exist. Creating..."
    mkdir -p "$ecosim_build_dir"
else
    echo "Directory $ecosim_build_dir already exists."
fi

# Configure the EcoSIM build
# Note: many of the options that used to be in the CMakeLists.txt file
# have been moved here to remove redundancies

echo "cmd_configure: $cmd_configure"
echo "building in: $ecosim_build_dir"

#run the configure command
cd "${ecosim_build_dir}" || exit 1
#${cmd_configure}
pwd

timed_step configure "$cmake_binary" ${CONFIG_FLAGS} "$ecosim_source_dir"

if [ $? -ne 0 ]; then
  echo "Configuration failed, aborting."
  exit 1
fi

#This does the build
timed_step build "$cmake_binary" --build . -j "${parallel_jobs}"

if [ $? -ne 0 ]; then
  echo "Build failed, aborting."
  exit 1
fi

#Does the install
timed_step install "$cmake_binary" --build . --target install -j "${parallel_jobs}"

if [ $? -ne 0 ]; then
  echo "Install failed, aborting."
  exit 1
fi

cd ../../

mkdir -p ./local
if [ -L ./local/bin ]; then
  rm -rf ./local/bin
fi
mkdir ./local/bin

build_path=$(realpath "$ecosim_build_dir/local/bin")

for file in $(find "$build_path" -type f); do
    basename=$(basename "$file")
    link_path="./local/bin/$basename"

    # Create the symlink only if it doesn't exist
    if [ ! -e "$link_path" ]; then
        ln -s "$file" "$link_path"
        echo "Created symlink: $link_path"
    else
        echo "Symlink for $basename already exists."
    fi
done

if [ "$regression_test" -eq 1 ]; then
  timed_step regression make -C ./regression-tests test --no-print-directory ${MAKEFLAGS} compiler=gnu;
fi
