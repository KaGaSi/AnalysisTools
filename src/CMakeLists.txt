cmake_minimum_required(VERSION 3.16...4.0)

project(project)
enable_language(C Fortran)

# Set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type - Release unless the user says otherwise
# (e.g. cmake ../src -DCMAKE_BUILD_TYPE=Debug)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# compiler flags (appended so user-supplied CFLAGS/CMAKE_C_FLAGS survive)
string(APPEND CMAKE_C_FLAGS " -std=c23 -Wall")
set(CMAKE_C_FLAGS_DEBUG "-ggdb -g")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")

# Optional sanitizers (appended so the hardcoded flags above are preserved):
#   cmake ../src -DSANITIZE=ON
option(SANITIZE "Build with Address + UndefinedBehaviour sanitizers" OFF)
if(SANITIZE)
  # -fno-sanitize-recover=all makes UBSan abort (non-zero exit) on the first
  # violation instead of just printing and continuing, so CI actually fails
  string(APPEND CMAKE_C_FLAGS
         " -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fsanitize=address,undefined")
  message(STATUS "Sanitizers: address, undefined")
endif()

# Link-time optimisation, so the ArrND accessors in Arrays.c can be inlined
# into the utilities that use them (they're in a separate translation unit,
# and some are called in per-bead loops). Optimised configurations only -
# it does nothing useful for Debug and interferes with the sanitizers.
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
if(ipo_supported AND NOT SANITIZE)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
  message(STATUS "IPO/LTO: enabled for Release/RelWithDebInfo")
elseif(NOT ipo_supported)
  message(STATUS "IPO/LTO: not supported (${ipo_error})")
endif()

# Find GSL
find_package(GSL REQUIRED)
find_library(GSL_LIBRARIES NAMES gsl)
find_path(GSL_INCLUDE_DIRS NAMES gsl/gsl_vector.h)

if(NOT GSL_LIBRARIES OR NOT GSL_INCLUDE_DIRS)
  message(FATAL_ERROR "GSL not found!")
endif()

message(STATUS "Found GSL: ${GSL_LIBRARIES}")
message(STATUS "GSL include dirs: ${GSL_INCLUDE_DIRS}")

# Find BLAS (try GSL's linked BLAS first)
execute_process(
  COMMAND ldd ${GSL_LIBRARIES} OUTPUT_VARIABLE LDD_OUTPUT ERROR_QUIET
)
string(REGEX MATCH ".*lib(blas|openblas).so[^ ]*" BLAS_PATH "${LDD_OUTPUT}")

if(BLAS_PATH)
  string(REGEX REPLACE ".* => ([^ ]*).*" "\\1" BLAS_LIB ${BLAS_PATH})
  message(STATUS "Detected BLAS used by GSL: ${BLAS_LIB}")
else()
  message(WARNING "Could not auto-detect BLAS via ldd. Searching manually.")
  find_library(BLAS_LIB NAMES openblas blas)
  if(NOT BLAS_LIB)
    message(WARNING "No BLAS library found, skipping BLAS linking.")
  else()
    message(STATUS "Found BLAS manually: ${BLAS_LIB}")
  endif()
endif()

# If no BLAS was found, we won’t link it
if(NOT BLAS_LIB)
  message(STATUS "Skipping BLAS linking.")
endif()

# Define library
add_library(AnalysisTools
  Aggregates.c
  Arrays.c
  Coordinates.c
  Debug.c
  Errors.c
  General.c
  Globals.c
  MainLoop.c
  MathUtils.h
  Options.c
  PairHash.c
  Pairs.c
  pcg32.c
  ReadLibrary.c
  ReadWrite.c
  ReadWriteConfig.c
  ReadWriteField.c
  ReadWriteGromacs.c
  ReadWriteLdata.c
  ReadWriteLtrj.c
  ReadWriteVtf.c
  ReadWriteXyz.c
  Structs.c
  ReadWriteSysInfo.c
  System.c)

# List of utilities
set(UTILS
  AddToSystem
  AggNumber
  BilayerInPlane
  BilayerRigidity
  BilayerSq
  Aggregates
  AggResidence
  AngleMolecules
  Average
  BondLength
  DensityAggregates
  DensityBox
  DiamondNetwork
  GyrationAggregates
  GyrationMolecules
  Info
  JoinSystems
  MolOverlap
  OrientOrder
  PairCorrel
  PersistenceLength
  Selected
  Surface
  VisualizeAgg)

# Conditionally link BLAS
target_link_libraries(AnalysisTools PRIVATE ${GSL_LIBRARIES} m ${BLAS_LIB})
target_include_directories(AnalysisTools PRIVATE ${GSL_INCLUDE_DIRS})

# Create utility executables
foreach(utility ${UTILS})
  add_executable(${utility} ../Utility/${utility}.c)
  target_link_libraries(${utility} PRIVATE ${GSL_LIBRARIES} m ${BLAS_LIB} AnalysisTools)
  target_include_directories(${utility} PRIVATE ${GSL_INCLUDE_DIRS})
endforeach()

# Fortran utilities
set(F_UTILS traject-v2_7)
foreach(utility ${F_UTILS})
  add_executable(${utility} ../Utility/${utility}.f90)
endforeach()

# ---- unit tests (CTest) ----
# Run with: cmake ../src && make && ctest --output-on-failure
# or under sanitizers: cmake ../src -DSANITIZE=ON && make && ctest --output-on-failure
enable_testing()
set(TESTS
  t_geometry
  t_traversal
  t_readers
  t_writers
  t_aggregates
  t_arrays
  t_general
  t_system
  t_library
  t_coordinates)
foreach(t ${TESTS})
  add_executable(${t} ../tests/${t}.c)
  target_link_libraries(${t} PRIVATE ${GSL_LIBRARIES} m ${BLAS_LIB} AnalysisTools)
  target_include_directories(${t} PRIVATE ${GSL_INCLUDE_DIRS})
  # absolute paths to the committed reader fixtures and to the shipped
  # examples (t_library checks the example molecule library still parses)
  target_compile_definitions(${t} PRIVATE
    FIXTURES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../tests/fixtures"
    EXAMPLES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../Examples")
  add_test(NAME ${t} COMMAND ${t})
  # Make sanitizer findings fatal so a mutation that slips past validation
  # aborts (SIGABRT) instead of printing a recoverable "runtime error" and
  # continuing - GCC's UBSan array-bounds check recovers by default, which
  # would otherwise let out-of-bounds indexing pass the reader fuzz cases.
  set_tests_properties(${t} PROPERTIES ENVIRONMENT
    "UBSAN_OPTIONS=halt_on_error=1:abort_on_error=1;ASAN_OPTIONS=abort_on_error=1")
endforeach()

# Ensure Fortran compiler is set
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

