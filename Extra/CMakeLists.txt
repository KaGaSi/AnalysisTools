cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Cygwin compatibility setting

project(project)
enable_language (C)
find_package(GSL REQUIRED)  # Find GSL package
message(STATUS "GSL package found: ${GSL_INCLUDE_DIRS} ${GSL_LIBRARIES}")

# Set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../build/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type
# set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_BUILD_TYPE "Release")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# compiler flags
set(CMAKE_C_FLAGS "-std=c11 -Wall")
set(CMAKE_C_FLAGS_DEBUG "-ggdb -g")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")

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
  ../src/AnalysisTools.c
  ../src/Aggregates.c
  ../src/Arrays.c
  ../src/Coordinates.c
  ../src/Debug.c
  ../src/Errors.c
  ../src/General.c
  ../src/Globals.c
  ../src/MainLoop.c
  ../src/Options.c
  ../src/PairHash.c
  ../src/Pairs.c
  ../src/pcg32.c
  ../src/ReadWrite.c
  ../src/ReadWriteVtf.c
  ../src/ReadWriteXyz.c
  ../src/ReadWriteLtrj.c
  ../src/ReadWriteLdata.c
  ../src/ReadWriteField.c
  ../src/ReadWriteConfig.c
  ../src/ReadWriteGromacs.c
  ../src/System.c
  ../src/Structs.c)

# List of utilities
set(UTILS
  # CICount
  # Distribution
  # dbscan
  # Entropy
  CA2Middle
  FuzzyCoor
  GenGal
  GalContacts
  LargestAgg
  # NearestNeighbour
  # NearestNeighbour-wall
  RandomChain
  Structure
  # SurfaceTension
  )

# Conditionally link BLAS
target_link_libraries(AnalysisTools PRIVATE ${GSL_LIBRARIES} m ${BLAS_LIB})
target_include_directories(AnalysisTools PRIVATE ${GSL_INCLUDE_DIRS})

# Create utility executables
foreach(utility ${UTILS})
  add_executable(${utility} ../Extra/${utility}.c)
  target_link_libraries(${utility} PRIVATE ${GSL_LIBRARIES} m ${BLAS_LIB} AnalysisTools)
  target_include_directories(${utility} PRIVATE ${GSL_INCLUDE_DIRS})
endforeach()
