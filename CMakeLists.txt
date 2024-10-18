cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # something for cygwin

project(project)
enable_language (C Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_BUILD_TYPE "Debug")
# set(CMAKE_BUILD_TYPE "Release")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# c compiler flags
#set(CMAKE_C_FLAGS "-std=c11 -Wall") # only show errors
set(CMAKE_C_FLAGS "-std=c11 -Wall")
set(CMAKE_C_FLAGS_DEBUG "-ggdb -g") # for -DCMAKE_BUILD_TYPE=Debug
# set(CMAKE_C_FLAGS_RELEASE "-Werror -O3") # for -DCMAKE_BUILD_TYPE=Release
set(CMAKE_C_FLAGS_RELEASE "-O3") # for -DCMAKE_BUILD_TYPE=Release

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(AnalysisTools
            AnalysisTools.c
            Read.c
            Write.c
            Options.c
            System.c
            Errors.c
            General.c
            Structs.c)

set(UTILS AddToSystem
          Aggregates
          Aggregates-NotSameBeads
          AngleMolecules
          Average
          BondLength
          DensityAggregates
          DensityBox
          DistrAgg
          GenFIELD
          Info
          JoinSystems
          PersistenceLength
          Selected
          Surface
          VisualizeAgg)
foreach(utility ${UTILS})
  add_executable(${utility} Utility/${utility}.c)
  # target_link_libraries(${utility} ${LIBRARIES} m)
  target_link_libraries(${utility} AnalysisTools m)
endforeach(utility)

# this somehow specifies fortran compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(F_UTILS traject-v2_7)
foreach(utility ${F_UTILS})
  add_executable(${utility} Utility/${utility}.f90)
endforeach(utility)
