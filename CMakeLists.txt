cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # something for cygwin

project(project)
enable_language (C Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_BUILD_TYPE "Release")
# set(CMAKE_BUILD_TYPE "Debug")
# c compiler flags
#set(CMAKE_C_FLAGS "-std=c11 -Wall -Werror") # always stop on errors
set(CMAKE_C_FLAGS "-std=c11 -Wall") # only show errors
set(CMAKE_C_FLAGS_DEBUG "-ggdb -g") # for -DCMAKE_BUILD_TYPE=Debug
set(CMAKE_C_FLAGS_RELEASE "-O3") # for -DCMAKE_BUILD_TYPE=Release

set(LIBRARIES Read
              Write
              Options
              AnalysisTools
              Errors
              General)
foreach(library ${LIBRARIES})
  add_library(${library} ${library}.c)
endforeach(library)

set(UTILS AddToSystem
          Aggregates
          Aggregates-NotSameBeads
          Aggregates-wall
          AngleMolecules
          Average
          BondLength
          Config
          Config_from_xyz
          Density2D
          DensityAggregates
          DensityBox
          DensityMolecules
          DihedralMolecules
          DistrAgg
          GenSystem
          GenLayers
          GyrationAggregates
          GyrationMolecules
          Info
          JoinAggregates
          JoinRuns
          lmp_data
          lmp_vtf
          Orientation
          PairCorrel
          PotentialAggregates
          SelectedVcf
          Surface
          TransformVsf
          VisualizeAgg)
        foreach(utility ${UTILS})
  add_executable(${utility} Utility/${utility}.c)
  target_link_libraries(${utility} ${LIBRARIES} m)
endforeach(utility)

# this somehow specifies fortran compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(F_UTILS traject-v2_5
            traject-v2_6)
foreach(utility ${F_UTILS})
  add_executable(${utility} Utility/${utility}.f90)
endforeach(utility)
