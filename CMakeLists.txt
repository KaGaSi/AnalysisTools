cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # something for cygwin

project(project)
enable_language (C Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_BUILD_TYPE "Debug")
#set(CMAKE_BUILD_TYPE "Release")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# c compiler flags
#set(CMAKE_C_FLAGS "-std=c11 -Wall") # only show errors
set(CMAKE_C_FLAGS "-std=c11 -Wall -Wno-format-truncation -Werror")
set(CMAKE_C_FLAGS_DEBUG "-ggdb -g") # for -DCMAKE_BUILD_TYPE=Debug
set(CMAKE_C_FLAGS_RELEASE "-O3") # for -DCMAKE_BUILD_TYPE=Release

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(LIBRARIES Read
              Write
              Options
              AnalysisTools
              Errors
              General
              Structs)
foreach(library ${LIBRARIES})
  add_library(${library} ${library}.c)
endforeach(library)

set(UTILS AddToSystem
          AddSystemToSystem
          AddRandomToSystem
          Aggregates
          Aggregates-NotSameBeads
          AngleMolecules
          Average
          BondLength
          Config
          Config_from_xyz
          DensityAggregates
          DensityBox
          DensityMolecules
          DihedralMolecules
          DistrAgg
          GenLayers
          GenAnglesDihedrals
          GyrationAggregates
          GyrationMolecules
          Info
          JoinAggregates
#         JoinRuns
          lmp_data
          lmp_vtf
          NumDeriv
          OrientationMolecules
#         Orientation-2
          OrientationAggregates
          PairCorrel
          PotentialAggregates
          SelectedVcf
          SuttonChenParameters
          Surface
          Velocity
          VisualizeAgg)
foreach(utility ${UTILS})
  add_executable(${utility} Utility/${utility}.c)
  target_link_libraries(${utility} ${LIBRARIES} m)
endforeach(utility)

#set(TESTING Coordinates)
#          foreach(utility ${TESTING})
#  add_executable(${utility} Testing/${utility}.c)
#  target_link_libraries(${utility} m ${LIBRARIES})
#endforeach(utility)

# this somehow specifies fortran compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(F_UTILS traject-v2_7)
foreach(utility ${F_UTILS})
  add_executable(${utility} Utility/${utility}.f90)
endforeach(utility)
