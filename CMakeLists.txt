cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # something for cygwin

project(project)
enable_language (C Fortran)
find_package(GSL REQUIRED)  # Find GSL package
message(STATUS "GSL package found: ${GSL_INCLUDE_DIRS} ${GSL_LIBRARIES}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_BUILD_TYPE "Release")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# c compiler flags
#set(CMAKE_C_FLAGS "-std=c11 -Wall") # only show errors
set(CMAKE_C_FLAGS "-std=c11 -Wall")
set(CMAKE_C_FLAGS_DEBUG "-ggdb -g") # for -DCMAKE_BUILD_TYPE=Debug
# set(CMAKE_C_FLAGS_RELEASE "-Werror -O3") # for -DCMAKE_BUILD_TYPE=Release
set(CMAKE_C_FLAGS_RELEASE "-O3") # for -DCMAKE_BUILD_TYPE=Release

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# # Ensure PIC is used (this is typically unnecessary for shared libraries)
# set_target_properties(mylib PROPERTIES POSITION_INDEPENDENT_CODE ON)
# add_library(AnalysisTools SHARED
add_library(AnalysisTools
            AnalysisTools.c
            Globals.c
            MathUtils.h
            ReadWrite.c
            ReadWriteVtf.c
            ReadWriteXyz.c
            ReadWriteLtrj.c
            ReadWriteLdata.c
            ReadWriteField.c
            ReadWriteConfig.c
            ReadWriteGromacs.c
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
          BOOP
          DensityAggregates
          DensityBox
          DistrAgg
          GenFIELD
          GyrationAggregates
          Info
          JoinSystems
          PairCorrel
          PersistenceLength
          Selected
          Surface
          VisualizeAgg)
foreach(utility ${UTILS})
  add_executable(${utility} Utility/${utility}.c)
  # target_link_libraries(${utility} ${LIBRARIES} m)
  # target_link_libraries(${utility} PRIVATE AnalysisTools m)
  target_link_libraries(${utility} PRIVATE ${GSL_LIBRARIES} AnalysisTools m)
  target_include_directories(${utility} PRIVATE ${GSL_INCLUDE_DIRS})
endforeach(utility)
# # Create the executable for BOOP before linking libraries
# add_executable(BOOP Utility/BOOP.c)  # Ensure the correct source file is specified
# # target_link_libraries(BOOP PRIVATE GSL::GSL)  # Link GSL to BOOP
# target_link_libraries(BOOP PRIVATE ${GSL_LIBRARIES} AnalysisTools m)
# target_include_directories(BOOP PRIVATE ${GSL_INCLUDE_DIRS})


# this somehow specifies fortran compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(F_UTILS traject-v2_7)
foreach(utility ${F_UTILS})
  add_executable(${utility} Utility/${utility}.f90)
endforeach(utility)
