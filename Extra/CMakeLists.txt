cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # something for cygwin

project(project)
enable_language (C)
find_package(GSL REQUIRED)  # Find GSL package
message(STATUS "GSL package found: ${GSL_INCLUDE_DIRS} ${GSL_LIBRARIES}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../build/bin)
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

add_library(AnalysisTools
            ../AnalysisTools.c
            ../Globals.c
            ../MathUtils.h
            ../ReadWrite.c
            ../ReadWriteVtf.c
            ../ReadWriteXyz.c
            ../ReadWriteLtrj.c
            ../ReadWriteLdata.c
            ../ReadWriteField.c
            ../ReadWriteConfig.c
            ../ReadWriteGromacs.c
            ../Options.c
            ../System.c
            ../Errors.c
            ../General.c
            ../Structs.c)

set(UTILS GenGal
          GalContacts
          LargestAgg
          Template
          test)
foreach(utility ${UTILS})
  add_executable(${utility} ../Extra/${utility}.c)
  target_link_libraries(${utility} PRIVATE ${GSL_LIBRARIES} AnalysisTools m)
  target_include_directories(${utility} PRIVATE ${GSL_INCLUDE_DIRS})
endforeach(utility)
