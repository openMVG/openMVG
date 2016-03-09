cmake_minimum_required(VERSION 2.6)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

project(flann)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)

include(${PROJECT_SOURCE_DIR}/cmake/flann_utils.cmake)
set(FLANN_VERSION 1.8.4)
DISSECT_VERSION()
GET_OS_INFO()

# CMAKE_C_COMPILER_ID and CMAKE_CXX_COMPILER_ID are documented as internal variables subject to
# change at: http://www.cmake.org/cmake/help/v2.8.8/cmake.html#variable:CMAKE_LANG_COMPILER_ID
# however currently this seems the most robust way to detect the presence of the Clang compiler
if(CMAKE_C_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_COMPILER_IS_CLANG 1)
endif ()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_COMPILER_IS_CLANGXX 1)
endif ()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Add an "uninstall" target
#CONFIGURE_FILE ("${PROJECT_SOURCE_DIR}/cmake/uninstall_target.cmake.in"
#    "${PROJECT_BINARY_DIR}/uninstall_target.cmake" IMMEDIATE @ONLY)
#ADD_CUSTOM_TARGET (uninstall "${CMAKE_COMMAND}" -P
#    "${PROJECT_BINARY_DIR}/uninstall_target.cmake")

# Set the build type.  Options are:
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    #set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
    #set(CMAKE_BUILD_TYPE Debug)
endif()

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
# set output path for tests
set(TEST_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/test)

option(BUILD_C_BINDINGS "Build C bindings" OFF)
#option(BUILD_PYTHON_BINDINGS "Build Python bindings" OFF)
#option(BUILD_MATLAB_BINDINGS "Build Matlab bindings" OFF)
option(BUILD_CUDA_LIB "Build CUDA library" OFF)
option(USE_OPENMP "Use OpenMP multi-threading" ON)
#option(USE_MPI "Use MPI" OFF)

set(NVCC_COMPILER_BINDIR "" CACHE PATH  "Directory where nvcc should look for C++ compiler. This is passed to nvcc through the --compiler-bindir option.")

#if (NOT BUILD_C_BINDINGS)
#    set(BUILD_PYTHON_BINDINGS OFF)
#    set(BUILD_MATLAB_BINDINGS OFF)
#endif()


if (USE_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    else()
        message(WARNING "OpenMP NOT found")
        set(USE_OPENMP OFF)
    endif()
endif()


#set the C/C++ include path to the "include" directory
include_directories(${PROJECT_SOURCE_DIR}/src/cpp)

# require proper c++
#add_definitions( "-Wall -ansi -pedantic" )
# HDF5 uses long long which is not ansi
if (WIN32)
    # lots of warnings with cl.exe right now, use /W1
    add_definitions("/W1 -D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS /bigobj")
else(WIN32)
    add_definitions( "-Wall -Wno-unknown-pragmas -Wno-unused-function" )
endif(WIN32)

add_subdirectory( cmake )
add_subdirectory( src )

