# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2010, 2011, 2012 Google Inc. All rights reserved.
# http://code.google.com/p/ceres-solver/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Google Inc. nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Author: keir@google.com (Keir Mierle)
# Modified for openMVG purpose:
# - Use CXSparse by default
# - Use minilog as log library

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)

IF (COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF (COMMAND cmake_policy)

PROJECT(CERES C CXX)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Important: Always bump the second number (e.g. 1.3.x to 1.4.0) for any
# release that changes the ABI. The ABI changes for almost any modification to
# include/ceres (e.g. the public API). If you are unsure about whether
# something is an ABI change, please ask on the list.
#
# For versions without ABI changes, bump the smallest number in CERES_VERSION,
# but leave the CERES_ABI_VERSION unchanged.
SET(CERES_VERSION 1.4.0)
SET(CERES_ABI_VERSION 1.4.0)

#ENABLE_TESTING()

OPTION(BUILD_TESTING
       "Enable tests"
       OFF)

OPTION(BUILD_ANDROID
       "Build for Android. Use build_android.sh instead of setting this."
       OFF)

# To get a more static build, try the following line on Mac and Linux:
# SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

# Default locations to search for on various platforms.
# Locations to search for Eigen
SET(EIGEN_SEARCH_HEADERS ../eigen)

# Locations to search for SuiteSparse
SET(SUITESPARSE_SEARCH_LIBS ${SEARCH_LIBS})

SET(SUITESPARSE_SEARCH_HEADERS ${SEARCH_HEADERS})

SET(CXSPARSE_SEARCH_LIBS ${SEARCH_LIBS})
SET(CXSPARSE_SEARCH_HEADERS ../cxsparse/Include)

# Check for SuiteSparse dependencies
SET(AMD_FOUND FALSE)
SET(CAMD_FOUND FALSE)
SET(COLAMD_FOUND FALSE)
SET(CCOLAMD_FOUND FALSE)
SET(CHOLMOD_FOUND FALSE)
SET(SUITESPARSE_CONFIG_FOUND FALSE)
SET(UFCONFIG_FOUND FALSE)
SET(BLAS_AND_LAPACK_FOUND FALSE)

SET(SUITESPARSE_FOUND
    ${AMD_FOUND} AND
    ${CAMD_FOUND} AND
    ${COLAMD_FOUND} AND
    ${CCOLAMD_FOUND} AND
    ${CHOLMOD_FOUND} AND
    (${SUITESPARSE_CONFIG_FOUND} OR ${UFCONFIG_FOUND}) AND
    ${BLAS_AND_LAPACK_FOUND})

SET(SUITESPARSE OFF)
ADD_DEFINITIONS(-DCERES_NO_SUITESPARSE)

# By default, if all of CXSparse's dependencies are found, Ceres is
# built with CXSparse support. -DCXSPARSE=ON/OFF can be used to
# enable/disable CXSparse explicitly.
MESSAGE("-- Check for CXSparse")
SET(CXSPARSE_FOUND ON)

FIND_PATH(CXSPARSE_INCLUDE NAMES cs.h PATHS ${CXSPARSE_SEARCH_HEADERS})
IF (EXISTS ${CXSPARSE_INCLUDE})
  MESSAGE("-- Found CXSparse header in: ${CXSPARSE_INCLUDE}")
ELSE (EXISTS ${CXSPARSE_INCLUDE})
  MESSAGE("-- Did not find CXSparse header")
  SET(CXSPARSE_FOUND FALSE)
ENDIF (EXISTS ${CXSPARSE_INCLUDE})

IF (DEFINED CXSPARSE)
  IF (${CXSPARSE})
    IF (NOT ${CXSPARSE_FOUND})
      MESSAGE(FATAL_ERROR "-- CXSparse not found.")
    ENDIF (NOT ${CXSPARSE_FOUND})
  ELSE (${CXSPARSE})
    ADD_DEFINITIONS(-DCERES_NO_CXSPARSE)
  ENDIF (${CXSPARSE})
ELSE (DEFINED CXSPARSE)
  IF (${CXSPARSE_FOUND})
    MESSAGE("-- Building with CXSparse support.")
    SET(CXSPARSE ON)
  ELSE (${CXSPARSE_FOUND})
    MESSAGE("-- Building without CXSparse.")
    SET(CXSPARSE OFF)
    ADD_DEFINITIONS(-DCERES_NO_CXSPARSE)
  ENDIF (${CXSPARSE_FOUND})
ENDIF (DEFINED CXSPARSE)

ADD_DEFINITIONS(-DCERES_NO_GFLAGS)
OPTION(MINIGLOG
       "Enable MiniGlog."
       ON)

# Google Logging
IF (${MINIGLOG})
  SET(GLOG_LIB miniglog)
  MESSAGE("-- Using minimal Glog substitute for Android (library): ${GLOG_LIB}")
  SET(GLOG_INCLUDE internal/ceres/miniglog)
  MESSAGE("-- Using minimal Glog substitute for Android (include): ${GLOG_INCLUDE}")
ENDIF (${MINIGLOG})

# Eigen
MESSAGE("-- Check for Eigen 3.0")
FIND_PATH(EIGEN_INCLUDE NAMES Eigen/Core PATHS ${EIGEN_SEARCH_HEADERS})
IF (NOT EXISTS ${EIGEN_INCLUDE})
  MESSAGE(FATAL_ERROR "Can't find Eigen. Try passing -DEIGEN_INCLUDE=...")
ENDIF (NOT EXISTS ${EIGEN_INCLUDE})
MESSAGE("-- Found Eigen 3.0: ${EIGEN_INCLUDE}")

# Template specializations for the Schur complement based solvers. If
# compile time, binary size or compiler performance is an issue, you
# may consider disabling this.
OPTION(SCHUR_SPECIALIZATIONS
       "Enable fixed-size schur specializations."
       OFF)

IF (NOT ${SCHUR_SPECIALIZATIONS})
  ADD_DEFINITIONS(-DCERES_RESTRICT_SCHUR_SPECIALIZATION)
  MESSAGE("-- Disabling Schur specializations (faster compiles)")
ENDIF (NOT ${SCHUR_SPECIALIZATIONS})

# Multithreading using OpenMP
OPTION(OPENMP
       "Enable threaded solving in Ceres (requires OpenMP)"
       ON)

IF (${OPENMP})
  FIND_PACKAGE(OpenMP)
  IF(${OPENMP_FOUND})
    MESSAGE("-- Found OpenMP.")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    ADD_DEFINITIONS(-DCERES_USE_OPENMP)
  ELSE(${OPENMP_FOUND})
    MESSAGE("-- Can't find OpenMP. Continuing without it.")
    SET(OPENMP OFF)
  ENDIF(${OPENMP_FOUND})
ENDIF (${OPENMP})

# Protocol buffers
OPTION(PROTOBUF
       "Enable protocol buffers support."
       OFF)
  ADD_DEFINITIONS(-DCERES_NO_PROTOCOL_BUFFERS)

# Use threads but not on Android, where there is no support for OpenMP yet.
IF ("${UNIX}" AND NOT ${BUILD_ANDROID})
  # At least on Linux, we need pthreads to be enabled for mutex to compile.
  # This may not work on windows or android.
  FIND_PACKAGE(Threads REQUIRED)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")
  ADD_DEFINITIONS(-DCERES_HAVE_PTHREAD)
  ADD_DEFINITIONS(-DCERES_HAVE_RWLOCK)
ENDIF ("${UNIX}" AND NOT ${BUILD_ANDROID})

# Disable threads in mutex.h. Someday, after there is OpenMP support in
# Android, this can get removed. Also turn on a workaround for an NDK bug.
IF (${BUILD_ANDROID})
  ADD_DEFINITIONS(-DCERES_NO_THREADS)
  ADD_DEFINITIONS(-DCERES_WORK_AROUND_ANDROID_NDK_COMPILER_BUG)
ENDIF (${BUILD_ANDROID})

OPTION(DISABLE_TR1
       "Don't use TR1. This replaces some hash tables with sets. Slower."
       OFF)

IF (${DISABLE_TR1})
  MESSAGE("-- Replacing unordered_map/set with map/set (warning: slower!)")
  ADD_DEFINITIONS(-DCERES_NO_TR1)
ELSE (${DISABLE_TR1})
  MESSAGE("-- Using normal TR1 unordered_map/set")
  # Use the std namespace for the hash<> and related templates. This may vary by
  # system.
  IF (MSVC)
    # This is known to work with Visual Studio 2010 Express.
    ADD_DEFINITIONS("\"-DCERES_HASH_NAMESPACE_START=namespace std {\"")
    ADD_DEFINITIONS("\"-DCERES_HASH_NAMESPACE_END=}\"")
  ELSE (MSVC)
    # This is known to work with recent versions of Linux and Mac OS X.
    ADD_DEFINITIONS("\"-DCERES_HASH_NAMESPACE_START=namespace std { namespace tr1 {\"")
    ADD_DEFINITIONS("\"-DCERES_HASH_NAMESPACE_END=}}\"")
  ENDIF (MSVC)
ENDIF (${DISABLE_TR1})

INCLUDE_DIRECTORIES(
  include
  internal
  internal/ceres
  ${GLOG_INCLUDE}
  ${EIGEN_INCLUDE}
  )

IF (${CXSPARSE})
  INCLUDE_DIRECTORIES(${CXSPARSE_INCLUDE})
ENDIF(${CXSPARSE})

# Change the default build type from Debug to Release, while still
# supporting overriding the build type.
#
# The CACHE STRING logic here and elsewhere is needed to force CMake
# to pay attention to the value of these variables.
IF (NOT CMAKE_BUILD_TYPE)
  MESSAGE("-- No build type specified; defaulting to CMAKE_BUILD_TYPE=Release.")
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ELSE (NOT CMAKE_BUILD_TYPE)
  IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
    MESSAGE("\n=================================================================================")
    MESSAGE("\n-- Build type: Debug. Performance will be terrible!")
    MESSAGE("-- Add -DCMAKE_BUILD_TYPE=Release to the CMake command line to get an optimized build.")
    MESSAGE("\n=================================================================================")
  ENDIF (CMAKE_BUILD_TYPE STREQUAL "Debug")
ENDIF (NOT CMAKE_BUILD_TYPE)

# Set the default Ceres flags to an empty string.
SET (CERES_CXX_FLAGS)

IF (CMAKE_BUILD_TYPE STREQUAL "Release")
  IF (CMAKE_COMPILER_IS_GNUCXX)
    IF (${BUILD_ANDROID})
      # TODO(keir): Figure out what flags should go here to make an optimized
      # native ARM binary for Android.
    ELSE (${BUILD_ANDROID})
      # Linux
      IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        SET (CERES_CXX_FLAGS "${CERES_CXX_FLAGS} -march=native -mtune=native")
      ENDIF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
      # Mac OS X
      IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        SET (CERES_CXX_FLAGS "${CERES_CXX_FLAGS} -fast -msse3")
      ENDIF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    ENDIF (${BUILD_ANDROID})
  ENDIF (CMAKE_COMPILER_IS_GNUCXX)
  SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CERES_CXX_FLAGS}"
       CACHE STRING "Release mode flags to the C++ Compiler" FORCE)
ENDIF (CMAKE_BUILD_TYPE STREQUAL "Release")

# After the tweaks for the compile settings, disable some warnings on MSVC.
IF (MSVC)
  # Disable signed/unsigned int conversion warnings.
  ADD_DEFINITIONS("/wd4018")
  # Disable warning about using struct/class for the same symobl.
  ADD_DEFINITIONS("/wd4099")
  # Disable warning about the insecurity of using "std::copy".
  ADD_DEFINITIONS("/wd4996")
  # Disable performance warning about int-to-bool conversion.
  ADD_DEFINITIONS("/wd4800")
  # Disable performance warning about fopen insecurity.
  ADD_DEFINITIONS("/wd4996")
  # Disable warning about int64 to int32 conversion. Disabling
  # this warning may not be correct; needs investigation.
  # TODO(keir): Investigate these warnings in more detail.
  ADD_DEFINITIONS("/wd4244")
  # It's not possible to use STL types in DLL interfaces in a portable and
  # reliable way. However, that's what happens with Google Log and Google Flags
  # on Windows. MSVC gets upset about this and throws warnings that we can't do
  # much about. The real solution is to link static versions of Google Log and
  # Google Test, but that seems tricky on Windows. So, disable the warning.
  ADD_DEFINITIONS("/wd4251")

  # Google Flags doesn't have their DLL import/export stuff set up correctly,
  # which results in linker warnings. This is irrelevant for Ceres, so ignore
  # the warnings.
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ignore:4049")
ENDIF (MSVC)

# GCC is not strict enough by default, so enable most of the warnings.
IF ("${UNIX}")
  SET(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -Wno-unused -Wno-unused-parameter")
ENDIF ("${UNIX}")

ADD_SUBDIRECTORY(internal/ceres)

OPTION(BUILD_EXAMPLES "Build examples" ON)

IF (${BUILD_EXAMPLES})
  MESSAGE("-- Build the examples.")
  ADD_SUBDIRECTORY(examples)
ELSE (${BUILD_EXAMPLES})
  MESSAGE("-- Do not build any example.")
ENDIF (${BUILD_EXAMPLES})

