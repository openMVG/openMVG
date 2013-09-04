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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)
CMAKE_POLICY(VERSION 2.8)

IF (COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF (COMMAND cmake_policy)

PROJECT(CERES C CXX)

# Set up the git hook to make Gerrit Change-Id: lines in commit messages.
SET (LOCAL_GIT_DIRECTORY)
IF (EXISTS ${CMAKE_SOURCE_DIR}/.git)
  # .git directory can be found on Unix based system, or on Windows with
  # Git Bash (shipped with msysgit)
  SET (LOCAL_GIT_DIRECTORY ${CMAKE_SOURCE_DIR}/.git)
ELSE (EXISTS ${CMAKE_SOURCE_DIR}/.git)
  # TODO(keir) Add proper Windows support
ENDIF (EXISTS ${CMAKE_SOURCE_DIR}/.git)

IF (EXISTS ${LOCAL_GIT_DIRECTORY})
  IF (NOT EXISTS ${LOCAL_GIT_DIRECTORY}/hooks/commit-msg)
    # Download the hook only if it is not already present
    FILE(DOWNLOAD https://ceres-solver-review.googlesource.com/tools/hooks/commit-msg
         ${CMAKE_BINARY_DIR}/commit-msg)

    # Make the downloaded file executable, since it is not by default.
    FILE(COPY ${CMAKE_BINARY_DIR}/commit-msg
         DESTINATION ${LOCAL_GIT_DIRECTORY}/hooks/
         FILE_PERMISSIONS
           OWNER_READ OWNER_WRITE OWNER_EXECUTE
           GROUP_READ GROUP_WRITE GROUP_EXECUTE
           WORLD_READ WORLD_EXECUTE)
  ENDIF (NOT EXISTS ${LOCAL_GIT_DIRECTORY}/hooks/commit-msg)
ENDIF (EXISTS ${LOCAL_GIT_DIRECTORY})

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
SET(CERES_VERSION_MAJOR 1)
SET(CERES_VERSION_MINOR 7)
SET(CERES_VERSION_PATCH 0)
SET(CERES_VERSION
    ${CERES_VERSION_MAJOR}.${CERES_VERSION_MINOR}.${CERES_VERSION_PATCH})
SET(CERES_ABI_VERSION 1.7.0)

#ENABLE_TESTING()

OPTION(MINIGLOG "Use a stripped down version of glog" ON)
OPTION(GFLAGS "Enable Google Flags." OFF)
# Template specializations for the Schur complement based solvers. If
# compile time, binary size or compiler performance is an issue, you
# may consider disabling this.
OPTION(SCHUR_SPECIALIZATIONS "Enable fixed-size schur specializations." ON)
OPTION(CUSTOM_BLAS
       "Use handcoded BLAS routines (usually faster) instead of Eigen."
       ON)
# Multithreading using OpenMP
OPTION(OPENMP "Enable threaded solving in Ceres (requires OpenMP)" ON)
# TODO(sameeragarwal): Replace this with a positive option instead?
OPTION(DISABLE_TR1
       "Don't use TR1. This replaces some hash tables with sets. Slower."
       OFF)
# Line search minimizer is useful for large scale problems or when
# sparse linear algebra libraries are not available. If compile time,
# binary size or compiler performance is an issue, consider disabling
# this.
OPTION(LINE_SEARCH_MINIMIZER "Enable the line search minimizer." ON)
OPTION(BUILD_TESTING "Enable tests" OFF)
OPTION(BUILD_DOCUMENTATION "Build User's Guide (html)" OFF)
OPTION(BUILD_EXAMPLES "Build examples" OFF)

# Default locations to search for on various platforms.

# Libraries
LIST(APPEND CMAKE_LIBRARY_PATH /opt/local/lib)
LIST(APPEND CMAKE_LIBRARY_PATH /opt/local/lib/ufsparse) # Mac OS X
LIST(APPEND CMAKE_LIBRARY_PATH /usr/lib)
LIST(APPEND CMAKE_LIBRARY_PATH /usr/lib/atlas)
LIST(APPEND CMAKE_LIBRARY_PATH /usr/lib/suitesparse) # Ubuntu
LIST(APPEND CMAKE_LIBRARY_PATH /usr/lib64/atlas)
LIST(APPEND CMAKE_LIBRARY_PATH /usr/local/homebrew/lib) # Mac OS X
LIST(APPEND CMAKE_LIBRARY_PATH /usr/local/lib)
LIST(APPEND CMAKE_LIBRARY_PATH /usr/local/lib/suitesparse)

# Headers
LIST(APPEND CMAKE_INCLUDE_PATH ../eigen) # OpenMVG
LIST(APPEND CMAKE_INCLUDE_PATH ../cxsparse/Include) # OpenMVG
#LIST(APPEND CMAKE_INCLUDE_PATH /opt/local/include)
#LIST(APPEND CMAKE_INCLUDE_PATH /opt/local/include/ufsparse) # Mac OS X
#LIST(APPEND CMAKE_INCLUDE_PATH /opt/local/var/macports/software/eigen3/opt/local/include/eigen3) # Mac OS X
#LIST(APPEND CMAKE_INCLUDE_PATH /usr/include)
#LIST(APPEND CMAKE_INCLUDE_PATH /usr/include/eigen3) # Ubuntu 10.04's default location.
#LIST(APPEND CMAKE_INCLUDE_PATH /usr/include/suitesparse) # Ubuntu
#LIST(APPEND CMAKE_INCLUDE_PATH /usr/local/homebrew/include) # Mac OS X
#LIST(APPEND CMAKE_INCLUDE_PATH /usr/local/homebrew/include/eigen3)  # Mac OS X
#LIST(APPEND CMAKE_INCLUDE_PATH /usr/local/include)
#LIST(APPEND CMAKE_INCLUDE_PATH /usr/local/include/eigen3)
#LIST(APPEND CMAKE_INCLUDE_PATH /usr/local/include/suitesparse)

# Eigen
FIND_PATH(EIGEN_INCLUDE NAMES Eigen/Core)
IF (NOT EXISTS ${EIGEN_INCLUDE})
  MESSAGE(FATAL_ERROR "Can't find Eigen. Try passing -DEIGEN_INCLUDE=...")
ELSE (NOT EXISTS ${EIGEN_INCLUDE})
  MESSAGE("-- Found Eigen 3.x: ${EIGEN_INCLUDE}")
ENDIF (NOT EXISTS ${EIGEN_INCLUDE})

SET(BLAS_AND_LAPACK_FOUND TRUE)
IF ((NOT DEFINED LAPACK) OR (DEFINED LAPACK AND LAPACK))
  FIND_PACKAGE(LAPACK)
  IF (LAPACK_FOUND)
    MESSAGE("-- Found LAPACK library: ${LAPACK_LIBRARIES}")
  ELSE (LAPACK_FOUND)
    MESSAGE("-- Did not find LAPACK library")
    SET(BLAS_AND_LAPACK_FOUND FALSE)
  ENDIF (LAPACK_FOUND)

  FIND_PACKAGE(BLAS)
  IF (BLAS_FOUND)
    MESSAGE("-- Found BLAS library: ${BLAS_LIBRARIES}")
  ELSE (BLAS_FOUND)
    MESSAGE("-- Did not find BLAS library")
    SET(BLAS_AND_BLAS_FOUND FALSE)
  ENDIF (BLAS_FOUND)

ELSE ((NOT DEFINED LAPACK) OR (DEFINED LAPACK AND LAPACK))
  SET(BLAS_AND_LAPACK_FOUND FALSE)
ENDIF ((NOT DEFINED LAPACK) OR (DEFINED LAPACK AND LAPACK))

IF (NOT BLAS_AND_LAPACK_FOUND)
  ADD_DEFINITIONS(-DCERES_NO_LAPACK)
ENDIF (NOT BLAS_AND_LAPACK_FOUND)

IF ((NOT DEFINED SUITESPARSE) OR (DEFINED SUITESPARSE AND SUITESPARSE))
# Check for SuiteSparse dependencies

SET(AMD_FOUND TRUE)
FIND_LIBRARY(AMD_LIB NAMES amd)
IF (EXISTS ${AMD_LIB})
  MESSAGE("-- Found AMD library: ${AMD_LIB}")
ELSE (EXISTS ${AMD_LIB})
  MESSAGE("-- Did not find AMD library")
  SET(AMD_FOUND FALSE)
ENDIF (EXISTS ${AMD_LIB})

FIND_PATH(AMD_INCLUDE NAMES amd.h)
IF (EXISTS ${AMD_INCLUDE})
  MESSAGE("-- Found AMD header in: ${AMD_INCLUDE}")
ELSE (EXISTS ${AMD_INCLUDE})
  MESSAGE("-- Did not find AMD header")
  SET(AMD_FOUND FALSE)
ENDIF (EXISTS ${AMD_INCLUDE})

SET(CAMD_FOUND TRUE)
FIND_LIBRARY(CAMD_LIB NAMES camd)
IF (EXISTS ${CAMD_LIB})
  MESSAGE("-- Found CAMD library: ${CAMD_LIB}")
ELSE (EXISTS ${CAMD_LIB})
  MESSAGE("-- Did not find CAMD library")
  SET(CAMD_FOUND FALSE)
ENDIF (EXISTS ${CAMD_LIB})

FIND_PATH(CAMD_INCLUDE NAMES camd.h)
IF (EXISTS ${CAMD_INCLUDE})
  MESSAGE("-- Found CAMD header in: ${CAMD_INCLUDE}")
ELSE (EXISTS ${CAMD_INCLUDE})
  MESSAGE("-- Did not find CAMD header")
  SET(CAMD_FOUND FALSE)
ENDIF (EXISTS ${CAMD_INCLUDE})

SET(COLAMD_FOUND TRUE)
FIND_LIBRARY(COLAMD_LIB NAMES colamd)
IF (EXISTS ${COLAMD_LIB})
  MESSAGE("-- Found COLAMD library: ${COLAMD_LIB}")
ELSE (EXISTS ${COLAMD_LIB})
  MESSAGE("-- Did not find COLAMD library")
  SET(COLAMD_FOUND FALSE)
ENDIF (EXISTS ${COLAMD_LIB})

FIND_PATH(COLAMD_INCLUDE NAMES colamd.h)
IF (EXISTS ${COLAMD_INCLUDE})
  MESSAGE("-- Found COLAMD header in: ${COLAMD_INCLUDE}")
ELSE (EXISTS ${COLAMD_INCLUDE})
  MESSAGE("-- Did not find COLAMD header")
  SET(COLAMD_FOUND FALSE)
ENDIF (EXISTS ${COLAMD_INCLUDE})

SET(CCOLAMD_FOUND TRUE)
FIND_LIBRARY(CCOLAMD_LIB NAMES ccolamd)
IF (EXISTS ${CCOLAMD_LIB})
  MESSAGE("-- Found CCOLAMD library: ${CCOLAMD_LIB}")
ELSE (EXISTS ${CCOLAMD_LIB})
  MESSAGE("-- Did not find CCOLAMD library")
  SET(CCOLAMD_FOUND FALSE)
ENDIF (EXISTS ${CCOLAMD_LIB})

FIND_PATH(CCOLAMD_INCLUDE NAMES ccolamd.h)
IF (EXISTS ${CCOLAMD_INCLUDE})
  MESSAGE("-- Found CCOLAMD header in: ${CCOLAMD_INCLUDE}")
ELSE (EXISTS ${CCOLAMD_INCLUDE})
  MESSAGE("-- Did not find CCOLAMD header")
  SET(CCOLAMD_FOUND FALSE)
ENDIF (EXISTS ${CCOLAMD_INCLUDE})

SET(CHOLMOD_FOUND TRUE)
FIND_LIBRARY(CHOLMOD_LIB NAMES cholmod)
IF (EXISTS ${CHOLMOD_LIB})
  MESSAGE("-- Found CHOLMOD library: ${CHOLMOD_LIB}")
ELSE (EXISTS ${CHOLMOD_LIB})
  MESSAGE("-- Did not find CHOLMOD library")
  SET(CHOLMOD_FOUND FALSE)
ENDIF (EXISTS ${CHOLMOD_LIB})

FIND_PATH(CHOLMOD_INCLUDE NAMES cholmod.h)
IF (EXISTS ${CHOLMOD_INCLUDE})
  MESSAGE("-- Found CHOLMOD header in: ${CHOLMOD_INCLUDE}")
ELSE (EXISTS ${CHOLMOD_INCLUDE})
  MESSAGE("-- Did not find CHOLMOD header")
  SET(CHOLMOD_FOUND FALSE)
ENDIF (EXISTS ${CHOLMOD_INCLUDE})

SET(SUITESPARSEQR_FOUND TRUE)
FIND_LIBRARY(SUITESPARSEQR_LIB NAMES spqr)
IF (EXISTS ${SUITESPARSEQR_LIB})
  MESSAGE("-- Found SUITESPARSEQR library: ${SUITESPARSEQR_LIB}")
ELSE (EXISTS ${SUITESPARSEQR_LIB})
  MESSAGE("-- Did not find SUITESPARSEQR library")
  SET(SUITESPARSEQR_FOUND FALSE)
ENDIF (EXISTS ${SUITESPARSEQR_LIB})

FIND_PATH(SUITESPARSEQR_INCLUDE NAMES SuiteSparseQR.hpp)
IF (EXISTS ${SUITESPARSEQR_INCLUDE})
  MESSAGE("-- Found SUITESPARSEQR header in: ${SUITESPARSEQR_INCLUDE}")
ELSE (EXISTS ${SUITESPARSEQR_INCLUDE})
  MESSAGE("-- Did not find SUITESPARSEQR header")
  SET(SUITESPARSEQR_FOUND FALSE)
ENDIF (EXISTS ${SUITESPARSEQR_INCLUDE})

# If SuiteSparse version is >= 4 then SuiteSparse_config is required.
# For SuiteSparse 3, UFconfig.h is required.
SET(SUITESPARSE_CONFIG_FOUND TRUE)
SET(UFCONFIG_FOUND TRUE)

FIND_LIBRARY(SUITESPARSE_CONFIG_LIB NAMES suitesparseconfig)
IF (EXISTS ${SUITESPARSE_CONFIG_LIB})
  MESSAGE("-- Found SuiteSparse_config library: ${SUITESPARSE_CONFIG_LIB}")
ELSE (EXISTS ${SUITESPARSE_CONFIG_LIB})
  MESSAGE("-- Did not find SuiteSparse_config library")
ENDIF (EXISTS ${SUITESPARSE_CONFIG_LIB})

FIND_PATH(SUITESPARSE_CONFIG_INCLUDE NAMES SuiteSparse_config.h)
IF (EXISTS ${SUITESPARSE_CONFIG_INCLUDE})
  MESSAGE("-- Found SuiteSparse_config header in: ${SUITESPARSE_CONFIG_INCLUDE}")
  SET(UFCONFIG_FOUND FALSE)
ELSE (EXISTS ${SUITESPARSE_CONFIG_INCLUDE})
  MESSAGE("-- Did not find SuiteSparse_config header")
ENDIF (EXISTS ${SUITESPARSE_CONFIG_INCLUDE})

IF (NOT EXISTS ${SUITESPARSE_CONFIG_LIB} OR
    NOT EXISTS ${SUITESPARSE_CONFIG_INCLUDE})
  SET(SUITESPARSE_CONFIG_FOUND FALSE)
  FIND_PATH(UFCONFIG_INCLUDE NAMES UFconfig.h)
  IF (EXISTS ${UFCONFIG_INCLUDE})
    MESSAGE("-- Found UFconfig header in: ${UFCONFIG_INCLUDE}")
  ELSE (EXISTS ${UFCONFIG_INCLUDE})
    MESSAGE("-- Did not find UFconfig header")
    SET(UFCONFIG_FOUND FALSE)
  ENDIF (EXISTS ${UFCONFIG_INCLUDE})
ENDIF (NOT EXISTS ${SUITESPARSE_CONFIG_LIB} OR
       NOT EXISTS ${SUITESPARSE_CONFIG_INCLUDE})

FIND_LIBRARY(METIS_LIB NAMES metis)
IF (EXISTS ${METIS_LIB})
  MESSAGE("-- Found METIS library: ${METIS_LIB}")
ELSE (EXISTS ${METIS_LIB})
  MESSAGE("-- Did not find METIS library")
ENDIF (EXISTS ${METIS_LIB})

# SuiteSparseQR may be compiled with Intel Threading Building Blocks.
SET(TBB_FOUND TRUE)
FIND_LIBRARY(TBB_LIB NAMES tbb)
IF (EXISTS ${TBB_LIB})
  MESSAGE("-- Found TBB library: ${TBB_LIB}")
ELSE (EXISTS ${TBB_LIB})
  MESSAGE("-- Did not find TBB library")
  SET(TBB_FOUND FALSE)
ENDIF (EXISTS ${TBB_LIB})

FIND_LIBRARY(TBB_MALLOC_LIB NAMES tbbmalloc)
IF (EXISTS ${TBB_MALLOC_LIB})
  MESSAGE("-- Found TBB Malloc library: ${TBB_MALLOC_LIB}")
ELSE (EXISTS ${TBB_MALLOC_LIB})
  MESSAGE("-- Did not find TBB library")
  SET(TBB_FOUND FALSE)
ENDIF (EXISTS ${TBB_MALLOC_LIB})

# We don't use SET(SUITESPARSE_FOUND ${AMD_FOUND} ...) in order to be
# able to check whether SuiteSparse is available without expanding
# SUITESPARSE_FOUND with ${}. This means further checks could be:
#
#   IF (SUITESPARSE_FOUND)
#
# and not:
#
#   IF (${SUITESPARSE_FOUND})
#
IF (${AMD_FOUND} AND
    ${CAMD_FOUND} AND
    ${COLAMD_FOUND} AND
    ${CCOLAMD_FOUND} AND
    ${CHOLMOD_FOUND} AND
    (${SUITESPARSE_CONFIG_FOUND} OR ${UFCONFIG_FOUND}) AND
    ${BLAS_AND_LAPACK_FOUND})
  SET(SUITESPARSE_FOUND TRUE)
ELSE ()
  SET(SUITESPARSE_FOUND FALSE)
ENDIF ()

ENDIF ((NOT DEFINED SUITESPARSE) OR (DEFINED SUITESPARSE AND SUITESPARSE))
# By default, if all of SuiteSparse's dependencies are found, Ceres is
# built with SuiteSparse support. -DSUITESPARSE=ON/OFF can be used to
# enable/disable SuiteSparse explicitly.
IF (DEFINED SUITESPARSE)
  IF (SUITESPARSE)
    IF (NOT SUITESPARSE_FOUND)
      MESSAGE(FATAL_ERROR "One or more of SuiteSparse's dependencies was not found")
    ENDIF (NOT SUITESPARSE_FOUND)
  ELSE (SUITESPARSE)
    ADD_DEFINITIONS(-DCERES_NO_SUITESPARSE)
  ENDIF (SUITESPARSE)
ELSE (DEFINED SUITESPARSE)
  IF (SUITESPARSE_FOUND)
    MESSAGE("-- Found all SuiteSparse dependencies. Building with SuiteSparse")
    SET(SUITESPARSE ON)
  ELSE (SUITESPARSE_FOUND)
    MESSAGE("-- Did not find all SuiteSparse dependencies. Building without SuiteSparse")
    SET(SUITESPARSE OFF)
    ADD_DEFINITIONS(-DCERES_NO_SUITESPARSE)
  ENDIF (SUITESPARSE_FOUND)
ENDIF (DEFINED SUITESPARSE)

# By default, if all of CXSparse's dependencies are found, Ceres is
# built with CXSparse support. -DCXSPARSE=ON/OFF can be used to
# enable/disable CXSparse explicitly.
IF ((NOT DEFINED CXSPARSE) OR (DEFINED CXSPARSE AND CXSPARSE))

SET(CXSPARSE_FOUND ON)
# openMVG: cxsparse already exists
SET(CXSPARSE_LIB cxsparse)
#FIND_LIBRARY(CXSPARSE_LIB NAMES cxsparse)
#IF (EXISTS ${CXSPARSE_LIB})
#  MESSAGE("-- Found CXSparse library in: ${CXSPARSE_LIB}")
#ELSE (EXISTS ${CXSPARSE_LIB})
#  MESSAGE("-- Did not find CXSparse header")
#  SET(CXSPARSE_FOUND FALSE)
#ENDIF (EXISTS ${CXSPARSE_LIB})

FIND_PATH(CXSPARSE_INCLUDE NAMES cs.h)
IF (EXISTS ${CXSPARSE_INCLUDE})
  MESSAGE("-- Found CXSparse header in: ${CXSPARSE_INCLUDE}")
ELSE (EXISTS ${CXSPARSE_INCLUDE})
  MESSAGE("-- Did not find CXSparse header")
  SET(CXSPARSE_FOUND FALSE)
ENDIF (EXISTS ${CXSPARSE_INCLUDE})
ENDIF ((NOT DEFINED CXSPARSE) OR (DEFINED CXSPARSE AND CXSPARSE))

IF (DEFINED CXSPARSE)
  IF (CXSPARSE)
    IF (NOT CXSPARSE_FOUND)
      MESSAGE(FATAL_ERROR "-- CXSparse not found.")
    ENDIF (NOT CXSPARSE_FOUND)
  ELSE (CXSPARSE)
    ADD_DEFINITIONS(-DCERES_NO_CXSPARSE)
  ENDIF (CXSPARSE)
ELSE (DEFINED CXSPARSE)
  IF (CXSPARSE_FOUND)
    MESSAGE("-- Building with CXSparse support.")
    SET(CXSPARSE ON)
  ELSE (CXSPARSE_FOUND)
    MESSAGE("-- Building without CXSparse.")
    SET(CXSPARSE OFF)
    ADD_DEFINITIONS(-DCERES_NO_CXSPARSE)
  ENDIF (CXSPARSE_FOUND)
ENDIF (DEFINED CXSPARSE)

IF (GFLAGS)
  FIND_LIBRARY(GFLAGS_LIB NAMES gflags)
  IF (NOT EXISTS ${GFLAGS_LIB})
    MESSAGE(FATAL_ERROR
            "Can't find Google Flags. Please specify: "
            "-DGFLAGS_LIB=...")
  ENDIF (NOT EXISTS ${GFLAGS_LIB})
  MESSAGE("-- Found Google Flags library: ${GFLAGS_LIB}")
  FIND_PATH(GFLAGS_INCLUDE NAMES gflags/gflags.h)
  IF (NOT EXISTS ${GFLAGS_INCLUDE})
    MESSAGE(FATAL_ERROR
            "Can't find Google Flags. Please specify: "
            "-DGFLAGS_INCLUDE=...")
  ENDIF (NOT EXISTS ${GFLAGS_INCLUDE})
  MESSAGE("-- Found Google Flags header in: ${GFLAGS_INCLUDE}")
ELSE (GFLAGS)
  MESSAGE("-- Google Flags disabled; no tests or tools will be built!")
  ADD_DEFINITIONS(-DCERES_NO_GFLAGS)
ENDIF (GFLAGS)

IF (MINIGLOG)
  SET(GLOG_LIB miniglog)
  MESSAGE("-- Using minimal Glog substitute (library): ${GLOG_LIB}")
  SET(GLOG_INCLUDE internal/ceres/miniglog)
  MESSAGE("-- Using minimal Glog substitute (include): ${GLOG_INCLUDE}")
ELSE (MINIGLOG)
  FIND_LIBRARY(GLOG_LIB NAMES glog)
  IF (EXISTS ${GLOG_LIB})
    MESSAGE("-- Found Google Log library: ${GLOG_LIB}")
  ELSE (EXISTS ${GLOG_LIB})
    MESSAGE(FATAL_ERROR
            "Can't find Google Log. Please specify: -DGLOG_LIB=...")
  ENDIF (EXISTS ${GLOG_LIB})

  FIND_PATH(GLOG_INCLUDE NAMES glog/logging.h)
  IF (EXISTS ${GLOG_INCLUDE})
    MESSAGE("-- Found Google Log header in: ${GLOG_INCLUDE}")
  ELSE (EXISTS ${GLOG_INCLUDE})
    MESSAGE(FATAL_ERROR
            "Can't find Google Log. Please specify: -DGLOG_INCLUDE=...")
  ENDIF (EXISTS ${GLOG_INCLUDE})
ENDIF (MINIGLOG)

IF (NOT SCHUR_SPECIALIZATIONS)
  ADD_DEFINITIONS(-DCERES_RESTRICT_SCHUR_SPECIALIZATION)
  MESSAGE("-- Disabling Schur specializations (faster compiles)")
ENDIF (NOT SCHUR_SPECIALIZATIONS)

IF (NOT LINE_SEARCH_MINIMIZER)
  ADD_DEFINITIONS(-DCERES_NO_LINE_SEARCH_MINIMIZER)
  MESSAGE("-- Disabling line search minimizer")
ENDIF (NOT LINE_SEARCH_MINIMIZER)

IF (NOT CUSTOM_BLAS)
  ADD_DEFINITIONS(-DCERES_NO_CUSTOM_BLAS)
  MESSAGE("-- Disabling custom blas")
ENDIF (NOT CUSTOM_BLAS)

IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET(OPENMP_FOUND FALSE)
ELSE (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  IF (OPENMP)
    FIND_PACKAGE(OpenMP)
  ENDIF (OPENMP)
ENDIF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

IF (OPENMP_FOUND)
  MESSAGE("-- Found OpenMP.")
  ADD_DEFINITIONS(-DCERES_USE_OPENMP)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  IF (UNIX)
    # At least on Linux, we need pthreads to be enabled for mutex to
    # compile.  This may not work on Windows or Android.
    FIND_PACKAGE(Threads REQUIRED)
    SET(STATIC_LIBRARY_FLAGS
        "${STATIC_LIBRARY_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")
    SET(CMAKE_SHARED_LINKER_FLAGS
        "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")
    ADD_DEFINITIONS(-DCERES_HAVE_PTHREAD)
    ADD_DEFINITIONS(-DCERES_HAVE_RWLOCK)
  ENDIF (UNIX)
ELSE (OPENMP_FOUND)
  MESSAGE("-- Can't find OpenMP. Disabling multithreading.")
  ADD_DEFINITIONS(-DCERES_NO_THREADS)
ENDIF (OPENMP_FOUND)

IF (DISABLE_TR1)
  MESSAGE("-- Replacing unordered_map/set with map/set (warning: slower!)")
  ADD_DEFINITIONS(-DCERES_NO_TR1)
ELSE (DISABLE_TR1)
  MESSAGE("-- Using normal TR1 unordered_map/set")
  # Use the std namespace for the hash<> and related templates. This may vary by
  # system.
  IF (MSVC)
    IF (MSVC90)
      # Special case for Visual Studio 2008.
      # Newer versions have got tr1 symbols in another namespace,
      # and this is being handled in Else branch of this condition.
      # Probably Visual studio 2003 and 2005 also shall be handled here,
      # but don't have by hand to verify and most likely they're not
      # used by Ceres users anyway.
      ADD_DEFINITIONS("\"-DCERES_HASH_NAMESPACE_START=namespace std { namespace tr1 {\"")
      ADD_DEFINITIONS("\"-DCERES_HASH_NAMESPACE_END=}}\"")
    ELSE (MSVC90)
      # This is known to work with Visual Studio 2010 Express.
      # Further, for as long Visual Studio 2012 didn't move tr1 to
      # just another namespace, the same define will work for it as well.
      # Hopefully all further versions will also keep working with this define.
      ADD_DEFINITIONS("\"-DCERES_HASH_NAMESPACE_START=namespace std {\"")
      ADD_DEFINITIONS("\"-DCERES_HASH_NAMESPACE_END=}\"")
    ENDIF(MSVC90)
  ELSE (MSVC)
    # This is known to work with recent versions of Linux and Mac OS X.
    ADD_DEFINITIONS("\"-DCERES_HASH_NAMESPACE_START=namespace std { namespace tr1 {\"")
    ADD_DEFINITIONS("\"-DCERES_HASH_NAMESPACE_END=}}\"")
  ENDIF (MSVC)
ENDIF (DISABLE_TR1)

INCLUDE_DIRECTORIES(
  include
  internal
  internal/ceres
  ${GLOG_INCLUDE}
  ${EIGEN_INCLUDE}
  )

FILE(GLOB CERES_HDRS ${CMAKE_SOURCE_DIR}/include/ceres/*.h)
INSTALL(FILES ${CERES_HDRS} DESTINATION include/ceres)

FILE(GLOB CERES_PUBLIC_INTERNAL_HDRS ${CMAKE_SOURCE_DIR}/include/ceres/internal/*.h)
INSTALL(FILES ${CERES_PUBLIC_INTERNAL_HDRS} DESTINATION include/ceres/internal)

IF (SUITESPARSE)
  INCLUDE_DIRECTORIES(${AMD_INCLUDE})
  INCLUDE_DIRECTORIES(${CAMD_INCLUDE})
  INCLUDE_DIRECTORIES(${COLAMD_INCLUDE})
  INCLUDE_DIRECTORIES(${CCOLAMD_INCLUDE})
  INCLUDE_DIRECTORIES(${CHOLMOD_INCLUDE})
  INCLUDE_DIRECTORIES(${SUITESPARSEQR_INCLUDE})
  IF (SUITESPARSE_CONFIG_FOUND)
    INCLUDE_DIRECTORIES(${SUITESPARSE_CONFIG_INCLUDE})
  ENDIF (SUITESPARSE_CONFIG_FOUND)
  IF (UFCONFIG_FOUND)
    INCLUDE_DIRECTORIES(${UFCONFIG_INCLUDE})
  ENDIF (UFCONFIG_FOUND)
ENDIF (SUITESPARSE)

IF (CXSPARSE)
  INCLUDE_DIRECTORIES(${CXSPARSE_INCLUDE})
ENDIF (CXSPARSE)

IF (GFLAGS)
  INCLUDE_DIRECTORIES(${GFLAGS_INCLUDE})
ENDIF (GFLAGS)

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
    # Linux
    IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
      SET (CERES_CXX_FLAGS "${CERES_CXX_FLAGS} -march=native -mtune=native")
    ENDIF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    # Mac OS X
    IF (CMAKE_SYSTEM_NAME MATCHES "Darwin")
      SET (CERES_CXX_FLAGS "${CERES_CXX_FLAGS} -msse3")
      # Use of -fast only applicable for Apple's GCC
      # Assume this is being used if GCC version < 4.3 on OSX
      EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
      IF (GCC_VERSION VERSION_LESS 4.3)
        SET (CERES_CXX_FLAGS "${CERES_CXX_FLAGS} -fast")
      ENDIF (GCC_VERSION VERSION_LESS 4.3)
    ENDIF (CMAKE_SYSTEM_NAME MATCHES "Darwin")
  ENDIF (CMAKE_COMPILER_IS_GNUCXX)
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Use of -O4 requires use of gold linker & LLVM-gold plugin, which might
    # well not be present / in use and without which files will compile, but
    # not link ('file not recognized') so explicitly check for support
    INCLUDE(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-O4" HAVE_LTO_SUPPORT)
    IF (HAVE_LTO_SUPPORT)
      MESSAGE(STATUS "Enabling link-time optimization (-O4)")
      SET(CERES_CXX_FLAGS "${CERES_CXX_FLAGS} -O4")
    ELSE ()
      MESSAGE(STATUS "Compiler/linker does not support link-time optimization (-O4), disabling.")
    ENDIF (HAVE_LTO_SUPPORT)
  ENDIF ()
ENDIF (CMAKE_BUILD_TYPE STREQUAL "Release")

SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CERES_CXX_FLAGS}")

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

  # Tuple sizes of 10 are used by Gtest.
  ADD_DEFINITIONS("-D_VARIADIC_MAX=10")
ENDIF (MSVC)

IF (UNIX)
  # GCC is not strict enough by default, so enable most of the warnings.
  SET(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -Wno-unused-parameter -Wno-missing-field-initializers")
ENDIF (UNIX)

# Use a larger inlining threshold for Clang, since it hobbles Eigen,
# resulting in an unreasonably slow version of the blas routines. The
# -Qunused-arguments is needed because CMake passes the inline
# threshold to the linker and clang complains about it and dies.
IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Qunused-arguments -mllvm -inline-threshold=600 -Wno-return-type-c-linkage")
ENDIF ()

ADD_SUBDIRECTORY(internal/ceres)

IF (BUILD_DOCUMENTATION)
  MESSAGE("-- Documentation building is enabled")

  # Make CMake aware of the cmake folder, in order to find 'FindSphinx.cmake'
  SET (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

  # Generate the User's Guide (html).
  # The corresponding target is UserGuide, but is included in ALL.
  ADD_SUBDIRECTORY(docs)
ENDIF (BUILD_DOCUMENTATION)

IF (BUILD_EXAMPLES)
  MESSAGE("-- Build the examples.")
  ADD_SUBDIRECTORY(examples)
ELSE (BUILD_EXAMPLES)
  MESSAGE("-- Do not build any example.")
ENDIF (BUILD_EXAMPLES)


