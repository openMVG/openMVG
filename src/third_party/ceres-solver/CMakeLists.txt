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
# Authors: keir@google.com (Keir Mierle)
#          alexs.mac@gmail.com (Alex Stewart)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)
CMAKE_POLICY(VERSION 2.8)

IF (COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF (COMMAND cmake_policy)

PROJECT(CERES C CXX)

# Make CMake aware of the cmake folder for local FindXXX scripts,
# append rather than set in case the user has passed their own
# additional paths via -D.
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
INCLUDE(UpdateCacheVariable)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# Set postfixes for generated libraries based on buildtype.
SET(CMAKE_RELEASE_POSTFIX "")
SET(CMAKE_DEBUG_POSTFIX "-debug")

# Important: Always bump the second number (e.g. 1.3.x to 1.4.0) for any
# release that changes the ABI. The ABI changes for almost any modification to
# include/ceres (e.g. the public API). If you are unsure about whether
# something is an ABI change, please ask on the list.
#
# For versions without ABI changes, bump the smallest number in CERES_VERSION,
# but leave the CERES_ABI_VERSION unchanged.
SET(CERES_VERSION_MAJOR 1)
SET(CERES_VERSION_MINOR 10)
SET(CERES_VERSION_PATCH 0)
SET(CERES_VERSION
    ${CERES_VERSION_MAJOR}.${CERES_VERSION_MINOR}.${CERES_VERSION_PATCH})
SET(CERES_ABI_VERSION 1.10.0)

ENABLE_TESTING()

OPTION(MINIGLOG "Use a stripped down version of glog." ON)
OPTION(GFLAGS "Enable Google Flags." OFF)
OPTION(SUITESPARSE "Enable SuiteSparse." OFF)
OPTION(CXSPARSE "Enable CXSparse." ON)
OPTION(LAPACK "Enable use of LAPACK." ON)
# Template specializations for the Schur complement based solvers. If
# compile time, binary size or compiler performance is an issue, you
# may consider disabling this.
OPTION(SCHUR_SPECIALIZATIONS "Enable fixed-size schur specializations." ON)
OPTION(CUSTOM_BLAS
       "Use handcoded BLAS routines (usually faster) instead of Eigen."
       ON)
# Multithreading using OpenMP
OPTION(OPENMP "Enable threaded solving in Ceres (requires OpenMP)" ON)
OPTION(EIGENSPARSE
  "Enable the use of Eigen as a sparse linear algebra library for
   solving the nonlinear least squares problems. Enabling this
   option will result in an LGPL licensed version of Ceres Solver
   as the Simplicial Cholesky factorization in Eigen is licensed under the LGPL.
   This does not affect the covariance estimation algorithm, as it
   depends on the sparse QR factorization algorithm, which is licensed
   under the MPL."
  OFF)
OPTION(BUILD_TESTING "Enable tests" OFF)
OPTION(BUILD_DOCUMENTATION "Build User's Guide (html)" OFF)
OPTION(BUILD_EXAMPLES "Build examples" OFF)
OPTION(BUILD_SHARED_LIBS "Build Ceres as a shared library." OFF)
IF (MSVC)
  OPTION(MSVC_USE_STATIC_CRT
    "MS Visual Studio: Use static C-Run Time Library in place of shared." OFF)

  IF (BUILD_TESTING AND BUILD_SHARED_LIBS)
    MESSAGE(
      "-- Disabling tests. The flags BUILD_TESTING and BUILD_SHARED_LIBS"
      " are incompatible with MSVC."
      )
    UPDATE_CACHE_VARIABLE(BUILD_TESTING OFF)
  ENDIF (BUILD_TESTING AND BUILD_SHARED_LIBS)
ENDIF (MSVC)

# Use ios-cmake to build a static library for iOS
#
# We need to add isysroot to force cmake to find the toolchains from the iOS SDK
# instead of using the standard ones. And add flag mios-simulator-version so clang
# knows we are building for ios simulator but not mac.
#
# You can build for OS (armv7, armv7s, arm64), SIMULATOR (i386) or SIMULATOR64 (x86_64)
# separately and use lipo to merge them into one static library.
#
# There are some features/algorithms are not available in iOS version and the
# minimum supported iOS version is 6.0 now.
#
# Use cmake ../ceres-solver -DCMAKE_TOOLCHAIN_FILE=../ceres-solver/cmake/iOS.cmake \
# -DIOS_PLATFORM=PLATFORM -DEIGEN_INCLUDE_DIR=/path/to/eigen/header
# to config the cmake. The PLATFORM can be one of OS, SIMULATOR and SIMULATOR64.
# Check the documentation in iOS.cmake to find more options.
#
# After building, you will get a single library: libceres.a, which
# you need to add to your Xcode project.
#
# If you use the lapack and blas, then you also need to add Accelerate.framework
# to your Xcode project's linking dependency.
IF (IOS)
  MESSAGE(STATUS "Building Ceres for iOS platform: ${IOS_PLATFORM}")

  UPDATE_CACHE_VARIABLE(MINIGLOG ON)
  MESSAGE(STATUS "Building for iOS, forcing use of miniglog instead of glog.")

  UPDATE_CACHE_VARIABLE(SUITESPARSE OFF)
  UPDATE_CACHE_VARIABLE(CXSPARSE ON)
  UPDATE_CACHE_VARIABLE(GFLAGS OFF)
  UPDATE_CACHE_VARIABLE(OPENMP OFF)

  MESSAGE(STATUS "Building for iOS: SuiteSparse, CXSparse, gflags and OpenMP are not available.")

  UPDATE_CACHE_VARIABLE(BUILD_EXAMPLES OFF)
  MESSAGE(STATUS "Building for iOS, will not build examples.")

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fobjc-abi-version=2 -fobjc-arc -isysroot ${CMAKE_OSX_SYSROOT}")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fobjc-abi-version=2 -fobjc-arc -isysroot ${CMAKE_OSX_SYSROOT}")

  IF (${IOS_PLATFORM} STREQUAL "SIMULATOR" OR ${IOS_PLATFORM} STREQUAL "SIMULATOR64")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mios-simulator-version-min=6.0")
  ENDIF()
ENDIF (IOS)

# Prior to October 2013, Ceres used some non-CMake standardised variables to
# hold user-specified (as opposed to FindPackage found) include directory and
# library paths for Ceres dependencies.  These were were of the form:
# <DEPENDENCY>_LIB / <DEPENDENCY>_INCLUDE.  Since then, Ceres now has
# FindPackage() scripts for all of its dependencies which obey the standard
# CMake variables: <DEPENDENCY>_LIBRARIES & <DEPENDENCY>_INCLUDE_DIRS.  In order
# to ensure backwards compatibility, we use convert any legacy variables to
# _directory_ hints for the FindPackage() scripts.
MACRO(HANDLE_LEGACY_INCLUDE_DEPENDENCY_HINT
    LEGACY_VAR DIRECTORY_HINT_VAR)
  IF (DEFINED ${LEGACY_VAR})
    # Get the dependency name (all caps) from the hint directory variable
    # for the warning message.
    STRING(REGEX MATCH "^[^_]*" DEPENDENCY_NAME ${DIRECTORY_HINT_VAR})
    MESSAGE(WARNING "You are defining a legacy variable ${LEGACY_VAR} "
      "to specify the include directory for ${DEPENDENCY_NAME}.  This is "
      "deprecated and support for it will be removed in a future release. "
      "Please use either the search directory hints variable: "
      "${DIRECTORY_HINT_VAR} or ${DEPENDENCY_NAME}_INCLUDE_DIR to specify "
      "exactly the directory used (no search performed), see: "
      "http://homes.cs.washington.edu/~sagarwal/ceres-solver/dev/building.html "
      "for more information.")
    LIST(APPEND ${DIRECTORY_HINT_VAR} ${${LEGACY_VAR}})
  ENDIF (DEFINED ${LEGACY_VAR})
ENDMACRO(HANDLE_LEGACY_INCLUDE_DEPENDENCY_HINT)

MACRO(HANDLE_LEGACY_LIBRARY_DEPENDENCY_HINT
    LEGACY_VAR DIRECTORY_HINT_VAR)
  IF (DEFINED ${LEGACY_VAR})
    # Get the dependency name (all caps) from the hint directory variable
    # for the warning message.
    STRING(REGEX MATCH "^[^_]*" DEPENDENCY_NAME ${DIRECTORY_HINT_VAR})
    MESSAGE(WARNING "You are defining a legacy variable ${LEGACY_VAR} "
      "to specify the library for ${DEPENDENCY_NAME}.  This is "
      "deprecated and support for it will be removed in a future release. "
      "Please use either the search directory hints variable: "
      "${DIRECTORY_HINT_VAR} or ${DEPENDENCY_NAME}_LIBRARY to specify "
      "exactly the library used (no search performed), see: "
      "http://homes.cs.washington.edu/~sagarwal/ceres-solver/dev/building.html "
      "for more information.")
    IF (EXISTS ${${LEGACY_VAR}} AND
        NOT IS_DIRECTORY ${${LEGACY_VAR}})
      # User specified an explicit (library) file using the legacy variable
      # interface, hints to FindPackage() scripts are directories so add the
      # parent directory of the specified file.
      GET_FILENAME_COMPONENT(DIR_HINT ${${LEGACY_VAR}} PATH)
      LIST(APPEND ${DIRECTORY_HINT_VAR} ${DIR_HINT})
    ELSEIF (EXISTS ${${LEGACY_VAR}} AND
            IS_DIRECTORY ${${LEGACY_VAR}})
      # User specified a directory hint using the legacy variable, use it.
      LIST(APPEND ${DIRECTORY_HINT_VAR} ${${LEGACY_VAR}})
    ENDIF()
  ENDIF (DEFINED ${LEGACY_VAR})
ENDMACRO(HANDLE_LEGACY_LIBRARY_DEPENDENCY_HINT)

UNSET(CERES_COMPILE_OPTIONS)

# Eigen.
#HANDLE_LEGACY_INCLUDE_DEPENDENCY_HINT(EIGEN_INCLUDE EIGEN_INCLUDE_DIR_HINTS)
LIST(APPEND CMAKE_INCLUDE_PATH ../eigen) # OpenMVG
SET(EIGEN_INCLUDE_DIRS ../eigen) # OpenMVG
FIND_PATH(EIGEN_INCLUDE NAMES Eigen/Core) # OpenMVG
FIND_PACKAGE(Eigen REQUIRED)
IF (EIGEN_FOUND)
  MESSAGE("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
  # Ensure that only MPL2 licensed code is part of the default build.
  MESSAGE("")
  MESSAGE("   ===============================================================")
  IF (EIGENSPARSE)
    LIST(APPEND CERES_COMPILE_OPTIONS CERES_USE_EIGEN_SPARSE)
    MESSAGE("   Enabling the use of Eigen as a sparse linear algebra library ")
    MESSAGE("   for solving the nonlinear least squares problems. Enabling ")
    MESSAGE("   this option will result in an LGPL licensed version of ")
    MESSAGE("   Ceres Solver as the Simplicial Cholesky factorization in Eigen")
    MESSAGE("   is licensed under the LGPL. ")

    IF (EIGEN_VERSION VERSION_LESS 3.2.2)
      MESSAGE("   WARNING:")
      MESSAGE("")
      MESSAGE("   Your version of Eigen is older than version 3.2.2.")
      MESSAGE("   The performance of SPARSE_NORMAL_CHOLESKY and SPARSE_SCHUR")
      MESSAGE("   linear solvers will suffer. ")
    ENDIF (EIGEN_VERSION VERSION_LESS 3.2.2)

  ELSE (EIGENSPARSE)
    MESSAGE("   Disabling the use of Eigen as a sparse linear algebra library.")
    MESSAGE("   This does not affect the covariance estimation algorithm ")
    MESSAGE("   which can still use the EIGEN_SPARSE_QR algorithm.")
    ADD_DEFINITIONS(-DEIGEN_MPL2_ONLY)
  ENDIF (EIGENSPARSE)
    MESSAGE("   ===============================================================")
    MESSAGE("")
ENDIF (EIGEN_FOUND)

# LAPACK (& BLAS).
IF (LAPACK)
  FIND_PACKAGE(LAPACK QUIET)
  IF (LAPACK_FOUND)
    MESSAGE("-- Found LAPACK library: ${LAPACK_LIBRARIES}")
  ELSE (LAPACK_FOUND)
    MESSAGE("-- Did not find LAPACK library, disabling LAPACK support.")
  ENDIF (LAPACK_FOUND)

  FIND_PACKAGE(BLAS QUIET)
  IF (BLAS_FOUND)
    MESSAGE("-- Found BLAS library: ${BLAS_LIBRARIES}")
  ELSE (BLAS_FOUND)
    MESSAGE("-- Did not find BLAS library, disabling LAPACK support.")
  ENDIF (BLAS_FOUND)

  IF (NOT (LAPACK_FOUND AND BLAS_FOUND))
    UPDATE_CACHE_VARIABLE(LAPACK OFF)
    LIST(APPEND CERES_COMPILE_OPTIONS CERES_NO_LAPACK)
  ENDIF (NOT (LAPACK_FOUND AND BLAS_FOUND))
ELSE (LAPACK)
  MESSAGE("-- Building without LAPACK.")
  LIST(APPEND CERES_COMPILE_OPTIONS CERES_NO_LAPACK)
ENDIF (LAPACK)

# SuiteSparse.
IF (SUITESPARSE AND NOT LAPACK)
  # If user has disabled LAPACK, but left SUITESPARSE ON, turn it OFF,
  # LAPACK controls whether Ceres will be linked, directly or indirectly
  # via SuiteSparse to LAPACK.
  MESSAGE("-- Disabling SuiteSparse as use of LAPACK has been disabled, "
    "turn ON LAPACK to enable (optional) building with SuiteSparse.")
  UPDATE_CACHE_VARIABLE(SUITESPARSE OFF)
ENDIF (SUITESPARSE AND NOT LAPACK)
IF (SUITESPARSE)
  # By default, if SuiteSparse and all dependencies are found, Ceres is
  # built with SuiteSparse support.

  # Check for SuiteSparse and dependencies.
  FIND_PACKAGE(SuiteSparse)
  IF (SUITESPARSE_FOUND)
    # On Ubuntu the system install of SuiteSparse (v3.4.0) up to at least
    # Ubuntu 13.10 cannot be used to link shared libraries.
    IF (BUILD_SHARED_LIBS AND
        SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION)
      MESSAGE(FATAL_ERROR "You are attempting to build Ceres as a shared "
        "library on Ubuntu using a system package install of SuiteSparse "
        "3.4.0. This package is broken and does not support the "
        "construction of shared libraries (you can still build Ceres as "
        "a static library).  If you wish to build a shared version of Ceres "
        "you should uninstall the system install of SuiteSparse "
        "(libsuitesparse-dev) and perform a source install of SuiteSparse "
        "(we recommend that you use the latest version), "
        "see: http://homes.cs.washington.edu/~sagarwal"
        "/ceres-solver/dev/building.html for more information.")
    ENDIF (BUILD_SHARED_LIBS AND
      SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION)

    # By default, if all of SuiteSparse's dependencies are found, Ceres is
    # built with SuiteSparse support.
    MESSAGE("-- Found SuiteSparse ${SUITESPARSE_VERSION}, "
            "building with SuiteSparse.")
  ELSE (SUITESPARSE_FOUND)
    # Disable use of SuiteSparse if it cannot be found and continue.
    MESSAGE("-- Did not find all SuiteSparse dependencies, disabling "
      "SuiteSparse support.")
    UPDATE_CACHE_VARIABLE(SUITESPARSE OFF)
    LIST(APPEND CERES_COMPILE_OPTIONS CERES_NO_SUITESPARSE)
  ENDIF (SUITESPARSE_FOUND)
ELSE (SUITESPARSE)
  MESSAGE("-- Building without SuiteSparse.")
  LIST(APPEND CERES_COMPILE_OPTIONS CERES_NO_SUITESPARSE)
ENDIF (SUITESPARSE)

# CXSparse.
IF (CXSPARSE)
  SET(CXSPARSE_INCLUDE_DIRS ../cxsparse/Include) # OpenMVG
  SET(CXSPARSE_FOUND ON)
  # openMVG: cxsparse already exists
  SET(CXSPARSE_LIBRARIES cxsparse)
  # Don't search with REQUIRED as we can continue without CXSparse.
#  FIND_PACKAGE(CXSparse)
#  IF (CXSPARSE_FOUND)
#    # By default, if CXSparse and all dependencies are found, Ceres is
#    # built with CXSparse support.
#    MESSAGE("-- Found CXSparse version: ${CXSPARSE_VERSION}, "
#      "building with CXSparse.")
#  ELSE (CXSPARSE_FOUND)
#    # Disable use of CXSparse if it cannot be found and continue.
#    MESSAGE("-- Did not find CXSparse, Building without CXSparse.")
#    UPDATE_CACHE_VARIABLE(CXSPARSE OFF)
#    LIST(APPEND CERES_COMPILE_OPTIONS CERES_NO_CXSPARSE)
#  ENDIF (CXSPARSE_FOUND)
ELSE (CXSPARSE)
  MESSAGE("-- Building without CXSparse.")
  LIST(APPEND CERES_COMPILE_OPTIONS CERES_NO_CXSPARSE)
  # Mark as advanced (remove from default GUI view) the CXSparse search
  # variables in case user enabled CXSPARSE, FindCXSparse did not find it, so
  # made search variables visible in GUI for user to set, but then user disables
  # CXSPARSE instead of setting them.
  MARK_AS_ADVANCED(FORCE CXSPARSE_INCLUDE_DIR
                         CXSPARSE_LIBRARY)
ENDIF (CXSPARSE)

# GFlags.
IF (GFLAGS)
  HANDLE_LEGACY_INCLUDE_DEPENDENCY_HINT(GFLAGS_INCLUDE GFLAGS_INCLUDE_DIR_HINTS)
  HANDLE_LEGACY_LIBRARY_DEPENDENCY_HINT(GFLAGS_LIB GFLAGS_LIBRARY_DIR_HINTS)

  # Don't search with REQUIRED as we can continue without gflags.
  FIND_PACKAGE(Gflags)
  IF (GFLAGS_FOUND)
    MESSAGE("-- Found Google Flags header in: ${GFLAGS_INCLUDE_DIRS}, "
      "in namespace: ${GFLAGS_NAMESPACE}")
  ELSE (GFLAGS_FOUND)
    MESSAGE("-- Did not find Google Flags (gflags), Building without gflags "
      "- no tests or tools will be built!")
    UPDATE_CACHE_VARIABLE(GFLAGS OFF)
  ENDIF (GFLAGS_FOUND)
ELSE (GFLAGS)
  MESSAGE("-- Google Flags disabled; no tests or tools will be built!")
  # Mark as advanced (remove from default GUI view) the gflags search
  # variables in case user enabled GFLAGS, FindGflags did not find it, so
  # made search variables visible in GUI for user to set, but then user disables
  # GFLAGS instead of setting them.
  MARK_AS_ADVANCED(FORCE GFLAGS_INCLUDE_DIR
                         GFLAGS_LIBRARY)
ENDIF (GFLAGS)

# MiniGLog.
IF (MINIGLOG)
  MESSAGE("-- Compiling minimal glog substitute into Ceres.")
  SET(GLOG_INCLUDE_DIRS internal/ceres/miniglog)
  MESSAGE("-- Using minimal glog substitute (include): ${GLOG_INCLUDE_DIRS}")

  # Mark as advanced (remove from default GUI view) the glog search
  # variables in case user disables MINIGLOG, FindGlog did not find it, so
  # made search variables visible in GUI for user to set, but then user enables
  # MINIGLOG instead of setting them.
  MARK_AS_ADVANCED(FORCE GLOG_INCLUDE_DIR
                         GLOG_LIBRARY)
ELSE (MINIGLOG)
  HANDLE_LEGACY_INCLUDE_DEPENDENCY_HINT(GLOG_INCLUDE GLOG_INCLUDE_DIR_HINTS)
  HANDLE_LEGACY_LIBRARY_DEPENDENCY_HINT(GLOG_LIB GLOG_LIBRARY_DIR_HINTS)

  # Don't search with REQUIRED so that configuration continues if not found and
  # we can output an error messages explaining MINIGLOG option.
  FIND_PACKAGE(Glog)
  IF (GLOG_FOUND)
    MESSAGE("-- Found Google Log header in: ${GLOG_INCLUDE_DIRS}")
  ELSE (GLOG_FOUND)
    MESSAGE(FATAL_ERROR "Can't find Google Log. Please set GLOG_INCLUDE_DIR & "
      "GLOG_LIBRARY or enable MINIGLOG option to use minimal glog "
      "implementation.")
  ENDIF (GLOG_FOUND)
ENDIF (MINIGLOG)

IF (NOT SCHUR_SPECIALIZATIONS)
  LIST(APPEND CERES_COMPILE_OPTIONS CERES_RESTRICT_SCHUR_SPECIALIZATION)
  MESSAGE("-- Disabling Schur specializations (faster compiles)")
ENDIF (NOT SCHUR_SPECIALIZATIONS)

IF (NOT CUSTOM_BLAS)
  LIST(APPEND CERES_COMPILE_OPTIONS CERES_NO_CUSTOM_BLAS)
  MESSAGE("-- Disabling custom blas")
ENDIF (NOT CUSTOM_BLAS)

IF (OPENMP)
  # Clang does not (yet) support OpenMP.
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    UPDATE_CACHE_VARIABLE(OPENMP OFF)
    MESSAGE("-- Compiler is Clang, disabling OpenMP.")
    LIST(APPEND CERES_COMPILE_OPTIONS CERES_NO_THREADS)
  ELSE (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Find quietly s/t as we can continue without OpenMP if it is not found.
    FIND_PACKAGE(OpenMP QUIET)
    IF (OPENMP_FOUND)
      MESSAGE("-- Building with OpenMP.")
      LIST(APPEND CERES_COMPILE_OPTIONS CERES_USE_OPENMP)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      IF (UNIX)
        # At least on Linux, we need pthreads to be enabled for mutex to
        # compile.  This may not work on Windows or Android.
        FIND_PACKAGE(Threads REQUIRED)
        LIST(APPEND CERES_COMPILE_OPTIONS CERES_HAVE_PTHREAD)
        LIST(APPEND CERES_COMPILE_OPTIONS CERES_HAVE_RWLOCK)
      ENDIF (UNIX)
    ELSE (OPENMP_FOUND)
      MESSAGE("-- Failed to find OpenMP, disabling.")
      UPDATE_CACHE_VARIABLE(OPENMP OFF)
      LIST(APPEND CERES_COMPILE_OPTIONS CERES_NO_THREADS)
    ENDIF (OPENMP_FOUND)
  ENDIF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
ELSE (OPENMP)
  MESSAGE("-- Building without OpenMP (disabling multithreading).")
  LIST(APPEND CERES_COMPILE_OPTIONS CERES_NO_THREADS)
ENDIF (OPENMP)

INCLUDE(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(unordered_map HAVE_STD_UNORDERED_MAP_HEADER)
IF (HAVE_STD_UNORDERED_MAP_HEADER)
  # Finding the unordered_map header doesn't mean that unordered_map
  # is in std namespace.
  #
  # In particular, MSVC 2008 has unordered_map declared in std::tr1.
  # In order to support this, we do an extra check to see which
  # namespace should be used.
  INCLUDE(CheckCXXSourceCompiles)
  CHECK_CXX_SOURCE_COMPILES("#include <unordered_map>
                             int main() {
                               std::unordered_map<int, int> map;
                               return 0;
                             }"
                            HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
  IF (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
    LIST(APPEND CERES_COMPILE_OPTIONS CERES_STD_UNORDERED_MAP)
    MESSAGE("-- Found unordered_map/set in std namespace.")
  ELSE (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
    CHECK_CXX_SOURCE_COMPILES("#include <unordered_map>
                               int main() {
                                 std::tr1::unordered_map<int, int> map;
                                 return 0;
                               }"
                              HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
    IF (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
      LIST(APPEND CERES_COMPILE_OPTIONS CERES_STD_UNORDERED_MAP_IN_TR1_NAMESPACE)
      MESSAGE("-- Found unordered_map/set in std::tr1 namespace.")
    ELSE (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
      MESSAGE("-- Found <unordered_map> but cannot find either std::unordered_map "
              "or std::tr1::unordered_map.")
      MESSAGE("-- Replacing unordered_map/set with map/set (warning: slower!)")
      LIST(APPEND CERES_COMPILE_OPTIONS CERES_NO_UNORDERED_MAP)
    ENDIF (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
  ENDIF (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
ELSE (HAVE_STD_UNORDERED_MAP_HEADER)
  CHECK_INCLUDE_FILE_CXX("tr1/unordered_map" HAVE_TR1_UNORDERED_MAP_HEADER)
  IF (HAVE_TR1_UNORDERED_MAP_HEADER)
    LIST(APPEND CERES_COMPILE_OPTIONS CERES_TR1_UNORDERED_MAP)
    MESSAGE("-- Found tr1/unordered_map/set in std::tr1 namespace.")
  ELSE (HAVE_TR1_UNORDERED_MAP_HEADE)
    MESSAGE("-- Unable to find <unordered_map> or <tr1/unordered_map>. ")
    MESSAGE("-- Replacing unordered_map/set with map/set (warning: slower!)")
    LIST(APPEND CERES_COMPILE_OPTIONS CERES_NO_UNORDERED_MAP)
  ENDIF (HAVE_TR1_UNORDERED_MAP_HEADER)
ENDIF (HAVE_STD_UNORDERED_MAP_HEADER)

INCLUDE(FindSharedPtr)
FIND_SHARED_PTR()
IF (SHARED_PTR_FOUND)
  IF (SHARED_PTR_TR1_MEMORY_HEADER)
    LIST(APPEND CERES_COMPILE_OPTIONS CERES_TR1_MEMORY_HEADER)
  ENDIF (SHARED_PTR_TR1_MEMORY_HEADER)
  IF (SHARED_PTR_TR1_NAMESPACE)
    LIST(APPEND CERES_COMPILE_OPTIONS CERES_TR1_SHARED_PTR)
  ENDIF (SHARED_PTR_TR1_NAMESPACE)
ELSE (SHARED_PTR_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find shared_ptr.")
ENDIF (SHARED_PTR_FOUND)

INCLUDE_DIRECTORIES(
  include
  internal
  internal/ceres
  ${GLOG_INCLUDE_DIRS})
# Eigen SparseQR generates various compiler warnings related to unused and
# uninitialised local variables, which prevents Ceres compilation as we use
# -Werror.  To avoid having to individually suppress these warnings around
# the #include statments for Eigen headers across all GCC/Clang versions, we
# tell CMake to treat Eigen headers as system headers.  This results in all
# compiler warnings from them being suppressed.
#
# Note that this is *not* propagated to clients, ie CERES_INCLUDE_DIRS
# used by clients after find_package(Ceres) does not identify Eigen as
# as system headers.
INCLUDE_DIRECTORIES(SYSTEM ${EIGEN_INCLUDE_DIRS})

IF (SUITESPARSE)
  INCLUDE_DIRECTORIES(${SUITESPARSE_INCLUDE_DIRS})
ENDIF (SUITESPARSE)

IF (CXSPARSE)
  INCLUDE_DIRECTORIES(${CXSPARSE_INCLUDE_DIRS})
ENDIF (CXSPARSE)

IF (GFLAGS)
  INCLUDE_DIRECTORIES(${GFLAGS_INCLUDE_DIRS})
ENDIF (GFLAGS)

IF (BUILD_SHARED_LIBS)
  MESSAGE("-- Building Ceres as a shared library.")
  # The CERES_BUILDING_SHARED_LIBRARY compile definition is NOT stored in
  # CERES_COMPILE_OPTIONS as it must only be defined when Ceres is compiled
  # not when it is used as it controls the CERES_EXPORT macro which provides
  # dllimport/export support in MSVC.
  ADD_DEFINITIONS(-DCERES_BUILDING_SHARED_LIBRARY)
  LIST(APPEND CERES_COMPILE_OPTIONS CERES_USING_SHARED_LIBRARY)
ELSE (BUILD_SHARED_LIBS)
  MESSAGE("-- Building Ceres as a static library.")
ENDIF (BUILD_SHARED_LIBS)

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
      IF (NOT GCC_VERSION VERSION_LESS 4.2)
        SET (CERES_CXX_FLAGS "${CERES_CXX_FLAGS} -march=native -mtune=native")
      ENDIF (NOT GCC_VERSION VERSION_LESS 4.2)
    ENDIF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    # Mac OS X
    IF (CMAKE_SYSTEM_NAME MATCHES "Darwin")
      SET (CERES_CXX_FLAGS "${CERES_CXX_FLAGS} -msse3")
      # Use of -fast only applicable for Apple's GCC
      # Assume this is being used if GCC version < 4.3 on OSX
      EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER}
        ARGS ${CMAKE_CXX_COMPILER_ARG1} -dumpversion
        OUTPUT_VARIABLE GCC_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      IF (GCC_VERSION VERSION_LESS 4.3)
        SET (CERES_CXX_FLAGS "${CERES_CXX_FLAGS} -fast")
      ENDIF (GCC_VERSION VERSION_LESS 4.3)
    ENDIF (CMAKE_SYSTEM_NAME MATCHES "Darwin")
  ENDIF (CMAKE_COMPILER_IS_GNUCXX)
#  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#    # Use of -flto requires use of gold linker & LLVM-gold plugin, which might
#    # well not be present / in use and without which files will compile, but
#    # not link ('file not recognized') so explicitly check for support
#    INCLUDE(CheckCXXCompilerFlag)
#    CHECK_CXX_COMPILER_FLAG("-flto" HAVE_LTO_SUPPORT)
#    IF (HAVE_LTO_SUPPORT)
#      MESSAGE(STATUS "Enabling link-time optimization (-flto)")
#      SET(CERES_CXX_FLAGS "${CERES_CXX_FLAGS} -flto")
#    ELSE ()
#      MESSAGE(STATUS "Compiler/linker does not support link-time optimization (-flto), disabling.")
#    ENDIF (HAVE_LTO_SUPPORT)
#  ENDIF ()
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

  # Update the C/CXX flags for MSVC to use either the static or shared
  # C-Run Time (CRT) library based on the user option: MSVC_USE_STATIC_CRT.
  LIST(APPEND C_CXX_FLAGS
    CMAKE_CXX_FLAGS
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS_RELWITHDEBINFO)

  FOREACH(FLAG_VAR ${C_CXX_FLAGS})
    IF (MSVC_USE_STATIC_CRT)
      # Use static CRT.
      IF (${FLAG_VAR} MATCHES "/MD")
        STRING(REGEX REPLACE "/MD" "/MT" ${FLAG_VAR} "${${FLAG_VAR}}")
      ENDIF (${FLAG_VAR} MATCHES "/MD")
    ELSE (MSVC_USE_STATIC_CRT)
      # Use shared, not static, CRT.
      IF (${FLAG_VAR} MATCHES "/MT")
        STRING(REGEX REPLACE "/MT" "/MD" ${FLAG_VAR} "${${FLAG_VAR}}")
      ENDIF (${FLAG_VAR} MATCHES "/MT")
    ENDIF (MSVC_USE_STATIC_CRT)
  ENDFOREACH()

  # Tuple sizes of 10 are used by Gtest.
  ADD_DEFINITIONS("-D_VARIADIC_MAX=10")
ENDIF (MSVC)

IF (UNIX)
  # GCC is not strict enough by default, so enable most of the warnings.
  SET(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Werror=all -Werror=extra -Wno-unknown-pragmas -Wno-sign-compare -Wno-unused-parameter -Wno-missing-field-initializers")
ENDIF (UNIX)

# Use a larger inlining threshold for Clang, since it hobbles Eigen,
# resulting in an unreasonably slow version of the blas routines. The
# -Qunused-arguments is needed because CMake passes the inline
# threshold to the linker and clang complains about it and dies.
IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Qunused-arguments -mllvm -inline-threshold=600")
  # Older versions of Clang (<= 2.9) do not support the 'return-type-c-linkage'
  # option, so check for its presence before adding it to the default flags set.
  INCLUDE(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-Wno-return-type-c-linkage"
                          HAVE_RETURN_TYPE_C_LINKAGE)
  IF (HAVE_RETURN_TYPE_C_LINKAGE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-return-type-c-linkage")
  ENDIF(HAVE_RETURN_TYPE_C_LINKAGE)
ENDIF ()

# Xcode 4.5.x used Clang 4.1 (Apple version), this has a bug that prevents
# compilation of Ceres.
IF (APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER}
    ARGS ${CMAKE_CXX_COMPILER_ARG1} -dumpversion
    OUTPUT_VARIABLE CLANG_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  # Use version > 4.0 & < 4.2 to catch all 4.1(.x) versions.
  IF (CLANG_VERSION VERSION_GREATER 4.0 AND
      CLANG_VERSION VERSION_LESS 4.2)
    MESSAGE(FATAL_ERROR "You are attempting to build Ceres on OS X using Xcode "
      "4.5.x (Clang version: ${CLANG_VERSION}). This version of Clang has a "
      "bug that prevents compilation of Ceres, please update to "
      "Xcode >= 4.6.3.")
  ENDIF (CLANG_VERSION VERSION_GREATER 4.0 AND
    CLANG_VERSION VERSION_LESS 4.2)
ENDIF (APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

# Configure the Ceres config.h compile options header using the current
# compile options and put the configured header into the Ceres build
# directory.  Note that the ceres/internal subdir in <build>/config where
# the configured config.h is placed is important, because Ceres will be
# built against this configured header, it needs to have the same relative
# include path as it would if it were in the source tree (or installed).
LIST(REMOVE_DUPLICATES CERES_COMPILE_OPTIONS)
INCLUDE(CreateCeresConfig)
CREATE_CERES_CONFIG("${CERES_COMPILE_OPTIONS}"
  ${CMAKE_CURRENT_BINARY_DIR}/config/ceres/internal)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/config)

ADD_SUBDIRECTORY(internal/ceres)

IF (BUILD_DOCUMENTATION)
  MESSAGE("-- Documentation building is enabled")

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
