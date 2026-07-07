# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2023 Google Inc. All rights reserved.
# http://ceres-solver.org/
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
# Author: alexs.mac@gmail.com (Alex Stewart)
#

#[=======================================================================[.rst:
FindSuiteSparse
===============

Module for locating SuiteSparse libraries and its dependencies.

This module defines the following variables:

``SuiteSparse_FOUND``
   ``TRUE`` iff SuiteSparse and all dependencies have been found.

``SuiteSparse_VERSION``
   Extracted from ``SuiteSparse_config.h`` (>= v4).

``SuiteSparse_VERSION_MAJOR``
    Equal to 4 if ``SuiteSparse_VERSION`` = 4.2.1

``SuiteSparse_VERSION_MINOR``
    Equal to 2 if ``SuiteSparse_VERSION`` = 4.2.1

``SuiteSparse_VERSION_PATCH``
    Equal to 1 if ``SuiteSparse_VERSION`` = 4.2.1

The following variables control the behaviour of this module:

``SuiteSparse_NO_CMAKE``
  Do not attempt to use the native SuiteSparse CMake package configuration.


Targets
-------

The following targets define the SuiteSparse components searched for.

``SuiteSparse::AMD``
    Symmetric Approximate Minimum Degree (AMD)

``SuiteSparse::CAMD``
    Constrained Approximate Minimum Degree (CAMD)

``SuiteSparse::COLAMD``
    Column Approximate Minimum Degree (COLAMD)

``SuiteSparse::CCOLAMD``
    Constrained Column Approximate Minimum Degree (CCOLAMD)

``SuiteSparse::CHOLMOD``
    Sparse Supernodal Cholesky Factorization and Update/Downdate (CHOLMOD)

``SuiteSparse::Partition``
    CHOLMOD with METIS support

``SuiteSparse::SPQR``
    Multifrontal Sparse QR (SuiteSparseQR)

``SuiteSparse::Config``
    Common configuration for all but CSparse (SuiteSparse version >= 4).

Optional SuiteSparse dependencies:

``METIS::METIS``
    Serial Graph Partitioning and Fill-reducing Matrix Ordering (METIS)
]=======================================================================]

if (NOT SuiteSparse_NO_CMAKE)
  find_package (SuiteSparse NO_MODULE QUIET)
endif (NOT SuiteSparse_NO_CMAKE)

if (SuiteSparse_FOUND)
  return ()
endif (SuiteSparse_FOUND)

# Push CMP0057 to enable support for IN_LIST, when cmake_minimum_required is
# set to <3.3.
cmake_policy (PUSH)
cmake_policy (SET CMP0057 NEW)

if (NOT SuiteSparse_FIND_COMPONENTS)
  set (SuiteSparse_FIND_COMPONENTS
    AMD
    CAMD
    CCOLAMD
    CHOLMOD
    COLAMD
    SPQR
  )

  foreach (component IN LISTS SuiteSparse_FIND_COMPONENTS)
    set (SuiteSparse_FIND_REQUIRED_${component} TRUE)
  endforeach (component IN LISTS SuiteSparse_FIND_COMPONENTS)
endif (NOT SuiteSparse_FIND_COMPONENTS)

# Assume SuiteSparse was found and set it to false only if third-party
# dependencies could not be located. SuiteSparse components are handled by
# FindPackageHandleStandardArgs HANDLE_COMPONENTS option.
set (SuiteSparse_FOUND TRUE)

include (CheckLibraryExists)
include (CheckSymbolExists)
include (CMakePushCheckState)

# Config is a base component and thus always required
set (SuiteSparse_IMPLICIT_COMPONENTS Config)

# CHOLMOD depends on AMD, CAMD, CCOLAMD, and COLAMD.
if (CHOLMOD IN_LIST SuiteSparse_FIND_COMPONENTS)
  list (APPEND SuiteSparse_IMPLICIT_COMPONENTS AMD CAMD CCOLAMD COLAMD)
endif (CHOLMOD IN_LIST SuiteSparse_FIND_COMPONENTS)

# SPQR depends on CHOLMOD.
if (SPQR IN_LIST SuiteSparse_FIND_COMPONENTS)
  list (APPEND SuiteSparse_IMPLICIT_COMPONENTS CHOLMOD)
endif (SPQR IN_LIST SuiteSparse_FIND_COMPONENTS)

# Implicit components are always required
foreach (component IN LISTS SuiteSparse_IMPLICIT_COMPONENTS)
  set (SuiteSparse_FIND_REQUIRED_${component} TRUE)
endforeach (component IN LISTS SuiteSparse_IMPLICIT_COMPONENTS)

list (APPEND SuiteSparse_FIND_COMPONENTS ${SuiteSparse_IMPLICIT_COMPONENTS})

# Do not list components multiple times.
list (REMOVE_DUPLICATES SuiteSparse_FIND_COMPONENTS)

# Reset CALLERS_CMAKE_FIND_LIBRARY_PREFIXES to its value when
# FindSuiteSparse was invoked.
macro(SuiteSparse_RESET_FIND_LIBRARY_PREFIX)
  if (MSVC)
    set(CMAKE_FIND_LIBRARY_PREFIXES "${CALLERS_CMAKE_FIND_LIBRARY_PREFIXES}")
  endif (MSVC)
endmacro(SuiteSparse_RESET_FIND_LIBRARY_PREFIX)

# Called if we failed to find SuiteSparse or any of it's required dependencies,
# unsets all public (designed to be used externally) variables and reports
# error message at priority depending upon [REQUIRED/QUIET/<NONE>] argument.
macro(SuiteSparse_REPORT_NOT_FOUND REASON_MSG)
  # Will be set to FALSE by find_package_handle_standard_args
  unset (SuiteSparse_FOUND)

  # Do NOT unset SuiteSparse_REQUIRED_VARS here, as it is used by
  # FindPackageHandleStandardArgs() to generate the automatic error message on
  # failure which highlights which components are missing.

  suitesparse_reset_find_library_prefix()

  # Note <package>_FIND_[REQUIRED/QUIETLY] variables defined by FindPackage()
  # use the camelcase library name, not uppercase.
  if (SuiteSparse_FIND_QUIETLY)
    message(STATUS "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
  elseif (SuiteSparse_FIND_REQUIRED)
    message(FATAL_ERROR "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
  else()
    # Neither QUIETLY nor REQUIRED, use no priority which emits a message
    # but continues configuration and allows generation.
    message("-- Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
  endif (SuiteSparse_FIND_QUIETLY)

  # Do not call return(), s/t we keep processing if not called with REQUIRED
  # and report all missing components, rather than bailing after failing to find
  # the first.
endmacro(SuiteSparse_REPORT_NOT_FOUND)

# Handle possible presence of lib prefix for libraries on MSVC, see
# also SuiteSparse_RESET_FIND_LIBRARY_PREFIX().
if (MSVC)
  # Preserve the caller's original values for CMAKE_FIND_LIBRARY_PREFIXES
  # s/t we can set it back before returning.
  set(CALLERS_CMAKE_FIND_LIBRARY_PREFIXES "${CMAKE_FIND_LIBRARY_PREFIXES}")
  # The empty string in this list is important, it represents the case when
  # the libraries have no prefix (shared libraries / DLLs).
  set(CMAKE_FIND_LIBRARY_PREFIXES "lib" "" "${CMAKE_FIND_LIBRARY_PREFIXES}")
endif (MSVC)

# Additional suffixes to try appending to each search path.
list(APPEND SuiteSparse_CHECK_PATH_SUFFIXES
  suitesparse) # Windows/Ubuntu

# Wrappers to find_path/library that pass the SuiteSparse search hints/paths.
#
# suitesparse_find_component(<component> [FILES name1 [name2 ...]]
#                                        [LIBRARIES name1 [name2 ...]])
macro(suitesparse_find_component COMPONENT)
  include(CMakeParseArguments)
  set(MULTI_VALUE_ARGS FILES LIBRARIES)
  cmake_parse_arguments(SuiteSparse_FIND_COMPONENT_${COMPONENT}
    "" "" "${MULTI_VALUE_ARGS}" ${ARGN})

  set(SuiteSparse_${COMPONENT}_FOUND TRUE)
  if (SuiteSparse_FIND_COMPONENT_${COMPONENT}_FILES)
    find_path(SuiteSparse_${COMPONENT}_INCLUDE_DIR
      NAMES ${SuiteSparse_FIND_COMPONENT_${COMPONENT}_FILES}
      PATH_SUFFIXES ${SuiteSparse_CHECK_PATH_SUFFIXES})
    if (SuiteSparse_${COMPONENT}_INCLUDE_DIR)
      message(STATUS "Found ${COMPONENT} headers in: "
        "${SuiteSparse_${COMPONENT}_INCLUDE_DIR}")
      mark_as_advanced(SuiteSparse_${COMPONENT}_INCLUDE_DIR)
    else()
      # Specified headers not found.
      set(SuiteSparse_${COMPONENT}_FOUND FALSE)
      if (SuiteSparse_FIND_REQUIRED_${COMPONENT})
        suitesparse_report_not_found(
          "Did not find ${COMPONENT} header (required SuiteSparse component).")
      else()
        message(STATUS "Did not find ${COMPONENT} header (optional "
          "SuiteSparse component).")
        # Hide optional vars from CMake GUI even if not found.
        mark_as_advanced(SuiteSparse_${COMPONENT}_INCLUDE_DIR)
      endif()
    endif()
  endif()

  if (SuiteSparse_FIND_COMPONENT_${COMPONENT}_LIBRARIES)
    find_library(SuiteSparse_${COMPONENT}_LIBRARY
      NAMES ${SuiteSparse_FIND_COMPONENT_${COMPONENT}_LIBRARIES}
      PATH_SUFFIXES ${SuiteSparse_CHECK_PATH_SUFFIXES})
    if (SuiteSparse_${COMPONENT}_LIBRARY)
      message(STATUS "Found ${COMPONENT} library: ${SuiteSparse_${COMPONENT}_LIBRARY}")
      mark_as_advanced(SuiteSparse_${COMPONENT}_LIBRARY)
    else ()
      # Specified libraries not found.
      set(SuiteSparse_${COMPONENT}_FOUND FALSE)
      if (SuiteSparse_FIND_REQUIRED_${COMPONENT})
        suitesparse_report_not_found(
          "Did not find ${COMPONENT} library (required SuiteSparse component).")
      else()
        message(STATUS "Did not find ${COMPONENT} library (optional SuiteSparse "
          "dependency)")
        # Hide optional vars from CMake GUI even if not found.
        mark_as_advanced(SuiteSparse_${COMPONENT}_LIBRARY)
      endif()
    endif()
  endif()

  # A component can be optional (given to OPTIONAL_COMPONENTS). However, if the
  # component is implicit (must be always present, such as the Config component)
  # assume it be required as well.
  if (SuiteSparse_FIND_REQUIRED_${COMPONENT})
    list (APPEND SuiteSparse_REQUIRED_VARS SuiteSparse_${COMPONENT}_INCLUDE_DIR)
    list (APPEND SuiteSparse_REQUIRED_VARS SuiteSparse_${COMPONENT}_LIBRARY)
  endif (SuiteSparse_FIND_REQUIRED_${COMPONENT})

  # Define the target only if the include directory and the library were found
  if (SuiteSparse_${COMPONENT}_INCLUDE_DIR AND SuiteSparse_${COMPONENT}_LIBRARY)
    if (NOT TARGET SuiteSparse::${COMPONENT})
      add_library(SuiteSparse::${COMPONENT} IMPORTED UNKNOWN)
    endif (NOT TARGET SuiteSparse::${COMPONENT})

    set_property(TARGET SuiteSparse::${COMPONENT} PROPERTY
      INTERFACE_INCLUDE_DIRECTORIES ${SuiteSparse_${COMPONENT}_INCLUDE_DIR})
    set_property(TARGET SuiteSparse::${COMPONENT} PROPERTY
      IMPORTED_LOCATION ${SuiteSparse_${COMPONENT}_LIBRARY})
  endif (SuiteSparse_${COMPONENT}_INCLUDE_DIR AND SuiteSparse_${COMPONENT}_LIBRARY)
endmacro()

# Given the number of components of SuiteSparse, and to ensure that the
# automatic failure message generated by FindPackageHandleStandardArgs()
# when not all required components are found is helpful, we maintain a list
# of all variables that must be defined for SuiteSparse to be considered found.
unset(SuiteSparse_REQUIRED_VARS)

# BLAS.
find_package(BLAS QUIET)
if (NOT BLAS_FOUND)
  suitesparse_report_not_found(
    "Did not find BLAS library (required for SuiteSparse).")
endif (NOT BLAS_FOUND)

# LAPACK.
find_package(LAPACK QUIET)
if (NOT LAPACK_FOUND)
  suitesparse_report_not_found(
    "Did not find LAPACK library (required for SuiteSparse).")
endif (NOT LAPACK_FOUND)

foreach (component IN LISTS SuiteSparse_FIND_COMPONENTS)
  if (component STREQUAL Partition)
    # Partition is a meta component that neither provides additional headers nor
    # a separate library. It is strictly part of CHOLMOD.
    continue ()
  endif (component STREQUAL Partition)
  string (TOLOWER ${component} component_library)

  if (component STREQUAL "Config")
    set (component_header SuiteSparse_config.h)
    set (component_library suitesparseconfig)
  elseif (component STREQUAL "SPQR")
    set (component_header SuiteSparseQR.hpp)
  else (component STREQUAL "SPQR")
    set (component_header ${component_library}.h)
  endif (component STREQUAL "Config")

  suitesparse_find_component(${component}
    FILES ${component_header}
    LIBRARIES ${component_library})
endforeach (component IN LISTS SuiteSparse_FIND_COMPONENTS)

if (TARGET SuiteSparse::SPQR)
  # SuiteSparseQR may be compiled with Intel Threading Building Blocks,
  # we assume that if TBB is installed, SuiteSparseQR was compiled with
  # support for it, this will do no harm if it wasn't.
  find_package(TBB QUIET)
  if (TBB_FOUND)
    message(STATUS "Found Intel Thread Building Blocks (TBB) library "
      "(${TBB_VERSION_MAJOR}.${TBB_VERSION_MINOR} / ${TBB_INTERFACE_VERSION}) "
      "include location: ${TBB_INCLUDE_DIRS}. Assuming SuiteSparseQR was "
      "compiled with TBB.")
    # Add the TBB libraries to the SuiteSparseQR libraries (the only
    # libraries to optionally depend on TBB).
    if (TARGET TBB::tbb)
      # Native TBB package configuration provides an imported target. Use it if
      # available.
      set_property (TARGET SuiteSparse::SPQR APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES TBB::tbb)
    else (TARGET TBB::tbb)
      set_property (TARGET SuiteSparse::SPQR APPEND PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES ${TBB_INCLUDE_DIRS})
      set_property (TARGET SuiteSparse::SPQR APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES ${TBB_LIBRARIES})
    endif (TARGET TBB::tbb)
  else (TBB_FOUND)
    message(STATUS "Did not find Intel TBB library, assuming SuiteSparseQR was "
      "not compiled with TBB.")
  endif (TBB_FOUND)
endif (TARGET SuiteSparse::SPQR)

check_library_exists(rt shm_open "" HAVE_LIBRT)

if (TARGET SuiteSparse::Config)
  # SuiteSparse_config (SuiteSparse version >= 4) requires librt library for
  # timing by default when compiled on Linux or Unix, but not on OSX (which
  # does not have librt).
  if (HAVE_LIBRT)
    message(STATUS "Adding librt to "
      "SuiteSparse_config libraries (required on Linux & Unix [not OSX] if "
      "SuiteSparse is compiled with timing).")
    set_property (TARGET SuiteSparse::Config APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES $<LINK_ONLY:rt>)
  else (HAVE_LIBRT)
    message(STATUS "Could not find librt, but found SuiteSparse_config, "
      "assuming that SuiteSparse was compiled without timing.")
  endif (HAVE_LIBRT)

  # Add BLAS and LAPACK as dependencies of SuiteSparse::Config for convenience
  # given that all components depend on it.
  if (BLAS_FOUND)
    if (TARGET BLAS::BLAS)
      set_property (TARGET SuiteSparse::Config APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES $<LINK_ONLY:BLAS::BLAS>)
    else (TARGET BLAS::BLAS)
      set_property (TARGET SuiteSparse::Config APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES ${BLAS_LIBRARIES})
    endif (TARGET BLAS::BLAS)
  endif (BLAS_FOUND)

  if (LAPACK_FOUND)
    if (TARGET LAPACK::LAPACK)
      set_property (TARGET SuiteSparse::Config APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES $<LINK_ONLY:LAPACK::LAPACK>)
    else (TARGET LAPACK::LAPACK)
      set_property (TARGET SuiteSparse::Config APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})
    endif (TARGET LAPACK::LAPACK)
  endif (LAPACK_FOUND)

  # SuiteSparse version >= 4.
  set(SuiteSparse_VERSION_FILE
    ${SuiteSparse_Config_INCLUDE_DIR}/SuiteSparse_config.h)
  if (NOT EXISTS ${SuiteSparse_VERSION_FILE})
    suitesparse_report_not_found(
      "Could not find file: ${SuiteSparse_VERSION_FILE} containing version "
      "information for >= v4 SuiteSparse installs, but SuiteSparse_config was "
      "found (only present in >= v4 installs).")
  else (NOT EXISTS ${SuiteSparse_VERSION_FILE})
    file(READ ${SuiteSparse_VERSION_FILE} Config_CONTENTS)

    string(REGEX MATCH "#define SUITESPARSE_MAIN_VERSION[ \t]+([0-9]+)"
      SuiteSparse_VERSION_LINE "${Config_CONTENTS}")
    set (SuiteSparse_VERSION_MAJOR ${CMAKE_MATCH_1})

    string(REGEX MATCH "#define SUITESPARSE_SUB_VERSION[ \t]+([0-9]+)"
      SuiteSparse_VERSION_LINE "${Config_CONTENTS}")
    set (SuiteSparse_VERSION_MINOR ${CMAKE_MATCH_1})

    string(REGEX MATCH "#define SUITESPARSE_SUBSUB_VERSION[ \t]+([0-9]+)"
      SuiteSparse_VERSION_LINE "${Config_CONTENTS}")
    set (SuiteSparse_VERSION_PATCH ${CMAKE_MATCH_1})

    unset (SuiteSparse_VERSION_LINE)

    # This is on a single line s/t CMake does not interpret it as a list of
    # elements and insert ';' separators which would result in 4.;2.;1 nonsense.
    set(SuiteSparse_VERSION
      "${SuiteSparse_VERSION_MAJOR}.${SuiteSparse_VERSION_MINOR}.${SuiteSparse_VERSION_PATCH}")

    if (SuiteSparse_VERSION MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+")
      set(SuiteSparse_VERSION_COMPONENTS 3)
    else (SuiteSparse_VERSION MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+")
      message (WARNING "Could not parse SuiteSparse_config.h: SuiteSparse "
        "version will not be available")

      unset (SuiteSparse_VERSION)
      unset (SuiteSparse_VERSION_MAJOR)
      unset (SuiteSparse_VERSION_MINOR)
      unset (SuiteSparse_VERSION_PATCH)
    endif (SuiteSparse_VERSION MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+")
  endif (NOT EXISTS ${SuiteSparse_VERSION_FILE})
endif (TARGET SuiteSparse::Config)

# CHOLMOD requires AMD CAMD CCOLAMD COLAMD
if (TARGET SuiteSparse::CHOLMOD)
  foreach (component IN ITEMS AMD CAMD CCOLAMD COLAMD)
    if (TARGET SuiteSparse::${component})
      set_property (TARGET SuiteSparse::CHOLMOD APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES SuiteSparse::${component})
    else (TARGET SuiteSparse::${component})
      # Consider CHOLMOD not found if COLAMD cannot be found
      set (SuiteSparse_CHOLMOD_FOUND FALSE)
    endif (TARGET SuiteSparse::${component})
  endforeach (component IN ITEMS AMD CAMD CCOLAMD COLAMD)
endif (TARGET SuiteSparse::CHOLMOD)

# SPQR requires CHOLMOD
if (TARGET SuiteSparse::SPQR)
  if (TARGET SuiteSparse::CHOLMOD)
    set_property (TARGET SuiteSparse::SPQR APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES SuiteSparse::CHOLMOD)
  else (TARGET SuiteSparse::CHOLMOD)
    # Consider SPQR not found if CHOLMOD cannot be found
    set (SuiteSparse_SQPR_FOUND FALSE)
  endif (TARGET SuiteSparse::CHOLMOD)
endif (TARGET SuiteSparse::SPQR)

# Add SuiteSparse::Config as dependency to all components
if (TARGET SuiteSparse::Config)
  foreach (component IN LISTS SuiteSparse_FIND_COMPONENTS)
    if (component STREQUAL Config)
      continue ()
    endif (component STREQUAL Config)

    if (TARGET SuiteSparse::${component})
      set_property (TARGET SuiteSparse::${component} APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES SuiteSparse::Config)
    endif (TARGET SuiteSparse::${component})
  endforeach (component IN LISTS SuiteSparse_FIND_COMPONENTS)
endif (TARGET SuiteSparse::Config)

# Check whether CHOLMOD was compiled with METIS support. The check can be
# performed only after the main components have been set up.
if (TARGET SuiteSparse::CHOLMOD)
  # NOTE If SuiteSparse was compiled as a static library we'll need to link
  # against METIS already during the check. Otherwise, the check can fail due to
  # undefined references even though SuiteSparse was compiled with METIS.
  find_package (METIS)

  if (TARGET METIS::METIS)
    cmake_push_check_state (RESET)
    set (CMAKE_REQUIRED_LIBRARIES SuiteSparse::CHOLMOD METIS::METIS)
    check_symbol_exists (cholmod_metis cholmod.h SuiteSparse_CHOLMOD_USES_METIS)
    cmake_pop_check_state ()

    if (SuiteSparse_CHOLMOD_USES_METIS)
      set_property (TARGET SuiteSparse::CHOLMOD APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES $<LINK_ONLY:METIS::METIS>)

      # Provide the SuiteSparse::Partition component whose availability indicates
      # that CHOLMOD was compiled with the Partition module.
      if (NOT TARGET SuiteSparse::Partition)
        add_library (SuiteSparse::Partition IMPORTED INTERFACE)
      endif (NOT TARGET SuiteSparse::Partition)

      set_property (TARGET SuiteSparse::Partition APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES SuiteSparse::CHOLMOD)
    endif (SuiteSparse_CHOLMOD_USES_METIS)
  endif (TARGET METIS::METIS)
endif (TARGET SuiteSparse::CHOLMOD)

# We do not use suitesparse_find_component to find Partition and therefore must
# handle the availability in an extra step.
if (TARGET SuiteSparse::Partition)
  set (SuiteSparse_Partition_FOUND TRUE)
else (TARGET SuiteSparse::Partition)
  set (SuiteSparse_Partition_FOUND FALSE)
endif (TARGET SuiteSparse::Partition)

suitesparse_reset_find_library_prefix()

# Handle REQUIRED and QUIET arguments to FIND_PACKAGE
include(FindPackageHandleStandardArgs)
if (SuiteSparse_FOUND)
  find_package_handle_standard_args(SuiteSparse
    REQUIRED_VARS ${SuiteSparse_REQUIRED_VARS}
    VERSION_VAR SuiteSparse_VERSION
    FAIL_MESSAGE "Failed to find some/all required components of SuiteSparse."
    HANDLE_COMPONENTS)
else (SuiteSparse_FOUND)
  # Do not pass VERSION_VAR to FindPackageHandleStandardArgs() if we failed to
  # find SuiteSparse to avoid a confusing autogenerated failure message
  # that states 'not found (missing: FOO) (found version: x.y.z)'.
  find_package_handle_standard_args(SuiteSparse
    REQUIRED_VARS ${SuiteSparse_REQUIRED_VARS}
    FAIL_MESSAGE "Failed to find some/all required components of SuiteSparse."
    HANDLE_COMPONENTS)
endif (SuiteSparse_FOUND)

# Pop CMP0057.
cmake_policy (POP)
