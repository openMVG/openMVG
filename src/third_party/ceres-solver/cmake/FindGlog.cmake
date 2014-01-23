# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2013 Google Inc. All rights reserved.
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
# Author: alexs.mac@gmail.com (Alex Stewart)
#

# FindGlog.cmake - Find Google glog logging library.
#
# This module defines the following variables:
#
# GLOG_FOUND: TRUE iff glog is found.
# GLOG_INCLUDE_DIRS: Include directories for glog.
# GLOG_LIBRARIES: Libraries required to link glog.
#
# The following variables control the behaviour of this module:
#
# GLOG_INCLUDE_DIRS_HINTS: List of additional directories in which to
#                          search for glog includes, e.g: /timbuktu/include.
# GLOG_LIBRARY_DIRS_HINTS: List of additional directories in which to
#                          search for glog libraries, e.g: /timbuktu/lib.
#
# The following variables are also defined by this module, but in line with
# CMake recommended FindPackage() module style should NOT be referenced directly
# by callers (use the plural variables detailed above instead).  These variables
# do however affect the behaviour of the module via FIND_[PATH/LIBRARY]() which
# are NOT re-called (i.e. search for library is not repeated) if these variables
# are set with valid values _in the CMake cache_. This means that if these
# variables are set directly in the cache, either by the user in the CMake GUI,
# or by the user passing -DVAR=VALUE directives to CMake when called (which
# explicitly defines a cache variable), then they will be used verbatim,
# bypassing the HINTS variables and other hard-coded search locations.
#
# GLOG_INCLUDE_DIR: Include directory for glog, not including the
#                   include directory of any dependencies.
# GLOG_LIBRARY: glog library, not including the libraries of any
#               dependencies.

# Called if we failed to find glog or any of it's required dependencies,
# unsets all public (designed to be used externally) variables and reports
# error message at priority depending upon [REQUIRED/QUIET/<NONE>] argument.
MACRO(GLOG_REPORT_NOT_FOUND REASON_MSG)
  UNSET(GLOG_FOUND)
  UNSET(GLOG_INCLUDE_DIRS)
  UNSET(GLOG_LIBRARIES)
  # Make results of search visible in the CMake GUI if glog has not
  # been found so that user does not have to toggle to advanced view.
  MARK_AS_ADVANCED(CLEAR GLOG_INCLUDE_DIR
                         GLOG_LIBRARY)
  # Note <package>_FIND_[REQUIRED/QUIETLY] variables defined by FindPackage()
  # use the camelcase library name, not uppercase.
  IF (Glog_FIND_QUIETLY)
    MESSAGE(STATUS "Failed to find glog - " ${REASON_MSG} ${ARGN})
  ELSEIF (Glog_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Failed to find glog - " ${REASON_MSG} ${ARGN})
  ELSE()
    # Neither QUIETLY nor REQUIRED, use no priority which emits a message
    # but continues configuration and allows generation.
    MESSAGE("-- Failed to find glog - " ${REASON_MSG} ${ARGN})
  ENDIF ()
ENDMACRO(GLOG_REPORT_NOT_FOUND)

# Search user-installed locations first, so that we prefer user installs
# to system installs where both exist.
#
# TODO: Add standard Windows search locations for glog.
LIST(APPEND GLOG_CHECK_INCLUDE_DIRS
  /usr/local/include
  /usr/local/homebrew/include # Mac OS X
  /opt/local/var/macports/software # Mac OS X.
  /opt/local/include
  /usr/include)
LIST(APPEND GLOG_CHECK_LIBRARY_DIRS
  /usr/local/lib
  /usr/local/homebrew/lib # Mac OS X.
  /opt/local/lib
  /usr/lib)

# Search supplied hint directories first if supplied.
FIND_PATH(GLOG_INCLUDE_DIR
  NAMES glog/logging.h
  PATHS ${GLOG_INCLUDE_DIR_HINTS}
  ${GLOG_CHECK_INCLUDE_DIRS})
IF (NOT GLOG_INCLUDE_DIR OR
    NOT EXISTS ${GLOG_INCLUDE_DIR})
  GLOG_REPORT_NOT_FOUND(
    "Could not find glog include directory, set GLOG_INCLUDE_DIR "
    "to directory containing glog/logging.h")
ENDIF (NOT GLOG_INCLUDE_DIR OR
       NOT EXISTS ${GLOG_INCLUDE_DIR})

FIND_LIBRARY(GLOG_LIBRARY NAMES glog
  PATHS ${GLOG_LIBRARY_DIR_HINTS}
  ${GLOG_CHECK_LIBRARY_DIRS})
IF (NOT GLOG_LIBRARY OR
    NOT EXISTS ${GLOG_LIBRARY})
  GLOG_REPORT_NOT_FOUND(
    "Could not find glog library, set GLOG_LIBRARY "
    "to full path to libglog.")
ENDIF (NOT GLOG_LIBRARY OR
       NOT EXISTS ${GLOG_LIBRARY})

# Mark internally as found, then verify. GLOG_REPORT_NOT_FOUND() unsets
# if called.
SET(GLOG_FOUND TRUE)

# Glog does not seem to provide any record of the version in its
# source tree, thus cannot extract version.

# Catch case when caller has set GLOG_INCLUDE_DIR in the cache / GUI and
# thus FIND_[PATH/LIBRARY] are not called, but specified locations are
# invalid, otherwise we would report the library as found.
IF (GLOG_INCLUDE_DIR AND
    NOT EXISTS ${GLOG_INCLUDE_DIR}/glog/logging.h)
  GLOG_REPORT_NOT_FOUND(
    "Caller defined GLOG_INCLUDE_DIR:"
    " ${GLOG_INCLUDE_DIR} does not contain glog/logging.h header.")
ENDIF (GLOG_INCLUDE_DIR AND
       NOT EXISTS ${GLOG_INCLUDE_DIR}/glog/logging.h)
# TODO: This regex for glog library is pretty primitive, we use lowercase
#       for comparison to handle Windows using CamelCase library names, could
#       this check be better?
STRING(TOLOWER "${GLOG_LIBRARY}" LOWERCASE_GLOG_LIBRARY)
IF (GLOG_LIBRARY AND
    NOT "${LOWERCASE_GLOG_LIBRARY}" MATCHES ".*glog[^/]*")
  GLOG_REPORT_NOT_FOUND(
    "Caller defined GLOG_LIBRARY: "
    "${GLOG_LIBRARY} does not match glog.")
ENDIF (GLOG_LIBRARY AND
       NOT "${LOWERCASE_GLOG_LIBRARY}" MATCHES ".*glog[^/]*")

# Set standard CMake FindPackage variables if found.
IF (GLOG_FOUND)
  SET(GLOG_INCLUDE_DIRS ${GLOG_INCLUDE_DIR})
  SET(GLOG_LIBRARIES ${GLOG_LIBRARY})
ENDIF (GLOG_FOUND)

# Handle REQUIRED / QUIET optional arguments.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Glog DEFAULT_MSG
  GLOG_INCLUDE_DIRS GLOG_LIBRARIES)

# Only mark internal variables as advanced if we found glog, otherwise
# leave them visible in the standard GUI for the user to set manually.
IF (GLOG_FOUND)
  MARK_AS_ADVANCED(FORCE GLOG_INCLUDE_DIR
                         GLOG_LIBRARY)
ENDIF (GLOG_FOUND)
