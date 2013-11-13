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

# FindCXSparse.cmake - Find CXSparse libraries & dependencies.
#
# This module defines the following variables which should be referenced
# by the caller to use the library.
#
# CXSPARSE_FOUND: TRUE iff CXSparse and all dependencies have been found.
# CXSPARSE_INCLUDE_DIRS: Include directories for CXSparse.
# CXSPARSE_LIBRARIES: Libraries for CXSparse and all dependencies.
#
# CXSPARSE_VERSION: Extracted from cs.h.
# CXSPARSE_MAIN_VERSION: Equal to 3 if CXSPARSE_VERSION = 3.1.2
# CXSPARSE_SUB_VERSION: Equal to 1 if CXSPARSE_VERSION = 3.1.2
# CXSPARSE_SUBSUB_VERSION: Equal to 2 if CXSPARSE_VERSION = 3.1.2
#
# The following variables control the behaviour of this module:
#
# CXSPARSE_INCLUDE_DIR_HINTS: List of additional directories in which to
#                             search for CXSparse includes,
#                             e.g: /timbuktu/include.
# CXSPARSE_LIBRARY_DIR_HINTS: List of additional directories in which to
#                             search for CXSparse libraries, e.g: /timbuktu/lib.
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
# CXSPARSE_INCLUDE_DIR: Include directory for CXSparse, not including the
#                       include directory of any dependencies.
# CXSPARSE_LIBRARY: CXSparse library, not including the libraries of any
#                   dependencies.

# Called if we failed to find CXSparse or any of it's required dependencies,
# unsets all public (designed to be used externally) variables and reports
# error message at priority depending upon [REQUIRED/QUIET/<NONE>] argument.
MACRO(CXSPARSE_REPORT_NOT_FOUND REASON_MSG)
  UNSET(CXSPARSE_FOUND)
  UNSET(CXSPARSE_INCLUDE_DIRS)
  UNSET(CXSPARSE_LIBRARIES)
  # Make results of search visible in the CMake GUI if CXSparse has not
  # been found so that user does not have to toggle to advanced view.
  MARK_AS_ADVANCED(CLEAR CXSPARSE_INCLUDE_DIR
                         CXSPARSE_LIBRARY)
  # Note <package>_FIND_[REQUIRED/QUIETLY] variables defined by FindPackage()
  # use the camelcase library name, not uppercase.
  IF (CXSparse_FIND_QUIETLY)
    MESSAGE(STATUS "Failed to find CXSparse - " ${REASON_MSG} ${ARGN})
  ELSEIF (CXSparse_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Failed to find CXSparse - " ${REASON_MSG} ${ARGN})
  ELSE()
    # Neither QUIETLY nor REQUIRED, use no priority which emits a message
    # but continues configuration and allows generation.
    MESSAGE("-- Failed to find CXSparse - " ${REASON_MSG} ${ARGN})
  ENDIF ()
ENDMACRO(CXSPARSE_REPORT_NOT_FOUND)

# Search user-installed locations first, so that we prefer user installs
# to system installs where both exist.
#
# TODO: Add standard Windows search locations for CXSparse.
LIST(APPEND CXSPARSE_CHECK_INCLUDE_DIRS
  /usr/local/include
  /usr/local/homebrew/include # Mac OS X
  /opt/local/var/macports/software # Mac OS X.
  /opt/local/include
  /usr/include)
LIST(APPEND CXSPARSE_CHECK_LIBRARY_DIRS
  /usr/local/lib
  /usr/local/homebrew/lib # Mac OS X.
  /opt/local/lib
  /usr/lib)

# Search supplied hint directories first if supplied.
FIND_PATH(CXSPARSE_INCLUDE_DIR
  NAMES cs.h
  PATHS ${CXSPARSE_INCLUDE_DIR_HINTS}
  ${CXSPARSE_CHECK_INCLUDE_DIRS})
IF (NOT CXSPARSE_INCLUDE_DIR OR
    NOT EXISTS ${CXSPARSE_INCLUDE_DIR})
  CXSPARSE_REPORT_NOT_FOUND(
    "Could not find CXSparse include directory, set CXSPARSE_INCLUDE_DIR "
    "to directory containing cs.h")
ENDIF (NOT CXSPARSE_INCLUDE_DIR OR
       NOT EXISTS ${CXSPARSE_INCLUDE_DIR})

FIND_LIBRARY(CXSPARSE_LIBRARY NAMES cxsparse
  PATHS ${CXSPARSE_LIBRARY_DIR_HINTS}
  ${CXSPARSE_CHECK_LIBRARY_DIRS})
IF (NOT CXSPARSE_LIBRARY OR
    NOT EXISTS ${CXSPARSE_LIBRARY})
  CXSPARSE_REPORT_NOT_FOUND(
    "Could not find CXSparse library, set CXSPARSE_LIBRARY "
    "to full path to libcxsparse.")
ENDIF (NOT CXSPARSE_LIBRARY OR
       NOT EXISTS ${CXSPARSE_LIBRARY})

# Mark internally as found, then verify. CXSPARSE_REPORT_NOT_FOUND() unsets
# if called.
SET(CXSPARSE_FOUND TRUE)

# Extract CXSparse version from cs.h
IF (CXSPARSE_INCLUDE_DIR)
  SET(CXSPARSE_VERSION_FILE ${CXSPARSE_INCLUDE_DIR}/cs.h)
  IF (NOT EXISTS ${CXSPARSE_VERSION_FILE})
    CXSPARSE_REPORT_NOT_FOUND(
      "Could not find file: ${CXSPARSE_VERSION_FILE} "
      "containing version information in CXSparse install located at: "
      "${CXSPARSE_INCLUDE_DIR}.")
  ELSE (NOT EXISTS ${CXSPARSE_VERSION_FILE})
    FILE(READ ${CXSPARSE_INCLUDE_DIR}/cs.h CXSPARSE_VERSION_FILE_CONTENTS)

    STRING(REGEX MATCH "#define CS_VER [0-9]+"
      CXSPARSE_MAIN_VERSION "${CXSPARSE_VERSION_FILE_CONTENTS}")
    STRING(REGEX REPLACE "#define CS_VER ([0-9]+)" "\\1"
      CXSPARSE_MAIN_VERSION "${CXSPARSE_MAIN_VERSION}")

    STRING(REGEX MATCH "#define CS_SUBVER [0-9]+"
      CXSPARSE_SUB_VERSION "${CXSPARSE_VERSION_FILE_CONTENTS}")
    STRING(REGEX REPLACE "#define CS_SUBVER ([0-9]+)" "\\1"
      CXSPARSE_SUB_VERSION "${CXSPARSE_SUB_VERSION}")

    STRING(REGEX MATCH "#define CS_SUBSUB [0-9]+"
      CXSPARSE_SUBSUB_VERSION "${CXSPARSE_VERSION_FILE_CONTENTS}")
    STRING(REGEX REPLACE "#define CS_SUBSUB ([0-9]+)" "\\1"
      CXSPARSE_SUBSUB_VERSION "${CXSPARSE_SUBSUB_VERSION}")

    # This is on a single line s/t CMake does not interpret it as a list of
    # elements and insert ';' separators which would result in 3.;1.;2 nonsense.
    SET(CXSPARSE_VERSION "${CXSPARSE_MAIN_VERSION}.${CXSPARSE_SUB_VERSION}.${CXSPARSE_SUBSUB_VERSION}")
  ENDIF (NOT EXISTS ${CXSPARSE_VERSION_FILE})
ENDIF (CXSPARSE_INCLUDE_DIR)

# Catch the case when the caller has set CXSPARSE_LIBRARY in the cache / GUI and
# thus FIND_LIBRARY was not called, but specified library is invalid, otherwise
# we would report CXSparse as found.
# TODO: This regex for CXSparse library is pretty primitive, we use lowercase
#       for comparison to handle Windows using CamelCase library names, could
#       this check be better?
STRING(TOLOWER "${CXSPARSE_LIBRARY}" LOWERCASE_CXSPARSE_LIBRARY)
IF (CXSPARSE_LIBRARY AND
    EXISTS ${CXSPARSE_LIBRARY} AND
    NOT "${LOWERCASE_CXSPARSE_LIBRARY}" MATCHES ".*cxsparse[^/]*")
  CXSPARSE_REPORT_NOT_FOUND(
    "Caller defined CXSPARSE_LIBRARY: "
    "${CXSPARSE_LIBRARY} does not match CXSparse.")
ENDIF (CXSPARSE_LIBRARY AND
       EXISTS ${CXSPARSE_LIBRARY} AND
       NOT "${LOWERCASE_CXSPARSE_LIBRARY}" MATCHES ".*cxsparse[^/]*")

# Set standard CMake FindPackage variables if found.
IF (CXSPARSE_FOUND)
  SET(CXSPARSE_INCLUDE_DIRS ${CXSPARSE_INCLUDE_DIR})
  SET(CXSPARSE_LIBRARIES ${CXSPARSE_LIBRARY})
ENDIF (CXSPARSE_FOUND)

# Handle REQUIRED / QUIET optional arguments and version.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CXSparse
  REQUIRED_VARS CXSPARSE_INCLUDE_DIRS CXSPARSE_LIBRARIES
  VERSION_VAR CXSPARSE_VERSION)

# Only mark internal variables as advanced if we found CXSparse, otherwise
# leave them visible in the standard GUI for the user to set manually.
IF (CXSPARSE_FOUND)
  MARK_AS_ADVANCED(FORCE CXSPARSE_INCLUDE_DIR
                         CXSPARSE_LIBRARY)
ENDIF (CXSPARSE_FOUND)
