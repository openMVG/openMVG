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

# FindGflags.cmake - Find Google gflags logging library.
#
# This module defines the following variables:
#
# GFLAGS_FOUND: TRUE iff gflags is found.
# GFLAGS_INCLUDE_DIRS: Include directories for gflags.
# GFLAGS_LIBRARIES: Libraries required to link gflags.
# GFLAGS_NAMESPACE: The namespace in which gflags is defined.  In versions of
#                   gflags < 2.1, this was google, for versions >= 2.1 it is
#                   by default gflags, although can be configured when building
#                   gflags to be something else (i.e. google for legacy
#                   compatibility).
#
# The following variables control the behaviour of this module:
#
# GFLAGS_INCLUDE_DIR_HINTS: List of additional directories in which to
#                           search for gflags includes, e.g: /timbuktu/include.
# GFLAGS_LIBRARY_DIR_HINTS: List of additional directories in which to
#                           search for gflags libraries, e.g: /timbuktu/lib.
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
# GFLAGS_INCLUDE_DIR: Include directory for gflags, not including the
#                     include directory of any dependencies.
# GFLAGS_LIBRARY: gflags library, not including the libraries of any
#                 dependencies.

# Called if we failed to find gflags or any of it's required dependencies,
# unsets all public (designed to be used externally) variables and reports
# error message at priority depending upon [REQUIRED/QUIET/<NONE>] argument.
MACRO(GFLAGS_REPORT_NOT_FOUND REASON_MSG)
  UNSET(GFLAGS_FOUND)
  UNSET(GFLAGS_INCLUDE_DIRS)
  UNSET(GFLAGS_LIBRARIES)
  UNSET(GFLAGS_NAMESPACE)
  # Make results of search visible in the CMake GUI if gflags has not
  # been found so that user does not have to toggle to advanced view.
  MARK_AS_ADVANCED(CLEAR GFLAGS_INCLUDE_DIR
                         GFLAGS_LIBRARY)
  # Note <package>_FIND_[REQUIRED/QUIETLY] variables defined by FindPackage()
  # use the camelcase library name, not uppercase.
  IF (Gflags_FIND_QUIETLY)
    MESSAGE(STATUS "Failed to find gflags - " ${REASON_MSG} ${ARGN})
  ELSEIF (Gflags_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Failed to find gflags - " ${REASON_MSG} ${ARGN})
  ELSE()
    # Neither QUIETLY nor REQUIRED, use no priority which emits a message
    # but continues configuration and allows generation.
    MESSAGE("-- Failed to find gflags - " ${REASON_MSG} ${ARGN})
  ENDIF ()
ENDMACRO(GFLAGS_REPORT_NOT_FOUND)

# Search user-installed locations first, so that we prefer user installs
# to system installs where both exist.
#
# TODO: Add standard Windows search locations for gflags.
LIST(APPEND GFLAGS_CHECK_INCLUDE_DIRS
  /usr/local/include
  /usr/local/homebrew/include # Mac OS X
  /opt/local/var/macports/software # Mac OS X.
  /opt/local/include
  /usr/include)
LIST(APPEND GFLAGS_CHECK_LIBRARY_DIRS
  /usr/local/lib
  /usr/local/homebrew/lib # Mac OS X.
  /opt/local/lib
  /usr/lib)

# Search supplied hint directories first if supplied.
FIND_PATH(GFLAGS_INCLUDE_DIR
  NAMES gflags/gflags.h
  PATHS ${GFLAGS_INCLUDE_DIR_HINTS}
  ${GFLAGS_CHECK_INCLUDE_DIRS})
IF (NOT GFLAGS_INCLUDE_DIR OR
    NOT EXISTS ${GFLAGS_INCLUDE_DIR})
  GFLAGS_REPORT_NOT_FOUND(
    "Could not find gflags include directory, set GFLAGS_INCLUDE_DIR "
    "to directory containing gflags/gflags.h")
ENDIF (NOT GFLAGS_INCLUDE_DIR OR
       NOT EXISTS ${GFLAGS_INCLUDE_DIR})

FIND_LIBRARY(GFLAGS_LIBRARY NAMES gflags
  PATHS ${GFLAGS_LIBRARY_DIR_HINTS}
  ${GFLAGS_CHECK_LIBRARY_DIRS})
IF (NOT GFLAGS_LIBRARY OR
    NOT EXISTS ${GFLAGS_LIBRARY})
  GFLAGS_REPORT_NOT_FOUND(
    "Could not find gflags library, set GFLAGS_LIBRARY "
    "to full path to libgflags.")
ENDIF (NOT GFLAGS_LIBRARY OR
       NOT EXISTS ${GFLAGS_LIBRARY})

# Mark internally as found, then verify. GFLAGS_REPORT_NOT_FOUND() unsets
# if called.
SET(GFLAGS_FOUND TRUE)

# Identify what namespace gflags was built with.
IF (GFLAGS_INCLUDE_DIR)
  # First try the (older) google namespace.
  INCLUDE(CheckCXXSourceCompiles)
  # Setup include path & link library for gflags for CHECK_CXX_SOURCE_COMPILES
  SET(CMAKE_REQUIRED_INCLUDES ${GFLAGS_INCLUDE_DIR})
  SET(CMAKE_REQUIRED_LIBRARIES ${GFLAGS_LIBRARY})
  CHECK_CXX_SOURCE_COMPILES(
    "#include <gflags/gflags.h>
     int main(int argc, char * argv[]) {
       google::ParseCommandLineFlags(&argc, &argv, true);
       return 0;
     }"
     GFLAGS_IN_GOOGLE_NAMESPACE)
  IF (GFLAGS_IN_GOOGLE_NAMESPACE)
    SET(GFLAGS_NAMESPACE google)
  ELSE (GFLAGS_IN_GOOGLE_NAMESPACE)
    # Try (newer) gflags namespace instead.
    #
    # Setup include path & link library for gflags for CHECK_CXX_SOURCE_COMPILES
    SET(CMAKE_REQUIRED_INCLUDES ${GFLAGS_INCLUDE_DIR})
    SET(CMAKE_REQUIRED_LIBRARIES ${GFLAGS_LIBRARY})
    CHECK_CXX_SOURCE_COMPILES(
      "#include <gflags/gflags.h>
       int main(int argc, char * argv[]) {
         gflags::ParseCommandLineFlags(&argc, &argv, true);
         return 0;
       }"
       GFLAGS_IN_GFLAGS_NAMESPACE)
    IF (GFLAGS_IN_GFLAGS_NAMESPACE)
      SET(GFLAGS_NAMESPACE gflags)
    ENDIF (GFLAGS_IN_GFLAGS_NAMESPACE)
  ENDIF (GFLAGS_IN_GOOGLE_NAMESPACE)

  IF (NOT GFLAGS_NAMESPACE)
    GFLAGS_REPORT_NOT_FOUND(
      "Failed to determine gflags namespace, it is not google or gflags.")
  ENDIF (NOT GFLAGS_NAMESPACE)
ENDIF (GFLAGS_INCLUDE_DIR)

# gflags does not seem to provide any record of the version in its
# source tree, thus cannot extract version.

# Catch case when caller has set GFLAGS_INCLUDE_DIR in the cache / GUI and
# thus FIND_[PATH/LIBRARY] are not called, but specified locations are
# invalid, otherwise we would report the library as found.
IF (GFLAGS_INCLUDE_DIR AND
    NOT EXISTS ${GFLAGS_INCLUDE_DIR}/gflags/gflags.h)
  GFLAGS_REPORT_NOT_FOUND(
    "Caller defined GFLAGS_INCLUDE_DIR:"
    " ${GFLAGS_INCLUDE_DIR} does not contain gflags/gflags.h header.")
ENDIF (GFLAGS_INCLUDE_DIR AND
       NOT EXISTS ${GFLAGS_INCLUDE_DIR}/gflags/gflags.h)
# TODO: This regex for gflags library is pretty primitive, we use lowercase
#       for comparison to handle Windows using CamelCase library names, could
#       this check be better?
STRING(TOLOWER "${GFLAGS_LIBRARY}" LOWERCASE_GFLAGS_LIBRARY)
IF (GFLAGS_LIBRARY AND
    NOT "${LOWERCASE_GFLAGS_LIBRARY}" MATCHES ".*gflags[^/]*")
  GFLAGS_REPORT_NOT_FOUND(
    "Caller defined GFLAGS_LIBRARY: "
    "${GFLAGS_LIBRARY} does not match gflags.")
ENDIF (GFLAGS_LIBRARY AND
       NOT "${LOWERCASE_GFLAGS_LIBRARY}" MATCHES ".*gflags[^/]*")

# Set standard CMake FindPackage variables if found.
IF (GFLAGS_FOUND)
  SET(GFLAGS_INCLUDE_DIRS ${GFLAGS_INCLUDE_DIR})
  SET(GFLAGS_LIBRARIES ${GFLAGS_LIBRARY})
ENDIF (GFLAGS_FOUND)

# Handle REQUIRED / QUIET optional arguments.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Gflags DEFAULT_MSG
  GFLAGS_INCLUDE_DIRS GFLAGS_LIBRARIES GFLAGS_NAMESPACE)

# Only mark internal variables as advanced if we found gflags, otherwise
# leave them visible in the standard GUI for the user to set manually.
IF (GFLAGS_FOUND)
  MARK_AS_ADVANCED(FORCE GFLAGS_INCLUDE_DIR
                         GFLAGS_LIBRARY
                         GFLAGS_NAMESPACE)
ENDIF (GFLAGS_FOUND)
