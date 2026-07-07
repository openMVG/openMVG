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

# This must take place outside of CONFIGURE_CERES_CONFIG() in order that
# we can determine where *this* file is, and thus the relative path to
# config.h.in.  Inside of CONFIGURE_CERES_CONFIG(), CMAKE_CURRENT_LIST_DIR
# refers to the caller of CONFIGURE_CERES_CONFIG(), not this file.
set(CERES_CONFIG_IN_FILE "${CMAKE_CURRENT_LIST_DIR}/config.h.in")

# CreateCeresConfig.cmake - Create the config.h for Ceres.
#
# This function configures the Ceres config.h file based on the current
# compile options and copies it into the specified location.  It should be
# called before Ceres is built so that the correct config.h is used when
# Ceres is compiled.
#
# INPUTS:
#   CURRENT_CERES_COMPILE_OPTIONS: List of currently enabled Ceres compile
#                                  options. These are compared against the
#                                  full list of valid options, which are read
#                                  from config.h.in.  Any options present
#                                  which are not part of the valid set will
#                                  invoke an error.  Any valid option present
#                                  will be enabled in the resulting config.h,
#                                  all other options will be disabled.
#
#   CERES_CONFIG_OUTPUT_DIRECTORY: Path to output directory in which to save
#                                  the configured config.h.  Typically this
#                                  will be <src>/include/ceres/internal.

function(CREATE_CERES_CONFIG CURRENT_CERES_COMPILE_OPTIONS CERES_CONFIG_OUTPUT_DIRECTORY)
  # Create the specified output directory if it does not exist.
  if (NOT EXISTS "${CERES_CONFIG_OUTPUT_DIRECTORY}")
    message(STATUS "Creating configured Ceres config.h output directory: "
      "${CERES_CONFIG_OUTPUT_DIRECTORY}")
    file(MAKE_DIRECTORY "${CERES_CONFIG_OUTPUT_DIRECTORY}")
  endif()
  if (EXISTS "${CERES_CONFIG_OUTPUT_DIRECTORY}" AND
      NOT IS_DIRECTORY "${CERES_CONFIG_OUTPUT_DIRECTORY}")
    message(FATAL_ERROR "Ceres Bug: Specified CERES_CONFIG_OUTPUT_DIRECTORY: "
      "${CERES_CONFIG_OUTPUT_DIRECTORY} exists, but is not a directory.")
  endif()

  # Read all possible configurable compile options from config.h.in, this avoids
  # us having to hard-code in this file what the valid options are.
  file(READ ${CERES_CONFIG_IN_FILE} CERES_CONFIG_IN_CONTENTS)
  string(REGEX MATCHALL "@[^@ $]+@"
    ALL_CONFIGURABLE_CERES_OPTIONS "${CERES_CONFIG_IN_CONTENTS}")
  # Removing @ symbols at beginning and end of each option.
  string(REPLACE "@" ""
    ALL_CONFIGURABLE_CERES_OPTIONS "${ALL_CONFIGURABLE_CERES_OPTIONS}")

  # Ensure that there are no repetitions in the current compile options.
  list(REMOVE_DUPLICATES CURRENT_CERES_COMPILE_OPTIONS)

  foreach (CERES_OPTION ${ALL_CONFIGURABLE_CERES_OPTIONS})
    # Try and find the option in the list of current compile options, if it
    # is present, then the option is enabled, otherwise it is disabled.
    list(FIND CURRENT_CERES_COMPILE_OPTIONS ${CERES_OPTION} OPTION_ENABLED)

    # list(FIND ..) returns -1 if the element was not in the list, but CMake
    # interprets if (VAR) to be true if VAR is any non-zero number, even
    # negative ones, hence we have to explicitly check for >= 0.
    if (OPTION_ENABLED GREATER -1)
      message(STATUS "Enabling ${CERES_OPTION} in Ceres config.h")
      set(${CERES_OPTION} "#define ${CERES_OPTION}")

      # Remove the item from the list of current options so that we can identify
      # any options that were in CURRENT_CERES_COMPILE_OPTIONS, but not in
      # ALL_CONFIGURABLE_CERES_OPTIONS (which is an error).
      list(REMOVE_ITEM CURRENT_CERES_COMPILE_OPTIONS ${CERES_OPTION})
    else()
      set(${CERES_OPTION} "// #define ${CERES_OPTION}")
    endif()
  endforeach()

  # CURRENT_CERES_COMPILE_OPTIONS should now be an empty list, any elements
  # remaining were not present in ALL_CONFIGURABLE_CERES_OPTIONS read from
  # config.h.in.
  if (CURRENT_CERES_COMPILE_OPTIONS)
    message(FATAL_ERROR "Ceres Bug: CURRENT_CERES_COMPILE_OPTIONS contained "
      "the following options which were not present in config.h.in: "
      "${CURRENT_CERES_COMPILE_OPTIONS}")
  endif()

  configure_file(${CERES_CONFIG_IN_FILE}
    "${CERES_CONFIG_OUTPUT_DIRECTORY}/config.h" @ONLY)
endfunction()
