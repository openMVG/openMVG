###########################################################
#                  Find Flann Library
#----------------------------------------------------------

# This sets the following variables:
# FLANN_FOUND - True if FLANN was found.
# FLANN_INCLUDE_DIRS - Directories containing the FLANN include files.
# FLANN_LIBRARIES - Libraries needed to use FLANN.
# FLANN_DEFINITIONS - Compiler flags for FLANN.

##====================================================
## Use PKG-CONFIG for local FLANN instance
##----------------------------------------------------
find_package(PkgConfig)
pkg_check_modules(PC_FLANN flann)
IF(PC_FLANN_FOUND)
  set(FLANN_FOUND ${PC_FLANN_FOUND})
  set(FLANN_DEFINITIONS ${PC_FLANN_CFLAGS_OTHER})

  find_path(FLANN_INCLUDE_DIR flann/flann.hpp
      HINTS ${PC_FLANN_INCLUDEDIR} ${PC_FLANN_INCLUDE_DIRS})

  find_library(FLANN_LIBRARY flann
      HINTS ${PC_FLANN_LIBDIR} ${PC_FLANN_LIBRARY_DIRS})

  set(FLANN_INCLUDE_DIRS ${FLANN_INCLUDE_DIR})
  set(FLANN_LIBRARIES  ${PC_FLANN_LIBRARIES})
  set(FLANN_VERSION ${PC_FLANN_VERSION})

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Flann DEFAULT_MSG
      FLANN_LIBRARY FLANN_INCLUDE_DIR)

  mark_as_advanced(FLANN_LIBRARY FLANN_INCLUDE_DIR)
ELSE()
  FIND_PATH(FLANN_DIR flann.hpp
      HINTS "${FLANN_ROOT}" "$ENV{FLANN_ROOT}" "${FLANN_INCLUDE_DIR_HINTS}"
      PATHS "$ENV{PROGRAMFILES}/flann" "$ENV{PROGRAMW6432}/flann"
      PATH_SUFFIXES flann
      DOC "Root directory of FLANN includes")

  ##====================================================
  ## Include FLANN library
  ##----------------------------------------------------
  IF(EXISTS "${FLANN_DIR}" AND NOT "${FLANN_DIR}" STREQUAL "")
          SET(FLANN_FOUND TRUE)
          SET(FLANN_INCLUDE_DIRS ${FLANN_DIR})
          SET(FLANN_DIR "${FLANN_DIR}" CACHE PATH "" FORCE)
          MARK_AS_ADVANCED(FLANN_DIR)

          # Extract Flann version from config.h
          SET(FLANN_VERSION_FILE ${FLANN_INCLUDE_DIRS}/config.h)
          IF (NOT EXISTS ${FLANN_VERSION_FILE})
                  FLANN_REPORT_NOT_FOUND(
                    "Could not find file: ${FLANN_VERSION_FILE} "
                    "containing version information in Flann install located at: "
                    "${FLANN_INCLUDE_DIRS}.")
          ELSE (NOT EXISTS ${FLANN_VERSION_FILE})
              FILE(READ ${FLANN_VERSION_FILE} FLANN_VERSION_FILE_CONTENTS)
              STRING(REGEX MATCH "#define FLANN_VERSION_ \"([0-9.]+)\""
                FLANN_VERSION "${FLANN_VERSION_FILE_CONTENTS}")
              STRING(REGEX REPLACE "#define FLANN_VERSION_ \"([0-9.]+)\"" "\\1"
                FLANN_VERSION "${FLANN_VERSION}")
          ENDIF (NOT EXISTS ${FLANN_VERSION_FILE})
          SET(FLANN_INCLUDE_DIR ${FLANN_DIR})

          FIND_LIBRARY(FLANN_LIBRARY NAMES flann)

          # locate Flann libraries
          IF(DEFINED FLANN_LIBRARY)
            SET(FLANN_LIBRARIES ${FLANN_LIBRARY})
          ENDIF()
          SET(OpenMVG_USE_INTERNAL_FLANN TRUE)
          MESSAGE(STATUS "Flann ${FLANN_VERSION} found (include: ${FLANN_INCLUDE_DIRS})")
  ELSE()
    MESSAGE(FATAL_ERROR "You are attempting to build without Flann. "
            "Please use cmake variable -DFLANN_INCLUDE_DIR_HINTS:STRING=\"PATH\" "
            "or FLANN_INCLUDE_DIR_HINTS env. variable to a valid Flann path. "
            "Or install last Flann version.")
    package_report_not_found(FLANN "Flann cannot be found")
  ENDIF()
ENDIF()

##====================================================
