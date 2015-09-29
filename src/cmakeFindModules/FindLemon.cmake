###########################################################
#                  Find Lemon Library
#----------------------------------------------------------

FIND_PATH(LEMON_DIR list_graph.h
    HINTS "${LEMON_ROOT}" "$ENV{LEMON_ROOT}" "${LEMON_INCLUDE_DIR_HINTS}"
    PATHS "$ENV{PROGRAMFILES}/lemon" "$ENV{PROGRAMW6432}/lemon"
    PATH_SUFFIXES lemon
    DOC "Root directory of LEMON includes")

##====================================================
## Include LEMON library
##----------------------------------------------------
IF(EXISTS "${LEMON_DIR}" AND NOT "${LEMON_DIR}" STREQUAL "")
  SET(LEMON_FOUND TRUE)
  # Remove /lemon from path (math.h cannot be exposed all time)
  GET_FILENAME_COMPONENT(LEMON_INCLUDE_DIRS "${LEMON_DIR}" PATH)
  SET(LEMON_DIR "${LEMON_DIR}" CACHE PATH "" FORCE)
  MARK_AS_ADVANCED(LEMON_DIR)
  # Extract Lemon version from config.h
  SET(LEMON_VERSION_FILE ${LEMON_INCLUDE_DIRS}/lemon/config.h)
  IF (NOT EXISTS ${LEMON_VERSION_FILE})
    IF (Lemon_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR
        "Could not find file: ${LEMON_VERSION_FILE} "
        "containing version information in Lemon install located at: "
        "${LEMON_INCLUDE_DIRS}.")
    ELSEIF(NOT Lemon_FIND_QUIETLY)
      MESSAGE(SEND_ERROR
      "Could not find file: ${LEMON_VERSION_FILE} "
      "containing version information in Lemon install located at: "
      "${LEMON_INCLUDE_DIRS}.")
    ENDIF()
  ELSE (NOT EXISTS ${LEMON_VERSION_FILE})
    FILE(READ ${LEMON_VERSION_FILE} LEMON_VERSION_FILE_CONTENTS)
    STRING(REGEX MATCH "#define LEMON_VERSION \"([0-9.]+)\""
    LEMON_VERSION "${LEMON_VERSION_FILE_CONTENTS}")
    STRING(REGEX REPLACE "#define LEMON_VERSION \"([0-9.]+)\"" "\\1"
    LEMON_VERSION "${LEMON_VERSION}")
  ENDIF (NOT EXISTS ${LEMON_VERSION_FILE})
  
  SET(LEMON_INCLUDE_DIR ${LEMON_DIR})

  FIND_LIBRARY(LEMON_LIBRARY NAMES Lemon lemon emon)

  # locate Lemon libraries
  IF(DEFINED LEMON_LIBRARY)
    SET(LEMON_LIBRARIES ${LEMON_LIBRARY})
  ENDIF()

  MESSAGE(STATUS "Lemon ${LEMON_VERSION} found (include: ${LEMON_DIR})")
ELSE()
  MESSAGE(FATAL_ERROR "You are attempting to build without Lemon. "
          "Please use cmake variable -DLEMON_INCLUDE_DIR_HINTS:STRING=\"PATH\" "
          "or LEMON_INCLUDE_DIR_HINTS env. variable to a valid Lemon path. "
          "Or install last Lemon version.")
  package_report_not_found(LEMON "Lemon cannot be found")
ENDIF()
##====================================================
