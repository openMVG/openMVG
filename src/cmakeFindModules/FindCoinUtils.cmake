###########################################################
#                  Find CoinUtils Library
#----------------------------------------------------------

FIND_PATH(COINUTILS_DIR CoinUtilsConfig.h
    HINTS "${COINUTILS_ROOT}" "$ENV{COINUTILS_ROOT}" "${COINUTILS_INCLUDE_DIR_HINTS}"
    PATHS "$ENV{PROGRAMFILES}/CoinUtils" "$ENV{PROGRAMW6432}/CoinUtils" "/usr" "/usr/local"
    PATH_SUFFIXES CoinUtils
    DOC "Root directory of COINUTILS includes")

##====================================================
## Include COINUTILS library
##----------------------------------------------------
IF(EXISTS "${COINUTILS_DIR}" AND NOT "${COINUTILS_DIR}" STREQUAL "")
        SET(COINUTILS_FOUND TRUE)
        SET(COINUTILS_INCLUDE_DIRS ${COINUTILS_DIR})
        SET(COINUTILS_DIR "${COINUTILS_DIR}" CACHE PATH "" FORCE)
        MARK_AS_ADVANCED(COINUTILS_DIR)

        # Extract CoinUtils version from CoinUtilsConfig.h
        SET(COINUTILS_VERSION_FILE ${COINUTILS_INCLUDE_DIRS}/CoinUtilsConfig.h)
        # Extract CoinUtils version from alternative config_coinutils_default.h
        IF (EXISTS ${COINUTILS_INCLUDE_DIRS}/config_coinutils_default.h)
          SET(COINUTILS_VERSION_FILE ${COINUTILS_INCLUDE_DIRS}/config_coinutils_default.h)
        ENDIF()
        IF (NOT EXISTS ${COINUTILS_VERSION_FILE})
                COINUTILS_REPORT_NOT_FOUND(
                  "Could not find file: ${COINUTILS_VERSION_FILE} "
                  "containing version information in CoinUtils install located at: "
                  "${COINUTILS_INCLUDE_DIRS}.")
        ELSE (NOT EXISTS ${COINUTILS_VERSION_FILE})
            FILE(READ ${COINUTILS_VERSION_FILE} COINUTILS_VERSION_FILE_CONTENTS)

                STRING(REGEX MATCH "#define COINUTILS_VERSION_MAJOR [0-9]+"
                  COINUTILS_VERSION_MAJOR "${COINUTILS_VERSION_FILE_CONTENTS}")
                STRING(REGEX REPLACE "#define COINUTILS_VERSION_MAJOR ([0-9]+)" "\\1"
                  COINUTILS_VERSION_MAJOR "${COINUTILS_VERSION_MAJOR}")

                STRING(REGEX MATCH "#define COINUTILS_VERSION_MINOR [0-9]+"
                  COINUTILS_VERSION_MINOR "${COINUTILS_VERSION_FILE_CONTENTS}")
                STRING(REGEX REPLACE "#define COINUTILS_VERSION_MINOR ([0-9]+)" "\\1"
                  COINUTILS_VERSION_MINOR "${COINUTILS_VERSION_MINOR}")

                STRING(REGEX MATCH "#define COINUTILS_VERSION_RELEASE [0-9]+"
                  COINUTILS_VERSION_RELEASE "${COINUTILS_VERSION_FILE_CONTENTS}")
                STRING(REGEX REPLACE "#define COINUTILS_VERSION_RELEASE ([0-9]+)" "\\1"
                  COINUTILS_VERSION_RELEASE "${COINUTILS_VERSION_RELEASE}")

                SET(COINUTILS_VERSION "${COINUTILS_VERSION_MAJOR}.${COINUTILS_VERSION_MINOR}.${COINUTILS_VERSION_RELEASE}")
        ENDIF (NOT EXISTS ${COINUTILS_VERSION_FILE})
        SET(COINUTILS_INCLUDE_DIR ${COINUTILS_DIR})

        FIND_LIBRARY(COINUTILS_LIBRARY NAMES CoinUtils)

        # locate CoinUtils libraries
        IF(DEFINED COINUTILS_LIBRARY)
          SET(COINUTILS_LIBRARIES ${COINUTILS_LIBRARY})
        ENDIF()

        MESSAGE(STATUS "CoinUtils ${COINUTILS_VERSION} found (include: ${COINUTILS_INCLUDE_DIRS})")
ELSE()
  MESSAGE(FATAL_ERROR "You are attempting to build without CoinUtils. "
          "Please use cmake variable -DCOINUTILS_INCLUDE_DIR_HINTS:STRING=\"PATH\" "
          "or COINUTILS_INCLUDE_DIR_HINTS env. variable to a valid CoinUtils path. "
          "Or install last CoinUtils version.")
  package_report_not_found(COINUTILS "CoinUtils cannot be found")
ENDIF()
##====================================================
