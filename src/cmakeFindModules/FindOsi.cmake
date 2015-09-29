###########################################################
#                  Find Osi Library
#----------------------------------------------------------

FIND_PATH(OSI_DIR OsiConfig.h
    HINTS "${OSI_ROOT}" "$ENV{OSI_ROOT}" "${OSI_INCLUDE_DIR_HINTS}"
    PATHS "$ENV{PROGRAMFILES}/Osi" "$ENV{PROGRAMW6432}/Osi" "/usr" "/usr/local"
    PATH_SUFFIXES Osi
    DOC "Root directory of OSI includes")

##====================================================
## Include OSI library
##----------------------------------------------------
IF(EXISTS "${OSI_DIR}" AND NOT "${OSI_DIR}" STREQUAL "")
        SET(OSI_FOUND TRUE)
        SET(OSI_INCLUDE_DIRS ${OSI_DIR})
        SET(OSI_DIR "${OSI_DIR}" CACHE PATH "" FORCE)
        MARK_AS_ADVANCED(OSI_DIR)

        # Extract Osi version from OsiConfig.h
        SET(OSI_VERSION_FILE ${OSI_INCLUDE_DIRS}/OsiConfig.h)
        # Extract Osi version from alternative config_osi_default.h
        IF (EXISTS ${OSI_INCLUDE_DIRS}/config_osi_default.h)
          SET(OSI_VERSION_FILE ${OSI_INCLUDE_DIRS}/config_osi_default.h)
        ENDIF()
        IF (NOT EXISTS ${OSI_VERSION_FILE})
                OSI_REPORT_NOT_FOUND(
                  "Could not find file: ${OSI_VERSION_FILE} "
                  "containing version information in Osi install located at: "
                  "${OSI_INCLUDE_DIRS}.")
        ELSE (NOT EXISTS ${OSI_VERSION_FILE})
            FILE(READ ${OSI_VERSION_FILE} OSI_VERSION_FILE_CONTENTS)

                STRING(REGEX MATCH "#define OSI_VERSION_MAJOR [0-9]+"
                  OSI_VERSION_MAJOR "${OSI_VERSION_FILE_CONTENTS}")
                STRING(REGEX REPLACE "#define OSI_VERSION_MAJOR ([0-9]+)" "\\1"
                  OSI_VERSION_MAJOR "${OSI_VERSION_MAJOR}")

                STRING(REGEX MATCH "#define OSI_VERSION_MINOR [0-9]+"
                  OSI_VERSION_MINOR "${OSI_VERSION_FILE_CONTENTS}")
                STRING(REGEX REPLACE "#define OSI_VERSION_MINOR ([0-9]+)" "\\1"
                  OSI_VERSION_MINOR "${OSI_VERSION_MINOR}")

                STRING(REGEX MATCH "#define OSI_VERSION_RELEASE [0-9]+"
                  OSI_VERSION_RELEASE "${OSI_VERSION_FILE_CONTENTS}")
                STRING(REGEX REPLACE "#define OSI_VERSION_RELEASE ([0-9]+)" "\\1"
                  OSI_VERSION_RELEASE "${OSI_VERSION_RELEASE}")

                SET(OSI_VERSION "${OSI_VERSION_MAJOR}.${OSI_VERSION_MINOR}.${OSI_VERSION_RELEASE}")
        ENDIF (NOT EXISTS ${OSI_VERSION_FILE})
        SET(OSI_INCLUDE_DIR ${OSI_DIR})

        FIND_LIBRARY(OSI_LIBRARY NAMES Osi)

        # locate Osi libraries
        IF(DEFINED OSI_LIBRARY)
          SET(OSI_LIBRARIES ${OSI_LIBRARY})
        ENDIF()

        MESSAGE(STATUS "Osi ${OSI_VERSION} found (include: ${OSI_INCLUDE_DIRS})")
ELSE()
  MESSAGE(FATAL_ERROR "You are attempting to build without Osi. "
          "Please use cmake variable -DOSI_INCLUDE_DIR_HINTS:STRING=\"PATH\" "
          "or OSI_INCLUDE_DIR_HINTS env. variable to a valid Osi path. "
          "Or install last Osi version.")
  package_report_not_found(OSI "Osi cannot be found")
ENDIF()
##====================================================
