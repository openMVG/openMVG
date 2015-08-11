###########################################################
#                  Find Clp Library
#----------------------------------------------------------

FIND_PATH(CLP_DIR ClpConfig.h
    HINTS "${CLP_ROOT}" "$ENV{CLP_ROOT}" "${CLP_INCLUDE_DIR_HINTS}"
    PATHS "$ENV{PROGRAMFILES}/Clp" "$ENV{PROGRAMW6432}/Clp" "/usr" "/usr/local"
    PATH_SUFFIXES Clp
    DOC "Root directory of CLP includes")

##====================================================
## Include CLP library
##----------------------------------------------------
IF(EXISTS "${CLP_DIR}" AND NOT "${CLP_DIR}" STREQUAL "")
        SET(CLP_FOUND TRUE)
        SET(CLP_INCLUDE_DIRS ${CLP_DIR})
        SET(CLP_DIR "${CLP_DIR}" CACHE PATH "" FORCE)
        MARK_AS_ADVANCED(CLP_DIR)

        # Extract Clp version from ClpConfig.h
        SET(CLP_VERSION_FILE ${CLP_INCLUDE_DIRS}/ClpConfig.h)
        # Extract Clp version from alternative config_clp_default.h
        IF (EXISTS ${CLP_INCLUDE_DIRS}/config_clp_default.h)
          SET(CLP_VERSION_FILE ${CLP_INCLUDE_DIRS}/config_clp_default.h)
        ENDIF()
        IF (NOT EXISTS ${CLP_VERSION_FILE})
                CLP_REPORT_NOT_FOUND(
                  "Could not find file: ${CLP_VERSION_FILE} "
                  "containing version information in Clp install located at: "
                  "${CLP_INCLUDE_DIRS}.")
        ELSE (NOT EXISTS ${CLP_VERSION_FILE})
            FILE(READ ${CLP_VERSION_FILE} CLP_VERSION_FILE_CONTENTS)

                STRING(REGEX MATCH "#define CLP_VERSION_MAJOR [0-9]+"
                  CLP_VERSION_MAJOR "${CLP_VERSION_FILE_CONTENTS}")
                STRING(REGEX REPLACE "#define CLP_VERSION_MAJOR ([0-9]+)" "\\1"
                  CLP_VERSION_MAJOR "${CLP_VERSION_MAJOR}")

                STRING(REGEX MATCH "#define CLP_VERSION_MINOR [0-9]+"
                  CLP_VERSION_MINOR "${CLP_VERSION_FILE_CONTENTS}")
                STRING(REGEX REPLACE "#define CLP_VERSION_MINOR ([0-9]+)" "\\1"
                  CLP_VERSION_MINOR "${CLP_VERSION_MINOR}")

                STRING(REGEX MATCH "#define CLP_VERSION_RELEASE [0-9]+"
                  CLP_VERSION_RELEASE "${CLP_VERSION_FILE_CONTENTS}")
                STRING(REGEX REPLACE "#define CLP_VERSION_RELEASE ([0-9]+)" "\\1"
                  CLP_VERSION_RELEASE "${CLP_VERSION_RELEASE}")

                SET(CLP_VERSION "${CLP_VERSION_MAJOR}.${CLP_VERSION_MINOR}.${CLP_VERSION_RELEASE}")
        ENDIF (NOT EXISTS ${CLP_VERSION_FILE})
        SET(CLP_INCLUDE_DIR ${CLP_DIR})

        FIND_LIBRARY(CLP_LIBRARY NAMES Clp)
        FIND_LIBRARY(CLPSOLVER_LIBRARY NAMES ClpSolver)
        FIND_LIBRARY(OSICLP_LIBRARY NAMES OsiClp)

        # locate Clp libraries
        IF(DEFINED CLP_LIBRARY AND DEFINED CLPSOLVER_LIBRARY AND DEFINED OSICLP_LIBRARY)
          SET(CLP_LIBRARIES ${CLP_LIBRARY} ${CLPSOLVER_LIBRARY} ${OSICLP_LIBRARY})
        ENDIF()

        MESSAGE(STATUS "Clp ${CLP_VERSION} found (include: ${CLP_INCLUDE_DIRS})")
ELSE()
  MESSAGE(FATAL_ERROR "You are attempting to build without Clp. "
          "Please use cmake variable -DCLP_INCLUDE_DIR_HINTS:STRING=\"PATH\" "
          "or CLP_INCLUDE_DIR_HINTS env. variable to a valid Clp path. "
          "Or install last Clp version.")
  package_report_not_found(CLP "Clp cannot be found")
ENDIF()
##====================================================
