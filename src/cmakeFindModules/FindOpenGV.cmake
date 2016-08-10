# Locate the openGV libraries.
#
# Defines the following variables:
#
#   OPENGV_FOUND        - TRUE if the openGV headers and libs are found
#   OPENGV_INCLUDE_DIRS - The path to openGV headers
#
#   OPENGV_LIBRARY      - The opengv library
#   OPENGV_LIBRARY_DIR  - The directory where the libraries are located
#
# Accepts the following variables as input:
#
#   OPENGV_DIR - (as a CMake or environment variable)
#                The root directory of the openGV install prefix

MESSAGE(STATUS "Looking for OpenGV.")

FIND_PATH(OPENGV_INCLUDE_DIR opengv/types.hpp
  HINTS
  $ENV{OPENGV_DIR}/include
  ${OPENGV_DIR}/include
  PATH_SUFFIXES
  openGV
)

IF(OPENGV_INCLUDE_DIR)
  MESSAGE(STATUS "OpenGV headers found in ${OPENGV_INCLUDE_DIR}")
  IF(NOT EIGEN_FOUND)
    MESSAGE(STATUS "Looking for Eigen dependency...")
    FIND_PACKAGE(Eigen QUIET)
    IF(EIGEN_FOUND)
        INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
    ELSE()
      MESSAGE(WARNING "Couldn't find Eigen, this is needed for compiling with openGV")
      # this is to make the find_package_handle_standard_args  fail
      SET(OPENGV_INCLUDE_DIR "OPENGV_INCLUDE_DIR-NOTFOUND")
    ENDIF(EIGEN_FOUND)
  ELSE(NOT EIGEN_FOUND)
    MESSAGE(STATUS "Eigen already found")
  ENDIF(NOT EIGEN_FOUND)
ELSE(OPENGV_INCLUDE_DIR)
  MESSAGE(STATUS "OpenGV headers not found!")
ENDIF(OPENGV_INCLUDE_DIR)

FIND_LIBRARY(OPENGV_LIBRARY NAMES opengv
  HINTS
  $ENV{OPENGV_DIR}/lib
  ${OPENGV_DIR}/lib
  PATH_SUFFIXES
  openGV
)

GET_FILENAME_COMPONENT(OPENGV_LIBRARY_DIR "${OPENGV_LIBRARY}" PATH)

IF(OPENGV_LIBRARY)
  MESSAGE(STATUS "OpenGV libraries found: ${OPENGV_LIBRARY}")
  MESSAGE(STATUS "OpenGV libraries directories: ${OPENGV_LIBRARY_DIR}")
ENDIF (OPENGV_LIBRARY)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set OPENGV_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(OpenGV  DEFAULT_MSG
                                  OPENGV_LIBRARY OPENGV_INCLUDE_DIR)

MARK_AS_ADVANCED(OPENGV_INCLUDE_DIR OPENGV_LIBRARY)