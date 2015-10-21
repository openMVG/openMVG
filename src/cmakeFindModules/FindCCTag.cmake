# Locate the openMVG libraries.
#
# Defines the following variables:
#
#   CCTAG_FOUND        - TRUE if the CCTag headers and libs are found
#   CCTAG_INCLUDE_DIRS - The path to CCTag headers
#
#   CCTAG_LIBRARIES    - Libraries to link against to use CCTag.
#   CCTAG_LIBRARY_DIR  - The base directory to search for CCTag.
#
# Accepts the following variables as input:
#
#   CCTAG_DIR - (as a CMake or environment variable)
#                The root directory of the CCTag install prefix

MESSAGE(STATUS "Looking for CCTag.")

FIND_PATH(CCTAG_INCLUDE_DIR cctag/ICCTag.hpp
  HINTS
  $ENV{CCTAG_DIR}/include
  ${CCTAG_DIR}/include
  PATH_SUFFIXES
  CCTag
)

find_package(OpenCV QUIET)
find_package(OPTPP QUIET)
find_package(Ceres QUIET)
find_package(Glog QUIET)

IF(CCTAG_INCLUDE_DIR)
  MESSAGE(STATUS "CCTag headers found in ${CCTAG_INCLUDE_DIRS}")
ELSE()
  MESSAGE(STATUS "NOT FOUND")
ENDIF (CCTAG_INCLUDE_DIR)

SET(CCTAG_LIBRARIES_NAMES  
  CCTag
  #third_party libraries
  boost_filesystem
  boost_system
  boost_serialization
  dl
  ${OpenCV_LIBS}
  ${OPTPP_LIBRARIES}
  ${Ceres_LIBRARIES}
  ${GLOG_LIBRARIES}
)

FIND_LIBRARY(CCTAG_LIBRARY NAMES ${CCTAG_LIBRARIES_NAMES}
  HINTS
  $ENV{CCTAG_DIR}/lib
  ${CCTAG_DIR}/lib
  PATH_SUFFIXES
  CCTag
  #third_party libraries
  ${OpenCV_LIB_DIR}
  ${OPTPP_LIBRARY_DIRS}
  ${Ceres_LIBRARY}
  ${GLOG_LIBRARY_DIR_HINTS}
)
GET_FILENAME_COMPONENT(CCTAG_LIBRARY_DIR "${CCTAG_LIBRARY}" PATH)

SET(CCTAG_LIBRARY "")
FOREACH(lib ${CCTAG_LIBRARIES_NAMES})
 LIST(APPEND CCTAG_LIBRARY ${lib})  
ENDFOREACH()

SET(CCTAG_LIBRARIES ${CCTAG_LIBRARY})
SET(CCTAG_INCLUDE_DIRS ${CCTAG_INCLUDE_DIR})

IF(CCTAG_LIBRARY)
  MESSAGE(STATUS "CCTag libraries found: ${CCTAG_LIBRARY}")
  MESSAGE(STATUS "CCTag libraries directories: ${CCTAG_LIBRARY_DIR}")
ENDIF (CCTAG_LIBRARY)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CCTAG_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CCTag  DEFAULT_MSG
                                  CCTAG_LIBRARY CCTAG_INCLUDE_DIR)

MARK_AS_ADVANCED(CCTAG_INCLUDE_DIR CCTAG_LIBRARY)

#Third parties:
# - include directories

IF(CCTAG_FOUND)
  SET(CCTAG_INCLUDE_DIRS
    ${CCTAG_INCLUDE_DIR}
    ${OpenCV_INCLUDE_DIRS}
  )
ENDIF(CCTAG_FOUND)
