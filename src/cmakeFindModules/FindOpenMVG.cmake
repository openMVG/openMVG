# Locate the openMVG libraries.
#
# Defines the following variables:
#
#   OPENMVG_FOUND        - TRUE if the openMVG headers and libs are found
#   OPENMVG_INCLUDE_DIRS - The path to openMVG headers
#
#   OPENMVG_LIBRARIES    - All openMVG libraries
#   OPENMVG_LIBRARY_DIR  - The directory where the libraries are located
#
# Accepts the following variables as input:
#
#   OPENMVG_DIR - (as a CMake or environment variable)
#                The root directory of the openMVG install prefix

MESSAGE(STATUS "Looking for OpenMVG.")

FIND_PATH(OPENMVG_INCLUDE_DIR openMVG/version.hpp
  HINTS
  $ENV{OPENMVG_DIR}/include
  ${OPENMVG_DIR}/include
  PATH_SUFFIXES
  openMVG
)

IF(OPENMVG_INCLUDE_DIR)
  MESSAGE(STATUS "OpenMVG headers found in ${OPENMVG_INCLUDE_DIR}")
ELSE()
  MESSAGE(STATUS "NOT FOUND")
ENDIF (OPENMVG_INCLUDE_DIR)

SET(OPENMVG_LIBRARIES_NAMES  
  openMVG_numeric
  openMVG_system
  openMVG_image
  openMVG_kvld
  openMVG_lInftyComputerVision
  openMVG_multiview
  #third_party libraries
  ceres
  lemon
  stlplus
  easyexif
  #optional third_party
  vlsift
  jpeg
  png
  tiff
  zlib)

FIND_LIBRARY(OPENMVG_LIBRARY NAMES ${OPENMVG_LIBRARIES_NAMES}
  HINTS
  $ENV{OPENMVG_DIR}/lib
  ${OPENMVG_DIR}/lib
  PATH_SUFFIXES
  openMVG
)
GET_FILENAME_COMPONENT(OPENMVG_LIBRARY_DIR "${OPENMVG_LIBRARY}" PATH)

SET(OPENMVG_LIBRARY "")
FOREACH(lib ${OPENMVG_LIBRARIES_NAMES})
 LIST(APPEND OPENMVG_LIBRARY ${lib})  
ENDFOREACH()

SET(OPENMVG_LIBRARIES ${OPENMVG_LIBRARY})
SET(OPENMVG_INCLUDE_DIRS ${OPENMVG_INCLUDE_DIR})

IF(OPENMVG_LIBRARY)
  MESSAGE(STATUS "OpenMVG libraries found: ${OPENMVG_LIBRARY}")
  MESSAGE(STATUS "OpenMVG libraries directories: ${OPENMVG_LIBRARY_DIR}")
ENDIF (OPENMVG_LIBRARY)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set OPENMVG_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(OpenMVG  DEFAULT_MSG
                                  OPENMVG_LIBRARY OPENMVG_INCLUDE_DIR)

MARK_AS_ADVANCED(OPENMVG_INCLUDE_DIR OPENMVG_LIBRARY)

#Third parties:
# - include directories

IF(OPENMVG_FOUND)
  SET(OPENMVG_INCLUDE_DIRS
    ${OPENMVG_INCLUDE_DIR}
    ${OPENMVG_INCLUDE_DIR}/openMVG_third_party
    ${OPENMVG_INCLUDE_DIR}/openMVG_third_party/eigen
    #${OPENMVG_INCLUDE_DIR}/openMVG_third_party/lemon
    #${OPENMVG_INCLUDE_DIR}/openMVG_third_party/ceres-solver/include
    #${OPENMVG_INCLUDE_DIR}/openMVG_third_party/ceres-solver/internal/ceres/miniglog
    #${OPENMVG_INCLUDE_DIR}/openMVG_third_party/ceres-solver/config
    #${OPENMVG_INCLUDE_DIR}/openMVG_third_party/flann/src/cpp
  )
ENDIF(OPENMVG_FOUND)
