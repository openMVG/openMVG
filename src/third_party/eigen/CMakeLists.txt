project(Eigen)

cmake_minimum_required(VERSION 2.8.2)

# guard against in-source builds

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

IF(NOT CMAKE_BUILD_TYPE AND NOT MSVC)
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF(NOT CMAKE_BUILD_TYPE AND NOT MSVC)

#############################################################################
# retrieve version information                                               #
#############################################################################

# automatically parse the version number
file(READ "${PROJECT_SOURCE_DIR}/Eigen/src/Core/util/Macros.h" _eigen_version_header)
string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" _eigen_world_version_match "${_eigen_version_header}")
set(EIGEN_WORLD_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _eigen_major_version_match "${_eigen_version_header}")
set(EIGEN_MAJOR_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _eigen_minor_version_match "${_eigen_version_header}")
set(EIGEN_MINOR_VERSION "${CMAKE_MATCH_1}")
set(EIGEN_VERSION_NUMBER ${EIGEN_WORLD_VERSION}.${EIGEN_MAJOR_VERSION}.${EIGEN_MINOR_VERSION})

# if the mercurial program is absent, this will leave the EIGEN_HG_CHANGESET string empty,
# but won't stop CMake.
execute_process(COMMAND hg tip -R ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE EIGEN_HGTIP_OUTPUT)
execute_process(COMMAND hg branch -R ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE EIGEN_BRANCH_OUTPUT)

# if this is the default (aka development) branch, extract the mercurial changeset number from the hg tip output...
if(EIGEN_BRANCH_OUTPUT MATCHES "default")
string(REGEX MATCH "^changeset: *[0-9]*:([0-9;a-f]+).*" EIGEN_HG_CHANGESET_MATCH "${EIGEN_HGTIP_OUTPUT}")
set(EIGEN_HG_CHANGESET "${CMAKE_MATCH_1}")
endif(EIGEN_BRANCH_OUTPUT MATCHES "default")
#...and show it next to the version number
if(EIGEN_HG_CHANGESET)
  set(EIGEN_VERSION "${EIGEN_VERSION_NUMBER} (mercurial changeset ${EIGEN_HG_CHANGESET})")
else(EIGEN_HG_CHANGESET)
  set(EIGEN_VERSION "${EIGEN_VERSION_NUMBER}")
endif(EIGEN_HG_CHANGESET)

# the user modifiable install path for header files
set(EIGEN_INCLUDE_INSTALL_DIR ${EIGEN_INCLUDE_INSTALL_DIR} CACHE PATH "The directory where we install the header files (optional)")

# set the internal install path for header files which depends on whether the user modifiable
# EIGEN_INCLUDE_INSTALL_DIR has been set by the user or not.
if(EIGEN_INCLUDE_INSTALL_DIR)
  set(INCLUDE_INSTALL_DIR
    ${EIGEN_INCLUDE_INSTALL_DIR}
    CACHE INTERNAL
    "The directory where we install the header files (internal)"
  )
else()
  set(INCLUDE_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/include/eigen3"
    CACHE INTERNAL
    "The directory where we install the header files (internal)"
  )
endif()

SET(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_subdirectory(Eigen)

add_subdirectory(unsupported)

message(STATUS "")
message(STATUS "Configured Eigen ${EIGEN_VERSION_NUMBER}")
message(STATUS "")
