project(Eigen)
cmake_minimum_required(VERSION 2.8.5)

# guard against in-source builds

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# guard against bad build-type strings

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(    NOT cmake_build_type_tolower STREQUAL "debug"
   AND NOT cmake_build_type_tolower STREQUAL "release"
   AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()


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

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Backward compatibility support for EIGEN_INCLUDE_INSTALL_DIR
if(EIGEN_INCLUDE_INSTALL_DIR AND NOT INCLUDE_INSTALL_DIR)
  set(INCLUDE_INSTALL_DIR ${EIGEN_INCLUDE_INSTALL_DIR}
      CACHE PATH "The directory relative to CMAKE_PREFIX_PATH where Eigen header files are installed")
else()
  set(INCLUDE_INSTALL_DIR
      "${CMAKE_INSTALL_INCLUDEDIR}/eigen3"
      CACHE PATH "The directory relative to CMAKE_PREFIX_PATH where Eigen header files are installed"
      )
endif()

set(CMAKEPACKAGE_INSTALL_DIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/eigen3"
    CACHE PATH "The directory relative to CMAKE_PREFIX_PATH where Eigen3Config.cmake is installed"
    )
set(PKGCONFIG_INSTALL_DIR
    "${CMAKE_INSTALL_DATADIR}/pkgconfig"
    CACHE PATH "The directory relative to CMAKE_PREFIX_PATH where eigen3.pc is installed"
    )

add_subdirectory(Eigen)

add_subdirectory(unsupported)

message(STATUS "")
message(STATUS "Configured Eigen ${EIGEN_VERSION_NUMBER}")
message(STATUS "")
