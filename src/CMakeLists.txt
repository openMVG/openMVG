
# Copyright (c) 2012, 2013 openMVG authors.

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(openMVG C CXX)

# guard against in-source builds
IF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed.")
ENDIF()

# ==============================================================================
# OpenMVG build options
# ==============================================================================
OPTION(OpenMVG_BUILD_SHARED "Build OpenMVG shared libs" OFF)
OPTION(OpenMVG_BUILD_TESTS "Build OpenMVG tests" OFF)
OPTION(OpenMVG_BUILD_DOC "Build OpenMVG documentation" ON)
OPTION(OpenMVG_BUILD_EXAMPLES "Build OpenMVG samples applications." ON)
OPTION(OpenMVG_BUILD_OPENGL_EXAMPLES "Build OpenMVG openGL examples" OFF)
OPTION(OpenMVG_BUILD_COVERAGE "Enable code coverage generation (gcc only)" OFF)
OPTION(OpenMVG_USE_OPENMP "Enable OpenMP parallelization" ON)
# ==============================================================================
# Opencv is not used by openMVG but some samples show how to use openCV
#  and openMVG simultaneously
# ==============================================================================
OPTION(OpenMVG_USE_OPENCV "Build or not opencv+openMVG samples programs" OFF)

# ==============================================================================
# Since OpenCV 3, SIFT is no longer in the default modules. See
# https://github.com/itseez/opencv_contrib for more informations.
# Enable this to be able to use OpenCV SIFT in main_ComputeFeatures_OpenCV.
# ==============================================================================
OPTION(OpenMVG_USE_OCVSIFT "Add or not OpenCV SIFT in available features" OFF)

# ==============================================================================
# OpenMVG version
# ==============================================================================
file(STRINGS "openMVG/version.hpp" _OPENMVG_VERSION_HPP_CONTENTS REGEX "#define OPENMVG_VERSION_")
foreach(v MAJOR MINOR REVISION)
  if("${_OPENMVG_VERSION_HPP_CONTENTS}" MATCHES "#define OPENMVG_VERSION_${v} ([0-9]+)")
    set(OPENMVG_VERSION_${v} "${CMAKE_MATCH_1}")
 else()
    message(FATAL_ERROR "Failed to retrieve the OpenMVG version from the source code. Missing OPENMVG_VERSION_${v}.")
  endif()
endforeach()
SET(OPENMVG_VERSION ${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}.${OPENMVG_VERSION_REVISION})

# Default build is in Release mode
IF(NOT CMAKE_BUILD_TYPE AND NOT MSVC)
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF(NOT CMAKE_BUILD_TYPE AND NOT MSVC)

# Set build path
SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_BUILD_TYPE}")
SET(LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_BUILD_TYPE}")

IF (WIN32)
  ADD_DEFINITIONS(-DNOMINMAX)
  ADD_DEFINITIONS(-DVC_EXTRALEAN)
  IF (MSVC)
    ADD_DEFINITIONS(/bigobj)
  ENDIF (MSVC)
ENDIF (WIN32)

# Set shared or static mode
IF(OpenMVG_BUILD_SHARED)
  SET(BUILD_SHARED_LIBS ON)
ELSEIF()
  SET(BUILD_SHARED_LIBS OFF)
ENDIF()

# ==============================================================================
# Check that submodule have been initialized and updated
# ==============================================================================
if(NOT EXISTS ${PROJECT_SOURCE_DIR}/dependencies/cereal/include)
  message(FATAL_ERROR
    "\n submodule(s) are missing, please update your repository:\n"
    "  > git submodule update -i\n")
endif()

# ==============================================================================
# Additional cmake find modules
# ==============================================================================
SET(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmakeFindModules)
INCLUDE(OptimizeForArchitecture)
OptimizeForArchitecture()
IF (SSE2_FOUND OR TARGET_ARCHITECTURE STREQUAL "native")
  IF (MSVC AND NOT ${CMAKE_CL_64})
    ADD_DEFINITIONS(/arch:SSE2)
  ENDIF (MSVC AND NOT ${CMAKE_CL_64})
  ADD_DEFINITIONS(-DOPENMVG_USE_SSE)
ENDIF ()
IF (UNIX AND NOT OpenMVG_BUILD_COVERAGE)
  SET(CMAKE_C_FLAGS_RELEASE "-O3")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
ENDIF ()

# ==============================================================================
# Check C++11 support
# ==============================================================================
INCLUDE(CXX11)
check_for_cxx11_compiler(CXX11_COMPILER)
# If a C++11 compiler is available, then set the appropriate flags
IF(CXX11_COMPILER)
  ADD_DEFINITIONS(-DOPENMVG_USE_CXX11)
  enable_cxx11()
  INCLUDE(CheckIncludeFileCXX)
  CHECK_INCLUDE_FILE_CXX(thread HAVE_CXX11_THREAD)
  IF (HAVE_CXX11_THREAD)
    ADD_DEFINITIONS(-DHAVE_CXX11_THREAD)
  ENDIF(HAVE_CXX11_THREAD)
    CHECK_INCLUDE_FILE_CXX(chrono HAVE_CXX11_CHRONO)
  IF (HAVE_CXX11_CHRONO)
    ADD_DEFINITIONS(-DHAVE_CXX11_CHRONO)
  ENDIF(HAVE_CXX11_CHRONO)
ENDIF(CXX11_COMPILER)

# ==============================================================================
# OpenMP detection
# ==============================================================================
IF(OpenMVG_USE_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    OPTION(OpenMVG_USE_OPENMP "Use OpenMP for parallelization" ON)
    ADD_DEFINITIONS(-DOPENMVG_USE_OPENMP)
    IF (NOT MSVC)
      IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        # for those using the clang with OpenMP support
        LIST(APPEND OPENMVG_LIBRARY_DEPENDENCIES omp)
      ELSE()
        LIST(APPEND OPENMVG_LIBRARY_DEPENDENCIES gomp)
      ENDIF()
    ENDIF(NOT MSVC)
  ENDIF(OPENMP_FOUND)
ELSE(OpenMVG_USE_OPENMP)
    OPTION(OpenMVG_USE_OPENMP "Use OpenMP for parallelization" OFF)
    UPDATE_CACHE_VARIABLE(OpenMVG_USE_OPENMP OFF)
    REMOVE_DEFINITIONS(-DOPENMVG_USE_OPENMP)
ENDIF(OpenMVG_USE_OPENMP)

# ==============================================================================
# enable code coverage generation (only with GCC)
# ==============================================================================
IF(OpenMVG_BUILD_COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
  MESSAGE("OpenMVG_BUILD_COVERAGE enabled")
  SET(CMAKE_BUILD_TYPE "Debug")
  ADD_DEFINITIONS(--coverage -fprofile-arcs -ftest-coverage)
  SET(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
ENDIF()

# ==============================================================================
# IMAGE IO detection
# ==============================================================================
IF (NOT APPLE)
  FIND_PACKAGE(JPEG QUIET)
  FIND_PACKAGE(PNG QUIET)
  FIND_PACKAGE(TIFF QUIET)
ENDIF (NOT APPLE)

# Folders
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

# ==============================================================================
# SUBMODULE CONFIGURATION
# ==============================================================================
#- glfw
# ==============================================================================
IF (OpenMVG_BUILD_OPENGL_EXAMPLES)
  SET(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "Do not build the GLFW example programs")
  SET(GLFW_BUILD_TESTS    OFF CACHE BOOL "Do not build the GLFW tests programs")
  SET(GLFW_BUILD_DOCS     OFF CACHE BOOL "Do not build the GLFW documentation")
  SET(GLFW_INSTALL        OFF CACHE BOOL "Do not generate the GLFW installation target")
  ADD_SUBDIRECTORY(dependencies/glfw)
  INCLUDE_DIRECTORIES(SYSTEM dependencies/glfw/include)
  SET_PROPERTY(TARGET glfw PROPERTY FOLDER OpenMVG/3rdParty/glfw)
ENDIF (OpenMVG_BUILD_OPENGL_EXAMPLES)

#- Mosek (linear programming interface)
#===============================================================================
FIND_PACKAGE(Mosek)
IF(MOSEK_FOUND)
  ADD_DEFINITIONS(-DOPENMVG_HAVE_MOSEK)
  SET(LP_INCLUDE_DIRS
    ${MOSEK_INCLUDE}
    ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Osi/src/OsiMsk/
  )
  INCLUDE_DIRECTORIES(${LP_INCLUDE_DIRS})
ENDIF(MOSEK_FOUND)

#INSTALL RULES
INSTALL(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/
  DESTINATION include/openMVG_dependencies
  COMPONENT headers
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)

# ==============================================================================
# --END-- SUBMODULE CONFIGURATION
# ==============================================================================

# ==============================================================================
# Enable cmake UNIT TEST framework
# ==============================================================================
IF(OpenMVG_BUILD_TESTS)
  ENABLE_TESTING()
ENDIF(OpenMVG_BUILD_TESTS)
# MACRO to ease UnitTesting
MACRO(UNIT_TEST NAMESPACE NAME EXTRA_LIBS)
  IF(OpenMVG_BUILD_TESTS)
    ADD_EXECUTABLE(${NAMESPACE}_test_${NAME} ${NAME}_test.cpp)

    SET_PROPERTY(TARGET ${NAMESPACE}_test_${NAME} PROPERTY FOLDER OpenMVG/test)

    TARGET_LINK_LIBRARIES(${NAMESPACE}_test_${NAME}
                          ${EXTRA_LIBS} # Extra libs MUST be first.
                          CppUnitLite ${OPENMVG_LIBRARY_DEPENDENCIES})
    ADD_TEST(NAME ${NAMESPACE}_test_${NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
             COMMAND $<TARGET_FILE:${NAMESPACE}_test_${NAME}>)
  ENDIF(OpenMVG_BUILD_TESTS)
ENDMACRO(UNIT_TEST)

# ==============================================================================
# Eigen
# ==============================================================================
# - internal by default,
# - external if EIGEN_INCLUDE_DIR_HINTS is defined
# ==============================================================================
IF(NOT DEFINED EIGEN_INCLUDE_DIR_HINTS)
  SET(EIGEN_INCLUDE_DIR_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen)
  SET(OpenMVG_USE_INTERNAL_EIGEN ON)
ENDIF()
FIND_PACKAGE(Eigen QUIET)
IF(EIGEN_FOUND)
  INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
ENDIF(EIGEN_FOUND)

# ==============================================================================
# Ceres
# ==============================================================================
# - internal by default (ceres-solver+cxsparse+miniglog),
# - external if CERES_DIR_HINTS and FIND_PACKAGE return a valid Ceres setup
# ==============================================================================
FIND_PACKAGE(Ceres QUIET HINTS ${CERES_DIR_HINTS})
IF(NOT Ceres_FOUND)
  SET(OpenMVG_USE_INTERNAL_CERES ON)
  SET(CERES_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/ceres-solver/include
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/ceres-solver/internal/ceres/miniglog
    ${PROJECT_BINARY_DIR}/third_party/ceres-solver/config)
    FILE(READ "${CMAKE_CURRENT_SOURCE_DIR}/third_party/ceres-solver/VERSION" CERES_CONFIG)
    STRING(REGEX REPLACE "version ([0-9.]+).*" "\\1" CERES_VERSION ${CERES_CONFIG})
  SET(CERES_LIBRARIES ceres cxsparse)
ENDIF()

# ==============================================================================
# Flann
# ==============================================================================
# - internal by default (flann),
# - external if FLANN_INCLUDE_DIR_HINTS and a valid Flann setup is found
# ==============================================================================
IF(NOT DEFINED FLANN_INCLUDE_DIR_HINTS)
  SET(FLANN_INCLUDE_DIR_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/flann/src/cpp)
  SET(OpenMVG_USE_INTERNAL_FLANN ON)
ENDIF()
FIND_PACKAGE(Flann QUIET)
IF(NOT FLANN_FOUND OR OpenMVG_USE_INTERNAL_FLANN)
  SET(FLANN_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/flann/src/cpp)
  SET(FLANN_LIBRARY flann_cpp_s)
ENDIF()

# ==============================================================================
# CoinUtils
# ==============================================================================
# - internal by default (CoinUtils),
# - external if COINUTILS_INCLUDE_DIR_HINTS and a valid CoinUtils setup is found
# ==============================================================================
IF(NOT DEFINED COINUTILS_INCLUDE_DIR_HINTS)
  SET(COINUTILS_INCLUDE_DIR_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/CoinUtils/src/)
  SET(OpenMVG_USE_INTERNAL_COINUTILS ON)
ENDIF()
FIND_PACKAGE(CoinUtils QUIET)
IF(NOT COINUTILS_FOUND OR OpenMVG_USE_INTERNAL_COINUTILS)
  SET(COINUTILS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/CoinUtils/src/)
  SET(COINUTILS_LIBRARY lib_CoinUtils)
ENDIF()

## ==============================================================================
## Clp
## ==============================================================================
## - internal by default (Clp),
## - external if CLP_INCLUDE_DIR_HINTS and a valid Clp setup is found
## ==============================================================================
IF(NOT DEFINED CLP_INCLUDE_DIR_HINTS)
  SET(CLP_INCLUDE_DIR_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Clp/src/)
  SET(OpenMVG_USE_INTERNAL_CLP ON)
ENDIF()
FIND_PACKAGE(Clp QUIET)
IF(NOT CLP_FOUND OR OpenMVG_USE_INTERNAL_CLP)
  SET(CLP_INCLUDE_DIRS
     ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Clp/src/
     ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Clp/src/OsiClp/)
  SET(CLP_LIBRARIES lib_clp lib_OsiClpSolver)
ENDIF()

# ==============================================================================
# Osi
# ==============================================================================
# - internal by default (Osi),
# - external if OSI_INCLUDE_DIR_HINTS and a valid Osi setup is found
# ==============================================================================
IF(NOT DEFINED OSI_INCLUDE_DIR_HINTS)
  SET(OSI_INCLUDE_DIR_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Osi/src/)
  SET(OpenMVG_USE_INTERNAL_OSI ON)
ENDIF()
FIND_PACKAGE(Osi QUIET)
IF(NOT OSI_FOUND OR OpenMVG_USE_INTERNAL_OSI)
  SET(OSI_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Osi/src/Osi/)
  SET(OSI_LIBRARY lib_Osi)
ENDIF()

# ==============================================================================
# Internal CLP/OSI/COINUTILS libraries:
# ==============================================================================
IF (OpenMVG_USE_INTERNAL_OSI AND OpenMVG_USE_INTERNAL_CLP AND OpenMVG_USE_INTERNAL_COINUTILS)
  ADD_SUBDIRECTORY(dependencies/osi_clp/)
  SET_PROPERTY(TARGET lib_clp PROPERTY FOLDER OpenMVG/3rdParty/lib_clp)
  SET_PROPERTY(TARGET lib_CoinUtils PROPERTY FOLDER OpenMVG/3rdParty/lib_CoinUtils)
  SET_PROPERTY(TARGET lib_Osi PROPERTY FOLDER OpenMVG/3rdParty/lib_Osi)
  SET_PROPERTY(TARGET lib_OsiClpSolver PROPERTY FOLDER OpenMVG/3rdParty/lib_OsiClpSolver)
ENDIF()

# ==============================================================================
# Lemon
# ==============================================================================
# - internal by default (Lemon),
# - external if LEMON_INCLUDE_DIR_HINTS and a valid Lemon setup is found
# ==============================================================================
IF(NOT DEFINED LEMON_INCLUDE_DIR_HINTS)
  SET(LEMON_INCLUDE_DIR_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lemon)
  SET(OpenMVG_USE_INTERNAL_LEMON ON)
ENDIF()
FIND_PACKAGE(Lemon QUIET)
IF(NOT LEMON_FOUND OR OpenMVG_USE_INTERNAL_LEMON)
  SET(LEMON_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lemon
    ${PROJECT_BINARY_DIR}/third_party/lemon)
  SET(LEMON_LIBRARY lemon)
ENDIF()

# ==============================================================================
# OpenCV
# ==============================================================================
# - only external and enabled only if OpenMVG_USE_OPENCV is set to ON
# ==============================================================================
IF(OpenMVG_USE_OPENCV)
  FIND_PACKAGE( OpenCV )
  IF(NOT OpenCV_FOUND)
    MESSAGE(STATUS "OpenCV was not found. Disabling OpenCV support.")
    UPDATE_CACHE_VARIABLE(OpenMVG_USE_OPENCV OFF)
    UPDATE_CACHE_VARIABLE(OpenMVG_USE_OCVSIFT OFF)
  ELSE(NOT OpenCV_FOUND)
    INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
    IF(OpenMVG_USE_OCVSIFT)
      ADD_DEFINITIONS(-DUSE_OCVSIFT)
    ENDIF(OpenMVG_USE_OCVSIFT)
    # add a definition that allows the conditional compiling
    ADD_DEFINITIONS(-DHAVE_OPENCV)
  ENDIF(NOT OpenCV_FOUND)
ENDIF(OpenMVG_USE_OPENCV)


# ==============================================================================
# Third-party libraries:
# ==============================================================================
ADD_SUBDIRECTORY(third_party)

# ==============================================================================
# Include directories
# ==============================================================================
SET(OpenMVG_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${JPEG_INCLUDE_DIR}
  ${PNG_INCLUDE_DIRS}
  ${TIFF_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party
  ${CMAKE_CURRENT_SOURCE_DIR}/dependencies
  ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/cereal/include
  ${LEMON_INCLUDE_DIRS}
  ${EIGEN_INCLUDE_DIRS}
  ${CERES_INCLUDE_DIRS}
  ${FLANN_INCLUDE_DIRS}
  ${LP_INCLUDE_DIRS}
  ${COINUTILS_INCLUDE_DIRS}
  ${CLP_INCLUDE_DIRS}
  ${OSI_INCLUDE_DIRS}
  CACHE PATH "include directories for openMVG and its dependencies"
)

INCLUDE_DIRECTORIES(${OpenMVG_INCLUDE_DIRS})

# ==============================================================================
# openMVG modules
# ==============================================================================
# The openMVG library itself
ADD_SUBDIRECTORY(openMVG)

# ==============================================================================
# openMVG libraries
# ==============================================================================
SET(OpenMVG_LIBRARIES
  openMVG_image
  openMVG_features
  openMVG_matching_image_collection
  openMVG_kvld
  openMVG_multiview
  openMVG_lInftyComputerVision
  openMVG_system
  openMVG_sfm)

FOREACH(omvglib ${OpenMVG_LIBRARIES})
  SET_PROPERTY(TARGET ${omvglib} PROPERTY FOLDER OpenMVG/OpenMVG)
ENDFOREACH()

# openMVG tutorial examples
IF (OpenMVG_BUILD_EXAMPLES)
  ADD_SUBDIRECTORY(openMVG_Samples)
ENDIF (OpenMVG_BUILD_EXAMPLES)

# Complete software(s) build on openMVG libraries
ADD_SUBDIRECTORY(software)

# Included for research purpose only
ADD_SUBDIRECTORY(nonFree)

# ==============================================================================
# Documentation
# --------------------------
# Sphinx detection
# ==============================================================================
IF (OpenMVG_BUILD_DOC)
  FIND_PACKAGE(Sphinx)
  IF (EXISTS ${SPHINX_EXECUTABLE})
    SET(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/htmlDoc")

    CONFIGURE_FILE(
      "${CMAKE_CURRENT_SOURCE_DIR}/../docs/sphinx/rst/conf.py"
      "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
      @ONLY)

    ADD_CUSTOM_TARGET(doc ALL
      ${SPHINX_EXECUTABLE}
      -b html
      "${CMAKE_CURRENT_SOURCE_DIR}/../docs/sphinx/rst"
      "${SPHINX_HTML_DIR}"
      COMMENT "Building HTML documentation with Sphinx")

    SET_PROPERTY(TARGET doc PROPERTY FOLDER OpenMVG)

  ELSE (EXISTS ${SPHINX_EXECUTABLE})
    MESSAGE("Sphinx need to be installed to generate the sphinx documentation")
  ENDIF (EXISTS ${SPHINX_EXECUTABLE})
  
  
  # ------
  # Doxygen detection
  # ======
  FIND_PACKAGE(Doxygen)
  IF (DOXYGEN_FOUND)
  
    CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/../docs/doxygen/Doxyfile.in
                    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY )

    ADD_CUSTOM_TARGET( doc_doxygen ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMENT "Generating API documentation with Doxygen" VERBATIM )
  
  ELSE (DOXYGEN_FOUND)
    MESSAGE("Doxygen need to be installed to generate the doxygen documentation")
  ENDIF (DOXYGEN_FOUND)
ENDIF()



MESSAGE("\n")

MESSAGE("** OpenMVG version: " ${OPENMVG_VERSION})
MESSAGE("** Build Shared libs: " ${OpenMVG_BUILD_SHARED})
MESSAGE("** Build OpenMVG tests: " ${OpenMVG_BUILD_TESTS})
MESSAGE("** Build OpenMVG documentation: " ${OpenMVG_BUILD_DOC})
MESSAGE("** Build OpenMVG samples applications: " ${OpenMVG_BUILD_EXAMPLES})
MESSAGE("** Build OpenMVG openGL examples: " ${OpenMVG_BUILD_OPENGL_EXAMPLES})
MESSAGE("** Enable code coverage generation: " ${OpenMVG_BUILD_COVERAGE})
MESSAGE("** Enable OpenMP parallelization: " ${OpenMVG_USE_OPENMP})
MESSAGE("** Build OpenCV+OpenMVG samples programs: " ${OpenMVG_USE_OPENCV})
MESSAGE("** Use OpenCV SIFT features: " ${OpenMVG_USE_OCVSIFT})

MESSAGE("\n")

IF(DEFINED OpenMVG_USE_INTERNAL_EIGEN)
  MESSAGE(STATUS "EIGEN: " ${EIGEN_VERSION} " (internal)")
ELSE()
  MESSAGE(STATUS "EIGEN: " ${EIGEN_VERSION} " (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_CERES)
  MESSAGE(STATUS "CERES: " ${CERES_VERSION} " (internal)")
ELSE()
  MESSAGE(STATUS "CERES: " ${CERES_VERSION} " (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_FLANN)
  MESSAGE(STATUS "FLANN: " ${FLANN_VERSION} " (internal)")
ELSE()
  MESSAGE(STATUS "FLANN: " ${FLANN_VERSION} " (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_TIFF)
  MESSAGE(STATUS "LIBTIFF: " ${TIFF_VERSION_STRING} " (internal)")
ELSE()
  MESSAGE(STATUS "LIBTIFF: " ${TIFF_VERSION_STRING} " (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_PNG)
  MESSAGE(STATUS "LIBPNG: " ${PNG_VERSION_STRING} " (internal)")
ELSE()
  MESSAGE(STATUS "LIBPNG: " ${PNG_VERSION_STRING} " (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_JPEG)
  MESSAGE(STATUS "LIBJPEG (internal)")
ELSE()
  MESSAGE(STATUS "LIBJPEG (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_CLP)
  MESSAGE(STATUS "CLP: " ${CLP_VERSION} " (internal)")
ELSE()
  MESSAGE(STATUS "CLP: " ${CLP_VERSION} " (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_COINUTILS)
  MESSAGE(STATUS "COINUTILS: " ${COINUTILS_VERSION} " (internal)")
ELSE()
  MESSAGE(STATUS "COINUTILS: " ${COINUTILS_VERSION} " (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_OSI)
  MESSAGE(STATUS "OSI: " ${OSI_VERSION} " (internal)")
ELSE()
  MESSAGE(STATUS "OSI: " ${OSI_VERSION} " (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_LEMON)
  MESSAGE(STATUS "LEMON: " ${LEMON_VERSION} " (internal)")
ELSE()
  MESSAGE(STATUS "LEMON: " ${LEMON_VERSION} " (external)")
ENDIF()

MESSAGE("\n")

# ==============================================================================
# INSTALL RULES
# ==============================================================================

INSTALL(EXPORT openMVG-targets DESTINATION lib)
INSTALL(EXPORT openMVG-targets
        DESTINATION share/openMVG/cmake FILE OpenMVGTargets.cmake)

#Adapt build include paths to install path
SET(OpenMVG_INCLUDE_DIRS
  "${OpenMVG_INCLUDE_DIRS}"
  "${CMAKE_INSTALL_PREFIX}/include/openMVG")

STRING(REGEX REPLACE
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_INSTALL_PREFIX}/include"
  OpenMVG_INCLUDE_DIRS
  "${OpenMVG_INCLUDE_DIRS}"
)
STRING(REGEX REPLACE
  "third_party"
  "openMVG/third_party"
  OpenMVG_INCLUDE_DIRS
  "${OpenMVG_INCLUDE_DIRS}"
)
STRING(REGEX REPLACE
  "dependencies"
  "openMVG_dependencies"
  OpenMVG_INCLUDE_DIRS
  "${OpenMVG_INCLUDE_DIRS}"
)
IF (OpenMVG_USE_INTERNAL_CERES)
  STRING(REGEX REPLACE
    "ceres-solver"
    "ceres"
    OpenMVG_INCLUDE_DIRS
    "${OpenMVG_INCLUDE_DIRS}"
  )
  STRING(REGEX REPLACE
    "internal/ceres/"
    ""
    OpenMVG_INCLUDE_DIRS
    "${OpenMVG_INCLUDE_DIRS}"
  )
ENDIF(OpenMVG_USE_INTERNAL_CERES)

# Create a OpenMVGConfig.cmake file. <name>Config.cmake files are searched by
# FIND_PACKAGE() automatically. We configure that file so that we can put any
# information we want in it, e.g. version numbers, include directories, etc.
CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/cmakeFindModules/OpenMVGConfig.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/OpenMVGConfig.cmake" @ONLY)

INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/OpenMVGConfig.cmake"
        DESTINATION share/openMVG/cmake)

