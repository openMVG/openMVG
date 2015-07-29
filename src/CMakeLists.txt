
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

# build path
SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_BUILD_TYPE}")
SET(LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_BUILD_TYPE}")

# ==============================================================================
# OpenMVG build options
# ==============================================================================
OPTION(OpenMVG_BUILD_SHARED "Build OpenMVG shared libs" OFF)
OPTION(OpenMVG_BUILD_TESTS "Build OpenMVG tests" OFF)
OPTION(OpenMVG_BUILD_EXAMPLES "Build OpenMVG samples applications." ON)
OPTION(OpenMVG_BUILD_OPENGL_EXAMPLES "Build OpenMVG openGL examples" OFF)
OPTION(OpenMVG_BUILD_COVERAGE "Enable code coverage generation (gcc only)" OFF)
OPTION(OpenMVG_USE_OPENMP "Enable OpenMP parallelization" ON)

SET(OPENMVG_VERSION_MAJOR 0)
SET(OPENMVG_VERSION_MINOR 8)
SET(OPENMVG_VERSION_PATCH 1)
SET(OPENMVG_VERSION
    ${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}.${OPENMVG_VERSION_PATCH})

# Default build is in Release mode
IF(NOT CMAKE_BUILD_TYPE AND NOT MSVC)
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF(NOT CMAKE_BUILD_TYPE AND NOT MSVC)

IF (WIN32)
  ADD_DEFINITIONS(-DNOMINMAX)
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
# Additional cmake find modules
# ==============================================================================
SET(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmakeFindModules)
INCLUDE(OptimizeForArchitecture)
OptimizeForArchitecture()
IF (SSE2_FOUND)
  IF (MSVC AND NOT ${CMAKE_CL_64})
    ADD_DEFINITIONS(/arch:SSE2)
  ENDIF (MSVC AND NOT ${CMAKE_CL_64})
  ADD_DEFINITIONS(-DOPENMVG_USE_SSE)
ENDIF (SSE2_FOUND)
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
IF (OpenMVG_USE_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF (OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    OPTION(OpenMVG_USE_OPENMP "Use OpenMP for parallelization" ON)
    ADD_DEFINITIONS(-DOPENMVG_USE_OPENMP)
    IF (NOT MSVC)
      LIST(APPEND OPENMVG_LIBRARY_DEPENDENCIES gomp)
    ENDIF (NOT MSVC)
ELSE (OpenMVG_USE_OPENMP)
    OPTION(OpenMVG_USE_OPENMP "Use OpenMP for parallelization" OFF)
    SET(OpenMVG_USE_OPENMP OFF CACHE BOOL "Use OpenMP for parallelization")
    REMOVE_DEFINITIONS(-DOPENMVG_USE_OPENMP)
  ENDIF (OPENMP_FOUND)
ENDIF (OpenMVG_USE_OPENMP)

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
# Enable cmake UNIT TEST framework
# ==============================================================================
IF (OpenMVG_BUILD_TESTS)
  ENABLE_TESTING()
ENDIF (OpenMVG_BUILD_TESTS)
# MACRO to ease UnitTesting
MACRO (UNIT_TEST NAMESPACE NAME EXTRA_LIBS)
  IF (OpenMVG_BUILD_TESTS)
    ADD_EXECUTABLE(${NAMESPACE}_test_${NAME} ${NAME}_test.cpp)

    SET_PROPERTY(TARGET ${NAMESPACE}_test_${NAME} PROPERTY FOLDER OpenMVG/test)

    TARGET_LINK_LIBRARIES(${NAMESPACE}_test_${NAME}
                          ${EXTRA_LIBS} # Extra libs MUST be first.
                          CppUnitLite ${OPENMVG_LIBRARY_DEPENDENCIES})
    ADD_TEST(NAME ${NAMESPACE}_test_${NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
             COMMAND $<TARGET_FILE:${NAMESPACE}_test_${NAME}>)
  ENDIF (OpenMVG_BUILD_TESTS)
ENDMACRO (UNIT_TEST)

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
# Configure Eigen to use only MPL2 licensed code
ADD_DEFINITIONS(-DEIGEN_MPL2_ONLY)

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
# - external if FLANN_DIR_HINTS and a valid Flann setup is found
# ==============================================================================
FIND_PATH(FLANN_INCLUDE_DIRS flann.hpp ${FLANN_DIR_HINTS})
IF(FLANN_INCLUDE_DIRS)
  FILE(READ "${FLANN_INCLUDE_DIRS}/config.h" FLANN_CONFIG)
  STRING(REGEX REPLACE ".*FLANN_VERSION_.+\"([0-9.]+)\".*" "\\1" FLANN_VERSION ${FLANN_CONFIG})
  FIND_LIBRARY(FLANN_LIBRARY NAMES flann_cpp)
  IF(NOT DEFINED FLANN_LIBRARY)
    SET(OpenMVG_USE_INTERNAL_FLANN ON)
  ENDIF()
ELSE(FLANN_INCLUDE_DIRS)
  SET(OpenMVG_USE_INTERNAL_FLANN ON)
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_FLANN)
  SET(FLANN_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/flann/src/cpp)
  FILE(READ "${FLANN_INCLUDE_DIRS}/flann/config.h" FLANN_CONFIG)
  STRING(REGEX REPLACE ".*FLANN_VERSION_.+\"([0-9.]+)\".*" "\\1" FLANN_VERSION ${FLANN_CONFIG})
  SET(FLANN_LIBRARY flann_cpp_s)
ENDIF()

# ==============================================================================
# CoinUtils
# ==============================================================================
# - internal by default (CoinUtils),
# - external if COINUTILS_DIR_HINTS and a valid CoinUtils setup is found
# ==============================================================================
FIND_PATH(COINUTILS_INCLUDE_DIRS CoinUtilsConfig.h ${COINUTILS_DIR_HINTS})
IF(COINUTILS_INCLUDE_DIRS)
  FILE(READ "${COINUTILS_INCLUDE_DIRS}/CoinUtilsConfig.h" COINUTILS_CONFIG)
  STRING(REGEX REPLACE ".*COINUTILS_VERSION \"([0-9.]+)\".*" "\\1" COINUTILS_VERSION ${COINUTILS_CONFIG})
  FIND_LIBRARY(COINUTILS_LIBRARY NAMES CoinUtils)
  IF(NOT DEFINED COINUTILS_LIBRARY)
    SET(OpenMVG_USE_INTERNAL_COINUTILS ON)
  ENDIF()
ELSE(COINUTILS_INCLUDE_DIRS)
  SET(OpenMVG_USE_INTERNAL_COINUTILS ON)
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_COINUTILS)
  SET(COINUTILS_INCLUDE_DIRS
     ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/CoinUtils/src/)
  SET(COINUTILS_LIBRARY lib_CoinUtils)
ENDIF()

# ==============================================================================
# Clp
# ==============================================================================
# - internal by default (Clp),
# - external if CLP_DIR_HINTS and a valid Clp setup is found
# ==============================================================================
FIND_PATH(CLP_INCLUDE_DIRS ClpConfig.h ${CLP_DIR_HINTS})
IF(CLP_INCLUDE_DIRS)
  FILE(READ "${CLP_INCLUDE_DIRS}/ClpConfig.h" CLP_CONFIG)
  STRING(REGEX REPLACE ".*CLP_VERSION \"([0-9.]+)\".*" "\\1" CLP_VERSION ${CLP_CONFIG})
  FIND_LIBRARY(CLP_LIBRARY NAMES Clp)
  FIND_LIBRARY(CLPSOLVER_LIBRARY NAMES ClpSolver)
  FIND_LIBRARY(OSICLP_LIBRARY NAMES OsiClp)
  IF(NOT DEFINED CLP_LIBRARY AND NOT DEFINED CLP_SOLVER_LIBRARY AND NOT DEFINED OSICLP_LIBRARY)
    SET(OpenMVG_USE_INTERNAL_CLP ON)
  ELSE()
    SET(CLP_LIBRARIES ${CLP_LIBRARY} ${CLPSOLVER_LIBRARY} ${OSICLP_LIBRARY})
  ENDIF()
ELSE(CLP_INCLUDE_DIRS)
  SET(OpenMVG_USE_INTERNAL_CLP ON)
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_CLP)
  SET(CLP_INCLUDE_DIRS
     ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Clp/src/
     ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Clp/src/OsiClp/)
  SET(CLP_LIBRARIES lib_clp lib_OsiClpSolver)
ENDIF()

# ==============================================================================
# Osi
# ==============================================================================
# - internal by default (Osi),
# - external if OSI_DIR_HINTS and a valid Osi setup is found
# ==============================================================================
FIND_PATH(OSI_INCLUDE_DIRS OsiConfig.h ${OSI_DIR_HINTS})
IF(OSI_INCLUDE_DIRS)
  FILE(READ "${OSI_INCLUDE_DIRS}/OsiConfig.h" OSI_CONFIG)
  STRING(REGEX REPLACE ".*OSI_VERSION \"([0-9.]+)\".*" "\\1" OSI_VERSION ${OSI_CONFIG})
  FIND_LIBRARY(OSI_LIBRARY NAMES Osi)
  IF(NOT DEFINED OSI_LIBRARY)
    SET(OpenMVG_USE_INTERNAL_OSI ON)
  ENDIF()
ELSE(OSI_INCLUDE_DIRS)
  SET(OpenMVG_USE_INTERNAL_OSI ON)
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_OSI)
  SET(OSI_INCLUDE_DIRS
     ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Osi/src/Osi/)
  SET(OSI_LIBRARY lib_Osi)
ENDIF()

# ==============================================================================
# Lemon
# ==============================================================================
# - internal by default (Lemon),
# - external if LEMON_DIR_HINTS and a valid Lemon setup is found
# ==============================================================================
FIND_PATH(LEMON_INCLUDE_DIR config.h ${LEMON_DIR_HINTS})
IF(LEMON_INCLUDE_DIR)
  FILE(READ "${LEMON_INCLUDE_DIR}/config.h" LEMON_CONFIG)
  STRING(REGEX REPLACE ".*LEMON_VERSION \"([0-9.]+)\".*" "\\1" LEMON_VERSION ${LEMON_CONFIG})
  FIND_LIBRARY(LEMON_LIBRARY NAMES emon)
  IF(NOT DEFINED LEMON_LIBRARY)
    SET(OpenMVG_USE_INTERNAL_LEMON ON)
  ELSE()
    GET_FILENAME_COMPONENT(LEMON_INCLUDE_DIRS ${LEMON_INCLUDE_DIR} PATH)
  ENDIF()
ELSE(LEMON_INCLUDE_DIR)
  SET(OpenMVG_USE_INTERNAL_LEMON ON)
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_LEMON)
  SET(LEMON_INCLUDE_DIRS
     ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lemon
     ${PROJECT_BINARY_DIR}/third_party/lemon)
  SET(LEMON_LIBRARY lemon)
ENDIF()

# ==============================================================================
# Internal CLP/OSI/COINUTILS libraries:
# ==============================================================================
ADD_SUBDIRECTORY(dependencies/osi_clp/)


# ==============================================================================
# Third-party libraries:
# ==============================================================================
ADD_SUBDIRECTORY(third_party)

# ==============================================================================
# Include directories
# ==============================================================================
SET(OpenMVG_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party
  ${JPEG_INCLUDE_DIR}
  ${PNG_INCLUDE_DIRS}
  ${TIFF_INCLUDE_DIR}
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

# software(s) under patent or commercial licence
# Included for research purpose only
ADD_SUBDIRECTORY(nonFree)

# ==============================================================================
# Documentation
# --------------------------
# Sphinx detection
# ==============================================================================
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
  MESSAGE("Sphinx need to be installed to generate the documentation")
ENDIF (EXISTS ${SPHINX_EXECUTABLE})


MESSAGE("\n")

MESSAGE("** OpenMVG version: " ${OPENMVG_VERSION})
MESSAGE("** Build Shared libs: " ${OpenMVG_BUILD_SHARED})
MESSAGE("** Build OpenMVG tests: " ${OpenMVG_BUILD_TESTS})
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
  MESSAGE(STATUS "LIBPNG (internal)")
ELSE()
  MESSAGE(STATUS "LIBPNG (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_CLP)
  MESSAGE(STATUS "CLP (internal)")
ELSE()
  MESSAGE(STATUS "CLP: " ${CLP_VERSION} " (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_COINUTILS)
  MESSAGE(STATUS "COINUTILS (internal)")
ELSE()
  MESSAGE(STATUS "COINUTILS: " ${COINUTILS_VERSION} " (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_OSI)
  MESSAGE(STATUS "OSI (internal)")
ELSE()
  MESSAGE(STATUS "OSI: " ${OSI_VERSION} " (external)")
ENDIF()

IF(DEFINED OpenMVG_USE_INTERNAL_LEMON)
  MESSAGE(STATUS "LEMON (internal)")
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

