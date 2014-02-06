
# Copyright (c) 2012, 2013 openMVG authors.

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(openMVG C CXX)

# By default build in Release mode
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF(NOT CMAKE_BUILD_TYPE)

IF (WIN32)
  ADD_DEFINITIONS(-DNOMINMAX)
ENDIF (WIN32)

# ==============================================================================
# Additional cmake find modules
# ==============================================================================
SET(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmakeFindModules)
INCLUDE(OptimizeForArchitecture)
OptimizeForArchitecture()
IF (SSE2_FOUND)
  ADD_DEFINITIONS(-DUSE_SSE)
ENDIF (SSE2_FOUND)

# ==============================================================================
# OpenMP detection
# ==============================================================================
FIND_PACKAGE(OpenMP)
IF (OPENMP_FOUND)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)
  ADD_DEFINITIONS(-DUSE_OPENMP)
  IF (NOT MSVC)
    LIST(APPEND OPENMVG_LIBRARY_DEPENDENCIES gomp)
  ENDIF (NOT MSVC)
ELSE (OPENMP_FOUND)
  OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)
ENDIF (OPENMP_FOUND)

# ==============================================================================
# IMAGE IO detection
# ==============================================================================
IF (NOT APPLE)
  FIND_PACKAGE(JPEG QUIET)
  FIND_PACKAGE(PNG QUIET)
ENDIF (NOT APPLE)

# ==============================================================================
# SUBMODULE CONFIGURATION
# ==============================================================================
#- glfw
# ==============================================================================
SET(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "Do not build the GLFW example programs" FORCE)
SET(GLFW_BUILD_TESTS    OFF CACHE BOOL "Do not build the GLFW tests programs"   FORCE)
ADD_SUBDIRECTORY(dependencies/glfw)
INCLUDE_DIRECTORIES(SYSTEM dependencies/glfw/include)

#- OpenExif
# ==============================================================================
ADD_SUBDIRECTORY(dependencies/OpenExif)

#- EXIV 2
# ==============================================================================
IF (UNIX)
  FIND_PACKAGE( Exiv2 QUIET)
  IF (EXIV2_FOUND)
    ADD_DEFINITIONS(-DUSE_EXIV2)
  ENDIF (EXIV2_FOUND)
ENDIF (UNIX)

# ==============================================================================
# Opencv is not used by openMVG but some samples show how to use openCV
#  and openMVG simultaneously
# ==============================================================================
SET(USE_OPENCV OFF CACHE BOOL "Do not build opencv+openMVG samples programs" FORCE)
IF (USE_OPENCV)
# In order use openMVG+openCV samples set USE_OPENCV to ON and set the OpenCV_DIR
#  to the build directory of openCV. Prefer use static version of openCV libraries.
  SET(OpenCV_DIR "/home/pierre/Documents/Dev_OpenCV/opencv_Build")
ENDIF (USE_OPENCV)



# ==============================================================================
# Enable cmake UNIT TEST framework
# ==============================================================================
ENABLE_TESTING()
# MACRO to ease UnitTesting
MACRO (UNIT_TEST NAMESPACE NAME EXTRA_LIBS)
  ADD_EXECUTABLE(${NAMESPACE}_test_${NAME} ${NAME}_test.cpp)
  TARGET_LINK_LIBRARIES(${NAMESPACE}_test_${NAME}
                        ${EXTRA_LIBS} # Extra libs MUST be first.
                        CppUnitLite ${OPENMVG_LIBRARY_DEPENDENCIES})
  ADD_TEST(${NAMESPACE}_test_${NAME} ${NAMESPACE}_test_${NAME})
ENDMACRO (UNIT_TEST)

# Configure Eigen to use only MPL2 licensed code
ADD_DEFINITIONS(-DEIGEN_MPL2_ONLY)

# ==============================================================================
# Third-party libraries:
# ==============================================================================
ADD_SUBDIRECTORY(third_party)

# ==============================================================================
# Include directories
# ==============================================================================
INCLUDE_DIRECTORIES(./
  ./third_party/
  ${JPEG_INCLUDE_DIR}
  ${PNG_INCLUDE_DIRS}
  ./third_party/eigen
  ./third_party/lemon
  ${PROJECT_BINARY_DIR}/third_party/lemon
  ./third_party/ceres-solver/include
  ./third_party/ceres-solver/internal/ceres/miniglog
  ./third_party/flann/src/cpp
  ./third_party/JsonBox/include
)

# ==============================================================================
# openMVG modules
# ==============================================================================
# The openMVG library itself
ADD_SUBDIRECTORY(openMVG)

# software under patent
# Included for research purpose only
ADD_SUBDIRECTORY(patented)

# openMVG tutorial examples
ADD_SUBDIRECTORY(openMVG_Samples)

# Complete software build on openMVG
ADD_SUBDIRECTORY(software)

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

ELSE (EXISTS ${SPHINX_EXECUTABLE})
  MESSAGE("Sphinx need to be installed to generate the documentation")
ENDIF (EXISTS ${SPHINX_EXECUTABLE})

