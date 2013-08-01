
# Copyright (c) 2012, 2013 openMVG authors.

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(openMVG C CXX)

# By default build in Release mode
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -msse -msse2 -msse3")

# ==============================================================================
# Additional cmake find modules
# ==============================================================================
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmakeFindModules)
include (OptimizeForArchitecture)
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
  option(USE_OPENMP "Use OpenMP for parallelization" ON)
  ADD_DEFINITIONS(-DUSE_OPENMP)
  IF (NOT MSVC)
    LIST(APPEND OPENMVG_LIBRARY_DEPENDENCIES gomp)
  ENDIF (NOT MSVC)
ELSE (OPENMP_FOUND)
  option(USE_OPENMP "Use OpenMP for parallelization" OFF)
ENDIF (OPENMP_FOUND)

# ==============================================================================
# IMAGE IO detection
# ==============================================================================
FIND_PACKAGE(JPEG QUIET)
FIND_PACKAGE(PNG QUIET)

# ==============================================================================
# glfw
# ==============================================================================
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "Do not build the GLFW example programs" FORCE)
set(GLFW_BUILD_TESTS    OFF CACHE BOOL "Do not build the GLFW tests programs"   FORCE)
add_subdirectory(dependencies/glfw)

include_directories(SYSTEM dependencies/glfw/include)
include_directories(SYSTEM dependencies/glm)

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

