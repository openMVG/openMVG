
# Copyright (c) 2012, 2013 openMVG authors.

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

project(openMVG C CXX)

# guard against in-source builds
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed.")
endif()

# ==============================================================================
# OpenMVG build options
# ==============================================================================
option(OpenMVG_BUILD_SHARED "Build OpenMVG shared libs" OFF)
option(OpenMVG_BUILD_TESTS "Build OpenMVG tests" OFF)
option(OpenMVG_BUILD_DOC "Build OpenMVG documentation" ON)
option(OpenMVG_BUILD_EXAMPLES "Build OpenMVG samples applications." ON)
option(OpenMVG_BUILD_OPENGL_EXAMPLES "Build OpenMVG openGL examples" OFF)
option(OpenMVG_BUILD_SOFTWARES "Build OpenMVG softwares" ON)
option(OpenMVG_BUILD_GUI_SOFTWARES "Build OpenMVG GUI softwares (QT5)" ON)
option(OpenMVG_BUILD_COVERAGE "Enable code coverage generation (gcc only)" OFF)
option(OpenMVG_USE_OPENMP "Enable OpenMP parallelization" ON)
# ==============================================================================
# Opencv is not used by openMVG but some samples show how to use openCV
#  and openMVG simultaneously
# ==============================================================================
option(OpenMVG_USE_OPENCV "Build or not opencv+openMVG samples programs" OFF)

# ==============================================================================
# Since OpenCV 3, SIFT is no longer in the default modules. See
# https://github.com/itseez/opencv_contrib for more informations.
# Enable this to be able to use OpenCV SIFT in main_ComputeFeatures_OpenCV.
# ==============================================================================
option(OpenMVG_USE_OCVSIFT "Add or not OpenCV SIFT in available features" OFF)

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
set(OPENMVG_VERSION ${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}.${OPENMVG_VERSION_REVISION})

# Default build is in Release mode
if (NOT CMAKE_BUILD_TYPE AND NOT MSVC)
  set(CMAKE_BUILD_TYPE "Release")
endif (NOT CMAKE_BUILD_TYPE AND NOT MSVC)

# Set build path
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_BUILD_TYPE}")
set(LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_BUILD_TYPE}")

# ==============================================================================
# Macro used to forward compiler preprocessors to OpenMVGConfig.cmake targets
# ==============================================================================
set(OpenMVG_DEFINITIONS "")
MACRO(register_definitions DEF)
  add_definitions(${DEF})
  string(REPLACE "-D" "" STRIPPED_DEF ${DEF})
  list(APPEND OpenMVG_DEFINITIONS ${STRIPPED_DEF})
ENDMACRO()

if (WIN32)
  register_definitions(-DNOMINMAX)
  register_definitions(-DVC_EXTRALEAN)
  register_definitions(-D_USE_MATH_DEFINES)
  if (MSVC)
    register_definitions(/bigobj)
  endif (MSVC)
endif (WIN32)

# Set shared or static mode
if (OpenMVG_BUILD_SHARED)
  set(BUILD_SHARED_LIBS ON)
else (OpenMVG_BUILD_SHARED)
  set(BUILD_SHARED_LIBS OFF)
endif()

# For both regular Clang and AppleClang
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  message("clang detected: using std::unordered_map for Hash_Map")
  register_definitions(-DOPENMVG_STD_UNORDERED_MAP)
endif()

# ==============================================================================
# Check that submodule have been initialized and updated
# ==============================================================================
if (NOT EXISTS ${PROJECT_SOURCE_DIR}/dependencies/cereal/include)
  message(FATAL_ERROR
    "\n submodule(s) are missing, please update your repository:\n"
    "  > git submodule update -i\n")
endif()

# ==============================================================================
# Additional cmake find modules
# ==============================================================================
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmakeFindModules)
include(OptimizeForArchitecture)
AutodetectHostArchitecture()
OptimizeForArchitecture()

if (USE_AVX2)
  register_definitions(-DOPENMVG_USE_AVX2)
endif ()
if (UNIX AND NOT OpenMVG_BUILD_COVERAGE)
  set(CMAKE_C_FLAGS_RELEASE "-O3")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif ()

# ==============================================================================
# Check C++11 support
# ==============================================================================
include(CXX11)
check_for_cxx11_compiler(CXX11_COMPILER)
# If a C++11 compiler is available, then set the appropriate flags
if (CXX11_COMPILER)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
else (CXX11_COMPILER)
  message(FATAL_ERROR "A CXX11 compiler is required.")
endif (CXX11_COMPILER)

# ==============================================================================
# OpenMP detection
# ==============================================================================
if (OpenMVG_USE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    option(OpenMVG_USE_OPENMP "Use OpenMP for parallelization" ON)
    register_definitions(-DOPENMVG_USE_OPENMP)
    if (NOT MSVC)
      IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        # for those using the clang with OpenMP support
        list(APPEND OPENMVG_LIBRARY_DEPENDENCIES omp)
      ELSE()
        list(APPEND OPENMVG_LIBRARY_DEPENDENCIES gomp)
      endif()
    endif (NOT MSVC)
  endif (OPENMP_FOUND)
else (OpenMVG_USE_OPENMP)
    option(OpenMVG_USE_OPENMP "Use OpenMP for parallelization" OFF)
    include(UpdateCacheVariable)
    UPDATE_CACHE_VARIABLE(OpenMVG_USE_OPENMP OFF)
    remove_definitions(-DOPENMVG_USE_OPENMP)
endif (OpenMVG_USE_OPENMP)

# ==============================================================================
# enable code coverage generation (only with GCC)
# ==============================================================================
if (OpenMVG_BUILD_COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
  message("OpenMVG_BUILD_COVERAGE enabled")
  set(CMAKE_BUILD_TYPE "Debug")
  add_definitions(--coverage -fprofile-arcs -ftest-coverage)
  set(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
endif()

# ==============================================================================
# IMAGE IO detection
# ==============================================================================
find_package(JPEG QUIET)
find_package(PNG QUIET)
find_package(TIFF QUIET)

# Folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# ==============================================================================
# SUBMODULE CONFIGURATION
# ==============================================================================
#- glfw
# ==============================================================================
if (OpenMVG_BUILD_OPENGL_EXAMPLES)
  set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "Do not build the GLFW example programs")
  set(GLFW_BUILD_TESTS    OFF CACHE BOOL "Do not build the GLFW tests programs")
  set(GLFW_BUILD_DOCS     OFF CACHE BOOL "Do not build the GLFW documentation")
  set(GLFW_INSTALL        OFF CACHE BOOL "Do not generate the GLFW installation target")
  add_subdirectory(dependencies/glfw)
  include_directories(SYSTEM dependencies/glfw/include)
  set_property(TARGET glfw PROPERTY FOLDER OpenMVG/3rdParty/glfw)
endif (OpenMVG_BUILD_OPENGL_EXAMPLES)

#INSTALL RULES
install(
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
if (OpenMVG_BUILD_TESTS)
  enable_testing()
endif (OpenMVG_BUILD_TESTS)
# MACRO to ease UnitTesting
MACRO(UNIT_TEST NAMESPACE NAME EXTRA_LIBS)
  if (OpenMVG_BUILD_TESTS)
    add_executable(${NAMESPACE}_test_${NAME} ${NAME}_test.cpp)

    set_property(TARGET ${NAMESPACE}_test_${NAME} PROPERTY FOLDER OpenMVG/test)
    target_include_directories(${NAMESPACE}_test_${NAME}
                               PRIVATE
                               ${CMAKE_SOURCE_DIR}/third_party)
    target_link_libraries(${NAMESPACE}_test_${NAME}
                          ${EXTRA_LIBS} # Extra libs MUST be first.
                          CppUnitLite ${OPENMVG_LIBRARY_DEPENDENCIES})
    add_test(NAME ${NAMESPACE}_test_${NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
             COMMAND $<TARGET_FILE:${NAMESPACE}_test_${NAME}>)
  target_compile_definitions(${NAMESPACE}_test_${NAME}
    PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  endif (OpenMVG_BUILD_TESTS)
ENDMACRO(UNIT_TEST)

# ==============================================================================
# Eigen
# ==============================================================================
# - internal by default,
# - external if EIGEN_INCLUDE_DIR_HINTS is defined
# ==============================================================================
if (NOT DEFINED EIGEN_INCLUDE_DIR_HINTS)
  set(EIGEN_INCLUDE_DIR_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen)
  set(OpenMVG_USE_INTERNAL_EIGEN ON)
endif()
find_package(Eigen QUIET)
if (EIGEN_FOUND)
  include_directories(${EIGEN_INCLUDE_DIRS})
endif (EIGEN_FOUND)

# ==============================================================================
# Ceres
# ==============================================================================
# - internal by default (ceres-solver+cxsparse+miniglog),
# - external if CERES_DIR_HINTS and find_package return a valid Ceres setup
# ==============================================================================
find_package(Ceres QUIET HINTS ${CERES_DIR_HINTS})
if (NOT Ceres_FOUND)
  set(OpenMVG_USE_INTERNAL_CERES ON)
  set(CERES_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/ceres-solver/include
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/ceres-solver/internal/ceres/miniglog
    ${PROJECT_BINARY_DIR}/third_party/ceres-solver/config)
    FILE(READ "${CMAKE_CURRENT_SOURCE_DIR}/third_party/ceres-solver/VERSION" CERES_CONFIG)
    STRING(REGEX REPLACE "version ([0-9.]+).*" "\\1" CERES_VERSION ${CERES_CONFIG})
  set(CERES_LIBRARIES ceres cxsparse)
endif()

# ==============================================================================
# Flann
# ==============================================================================
# - internal by default (flann),
# - external if FLANN_INCLUDE_DIR_HINTS and a valid Flann setup is found
# ==============================================================================
if (NOT DEFINED FLANN_INCLUDE_DIR_HINTS)
  set(FLANN_INCLUDE_DIR_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/flann/src/cpp)
  set(OpenMVG_USE_INTERNAL_FLANN ON)
endif()
find_package(Flann QUIET)
if (NOT FLANN_FOUND OR OpenMVG_USE_INTERNAL_FLANN)
  set(FLANN_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/flann/src/cpp)
endif()

# ==============================================================================
# CoinUtils
# ==============================================================================
# - internal by default (CoinUtils),
# - external if COINUTILS_INCLUDE_DIR_HINTS and a valid CoinUtils setup is found
# ==============================================================================
if (NOT DEFINED COINUTILS_INCLUDE_DIR_HINTS)
  set(COINUTILS_INCLUDE_DIR_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/CoinUtils/src/)
  set(OpenMVG_USE_INTERNAL_COINUTILS ON)
endif()
find_package(CoinUtils QUIET)
if (NOT COINUTILS_FOUND OR OpenMVG_USE_INTERNAL_COINUTILS)
  set(COINUTILS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/CoinUtils/src/)
  set(COINUTILS_LIBRARY lib_CoinUtils)
endif()

## ==============================================================================
## Clp
## ==============================================================================
## - internal by default (Clp),
## - external if CLP_INCLUDE_DIR_HINTS and a valid Clp setup is found
## ==============================================================================
if (NOT DEFINED CLP_INCLUDE_DIR_HINTS)
  set(CLP_INCLUDE_DIR_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Clp/src/)
  set(OpenMVG_USE_INTERNAL_CLP ON)
endif()
find_package(Clp QUIET)
if (NOT CLP_FOUND OR OpenMVG_USE_INTERNAL_CLP)
  set(CLP_INCLUDE_DIRS
     ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Clp/src/
     ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Clp/src/OsiClp/)
  set(CLP_LIBRARIES lib_clp lib_OsiClpSolver)
endif()

# ==============================================================================
# Osi
# ==============================================================================
# - internal by default (Osi),
# - external if OSI_INCLUDE_DIR_HINTS and a valid Osi setup is found
# ==============================================================================
if (NOT DEFINED OSI_INCLUDE_DIR_HINTS)
  set(OSI_INCLUDE_DIR_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Osi/src/)
  set(OpenMVG_USE_INTERNAL_OSI ON)
endif()
find_package(Osi QUIET)
if (NOT OSI_FOUND OR OpenMVG_USE_INTERNAL_OSI)
  set(OSI_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/osi_clp/Osi/src/Osi/)
  set(OSI_LIBRARY lib_Osi)
endif()

# ==============================================================================
# Internal CLP/OSI/COINUTILS libraries:
# ==============================================================================
if (OpenMVG_USE_INTERNAL_OSI AND OpenMVG_USE_INTERNAL_CLP AND OpenMVG_USE_INTERNAL_COINUTILS)
  add_subdirectory(dependencies/osi_clp/)
  set_property(TARGET lib_clp PROPERTY FOLDER OpenMVG/3rdParty/lib_clp)
  set_property(TARGET lib_CoinUtils PROPERTY FOLDER OpenMVG/3rdParty/lib_CoinUtils)
  set_property(TARGET lib_Osi PROPERTY FOLDER OpenMVG/3rdParty/lib_Osi)
  set_property(TARGET lib_OsiClpSolver PROPERTY FOLDER OpenMVG/3rdParty/lib_OsiClpSolver)
endif()

# ==============================================================================
# Lemon
# ==============================================================================
# - internal by default (Lemon),
# - external if LEMON_INCLUDE_DIR_HINTS and a valid Lemon setup is found
# ==============================================================================
if (NOT DEFINED LEMON_INCLUDE_DIR_HINTS)
  set(LEMON_INCLUDE_DIR_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lemon)
  set(OpenMVG_USE_INTERNAL_LEMON ON)
endif()
find_package(Lemon QUIET)
if (NOT LEMON_FOUND OR OpenMVG_USE_INTERNAL_LEMON)
  set(LEMON_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lemon
    ${PROJECT_BINARY_DIR}/third_party/lemon)
  set(LEMON_LIBRARY lemon)
endif()

# ==============================================================================
# OpenCV
# ==============================================================================
# - only external and enabled only if OpenMVG_USE_OPENCV is set to ON
# ==============================================================================
if (OpenMVG_USE_OPENCV)
  find_package( OpenCV 3.0 QUIET )
  if (NOT OpenCV_FOUND)
    message(STATUS "OpenCV was not found. Disabling OpenCV support.")
    UPDATE_CACHE_VARIABLE(OpenMVG_USE_OPENCV OFF)
    UPDATE_CACHE_VARIABLE(OpenMVG_USE_OCVSIFT OFF)
  endif (NOT OpenCV_FOUND)
endif (OpenMVG_USE_OPENCV)


# ==============================================================================
# Third-party libraries:
# ==============================================================================
add_subdirectory(third_party)

# ==============================================================================
# Include directories
# ==============================================================================
set(OpenMVG_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/dependencies
  ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/cereal/include
  ${LEMON_INCLUDE_DIRS}
  ${EIGEN_INCLUDE_DIRS}
  ${FLANN_INCLUDE_DIRS}
  CACHE PATH "include directories for openMVG and its dependencies"
)

include_directories(${OpenMVG_INCLUDE_DIRS})

# ==============================================================================
# openMVG modules
# ==============================================================================
# The openMVG library itself
add_subdirectory(openMVG)

# ==============================================================================
# openMVG libraries
# ==============================================================================
set(OpenMVG_LIBRARIES
  openMVG_image
  openMVG_features
  openMVG_geometry
  openMVG_matching
  openMVG_matching_image_collection
  openMVG_kvld
  openMVG_multiview
  openMVG_lInftyComputerVision
  openMVG_system
  openMVG_sfm)

foreach(omvglib ${OpenMVG_LIBRARIES})
  set_property(TARGET ${omvglib} PROPERTY FOLDER OpenMVG/OpenMVG)
  set_property(TARGET ${omvglib} APPEND PROPERTY
    INTERFACE_COMPILE_DEFINITIONS ${OpenMVG_DEFINITIONS})
endforeach()

# openMVG tutorial examples
if (OpenMVG_BUILD_EXAMPLES)
  add_subdirectory(openMVG_Samples)
endif (OpenMVG_BUILD_EXAMPLES)

# Complete software(s) build on openMVG libraries
if (OpenMVG_BUILD_SOFTWARES)
  add_subdirectory(software)
endif (OpenMVG_BUILD_SOFTWARES)

# Included for research purpose only
add_subdirectory(nonFree)

# ==============================================================================
# Documentation
# --------------------------
# Sphinx detection
# ==============================================================================
if (OpenMVG_BUILD_DOC)
  find_package(Sphinx)
  if (EXISTS ${SPHINX_EXECUTABLE})
    set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/htmlDoc")

    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/../docs/sphinx/rst/conf.py"
      "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
      @ONLY)

    add_custom_target(doc ALL
      ${SPHINX_EXECUTABLE}
      -b html
      "${CMAKE_CURRENT_SOURCE_DIR}/../docs/sphinx/rst"
      "${SPHINX_HTML_DIR}"
      COMMENT "Building HTML documentation with Sphinx")

    set_property(TARGET doc PROPERTY FOLDER OpenMVG)

  else (EXISTS ${SPHINX_EXECUTABLE})
    message("Sphinx need to be installed to generate the sphinx documentation")
  endif (EXISTS ${SPHINX_EXECUTABLE})


  # ------
  # Doxygen detection
  # ======
  find_package(Doxygen)
  if (DOXYGEN_FOUND)

    configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/../docs/doxygen/Doxyfile.in
                    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY )

    add_custom_target( doc_doxygen ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMENT "Generating API documentation with Doxygen" VERBATIM )

  else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
  endif (DOXYGEN_FOUND)
endif()



message("\n")

message("** OpenMVG version: " ${OPENMVG_VERSION})
message("** Build Shared libs: " ${OpenMVG_BUILD_SHARED})
message("** Build OpenMVG tests: " ${OpenMVG_BUILD_TESTS})
message("** Build OpenMVG softwares: " ${OpenMVG_BUILD_SOFTWARES})
message("** Build OpenMVG documentation: " ${OpenMVG_BUILD_DOC})
message("** Build OpenMVG samples applications: " ${OpenMVG_BUILD_EXAMPLES})
message("** Build OpenMVG openGL examples: " ${OpenMVG_BUILD_OPENGL_EXAMPLES})
message("** Enable code coverage generation: " ${OpenMVG_BUILD_COVERAGE})
message("** Enable OpenMP parallelization: " ${OpenMVG_USE_OPENMP})
message("** Build OpenCV+OpenMVG samples programs: " ${OpenMVG_USE_OPENCV})
message("** Use OpenCV SIFT features: " ${OpenMVG_USE_OCVSIFT})

message("\n")

if (DEFINED OpenMVG_USE_INTERNAL_EIGEN)
  message(STATUS "EIGEN: " ${EIGEN_VERSION} " (internal)")
else()
  message(STATUS "EIGEN: " ${EIGEN_VERSION} " (external)")
endif()

if (DEFINED OpenMVG_USE_INTERNAL_CERES)
  message(STATUS "CERES: " ${CERES_VERSION} " (internal)")
else()
  message(STATUS "CERES: " ${CERES_VERSION} " (external)")
endif()

if (DEFINED OpenMVG_USE_INTERNAL_FLANN)
  message(STATUS "FLANN: " ${FLANN_VERSION} " (internal)")
else()
  message(STATUS "FLANN: " ${FLANN_VERSION} " (external)")
endif()

if (DEFINED OpenMVG_USE_INTERNAL_TIFF)
  message(STATUS "LIBTIFF: " ${TIFF_VERSION_STRING} " (internal)")
else()
  message(STATUS "LIBTIFF: " ${TIFF_VERSION_STRING} " (external)")
endif()

if (DEFINED OpenMVG_USE_INTERNAL_PNG)
  message(STATUS "LIBPNG: " ${PNG_VERSION_STRING} " (internal)")
else()
  message(STATUS "LIBPNG: " ${PNG_VERSION_STRING} " (external)")
endif()

if (DEFINED OpenMVG_USE_INTERNAL_JPEG)
  message(STATUS "LIBJPEG (internal)")
else()
  message(STATUS "LIBJPEG (external)")
endif()

if (DEFINED OpenMVG_USE_INTERNAL_CLP)
  message(STATUS "CLP: " ${CLP_VERSION} " (internal)")
else()
  message(STATUS "CLP: " ${CLP_VERSION} " (external)")
endif()

if (DEFINED OpenMVG_USE_INTERNAL_COINUTILS)
  message(STATUS "COINUTILS: " ${COINUTILS_VERSION} " (internal)")
else()
  message(STATUS "COINUTILS: " ${COINUTILS_VERSION} " (external)")
endif()

if (DEFINED OpenMVG_USE_INTERNAL_OSI)
  message(STATUS "OSI: " ${OSI_VERSION} " (internal)")
else()
  message(STATUS "OSI: " ${OSI_VERSION} " (external)")
endif()

if (DEFINED OpenMVG_USE_INTERNAL_LEMON)
  message(STATUS "LEMON: " ${LEMON_VERSION} " (internal)")
else()
  message(STATUS "LEMON: " ${LEMON_VERSION} " (external)")
endif()

message("\n")

# ==============================================================================
# INSTALL RULES
# ==============================================================================

install(EXPORT openMVG-targets DESTINATION lib)
install(EXPORT openMVG-targets
        DESTINATION share/openMVG/cmake FILE OpenMVGTargets.cmake)

#Adapt build include paths to install path
set(OpenMVG_INCLUDE_DIRS
  "${OpenMVG_INCLUDE_DIRS}"
  "${CMAKE_INSTALL_PREFIX}/include/openMVG")

string(REGEX REPLACE
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_INSTALL_PREFIX}/include"
  OpenMVG_INCLUDE_DIRS
  "${OpenMVG_INCLUDE_DIRS}"
)
string(REGEX REPLACE
  "third_party"
  "openMVG/third_party"
  OpenMVG_INCLUDE_DIRS
  "${OpenMVG_INCLUDE_DIRS}"
)
string(REGEX REPLACE
  "dependencies"
  "openMVG_dependencies"
  OpenMVG_INCLUDE_DIRS
  "${OpenMVG_INCLUDE_DIRS}"
)
if (OpenMVG_USE_INTERNAL_CERES)
  string(REGEX REPLACE
    "ceres-solver"
    "ceres"
    OpenMVG_INCLUDE_DIRS
    "${OpenMVG_INCLUDE_DIRS}"
  )
  string(REGEX REPLACE
    "internal/ceres/"
    ""
    OpenMVG_INCLUDE_DIRS
    "${OpenMVG_INCLUDE_DIRS}"
  )
endif (OpenMVG_USE_INTERNAL_CERES)

# Build OpenMVG_CFLAGS variable for export
set(OpenMVG_CFLAGS "")
foreach(d ${OpenMVG_DEFINITIONS})
  list(APPEND OpenMVG_CFLAGS "-D${d}")
endforeach()

# Create a OpenMVGConfig.cmake file. <name>Config.cmake files are searched by
# find_package() automatically. We configure that file so that we can put any
# information we want in it, e.g. version numbers, include directories, etc.
configure_file("${CMAKE_SOURCE_DIR}/cmakeFindModules/OpenMVGConfig.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/OpenMVGConfig.cmake" @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/OpenMVGConfig.cmake"
        DESTINATION share/openMVG/cmake)

# ==============================================================================
# UNINSTALL RULES
# ==============================================================================

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
