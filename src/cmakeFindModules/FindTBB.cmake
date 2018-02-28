# - Find ThreadingBuildingBlocks include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(TBB
#    [REQUIRED]             # Fail with error if TBB is not found
#    )                      #
# Once done, this will define
#
#  TBB_FOUND - system has TBB
#  TBB_INCLUDE_DIRS - the TBB include directories
#  TBB_LIBRARIES - TBB libraries to be lined, doesn't include malloc or
#                  malloc proxy
#
#  TBB_VERSION_MAJOR - Major Product Version Number
#  TBB_VERSION_MINOR - Minor Product Version Number
#  TBB_INTERFACE_VERSION - Engineering Focused Version Number
#  TBB_COMPATIBLE_INTERFACE_VERSION - The oldest major interface version
#                                     still supported. This uses the engineering
#                                     focused interface version numbers.
#
#  TBB_MALLOC_FOUND - system has TBB malloc library
#  TBB_MALLOC_INCLUDE_DIRS - the TBB malloc include directories
#  TBB_MALLOC_LIBRARIES - The TBB malloc libraries to be lined
#
#  TBB_MALLOC_PROXY_FOUND - system has TBB malloc proxy library
#  TBB_MALLOC_PROXY_INCLUDE_DIRS = the TBB malloc proxy include directories
#  TBB_MALLOC_PROXY_LIBRARIES - The TBB malloc proxy libraries to be lined
#
#
# This module reads hints about search locations from variables:
#  ENV TBB_ARCH_PLATFORM - for eg. set it to "mic" for Xeon Phi builds
#  ENV TBB_ROOT or just TBB_ROOT - root directory of tbb installation
#  ENV TBB_BUILD_PREFIX - specifies the build prefix for user built tbb
#                         libraries. Should be specified with ENV TBB_ROOT
#                         and optionally...
#  ENV TBB_BUILD_DIR - if build directory is different than ${TBB_ROOT}/build
#
#
# Modified by Robert Maynard from the original OGRE source
#
#-------------------------------------------------------------------
# This file is part of the CMake build system for OGRE
#     (Object-oriented Graphics Rendering Engine)
# For the latest info, see http://www.ogre3d.org/
#
# The contents of this file are placed in the public domain. Feel
# free to make use of it in any way you like.
#-------------------------------------------------------------------
#
#=============================================================================
# Copyright 2010-2012 Kitware, Inc.
# Copyright 2012      Rolf Eike Beer <eike@sf-mail.de>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)


#=============================================================================
#  FindTBB helper functions and macros
#

#===============================================
# Do the final processing for the package find.
#===============================================
macro(findpkg_finish PREFIX)
    # skip if already processed during this run
    if (NOT ${PREFIX}_FOUND)
        if (${PREFIX}_INCLUDE_DIR AND ${PREFIX}_LIBRARY)
            set(${PREFIX}_FOUND TRUE)
            set (${PREFIX}_INCLUDE_DIRS ${${PREFIX}_INCLUDE_DIR})
            set (${PREFIX}_LIBRARIES ${${PREFIX}_LIBRARY})
        else ()
            if (${PREFIX}_FIND_REQUIRED AND NOT ${PREFIX}_FIND_QUIETLY)
                message(FATAL_ERROR "Required library ${PREFIX} not found.")
            endif ()
        endif ()

        #mark the following variables as internal variables
        mark_as_advanced(${PREFIX}_INCLUDE_DIR
                ${PREFIX}_LIBRARY
                ${PREFIX}_LIBRARY_DEBUG
                ${PREFIX}_LIBRARY_RELEASE)
    endif ()
endmacro(findpkg_finish)

#===============================================
# Generate debug names from given RELEASEease names
#===============================================
macro(get_debug_names PREFIX)
    foreach(i ${${PREFIX}})
        set(${PREFIX}_DEBUG ${${PREFIX}_DEBUG} ${i}d ${i}D ${i}_d ${i}_D ${i}_debug ${i})
    endforeach(i)
endmacro(get_debug_names)

#===============================================
# See if we have env vars to help us find tbb
#===============================================
macro(getenv_path VAR)
    set(ENV_${VAR} $ENV{${VAR}})
    # replace won't work if var is blank
    if (ENV_${VAR})
        string( REGEX REPLACE "\\\\" "/" ENV_${VAR} ${ENV_${VAR}} )
    endif ()
endmacro(getenv_path)

#===============================================
# Couple a set of RELEASEease AND debug libraries
#===============================================
macro(make_library_set PREFIX)
    if (${PREFIX}_RELEASE AND ${PREFIX}_DEBUG)
        set(${PREFIX} optimized ${${PREFIX}_RELEASE} debug ${${PREFIX}_DEBUG})
    elseif (${PREFIX}_RELEASE)
        set(${PREFIX} ${${PREFIX}_RELEASE})
    elseif (${PREFIX}_DEBUG)
        set(${PREFIX} ${${PREFIX}_DEBUG})
    endif ()
endmacro(make_library_set)


#=============================================================================
#  Now to actually find TBB
#

# Get path, convert backslashes as ${ENV_${var}}
getenv_path(TBB_ROOT)

# initialize search paths
set(TBB_PREFIX_PATH ${TBB_ROOT} ${ENV_TBB_ROOT})
set(TBB_INC_SEARCH_PATH "")
set(TBB_LIB_SEARCH_PATH "")


# If user built from sources
set(TBB_BUILD_PREFIX $ENV{TBB_BUILD_PREFIX})
if (TBB_BUILD_PREFIX AND ENV_TBB_ROOT)
    getenv_path(TBB_BUILD_DIR)
    if (NOT ENV_TBB_BUILD_DIR)
        set(ENV_TBB_BUILD_DIR ${ENV_TBB_ROOT}/build)
    endif ()

    # include directory under ${ENV_TBB_ROOT}/include
    list(APPEND TBB_LIB_SEARCH_PATH
            ${ENV_TBB_BUILD_DIR}/${TBB_BUILD_PREFIX}_release
            ${ENV_TBB_BUILD_DIR}/${TBB_BUILD_PREFIX}_debug)
endif ()


# For Windows, let's assume that the user might be using the precompiled
# TBB packages from the main website. These use a rather awkward directory
# structure (at least for automatically finding the right files) depending
# on platform and compiler, but we'll do our best to accommodate it.
# Not adding the same effort for the precompiled linux builds, though. Those
# have different versions for CC compiler versions and linux kernels which
# will never adequately match the user's setup, so there is no feasible way
# to detect the "best" version to use. The user will have to manually
# select the right files. (Chances are the distributions are shipping their
# custom version of tbb, anyway, so the problem is probably nonexistant.)
if (WIN32 AND MSVC)
    set(COMPILER_PREFIX "vc7.1")
    if (MSVC_VERSION EQUAL 1400)
        set(COMPILER_PREFIX "vc8")
    elseif(MSVC_VERSION EQUAL 1500)
        set(COMPILER_PREFIX "vc9")
    elseif(MSVC_VERSION EQUAL 1600)
        set(COMPILER_PREFIX "vc10")
    elseif(MSVC_VERSION EQUAL 1700)
        set(COMPILER_PREFIX "vc11")
    elseif(MSVC_VERSION EQUAL 1800)
        set(COMPILER_PREFIX "vc12")
    elseif(MSVC_VERSION EQUAL 1900)
        set(COMPILER_PREFIX "vc14")
    endif ()

    # for each prefix path, add ia32/64\${COMPILER_PREFIX}\lib to the lib search path
    foreach (dir ${TBB_PREFIX_PATH})
        if (CMAKE_CL_64)
            list(APPEND TBB_LIB_SEARCH_PATH ${dir}/ia64/${COMPILER_PREFIX}/lib)
            list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/ia64/${COMPILER_PREFIX})
            list(APPEND TBB_LIB_SEARCH_PATH ${dir}/intel64/${COMPILER_PREFIX}/lib)
            list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/intel64/${COMPILER_PREFIX})
        else ()
            list(APPEND TBB_LIB_SEARCH_PATH ${dir}/ia32/${COMPILER_PREFIX}/lib)
            list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/ia32/${COMPILER_PREFIX})
        endif ()
    endforeach ()
endif ()

# For OS X binary distribution, choose libc++ based libraries for Maverics and
# above and AppleClang
if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin" AND
        NOT ${CMAKE_SYSTEM_VERSION} LESS 13.0)
    set (USE_LIBCXX OFF)
    cmake_policy(GET CMP0025 POLICY_VAR)

    if ("${POLICY_VAR}" STREQUAL "NEW")
        if (${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
            set (USE_LIBCXX ON)
        endif ()
    else ()
        if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
            set (USE_LIBCXX ON)
        endif ()
    endif ()

    if (${USE_LIBCXX})
        foreach (dir ${TBB_PREFIX_PATH})
            list (APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/libc++ ${dir}/libc++/lib)
        endforeach ()
    endif ()
endif ()

# check compiler ABI
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
        CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.4)
    set(COMPILER_PREFIX "gcc4.1")
else () # Assume compatibility with 4.4 for other compilers
    set(COMPILER_PREFIX "gcc4.4")
endif ()

# if platform architecture is explicitly specified
set(TBB_ARCH_PLATFORM $ENV{TBB_ARCH_PLATFORM})
if (TBB_ARCH_PLATFORM)
    foreach (dir ${TBB_PREFIX_PATH})
        list(APPEND TBB_LIB_SEARCH_PATH ${dir}/${TBB_ARCH_PLATFORM}/lib)
        list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/${TBB_ARCH_PLATFORM})
    endforeach ()
endif ()

foreach (dir ${TBB_PREFIX_PATH})
    if (CMAKE_SIZEOF_VOID_P EQUAL 8)
        list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/intel64)
        list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/intel64/${COMPILER_PREFIX})
        list(APPEND TBB_LIB_SEARCH_PATH ${dir}/intel64/lib)
        list(APPEND TBB_LIB_SEARCH_PATH ${dir}/intel64/${COMPILER_PREFIX}/lib)
    else ()
        list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/ia32)
        list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib/ia32/${COMPILER_PREFIX})
        list(APPEND TBB_LIB_SEARCH_PATH ${dir}/ia32/lib)
        list(APPEND TBB_LIB_SEARCH_PATH ${dir}/ia32/${COMPILER_PREFIX}/lib)
    endif ()
endforeach ()

# add general search paths
foreach (dir ${TBB_PREFIX_PATH})
    list(APPEND TBB_LIB_SEARCH_PATH ${dir}/lib ${dir}/Lib ${dir}/lib/tbb
            ${dir}/Libs)
    list(APPEND TBB_INC_SEARCH_PATH ${dir}/include ${dir}/Include
            ${dir}/include/tbb)
endforeach ()

set(TBB_LIBRARY_NAMES tbb)
get_debug_names(TBB_LIBRARY_NAMES)


find_path(TBB_INCLUDE_DIR
        NAMES tbb/tbb.h
        PATHS ${TBB_INC_SEARCH_PATH})

find_library(TBB_LIBRARY_RELEASE
        NAMES ${TBB_LIBRARY_NAMES}
        PATHS ${TBB_LIB_SEARCH_PATH})
find_library(TBB_LIBRARY_DEBUG
        NAMES ${TBB_LIBRARY_NAMES_DEBUG}
        PATHS ${TBB_LIB_SEARCH_PATH})
make_library_set(TBB_LIBRARY)

findpkg_finish(TBB)

#if we haven't found TBB no point on going any further
if (NOT TBB_FOUND)
    return()
endif ()

#=============================================================================
# Look for TBB's malloc package
set(TBB_MALLOC_LIBRARY_NAMES tbbmalloc)
get_debug_names(TBB_MALLOC_LIBRARY_NAMES)

find_path(TBB_MALLOC_INCLUDE_DIR
        NAMES tbb/tbb.h
        PATHS ${TBB_INC_SEARCH_PATH})

find_library(TBB_MALLOC_LIBRARY_RELEASE
        NAMES ${TBB_MALLOC_LIBRARY_NAMES}
        PATHS ${TBB_LIB_SEARCH_PATH})
find_library(TBB_MALLOC_LIBRARY_DEBUG
        NAMES ${TBB_MALLOC_LIBRARY_NAMES_DEBUG}
        PATHS ${TBB_LIB_SEARCH_PATH})
make_library_set(TBB_MALLOC_LIBRARY)

findpkg_finish(TBB_MALLOC)

#=============================================================================
# Look for TBB's malloc proxy package
set(TBB_MALLOC_PROXY_LIBRARY_NAMES tbbmalloc_proxy)
get_debug_names(TBB_MALLOC_PROXY_LIBRARY_NAMES)

find_path(TBB_MALLOC_PROXY_INCLUDE_DIR
        NAMES tbb/tbbmalloc_proxy.h
        PATHS ${TBB_INC_SEARCH_PATH})

find_library(TBB_MALLOC_PROXY_LIBRARY_RELEASE
        NAMES ${TBB_MALLOC_PROXY_LIBRARY_NAMES}
        PATHS ${TBB_LIB_SEARCH_PATH})
find_library(TBB_MALLOC_PROXY_LIBRARY_DEBUG
        NAMES ${TBB_MALLOC_PROXY_LIBRARY_NAMES_DEBUG}
        PATHS ${TBB_LIB_SEARCH_PATH})
make_library_set(TBB_MALLOC_PROXY_LIBRARY)

findpkg_finish(TBB_MALLOC_PROXY)


#=============================================================================
#parse all the version numbers from tbb
if(NOT TBB_VERSION)

    #only read the start of the file
    file(READ
            "${TBB_INCLUDE_DIR}/tbb/tbb_stddef.h"
            TBB_VERSION_CONTENTS
            LIMIT 2048)

    string(REGEX REPLACE
            ".*#define TBB_VERSION_MAJOR ([0-9]+).*" "\\1"
            TBB_VERSION_MAJOR "${TBB_VERSION_CONTENTS}")

    string(REGEX REPLACE
            ".*#define TBB_VERSION_MINOR ([0-9]+).*" "\\1"
            TBB_VERSION_MINOR "${TBB_VERSION_CONTENTS}")

    string(REGEX REPLACE
            ".*#define TBB_INTERFACE_VERSION ([0-9]+).*" "\\1"
            TBB_INTERFACE_VERSION "${TBB_VERSION_CONTENTS}")

    string(REGEX REPLACE
            ".*#define TBB_COMPATIBLE_INTERFACE_VERSION ([0-9]+).*" "\\1"
            TBB_COMPATIBLE_INTERFACE_VERSION "${TBB_VERSION_CONTENTS}")

endif()
