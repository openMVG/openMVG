# CMakeLists.txt

# Copyright (C) 2007-2013 Glenn Randers-Pehrson

# This code is released under the libpng license.
# For conditions of distribution and use, see the disclaimer
# and license in png.h

cmake_minimum_required(VERSION 2.4.4)

project(libpng C)
enable_testing()

set(PNGLIB_MAJOR 1)
set(PNGLIB_MINOR 6)
set(PNGLIB_RELEASE 6)
set(PNGLIB_NAME libpng${PNGLIB_MAJOR}${PNGLIB_MINOR})
set(PNGLIB_VERSION ${PNGLIB_MAJOR}.${PNGLIB_MINOR}.${PNGLIB_RELEASE})

if(NOT WIN32)
  find_library(M_LIBRARY NAMES m PATHS /usr/lib /usr/local/lib)
  if(NOT M_LIBRARY)
    message(STATUS
      "math library 'libm' not found - floating point support disabled")
  endif()
else()
  # not needed on windows
  set(M_LIBRARY "")
endif()


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pnglibconf.h.prebuilt
               ${CMAKE_CURRENT_BINARY_DIR}/config/pnglibconf.h)

set(libpng_public_hdrs
  png.h
  pngconf.h
  ${CMAKE_CURRENT_BINARY_DIR}/config/pnglibconf.h
)
set(libpng_sources
  ${libpng_public_hdrs}
  pngdebug.h
  pnginfo.h
  pngpriv.h
  pngstruct.h
  png.c
  pngerror.c
  pngget.c
  pngmem.c
  pngpread.c
  pngread.c
  pngrio.c
  pngrtran.c
  pngrutil.c
  pngset.c
  pngtrans.c
  pngwio.c
  pngwrite.c
  pngwtran.c
  pngwutil.c
)

if(MSVC)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
endif(MSVC)
add_library(openMVG_png STATIC ${libpng_sources})
target_link_libraries(openMVG_png PUBLIC openMVG_zlib ${M_LIBRARY})
set_property(TARGET openMVG_png PROPERTY FOLDER OpenMVG/3rdParty)

target_include_directories(openMVG_png PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/config/>
)

install(
  TARGETS openMVG_png
  DESTINATION lib
  EXPORT openMVG-targets
)
