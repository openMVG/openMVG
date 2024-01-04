PROJECT(jpeg C)
cmake_minimum_required (VERSION 2.8)

include(CheckIncludeFile)
check_include_file (stddef.h HAVE_STDDEF_H)
check_include_file (stdlib.h HAVE_STDLIB_H)
if (WIN32 AND NOT CYGWIN)
  set (TWO_FILE_COMMANDLINE true)
endif ()
configure_file (jconfig.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config/jconfig.h)

set (HEADERS jerror.h jmorecfg.h jpeglib.h ${CMAKE_CURRENT_BINARY_DIR}/config/jconfig.h)

set (SRC jmemnobs.c jaricom.c jcapimin.c jcapistd.c jcarith.c jccoefct.c jccolor.c
  jcdctmgr.c jchuff.c jcinit.c jcmainct.c jcmarker.c jcmaster.c jcomapi.c jcparam.c
  jcprepct.c jcsample.c jctrans.c jdapimin.c jdapistd.c jdarith.c jdatadst.c jdatasrc.c
  jdcoefct.c jdcolor.c jddctmgr.c jdhuff.c jdinput.c jdmainct.c jdmarker.c jdmaster.c
  jdmerge.c jdpostct.c jdsample.c jdtrans.c jerror.c jfdctflt.c jfdctfst.c jfdctint.c
  jidctflt.c jidctfst.c jidctint.c jquant1.c jquant2.c jutils.c jmemmgr.c jdct.h
  jinclude.h jmemsys.h jpegint.h jversion.h)

if(WIN32)
  add_definitions( -DJPEGSTATIC -D_CRT_SECURE_NO_WARNINGS)
endif(WIN32)

add_library (openMVG_jpeg STATIC ${SRC} ${HEADERS})
set_property(TARGET openMVG_jpeg PROPERTY FOLDER OpenMVG/3rdParty)

target_include_directories(openMVG_jpeg PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/config>
)

install(
  TARGETS openMVG_jpeg
  DESTINATION lib
  EXPORT openMVG-targets
)
