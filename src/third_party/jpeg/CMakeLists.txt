PROJECT(jpeg C)
cmake_minimum_required (VERSION 2.8)

INCLUDE (CheckIncludeFile)
check_include_file (stddef.h HAVE_STDDEF_H)
check_include_file (stdlib.h HAVE_STDLIB_H)
if (WIN32 AND NOT CYGWIN)
  SET (TWO_FILE_COMMANDLINE true)
endif ()
CONFIGURE_FILE (jconfig.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config/jconfig.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/config/)

# LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/lib)

SET (HEADERS jerror.h jmorecfg.h jpeglib.h ${CMAKE_CURRENT_BINARY_DIR}/config/jconfig.h)

SET (SRC jmemnobs.c jaricom.c jcapimin.c jcapistd.c jcarith.c jccoefct.c jccolor.c 
  jcdctmgr.c jchuff.c jcinit.c jcmainct.c jcmarker.c jcmaster.c jcomapi.c jcparam.c 
  jcprepct.c jcsample.c jctrans.c jdapimin.c jdapistd.c jdarith.c jdatadst.c jdatasrc.c 
  jdcoefct.c jdcolor.c jddctmgr.c jdhuff.c jdinput.c jdmainct.c jdmarker.c jdmaster.c 
  jdmerge.c jdpostct.c jdsample.c jdtrans.c jerror.c jfdctflt.c jfdctfst.c jfdctint.c 
  jidctflt.c jidctfst.c jidctint.c jquant1.c jquant2.c jutils.c jmemmgr.c jdct.h 
  jinclude.h jmemsys.h jpegint.h jversion.h)

IF(WIN32)
  ADD_DEFINITIONS( -DJPEGSTATIC -D_CRT_SECURE_NO_WARNINGS)
ENDIF(WIN32)

ADD_LIBRARY (jpeg STATIC ${SRC} ${HEADERS})
SET_PROPERTY(TARGET jpeg PROPERTY FOLDER OpenMVG/3rdParty)

if(JPEG_INCLUDE_INSTALL_DIR)
  set(INCLUDE_INSTALL_DIR ${JPEG_INCLUDE_INSTALL_DIR})
else()
  set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/jpeg")
endif()

install(FILES
  ${HEADERS}
  DESTINATION ${INCLUDE_INSTALL_DIR} COMPONENT Devel
)
INSTALL(TARGETS jpeg DESTINATION lib EXPORT openMVG-targets)
