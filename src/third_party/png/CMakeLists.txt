PROJECT(png C)

SET(PNG_SRCS
  png.h  pnglibconf.h pngconf.h  pngdebug.h  pnginfo.h  pngpriv.h  pngstruct.h
  png.c pngerror.c  pngget.c
  pngmem.c  pngpread.c  pngread.c  pngrio.c  pngrtran.c  pngrutil.c  pngset.c
  pngtrans.c  pngwio.c  pngwrite.c  pngwtran.c  pngwutil.c	   
  )

IF (WIN32)
  ADD_DEFINITIONS(-DPNG_STATIC -D_CRT_SECURE_NO_WARNINGS)
#  # TODO(keir): I don't know why, but CMake *refuses* to add this to the
#  # command line flags. After much frustration, I am giving up and leaving the
#  # compile errors. Please, someone who knows CMake, try to fix this.
#  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4996")
ENDIF (WIN32)

INCLUDE_DIRECTORIES(../zlib)
ADD_LIBRARY(png ${PNG_SRCS})
TARGET_LINK_LIBRARIES(png zlib)
