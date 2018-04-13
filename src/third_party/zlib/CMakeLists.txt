PROJECT(zlib)

ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)

ADD_LIBRARY(openMVG_zlib STATIC
  adler32.c    compress.c   crc32.c      deflate.c    gzio.c       infback.c
  inffast.c    inflate.c    inftrees.c   minigzip.c   trees.c      uncompr.c
  zutil.c      crc32.h      deflate.h    inffast.h    inffixed.h   inflate.h
  inftrees.h   trees.h      zconf.h      zconf.in.h   zlib.h       zutil.h
)
SET_PROPERTY(TARGET openMVG_zlib PROPERTY FOLDER OpenMVG/3rdParty)

target_include_directories(openMVG_zlib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

install(
  TARGETS openMVG_zlib
  DESTINATION lib
  EXPORT openMVG-targets
)
