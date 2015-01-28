PROJECT(zlib)

ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)

ADD_LIBRARY(zlib STATIC
  adler32.c    compress.c   crc32.c      deflate.c    gzio.c       infback.c
  inffast.c    inflate.c    inftrees.c   minigzip.c   trees.c      uncompr.c
  zutil.c      crc32.h      deflate.h    inffast.h    inffixed.h   inflate.h
  inftrees.h   trees.h      zconf.h      zconf.in.h   zlib.h       zutil.h
)
SET_PROPERTY(TARGET zlib PROPERTY FOLDER OpenMVG/3rdParty)

if(ZLIB_INCLUDE_INSTALL_DIR)
  set(INCLUDE_INSTALL_DIR ${ZLIB_INCLUDE_INSTALL_DIR})
else()
  set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/zlib")
endif()

install(FILES
  crc32.h      deflate.h    inffast.h    inffixed.h   inflate.h
  inftrees.h   trees.h      zconf.h      zconf.in.h   zlib.h       zutil.h
  DESTINATION ${INCLUDE_INSTALL_DIR} COMPONENT Devel
)
INSTALL(TARGETS zlib DESTINATION lib/)
