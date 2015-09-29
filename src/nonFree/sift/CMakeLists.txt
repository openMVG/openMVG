
# libs should be static
SET(BUILD_SHARED_LIBS OFF)

# use PIC code for link into shared lib
IF(UNIX)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
ENDIF(UNIX)

INCLUDE_DIRECTORIES(./vl)
SET(FEATS 
  vl/generic.c
  vl/imopv_sse2.c
  vl/sift.c
  vl/imopv.c
  vl/mathop_sse2.c
  vl/sift.c
  vl/host.c
  vl/mathop.c
  vl/random.c)
SET_SOURCE_FILES_PROPERTIES(${FEATS} PROPERTIES LANGUAGE C)
ADD_LIBRARY(vlsift ${FEATS})
INSTALL(TARGETS vlsift DESTINATION lib EXPORT openMVG-targets)
SET_PROPERTY(TARGET vlsift PROPERTY FOLDER OpenMVG/nonFree)
INSTALL(
  DIRECTORY .
  DESTINATION include/openMVG_dependencies/nonFree/sift
  COMPONENT headers
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)
