
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
SET_PROPERTY(TARGET vlsift PROPERTY FOLDER OpenMVG/patented)
