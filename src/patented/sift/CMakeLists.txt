
INCLUDE_DIRECTORIES(./vl)
SET(FEATS 
  vl/generic.c
  vl/generic.h
  vl/host.h
  vl/imopv_sse2.c
  vl/mathop.h
  vl/sift.c
  vl/generic.h
  vl/imopv.c
  vl/imopv_sse2.h
  vl/mathop_sse2.c
  vl/sift.h
  vl/sift.c
  vl/host.c
  vl/imopv.h
  vl/mathop.c
  vl/mathop_sse2.h
  vl/random.h
  vl/random.c)
SET_SOURCE_FILES_PROPERTIES(${FEATS} PROPERTIES LANGUAGE C)
ADD_LIBRARY(vlsift ${FEATS})

