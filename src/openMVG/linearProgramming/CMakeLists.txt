FILE(
	GLOB
	linearProgramming_headers
	*.hpp
)
FILE(
	GLOB
	linearProgramming_cpp
	*.cpp
)
FILE(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)

SET_SOURCE_FILES_PROPERTIES(${linearProgramming_headers} PROPERTIES LANGUAGE CXX)

#Remove the future main files
LIST(REMOVE_ITEM linearProgramming_cpp ${REMOVEFILESUNITTEST})

ADD_LIBRARY(openMVG_linearProgramming ${linearProgramming_headers} ${linearProgramming_cpp})
TARGET_LINK_LIBRARIES(openMVG_linearProgramming
    lib_OsiClpSolver # solver wrapper
    lib_CoinUtils # container tools
    lib_Osi # generic LP
    lib_clp # Solver
)
IF (MOSEK_FOUND)
  TARGET_LINK_LIBRARIES(openMVG_linearProgramming
    pthread
    ${MOSEK_LIB}
    lib_Osi_Msk # OSI solver wrapper for the Mosek backend
  )
ENDIF (MOSEK_FOUND)

UNIT_TEST(openMVG_linearProgramming linearProgramming "openMVG_linearProgramming")

ADD_SUBDIRECTORY(lInfinitycomputervision)
