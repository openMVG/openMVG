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

IF (NOT MOSEK_FOUND)

  FILE(GLOB_RECURSE REMOVEFILESMOSEK *MOSEK*)
  LIST(REMOVE_ITEM linearProgramming_cpp ${REMOVEFILESMOSEK})
  LIST(REMOVE_ITEM linearProgramming_headers ${REMOVEFILESMOSEK})
  
  SET(openMVG_linearProgramming
    lib_OsiClpSolver # solver wrapper
    lib_CoinUtils    # container tools
    lib_Osi          # generic LP
    lib_clp          # Solver
  )
  
ELSE (NOT MOSEK_FOUND)
  # MOSEK_FOUND is true
  ADD_LIBRARY(openMVG_linearProgrammingMSK ${linearProgramming_headers} ${linearProgramming_cpp})
  SET_PROPERTY(TARGET openMVG_linearProgrammingMSK PROPERTY FOLDER OpenMVG/OpenMVG)
  INSTALL(TARGETS openMVG_linearProgrammingMSK DESTINATION lib/)

  SET(openMVG_linearProgramming
    lib_Osi_Msk      # OSI solver wrapper for the Mosek backend
    lib_OsiClpSolver # solver wrapper
    lib_CoinUtils    # container tools
    lib_Osi          # generic LP
    lib_clp          # Solver

    openMVG_linearProgrammingMSK # openMVG Mosek wrapper
    ${MOSEK_LIB}    # Mosek lib
  )
ENDIF (NOT MOSEK_FOUND)

# set the variable in the parent_scope
SET(openMVG_linearProgramming ${openMVG_linearProgramming} PARENT_SCOPE)

UNIT_TEST(openMVG_linearProgramming linearProgramming "${openMVG_linearProgramming}")

ADD_SUBDIRECTORY(lInfinityCV)
