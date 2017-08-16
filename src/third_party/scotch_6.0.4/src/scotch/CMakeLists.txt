## Copyright 2014,2015 IPB, Universite de Bordeaux, INRIA & CNRS
##
## This file is part of the Scotch software package for static mapping,
## graph partitioning and sparse matrix ordering.
##
## This software is governed by the CeCILL-C license under French law
## and abiding by the rules of distribution of free software. You can
## use, modify and/or redistribute the software under the terms of the
## CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
## URL: "http://www.cecill.info".
## 
## As a counterpart to the access to the source code and rights to copy,
## modify and redistribute granted by the license, users are provided
## only with a limited warranty and the software's author, the holder of
## the economic rights, and the successive licensors have only limited
## liability.
## 
## In this respect, the user's attention is drawn to the risks associated
## with loading, using, modifying and/or developing or reproducing the
## software by the user in light of its specific status of free software,
## that may mean that it is complicated to manipulate, and that also
## therefore means that it is reserved for developers and experienced
## professionals having in-depth computer knowledge. Users are therefore
## encouraged to load and test the software's suitability as regards
## their requirements in conditions enabling the security of their
## systems and/or data to be ensured and, more generally, to use and
## operate it in the same conditions as regards security.
## 
## The fact that you are presently reading this means that you have had
## knowledge of the CeCILL-C license and that you accept its terms.
##
##############################################################
###                                                        ###
###   AUTHOR     : Cedric LACHAT                           ###
###                                                        ###
###   FUNCTION   : Secondary configuration file for CMake  ###
###                                                        ###
###   DATES      : # Version 6.0  : from : 01 sep 2014     ###
###                                 to     01 mar 2015     ###
###                                                        ###
##############################################################

include_directories(include ${SCOTCH_BINARY_DIR}/src/libscotch)
include_directories(${SCOTCH_SOURCE_DIR}/src/libscotch)

###########
#  Execs  #
###########

# Search source files with components
set (SCOTCH_EXECS
acpl
amk_ccc
amk_fft2
amk_grf
amk_hy
amk_m2
amk_p2
atst
gbase
gcv
gmap
gmk_hy
gmk_m2
gmk_m3
gmk_msh
gmk_ub2
gmtst
gord
gotst
gout
gpart
gscat
gtst
mcv
mmk_m2
mmk_m3
mord
mtst
)

if (PTSCOTCH)
  set (PTSCOTCH_EXECS
    dggath
    dgmap
    dgord
    dgpart
    dgscat
    dgtst
    )
endif (PTSCOTCH)

####################
#  Parallel execs  #
####################

# Combine parallel source file with execution file
foreach(EXEC_NAME ${PTSCOTCH_EXECS})
  if (NOT ${EXEC_NAME} STREQUAL "dgpart")
    add_executable(${EXEC_NAME} ${EXEC_NAME}.c)
  else ()
    add_executable(${EXEC_NAME} dgmap.c)
    set_target_properties(${EXEC_NAME} PROPERTIES COMPILE_FLAGS "-DSCOTCH_COMPILE_PART")
  endif ()
  target_link_libraries(${EXEC_NAME} scotcherr ptscotch scotch ${SCOTCH_DEPENDENCIES})

  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME} DESTINATION bin
    PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE)
  list(APPEND C_EXECS ${EXEC_NAME})

  ## Add test
  #if (TEST_FILES MATCHES ".*${EXEC_NAME}.*")
  #  if (OLD_CMAKE)
  #    add_test( ${EXEC_NAME}
  #      ${CMAKE_MODULE_PATH}/run-partests.sh ${MPIEXEC} ${NB_PROC}
  #      ${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME}-out ${TEST_DIR}/${EXEC_NAME})
  #  else ()
  #    add_test(NAME ${EXEC_NAME}
  #      COMMAND ${CMAKE_MODULE_PATH}/run-partests.sh ${MPIEXEC} ${NB_PROC} ${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME}
  #      ${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME}-out ${TEST_DIR}/${EXEC_NAME})
  #    set_tests_properties (${EXEC_NAME} PROPERTIES
  #      TIMEOUT 60
  #      )
  #  endif ()
  #endif ()
  ## Add memory test
  #    if (OLD_CMAKE)
  #	add_test( MEM_${EXEC_NAME}
  #     ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NB_PROC} ${VALGRIND_EXE} ${VALGRIND_ARGS}
  #    ${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME})
  #    else ()
  #	add_test(NAME MEM_${EXEC_NAME}
  #    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NB_PROC} ${VALGRIND_EXE} ${VALGRIND_ARGS}
  #    ${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME})
  #    endif ()
  #set_tests_properties (MEM_${EXEC_NAME} PROPERTIES
  #  FAIL_REGULAR_EXPRESSION "problem with execution of;mpiexec was unable to launch;definitely lost: [^0][0-9,]* bytes in [^0][0-9,]* blocks;indirectly lost: [^0][0-9,]* bytes in [^0][0-9,]* blocks;possibly lost: [^0][0-9,]* bytes in [^0][0-9,]* blocks;ERROR SUMMARY: [^0][0-9,]* errors from [^0][0-9,]* contexts;Conditional jump or move depends on uninitialised value;Invalid write;Invalid read"
  #  TIMEOUT 60
  #  PASS_REGULAR_EXPRESSION "fin du prog")
endforeach()

######################
#  Sequential execs  #
######################

# Combine sequential source file with execution file

message (STATUS "depend: ${SCOTCH_DEPENDENCIES}")
foreach (EXEC_NAME ${SCOTCH_EXECS})
  if (NOT ${EXEC_NAME} STREQUAL "gout" AND NOT ${EXEC_NAME} STREQUAL "gpart")
    add_executable(${EXEC_NAME} ${EXEC_NAME}.c)
  else ()
    if (${EXEC_NAME} STREQUAL "gout")
      add_executable (${EXEC_NAME} ${EXEC_NAME}_c.c ${EXEC_NAME}_o.c)
    else ()
      add_executable (${EXEC_NAME} gmap.c)
      set_target_properties (${EXEC_NAME} PROPERTIES COMPILE_FLAGS "-DSCOTCH_COMPILE_PART")
    endif ()
  endif ()
  target_link_libraries (${EXEC_NAME} scotcherr scotch ${SCOTCH_DEPENDENCIES})

  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME} DESTINATION bin
    PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE)
  list (APPEND C_EXECS ${EXEC_NAME})

  ## Add test
  #if (TEST_FILES MATCHES ".*${EXEC_NAME}.*")
  #  	if (OLD_CMAKE)
  #    add_test (${EXEC_NAME}
  #    	 ${CMAKE_MODULE_PATH}/run-seqtests.sh
  #    	${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME}-out ${TEST_DIR}/${EXEC_NAME})
  #    	else ()
  #    add_test (NAME ${EXEC_NAME}
  #    	COMMAND ${CMAKE_MODULE_PATH}/run-seqtests.sh
  #    	${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME}-out ${TEST_DIR}/${EXEC_NAME})
  #    	endif ()
  #endif()
  ## Add memory test
  #    if (OLD_CMAKE)
  #	add_test (MEM_${EXEC_NAME}
  #     ${VALGRIND_EXE} ${VALGRIND_ARGS}
  #    ${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME})
  #    else ()
  #	add_test (NAME MEM_${EXEC_NAME}
  #    COMMAND ${VALGRIND_EXE} ${VALGRIND_ARGS}
  #    ${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME})
  #    endif ()
  #set_tests_properties (MEM_${EXEC_NAME} PROPERTIES
  #  FAIL_REGULAR_EXPRESSION "definitely lost: [^0][0-9,]* bytes in [^0][0-9,]* blocks;indirectly lost: [^0][0-9,]* bytes in [^0][0-9,]* blocks;possibly lost: [^0][0-9,]* bytes in [^0][0-9,]* blocks;ERROR SUMMARY: [^0][0-9,]* errors from [^0][0-9,]* contexts;Conditional jump or move depends on uninitialised value;Invalid write;Invalid read")
endforeach()
