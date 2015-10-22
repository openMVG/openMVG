# Searches for an installation of the OPTPP library (http://software.sandia.gov/opt++/)
#
# Defines:
#
#   OPTPP_FOUND          True if OPTPP was found, else false
#   OPTPP_LIBRARIES      Libraries to link
#   OPTPP_INCLUDE_DIRS   The directories containing the header files
#   OPTPP_LDFLAGS        Extra linker flags
#
# To specify an additional directory to search, set OPTPP_ROOT.
#
# Author: Siddhartha Chaudhuri, 2011
#

SET(_optpp_SEARCH_DIRS
  ${OPTPP_ROOT}
  /usr/local
  /sw
  /opt/local
  /opt/local/opt++
)

set(OPTPP_FOUND FALSE)

if(NOT OPTPP_INCLUDE_DIRS)
  # Look for the OPTPP header, first in the user-specified location and then in the system locations
  set(OPTPP_INCLUDE_DOC "The directory containing the OPTPP include file Opt.h")

  find_path(OPTPP_INCLUDE_DIRS
            NAMES Opt.h
            PATHS ${OPTPP_ROOT}
            HINTS ${_optpp_SEARCH_DIRS}
            PATH_SUFFIXES "include" "include/OPTPP" "include/OPTPP" "include/opt++" "include/optpp"
            DOC ${OPTPP_INCLUDE_DOC}
            NO_DEFAULT_PATH)

  IF(OPTPP_INCLUDE_DIRS)
    MESSAGE( STATUS "opt++ header files found at ${OPTPP_INCLUDE_DIRS}" )
  ELSE(OPTPP_INCLUDE_DIRS)
    # now look in system locations
    find_path( OPTPP_INCLUDE_DIRS
               NAMES Opt.h
               PATH_SUFFIXES "OPTPP" "OPTPP" "opt++" "optpp"
	       DOC ${OPTPP_INCLUDE_DOC})
    IF(OPTPP_INCLUDE_DIRS)
      MESSAGE( STATUS "opt++ header files found at ${OPTPP_INCLUDE_DIRS}" )
    ELSE(OPTPP_INCLUDE_DIRS)
      MESSAGE( FATAL_ERROR "opt++ header files not found" )
    ENDIF(OPTPP_INCLUDE_DIRS)
  endif()
endif()

# Only look for the library file in the immediate neighbourhood of the include directory
set(OPTPP_LIBRARY_DIRS ${OPTPP_INCLUDE_DIRS})

if("${OPTPP_LIBRARY_DIRS}" MATCHES "/(OPT|opt)([+][+]|pp)$")
    # Strip off the trailing "/OPTPP" from the path
    GET_FILENAME_COMPONENT(OPTPP_LIBRARY_DIRS ${OPTPP_LIBRARY_DIRS} PATH)
endif()

if("${OPTPP_LIBRARY_DIRS}" MATCHES "/include$")
    # Strip off the trailing "/include" from the path
    GET_FILENAME_COMPONENT(OPTPP_LIBRARY_DIRS ${OPTPP_LIBRARY_DIRS} PATH)
endif()

# Look for libopt
find_library(OPTPP_DEBUG_LIBRARY NAMES opt_d libopt_d optd liboptd PATH_SUFFIXES "" Debug
             PATHS ${OPTPP_LIBRARY_DIRS} ${OPTPP_LIBRARY_DIRS}/lib NO_DEFAULT_PATH)

find_library(OPTPP_RELEASE_LIBRARY NAMES opt libopt PATH_SUFFIXES "" Release
             PATHS ${OPTPP_LIBRARY_DIRS} ${OPTPP_LIBRARY_DIRS}/lib NO_DEFAULT_PATH)

if(OPTPP_DEBUG_LIBRARY AND OPTPP_RELEASE_LIBRARY)
  message(STATUS "opt++ release library found at ${OPTPP_RELEASE_LIBRARY}" )
  message(STATUS "opt++ debug library found at   ${OPTPP_DEBUG_LIBRARY}" )
  set(OPTPP_LIBRARIES debug ${OPTPP_DEBUG_LIBRARY} optimized ${OPTPP_RELEASE_LIBRARY})
elseif(OPTPP_DEBUG_LIBRARY)
  message(STATUS "opt++ debug library found at   ${OPTPP_DEBUG_LIBRARY}" )
  set(OPTPP_LIBRARIES ${OPTPP_DEBUG_LIBRARY})
elseif(OPTPP_RELEASE_LIBRARY)
  message(STATUS "opt++ release library found at ${OPTPP_RELEASE_LIBRARY}" )
  set(OPTPP_LIBRARIES ${OPTPP_RELEASE_LIBRARY})
else()
  message(FATAL_ERROR "did not found any opt++ libraries")
endif()

# Look for libnewmat
find_path(OPTPP_NEWMAT_INCLUDE_DIRS
          NAMES include.h
          PATHS ${OPTPP_ROOT}
          HINTS ${_optpp_SEARCH_DIRS}
          PATH_SUFFIXES "include" "include/newmat11" "include/newmat" "newmat11" "newmat"
          NO_DEFAULT_PATH)

IF(OPTPP_NEWMAT_INCLUDE_DIRS)
  MESSAGE(STATUS "newmat include dirs found at ${OPTPP_NEWMAT_INCLUDE_DIRS}")
ELSE(OPTPP_NEWMAT_INCLUDE_DIRS)
  FIND_PATH( OPTPP_NEWMAT_INCLUDE_DIRS
             NAMES include.h
             PATH_SUFFIXES "include" "include/newmat11" "include/newmat" "newmat11" "newmat"
	     DOC ${OPTPP_INCLUDE_DOC})
  IF(NOT OPTPP_NEWMAT_INCLUDE_DIRS)
    MESSAGE( WARNING "newmat include dirs not found")
  ENDIF(NOT OPTPP_NEWMAT_INCLUDE_DIRS)
ENDIF(OPTPP_NEWMAT_INCLUDE_DIRS)

find_library(OPTPP_NEWMAT_DEBUG_LIBRARY
             NAMES newmat_d libnewmat_d newmatd libnewmatd PATH_SUFFIXES "" Debug
             PATHS ${OPTPP_LIBRARY_DIRS} ${OPTPP_LIBRARY_DIRS}/lib
	     NO_DEFAULT_PATH)
  
find_library(OPTPP_NEWMAT_RELEASE_LIBRARY
             NAMES newmat libnewmat
	     PATH_SUFFIXES "" Release
             PATHS ${OPTPP_LIBRARY_DIRS} ${OPTPP_LIBRARY_DIRS}/lib
	     NO_DEFAULT_PATH)

IF(NOT OPTPP_INCLUDE_DIRS STREQUAL OPTPP_NEWMAT_INCLUDE_DIRS)
  SET(OPTPP_INCLUDE_DIRS ${OPTPP_INCLUDE_DIRS} ${OPTPP_NEWMAT_INCLUDE_DIRS})
ENDIF(NOT OPTPP_INCLUDE_DIRS STREQUAL OPTPP_NEWMAT_INCLUDE_DIRS)

if(OPTPP_NEWMAT_DEBUG_LIBRARY AND OPTPP_NEWMAT_RELEASE_LIBRARY)
  message(STATUS "newmat release library found at ${OPTPP_NEWMAT_RELEASE_LIBRARY}" )
  message(STATUS "newmat debug library found at   ${OPTPP_NEWMAT_DEBUG_LIBRARY}" )
  set(OPTPP_LIBRARIES ${OPTPP_LIBRARIES} debug ${OPTPP_NEWMAT_DEBUG_LIBRARY} optimized ${OPTPP_NEWMAT_RELEASE_LIBRARY})
elseif(OPTPP_NEWMAT_DEBUG_LIBRARY)
  message(STATUS "newmat debug library found at   ${OPTPP_NEWMAT_DEBUG_LIBRARY}" )
  set(OPTPP_LIBRARIES ${OPTPP_LIBRARIES} ${OPTPP_NEWMAT_DEBUG_LIBRARY})
elseif(OPTPP_NEWMAT_RELEASE_LIBRARY)
  message(STATUS "newmat release library found at ${OPTPP_NEWMAT_RELEASE_LIBRARY}" )
  set(OPTPP_LIBRARIES ${OPTPP_LIBRARIES} ${OPTPP_NEWMAT_RELEASE_LIBRARY})
else()
  message(FATAL_ERROR "did not found any newmat libraries")
endif()

# Look for BLAS
FIND_PACKAGE(BLAS)
if(BLAS_FOUND)
  set(OPTPP_LIBRARIES ${OPTPP_LIBRARIES} ${BLAS_LIBRARIES})
  set(OPTPP_LDFLAGS ${BLAS_LINKER_FLAGS})
else()
  message(FATAL_ERROR "OPTPP: BLAS library not found")
  set(OPTPP_LIBRARIES )
endif()

# now we have opt++, newmat and BLAS, opt++ is probably going to work
set(OPTPP_FOUND TRUE)

if(OPTPP_FOUND)
  if(NOT OPTPP_FIND_QUIETLY)
    message(STATUS "Found OPTPP: headers at ${OPTPP_INCLUDE_DIRS}, libraries at ${OPTPP_LIBRARIES}")
  endif()
else()
  if(OPTPP_FIND_REQUIRED)
    message(FATAL_ERROR "OPTPP not found")
  endif()
endif()

