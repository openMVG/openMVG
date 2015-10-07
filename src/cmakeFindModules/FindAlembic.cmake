#
# This module tries to find and setup an Alembic configuration for openMVG.
# You can help the search by providing several environment variable or cmake 
# variable:
# ALEMBIC_ROOT
# HDF5_ROOT
# ILMBASE_ROOT
# OPENEXR_ROOT
#
# HDF5 and ILMBASE should point to the root dir used to compile alembic
#
# This module provides variables prefixed with ABC
# It also sets ALEMBIC_FOUND if all the libraries and include dirs were found
#

MESSAGE(STATUS "Looking for Alembic. 1.5.8")

################################################################################
# IlmBase include dir for half and ilm libraries used in alembic
################################################################################

# Alembic includes half.h for a single function "half to float", this is unfortunate
FIND_PATH(ABC_HALF_INCLUDE_DIR half.h
    HINTS
    ${ILMBASE_ROOT}/include/OpenEXR
    $ENV{ILMBASE_ROOT}/include/OpenEXR)

FIND_PATH(ABC_ILMBASE_LIBS_PATH NAMES libIex.so libIex.a 
    PATHS
        ${ILMBASE_ROOT}/lib 
        ${ILMBASE_ROOT}/lib64
        $ENV{ILMBASE_ROOT}/lib 
        $ENV{ILMBASE_ROOT}/lib64
    NO_DEFAULT_PATH)

FIND_LIBRARY(ABC_ILMBASE_IEX Iex PATHS ${ABC_ILMBASE_LIBS_PATH} NO_DEFAULT_PATH)
FIND_LIBRARY(ABC_ILMBASE_IEXMATH IexMath PATHS ${ABC_ILMBASE_LIBS_PATH} NO_DEFAULT_PATH)
SET(ABC_ILMBASE_LIBS ${ABC_ILMBASE_IEX} ${ABC_ILMBASE_IEXMATH}) 

# OpenEXR
FIND_LIBRARY(ABC_OPENEXR_LIBS IlmImf 
    PATHS 
        ${OPENEXR_ROOT}/lib
        ${OPENEXR_ROOT}/lib64
        $ENV{OPENEXR_ROOT}/lib 
        $ENV{OPENEXR_ROOT}/lib64
    NO_DEFAULT_PATH)

################################################################################
# HDF5 libraries used in alembic
################################################################################

# FIXME: hdf5 should be handled by a specialized module
FIND_PATH(ABC_HDF5_LIBS_PATH NAMES libhdf5.so libhdf5.a
       PATHS 
        ${HDF5_ROOT}/lib 
        ${HDF5_ROOT}/lib64
        $ENV{HDF5_ROOT}/lib 
        $ENV{HDF5_ROOT}/lib64
       NO_DEFAULT_PATH)
FIND_LIBRARY(ABC_HDF5 hdf5 PATHS ${ABC_HDF5_LIBS_PATH})
FIND_LIBRARY(ABC_HDF5_HL hdf5_hl PATHS ${ABC_HDF5_LIBS_PATH})
SET(ABC_HDF5_LIBS ${ABC_HDF5} ${ABC_HDF5_HL})

################################################################################
# ALEMBIC include and library dir
################################################################################

FIND_PATH(ABC_INCLUDE_DIR Alembic/Abc/All.h
    PATHS
        $ENV{ALEMBIC_ROOT}/include
        ${ALEMBIC_ROOT}/include
    PATH_SUFFIXES
        alembic
    NO_DEFAULT_PATH
)

#
# We force the use of dynamic libraries because we had initialization problems
# on linux when we used the static version of the libraries
#
FIND_PATH(ABC_LIBRARY_DIR NAMES libAlembicAbc.so libAlembicAbc.dylib
    PATHS
        ${ALEMBIC_ROOT}/lib
        ${ALEMBIC_ROOT}/lib64
        $ENV{ALEMBIC_ROOT}/lib
        $ENV{ALEMBIC_ROOT}/lib64
    NO_DEFAULT_PATH)
FIND_LIBRARY(ABC AlembicAbc PATHS ${ABC_LIBRARY_DIR})
FIND_LIBRARY(ABC_COREABSTRACT AlembicAbcCoreAbstract PATHS ${ABC_LIBRARY_DIR})
FIND_LIBRARY(ABC_COREHDF5 AlembicAbcCoreHDF5 PATHS ${ABC_LIBRARY_DIR})
FIND_LIBRARY(ABC_GEOM AlembicAbcGeom PATHS ${ABC_LIBRARY_DIR})
FIND_LIBRARY(ABC_UTIL AlembicUtil PATHS ${ABC_LIBRARY_DIR})
SET(ABC_CORE_LIBS ${ABC_GEOM} ${ABC} ${ABC_COREHDF5} ${ABC_COREABSTRACT} ${ABC_UTIL})

SET(ABC_LIBRARIES ${ABC_CORE_LIBS} ${ABC_HDF5_LIBS} "-ldl" ${ABC_OPENEXR_LIBS} ${ABC_ILMBASE_LIBS})
SET(ABC_INCLUDE_DIR ${ABC_INCLUDE_DIR} ${ABC_HALF_INCLUDE_DIR})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS("Alembic" DEFAULT_MSG ABC_LIBRARIES ABC_LIBRARY_DIR ABC_INCLUDE_DIR ABC_HDF5_LIBS)


mark_as_advanced(ABC_LIBRARY_DIR ABC_HDF5_LIBS_PATH ABC_ILMBASE_LIBS_PATH)
if (ALEMBIC_FOUND)
    message("Found Alembic - will build openmvg to abc exporter")
else()
    message("Alembic NOT FOUND")   
endif()

