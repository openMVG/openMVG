
FILE(
	GLOB
	lInftycomputervision_headers
	*.hpp
)

FILE(
	GLOB
	lInftycomputervision_cpp
	*.cpp
)

FILE(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)

# Remove the future main files
list(REMOVE_ITEM lInftycomputervision_cpp ${REMOVEFILESUNITTEST})

# Disable some warnings on MSVC.
if (MSVC)
  # Disable Warning	C4267: conversion from 'size_t' to 'int', possible loss of data
  ADD_DEFINITIONS("/wd4267")
endif (MSVC)

add_library(openMVG_lInftyComputerVision  ${lInftycomputervision_headers} ${lInftycomputervision_cpp})
set(LINEAR_PROGRAMMING_LIBRAIRIES
  ${CLP_LIBRARIES}     # clp + solver wrapper
  ${COINUTILS_LIBRARY} # container tools
  ${OSI_LIBRARY}       # generic LP
)
if (MOSEK_FOUND)
  list(APPEND LINEAR_PROGRAMMING_LIBRAIRIES
    openMVG_linearProgrammingMSK # openMVG Mosek wrapper
    ${MOSEK_LIB}    # Mosek lib
  )
endif (MOSEK_FOUND)
target_link_libraries(openMVG_lInftyComputerVision ${LINEAR_PROGRAMMING_LIBRAIRIES})
set_target_properties(openMVG_lInftyComputerVision PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
install(TARGETS openMVG_lInftyComputerVision DESTINATION lib EXPORT openMVG-targets)

set(LOCAL_LIBS
  openMVG_lInftyComputerVision
  openMVG_multiview_test_data)

# Linfinity Triangulation
UNIT_TEST(openMVG_lInftyCV triangulation "${LOCAL_LIBS}")

# Linfinity Pose/Resection estimation
UNIT_TEST(openMVG_lInftyCV resection "${LOCAL_LIBS}")
UNIT_TEST(openMVG_lInftyCV resection_robust "${LOCAL_LIBS}")

# Linfinity Global translations from triplets of heading directions
UNIT_TEST(openMVG_lInftyCV global_translations_fromTriplets "${LOCAL_LIBS}")
# Linfinity Global translations from heading directions
UNIT_TEST(openMVG_lInftyCV global_translations_fromTij "${LOCAL_LIBS}")

# Linfinity Translation and Structure
UNIT_TEST(openMVG_lInftyCV tijsAndXis "${LOCAL_LIBS}")
# Linfinity Translation and Structure with noise handling
UNIT_TEST(openMVG_lInftyCV tijsAndXis_outlier "${LOCAL_LIBS}")
