
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
LIST(REMOVE_ITEM lInftycomputervision_cpp ${REMOVEFILESUNITTEST})

# Disable some warnings on MSVC.
IF (MSVC)
  # Disable Warning	C4267: conversion from 'size_t' to 'int', possible loss of data
  ADD_DEFINITIONS("/wd4267")
ENDIF (MSVC)

ADD_LIBRARY(openMVG_lInftyComputerVision  ${lInftycomputervision_headers} ${lInftycomputervision_cpp})
INSTALL(TARGETS openMVG_lInftyComputerVision DESTINATION lib/)

SET(LOCAL_LIBS
  openMVG_multiview_test_data
  openMVG_multiview
  openMVG_lInftyComputerVision
  ${openMVG_linearProgramming})

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

