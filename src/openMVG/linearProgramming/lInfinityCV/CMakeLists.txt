
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
SET_SOURCE_FILES_PROPERTIES(${lInftycomputervision_headers} PROPERTIES LANGUAGE CXX)

FILE(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)

# Remove the future main files
LIST(REMOVE_ITEM lInftycomputervision_cpp ${REMOVEFILESUNITTEST})

ADD_LIBRARY(openMVG_lInftyComputerVision  ${lInftycomputervision_headers} ${lInftycomputervision_cpp})

# Linfinity Triangulation
UNIT_TEST(openMVG_lInftyCV triangulation "openMVG_multiview;openMVG_multiview_test_data;openMVG_linearProgramming;openMVG_lInftyComputerVision")

# Linfinity Pose/Resection estimation
UNIT_TEST(openMVG_lInftyCV resection 
 "openMVG_multiview_test_data;openMVG_multiview;openMVG_linearProgramming;openMVG_lInftyComputerVision")
UNIT_TEST(openMVG_lInftyCV resection_robust
 "openMVG_multiview;openMVG_multiview_test_data;openMVG_linearProgramming;openMVG_lInftyComputerVision")

# Linfinity Global translations from triplets of heading directions
UNIT_TEST(openMVG_lInftyCV global_translations_fromTriplets "openMVG_multiview_test_data;openMVG_multiview;openMVG_linearProgramming;openMVG_lInftyComputerVision")

# Linfinity Translation and Structure 
UNIT_TEST(openMVG_lInftyCV tijsAndXis "openMVG_multiview;openMVG_multiview_test_data;openMVG_linearProgramming;openMVG_lInftyComputerVision")
# Linfinity Translation and Structure with noise handling
UNIT_TEST(openMVG_lInftyCV tijsAndXis_outlier "openMVG_multiview;openMVG_multiview_test_data;openMVG_linearProgramming;openMVG_lInftyComputerVision")

