
FILE(
	GLOB_RECURSE
	sfm_files_header
	*.hpp
)
FILE(
	GLOB_RECURSE
	sfm_files_cpp
	*.cpp
)
FILE(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM sfm_files_cpp ${REMOVEFILESUNITTEST})
#Remove the test header main files
FILE(GLOB_RECURSE REMOVEFILESUNITTEST *_test.hpp)
list(REMOVE_ITEM sfm_files_header ${REMOVEFILESUNITTEST})

ADD_LIBRARY(openMVG_sfm ${sfm_files_header} ${sfm_files_cpp})
TARGET_LINK_LIBRARIES(openMVG_sfm ceres openMVG_multiview)
INSTALL(TARGETS openMVG_sfm DESTINATION lib/)

UNIT_TEST(openMVG sfm_data_io "stlplus;openMVG_features;openMVG_sfm")
UNIT_TEST(openMVG sfm_data_BA
  "openMVG_multiview_test_data;openMVG_features;openMVG_multiview;openMVG_sfm;stlplus")

ADD_SUBDIRECTORY(pipelines)

