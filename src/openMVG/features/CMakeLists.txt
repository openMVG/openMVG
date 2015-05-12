
project(features)

file(
	GLOB_RECURSE
	features_files_headers
	*.hpp
)

file(
	GLOB_RECURSE
	features_files_sources
	*.cpp
)

file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM features_files_sources ${REMOVEFILESUNITTEST})
			
set_source_files_properties(${features_files_sources} PROPERTIES LANGUAGE CXX)
ADD_LIBRARY(openMVG_features ${features_files_sources} ${features_files_headers})
SET_PROPERTY(TARGET openMVG_features PROPERTY FOLDER OpenMVG/OpenMVG)

UNIT_TEST(openMVG features "openMVG_features")

