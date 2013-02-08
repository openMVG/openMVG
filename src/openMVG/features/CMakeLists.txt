
project(features)

file(
	GLOB_RECURSE
	features_files
	*.hpp
)
			
set_source_files_properties(${features_files} PROPERTIES LANGUAGE CXX)

UNIT_TEST(openMVG features "")

