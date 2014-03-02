
project(matching)

file(
	GLOB_RECURSE
	matching_files_header
	*.hpp
)
file(
	GLOB_RECURSE
	matching_files_cpp
	*.cpp
)
file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)

#Remove the future main files
list(REMOVE_ITEM matching_files_cpp ${REMOVEFILESUNITTEST})

set_source_files_properties(${matching_files_header} PROPERTIES LANGUAGE CXX)
set_source_files_properties(${matching_files_cpp} PROPERTIES LANGUAGE CXX)

UNIT_TEST(openMVG matching "flann_cpp_s")
UNIT_TEST(openMVG matching_filters "")
UNIT_TEST(openMVG indMatch "")

ADD_SUBDIRECTORY(kvld)
