
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

file(GLOB_RECURSE DONOTCOMPILE *.cc)
list(REMOVE_ITEM matching_files_cpp ${DONOTCOMPILE})

set_source_files_properties(${matching_files_header} PROPERTIES LANGUAGE CXX)
set_source_files_properties(${matching_files_cpp} PROPERTIES LANGUAGE CXX)

UNIT_TEST(openMVG matching "flann_cpp_s")
UNIT_TEST(openMVG matching_filters "")
UNIT_TEST(openMVG indMatch "")
UNIT_TEST(openMVG metric "")

ADD_SUBDIRECTORY(kvld)
