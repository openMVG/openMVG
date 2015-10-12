
project(matching)

file(
	GLOB
	matching_files_header
	*.hpp
)
file(
	GLOB
	matching_files_cpp
	*.cpp
)
file(GLOB REMOVEFILESUNITTEST *_test.cpp)

#Remove the future main files
list(REMOVE_ITEM matching_files_cpp ${REMOVEFILESUNITTEST})

set_source_files_properties(${matching_files_cpp} PROPERTIES LANGUAGE CXX)
add_library(openMVG_matching
  ${matching_files_header}
  ${matching_files_cpp})
set_target_properties(openMVG_matching PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
target_link_libraries(openMVG_matching ${FLANN_LIBRARY})
set_property(TARGET openMVG_matching PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_matching DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG matching "stlplus")
UNIT_TEST(openMVG matching_filters "stlplus")
UNIT_TEST(openMVG indMatch "stlplus")
UNIT_TEST(openMVG metric "")

add_subdirectory(kvld)
