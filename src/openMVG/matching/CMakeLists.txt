
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

ADD_LIBRARY(openMVG_matching
  ${matching_files_header}
  ${matching_files_cpp})
SET_TARGET_PROPERTIES(openMVG_matching PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
TARGET_LINK_LIBRARIES(openMVG_matching ${FLANN_LIBRARY})
SET_PROPERTY(TARGET openMVG_matching PROPERTY FOLDER OpenMVG)
INSTALL(TARGETS openMVG_matching DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG matching "")
UNIT_TEST(openMVG matching_filters "")
UNIT_TEST(openMVG indMatch "")
UNIT_TEST(openMVG metric "")

ADD_SUBDIRECTORY(kvld)
