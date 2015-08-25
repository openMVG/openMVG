ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

project(exif)

file(
	GLOB_RECURSE
	exif_files_headers
	*.hpp
)

file(
	GLOB_RECURSE
	exif_files_sources
	*.cpp
)

file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM exif_files_sources ${REMOVEFILESUNITTEST})
			
set_source_files_properties(${exif_files_sources} PROPERTIES LANGUAGE CXX)
ADD_LIBRARY(openMVG_exif ${exif_files_sources} ${exif_files_headers})
TARGET_LINK_LIBRARIES(openMVG_exif stlplus easyexif)

SET_TARGET_PROPERTIES(openMVG_exif PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
INSTALL(TARGETS openMVG_exif DESTINATION lib EXPORT openMVG-targets)
SET_PROPERTY(TARGET openMVG_exif PROPERTY FOLDER OpenMVG/OpenMVG)

UNIT_TEST(openMVG exif_IO "openMVG_exif")

ADD_SUBDIRECTORY(sensor_width_database)
