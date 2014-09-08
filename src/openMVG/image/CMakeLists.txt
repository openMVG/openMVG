
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

file(
	GLOB_RECURSE
	image_files_header
	*.hpp
)
file(
	GLOB_RECURSE
	image_files_cpp
	*.cpp
)
file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)

#Remove the future main files
list(REMOVE_ITEM image_files_cpp ${REMOVEFILESUNITTEST})

ADD_LIBRARY(openMVG_image ${image_files_header} ${image_files_cpp})
TARGET_LINK_LIBRARIES(openMVG_image ${PNG_LIBRARIES} ${JPEG_LIBRARY} ${TIFF_LIBRARIES})

UNIT_TEST(openMVG image "openMVG_numeric;openMVG_image")
UNIT_TEST(openMVG image_drawing "openMVG_numeric;openMVG_image")
UNIT_TEST(openMVG image_io "openMVG_numeric;openMVG_image")

