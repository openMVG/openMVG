add_definitions(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

UNIT_TEST(openMVG sift "openMVG_image;openMVG_features;openMVG_system")
