add_definitions(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
UNIT_TEST(openMVG exif_IO "stlplus;easyexif")

add_subdirectory(sensor_width_database)
