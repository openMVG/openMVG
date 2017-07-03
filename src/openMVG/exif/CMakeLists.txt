
UNIT_TEST(openMVG exif_IO "stlplus;easyexif")
target_compile_definitions(openMVG_test_exif_IO
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(sensor_width_database)
