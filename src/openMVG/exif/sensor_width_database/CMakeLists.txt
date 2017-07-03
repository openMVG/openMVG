
UNIT_TEST(openMVG ParseDatabase "stlplus")
target_compile_definitions(openMVG_test_ParseDatabase
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

# Installation rules
install(FILES sensor_width_camera_database.txt DESTINATION share/openMVG)
