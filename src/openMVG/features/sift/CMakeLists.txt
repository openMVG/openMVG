
UNIT_TEST(openMVG sift "openMVG_image;openMVG_features;openMVG_system")
target_compile_definitions(openMVG_test_sift
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
