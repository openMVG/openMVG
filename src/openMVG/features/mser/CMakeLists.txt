
UNIT_TEST(openMVG mser "openMVG_image;openMVG_features")
if (OpenMVG_BUILD_TESTS)
  target_compile_definitions(openMVG_test_mser
      PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../image/")
endif (OpenMVG_BUILD_TESTS)