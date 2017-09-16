
add_executable(openMVG_sample_features_kvld kvld_filter.cpp)
target_link_libraries(openMVG_sample_features_kvld
  openMVG_image
  openMVG_multiview
  openMVG_kvld
  openMVG_features
  openMVG_matching
  stlplus)
target_compile_definitions(openMVG_sample_features_kvld
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET openMVG_sample_features_kvld PROPERTY FOLDER OpenMVG/Samples)
