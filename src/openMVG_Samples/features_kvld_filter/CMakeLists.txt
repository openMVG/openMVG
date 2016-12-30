
add_definitions(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(openMVG_sample_features_kvld kvld_filter.cpp)
target_link_libraries(openMVG_sample_features_kvld
  openMVG_image
  openMVG_multiview
  openMVG_kvld
  openMVG_features
  openMVG_matching
  stlplus
  vlsift)

set_property(TARGET openMVG_sample_features_kvld PROPERTY FOLDER OpenMVG/Samples)
