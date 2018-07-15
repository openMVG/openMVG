
add_executable(openMVG_sample_image_matching describe_and_match.cpp)
target_link_libraries(openMVG_sample_image_matching
  openMVG_image
  openMVG_features
  openMVG_matching
  ${STLPLUS_LIBRARY}
  vlsift)
target_compile_definitions(openMVG_sample_image_matching
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

set_property(TARGET openMVG_sample_image_matching PROPERTY FOLDER OpenMVG/Samples)
