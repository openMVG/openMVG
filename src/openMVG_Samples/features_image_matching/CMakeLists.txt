
add_definitions(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(openMVG_sample_image_matching describe_and_match.cpp)
target_link_libraries(openMVG_sample_image_matching
  openMVG_image
  openMVG_features
  openMVG_matching
  stlplus
  vlsift)

set_property(TARGET openMVG_sample_image_matching PROPERTY FOLDER OpenMVG/Samples)
