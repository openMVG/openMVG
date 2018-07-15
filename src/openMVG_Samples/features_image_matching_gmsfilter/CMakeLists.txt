
add_executable(openMVG_sample_image_matching_gmsfilter describe_and_match_gmsfilter.cpp)
target_link_libraries(openMVG_sample_image_matching_gmsfilter
  openMVG_image
  openMVG_features
  openMVG_matching
  openMVG_multiview
  openMVG_robust_estimation
  ${STLPLUS_LIBRARY}
  )
target_compile_definitions(openMVG_sample_image_matching_gmsfilter
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

set_property(TARGET openMVG_sample_image_matching_gmsfilter PROPERTY FOLDER OpenMVG/Samples)
