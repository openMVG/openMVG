
add_executable(openMVG_sample_multiview_robustHomography robust_homography.cpp)
target_link_libraries(openMVG_sample_multiview_robustHomography
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  ${STLPLUS_LIBRARY})
target_compile_definitions(openMVG_sample_multiview_robustHomography
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET openMVG_sample_multiview_robustHomography PROPERTY FOLDER OpenMVG/Samples)
