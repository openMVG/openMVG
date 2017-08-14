
add_executable(openMVG_sample_multiview_robustEssential_ba robust_essential_ba.cpp)
target_link_libraries(openMVG_sample_multiview_robustEssential_ba
  openMVG_image
  openMVG_features
  openMVG_matching
  openMVG_system
  openMVG_sfm
  stlplus)
target_compile_definitions(openMVG_sample_multiview_robustEssential_ba
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET openMVG_sample_multiview_robustEssential_ba PROPERTY FOLDER OpenMVG/Samples)
