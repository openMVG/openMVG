
add_executable(openMVG_sample_multiview_robustEssential robust_essential.cpp)
target_link_libraries(openMVG_sample_multiview_robustEssential
  openMVG_image
  openMVG_multiview
  openMVG_system
  openMVG_features
  openMVG_geometry
  openMVG_matching
  openMVG_sfm)
target_compile_definitions(openMVG_sample_multiview_robustEssential
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET openMVG_sample_multiview_robustEssential PROPERTY FOLDER OpenMVG/Samples)
