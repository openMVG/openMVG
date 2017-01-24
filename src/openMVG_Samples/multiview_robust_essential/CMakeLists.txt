
add_definitions(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(openMVG_sample_multiview_robustEssential robust_essential.cpp)
target_link_libraries(openMVG_sample_multiview_robustEssential
  openMVG_image
  openMVG_multiview
  openMVG_system
  openMVG_features
  openMVG_matching
  openMVG_sfm
  vlsift)

set_property(TARGET openMVG_sample_multiview_robustEssential PROPERTY FOLDER OpenMVG/Samples)
