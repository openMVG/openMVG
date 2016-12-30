
add_definitions(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(openMVG_sample_robust_estimation_Homography robust_estimation_tutorial_homography.cpp)
target_link_libraries(openMVG_sample_robust_estimation_Homography
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  vlsift
  stlplus)
set_property(TARGET openMVG_sample_robust_estimation_Homography PROPERTY FOLDER OpenMVG/Samples)
