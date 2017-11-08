
add_executable(OMVGsmpl_robustEstHomography robust_estimation_tutorial_homography.cpp)
target_link_libraries(OMVGsmpl_robustEstHomography
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  stlplus)
target_compile_definitions(OMVGsmpl_robustEstHomography
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET OMVGsmpl_robustEstHomography PROPERTY FOLDER OpenMVG/Samples)
