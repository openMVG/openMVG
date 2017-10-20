
add_executable(OMVGsmpl_mvRobustHomographyGuided robust_homography_guided.cpp)
target_link_libraries(OMVGsmpl_mvRobustHomographyGuided
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  stlplus)
target_compile_definitions(OMVGsmpl_mvRobustHomographyGuided
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET OMVGsmpl_mvRobustHomographyGuided PROPERTY FOLDER OpenMVG/Samples)
