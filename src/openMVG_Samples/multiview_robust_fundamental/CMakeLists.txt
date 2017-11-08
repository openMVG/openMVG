
add_executable(OMVGsmpl_mvRobustFundamental robust_fundamental.cpp)
target_link_libraries(OMVGsmpl_mvRobustFundamental
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  stlplus)
target_compile_definitions(OMVGsmpl_mvRobustFundamental
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET OMVGsmpl_mvRobustFundamental PROPERTY FOLDER OpenMVG/Samples)
