
add_executable(OMVGsmpl_mvRobustFundamentalGuided robust_fundamental_guided.cpp)
target_link_libraries(OMVGsmpl_mvRobustFundamentalGuided
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  stlplus)
target_compile_definitions(OMVGsmpl_mvRobustFundamentalGuided
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET OMVGsmpl_mvRobustFundamentalGuided PROPERTY FOLDER OpenMVG/Samples)
