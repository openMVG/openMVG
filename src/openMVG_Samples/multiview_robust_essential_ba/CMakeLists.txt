
add_executable(OMVGsmpl_mvRobustEssentialBA robust_essential_ba.cpp)
target_link_libraries(OMVGsmpl_mvRobustEssentialBA
  openMVG_image
  openMVG_features
  openMVG_matching
  openMVG_system
  openMVG_sfm
  stlplus)
target_compile_definitions(OMVGsmpl_mvRobustEssentialBA
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET OMVGsmpl_mvRobustEssentialBA PROPERTY FOLDER OpenMVG/Samples)
