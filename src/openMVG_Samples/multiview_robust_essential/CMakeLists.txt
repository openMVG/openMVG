
add_executable(OMVGsmpl_mvRobustEssential robust_essential.cpp)
target_link_libraries(OMVGsmpl_mvRobustEssential
  openMVG_image
  openMVG_multiview
  openMVG_system
  openMVG_features
  openMVG_geometry
  openMVG_matching
  openMVG_sfm)
target_compile_definitions(OMVGsmpl_mvRobustEssential
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET OMVGsmpl_mvRobustEssential PROPERTY FOLDER OpenMVG/Samples)
