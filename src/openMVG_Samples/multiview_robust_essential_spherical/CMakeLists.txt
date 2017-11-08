
add_executable(OMVGsmpl_mvRobustEssentialSpherical robust_essential_spherical.cpp)
target_link_libraries(OMVGsmpl_mvRobustEssentialSpherical
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  openMVG_system
  openMVG_sfm
  stlplus)
target_compile_definitions(OMVGsmpl_mvRobustEssentialSpherical
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET OMVGsmpl_mvRobustEssentialSpherical PROPERTY FOLDER OpenMVG/Samples)
