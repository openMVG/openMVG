
add_executable(openMVG_sample_cameras_undistoBrown undistoBrown.cpp)
target_link_libraries(openMVG_sample_cameras_undistoBrown
  openMVG_geometry
  openMVG_image
  openMVG_multiview
  stlplus
  ${CERES_LIBRARIES})
set_property(TARGET openMVG_sample_cameras_undistoBrown PROPERTY FOLDER OpenMVG/Samples)
