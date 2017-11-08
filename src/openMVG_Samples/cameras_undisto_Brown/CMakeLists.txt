
add_executable(OMVGsmpl_camUndistBrown undistoBrown.cpp)
target_link_libraries(OMVGsmpl_camUndistBrown
  openMVG_geometry
  openMVG_image
  openMVG_multiview
  stlplus
  ${CERES_LIBRARIES})
set_property(TARGET OMVGsmpl_camUndistBrown PROPERTY FOLDER OpenMVG/Samples)
