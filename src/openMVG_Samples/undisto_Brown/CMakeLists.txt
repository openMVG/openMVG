
ADD_EXECUTABLE(openMVG_sample_undistoBrown undistoBrown.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_undistoBrown
  openMVG_image
  openMVG_multiview
  stlplus
  ${CERES_LIBRARIES})
SET_PROPERTY(TARGET openMVG_sample_undistoBrown PROPERTY FOLDER OpenMVG/Samples)
