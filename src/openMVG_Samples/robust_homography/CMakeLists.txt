
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_robustHomography robust_homography.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_robustHomography
  openMVG_image
  openMVG_multiview
  openMVG_features
  ${OPENMVG_LIBRARY_DEPENDENCIES}
  vlsift
  stlplus
  ${CLP_LIBRARIES} ${COINUTILS_LIBRARY}
  ${CERES_LIBRARIES})
SET_PROPERTY(TARGET openMVG_sample_robustHomography PROPERTY FOLDER OpenMVG/Samples)
