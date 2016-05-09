
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_multiview_robustHomography robust_homography.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_multiview_robustHomography
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  vlsift
  stlplus)
SET_PROPERTY(TARGET openMVG_sample_multiview_robustHomography PROPERTY FOLDER OpenMVG/Samples)
