
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_robust_estimation_Homography robust_estimation_tutorial_homography.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_robust_estimation_Homography
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  vlsift
  stlplus)
SET_PROPERTY(TARGET openMVG_sample_robust_estimation_Homography PROPERTY FOLDER OpenMVG/Samples)
