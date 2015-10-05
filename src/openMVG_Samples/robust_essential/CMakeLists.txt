
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_robustEssential robust_essential.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_robustEssential
  openMVG_image
  openMVG_multiview
  openMVG_system
  openMVG_features
  openMVG_matching
  openMVG_sfm
  vlsift)

SET_PROPERTY(TARGET openMVG_sample_robustEssential PROPERTY FOLDER OpenMVG/Samples)
