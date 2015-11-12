
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_robustFundamental robust_fundamental.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_robustFundamental
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  stlplus
  vlsift)

SET_PROPERTY(TARGET openMVG_sample_robustFundamental PROPERTY FOLDER OpenMVG/Samples)
