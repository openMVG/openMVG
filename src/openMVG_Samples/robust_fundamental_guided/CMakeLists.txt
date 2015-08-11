
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_robustFundamental_guided robust_fundamental_guided.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_robustFundamental_guided 
  openMVG_image
  openMVG_multiview
  openMVG_features
  stlplus
  vlsift
  ${CERES_LIBRARIES})

SET_PROPERTY(TARGET openMVG_sample_robustFundamental_guided PROPERTY FOLDER OpenMVG/Samples)
