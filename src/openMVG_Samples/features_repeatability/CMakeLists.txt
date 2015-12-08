
ADD_EXECUTABLE(openMVG_sample_main_repeatability_dataset main_repeatability_dataset.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_main_repeatability_dataset
  openMVG_image
  openMVG_features
  openMVG_matching
  openMVG_system
  openMVG_multiview
  vlsift
  stlplus)

SET_PROPERTY(TARGET openMVG_sample_main_repeatability_dataset PROPERTY FOLDER OpenMVG/Samples)
