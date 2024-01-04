
add_executable(openMVG_sample_main_features_repeatability_dataset main_repeatability_dataset.cpp)
target_link_libraries(openMVG_sample_main_features_repeatability_dataset
  openMVG_image
  openMVG_features
  openMVG_matching
  openMVG_system
  openMVG_multiview
  vlsift
  ${STLPLUS_LIBRARY})

set_property(TARGET openMVG_sample_main_features_repeatability_dataset PROPERTY FOLDER OpenMVG/Samples)
