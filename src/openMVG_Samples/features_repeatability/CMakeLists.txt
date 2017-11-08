
add_executable(OMVGsmpl_mainFeatRepeatDataset main_repeatability_dataset.cpp)
target_link_libraries(OMVGsmpl_mainFeatRepeatDataset
  openMVG_image
  openMVG_features
  openMVG_matching
  openMVG_system
  openMVG_multiview
  vlsift
  stlplus)

set_property(TARGET OMVGsmpl_mainFeatRepeatDataset PROPERTY FOLDER OpenMVG/Samples)
