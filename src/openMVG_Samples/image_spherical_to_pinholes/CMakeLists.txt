
add_executable(openMVG_sample_pano_converter main_pano_converter.cpp)
target_link_libraries(openMVG_sample_pano_converter
  openMVG_image
  openMVG_multiview
  openMVG_numeric
  stlplus)

set_property(TARGET openMVG_sample_pano_converter PROPERTY FOLDER OpenMVG/Samples)
