
add_executable(openMVG_sample_pano_converter main_pano_converter.cpp)
target_link_libraries(openMVG_sample_pano_converter
  openMVG_numeric
  openMVG_image
  stlplus)

set_property(TARGET openMVG_sample_pano_converter PROPERTY FOLDER OpenMVG/Samples)
