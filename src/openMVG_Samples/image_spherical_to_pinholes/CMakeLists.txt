
add_executable(openMVG_sample_pano_converter main_pano_converter.cpp)
target_link_libraries(openMVG_sample_pano_converter
  openMVG_camera
  openMVG_image
  ${STLPLUS_LIBRARY})

set_property(TARGET openMVG_sample_pano_converter PROPERTY FOLDER OpenMVG/Samples)
