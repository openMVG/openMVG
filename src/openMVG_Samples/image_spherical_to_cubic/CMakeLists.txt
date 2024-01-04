
add_executable(openMVG_sample_pano_spherical_to_cubic main_pano_spherical_to_cubic.cpp)
target_link_libraries(openMVG_sample_pano_spherical_to_cubic
  openMVG_camera
  openMVG_image
  ${STLPLUS_LIBRARY})

set_property(TARGET openMVG_sample_pano_spherical_to_cubic PROPERTY FOLDER OpenMVG/Samples)
