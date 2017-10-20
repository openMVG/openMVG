
add_executable(OMVGsampl_panoConv main_pano_converter.cpp)
target_link_libraries(OMVGsampl_panoConv
  openMVG_numeric
  openMVG_image
  stlplus)

set_property(TARGET OMVGsampl_panoConv PROPERTY FOLDER OpenMVG/Samples)
