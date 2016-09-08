
ADD_EXECUTABLE(openMVG_sample_pano_converter main_pano_converter.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_pano_converter
  openMVG_numeric
  openMVG_image
  stlplus)

