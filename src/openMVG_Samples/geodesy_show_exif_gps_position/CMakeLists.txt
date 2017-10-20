
add_executable(OMVGsmpl_geoShowExifGps show_exif_gps_position_demo.cpp)
target_link_libraries(OMVGsmpl_geoShowExifGps
  easyexif
  stlplus)
target_compile_definitions(OMVGsmpl_geoShowExifGps
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET OMVGsmpl_geoShowExifGps PROPERTY FOLDER OpenMVG/Samples)
