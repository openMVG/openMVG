
add_executable(openMVG_sample_geodesy_show_exif_gps_position show_exif_gps_position_demo.cpp)
target_link_libraries(openMVG_sample_geodesy_show_exif_gps_position
  PRIVATE
    openMVG_exif
    openMVG_geodesy
    ${STLPLUS_LIBRARY})
target_compile_definitions(openMVG_sample_geodesy_show_exif_gps_position
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET openMVG_sample_geodesy_show_exif_gps_position PROPERTY FOLDER OpenMVG/Samples)
