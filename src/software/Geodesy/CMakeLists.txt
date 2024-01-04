
add_executable(openMVG_main_geodesy_registration_to_gps_position registration_to_exif_gps_position.cpp)
target_link_libraries(openMVG_main_geodesy_registration_to_gps_position
  PRIVATE
    openMVG_exif
    openMVG_features
    openMVG_geodesy
    openMVG_sfm
    openMVG_system
    ${STLPLUS_LIBRARY}
)

set_property(TARGET openMVG_main_geodesy_registration_to_gps_position PROPERTY FOLDER OpenMVG/software)
install(TARGETS openMVG_main_geodesy_registration_to_gps_position DESTINATION bin/)
