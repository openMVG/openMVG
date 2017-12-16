
add_library(openMVG_exif INTERFACE)
target_link_libraries(openMVG_exif INTERFACE easyexif)
install(TARGETS openMVG_exif DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG exif_IO "openMVG_testing;stlplus;openMVG_exif")

add_subdirectory(sensor_width_database)
