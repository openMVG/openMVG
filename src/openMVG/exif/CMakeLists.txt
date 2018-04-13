
add_library(openMVG_exif exif_IO_EasyExif.cpp)
target_link_libraries(openMVG_exif PRIVATE easyexif)
target_include_directories(openMVG_exif PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/third_party/include/easyexif)
install(TARGETS openMVG_exif DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG exif_IO "openMVG_testing;stlplus;openMVG_exif")

add_subdirectory(sensor_width_database)
