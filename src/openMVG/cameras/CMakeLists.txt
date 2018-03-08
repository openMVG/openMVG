
add_library(openMVG_camera INTERFACE)
target_link_libraries(openMVG_camera INTERFACE openMVG_numeric openMVG_geometry cereal)
install(TARGETS openMVG_camera DESTINATION lib EXPORT openMVG-targets)

add_library(openMVG_camera_test INTERFACE)
target_link_libraries(openMVG_camera_test INTERFACE openMVG_camera openMVG_multiview cereal)

UNIT_TEST(openMVG Camera_IO "openMVG_camera_test;stlplus")

UNIT_TEST(openMVG Camera_Pinhole "openMVG_camera_test")

UNIT_TEST(openMVG Camera_Pinhole_Radial "openMVG_camera_test")

UNIT_TEST(openMVG Camera_Pinhole_Brown "openMVG_camera_test")

UNIT_TEST(openMVG Camera_Pinhole_Fisheye "openMVG_camera_test")

UNIT_TEST(openMVG Camera_Spherical "openMVG_camera_test")

UNIT_TEST(openMVG Camera_Subset_Parametrization "openMVG_camera_test")
