
add_library(openMVG_camera INTERFACE)
set_property(TARGET openMVG_camera APPEND PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES "$<INSTALL_INTERFACE:include>")

target_compile_features(openMVG_camera INTERFACE ${CXX11_FEATURES})
target_link_libraries(openMVG_camera INTERFACE openMVG_numeric ${OPENMVG_LIBRARY_DEPENDENCIES})
install(TARGETS openMVG_camera DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG Camera_Pinhole openMVG_camera)

UNIT_TEST(openMVG Camera_Pinhole_Radial openMVG_camera)

UNIT_TEST(openMVG Camera_Pinhole_Brown openMVG_camera)

UNIT_TEST(openMVG Camera_Pinhole_Fisheye openMVG_camera)

UNIT_TEST(openMVG Camera_Spherical openMVG_camera)

UNIT_TEST(openMVG Camera_Subset_Parametrization openMVG_camera)

add_library(openMVG_camera_test INTERFACE)
target_link_libraries(openMVG_camera_test INTERFACE openMVG_camera)
target_include_directories(openMVG_camera_test
  INTERFACE $<BUILD_INTERFACE:${CEREAL_INCLUDE_DIRS}>)

UNIT_TEST(openMVG Camera_IO "openMVG_camera_test;${STLPLUS_LIBRARY}")
