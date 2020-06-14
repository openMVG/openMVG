
add_library(openMVG_spherical INTERFACE)
target_compile_features(openMVG_spherical INTERFACE ${CXX11_FEATURES})
target_link_libraries(openMVG_spherical INTERFACE openMVG_image)
install(TARGETS openMVG_spherical DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG_spherical spherical "openMVG_spherical")
UNIT_TEST(openMVG convert "openMVG_testing;openMVG_image;openMVG_camera")
