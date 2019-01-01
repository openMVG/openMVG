
add_library(openMVG_geodesy INTERFACE)
target_compile_features(openMVG_geodesy INTERFACE ${CXX11_FEATURES})
set_property(TARGET openMVG_camera APPEND PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES "$<INSTALL_INTERFACE:include>")
target_link_libraries(openMVG_geodesy
  INTERFACE
    openMVG_numeric
)
install(TARGETS openMVG_geodesy DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG_geodesy geodesy "openMVG_geodesy")
