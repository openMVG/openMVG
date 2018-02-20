
add_library(openMVG_geodesy INTERFACE)
target_link_libraries(openMVG_geodesy
  INTERFACE
    openMVG_numeric
)
install(TARGETS openMVG_geodesy DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG_geodesy geodesy "openMVG_geodesy")
