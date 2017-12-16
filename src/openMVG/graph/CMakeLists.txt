
add_library(openMVG_graph INTERFACE)
target_link_libraries(openMVG_graph INTERFACE ${LEMON_LIBRARY})
install(TARGETS openMVG_graph DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG_graph connectedComponent "openMVG_testing;openMVG_graph")
UNIT_TEST(openMVG_graph triplet_finder "openMVG_testing;openMVG_graph")
