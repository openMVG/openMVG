

ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

UNIT_TEST(openMVG ParseDatabase "stlplus")

# Installation rules
INSTALL(FILES sensor_width_camera_database.txt DESTINATION share/openMVG)
