
file(
  GLOB_RECURSE
  cameras_files_header
  *.hpp
)
file(
  GLOB_RECURSE
  cameras_files_cpp
  *.cpp
)
file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM cameras_files_cpp ${REMOVEFILESUNITTEST})

add_library(openMVG_cameras ${cameras_files_header} ${cameras_files_cpp})
set_target_properties(openMVG_cameras PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

target_link_libraries(openMVG_cameras openMVG_multiview openMVG_geometry)

install(TARGETS openMVG_cameras DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG Camera_IO "stlplus;openMVG_cameras")

UNIT_TEST(openMVG Camera_Pinhole "openMVG_cameras")

UNIT_TEST(openMVG Camera_Pinhole_Radial "openMVG_cameras")

UNIT_TEST(openMVG Camera_Pinhole_Brown "openMVG_cameras")

UNIT_TEST(openMVG Camera_Pinhole_Fisheye "openMVG_cameras")

UNIT_TEST(openMVG Camera_Spherical "openMVG_cameras")

UNIT_TEST(openMVG Camera_Subset_Parametrization "openMVG_cameras")
