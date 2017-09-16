
file(
  GLOB_RECURSE
  geometry_files_header
  *.hpp
)
file(
  GLOB_RECURSE
  geometry_files_cpp
  *.cpp
)
file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM geometry_files_cpp ${REMOVEFILESUNITTEST})

add_library(openMVG_geometry ${geometry_files_header} ${geometry_files_cpp})
set_target_properties(openMVG_geometry PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

target_link_libraries(openMVG_geometry openMVG_linearProgramming)

install(TARGETS openMVG_geometry DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG rigid_transformation3D_srt "openMVG_numeric;openMVG_geometry")

UNIT_TEST(openMVG half_space_intersection "openMVG_numeric;openMVG_linearProgramming;openMVG_geometry")

UNIT_TEST(openMVG frustum_intersection "openMVG_numeric;openMVG_multiview_test_data;openMVG_multiview;openMVG_linearProgramming;openMVG_geometry")

UNIT_TEST(openMVG frustum_box_intersection "openMVG_numeric;openMVG_multiview_test_data;openMVG_multiview;openMVG_linearProgramming;openMVG_geometry")
