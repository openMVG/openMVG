
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
target_compile_features(openMVG_geometry INTERFACE ${CXX11_FEATURES})
set_target_properties(openMVG_geometry PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_geometry PROPERTY FOLDER OpenMVG/OpenMVG)

target_link_libraries(openMVG_geometry
  PUBLIC
    openMVG_numeric
  PRIVATE
    openMVG_linearProgramming
)
target_include_directories(openMVG_geometry
  PUBLIC
    $<INSTALL_INTERFACE:include>
)

install(TARGETS openMVG_geometry DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG convex_hull "openMVG_geometry")
UNIT_TEST(openMVG rigid_transformation3D_srt "openMVG_geometry")
UNIT_TEST(openMVG plane_estimation_kernel "openMVG_geometry")
UNIT_TEST(openMVG half_space_intersection "openMVG_geometry")
UNIT_TEST(openMVG frustum_intersection "openMVG_multiview_test_data;openMVG_multiview;openMVG_geometry")
UNIT_TEST(openMVG frustum_box_intersection "openMVG_multiview_test_data;openMVG_multiview;openMVG_geometry")
