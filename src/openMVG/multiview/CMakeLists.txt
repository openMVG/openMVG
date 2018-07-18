
file(
  GLOB_RECURSE
  multiview_files_header
  *.hpp
)
file(
  GLOB_RECURSE
  multiview_files_cpp
  *.cpp
)
file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
file(GLOB_RECURSE REMOVEFILESUNITTEST_HEADER *_test.hpp)
file(GLOB_RECURSE MULTIVIEWTESTDATA *_data_sets*)

#Remove the unit test files (not been used by the library)
list(REMOVE_ITEM multiview_files_cpp ${REMOVEFILESUNITTEST})
list(REMOVE_ITEM multiview_files_header ${REMOVEFILESUNITTEST_HEADER})
list(REMOVE_ITEM multiview_files_cpp ${MULTIVIEWTESTDATA})

add_library(openMVG_multiview
  ${multiview_files_header}
  ${multiview_files_cpp})
target_compile_features(openMVG_multiview INTERFACE ${CXX11_FEATURES})
set_property(TARGET openMVG_multiview PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_multiview DESTINATION lib EXPORT openMVG-targets)
target_link_libraries(openMVG_multiview
  PUBLIC
    openMVG_numeric
  PRIVATE
    openMVG_graph
    ${CERES_LIBRARIES}
)
target_include_directories(openMVG_multiview
  PRIVATE
    ${CERES_INCLUDE_DIRS}
  PUBLIC
    $<INSTALL_INTERFACE:include>
)
set_target_properties(openMVG_multiview PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

add_library(openMVG_multiview_test_data ${MULTIVIEWTESTDATA})
target_link_libraries(openMVG_multiview_test_data PRIVATE openMVG_numeric openMVG_multiview)
set_property(TARGET openMVG_multiview_test_data PROPERTY FOLDER OpenMVG/OpenMVG)

# Triangulation routines
UNIT_TEST(openMVG triangulation "openMVG_multiview_test_data;openMVG_multiview")
UNIT_TEST(openMVG triangulation_nview "openMVG_multiview_test_data;openMVG_multiview")

# Solvers
UNIT_TEST(openMVG solver_affine "openMVG_multiview")
UNIT_TEST(openMVG solver_essential_eight_point "openMVG_multiview_test_data;openMVG_multiview")
UNIT_TEST(openMVG solver_essential_five_point "openMVG_multiview_test_data;openMVG_multiview")
UNIT_TEST(openMVG solver_essential_three_point "openMVG_multiview")
UNIT_TEST(openMVG solver_fundamental_kernel "openMVG_multiview")
UNIT_TEST(openMVG solver_homography_kernel "openMVG_multiview_test_data;openMVG_multiview")
UNIT_TEST(openMVG solver_resection_kernel "openMVG_multiview_test_data;openMVG_multiview")
UNIT_TEST(openMVG solver_translation_knownRotation_kernel "openMVG_multiview_test_data;openMVG_multiview")

# Averaging routines
UNIT_TEST(openMVG rotation_averaging "openMVG_multiview_test_data;openMVG_multiview")
UNIT_TEST(openMVG translation_averaging "openMVG_multiview_test_data;openMVG_multiview")
