
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
file(GLOB_RECURSE MULTIVIEWTESTDATA *_data_sets*)

#Remove the unit test files (not been used by the library)
list(REMOVE_ITEM multiview_files_cpp ${REMOVEFILESUNITTEST})
list(REMOVE_ITEM multiview_files_cpp ${MULTIVIEWTESTDATA})

ADD_LIBRARY(openMVG_multiview
  ${multiview_files_header}
  ${multiview_files_cpp})
SET_PROPERTY(TARGET openMVG_multiview PROPERTY FOLDER OpenMVG)
INSTALL(TARGETS openMVG_multiview DESTINATION lib/)

#Make multiview library integrate numeric library
ADD_LIBRARY(openMVG_multiview_test_data ${MULTIVIEWTESTDATA})
SET_PROPERTY(TARGET openMVG_multiview_test_data PROPERTY FOLDER OpenMVG/OpenMVG)
TARGET_LINK_LIBRARIES(openMVG_multiview openMVG_numeric)

UNIT_TEST(openMVG triangulation "openMVG_multiview;openMVG_multiview_test_data")
UNIT_TEST(openMVG triangulation_nview "openMVG_multiview;openMVG_multiview_test_data")

UNIT_TEST(openMVG solver_affine "openMVG_multiview")
UNIT_TEST(openMVG solver_fundamental_kernel "openMVG_multiview")
UNIT_TEST(openMVG solver_essential_kernel "openMVG_multiview;openMVG_multiview_test_data")
UNIT_TEST(openMVG solver_homography_kernel "openMVG_multiview;openMVG_multiview_test_data")
UNIT_TEST(openMVG solver_essential_five_point "openMVG_multiview;openMVG_multiview_test_data")
UNIT_TEST(openMVG solver_resection_kernel "openMVG_multiview;openMVG_multiview_test_data")
UNIT_TEST(openMVG rotation_averaging "lemon;openMVG_multiview;openMVG_multiview_test_data")
UNIT_TEST(openMVG translation_averaging "lemon;openMVG_multiview;openMVG_multiview_test_data;ceres")
