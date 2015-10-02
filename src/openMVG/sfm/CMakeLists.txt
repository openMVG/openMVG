
file(
	GLOB_RECURSE
	sfm_files_header
	*.hpp
)
file(
	GLOB_RECURSE
	sfm_files_cpp
	*.cpp
)
file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM sfm_files_cpp ${REMOVEFILESUNITTEST})
#Remove the test header main files
file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.hpp)
list(REMOVE_ITEM sfm_files_header ${REMOVEFILESUNITTEST})

add_library(openMVG_sfm ${sfm_files_header} ${sfm_files_cpp})
target_link_libraries(
  openMVG_sfm
  openMVG_multiview
  stlplus
  ${CERES_LIBRARIES}
  openMVG_lInftyComputerVision
  openMVG_system
  openMVG_matching
)
set_target_properties(openMVG_sfm PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

install(TARGETS openMVG_sfm DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG sfm_data_io
  "openMVG_features;openMVG_sfm;openMVG_system;stlplus")
UNIT_TEST(openMVG sfm_data_BA
  "openMVG_multiview_test_data;openMVG_features;openMVG_multiview;openMVG_sfm;openMVG_system;stlplus")
UNIT_TEST(openMVG sfm_data_utils
  "openMVG_features;openMVG_multiview;openMVG_system;openMVG_sfm;stlplus")

add_subdirectory(pipelines)

