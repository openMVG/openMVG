file(
  GLOB
  global_headers
  *.hpp)

file(
  GLOB
  global_cpp
  *.cpp)

file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM global_cpp ${REMOVEFILESUNITTEST})

add_library(openMVG_pipeline_global ${global_headers} ${global_cpp})

target_link_libraries(
  openMVG_pipeline_global
  PRIVATE
    openMVG_sfm)

target_include_directories(
  openMVG_pipeline_global
  PRIVATE
  ${CERES_INCLUDE_DIRS})

set_target_properties(openMVG_pipeline_global PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

install(TARGETS openMVG_pipeline_global DESTINATION lib EXPORT openMVG-targets)
UNIT_TEST(openMVG global_SfM
  "openMVG_multiview_test_data;openMVG_features;openMVG_multiview;openMVG_pipeline_global;openMVG_system;stlplus")

