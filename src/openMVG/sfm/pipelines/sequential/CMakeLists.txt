file(
  GLOB
  sequential_headers
  *.hpp)

file(
  GLOB
  sequential_cpp
  *.cpp)

file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM sequential_cpp ${REMOVEFILESUNITTEST})

add_library(openMVG_pipeline_sequential ${sequential_headers} ${sequential_cpp})

target_link_libraries(
  openMVG_pipeline_sequential
  PRIVATE
    openMVG_sfm)

target_include_directories(
  openMVG_pipeline_sequential
  PRIVATE
  ${CERES_INCLUDE_DIRS})

set_target_properties(openMVG_pipeline_sequential PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

install(TARGETS openMVG_pipeline_sequential DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG sequential_SfM
  "openMVG_multiview_test_data;openMVG_features;openMVG_multiview;openMVG_pipeline_sequential;openMVG_system;stlplus")

