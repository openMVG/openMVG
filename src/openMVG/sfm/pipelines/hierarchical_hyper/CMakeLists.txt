file(
  GLOB
  hierarchical_hyper_headers
  *.hpp)

file(
  GLOB
  hierarchical_hyper_cpp
  *.cpp)

file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM hierarchical_hyper_cpp ${REMOVEFILESUNITTEST})

add_library(openMVG_pipeline_hierarchical_hyper ${hierarchical_hyper_headers} ${hierarchical_hyper_cpp})

target_link_libraries(
  openMVG_pipeline_hierarchical_hyper
  PRIVATE
    openMVG_sfm
    openMVG_pipeline_sequential
    scotch
    scotcherr
    )

target_include_directories(
  openMVG_pipeline_hierarchical_hyper
  PRIVATE
  ${CERES_INCLUDE_DIRS})

set_target_properties(openMVG_pipeline_hierarchical_hyper PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

install(TARGETS openMVG_pipeline_hierarchical_hyper DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG scene_aligner
  "openMVG_features;openMVG_multiview_test_data;openMVG_pipeline_hierarchical_hyper")
UNIT_TEST(openMVG submap_BA_fixed_separators
  "openMVG_features;openMVG_multiview_test_data;openMVG_pipeline_hierarchical_hyper")
UNIT_TEST(openMVG hypercluster "openMVG_pipeline_hierarchical_hyper")
