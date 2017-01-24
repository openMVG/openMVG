
add_definitions(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(openMVG_sample_features_affine features_affine_demo.cpp)
target_link_libraries(openMVG_sample_features_affine
  openMVG_image
  openMVG_features
  openMVG_system
  stlplus
  ${OPENMVG_LIBRARY_DEPENDENCIES})

set_property(TARGET openMVG_sample_features_affine PROPERTY FOLDER OpenMVG/Samples)
