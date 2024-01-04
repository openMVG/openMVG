
add_executable(openMVG_sample_features_siftPutative siftmatch.cpp)
target_link_libraries(openMVG_sample_features_siftPutative
  openMVG_image
  openMVG_features
  openMVG_matching
  ${STLPLUS_LIBRARY})
target_compile_definitions(openMVG_sample_features_siftPutative
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET openMVG_sample_features_siftPutative PROPERTY FOLDER OpenMVG/Samples)
