
add_executable(OMVGsmpl_featSiftPutative siftmatch.cpp)
target_link_libraries(OMVGsmpl_featSiftPutative
  openMVG_image
  openMVG_features
  openMVG_matching
  stlplus)
target_compile_definitions(OMVGsmpl_featSiftPutative
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET OMVGsmpl_featSiftPutative PROPERTY FOLDER OpenMVG/Samples)
