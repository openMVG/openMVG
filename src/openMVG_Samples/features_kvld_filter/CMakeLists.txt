
add_executable(OMVGsmpl_featKvld kvld_filter.cpp)
target_link_libraries(OMVGsmpl_featKvld
  openMVG_image
  openMVG_multiview
  openMVG_kvld
  openMVG_features
  openMVG_matching
  stlplus)
target_compile_definitions(OMVGsmpl_featKvld
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET OMVGsmpl_featKvld PROPERTY FOLDER OpenMVG/Samples)
