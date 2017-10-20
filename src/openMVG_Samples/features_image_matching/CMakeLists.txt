
add_executable(OMVGsmpl_imgMatching describe_and_match.cpp)
target_link_libraries(OMVGsmpl_imgMatching
  openMVG_image
  openMVG_features
  openMVG_matching
  stlplus
  vlsift)
target_compile_definitions(OMVGsmpl_imgMatching
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

set_property(TARGET OMVGsmpl_imgMatching PROPERTY FOLDER OpenMVG/Samples)
