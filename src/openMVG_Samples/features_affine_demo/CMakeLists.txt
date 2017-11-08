
add_executable(OMVGsmpl_featAffine features_affine_demo.cpp)
target_link_libraries(OMVGsmpl_featAffine
  openMVG_image
  openMVG_features
  openMVG_system
  stlplus)
target_compile_definitions(OMVGsmpl_featAffine
    PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET OMVGsmpl_featAffine PROPERTY FOLDER OpenMVG/Samples)
