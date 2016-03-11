
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_features_affine features_affine_demo.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_features_affine
  openMVG_image
  openMVG_features
  openMVG_system
  stlplus
  ${OPENMVG_LIBRARY_DEPENDENCIES})

SET_PROPERTY(TARGET openMVG_sample_features_affine PROPERTY FOLDER OpenMVG/Samples)
