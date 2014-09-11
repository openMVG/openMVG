
ADD_EXECUTABLE(openMVG_sample_undistoBrown undistoBrown.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_undistoBrown
  openMVG_numeric
  openMVG_image
  stlplus
  ${OPENMVG_LIBRARY_DEPENDENCIES})

SET_PROPERTY(TARGET openMVG_sample_undistoBrown PROPERTY FOLDER OpenMVG/Samples)
