
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_akaze akaze_demo.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_akaze
  openMVG_image
  openMVG_features
  stlplus
  ${OPENMVG_LIBRARY_DEPENDENCIES})
SET_PROPERTY(TARGET openMVG_sample_akaze PROPERTY FOLDER OpenMVG/Samples)

ADD_EXECUTABLE( openMVG_main_computeAKAZE main_computeAKAZE.cpp )
TARGET_LINK_LIBRARIES( openMVG_main_computeAKAZE
  openMVG_image
  openMVG_features
  openMVG_system
  stlplus
  ${OPENMVG_LIBRARY_DEPENDENCIES})
SET_PROPERTY(TARGET openMVG_main_computeAKAZE PROPERTY FOLDER OpenMVG/Samples)

