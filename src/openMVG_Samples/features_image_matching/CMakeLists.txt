
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_image_matching describe_and_match.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_image_matching
  openMVG_image
  openMVG_features
  openMVG_matching
  stlplus
  vlsift)

SET_PROPERTY(TARGET openMVG_sample_image_matching PROPERTY FOLDER OpenMVG/Samples)
