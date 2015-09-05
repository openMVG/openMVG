
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_siftPutative siftmatch.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_siftPutative
  openMVG_image
  openMVG_features
  openMVG_matching
  vlsift
  stlplus)

SET_PROPERTY(TARGET openMVG_sample_siftPutative PROPERTY FOLDER OpenMVG/Samples)
