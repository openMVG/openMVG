
ADD_EXECUTABLE(openMVG_sample_undistoBrown undistoBrown.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_undistoBrown
  ${OpenMVG_LIBRARIES})
SET_PROPERTY(TARGET openMVG_sample_undistoBrown PROPERTY FOLDER OpenMVG/Samples)
