
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_describe_and_match describe_and_match.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_describe_and_match
  openMVG_features
  ${OpenMVG_LIBS}
  vlsift
  stlplus
  flann_cpp_s  
  ${OPENMVG_LIBRARY_DEPENDENCIES})

SET_PROPERTY(TARGET openMVG_sample_describe_and_match PROPERTY FOLDER OpenMVG/Samples)
