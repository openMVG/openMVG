
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
ADD_EXECUTABLE(openMVG_sample_kvld kvld_filter.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_kvld
  openMVG_image
  openMVG_multiview
  vlsift
  stlplus
  openMVG_kvld
  )

SET_PROPERTY(TARGET openMVG_sample_kvld PROPERTY FOLDER OpenMVG/Samples)
