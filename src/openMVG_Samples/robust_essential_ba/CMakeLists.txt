
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_robustEssential_ba robust_essential_ba.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_robustEssential_ba 
  ${OpenMVG_LIBRARIES})

SET_PROPERTY(TARGET openMVG_sample_robustEssential_ba PROPERTY FOLDER OpenMVG/Samples)
