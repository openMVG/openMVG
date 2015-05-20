
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_robustEssential robust_essential.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_robustEssential 
  ${OpenMVG_LIBS})

SET_PROPERTY(TARGET openMVG_sample_robustEssential PROPERTY FOLDER OpenMVG/Samples)
