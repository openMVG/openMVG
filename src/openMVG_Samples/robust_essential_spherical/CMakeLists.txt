
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_robustEssential_spherical robust_essential_spherical.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_robustEssential_spherical
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  openMVG_system
  openMVG_sfm
  stlplus
  vlsift
  )

SET_PROPERTY(TARGET openMVG_sample_robustEssential_spherical PROPERTY FOLDER OpenMVG/Samples)

