
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_robustEssential_spherical robust_essential_spherical.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_robustEssential_spherical
  openMVG_image
  openMVG_multiview
  vlsift
  stlplus
  flann_cpp_s
  ${OPENMVG_LIBRARY_DEPENDENCIES})

  SET_PROPERTY(TARGET openMVG_sample_robustEssential_spherical PROPERTY FOLDER OpenMVG/Samples)

