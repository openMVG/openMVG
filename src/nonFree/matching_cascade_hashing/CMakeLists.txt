#Add executable that computes:
# - features and descriptors
# - putatives matches (thanks to Cascade Hashing CVPR14)
# - geometric filtered matches
#
ADD_EXECUTABLE(openMVG_main_computeMatchesCasHas
  main_computeMatchesCasHas.cpp
  CasHash.cpp
  CasHash.hpp
  Matcher_CascadeHashing.hpp
)
TARGET_LINK_LIBRARIES(openMVG_main_computeMatchesCasHas
  openMVG_image
  openMVG_multiview
  openMVG_system
  vlsift
  stlplus
  flann_cpp_s)
SET_PROPERTY(TARGET openMVG_main_computeMatchesCasHas PROPERTY FOLDER OpenMVG/nonFree/software)

