IF(OpenMVG_USE_BOOST)
file(
  GLOB_RECURSE
  voctree_header
  *.hpp
)
file(
  GLOB_RECURSE
  voctree_cpp
  *.cpp
)

ADD_LIBRARY(openMVG_voctree
  ${voctree_header}
  ${voctree_cpp})
SET_PROPERTY(TARGET openMVG_voctree PROPERTY FOLDER OpenMVG/OpenMVG)
TARGET_LINK_LIBRARIES( openMVG_voctree openMVG_features openMVG_sfm ${Boost_LIBRARIES})
INSTALL(TARGETS openMVG_voctree DESTINATION lib/ EXPORT openMVG-targets)

ENDIF(OpenMVG_USE_BOOST)
