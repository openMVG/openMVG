IF(OpenMVG_BUILD_VOCTREE)
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

file(GLOB_RECURSE REMOVE_FILES_UNIT_TEST *_test.cpp)
list(REMOVE_ITEM voctree_cpp ${REMOVE_FILES_UNIT_TEST})

ADD_LIBRARY(openMVG_voctree
  ${voctree_header}
  ${voctree_cpp})
SET_PROPERTY(TARGET openMVG_voctree PROPERTY FOLDER OpenMVG/OpenMVG)
TARGET_LINK_LIBRARIES( openMVG_voctree openMVG_features openMVG_sfm ${Boost_LIBRARIES})
INSTALL(TARGETS openMVG_voctree DESTINATION lib/ EXPORT openMVG-targets)

UNIT_TEST(openMVG voctree "openMVG_voctree")

ENDIF(OpenMVG_BUILD_VOCTREE)
