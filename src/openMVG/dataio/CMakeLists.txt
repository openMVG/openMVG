file(
  GLOB_RECURSE
  dataio_header
  *.hpp
)
file(
  GLOB_RECURSE
  dataio_cpp
  *.cpp
)

#file(GLOB_RECURSE REMOVE_FILES_UNIT_TEST *_test.cpp)
#list(REMOVE_ITEM dataio_cpp ${REMOVE_FILES_UNIT_TEST})

set_source_files_properties(${dataio_cpp} PROPERTIES LANGUAGE CXX)
ADD_LIBRARY(openMVG_dataio
  ${dataio_header}
  ${dataio_cpp})
SET_PROPERTY(TARGET openMVG_dataio PROPERTY FOLDER OpenMVG/OpenMVG)
TARGET_LINK_LIBRARIES( openMVG_dataio openMVG_sfm ${ABC_LIBRARIES})
INSTALL(TARGETS openMVG_dataio DESTINATION lib/ EXPORT openMVG-targets)

#UNIT_TEST(openMVG dataio "openMVG_dataio")