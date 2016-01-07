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
add_library(openMVG_dataio
  ${dataio_header}
  ${dataio_cpp})
set_property(TARGET openMVG_dataio PROPERTY FOLDER OpenMVG/OpenMVG)
target_link_libraries( openMVG_dataio openMVG_sfm ${Boost_LIBRARIES})
if(OpenMVG_USE_OPENCV)
  target_link_libraries( openMVG_dataio ${OpenCV_LIBS})
endif(OpenMVG_USE_OPENCV)
if(OpenMVG_USE_ALEMBIC)
  target_link_libraries( openMVG_dataio ${ABC_LIBRARIES})
endif(OpenMVG_USE_ALEMBIC)
install(TARGETS openMVG_dataio DESTINATION lib/ EXPORT openMVG-targets)

UNIT_TEST(openMVG AlembicImporter "openMVG_dataio")