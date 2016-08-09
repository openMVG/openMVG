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

set_source_files_properties(${dataio_cpp} PROPERTIES LANGUAGE CXX)
add_library(openMVG_dataio
  ${dataio_header}
  ${dataio_cpp})
set_property(TARGET openMVG_dataio PROPERTY FOLDER OpenMVG/OpenMVG)
target_link_libraries( openMVG_dataio openMVG_sfm ${Boost_LIBRARIES})
if(OpenMVG_USE_OPENCV)
  target_link_libraries( openMVG_dataio ${OpenCV_LIBS})
endif(OpenMVG_USE_OPENCV)
set_target_properties(openMVG_dataio PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

install(TARGETS openMVG_dataio DESTINATION lib/ EXPORT openMVG-targets)