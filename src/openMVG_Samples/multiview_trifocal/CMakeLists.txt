set(trifocal_sources
  trifocal.cpp trifocal.h  # TODO(future): make into .hpp header-only lib
  trifocal-app.cpp trifocal-app.h 
  trifocal-util.cpp trifocal-util.h

)
add_library(openMVG_trifocal ${trifocal_sources})
target_link_libraries(openMVG_trifocal
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  ${STLPLUS_LIBRARY})
UNIT_TEST(openMVG trifocal_stub "openMVG_trifocal")
 
 

add_executable(openMVG_sample_multiview_trifocal ${trifocal_sources} )
target_link_libraries(openMVG_sample_multiview_trifocal openMVG_trifocal)
target_compile_definitions(openMVG_sample_multiview_trifocal
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET openMVG_sample_multiview_trifocal PROPERTY FOLDER OpenMVG/Samples)

