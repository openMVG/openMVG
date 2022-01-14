set(trifocal_sources
  trifocal.cpp trifocal.h  # TODO(future): make into .hpp header-only lib
  trifocal-util.cpp trifocal-util.h
)
set(trifocal_app_sources
  trifocal-app.cpp trifocal-app.h
) 
set(trifocal_sample
  trifocal-sample.cpp 
) 
add_library(openMVG_trifocal ${trifocal_sources})
add_library(openMVG_trifocal_app ${trifocal_app_sources})
target_link_libraries(openMVG_trifocal
  openMVG_image
  openMVG_multiview
  openMVG_features
  openMVG_matching
  ${STLPLUS_LIBRARY})
UNIT_TEST(openMVG trifocal_stub "openMVG_trifocal")
UNIT_TEST(openMVG trifocal "openMVG_trifocal_app")

target_link_libraries(openMVG_trifocal_app openMVG_trifocal)
add_executable(openMVG_sample_multiview_trifocal ${trifocal_sample})
target_link_libraries(openMVG_sample_multiview_trifocal openMVG_trifocal_app)
target_compile_definitions(openMVG_sample_multiview_trifocal
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET openMVG_sample_multiview_trifocal PROPERTY FOLDER OpenMVG/Samples)

