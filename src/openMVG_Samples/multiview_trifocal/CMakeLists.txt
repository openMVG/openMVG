set(trifocal_app_sources
  trifocal_app.cpp trifocal_app.hpp
  trifocal_util.cpp trifocal_util.hpp
)
set(trifocal_sample_sources
  trifocal_sample.cpp
)

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${MINUS_INCLUDE_DIRS})

# ------------------------------------------------------------------------------

add_library(openMVG_trifocal_app ${trifocal_app_sources})

target_link_libraries(openMVG_trifocal_app
  openMVG_image
  openMVG_multiview
  openMVG_multiview_trifocal
  openMVG_features
  openMVG_matching
  ${STLPLUS_LIBRARY})

UNIT_TEST(openMVG trifocal_stub "openMVG_multiview_trifocal")
UNIT_TEST(openMVG trifocal "openMVG_trifocal_app")

#set_target_properties(openMVG_test_trifocal PROPERTIES COMPILE_FLAGS "-fsanitize=address -fno-omit-frame-pointer")
set_target_properties(openMVG_trifocal_app PROPERTIES COMPILE_FLAGS "-g")

add_executable(openMVG_sample_trifocal ${trifocal_sample_sources})
target_link_libraries(openMVG_sample_trifocal openMVG_trifocal_app)
target_compile_definitions(openMVG_sample_trifocal
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET openMVG_sample_trifocal PROPERTY FOLDER OpenMVG/Samples)
