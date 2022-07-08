set(trifocal_app_sources
  trifocal_app.cpp trifocal_app.hpp
) 
set(trifocal_sample_sources
  trifocal_sample.cpp 
) 

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# always build minus with specially tweaked optimizations.
# these were copied directly from minus original flags
set_target_properties(minus PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(minus PROPERTIES COMPILE_FLAGS "-O3 -march=native -ffast-math")
set_target_properties(minus PROPERTIES COMPILE_DEFINITIONS "NDEBUG")
add_library(openMVG_trifocal_app ${trifocal_app_sources})

target_link_libraries(openMVG_trifocal_app
  openMVG_image
  openMVG_multiview
  openMVG_multiview_trifocal
  openMVG_features
  openMVG_matching
  ${STLPLUS_LIBRARY})


# set(CMAKE_CXX_FLAGS "-g") # XXX temporary for build cycle
# set(CMAKE_C_FLAGS "-g") # XXX temporary for build cycle
# These work better: they add -g but they won't override previous flags
#set_target_properties(openMVG_test_trifocal PROPERTIES COMPILE_FLAGS "-fsanitize=address -fno-omit-frame-pointer")
set_target_properties(openMVG_trifocal_app PROPERTIES COMPILE_FLAGS "-g")

add_executable(openMVG_sample_trifocal ${trifocal_sample_sources})
target_link_libraries(openMVG_sample_trifocal openMVG_trifocal_app)
target_compile_definitions(openMVG_sample_trifocal
  PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET openMVG_sample_trifocal PROPERTY FOLDER OpenMVG/Samples)

