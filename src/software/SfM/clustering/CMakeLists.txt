add_subdirectory(domSetLibrary)

# compute clusters for a given sfm_data.bin
add_executable(openMVG_main_ComputeClusters main_ComputeClusters.cpp)
target_link_libraries(openMVG_main_ComputeClusters 
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_sfm
  domset)

set_property(TARGET openMVG_main_ComputeClusters PROPERTY FOLDER OpenMVG/software/clustering)
install(TARGETS openMVG_main_ComputeClusters DESTINATION bin/)

