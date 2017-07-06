
add_executable(openMVG_main_ColHarmonize main_ColHarmonize.cpp
  colorHarmonizeEngineGlobal.cpp)
target_link_libraries(openMVG_main_ColHarmonize
  openMVG_system
  openMVG_image
  openMVG_kvld
  openMVG_features
  openMVG_sfm
  stlplus
  ${CLP_LIBRARIES} ${COINUTILS_LIBRARY}
  ${OSI_LIBRARY})

set_property(TARGET openMVG_main_ColHarmonize PROPERTY FOLDER OpenMVG/software)

install(TARGETS openMVG_main_ColHarmonize DESTINATION bin/)
