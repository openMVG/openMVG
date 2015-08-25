
ADD_EXECUTABLE(openMVG_main_ColHarmonize main_ColHarmonize.cpp
  colorHarmonizeEngineGlobal.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_ColHarmonize
  openMVG_system
  openMVG_image
  openMVG_kvld
  openMVG_features
  openMVG_sfm
  stlplus
  ${CLP_LIBRARIES} ${COINUTILS_LIBRARY}
  ${OSI_LIBRARY})

SET_PROPERTY(TARGET openMVG_main_ColHarmonize PROPERTY FOLDER OpenMVG/software)

INSTALL(TARGETS openMVG_main_ColHarmonize DESTINATION bin/)
