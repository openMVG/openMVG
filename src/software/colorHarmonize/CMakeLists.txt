
ADD_EXECUTABLE(openMVG_main_ColHarmonize main_ColHarmonize.cpp
  colorHarmonizeEngineGlobal.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_ColHarmonize
  ${OpenMVG_LIBS}
  openMVG_kvld)

SET_PROPERTY(TARGET openMVG_main_ColHarmonize PROPERTY FOLDER OpenMVG/software)

INSTALL(TARGETS openMVG_main_ColHarmonize DESTINATION bin/)
