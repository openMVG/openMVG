
ADD_EXECUTABLE(openMVG_main_GlobalSfM main_GlobalSfM.cpp)
TARGET_LINK_LIBRARIES(
  openMVG_main_GlobalSfM
  ${OpenMVG_LIBS})

SET_PROPERTY(TARGET openMVG_main_GlobalSfM PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_GlobalSfM DESTINATION bin/)
