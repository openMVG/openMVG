
ADD_LIBRARY(openMVG_SfM_Global SfMGlobalEngine.cpp)

ADD_EXECUTABLE(openMVG_main_GlobalSfM main_GlobalSfM.cpp)
TARGET_LINK_LIBRARIES(
  openMVG_main_GlobalSfM
  openMVG_SfM_Global
  ${OpenMVG_LIBS})

set_property(TARGET openMVG_main_GlobalSfM PROPERTY FOLDER OpenMVG/software)
set_property(TARGET openMVG_SfM_Global PROPERTY FOLDER OpenMVG/OpenMVG)

