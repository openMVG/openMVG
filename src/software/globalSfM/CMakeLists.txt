
SET(LIBS ${EXTRA_LIBS}
					  openMVG_multiview
					  openMVG_image
					  vlsift
					  lemon
			      stlplus
            ceres)

ADD_LIBRARY(openMVG_GlobalSfM SfMGlobalEngine.cpp)

ADD_EXECUTABLE(openMVG_main_GlobalSfM main_GlobalSfM.cpp)
TARGET_LINK_LIBRARIES(
  openMVG_main_GlobalSfM
  openMVG_GlobalSfM
  ${LIBS}
  openMVG_linearProgramming
  openMVG_lInftyComputerVision)
