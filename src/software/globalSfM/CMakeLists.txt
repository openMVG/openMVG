
SET(LIBS ${EXTRA_LIBS}
					  openMVG_multiview
					  openMVG_image
					  vlsift
					  lemon
			      stlplus
            ceres)

ADD_LIBRARY(openMVG_SfM_Global SfMGlobalEngine.cpp)

ADD_EXECUTABLE(openMVG_main_GlobalSfM main_GlobalSfM.cpp)
TARGET_LINK_LIBRARIES(
  openMVG_main_GlobalSfM
  openMVG_SfM_Global
  ${LIBS}
  openMVG_linearProgramming
  openMVG_lInftyComputerVision)
