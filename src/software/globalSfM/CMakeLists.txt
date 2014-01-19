
SET(LIBS ${EXTRA_LIBS}
					  openMVG_multiview
					  openMVG_image
					  vlsift
					  lemon
			      stlplus
            ceres)

ADD_EXECUTABLE(
  openMVG_main_GlobalSfM
  main_GlobalSfM.cpp
  SfMGlobalEngine.cpp)
TARGET_LINK_LIBRARIES(
  openMVG_main_GlobalSfM
    ${LIBS}
    openMVG_linearProgramming
    openMVG_lInftyComputerVision)

ADD_EXECUTABLE(
  openMVG_main_computeMatchesEssential
  main_computeMatchesEssential.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_computeMatchesEssential
                     ${LIBS})

