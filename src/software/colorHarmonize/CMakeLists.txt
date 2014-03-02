
SET(LIBS ${EXTRA_LIBS}
			      openMVG_multiview
			      openMVG_image
			      vlsift
			      lemon
			      stlplus
            ${OPENMVG_LIBRARY_DEPENDENCIES}
            openMVG_linearProgramming)

#Add executable that computes:

ADD_EXECUTABLE(openMVG_main_ColHarmonize main_ColHarmonize.cpp
  colorHarmonizeEngineGlobal.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_ColHarmonize
  ${LIBS}
  openMVG_kvld)

