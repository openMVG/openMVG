add_subdirectory(SfM)
add_subdirectory(Localization)
add_subdirectory(Geodesy)
#add_subdirectory(opencv)

if (OpenMVG_BUILD_OPENGL_EXAMPLES)
  add_subdirectory(SfMViewer)
endif (OpenMVG_BUILD_OPENGL_EXAMPLES)

add_subdirectory(colorHarmonize)
add_subdirectory(VO)

find_package(Qt4 QUIET)
if (QT_FOUND)
  add_subdirectory(ui)
endif (QT_FOUND)
