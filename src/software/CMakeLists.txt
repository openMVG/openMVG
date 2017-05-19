add_subdirectory(SfM)
add_subdirectory(Localization)
add_subdirectory(Geodesy)
#add_subdirectory(opencv)

if (OpenMVG_BUILD_OPENGL_EXAMPLES)
  add_subdirectory(SfMViewer)
endif (OpenMVG_BUILD_OPENGL_EXAMPLES)

add_subdirectory(colorHarmonize)
add_subdirectory(VO)

find_package(Qt5 COMPONENTS Core Widgets Svg REQUIRED)
if (Qt5_FOUND)
  add_subdirectory(ui)
  message("** Found Qt5, ui part will be build")
endif (Qt5_FOUND)

add_subdirectory( SfMWebGLViewer )
