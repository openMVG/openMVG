add_subdirectory(SfM)
add_subdirectory(Localization)
add_subdirectory(Geodesy)
#add_subdirectory(opencv)

if (OpenMVG_BUILD_OPENGL_EXAMPLES)
  add_subdirectory(SfMViewer)
endif (OpenMVG_BUILD_OPENGL_EXAMPLES)

add_subdirectory(colorHarmonize)
add_subdirectory(VO)

if (OpenMVG_BUILD_GUI_SOFTWARES)
  find_package(Qt5 COMPONENTS Core Widgets Svg)
  if (Qt5_FOUND)
    message("** Found Qt5, GUI modules will be build")
    add_subdirectory(ui)
  else (Qt5_FOUND)
    message("** Qt5 not found, GUI modules will not be build")
  endif (Qt5_FOUND)
endif (OpenMVG_BUILD_GUI_SOFTWARES)

add_subdirectory( SfMWebGLViewer )
