# This file is part of OpenMVG (Open Multiple View Geometry) C++ library.

# Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

find_package( Qt5Widgets )
find_package( Qt5Svg )

if (Qt5Widgets_FOUND AND Qt5Svg_FOUND )

  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  set(CMAKE_AUTOMOC ON)

  set( VIEWER_SOURCES main.cpp mainwindow.cpp mainframe.cpp pairgraphicsitem.cpp matchingpairgraphicsview.cpp )
  set( VIEWER_HEADERS pairgraphicsitem.h mainframe.h mainwindow.h matchingpairgraphicsview.h )

  if (APPLE )
    add_executable( ui_openMVG_MatchesViewer MACOSX_BUNDLE ${VIEWER_SOURCES} ${ALTERNATIVE_VO_HEADERS} )
  else( APPLE )
    add_executable( ui_openMVG_MatchesViewer ${VIEWER_SOURCES} ${VIEWER_HEADERS} )
  endif( APPLE )

  target_include_directories( ui_openMVG_MatchesViewer PUBLIC ${OPENMVG_INCLUDE_DIRS} )
  target_link_libraries( ui_openMVG_MatchesViewer Qt5::Widgets Qt5::Svg ${OpenMVG_LIBRARIES} )

  find_package( Qt5OpenGL )
  if (Qt5OpenGL_FOUND )
    target_link_libraries( ui_openMVG_MatchesViewer Qt5::OpenGL)
  else()
    add_definitions(-DQT_NO_OPENGL)
  endif()

  set_property(TARGET ui_openMVG_MatchesViewer PROPERTY FOLDER OpenMVG/software/ui)
  install(TARGETS ui_openMVG_MatchesViewer DESTINATION bin/)
endif(Qt5Widgets_FOUND AND Qt5Svg_FOUND )
