
if ( OpenMVG_BUILD_GUI_SOFTWARES )

  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  set(CMAKE_AUTOMOC ON)

  find_package( Qt5Widgets )

  if ( Qt5Widgets_FOUND )

    set( ALTERNATIVE_VO_SOURCES main.cpp ControlButtonsPanel.cpp MainWindow.cpp VOFolderProcessor.cpp VOViewerPanel.cpp VOViewerDrawableElements.cpp )
    set( ALTERNATIVE_VO_HEADERS ControlButtonsPanel.hpp MainWindow.hpp VOFolderProcessor.hpp VOViewerPanel.hpp VOViewerDrawableElements.hpp )

    if ( APPLE )
      add_executable( openMVG_main_AlternativeVO MACOSX_BUNDLE ${ALTERNATIVE_VO_SOURCES} ${ALTERNATIVE_VO_HEADERS} )
    else( APPLE )
      if ( WIN32 )
        add_executable( openMVG_main_AlternativeVO WIN32 ${ALTERNATIVE_VO_SOURCES} ${ALTERNATIVE_VO_HEADERS}  )
      else ( WIN32 )
        add_executable( openMVG_main_AlternativeVO ${ALTERNATIVE_VO_SOURCES} ${ALTERNATIVE_VO_HEADERS}  )
      endif( WIN32 )
    endif( APPLE )

    target_link_libraries( openMVG_main_AlternativeVO Qt5::Widgets openMVG_image openMVG_features ${STLPLUS_LIBRARY} )

    if ( OpenMVG_USE_OPENCV )
      target_link_libraries( openMVG_main_AlternativeVO ${OpenCV_LIBS} )
      target_include_directories( openMVG_main_AlternativeVO PRIVATE ${OpenCV_INCLUDE_DIRS} )
      target_compile_definitions( openMVG_main_AlternativeVO PRIVATE HAVE_OPENCV )
    endif( OpenMVG_USE_OPENCV )

    set_property( TARGET openMVG_main_AlternativeVO PROPERTY CXX_STANDARD 11 )
    set_property(TARGET openMVG_main_AlternativeVO PROPERTY FOLDER OpenMVG/software)
    install(TARGETS openMVG_main_AlternativeVO DESTINATION bin/)

  endif( Qt5Widgets_FOUND )

endif(OpenMVG_BUILD_GUI_SOFTWARES)
