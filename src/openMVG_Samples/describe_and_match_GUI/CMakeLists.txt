
if (OpenMVG_BUILD_GUI_SOFTWARES)

  find_package( Qt5Widgets )

  if( Qt5Widgets_FOUND )

    cmake_policy( SET CMP0025 OLD )

    set(CMAKE_INCLUDE_CURRENT_DIR ON)
    set(CMAKE_AUTOMOC ON)

    set( FEATURES_PAIR_DEMO_SRCS main.cpp ImageView.cpp MainWindow.cpp )
    set( FEATURES_PAIR_DEMO_HDRS ImageView.hpp MainWindow.hpp )

    if( APPLE )
      add_executable( openMVG_sample_describe_and_match_GUI MACOSX_BUNDLE ${FEATURES_PAIR_DEMO_SRCS} ${FEATURES_PAIR_DEMO_HDRS} )
      set_target_properties(openMVG_sample_describe_and_match_GUI PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist.in)
    else( APPLE )
      add_executable( openMVG_sample_describe_and_match_GUI WIN32 ${FEATURES_PAIR_DEMO_SRCS} ${FEATURES_PAIR_DEMO_HDRS} )
    endif( APPLE )

    target_link_libraries( openMVG_sample_describe_and_match_GUI Qt5::Widgets openMVG_features openMVG_image openMVG_matching vlsift )

    set_target_properties( openMVG_sample_describe_and_match_GUI PROPERTIES CXX_STANDARD 11)

    set_property(TARGET openMVG_sample_describe_and_match_GUI PROPERTY FOLDER OpenMVG/Samples)

  endif(Qt5Widgets_FOUND)

endif(OpenMVG_BUILD_GUI_SOFTWARES)
