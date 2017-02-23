find_package( Qt5Widgets )

if( Qt5Widgets_FOUND )

cmake_policy( SET CMP0025 OLD )

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

set( FEATURES_PAIR_DEMO_SRCS main.cpp ImageView.cpp MainWindow.cpp )
set( FEATURES_PAIR_DEMO_HDRS ImageView.hpp MainWindow.hpp )

add_executable( openMVG_sample_describe_and_match_GUI ${FEATURES_PAIR_DEMO_SRCS} ${FEATURES_PAIR_DEMO_HDRS} )
target_link_libraries( openMVG_sample_describe_and_match_GUI Qt5::Widgets openMVG_features openMVG_image openMVG_matching vlsift )

set_target_properties( openMVG_sample_describe_and_match_GUI PROPERTIES CXX_STANDARD 11)

endif()