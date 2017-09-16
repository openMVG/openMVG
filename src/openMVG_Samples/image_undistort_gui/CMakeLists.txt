find_package( Qt5Widgets QUIET )

if( Qt5Widgets_FOUND )

  set(CMAKE_INCLUDE_CURRENT_DIR ON)

  set(CMAKE_AUTOMOC ON)

  set( OPENMVG_SAMPLE_IMAGE_UNDISTORT_SRCS main.cc MainWindow.cc QImageInterface.cc )
  set( OPENMVG_SAMPLE_IMAGE_UNDISTORT_HDRS MainWindow.hh QImageInterface.hh )

  add_executable( openMVG_sample_image_undistort ${OPENMVG_SAMPLE_IMAGE_UNDISTORT_SRCS} ${OPENMVG_SAMPLE_IMAGE_UNDISTORT_HDRS})
  target_link_libraries( openMVG_sample_image_undistort Qt5::Widgets openMVG_image openMVG_multiview )

   set_property(TARGET openMVG_sample_image_undistort PROPERTY CXX_STANDARD 11)
endif( Qt5Widgets_FOUND)