# ==============================================================================
# Generic CMake File for QTBased project
# ==============================================================================

find_package( Qt5Widgets )

if(Qt5Widgets_FOUND)

  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  set(CMAKE_AUTOMOC ON)
  # The name of my project.
  SET(ProjectName ui_openMVG_control_points_registration)
  PROJECT(${ProjectName})

  # All Source Files here
  SET(SourceFiles
  main.cpp
  MainWindow.cpp
  GraphicsView.cpp
  ControlPointTableView.cpp
  ControlPoint2DNode.cpp)

  # All non moc header files here
  SET(HeaderFiles
  MainWindow.hpp
  GraphicsView.hpp
  ControlPoint2DNode.hpp
  document.hpp
  ControlPointTableView.hpp)

  # ==============================================================================
  # Add the following files, modules and libraries to the project.
  # ==============================================================================
  # Files
  ADD_EXECUTABLE(${ProjectName}
  ${SourceFiles}
  ${HeaderFiles})

  # Linking with Qt Libraries.
  TARGET_LINK_LIBRARIES(${ProjectName}
  Qt5::Widgets Qt5::Gui
  openMVG_sfm
  ${STLPLUS_LIBRARY})

  SET_PROPERTY(TARGET ${ProjectName} PROPERTY FOLDER OpenMVG/software/ui)
  INSTALL(TARGETS ${ProjectName} DESTINATION bin/)
endif(Qt5Widgets_FOUND)
