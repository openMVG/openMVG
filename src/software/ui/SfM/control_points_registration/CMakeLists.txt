# ==============================================================================
# Generic CMake File for QTBased project
# ==============================================================================

# The name of my project.
SET(ProjectName ui_openMVG_control_points_registration)
PROJECT(${ProjectName})

# CMake version control
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# ==============================================================================
# Libraries setup
# ==============================================================================
# Find and setup Qt4.
FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

# ==============================================================================
# List all source files here.
# ==============================================================================

# All Source Files here
SET(SourceFiles
  main.cpp
  mainLayout.cpp
  GraphicsMainWindow.cpp
  ControlPointTableView.cpp
  node.cpp)

# All non moc header files here
SET(HeaderFiles
  mainLayout.hpp
  GraphicsMainWindow.hpp
  node.hpp
  document.hpp)

# ==============================================================================
# Qt4 specific moc source files generation files
# ==============================================================================
# After this call, moc source files will be generated.
QT4_AUTOMOC(${SourceFiles})

# Call this instruction if you don't include moc files.
QT4_WRAP_CPP(MocSourceFiles ${HeaderFiles})

# ==============================================================================
# Add the following files, modules and libraries to the project.
# ==============================================================================
# Files
ADD_EXECUTABLE(${ProjectName} 
  ${SourceFiles}
  ${HeaderFiles}
  ${MocSourceFiles})

# Linking with Qt Libraries.
TARGET_LINK_LIBRARIES(${ProjectName}
  ${QT_LIBRARIES}
  ${OpenMVG_LIBRARIES}
)

SET_PROPERTY(TARGET ${ProjectName} PROPERTY FOLDER OpenMVG/software/ui)
INSTALL(TARGETS ${ProjectName} DESTINATION bin/)

