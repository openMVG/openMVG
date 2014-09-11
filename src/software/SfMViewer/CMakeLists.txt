
#
# OpenGL
#
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})
add_definitions(-DGL_GLEXT_PROTOTYPES -DGL3_PROTOTYPES)

if(APPLE)
	add_definitions(-DGLFW_INCLUDE_GLCOREARB)
endif()

#
# GLFW
#
if(!APPLE)
add_definitions(-DGLFW_NO_GLU)
endif()

FIND_FILE(GL3HEADER gl3.h)
if(NOT ${GL3HEADER} MATCHES NOTFOUND)
add_definitions(-DGLFW_INCLUDE_GL3)
endif()

INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_BINARY_DIR}
)

SET(SOURCES
  main.cpp
)

ADD_EXECUTABLE(openMVG_main_sfmViewer ${SOURCES} )
TARGET_LINK_LIBRARIES(openMVG_main_sfmViewer
     ${OPENGL_gl_LIBRARY} glfw ${GLFW_LIBRARIES}
  stlplus
  openMVG_numeric
  openMVG_multiview
  openMVG_image
)
SET_PROPERTY(TARGET openMVG_main_sfmViewer PROPERTY FOLDER OpenMVG/software)
