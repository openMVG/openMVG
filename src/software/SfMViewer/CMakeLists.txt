
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

find_file(GL3HEADER gl3.h)
if(NOT ${GL3HEADER} MATCHES NOTFOUND)
  add_definitions(-DGLFW_INCLUDE_GL3)
endif()

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
)

set(SOURCES
  main.cpp
)

add_executable(openMVG_main_sfmViewer ${SOURCES} )
target_link_libraries(openMVG_main_sfmViewer
  openMVG_image
  openMVG_features
  openMVG_multiview
  openMVG_sfm
  ${STLPLUS_LIBRARY}
  ${OPENGL_gl_LIBRARY}
  glfw
  ${GLFW_LIBRARIES})

set_property(TARGET openMVG_main_sfmViewer PROPERTY FOLDER OpenMVG/software)
install(TARGETS openMVG_main_sfmViewer DESTINATION bin/)
