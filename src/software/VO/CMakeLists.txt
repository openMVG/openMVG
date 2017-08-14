
add_subdirectory( AlternativeVO )

if (OpenMVG_BUILD_OPENGL_EXAMPLES)

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

  #
  # END - OpenGL
  #

  # - VO (WIP)
  #

  add_executable(openMVG_main_VO main_VO.cpp Monocular_VO.hpp CGlWindow.hpp)
  target_link_libraries(openMVG_main_VO
    ${OPENGL_gl_LIBRARY}
    glfw
    ${GLFW_LIBRARIES}
    ${OpenMVG_LIBRARIES}
  )

  if(OpenMVG_USE_OPENCV)
    target_link_libraries(openMVG_main_VO
        ${OpenCV_LIBS})
    target_include_directories(openMVG_main_VO PRIVATE ${OpenCV_INCLUDE_DIRS})
    target_compile_definitions(openMVG_main_VO PRIVATE HAVE_OPENCV)
  endif(OpenMVG_USE_OPENCV)

  set_property(TARGET openMVG_main_VO PROPERTY FOLDER OpenMVG/software)

endif()
