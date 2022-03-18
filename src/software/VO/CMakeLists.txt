
add_subdirectory( AlternativeVO )

# - VO (WIP)
#
add_executable(openMVG_main_VO main_VO.cpp Monocular_VO.hpp CGlWindow.hpp)
target_link_libraries(openMVG_main_VO
  openMVG_sfm
  ${STLPLUS_LIBRARY}
  openMVG_image
)
set_property(TARGET openMVG_main_VO PROPERTY FOLDER OpenMVG/software)

if (OpenMVG_BUILD_OPENGL_EXAMPLES)
  # deal with platform specifics
  if(OpenMVG_USE_OPENCV)
    target_link_libraries(openMVG_main_VO
        ${OpenCV_LIBS})
    target_include_directories(openMVG_main_VO PRIVATE ${OpenCV_INCLUDE_DIRS})
    target_compile_definitions(openMVG_main_VO PRIVATE HAVE_OPENCV)
  endif(OpenMVG_USE_OPENCV)

  #
  # OpenGL
  #
  find_package(OpenGL QUIET)
  if(${OPENGL_FOUND})

    target_compile_definitions(openMVG_main_VO PRIVATE HAVE_OPENGL)
    target_include_directories(openMVG_main_VO PRIVATE ${OPENGL_INCLUDE_DIR})
    target_compile_definitions(openMVG_main_VO PRIVATE GL_GLEXT_PROTOTYPES GL3_PROTOTYPES)

    if(APPLE)
      target_compile_definitions(openMVG_main_VO PRIVATE GLFW_INCLUDE_GLCOREARB)
    else()
      target_compile_definitions(openMVG_main_VO PRIVATE DGLFW_NO_GLU)
    endif()

    find_file(GL3HEADER gl3.h)
    if(NOT ${GL3HEADER} MATCHES NOTFOUND)
      target_compile_definitions(openMVG_main_VO PRIVATE GLFW_INCLUDE_GL3)
    endif()

    target_link_libraries(openMVG_main_VO
      ${OPENGL_gl_LIBRARY}
      glfw
      ${GLFW_LIBRARIES})
  endif()
endif()
