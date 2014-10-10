
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

#
# END - OpenGL
#

# - VO (WIP)
#

IF(USE_OPENCV)
  FIND_PACKAGE( OpenCV )

  IF(NOT OpenCV_FOUND)
    MESSAGE(STATUS "WARNING: OpenCV was not found. ")
  ELSE(NOT OpenCV_FOUND)
    ADD_EXECUTABLE(openMVG_main_VO main_VO.cpp)
    TARGET_LINK_LIBRARIES(openMVG_main_VO
      ${OPENGL_gl_LIBRARY}
      glfw
      ${GLFW_LIBRARIES}
      ${OpenMVG_LIBS}
      fast
      ${OpenCV_LIBS})
    SET_PROPERTY(TARGET openMVG_main_VO PROPERTY FOLDER OpenMVG/software)
    INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
  ENDIF(NOT OpenCV_FOUND)
ENDIF(USE_OPENCV)

