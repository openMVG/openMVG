set( OPENMVG_SOFTWARE_SFMWEBGLVIEWER_SRC_DIR ${CMAKE_CURRENT_LIST_DIR} )
set( OPENMVG_SOFTWARE_SFMWEBGLVIEWER_INSTALL_SHARED_DIR ${CMAKE_INSTALL_PREFIX}/share/openMVG/webgl )
configure_file( config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h )

add_executable(openMVG_main_openMVG2WebGL openMVG_main_openMVG2WebGL.cpp)
target_link_libraries(openMVG_main_openMVG2WebGL
  openMVG_system
  openMVG_features
  openMVG_image
  openMVG_sfm
  ${STLPLUS_LIBRARY}
  )

target_include_directories( openMVG_main_openMVG2WebGL PRIVATE ${CMAKE_CURRENT_BINARY_DIR} )

install(TARGETS openMVG_main_openMVG2WebGL DESTINATION bin/)
set_property(TARGET openMVG_main_openMVG2WebGL PROPERTY FOLDER OpenMVG/software)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/camera.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/camera_gizmo.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/common.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/dual_quaternion.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/index.html"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/main.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/matrix3.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/matrix4.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/plane.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/point_cloud.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/quaternion.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/render_context.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/shader.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/style.css"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/trackball.js"  DESTINATION share/openMVG/webgl/common/ )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/common/vector.js"  DESTINATION share/openMVG/webgl/common/ )
