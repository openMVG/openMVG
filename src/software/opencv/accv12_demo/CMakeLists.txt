if(OpenMVG_USE_OPENCV)
    add_executable(opencv_fitting main.cpp)
    target_link_libraries(opencv_fitting opencv::opencv openMVG_multiview)
    target_include_directories(opencv_fitting PRIVATE ${opencv_INCLUDE_DIRS})
    target_compile_definitions(opencv_fitting
      PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    install(TARGETS opencv_fitting DESTINATION bin/)
endif(OpenMVG_USE_OPENCV)
