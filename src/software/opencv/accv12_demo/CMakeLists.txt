if(OpenMVG_USE_OPENCV)
    add_definitions(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    add_executable( opencv_fitting main.cpp )
    target_link_libraries( opencv_fitting ${OpenCV_LIBS} openMVG_multiview) 
    install(TARGETS opencv_fitting DESTINATION bin/)
endif(OpenMVG_USE_OPENCV)


