
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

find_package( OpenCV )

if(NOT OpenCV_FOUND)
  message(STATUS "WARNING: OpenCV was not found. If you want compile this sample, configure the OpenCV_DIR variable in src/software/opencv/accv12_demo/CMakelists.txt")
else(NOT OpenCV_FOUND)

add_executable( opencv_fitting main.cpp )
target_link_libraries( opencv_fitting ${OpenCV_LIBS} openMVG_multiview) 

endif(NOT OpenCV_FOUND)

