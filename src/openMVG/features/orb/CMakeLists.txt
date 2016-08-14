PROJECT(orb_cuda)

FIND_PACKAGE(CUDA REQUIRED)
FIND_PACKAGE(OpenCV 3 REQUIRED)

FILE(GLOB 
	SOURCES "OrbClassifierOpenMVG.cpp" "*.c" "*.h" ".hpp")

INCLUDE(FindCUDA)

#LIST(APPEND CMAKE_CXX_FLAGS "-std=c++11 -O3 -I/opt/cuda/include -L/opt/cuda/lib64 -lcuda -lcudart -gencode arch=compute_30,code=sm_30")
#LIST(APPEND CMAKE_C_FLAGS "-std=c++11 -O3 -I/opt/cuda/include -L/opt/cuda/lib64 -lcuda -lcudart")

add_library(orb_cuda SHARED ${SOURCES})

install(TARGETS orb_cuda DESTINATION lib EXPORT openMVG-targets)

include_directories(${CUDA_INCLUDE_DIRS})
target_link_libraries(orb_cuda ${CUDA_LIBRARIES})

include_directories(${OpenCV_INCLUDE_DIRS})
target_link_libraries(orb_cuda ${OpenCV_LIBS})

add_definitions(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../image/")
