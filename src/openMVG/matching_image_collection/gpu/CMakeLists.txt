PROJECT(bitMatcher_cuda)

FIND_PACKAGE(CUDA REQUIRED)
FIND_PACKAGE(OpenCV 3 REQUIRED)

FILE(GLOB 
    SOURCES "*.h" ".hpp" ".cpp")

INCLUDE(FindCUDA)

link_directories(${cuda_library_dirs})

SET(CUDASRCS
    bitMatcher.cu
)

LIST(APPEND CUDA_NVCC_FLAGS "-D_MWAITXINTRIN_H_INCLUDED -lineinfo -Xptxas -v -use_fast_math -O3 -gencode arch=compute_60,code=sm_60 -std=c++11")
message("CUDA_NVCC_FLAGS: ${CUDA_NVCC_FLAGS}")

SET(CUDA_PROPAGATE_HOST_FLAGS OFF)

cuda_add_library(bitMatcher_cuda ${CUDASRCS} ${SOURCES} SHARED)
target_link_libraries(bitMatcher_cuda ${OpenCV_LIBS} ${CUDA_LIBRARIES})

install(TARGETS bitMatcher_cuda DESTINATION lib EXPORT openMVG-targets)

add_definitions(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../image/")
