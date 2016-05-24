PROJECT(bitMatcher_cuda)

FIND_PACKAGE(CUDA REQUIRED)

FILE(GLOB 
    SOURCES "*.h" ".hpp" ".cpp")

INCLUDE(FindCUDA)

SET(CUDASRCS
    bitMatcher.cu
)

LIST(APPEND CUDA_NVCC_FLAGS "-D_MWAITXINTRIN_H_INCLUDED -lineinfo -Xptxas -v -use_fast_math -O3 -gencode arch=compute_30,code=sm_30 --default-stream per-thread")
message("CUDA_NVCC_FLAGS: ${CUDA_NVCC_FLAGS}")

SET(CUDA_PROPAGATE_HOST_FLAGS OFF)

cuda_add_library(bitMatcher_cuda ${CUDASRCS} ${SOURCES} SHARED)

install(TARGETS bitMatcher_cuda DESTINATION lib EXPORT openMVG-targets)

add_definitions(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../image/")
