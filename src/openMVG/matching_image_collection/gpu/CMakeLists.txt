project(openMVG_gpu_matching_image_collection)

if (OpenMVG_USE_CUDA)
  file(GLOB 
      SOURCES "*.h" ".hpp" ".cpp")

  include(FindCUDA)

  link_directories(${cuda_library_dirs})

  set(CUDASRCS
    CudaBruteForceMatcher.cu
  )

  list(APPEND CUDA_NVCC_FLAGS "-D_MWAITXINTRIN_H_INCLUDED -lineinfo -Xptxas -v -use_fast_math -O3 -gencode arch=compute_35,code=sm_35 -std=c++11")
  message("CUDA_NVCC_FLAGS: ${CUDA_NVCC_FLAGS}")

  set(CUDA_PROPAGATE_HOST_FLAGS OFF)

  cuda_add_library(openMVG_gpu_matching_image_collection ${CUDASRCS} ${SOURCES} STATIC)
  target_link_libraries(openMVG_gpu_matching_image_collection ${CUDA_LIBRARIES})

  install(TARGETS openMVG_gpu_matching_image_collection DESTINATION lib EXPORT openMVG-targets)
endif()
