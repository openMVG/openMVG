
if (OpenMVG_USE_CUDA)
  find_package(CUDA REQUIRED)

  set(SOURCES CudaBruteForceMatcher.cu LatchBitMatcher.cpp)
  set(HEADERS CudaBruteForceMatcher.h LatchBitMatcher.hpp)

  list(APPEND CUDA_NVCC_FLAGS "-D_MWAITXINTRIN_H_INCLUDED -lineinfo -Xptxas -v -Xcompiler -use_fast_math -O3 -gencode arch=compute_60,code=sm_60 -std=c++11")

  set(CUDA_PROPAGATE_HOST_FLAGS ON)

  cuda_add_library(openMVG_gpu_matching_image_collection STATIC ${SOURCES} ${HEADERS})
  target_link_libraries(openMVG_gpu_matching_image_collection ${CUDA_LIBRARIES})
  target_include_directories(openMVG_gpu_matching_image_collection PUBLIC ${CUDA_INCLUDE_DIRS})
  set_property(TARGET openMVG_gpu_matching_image_collection PROPERTY FOLDER OpenMVG/OpenMVG)

  install(TARGETS openMVG_gpu_matching_image_collection DESTINATION lib EXPORT openMVG-targets)
endif()
