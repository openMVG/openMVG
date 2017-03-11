cmake_minimum_required(VERSION 2.8)

project(cudasift)
set(cudaSift_VERSION_MAJOR 2)
set(cudaSift_VERSION_MINOR 0)
set(cudaSift_VERSION_PATCH 0)

set(CPACK_PACKAGE_VERSION_MAJOR "${cudaSift_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${cudaSift_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${cudaSift_VERSION_PATCH}")
set(CPACK_GENERATOR "ZIP")
include(CPack)

find_package(CUDA)
if (NOT CUDA_FOUND)
  message(STATUS "CUDA not found. CUDA SIFT will not be built.")
endif(NOT CUDA_FOUND)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -msse2 ")
list(APPEND CUDA_NVCC_FLAGS "-arch=sm_35;-lineinfo;--compiler-options;-O2;-DVERBOSE")

set(cuda_sources
  # dynamic.cu
  cudaImage.cu
  cudaImage.h
  cudaSiftH.cu
  cudaSiftH.h
  matching.cu
  cudaSiftD.h
  cudaSift.h
  cudautils.h
)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)

set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/cudasift")

set(
    CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -O3 -gencode arch=compute_35,code=sm_35;
    )

cuda_add_library(cudasift SHARED
  cudaImage.cu
  cudaImage.h
  cudaSiftH.cu
  cudaSiftH.h
  matching.cu
  cudaSiftD.h
  cudaSift.h
  cudautils.h
)

SET(CUDA_SEPARABLE_COMPILATION ON)


set_property(TARGET cudasift PROPERTY FOLDER OpenMVG/3rdParty)

install(
  DIRECTORY .
  DESTINATION ${INCLUDE_INSTALL_DIR}
  COMPONENT headers
  FILES_MATCHING PATTERN "*.h *.hpp"
)

install(TARGETS cudasift DESTINATION lib EXPORT openMVG-targets)
