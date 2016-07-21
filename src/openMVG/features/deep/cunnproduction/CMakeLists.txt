CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
CMAKE_POLICY(VERSION 2.8)

PROJECT(cunnproduction)

FIND_PACKAGE(Torch)
FIND_PACKAGE(CUDA 5.5 REQUIRED) # required for std::shared_ptr

add_definitions(-std=c++11 -fPIC)
LIST(APPEND CUDA_NVCC_FLAGS "-arch=sm_60 -Xcompiler -fPIC")

IF(DEFINED Torch_INSTALL_INCLUDE)
  SET(THC_INSTALL_INCLUDE ${Torch_INSTALL_INCLUDE})
  SET(TH_INSTALL_INCLUDE ${Torch_INSTALL_INCLUDE})
  SET(TH_INSTALL_LIB ${Torch_INSTALL_LIB})
  SET(THC_INSTALL_LIB ${Torch_INSTALL_LIB})
ENDIF()

INCLUDE_DIRECTORIES(
  ${THC_INSTALL_INCLUDE}/THC
  ${THC_INSTALL_INCLUDE}
  ${TH_INSTALL_INCLUDE}/TH
  ${TH_INSTALL_INCLUDE}
  ${CUDA_INCLUDE_DIRS}
)

LINK_DIRECTORIES(${TH_INSTALL_LIB} ${THC_INSTALL_LIB})

SET(src-cuda SpatialMaxPooling.cu SpatialAveragePooling.cu SpatialConvolution.cu Linear.cu ReLU.cu)
LIST(APPEND CUDA_NVCC_FLAGS "-D_MWAITXINTRIN_H_INCLUDED -lineinfo -Xptxas -v -Xcompiler -fopenmp -use_fast_math -O3 -gencode arch=compute_60,code=sm_60 --default-stream per-thread")

SET(CUDA_PROPAGATE_HOST_FLAGS OFF)

CUDA_ADD_LIBRARY(cunn STATIC ${src-cuda})

ADD_LIBRARY(cunnproduction SHARED cunn.cpp)
ADD_LIBRARY(cunnproduction_static STATIC cunn.cpp)

install(TARGETS cunnproduction_static DESTINATION lib EXPORT openMVG-targets)
install(TARGETS cunn DESTINATION lib EXPORT openMVG-targets)

TARGET_LINK_LIBRARIES(cunnproduction TH THC cunn)
TARGET_LINK_LIBRARIES(cunnproduction_static TH THC cunn)
