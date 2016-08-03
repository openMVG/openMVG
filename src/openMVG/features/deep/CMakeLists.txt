CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
CMAKE_POLICY(VERSION 2.8)

PROJECT(OpenMVG_DeepDescriptors)

ADD_SUBDIRECTORY(cunnproduction)
ADD_SUBDIRECTORY(thnets)

FIND_PACKAGE(CUDA REQUIRED)
FIND_PACKAGE(Torch REQUIRED)
FIND_PACKAGE(OpenCV REQUIRED)

ADD_DEFINITIONS(-std=c++11 -fPIC)

IF(DEFINED Torch_INSTALL_INCLUDE)
  SET(THC_INSTALL_INCLUDE ${Torch_INSTALL_INCLUDE})
  SET(TH_INSTALL_INCLUDE ${Torch_INSTALL_INCLUDE})
  SET(TH_INSTALL_LIB ${Torch_INSTALL_LIB})
  SET(THC_INSTALL_LIB ${Torch_INSTALL_LIB})
ENDIF()

ADD_LIBRARY(loader_static STATIC src/loader.cpp)
ADD_LIBRARY(wrapper_static STATIC src/wrapper.cpp)
ADD_LIBRARY(deepDescriptor_static STATIC 
	src/DeepClassifier.hpp 
	src/DeepClassifierKeypoint.hpp 
	src/DeepClassifierMatch.hpp 
	src/DeepClassifierTorch.cpp 
	src/DeepClassifierTorch.hpp 
	src/loader.h src/wrapper.h
)

TARGET_INCLUDE_DIRECTORIES (loader_static PRIVATE
  ${THC_INSTALL_INCLUDE}/THC
  ${THC_INSTALL_INCLUDE}
  ${TH_INSTALL_INCLUDE}/TH
  ${TH_INSTALL_INCLUDE}
  ${CUDA_INCLUDE_DIRS}
	cunnproduction
)

TARGET_INCLUDE_DIRECTORIES (wrapper_static PRIVATE
  ${THC_INSTALL_INCLUDE}/THC
  ${THC_INSTALL_INCLUDE}
  ${TH_INSTALL_INCLUDE}/TH
  ${TH_INSTALL_INCLUDE}
  ${CUDA_INCLUDE_DIRS}
	cunnproduction
)

TARGET_INCLUDE_DIRECTORIES (deepDescriptor_static PRIVATE
  ${THC_INSTALL_INCLUDE}/THC
  ${THC_INSTALL_INCLUDE}
  ${TH_INSTALL_INCLUDE}/TH
  ${TH_INSTALL_INCLUDE}
  ${CUDA_INCLUDE_DIRS}
	cunnproduction
)

ADD_LIBRARY(deepDescriptorTHNets_static STATIC 
	src/DeepClassifier.hpp 
	src/DeepClassifierKeypoint.hpp 
	src/DeepClassifierMatch.hpp 
	src/DeepClassifierTHNets.hpp 
	src/DeepClassifierTHNets.cpp
)

TARGET_INCLUDE_DIRECTORIES(deepDescriptorTHNets_static PRIVATE thnets)

INSTALL(TARGETS loader_static DESTINATION lib EXPORT openMVG-targets)
INSTALL(TARGETS wrapper_static DESTINATION lib EXPORT openMVG-targets)
INSTALL(TARGETS deepDescriptor_static DESTINATION lib EXPORT openMVG-targets)
INSTALL(TARGETS deepDescriptorTHNets_static DESTINATION lib EXPORT openMVG-targets)

ADD_DEPENDENCIES(deepDescriptor_static cunnproduction_static loader_static wrapper_static)
ADD_DEPENDENCIES(deepDescriptorTHNets_static THNets_CPU)

TARGET_LINK_LIBRARIES(loader_static TH THC cunnproduction_static)
TARGET_LINK_LIBRARIES(wrapper_static TH THC cunnproduction_static loader_static)
TARGET_LINK_LIBRARIES(deepDescriptor_static TH THC cunnproduction_static loader_static wrapper_static ${OpenCV_LIBS})
TARGET_LINK_LIBRARIES(deepDescriptorTHNets_static THNets_CPU jpeg png ${OpenCV_LIBS})
