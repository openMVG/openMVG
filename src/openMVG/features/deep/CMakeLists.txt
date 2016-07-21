CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
CMAKE_POLICY(VERSION 2.8)

PROJECT(OpenMVG_DeepDescriptors)

ADD_SUBDIRECTORY(cunnproduction)

FIND_PACKAGE(CUDA 5.5 REQUIRED)
FIND_PACKAGE(Torch REQUIRED)
FIND_PACKAGE(OpenCV REQUIRED)

ADD_DEFINITIONS(-std=c++11 -fPIC)

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
  cunnproduction
)

LINK_DIRECTORIES(${TH_INSTALL_LIB} ${THC_INSTALL_LIB})
ADD_LIBRARY(loader_static STATIC src/loader.cpp)
ADD_LIBRARY(wrapper_static STATIC src/wrapper.cpp)
ADD_LIBRARY(deepDescriptor_static STATIC src/DeepClassifier.hpp 
										src/DeepClassifierKeypoint.hpp 
										src/DeepClassifierMatch.hpp 
										src/DeepClassifierTorch.cpp 
										src/DeepClassifierTorch.hpp 
										src/loader.h src/wrapper.h)

install(TARGETS loader_static DESTINATION lib EXPORT openMVG-targets)
install(TARGETS wrapper_static DESTINATION lib EXPORT openMVG-targets)
install(TARGETS deepDescriptor_static DESTINATION lib EXPORT openMVG-targets)


TARGET_LINK_LIBRARIES(loader_static TH THC cunnproduction_static)
TARGET_LINK_LIBRARIES(wrapper_static TH THC cunnproduction_static loader_static)
TARGET_LINK_LIBRARIES(deepDescriptor_static TH THC cunnproduction_static loader_static wrapper_static ${OpenCV_LIBS})
