CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
CMAKE_POLICY(VERSION 2.8)

PROJECT(OpenMVG_DeepDescriptors)

ADD_SUBDIRECTORY(thnets)

FIND_PACKAGE(CUDA REQUIRED)
# FIND_PACKAGE(Torch REQUIRED)
FIND_PACKAGE(OpenCV REQUIRED)

ADD_DEFINITIONS(-std=c++11 -fPIC)

link_directories(${cuda_library_dirs})

ADD_LIBRARY(deepDescriptorTHNets_static STATIC 
	src/DeepClassifier.hpp 
	src/DeepClassifierKeypoint.hpp 
	src/DeepClassifierMatch.hpp 
	src/DeepClassifierTHNets.hpp 
	src/DeepClassifierTHNets.cpp
)

TARGET_INCLUDE_DIRECTORIES(deepDescriptorTHNets_static PRIVATE thnets)

INSTALL(TARGETS deepDescriptorTHNets_static DESTINATION lib EXPORT openMVG-targets)

ADD_DEPENDENCIES(deepDescriptorTHNets_static THNets_CPU)

TARGET_LINK_LIBRARIES(deepDescriptorTHNets_static THNets_CPU jpeg png ${OpenCV_LIBS})
