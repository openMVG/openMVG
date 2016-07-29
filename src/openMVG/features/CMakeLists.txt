
project(features)

add_subdirectory(latch)
add_subdirectory(deep)

file(GLOB features_files_headers
       *.hpp
       akaze/*.hpp
       liop/*.hpp
       mser/*.hpp
       tbmr/*.hpp)

file(GLOB features_files_sources
       *.cpp
       akaze/*.cpp
       liop/*.cpp
       mser/*.cpp
       tbmr/*.cpp)

file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM features_files_sources ${REMOVEFILESUNITTEST})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
	/usr/local/cuda/include
)

set_source_files_properties(${features_files_sources} PROPERTIES LANGUAGE CXX)
add_library(openMVG_features ${features_files_sources} ${features_files_headers})
set_target_properties(openMVG_features PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
install(TARGETS openMVG_features DESTINATION lib EXPORT openMVG-targets)
set_property(TARGET openMVG_features PROPERTY FOLDER OpenMVG/OpenMVG)

add_dependencies(openMVG_features latch_cuda)
target_link_libraries(openMVG_features latch_cuda)

add_dependencies(openMVG_features deepDescriptor_static)
target_link_libraries(openMVG_features deepDescriptor_static)

add_dependencies(openMVG_features deepDescriptorTHNets_static)
target_link_libraries(openMVG_features deepDescriptorTHNets_static)

UNIT_TEST(openMVG features "openMVG_features;stlplus")

add_subdirectory(mser)
add_subdirectory(tbmr)
