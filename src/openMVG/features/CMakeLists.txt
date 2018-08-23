
project(features)

file(
  GLOB_RECURSE
  features_files_headers
  *.hpp
)

file(
  GLOB_RECURSE
  features_files_sources
  *.cpp
)

file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM features_files_sources ${REMOVEFILESUNITTEST})

set_source_files_properties(${features_files_sources} PROPERTIES LANGUAGE CXX)
add_library(openMVG_features ${features_files_sources} ${features_files_headers})
target_compile_features(openMVG_features INTERFACE ${CXX11_FEATURES})
string(REPLACE " " ";" REPLACED_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
target_compile_options(openMVG_features INTERFACE ${REPLACED_CMAKE_CXX_FLAGS})
target_include_directories(openMVG_features
  PUBLIC
    $<BUILD_INTERFACE:${EIGEN_INCLUDE_DIRS}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CEREAL_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)
target_link_libraries(openMVG_features
  PRIVATE openMVG_fast ${STLPLUS_LIBRARY}
  PUBLIC ${OPENMVG_LIBRARY_DEPENDENCIES})
if (MSVC)
  set_target_properties(openMVG_features PROPERTIES COMPILE_FLAGS "/bigobj")
  target_compile_options(openMVG_features PUBLIC "-D_USE_MATH_DEFINES")
endif (MSVC)
set_target_properties(openMVG_features PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
install(TARGETS openMVG_features DESTINATION lib EXPORT openMVG-targets)
set_property(TARGET openMVG_features PROPERTY FOLDER OpenMVG/OpenMVG)

UNIT_TEST(openMVG features "openMVG_features")
UNIT_TEST(openMVG image_describer "openMVG_features;${STLPLUS_LIBRARY}")

add_subdirectory(akaze)
add_subdirectory(mser)
add_subdirectory(sift)
add_subdirectory(tbmr)
