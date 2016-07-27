
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

file(GLOB_RECURSE REMOVEFILESFASTCPP fast/*.cpp)
list(REMOVE_ITEM features_files_sources ${REMOVEFILESFASTCPP})

file(GLOB_RECURSE REMOVEFILESFASTHPP fast/*.hpp)
list(REMOVE_ITEM features_files_sources ${REMOVEFILESFASTHPP})

set_source_files_properties(${features_files_sources} PROPERTIES LANGUAGE CXX)
add_library(openMVG_features ${features_files_sources} ${features_files_headers})
set_target_properties(openMVG_features PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
install(TARGETS openMVG_features DESTINATION lib EXPORT openMVG-targets)
set_property(TARGET openMVG_features PROPERTY FOLDER OpenMVG/OpenMVG)
add_dependencies(openMVG_features fast)
target_link_libraries(openMVG_features fast)

UNIT_TEST(openMVG features "openMVG_features;stlplus")

add_subdirectory(fast)
add_subdirectory(mser)
add_subdirectory(tbmr)

