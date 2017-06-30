
project(matching)

file(
  GLOB
  matching_files_header
  *.hpp
)
file(
  GLOB
  matching_files_cpp
  *.cpp
)
file(GLOB REMOVEFILESUNITTEST *_test.cpp)

#Remove the future main files
list(REMOVE_ITEM matching_files_cpp ${REMOVEFILESUNITTEST})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set_source_files_properties(${matching_files_cpp} PROPERTIES LANGUAGE CXX)
add_library(openMVG_matching
  ${matching_files_header}
  ${matching_files_cpp})
target_link_libraries(openMVG_matching PRIVATE openMVG_features stlplus)
target_link_libraries(openMVG_matching PUBLIC Threads::Threads)
set_target_properties(openMVG_matching PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_matching PROPERTY FOLDER OpenMVG/OpenMVG)
if (UNIX AND AVX2_FOUND)
  target_compile_options(openMVG_matching PUBLIC "-mavx2")
endif(UNIX AND AVX2_FOUND)
install(TARGETS openMVG_matching DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG matching "openMVG_matching")
UNIT_TEST(openMVG matching_filters "openMVG_matching")
UNIT_TEST(openMVG indMatch "openMVG_matching")
UNIT_TEST(openMVG metric "openMVG_matching")

add_subdirectory(kvld)
