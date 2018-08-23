
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
target_compile_features(openMVG_matching INTERFACE ${CXX11_FEATURES})
target_include_directories(openMVG_matching
  PUBLIC
    $<BUILD_INTERFACE:${FLANN_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)
target_link_libraries(openMVG_matching
  PRIVATE
    ${STLPLUS_LIBRARY}
    ${OPENMVG_LIBRARY_DEPENDENCIES}
  PUBLIC
    openMVG_features
    Threads::Threads
)
if (NOT DEFINED OpenMVG_USE_INTERNAL_FLANN)
target_link_libraries(openMVG_matching
  PUBLIC
    ${FLANN_LIBRARIES}
)
endif()
set_target_properties(openMVG_matching PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_matching PROPERTY FOLDER OpenMVG/OpenMVG)

if (USE_AVX2)
  target_compile_options(openMVG_matching PUBLIC "-DOPENMVG_USE_AVX2")
  if (UNIX)
    target_compile_options(openMVG_matching PUBLIC "-mavx2")
  endif (UNIX)
endif (USE_AVX2)

install(TARGETS openMVG_matching DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG matching "openMVG_matching")
UNIT_TEST(openMVG matching_filters "openMVG_matching")
UNIT_TEST(openMVG indMatch "openMVG_matching")
UNIT_TEST(openMVG metric "openMVG_matching")

add_subdirectory(kvld)
