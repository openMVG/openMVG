file(
  GLOB_RECURSE
  spherical_files_header
  *.hpp
)
file(
  GLOB_RECURSE
  spherical_files_cpp
  *.cpp
)

file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
list(REMOVE_ITEM spherical_files_cpp ${REMOVEFILESUNITTEST})

add_library(openMVG_spherical 
    ${spherical_files_header}
    ${spherical_files_cpp})
target_compile_features(openMVG_spherical INTERFACE ${CXX11_FEATURES})
target_include_directories(openMVG_spherical
  PUBLIC
    $<BUILD_INTERFACE:${EIGEN_INCLUDE_DIRS}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
)

target_link_libraries(openMVG_spherical INTERFACE openMVG_image)
install(TARGETS openMVG_spherical DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG tangent_images "openMVG_testing;openMVG_image;openMVG_spherical;openMVG_features;vlsift")
UNIT_TEST(openMVG convert "openMVG_testing;openMVG_image;openMVG_camera")
