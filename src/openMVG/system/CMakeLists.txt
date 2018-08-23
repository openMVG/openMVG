
add_library(openMVG_system
  timer.hpp
  timer.cpp)
target_include_directories(openMVG_system PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>)
target_compile_features(openMVG_system INTERFACE ${CXX11_FEATURES})
set_target_properties(openMVG_system PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_system PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_system DESTINATION lib/ EXPORT openMVG-targets)

add_library(openMVG_progress_test INTERFACE)
target_link_libraries(openMVG_progress_test INTERFACE ${OPENMVG_LIBRARY_DEPENDENCIES})
target_include_directories(openMVG_progress_test INTERFACE ${EIGEN_INCLUDE_DIRS})

UNIT_TEST(openMVG progress "openMVG_system;openMVG_progress_test;openMVG_testing")
