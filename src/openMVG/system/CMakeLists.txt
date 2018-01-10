
add_library(openMVG_system
  timer.hpp
  timer.cpp)
target_include_directories(openMVG_system PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>)
set_target_properties(openMVG_system PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_system PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_system DESTINATION lib/ EXPORT openMVG-targets)

add_library(openMVG_progress_test INTERFACE)
target_include_directories(openMVG_progress_test INTERFACE ${EIGEN_INCLUDE_DIRS})

UNIT_TEST(openMVG progress "openMVG_system;openMVG_progress_test")
