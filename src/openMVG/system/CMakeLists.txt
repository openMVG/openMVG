
add_library(openMVG_system
  timer.hpp
  timer.cpp)
set_target_properties(openMVG_system PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_system PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_system DESTINATION lib/ EXPORT openMVG-targets)

UNIT_TEST(openMVG progress "")