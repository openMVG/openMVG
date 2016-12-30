add_library(openMVG_numeric
  numeric.cpp
  numeric.h)
set_target_properties(openMVG_numeric PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_numeric PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_numeric DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG numeric openMVG_numeric)
UNIT_TEST(openMVG poly openMVG_numeric)
UNIT_TEST(openMVG lm openMVG_numeric)
UNIT_TEST(openMVG l1_solver openMVG_numeric)
