ADD_LIBRARY(openMVG_numeric
  numeric.cpp
  numeric.h)
SET_TARGET_PROPERTIES(openMVG_numeric PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
SET_PROPERTY(TARGET openMVG_numeric PROPERTY FOLDER OpenMVG/OpenMVG)
INSTALL(TARGETS openMVG_numeric DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG numeric openMVG_numeric)
UNIT_TEST(openMVG poly openMVG_numeric)
UNIT_TEST(openMVG lm openMVG_numeric)
UNIT_TEST(openMVG l1_solver openMVG_numeric)
