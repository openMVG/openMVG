
add_library(openMVG_kvld kvld.cpp kvld.h algorithm.cpp algorithm.h)
set_target_properties(openMVG_kvld PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_kvld PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_kvld DESTINATION lib EXPORT openMVG-targets)
