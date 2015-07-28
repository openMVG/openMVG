
file(
  GLOB_RECURSE
  system_files_header
  *.hpp
)
file(
  GLOB_RECURSE
  system_files_cpp
  *.cpp
)

ADD_LIBRARY(openMVG_system
  ${sytem_files_header}
  ${system_files_cpp})
SET_TARGET_PROPERTIES(openMVG_system PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
SET_PROPERTY(TARGET openMVG_system PROPERTY FOLDER OpenMVG/OpenMVG)
INSTALL(TARGETS openMVG_system DESTINATION lib/ EXPORT openMVG-targets)

