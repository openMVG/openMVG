
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

add_library(openMVG_system
  ${sytem_files_header}
  ${system_files_cpp})
set_target_properties(openMVG_system PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_system PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_system DESTINATION lib/ EXPORT openMVG-targets)

