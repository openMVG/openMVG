set( OPENMVG_GEOMETRY_IO_SRCS ply.cpp )
set( OPENMVG_GEOMETRY_IO_HDRS ply.hpp )

add_library(openMVG_geometry_io ${OPENMVG_GEOMETRY_IO_SRCS} ${OPENMVG_GEOMETRY_IO_HDRS} )
set_target_properties(openMVG_geometry_io PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
install(TARGETS openMVG_geometry_io DESTINATION lib EXPORT openMVG-targets)
set_property(TARGET openMVG_geometry_io PROPERTY FOLDER OpenMVG/OpenMVG)


UNIT_TEST(openMVG ply "openMVG_numeric;openMVG_geometry_io")

