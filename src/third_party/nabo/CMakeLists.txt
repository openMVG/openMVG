set( NABO_SRCS nabo.cpp brute_force_cpu.cpp kdtree_cpu.cpp )
set( NABO_HDRS nabo.h index_heap.h nabo_private.h )

add_library( libnabo STATIC ${NABO_SRCS} )
SET_PROPERTY(TARGET libnabo PROPERTY FOLDER OpenMVG/3rdParty)

install(FILES ${NABO_HDRS}
  DESTINATION ${INCLUDE_INSTALL_DIR} COMPONENT Devel
)

INSTALL(TARGETS libnabo DESTINATION lib EXPORT openMVG-targets)
