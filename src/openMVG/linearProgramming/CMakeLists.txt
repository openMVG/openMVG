
add_library(openMVG_linearProgramming
  bisectionLP.hpp
  linearProgrammingInterface.hpp
  linearProgrammingOSI_X.cpp
  linearProgrammingOSI_X.hpp
  linearProgramming.hpp)

set_property(TARGET openMVG_linearProgramming PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_linearProgramming DESTINATION lib EXPORT openMVG-targets)
set_target_properties(openMVG_linearProgramming PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

target_link_libraries(openMVG_linearProgramming
  PRIVATE
    ${CLP_LIBRARIES}     # clp + solver wrapper
    ${COINUTILS_LIBRARY} # container tools
    ${OSI_LIBRARY}       # generic LP
)

target_include_directories(openMVG_linearProgramming
  PRIVATE
    ${CLP_INCLUDE_DIRS}
    ${COINUTILS_INCLUDE_DIRS}
    ${OSI_INCLUDE_DIRS}
)

set(UNIT_TEST_LIBS openMVG_linearProgramming)
UNIT_TEST(openMVG linearProgramming "${UNIT_TEST_LIBS}")

add_subdirectory(lInfinityCV)
