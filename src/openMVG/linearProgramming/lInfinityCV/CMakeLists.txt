
file(
  GLOB
  lInftycomputervision_headers
  *.hpp
)

file(
  GLOB
  lInftycomputervision_cpp
  *.cpp
)

file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)

# Remove the future main files
list(REMOVE_ITEM lInftycomputervision_cpp ${REMOVEFILESUNITTEST})

add_library(openMVG_lInftyComputerVision ${lInftycomputervision_headers} ${lInftycomputervision_cpp})
target_compile_features(openMVG_lInftyComputerVision INTERFACE ${CXX11_FEATURES})
target_link_libraries(openMVG_lInftyComputerVision PUBLIC openMVG_linearProgramming openMVG_multiview)
set_target_properties(openMVG_lInftyComputerVision PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_lInftyComputerVision PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_lInftyComputerVision DESTINATION lib EXPORT openMVG-targets)

set(LOCAL_LIBS
  openMVG_lInftyComputerVision
  openMVG_multiview_test_data
  openMVG_multiview)

# Linfinity Triangulation
UNIT_TEST(openMVG_lInftyCV triangulation "${LOCAL_LIBS}")

# Linfinity Pose/Resection estimation
UNIT_TEST(openMVG_lInftyCV resection "${LOCAL_LIBS}")
UNIT_TEST(openMVG_lInftyCV resection_robust "${LOCAL_LIBS}")

# Linfinity Global translations from triplets of heading directions
UNIT_TEST(openMVG_lInftyCV global_translations_fromTriplets "${LOCAL_LIBS}")
# Linfinity Global translations from heading directions
UNIT_TEST(openMVG_lInftyCV global_translations_fromTij "${LOCAL_LIBS}")

# Linfinity Translation and Structure
UNIT_TEST(openMVG_lInftyCV tijsAndXis "${LOCAL_LIBS}")
# Linfinity Translation and Structure with noise handling
UNIT_TEST(openMVG_lInftyCV tijsAndXis_outlier "${LOCAL_LIBS}")
