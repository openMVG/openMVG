
UNIT_TEST(openMVG rand_sampling "openMVG_testing")
UNIT_TEST(openMVG robust_estimator_lineKernel "openMVG_testing")
UNIT_TEST(openMVG robust_estimator_MaxConsensus "openMVG_testing")
UNIT_TEST(openMVG robust_estimator_Ransac "openMVG_testing")
#UNIT_TEST(openMVG robust_estimator_LMeds "openMVG_testing")
UNIT_TEST(openMVG robust_estimator_ACRansac "openMVG_testing")

add_library(openMVG_robust_estimation
  gms_filter.hpp gms_filter.cpp
)
target_compile_features(openMVG_robust_estimation INTERFACE ${CXX11_FEATURES})
target_include_directories(openMVG_robust_estimation
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)
target_link_libraries(openMVG_robust_estimation
  PUBLIC
    openMVG_numeric
)
set_property(TARGET openMVG_robust_estimation PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_robust_estimation DESTINATION lib EXPORT openMVG-targets)
set_target_properties(
  openMVG_robust_estimation
  PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR}
  VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

UNIT_TEST(openMVG gms_filter "openMVG_robust_estimation")
