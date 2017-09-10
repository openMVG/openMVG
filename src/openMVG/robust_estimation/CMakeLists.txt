
UNIT_TEST(openMVG rand_sampling "")
UNIT_TEST(openMVG robust_estimator_lineKernel "")
UNIT_TEST(openMVG robust_estimator_MaxConsensus "")
UNIT_TEST(openMVG robust_estimator_Ransac "")
#UNIT_TEST(openMVG robust_estimator_LMeds "")
UNIT_TEST(openMVG robust_estimator_ACRansac "")

add_library(openMVG_robust_estimation
  gms_filter.hpp gms_filter.cpp
)
set_property(TARGET openMVG_robust_estimation PROPERTY FOLDER OpenMVG)
install(TARGETS openMVG_robust_estimation DESTINATION lib EXPORT openMVG-targets)
set_target_properties(
  openMVG_robust_estimation
  PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR}
  VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

UNIT_TEST(openMVG gms_filter "openMVG_robust_estimation")
