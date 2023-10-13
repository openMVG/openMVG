set(sfm_sequential_files_header
  SfmSceneInitializerStellar.hpp
  sequential_SfM_util.hpp
  sequential_SfM_base.hpp
  SfmSceneInitializerMaxPair.hpp
  SfmSceneInitializer.hpp
  sequential_SfM.hpp
  sequential_SfM2.hpp
  sfm_robust_model_estimation_trifocal.hpp
)

set(sfm_sequential_files_cpp
  sequential_SfM.cpp
  sequential_SfM2.cpp
  sequential_SfM_base.cpp
  SfmSceneInitializerMaxPair.cpp
  SfmSceneInitializerStellar.cpp
  sequential_SfM_util.cpp
  # if stable, move to base or base_extras:
  sfm_robust_model_estimation_trifocal.cpp
)

add_library(openMVG_sfm_sequential ${sfm_sequential_files_header} ${sfm_sequential_files_cpp})
target_compile_features(openMVG_sfm_sequential INTERFACE ${CXX11_FEATURES})

target_link_libraries(openMVG_sfm_sequential
  PUBLIC
    openMVG_geometry
    openMVG_features
    openMVG_graph
    openMVG_matching
    openMVG_multiview
    openMVG_multiview_trifocal
    openMVG_sfm_base
    openMVG_sfm_stellar
    ${cereal_TARGET}
    ${OPENMVG_LIBRARY_DEPENDENCIES}
    ${STLPLUS_LIBRARY}
)

target_include_directories(openMVG_sfm_sequential
  PRIVATE
    ${CERES_INCLUDE_DIRS}
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)


if (MSVC)
  set_target_properties(openMVG_sfm_sequential PROPERTIES COMPILE_FLAGS "/bigobj")
endif (MSVC)

set_target_properties(openMVG_sfm_sequential PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_sfm_sequential PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_sfm_sequential DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG sequential_SfM
  "openMVG_multiview_test_data;openMVG_sfm_sequential;openMVG_multiview_trifocal;${STLPLUS_LIBRARY}")

UNIT_TEST(openMVG sequential_SfM2
  "openMVG_multiview_test_data;openMVG_sfm_sequential;${STLPLUS_LIBRARY}")
