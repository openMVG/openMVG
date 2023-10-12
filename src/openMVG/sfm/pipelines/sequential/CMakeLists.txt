set(sequential_sfm_files_header
  SfmSceneInitializerStellar.hpp
  sequential_SfM_util.hpp
  sequential_SfM_base.hpp
  SfmSceneInitializerMaxPair.hpp
  SfmSceneInitializer.hpp
  sequential_SfM.hpp
  sequential_SfM2.hpp
)

set(sequential_sfm_files_cpp
  sequential_SfM.cpp
  sequential_SfM2.cpp
  sequential_SfM_base.cpp
  SfmSceneInitializerMaxPair.cpp
  SfmSceneInitializerStellar.cpp
  sequential_SfM_util.cpp
)

add_library(openMVG_sequential_sfm ${sequential_sfm_files_header} ${sequential_sfm_files_cpp})
target_compile_features(openMVG_sequential_sfm INTERFACE ${CXX11_FEATURES})

target_link_libraries(openMVG_sequential_sfm
  PUBLIC
    openMVG_geometry
    openMVG_features
    openMVG_graph
    openMVG_matching
    openMVG_multiview
    openMVG_multiview_trifocal
    ${cereal_TARGET}
    ${OPENMVG_LIBRARY_DEPENDENCIES}
)

target_include_directories(openMVG_sequential_sfm
  PRIVATE
    ${CERES_INCLUDE_DIRS}
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)


if (MSVC)
  set_target_properties(openMVG_sequential_sfm PROPERTIES COMPILE_FLAGS "/bigobj")
endif (MSVC)

set_target_properties(openMVG_sequential_sfm PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_sequential_sfm PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_sequential_sfm DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG sequential_SfM
  "openMVG_multiview_test_data;openMVG_sfm;openMVG_multiview_trifocal;${STLPLUS_LIBRARY}")

UNIT_TEST(openMVG sequential_SfM2
  "openMVG_multiview_test_data;openMVG_sfm;${STLPLUS_LIBRARY}")
