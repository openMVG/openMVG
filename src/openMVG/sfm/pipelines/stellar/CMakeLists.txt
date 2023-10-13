set( # TODO: each pipelines/subfolder its own library/build module
  sfm_stellar_files_header
  stellar_definitions.hpp
  relative_scale.hpp
  stellar_solver.hpp
  #   pipelines/sequential/SfmSceneInitializerStellar.hpp
)

set(
  sfm_stellar_files_cpp
  stellar_solver.cpp
  #   pipelines/sequential/SfmSceneInitializerStellar.cpp
)

add_library(openMVG_sfm_stellar ${sfm_stellar_files_header} ${sfm_stellar_files_cpp})
target_compile_features(openMVG_sfm_stellar INTERFACE ${CXX11_FEATURES})
target_link_libraries(openMVG_sfm_stellar
  PUBLIC
    openMVG_sfm_base
    ${cereal_TARGET}
    ${OPENMVG_LIBRARY_DEPENDENCIES}
)

target_link_libraries(openMVG_sfm_stellar
  LINK_PRIVATE
)
target_include_directories(openMVG_sfm_stellar
  PRIVATE
    ${CERES_INCLUDE_DIRS}
    ${LEMON_INCLUDE_DIRS}
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)
set_target_properties(openMVG_sfm_stellar PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_sfm_stellar PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_sfm_stellar DESTINATION lib EXPORT openMVG-targets)


UNIT_TEST(openMVG stellar_SfM
  "openMVG_multiview_test_data;openMVG_sfm_stellar;${STLPLUS_LIBRARY}")

add_subdirectory(engine)

