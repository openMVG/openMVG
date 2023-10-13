set( # TODO: each pipelines/subfolder its own library/build module
  sfm_stellar_engine_files_header
  sfm_stellar_engine.hpp
)

set(
  sfm_stellar_engine_files_cpp
  sfm_stellar_engine.cpp
)

add_library(openMVG_sfm_stellar_engine ${sfm_stellar_engine_files_header} ${sfm_stellar_engine_files_cpp})
target_compile_features(openMVG_sfm_stellar_engine INTERFACE ${CXX11_FEATURES})
target_link_libraries(openMVG_sfm_stellar_engine
  PUBLIC
    openMVG_sfm
    ${cereal_TARGET}
    ${OPENMVG_LIBRARY_DEPENDENCIES}
)

target_link_libraries(openMVG_sfm_stellar_engine
  LINK_PRIVATE
)
target_include_directories(openMVG_sfm_stellar_engine
  PRIVATE
    ${CERES_INCLUDE_DIRS}
    ${LEMON_INCLUDE_DIRS}
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)
set_target_properties(openMVG_sfm_stellar_engine PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_sfm_stellar_engine PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_sfm_stellar_engine DESTINATION lib EXPORT openMVG-targets)
