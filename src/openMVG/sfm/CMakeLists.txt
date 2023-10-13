add_subdirectory(base)
add_subdirectory(pipelines)

set( # TODO: each pipelines/subfolder its own library/build module
  sfm_files_header
  sfm.hpp
  # sequential is built into its own lib
  #   pipelines/sequential/SfmSceneInitializerStellar.hpp
  #   pipelines/sequential/sequential_SfM_util.hpp
  #   pipelines/sequential/sequential_SfM_base.hpp
  #   pipelines/sequential/SfmSceneInitializerMaxPair.hpp
  #   pipelines/sequential/SfmSceneInitializer.hpp
  #   pipelines/sequential/sequential_SfM.hpp
  #   pipelines/sequential/sequential_SfM2.hpp
  pipelines/global/mutexSet.hpp
  pipelines/global/GlobalSfM_rotation_averaging.hpp
  pipelines/global/sfm_global_engine_relative_motions.hpp
  pipelines/global/GlobalSfM_translation_averaging.hpp
  pipelines/global/triplet_t_ACRansac_kernelAdaptator.hpp
  pipelines/global/sfm_global_reindex.hpp
  pipelines/structure_from_known_poses/structure_estimator.hpp
)

set(
  sfm_files_cpp
  # GLOB_RECURSE IST EVILN
  #  GLOB_RECURSE 
  # sequential is built into its own lib
  #   pipelines/sequential/sequential_SfM.cpp
  #   pipelines/sequential/sequential_SfM2.cpp
  #   pipelines/sequential/sequential_SfM_base.cpp
  #   pipelines/sequential/SfmSceneInitializerMaxPair.cpp
  #   pipelines/sequential/SfmSceneInitializerStellar.cpp
  #   pipelines/sequential/sequential_SfM_util.cpp
  pipelines/global/GlobalSfM_translation_averaging.cpp
  pipelines/global/sfm_global_engine_relative_motions.cpp
  pipelines/global/GlobalSfM_rotation_averaging.cpp
  pipelines/structure_from_known_poses/structure_estimator.cpp
)

add_library(openMVG_sfm ${sfm_files_header} ${sfm_files_cpp})
target_compile_features(openMVG_sfm INTERFACE ${CXX11_FEATURES})
target_link_libraries(openMVG_sfm
  PUBLIC
    openMVG_geometry
    openMVG_features
    openMVG_graph
    openMVG_matching
    openMVG_sfm_base
    openMVG_sfm_stellar
    openMVG_sfm_sequential
    ${cereal_TARGET}
    ${OPENMVG_LIBRARY_DEPENDENCIES}
)

target_link_libraries(openMVG_sfm
  LINK_PRIVATE
    openMVG_image
    openMVG_lInftyComputerVision
    openMVG_system
    ${CERES_LIBRARIES}
    ${STLPLUS_LIBRARY}
)
target_include_directories(openMVG_sfm
  PRIVATE
    ${CERES_INCLUDE_DIRS}
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)
if (MSVC)
  set_target_properties(openMVG_sfm PROPERTIES COMPILE_FLAGS "/bigobj")
endif (MSVC)
set_target_properties(openMVG_sfm PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_sfm PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_sfm DESTINATION ${CMAKE_INSTALL_LIBDIR} EXPORT openMVG-targets)

