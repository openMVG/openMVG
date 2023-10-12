set(sfm_localization_files_header
  SfM_Localizer.hpp
  SfM_Localizer_Single_3DTrackObservation_Database.hpp
)

set(sfm_localization_files_cpp
  SfM_Localizer.cpp
  SfM_Localizer_Single_3DTrackObservation_Database.cpp
)

add_library(openMVG_sfm_localization ${sfm_localization_files_header} ${sfm_localization_files_cpp})
target_compile_features(openMVG_sfm_localization INTERFACE ${CXX11_FEATURES})

target_link_libraries(openMVG_sfm_localization
  PUBLIC
#     openMVG_geometry
#     openMVG_multiview
    openmvg_sfm_base
#     openMVG_graph
#     openMVG_matching
#     ${cereal_TARGET}
)

target_include_directories(openMVG_sfm_localization
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)

set_target_properties(openMVG_sfm_localization PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_sfm_localization PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_sfm_localization DESTINATION lib EXPORT openMVG-targets)
