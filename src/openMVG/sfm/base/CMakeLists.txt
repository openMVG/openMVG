set(
  sfm_base_files_header
  sfm_landmark.hpp
  sfm_report.hpp
  sfm_data_io_baf.hpp
  sfm_data_BA_ceres.hpp
  sfm_data_utils.hpp
  sfm_view_priors_io.hpp
  sfm_data_BA.hpp
  sfm_data_io_ply.hpp
  sfm_view.hpp
  sfm_data_graph_utils.hpp
  sfm_data.hpp
  sfm_data_transform.hpp
  sfm_data_filters.hpp
  sfm_view_priors.hpp
  sfm_landmark_io.hpp
  sfm_data_io.hpp
  sfm_data_filters_frustum.hpp
  sfm_data_BA_ceres_camera_functor.hpp
  sfm_data_triangulation.hpp
  sfm_data_colorization.hpp
  sfm_filters.hpp
  sfm_data_io_cereal.hpp
  sfm_view_io.hpp
)

set(
  sfm_base_files_cpp
  sfm_data_io.cpp
  sfm_data_filters_frustum.cpp
  sfm_data_colorization.cpp
  sfm_data_io_cereal.cpp
  sfm_data_triangulation.cpp
  sfm_data_transform.cpp
  sfm_data_filters.cpp
  sfm_data_graph_utils.cpp
  sfm_report.cpp
  sfm_data_BA_ceres.cpp
  sfm_data_utils.cpp
)

add_library(openMVG_sfm_base ${sfm_base_files_header} ${sfm_base_files_cpp})
target_compile_features(openMVG_sfm_base INTERFACE ${CXX11_FEATURES})
target_link_libraries(openMVG_sfm_base
  PUBLIC
    openMVG_geometry
    openMVG_multiview
#     openMVG_features
#     openMVG_graph
#     openMVG_matching
#     ${cereal_TARGET}
#     ${OPENMVG_LIBRARY_DEPENDENCIES}
)

target_link_libraries(openMVG_sfm_base
  LINK_PRIVATE
    openMVG_image
#     openMVG_lInftyComputerVision
#     openMVG_system
    ${CERES_LIBRARIES}
    ${STLPLUS_LIBRARY}
)
target_include_directories(openMVG_sfm_base
  PRIVATE
    ${CERES_INCLUDE_DIRS}
    ${LEMON_INCLUDE_DIRS}
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)

set_target_properties(openMVG_sfm_base PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_sfm_base PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_sfm_base DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG sfm_data_io "openMVG_sfm_base;${STLPLUS_LIBRARY}")
UNIT_TEST(openMVG sfm_data_filters "openMVG_sfm_base")
UNIT_TEST(openMVG sfm_data_graph_utils "openMVG_sfm_base")
