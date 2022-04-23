
file(
  GLOB
  multiview_trifocal_files_header
  *.hpp
)
file(
  GLOB
  multiview_trifocal_files_cpp
  *.cpp
)
file(GLOB REMOVEFILESUNITTEST *_test.cpp)
# Remove the unit test files (not been used by the library)
list(REMOVE_ITEM multiview_trifocal_files_cpp ${REMOVEFILESUNITTEST})

add_library(openMVG_multiview_trifocal
  ${multiview_trifocal_files_header}
  ${multiview_trifocal_files_cpp})
target_compile_features(openMVG_multiview_trifocal INTERFACE ${CXX11_FEATURES})
set_property(TARGET openMVG_multiview_trifocal PROPERTY FOLDER OpenMVG/OpenMVG)
# install(TARGETS openMVG_multiview_trifocal DESTINATION lib EXPORT openMVG-targets)

target_link_libraries(openMVG_multiview_trifocal
  PUBLIC
    openMVG_multiview
)
target_link_libraries(openMVG_multiview_trifocal
  LINK_PRIVATE
    minus
)
target_include_directories(openMVG_multiview_trifocal
  PRIVATE
    ${MINUS_INCLUDE_DIRS}
  PUBLIC
    $<INSTALL_INTERFACE:include>
)

set_target_properties(openMVG_multiview_trifocal PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

# Solver
UNIT_TEST(openMVG solver_trifocal_three_point "openMVG_multiview_trifocal;openMVG_multiview;minus")

# Error function
UNIT_TEST(openMVG solver_trifocal_metrics "openMVG_multiview_trifocal;openMVG_multiview;minus")

# Kernel / Robust solve
UNIT_TEST(openMVG three_view_kernel "openMVG_multiview_trifocal;openMVG_multiview;minus")
