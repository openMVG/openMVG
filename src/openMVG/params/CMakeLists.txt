
file(
	GLOB_RECURSE
	params_files_header
	*.hpp
)
file(
	GLOB_RECURSE
	params_files_cpp
	*.cpp
)
file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM params_files_cpp ${REMOVEFILESUNITTEST})
#Remove the test header main files
file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.hpp)
#list(REMOVE_ITEM params_files_header ${REMOVEFILESUNITTEST})

add_library(openMVG_params ${params_files_header} ${params_files_cpp})
target_link_libraries(
  openMVG_params
  openMVG_multiview
  stlplus
)
set_target_properties(openMVG_params PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")

install(TARGETS openMVG_params DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG params_io
  "openMVG_params;openMVG_features;openMVG_sfm;openMVG_system;stlplus")
