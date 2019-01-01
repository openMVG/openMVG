
file(
  GLOB_RECURSE
  matching_collection_images_files_header
  *.hpp
)
file(
  GLOB_RECURSE
  matching_collection_images_files_cpp
  *.cpp
)
file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)

#Remove the unit test files (not been used by the library)
list(REMOVE_ITEM matching_collection_images_files_cpp ${REMOVEFILESUNITTEST})

add_library(openMVG_matching_image_collection
  ${matching_collection_images_files_header}
  ${matching_collection_images_files_cpp})
target_compile_features(openMVG_matching_image_collection INTERFACE ${CXX11_FEATURES})

target_link_libraries(openMVG_matching_image_collection
  PUBLIC
    openMVG_matching
    openMVG_multiview
    ${OPENMVG_LIBRARY_DEPENDENCIES})
target_include_directories(openMVG_matching_image_collection
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/openMVG>
)
set_target_properties(openMVG_matching_image_collection PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_matching_image_collection PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_matching_image_collection DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG Pair_Builder "openMVG_matching_image_collection")
