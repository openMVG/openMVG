
file(
  GLOB
  matching_collection_images_files_header
  *.hpp
)
file(
  GLOB
  matching_collection_images_files_cpp
  *.cpp
)
file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)

#Remove the unit test files (not been used by the library)
list(REMOVE_ITEM matching_collection_images_files_cpp ${REMOVEFILESUNITTEST})

add_library(openMVG_matching_image_collection
  ${matching_collection_images_files_header}
  ${matching_collection_images_files_cpp})
if (OpenMVG_USE_CUDA)
  target_link_libraries(openMVG_matching_image_collection openMVG_gpu_matching_image_collection)
endif()
target_link_libraries(openMVG_matching_image_collection openMVG_matching openMVG_multiview)
set_target_properties(openMVG_matching_image_collection PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_matching_image_collection PROPERTY FOLDER OpenMVG)
install(TARGETS openMVG_matching_image_collection DESTINATION lib EXPORT openMVG-targets)

UNIT_TEST(openMVG Pair_Builder "")

if (OpenMVG_USE_CUDA) 
  add_subdirectory(gpu)
endif()
