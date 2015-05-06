
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

ADD_LIBRARY(openMVG_matching_image_collection
  ${matching_collection_images_files_header}
  ${matching_collection_images_files_cpp})
SET_PROPERTY(TARGET openMVG_matching_image_collection PROPERTY FOLDER OpenMVG)
INSTALL(TARGETS openMVG_matching_image_collection DESTINATION lib/)

UNIT_TEST(openMVG Pair_Builder "")
