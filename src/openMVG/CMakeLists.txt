ADD_SUBDIRECTORY(bundle_adjustment)
ADD_SUBDIRECTORY(cameras)
ADD_SUBDIRECTORY(exif_IO)
ADD_SUBDIRECTORY(features)
ADD_SUBDIRECTORY(geometry)
ADD_SUBDIRECTORY(graph)
ADD_SUBDIRECTORY(image)
ADD_SUBDIRECTORY(linearProgramming)
ADD_SUBDIRECTORY(matching)
ADD_SUBDIRECTORY(matching_image_collection)
ADD_SUBDIRECTORY(multiview)
ADD_SUBDIRECTORY(numeric)
ADD_SUBDIRECTORY(robust_estimation)
ADD_SUBDIRECTORY(split)
ADD_SUBDIRECTORY(tracks)
ADD_SUBDIRECTORY(color_harmonization)
ADD_SUBDIRECTORY(system)

# set the variable in the parent_scope
SET(openMVG_linearProgramming ${openMVG_linearProgramming} PARENT_SCOPE)

#INSTALL RULES
INSTALL(
  DIRECTORY .
  DESTINATION include/openMVG
  COMPONENT headers
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)
