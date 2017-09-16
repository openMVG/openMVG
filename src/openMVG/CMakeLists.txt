add_subdirectory(cameras)
add_subdirectory(exif)
add_subdirectory(features)
add_subdirectory(graph)
add_subdirectory(image)
add_subdirectory(linearProgramming)
add_subdirectory(geodesy)
add_subdirectory(geometry)
add_subdirectory(matching)
add_subdirectory(matching_image_collection)
add_subdirectory(multiview)
add_subdirectory(numeric)
add_subdirectory(robust_estimation)
add_subdirectory(tracks)
add_subdirectory(color_harmonization)
add_subdirectory(system)
add_subdirectory(sfm)
add_subdirectory(stl)

#INSTALL RULES
install(
  DIRECTORY .
  DESTINATION include/openMVG
  COMPONENT headers
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)
