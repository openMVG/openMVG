
###
# Localize images in an existing SfM reconstruction
###
add_executable(openMVG_main_SfM_Localization main_SfM_Localization.cpp)
target_link_libraries(openMVG_main_SfM_Localization
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_matching_image_collection
  openMVG_sfm
  easyexif
  vlsift
  )

# Installation rules
set_property(TARGET openMVG_main_SfM_Localization PROPERTY FOLDER OpenMVG/software)
install(TARGETS openMVG_main_SfM_Localization DESTINATION bin/)
