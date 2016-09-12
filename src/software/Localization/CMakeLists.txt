
###
# Localize images in an existing SfM reconstruction
###
ADD_EXECUTABLE(openMVG_main_SfM_Localization main_SfM_Localization.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_SfM_Localization
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_matching_image_collection
  openMVG_sfm
  easyexif
  vlsift
  )

# Installation rules
SET_PROPERTY(TARGET openMVG_main_SfM_Localization PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_SfM_Localization DESTINATION bin/)

###
# Localize images in an existing SfM reconstruction (Multiple images)
###
ADD_EXECUTABLE(openMVG_main_SfM_Localization_MultipleImages main_SfM_Localization_MultipleImages.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_SfM_Localization_MultipleImages
  openMVG_system
  openMVG_image
  openMVG_features
  openMVG_matching_image_collection
  openMVG_sfm
  easyexif
  vlsift
  )

# Installation rules
SET_PROPERTY(TARGET openMVG_main_SfM_Localization_MultipleImages PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_SfM_Localization_MultipleImages DESTINATION bin/)
