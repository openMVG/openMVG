IF(OpenMVG_BUILD_VOCTREE)
    ###
    # Test the localization script
    ###
    ADD_EXECUTABLE(openMVG_main_localization main_localization.cpp)
    INSTALL(TARGETS openMVG_main_localization DESTINATION bin/)
    TARGET_LINK_LIBRARIES(openMVG_main_localization openMVG_localization openMVG_image openMVG_features vlsift ${BOOST_LIBRARIES})

    ###
    # Test the localization class
    ###
    ADD_EXECUTABLE(openMVG_main_voctreLocalizer main_voctreeLocalizer.cpp)
    INSTALL(TARGETS openMVG_main_voctreLocalizer DESTINATION bin/)
    TARGET_LINK_LIBRARIES(openMVG_main_voctreLocalizer openMVG_localization openMVG_dataio openMVG_image openMVG_features vlsift ${BOOST_LIBRARIES})
    if(OpenMVG_USE_ALEMBIC)
        target_link_libraries(openMVG_main_voctreLocalizer ${ABC_LIBRARIES})
    endif()
    ###
    # Create the vocabulary tree
    ###
    ADD_EXECUTABLE(openMVG_main_createVoctree createVoctree.cpp)
    INSTALL(TARGETS openMVG_main_createVoctree DESTINATION bin/)
    TARGET_LINK_LIBRARIES(openMVG_main_createVoctree openMVG_voctree openMVG_features ${BOOST_LIBRARIES})

    ###
    # generate the pair lists
    ###
    ADD_EXECUTABLE(openMVG_main_generatePairList pairListFromVoctree.cpp)
    INSTALL(TARGETS openMVG_main_generatePairList DESTINATION bin/)
    TARGET_LINK_LIBRARIES(openMVG_main_generatePairList openMVG_voctree ${BOOST_LIBRARIES})


    #@TODO decide where to put the binaries

    # Installation rules
    #SET_PROPERTY(TARGET openMVG_main_createVoctree PROPERTY FOLDER OpenMVG/software)
    #INSTALL(TARGETS openMVG_main_createVoctree DESTINATION bin/)
    #SET_PROPERTY(TARGET openMVG_main_generatePairList PROPERTY FOLDER OpenMVG/software)
    #INSTALL(TARGETS openMVG_main_generatePairList DESTINATION bin/)

ENDIF(OpenMVG_BUILD_VOCTREE)

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
