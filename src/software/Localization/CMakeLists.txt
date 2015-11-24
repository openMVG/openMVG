IF(OpenMVG_BUILD_VOCTREE)
    ###
    # Test the localization class
    ###
    ADD_EXECUTABLE(openMVG_main_voctreeLocalizer main_voctreeLocalizer.cpp)
    INSTALL(TARGETS openMVG_main_voctreeLocalizer DESTINATION bin/)
    TARGET_LINK_LIBRARIES(openMVG_main_voctreeLocalizer openMVG_localization openMVG_dataio openMVG_image openMVG_features vlsift ${BOOST_LIBRARIES})
    if(OpenMVG_USE_ALEMBIC)
        target_link_libraries(openMVG_main_voctreeLocalizer ${ABC_LIBRARIES})
    endif()
    ###
    # Create the vocabulary tree
    ###
    ADD_EXECUTABLE(openMVG_main_createVoctree createVoctree.cpp)
    INSTALL(TARGETS openMVG_main_createVoctree DESTINATION bin/)
    TARGET_LINK_LIBRARIES(openMVG_main_createVoctree openMVG_voctree openMVG_features ${BOOST_LIBRARIES})

    IF(OpenMVG_USE_CCTAG)
      ###
      # Localize images in an existing SfM reconstruction based on CCTag features
      ###
      ADD_EXECUTABLE(openMVG_main_cctagLocalizer main_cctagLocalizer.cpp)
      INSTALL(TARGETS openMVG_main_cctagLocalizer DESTINATION bin/)
      TARGET_LINK_LIBRARIES(openMVG_main_cctagLocalizer openMVG_localization openMVG_dataio openMVG_image openMVG_features ${CCTAG_LIBRARIES} ${BOOST_LIBRARIES})
      if(OpenMVG_USE_ALEMBIC)
        target_link_libraries(openMVG_main_cctagLocalizer ${ABC_LIBRARIES})
      endif()
    ENDIF(OpenMVG_USE_CCTAG)

    ###
    # generate the pair lists
    ###
    ADD_EXECUTABLE(openMVG_main_generatePairList pairListFromVoctree.cpp)
    INSTALL(TARGETS openMVG_main_generatePairList DESTINATION bin/)
    TARGET_LINK_LIBRARIES(openMVG_main_generatePairList openMVG_voctree ${BOOST_LIBRARIES})

    ###
    # voctree query utility
    ###
    ADD_EXECUTABLE(openMVG_main_voctreeQueryUtility voctreeQueryUtility.cpp)
    INSTALL(TARGETS openMVG_main_voctreeQueryUtility DESTINATION bin/)
    TARGET_LINK_LIBRARIES(openMVG_main_voctreeQueryUtility openMVG_voctree ${BOOST_LIBRARIES})


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
