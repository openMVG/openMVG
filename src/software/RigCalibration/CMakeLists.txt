IF(OpenMVG_BUILD_VOCTREE)

    IF(OpenMVG_USE_CCTAG)
      ###
      # Calibrate a rig based on CCTag localizers
      ###
      ADD_EXECUTABLE(openMVG_main_rigCalibration main_rigCalibration.cpp)
      INSTALL(TARGETS openMVG_main_rigCalibration DESTINATION bin/)
      TARGET_LINK_LIBRARIES(openMVG_main_rigCalibration openMVG_rig openMVG_localization  openMVG_dataio openMVG_image openMVG_features vlsift ${CCTAG_LIBRARIES} ${BOOST_LIBRARIES}) #todo: delete vlsift
    ENDIF(OpenMVG_USE_CCTAG)

ENDIF(OpenMVG_BUILD_VOCTREE)
