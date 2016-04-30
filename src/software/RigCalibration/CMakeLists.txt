# Rig calibration
IF(OpenMVG_BUILD_VOCTREE)

      ###
      # Calibrate a rig based on CCTag localizers
      ###
      ADD_EXECUTABLE(openMVG_main_rigCalibration main_rigCalibration.cpp)
      INSTALL(TARGETS openMVG_main_rigCalibration DESTINATION bin/)
      TARGET_LINK_LIBRARIES(openMVG_main_rigCalibration openMVG_rig openMVG_localization  openMVG_dataio openMVG_image openMVG_features vlsift ${CCTAG_LIBRARIES} ${BOOST_LIBRARIES}) #todo: delete vlsift
  if(OpenMVG_USE_CCTAG)
    target_link_libraries(openMVG_main_rigCalibration
      ${CCTAG_LIBRARIES}
      )
  endif(OpenMVG_USE_CCTAG)

ENDIF(OpenMVG_BUILD_VOCTREE)

# Camera calibration
if(OpenCV_FOUND)
  message(STATUS "OpenCV found, will build openMVG_main_cameraCalibration")
  add_executable(openMVG_main_cameraCalibration main_cameraCalibration.cpp)
  target_link_libraries(openMVG_main_cameraCalibration ${OpenCV_LIBRARIES})
  install(TARGETS openMVG_main_cameraCalibration DESTINATION bin/)
endif(OpenCV_FOUND)
