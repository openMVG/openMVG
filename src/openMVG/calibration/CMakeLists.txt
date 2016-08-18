set(calibration_files_headers
  calibration.hpp
  bestImages.hpp
  patternDetect.hpp
  exportData.hpp
)
set(calibration_files_sources
  calibration.cpp
  bestImages.cpp
  patternDetect.cpp
  exportData.cpp
)

ADD_LIBRARY(openMVG_calibration
  ${calibration_files_headers}
  ${calibration_files_sources})
TARGET_LINK_LIBRARIES(openMVG_calibration openMVG_image openMVG_dataio ${BOOST_LIBRARIES}  ${OpenCV_LIBS})
IF(OpenMVG_USE_CCTAG)
  TARGET_LINK_LIBRARIES(openMVG_calibration ${CCTAG_LIBRARIES})
ENDIF(OpenMVG_USE_CCTAG)

SET_PROPERTY(TARGET openMVG_calibration PROPERTY FOLDER OpenMVG)

INSTALL(TARGETS openMVG_calibration DESTINATION lib EXPORT openMVG-targets)