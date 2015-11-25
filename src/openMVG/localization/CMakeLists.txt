
file(
  GLOB_RECURSE
  localization_files_header
  *.hpp
)
file(
  GLOB_RECURSE
  localization_files_cpp
  *.cpp
)
#file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#
##Remove the unit test files (not been used by the library)
#list(REMOVE_ITEM localization_files_cpp ${REMOVEFILESUNITTEST})

ADD_LIBRARY(openMVG_localization
  ${localization_files_header}
  ${localization_files_cpp})
TARGET_LINK_LIBRARIES(openMVG_localization openMVG_sfm openMVG_voctree openMVG_numeric ${BOOST_LIBRARIES})
IF(OpenMVG_USE_CCTAG)
  TARGET_LINK_LIBRARIES(openMVG_localization ${CCTAG_LIBRARIES})
ENDIF(OpenMVG_USE_CCTAG)
SET_PROPERTY(TARGET openMVG_localization PROPERTY FOLDER OpenMVG)
INSTALL(TARGETS openMVG_localization DESTINATION lib EXPORT openMVG-targets)
