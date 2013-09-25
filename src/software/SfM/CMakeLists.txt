
SET(LIBS ${EXTRA_LIBS}
  openMVG_multiview
  openMVG_image
  vlsift
  lemon
  stlplus
  ${OPENMVG_LIBRARY_DEPENDENCIES})

#Add executable that computes:
# - features and descriptors
# - putatives matches
# - geometric filtered matches
#
ADD_EXECUTABLE(openMVG_main_computeMatches main_computeMatches.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_computeMatches
                     ${LIBS} flann_cpp_s)

#Add executable that perform SfM:
# - the required input is matches computed by openMVG_main_computeMatches
#
ADD_EXECUTABLE(openMVG_main_IncrementalSfM main_IncrementalSfM.cpp
	SfMEngine.hpp SfMIncrementalEngine.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_IncrementalSfM
                      ${LIBS} ceres)

#---------
# Tools to visualize feature tracking data
# - View of the extracted Keypoints
# - View of the computed matches (putatives, geometrics) per image pair
# - View of the tracks per image pair
#--------

# - View of the extracted Keypoints
#
ADD_EXECUTABLE(openMVG_main_exportKeypoints main_exportKeypoints.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_exportKeypoints
                     ${LIBS})

# - View of the computed matches (putatives, geometrics) per image pair
#
ADD_EXECUTABLE(openMVG_main_exportMatches main_exportMatches.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_exportMatches
                     ${LIBS})

# - View of the tracks per image pair
#
ADD_EXECUTABLE(openMVG_main_exportTracks main_exportTracks.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_exportTracks
                     ${LIBS})

# - Export a computed SfM_Ouput scene to pmvs format
#
ADD_EXECUTABLE(openMVG_main_openMVG2PMVS main_openMVG2PMVS.cpp ../SfMViewer/Ply.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_openMVG2PMVS
                     ${LIBS})

# -
# Use openCV Feature detector with the openMVG matching framework
# -
IF(USE_OPENCV)
  FIND_PACKAGE( OpenCV )

  if(NOT OpenCV_FOUND)
    MESSAGE(STATUS "WARNING: OpenCV was not found. ")
  else(NOT OpenCV_FOUND)

  ADD_EXECUTABLE(openMVG_main_computeMatches_OpenCV main_computeMatches_OpenCV.cpp)
  TARGET_LINK_LIBRARIES(openMVG_main_computeMatches_OpenCV
  ${LIBS} flann_cpp_s ${OpenCV_LIBS})

  ENDIF(NOT OpenCV_FOUND)
ENDIF(USE_OPENCV)

# -
# Export a SfM openMVG scene to meshlab scene with rasters
# -
ADD_EXECUTABLE(openMVG_main_openMVG2MESHLAB main_openMVG2MESHLAB.cpp ../SfMViewer/Ply.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_openMVG2MESHLAB
                     ${LIBS})
