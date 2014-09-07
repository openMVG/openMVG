
ADD_EXECUTABLE(openMVG_main_CreateList main_CreateList.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_CreateList
                     ${OpenMVG_LIBS} easyexif)

#Add executable that computes:
# - features and descriptors
# - putatives matches
# - geometric filtered matches
#
ADD_EXECUTABLE(openMVG_main_computeMatches main_computeMatches.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_computeMatches
                     ${OpenMVG_LIBS})

#Add library and executable that perform SfM:
# - the required input is matches computed by openMVG_main_computeMatches
#
ADD_LIBRARY(openMVG_SfM_Sequential SfMIncrementalEngine.cpp)
ADD_EXECUTABLE(openMVG_main_IncrementalSfM main_IncrementalSfM.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_IncrementalSfM
                      openMVG_SfM_Sequential
                      ${OpenMVG_LIBS})

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
                     ${OpenMVG_LIBS})

# - View of the computed matches (putatives, geometrics) per image pair
#
ADD_EXECUTABLE(openMVG_main_exportMatches main_exportMatches.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_exportMatches
                     ${OpenMVG_LIBS})

# - View of the tracks per image pair
#
ADD_EXECUTABLE(openMVG_main_exportTracks main_exportTracks.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_exportTracks
                     ${OpenMVG_LIBS})

# - Export a computed SfM_Ouput scene to PMVS format
#
ADD_EXECUTABLE(openMVG_main_openMVG2PMVS main_openMVG2PMVS.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_openMVG2PMVS
                     ${OpenMVG_LIBS})

# - Export a computed SfM_Ouput scene to CMPMVS format
#
ADD_EXECUTABLE(openMVG_main_openMVG2CMPMVS main_openMVG2CMPMVS.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_openMVG2CMPMVS
                     ${OpenMVG_LIBS})

# -
# Use openCV Feature detector with the openMVG matching framework
# -
IF(USE_OPENCV)
  FIND_PACKAGE( OpenCV )

  IF(NOT OpenCV_FOUND)
    MESSAGE(STATUS "WARNING: OpenCV was not found. ")
  ELSE(NOT OpenCV_FOUND)

  ADD_EXECUTABLE(openMVG_main_computeMatches_OpenCV main_computeMatches_OpenCV.cpp)
  TARGET_LINK_LIBRARIES(openMVG_main_computeMatches_OpenCV
                        ${OpenMVG_LIBS}
                        ${OpenCV_LIBS})

  ENDIF(NOT OpenCV_FOUND)
ENDIF(USE_OPENCV)

# -
# Export a SfM openMVG scene to meshlab scene with rasters
# -
ADD_EXECUTABLE(openMVG_main_openMVG2MESHLAB main_openMVG2MESHLAB.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_openMVG2MESHLAB
                     ${OpenMVG_LIBS})

# - Quality comparison against a GT camera path (MultiView Evaluation dataset)
#
ADD_EXECUTABLE(openMVG_main_evalQuality main_evalQuality.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_evalQuality
                     ${OpenMVG_LIBS})
                     
UNIT_TEST(openMVG SfMIOHelper "stlplus")
