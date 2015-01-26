
ADD_EXECUTABLE(openMVG_main_CreateList main_CreateList.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_CreateList
                     ${OpenMVG_LIBS})
SET_PROPERTY(TARGET openMVG_main_CreateList PROPERTY FOLDER OpenMVG/software)

INSTALL(TARGETS openMVG_main_CreateList DESTINATION bin/)

#Add executable that computes:
# - features and descriptors
# - putatives matches
# - geometric filtered matches
#
ADD_EXECUTABLE(openMVG_main_computeMatches main_computeMatches.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_computeMatches
                     ${OpenMVG_LIBS})
SET_PROPERTY(TARGET openMVG_main_computeMatches PROPERTY FOLDER OpenMVG/software)

INSTALL(TARGETS openMVG_main_computeMatches DESTINATION bin/)

#Add library and executable that perform SfM:
# - the required input is matches computed by openMVG_main_computeMatches
#
ADD_LIBRARY(openMVG_SfM_Sequential SfMIncrementalEngine.cpp)
ADD_EXECUTABLE(openMVG_main_IncrementalSfM main_IncrementalSfM.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_IncrementalSfM
                      openMVG_SfM_Sequential
                      ${OpenMVG_LIBS})
SET_PROPERTY(TARGET openMVG_main_IncrementalSfM PROPERTY FOLDER OpenMVG/software)
SET_PROPERTY(TARGET openMVG_SfM_Sequential PROPERTY FOLDER OpenMVG/OpenMVG)

INSTALL(TARGETS openMVG_SfM_Sequential DESTINATION lib/)
INSTALL(TARGETS openMVG_main_IncrementalSfM DESTINATION bin/)

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
SET_PROPERTY(TARGET openMVG_main_exportKeypoints PROPERTY FOLDER OpenMVG/software)

INSTALL(TARGETS openMVG_main_exportKeypoints DESTINATION bin/)

# - View of the computed matches (putatives, geometrics) per image pair
#
ADD_EXECUTABLE(openMVG_main_exportMatches main_exportMatches.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_exportMatches
                     ${OpenMVG_LIBS})
SET_PROPERTY(TARGET openMVG_main_exportMatches PROPERTY FOLDER OpenMVG/software)

INSTALL(TARGETS openMVG_main_exportMatches DESTINATION bin/)

# - View of the tracks per image pair
#
ADD_EXECUTABLE(openMVG_main_exportTracks main_exportTracks.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_exportTracks
                     ${OpenMVG_LIBS})
SET_PROPERTY(TARGET openMVG_main_exportTracks PROPERTY FOLDER OpenMVG/software)

INSTALL(TARGETS openMVG_main_exportTracks DESTINATION bin/)

# - Export a computed SfM_Ouput scene to PMVS format
#
ADD_EXECUTABLE(openMVG_main_openMVG2PMVS main_openMVG2PMVS.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_openMVG2PMVS
                     ${OpenMVG_LIBS})
SET_PROPERTY(TARGET openMVG_main_openMVG2PMVS PROPERTY FOLDER OpenMVG/software)

INSTALL(TARGETS openMVG_main_openMVG2PMVS DESTINATION bin/)

# - Export a computed SfM_Ouput scene to CMPMVS format
#
ADD_EXECUTABLE(openMVG_main_openMVG2CMPMVS main_openMVG2CMPMVS.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_openMVG2CMPMVS
                     ${OpenMVG_LIBS})
SET_PROPERTY(TARGET openMVG_main_openMVG2CMPMVS PROPERTY FOLDER OpenMVG/software)

INSTALL(TARGETS openMVG_main_openMVG2CMPMVS DESTINATION bin/)

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

    SET_PROPERTY(TARGET openMVG_main_computeMatches_OpenCV PROPERTY FOLDER OpenMVG/software)

    INSTALL(TARGETS openMVG_main_computeMatches_OpenCV DESTINATION bin/)

  ENDIF(NOT OpenCV_FOUND)
ENDIF(USE_OPENCV)

# -
# Export a SfM openMVG scene to meshlab scene with rasters
# -
ADD_EXECUTABLE(openMVG_main_openMVG2MESHLAB main_openMVG2MESHLAB.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_openMVG2MESHLAB
                     ${OpenMVG_LIBS})
SET_PROPERTY(TARGET openMVG_main_openMVG2MESHLAB PROPERTY FOLDER OpenMVG/software)

INSTALL(TARGETS openMVG_main_openMVG2MESHLAB DESTINATION bin/)

# - Quality comparison against a GT camera path (MultiView Evaluation dataset)
#
ADD_EXECUTABLE(openMVG_main_evalQuality main_evalQuality.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_evalQuality
                     ${OpenMVG_LIBS})
SET_PROPERTY(TARGET openMVG_main_evalQuality PROPERTY FOLDER OpenMVG/software)

INSTALL(TARGETS openMVG_main_evalQuality DESTINATION bin/)

                     
UNIT_TEST(openMVG SfMIOHelper "stlplus")

INSTALL(FILES cameraSensorWidth/cameraGenerated.txt DESTINATION resources/)

