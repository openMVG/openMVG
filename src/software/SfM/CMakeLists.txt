
###
# Intrinsic image analysis and SfM_Data container initialization
###
ADD_EXECUTABLE(openMVG_main_SfMInit_ImageListing main_SfMInit_ImageListing.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_SfMInit_ImageListing ${OpenMVG_LIBRARIES})

#convert a v0.6 lists.txt file to the new sfm_data.X format
ADD_EXECUTABLE(openMVG_main_ConvertList main_ConvertList.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_ConvertList ${OpenMVG_LIBRARIES})

# Installation rules
SET_PROPERTY(TARGET openMVG_main_SfMInit_ImageListing PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_SfMInit_ImageListing DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_ConvertList PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_ConvertList DESTINATION bin/)

###
# Add executable that computes:
# - openMVG_main_ComputeFeatures: features and descriptors
# - openMVG_main_ComputeMatches: putatives matches + geometric filtered matches
###

ADD_EXECUTABLE(openMVG_main_ComputeFeatures main_ComputeFeatures.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_ComputeFeatures ${OpenMVG_LIBRARIES})

ADD_EXECUTABLE(openMVG_main_ComputeMatches main_ComputeMatches.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_ComputeMatches ${OpenMVG_LIBRARIES})

# Installation rules
SET_PROPERTY(TARGET openMVG_main_ComputeFeatures PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_ComputeFeatures DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_ComputeMatches PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_ComputeMatches DESTINATION bin/)

###
# SfM Pipelines
# - Incremental/Sequential
# - Global
###
# - Convert SfM_Data format (from one to another)
# - Pair filtering (from frustum)
# - Compute structure from known camera poses
# - Compute structure color
###
ADD_EXECUTABLE(openMVG_main_IncrementalSfM main_IncrementalSfM.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_IncrementalSfM ${OpenMVG_LIBRARIES})

ADD_EXECUTABLE(openMVG_main_GlobalSfM main_GlobalSfM.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_GlobalSfM ${OpenMVG_LIBRARIES})

ADD_EXECUTABLE(openMVG_main_ConvertSfM_DataFormat main_ConvertSfM_DataFormat.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_ConvertSfM_DataFormat ${OpenMVG_LIBRARIES})

ADD_EXECUTABLE(openMVG_main_FrustumFiltering main_FrustumFiltering.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_FrustumFiltering ${OpenMVG_LIBRARIES})

ADD_EXECUTABLE(openMVG_main_ComputeStructureFromKnownPoses main_ComputeStructureFromKnownPoses.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_ComputeStructureFromKnownPoses ${OpenMVG_LIBRARIES})

ADD_EXECUTABLE(openMVG_main_ComputeSfM_DataColor main_ComputeSfM_DataColor.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_ComputeSfM_DataColor ${OpenMVG_LIBRARIES})

# Installation rules
INSTALL(FILES cameraSensorWidth/cameraGenerated.txt DESTINATION share/openMVG)
SET_PROPERTY(TARGET openMVG_main_IncrementalSfM PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_IncrementalSfM DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_GlobalSfM PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_GlobalSfM DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_ConvertSfM_DataFormat PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_ConvertSfM_DataFormat DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_FrustumFiltering PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_FrustumFiltering DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_ComputeStructureFromKnownPoses PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_ComputeStructureFromKnownPoses DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_ComputeSfM_DataColor PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_ComputeSfM_DataColor DESTINATION bin/)

###
# SfM tools to visualize feature tracking data
###

# - View extracted Keypoints
#
ADD_EXECUTABLE(openMVG_main_exportKeypoints main_exportKeypoints.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_exportKeypoints ${OpenMVG_LIBRARIES})


# - View computed matches (putatives, geometrics) per image pair
#
ADD_EXECUTABLE(openMVG_main_exportMatches main_exportMatches.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_exportMatches ${OpenMVG_LIBRARIES})

# - View tracks per image pair
#
ADD_EXECUTABLE(openMVG_main_exportTracks main_exportTracks.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_exportTracks ${OpenMVG_LIBRARIES})

# - Export undistorted images related to a sfm_data file
#
ADD_EXECUTABLE(openMVG_main_ExportUndistortedImages main_ExportUndistortedImages.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_ExportUndistortedImages ${OpenMVG_LIBRARIES})

# installation rules
SET_PROPERTY(TARGET openMVG_main_exportKeypoints PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_exportKeypoints DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_exportMatches PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_exportMatches DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_exportTracks PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_exportTracks DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_ExportUndistortedImages PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_ExportUndistortedImages DESTINATION bin/)

###
# SfM export to X
# - PMVS + bundler data format
# - CMPMVS
# - Meshlab
###

# - Export a SfM openMVG scene to PMVS format
#
ADD_EXECUTABLE(openMVG_main_openMVG2PMVS main_openMVG2PMVS.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_openMVG2PMVS ${OpenMVG_LIBRARIES})

# - Export a SfM openMVG scene to CMPMVS format
#
ADD_EXECUTABLE(openMVG_main_openMVG2CMPMVS main_openMVG2CMPMVS.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_openMVG2CMPMVS ${OpenMVG_LIBRARIES})

# - Export a SfM openMVG scene to Alembic format
#
FIND_PACKAGE(Alembic)
if (ALEMBIC_FOUND)
  ADD_EXECUTABLE(openMVG_main_openMVG2ABC main_openMVG2ABC.cpp)
  INCLUDE_DIRECTORIES(${ABC_INCLUDE_DIR})
  TARGET_LINK_LIBRARIES(openMVG_main_openMVG2ABC ${ABC_LIBRARIES} ${OpenMVG_LIBRARIES})
  SET_PROPERTY(TARGET openMVG_main_openMVG2ABC PROPERTY FOLDER OpenMVG/software)
  INSTALL(TARGETS openMVG_main_openMVG2ABC DESTINATION bin/)
endif(ALEMBIC_FOUND)

# - Export a SfM openMVG scene to meshlab scene with rasters
# -
ADD_EXECUTABLE(openMVG_main_openMVG2MESHLAB main_openMVG2MESHLAB.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_openMVG2MESHLAB ${OpenMVG_LIBRARIES})

# - Export a SfM openMVG scene to mvs-texturing scene folder
# -
ADD_EXECUTABLE(openMVG_main_openMVG2MVSTEXTURING main_openMVG2MVSTEXTURING.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_openMVG2MVSTEXTURING ${OpenMVG_LIBRARIES})

# - Export SfM openMVG camera scene as triangle meshes
# -
ADD_EXECUTABLE(openMVG_main_ExportCameraFrustums main_ExportCameraFrustums.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_ExportCameraFrustums ${OpenMVG_LIBRARIES})

# installation rules
SET_PROPERTY(TARGET openMVG_main_openMVG2PMVS PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_openMVG2PMVS DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_openMVG2CMPMVS PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_openMVG2CMPMVS DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_openMVG2MESHLAB PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_openMVG2MESHLAB DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_openMVG2MVSTEXTURING PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_openMVG2MVSTEXTURING DESTINATION bin/)
SET_PROPERTY(TARGET openMVG_main_ExportCameraFrustums PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_ExportCameraFrustums DESTINATION bin/)

# -
# Use openCV Feature detector to describe the images
# -
IF(OpenMVG_USE_OPENCV)
  FIND_PACKAGE( OpenCV )
  INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
  IF(NOT OpenCV_FOUND)
    MESSAGE(STATUS "WARNING: OpenCV was not found. ")
  ELSE(NOT OpenCV_FOUND)

    ADD_EXECUTABLE(openMVG_main_ComputeFeatures_OpenCV
      main_ComputeFeatures_OpenCV.cpp)
    TARGET_LINK_LIBRARIES(openMVG_main_ComputeFeatures_OpenCV
                          ${OpenMVG_LIBRARIES}
                          ${OpenCV_LIBS})

    SET_PROPERTY(TARGET openMVG_main_ComputeFeatures_OpenCV PROPERTY FOLDER OpenMVG/software)
    INSTALL(TARGETS openMVG_main_ComputeFeatures_OpenCV DESTINATION bin/)

  ENDIF(NOT OpenCV_FOUND)
ENDIF(OpenMVG_USE_OPENCV)


###
# SfM Research tools
#
# - Quality comparison against a GT camera path (MultiView Evaluation dataset)
###
ADD_EXECUTABLE(openMVG_main_evalQuality main_evalQuality.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_evalQuality ${OpenMVG_LIBRARIES})

#installation rules
SET_PROPERTY(TARGET openMVG_main_evalQuality PROPERTY FOLDER OpenMVG/software)
INSTALL(TARGETS openMVG_main_evalQuality DESTINATION bin/)

##
# Export tutorial with valid path thanks to a configuration file
##
SET(OPENMVG_SOFTWARE_SFM_SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
IF(MSVC OR APPLE)
  SET(OPENMVG_SOFTWARE_SFM_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE}")
ELSE()
  SET(OPENMVG_SOFTWARE_SFM_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
ENDIF()
SET(OPENMVG_TUTORIAL_IN_FILE "${CMAKE_CURRENT_LIST_DIR}/tutorial_demo.py.in")
CONFIGURE_FILE(${OPENMVG_TUTORIAL_IN_FILE}
    "${CMAKE_CURRENT_BINARY_DIR}/tutorial_demo.py" @ONLY)

SET(OPENMVG_TUTORIAL_IN_FILE "${CMAKE_CURRENT_LIST_DIR}/SfM_GlobalPipeline.py.in")
CONFIGURE_FILE(${OPENMVG_TUTORIAL_IN_FILE}
    "${CMAKE_CURRENT_BINARY_DIR}/SfM_GlobalPipeline.py" @ONLY)

SET(OPENMVG_TUTORIAL_IN_FILE "${CMAKE_CURRENT_LIST_DIR}/SfM_SequentialPipeline.py.in")
CONFIGURE_FILE(${OPENMVG_TUTORIAL_IN_FILE}
    "${CMAKE_CURRENT_BINARY_DIR}/SfM_SequentialPipeline.py" @ONLY)

