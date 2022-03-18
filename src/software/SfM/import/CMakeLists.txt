
add_executable(openMVG_main_SfMInit_ImageListingFromKnownPoses
  main_SfMInit_ImageListingFromKnownPoses.cpp
    io_readGTInterface.hpp
    io_readGTStrecha.hpp
    io_readGTMiddleBury.hpp
    io_readGTDTUMVS.hpp
    io_readGTETH3D.hpp
    io_readGTKitti.hpp
)
target_link_libraries(openMVG_main_SfMInit_ImageListingFromKnownPoses
  PRIVATE
    openMVG_system
    openMVG_image
    openMVG_features
    openMVG_sfm
    openMVG_exif
    ${STLPLUS_LIBRARY}
)

# Installation rules
set_property(TARGET openMVG_main_SfMInit_ImageListingFromKnownPoses PROPERTY FOLDER OpenMVG/software)
install(TARGETS openMVG_main_SfMInit_ImageListingFromKnownPoses DESTINATION bin/)

SET(OPENMVG_TUTORIAL_IN_FILE "${CMAKE_CURRENT_LIST_DIR}/SfM_StructurePipeline.py.in")
CONFIGURE_FILE(${OPENMVG_TUTORIAL_IN_FILE}
    "${CMAKE_CURRENT_BINARY_DIR}/SfM_StructurePipeline.py" @ONLY)
