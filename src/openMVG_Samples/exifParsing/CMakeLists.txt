
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE( main_openmvg_exifParsing exifParsing.cpp)
TARGET_LINK_LIBRARIES(main_openmvg_exifParsing
                     ${LIBS} easyexif)

SET_PROPERTY(TARGET main_openmvg_exifParsing PROPERTY FOLDER OpenMVG/Samples)
