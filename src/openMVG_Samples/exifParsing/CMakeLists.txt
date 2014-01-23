
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
SET(LIBS 
    ${EXIV2_LIBRARIES}
    OpenExif)

ADD_EXECUTABLE( main_openmvg_exifParsing exifParsing.cpp)
TARGET_LINK_LIBRARIES(main_openmvg_exifParsing
                     ${LIBS})

