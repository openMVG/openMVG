
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_main_exifParsing exifParsing.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_exifParsing easyexif)

SET_PROPERTY(TARGET openMVG_main_exifParsing PROPERTY FOLDER OpenMVG/Samples)
