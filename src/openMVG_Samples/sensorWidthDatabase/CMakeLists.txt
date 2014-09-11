
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../software/SfM/cameraSensorWidth")

# - Parse database
ADD_EXECUTABLE(openMVG_main_ParseDatabase main_ParseDatabase.cpp )
TARGET_LINK_LIBRARIES(openMVG_main_ParseDatabase
                     stlplus)

UNIT_TEST(openMVG ParseDatabase "stlplus")

SET_PROPERTY(TARGET openMVG_main_ParseDatabase PROPERTY FOLDER OpenMVG/Samples)

