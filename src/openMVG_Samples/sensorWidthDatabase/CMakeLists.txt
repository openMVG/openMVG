
# - Parse database
ADD_EXECUTABLE(openMVG_main_ParseDatabase main_ParseDatabase.cpp )
TARGET_LINK_LIBRARIES(openMVG_main_ParseDatabase
                     stlplus)
SET_PROPERTY(TARGET openMVG_main_ParseDatabase PROPERTY FOLDER OpenMVG/Samples)

