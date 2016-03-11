
# - Parse database
ADD_EXECUTABLE(openMVG_main_exif_ParseDatabase main_ParseDatabase.cpp )
TARGET_LINK_LIBRARIES(openMVG_main_exif_ParseDatabase
                     stlplus)
SET_PROPERTY(TARGET openMVG_main_exif_ParseDatabase PROPERTY FOLDER OpenMVG/Samples)

