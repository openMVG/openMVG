project (openMVG_easyexif CXX)

add_library(openMVG_easyexif STATIC exif.h exif.cpp)
set_property(TARGET openMVG_easyexif PROPERTY FOLDER OpenMVG/3rdParty)

install(TARGETS openMVG_easyexif DESTINATION lib EXPORT openMVG-targets)
