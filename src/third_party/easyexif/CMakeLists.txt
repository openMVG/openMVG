PROJECT (easyexif CXX)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

ADD_LIBRARY (easyexif exif.h exif.cpp)
set_property(TARGET easyexif PROPERTY FOLDER OpenMVG/3rdParty)
