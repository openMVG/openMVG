#---------------------
#Third party libraries
#---------------------

# Unit test library
ADD_SUBDIRECTORY(CppUnitLite)

# Basic filesystem utils
ADD_SUBDIRECTORY(stlplus3)

# Add linear algebra library
#  Eigen is used as header only

# Add graph library
ADD_SUBDIRECTORY(lemon)

# Image I/O
IF(NOT JPEG_FOUND)
  ADD_SUBDIRECTORY(jpeg)
  SET(JPEG_LIBRARY jpeg PARENT_SCOPE)
  SET(JPEG_LIBRARY jpeg CACHE STRING "Bundled JPEG" FORCE)
  LIST(APPEND JPEG_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/jpeg ${CMAKE_CURRENT_BINARY_DIR}/jpeg/config)
  SET(JPEG_INCLUDE_DIR ${JPEG_INCLUDE_DIRECTORIES} PARENT_SCOPE)
  SET(JPEG_INCLUDE_DIR ${JPEG_INCLUDE_DIRECTORIES} CACHE PATH "Bundled JPEG" FORCE)
ENDIF(NOT JPEG_FOUND)

IF (NOT PNG_FOUND)
  ADD_SUBDIRECTORY(zlib)
  ADD_SUBDIRECTORY(png)
  SET(PNG_LIBRARIES png zlib PARENT_SCOPE)
  LIST(APPEND PNG_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/png ${CMAKE_CURRENT_SOURCE_DIR}/zlib ${CMAKE_CURRENT_BINARY_DIR}/png/config)
  SET(PNG_INCLUDE_DIRS ${PNG_INCLUDE_DIRECTORIES} PARENT_SCOPE)
ENDIF (NOT PNG_FOUND)

# Add svg library
ADD_SUBDIRECTORY(vectorGraphics)

# Add ceres-solver (A Nonlinear Least Squares Minimizer)
ADD_SUBDIRECTORY(cxsparse)
ADD_SUBDIRECTORY(ceres-solver)

# Add an Approximate Nearest Neighbor library
ADD_SUBDIRECTORY(flann)

# Exif data parsing
ADD_SUBDIRECTORY(easyexif)

