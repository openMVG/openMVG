#---------------------
#Third party libraries
#---------------------

# libs should be static
set(BUILD_SHARED_LIBS OFF)

# use PIC code for link into shared lib
if(UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif(UNIX)

if (MSVC)
  # Disable warning about size_t conversions for all third_party libraries.
  add_compile_options(/wd4267)
endif()

# Unit test library
if (OpenMVG_BUILD_TESTS)
  add_subdirectory(CppUnitLite)
endif()

# Basic filesystem utils
set(STLPLUS_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/openMVG/third_party/stlplus3)
add_subdirectory(stlplus3)

# Add graph library
if(DEFINED OpenMVG_USE_INTERNAL_LEMON)
  add_subdirectory(lemon)
  set_property(TARGET lemon PROPERTY FOLDER OpenMVG/3rdParty/lemon)
  set_property(TARGET check PROPERTY FOLDER OpenMVG/3rdParty/lemon)
endif()

##
# Image I/O
##
if(NOT JPEG_FOUND)
  set(OpenMVG_USE_INTERNAL_JPEG ON PARENT_SCOPE)
  set(JPEG_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/openMVG/third_party/jpeg)
  add_subdirectory(jpeg)
  list(APPEND JPEG_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/jpeg ${CMAKE_CURRENT_BINARY_DIR}/jpeg/config)
  set(JPEG_INCLUDE_DIR ${JPEG_INCLUDE_DIRECTORIES})
  set(JPEG_LIBRARIES jpeg PARENT_SCOPE)
  set(JPEG_INCLUDE_DIR ${JPEG_INCLUDE_DIR} PARENT_SCOPE)
endif(NOT JPEG_FOUND)

# TIFF and PNG depend on zlib, if one of them is not found add the internal zlib
if(NOT PNG_FOUND OR NOT TIFF_FOUND)
  set(ZLIB_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/openMVG/third_party/zlib)
  add_subdirectory(zlib)
endif(NOT PNG_FOUND OR NOT TIFF_FOUND)

if (NOT PNG_FOUND)
  set(OpenMVG_USE_INTERNAL_PNG ON PARENT_SCOPE)
  set(PNG_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/openMVG/third_party/png)
  add_subdirectory(png)
  set(PNG_LIBRARIES png zlib)
  list(APPEND PNG_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/png ${CMAKE_CURRENT_SOURCE_DIR}/zlib ${CMAKE_CURRENT_BINARY_DIR}/png/config)
  set(PNG_INCLUDE_DIRS ${PNG_INCLUDE_DIRECTORIES})
  set(PNG_LIBRARIES ${PNG_LIBRARIES} PARENT_SCOPE)
  set(PNG_INCLUDE_DIRS ${PNG_INCLUDE_DIRS} PARENT_SCOPE)
endif (NOT PNG_FOUND)

if (NOT TIFF_FOUND)
  set(OpenMVG_USE_INTERNAL_TIFF ON PARENT_SCOPE)
  set(TIFF_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/openMVG/third_party/tiff)
  add_subdirectory(tiff)
  list(APPEND TIFF_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/tiff ${CMAKE_CURRENT_BINARY_DIR}/tiff)
  set(TIFF_INCLUDE_DIR ${TIFF_INCLUDE_DIRECTORIES})
  set(TIFF_LIBRARIES tiff)
  set(TIFF_LIBRARIES ${TIFF_LIBRARIES} PARENT_SCOPE)
  set(TIFF_INCLUDE_DIR ${TIFF_INCLUDE_DIR} PARENT_SCOPE)
endif (NOT TIFF_FOUND)
##
# End - Image I/O
##

# Add svg library
add_subdirectory(vectorGraphics)

# Add ceres-solver (A Nonlinear Least Squares Minimizer)
if (DEFINED OpenMVG_USE_INTERNAL_CERES)
  set(CXSPARSE_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/openMVG/third_party/cxsparse)
  add_subdirectory(cxsparse)
  set(CERES_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/openMVG/third_party/ceres)
  add_subdirectory(ceres-solver)
  set_property(TARGET cxsparse PROPERTY FOLDER OpenMVG/3rdParty/ceres)
  set_property(TARGET ceres PROPERTY FOLDER OpenMVG/3rdParty/ceres)
endif()

# Add an Approximate Nearest Neighbor library
if (DEFINED OpenMVG_USE_INTERNAL_FLANN)
  set(FLANN_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/openMVG/third_party/flann/src/cpp)
  add_subdirectory(flann)
endif()

# Exif data parsing
set(EASYEXIF_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/openMVG/third_party/easyexif)
add_subdirectory(easyexif)

# Fast corner detector
set(FAST_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/openMVG/third_party/fast)
add_subdirectory(fast)

##
# Install Header only libraries if necessary
##
if (DEFINED OpenMVG_USE_INTERNAL_EIGEN)
  #Configure Eigen install
  set(EIGEN_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/openMVG/third_party/eigen)
  add_subdirectory(eigen)
endif()

list(APPEND directories cmdLine histogram htmlDoc progress vectorGraphics)
foreach(inDirectory ${directories})
install(
  DIRECTORY ./${inDirectory}
  DESTINATION include/openMVG/third_party/
  COMPONENT headers
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)
endforeach(inDirectory)
