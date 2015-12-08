
project(features)

set(features_files_headers
  regions.hpp
  akaze/mldb_descriptor.hpp
  akaze/msurf_descriptor.hpp
  akaze/AKAZE.hpp
  image_describer.hpp
  liop/liop_descriptor.hpp
  descriptor.hpp
  io_regions_type.hpp
  image_describer_akaze.hpp
  regions_factory.hpp
  feature.hpp
  features.hpp
  keypointSet.hpp)

set( features_files_sources
  akaze/AKAZE.cpp
  liop/liop_descriptor.cpp
  features_test.cpp)

if (CCTAG_FOUND)
  list(APPEND features_files_headers 
        cctag/CCTAG_describer.hpp 
        cctag/SIFT_CCTAG_describer.hpp)
  list(APPEND features_files_sources 
        cctag/SIFT_CCTAG_describer.cpp
        cctag/CCTAG_describer.cpp)
endif()

file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#Remove the future main files
list(REMOVE_ITEM features_files_sources ${REMOVEFILESUNITTEST})

set_source_files_properties(${features_files_sources} PROPERTIES LANGUAGE CXX)
ADD_LIBRARY(openMVG_features ${features_files_sources} ${features_files_headers})

target_link_libraries( openMVG_features openMVG_numeric openMVG_image vlsift stlplus)
if(OpenMVG_USE_CCTAG)
  target_link_libraries( openMVG_features ${CCTAG_LIBRARIES})
endif()

SET_TARGET_PROPERTIES(openMVG_features PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
INSTALL(TARGETS openMVG_features DESTINATION lib EXPORT openMVG-targets)
SET_PROPERTY(TARGET openMVG_features PROPERTY FOLDER OpenMVG/OpenMVG)

UNIT_TEST(openMVG features "openMVG_features;stlplus")

