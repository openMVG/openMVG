
set(localization_files_headers
  LocalizationResult.hpp
  VoctreeLocalizer.hpp
  optimization.hpp
  reconstructed_regions.hpp
  svgVisualization.hpp
  ILocalizer.hpp
)
set(localization_files_sources
  LocalizationResult.cpp
  VoctreeLocalizer.cpp
  optimization.cpp
  svgVisualization.cpp
)

if (CCTAG_FOUND)
  list(APPEND localization_files_headers 
        CCTagLocalizer.hpp)
  list(APPEND localization_files_sources 
        CCTagLocalizer.cpp)
endif()

ADD_LIBRARY(openMVG_localization
  ${localization_files_headers}
  ${localization_files_sources})
TARGET_LINK_LIBRARIES(openMVG_localization openMVG_sfm openMVG_voctree openMVG_numeric ${BOOST_LIBRARIES})
IF(OpenMVG_USE_CCTAG)
  TARGET_LINK_LIBRARIES(openMVG_localization ${CCTAG_LIBRARIES})
ENDIF(OpenMVG_USE_CCTAG)
SET_PROPERTY(TARGET openMVG_localization PROPERTY FOLDER OpenMVG)
INSTALL(TARGETS openMVG_localization DESTINATION lib EXPORT openMVG-targets)
