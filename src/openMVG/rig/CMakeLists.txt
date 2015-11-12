file(
  GLOB_RECURSE
  rig_files_header
  *.hpp
)
file(
  GLOB_RECURSE
  rig_files_cpp
  *.cpp
)
#file(GLOB_RECURSE REMOVEFILESUNITTEST *_test.cpp)
#
##Remove the unit test files (not been used by the library)
#list(REMOVE_ITEM localization_files_cpp ${REMOVEFILESUNITTEST})

ADD_LIBRARY(openMVG_rig
  ${rig_files_header}
  ${rig_files_cpp})
TARGET_LINK_LIBRARIES(openMVG_rig openMVG_localization)
SET_PROPERTY(TARGET openMVG_rig PROPERTY FOLDER OpenMVG)
INSTALL(TARGETS openMVG_rig DESTINATION lib EXPORT openMVG-targets)
