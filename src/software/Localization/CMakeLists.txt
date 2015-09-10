IF(OpenMVG_USE_BOOST)
###
# Create the vocabulary tree
###
ADD_EXECUTABLE(openMVG_main_localization main_localization.cpp)
INSTALL(TARGETS openMVG_main_localization DESTINATION bin/)
TARGET_LINK_LIBRARIES(openMVG_main_localization openMVG_localization openMVG_image openMVG_features vlsift ${BOOST_LIBRARIES})

###
# Create the vocabulary tree
###
ADD_EXECUTABLE(openMVG_main_createVoctree createVoctree.cpp)
INSTALL(TARGETS openMVG_main_createVoctree DESTINATION bin/)
TARGET_LINK_LIBRARIES(openMVG_main_createVoctree openMVG_voctree openMVG_features ${BOOST_LIBRARIES})

###
# generate the pair lists
###
ADD_EXECUTABLE(openMVG_main_generatePairList pairListFromVoctree.cpp)
INSTALL(TARGETS openMVG_main_generatePairList DESTINATION bin/)
TARGET_LINK_LIBRARIES(openMVG_main_generatePairList openMVG_voctree ${BOOST_LIBRARIES})


#@TODO decide where to put the binaries

# Installation rules
#SET_PROPERTY(TARGET openMVG_main_createVoctree PROPERTY FOLDER OpenMVG/software)
#INSTALL(TARGETS openMVG_main_createVoctree DESTINATION bin/)
#SET_PROPERTY(TARGET openMVG_main_generatePairList PROPERTY FOLDER OpenMVG/software)
#INSTALL(TARGETS openMVG_main_generatePairList DESTINATION bin/)

ENDIF(OpenMVG_USE_BOOST)
