
IF(OpenMVG_USE_BOOST)
###
# Create the vocabulary tree
###
ADD_EXECUTABLE(openMVG_main_createVoctree createVoctree.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_createVoctree openMVG_voctree ${BOOST_LIBRARIES})
ENDIF(OpenMVG_USE_BOOST)