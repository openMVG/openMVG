
ADD_LIBRARY(openMVG_kvld kvld.cpp kvld.h algorithm.cpp algorithm.h)
TARGET_LINK_LIBRARIES(openMVG_kvld openMVG_image)
SET_PROPERTY(TARGET openMVG_kvld PROPERTY FOLDER OpenMVG/OpenMVG)

