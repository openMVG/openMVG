
ADD_LIBRARY(openMVG_kvld kvld.cpp kvld.h algorithm.cpp algorithm.h)
SET_PROPERTY(TARGET openMVG_kvld PROPERTY FOLDER OpenMVG/OpenMVG)
INSTALL(TARGETS openMVG_kvld DESTINATION lib/)