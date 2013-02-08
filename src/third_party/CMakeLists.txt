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
IF(WIN32 OR APPLE)
  ADD_SUBDIRECTORY(jpeg)
  ADD_SUBDIRECTORY(zlib)
  ADD_SUBDIRECTORY(png)
ENDIF(WIN32 OR APPLE)

# Add svg library
ADD_SUBDIRECTORY(vectorGraphics)

# Add ceres-solver (A Nonlinear Least Squares Minimizer)
ADD_SUBDIRECTORY(cxsparse)
ADD_SUBDIRECTORY(ceres-solver)

# Add an Approximate Nearest Neighbor library
ADD_SUBDIRECTORY(flann)

