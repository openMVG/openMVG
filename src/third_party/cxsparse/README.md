#[CXSparse](http://www.cise.ufl.edu/research/sparse/CXSparse/) Multi-platform compilation

CXSparse version 3.1.1.

##Compilation
- Unzip CXSparse
- mkdir CXSparseBuild
- cd CXSparseBuild
- Run cmake . ../CXSparse
- make

##Test
- cd Demo
- ./cs_demo3.exe < ../../CXSparse/Matrix/bcsstk16

 > residual must be very small (1.e-023)

#Modification:
- Disabled the complex support

#Authors:
- TheFrenchLeaf (cmake based build)

