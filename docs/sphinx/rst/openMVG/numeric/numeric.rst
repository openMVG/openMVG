*******************
numeric
*******************

This module provides math and linear algebra utils that relies on [Eigen]_ library.
Eigen is a C++ template library for linear algebra.

Basic idea is to provide to openMVG :

- a high level memory container for matrices and vectors,
- an easy matrices and vectors manipulation,
- a collection of numeric solvers and related algorithms.

Vector, Matrix containers
===========================

OpenMVG redefines some Eigen basis type (points, vectors, matrices) for code consitency and clarity:

* ``Vec2`` a single 2d point stored as a column matrix (x,y),
* ``Vec3`` a single 3d point stored as a column matrix (x,y,z),
* ``Vec2f, Vec3f`` float version.

|

* ``Vec`` a vector of value (double precision),
* ``Vecf`` a vector of floating point value,

|

* ``Mat`` the generic matrix container,
* ``Mat2X`` a collection of 2d points stored by column,
* ``Mat3X`` a collection of 3d points stored as column.

Note: Default memory alignment is column major.

.. code-block:: c++ 

  // Create a set of 2D points store as column
  Mat2X A(2, 5);
  A << 1, 2, 3, 4, 5,
       6, 7, 8, 9, 10;
  A.col(); // return a column vector : (1,6)^T
  A.row(); // return a row vector : (1,2,3,4,5)

Linear algebra 
===============

* SVD/QR/LU decomposition.
  
To know more
--------------

Please visit: http://eigen.tuxfamily.org/dox/group__QuickRefPage.html


