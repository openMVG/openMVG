.. _chapter-introduction:

============
Introduction
============

Solving nonlinear least squares problems [#f1]_ comes up in a broad
range of areas across science and engineering - from fitting curves in
statistics, to constructing 3D models from photographs in computer
vision. Ceres Solver [#f2]_ [#f3]_ is a portable C++ library for
solving non-linear least squares problems accurately and efficiently.

**Features**

#. A friendly :ref:`chapter-modeling` API.

#. Automatic and numeric differentiation.

#. Robust loss functions and local parameterizations.

#. Multithreading.

#. Trust-Region (Levenberg-Marquardt and Dogleg) and Line Search
   (Nonlinear CG and L-BFGS) solvers.

#. Variety of linear solvers.

   a. Dense QR and Cholesky factorization (using `Eigen
      <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_) for
      small problems.

   b. Sparse Cholesky factorization (using `SuiteSparse
      <http://www.cise.ufl.edu/research/sparse/SuiteSparse/>`_ and
      `CXSparse <http://www.cise.ufl.edu/research/sparse/CSparse/>`_) for
      large sparse problems.

   c. Specialized solvers for bundle adjustment problems in computer
      vision.

   d. Iterative linear solvers with preconditioners for general sparse
      and bundle adjustment problems.

#. Portable: Runs on Linux, Windows, Mac OS X and Android.


At Google, Ceres Solver has been used for solving a variety of
problems in computer vision and machine learning. e.g., it is used to
to estimate the pose of Street View cars, aircrafts, and satellites;
to build 3D models for PhotoTours; to estimate satellite image sensor
characteristics, and more.

`Blender <http://www.blender.org>`_ uses Ceres for `motion tracking
<http://mango.blender.org/development/planar-tracking-preview/>`_ and
`bundle adjustment
<http://wiki.blender.org/index.php/Dev:Ref/Release_Notes/2.67/Motion_Tracker>`_.


.. rubric:: Footnotes

.. [#f1] For a gentle but brief introduction to non-linear least
         squares problems, please start by reading the
         :ref:`chapter-tutorial`.

.. [#f2] While there is some debate as to who invented the method of
         Least Squares [Stigler]_, there is no debate that it was
         `Carl Friedrich Gauss
         <http://en.wikipedia.org/wiki/Carl_Friedrich_Gauss>`_ who
         brought it to the attention of the world. Using just 22
         observations of the newly discovered asteroid `Ceres
         <http://en.wikipedia.org/wiki/Ceres_(dwarf_planet)>`_, Gauss
         used the method of least squares to correctly predict when
         and where the asteroid will emerge from behind the Sun
         [TenenbaumDirector]_. We named our solver after Ceres to
         celebrate this seminal event in the history of astronomy,
         statistics and optimization.

.. [#f3] For brevity, in the rest of this document we will just use
         the term Ceres.



