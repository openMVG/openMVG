.. _chapter-building:

=====================
Building Ceres Solver
=====================

Stable Ceres Solver releases are available for download at
`code.google.com <http://code.google.com/p/ceres-solver/>`_. For the
more adventurous, the git repository is hosted on `Gerrit
<https://ceres-solver-review.googlesource.com/>`_.

.. _section-dependencies:

Dependencies
============

Ceres relies on a number of open source libraries, some of which are
optional. For details on customizing the build process, see
:ref:`section-customizing` .

1. `CMake <http://www.cmake.org>`_ is a cross platform build
system. Ceres needs a relatively recent version of CMake (version
2.8.0 or better).

2. `eigen3 <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ is
used for doing all the low level matrix and linear algebra operations.

3. `google-glog <http://code.google.com/p/google-glog>`_ is
used for error checking and logging. Ceres needs glog version 0.3.1 or
later. Version 0.3 (which ships with Fedora 16) has a namespace bug
which prevents Ceres from building.

4. `gflags <http://code.google.com/p/gflags>`_ is a library for
processing command line flags. It is used by some of the examples and
tests. While it is not strictly necessary to build the library, we
strongly recommend building the library with gflags.


5. `SuiteSparse
<http://www.cise.ufl.edu/research/sparse/SuiteSparse/>`_ is used for
sparse matrix analysis, ordering and factorization. In particular
Ceres uses the AMD, CAMD, COLAMD and CHOLMOD libraries. This is an optional
dependency.

6. `CXSparse <http://www.cise.ufl.edu/research/sparse/CXSparse/>`_ is
a sparse matrix library similar in scope to ``SuiteSparse`` but with
no dependencies on ``LAPACK`` and ``BLAS``. This makes for a simpler
build process and a smaller binary.  The simplicity comes at a cost --
for all but the most trivial matrices, ``SuiteSparse`` is
significantly faster than ``CXSparse``.

7. `BLAS <http://www.netlib.org/blas/>`_ and `LAPACK
<http://www.netlib.org/lapack/>`_ routines are needed by
SuiteSparse. We recommend `ATLAS
<http://math-atlas.sourceforge.net/>`_, which includes BLAS and LAPACK
routines. It is also possible to use `OpenBLAS
<https://github.com/xianyi/OpenBLAS>`_ . However, one needs to be
careful to `turn off the threading
<https://github.com/xianyi/OpenBLAS/wiki/faq#wiki-multi-threaded>`_
inside ``OpenBLAS`` as it conflicts with use of threads in Ceres.

.. _section-linux:

Building on Linux
=================
We will use `Ubuntu <http://www.ubuntu.com>`_ as our example
platform. Start by installing all the dependencies.

.. code-block:: bash

     # CMake
     sudo apt-get install cmake
     # gflags
     tar -xvzf gflags-2.0.tar.gz
     cd gflags-2.0
     ./configure --prefix=/usr/local
     make
     sudo make install.
     # google-glog must be configured to use the previously installed gflags
     tar -xvzf glog-0.3.2.tar.gz
     cd glog-0.3.2
     ./configure --with-gflags=/usr/local/
     make
     sudo make install
     # BLAS & LAPACK
     sudo apt-get install libatlas-base-dev
     # Eigen3
     sudo apt-get install libeigen3-dev
     # SuiteSparse and CXSparse
     sudo apt-get install libsuitesparse-dev

We are now ready to build and test Ceres.

.. code-block:: bash

 tar zxf ceres-solver-1.7.0.tar.gz
 mkdir ceres-bin
 cd ceres-bin
 cmake ../ceres-solver-1.7.0
 make -j3
 make test

You can also try running the command line bundling application with one of the
included problems, which comes from the University of Washington's BAL
dataset [Agarwal]_.

.. code-block:: bash

 bin/simple_bundle_adjuster ../ceres-solver-1.7.0/data/problem-16-22106-pre.txt

This runs Ceres for a maximum of 10 iterations using the
``DENSE_SCHUR`` linear solver. The output should look something like
this.

.. code-block:: bash

    0: f: 4.185660e+06 d: 0.00e+00 g: 1.09e+08 h: 0.00e+00 rho: 0.00e+00 mu: 1.00e+04 li:  0 it: 1.16e-01 tt: 3.39e-01
    1: f: 1.062590e+05 d: 4.08e+06 g: 8.99e+06 h: 5.36e+02 rho: 9.82e-01 mu: 3.00e+04 li:  1 it: 3.90e-01 tt: 7.29e-01
    2: f: 4.992817e+04 d: 5.63e+04 g: 8.32e+06 h: 3.19e+02 rho: 6.52e-01 mu: 3.09e+04 li:  1 it: 3.52e-01 tt: 1.08e+00
    3: f: 1.899774e+04 d: 3.09e+04 g: 1.60e+06 h: 1.24e+02 rho: 9.77e-01 mu: 9.26e+04 li:  1 it: 3.60e-01 tt: 1.44e+00
    4: f: 1.808729e+04 d: 9.10e+02 g: 3.97e+05 h: 6.39e+01 rho: 9.51e-01 mu: 2.78e+05 li:  1 it: 3.62e-01 tt: 1.80e+00
    5: f: 1.803399e+04 d: 5.33e+01 g: 1.48e+04 h: 1.23e+01 rho: 9.99e-01 mu: 8.33e+05 li:  1 it: 3.54e-01 tt: 2.16e+00
    6: f: 1.803390e+04 d: 9.02e-02 g: 6.35e+01 h: 8.00e-01 rho: 1.00e+00 mu: 2.50e+06 li:  1 it: 3.59e-01 tt: 2.52e+00

 Ceres Solver Report
 -------------------
                                      Original                  Reduced
 Parameter blocks                        22122                    22122
 Parameters                              66462                    66462
 Residual blocks                         83718                    83718
 Residual                               167436                   167436
 Trust Region Strategy     LEVENBERG_MARQUARDT

                                         Given                     Used
 Linear solver                     DENSE_SCHUR              DENSE_SCHUR
 Preconditioner                            N/A                      N/A
 Threads:                                    1                        1
 Linear solver threads                       1                        1
 Linear solver ordering              AUTOMATIC                 22106,16

 Cost:
 Initial                          4.185660e+06
 Final                            1.803390e+04
 Change                           4.167626e+06

 Number of iterations:
 Successful                                  6
 Unsuccessful                                0
 Total                                       6

 Time (in seconds):
 Preprocessor                        2.229e-01

   Evaluator::Residuals              7.438e-02
   Evaluator::Jacobians              6.790e-01
   Linear Solver                     1.681e+00
 Minimizer                           2.547e+00

 Postprocessor                       1.920e-02
 Total                               2.823e+00

 Termination:               FUNCTION_TOLERANCE

.. section-osx:

Building on Mac OS X
====================

On OS X, we recommend using the `homebrew
<http://mxcl.github.com/homebrew/>`_ package manager to install the
dependencies. There is no need to install ``BLAS`` or ``LAPACK``
separately as OS X ships with optimized ``BLAS`` and ``LAPACK``
routines as part of the `vecLib
<https://developer.apple.com/library/mac/#documentation/Performance/Conceptual/vecLib/Reference/reference.html>`_
framework.

.. code-block:: bash

      # CMake
      brew install cmake
      # google-glog and gflags
      brew install glog
      # Eigen3
      brew install eigen
      # SuiteSparse and CXSparse
      brew install suite-sparse


We are now ready to build and test Ceres.

.. code-block:: bash

   tar zxf ceres-solver-1.7.0.tar.gz
   mkdir ceres-bin
   cd ceres-bin
   cmake ../ceres-solver-1.7.0
   make -j3
   make test


Like the Linux build, you should now be able to run
``bin/simple_bundle_adjuster``.

.. _section-windows:

Building on Windows with Visual Studio
======================================

On Windows, we support building with Visual Studio 2010 or newer. Note
that the Windows port is less featureful and less tested than the
Linux or Mac OS X versions due to the unavailability of SuiteSparse
and ``CXSparse``. Building is also more involved since there is no
automated way to install the dependencies.

#. Make a toplevel directory for deps & build & src somewhere: ``ceres/``
#. Get dependencies; unpack them as subdirectories in ``ceres/``
   (``ceres/eigen``, ``ceres/glog``, etc)

   #. ``Eigen`` 3.1 (needed on Windows; 3.0.x will not work). There is
      no need to build anything; just unpack the source tarball.

   #. ``google-glog`` Open up the Visual Studio solution and build it.
   #. ``gflags`` Open up the Visual Studio solution and build it.

#. Unpack the Ceres tarball into ``ceres``. For the tarball, you
   should get a directory inside ``ceres`` similar to
   ``ceres-solver-1.3.0``. Alternately, checkout Ceres via ``git`` to
   get ``ceres-solver.git`` inside ``ceres``.

#. Install ``CMake``,

#. Make a dir ``ceres/ceres-bin`` (for an out-of-tree build)

#. Run ``CMake``; select the ``ceres-solver-X.Y.Z`` or
   ``ceres-solver.git`` directory for the CMake file. Then select the
   ``ceres-bin`` for the build dir.

#. Try running ``Configure``. It won't work. It'll show a bunch of options.
   You'll need to set:

   #. ``GLOG_INCLUDE``
   #. ``GLOG_LIB``
   #. ``GFLAGS_LIB``
   #. ``GFLAGS_INCLUDE``

   to the appropriate place where you unpacked/built them.

#. You may have to tweak some more settings to generate a MSVC
   project.  After each adjustment, try pressing Configure & Generate
   until it generates successfully.

#. Open the solution and build it in MSVC


To run the tests, select the ``RUN_TESTS`` target and hit **Build
RUN_TESTS** from the build menu.

Like the Linux build, you should now be able to run ``bin/simple_bundle_adjuster``.

Notes:

#. The default build is Debug; consider switching it to release mode.
#. Currently ``system_test`` is not working properly.
#. Building Ceres as a DLL is not supported; patches welcome.
#. CMake puts the resulting test binaries in ``ceres-bin/examples/Debug``
   by default.
#. The solvers supported on Windows are ``DENSE_QR``, ``DENSE_SCHUR``,
   ``CGNR``, and ``ITERATIVE_SCHUR``.
#. We're looking for someone to work with upstream ``SuiteSparse`` to
   port their build system to something sane like ``CMake``, and get a
   supported Windows port.


.. _section-android:

Building on Android
===================


Download the ``Android NDK``. Run ``ndk-build`` from inside the
``jni`` directory. Use the ``libceres.a`` that gets created.

.. _section-customizing:

Customizing the build
=====================

It is possible to reduce the libraries needed to build Ceres and
customize the build process by passing appropriate flags to
``CMake``. Use these flags only if you really know what you are doing.

#. ``-DSUITESPARSE=OFF``: By default, Ceres will link to
   ``SuiteSparse`` if all its dependencies are present. Use this flag
   to build Ceres without ``SuiteSparse``. This will also disable
   dependency checking for ``LAPACK`` and ``BLAS``. This will reduce
   Ceres' dependencies down to ``Eigen``, ``gflags`` and
   ``google-glog``.

#. ``-DCXSPARSE=OFF``: By default, Ceres will link to ``CXSparse`` if
   all its dependencies are present. Use this flag to builds Ceres
   without ``CXSparse``. This will reduce Ceres' dependencies down to
   ``Eigen``, ``gflags`` and ``google-glog``.

#. ``-DGFLAGS=OFF``: Use this flag to build Ceres without
   ``gflags``. This will also prevent some of the example code from
   building.

#. ``-DSCHUR_SPECIALIZATIONS=OFF``: If you are concerned about binary
   size/compilation time over some small (10-20%) performance gains in
   the ``SPARSE_SCHUR`` solver, you can disable some of the template
   specializations by using this flag.

#. ``-DLINE_SEARCH_MINIMIZER=OFF``: The line search based minimizer is
   mostly suitable for large scale optimization problems, or when sparse
   linear algebra libraries are not available. You can further save on
   some compile time and binary size by using this flag.

#. ``-DOPENMP=OFF``: On certain platforms like Android,
   multi-threading with ``OpenMP`` is not supported. Use this flag to
   disable multithreading.

#. ``-DBUILD_DOCUMENTATION=ON``: Use this flag to enable building the
   documentation. In addition, ``make ceres_docs`` can be used to
   build only the documentation.

.. _section-using-ceres:

Using Ceres with CMake
======================

Once the library is installed with ``make install``, it is possible to
use CMake with `FIND_PACKAGE()
<http://www.cmake.org/cmake/help/v2.8.10/cmake.html#command:find_package>`_
in order to compile **user code** against Ceres. For example, for
`examples/helloworld.cc
<https://ceres-solver.googlesource.com/ceres-solver/+/master/examples/helloworld.cc>`_
the following CMakeList.txt can be used:

.. code-block:: cmake

    CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

    PROJECT(helloworld)

    FIND_PACKAGE(Ceres REQUIRED)
    INCLUDE_DIRECTORIES(${CERES_INCLUDES})

    # helloworld
    ADD_EXECUTABLE(helloworld helloworld.cc)
    TARGET_LINK_LIBRARIES(helloworld ${CERES_LIBRARIES})

Specify Ceres version
---------------------

Additionally, when CMake has found Ceres it can check the package
version, if it has been specified in the `FIND_PACKAGE()
<http://www.cmake.org/cmake/help/v2.8.10/cmake.html#command:find_package>`_
call.  For example:

.. code-block:: cmake

    FIND_PACKAGE(Ceres 1.2.3 REQUIRED)

The version is an optional argument.

Local installations
-------------------

If Ceres was installed in a non-standard path by specifying
-DCMAKE_INSTALL_PREFIX="/some/where/local", then the user should add
the **PATHS** option to the ``FIND_PACKAGE()`` command. e.g.,

.. code-block:: cmake

   FIND_PACKAGE(Ceres REQUIRED PATHS "/some/where/local/")

Note that this can be used to have multiple versions of Ceres installed.

Compiling against static or shared library
------------------------------------------

.. code-block:: cmake

    TARGET_LINK_LIBRARIES(helloworld ${CERES_LIBRARIES})

will result in a statically linked binary. Changing this line to

.. code-block:: cmake

    TARGET_LINK_LIBRARIES(helloworld ${CERES_LIBRARIES_SHARED})

will result in a dynamically linked binary.
