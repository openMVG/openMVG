.. _chapter-building:

============
Installation
============

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
which prevents Ceres from building. Ceres contains a stripped-down,
minimal version of ``glog`` called ``miniglog``, which can be enabled
with the ``MINIGLOG`` build option. If enabled, it replaces the
requirement for ``glog``. However, in general it is recommended that
you use the full ``glog``.

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
significantly faster than ``CXSparse``. This is an optional dependency.

7. `BLAS <http://www.netlib.org/blas/>`_ and `LAPACK
<http://www.netlib.org/lapack/>`_ routines are needed by
SuiteSparse, and optionally used by Ceres directly for some operations.
We recommend `ATLAS <http://math-atlas.sourceforge.net/>`_,
which includes BLAS and LAPACK routines. It is also possible to use
`OpenBLAS <https://github.com/xianyi/OpenBLAS>`_ . However, one needs
to be careful to `turn off the threading
<https://github.com/xianyi/OpenBLAS/wiki/faq#wiki-multi-threaded>`_
inside ``OpenBLAS`` as it conflicts with use of threads in Ceres.

.. _section-linux:

Building on Linux
=================
We will use `Ubuntu <http://www.ubuntu.com>`_ as our example
platform. Start by installing all the dependencies.

.. NOTE::

 Up to at least Ubuntu 13.10, the SuiteSparse package in the official
 package repository (built from SuiteSparse v3.4.0) **cannot** be used
 to build Ceres as a *shared* library.  Thus if you want to build
 Ceres as a shared library using SuiteSparse, you must perform a
 source install of SuiteSparse.  It is recommended that you use the
 current version of SuiteSparse (4.2.1 at the time of writing).

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
     # SuiteSparse and CXSparse (optional)
     # - If you want to build Ceres as a *static* library (the default)
     #   you can use the SuiteSparse package in the main Ubuntu package
     #   repository:
     sudo apt-get install libsuitesparse-dev
     # - However, if you want to build Ceres as a *shared* library, you must
     #   perform a source install of SuiteSparse (and uninstall the Ubuntu
     #   package if it is currently installed.

We are now ready to build and test Ceres.

.. code-block:: bash

 tar zxf ceres-solver-1.8.0.tar.gz
 mkdir ceres-bin
 cd ceres-bin
 cmake ../ceres-solver-1.8.0
 make -j3
 make test

You can also try running the command line bundling application with one of the
included problems, which comes from the University of Washington's BAL
dataset [Agarwal]_.

.. code-block:: bash

 bin/simple_bundle_adjuster ../ceres-solver-1.8.0/data/problem-16-22106-pre.txt

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

.. NOTE::

 Ceres will not compile using Xcode 4.5.x (Clang version 4.1) due to a bug in that version of
 Clang.  If you are running Xcode 4.5.x, please update to Xcode >= 4.6.x before attempting to
 build Ceres.

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

   tar zxf ceres-solver-1.8.0.tar.gz
   mkdir ceres-bin
   cd ceres-bin
   cmake ../ceres-solver-1.8.0
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

   #. ``EIGEN_INCLUDE_DIR``
   #. ``GLOG_INCLUDE_DIR``
   #. ``GLOG_LIBRARY``
   #. ``GFLAGS_INCLUDE_DIR``
   #. ``GFLAGS_LIBRARY``

   to the appropriate place where you unpacked/built them. If any of the
   variables are not visible in the ``CMake`` GUI, toggle to the
   *Advanced View* with ``<t>``.

#. You may have to tweak some more settings to generate a MSVC
   project.  After each adjustment, try pressing Configure & Generate
   until it generates successfully.

#. Open the solution and build it in MSVC


To run the tests, select the ``RUN_TESTS`` target and hit **Build
RUN_TESTS** from the build menu.

Like the Linux build, you should now be able to run
``bin/simple_bundle_adjuster``.

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
customize the build process by setting the appropriate options in
``CMake``.  These options can either be set in the ``CMake`` GUI,
or via ``-D<OPTION>=<ON/OFF>`` when running ``CMake`` from the
command line.  In general, you should only modify these options from
their defaults if you know what you are doing.

.. NOTE::

 If you are setting variables via ``-D<VARIABLE>=<VALUE>`` when calling
 ``CMake``, it is important to understand that this forcibly **overwrites** the
 variable ``<VARIABLE>`` in the ``CMake`` cache at the start of *every configure*.

 This can lead to confusion if you are invoking the ``CMake``
 `curses <http://www.gnu.org/software/ncurses/ncurses.html>`_ terminal GUI
 (via ``ccmake``, e.g. ```ccmake -D<VARIABLE>=<VALUE> <PATH_TO_SRC>``).
 In this case, even if you change the value of ``<VARIABLE>`` in the ``CMake``
 GUI, your changes will be **overwritten** with the value passed via
 ``-D<VARIABLE>=<VALUE>`` (if one exists) at the start of each configure.

 As such, it is generally easier not to pass values to ``CMake`` via ``-D``
 and instead interactively experiment with their values in the ``CMake`` GUI.
 If they are not present in the *Standard View*, toggle to the *Advanced View*
 with ``<t>``.

Options controlling Ceres configuration
---------------------------------------

#. ``LAPACK [Default: ON]``: By default Ceres will use ``LAPACK`` (&
   ``BLAS``) if they are found.  Turn this ``OFF`` to build Ceres
   without ``LAPACK``. Turning this ``OFF`` also disables
   ``SUITESPARSE`` as it depends on ``LAPACK``.

#. ``SUITESPARSE [Default: ON]``: By default, Ceres will link to
   ``SuiteSparse`` if it and all of its dependencies are present. Turn
   this ``OFF`` to build Ceres without ``SuiteSparse``. Note that
   ``LAPACK`` must be ``ON`` in order to build with ``SuiteSparse``.

#. ``CXSPARSE [Default: ON]``: By default, Ceres will link to
   ``CXSparse`` if all its dependencies are present. Turn this ``OFF``
   to build Ceres without ``CXSparse``.

#. ``GFLAGS [Default: ON]``: Turn this ``OFF`` to build Ceres without
   ``gflags``. This will also prevent some of the example code from
   building.

#. ``MINIGLOG [Default: OFF]``: Ceres includes a stripped-down,
   minimal implementation of ``glog`` which can optionally be used as
   a substitute for ``glog``, thus removing ``glog`` as a required
   dependency. Turn this ``ON`` to use this minimal ``glog``
   implementation.

#. ``SCHUR_SPECIALIZATIONS [Default: ON]``: If you are concerned about
   binary size/compilation time over some small (10-20%) performance
   gains in the ``SPARSE_SCHUR`` solver, you can disable some of the
   template specializations by turning this ``OFF``.

#. ``LINE_SEARCH_MINIMIZER [Default: ON]``: The line search based
   minimizer is mostly suitable for large scale optimization problems,
   or when sparse linear algebra libraries are not available. You can
   further save on some compile time and binary size by turning this
   ``OFF``.

#. ``OPENMP [Default: ON]``: On certain platforms like Android,
   multi-threading with ``OpenMP`` is not supported. Turn this ``OFF``
   to disable multithreading.

#. ``BUILD_SHARED_LIBS [Default: OFF]``: By default Ceres is built as
   a static library, turn this ``ON`` to instead build Ceres as a
   shared library.

#. ``BUILD_DOCUMENTATION [Default: OFF]``: Use this to enable building
   the documentation, requires `Sphinx <http://sphinx-doc.org/>`_. In
   addition, ``make ceres_docs`` can be used to build only the
   documentation.

#. ``MSVC_USE_STATIC_CRT [Default: OFF]`` *Windows Only*: By default
   Ceres will use the Visual Studio default, *shared* C-Run Time (CRT) library.
   Turn this ``ON`` to use the *static* C-Run Time library instead.


Options controlling Ceres dependency locations
----------------------------------------------

Ceres uses the ``CMake``
`find_package <http://www.cmake.org/cmake/help/v2.8.12/cmake.html#command:find_package>`_
function to find all of its dependencies using
``Find<DEPENDENCY_NAME>.cmake`` scripts which are either included in Ceres
(for most dependencies) or are shipped as standard with ``CMake``
(for ``LAPACK`` & ``BLAS``).  These scripts will search all of the "standard"
install locations for various OSs for each dependency.  However, particularly
for Windows, they may fail to find the library, in this case you will have to
manually specify its installed location.  The ``Find<DEPENDENCY_NAME>.cmake``
scripts shipped with Ceres support two ways for you to do this:

#. Set the *hints* variables specifying the *directories* to search in
   preference, but in addition, to the search directories in the
   ``Find<DEPENDENCY_NAME>.cmake`` script:

   - ``<DEPENDENCY_NAME (CAPS)>_INCLUDE_DIR_HINTS``
   - ``<DEPENDENCY_NAME (CAPS)>_LIBRARY_DIR_HINTS``

   These variables should be set via ``-D<VAR>=<VALUE>``
   ``CMake`` arguments as they are not visible in the GUI.

#. Set the variables specifying the *explicit* include directory
   and library file to use:

   - ``<DEPENDENCY_NAME (CAPS)>_INCLUDE_DIR``
   - ``<DEPENDENCY_NAME (CAPS)>_LIBRARY``

   This bypasses *all* searching in the
   ``Find<DEPENDENCY_NAME>.cmake`` script, but validation is still
   performed.

   These variables are available to set in the ``CMake`` GUI. They
   are visible in the *Standard View* if the library has not been
   found (but the current Ceres configuration requires it), but
   are always visible in the *Advanced View*.  They can also be
   set directly via ``-D<VAR>=<VALUE>`` arguments to ``CMake``.

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
    INCLUDE_DIRECTORIES(${CERES_INCLUDE_DIRS})

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

Note that this can be used to have multiple versions of Ceres
installed.

