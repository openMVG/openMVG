.. _chapter-building:

=======================
Building & Installation
=======================

Getting the source code
=======================
.. _section-source:

You can start with the `latest stable release
<http://ceres-solver.org/ceres-solver-1.9.0.tar.gz>`_ . Or if you want
the latest version, you can clone the git repository

.. code-block:: bash

       git clone https://ceres-solver.googlesource.com/ceres-solver

.. _section-dependencies:

Dependencies
============

Ceres relies on a number of open source libraries, some of which are
optional. For details on customizing the build process, see
:ref:`section-customizing` .

- `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ 3.0 or later.
  **Required**

- `CMake <http://www.cmake.org>`_ 2.8.0 or later.
  **Required on all platforms except for Android.**

- `Google Log <http://code.google.com/p/google-glog>`_ 0.3.1 or
  later. **Recommended**

  Ceres has a minimal replacement of ``glog`` called ``miniglog``,
  enabled with the ``MINIGLOG`` build option. ``miniglog`` replaces
  the requirement for ``glog``. We advise using full ``glog`` due to
  performance compromises in ``miniglog``. ``miniglog`` is needed on
  Android.

- `Google Flags <http://code.google.com/p/gflags>`_. Needed to build
  examples and tests.

- `SuiteSparse
  <http://www.cise.ufl.edu/research/sparse/SuiteSparse/>`_. Needed for
  analyzing and solving sparse systems. Ceres useses the AMD, CAMD,
  COLAMD and CHOLMOD libraries.
  **Optional; strongly recomended for bundle adjustment**

- `CXSparse <http://www.cise.ufl.edu/research/sparse/CXSparse/>`_.
  Similar to ``SuiteSparse`` but simpler and slower. CXSparse has
  no dependencies on ``LAPACK`` and ``BLAS``. This makes for a simpler
  build process and a smaller binary. **Optional**

- `BLAS <http://www.netlib.org/blas/>`_ and `LAPACK
  <http://www.netlib.org/lapack/>`_ routines are needed by
  SuiteSparse, and optionally used by Ceres directly for some operations.
  We recommend `ATLAS <http://math-atlas.sourceforge.net/>`_,
  which includes BLAS and LAPACK routines. It is also possible to use
  `OpenBLAS <https://github.com/xianyi/OpenBLAS>`_ . However, one needs
  to be careful to `turn off the threading
  <https://github.com/xianyi/OpenBLAS/wiki/faq#wiki-multi-threaded>`_
  inside ``OpenBLAS`` as it conflicts with use of threads in Ceres.
  **Optional but required for SuiteSparse**

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

 tar zxf ceres-solver-1.9.0.tar.gz
 mkdir ceres-bin
 cd ceres-bin
 cmake ../ceres-solver-1.9.0
 make -j3
 make test

You can also try running the command line bundling application with one of the
included problems, which comes from the University of Washington's BAL
dataset [Agarwal]_.

.. code-block:: bash

 bin/simple_bundle_adjuster ../ceres-solver-1.9.0/data/problem-16-22106-pre.txt

This runs Ceres for a maximum of 10 iterations using the
``DENSE_SCHUR`` linear solver. The output should look something like
this.

.. code-block:: bash

   0: f: 4.185660e+06 d: 0.00e+00 g: 1.09e+08 h: 0.00e+00 rho: 0.00e+00 mu: 1.00e+04 li:  0 it: 8.73e-02 tt: 2.61e-01
   1: f: 1.062590e+05 d: 4.08e+06 g: 8.99e+06 h: 5.36e+02 rho: 9.82e-01 mu: 3.00e+04 li:  1 it: 1.85e-01 tt: 4.46e-01
   2: f: 4.992817e+04 d: 5.63e+04 g: 8.32e+06 h: 3.19e+02 rho: 6.52e-01 mu: 3.09e+04 li:  1 it: 1.74e-01 tt: 6.20e-01
   3: f: 1.899774e+04 d: 3.09e+04 g: 1.60e+06 h: 1.24e+02 rho: 9.77e-01 mu: 9.26e+04 li:  1 it: 1.74e-01 tt: 7.94e-01
   4: f: 1.808729e+04 d: 9.10e+02 g: 3.97e+05 h: 6.39e+01 rho: 9.51e-01 mu: 2.78e+05 li:  1 it: 1.73e-01 tt: 9.67e-01
   5: f: 1.803399e+04 d: 5.33e+01 g: 1.48e+04 h: 1.23e+01 rho: 9.99e-01 mu: 8.33e+05 li:  1 it: 1.75e-01 tt: 1.14e+00
   6: f: 1.803390e+04 d: 9.02e-02 g: 6.35e+01 h: 8.00e-01 rho: 1.00e+00 mu: 2.50e+06 li:  1 it: 1.75e-01 tt: 1.32e+00

   Ceres Solver Report
   -------------------
                                        Original                  Reduced
   Parameter blocks                        22122                    22122
   Parameters                              66462                    66462
   Residual blocks                         83718                    83718
   Residual                               167436                   167436

   Minimizer                        TRUST_REGION

   Dense linear algebra library            EIGEN
   Trust region strategy     LEVENBERG_MARQUARDT

                                           Given                     Used
   Linear solver                     DENSE_SCHUR              DENSE_SCHUR
   Threads                                     1                        1
   Linear solver threads                       1                        1
   Linear solver ordering              AUTOMATIC                22106, 16

   Cost:
   Initial                          4.185660e+06
   Final                            1.803390e+04
   Change                           4.167626e+06

   Minimizer iterations                        6
   Successful steps                            6
   Unsuccessful steps                          0

   Time (in seconds):
   Preprocessor                            0.173

     Residual evaluation                   0.115
     Jacobian evaluation                   0.498
     Linear solver                         0.517
   Minimizer                               1.242

   Postprocessor                           0.003
   Total                                   1.437

   Termination:                      CONVERGENCE (Function tolerance reached. |cost_change|/cost: 1.769750e-09 <= 1.000000e-06)

.. section-osx:

Building on Mac OS X
====================
.. NOTE::

 Ceres will not compile using Xcode 4.5.x (Clang version 4.1) due to a bug in that version of
 Clang.  If you are running Xcode 4.5.x, please update to Xcode >= 4.6.x before attempting to
 build Ceres.


On OS X, we recommend using the `homebrew
<http://mxcl.github.com/homebrew/>`_ package manager to install Ceres.

.. code-block:: bash

      brew install ceres-solver

will install the latest stable version along with all the required
dependencies and

.. code-block:: bash

      brew install ceres-solver --HEAD

will install the latest version in the git repo.

You can also install each of the dependencies by hand using `homebrew
<http://mxcl.github.com/homebrew/>`_. There is no need to install
``BLAS`` or ``LAPACK`` separately as OS X ships with optimized
``BLAS`` and ``LAPACK`` routines as part of the `vecLib
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

   tar zxf ceres-solver-1.9.0.tar.gz
   mkdir ceres-bin
   cd ceres-bin
   cmake ../ceres-solver-1.9.0
   make -j3
   make test

Like the Linux build, you should now be able to run
``bin/simple_bundle_adjuster``.

.. _section-windows:

Building on Windows with Visual Studio
======================================

On Windows, we support building with Visual Studio 2010 or newer. Note
that the Windows port is less featureful and less tested than the Linux or
Mac OS X versions due to the lack of an officially supported way of building
SuiteSparse and CXSparse.  There are however a number of unofficial ways of
building these libraries. Building on Windows also a bit more involved since
there is no automated way to install dependencies.

#. Make a toplevel directory for deps & build & src somewhere: ``ceres/``
#. Get dependencies; unpack them as subdirectories in ``ceres/``
   (``ceres/eigen``, ``ceres/glog``, etc)

   #. ``Eigen`` 3.1 (needed on Windows; 3.0.x will not work). There is
      no need to build anything; just unpack the source tarball.

   #. ``google-glog`` Open up the Visual Studio solution and build it.
   #. ``gflags`` Open up the Visual Studio solution and build it.

   #. (Experimental) ``SuiteSparse`` Previously SuiteSparse was not available
      on Windows, recently it has become possible to build it on Windows using
      the `suitesparse-metis-for-windows <https://github.com/jlblancoc/suitesparse-metis-for-windows>`_
      project.  If you wish to use ``SuiteSparse``, follow their instructions
      for obtaining and building it.

   #. (Experimental) ``CXSparse`` Previously CXSparse was not available on
      Windows, there are now several ports that enable it to be, including:
      `[1] <https://github.com/PetterS/CXSparse>`_ and
      `[2] <https://github.com/TheFrenchLeaf/CXSparse>`_.  If you wish to use
      ``CXSparse``, follow their instructions for obtaining and building it.

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

   #. ``EIGEN_INCLUDE_DIR_HINTS``
   #. ``GLOG_INCLUDE_DIR_HINTS``
   #. ``GLOG_LIBRARY_DIR_HINTS``
   #. ``GFLAGS_INCLUDE_DIR_HINTS``
   #. ``GFLAGS_LIBRARY_DIR_HINTS``
   #. (Optional) ``SUITESPARSE_INCLUDE_DIR_HINTS``
   #. (Optional) ``SUITESPARSE_LIBRARY_DIR_HINTS``
   #. (Optional) ``CXSPARSE_INCLUDE_DIR_HINTS``
   #. (Optional) ``CXSPARSE_LIBRARY_DIR_HINTS``

   to the appropriate directories where you unpacked/built them. If any of
   the variables are not visible in the ``CMake`` GUI, create a new entry
   for them.  We recommend using the ``<NAME>_(INCLUDE/LIBRARY)_DIR_HINTS``
   variables rather than setting the ``<NAME>_INCLUDE_DIR`` &
   ``<NAME>_LIBRARY`` variables directly to keep all of the validity
   checking, and to avoid having to specify the library files manually.

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
#. CMake puts the resulting test binaries in ``ceres-bin/examples/Debug``
   by default.
#. The solvers supported on Windows are ``DENSE_QR``, ``DENSE_SCHUR``,
   ``CGNR``, and ``ITERATIVE_SCHUR``.
#. We're looking for someone to work with upstream ``SuiteSparse`` to
   port their build system to something sane like ``CMake``, and get a
   fully supported Windows port.


.. _section-android:

Building on Android
===================

Download the ``Android NDK`` version ``r9d`` or later. Run
``ndk-build`` from inside the ``jni`` directory. Use the
``libceres.a`` that gets created.

.. _section-ios:

Building on iOS
===============
.. NOTE::

   You need iOS version 6.0 or higher to build Ceres Solver.

To build Ceres for iOS, we need to force ``CMake`` to find the toolchains from
the iOS SDK instead of using the standard ones. For example:

.. code-block:: bash

   cmake ../ceres-solver \
   -DCMAKE_TOOLCHAIN_FILE=../ceres-solver/cmake/iOS.cmake \
   -DEIGEN_INCLUDE_DIR=/path/to/eigen/header \
   -DIOS_PLATFORM=<PLATFORM>

``PLATFORM`` can be one of ``OS``, ``SIMULATOR`` and ``SIMULATOR64``. You can
build for ``OS`` (``armv7``, ``armv7s``, ``arm64``), ``SIMULATOR`` (``i386``) or
``SIMULATOR64`` (``x86_64``) separately and use ``LIPO`` to merge them into
one static library.  See ``cmake/iOS.cmake`` for more options.

After building, you will get ``libceres.a`` and ``libminiglog.a``
You need to add these two libraries into your XCode project.

The default CMake configuration builds a bare bones version of Ceres
Solver that only depends on Eigen and MINIGLOG, this should be
sufficient for solving small to moderate sized problems (No
``SPARSE_SCHUR``, ``SPARSE_NORMAL_CHOLESKY`` linear solvers and no
``CLUSTER_JACOBI`` and ``CLUSTER_TRIDIAGONAL`` preconditioners).

If you decide to use ``LAPACK`` and ``BLAS``, then you also need to add
``Accelerate.framework`` to your XCode project's linking dependency.

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

#. ``OPENMP [Default: ON]``: On certain platforms like Android,
   multi-threading with ``OpenMP`` is not supported. Turn this ``OFF``
   to disable multithreading.

#. ``BUILD_SHARED_LIBS [Default: OFF]``: By default Ceres is built as
   a static library, turn this ``ON`` to instead build Ceres as a
   shared library.

#. ``BUILD_DOCUMENTATION [Default: OFF]``: Use this to enable building
   the documentation, requires `Sphinx <http://sphinx-doc.org/>`_ and the
   `sphinx_rtd_theme <https://pypi.python.org/pypi/sphinx_rtd_theme>`_
   package available from the Python package index. In addition,
   ``make ceres_docs`` can be used to build only the documentation.

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

Building using custom BLAS & LAPACK installs
----------------------------------------------

If you are building on an exotic system, then the standard find package
scripts for ``BLAS`` & ``LAPACK`` which ship with ``CMake`` might not
work.  In this case, one option would be to write your own custom versions for
your environment and then set ``CMAKE_MODULE_PATH`` to the directory
containing these custom scripts when invoking ``CMake`` to build Ceres and they
will be used in preference to the default versions.  However, in order for this
to work, your scripts must provide the full set of variables provided by the
default scripts.  Also, if you are building Ceres with ``SuiteSparse``, the
versions of ``BLAS`` & ``LAPACK`` used by ``SuiteSparse`` and Ceres should be
the same.

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
