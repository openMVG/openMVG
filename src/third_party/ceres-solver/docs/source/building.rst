.. _chapter-building:

=======================
Building & Installation
=======================

Getting the source code
=======================
.. _section-source:

You can start with the `latest stable release
<http://ceres-solver.org/ceres-solver-1.10.0.tar.gz>`_ . Or if you want
the latest version, you can clone the git repository

.. code-block:: bash

       git clone https://ceres-solver.googlesource.com/ceres-solver

.. _section-dependencies:

Dependencies
============

Ceres relies on a number of open source libraries, some of which are
optional. For details on customizing the build process, see
:ref:`section-customizing` .

- `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_
  3.2.1 or later.  **Required**

  .. NOTE ::

    Ceres can also use Eigen as a sparse linear algebra
    library. Please see the documentation for ``-DEIGENSPARSE`` for`
    more details.

- `CMake <http://www.cmake.org>`_ 2.8.0 or later.
  **Required on all platforms except for Android.**

- `Google Log <http://code.google.com/p/google-glog>`_ 0.3.1 or
  later. **Recommended**

  .. NOTE::

    Ceres has a minimal replacement of ``glog`` called ``miniglog``
    that can be enabled with the ``MINIGLOG`` build
    option. ``miniglog`` is needed on Android as ``glog`` currently
    does not build using the NDK. It can however be used on other
    platforms too.

    **We do not advise using** ``miniglog`` **on platforms other than
    Android due to the various performance and functionality
    compromises in** ``miniglog``.

- `Google Flags <http://code.google.com/p/gflags>`_. Needed to build
  examples and tests.

- `SuiteSparse
  <http://www.cise.ufl.edu/research/sparse/SuiteSparse/>`_. Needed for
  solving large sparse linear systems. **Optional; strongly recomended
  for large scale bundle adjustment**

- `CXSparse <http://www.cise.ufl.edu/research/sparse/CXSparse/>`_.
  Similar to ``SuiteSparse`` but simpler and slower. CXSparse has
  no dependencies on ``LAPACK`` and ``BLAS``. This makes for a simpler
  build process and a smaller binary. **Optional**

- `BLAS <http://www.netlib.org/blas/>`_ and `LAPACK
  <http://www.netlib.org/lapack/>`_ routines are needed by
  ``SuiteSparse``, and optionally used by Ceres directly for some
  operations.

  On ``UNIX`` OSes other than Mac OS X we recommend `ATLAS
  <http://math-atlas.sourceforge.net/>`_, which includes ``BLAS`` and
  ``LAPACK`` routines. It is also possible to use `OpenBLAS
  <https://github.com/xianyi/OpenBLAS>`_ . However, one needs to be
  careful to `turn off the threading
  <https://github.com/xianyi/OpenBLAS/wiki/faq#wiki-multi-threaded>`_
  inside ``OpenBLAS`` as it conflicts with use of threads in Ceres.

  MAC OS X ships with an optimized ``LAPACK`` and ``BLAS``
  implementation as part of the ``Accelerate`` framework. The Ceres
  build system will automatically detect and use it.

  For Windows things are much more complicated. `LAPACK For
  Windows <http://icl.cs.utk.edu/lapack-for-windows/lapack/>`_
  has detailed instructions..

  **Optional but required for** ``SuiteSparse``.

.. _section-linux:

Linux
=====

We will use `Ubuntu <http://www.ubuntu.com>`_ as our example linux
distribution.

.. NOTE::

 Up to at least Ubuntu 13.10, the SuiteSparse package in the official
 package repository (built from SuiteSparse v3.4.0) **cannot** be used
 to build Ceres as a *shared* library.  Thus if you want to build
 Ceres as a shared library using SuiteSparse, you must perform a
 source install of SuiteSparse.  It is recommended that you use the
 current version of SuiteSparse (4.2.1 at the time of writing).


Start by installing all the dependencies.

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

We are now ready to build, test, and install Ceres.

.. code-block:: bash

 tar zxf ceres-solver-1.10.0.tar.gz
 mkdir ceres-bin
 cd ceres-bin
 cmake ../ceres-solver-1.10.0
 make -j3
 make test
 make install

You can also try running the command line bundling application with one of the
included problems, which comes from the University of Washington's BAL
dataset [Agarwal]_.

.. code-block:: bash

 bin/simple_bundle_adjuster ../ceres-solver-1.10.0/data/problem-16-22106-pre.txt

This runs Ceres for a maximum of 10 iterations using the
``DENSE_SCHUR`` linear solver. The output should look something like
this.

.. code-block:: bash

    iter      cost      cost_change  |gradient|   |step|    tr_ratio  tr_radius  ls_iter  iter_time  total_time
       0  4.185660e+06    0.00e+00    1.09e+08   0.00e+00   0.00e+00  1.00e+04       0    7.59e-02    3.37e-01
       1  1.062590e+05    4.08e+06    8.99e+06   5.36e+02   9.82e-01  3.00e+04       1    1.65e-01    5.03e-01
       2  4.992817e+04    5.63e+04    8.32e+06   3.19e+02   6.52e-01  3.09e+04       1    1.45e-01    6.48e-01
       3  1.899774e+04    3.09e+04    1.60e+06   1.24e+02   9.77e-01  9.26e+04       1    1.43e-01    7.92e-01
       4  1.808729e+04    9.10e+02    3.97e+05   6.39e+01   9.51e-01  2.78e+05       1    1.45e-01    9.36e-01
       5  1.803399e+04    5.33e+01    1.48e+04   1.23e+01   9.99e-01  8.33e+05       1    1.45e-01    1.08e+00
       6  1.803390e+04    9.02e-02    6.35e+01   8.00e-01   1.00e+00  2.50e+06       1    1.50e-01    1.23e+00

    Ceres Solver v1.10.0 Solve Report
    ----------------------------------
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
    Preprocessor                            0.261

      Residual evaluation                   0.082
      Jacobian evaluation                   0.412
      Linear solver                         0.442
    Minimizer                               1.051

    Postprocessor                           0.002
    Total                                   1.357

    Termination:                      CONVERGENCE (Function tolerance reached. |cost_change|/cost: 1.769766e-09 <= 1.000000e-06)

.. section-osx:

Mac OS X
========
.. NOTE::

 Ceres will not compile using Xcode 4.5.x (Clang version 4.1) due to a
 bug in that version of Clang.  If you are running Xcode 4.5.x, please
 update to Xcode >= 4.6.x before attempting to build Ceres.


On OS X, we recommend using the `homebrew
<http://mxcl.github.com/homebrew/>`_ package manager to install
Ceres. Assuming that you have the ``homebrew/science`` [#f1]_ tap
enabled, then

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

We are now ready to build, test, and install Ceres.

.. code-block:: bash

   tar zxf ceres-solver-1.10.0.tar.gz
   mkdir ceres-bin
   cd ceres-bin
   cmake ../ceres-solver-1.10.0
   make -j3
   make test
   make install

Like the Linux build, you should now be able to run
``bin/simple_bundle_adjuster``.


.. rubric:: Footnotes

.. [#f1] Ceres and many of its dependencies are in `homebrew/science
   <https://github.com/Homebrew/homebrew-science>`_ tap. So, if you
   don't have this tap enabled, then you will need to enable it as
   follows before executing any of the commands in this section.

   .. code-block:: bash

      brew tap homebrew/science


.. _section-windows:

Windows
=======

.. NOTE::

  If you find the following `CMake` difficult to set up, then you may
  be interested in a `Microsoft Visual Studio wrapper
  <https://github.com/tbennun/ceres-windows>`_ for Ceres Solver by Tal
  Ben-Nun.

On Windows, we support building with Visual Studio 2010 or newer. Note
that the Windows port is less featureful and less tested than the
Linux or Mac OS X versions due to the lack of an officially supported
way of building SuiteSparse and CXSparse.  There are however a number
of unofficial ways of building these libraries. Building on Windows
also a bit more involved since there is no automated way to install
dependencies.

.. NOTE:: Using ``google-glog`` & ``miniglog`` with windows.h.

 The windows.h header if used with GDI (Graphics Device Interface)
 defines ``ERROR``, which conflicts with the definition of ``ERROR``
 as a LogSeverity level in ``google-glog`` and ``miniglog``.  There
 are at least two possible fixes to this problem:

 #. Use ``google-glog`` and define ``GLOG_NO_ABBREVIATED_SEVERITIES``
    when building Ceres and your own project, as documented
    `here <http://google-glog.googlecode.com/svn/trunk/doc/glog.html>`__.
    Note that this fix will not work for ``miniglog``,
    but use of ``miniglog`` is strongly discouraged on any platform for which
    ``google-glog`` is available (which includes Windows).
 #. If you do not require GDI, then define ``NOGDI`` **before** including
    windows.h.  This solution should work for both ``google-glog`` and
    ``miniglog`` and is documented for ``google-glog``
    `here <https://code.google.com/p/google-glog/issues/detail?id=33>`__.

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

Android
=======

Download the ``Android NDK`` version ``r9d`` or later. Run
``ndk-build`` from inside the ``jni`` directory. Use the
``libceres.a`` that gets created.

.. _section-ios:

iOS
===

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

After building, you will get a ``libceres.a`` library, which you will need to
add to your Xcode project.

The default CMake configuration builds a bare bones version of Ceres
Solver that only depends on Eigen (``MINIGLOG`` is compiled into Ceres if it is
used), this should be sufficient for solving small to moderate sized problems
(No ``SPARSE_SCHUR``, ``SPARSE_NORMAL_CHOLESKY`` linear solvers and no
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

#. ``EIGENSPARSE [Default: OFF]``: By default, Ceres will not use
   Eigen's sparse Cholesky factorization. The is because this part of
   the code is licensed under the ``LGPL`` and since ``Eigen`` is a
   header only library, including this code will result in an ``LGPL``
   licensed version of Ceres.

   .. NOTE::

      For good performance, use Eigen version 3.2.2 or later.

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
   to disable multi-threading.

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

If the standard find package scripts for ``BLAS`` & ``LAPACK`` which ship with
``CMake`` fail to find the desired libraries on your system, try setting
``CMAKE_LIBRARY_PATH`` to the path(s) to the directories containing the
``BLAS`` & ``LAPACK`` libraries when invoking ``CMake`` to build Ceres via
``-D<VAR>=<VALUE>``.  This should result in the libraries being found for any
common variant of each.

If you are building on an exotic system, or setting ``CMAKE_LIBRARY_PATH``
does not work, or is not appropriate for some other reason, one option would be
to write your own custom versions of ``FindBLAS.cmake`` &
``FindLAPACK.cmake`` specific to your environment.  In this case you must set
``CMAKE_MODULE_PATH`` to the directory containing these custom scripts when
invoking ``CMake`` to build Ceres and they will be used in preference to the
default versions.  However, in order for this to work, your scripts must provide
the full set of variables provided by the default scripts.  Also, if you are
building Ceres with ``SuiteSparse``, the versions of ``BLAS`` & ``LAPACK``
used by ``SuiteSparse`` and Ceres should be the same.

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
the **PATHS** option to the ``FIND_PACKAGE()`` command, e.g.,

.. code-block:: cmake

   FIND_PACKAGE(Ceres REQUIRED PATHS "/some/where/local/")

Note that this can be used to have multiple versions of Ceres
installed.
