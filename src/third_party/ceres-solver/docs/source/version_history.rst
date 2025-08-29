.. default-domain:: cpp

.. highlight:: c++

.. cpp:namespace:: ceres


.. _chapter-version-history:

===============
Version History
===============

2.2.0
=====

New Features
------------

#. Substantial improvement to threading performance across the board
   (Dmitry Korchemkin)
#. Mixed precision solves + iterative refinement when using ``CUDA`` or
   CPU based dense linear solvers, or ``EIGEN_SPARSE`` as the sparse
   linear algebra library. (Sameer Agarwal & Joydeep Biswas)
#. Cuda based CGNR and preconditioner support (Joydeep Biswas & Sameer
   Agarwal)
#. Nested Dissection (``NESDIS``) is now supported as an ordering method
   in addition to ``AMD``. (Sameer Agarwal, Alex Stewart & Sergiu
   Deitsch)
#. **Power Bundle Adjustment** is available as a linear solver and as
   a preconditioner by the name of ``SCHUR POWER SERIES EXPANSION``
   (Mark Shachkov).
#. Generalized Euler Angle conversions (hs293go@)


Backward Incompatible API Changes
---------------------------------

#. :class:`LocalParameterization` has been removed, use
   :class:`Manifold` instead.
#. Ceres Solver now requires a C++17 compliant compiler.
#. Ceres Solver now requires CMake version 3.16 or later.
#. Ceres Solver now requires SuiteSparse version 4.5.6 or later.
#. OpenMP and NO_THREADING backends have been removed. C++ threads is
   how all threading is done.
#. Support for ``CX_SPARSE`` as a sparse linear algebra backend has
   been removed. Similar or better performance can be expected from
   ``Eigen`` as the sparse linear algebra library.


Bug Fixes & Minor Changes
-------------------------
#. Optimize the computation of the LM diagonal in TinySolver
#. Improvements to multi-threaded performance for small problems that
   had regressed due to changes to threading (Dmitrity Korchemkin)
#. Fix handling of M_PI for MSVC (Sergiu Deitsch)
#. Add a default value for Solver::Summary::linear_solver_ordering_type (Sameer Agarwal)
#. Make sure that the code compiles well with CUDA 11 (Dmitriy
   Korchemkin)
#. Rework MSVC warning suppression (Sergiu Deitsch)
#. Add an example for EvaluationCallback (Sameer Agarwal)
#. Add an example for IterationCallback (Sameer Agarwal)
#. Add end-to-end BA tests for SCHUR_POWER_SERIES_EXPANSION (Sameer Agarwal)
#. Update documentation for linear solvers (Sameer Agarwal)
#. Add an accessor for the CostFunctor in DynamicAutoDiffCostFunction (Sameer Agarwal)
#. Runtime check for cudaMallocAsync support (Dmitriy Korchemkin)
#. Remove cuda-memcheck based tests (Sameer Agarwal)
#. Modernize ``Sphinx`` related CMake handling as well the ``Sphinx``
   build process in the terminal. (Sergiu Deitsch)
#. Fix macos ``sprintf`` security related warnings (Sergiu Deitsch)
#. Lots of Cuda releated build system fixes (Sergiu Deitsch, Dmitriy
   Korchemkin, Jason Mak)
#. Improved windows build support (Sergiu Deitsch)
#. Various documentation fixes (Maxim Smolskiy, Evan Levine)
#. Improved handling of large Jacobians (Sameer Agarwal)
#. Improved handling of infinite initial cost (Sameer Agarwal)
#. Improved traits support for Jets (Sameer Agarwal)
#. Improved tests for Euler angle conversion routines (@Hs293Go)
#. Use a std::tuple to store ProductManifold for better efficiency
   (Sergiu Deitsch)
#. Allow default construction of ProductManifold when underlying
   manifolds have default constructors (Sergiu Deitsch)
#. Move LineManifold and SphereManifold into their own headers (Sameer
   Agarwal)
#. Fix a byte vs number of elements error when dealing with CUDA
   workspace computations (Joydeep Biswas)
#. Hide and prevent internal symbols from being exported (Sergiu
   Deitsch)
#. Switch to imported SuiteSparse, CXSparse & METIS targets.
#. Improve compilation on Ubuntu 20.04 (Sergiu Deitsch)
#. Update to using gtest 1.11.0 (Sameer Agarwal)
#. Fix Euler angle conversion code to not rely on constexpr
   constrctors for Jets. (Sameer Agarwal)
#. BlockRandomAccessSparseMatrix now uses a BlockSparseMatrix as
   storage instead of TripletSparseMatrix. (Dmitriy Korchemkin)
#. Deduction guide for DynamicAutoDiffCostFunction (Sergiu Deitsch)
#. Explicit conversions from long to ints (Alexander Ivanov)
#. Unused code deletion/commenting and code modernization (Alexander
   Ivanov)
#. Improve the bazel build & tests (Alexander Ivanov)
#. Fix a bug in QuaternionRotatePoint introduced by the use of hypot
   earlier in this release cycle (Jonathan Taylor & Sameer Agarwal)
#. Lots of GitHub CI improvements (Sergiu Deitsch & Dmitry Korchemkin)
#. Improve the robustness of the Cuda based dense linear algebra tests
   (Joydeep Biswas)
#. Refactor storage & threading support in BlockRandomAccessMatrix and
   its subclasses (Sameer Agarwal)
#. Fix a bug in CoordinateDescentMinimizer related to uninitialized
   variables (Sameer Agarwal)
#. Remove OpenMP and NO_THREADS backends. (Sameer Agarwal)
#. Fix version string parsing starting with SuiteSparse 6.0 (Sergiu
   Deitsch)
#. Use FindCUDAToolkit for CMake >= 3.17 (Alex Stewart)
#. Add a const accessor for the Problem::Options struct used by
   Problem. (Alex Stewart)
#. Fix a serious performance regression when using SuiteSparse
   introduced in `d09f7e9d5e
   <https://github.com/ceres-solver/ceres-solver/commit/d09f7e9d5e3bfab2d7ec7e81fd6a55786edca17a>`_. (Sameer
   Agarwal)
#. Fix the build on QNX (Alex Stewart)
#. Improve testing macros and documentation for Manifolds (Alex
   Stewart)
#. Improved code formatting (Tyler Hovanec)
#. Better use of std::unique_ptr in the code (Mike Vitus)
#. Fix a memory leak in ContextImpl (Sameer Agarwal)
#. Faster locking when num_thread = 1 (Sameer Agarwal)
#. Fix how x_norm is computed in TrustRegionMinimizer (Sameer Agarwal)
#. Faster JACOBI preconditioner for CGNR (Sameer Agarwal)
#. Convert internal enums to class enums (Sameer Agarwal)
#. Improve the code in small_blas to be more compiler friendly (Sameer
   Agarwal)
#. Add the ability to specify the pivot threshold in
   ::class::`Covariance::Options` (Sameer Agarwal)
#. Modernize the internals to use C++17 (Sameer Agarwal)
#. Choose SPMV algorithm based on the CUDA SDK Version (Joydeep
   Biswas)
#. Better defaults in ``bundle_adjuster.cc`` (Sameer Agarwal)
#. Use ``foo.data()`` instead of ``&foo[0]`` (Sameer Agarwal)
#. Fix GCC 12.1.1 LTO -Walloc-size-larger-than= warnings (Sergiu
   Deitsch)
#. Improved determinism in tests by re-using the same PRNG (Sergiu
   Deitsch)
#. Improved docs for ``vcpkg`` installation. (Sergiu Deitsch)
#. Update FindGlog.cmake to create glog::glog target (KrisThielemans@)
#. Improve consistency & correctness of Sphere & Line Manifolds
   (Julio L. Paneque)
#. Remove ``ceres/internal/random.h`` in favor of ``<random>``.
#. Fix a crash in ``InnerProductComputer`` (Sameer Agarwal)
#. Various fixes to improve compilation on windows using MinGW & MSVC
   (Sergiu Deitsch)
#. Fix fmin/fmax() to use Jet averaging on equality (Alex Stewart)
#. Fix use of conditional preprocessor checks within a macro in tests
   (Alex Stewart)
#. Better support for ``CUDA memcheck`` (Joydeep Biswas)
#. Improve the logic for linking to the platform specific threading
   library (Sergiu Deitsch)
#. Generate the version string at compile time (Sergiu Deitsch)
#. :class:`NumericDiffFirstOrderFunction` can now take a dynamically
   sized parameter vector. (Sameer Agarwal)
#. Fix compilation with SuiteSparse 7.2.0 (Mark Shackov)

2.1.0
=====

New Features
------------

#. Support for CUDA based dense solvers - ``DENSE_QR``,
   ``DENSE_NORMAL_CHOLESKY`` & ``DENSE_SCHUR`` (Joydeep Biswas, Sameer
   Agarwal)

#. :class:`Manifold` is the new
   :class:`LocalParameterization`. Version 2.1 is the transition
   release where users can use both :class:`LocalParameterization` as
   well as :class:`Manifold` objects as they transition from the
   former to the latter. :class:`LocalParameterization` will be
   removed in version 2.2. There should be no numerical change to the
   results as a result of this change. (Sameer Agarwal, Johannes Beck,
   Sergiu Deitsch)

#. A number of changes to :class:`Jet` s (Sergiu Deitsch)

   * :class:`Jet` gained support for, ``copysign``, ``fma`` (fused
     multiply-add), ``midpoint`` (C++20 and above), ``lerp`` (C++20
     and above), 3-argument ``hypot`` (C++17 and above), ``log10``,
     ``log1p``, ``exp1m``, ``norm`` (squared :math:`L^2` norm).

   * Quiet floating-point comparison: ``isless``, ``isgreater``,
     ``islessgreater``, ``islessequal``, ``isgreaterequal``,
     ``isunordered``, ``signbit``, ``fdim``

   * Categorization and comparison operations are applied exclusively
     and consistently to the scalar part of a Jet now: ``isnan``,
     ``isinf``, ``isnormal``, ``isfinite``, ``fpclassify`` (new),
     ``fmin``, ``fmax``

   * It is now possible to safely compare a :class:`Jet` against a scalar
     (or literal) without constructing a :class:`Jet` first (even if it's
     nested):

     .. code-block:: c++

        Jet<Jet<Jet<T, N>, M>, O> x;
        if (x == 2) { } // equivalent to x.a.a.a == 2


     This enables interaction with various arithmetic functions that
     expect a scalar like instance, such as ``boost::math::pow<-N>``
     for reciprocal computation.

#. Add :class:`NumericDiffFirstOrderFunction` (Sameer Agarwal)


Backward Incompatible API Changes
---------------------------------

#. :class:`LocalParameterization` is deprecated. It will be removed in
   version 2.2. Use :class:`Manifold` instead.
#. Classification functions like ``IsFinite`` are deprecated. Use the
   ``C++11`` functions (``isfinite``, ``isnan`` etc) going
   forward. However to maintain consistent behaviour with comparison
   operators, these functions only inspect the scalar part of the
   :class:`Jet`.

Bug Fixes & Minor Changes
-------------------------

#. Worked around an MSVC ordering bug when using C++17/20 (Sergiu
   Deitsch)
#. Added a CITATION.cff file. (Sergiu Deitsch)
#. Updated included gtest version to 1.11.0. This should fix some
   ``C++20`` compilation problems. (Sameer Agarwal).
#. Workaround ``MSVC`` ``STL`` deficiency in ``C++17`` mode (Sergiu
   Deitsch)
#. Fix ``Jet`` test failures on ``ARMv8`` with recent ``Xcode``
   (Sergiu Deitsch)
#. Fix unused arguments of ``Make1stOrderPerturbation`` (Dmitriy
   Korchemkin)
#. Fix ``SuiteSparse`` path and version reporting (Sergiu Deitsch)
#. Enable `GitHub` workflows and deprecate ``TravisCI`` (Sergiu
   Deitsch)
#. Add missing includes (Sergiu Deitsch, Sameer Agarwal)
#. Fix path for ``cuda-memcheck`` tests (Joydeep Biswas)
#. ClangFormat cleanup (Sameer Agarwal)
#. Set ``CMP0057`` policy for ``IN_LIST`` operator in
   ``FindSuiteSparse.cmake`` (Brent Yi)
#. Do not define unusable import targets (Sergiu Deitsch)
#. Fix Ubuntu 18.04 shared library build (Sergiu Deitsch)
#. Force ``C++`` linker when building the ``C`` API (Sergiu Deitsch)
#. Modernize the code to be inline with ``C++14`` (Sergiu Deitsch,
   Sameer Agarwal)
#. Lots of fixes to make Ceres compile out of the box on Windows
   (Sergiu Deitsch)
#. Standardize path handling using ``GNUImstallDirs`` (Sergiu Deitsch)
#. Add final specifier to classes to help the compiler with
   devirtualization (Sameer Agarwal)
#. LOTs of clean & modernization of the CMake build files (Sergiu
   Deitsch & Alex Stewart)
#. Simplification to the symbol export logic (Sergiu Deitsch)
#. Add cmake option ``ENABLE_BITCODE`` for iOS builds (John Harrison)
#. Add const accessor for functor wrapped by auto/numeric-diff objects
   (Alex Stewart)
#. Cleanup & refactor ``jet_test.cc``. (Sameer Agarwal)
#. Fix docs of supported sparse backends for mixed precision solvers
   (Alex Stewart)
#. Fix C++20 compilation (Sergiu Deitsch)
#. Add an example for ``BiCubicInterpolator`` (Dmitriy Korcchemkin)
#. Add a section to the documentation on implicit and inverse function
   theorems (Sameer Agarwal)
#. Add a note about Trigg's correction (Sameer Agarwal)
#. Fix the docs for ``Problem::RemoveResidualBlock`` &
   ``Problem::RemoveParameterBlock`` (Sameer Agarwal)
#. Fix an incorrect check in ``reorder_program.cc`` (William Gandler)
#. Add ``function_tolerance`` based convergence testing to ``TinySolver``
   (Sameer Agarwal).
#. Fix a number of typos in ``rotation.h`` (@yiping)
#. Fix a typo in ``interfacing_with_autodiff.rst`` (@tangobravo)
#. Fix a matrix sizing bug in covariance_impl.cc (William Gandler)
#. Fix a bug in ``system_test.cc`` (William Gandler)
#. Fix the Jacobian computation in ``trust_region_minimizer_test.cc``
   (William Gandler)
#. Fix a bug in ``local_parameterization_test.cc`` (William Gandler)
#. Add accessors to ``GradientProblem`` (Sameer Agarwal)
#. Refactor ``small_blas_gemm_benchmark`` (Ahmed Taei)
#. Refactor ``small_blas_test`` (Ahmed Taei)
#. Fix dependency check for building documentation (Sumit Dey)
#. Fix an errant double link in the docs (Timon Knigge)
#. Fix a typo in the version history (Noah Snavely)
#. Fix typo in LossFunctionWrapper sample code (Dmitriy Korchemkin)
#. Add fmax/fmin overloads for scalars (Alex Karatarakis)
#. Introduce benchmarks for ``Jet`` operations (Alexander Karatarakis)
#. Fix typos in documentation and fix the documentation for
   ``IterationSummary`` (Alexander Karatarakis)
#. Do not check MaxNumThreadsAvailable if the thread number is set
   to 1. (Fuhao Shi)
#. Add a macro ``CERES_GET_FLAG``. (Sameer Agarwal)
#. Reduce log spam in ``covariance_impl.cc`` (Daniel Henell)
#. Fix FindTBB version detection with TBB >= 2021.1.1 (Alex Stewart)
#. Fix Eigen3_VERSION (Florian Berchtold)
#. Allow Unity Build (Tobias Schluter)
#. Make miniglog's InitGoogleLogging argument const (Tobias Schluter)
#. Use portable expression for constant 2/sqrt(pi) (Tobias Schluter)
#. Fix a number of compile errors related (Austin Schuch)

   * ``format not a string literal``
   * ``-Wno-maybe-uninitialized error``
   * ``nonnull arg compared to NULL``
   * ``-Wno-format-nonliteral``
   * ``-Wmissing-field-initializers``
   * ``-Werror``

#. Fix ``cc_binary`` includes so examples build as an external repo
   (Austin Schuh)
#. Fix an explicit double in TinySolver (Bogdan Burlacu)
#. Fix unit quaternion rotation (Mykyta Kozlov)


2.0.0
=====

New Features
------------
#. Ceres Solver now requires a C++14 compatible compiler, Eigen
   version >= 3.3 & CMake version >= 3.5, XCode version >= 11.2 (Sameer
   Agarwal, Alex Stewart & Keir Mierle)
#. C++ threading based multi-threading support. (Mike Vitus)
#. :func:`Problem::AddResidualBlock`, :class:`SizedFunction`,
   :class:`AutoDiffCostFunction`, :class:`NumericDiffCostFunction`
   support an arbitrary number of parameter blocks using variadic
   templates (Johannes Beck)
#. On Apple platforms, support for Apple's Accelerate framework as a
   sparse linear algebra library. (Alex Stewart)
#. Significantly faster AutoDiff (Darius Rueckert)
#. Mixed precision solves when using
   ``SPARSE_NORMAL_CHOLESKY``. (Sameer Agarwal)
#. ``LocalParameterization`` objects can have a zero sized tangent
   size, which effectively makes the parameter block constant. In
   particular, this allows for a ``SubsetParameterization`` that holds
   all the coordinates of a parameter block constant. (Sameer Agarwal
   & Emil Ernerfeldt)
#. Visibility based preconditioning now works with ``Eigen`` and
   ``CXSparse``. (Sameer Agarwal)
#. Added :func:`Problem::EvaluateResidualBlock` and
   :func:`Problem::EvaluateResidualBlockAssumingParametersUnchanged`. (Sameer
   Agarwal)
#. ``GradientChecker`` now uses ``RIDDERS`` method for more accurate
   numerical derivatives. (Sameer Agarwal)
#. Covariance computation uses a faster SVD algorithm (Johannes Beck)
#. A new local parameterization for lines (Johannes Beck)
#. A new (``SUBSET``) preconditioner for problems with general
   sparsity. (Sameer Agarwal)
#. Faster Schur elimination using faster custom BLAS routines for
   small matrices. (yangfan)
#. Automatic differentiation for ``FirstOrderFunction`` in the form of
   :class:`AutoDiffFirstOrderFunction`. (Sameer Agarwal)
#. ``TinySolverAutoDiffFunction`` now supports dynamic number of residuals
   just like ``AutoDiffCostFunction``. (Johannes Graeter)

Backward Incompatible API Changes
---------------------------------

#. ``EvaluationCallback`` has been moved from ``Solver::Options`` to
   ``Problem::Options`` for a more correct API.
#. Removed ``Android.mk`` based build.
#. ``Solver::Options::num_linear_solver_threads`` is no more.

Bug Fixes & Minor Changes
-------------------------
#. Use CMAKE_PREFIX_PATH to pass Homebrew install location (Alex Stewart)
#. Add automatic differentiation support for ``Erf`` and ``Erfc``. (Morten Hennemose)
#. Add a move constructor to ``AutoDiffCostFunction``, ``NumericDiffCostFunction``, ``DynamicAutoDiffCostFunction`` and ``DynamicNumericDiffCostFunction``. (Julian Kent & Sameer Agarwal)
#. Fix potential for mismatched release/debug TBB libraries (Alex Stewart)
#. Trust region minimizer now reports the gradient of the current state, rather than zero when it encounters an unsuccessful step (Sameer Agarwal & Alex Stewart)
#. Unify symbol visibility configuration for all compilers (Taylor Braun-Jones)
#. Fix the Bazel build so that it points GitLab instead of the old BitBucket repo for Eigen (Sameer Agarwal)
#. Reformat source to be clang-format clean and add a script to format the repo using clang-format. (Nikolaus Demmel)
#. Various documentation improvements (Sameer Agarwal, Carl Dehlin,
   Bayes Nie, Chris Choi, Frank, Kuang Fangjun, Dmitriy Korchemkin,
   huangqinjin, Patrik Huber, Nikolaus Demmel, Lorenzo Lamia)
#. Huge number of build system simplification & cleanups (Alex
   Stewart, NeroBurner, Alastair Harrison, Linus Mårtensson, Nikolaus Demmel)
#. Intel TBB based threading removed (Mike Vitus)
#. Allow :class:`SubsetParameterization` to accept an empty vector of
   constant parameters. (Sameer Agarwal & Frédéric Devernay)
#. Fix a bug in DynamicAutoDiffCostFunction when all parameters are
   constant (Ky Waegel & Sameer Agarwal)
#. Fixed incorrect argument name in ``RotationMatrixToQuaternion``
   (Alex Stewart & Frank Dellaert)
#. Do not export class template LineParameterization (huangqinjin)
#. Change the type of parameter index/offset to match their getter/setter (huangqinjin)
#. Initialize integer variables with integer instead of double (huangqinjin)
#. Add std::numeric_limit specialization for Jets (Sameer Agarwal)
#. Fix a MSVC type deduction bug in ComputeHouseholderVector (Sameer Agarwal)
#. Allow LocalParameterizations to have zero local size. (Sameer Agarwal)
#. Add photometric and relative-pose residuals to autodiff benchmarks (Nikolaus Demmel)
#. Add a constant cost function to the autodiff benchmarks (Darius Rueckert)
#. Add const to GetCovarianceMatrix#. (Johannes Beck)
#. Fix Tukey loss function (Enrique Fernandez)
#. Fix 3+ nested Jet constructor (Julian Kent)
#. Fix windows MSVC build. (Johannes Beck)
#. Fix invert PSD matrix. (Johannes Beck)
#. Remove not used using declaration (Johannes Beck)
#. Let Problem::SetParameterization be called more than once. (Sameer Agarwal)
#. Make Problem movable. (Sameer Agarwal)
#. Make EventLogger more efficient. (Sameer Agarwal)
#. Remove a CHECK failure from covariance_impl.cc (Sameer Agarwal)
#. Add a missing cast in rotation.h (Sameer Agarwal)
#. Add a specialized SchurEliminator and integrate it for the case <2,3,6> (Sameer Agarwal)
#. Remove use of SetUsage as it creates compilation problems. (Sameer Agarwal)
#. Protect declarations of lapack functions under CERES_NO_LAPACK (Sameer Agarwal)
#. Drop ROS dependency on catkin (Scott K Logan)
#. Explicitly delete the copy constructor and copy assignment operator (huangqinjin)
#. Use selfAdjoingView<Upper> in InvertPSDMatrix. (Sameer Agarwal)
#. Speed up InvertPSDMatrix (Sameer Agarwal)
#. Allow Solver::Options::max_num_line_search_step_size_iterations = 0. (Sameer Agarwal)
#. Make LineSearchMinizer work correctly with negative valued functions. (Sameer Agarwal)
#. Fix missing declaration warnings in Ceres code (Sergey Sharybin)
#. Modernize ProductParameterization. (Johannes Beck)
#.  Add some missing string-to-enum-to-string convertors. (Sameer Agarwal)
#. Add checks in rotation.h for inplace operations. (Johannes Beck)
#. Update Bazel WORKSPACE for newest Bazel (Keir Mierle)
#. TripletSparseMatrix: guard against self-assignment (ngoclinhng)
#. Fix Eigen alignment issues. (Johannes Beck)
#. Add the missing <array> header to fixed_array.h (Sameer Agarwal)
#. Switch to FixedArray implementation from abseil. (Johannes Beck)
#. IdentityTransformation -> IdentityParameterization (Sameer Agarwal)
#. Reorder initializer list to make -Wreorder happy (Sam Hasinoff)
#. Reduce machoness of macro definition in cost_functor_to_function_test.cc (Sameer Agarwal)
#. Enable optional use of sanitizers (Alex Stewart)
#. Fix a typo in cubic_interpolation.h (Sameer Agarwal)
#. Update googletest/googlemock to db9b85e2. (Sameer Agarwal)
#. Fix Jacobian evaluation for constant parameter (Johannes Beck)
#. AutoDiffCostFunction: use static_assert to check if the correct overload of the constructor is used. (Christopher Wecht)
#. Avoid additional memory allocation in gradient checker (Justin Carpentier)
#. Swap the order of definition of IsValidParameterDimensionSequence. (Sameer Agarwal)
#. Add ParameterBlock::IsSetConstantByUser() (Sameer Agarwal)
#. Add parameter dims for variadic sized cost function (Johannes Beck)
#. Remove trailing zero parameter block sizes (Johannes Beck)
#. Adding integer sequence and algorithms (Johannes Beck)
#. Improve readability of LocalParameterization code. (Sameer Agarwal)
#. Simplifying Init in manual contructor (Johannes Beck)
#. Fix typo in NIST url. (Alessandro Gentilini)
#. Add a .clang-format file. (Sameer Agarwal)
#. Make ConditionedCostFunction compatible with repeated CostFunction. (Sameer Agarwal)
#. Remove conversions from a double to a Jet. (Kuang Fangjun)
#. close the file on return. (Kuang Fangjun)
#. Fix an error in the demo code for ceres::Jet. (Kuang Fangjun)
#. Recheck the residual after a new call. (Kuang Fangjun)
#. avoid recomputation. (Kuang Fangjun)
#. Fix calculation of Solver::Summary::num_threads_used. (Alex Stewart)
#. Convert calls to CHECK_NOTNULL to CHECK. (Sameer Agarwal)
#. Add a missing <cstdint> to block_structure.h (Sameer Agarwal)
#. Fix an uninitialized memory error in EvaluationCallbackTest (Sameer Agarwal)
#. Respect bounds when using Solver::Options::check_gradients (Sameer Agarwal)
#. Relax the limitation that SchurEliminator::Eliminate requires a rhs. (Sameer Agarwal)
#. Fix three out of bounds errors in CompressedRowSparseMatrix. (Sameer Agarwal)
#. Add Travis CI support. (Alex Stewart)
#. Refactor Ceres threading option configuration. (Alex Stewart)
#. Handle NULL permutation from SuiteSparseQR (Pau Gargallo)
#. Remove chunk shuffle in multithreaded SchurEliminator (Norbert Wenzel)
#. Add /bigobj to nist on MSVC. (Alex Stewart)
#. Fix 'xxx.cc has no symbols' warnings. (Alex Stewart)
#. Add a typedef to expose the scalar type used in a Jet. (Sameer Agarwal)
#. Fix a use after free bug in the tests. (Sameer Agarwal)
#. Simplify integration tests. (Sameer Agarwal)
#. Converts std::unique_lock to std::lock_guard. (Mike Vitus)
#. Bring the Bazel build in sync with the CMake build. (Sameer Agarwal)
#. Adds a ParallelFor wrapper for no threads and OpenMP. (Mike Vitus)
#. Improve the test coverage in small_blas_test (Sameer Agarwal)
#. Handle possible overflow in TrustRegionStepEvaluator. (Sameer Agarwal)
#. Fix lower-bound on result of minimising step-size polynomial. (Alex Stewart)
#. Adds missing functional include in thread_pool.h (Mike Vitus)


1.14.0
======

New Features
------------

#. New ``EvaluationCallback`` API. (Keir Mierle)
#. TBB based threading (Yury Prokazov & Mike Vitus)
#. C++11 threads based threading (Mike Vitus)
#. A ``ceres::Context`` object to cache and keep track of global
   state. (Mike Vitus)
#. TinySolver - A small dense solver meant for solving small problems
   really fast. [EXPERIMENTAL] (Keir Mierle & Sameer Agarwal)
#. Bazel Build. (Keir Mierle & Rodrigo Queiro)


Backward Incompatible API Changes
---------------------------------

#. ``Solver::Options::num_linear_solver_threads`` is deprecated,
   ``Solver::Options::num_threads`` controls all parallelism in Ceres
   Solver now. Similarly,
   ``Solver::Summary::num_linear_solver_threads_given`` and
   ``Solver::Summary::num_linear_solver_threads_used`` are also
   deprecated.


Bug Fixes & Minor Changes
-------------------------

#. Remove armv7 from target architectures when building for iOS >= 11. (Alex Stewart)
#. Corrects the documentation of Problem::AddResidualBlock. (Mike Vitus)
#. Fixes the configuration check in port.h. (Mike Vitus)
#. Add small_blas_gemm_benchmark. (Sameer Agarwal)
#. Implement some C++11 math functions for Jet (Emil Ernerfeldt)
#. Fix integer conversion warning in MSVC. (Alex Stewart)
#. Improve NDK build error handling (Keir Mierle)
#. Fix build: -Wreorder, test fail (Keir Mierle)
#. An implementation of SubsetPreconditioner. (Sameer Agarwal)
#. Split bundle adjustment tests into individual binaries (Keir Mierle)
#. Require Eigen >= 3.3.4 on aarch64. (Alex Stewart)
#. Fix TBB detection on Windows. (Alex Stewart)
#. Improve ExecutionSummary (Sameer Agarwal)
#. Remove as typo from callbacks.h (Sameer Agarwal)
#. Removes two unimplemented class functions. (Mike Vitus)
#. Update EigenTypes to deal with 1 column matrices (Sameer Agarwal)
#. Add GradientProblemSolver::Options::update_state_every_iteration (Sameer Agarwal)
#. Fixes the pose graph example documentation. (Mike Vitus)
#. Fix Eigen >= 3.3 compilation if EIGEN_DONT_VECTORIZE set (Janick Martinez Esturo)
#. Add an optional dependency on the Google Benchmark library. (Sameer Agarwal)
#. Fix the documentation for CostFunction::Evaluate. (Sameer Agarwal)
#. Fix a mathematical typo. (Sameer Agarwal)
#. Add TBB information to Ceres version string. (Alex Stewart)
#. Move discussion of dependency licensing to Sphinx docs. (Alex Stewart)
#. Fix an erroneous namespace comment (Sameer Agarwal)
#. Fix use of unnamed type as template argument warnings on Clang. (Alex Stewart)
#. Add link for CLA in docs; minor fixes (Keir Mierle)
#. Fix tiny_solver_test (Sameer Agarwal)
#. Improve compatibility with ceres::Solver (Sameer Agarwal)
#. Refactor nist.cc to be compatible with TinySolver (Sameer Agarwal)
#. Report timings with microsecond resolution (Thomas Gamper)
#. Add missing Eigen traits to Jets (Sameer Agarwal)
#. Use high-resolution timer on Windows (Thomas Gamper)
#. Add a comment about default constructed reference counts= (Keir Mierle)
#. Delete cost and loss functions when not in use. (Sameer Agarwal)
#. Fix assert_ndk_version for >= r11. (Alex Stewart)
#. Add docs explaining how to build Ceres with OpenMP on OS X. (Alex Stewart)
#. Update LAPACK option to refer to direct use by Ceres only. (Alex Stewart)
#. Hide optional SuiteSparse vars in CMake GUI by default. (Alex Stewart)
#. Always hide TBB_LIBRARY in CMake GUI by default. (Alex Stewart)
#. Fix typo in definition of f3 in powell example (x4 -> x3). (Alex Stewart)
#. Fix suppression of C++11 propagation warning. (Alex Stewart)
#. Add new Schur specialization for 2, 4, 6. (Chris Sweeney)
#. Use const keyword for 'int thread_id' variables. (pmoulon)


1.13.0
======

New Features
------------
#. ``LineSearchMinimizer`` and ``GradientProblemSolver`` are up to 2x
   faster due to fewer function evaluations. (Sameer Agarwal)
#. ``SPARSE_NORMAL_CHOLESKY`` is significantly faster because Ceres
   now computes the normal equations exploiting the static block
   sparsity structure. (Cheng Wang & Sameer Agarwal)
#. Add compound with scalar operators for Jets. (Alex Stewart)
#. Enable support for AVX instructions for Jets. (Alex Stewart)

Backward Incompatible API Changes
---------------------------------
The enum ``CovarianceAlgorithmType`` which controls the linear algebra
algorithm used to compute the covariance used to combine the choice of
the algorithm and the choice of the sparse linear algebra library into
the enum name. So we had ``SUITE_SPARSE_QR`` and
``EIGEN_SPARSE_QR``. ``Covariance::Options`` now has a separate member
allowing the user to choose the sparse linear algebra library, just
like the solver and ``CovarianceAlgorithmType`` now takes values
``DENSE_SVD`` and ``SPARSE_QR``. This is a forward looking change that
will allow us to develop more flexible covariance estimation
algorithms with multiple linear algebra backends.

Bug Fixes & Minor Changes
-------------------------
#. Fix ``InvertPSDMatrix`` as it was triggering an Eigen assert in
   Debug mode. (Philipp Hubner)
#. Fix cmake error from CeresConfig.cmake when Ceres not found (Taylor
   Braun-Jones)
#. Completely refactored ``SparseNormalCholeskySolver``. (Sameer
   Agarwal)
#. Fixed time reporting in ``Summary::FullReport`` when
   ``LineSearchMinimizer`` is used. (Sameer Agarwal)
#. Remove unused file: collections_port.cc. (Sameer Agarwal)
#. ``SPARSE_SCHUR`` + ``CX_SPARSE`` = Faster (Sameer Agarwal)
#. Refactored a number of linear solver tests to be more thorough and
   informative. (Sameer Agarwal)
#. Pass user-specified search hints as HINTS not PATHS. (Alex Stewart)
#. Prefer Eigen installs over exported build directories. (Alex
   Stewart)
#. Add OpenMP flags when compiling for C if enabled. (Alex Stewart)
#. Add a missing ``CERES_EXPORT`` to GradientChecker (Sameer Agarwal)
#. Use target_compile_features() to specify C++11 requirement if
   available. (Alex Stewart)
#. Update docs: .netrc --> .gitcookies (Keir Mierle)
#. Fix implicit precision loss warning on 64-bit archs (Ricardo
   Sanchez-Saez)
#. Optionally use exported Eigen CMake configuration if
   available. (Alex Stewart)
#. Use ``Ceres_[SOURCE/BINARY]_DIR`` not ``CMAKE_XXX_DIR`` to support
   nesting. (Alex Stewart)
#. Update ``Problem::EvaluateOptions`` documentation. (Sameer Agarwal)
#. Add public headers to CMake target for IDEs. (Devin Lane)
#. Add an article on interfacing with automatic
   differentiation. (Sameer Agarwal)
#. Add default Fedora/Debian locations for CXSparse to search
   paths. (Alex Stewart)
#. Add a test for ``LineSearchMinimizer`` (Sameer Agarwal)
#. Flatten the table of contents. (Sameer Agarwal)
#. Fix when ``LineSearchMinimizer`` adds the ``IterationSummary``` to
   ``Solver::Summary`` (Sameer Agarwal)
#. Fix search path for miniglog headers when Ceres is exported. (Alex
   Stewart)
#. Fix ambiguous reference to ``WARNING`` when using miniglog. (Alex
   Stewart)
#. Fix Jet/Eigen compatibility for Eigen > 3.3 (Julien Pilet)
#. Add max severity option when ``MINIGLOG`` is enabled (Taylor
   Braun-Jones)
#. Improvements to Schur template specializations (Sameer Agarwal)
#. Added an article on derivatives (Sameer Agarwal)
#. Require Eigen >= 3.3 to define ScalarBinaryOpTraits in Jet. (Alex
   Stewart)
#. A hacky fix for the Eigen::FullPivLU changes. (Sameer Agarwal)
#. Specify ``ScalarBinaryOpTraits`` for Jet types. (Chris Sweeney)
#. Remove spurious conversion from doubles to Jets. (Sameer Agarwal)
#. Fix an error in the tutorial code for ``NumericDiffCostFunction``
   (Sameer Agarwal)
#. ``CERES_EXPORT`` fix to compile Ceres as DLL (Je Hyeong Hong)
#. Fix detection of deprecated Bessel function names on MSVC. (Alex
   Stewart)
#. Ensure that partial evaluation of residuals triggers an error
   (Sameer Agarwal)
#. Fix detection of CMake-built glog on Windows. (Alex Stewart)
#. Add additional search paths for glog & Eigen on Windows. (Alex
   Stewart)
#. Various minor grammar and bug fixes to the documentation (Sameer
   Agarwal, Alex Stewart, William Rucklidge)


1.12.0
======

New Features
------------
#. Aligned ``Jet`` matrices for improved automatic differentiation
   performance. (Andrew Hunter)
#. Auto-differentiable implementations of Bessel functions, ``floor``,
   and ``ceil`` (Alessandro Gentilini & Michael Vitus)
#. New 2D and 3D SLAM examples. (Michael Vitus)
#. Added ``EigenQuaternionParameterization``. (Michael Vitus)
#. Added ``Problem::IsParameterBlockConstant`` (Thomas Schneider)
#. A complete refactoring of ``TrustRegionMinimizer``. (Sameer Agarwal)
#. Gradient checking cleanup and local parameterization bugfix (David
   Gossow)


Backward Incompatible API Changes
---------------------------------
#. ``Solver::Options::numeric_derivative_relative_step_size`` has been
   renamed to
   ``Solver::Options::gradient_check_numeric_derivative_relative_step_size``. (Sameer
   Agarwal)

Bug Fixes & Minor Changes
-------------------------
#. Clear XXX_FOUND in Find<XXX>.cmake prior to searching. (Alex
   Stewart)
#. Fix versioning in the documentation (Sameer Agarwal)
#. Fix missing gflags imported target definition in
   CeresConfig.cmake. (Alex Stewart)
#. Make gflags a public dependency of Ceres if it and glog are
   found. (Alex Stewart)
#. Add support for glog exported CMake target. (Alex Stewart)
#. Use ``google::GLOG_WARNING`` instead of ``WARNING`` in tests to
   support MSVC. (Alex Stewart)
#. Update gtest and gmock to
   ``a2b8a8e07628e5fd60644b6dd99c1b5e7d7f1f47`` (Sameer Agarwal)
#. Add MSVC-specific ``#define`` to expose math constants in
   ``<cmath>``. (Alex Stewart)
#. Fix typo. indepdendent -> independent (Hung Lun)
#. Fix potential invalid reset of CMAKE_FIND_LIBRARY_PREFIXES on MSVC
   (Alex Stewart)
#. Fix use of alignas(0) which is not ignored on GCC (Alex Stewart)
#. Use default alignment if alignof(std::max_align_t) < 16 with C++11
   (Alex Stewart)
#. Introduce a common base class for DynamicAutoDiffCostFunction and
   DynamicNumericDiffCostFunction. (Sameer Agarwal)
#. Fix an exact equality test causing breakage in
   gradient_checker_test. (Sameer Agarwal)
#. Add GradientProblemSolver::Options::parameter_tolerance. (Sameer
   Agarwal)
#. Add missing T() wrappers for constants. (Rob Carroll)
#. Remove two checks from rotation.h (Sameer Agarwal)
#. Relax the tolerance in QuaternionParameterizationTestHelper. (Je
   Hyeong Hong)
#. Occured -> Occurred. (Sameer Agarwal)
#. Fix a test error in autodiff_test.cc. (Je Hyeong Hong)
#. Fix documentation source for templated function in ``rotation.h``.
#. Add ``package.xml`` to enable Catkin builds. (Damon Kohler)
#. Relaxing Jacobian matching in Gradient Checker test. (David Gossow)
#. Allow SubsetParameterization to hold all parameters constant
   (Sameer Agarwal)
#. Fix an Intel compiler error in covariance_impl.cc (Je Hyeong Hong)
#. Removing duplicate include directive. (David Gossow)
#. Remove two DCHECKs from CubicHermiteSpline. (Sameer Agarwal)
#. Fix some compiler warnings. (Richard Trieu)
#. Update ExpectArraysClose to use ExpectClose instead of
   EXPECT_NEAR. (Phillip Hubner)
#. FindWithDefault returns by value rather than reference. (@aradval)
#. Fix compiler errors on some systems. (David Gossow)
#. Note that Problem::Evaluate cannot be called from an
   IterationCallback. (Sameer Agarwal)
#. Use ProductParameterization in bundle_adjuster.cc (Sameer Agarwal)
#. Enable support for OpenMP in Clang if detected. (Alex Stewart)
#. Remove duplicate entry for the NIST example in the docs. (Michael
   Vitus)
#. Add additional logging for analyzing orderings (Sameer Agarwal)
#. Add readme for the sampled_function example. (Michael Vitus)
#. Use _j[0,1,n]() Bessel functions on MSVC to avoid deprecation
   errors. (Alex Stewart & Kichang Kim)
#. Fix: Copy minimizer option ``is_silent`` to
   ``LineSearchDirection::Options`` (Nicolai Wojke)
#. Fix typos in ``users.rst`` (Sameer Agarwal)
#. Make some Jet comparisons exact. (Sameer Agarwal)
#. Add colmap to users.rst (Sameer Agarwal)
#. Fix step norm evaluation in LineSearchMinimizer (Sameer Agarwal)
#. Remove use of -Werror when compiling Ceres. (Alex Stewart)
#. Report Ceres compile options as components in find_package(). (Alex
   Stewart)
#. Fix a spelling error in nnls_modeling.rst (Timer)
#. Only use collapse() directive with OpenMP 3.0 or higher. (Keir
   Mierle)
#. Fix install path for CeresConfig.cmake to be architecture-aware.
#. Fix double conversion to degrees in rotation_test (Keir Mierle)
#. Make Jet string output more readable (Keir Mierle)
#. Fix rotation_test IsClose() and related tests (Keir Mierle)
#. Loosen an exact equality in local_parameterization_test (Sameer
   Agarwal)
#. make_docs: Pass the file encoding to open() (Niels Ole Salscheider)
#. Fix error message returned when using SUITE_SPARSE_QR in covariance
   estimation on a ceres built without SuiteSparse support. (Simon
   Rutishauser)
#. Fix CXX11 option to be available on MinGW & CygWin, but not
   MSVC. (Alex Stewart)
#. Fix missing early return() in xxx_not_found() dependency
   macros. (Alex Stewart)
#. Initialize ``inner_iterations_were_useful_`` correctly. (Sameer
   Agarwal)
#. Add an implementation for GradientProblemSolver::Options::IsValid
   (Sameer Agarwal)
#. Fix use of va_copy() if compiling with explicit C++ version <
   C++11. (Alex Stewart)
#. Install CMake files to lib/cmake/Ceres (Niels Ole Salscheider)
#. Allow users to override the documentation install directory. (Niels
   Ole Salscheider)
#. Add covariance matrix for a vector of parameters (Wannes Van Loock)
#. Saner tolerances & stricter LRE test. (Sameer Agarwal)
#. Fix a malformed sentence in the tutorial. (Sameer Agarwal)
#. Add logging for sparse Cholesky factorization using Eigen. (Sameer
   Agarwal)
#. Use std::adjacent_find instead of std::unique. (Sameer Agarwal)
#. Improve logging in CompressedRowJacobianWriter on crash. (Sameer
   Agarwal)
#. Fix free parameter block handling in covariance computation (Wannes
   Van Loock)
#. Report the number of line search steps in FullReport. (Sameer
   Agarwal)
#. Make CMake read Ceres version directly from
   include/ceres/version.h. (Alex Stewart)
#. Lots of code style/lint changes. (William Rucklidge)
#. Fix covariance computation for constant blocks (Wannes Van Loock)
#. Add IOS_DEPLOYMENT_TARGET variable to iOS.cmake (Eduard Feicho)
#. Make miniglog threadsafe on non-windows system by using
   localtime_r() instead of localtime() for time formatting (Simon
   Rutishauser)

1.11.0
======

New Features
------------
#. Adaptive numeric differentiation using Ridders' method. (Tal
   Ben-Nun)
#. Add ``CubicInterpolator`` and ``BiCubicInterpolator`` to allow
   smooth interpolation of sampled functions and integration with
   automatic differentiation.
#. Add method to return covariance in tangent space. (Michael Vitus &
   Steve Hsu)
#. Add Homogeneous vector parameterization. (Michael Vitus)
#. Add a ``ProductParameterization``, a local parameterization that
   can be constructed as a cartesian product of other local
   parameterization.
#. Add DynamicCostFunctionToFunctor. (David Gossow)
#. Optionally export Ceres build directory into local CMake package
   registry.
#. Faster ``SPARSE_NORMAL_CHOLESKY`` in the presence of dynamic
   sparsity.

Bug Fixes & Minor Changes
-------------------------
#. Remove use of link-time optimisation (LTO) for all compilers due to
   portability issues with gtest / type_info::operator== & Eigen with
   Clang on OS X vs GCC 4.9+ on Linux requiring contradictory 'fixes'.
#. Use link-time optimisation (LTO) only when compiling Ceres itself,
   not tests or examples, to bypass gtest / type_info::operator==
   issue.
#. Use old minimum iOS version flags on Xcode < 7.0.
#. Add gtest-specific flags when building/using as a shared library.
#. Clean up iOS.cmake to use xcrun/xcodebuild & libtool.
#. Import the latest version of ``googletest``.
#. Refactored ``system_test`` into ``bundle_adjustment_test`` and
   ``system_test``, where each test case is its own test.
#. Fix invalid memory access bug in
   ``CompressedRowSparseMatrix::AppendRows`` when it was called with a
   matrix of size zero.
#. Build position independent code when compiling Ceres statically
   (Alexander Alekhin).
#. Fix a bug in DetectStructure (Johannes Schonberger).
#. Reduce memory footprint of SubsetParameterization (Johannes
   Schonberger).
#. Fix for reorder program unit test when built without suitesparse
   (Sergey Sharybin).
#. Fix a bug in the Schur eliminator (Werner Trobin).
#. Fix a bug in the reordering code (Bernhard Zeisl).
#. Add missing CERES_EXPORT to ComposedLoss (Simon Rutishauser).
#. Add the option to use numeric differentiation to ``nist`` and
   ``more_garbow_hillstrom``.
#. Fix EIGENSPARSE option help s/t it displays in CMake ncurses GUI.
#. Fix SparseNormalCholeskySolver with dynamic sparsity (Richie
   Stebbing).
#. Remove legacy dependency detection macros.
#. Fix failed if() condition expansion if gflags is not found.
#. Update all CMake to lowercase function name style.
#. Update minimum iOS version to 7.0 for shared_ptr/unordered_map.
#. Fix bug in gflags' <= 2.1.2 exported CMake configuration.
#. Remove the spec file needed for generating RPMs.
#. Fix a typo in small_blas.h (Werber Trobin).
#. Cleanup FindGflags & use installed gflags CMake config if present.
#. Add default glog install location on Windows to search paths
   (bvanevery).
#. Add default Eigen install location on Windows to search paths
   (bvanevery).
#. Fix explanation of config.h generation in bare config.h.
#. Fix unused parameter compiler warnings in numeric_diff.h.
#. Increase tolerance for a test in polynomial_test (Taylor Braun
   Jones).
#. Fix addition of Gerrit commit hook when Ceres is a git submodule
   (Chris Cooper).
#. Fix missing EIGEN_VERSION expansion typo.
#. Fix links to SuiteSparse & CXSparse (Henrique Mendonça).
#. Ensure Eigen is at least 3.1.0 for Eigen/SparseCore.
#. Add option to use C++11 (not TR1) shared_ptr & unordered_map
   (Norman Goldstein).
#. Fix an incorrect usage message in bundle_adjuster.cc
#. Gracefully disable docs if Sphinx is not found.
#. Explicitly use (new) default OS X rpath policy if present.
#. Add support of EIGEN_SPARSE type in
   IsSparseLinearAlgebraLibraryTypeAvailable function (Pierre Moulon).
#. Allow the LossFunction contained in a LossFunctionWrapper to be
   NULL. This is consistent with how NULL LossFunctions are treated
   everywhere else. (Simon Rutishauser).
#. Improve numeric differentation near zero.
#. Refactored DynamicNumericDiffCostFunction to use NumericDiff (Tal
   Ben-Nun).
#. Remove use of :caption tag in Sphinx.
#. Add a small test to make sure GradientProblemSolver works correctly
   (Petter Strandmark).
#. Add simple unit tests for GradientProblem (Petter Strandmark).
#. Make the robust curve fitting example robust.
#. Homogenize convergence operators in docs and code (Johannes
   Schonberger).
#. Add parameter_tolerance convergence to line search minimizer
   (Johannes Schonberger).
#. Fix bug where pow(JetA,JetB) returned wrong result for JetA==0
   (Russell Smith).
#. Remove duplicate step norm computation (Johannes Schonberger).
#. Enhance usability when encountering Eigen version mismatches
   (Andrew Hundt).
#. Add PLY file logger before and after BA in order to ease visual
   comparison (Pierre Moulon).
#. Fix CMake config file docs to include 2.8.x & 3.x styles.
#. Python3 fixes (Markus Moll).
#. Remove confusing code from DenseJacobianWriter (Michael Vitus).
#. Add documentation on CMake package installation process.
#. Revert a call to SolveUpperTriangularUsingCholesky.
#. Make CERES_EIGEN_VERSION macro independent of CMake.
#. Add versions of dependencies used to FullReport().
#. Ensure local config.h is used if Ceres is already installed.
#. Small messaging and comment updates in CMake
#. Handle possible presence of library prefixes in MSVC (Sylvain
   Duchêne).
#. Use -O2 not -O3 on MinGW to workaround issue with Eigen
   (s1m3mu3@gmail.com).
#. Increase tolerance in small_blas test for Cygwin
   (s1m3mu3@gmail.com).
#. Fix iOS cmake file for cmake 3.0 (Jack Feng)
#. Fix missing gflags shlwapi dependency on MinGW (s1m3mu3@gmail.com).
#. Add thread dependency & fix namespace detection on Windows for
   gflags (arrigo.benedetti@gmail.com).
#. Rename macros in the public API to have a ``CERES_`` prefix.
#. Fix ``OrderedGroup::Reverse()`` when it is empty (Chris Sweeney).
#. Update the code to point to ceres-solver.org.
#. Update documentation to point to the GitHub issue tracker.
#. Disable ``LAPACK`` for iOS builds. (Greg Coombe)
#. Force use of single-thread in ``Problem::Evaluate()`` without
   OpenMP.
#. Less strict check for multithreading. (Chris Sweeney)
#. Update tolerances in small_blas_test.cc (Philipp Hubner)
#. Documentation corrections (Steve Hsu)
#. Fixed ``sampled_function.cc`` (Pablo Speciale)
#. Fix example code in the documentation. (Rodney Hoskinson)
#. Improve the error handling in Conjugate Gradients.
#. Improve preconditioner documentation.
#. Remove dead code from fpclassify.h.
#. Make Android.mk threads sensitive.
#. Changed the ``CURRENT_CONFIG_INSTALL_DIR`` to be a variable local
   to Ceres. (Chris Sweeney)
#. Fix typo in the comments in ``Jet.h``. (Julius Ziegler)
#. Add the ASL at ETH Zurich, Theia & OpenPTrack to the list of users.
#. Fixed a typo in the documentation. (Richard Stebbing)
#. Fixed a boundary handling bug in the BiCubic interpolation
   code. (Bernhard Zeisl)
#. Fixed a ``MSVC`` compilation bug in the cubic interpolation code
   (Johannes Schönberger)
#. Add covariance related files to the Android build.
#. Update Ubuntu 14.04 installation instructions. (Filippo Basso)
#. Improved logging for linear solver failures.
#. Improved crash messages in ``Problem``.
#. Hide Homebrew related variables in CMake GUI.
#. Add SuiteSparse link dependency for
   compressed_col_sparse_matrix_utils_test.
#. Autodetect Homebrew install prefix on OSX.
#. Lint changes from William Rucklidge and Jim Roseborough.
#. Remove ``using namespace std:`` from ``port.h``
#. Add note about glog not currently compiling against gflags 2.1.
#. Add explicit no sparse linear algebra library available option.
#. Improve some wording in the FAQ. (Vasily Vylkov)
#. Delete Incomplete LQ Factorization.
#. Add a pointer to MacPorts. (Markus Moll)


1.10.0
======

New Features
------------
#. Ceres Solver can now be used to solve general unconstrained
   optimization problems. See the documentation for
   ``GradientProblem`` and ``GradientProblemSolver``.
#. ``Eigen`` can now be as a sparse linear algebra backend. This can
   be done by setting
   ``Solver::Options::sparse_linear_algebra_library_type`` to
   ``EIGEN_SPARSE``. Performance should be comparable to
   ``CX_SPARSE``.

   .. NOTE::

      Because ``Eigen`` is a header only library, and some of the code
      related to sparse Cholesky factorization is LGPL, building Ceres
      with support for Eigen's sparse linear algebra is disabled by
      default and should be enabled explicitly.

   .. NOTE::

      For good performance, use Eigen version 3.2.2 or later.

#. Added ``EIGEN_SPARSE_QR`` algorithm for covariance estimation using
   ``Eigen``'s sparse QR factorization. (Michael Vitus)
#. Faster inner iterations when using multiple threads.
#. Faster ``ITERATIVE_SCHUR`` + ``SCHUR_JACOBI`` for small to medium
   sized problems (see documentation for
   ``Solver::Options::use_explicit_schur_complement``).
#. Faster automatic Schur ordering.
#. Reduced memory usage when solving problems with dynamic sparsity.
#. ``CostFunctionToFunctor`` now supports dynamic number of residuals.
#. A complete re-write of the problem preprocessing phase.
#. ``Solver::Summary::FullReport`` now reports the build configuration
   for Ceres.
#. When building on Android, the ``NDK`` version detection logic has
   been improved.
#. The ``CERES_VERSION`` macro has been improved and replaced with the
   ``CERES_VERSION_STRING`` macro.
#. Added ``Solver::Options::IsValid`` which allows users to validate
   their solver configuration before calling ``Solve``.
#. Added ``Problem::GetCostFunctionForResidualBlock`` and
   ``Problem::GetLossFunctionForResidualBlock``.
#. Added Tukey's loss function. (Michael Vitus)
#. Added RotationMatrixToQuaternion
#. Compute & report timing information for line searches.
#. Autodetect gflags namespace.
#. Expanded ``more_garbow_hillstrom.cc``.
#. Added a pointer to Tal Ben-Nun's MSVC wrapper to the docs.
#. Added the ``<2,3,6>`` Schur template specialization. (Alessandro
   Dal Grande)

Backward Incompatible API Changes
---------------------------------
#. ``NumericDiffFunctor`` has been removed. It's API was broken, and
   the implementation was an unnecessary layer of abstraction over
   ``CostFunctionToFunctor``.
#. ``POLAK_RIBIRERE`` conjugate gradients direction type has been
   renamed to ``POLAK_RIBIERE``.
#. ``Solver::Options::solver_log`` has been removed. If needed this
   iteration callback can easily be implemented in user code.
#. The ``SPARSE_CHOLESKY`` algorithm for covariance estimation has
   been removed. It is not rank revealing and numerically poorly
   behaved. Sparse QR factorization is a much better way to do this.
#. The ``SPARSE_QR`` algorithm for covariance estimation has been
   renamed to ``SUITE_SPARSE_QR`` to be consistent with
   ``EIGEN_SPARSE_QR``.
#. ``Solver::Summary::preconditioner_type`` has been replaced with
   ``Solver::Summary::preconditioner_type_given`` and
   ``Solver::Summary::preconditioner_type_used`` to be more consistent
   with how information about the linear solver is communicated.
#. ``CERES_VERSION`` and ``CERES_ABI_VERSION`` macros were not
   terribly useful. They have been replaced with
   ``CERES_VERSION_MAJOR``, ``CERES_VERSION_MINOR`` ,
   ``CERES_VERSION_REVISION`` and ``CERES_VERSION_ABI`` macros. In
   particular the functionality of ``CERES_VERSION`` is provided by
   ``CERES_VERSION_STRING`` macro.

Bug Fixes
---------
#. Do not try the gradient step if TR step line search fails.
#. Fix missing include in libmv_bundle_adjuster on OSX.
#. Conditionally log evaluation failure warnings.
#. Runtime uses four digits after the decimal in Summary:FullReport.
#. Better options checking for TrustRegionMinimizer.
#. Fix RotationMatrixToAngleAxis when the angle of rotation is near
   PI. (Tobias Strauss)
#. Sometimes gradient norm based convergence would miss a step with a
   substantial solution quality improvement. (Rodney Hoskinson)
#. Ignore warnings from within Eigen/SparseQR (3.2.2).
#. Fix empty Cache HELPSTRING parsing error on OS X 10.10 Yosemite.
#. Fix a formatting error TrustRegionMinimizer logging.
#. Add an explicit include for local_parameterization.h (cooordz)
#. Fix a number of typos in the documentation (Martin Baeuml)
#. Made the logging in TrustRegionMinimizer consistent with
   LineSearchMinimizer.
#. Fix some obsolete documentation in CostFunction::Evaluate.
#. Fix CG solver options for ITERATIVE_SCHUR, which did not copy
   min_num_iterations (Johannes Schönberger)
#. Remove obsolete include of numeric_diff_functor.h. (Martin Baeuml)
#. Fix max. linear solver iterations in ConjugateGradientsSolver
   (Johannes Schönberger)
#. Expand check for lack of a sparse linear algebra library. (Michael
   Samples and Domink Reitzle)
#. Fix Eigen Row/ColMajor bug in NumericDiffCostFunction. (Dominik
   Reitzle)
#. Fix crash in Covariance if # threads > 1 requested without OpenMP.
#. Fixed Malformed regex. (Björn Piltz)
#. Fixed MSVC error C2124: divide or mod by zero. (Björn Piltz)
#. Add missing #include of <limits> for loss functions.
#. Make canned loss functions more robust.
#. Fix type of suppressed compiler warning for Eigen 3.2.0.
#. Suppress unused variable warning from Eigen 3.2.0.
#. Add "make install" to the install instructions.
#. Correct formula in documentation of
   Solver::Options::function_tolerance. (Alessandro Gentilini)
#. Add release flags to iOS toolchain.
#. Fix a broken hyperlink in the documentation. (Henrique Mendonca)
#. Add fixes for multiple definitions of ERROR on Windows to docs.
#. Compile miniglog into Ceres if enabled on all platforms.
#. Add two missing files to Android.mk (Greg Coombe)
#. Fix Cmake error when using miniglog. (Greg Coombe)
#. Don't build miniglog unconditionally as a static library (Björn
   Piltz)
#. Added a missing include. (Björn Piltz)
#. Conditionally disable SparseNormalCholesky.
#. Fix a memory leak in program_test.cc.


1.9.0
=====

New Features
------------
#. Bounds constraints: Support for upper and/or lower bounds on
   parameters when using the trust region minimizer.
#. Dynamic Sparsity: Problems in which the sparsity structure of the
   Jacobian changes over the course of the optimization can now be
   solved much more efficiently. (Richard Stebbing)
#. Improved support for Microsoft Visual C++ including the ability to
   build and ship DLLs. (Björn Piltz, Alex Stewart and Sergey
   Sharybin)
#. Support for building on iOS 6.0 or higher (Jack Feng).
#. Autogeneration of config.h that captures all the defines used to
   build and use Ceres Solver.
#. Simpler and more informative solver termination type
   reporting. (See below for more details)
#. New `website <http://www.ceres-solver.org>`_ based entirely on
   Sphinx.
#. ``AutoDiffLocalParameterization`` allows the use of automatic
   differentiation for defining ``LocalParameterization`` objects
   (Alex Stewart)
#. LBFGS is faster due to fewer memory copies.
#. Parameter blocks are not restricted to be less than 32k in size,
   they can be up to 2G in size.
#. Faster ``SPARSE_NORMAL_CHOLESKY`` solver when using ``CX_SPARSE``
   as the sparse linear algebra library.
#. Added ``Problem::IsParameterBlockPresent`` and
   ``Problem::GetParameterization``.
#. Added the (2,4,9) and (2,4,8) template specializations.
#. An example demonstrating the use of
   DynamicAutoDiffCostFunction. (Joydeep Biswas)
#. Homography estimation example from Blender demonstrating the use of
   a custom ``IterationCallback``. (Sergey Sharybin)
#. Support user passing a custom CMAKE_MODULE_PATH (for BLAS /
   LAPACK).

Backward Incompatible API Changes
---------------------------------
#. ``Solver::Options::linear_solver_ordering`` used to be a naked
   pointer that Ceres took ownership of. This is error prone behaviour
   which leads to problems when copying the ``Solver::Options`` struct
   around. This has been replaced with a ``shared_ptr`` to handle
   ownership correctly across copies.

#. The enum used for reporting the termination/convergence status of
   the solver has been renamed from ``SolverTerminationType`` to
   ``TerminationType``.

   The enum values have also changed. ``FUNCTION_TOLERANCE``,
   ``GRADIENT_TOLERANCE`` and ``PARAMETER_TOLERANCE`` have all been
   replaced by ``CONVERGENCE``.

   ``NUMERICAL_FAILURE`` has been replaced by ``FAILURE``.

   ``USER_ABORT`` has been renamed to ``USER_FAILURE``.

   Further ``Solver::Summary::error`` has been renamed to
   ``Solver::Summary::message``. It contains a more detailed
   explanation for why the solver terminated.

#. ``Solver::Options::gradient_tolerance`` used to be a relative
   gradient tolerance. i.e., The solver converged when

   .. math:: \|g(x)\|_\infty < \text{gradient_tolerance} *
      \|g(x_0)\|_\infty

   where :math:`g(x)` is the gradient of the objective function at
   :math:`x` and :math:`x_0` is the parmeter vector at the start of
   the optimization.

   This has changed to an absolute tolerance, i.e. the solver
   converges when

   .. math:: \|g(x)\|_\infty < \text{gradient_tolerance}

#. Ceres cannot be built without the line search minimizer
   anymore. Thus the preprocessor define
   ``CERES_NO_LINE_SEARCH_MINIMIZER`` has been removed.

Bug Fixes
---------
#. Disabled warning C4251. (Björn Piltz)
#. Do not propagate 3d party libs through
   `IMPORTED_LINK_INTERFACE_LIBRARIES_[DEBUG/RELEASE]` mechanism when
   building shared libraries. (Björn Piltz)
#. Fixed errant verbose levels (Björn Piltz)
#. Variety of code cleanups, optimizations and bug fixes to the line
   search minimizer code (Alex Stewart)
#. Fixed ``BlockSparseMatrix::Transpose`` when the matrix has row and
   column blocks. (Richard Bowen)
#. Better error checking when ``Problem::RemoveResidualBlock`` is
   called. (Alex Stewart)
#. Fixed a memory leak in ``SchurComplementSolver``.
#. Added ``epsilon()`` method to ``NumTraits<ceres::Jet<T, N>
   >``. (Filippo Basso)
#. Fixed a bug in `CompressedRowSparseMatrix::AppendRows`` and
   ``DeleteRows``.q
#. Handle empty problems consistently.
#. Restore the state of the ``Problem`` after a call to
   ``Problem::Evaluate``. (Stefan Leutenegger)
#. Better error checking and reporting for linear solvers.
#. Use explicit formula to solve quadratic polynomials instead of the
   eigenvalue solver.
#. Fix constant parameter handling in inner iterations (Mikael
   Persson).
#. SuiteSparse errors do not cause a fatal crash anymore.
#. Fix ``corrector_test.cc``.
#. Relax the requirements on loss function derivatives.
#. Minor bugfix to logging.h (Scott Ettinger)
#. Updated ``gmock`` and ``gtest`` to the latest upstream version.
#. Fix build breakage on old versions of SuiteSparse.
#. Fixed build issues related to Clang / LLVM 3.4 (Johannes
   Schönberger)
#. METIS_FOUND is never set. Changed the commit to fit the setting of
   the other #._FOUND definitions. (Andreas Franek)
#. Variety of bug fixes and cleanups to the ``CMake`` build system
   (Alex Stewart)
#. Removed fictitious shared library target from the NDK build.
#. Solver::Options now uses ``shared_ptr`` to handle ownership of
   ``Solver::Options::linear_solver_ordering`` and
   ``Solver::Options::inner_iteration_ordering``. As a consequence the
   ``NDK`` build now depends on ``libc++`` from the ``LLVM`` project.
#. Variety of lint cleanups (William Rucklidge & Jim Roseborough)
#. Various internal cleanups including dead code removal.


1.8.0
=====

New Features
------------
#. Significant improved ``CMake`` files with better robustness,
   dependency checking and GUI support. (Alex Stewart)
#. Added ``DynamicNumericDiffCostFunction`` for numerically
   differentiated cost functions whose sizing is determined at run
   time.
#. ``NumericDiffCostFunction`` now supports a dynamic number of
   residuals just like ``AutoDiffCostFunction``.
#. ``Problem`` exposes more of its structure in its API.
#. Faster automatic differentiation (Tim Langlois)
#. Added the commonly occurring ``2_d_d`` template specialization for
   the Schur Eliminator.
#. Faster ``ITERATIVE_SCHUR`` solver using template specializations.
#. Faster ``SCHUR_JACOBI`` preconditioner construction.
#. Faster ``AngleAxisRotatePoint``.
#. Faster Jacobian evaluation when a loss function is used.
#. Added support for multiple clustering algorithms in visibility
   based preconditioning, including a new fast single linkage
   clustering algorithm.

Bug Fixes
---------
#. Fix ordering of ParseCommandLineFlags() & InitGoogleTest() for
   Windows. (Alex Stewart)
#. Remove DCHECK_GE checks from fixed_array.h.
#. Fix build on MSVC 2013 (Petter Strandmark)
#. Fixed ``AngleAxisToRotationMatrix`` near zero.
#. Move ``CERES_HASH_NAMESPACE`` macros to ``collections_port.h``.
#. Fix handling of unordered_map/unordered_set on OSX 10.9.0.
#. Explicitly link to libm for ``curve_fitting_c.c``. (Alex Stewart)
#. Minor type conversion fix to autodiff.h
#. Remove RuntimeNumericDiffCostFunction.
#. Fix operator= ambiguity on some versions of Clang. (Alex Stewart)
#. Various Lint cleanups (William Rucklidge & Jim Roseborough)
#. Modified installation folders for Windows. (Pablo Speciale)
#. Added librt to link libraries for SuiteSparse_config on
   Linux. (Alex Stewart)
#. Check for presence of return-type-c-linkage option with
   Clang. (Alex Stewart)
#. Fix Problem::RemoveParameterBlock after calling solve. (Simon
   Lynen)
#. Fix a free/delete bug in covariance_impl.cc
#. Fix two build errors. (Dustin Lang)
#. Add RequireInitialization = 1 to NumTraits::Jet.
#. Update gmock/gtest to 1.7.0
#. Added IterationSummary::gradient_norm.
#. Reduced verbosity of the inner iteration minimizer.
#. Fixed a bug in TrustRegionMinimizer. (Michael Vitus)
#. Removed android/build_android.sh.


1.7.0
=====

Backward Incompatible API Changes
---------------------------------

#. ``Solver::Options::sparse_linear_algebra_library`` has been renamed
   to ``Solver::Options::sparse_linear_algebra_library_type``.

New Features
------------
#. Sparse and dense covariance estimation.
#. A new Wolfe line search. (Alex Stewart)
#. ``BFGS`` line search direction. (Alex Stewart)
#. C API
#. Speeded up the use of loss functions > 17x.
#. Faster ``DENSE_QR``, ``DENSE_NORMAL_CHOLESKY`` and ``DENSE_SCHUR``
   solvers.
#. Support for multiple dense linear algebra backends. In particular
   optimized ``BLAS`` and ``LAPACK`` implementations (e.g., Intel MKL,
   ACML, OpenBLAS etc) can now be used to do the dense linear algebra
   for ``DENSE_QR``, ``DENSE_NORMAL_CHOLESKY`` and ``DENSE_SCHUR``
#. Use of Inner iterations can now be adaptively stopped. Iteration
   and runtime statistics for inner iterations are not reported in
   ``Solver::Summary`` and ``Solver::Summary::FullReport``.
#. Improved inner iteration step acceptance criterion.
#. Add BlockRandomAccessCRSMatrix.
#. Speeded up automatic differentiation by 7\%.
#. Bundle adjustment example from libmv/Blender (Sergey Sharybin)
#. Shared library building is now controlled by CMake, rather than a
   custom solution. Previously, Ceres had a custom option, but this is
   now deprecated in favor of CMake's built in support for switching
   between static and shared. Turn on BUILD_SHARED_LIBS to get shared
   Ceres libraries.
#. No more dependence on Protocol Buffers.
#. Incomplete LQ factorization.
#. Ability to write trust region problems to disk.
#. Add sinh, cosh, tanh and tan functions to automatic differentiation
   (Johannes Schönberger)
#. Simplifications to the cmake build file.
#. ``miniglog`` can now be used as a replacement for ``google-glog``
   on non Android platforms. (This is NOT recommended).

Bug Fixes
---------
#. Fix ``ITERATIVE_SCHUR`` solver to work correctly when the schur
   complement is of size zero. (Soohyun Bae)
#. Fix the ``spec`` file for generating ``RPM`` packages (Brian Pitts
   and Taylor Braun-Jones).
#. Fix how ceres calls CAMD (Manas Jagadev)
#. Fix breakage on old versions of SuiteSparse. (Fisher Yu)
#. Fix warning C4373 in Visual Studio (Petter Strandmark)
#. Fix compilation error caused by missing suitesparse headers and
   reorganize them to be more robust. (Sergey Sharybin)
#. Check GCC Version before adding -fast compiler option on
   OSX. (Steven Lovegrove)
#. Add documentation for minimizer progress output.
#. Lint and other cleanups (William Rucklidge and James Roseborough)
#. Collections port fix for MSC 2008 (Sergey Sharybin)
#. Various corrections and cleanups in the documentation.
#. Change the path where CeresConfig.cmake is installed (Pablo
   Speciale)
#. Minor errors in documentation (Pablo Speciale)
#. Updated depend.cmake to follow CMake IF convention. (Joydeep
   Biswas)
#. Stabilize the schur ordering algorithm.
#. Update license header in split.h.
#. Enabling -O4 (link-time optimization) only if compiler/linker
   support it. (Alex Stewart)
#. Consistent glog path across files.
#. ceres-solver.spec: Use cleaner, more conventional Release string
   (Taylor Braun-Jones)
#. Fix compile bug on RHEL6 due to missing header (Taylor Braun-Jones)
#. CMake file is less verbose.
#. Use the latest upstream version of google-test and gmock.
#. Rationalize some of the variable names in ``Solver::Options``.
#. Improve Summary::FullReport when line search is used.
#. Expose line search parameters in ``Solver::Options``.
#. Fix update of L-BFGS history buffers after they become full. (Alex
   Stewart)
#. Fix configuration error on systems without SuiteSparse installed
   (Sergey Sharybin)
#. Enforce the read call returns correct value in
   ``curve_fitting_c.c`` (Arnaud Gelas)
#. Fix DynamicAutoDiffCostFunction (Richard Stebbing)
#. Fix Problem::RemoveParameterBlock documentation (Johannes
   Schönberger)
#. Fix a logging bug in parameter_block.h
#. Refactor the preconditioner class structure.
#. Fix an uninitialized variable warning when building with ``GCC``.
#. Fix a reallocation bug in
   ``CreateJacobianBlockSparsityTranspose``. (Yuliy Schwartzburg)
#. Add a define for O_BINARY.
#. Fix miniglog-based Android NDK build; now works with NDK r9. (Scott
   Ettinger)


1.6.0
=====

New Features
------------
#. Major Performance improvements.

   a. Schur type solvers (``SPARSE_SCHUR``, ``DENSE_SCHUR``,
      ``ITERATIVE_SCHUR``) are significantly faster due to custom BLAS
      routines and fewer heap allocations.

   b. ``SPARSE_SCHUR`` when used with ``CX_SPARSE`` now uses a block
      AMD for much improved factorization performance.

   c. The jacobian matrix is pre-ordered so that
      ``SPARSE_NORMAL_CHOLESKY`` and ``SPARSE_SCHUR`` do not have to
      make copies inside ``CHOLMOD``.

   d. Faster autodiff by replacing division by multplication by inverse.

   e. When compiled without threads, the schur eliminator does not pay
      the penalty for locking and unlocking mutexes.

#. Users can now use ``linear_solver_ordering`` to affect the
   fill-reducing ordering used by ``SUITE_SPARSE`` for
   ``SPARSE_NORMAL_CHOLESKY``.
#. ``Problem`` can now report the set of parameter blocks it knows about.
#. ``TrustRegionMinimizer`` uses the evaluator to compute the gradient
   instead of a matrix vector multiply.
#. On ``Mac OS``, whole program optimization is enabled.
#. Users can now use automatic differentiation to define new
   ``LocalParameterization`` objects. (Sergey Sharybin)
#. Enable larger tuple sizes for Visual Studio 2012. (Petter Strandmark)


Bug Fixes
---------

#. Update the documentation for ``CostFunction``.
#. Fixed a typo in the documentation. (Pablo Speciale)
#. Fix a typo in suitesparse.cc.
#. Bugfix in ``NumericDiffCostFunction``. (Nicolas Brodu)
#. Death to BlockSparseMatrixBase.
#. Change Minimizer::Options::min_trust_region_radius to double.
#. Update to compile with stricter gcc checks. (Joydeep Biswas)
#. Do not modify cached CMAKE_CXX_FLAGS_RELEASE. (Sergey Sharybin)
#. ``<iterator>`` needed for back_insert_iterator. (Petter Strandmark)
#. Lint cleanup. (William Rucklidge)
#. Documentation corrections. (Pablo Speciale)


1.5.0
=====

Backward Incompatible API Changes
---------------------------------
#. Added ``Problem::Evaluate``. Now you can evaluate a problem or any
   part of it without calling the solver.

   In light of this the following settings have been deprecated and
   removed from the API.

   - ``Solver::Options::return_initial_residuals``
   - ``Solver::Options::return_initial_gradient``
   - ``Solver::Options::return_initial_jacobian``
   - ``Solver::Options::return_final_residuals``
   - ``Solver::Options::return_final_gradient``
   - ``Solver::Options::return_final_jacobian``

   Instead we recommend using something like this.

   .. code-block:: c++

     Problem problem;
     // Build problem

     vector<double> initial_residuals;
     problem.Evaluate(Problem::EvaluateOptions(),
                      NULL, /* No cost */
                      &initial_residuals,
                      NULL, /* No gradient */
                      NULL  /* No jacobian */);

     Solver::Options options;
     Solver::Summary summary;
     Solver::Solve(options, &problem, &summary);

     vector<double> final_residuals;
     problem.Evaluate(Problem::EvaluateOptions(),
                      NULL, /* No cost */
                      &final_residuals,
                      NULL, /* No gradient */
                      NULL  /* No jacobian */);


New Features
------------
#. Problem now supports removal of ParameterBlocks and
   ResidualBlocks. There is a space/time tradeoff in doing this which
   is controlled by
   ``Problem::Options::enable_fast_parameter_block_removal``.

#. Ceres now supports Line search based optimization algorithms in
   addition to trust region algorithms. Currently there is support for
   gradient descent, non-linear conjugate gradient and LBFGS search
   directions.
#. Added ``Problem::Evaluate``. Now you can evaluate a problem or any
   part of it without calling the solver. In light of this the
   following settings have been deprecated and removed from the API.

   - ``Solver::Options::return_initial_residuals``
   - ``Solver::Options::return_initial_gradient``
   - ``Solver::Options::return_initial_jacobian``
   - ``Solver::Options::return_final_residuals``
   - ``Solver::Options::return_final_gradient``
   - ``Solver::Options::return_final_jacobian``

#. New, much improved HTML documentation using Sphinx.
#. Changed ``NumericDiffCostFunction`` to take functors like
   ``AutoDiffCostFunction``.
#. Added support for mixing automatic, analytic and numeric
   differentiation. This is done by adding ``CostFunctionToFunctor``
   and ``NumericDiffFunctor`` objects to the API.
#. Sped up the robust loss function correction logic when residual is
   one dimensional.
#. Sped up ``DenseQRSolver`` by changing the way dense jacobians are
   stored. This is a 200-500% improvement in linear solver performance
   depending on the size of the problem.
#. ``DENSE_SCHUR`` now supports multi-threading.
#. Greatly expanded ``Summary::FullReport``:

   - Report the ordering used by the ``LinearSolver``.
   - Report the ordering used by the inner iterations.
   - Execution timing breakdown into evaluations and linear solves.
   - Effective size of the problem solved by the solver, which now
     accounts for the size of the tangent space when using a
     ``LocalParameterization``.
#. Ceres when run at the ``VLOG`` level 3 or higher will report
   detailed timing information about its internals.
#. Remove extraneous initial and final residual evaluations. This
   speeds up the solver a bit.
#. Automatic differenatiation with a dynamic number of parameter
   blocks. (Based on an idea by Thad Hughes).
#. Sped up problem construction and destruction.
#. Added matrix adapters to ``rotation.h`` so that the rotation matrix
   routines can work with row and column major matrices. (Markus Moll)
#. ``SCHUR_JACOBI`` can now be used without ``SuiteSparse``.
#. A ``.spec`` file for producing RPMs. (Taylor Braun-Jones)
#. ``CMake`` can now build the sphinx documentation (Pablo Speciale)
#. Add support for creating a CMake config file during build to make
   embedding Ceres in other CMake-using projects easier. (Pablo
   Speciale).
#. Better error reporting in ``Problem`` for missing parameter blocks.
#. A more flexible ``Android.mk`` and a more modular build. If binary
   size and/or compile time is a concern, larger parts of the solver
   can be disabled at compile time.

Bug Fixes
---------
#. Compilation fixes for MSVC2010 (Sergey Sharybin)
#. Fixed "deprecated conversion from string constant to char*"
   warnings. (Pablo Speciale)
#. Correctly propagate ifdefs when building without Schur eliminator
   template specializations.
#. Correct handling of ``LIB_SUFFIX`` on Linux. (Yuliy Schwartzburg).
#. Code and signature cleanup in ``rotation.h``.
#. Make examples independent of internal code.
#. Disable unused member in ``gtest`` which results in build error on
   OS X with latest Xcode. (Taylor Braun-Jones)
#. Pass the correct flags to the linker when using
   ``pthreads``. (Taylor Braun-Jones)
#. Only use ``cmake28`` macro when building on RHEL6. (Taylor
   Braun-Jones)
#. Remove ``-Wno-return-type-c-linkage`` when compiling with
   GCC. (Taylor Braun-Jones)
#. Fix ``No previous prototype`` warnings. (Sergey Sharybin)
#. MinGW build fixes. (Sergey Sharybin)
#. Lots of minor code and lint fixes. (William Rucklidge)
#. Fixed a bug in ``solver_impl.cc`` residual evaluation. (Markus
   Moll)
#. Fixed variadic evaluation bug in ``AutoDiff``.
#. Fixed ``SolverImpl`` tests.
#. Fixed a bug in ``DenseSparseMatrix::ToDenseMatrix()``.
#. Fixed an initialization bug in ``ProgramEvaluator``.
#. Fixes to Android.mk paths (Carlos Hernandez)
#. Modify ``nist.cc`` to compute accuracy based on ground truth
   solution rather than the ground truth function value.
#. Fixed a memory leak in ``cxsparse.cc``. (Alexander Mordvintsev).
#. Fixed the install directory for libraries by correctly handling
   ``LIB_SUFFIX``. (Taylor Braun-Jones)

1.4.0
=====

Backward Incompatible API Changes
---------------------------------
The new ordering API breaks existing code. Here the common case fixes.

**Before**

.. code-block:: c++

 options.linear_solver_type = ceres::DENSE_SCHUR
 options.ordering_type = ceres::SCHUR

**After**


.. code-block:: c++

  options.linear_solver_type = ceres::DENSE_SCHUR


**Before**

.. code-block:: c++

 options.linear_solver_type = ceres::DENSE_SCHUR;
 options.ordering_type = ceres::USER;
 for (int i = 0; i < num_points; ++i) {
   options.ordering.push_back(my_points[i])
 }
 for (int i = 0; i < num_cameras; ++i) {
   options.ordering.push_back(my_cameras[i])
 }
 options.num_eliminate_blocks = num_points;


**After**

.. code-block:: c++

 options.linear_solver_type = ceres::DENSE_SCHUR;
 options.ordering = new ceres::ParameterBlockOrdering;
 for (int i = 0; i < num_points; ++i) {
   options.linear_solver_ordering->AddElementToGroup(my_points[i], 0);
 }
 for (int i = 0; i < num_cameras; ++i) {
   options.linear_solver_ordering->AddElementToGroup(my_cameras[i], 1);
 }


New Features
------------
#. A new richer, more expressive and consistent API for ordering
   parameter blocks.
#. A non-linear generalization of Ruhe & Wedin's Algorithm II. This
   allows the user to use variable projection on separable and
   non-separable non-linear least squares problems. With
   multithreading, this results in significant improvements to the
   convergence behavior of the solver at a small increase in run time.
#. An image denoising example using fields of experts. (Petter
   Strandmark)
#. Defines for Ceres version and ABI version.
#. Higher precision timer code where available. (Petter Strandmark)
#. Example Makefile for users of Ceres.
#. IterationSummary now informs the user when the step is a
   non-monotonic step.
#. Fewer memory allocations when using ``DenseQRSolver``.
#. GradientChecker for testing CostFunctions (William Rucklidge)
#. Add support for cost functions with 10 parameter blocks in
   ``Problem``. (Fisher)
#. Add support for 10 parameter blocks in ``AutoDiffCostFunction``.


Bug Fixes
---------

#. static cast to force Eigen::Index to long conversion
#. Change LOG(ERROR) to LOG(WARNING) in ``schur_complement_solver.cc``.
#. Remove verbose logging from ``DenseQRSolve``.
#. Fix the Android NDK build.
#. Better handling of empty and constant Problems.
#. Remove an internal header that was leaking into the public API.
#. Memory leak in ``trust_region_minimizer.cc``
#. Schur ordering was operating on the wrong object (Ricardo Martin)
#. MSVC fixes (Petter Strandmark)
#. Various fixes to ``nist.cc`` (Markus Moll)
#. Fixed a jacobian scaling bug.
#. Numerically robust computation of ``model_cost_change``.
#. Signed comparison compiler warning fixes (Ricardo Martin)
#. Various compiler warning fixes all over.
#. Inclusion guard fixes (Petter Strandmark)
#. Segfault in test code (Sergey Popov)
#. Replaced ``EXPECT/ASSERT_DEATH`` with the more portable
   ``EXPECT_DEATH_IF_SUPPORTED`` macros.
#. Fixed the camera projection model in Ceres' implementation of
   Snavely's camera model. (Ricardo Martin)


1.3.0
=====

New Features
------------
#. Android Port (Scott Ettinger also contributed to the port)
#. Windows port. (Changchang Wu and Pierre Moulon also contributed to the port)
#. New subspace Dogleg Solver. (Markus Moll)
#. Trust region algorithm now supports the option of non-monotonic steps.
#. New loss functions ``ArcTanLossFunction``, ``TolerantLossFunction``
   and ``ComposedLossFunction``. (James Roseborough).
#. New ``DENSE_NORMAL_CHOLESKY`` linear solver, which uses Eigen's
   LDLT factorization on the normal equations.
#. Cached symbolic factorization when using ``CXSparse``.
   (Petter Strandark)
#. New example ``nist.cc`` and data from the NIST non-linear
   regression test suite. (Thanks to Douglas Bates for suggesting this.)
#. The traditional Dogleg solver now uses an elliptical trust
   region (Markus Moll)
#. Support for returning initial and final gradients & Jacobians.
#. Gradient computation support in the evaluators, with an eye
   towards developing first order/gradient based solvers.
#. A better way to compute ``Solver::Summary::fixed_cost``. (Markus Moll)
#. ``CMake`` support for building documentation, separate examples,
   installing and uninstalling the library and Gerrit hooks (Arnaud
   Gelas)
#. ``SuiteSparse4`` support (Markus Moll)
#. Support for building Ceres without ``TR1`` (This leads to
   slightly slower ``DENSE_SCHUR`` and ``SPARSE_SCHUR`` solvers).
#. ``BALProblem`` can now write a problem back to disk.
#. ``bundle_adjuster`` now allows the user to normalize and perturb the
   problem before solving.
#. Solver progress logging to file.
#. Added ``Program::ToString`` and ``ParameterBlock::ToString`` to
   help with debugging.
#. Ability to build Ceres as a shared library (MacOS and Linux only),
   associated versioning and build release script changes.
#. Portable floating point classification API.


Bug Fixes
---------
#. Fix how invalid step evaluations are handled.
#. Change the slop handling around zero for model cost changes to use
   relative tolerances rather than absolute tolerances.
#. Fix an inadvertant integer to bool conversion. (Petter Strandmark)
#. Do not link to ``libgomp`` when building on
   windows. (Petter Strandmark)
#. Include ``gflags.h`` in ``test_utils.cc``. (Petter
   Strandmark)
#. Use standard random number generation routines. (Petter Strandmark)
#. ``TrustRegionMinimizer`` does not implicitly negate the
   steps that it takes. (Markus Moll)
#. Diagonal scaling allows for equal upper and lower bounds. (Markus Moll)
#. TrustRegionStrategy does not misuse LinearSolver:Summary anymore.
#. Fix Eigen3 Row/Column Major storage issue. (Lena Gieseke)
#. QuaternionToAngleAxis now guarantees an angle in $[-\pi, \pi]$. (Guoxuan Zhang)
#. Added a workaround for a compiler bug in the Android NDK to the
   Schur eliminator.
#. The sparse linear algebra library is only logged in
   Summary::FullReport if it is used.
#. Rename the macro ``CERES_DONT_HAVE_PROTOCOL_BUFFERS``
   to ``CERES_NO_PROTOCOL_BUFFERS`` for consistency.
#. Fix how static structure detection for the Schur eliminator logs
   its results.
#. Correct example code in the documentation. (Petter Strandmark)
#. Fix ``fpclassify.h`` to work with the Android NDK and STLport.
#. Fix a memory leak in the ``levenber_marquardt_strategy_test.cc``
#. Fix an early return bug in the Dogleg solver. (Markus Moll)
#. Zero initialize Jets.
#. Moved ``internal/ceres/mock_log.h`` to ``internal/ceres/gmock/mock-log.h``
#. Unified file path handling in tests.
#. ``data_fitting.cc`` includes ``gflags``
#. Renamed Ceres' Mutex class and associated macros to avoid
   namespace conflicts.
#. Close the BAL problem file after reading it (Markus Moll)
#. Fix IsInfinite on Jets.
#. Drop alignment requirements for Jets.
#. Fixed Jet to integer comparison. (Keith Leung)
#. Fix use of uninitialized arrays. (Sebastian Koch & Markus Moll)
#. Conditionally compile gflag dependencies.(Casey Goodlett)
#. Add ``data_fitting.cc`` to the examples ``CMake`` file.


1.2.3
=====

Bug Fixes
---------
#. ``suitesparse_test`` is enabled even when ``-DSUITESPARSE=OFF``.
#. ``FixedArray`` internal struct did not respect ``Eigen``
   alignment requirements (Koichi Akabe & Stephan Kassemeyer).
#. Fixed ``quadratic.cc`` documentation and code mismatch
   (Nick Lewycky).

1.2.2
=====

Bug Fixes
---------
#. Fix constant parameter blocks, and other minor fixes (Markus Moll)
#. Fix alignment issues when combining ``Jet`` and
   ``FixedArray`` in automatic differeniation.
#. Remove obsolete ``build_defs`` file.

1.2.1
=====

New Features
------------
#. Powell's Dogleg solver
#. Documentation now has a brief overview of Trust Region methods and
   how the Levenberg-Marquardt and Dogleg methods work.

Bug Fixes
---------
#. Destructor for ``TrustRegionStrategy`` was not virtual (Markus
   Moll)
#. Invalid ``DCHECK`` in ``suitesparse.cc`` (Markus Moll)
#. Iteration callbacks were not properly invoked (Luis Alberto
   Zarrabeiti)
#. Logging level changes in ConjugateGradientsSolver
#. VisibilityBasedPreconditioner setup does not account for skipped
   camera pairs. This was debugging code.
#. Enable SSE support on MacOS
#. ``system_test`` was taking too long and too much memory (Koichi
   Akabe)

1.2.0
=====

New Features
------------

#. ``CXSparse`` support.
#. Block oriented fill reducing orderings. This reduces the
   factorization time for sparse ``CHOLMOD`` significantly.
#. New Trust region loop with support for multiple trust region step
   strategies. Currently only Levenberg-Marquardt is supported, but
   this refactoring opens the door for Dog-leg, Stiehaug and others.
#. ``CMake`` file restructuring.  Builds in ``Release`` mode by   default, and now has platform specific tuning flags.
#. Re-organized documentation. No new content, but better
   organization.


Bug Fixes
---------

#. Fixed integer overflow bug in ``block_random_access_sparse_matrix.cc``.
#. Renamed some macros to prevent name conflicts.
#. Fixed incorrect input to ``StateUpdatingCallback``.
#. Fixes to AutoDiff tests.
#. Various internal cleanups.


1.1.1
=====

Bug Fixes
---------
#. Fix a bug in the handling of constant blocks. (Louis Simard)
#. Add an optional lower bound to the Levenberg-Marquardt regularizer
   to prevent oscillating between well and ill posed linear problems.
#. Some internal refactoring and test fixes.

1.1.0
=====

New Features
------------
#. New iterative linear solver for general sparse problems - ``CGNR``
   and a block Jacobi preconditioner for it.
#. Changed the semantics of how ``SuiteSparse`` dependencies are
   checked and used. Now ``SuiteSparse`` is built by default, only if
   all of its dependencies are present.
#. Automatic differentiation now supports dynamic number of residuals.
#. Support for writing the linear least squares problems to disk in
   text format so that they can loaded into ``MATLAB``.
#. Linear solver results are now checked for nan and infinities.
#. Added ``.gitignore`` file.
#. A better more robust build system.


Bug Fixes
---------
#. Fixed a strict weak ordering bug in the schur ordering.
#. Grammar and typos in the documents and code comments.
#. Fixed tests which depended on exact equality between floating point
   values.

1.0.0
=====
Initial open source release. Nathan Wiegand contributed to the Mac OSX
port.


Origins
=======

Ceres Solver grew out of the need for general least squares solving at
Google. In early 2010, Sameer Agarwal and Frederik Schaffalitzky
started the development of Ceres Solver. Frederik left Google shortly
thereafter and Keir Mierle stepped in to take his place. After two
years of on-and-off development, Ceres Solver was released as open
source in May of 2012.
