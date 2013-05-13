.. _chapter-version-history:

===============
Version History
===============

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
                      NULL  /* No jacobian */ );

     Solver::Options options;
     Solver::Summary summary;
     Solver::Solve(options, &problem, &summary);

     vector<double> final_residuals;
     problem.Evaluate(Problem::EvaluateOptions(),
                      NULL, /* No cost */
                      &final_residuals,
                      NULL, /* No gradient */
                      NULL  /* No jacobian */ );


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

#. Fixed varidic evaluation bug in ``AutoDiff``.

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

#. Destructor for ``TrustRegionStrategy`` was not virtual (Markus Moll)

#. Invalid ``DCHECK`` in ``suitesparse.cc`` (Markus Moll)

#. Iteration callbacks were not properly invoked (Luis Alberto Zarrabeiti)

#. Logging level changes in ConjugateGradientsSolver

#. VisibilityBasedPreconditioner setup does not account for skipped camera pairs. This was debugging code.

#. Enable SSE support on MacOS

#. ``system_test`` was taking too long and too much memory (Koichi Akabe)

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

#. ``CMake`` file restructuring.  Builds in ``Release`` mode by
   default, and now has platform specific tuning flags.

#. Re-organized documentation. No new content, but better
   organization.


Bug Fixes
---------

#. Fixed integer overflow bug in ``block_random_access_sparse_matrix.cc``.

#. Renamed some macros to prevent name conflicts.

#. Fixed incorrent input to ``StateUpdatingCallback``.

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

#. Fixed tests which depended on exact equality between floating point values.

1.0.0
=====

Initial Release. Nathan Wiegand contributed to the Mac OSX port.
