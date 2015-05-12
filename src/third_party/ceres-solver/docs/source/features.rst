========
Features
========
.. _chapter-features:

* **Code Quality** - Ceres Solver has been used in production at
  Google for more than four years now. It is clean, extensively tested
  and well documented code that is actively developed and supported.

* **Modeling API** - It is rarely the case that one starts with the
  exact and complete formulation of the problem that one is trying to
  solve. Ceres's modeling API has been designed so that the user can
  easily build and modify the objective function, one term at a
  time. And to do so without worrying about how the solver is going to
  deal with the resulting changes in the sparsity/structure of the
  underlying problem.

  - **Derivatives** Supplying derivatives is perhaps the most tedious
    and error prone part of using an optimization library.  Ceres
    ships with `automatic`_ and `numeric`_ differentiation. So you
    never have to compute derivatives by hand (unless you really want
    to). Not only this, Ceres allows you to mix automatic, numeric and
    analytical derivatives in any combination that you want.

  - **Robust Loss Functions** Most non-linear least squares problems
    involve data. If there is data, there will be outliers. Ceres
    allows the user to *shape* their residuals using a
    :class:`LossFunction` to reduce the influence of outliers.

  - **Local Parameterization** In many cases, some parameters lie on a
    manifold other than Euclidean space, e.g., rotation matrices. In
    such cases, the user can specify the geometry of the local tangent
    space by specifying a :class:`LocalParameterization` object.

* **Solver Choice** Depending on the size, sparsity structure, time &
  memory budgets, and solution quality requiremnts, different
  optimization algorithms will suit different needs. To this end,
  Ceres Solver comes with a variety of optimization algorithms:

  - **Trust Region Solvers** - Ceres supports Levenberg-Marquardt,
    Powell's Dogleg, and Subspace dogleg methods. The key
    computational cost in all of these methods is the solution of a
    linear system. To this end Ceres ships with a variety of linear
    solvers - dense QR and dense Cholesky factorization (using
    `Eigen`_ or `LAPACK`_) for dense problems, sparse Cholesky
    factorization (`SuiteSparse`_, `CXSparse`_ or `Eigen`_) for large
    sparse problems custom Schur complement based dense, sparse, and
    iterative linear solvers for `bundle adjustment`_ problems.

  - **Line Search Solvers** - When the problem size is so large that
    storing and factoring the Jacobian is not feasible or a low
    accuracy solution is required cheaply, Ceres offers a number of
    line search based algorithms. This includes a number of variants
    of Non-linear Conjugate Gradients, BFGS and LBFGS.

* **Speed** - Ceres Solver has been extensively optimized, with C++
  templating, hand written linear algebra routines and OpenMP based
  multithreading of the Jacobian evaluation and the linear solvers.

* **Solution Quality** Ceres is the `best performing`_ solver on the NIST
  problem set used by Mondragon and Borchers for benchmarking
  non-linear least squares solvers.

* **Covariance estimation** - Evaluate the sensitivity/uncertainty of
  the solution by evaluating all or part of the covariance
  matrix. Ceres is one of the few solvers that allows you to to do
  this analysis at scale.

* **Community** Since its release as an open source software, Ceres
  has developed an active developer community that contributes new
  features, bug fixes and support.

* **Portability** - Runs on *Linux*, *Windows*, *Mac OS X*, *Android*
  *and iOS*.

* **BSD Licensed** The BSD license offers the flexibility to ship your
  application

.. _best performing: https://groups.google.com/forum/#!topic/ceres-solver/UcicgMPgbXw
.. _bundle adjustment: http://en.wikipedia.org/wiki/Bundle_adjustment
.. _SuiteSparse: http://www.cise.ufl.edu/research/sparse/SuiteSparse/
.. _Eigen: http://eigen.tuxfamily.org/
.. _LAPACK: http://www.netlib.org/lapack/
.. _CXSparse: https://www.cise.ufl.edu/research/sparse/CXSparse/
.. _automatic: http://en.wikipedia.org/wiki/Automatic_differentiation
.. _numeric: http://en.wikipedia.org/wiki/Numerical_differentiation
