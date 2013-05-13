
.. default-domain:: cpp

.. cpp:namespace:: ceres

.. _chapter-solving:

==========
Solver API
==========


Introduction
============

Effective use of Ceres requires some familiarity with the basic
components of a nonlinear least squares solver, so before we describe
how to configure the solver, we will begin by taking a brief look at
how some of the core optimization algorithms in Ceres work and the
various linear solvers and preconditioners that power it.


Let :math:`x \in \mathbb{R}^n` be an :math:`n`-dimensional vector of
variables, and
:math:`F(x) = \left[f_1(x), ... ,  f_{m}(x) \right]^{\top}` be a
:math:`m`-dimensional function of :math:`x`.  We are interested in
solving the following optimization problem [#f1]_ .

.. math:: \arg \min_x \frac{1}{2}\|F(x)\|^2\ .
  :label: nonlinsq

Here, the Jacobian :math:`J(x)` of :math:`F(x)` is an :math:`m\times
n` matrix, where :math:`J_{ij}(x) = \partial_j f_i(x)` and the
gradient vector :math:`g(x) = \nabla \frac{1}{2}\|F(x)\|^2 = J(x)^\top
F(x)`. Since the efficient global minimization of :eq:`nonlinsq` for
general :math:`F(x)` is an intractable problem, we will have to settle
for finding a local minimum.

The general strategy when solving non-linear optimization problems is
to solve a sequence of approximations to the original problem
[NocedalWright]_. At each iteration, the approximation is solved to
determine a correction :math:`\Delta x` to the vector :math:`x`. For
non-linear least squares, an approximation can be constructed by using
the linearization :math:`F(x+\Delta x) \approx F(x) + J(x)\Delta x`,
which leads to the following linear least squares problem:

.. math:: \min_{\Delta x} \frac{1}{2}\|J(x)\Delta x + F(x)\|^2
   :label: linearapprox

Unfortunately, naively solving a sequence of these problems and
updating :math:`x \leftarrow x+ \Delta x` leads to an algorithm that
may not converge.  To get a convergent algorithm, we need to control
the size of the step :math:`\Delta x`. Depending on how the size of
the step :math:`\Delta x` is controlled, non-linear optimization
algorithms can be divided into two major categories [NocedalWright]_.

1. **Trust Region** The trust region approach approximates the
   objective function using using a model function (often a quadratic)
   over a subset of the search space known as the trust region. If the
   model function succeeds in minimizing the true objective function
   the trust region is expanded; conversely, otherwise it is
   contracted and the model optimization problem is solved again.

2. **Line Search** The line search approach first finds a descent
   direction along which the objective function will be reduced and
   then computes a step size that decides how far should move along
   that direction. The descent direction can be computed by various
   methods, such as gradient descent, Newton's method and Quasi-Newton
   method. The step size can be determined either exactly or
   inexactly.

Trust region methods are in some sense dual to line search methods:
trust region methods first choose a step size (the size of the trust
region) and then a step direction while line search methods first
choose a step direction and then a step size.

Ceres implements multiple algorithms in both categories.

.. _section-trust-region-methods:

Trust Region Methods
====================

The basic trust region algorithm looks something like this.

   1. Given an initial point :math:`x` and a trust region radius :math:`\mu`.
   2. :math:`\arg \min_{\Delta x} \frac{1}{2}\|J(x)\Delta
      x + F(x)\|^2` s.t. :math:`\|D(x)\Delta x\|^2 \le \mu`
   3. :math:`\rho = \frac{\displaystyle \|F(x + \Delta x)\|^2 -
      \|F(x)\|^2}{\displaystyle \|J(x)\Delta x + F(x)\|^2 -
      \|F(x)\|^2}`
   4. if :math:`\rho > \epsilon` then  :math:`x = x + \Delta x`.
   5. if :math:`\rho > \eta_1` then :math:`\rho = 2  \rho`
   6. else if :math:`\rho < \eta_2` then :math:`\rho = 0.5 * \rho`
   7. Goto 2.

Here, :math:`\mu` is the trust region radius, :math:`D(x)` is some
matrix used to define a metric on the domain of :math:`F(x)` and
:math:`\rho` measures the quality of the step :math:`\Delta x`, i.e.,
how well did the linear model predict the decrease in the value of the
non-linear objective. The idea is to increase or decrease the radius
of the trust region depending on how well the linearization predicts
the behavior of the non-linear objective, which in turn is reflected
in the value of :math:`\rho`.

The key computational step in a trust-region algorithm is the solution
of the constrained optimization problem

.. math:: \arg\min_{\Delta x} \frac{1}{2}\|J(x)\Delta x +  F(x)\|^2\quad \text{such that}\quad  \|D(x)\Delta x\|^2 \le \mu
   :label: trp

There are a number of different ways of solving this problem, each
giving rise to a different concrete trust-region algorithm. Currently
Ceres, implements two trust-region algorithms - Levenberg-Marquardt
and Dogleg. The user can choose between them by setting
:member:`Solver::Options::trust_region_strategy_type`.

.. rubric:: Footnotes

.. [#f1] At the level of the non-linear solver, the block and
         structure is not relevant, therefore our discussion here is
         in terms of an optimization problem defined over a state
         vector of size :math:`n`.


.. _section-levenberg-marquardt:

Levenberg-Marquardt
-------------------

The Levenberg-Marquardt algorithm [Levenberg]_ [Marquardt]_ is the
most popular algorithm for solving non-linear least squares problems.
It was also the first trust region algorithm to be developed
[Levenberg]_ [Marquardt]_. Ceres implements an exact step [Madsen]_
and an inexact step variant of the Levenberg-Marquardt algorithm
[WrightHolt]_ [NashSofer]_.

It can be shown, that the solution to :eq:`trp` can be obtained by
solving an unconstrained optimization of the form

.. math:: \arg\min_{\Delta x}& \frac{1}{2}\|J(x)\Delta x + F(x)\|^2 +\lambda  \|D(x)\Delta x\|^2

Where, :math:`\lambda` is a Lagrange multiplier that is inverse
related to :math:`\mu`. In Ceres, we solve for

.. math:: \arg\min_{\Delta x}& \frac{1}{2}\|J(x)\Delta x + F(x)\|^2 + \frac{1}{\mu} \|D(x)\Delta x\|^2
   :label: lsqr

The matrix :math:`D(x)` is a non-negative diagonal matrix, typically
the square root of the diagonal of the matrix :math:`J(x)^\top J(x)`.

Before going further, let us make some notational simplifications. We
will assume that the matrix :math:`\sqrt{\mu} D` has been concatenated
at the bottom of the matrix :math:`J` and similarly a vector of zeros
has been added to the bottom of the vector :math:`f` and the rest of
our discussion will be in terms of :math:`J` and :math:`f`, i.e, the
linear least squares problem.

.. math:: \min_{\Delta x} \frac{1}{2} \|J(x)\Delta x + f(x)\|^2 .
   :label: simple

For all but the smallest problems the solution of :eq:`simple` in
each iteration of the Levenberg-Marquardt algorithm is the dominant
computational cost in Ceres. Ceres provides a number of different
options for solving :eq:`simple`. There are two major classes of
methods - factorization and iterative.

The factorization methods are based on computing an exact solution of
:eq:`lsqr` using a Cholesky or a QR factorization and lead to an exact
step Levenberg-Marquardt algorithm. But it is not clear if an exact
solution of :eq:`lsqr` is necessary at each step of the LM algorithm
to solve :eq:`nonlinsq`. In fact, we have already seen evidence
that this may not be the case, as :eq:`lsqr` is itself a regularized
version of :eq:`linearapprox`. Indeed, it is possible to
construct non-linear optimization algorithms in which the linearized
problem is solved approximately. These algorithms are known as inexact
Newton or truncated Newton methods [NocedalWright]_.

An inexact Newton method requires two ingredients. First, a cheap
method for approximately solving systems of linear
equations. Typically an iterative linear solver like the Conjugate
Gradients method is used for this
purpose [NocedalWright]_. Second, a termination rule for
the iterative solver. A typical termination rule is of the form

.. math:: \|H(x) \Delta x + g(x)\| \leq \eta_k \|g(x)\|.
   :label: inexact

Here, :math:`k` indicates the Levenberg-Marquardt iteration number and
:math:`0 < \eta_k <1` is known as the forcing sequence.  [WrightHolt]_
prove that a truncated Levenberg-Marquardt algorithm that uses an
inexact Newton step based on :eq:`inexact` converges for any
sequence :math:`\eta_k \leq \eta_0 < 1` and the rate of convergence
depends on the choice of the forcing sequence :math:`\eta_k`.

Ceres supports both exact and inexact step solution strategies. When
the user chooses a factorization based linear solver, the exact step
Levenberg-Marquardt algorithm is used. When the user chooses an
iterative linear solver, the inexact step Levenberg-Marquardt
algorithm is used.

.. _section-dogleg:

Dogleg
------

Another strategy for solving the trust region problem :eq:`trp` was
introduced by M. J. D. Powell. The key idea there is to compute two
vectors

.. math::

        \Delta x^{\text{Gauss-Newton}} &= \arg \min_{\Delta x}\frac{1}{2} \|J(x)\Delta x + f(x)\|^2.\\
        \Delta x^{\text{Cauchy}} &= -\frac{\|g(x)\|^2}{\|J(x)g(x)\|^2}g(x).

Note that the vector :math:`\Delta x^{\text{Gauss-Newton}}` is the
solution to :eq:`linearapprox` and :math:`\Delta
x^{\text{Cauchy}}` is the vector that minimizes the linear
approximation if we restrict ourselves to moving along the direction
of the gradient. Dogleg methods finds a vector :math:`\Delta x`
defined by :math:`\Delta x^{\text{Gauss-Newton}}` and :math:`\Delta
x^{\text{Cauchy}}` that solves the trust region problem. Ceres
supports two variants that can be chose by setting
:member:`Solver::Options::dogleg_type`.

``TRADITIONAL_DOGLEG`` as described by Powell, constructs two line
segments using the Gauss-Newton and Cauchy vectors and finds the point
farthest along this line shaped like a dogleg (hence the name) that is
contained in the trust-region. For more details on the exact reasoning
and computations, please see Madsen et al [Madsen]_.

``SUBSPACE_DOGLEG`` is a more sophisticated method that considers the
entire two dimensional subspace spanned by these two vectors and finds
the point that minimizes the trust region problem in this subspace
[ByrdSchnabel]_.

The key advantage of the Dogleg over Levenberg Marquardt is that if
the step computation for a particular choice of :math:`\mu` does not
result in sufficient decrease in the value of the objective function,
Levenberg-Marquardt solves the linear approximation from scratch with
a smaller value of :math:`\mu`. Dogleg on the other hand, only needs
to compute the interpolation between the Gauss-Newton and the Cauchy
vectors, as neither of them depend on the value of :math:`\mu`.

The Dogleg method can only be used with the exact factorization based
linear solvers.

.. _section-inner-iterations:

Inner Iterations
----------------

Some non-linear least squares problems have additional structure in
the way the parameter blocks interact that it is beneficial to modify
the way the trust region step is computed. e.g., consider the
following regression problem

.. math::   y = a_1 e^{b_1 x} + a_2 e^{b_3 x^2 + c_1}


Given a set of pairs :math:`\{(x_i, y_i)\}`, the user wishes to estimate
:math:`a_1, a_2, b_1, b_2`, and :math:`c_1`.

Notice that the expression on the left is linear in :math:`a_1` and
:math:`a_2`, and given any value for :math:`b_1, b_2` and :math:`c_1`,
it is possible to use linear regression to estimate the optimal values
of :math:`a_1` and :math:`a_2`. It's possible to analytically
eliminate the variables :math:`a_1` and :math:`a_2` from the problem
entirely. Problems like these are known as separable least squares
problem and the most famous algorithm for solving them is the Variable
Projection algorithm invented by Golub & Pereyra [GolubPereyra]_.

Similar structure can be found in the matrix factorization with
missing data problem. There the corresponding algorithm is known as
Wiberg's algorithm [Wiberg]_.

Ruhe & Wedin present an analysis of various algorithms for solving
separable non-linear least squares problems and refer to *Variable
Projection* as Algorithm I in their paper [RuheWedin]_.

Implementing Variable Projection is tedious and expensive. Ruhe &
Wedin present a simpler algorithm with comparable convergence
properties, which they call Algorithm II.  Algorithm II performs an
additional optimization step to estimate :math:`a_1` and :math:`a_2`
exactly after computing a successful Newton step.


This idea can be generalized to cases where the residual is not
linear in :math:`a_1` and :math:`a_2`, i.e.,

.. math:: y = f_1(a_1, e^{b_1 x}) + f_2(a_2, e^{b_3 x^2 + c_1})

In this case, we solve for the trust region step for the full problem,
and then use it as the starting point to further optimize just `a_1`
and `a_2`. For the linear case, this amounts to doing a single linear
least squares solve. For non-linear problems, any method for solving
the `a_1` and `a_2` optimization problems will do. The only constraint
on `a_1` and `a_2` (if they are two different parameter block) is that
they do not co-occur in a residual block.

This idea can be further generalized, by not just optimizing
:math:`(a_1, a_2)`, but decomposing the graph corresponding to the
Hessian matrix's sparsity structure into a collection of
non-overlapping independent sets and optimizing each of them.

Setting :member:`Solver::Options::use_inner_iterations` to ``true``
enables the use of this non-linear generalization of Ruhe & Wedin's
Algorithm II.  This version of Ceres has a higher iteration
complexity, but also displays better convergence behavior per
iteration.

Setting :member:`Solver::Options::num_threads` to the maximum number
possible is highly recommended.

.. _section-non-monotonic-steps:

Non-monotonic Steps
-------------------

Note that the basic trust-region algorithm described in
Algorithm~\ref{alg:trust-region} is a descent algorithm in that they
only accepts a point if it strictly reduces the value of the objective
function.

Relaxing this requirement allows the algorithm to be more efficient in
the long term at the cost of some local increase in the value of the
objective function.

This is because allowing for non-decreasing objective function values
in a princpled manner allows the algorithm to *jump over boulders* as
the method is not restricted to move into narrow valleys while
preserving its convergence properties.

Setting :member:`Solver::Options::use_nonmonotonic_steps` to ``true``
enables the non-monotonic trust region algorithm as described by Conn,
Gould & Toint in [Conn]_.

Even though the value of the objective function may be larger
than the minimum value encountered over the course of the
optimization, the final parameters returned to the user are the
ones corresponding to the minimum cost over all iterations.

The option to take non-monotonic steps is available for all trust
region strategies.


.. _section-line-search-methods:

Line Search Methods
===================

**The implementation of line search algorithms in Ceres Solver is
fairly new and not very well tested, so for now this part of the
solver should be considered beta quality. We welcome reports of your
experiences both good and bad on the mailinglist.**

Line search algorithms

   1. Given an initial point :math:`x`
   2. :math:`\Delta x = -H^{-1}(x) g(x)`
   3. :math:`\arg \min_\mu \frac{1}{2} \| F(x + \mu \Delta x) \|^2`
   4. :math:`x = x + \mu \Delta x`
   5. Goto 2.

Here :math:`H(x)` is some approximation to the Hessian of the
objective function, and :math:`g(x)` is the gradient at
:math:`x`. Depending on the choice of :math:`H(x)` we get a variety of
different search directions -`\Delta x`.

Step 4, which is a one dimensional optimization or `Line Search` along
:math:`\Delta x` is what gives this class of methods its name.

Different line search algorithms differ in their choice of the search
direction :math:`\Delta x` and the method used for one dimensional
optimization along :math:`\Delta x`. The choice of :math:`H(x)` is the
primary source of computational complexity in these
methods. Currently, Ceres Solver supports three choices of search
directions, all aimed at large scale problems.

1. ``STEEPEST_DESCENT`` This corresponds to choosing :math:`H(x)` to
   be the identity matrix. This is not a good search direction for
   anything but the simplest of the problems. It is only included here
   for completeness.

2. ``NONLINEAR_CONJUGATE_GRADIENT`` A generalization of the Conjugate
   Gradient method to non-linear functions. The generalization can be
   performed in a number of different ways, resulting in a variety of
   search directions. Ceres Solver currently supports
   ``FLETCHER_REEVES``, ``POLAK_RIBIRERE`` and ``HESTENES_STIEFEL``
   directions.

3. ``LBFGS`` In this method, a limited memory approximation to the
   inverse Hessian is maintained and used to compute a quasi-Newton
   step [Nocedal]_, [ByrdNocedal]_.

Currently Ceres Solver uses a backtracking and interpolation based
Armijo line search algorithm.

.. _section-linear-solver:

LinearSolver
============

Recall that in both of the trust-region methods described above, the
key computational cost is the solution of a linear least squares
problem of the form

.. math:: \min_{\Delta x} \frac{1}{2} \|J(x)\Delta x + f(x)\|^2 .
   :label: simple2

Let :math:`H(x)= J(x)^\top J(x)` and :math:`g(x) = -J(x)^\top
f(x)`. For notational convenience let us also drop the dependence on
:math:`x`. Then it is easy to see that solving :eq:`simple2` is
equivalent to solving the *normal equations*.

.. math:: H \Delta x = g
   :label: normal

Ceres provides a number of different options for solving :eq:`normal`.

.. _section-qr:

``DENSE_QR``
------------

For small problems (a couple of hundred parameters and a few thousand
residuals) with relatively dense Jacobians, ``DENSE_QR`` is the method
of choice [Bjorck]_. Let :math:`J = QR` be the QR-decomposition of
:math:`J`, where :math:`Q` is an orthonormal matrix and :math:`R` is
an upper triangular matrix [TrefethenBau]_. Then it can be shown that
the solution to :eq:`normal` is given by

.. math:: \Delta x^* = -R^{-1}Q^\top f


Ceres uses ``Eigen`` 's dense QR factorization routines.

.. _section-cholesky:

``DENSE_NORMAL_CHOLESKY`` & ``SPARSE_NORMAL_CHOLESKY``
------------------------------------------------------

Large non-linear least square problems are usually sparse. In such
cases, using a dense QR factorization is inefficient. Let :math:`H =
R^\top R` be the Cholesky factorization of the normal equations, where
:math:`R` is an upper triangular matrix, then the solution to
:eq:`normal` is given by

.. math::

    \Delta x^* = R^{-1} R^{-\top} g.


The observant reader will note that the :math:`R` in the Cholesky
factorization of :math:`H` is the same upper triangular matrix
:math:`R` in the QR factorization of :math:`J`. Since :math:`Q` is an
orthonormal matrix, :math:`J=QR` implies that :math:`J^\top J = R^\top
Q^\top Q R = R^\top R`. There are two variants of Cholesky
factorization -- sparse and dense.

``DENSE_NORMAL_CHOLESKY``  as the name implies performs a dense
Cholesky factorization of the normal equations. Ceres uses
``Eigen`` 's dense LDLT factorization routines.

``SPARSE_NORMAL_CHOLESKY``, as the name implies performs a sparse
Cholesky factorization of the normal equations. This leads to
substantial savings in time and memory for large sparse
problems. Ceres uses the sparse Cholesky factorization routines in
Professor Tim Davis' ``SuiteSparse`` or ``CXSparse`` packages [Chen]_.

.. _section-schur:

``DENSE_SCHUR`` & ``SPARSE_SCHUR``
----------------------------------

While it is possible to use ``SPARSE_NORMAL_CHOLESKY`` to solve bundle
adjustment problems, bundle adjustment problem have a special
structure, and a more efficient scheme for solving :eq:`normal`
can be constructed.

Suppose that the SfM problem consists of :math:`p` cameras and
:math:`q` points and the variable vector :math:`x` has the block
structure :math:`x = [y_{1}, ... ,y_{p},z_{1}, ... ,z_{q}]`. Where,
:math:`y` and :math:`z` correspond to camera and point parameters,
respectively.  Further, let the camera blocks be of size :math:`c` and
the point blocks be of size :math:`s` (for most problems :math:`c` =
:math:`6`--`9` and :math:`s = 3`). Ceres does not impose any constancy
requirement on these block sizes, but choosing them to be constant
simplifies the exposition.

A key characteristic of the bundle adjustment problem is that there is
no term :math:`f_{i}` that includes two or more point blocks.  This in
turn implies that the matrix :math:`H` is of the form

.. math:: H = \left[ \begin{matrix} B & E\\ E^\top & C \end{matrix} \right]\ ,
   :label: hblock

where, :math:`B \in \mathbb{R}^{pc\times pc}` is a block sparse matrix
with :math:`p` blocks of size :math:`c\times c` and :math:`C \in
\mathbb{R}^{qs\times qs}` is a block diagonal matrix with :math:`q` blocks
of size :math:`s\times s`. :math:`E \in \mathbb{R}^{pc\times qs}` is a
general block sparse matrix, with a block of size :math:`c\times s`
for each observation. Let us now block partition :math:`\Delta x =
[\Delta y,\Delta z]` and :math:`g=[v,w]` to restate :eq:`normal`
as the block structured linear system

.. math:: \left[ \begin{matrix} B & E\\ E^\top & C \end{matrix}
                \right]\left[ \begin{matrix} \Delta y \\ \Delta z
            	    \end{matrix} \right] = \left[ \begin{matrix} v\\ w
                    \end{matrix} \right]\ ,
   :label: linear2

and apply Gaussian elimination to it. As we noted above, :math:`C` is
a block diagonal matrix, with small diagonal blocks of size
:math:`s\times s`.  Thus, calculating the inverse of :math:`C` by
inverting each of these blocks is cheap. This allows us to eliminate
:math:`\Delta z` by observing that :math:`\Delta z = C^{-1}(w - E^\top
\Delta y)`, giving us

.. math:: \left[B - EC^{-1}E^\top\right] \Delta y = v - EC^{-1}w\ .
   :label: schur

The matrix

.. math:: S = B - EC^{-1}E^\top

is the Schur complement of :math:`C` in :math:`H`. It is also known as
the *reduced camera matrix*, because the only variables
participating in :eq:`schur` are the ones corresponding to the
cameras. :math:`S \in \mathbb{R}^{pc\times pc}` is a block structured
symmetric positive definite matrix, with blocks of size :math:`c\times
c`. The block :math:`S_{ij}` corresponding to the pair of images
:math:`i` and :math:`j` is non-zero if and only if the two images
observe at least one common point.


Now, eq-linear2 can be solved by first forming :math:`S`, solving for
:math:`\Delta y`, and then back-substituting :math:`\Delta y` to
obtain the value of :math:`\Delta z`.  Thus, the solution of what was
an :math:`n\times n`, :math:`n=pc+qs` linear system is reduced to the
inversion of the block diagonal matrix :math:`C`, a few matrix-matrix
and matrix-vector multiplies, and the solution of block sparse
:math:`pc\times pc` linear system :eq:`schur`.  For almost all
problems, the number of cameras is much smaller than the number of
points, :math:`p \ll q`, thus solving :eq:`schur` is
significantly cheaper than solving :eq:`linear2`. This is the
*Schur complement trick* [Brown]_.

This still leaves open the question of solving :eq:`schur`. The
method of choice for solving symmetric positive definite systems
exactly is via the Cholesky factorization [TrefethenBau]_ and
depending upon the structure of the matrix, there are, in general, two
options. The first is direct factorization, where we store and factor
:math:`S` as a dense matrix [TrefethenBau]_. This method has
:math:`O(p^2)` space complexity and :math:`O(p^3)` time complexity and
is only practical for problems with up to a few hundred cameras. Ceres
implements this strategy as the ``DENSE_SCHUR`` solver.


But, :math:`S` is typically a fairly sparse matrix, as most images
only see a small fraction of the scene. This leads us to the second
option: Sparse Direct Methods. These methods store :math:`S` as a
sparse matrix, use row and column re-ordering algorithms to maximize
the sparsity of the Cholesky decomposition, and focus their compute
effort on the non-zero part of the factorization [Chen]_. Sparse
direct methods, depending on the exact sparsity structure of the Schur
complement, allow bundle adjustment algorithms to significantly scale
up over those based on dense factorization. Ceres implements this
strategy as the ``SPARSE_SCHUR`` solver.

.. _section-cgnr:

``CGNR``
--------

For general sparse problems, if the problem is too large for
``CHOLMOD`` or a sparse linear algebra library is not linked into
Ceres, another option is the ``CGNR`` solver. This solver uses the
Conjugate Gradients solver on the *normal equations*, but without
forming the normal equations explicitly. It exploits the relation

.. math::
    H x = J^\top J x = J^\top(J x)


When the user chooses ``ITERATIVE_SCHUR`` as the linear solver, Ceres
automatically switches from the exact step algorithm to an inexact
step algorithm.

.. _section-iterative_schur:

``ITERATIVE_SCHUR``
-------------------

Another option for bundle adjustment problems is to apply PCG to the
reduced camera matrix :math:`S` instead of :math:`H`. One reason to do
this is that :math:`S` is a much smaller matrix than :math:`H`, but
more importantly, it can be shown that :math:`\kappa(S)\leq
\kappa(H)`.  Cseres implements PCG on :math:`S` as the
``ITERATIVE_SCHUR`` solver. When the user chooses ``ITERATIVE_SCHUR``
as the linear solver, Ceres automatically switches from the exact step
algorithm to an inexact step algorithm.

The cost of forming and storing the Schur complement :math:`S` can be
prohibitive for large problems. Indeed, for an inexact Newton solver
that computes :math:`S` and runs PCG on it, almost all of its time is
spent in constructing :math:`S`; the time spent inside the PCG
algorithm is negligible in comparison. Because PCG only needs access
to :math:`S` via its product with a vector, one way to evaluate
:math:`Sx` is to observe that

.. math::  x_1 &= E^\top x
.. math::  x_2 &= C^{-1} x_1
.. math::  x_3 &= Ex_2\\
.. math::  x_4 &= Bx\\
.. math::   Sx &= x_4 - x_3
   :label: schurtrick1

Thus, we can run PCG on :math:`S` with the same computational effort
per iteration as PCG on :math:`H`, while reaping the benefits of a
more powerful preconditioner. In fact, we do not even need to compute
:math:`H`, :eq:`schurtrick1` can be implemented using just the columns
of :math:`J`.

Equation :eq:`schurtrick1` is closely related to *Domain
Decomposition methods* for solving large linear systems that arise in
structural engineering and partial differential equations. In the
language of Domain Decomposition, each point in a bundle adjustment
problem is a domain, and the cameras form the interface between these
domains. The iterative solution of the Schur complement then falls
within the sub-category of techniques known as Iterative
Sub-structuring [Saad]_ [Mathew]_.

.. _section-preconditioner:

Preconditioner
--------------

The convergence rate of Conjugate Gradients for
solving :eq:`normal` depends on the distribution of eigenvalues
of :math:`H` [Saad]_. A useful upper bound is
:math:`\sqrt{\kappa(H)}`, where, :math:`\kappa(H)` is the condition
number of the matrix :math:`H`. For most bundle adjustment problems,
:math:`\kappa(H)` is high and a direct application of Conjugate
Gradients to :eq:`normal` results in extremely poor performance.

The solution to this problem is to replace :eq:`normal` with a
*preconditioned* system.  Given a linear system, :math:`Ax =b` and a
preconditioner :math:`M` the preconditioned system is given by
:math:`M^{-1}Ax = M^{-1}b`. The resulting algorithm is known as
Preconditioned Conjugate Gradients algorithm (PCG) and its worst case
complexity now depends on the condition number of the *preconditioned*
matrix :math:`\kappa(M^{-1}A)`.

The computational cost of using a preconditioner :math:`M` is the cost
of computing :math:`M` and evaluating the product :math:`M^{-1}y` for
arbitrary vectors :math:`y`. Thus, there are two competing factors to
consider: How much of :math:`H`'s structure is captured by :math:`M`
so that the condition number :math:`\kappa(HM^{-1})` is low, and the
computational cost of constructing and using :math:`M`.  The ideal
preconditioner would be one for which :math:`\kappa(M^{-1}A)
=1`. :math:`M=A` achieves this, but it is not a practical choice, as
applying this preconditioner would require solving a linear system
equivalent to the unpreconditioned problem.  It is usually the case
that the more information :math:`M` has about :math:`H`, the more
expensive it is use. For example, Incomplete Cholesky factorization
based preconditioners have much better convergence behavior than the
Jacobi preconditioner, but are also much more expensive.


The simplest of all preconditioners is the diagonal or Jacobi
preconditioner, i.e., :math:`M=\operatorname{diag}(A)`, which for
block structured matrices like :math:`H` can be generalized to the
block Jacobi preconditioner.

For ``ITERATIVE_SCHUR`` there are two obvious choices for block
diagonal preconditioners for :math:`S`. The block diagonal of the
matrix :math:`B` [Mandel]_ and the block diagonal :math:`S`, i.e, the
block Jacobi preconditioner for :math:`S`. Ceres's implements both of
these preconditioners and refers to them as ``JACOBI`` and
``SCHUR_JACOBI`` respectively.

For bundle adjustment problems arising in reconstruction from
community photo collections, more effective preconditioners can be
constructed by analyzing and exploiting the camera-point visibility
structure of the scene [KushalAgarwal]. Ceres implements the two
visibility based preconditioners described by Kushal & Agarwal as
``CLUSTER_JACOBI`` and ``CLUSTER_TRIDIAGONAL``. These are fairly new
preconditioners and Ceres' implementation of them is in its early
stages and is not as mature as the other preconditioners described
above.

.. _section-ordering:

Ordering
--------

The order in which variables are eliminated in a linear solver can
have a significant of impact on the efficiency and accuracy of the
method. For example when doing sparse Cholesky factorization, there
are matrices for which a good ordering will give a Cholesky factor
with :math:`O(n)` storage, where as a bad ordering will result in an
completely dense factor.

Ceres allows the user to provide varying amounts of hints to the
solver about the variable elimination ordering to use. This can range
from no hints, where the solver is free to decide the best ordering
based on the user's choices like the linear solver being used, to an
exact order in which the variables should be eliminated, and a variety
of possibilities in between.

Instances of the :class:`ParameterBlockOrdering` class are used to
communicate this information to Ceres.

Formally an ordering is an ordered partitioning of the parameter
blocks. Each parameter block belongs to exactly one group, and each
group has a unique integer associated with it, that determines its
order in the set of groups. We call these groups *Elimination Groups*

Given such an ordering, Ceres ensures that the parameter blocks in the
lowest numbered elimination group are eliminated first, and then the
parameter blocks in the next lowest numbered elimination group and so
on. Within each elimination group, Ceres is free to order the
parameter blocks as it chooses. e.g. Consider the linear system

.. math::
  x + y &= 3\\
  2x + 3y &= 7

There are two ways in which it can be solved. First eliminating
:math:`x` from the two equations, solving for y and then back
substituting for :math:`x`, or first eliminating :math:`y`, solving
for :math:`x` and back substituting for :math:`y`. The user can
construct three orderings here.

1. :math:`\{0: x\}, \{1: y\}` : Eliminate :math:`x` first.
2. :math:`\{0: y\}, \{1: x\}` : Eliminate :math:`y` first.
3. :math:`\{0: x, y\}`        : Solver gets to decide the elimination order.

Thus, to have Ceres determine the ordering automatically using
heuristics, put all the variables in the same elimination group. The
identity of the group does not matter. This is the same as not
specifying an ordering at all. To control the ordering for every
variable, create an elimination group per variable, ordering them in
the desired order.

If the user is using one of the Schur solvers (``DENSE_SCHUR``,
``SPARSE_SCHUR``, ``ITERATIVE_SCHUR``) and chooses to specify an
ordering, it must have one important property. The lowest numbered
elimination group must form an independent set in the graph
corresponding to the Hessian, or in other words, no two parameter
blocks in in the first elimination group should co-occur in the same
residual block. For the best performance, this elimination group
should be as large as possible. For standard bundle adjustment
problems, this corresponds to the first elimination group containing
all the 3d points, and the second containing the all the cameras
parameter blocks.

If the user leaves the choice to Ceres, then the solver uses an
approximate maximum independent set algorithm to identify the first
elimination group [LiSaad]_.

.. _section-solver-options:

:class:`Solver::Options`
------------------------

.. class:: Solver::Options

  :class:`Solver::Options` controls the overall behavior of the
  solver. We list the various settings and their default values below.


.. member:: MinimizerType Solver::Options::minimizer_type

   Default: ``TRUST_REGION``

   Choose between ``LINE_SEARCH`` and ``TRUST_REGION`` algorithms. See
   :ref:`section-trust-region-methods` and
   :ref:`section-line-search-methods` for more details.

.. member:: LineSearchDirectionType Solver::Options::line_search_direction_type

   Default: ``LBFGS``

   Choices are ``STEEPEST_DESCENT``, ``NONLINEAR_CONJUGATE_GRADIENT``
   and ``LBFGS``.

.. member:: LineSearchType Solver::Options::line_search_type

   Default: ``ARMIJO``

   ``ARMIJO`` is the only choice right now.

.. member:: NonlinearConjugateGradientType Solver::Options::nonlinear_conjugate_gradient_type

   Default: ``FLETCHER_REEVES``

   Choices are ``FLETCHER_REEVES``, ``POLAK_RIBIRERE`` and
   ``HESTENES_STIEFEL``.

.. member:: int Solver::Options::max_lbfs_rank

   Default: 20

   The LBFGS hessian approximation is a low rank approximation to the
   inverse of the Hessian matrix. The rank of the approximation
   determines (linearly) the space and time complexity of using the
   approximation. Higher the rank, the better is the quality of the
   approximation. The increase in quality is however is bounded for a
   number of reasons.

     1. The method only uses secant information and not actual
        derivatives.

     2. The Hessian approximation is constrained to be positive
        definite.

   So increasing this rank to a large number will cost time and space
   complexity without the corresponding increase in solution
   quality. There are no hard and fast rules for choosing the maximum
   rank. The best choice usually requires some problem specific
   experimentation.

.. member:: TrustRegionStrategyType Solver::Options::trust_region_strategy_type

   Default: ``LEVENBERG_MARQUARDT``

   The trust region step computation algorithm used by
   Ceres. Currently ``LEVENBERG_MARQUARDT`` and ``DOGLEG`` are the two
   valid choices. See :ref:`section-levenberg-marquardt` and
   :ref:`section-dogleg` for more details.

.. member:: DoglegType Solver::Options::dogleg_type

   Default: ``TRADITIONAL_DOGLEG``

   Ceres supports two different dogleg strategies.
   ``TRADITIONAL_DOGLEG`` method by Powell and the ``SUBSPACE_DOGLEG``
   method described by [ByrdSchnabel]_ .  See :ref:`section-dogleg`
   for more details.

.. member:: bool Solver::Options::use_nonmonotonic_steps

   Default: ``false``

   Relax the requirement that the trust-region algorithm take strictly
   decreasing steps. See :ref:`section-non-monotonic-steps` for more
   details.

.. member:: int Solver::Options::max_consecutive_nonmonotonic_steps

   Default: ``5``

   The window size used by the step selection algorithm to accept
   non-monotonic steps.

.. member:: int Solver::Options::max_num_iterations

   Default: ``50``

   Maximum number of iterations for which the solver should run.

.. member:: double Solver::Options::max_solver_time_in_seconds

   Default: ``1e6``
   Maximum amount of time for which the solver should run.

.. member:: int Solver::Options::num_threads

   Default: ``1``

   Number of threads used by Ceres to evaluate the Jacobian.

.. member::  double Solver::Options::initial_trust_region_radius

   Default: ``1e4``

   The size of the initial trust region. When the
   ``LEVENBERG_MARQUARDT`` strategy is used, the reciprocal of this
   number is the initial regularization parameter.

.. member:: double Solver::Options::max_trust_region_radius

   Default: ``1e16``

   The trust region radius is not allowed to grow beyond this value.

.. member:: double Solver::Options::min_trust_region_radius

   Default: ``1e-32``

   The solver terminates, when the trust region becomes smaller than
   this value.

.. member:: double Solver::Options::min_relative_decrease

   Default: ``1e-3``

   Lower threshold for relative decrease before a trust-region step is
   acceped.

.. member:: double Solver::Options::lm_min_diagonal

   Default: ``1e6``

   The ``LEVENBERG_MARQUARDT`` strategy, uses a diagonal matrix to
   regularize the the trust region step. This is the lower bound on
   the values of this diagonal matrix.

.. member:: double Solver::Options::lm_max_diagonal

   Default:  ``1e32``

   The ``LEVENBERG_MARQUARDT`` strategy, uses a diagonal matrix to
   regularize the the trust region step. This is the upper bound on
   the values of this diagonal matrix.

.. member:: int Solver::Options::max_num_consecutive_invalid_steps

   Default: ``5``

   The step returned by a trust region strategy can sometimes be
   numerically invalid, usually because of conditioning
   issues. Instead of crashing or stopping the optimization, the
   optimizer can go ahead and try solving with a smaller trust
   region/better conditioned problem. This parameter sets the number
   of consecutive retries before the minimizer gives up.

.. member:: double Solver::Options::function_tolerance

   Default: ``1e-6``

   Solver terminates if

   .. math:: \frac{|\Delta \text{cost}|}{\text{cost} < \text{function_tolerance}}

   where, :math:`\Delta \text{cost}` is the change in objective function
   value (up or down) in the current iteration of Levenberg-Marquardt.

.. member:: double Solver::Options::gradient_tolerance

   Default: ``1e-10``

   Solver terminates if

   .. math:: \frac{\|g(x)\|_\infty}{\|g(x_0)\|_\infty} < \text{gradient_tolerance}

   where :math:`\|\cdot\|_\infty` refers to the max norm, and :math:`x_0` is
   the vector of initial parameter values.

.. member:: double Solver::Options::parameter_tolerance

   Default: ``1e-8``

   Solver terminates if

   .. math:: \|\Delta x\| < (\|x\| + \text{parameter_tolerance}) * \text{parameter_tolerance}

   where :math:`\Delta x` is the step computed by the linear solver in the
   current iteration of Levenberg-Marquardt.

.. member:: LinearSolverType Solver::Options::linear_solver_type

   Default: ``SPARSE_NORMAL_CHOLESKY`` / ``DENSE_QR``

   Type of linear solver used to compute the solution to the linear
   least squares problem in each iteration of the Levenberg-Marquardt
   algorithm. If Ceres is build with ``SuiteSparse`` linked in then
   the default is ``SPARSE_NORMAL_CHOLESKY``, it is ``DENSE_QR``
   otherwise.

.. member:: PreconditionerType Solver::Options::preconditioner_type

   Default: ``JACOBI``

   The preconditioner used by the iterative linear solver. The default
   is the block Jacobi preconditioner. Valid values are (in increasing
   order of complexity) ``IDENTITY``, ``JACOBI``, ``SCHUR_JACOBI``,
   ``CLUSTER_JACOBI`` and ``CLUSTER_TRIDIAGONAL``. See
   :ref:`section-preconditioner` for more details.

.. member:: SparseLinearAlgebraLibrary Solver::Options::sparse_linear_algebra_library

   Default:``SUITE_SPARSE``

   Ceres supports the use of two sparse linear algebra libraries,
   ``SuiteSparse``, which is enabled by setting this parameter to
   ``SUITE_SPARSE`` and ``CXSparse``, which can be selected by setting
   this parameter to ```CX_SPARSE``. ``SuiteSparse`` is a
   sophisticated and complex sparse linear algebra library and should
   be used in general. If your needs/platforms prevent you from using
   ``SuiteSparse``, consider using ``CXSparse``, which is a much
   smaller, easier to build library. As can be expected, its
   performance on large problems is not comparable to that of
   ``SuiteSparse``.

.. member:: int Solver::Options::num_linear_solver_threads

   Default: ``1``

   Number of threads used by the linear solver.

.. member:: bool Solver::Options::use_inner_iterations

   Default: ``false``

   Use a non-linear version of a simplified variable projection
   algorithm. Essentially this amounts to doing a further optimization
   on each Newton/Trust region step using a coordinate descent
   algorithm.  For more details, see :ref:`section-inner-iterations`.

.. member:: ParameterBlockOrdering*  Solver::Options::inner_iteration_ordering

   Default: ``NULL``

   If :member:`Solver::Options::use_inner_iterations` true, then the user has
   two choices.

   1. Let the solver heuristically decide which parameter blocks to
      optimize in each inner iteration. To do this, set
      :member:`Solver::Options::inner_iteration_ordering` to ``NULL``.

   2. Specify a collection of of ordered independent sets. The lower
      numbered groups are optimized before the higher number groups
      during the inner optimization phase. Each group must be an
      independent set.

   See :ref:`section-ordering` for more details.

.. member:: ParameterBlockOrdering* Solver::Options::linear_solver_ordering

   Default: ``NULL``

   An instance of the ordering object informs the solver about the
   desired order in which parameter blocks should be eliminated by the
   linear solvers. See section~\ref{sec:ordering`` for more details.

   If ``NULL``, the solver is free to choose an ordering that it
   thinks is best. Note: currently, this option only has an effect on
   the Schur type solvers, support for the ``SPARSE_NORMAL_CHOLESKY``
   solver is forth coming.

   See :ref:`section-ordering` for more details.

.. member:: bool Solver::Options::use_block_amd

   Default: ``true``

   By virtue of the modeling layer in Ceres being block oriented, all
   the matrices used by Ceres are also block oriented.  When doing
   sparse direct factorization of these matrices, the fill-reducing
   ordering algorithms can either be run on the block or the scalar
   form of these matrices. Running it on the block form exposes more
   of the super-nodal structure of the matrix to the Cholesky
   factorization routines. This leads to substantial gains in
   factorization performance. Setting this parameter to true, enables
   the use of a block oriented Approximate Minimum Degree ordering
   algorithm. Settings it to ``false``, uses a scalar AMD
   algorithm. This option only makes sense when using
   :member:`Solver::Options::sparse_linear_algebra_library` = ``SUITE_SPARSE``
   as it uses the ``AMD`` package that is part of ``SuiteSparse``.

.. member:: int Solver::Options::linear_solver_min_num_iterations

   Default: ``1``

   Minimum number of iterations used by the linear solver. This only
   makes sense when the linear solver is an iterative solver, e.g.,
   ``ITERATIVE_SCHUR`` or ``CGNR``.

.. member:: int Solver::Options::linear_solver_max_num_iterations

   Default: ``500``

   Minimum number of iterations used by the linear solver. This only
   makes sense when the linear solver is an iterative solver, e.g.,
   ``ITERATIVE_SCHUR`` or ``CGNR``.

.. member:: double Solver::Options::eta

   Default: ``1e-1``

   Forcing sequence parameter. The truncated Newton solver uses this
   number to control the relative accuracy with which the Newton step
   is computed. This constant is passed to
   ``ConjugateGradientsSolver`` which uses it to terminate the
   iterations when

   .. math:: \frac{Q_i - Q_{i-1}}{Q_i} < \frac{\eta}{i}

.. member:: bool Solver::Options::jacobi_scaling

   Default: ``true``

   ``true`` means that the Jacobian is scaled by the norm of its
   columns before being passed to the linear solver. This improves the
   numerical conditioning of the normal equations.

.. member:: LoggingType Solver::Options::logging_type

   Default: ``PER_MINIMIZER_ITERATION``

.. member:: bool Solver::Options::minimizer_progress_to_stdout

   Default: ``false``

   By default the :class:`Minimizer` progress is logged to ``STDERR``
   depending on the ``vlog`` level. If this flag is set to true, and
   :member:`Solver::Options::logging_type` is not ``SILENT``, the logging
   output is sent to ``STDOUT``.

.. member:: vector<int> Solver::Options::lsqp_iterations_to_dump

   Default: ``empty``

   List of iterations at which the optimizer should dump the linear
   least squares problem to disk. Useful for testing and
   benchmarking. If ``empty``, no problems are dumped.

.. member:: string Solver::Options::lsqp_dump_directory

   Default: ``/tmp``

   If :member:`Solver::Options::lsqp_iterations_to_dump` is non-empty, then
   this setting determines the directory to which the files containing
   the linear least squares problems are written to.

.. member:: DumpFormatType Solver::Options::lsqp_dump_format

   Default: ``TEXTFILE``

   The format in which linear least squares problems should be logged
   when :member:`Solver::Options::lsqp_iterations_to_dump` is non-empty.
   There are three options:

   * ``CONSOLE`` prints the linear least squares problem in a human
      readable format to ``stderr``. The Jacobian is printed as a
      dense matrix. The vectors :math:`D`, :math:`x` and :math:`f` are
      printed as dense vectors. This should only be used for small
      problems.

   * ``PROTOBUF`` Write out the linear least squares problem to the
     directory pointed to by :member:`Solver::Options::lsqp_dump_directory` as
     a protocol buffer. ``linear_least_squares_problems.h/cc``
     contains routines for loading these problems. For details on the
     on disk format used, see ``matrix.proto``. The files are named
     ``lm_iteration_???.lsqp``. This requires that ``protobuf`` be
     linked into Ceres Solver.

   * ``TEXTFILE`` Write out the linear least squares problem to the
     directory pointed to by member::`Solver::Options::lsqp_dump_directory` as
     text files which can be read into ``MATLAB/Octave``. The Jacobian
     is dumped as a text file containing :math:`(i,j,s)` triplets, the
     vectors :math:`D`, `x` and `f` are dumped as text files
     containing a list of their values.

   A ``MATLAB/Octave`` script called ``lm_iteration_???.m`` is also
   output, which can be used to parse and load the problem into memory.

.. member:: bool Solver::Options::check_gradients

   Default: ``false``

   Check all Jacobians computed by each residual block with finite
   differences. This is expensive since it involves computing the
   derivative by normal means (e.g. user specified, autodiff, etc),
   then also computing it using finite differences. The results are
   compared, and if they differ substantially, details are printed to
   the log.

.. member:: double Solver::Options::gradient_check_relative_precision

   Default: ``1e08``

   Precision to check for in the gradient checker. If the relative
   difference between an element in a Jacobian exceeds this number,
   then the Jacobian for that cost term is dumped.

.. member:: double Solver::Options::numeric_derivative_relative_step_size

   Default: ``1e-6``

   Relative shift used for taking numeric derivatives. For finite
   differencing, each dimension is evaluated at slightly shifted
   values, e.g., for forward differences, the numerical derivative is

   .. math::

     \delta &= numeric\_derivative\_relative\_step\_size\\
     \Delta f &= \frac{f((1 + \delta)  x) - f(x)}{\delta x}

   The finite differencing is done along each dimension. The reason to
   use a relative (rather than absolute) step size is that this way,
   numeric differentiation works for functions where the arguments are
   typically large (e.g. :math:`10^9`) and when the values are small
   (e.g. :math:`10^{-5}`). It is possible to construct *torture cases*
   which break this finite difference heuristic, but they do not come
   up often in practice.

.. member:: vector<IterationCallback> Solver::Options::callbacks

   Callbacks that are executed at the end of each iteration of the
   :class:`Minimizer`. They are executed in the order that they are
   specified in this vector. By default, parameter blocks are updated
   only at the end of the optimization, i.e when the
   :class:`Minimizer` terminates. This behavior is controlled by
   :member:`Solver::Options::update_state_every_variable`. If the user wishes
   to have access to the update parameter blocks when his/her
   callbacks are executed, then set
   :member:`Solver::Options::update_state_every_iteration` to true.

   The solver does NOT take ownership of these pointers.

.. member:: bool Solver::Options::update_state_every_iteration

   Default: ``false``

   Normally the parameter blocks are only updated when the solver
   terminates. Setting this to true update them in every
   iteration. This setting is useful when building an interactive
   application using Ceres and using an :class:`IterationCallback`.

.. member:: string Solver::Options::solver_log

   Default: ``empty``

   If non-empty, a summary of the execution of the solver is recorded
   to this file.  This file is used for recording and Ceres'
   performance. Currently, only the iteration number, total time and
   the objective function value are logged. The format of this file is
   expected to change over time as the performance evaluation
   framework is fleshed out.

:class:`ParameterBlockOrdering`
-------------------------------

.. class:: ParameterBlockOrdering

   TBD

:class:`IterationCallback`
--------------------------

.. class:: IterationCallback

   TBD

:class:`CRSMatrix`
------------------

.. class:: CRSMatrix

   TBD

:class:`Solver::Summary`
------------------------

.. class:: Solver::Summary

   TBD

:class:`GradientChecker`
------------------------

.. class:: GradientChecker





