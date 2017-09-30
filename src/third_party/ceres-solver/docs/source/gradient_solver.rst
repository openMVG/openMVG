.. highlight:: c++

.. default-domain:: cpp

.. _chapter-gradient_problem_solver:

==================================
General Unconstrained Minimization
==================================

Modeling
========

:class:`FirstOrderFunction`
---------------------------

.. class:: FirstOrderFunction

  Instances of :class:`FirstOrderFunction` implement the evaluation of
  a function and its gradient.

  .. code-block:: c++

   class FirstOrderFunction {
     public:
      virtual ~FirstOrderFunction() {}
      virtual bool Evaluate(const double* const parameters,
                            double* cost,
                            double* gradient) const = 0;
      virtual int NumParameters() const = 0;
   };

.. function:: bool FirstOrderFunction::Evaluate(const double* const parameters, double* cost, double* gradient) const

   Evaluate the cost/value of the function. If ``gradient`` is not
   ``NULL`` then evaluate the gradient too. If evaluation is
   successful return, ``true`` else return ``false``.

   ``cost`` guaranteed to be never ``NULL``, ``gradient`` can be ``NULL``.

.. function:: int FirstOrderFunction::NumParameters() const

   Number of parameters in the domain of the function.


:class:`GradientProblem`
------------------------

.. class:: GradientProblem

.. code-block:: c++

  class GradientProblem {
   public:
    explicit GradientProblem(FirstOrderFunction* function);
    GradientProblem(FirstOrderFunction* function,
                    LocalParameterization* parameterization);
    int NumParameters() const;
    int NumLocalParameters() const;
    bool Evaluate(const double* parameters, double* cost, double* gradient) const;
    bool Plus(const double* x, const double* delta, double* x_plus_delta) const;
  };

Instances of :class:`GradientProblem` represent general non-linear
optimization problems that must be solved using just the value of the
objective function and its gradient. Unlike the :class:`Problem`
class, which can only be used to model non-linear least squares
problems, instances of :class:`GradientProblem` not restricted in the
form of the objective function.

Structurally :class:`GradientProblem` is a composition of a
:class:`FirstOrderFunction` and optionally a
:class:`LocalParameterization`.

The :class:`FirstOrderFunction` is responsible for evaluating the cost
and gradient of the objective function.

The :class:`LocalParameterization` is responsible for going back and
forth between the ambient space and the local tangent space. When a
:class:`LocalParameterization` is not provided, then the tangent space
is assumed to coincide with the ambient Euclidean space that the
gradient vector lives in.

The constructor takes ownership of the :class:`FirstOrderFunction` and
:class:`LocalParamterization` objects passed to it.


.. function:: void Solve(const GradientProblemSolver::Options& options, const GradientProblem& problem, double* parameters, GradientProblemSolver::Summary* summary)

   Solve the given :class:`GradientProblem` using the values in
   ``parameters`` as the initial guess of the solution.


Solving
=======

:class:`GradientProblemSolver::Options`
---------------------------------------

.. class:: GradientProblemSolver::Options

   :class:`GradientProblemSolver::Options` controls the overall
   behavior of the solver. We list the various settings and their
   default values below.

.. function:: bool GradientProblemSolver::Options::IsValid(string* error) const

   Validate the values in the options struct and returns true on
   success. If there is a problem, the method returns false with
   ``error`` containing a textual description of the cause.

.. member:: LineSearchDirectionType GradientProblemSolver::Options::line_search_direction_type

   Default: ``LBFGS``

   Choices are ``STEEPEST_DESCENT``, ``NONLINEAR_CONJUGATE_GRADIENT``,
   ``BFGS`` and ``LBFGS``.

.. member:: LineSearchType GradientProblemSolver::Options::line_search_type

   Default: ``WOLFE``

   Choices are ``ARMIJO`` and ``WOLFE`` (strong Wolfe conditions).
   Note that in order for the assumptions underlying the ``BFGS`` and
   ``LBFGS`` line search direction algorithms to be guaranteed to be
   satisifed, the ``WOLFE`` line search should be used.

.. member:: NonlinearConjugateGradientType GradientProblemSolver::Options::nonlinear_conjugate_gradient_type

   Default: ``FLETCHER_REEVES``

   Choices are ``FLETCHER_REEVES``, ``POLAK_RIBIERE`` and
   ``HESTENES_STIEFEL``.

.. member:: int GradientProblemSolver::Options::max_lbfs_rank

   Default: 20

   The L-BFGS hessian approximation is a low rank approximation to the
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

.. member:: bool GradientProblemSolver::Options::use_approximate_eigenvalue_bfgs_scaling

   Default: ``false``

   As part of the ``BFGS`` update step / ``LBFGS`` right-multiply
   step, the initial inverse Hessian approximation is taken to be the
   Identity.  However, [Oren]_ showed that using instead :math:`I *
   \gamma`, where :math:`\gamma` is a scalar chosen to approximate an
   eigenvalue of the true inverse Hessian can result in improved
   convergence in a wide variety of cases.  Setting
   ``use_approximate_eigenvalue_bfgs_scaling`` to true enables this
   scaling in ``BFGS`` (before first iteration) and ``LBFGS`` (at each
   iteration).

   Precisely, approximate eigenvalue scaling equates to

   .. math:: \gamma = \frac{y_k' s_k}{y_k' y_k}

   With:

  .. math:: y_k = \nabla f_{k+1} - \nabla f_k
  .. math:: s_k = x_{k+1} - x_k

  Where :math:`f()` is the line search objective and :math:`x` the
  vector of parameter values [NocedalWright]_.

  It is important to note that approximate eigenvalue scaling does
  **not** *always* improve convergence, and that it can in fact
  *significantly* degrade performance for certain classes of problem,
  which is why it is disabled by default.  In particular it can
  degrade performance when the sensitivity of the problem to different
  parameters varies significantly, as in this case a single scalar
  factor fails to capture this variation and detrimentally downscales
  parts of the Jacobian approximation which correspond to
  low-sensitivity parameters. It can also reduce the robustness of the
  solution to errors in the Jacobians.

.. member:: LineSearchIterpolationType GradientProblemSolver::Options::line_search_interpolation_type

   Default: ``CUBIC``

   Degree of the polynomial used to approximate the objective
   function. Valid values are ``BISECTION``, ``QUADRATIC`` and
   ``CUBIC``.

.. member:: double GradientProblemSolver::Options::min_line_search_step_size

   The line search terminates if:

   .. math:: \|\Delta x_k\|_\infty < \text{min_line_search_step_size}

   where :math:`\|\cdot\|_\infty` refers to the max norm, and
   :math:`\Delta x_k` is the step change in the parameter values at
   the :math:`k`-th iteration.

.. member:: double GradientProblemSolver::Options::line_search_sufficient_function_decrease

   Default: ``1e-4``

   Solving the line search problem exactly is computationally
   prohibitive. Fortunately, line search based optimization algorithms
   can still guarantee convergence if instead of an exact solution,
   the line search algorithm returns a solution which decreases the
   value of the objective function sufficiently. More precisely, we
   are looking for a step size s.t.

   .. math:: f(\text{step_size}) \le f(0) + \text{sufficient_decrease} * [f'(0) * \text{step_size}]

   This condition is known as the Armijo condition.

.. member:: double GradientProblemSolver::Options::max_line_search_step_contraction

   Default: ``1e-3``

   In each iteration of the line search,

   .. math:: \text{new_step_size} \geq \text{max_line_search_step_contraction} * \text{step_size}

   Note that by definition, for contraction:

   .. math:: 0 < \text{max_step_contraction} < \text{min_step_contraction} < 1

.. member:: double GradientProblemSolver::Options::min_line_search_step_contraction

   Default: ``0.6``

   In each iteration of the line search,

   .. math:: \text{new_step_size} \leq \text{min_line_search_step_contraction} * \text{step_size}

   Note that by definition, for contraction:

   .. math:: 0 < \text{max_step_contraction} < \text{min_step_contraction} < 1

.. member:: int GradientProblemSolver::Options::max_num_line_search_step_size_iterations

   Default: ``20``

   Maximum number of trial step size iterations during each line
   search, if a step size satisfying the search conditions cannot be
   found within this number of trials, the line search will stop.

   As this is an 'artificial' constraint (one imposed by the user, not
   the underlying math), if ``WOLFE`` line search is being used, *and*
   points satisfying the Armijo sufficient (function) decrease
   condition have been found during the current search (in :math:`\leq`
   ``max_num_line_search_step_size_iterations``).  Then, the step size
   with the lowest function value which satisfies the Armijo condition
   will be returned as the new valid step, even though it does *not*
   satisfy the strong Wolfe conditions.  This behaviour protects
   against early termination of the optimizer at a sub-optimal point.

.. member:: int GradientProblemSolver::Options::max_num_line_search_direction_restarts

   Default: ``5``

   Maximum number of restarts of the line search direction algorithm
   before terminating the optimization. Restarts of the line search
   direction algorithm occur when the current algorithm fails to
   produce a new descent direction. This typically indicates a
   numerical failure, or a breakdown in the validity of the
   approximations used.

.. member:: double GradientProblemSolver::Options::line_search_sufficient_curvature_decrease

   Default: ``0.9``

   The strong Wolfe conditions consist of the Armijo sufficient
   decrease condition, and an additional requirement that the
   step size be chosen s.t. the *magnitude* ('strong' Wolfe
   conditions) of the gradient along the search direction
   decreases sufficiently. Precisely, this second condition
   is that we seek a step size s.t.

   .. math:: \|f'(\text{step_size})\| \leq \text{sufficient_curvature_decrease} * \|f'(0)\|

   Where :math:`f()` is the line search objective and :math:`f'()` is the derivative
   of :math:`f` with respect to the step size: :math:`\frac{d f}{d~\text{step size}}`.

.. member:: double GradientProblemSolver::Options::max_line_search_step_expansion

   Default: ``10.0``

   During the bracketing phase of a Wolfe line search, the step size
   is increased until either a point satisfying the Wolfe conditions
   is found, or an upper bound for a bracket containing a point
   satisfying the conditions is found.  Precisely, at each iteration
   of the expansion:

   .. math:: \text{new_step_size} \leq \text{max_step_expansion} * \text{step_size}

   By definition for expansion

   .. math:: \text{max_step_expansion} > 1.0

.. member:: int GradientProblemSolver::Options::max_num_iterations

   Default: ``50``

   Maximum number of iterations for which the solver should run.

.. member:: double GradientProblemSolver::Options::max_solver_time_in_seconds

   Default: ``1e6``
   Maximum amount of time for which the solver should run.

.. member:: double GradientProblemSolver::Options::function_tolerance

   Default: ``1e-6``

   Solver terminates if

   .. math:: \frac{|\Delta \text{cost}|}{\text{cost}} \leq \text{function_tolerance}

   where, :math:`\Delta \text{cost}` is the change in objective
   function value (up or down) in the current iteration of the line search.

.. member:: double GradientProblemSolver::Options::gradient_tolerance

   Default: ``1e-10``

   Solver terminates if

   .. math:: \|x - \Pi \boxplus(x, -g(x))\|_\infty \leq \text{gradient_tolerance}

   where :math:`\|\cdot\|_\infty` refers to the max norm, :math:`\Pi`
   is projection onto the bounds constraints and :math:`\boxplus` is
   Plus operation for the overall local parameterization associated
   with the parameter vector.

.. member:: double GradientProblemSolver::Options::parameter_tolerance

   Default: ``1e-8``

   Solver terminates if

   .. math:: \|\Delta x\| \leq (\|x\| + \text{parameter_tolerance}) * \text{parameter_tolerance}

   where :math:`\Delta x` is the step computed by the linear solver in
   the current iteration of the line search.

.. member:: LoggingType GradientProblemSolver::Options::logging_type

   Default: ``PER_MINIMIZER_ITERATION``

.. member:: bool GradientProblemSolver::Options::minimizer_progress_to_stdout

   Default: ``false``

   By default the :class:`Minimizer` progress is logged to ``STDERR``
   depending on the ``vlog`` level. If this flag is set to true, and
   :member:`GradientProblemSolver::Options::logging_type` is not
   ``SILENT``, the logging output is sent to ``STDOUT``.

   The progress display looks like

   .. code-block:: bash

      0: f: 2.317806e+05 d: 0.00e+00 g: 3.19e-01 h: 0.00e+00 s: 0.00e+00 e:  0 it: 2.98e-02 tt: 8.50e-02
      1: f: 2.312019e+05 d: 5.79e+02 g: 3.18e-01 h: 2.41e+01 s: 1.00e+00 e:  1 it: 4.54e-02 tt: 1.31e-01
      2: f: 2.300462e+05 d: 1.16e+03 g: 3.17e-01 h: 4.90e+01 s: 2.54e-03 e:  1 it: 4.96e-02 tt: 1.81e-01

   Here

   #. ``f`` is the value of the objective function.
   #. ``d`` is the change in the value of the objective function if
      the step computed in this iteration is accepted.
   #. ``g`` is the max norm of the gradient.
   #. ``h`` is the change in the parameter vector.
   #. ``s`` is the optimal step length computed by the line search.
   #. ``it`` is the time take by the current iteration.
   #. ``tt`` is the total time taken by the minimizer.

.. member:: vector<IterationCallback> GradientProblemSolver::Options::callbacks

   Callbacks that are executed at the end of each iteration of the
   :class:`Minimizer`. They are executed in the order that they are
   specified in this vector. See the documentation for
   :class:`IterationCallback` for more details.

   The solver does NOT take ownership of these pointers.


:class:`GradientProblemSolver::Summary`
---------------------------------------

.. class:: GradientProblemSolver::Summary

   Summary of the various stages of the solver after termination.

.. function:: string GradientProblemSolver::Summary::BriefReport() const

   A brief one line description of the state of the solver after
   termination.

.. function:: string GradientProblemSolver::Summary::FullReport() const

   A full multiline description of the state of the solver after
   termination.

.. function:: bool GradientProblemSolver::Summary::IsSolutionUsable() const

   Whether the solution returned by the optimization algorithm can be
   relied on to be numerically sane. This will be the case if
   `GradientProblemSolver::Summary:termination_type` is set to `CONVERGENCE`,
   `USER_SUCCESS` or `NO_CONVERGENCE`, i.e., either the solver
   converged by meeting one of the convergence tolerances or because
   the user indicated that it had converged or it ran to the maximum
   number of iterations or time.

.. member:: TerminationType GradientProblemSolver::Summary::termination_type

   The cause of the minimizer terminating.

.. member:: string GradientProblemSolver::Summary::message

   Reason why the solver terminated.

.. member:: double GradientProblemSolver::Summary::initial_cost

   Cost of the problem (value of the objective function) before the
   optimization.

.. member:: double GradientProblemSolver::Summary::final_cost

   Cost of the problem (value of the objective function) after the
   optimization.

.. member:: vector<IterationSummary> GradientProblemSolver::Summary::iterations

   :class:`IterationSummary` for each minimizer iteration in order.

.. member:: int num_cost_evaluations

   Number of times the cost (and not the gradient) was evaluated.

.. member:: int num_gradient_evaluations

   Number of times the gradient (and the cost) were evaluated.

.. member:: double GradientProblemSolver::Summary::total_time_in_seconds

   Time (in seconds) spent in the solver.

.. member:: double GradientProblemSolver::Summary::cost_evaluation_time_in_seconds

   Time (in seconds) spent evaluating the cost vector.

.. member:: double GradientProblemSolver::Summary::gradient_evaluation_time_in_seconds

   Time (in seconds) spent evaluating the gradient vector.

.. member:: int GradientProblemSolver::Summary::num_parameters

   Number of parameters in the problem.

.. member:: int GradientProblemSolver::Summary::num_local_parameters

   Dimension of the tangent space of the problem. This is different
   from :member:`GradientProblemSolver::Summary::num_parameters` if a
   :class:`LocalParameterization` object is used.

.. member:: LineSearchDirectionType GradientProblemSolver::Summary::line_search_direction_type

   Type of line search direction used.

.. member:: LineSearchType GradientProblemSolver::Summary::line_search_type

   Type of the line search algorithm used.

.. member:: LineSearchInterpolationType GradientProblemSolver::Summary::line_search_interpolation_type

   When performing line search, the degree of the polynomial used to
   approximate the objective function.

.. member:: NonlinearConjugateGradientType GradientProblemSolver::Summary::nonlinear_conjugate_gradient_type

   If the line search direction is `NONLINEAR_CONJUGATE_GRADIENT`,
   then this indicates the particular variant of non-linear conjugate
   gradient used.

.. member:: int GradientProblemSolver::Summary::max_lbfgs_rank

   If the type of the line search direction is `LBFGS`, then this
   indicates the rank of the Hessian approximation.
