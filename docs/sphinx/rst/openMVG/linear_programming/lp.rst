*************************
linear programming
*************************

Linear programming [LP]_ is a technique for the optimization of a linear objective function, subject to linear equality and linear inequality constraints such as:

.. math::

    \begin{align} & \text{maximize} && \mathbf{c}^\mathrm{T} \mathbf{x}\\
    & \text{subject to} && A \mathbf{x} \leq \mathbf{b} \\
    & \text{and} && \mathbf{x} \ge \mathbf{0} \end{align}

where ``x`` represents the vector of variables (to be determined), ``c`` and ``b`` are vectors of (known) coefficients, ``A`` is a (known) matrix of coefficients.

openMVG linear programming tools
---------------------------------

openMVG provides tools to:

- configure Linear programs (LP container),
- solve Linear Programs (convex or quasi convex ones).

It results in a collection of L infinity based solver for computer vision problems.

.. toctree::
  :maxdepth: 1
  
  linfinityCV.rst


openMVG linear program container
---------------------------------

openMVG provides a generic container for LP (Linear Programming problems) that can be dense of sparse.

.. code-block:: c++

  // Dense LP
  LP_Constraints
  // Sparse LP
  LP_Constraints_Sparse


It allows to embed:

- objective function ``c`` and the problem type (minimization or maximization),
- constraints (coefficients ``A``, Sign, objective value ``b``),
- bounds over ``x`` parameters (<=, =, >=).

openMVG linear program solvers
---------------------------------

openMVG provide access to different solvers (not exhaustive):

- OSI_CLP (COIN-OR) project,
- MOSEK commercial, free in a research context.

Those solver have been choosen due to the stability of their results and ability to handle large problems without numerical stability (LPSolve and GPLK have been discarded after extensive experiments).

I refer the reader to openMVG/src/openMVG/linearProgramming/linear_programming_test.cpp to know more.

openMVG linear programming module usage
-------------------------------------------

The linear programming module of openMVG can be used for:

- solve classical linear problem (optimization),
- test the feasibility of linear problem,
- optimize upper bound of feasible problem (quasi-convex linear programs).

**classical linear problem solving (optimization)**

Here an example of usage of the framework:

.. code-block:: c++
  
  // Setup the LP (fill A,b,c and the constraint over x)
  LP_Constraints cstraint;
  BuildLinearProblem(cstraint);

  // Solve the LP with the solver of your choice
  std::vector<double> vec_solution(2);
  #if OPENMVG_HAVE_MOSEK  
    MOSEK_SolveWrapper solver(2);
  #else
    OSI_CLP_SolverWrapper solver(2);
  #endif
  // Send constraint to the LP solver
  solver.setup(cstraint);

  // If LP have a solution
  if (solver.solve())
    // Get back estimated parameters
    solver.getSolution(vec_solution);

**Linear programming, feasible problem**

openMVG can be use also to test only the feasibility of a given LP


.. math::

    \begin{align} & \text{find} && \mathbf{x}\\
    & \text{subject to} && A \mathbf{x} \leq \mathbf{b} \\
    & \text{and} && \mathbf{x} \ge \mathbf{0} \end{align}

**Linear programming, quasi convex optimization**

openMVG used a lot of L infinity minimisation formulation.
Often the posed problems are quasi-convex and dependent of an external parameter that we are looking for (i.e the maximal re-projection error for which the set of contraint is still feasible).


Optimization of this upper bound parameter can be done by iterating over all the possible value or by using a bisection that reduce the search range at each iteration.

.. code-block:: c++

  Require: gammaLow, gammUp (Low and upper bound of the parameter to optimize)
  Require: the LP problem (cstraintBuilder)
  Ensure: the optimal gamma value, or return infeasibility of the contraints set.
  
  BisectionLP(
    LP_Solver & solver,
    ConstraintBuilder & cstraintBuilder,
    double gammaUp  = 1.0,  // Upper bound
    double gammaLow = 0.0,  // lower bound
    double eps      = 1e-8, // precision that stop dichotomy
    const int maxIteration = 20) // max number of iteration
  {
    ConstraintType constraint;
    do
    {
      ++k; // One more iteration

      double gamma = (gammaLow + gammaUp) / 2.0;

      //-- Setup constraint and solver
      cstraintBuilder.Build(gamma, constraint);
      solver.setup( constraint );
      
      //-- Solving
      bool bFeasible = solver.solve();

      //-- According feasibility update the corresponding bound
      //-> Feasible, update the upper bound
      //-> Not feasible, update the lower bound
      (bFeasible) ? gammaUp = gamma; : gammaLow = gamma;
      
    } while (k < maxIteration && gammaUp - gammaLow > eps);
  }



.. [LP] http://en.wikipedia.org/wiki/Linear_programming
