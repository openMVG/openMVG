.. highlight:: c++

.. default-domain:: cpp

.. cpp:namespace:: ceres

.. _chapter-gradient_tutorial:

==================================
General Unconstrained Minimization
==================================

Ceres Solver besides being able to solve non-linear least squares
problem can also solve general unconstrained problems using just their
objective function value and gradients. In this chapter we will see
how to do this.

Rosenbrock's Function
=====================

Consider minimizing the famous `Rosenbrock's function
<http://en.wikipedia.org/wiki/Rosenbrock_function>`_ [#f1]_.

The simplest way to minimize is to define a templated functor to
evaluate the objective value of this function and then use Ceres
Solver's automatic differentiation to compute its derivatives.

We begin by defining a templated functor and then using
``AutoDiffFirstOrderFunction`` to construct an instance of the
``FirstOrderFunction`` interface. This is the object that is
responsible for computing the objective function value and the
gradient (if required). This is the analog of the
:class:`CostFunction` when defining non-linear least squares problems
in Ceres.

.. code::

  // f(x,y) = (1-x)^2 + 100(y - x^2)^2;
  struct Rosenbrock {
    template <typename T>
    bool operator()(const T* parameters, T* cost) const {
      const T x = parameters[0];
      const T y = parameters[1];
      cost[0] = (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x);
      return true;
    }

    static ceres::FirstOrderFunction* Create() {
      constexpr int kNumParameters = 2;
      return new ceres::AutoDiffFirstOrderFunction<Rosenbrock, kNumParameters>(
          new Rosenbrock);
    }
  };


Minimizing it then is a straightforward matter of constructing a
:class:`GradientProblem` object and calling :func:`Solve` on it.

.. code::

    double parameters[2] = {-1.2, 1.0};

    ceres::GradientProblem problem(Rosenbrock::Create());

    ceres::GradientProblemSolver::Options options;
    options.minimizer_progress_to_stdout = true;
    ceres::GradientProblemSolver::Summary summary;
    ceres::Solve(options, problem, parameters, &summary);

    std::cout << summary.FullReport() << "\n";

Executing this code results, solve the problem using limited memory
`BFGS
<http://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm>`_
algorithm.

.. code-block:: bash

       0: f: 2.420000e+01 d: 0.00e+00 g: 2.16e+02 h: 0.00e+00 s: 0.00e+00 e:  0 it: 1.19e-05 tt: 1.19e-05
       1: f: 4.280493e+00 d: 1.99e+01 g: 1.52e+01 h: 2.01e-01 s: 8.62e-04 e:  2 it: 7.30e-05 tt: 1.72e-04
       2: f: 3.571154e+00 d: 7.09e-01 g: 1.35e+01 h: 3.78e-01 s: 1.34e-01 e:  3 it: 1.60e-05 tt: 1.93e-04
       3: f: 3.440869e+00 d: 1.30e-01 g: 1.73e+01 h: 1.36e-01 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 1.97e-04
       4: f: 3.213597e+00 d: 2.27e-01 g: 1.55e+01 h: 1.06e-01 s: 4.59e-01 e:  1 it: 1.19e-06 tt: 2.00e-04
       5: f: 2.839723e+00 d: 3.74e-01 g: 1.05e+01 h: 1.34e-01 s: 5.24e-01 e:  1 it: 9.54e-07 tt: 2.03e-04
       6: f: 2.448490e+00 d: 3.91e-01 g: 1.29e+01 h: 3.04e-01 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 2.05e-04
       7: f: 1.943019e+00 d: 5.05e-01 g: 4.00e+00 h: 8.81e-02 s: 7.43e-01 e:  1 it: 9.54e-07 tt: 2.08e-04
       8: f: 1.731469e+00 d: 2.12e-01 g: 7.36e+00 h: 1.71e-01 s: 4.60e-01 e:  2 it: 2.15e-06 tt: 2.11e-04
       9: f: 1.503267e+00 d: 2.28e-01 g: 6.47e+00 h: 8.66e-02 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 2.14e-04
      10: f: 1.228331e+00 d: 2.75e-01 g: 2.00e+00 h: 7.70e-02 s: 7.90e-01 e:  1 it: 0.00e+00 tt: 2.16e-04
      11: f: 1.016523e+00 d: 2.12e-01 g: 5.15e+00 h: 1.39e-01 s: 3.76e-01 e:  2 it: 1.91e-06 tt: 2.25e-04
      12: f: 9.145773e-01 d: 1.02e-01 g: 6.74e+00 h: 7.98e-02 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 2.28e-04
      13: f: 7.508302e-01 d: 1.64e-01 g: 3.88e+00 h: 5.76e-02 s: 4.93e-01 e:  1 it: 9.54e-07 tt: 2.30e-04
      14: f: 5.832378e-01 d: 1.68e-01 g: 5.56e+00 h: 1.42e-01 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 2.33e-04
      15: f: 3.969581e-01 d: 1.86e-01 g: 1.64e+00 h: 1.17e-01 s: 1.00e+00 e:  1 it: 1.19e-06 tt: 2.36e-04
      16: f: 3.171557e-01 d: 7.98e-02 g: 3.84e+00 h: 1.18e-01 s: 3.97e-01 e:  2 it: 1.91e-06 tt: 2.39e-04
      17: f: 2.641257e-01 d: 5.30e-02 g: 3.27e+00 h: 6.14e-02 s: 1.00e+00 e:  1 it: 1.19e-06 tt: 2.42e-04
      18: f: 1.909730e-01 d: 7.32e-02 g: 5.29e-01 h: 8.55e-02 s: 6.82e-01 e:  1 it: 9.54e-07 tt: 2.45e-04
      19: f: 1.472012e-01 d: 4.38e-02 g: 3.11e+00 h: 1.20e-01 s: 3.47e-01 e:  2 it: 1.91e-06 tt: 2.49e-04
      20: f: 1.093558e-01 d: 3.78e-02 g: 2.97e+00 h: 8.43e-02 s: 1.00e+00 e:  1 it: 2.15e-06 tt: 2.52e-04
      21: f: 6.710346e-02 d: 4.23e-02 g: 1.42e+00 h: 9.64e-02 s: 8.85e-01 e:  1 it: 8.82e-06 tt: 2.81e-04
      22: f: 3.993377e-02 d: 2.72e-02 g: 2.30e+00 h: 1.29e-01 s: 4.63e-01 e:  2 it: 7.87e-06 tt: 2.96e-04
      23: f: 2.911794e-02 d: 1.08e-02 g: 2.55e+00 h: 6.55e-02 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 3.00e-04
      24: f: 1.457683e-02 d: 1.45e-02 g: 2.77e-01 h: 6.37e-02 s: 6.14e-01 e:  1 it: 1.19e-06 tt: 3.03e-04
      25: f: 8.577515e-03 d: 6.00e-03 g: 2.86e+00 h: 1.40e-01 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 3.06e-04
      26: f: 3.486574e-03 d: 5.09e-03 g: 1.76e-01 h: 1.23e-02 s: 1.00e+00 e:  1 it: 1.19e-06 tt: 3.09e-04
      27: f: 1.257570e-03 d: 2.23e-03 g: 1.39e-01 h: 5.08e-02 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 3.12e-04
      28: f: 2.783568e-04 d: 9.79e-04 g: 6.20e-01 h: 6.47e-02 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 3.15e-04
      29: f: 2.533399e-05 d: 2.53e-04 g: 1.68e-02 h: 1.98e-03 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 3.17e-04
      30: f: 7.591572e-07 d: 2.46e-05 g: 5.40e-03 h: 9.27e-03 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 3.20e-04
      31: f: 1.902460e-09 d: 7.57e-07 g: 1.62e-03 h: 1.89e-03 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 3.23e-04
      32: f: 1.003030e-12 d: 1.90e-09 g: 3.50e-05 h: 3.52e-05 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 3.26e-04
      33: f: 4.835994e-17 d: 1.00e-12 g: 1.05e-07 h: 1.13e-06 s: 1.00e+00 e:  1 it: 1.19e-06 tt: 3.34e-04
      34: f: 1.885250e-22 d: 4.84e-17 g: 2.69e-10 h: 1.45e-08 s: 1.00e+00 e:  1 it: 9.54e-07 tt: 3.37e-04

    Solver Summary (v 2.2.0-eigen-(3.4.0)-lapack-suitesparse-(7.1.0)-metis-(5.1.0)-acceleratesparse-eigensparse)

    Parameters                                  2
    Line search direction              LBFGS (20)
    Line search type                  CUBIC WOLFE


    Cost:
    Initial                          2.420000e+01
    Final                            1.955192e-27
    Change                           2.420000e+01

    Minimizer iterations                       36

    Time (in seconds):

      Cost evaluation                    0.000000 (0)
      Gradient & cost evaluation         0.000000 (44)
      Polynomial minimization            0.000061
    Total                                0.000438

    Termination:                      CONVERGENCE (Parameter tolerance reached. Relative step_norm: 1.890726e-11 <= 1.000000e-08.)

    Initial x: -1.2 y: 1
    Final   x: 1 y: 1




If you are unable to use automatic differentiation for some reason
(say because you need to call an external library), then you can
use numeric differentiation. In that case the functor is defined as
follows [#f2]_.

.. code::

  // f(x,y) = (1-x)^2 + 100(y - x^2)^2;
  struct Rosenbrock {
    bool operator()(const double* parameters, double* cost) const {
      const double x = parameters[0];
      const double y = parameters[1];
      cost[0] = (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x);
      return true;
    }

    static ceres::FirstOrderFunction* Create() {
      constexpr int kNumParameters = 2;
      return new ceres::NumericDiffFirstOrderFunction<Rosenbrock,
                                                      ceres::CENTRAL,
                                                      kNumParameters>(
          new Rosenbrock);
    }
  };

And finally, if you would rather compute the derivatives by hand (say
because the size of the parameter vector is too large to be
automatically differentiated). Then you should define an instance of
`FirstOrderFunction`, which is the analog of :class:`CostFunction` for
non-linear least squares problems [#f3]_.

.. code::

  // f(x,y) = (1-x)^2 + 100(y - x^2)^2;
  class Rosenbrock final  : public ceres::FirstOrderFunction {
    public:
      ~Rosenbrock() override {}

      bool Evaluate(const double* parameters,
                             double* cost,
                             double* gradient) const override {
         const double x = parameters[0];
         const double y = parameters[1];

         cost[0] = (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x);
         if (gradient) {
           gradient[0] = -2.0 * (1.0 - x) - 200.0 * (y - x * x) * 2.0 * x;
           gradient[1] = 200.0 * (y - x * x);
         }
        return true;
     }

     int NumParameters() const override { return 2; }
  };

.. rubric:: Footnotes

.. [#f1] `examples/rosenbrock.cc
   <https://ceres-solver.googlesource.com/ceres-solver/+/master/examples/rosenbrock.cc>`_

.. [#f2] `examples/rosenbrock_numeric_diff.cc
   <https://ceres-solver.googlesource.com/ceres-solver/+/master/examples/rosenbrock_numeric_diff.cc>`_

.. [#f3] `examples/rosenbrock_analytic_diff.cc
   <https://ceres-solver.googlesource.com/ceres-solver/+/master/examples/rosenbrock_analytic_diff.cc>`_
