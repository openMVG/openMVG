.. default-domain:: cpp

.. cpp:namespace:: ceres

.. _chapter-analytical_derivatives:

====================
Analytic Derivatives
====================

Consider the problem of fitting the following curve (`Rat43
<http://www.itl.nist.gov/div898/strd/nls/data/ratkowsky3.shtml>`_) to
data:

.. math::
  y = \frac{b_1}{(1+e^{b_2-b_3x})^{1/b_4}}

That is, given some data :math:`\{x_i, y_i\},\ \forall i=1,... ,n`,
determine parameters :math:`b_1, b_2, b_3` and :math:`b_4` that best
fit this data.

Which can be stated as the problem of finding the
values of :math:`b_1, b_2, b_3` and :math:`b_4` are the ones that
minimize the following objective function [#f1]_:

.. math::
   \begin{align}
   E(b_1, b_2, b_3, b_4)
   &= \sum_i f^2(b_1, b_2, b_3, b_4 ; x_i, y_i)\\
   &= \sum_i \left(\frac{b_1}{(1+e^{b_2-b_3x_i})^{1/b_4}} - y_i\right)^2\\
   \end{align}

To solve this problem using Ceres Solver, we need to define a
:class:`CostFunction` that computes the residual :math:`f` for a given
:math:`x` and :math:`y` and its derivatives with respect to
:math:`b_1, b_2, b_3` and :math:`b_4`.

Using elementary differential calculus, we can see that:

.. math::
  \begin{align}
  D_1 f(b_1, b_2, b_3, b_4; x,y) &= \frac{1}{(1+e^{b_2-b_3x})^{1/b_4}}\\
  D_2 f(b_1, b_2, b_3, b_4; x,y) &=
  \frac{-b_1e^{b_2-b_3x}}{b_4(1+e^{b_2-b_3x})^{1/b_4 + 1}} \\
  D_3 f(b_1, b_2, b_3, b_4; x,y) &=
  \frac{b_1xe^{b_2-b_3x}}{b_4(1+e^{b_2-b_3x})^{1/b_4 + 1}} \\
  D_4 f(b_1, b_2, b_3, b_4; x,y) & = \frac{b_1  \log\left(1+e^{b_2-b_3x}\right) }{b_4^2(1+e^{b_2-b_3x})^{1/b_4}}
  \end{align}

With these derivatives in hand, we can now implement the
:class:`CostFunction` as:

.. code-block:: c++

  class Rat43Analytic : public SizedCostFunction<1,4> {
     public:
       Rat43Analytic(const double x, const double y) : x_(x), y_(y) {}
       virtual ~Rat43Analytic() {}
       virtual bool Evaluate(double const* const* parameters,
                             double* residuals,
                             double** jacobians) const {
         const double b1 = parameters[0][0];
         const double b2 = parameters[0][1];
         const double b3 = parameters[0][2];
         const double b4 = parameters[0][3];

         residuals[0] = b1 *  pow(1 + exp(b2 -  b3 * x_), -1.0 / b4) - y_;

         if (!jacobians) return true;
         double* jacobian = jacobians[0];
         if (!jacobian) return true;

         jacobian[0] = pow(1 + exp(b2 - b3 * x_), -1.0 / b4);
         jacobian[1] = -b1 * exp(b2 - b3 * x_) *
                       pow(1 + exp(b2 - b3 * x_), -1.0 / b4 - 1) / b4;
         jacobian[2] = x_ * b1 * exp(b2 - b3 * x_) *
                       pow(1 + exp(b2 - b3 * x_), -1.0 / b4 - 1) / b4;
         jacobian[3] = b1 * log(1 + exp(b2 - b3 * x_)) *
                       pow(1 + exp(b2 - b3 * x_), -1.0 / b4) / (b4 * b4);
         return true;
       }

      private:
       const double x_;
       const double y_;
   };

This is tedious code, hard to read and with a lot of
redundancy. So in practice we will cache some sub-expressions to
improve its efficiency, which would give us something like:

.. code-block:: c++

  class Rat43AnalyticOptimized : public SizedCostFunction<1,4> {
     public:
       Rat43AnalyticOptimized(const double x, const double y) : x_(x), y_(y) {}
       virtual ~Rat43AnalyticOptimized() {}
       virtual bool Evaluate(double const* const* parameters,
                             double* residuals,
                             double** jacobians) const {
         const double b1 = parameters[0][0];
         const double b2 = parameters[0][1];
         const double b3 = parameters[0][2];
         const double b4 = parameters[0][3];

         const double t1 = exp(b2 -  b3 * x_);
         const double t2 = 1 + t1;
         const double t3 = pow(t2, -1.0 / b4);
         residuals[0] = b1 * t3 - y_;

         if (!jacobians) return true;
         double* jacobian = jacobians[0];
         if (!jacobian) return true;

         const double t4 = pow(t2, -1.0 / b4 - 1);
         jacobian[0] = t3;
         jacobian[1] = -b1 * t1 * t4 / b4;
         jacobian[2] = -x_ * jacobian[1];
         jacobian[3] = b1 * log(t2) * t3 / (b4 * b4);
         return true;
       }

     private:
       const double x_;
       const double y_;
   };

What is the difference in performance of these two implementations?

==========================   =========
CostFunction                 Time (ns)
==========================   =========
Rat43Analytic                      255
Rat43AnalyticOptimized              92
==========================   =========

``Rat43AnalyticOptimized`` is :math:`2.8` times faster than
``Rat43Analytic``.  This difference in run-time is not uncommon. To
get the best performance out of analytically computed derivatives, one
usually needs to optimize the code to account for common
sub-expressions.


When should you use analytical derivatives?
===========================================

#. The expressions are simple, e.g. mostly linear.

#. A computer algebra system like `Maple
   <https://www.maplesoft.com/products/maple/>`_ , `Mathematica
   <https://www.wolfram.com/mathematica/>`_, or `SymPy
   <http://www.sympy.org/en/index.html>`_ can be used to symbolically
   differentiate the objective function and generate the C++ to
   evaluate them.

#. Performance is of utmost concern and there is algebraic structure
   in the terms that you can exploit to get better performance than
   automatic differentiation.

   That said, getting the best performance out of analytical
   derivatives requires a non-trivial amount of work.  Before going
   down this path, it is useful to measure the amount of time being
   spent evaluating the Jacobian as a fraction of the total solve time
   and remember `Amdahl's Law
   <https://en.wikipedia.org/wiki/Amdahl's_law>`_ is your friend.

#. There is no other way to compute the derivatives, e.g. you
   wish to compute the derivative of the root of a polynomial:

   .. math::
     a_3(x,y)z^3 + a_2(x,y)z^2 + a_1(x,y)z + a_0(x,y) = 0


   with respect to :math:`x` and :math:`y`. This requires the use of
   the `Inverse Function Theorem
   <https://en.wikipedia.org/wiki/Inverse_function_theorem>`_

#. You love the chain rule and actually enjoy doing all the algebra by
   hand.


.. rubric:: Footnotes

.. [#f1] The notion of best fit depends on the choice of the objective
         function used to measure the quality of fit, which in turn
         depends on the underlying noise process which generated the
         observations. Minimizing the sum of squared differences is
         the right thing to do when the noise is `Gaussian
         <https://en.wikipedia.org/wiki/Normal_distribution>`_. In
         that case the optimal value of the parameters is the `Maximum
         Likelihood Estimate
         <https://en.wikipedia.org/wiki/Maximum_likelihood_estimation>`_.
