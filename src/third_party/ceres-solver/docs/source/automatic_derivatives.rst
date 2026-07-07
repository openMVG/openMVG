.. default-domain:: cpp

.. cpp:namespace:: ceres

.. _chapter-automatic_derivatives:

=====================
Automatic Derivatives
=====================

We will now consider automatic differentiation. It is a technique that
can compute exact derivatives, fast, while requiring about the same
effort from the user as is needed to use numerical differentiation.

Don't believe me? Well here goes. The following code fragment
implements an automatically differentiated ``CostFunction`` for `Rat43
<http://www.itl.nist.gov/div898/strd/nls/data/ratkowsky3.shtml>`_.

.. code-block:: c++

  struct Rat43CostFunctor {
    Rat43CostFunctor(const double x, const double y) : x_(x), y_(y) {}

    template <typename T>
    bool operator()(const T* parameters, T* residuals) const {
      const T b1 = parameters[0];
      const T b2 = parameters[1];
      const T b3 = parameters[2];
      const T b4 = parameters[3];
      residuals[0] = b1 * pow(1.0 + exp(b2 -  b3 * x_), -1.0 / b4) - y_;
      return true;
    }

    private:
      const double x_;
      const double y_;
  };


  CostFunction* cost_function =
        new AutoDiffCostFunction<Rat43CostFunctor, 1, 4>(
          new Rat43CostFunctor(x, y));

Notice that compared to numeric differentiation, the only difference
when defining the functor for use with automatic differentiation is
the signature of the ``operator()``.

In the case of numeric differentiation it was

.. code-block:: c++

   bool operator()(const double* parameters, double* residuals) const;

and for automatic differentiation it is a templated function of the
form

.. code-block:: c++

   template <typename T> bool operator()(const T* parameters, T* residuals) const;


So what does this small change buy us? The following table compares
the time it takes to evaluate the residual and the Jacobian for
`Rat43` using various methods.

==========================   =========
CostFunction                 Time (ns)
==========================   =========
Rat43Analytic                      255
Rat43AnalyticOptimized              92
Rat43NumericDiffForward            262
Rat43NumericDiffCentral            517
Rat43NumericDiffRidders           3760
Rat43AutomaticDiff                 129
==========================   =========

We can get exact derivatives using automatic differentiation
(``Rat43AutomaticDiff``) with about the same effort that is required
to write the code for numeric differentiation but only :math:`40\%`
slower than hand optimized analytical derivatives.

So how does it work? For this we will have to learn about **Dual
Numbers** and **Jets** .


Dual Numbers & Jets
===================

.. NOTE::

   Reading this and the next section on implementing Jets is not
   necessary to use automatic differentiation in Ceres Solver. But
   knowing the basics of how Jets work is useful when debugging and
   reasoning about the performance of automatic differentiation.

Dual numbers are an extension of the real numbers analogous to complex
numbers: whereas complex numbers augment the reals by introducing an
imaginary unit :math:`\iota` such that :math:`\iota^2 = -1`, dual
numbers introduce an *infinitesimal* unit :math:`\epsilon` such that
:math:`\epsilon^2 = 0` . A dual number :math:`a + v\epsilon` has two
components, the *real* component :math:`a` and the *infinitesimal*
component :math:`v`.

Surprisingly, this simple change leads to a convenient method for
computing exact derivatives without needing to manipulate complicated
symbolic expressions.

For example, consider the function

.. math::

   f(x) = x^2 ,

Then,

.. math::

   \begin{align}
   f(10 + \epsilon) &= (10 + \epsilon)^2\\
            &= 100 + 20 \epsilon + \epsilon^2\\
            &= 100 + 20 \epsilon
   \end{align}

Observe that the coefficient of :math:`\epsilon` is :math:`Df(10) =
20`. Indeed this generalizes to functions which are not
polynomial. Consider an arbitrary differentiable function
:math:`f(x)`. Then we can evaluate :math:`f(x + \epsilon)` by
considering the Taylor expansion of :math:`f` near :math:`x`, which
gives us the infinite series

.. math::
   \begin{align}
   f(x + \epsilon) &= f(x) + Df(x) \epsilon + D^2f(x)
   \frac{\epsilon^2}{2} + D^3f(x) \frac{\epsilon^3}{6} + \cdots\\
   f(x + \epsilon) &= f(x) + Df(x) \epsilon
   \end{align}

Here we are using the fact that :math:`\epsilon^2 = 0`.

A `Jet <https://en.wikipedia.org/wiki/Jet_(mathematics)>`_ is a
:math:`n`-dimensional dual number, where we augment the real numbers
with :math:`n` infinitesimal units :math:`\epsilon_i,\ i=1,...,n` with
the property that :math:`\forall i, j\ :\epsilon_i\epsilon_j = 0`. Then
a Jet consists of a *real* part :math:`a` and a :math:`n`-dimensional
*infinitesimal* part :math:`\mathbf{v}`, i.e.,

.. math::
   x = a + \sum_j v_{j} \epsilon_j

The summation notation gets tedious, so we will also just write

.. math::
   x = a + \mathbf{v}.

where the :math:`\epsilon_i`'s are implicit. Then, using the same
Taylor series expansion used above, we can see that:

.. math::

  f(a + \mathbf{v}) = f(a) + Df(a) \mathbf{v}.

Similarly for a multivariate function
:math:`f:\mathbb{R}^{n}\rightarrow \mathbb{R}^m`, evaluated on
:math:`x_i = a_i + \mathbf{v}_i,\ \forall i = 1,...,n`:

.. math::
   f(x_1,..., x_n) = f(a_1, ..., a_n) + \sum_i D_i f(a_1, ..., a_n) \mathbf{v}_i

So if each :math:`\mathbf{v}_i = e_i` were the :math:`i^{\text{th}}`
standard basis vector, then, the above expression would simplify to

.. math::
   f(x_1,..., x_n) = f(a_1, ..., a_n) + \sum_i D_i f(a_1, ..., a_n) \epsilon_i

and we can extract the coordinates of the Jacobian by inspecting the
coefficients of :math:`\epsilon_i`.

Implementing Jets
-----------------

In order for the above to work in practice, we will need the ability
to evaluate an arbitrary function :math:`f` not just on real numbers
but also on dual numbers, but one does not usually evaluate functions
by evaluating their Taylor expansions,

This is where C++ templates and operator overloading comes into
play. The following code fragment has a simple implementation of a
``Jet`` and some operators/functions that operate on them.

.. code-block:: c++

   template<int N> struct Jet {
     double a;
     Eigen::Matrix<double, 1, N> v;
   };

   template<int N> Jet<N> operator+(const Jet<N>& f, const Jet<N>& g) {
     return Jet<N>(f.a + g.a, f.v + g.v);
   }

   template<int N> Jet<N> operator-(const Jet<N>& f, const Jet<N>& g) {
     return Jet<N>(f.a - g.a, f.v - g.v);
   }

   template<int N> Jet<N> operator*(const Jet<N>& f, const Jet<N>& g) {
     return Jet<N>(f.a * g.a, f.a * g.v + f.v * g.a);
   }

   template<int N> Jet<N> operator/(const Jet<N>& f, const Jet<N>& g) {
     return Jet<N>(f.a / g.a, f.v / g.a - f.a * g.v / (g.a * g.a));
   }

   template <int N> Jet<N> exp(const Jet<N>& f) {
     return Jet<T, N>(exp(f.a), exp(f.a) * f.v);
   }

   // This is a simple implementation for illustration purposes, the
   // actual implementation of pow requires careful handling of a number
   // of corner cases.
   template <int N>  Jet<N> pow(const Jet<N>& f, const Jet<N>& g) {
     return Jet<N>(pow(f.a, g.a),
                   g.a * pow(f.a, g.a - 1.0) * f.v +
                   pow(f.a, g.a) * log(f.a); * g.v);
   }


With these overloaded functions in hand, we can now call
``Rat43CostFunctor`` with an array of Jets instead of doubles. Putting
that together with appropriately initialized Jets allows us to compute
the Jacobian as follows:

.. code-block:: c++

  class Rat43Automatic : public ceres::SizedCostFunction<1,4> {
   public:
    Rat43Automatic(const Rat43CostFunctor* functor) : functor_(functor) {}
    virtual ~Rat43Automatic() {}
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const {
      // Just evaluate the residuals if Jacobians are not required.
      if (!jacobians) return (*functor_)(parameters[0], residuals);

      // Initialize the Jets
      ceres::Jet<4> jets[4];
      for (int i = 0; i < 4; ++i) {
        jets[i].a = parameters[0][i];
        jets[i].v.setZero();
        jets[i].v[i] = 1.0;
      }

      ceres::Jet<4> result;
      (*functor_)(jets, &result);

      // Copy the values out of the Jet.
      residuals[0] = result.a;
      for (int i = 0; i < 4; ++i) {
        jacobians[0][i] = result.v[i];
      }
      return true;
    }

   private:
    std::unique_ptr<const Rat43CostFunctor> functor_;
  };

Indeed, this is essentially how :class:`AutoDiffCostFunction` works.

Pitfalls
========

Automatic differentiation frees the user from the burden of computing
and reasoning about the symbolic expressions for the Jacobians, but
this freedom comes at a cost. For example consider the following
simple functor:

.. code-block:: c++

   struct Functor {
     template <typename T> bool operator()(const T* x, T* residual) const {
       residual[0] = 1.0 - sqrt(x[0] * x[0] + x[1] * x[1]);
       return true;
     }
   };

Looking at the code for the residual computation, one does not foresee
any problems. However, if we look at the analytical expressions for
the Jacobian:

.. math::

      y &= 1 - \sqrt{x_0^2 + x_1^2}\\
   D_1y &= -\frac{x_0}{\sqrt{x_0^2 + x_1^2}},\
   D_2y = -\frac{x_1}{\sqrt{x_0^2 + x_1^2}}

we find that it is an indeterminate form at :math:`x_0 = 0, x_1 =
0`.

There is no single solution to this problem. In some cases one needs
to reason explicitly about the points where indeterminacy may occur
and use alternate expressions using `L'Hôpital's rule
<https://en.wikipedia.org/wiki/L'H%C3%B4pital's_rule>`_ (see for
example some of the conversion routines in `rotation.h
<https://github.com/ceres-solver/ceres-solver/blob/master/include/ceres/rotation.h>`_. In
other cases, one may need to regularize the expressions to eliminate
these points.
