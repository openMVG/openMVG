.. default-domain:: cpp

.. cpp:namespace:: ceres

.. _chapter-interfacing_with_automatic_differentiation:

Interfacing with Automatic Differentiation
==========================================

Automatic differentiation is straightforward to use in cases where an
explicit expression for the cost function is available. But this is
not always possible. Often one has to interface with external routines
or data. In this chapter we will consider a number of different ways
of doing so.

To do this, we will consider the problem of finding parameters
:math:`\theta` and :math:`t` that solve an optimization problem of the
form:

.. math::
   \min & \quad \sum_i \left \|y_i - f\left (\|q_{i}\|^2\right) q_i
   \right \|^2\\
   \text{such that} & \quad q_i = R(\theta) x_i + t

Here, :math:`R` is a two dimensional rotation matrix parameterized
using the angle :math:`\theta` and :math:`t` is a two dimensional
vector. :math:`f` is an external distortion function.

We begin by considering the case, where we have a templated function
:code:`TemplatedComputeDistortion` that can compute the function
:math:`f`. Then the implementation of the corresponding residual
functor is straightforward and will look as follows:

.. code-block:: c++
   :emphasize-lines: 21

   template <typename T> T TemplatedComputeDistortion(const T r2) {
     const double k1 = 0.0082;
     const double k2 = 0.000023;
     return 1.0 + k1 * r2 + k2 * r2 * r2;
   }

   struct Affine2DWithDistortion {
     Affine2DWithDistortion(const double x_in[2], const double y_in[2]) {
       x[0] = x_in[0];
       x[1] = x_in[1];
       y[0] = y_in[0];
       y[1] = y_in[1];
     }

     template <typename T>
     bool operator()(const T* theta,
                     const T* t,
                     T* residuals) const {
       const T q_0 =  cos(theta[0]) * x[0] - sin(theta[0]) * x[1] + t[0];
       const T q_1 =  sin(theta[0]) * x[0] + cos(theta[0]) * x[1] + t[1];
       const T f = TemplatedComputeDistortion(q_0 * q_0 + q_1 * q_1);
       residuals[0] = y[0] - f * q_0;
       residuals[1] = y[1] - f * q_1;
       return true;
     }

     double x[2];
     double y[2];
   };

So far so good, but let us now consider three ways of defining
:math:`f` which are not directly amenable to being used with automatic
differentiation:

#. A non-templated function that evaluates its value.
#. A function that evaluates its value and derivative.
#. A function that is defined as a table of values to be interpolated.

We will consider them in turn below.

A function that returns its value
----------------------------------

Suppose we were given a function :code:`ComputeDistortionValue` with
the following signature

.. code-block:: c++

   double ComputeDistortionValue(double r2);

that computes the value of :math:`f`. The actual implementation of the
function does not matter. Interfacing this function with
:code:`Affine2DWithDistortion` is a three step process:

1. Wrap :code:`ComputeDistortionValue` into a functor
   :code:`ComputeDistortionValueFunctor`.
2. Numerically differentiate :code:`ComputeDistortionValueFunctor`
   using :class:`NumericDiffCostFunction` to create a
   :class:`CostFunction`.
3. Wrap the resulting :class:`CostFunction` object using
   :class:`CostFunctionToFunctor`. The resulting object is a functor
   with a templated :code:`operator()` method, which pipes the
   Jacobian computed by :class:`NumericDiffCostFunction` into the
   appropriate :code:`Jet` objects.

An implementation of the above three steps looks as follows:

.. code-block:: c++
   :emphasize-lines: 15,16,17,18,19,20, 29

   struct ComputeDistortionValueFunctor {
     bool operator()(const double* r2, double* value) const {
       *value = ComputeDistortionValue(r2[0]);
       return true;
     }
   };

   struct Affine2DWithDistortion {
     Affine2DWithDistortion(const double x_in[2], const double y_in[2]) {
       x[0] = x_in[0];
       x[1] = x_in[1];
       y[0] = y_in[0];
       y[1] = y_in[1];

       compute_distortion.reset(new ceres::CostFunctionToFunctor<1, 1>(
            new ceres::NumericDiffCostFunction<ComputeDistortionValueFunctor,
                                               ceres::CENTRAL,
                                               1,
                                               1>(
               new ComputeDistortionValueFunctor)));
     }

     template <typename T>
     bool operator()(const T* theta, const T* t, T* residuals) const {
       const T q_0 = cos(theta[0]) * x[0] - sin(theta[0]) * x[1] + t[0];
       const T q_1 = sin(theta[0]) * x[0] + cos(theta[0]) * x[1] + t[1];
       const T r2 = q_0 * q_0 + q_1 * q_1;
       T f;
       (*compute_distortion)(&r2, &f);
       residuals[0] = y[0] - f * q_0;
       residuals[1] = y[1] - f * q_1;
       return true;
     }

     double x[2];
     double y[2];
     std::unique_ptr<ceres::CostFunctionToFunctor<1, 1> > compute_distortion;
   };


A function that returns its value and derivative
------------------------------------------------

Now suppose we are given a function :code:`ComputeDistortionValue`
that is able to compute its value and optionally its Jacobian on demand
and has the following signature:

.. code-block:: c++

   void ComputeDistortionValueAndJacobian(double r2,
                                          double* value,
                                          double* jacobian);

Again, the actual implementation of the function does not
matter. Interfacing this function with :code:`Affine2DWithDistortion`
is a two step process:

1. Wrap :code:`ComputeDistortionValueAndJacobian` into a
   :class:`CostFunction` object which we call
   :code:`ComputeDistortionFunction`.
2. Wrap the resulting :class:`ComputeDistortionFunction` object using
   :class:`CostFunctionToFunctor`. The resulting object is a functor
   with a templated :code:`operator()` method, which pipes the
   Jacobian computed by :class:`NumericDiffCostFunction` into the
   appropriate :code:`Jet` objects.

The resulting code will look as follows:

.. code-block:: c++
   :emphasize-lines: 21,22, 33

   class ComputeDistortionFunction : public ceres::SizedCostFunction<1, 1> {
    public:
     virtual bool Evaluate(double const* const* parameters,
                           double* residuals,
                           double** jacobians) const {
       if (!jacobians) {
         ComputeDistortionValueAndJacobian(parameters[0][0], residuals, nullptr);
       } else {
         ComputeDistortionValueAndJacobian(parameters[0][0], residuals, jacobians[0]);
       }
       return true;
     }
   };

   struct Affine2DWithDistortion {
     Affine2DWithDistortion(const double x_in[2], const double y_in[2]) {
       x[0] = x_in[0];
       x[1] = x_in[1];
       y[0] = y_in[0];
       y[1] = y_in[1];
       compute_distortion.reset(
           new ceres::CostFunctionToFunctor<1, 1>(new ComputeDistortionFunction));
     }

     template <typename T>
     bool operator()(const T* theta,
                     const T* t,
                     T* residuals) const {
       const T q_0 =  cos(theta[0]) * x[0] - sin(theta[0]) * x[1] + t[0];
       const T q_1 =  sin(theta[0]) * x[0] + cos(theta[0]) * x[1] + t[1];
       const T r2 = q_0 * q_0 + q_1 * q_1;
       T f;
       (*compute_distortion)(&r2, &f);
       residuals[0] = y[0] - f * q_0;
       residuals[1] = y[1] - f * q_1;
       return true;
     }

     double x[2];
     double y[2];
     std::unique_ptr<ceres::CostFunctionToFunctor<1, 1> > compute_distortion;
   };


A function that is defined as a table of values
-----------------------------------------------

The third and final case we will consider is where the function
:math:`f` is defined as a table of values on the interval :math:`[0,
100)`, with a value for each integer.

.. code-block:: c++

   vector<double> distortion_values;

There are many ways of interpolating a table of values. Perhaps the
simplest and most common method is linear interpolation. But it is not
a great idea to use linear interpolation because the interpolating
function is not differentiable at the sample points.

A simple (well behaved) differentiable interpolation is the `Cubic
Hermite Spline
<http://en.wikipedia.org/wiki/Cubic_Hermite_spline>`_. Ceres Solver
ships with routines to perform Cubic & Bi-Cubic interpolation that is
automatic differentiation friendly.

Using Cubic interpolation requires first constructing a
:class:`Grid1D` object to wrap the table of values and then
constructing a :class:`CubicInterpolator` object using it.

The resulting code will look as follows:

.. code-block:: c++
   :emphasize-lines: 10,11,12,13, 24, 32,33

   struct Affine2DWithDistortion {
     Affine2DWithDistortion(const double x_in[2],
                            const double y_in[2],
                            const std::vector<double>& distortion_values) {
       x[0] = x_in[0];
       x[1] = x_in[1];
       y[0] = y_in[0];
       y[1] = y_in[1];

       grid.reset(new ceres::Grid1D<double, 1>(
           &distortion_values[0], 0, distortion_values.size()));
       compute_distortion.reset(
           new ceres::CubicInterpolator<ceres::Grid1D<double, 1> >(*grid));
     }

     template <typename T>
     bool operator()(const T* theta,
                     const T* t,
                     T* residuals) const {
       const T q_0 =  cos(theta[0]) * x[0] - sin(theta[0]) * x[1] + t[0];
       const T q_1 =  sin(theta[0]) * x[0] + cos(theta[0]) * x[1] + t[1];
       const T r2 = q_0 * q_0 + q_1 * q_1;
       T f;
       compute_distortion->Evaluate(r2, &f);
       residuals[0] = y[0] - f * q_0;
       residuals[1] = y[1] - f * q_1;
       return true;
     }

     double x[2];
     double y[2];
     std::unique_ptr<ceres::Grid1D<double, 1> > grid;
     std::unique_ptr<ceres::CubicInterpolator<ceres::Grid1D<double, 1> > > compute_distortion;
   };

In the above example we used :class:`Grid1D` and
:class:`CubicInterpolator` to interpolate a one dimensional table of
values. :class:`Grid2D` combined with :class:`CubicInterpolator` lets
the user to interpolate two dimensional tables of values. Note that
neither :class:`Grid1D` or :class:`Grid2D` are limited to scalar
valued functions, they also work with vector valued functions.
