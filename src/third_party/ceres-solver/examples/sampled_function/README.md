Sampled Functions
--

It is common to not have an analytical representation of the optimization
problem but rather a table of values at specific inputs. This commonly occurs
when working with images or when the functions in the problem are expensive to
evaluate. To use this data in an optimization problem we can use interpolation
to evaluate the function and derivatives at intermediate input values.

There are many libraries that implement a variety of interpolation schemes, but
it is difficult to use them in Ceres' automatic differentiation framework.
Instead, Ceres provides the ability to interpolate one and two dimensional data.

The one dimensional interpolation is based on the Cubic Hermite Spline. This
interpolation method requires knowledge of the function derivatives at the
control points, however we only know the function values. Consequently, we will
use the data to estimate derivatives at the control points. The choice of how to
compute the derivatives is not unique and Ceres uses the Catmull–Rom Spline
variant which uses `0.5 * (p_{k+1} - p_{k-1})` as the derivative for control
point `p_k.` This produces a first order differentiable interpolating
function. The two dimensional interpolation scheme is a generalization of the
one dimensional scheme where the interpolating function is assumed to be
separable in the two dimensions.

This example shows how to use interpolation schemes within the Ceres automatic
differentiation framework. This is a one dimensional example and the objective
function is to minimize `0.5 * f(x)^2` where `f(x) = (x - 4.5)^2`.

It is also possible to use analytical derivatives with the provided
interpolation schemes by using a `SizedCostFunction` and defining the
``Evaluate` function. For this example, the evaluate function would be:

```c++
bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
  if (jacobians == NULL || jacobians[0] == NULL)
    interpolator_.Evaluate(parameters[0][0], residuals);
  else
    interpolator_.Evaluate(parameters[0][0], residuals, jacobians[0]);

  return true;
}
```
