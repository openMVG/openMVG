.. default-domain:: cpp

.. cpp:namespace:: ceres

.. _chapter-on_derivatives:

==============
On Derivatives
==============

Ceres Solver, like all gradient based optimization algorithms, depends
on being able to evaluate the objective function and its derivatives
at arbitrary points in its domain. Indeed, defining the objective
function and its `Jacobian
<https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant>`_ is
the principal task that the user is required to perform when solving
an optimization problem using Ceres Solver. The correct and efficient
computation of the Jacobian is the key to good performance.

Ceres Solver offers considerable flexibility in how the user can
provide derivatives to the solver. She can use:

#. :ref:`chapter-analytical_derivatives`: The user figures out the
   derivatives herself, by hand or using a tool like `Maple
   <https://www.maplesoft.com/products/maple/>`_ or `Mathematica
   <https://www.wolfram.com/mathematica/>`_, and implements them in a
   :class:`CostFunction`.
#. :ref:`chapter-numerical_derivatives`: Ceres numerically computes
   the derivative using finite differences.
#. :ref:`chapter-automatic_derivatives`: Ceres automatically computes
   the analytic derivative using C++ templates and operator
   overloading.

Which of these three approaches (alone or in combination) should be
used depends on the situation and the tradeoffs the user is willing to
make. Unfortunately, numerical optimization textbooks rarely discuss
these issues in detail and the user is left to her own devices.

The aim of this article is to fill this gap and describe each of these
three approaches in the context of Ceres Solver with sufficient detail
that the user can make an informed choice.

For the impatient amongst you, here is some high level advice:

#. Use :ref:`chapter-automatic_derivatives`.
#. In some cases it maybe worth using
   :ref:`chapter-analytical_derivatives`.
#. Avoid :ref:`chapter-numerical_derivatives`. Use it as a measure of
   last resort, mostly to interface with external libraries.

For the rest, read on.

.. toctree::
   :maxdepth: 1

   spivak_notation
   analytical_derivatives
   numerical_derivatives
   automatic_derivatives
   interfacing_with_autodiff
