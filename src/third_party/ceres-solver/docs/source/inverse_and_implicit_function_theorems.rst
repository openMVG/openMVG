.. default-domain:: cpp

.. cpp:namespace:: ceres

.. _chapter-inverse_function_theorem:

==========================================
Using Inverse & Implicit Function Theorems
==========================================

Until now we have considered methods for computing derivatives that
work directly on the function being differentiated. However, this is
not always possible. For example, if the function can only be computed
via an iterative algorithm, or there is no explicit definition of the
function available.  In this section we will see how we can use two
basic results from calculus to get around these difficulties.


Inverse Function Theorem
========================

Suppose we wish to evaluate the derivative of a function :math:`f(x)`,
but evaluating :math:`f(x)` is not easy. Say it involves running an
iterative algorithm. You could try automatically differentiating the
iterative algorithm, but even if that is possible, it can become quite
expensive.

In some cases we get lucky, and computing the inverse of :math:`f(x)`
is an easy operation. In these cases, we can use the `Inverse Function
Theorem <http://en.wikipedia.org/wiki/Inverse_function_theorem>`_ to
compute the derivative exactly. Here is the key idea:

Assuming that :math:`y=f(x)` is continuously differentiable in a
neighborhood of a point :math:`x` and :math:`Df(x)` is the invertible
Jacobian of :math:`f` at :math:`x`, then by applying the chain rule to
the identity :math:`f^{-1}(f(x)) = x`, we have
:math:`Df^{-1}(f(x))Df(x) = I`, or :math:`Df^{-1}(y) = (Df(x))^{-1}`,
i.e., the Jacobian of :math:`f^{-1}` is the inverse of the Jacobian of
:math:`f`, or :math:`Df(x) = (Df^{-1}(y))^{-1}`.

For example, let :math:`f(x) = e^x`. Now of course we know that
:math:`Df(x) = e^x`, but let's try and compute it via the Inverse
Function Theorem. For :math:`x > 0`, we have :math:`f^{-1}(y) = \log
y`, so :math:`Df^{-1}(y) = \frac{1}{y}`, so :math:`Df(x) =
(Df^{-1}(y))^{-1} = y = e^x`.

You maybe wondering why the above is true. A smoothly differentiable
function in a small neighborhood is well approximated by a linear
function. Indeed this is a good way to think about the Jacobian, it is
the matrix that best approximates the function linearly. Once you do
that, it is straightforward to see that *locally* :math:`f^{-1}(y)` is
best approximated linearly by the inverse of the Jacobian of
:math:`f(x)`.

Let us now consider a more practical example.

Geodetic Coordinate System Conversion
-------------------------------------

When working with data related to the Earth, one can use two different
coordinate systems. The familiar (latitude, longitude, height)
Latitude-Longitude-Altitude coordinate system or the `ECEF
<http://en.wikipedia.org/wiki/ECEF>`_ coordinate systems. The former
is familiar but is not terribly convenient analytically. The latter is
a Cartesian system but not particularly intuitive. So systems that
process earth related data have to go back and forth between these
coordinate systems.

The conversion between the LLA and the ECEF coordinate system requires
a model of the Earth, the most commonly used one being `WGS84
<https://en.wikipedia.org/wiki/World_Geodetic_System#1984_version>`_.

Going from the spherical :math:`(\phi,\lambda,h)` to the ECEF
:math:`(x,y,z)` coordinates is easy.

.. math::

   \chi &= \sqrt{1 - e^2 \sin^2 \phi}

   X &= \left( \frac{a}{\chi} + h \right) \cos \phi \cos \lambda

   Y &= \left( \frac{a}{\chi} + h \right) \cos \phi \sin \lambda

   Z &= \left(\frac{a(1-e^2)}{\chi}  +h \right) \sin \phi

Here :math:`a` and :math:`e^2` are constants defined by `WGS84
<https://en.wikipedia.org/wiki/World_Geodetic_System#1984_version>`_.

Going from ECEF to LLA coordinates requires an iterative algorithm. So
to compute the derivative of the this transformation we invoke the
Inverse Function Theorem as follows:

.. code-block:: c++

   Eigen::Vector3d ecef; // Fill some values
   // Iterative computation.
   Eigen::Vector3d lla = ECEFToLLA(ecef);
   // Analytic derivatives
   Eigen::Matrix3d lla_to_ecef_jacobian = LLAToECEFJacobian(lla);
   bool invertible;
   Eigen::Matrix3d ecef_to_lla_jacobian;
   lla_to_ecef_jacobian.computeInverseWithCheck(ecef_to_lla_jacobian, invertible);


Implicit Function Theorem
=========================

Consider now the problem where we have two variables :math:`x \in
\mathbb{R}^m` and :math:`y \in \mathbb{R}^n` and a function
:math:`F:\mathbb{R}^m \times \mathbb{R}^n \rightarrow \mathbb{R}^n`
such that :math:`F(x,y) = 0` and we wish to calculate the Jacobian of
:math:`y` with respect to `x`. How do we do this?

If for a given value of :math:`(x,y)`, the partial Jacobian
:math:`D_2F(x,y)` is full rank, then the `Implicit Function Theorem
<https://en.wikipedia.org/wiki/Implicit_function_theorem>`_ tells us
that there exists a neighborhood of :math:`x` and a function :math:`G`
such :math:`y = G(x)` in this neighborhood. Differentiating
:math:`F(x,G(x)) = 0` gives us

.. math::

   D_1F(x,y) + D_2F(x,y)DG(x) &= 0

                        DG(x) &= -(D_2F(x,y))^{-1} D_1 F(x,y)

                        D y(x) &= -(D_2F(x,y))^{-1} D_1 F(x,y)

This means that we can compute the derivative of :math:`y` with
respect to :math:`x` by multiplying the Jacobian of :math:`F` w.r.t
:math:`x` by the inverse of the Jacobian of :math:`F` w.r.t :math:`y`.

Let's consider two examples.

Roots of a Polynomial
---------------------

The first example we consider is a classic. Let :math:`p(x) = a_0 +
a_1 x + \dots + a_n x^n` be a degree :math:`n` polynomial, and we wish
to compute the derivative of its roots with respect to its
coefficients. There is no closed form formula for computing the roots
of a general degree :math:`n` polynomial. `Galois
<https://en.wikipedia.org/wiki/%C3%89variste_Galois>`_ and `Abel
<https://en.wikipedia.org/wiki/Niels_Henrik_Abel>`_ proved that. There
are numerical algorithms like computing the eigenvalues of the
`Companion Matrix
<https://nhigham.com/2021/03/23/what-is-a-companion-matrix/>`_, but
differentiating an eigenvalue solver does not seem like fun. But the
Implicit Function Theorem offers us a simple path.

If :math:`x` is a root of :math:`p(x)`, then :math:`F(\mathbf{a}, x) =
a_0 + a_1 x + \dots + a_n x^n = 0`. So,

.. math::

   D_1 F(\mathbf{a}, x) &= [1, x, x^2, \dots, x^n]

   D_2 F(\mathbf{a}, x) &= \sum_{k=1}^n k a_k x^{k-1} = Dp(x)

        Dx(a) &= \frac{-1}{Dp(x)} [1, x, x^2, \dots, x^n]

Differentiating the Solution to an Optimization Problem
-------------------------------------------------------

Sometimes we are required to solve optimization problems inside
optimization problems, and this requires computing the derivative of
the optimal solution (or a fixed point) of an optimization problem
w.r.t its parameters.

Let :math:`\theta \in \mathbb{R}^m` be a vector, :math:`A(\theta) \in
\mathbb{R}^{k\times n}` be a matrix whose entries are a function of
:math:`\theta` with :math:`k \ge n` and let :math:`b \in \mathbb{R}^k`
be a constant vector, then consider the linear least squares problem:

.. math::

   x^* = \arg \min_x \|A(\theta) x - b\|_2^2

How do we compute :math:`D_\theta x^*(\theta)`?

One approach would be to observe that :math:`x^*(\theta) =
(A^\top(\theta)A(\theta))^{-1}A^\top(\theta)b` and then differentiate
this w.r.t :math:`\theta`. But this would require differentiating
through the inverse of the matrix
:math:`(A^\top(\theta)A(\theta))^{-1}`. Not exactly easy. Let's use
the Implicit Function Theorem instead.

The first step is to observe that :math:`x^*` satisfies the so called
*normal equations*.

.. math::

   A^\top(\theta)A(\theta)x^* - A^\top(\theta)b = 0

We will compute :math:`D_\theta x^*` column-wise, treating
:math:`A(\theta)` as a function of one coordinate (:math:`\theta_i`)
of :math:`\theta` at a time. So using the normal equations, let's
define :math:`F(\theta_i, x^*) = A^\top(\theta_i)A(\theta_i)x^* -
A^\top(\theta_i)b = 0`. Using which can now compute:

.. math::

   D_1F(\theta_i, x^*) &= D_{\theta_i}A^\top A + A^\top
   D_{\theta_i}Ax^* - D_{\theta_i} A^\top b = g_i

   D_2F(\theta_i, x^*) &= A^\top A

   Dx^*(\theta_i) & = -(A^\top A)^{-1} g_i

   Dx^*(\theta) & = -(A^\top A )^{-1} \left[g_1, \dots, g_m\right]

Observe that we only need to compute the inverse of :math:`A^\top A`,
to compute :math:`D x^*(\theta)`, which we needed anyways to compute
:math:`x^*`.
