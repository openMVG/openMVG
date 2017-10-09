.. _chapter-modeling_faqs:

.. default-domain:: cpp

.. cpp:namespace:: ceres

========
Modeling
========

#. Use analytical/automatic derivatives.

   This is the single most important piece of advice we can give to
   you. It is tempting to take the easy way out and use numeric
   differentiation. This is a bad idea. Numeric differentiation is
   slow, ill-behaved, hard to get right, and results in poor
   convergence behaviour.

   Ceres allows the user to define templated functors which will
   be automatically differentiated. For most situations this is enough
   and we recommend using this facility. In some cases the derivatives
   are simple enough or the performance considerations are such that
   the overhead of automatic differentiation is too much. In such
   cases, analytic derivatives are recommended.

   The use of numerical derivatives should be a measure of last
   resort, where it is simply not possible to write a templated
   implementation of the cost function.

   In many cases it is not possible to do analytic or automatic
   differentiation of the entire cost function, but it is generally
   the case that it is possible to decompose the cost function into
   parts that need to be numerically differentiated and parts that can
   be automatically or analytically differentiated.

   To this end, Ceres has extensive support for mixing analytic,
   automatic and numeric differentiation. See
   :class:`CostFunctionToFunctor`.

#. When using Quaternions,  consider using :class:`QuaternionParameterization`.

   `Quaternions <https://en.wikipedia.org/wiki/Quaternion>`_ are a
   four dimensional parameterization of the space of three dimensional
   rotations :math:`SO(3)`.  However, the :math:`SO(3)` is a three
   dimensional set, and so is the tangent space of a
   Quaternion. Therefore, it is sometimes (not always) benefecial to
   associate a local parameterization with parameter blocks
   representing a Quaternion. Assuming that the order of entries in
   your parameter block is :math:`w,x,y,z`, you can use
   :class:`QuaternionParameterization`.

   .. NOTE::

     If you are using `Eigen's Quaternion
     <http://eigen.tuxfamily.org/dox/classEigen_1_1Quaternion.html>`_
     object, whose layout is :math:`x,y,z,w`, then you should use
     :class:`EigenQuaternionParameterization`.


#. How do I solve problems with general linear & non-linear
   **inequality** constraints with Ceres Solver?

   Currently, Ceres Solver only supports upper and lower bounds
   constraints on the parameter blocks.

   A crude way of dealing with inequality constraints is have one or
   more of your cost functions check if the inequalities you are
   interested in are satisfied, and if not return false instead of
   true. This will prevent the solver from ever stepping into an
   infeasible region.

   This requires that the starting point for the optimization be a
   feasible point.  You also risk pre-mature convergence using this
   method.

#. How do I solve problems with general linear & non-linear **equality**
   constraints with Ceres Solver?

   There is no built in support in ceres for solving problems with
   equality constraints.  Currently, Ceres Solver only supports upper
   and lower bounds constraints on the parameter blocks.

   The trick described above for dealing with inequality
   constraints will **not** work for equality constraints.

#. How do I set one or more components of a parameter block constant?

   Using :class:`SubsetParameterization`.

#. Putting `Inverse Function Theorem
   <http://en.wikipedia.org/wiki/Inverse_function_theorem>`_ to use.

   Every now and then we have to deal with functions which cannot be
   evaluated analytically. Computing the Jacobian in such cases is
   tricky. A particularly interesting case is where the inverse of the
   function is easy to compute analytically. An example of such a
   function is the Coordinate transformation between the `ECEF
   <http://en.wikipedia.org/wiki/ECEF>`_ and the `WGS84
   <http://en.wikipedia.org/wiki/World_Geodetic_System>`_ where the
   conversion from WGS84 to ECEF is analytic, but the conversion
   back to WGS84 uses an iterative algorithm. So how do you compute the
   derivative of the ECEF to WGS84 transformation?

   One obvious approach would be to numerically
   differentiate the conversion function. This is not a good idea. For
   one, it will be slow, but it will also be numerically quite
   bad.

   Turns out you can use the `Inverse Function Theorem
   <http://en.wikipedia.org/wiki/Inverse_function_theorem>`_ in this
   case to compute the derivatives more or less analytically.

   The key result here is. If :math:`x = f^{-1}(y)`, and :math:`Df(x)`
   is the invertible Jacobian of :math:`f` at :math:`x`. Then the
   Jacobian :math:`Df^{-1}(y) = [Df(x)]^{-1}`, i.e., the Jacobian of
   the :math:`f^{-1}` is the inverse of the Jacobian of :math:`f`.

   Algorithmically this means that given :math:`y`, compute :math:`x =
   f^{-1}(y)` by whatever means you can. Evaluate the Jacobian of
   :math:`f` at :math:`x`. If the Jacobian matrix is invertible, then
   its inverse is the Jacobian of :math:`f^{-1}(y)` at  :math:`y`.

   One can put this into practice with the following code fragment.

   .. code-block:: c++

      Eigen::Vector3d ecef; // Fill some values
      // Iterative computation.
      Eigen::Vector3d lla = ECEFToLLA(ecef);
      // Analytic derivatives
      Eigen::Matrix3d lla_to_ecef_jacobian = LLAToECEFJacobian(lla);
      bool invertible;
      Eigen::Matrix3d ecef_to_lla_jacobian;
      lla_to_ecef_jacobian.computeInverseWithCheck(ecef_to_lla_jacobian, invertible);
