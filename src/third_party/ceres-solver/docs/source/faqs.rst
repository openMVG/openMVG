.. _chapter-tricks:

===================
FAQS, Tips & Tricks
===================

Answers to frequently asked questions, tricks of the trade and general
wisdom.

Building
========

#. Use `google-glog <http://code.google.com/p/google-glog>`_.

   Ceres has extensive support for logging detailed information about
   memory allocations and time consumed in various parts of the solve,
   internal error conditions etc. This is done logging using the
   `google-glog <http://code.google.com/p/google-glog>`_ library. We
   use it extensively to observe and analyze Ceres's
   performance. `google-glog <http://code.google.com/p/google-glog>`_
   allows you to control its behaviour from the command line `flags
   <http://google-glog.googlecode.com/svn/trunk/doc/glog.html>`_. Starting
   with ``-logtostdterr`` you can add ``-v=N`` for increasing values
   of ``N`` to get more and more verbose and detailed information
   about Ceres internals.

   In an attempt to reduce dependencies, it is tempting to use
   `miniglog` - a minimal implementation of the ``glog`` interface
   that ships with Ceres. This is a bad idea. ``miniglog`` was written
   primarily for building and using Ceres on Android because the
   current version of `google-glog
   <http://code.google.com/p/google-glog>`_ does not build using the
   NDK. It has worse performance than the full fledged glog library
   and is much harder to control and use.


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
   :class:`NumericDiffFunctor` and :class:`CostFunctionToFunctor`.

#. Putting `Inverse Function Theorem
   <http://en.wikipedia.org/wiki/Inverse_function_theorem>`_ to use.

   Every now and then we have to deal with functions which cannot be
   evaluated analytically. Computing the Jacobian in such cases is
   tricky. A particularly interesting case is where the inverse of the
   function is easy to compute analytically. An example of such a
   function is the Coordinate transformation between the `ECEF
   <http://en.wikipedia.org/wiki/ECEF>`_ and the `WGS84
   <http://en.wikipedia.org/wiki/World_Geodetic_System>`_ where the
   conversion from WGS84 from ECEF is analytic, but the conversion
   back to ECEF uses an iterative algorithm. So how do you compute the
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
   the inverse is the Jacobian of the inverse at :math:`y`.

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

#. When using Quaternions, use :class:`QuaternionParameterization`.

   TBD

#. How to choose a parameter block size?

   TBD

Solving
=======

#. Choosing a linear solver.

   When using the ``TRUST_REGION`` minimizer, the choice of linear
   solver is an important decision. It affects solution quality and
   runtime. Here is a simple way to reason about it.

   1. For small (a few hundred parameters) or dense problems use
      ``DENSE_QR``.

   2. For general sparse problems (i.e., the Jacobian matrix has a
      substantial number of zeros) use
      ``SPARSE_NORMAL_CHOLESKY``. This requires that you have
      ``SuiteSparse`` or ``CXSparse`` installed.

   3. For bundle adjustment problems with up to a hundred or so
      cameras, use ``DENSE_SCHUR``.

   4. For larger bundle adjustment problems with sparse Schur
      Complement/Reduced camera matrices use ``SPARSE_SCHUR``. This
      requires that you have ``SuiteSparse`` or ``CXSparse``
      installed.

   5. For large bundle adjustment problems (a few thousand cameras or
      more) use the ``ITERATIVE_SCHUR`` solver. There are a number of
      preconditioner choices here. ``SCHUR_JACOBI`` offers an
      excellent balance of speed and accuracy. This is also the
      recommended option if you are solving medium sized problems for
      which ``DENSE_SCHUR`` is too slow but ``SuiteSparse`` is not
      available.

      If you are not satisfied with ``SCHUR_JACOBI``'s performance try
      ``CLUSTER_JACOBI`` and ``CLUSTER_TRIDIAGONAL`` in that
      order. They require that you have ``SuiteSparse``
      installed. Both of these preconditioners use a clustering
      algorithm. Use ``SINGLE_LINKAGE`` before ``CANONICAL_VIEWS``.

#. Use `Solver::Summary::FullReport` to diagnose performance problems.

   When diagnosing Ceres performance issues - runtime and convergence,
   the first place to start is by looking at the output of
   ``Solver::Summary::FullReport``. Here is an example

   .. code-block:: bash

     ./bin/bundle_adjuster --input ../data/problem-16-22106-pre.txt


     0: f: 4.185660e+06 d: 0.00e+00 g: 2.16e+07 h: 0.00e+00 rho: 0.00e+00 mu: 1.00e+04 li:  0 it: 9.20e-02 tt: 3.35e-01
     1: f: 1.980525e+05 d: 3.99e+06 g: 5.34e+06 h: 2.40e+03 rho: 9.60e-01 mu: 3.00e+04 li:  1 it: 1.99e-01 tt: 5.34e-01
     2: f: 5.086543e+04 d: 1.47e+05 g: 2.11e+06 h: 1.01e+03 rho: 8.22e-01 mu: 4.09e+04 li:  1 it: 1.61e-01 tt: 6.95e-01
     3: f: 1.859667e+04 d: 3.23e+04 g: 2.87e+05 h: 2.64e+02 rho: 9.85e-01 mu: 1.23e+05 li:  1 it: 1.63e-01 tt: 8.58e-01
     4: f: 1.803857e+04 d: 5.58e+02 g: 2.69e+04 h: 8.66e+01 rho: 9.93e-01 mu: 3.69e+05 li:  1 it: 1.62e-01 tt: 1.02e+00
     5: f: 1.803391e+04 d: 4.66e+00 g: 3.11e+02 h: 1.02e+01 rho: 1.00e+00 mu: 1.11e+06 li:  1 it: 1.61e-01 tt: 1.18e+00

     Ceres Solver Report
     -------------------
                                          Original                  Reduced
     Parameter blocks                        22122                    22122
     Parameters                              66462                    66462
     Residual blocks                         83718                    83718
     Residual                               167436                   167436

     Minimizer                        TRUST_REGION

     Sparse linear algebra library    SUITE_SPARSE
     Trust region strategy     LEVENBERG_MARQUARDT

                                              Given                     Used
     Linear solver                    SPARSE_SCHUR             SPARSE_SCHUR
     Threads                                     1                        1
     Linear solver threads                       1                        1
     Linear solver ordering              AUTOMATIC                22106, 16

     Cost:
     Initial                          4.185660e+06
     Final                            1.803391e+04
     Change                           4.167626e+06

     Minimizer iterations                        5
     Successful steps                            5
     Unsuccessful steps                          0

     Time (in seconds):
     Preprocessor                            0.243

       Residual evaluation                   0.053
       Jacobian evaluation                   0.435
       Linear solver                         0.371
     Minimizer                               0.940

     Postprocessor                           0.002
     Total                                   1.221

     Termination:                   NO_CONVERGENCE (Maximum number of iterations reached.)

  Let us focus on run-time performance. The relevant lines to look at
  are


   .. code-block:: bash

     Time (in seconds):
     Preprocessor                            0.243

       Residual evaluation                   0.053
       Jacobian evaluation                   0.435
       Linear solver                         0.371
     Minimizer                               0.940

     Postprocessor                           0.002
     Total                                   1.221

  Which tell us that of the total 1.2 seconds, about .4 seconds was
  spent in the linear solver and the rest was mostly spent in
  preprocessing and jacobian evaluation.

  The preprocessing seems particularly expensive. Looking back at the
  report, we observe

   .. code-block:: bash

     Linear solver ordering              AUTOMATIC                22106, 16

  Which indicates that we are using automatic ordering for the
  ``SPARSE_SCHUR`` solver. This can be expensive at times. A straight
  forward way to deal with this is to give the ordering manually. For
  ``bundle_adjuster`` this can be done by passing the flag
  ``-ordering=user``. Doing so and looking at the timing block of the
  full report gives us

   .. code-block:: bash

     Time (in seconds):
     Preprocessor                            0.058

       Residual evaluation                   0.050
       Jacobian evaluation                   0.416
       Linear solver                         0.360
     Minimizer                               0.903

     Postprocessor                           0.002
     Total                                   0.998

  The preprocessor time has gone down by more than 4x!.

Further Reading
===============

For a short but informative introduction to the subject we recommend
the booklet by [Madsen]_ . For a general introduction to non-linear
optimization we recommend [NocedalWright]_. [Bjorck]_ remains the
seminal reference on least squares problems. [TrefethenBau]_ book is
our favorite text on introductory numerical linear algebra. [Triggs]_
provides a thorough coverage of the bundle adjustment problem.
