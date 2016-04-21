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
   with ``-logtostderr`` you can add ``-v=N`` for increasing values
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
   :class:`CostFunctionToFunctor`.

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
      requires that you build Ceres with support for ``SuiteSparse``,
      ``CXSparse`` or Eigen's sparse linear algebra libraries.

      If you do not have access to these libraries for whatever
      reason, ``ITERATIVE_SCHUR`` with ``SCHUR_JACOBI`` is an
      excellent alternative.

   5. For large bundle adjustment problems (a few thousand cameras or
      more) use the ``ITERATIVE_SCHUR`` solver. There are a number of
      preconditioner choices here. ``SCHUR_JACOBI`` offers an
      excellent balance of speed and accuracy. This is also the
      recommended option if you are solving medium sized problems for
      which ``DENSE_SCHUR`` is too slow but ``SuiteSparse`` is not
      available.

      .. NOTE::

        If you are solving small to medium sized problems, consider
        setting ``Solver::Options::use_explicit_schur_complement`` to
        ``true``, it can result in a substantial performance boost.

      If you are not satisfied with ``SCHUR_JACOBI``'s performance try
      ``CLUSTER_JACOBI`` and ``CLUSTER_TRIDIAGONAL`` in that
      order. They require that you have ``SuiteSparse``
      installed. Both of these preconditioners use a clustering
      algorithm. Use ``SINGLE_LINKAGE`` before ``CANONICAL_VIEWS``.

#. Use :func:`Solver::Summary::FullReport` to diagnose performance problems.

   When diagnosing Ceres performance issues - runtime and convergence,
   the first place to start is by looking at the output of
   ``Solver::Summary::FullReport``. Here is an example

   .. code-block:: bash

     ./bin/bundle_adjuster --input ../data/problem-16-22106-pre.txt

     iter      cost      cost_change  |gradient|   |step|    tr_ratio  tr_radius  ls_iter  iter_time  total_time
        0  4.185660e+06    0.00e+00    2.16e+07   0.00e+00   0.00e+00  1.00e+04       0    7.50e-02    3.58e-01
        1  1.980525e+05    3.99e+06    5.34e+06   2.40e+03   9.60e-01  3.00e+04       1    1.84e-01    5.42e-01
        2  5.086543e+04    1.47e+05    2.11e+06   1.01e+03   8.22e-01  4.09e+04       1    1.53e-01    6.95e-01
        3  1.859667e+04    3.23e+04    2.87e+05   2.64e+02   9.85e-01  1.23e+05       1    1.71e-01    8.66e-01
        4  1.803857e+04    5.58e+02    2.69e+04   8.66e+01   9.93e-01  3.69e+05       1    1.61e-01    1.03e+00
        5  1.803391e+04    4.66e+00    3.11e+02   1.02e+01   1.00e+00  1.11e+06       1    1.49e-01    1.18e+00

     Ceres Solver v1.11.0 Solve Report
     ----------------------------------
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
     Preprocessor                            0.283

       Residual evaluation                   0.061
       Jacobian evaluation                   0.361
       Linear solver                         0.382
     Minimizer                               0.895

     Postprocessor                           0.002
     Total                                   1.220

     Termination:                   NO_CONVERGENCE (Maximum number of iterations reached.)

  Let us focus on run-time performance. The relevant lines to look at
  are


   .. code-block:: bash

     Time (in seconds):
     Preprocessor                            0.283

       Residual evaluation                   0.061
       Jacobian evaluation                   0.361
       Linear solver                         0.382
     Minimizer                               0.895

     Postprocessor                           0.002
     Total                                   1.220


  Which tell us that of the total 1.2 seconds, about .3 seconds was
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
     Preprocessor                            0.051

       Residual evaluation                   0.053
       Jacobian evaluation                   0.344
       Linear solver                         0.372
     Minimizer                               0.854

     Postprocessor                           0.002
     Total                                   0.935



  The preprocessor time has gone down by more than 5.5x!.

Further Reading
===============

For a short but informative introduction to the subject we recommend
the booklet by [Madsen]_ . For a general introduction to non-linear
optimization we recommend [NocedalWright]_. [Bjorck]_ remains the
seminal reference on least squares problems. [TrefethenBau]_ book is
our favorite text on introductory numerical linear algebra. [Triggs]_
provides a thorough coverage of the bundle adjustment problem.
