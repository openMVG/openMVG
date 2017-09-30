.. _chapter-solving_faqs:

.. default-domain:: cpp

.. cpp:namespace:: ceres

=======
Solving
=======

#. How do I evaluate the Jacobian for a solved problem?

   Using :func:`Problem::Evaluate`.

#. How do I choose the right linear solver?

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

     Ceres Solver v1.12.0 Solve Report
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
