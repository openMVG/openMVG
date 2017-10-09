
.. default-domain:: cpp

.. cpp:namespace:: ceres

.. _chapter-nnls_covariance:

=====================
Covariance Estimation
=====================

Introduction
============

One way to assess the quality of the solution returned by a non-linear
least squares solver is to analyze the covariance of the solution.

Let us consider the non-linear regression problem

.. math::  y = f(x) + N(0, I)

i.e., the observation :math:`y` is a random non-linear function of the
independent variable :math:`x` with mean :math:`f(x)` and identity
covariance. Then the maximum likelihood estimate of :math:`x` given
observations :math:`y` is the solution to the non-linear least squares
problem:

.. math:: x^* = \arg \min_x \|f(x)\|^2

And the covariance of :math:`x^*` is given by

.. math:: C(x^*) = \left(J'(x^*)J(x^*)\right)^{-1}

Here :math:`J(x^*)` is the Jacobian of :math:`f` at :math:`x^*`. The
above formula assumes that :math:`J(x^*)` has full column rank.

If :math:`J(x^*)` is rank deficient, then the covariance matrix :math:`C(x^*)`
is also rank deficient and is given by the Moore-Penrose pseudo inverse.

.. math:: C(x^*) =  \left(J'(x^*)J(x^*)\right)^{\dagger}

Note that in the above, we assumed that the covariance matrix for
:math:`y` was identity. This is an important assumption. If this is
not the case and we have

.. math:: y = f(x) + N(0, S)

Where :math:`S` is a positive semi-definite matrix denoting the
covariance of :math:`y`, then the maximum likelihood problem to be
solved is

.. math:: x^* = \arg \min_x f'(x) S^{-1} f(x)

and the corresponding covariance estimate of :math:`x^*` is given by

.. math:: C(x^*) = \left(J'(x^*) S^{-1} J(x^*)\right)^{-1}

So, if it is the case that the observations being fitted to have a
covariance matrix not equal to identity, then it is the user's
responsibility that the corresponding cost functions are correctly
scaled, e.g. in the above case the cost function for this problem
should evaluate :math:`S^{-1/2} f(x)` instead of just :math:`f(x)`,
where :math:`S^{-1/2}` is the inverse square root of the covariance
matrix :math:`S`.

Gauge Invariance
================

In structure from motion (3D reconstruction) problems, the
reconstruction is ambiguous upto a similarity transform. This is
known as a *Gauge Ambiguity*. Handling Gauges correctly requires the
use of SVD or custom inversion algorithms. For small problems the
user can use the dense algorithm. For more details see the work of
Kanatani & Morris [KanataniMorris]_.


:class:`Covariance`
===================

:class:`Covariance` allows the user to evaluate the covariance for a
non-linear least squares problem and provides random access to its
blocks. The computation assumes that the cost functions compute
residuals such that their covariance is identity.

Since the computation of the covariance matrix requires computing the
inverse of a potentially large matrix, this can involve a rather large
amount of time and memory. However, it is usually the case that the
user is only interested in a small part of the covariance
matrix. Quite often just the block diagonal. :class:`Covariance`
allows the user to specify the parts of the covariance matrix that she
is interested in and then uses this information to only compute and
store those parts of the covariance matrix.

Rank of the Jacobian
====================

As we noted above, if the Jacobian is rank deficient, then the inverse
of :math:`J'J` is not defined and instead a pseudo inverse needs to be
computed.

The rank deficiency in :math:`J` can be *structural* -- columns
which are always known to be zero or *numerical* -- depending on the
exact values in the Jacobian.

Structural rank deficiency occurs when the problem contains parameter
blocks that are constant. This class correctly handles structural rank
deficiency like that.

Numerical rank deficiency, where the rank of the matrix cannot be
predicted by its sparsity structure and requires looking at its
numerical values is more complicated. Here again there are two
cases.

  a. The rank deficiency arises from overparameterization. e.g., a
     four dimensional quaternion used to parameterize :math:`SO(3)`,
     which is a three dimensional manifold. In cases like this, the
     user should use an appropriate
     :class:`LocalParameterization`. Not only will this lead to better
     numerical behaviour of the Solver, it will also expose the rank
     deficiency to the :class:`Covariance` object so that it can
     handle it correctly.

  b. More general numerical rank deficiency in the Jacobian requires
     the computation of the so called Singular Value Decomposition
     (SVD) of :math:`J'J`. We do not know how to do this for large
     sparse matrices efficiently. For small and moderate sized
     problems this is done using dense linear algebra.


:class:`Covariance::Options`

.. class:: Covariance::Options

.. member:: int Covariance::Options::num_threads

   Default: ``1``

   Number of threads to be used for evaluating the Jacobian and
   estimation of covariance.

.. member:: SparseLinearAlgebraLibraryType Covariance::Options::sparse_linear_algebra_library_type

   Default: ``SUITE_SPARSE`` Ceres Solver is built with support for
   `SuiteSparse <http://faculty.cse.tamu.edu/davis/suitesparse.html>`_
   and ``EIGEN_SPARSE`` otherwise. Note that ``EIGEN_SPARSE`` is
   always available.

.. member:: CovarianceAlgorithmType Covariance::Options::algorithm_type

   Default: ``SPARSE_QR``

   Ceres supports two different algorithms for covariance estimation,
   which represent different tradeoffs in speed, accuracy and
   reliability.

   1. ``SPARSE_QR`` uses the sparse QR factorization algorithm to
      compute the decomposition

       .. math::

          QR &= J\\
          \left(J^\top J\right)^{-1} &= \left(R^\top R\right)^{-1}

      The speed of this algorithm depends on the sparse linear algebra
      library being used. ``Eigen``'s sparse QR factorization is a
      moderately fast algorithm suitable for small to medium sized
      matrices. For best performance we recommend using
      ``SuiteSparseQR`` which is enabled by setting
      :member:`Covaraince::Options::sparse_linear_algebra_library_type`
      to ``SUITE_SPARSE``.

      Neither ``SPARSE_QR`` cannot compute the covariance if the
      Jacobian is rank deficient.


   2. ``DENSE_SVD`` uses ``Eigen``'s ``JacobiSVD`` to perform the
      computations. It computes the singular value decomposition

      .. math::   U S V^\top = J

      and then uses it to compute the pseudo inverse of J'J as

      .. math::   (J'J)^{\dagger} = V  S^{\dagger}  V^\top

      It is an accurate but slow method and should only be used for
      small to moderate sized problems. It can handle full-rank as
      well as rank deficient Jacobians.


.. member:: int Covariance::Options::min_reciprocal_condition_number

   Default: :math:`10^{-14}`

   If the Jacobian matrix is near singular, then inverting :math:`J'J`
   will result in unreliable results, e.g, if

   .. math::

     J = \begin{bmatrix}
         1.0& 1.0 \\
         1.0& 1.0000001
         \end{bmatrix}

   which is essentially a rank deficient matrix, we have

   .. math::

     (J'J)^{-1} = \begin{bmatrix}
                  2.0471e+14&  -2.0471e+14 \\
                  -2.0471e+14   2.0471e+14
                  \end{bmatrix}


   This is not a useful result. Therefore, by default
   :func:`Covariance::Compute` will return ``false`` if a rank
   deficient Jacobian is encountered. How rank deficiency is detected
   depends on the algorithm being used.

   1. ``DENSE_SVD``

      .. math:: \frac{\sigma_{\text{min}}}{\sigma_{\text{max}}}  < \sqrt{\text{min_reciprocal_condition_number}}

      where :math:`\sigma_{\text{min}}` and
      :math:`\sigma_{\text{max}}` are the minimum and maxiumum
      singular values of :math:`J` respectively.

   2. ``SPARSE_QR``

       .. math:: \operatorname{rank}(J) < \operatorname{num\_col}(J)

       Here :math:`\operatorname{rank}(J)` is the estimate of the rank
       of :math:`J` returned by the sparse QR factorization
       algorithm. It is a fairly reliable indication of rank
       deficiency.

.. member:: int Covariance::Options::null_space_rank

    When using ``DENSE_SVD``, the user has more control in dealing
    with singular and near singular covariance matrices.

    As mentioned above, when the covariance matrix is near singular,
    instead of computing the inverse of :math:`J'J`, the Moore-Penrose
    pseudoinverse of :math:`J'J` should be computed.

    If :math:`J'J` has the eigen decomposition :math:`(\lambda_i,
    e_i)`, where :math:`\lambda_i` is the :math:`i^\textrm{th}`
    eigenvalue and :math:`e_i` is the corresponding eigenvector, then
    the inverse of :math:`J'J` is

    .. math:: (J'J)^{-1} = \sum_i \frac{1}{\lambda_i} e_i e_i'

    and computing the pseudo inverse involves dropping terms from this
    sum that correspond to small eigenvalues.

    How terms are dropped is controlled by
    `min_reciprocal_condition_number` and `null_space_rank`.

    If `null_space_rank` is non-negative, then the smallest
    `null_space_rank` eigenvalue/eigenvectors are dropped irrespective
    of the magnitude of :math:`\lambda_i`. If the ratio of the
    smallest non-zero eigenvalue to the largest eigenvalue in the
    truncated matrix is still below min_reciprocal_condition_number,
    then the `Covariance::Compute()` will fail and return `false`.

    Setting `null_space_rank = -1` drops all terms for which

    .. math::  \frac{\lambda_i}{\lambda_{\textrm{max}}} < \textrm{min_reciprocal_condition_number}

    This option has no effect on ``SPARSE_QR``.

.. member:: bool Covariance::Options::apply_loss_function

   Default: `true`

   Even though the residual blocks in the problem may contain loss
   functions, setting ``apply_loss_function`` to false will turn off
   the application of the loss function to the output of the cost
   function and in turn its effect on the covariance.

.. class:: Covariance

   :class:`Covariance::Options` as the name implies is used to control
   the covariance estimation algorithm. Covariance estimation is a
   complicated and numerically sensitive procedure. Please read the
   entire documentation for :class:`Covariance::Options` before using
   :class:`Covariance`.

.. function:: bool Covariance::Compute(const vector<pair<const double*, const double*> >& covariance_blocks, Problem* problem)

   Compute a part of the covariance matrix.

   The vector ``covariance_blocks``, indexes into the covariance
   matrix block-wise using pairs of parameter blocks. This allows the
   covariance estimation algorithm to only compute and store these
   blocks.

   Since the covariance matrix is symmetric, if the user passes
   ``<block1, block2>``, then ``GetCovarianceBlock`` can be called with
   ``block1``, ``block2`` as well as ``block2``, ``block1``.

   ``covariance_blocks`` cannot contain duplicates. Bad things will
   happen if they do.

   Note that the list of ``covariance_blocks`` is only used to
   determine what parts of the covariance matrix are computed. The
   full Jacobian is used to do the computation, i.e. they do not have
   an impact on what part of the Jacobian is used for computation.

   The return value indicates the success or failure of the covariance
   computation. Please see the documentation for
   :class:`Covariance::Options` for more on the conditions under which
   this function returns ``false``.

.. function:: bool GetCovarianceBlock(const double* parameter_block1, const double* parameter_block2, double* covariance_block) const

   Return the block of the cross-covariance matrix corresponding to
   ``parameter_block1`` and ``parameter_block2``.

   Compute must be called before the first call to ``GetCovarianceBlock``
   and the pair ``<parameter_block1, parameter_block2>`` OR the pair
   ``<parameter_block2, parameter_block1>`` must have been present in the
   vector covariance_blocks when ``Compute`` was called. Otherwise
   ``GetCovarianceBlock`` will return false.

   ``covariance_block`` must point to a memory location that can store
   a ``parameter_block1_size x parameter_block2_size`` matrix. The
   returned covariance will be a row-major matrix.

.. function:: bool GetCovarianceBlockInTangentSpace(const double* parameter_block1, const double* parameter_block2, double* covariance_block) const

   Return the block of the cross-covariance matrix corresponding to
   ``parameter_block1`` and ``parameter_block2``.
   Returns cross-covariance in the tangent space if a local
   parameterization is associated with either parameter block;
   else returns cross-covariance in the ambient space.

   Compute must be called before the first call to ``GetCovarianceBlock``
   and the pair ``<parameter_block1, parameter_block2>`` OR the pair
   ``<parameter_block2, parameter_block1>`` must have been present in the
   vector covariance_blocks when ``Compute`` was called. Otherwise
   ``GetCovarianceBlock`` will return false.

   ``covariance_block`` must point to a memory location that can store
   a ``parameter_block1_local_size x parameter_block2_local_size`` matrix. The
   returned covariance will be a row-major matrix.

Example Usage
=============

.. code-block:: c++

 double x[3];
 double y[2];

 Problem problem;
 problem.AddParameterBlock(x, 3);
 problem.AddParameterBlock(y, 2);
 <Build Problem>
 <Solve Problem>

 Covariance::Options options;
 Covariance covariance(options);

 vector<pair<const double*, const double*> > covariance_blocks;
 covariance_blocks.push_back(make_pair(x, x));
 covariance_blocks.push_back(make_pair(y, y));
 covariance_blocks.push_back(make_pair(x, y));

 CHECK(covariance.Compute(covariance_blocks, &problem));

 double covariance_xx[3 * 3];
 double covariance_yy[2 * 2];
 double covariance_xy[3 * 2];
 covariance.GetCovarianceBlock(x, x, covariance_xx)
 covariance.GetCovarianceBlock(y, y, covariance_yy)
 covariance.GetCovarianceBlock(x, y, covariance_xy)
