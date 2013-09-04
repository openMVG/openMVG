.. default-domain:: cpp

.. cpp:namespace:: ceres

.. _`chapter-modeling`:

========
Modeling
========

Recall that Ceres solves robustified non-linear least squares problems
of the form

.. math:: \frac{1}{2}\sum_{i=1} \rho_i\left(\left\|f_i\left(x_{i_1}, ... ,x_{i_k}\right)\right\|^2\right).
   :label: ceresproblem

The expression
:math:`\rho_i\left(\left\|f_i\left(x_{i_1},...,x_{i_k}\right)\right\|^2\right)`
is known as a ``ResidualBlock``, where :math:`f_i(\cdot)` is a
:class:`CostFunction` that depends on the parameter blocks
:math:`\left[x_{i_1},... , x_{i_k}\right]`. In most optimization
problems small groups of scalars occur together. For example the three
components of a translation vector and the four components of the
quaternion that define the pose of a camera. We refer to such a group
of small scalars as a ``ParameterBlock``. Of course a
``ParameterBlock`` can just be a single parameter. :math:`\rho_i` is a
:class:`LossFunction`. A :class:`LossFunction` is a scalar function
that is used to reduce the influence of outliers on the solution of
non-linear least squares problems.

In this chapter we will describe the various classes that are part of
Ceres Solver's modeling API, and how they can be used to construct an
optimization problem. Once a problem has been constructed, various
methods for solving them will be discussed in
:ref:`chapter-solving`. It is by design that the modeling and the
solving APIs are orthogonal to each other. This enables
switching/tweaking of various solver parameters without having to
touch the problem once it has been successfully modeled.

:class:`CostFunction`
---------------------

The single biggest task when modeling a problem is specifying the
residuals and their derivatives. This is done using
:class:`CostFunction` objects.

.. class:: CostFunction

   .. code-block:: c++

    class CostFunction {
     public:
      virtual bool Evaluate(double const* const* parameters,
                            double* residuals,
                            double** jacobians) = 0;
      const vector<int16>& parameter_block_sizes();
      int num_residuals() const;

     protected:
      vector<int16>* mutable_parameter_block_sizes();
      void set_num_residuals(int num_residuals);
    };

   Given parameter blocks :math:`\left[x_{i_1}, ... , x_{i_k}\right]`,
   a :class:`CostFunction` is responsible for computing a vector of
   residuals and if asked a vector of Jacobian matrices, i.e., given
   :math:`\left[x_{i_1}, ... , x_{i_k}\right]`, compute the vector
   :math:`f_i\left(x_{i_1},...,x_{i_k}\right)` and the matrices

   .. math:: J_{ij} = \frac{\partial}{\partial x_{i_j}}f_i\left(x_{i_1},...,x_{i_k}\right),\quad \forall j \in \{i_1,..., i_k\}

   The signature of the :class:`CostFunction` (number and sizes of
   input parameter blocks and number of outputs) is stored in
   :member:`CostFunction::parameter_block_sizes_` and
   :member:`CostFunction::num_residuals_` respectively. User code
   inheriting from this class is expected to set these two members
   with the corresponding accessors. This information will be verified
   by the :class:`Problem` when added with
   :func:`Problem::AddResidualBlock`.

.. function:: bool CostFunction::Evaluate(double const* const* parameters, double* residuals, double** jacobians)

   Compute the residual vector and the Jacobian matrices.

   ``parameters`` is an array of pointers to arrays containing the
   various parameter blocks. ``parameters`` has the same number of
   elements as :member:`CostFunction::parameter_block_sizes_` and the
   parameter blocks are in the same order as
   :member:`CostFunction::parameter_block_sizes_`.

   ``residuals`` is an array of size ``num_residuals_``.

   ``jacobians`` is an array of size
   :member:`CostFunction::parameter_block_sizes_` containing pointers
   to storage for Jacobian matrices corresponding to each parameter
   block. The Jacobian matrices are in the same order as
   :member:`CostFunction::parameter_block_sizes_`. ``jacobians[i]`` is
   an array that contains :member:`CostFunction::num_residuals_` x
   :member:`CostFunction::parameter_block_sizes_` ``[i]``
   elements. Each Jacobian matrix is stored in row-major order, i.e.,
   ``jacobians[i][r * parameter_block_size_[i] + c]`` =
   :math:`\frac{\partial residual[r]}{\partial parameters[i][c]}`


   If ``jacobians`` is ``NULL``, then no derivatives are returned;
   this is the case when computing cost only. If ``jacobians[i]`` is
   ``NULL``, then the Jacobian matrix corresponding to the
   :math:`i^{\textrm{th}}` parameter block must not be returned, this
   is the case when a parameter block is marked constant.

   **NOTE** The return value indicates whether the computation of the
   residuals and/or jacobians was successful or not.

   This can be used to communicate numerical failures in Jacobian
   computations for instance.

   A more interesting and common use is to impose constraints on the
   parameters. If the initial values of the parameter blocks satisfy
   the constraints, then returning false whenever the constraints are
   not satisfied will prevent the solver from moving into the
   infeasible region. This is not a very sophisticated mechanism for
   enforcing constraints, but is often good enough for things like
   non-negativity constraints.

   Note that it is important that the initial values of the parameter
   block must be feasible, otherwise the solver will declare a
   numerical problem at iteration 0.


:class:`SizedCostFunction`
--------------------------

.. class:: SizedCostFunction

   If the size of the parameter blocks and the size of the residual
   vector is known at compile time (this is the common case),
   :class:`SizeCostFunction` can be used where these values can be
   specified as template parameters and the user only needs to
   implement :func:`CostFunction::Evaluate`.

   .. code-block:: c++

    template<int kNumResiduals,
             int N0 = 0, int N1 = 0, int N2 = 0, int N3 = 0, int N4 = 0,
             int N5 = 0, int N6 = 0, int N7 = 0, int N8 = 0, int N9 = 0>
    class SizedCostFunction : public CostFunction {
     public:
      virtual bool Evaluate(double const* const* parameters,
                            double* residuals,
                            double** jacobians) const = 0;
    };


:class:`AutoDiffCostFunction`
-----------------------------

.. class:: AutoDiffCostFunction

   Defining a :class:`CostFunction` or a :class:`SizedCostFunction`
   can be a tedious and error prone especially when computing
   derivatives.  To this end Ceres provides `automatic differentiation
   <http://en.wikipedia.org/wiki/Automatic_differentiation>`_.

   .. code-block:: c++

     template <typename CostFunctor,
            int M,        // Number of residuals, or ceres::DYNAMIC.
            int N0,       // Number of parameters in block 0.
            int N1 = 0,   // Number of parameters in block 1.
            int N2 = 0,   // Number of parameters in block 2.
            int N3 = 0,   // Number of parameters in block 3.
            int N4 = 0,   // Number of parameters in block 4.
            int N5 = 0,   // Number of parameters in block 5.
            int N6 = 0,   // Number of parameters in block 6.
            int N7 = 0,   // Number of parameters in block 7.
            int N8 = 0,   // Number of parameters in block 8.
            int N9 = 0>   // Number of parameters in block 9.
     class AutoDiffCostFunction : public
     SizedCostFunction<M, N0, N1, N2, N3, N4, N5, N6, N7, N8, N9> {
     };

   To get an auto differentiated cost function, you must define a
   class with a templated ``operator()`` (a functor) that computes the
   cost function in terms of the template parameter ``T``. The
   autodiff framework substitutes appropriate ``Jet`` objects for
   ``T`` in order to compute the derivative when necessary, but this
   is hidden, and you should write the function as if ``T`` were a
   scalar type (e.g. a double-precision floating point number).

   The function must write the computed value in the last argument
   (the only non-``const`` one) and return true to indicate success.
   Please see :class:`CostFunction` for details on how the return
   value may be used to impose simple constraints on the parameter
   block.

   For example, consider a scalar error :math:`e = k - x^\top y`,
   where both :math:`x` and :math:`y` are two-dimensional vector
   parameters and :math:`k` is a constant. The form of this error,
   which is the difference between a constant and an expression, is a
   common pattern in least squares problems. For example, the value
   :math:`x^\top y` might be the model expectation for a series of
   measurements, where there is an instance of the cost function for
   each measurement :math:`k`.

   The actual cost added to the total problem is :math:`e^2`, or
   :math:`(k - x^\top y)^2`; however, the squaring is implicitly done
   by the optimization framework.

   To write an auto-differentiable cost function for the above model,
   first define the object

   .. code-block:: c++

    class MyScalarCostFunctor {
      MyScalarCostFunctor(double k): k_(k) {}

      template <typename T>
      bool operator()(const T* const x , const T* const y, T* e) const {
        e[0] = T(k_) - x[0] * y[0] - x[1] * y[1];
        return true;
      }

     private:
      double k_;
    };


   Note that in the declaration of ``operator()`` the input parameters
   ``x`` and ``y`` come first, and are passed as const pointers to arrays
   of ``T``. If there were three input parameters, then the third input
   parameter would come after ``y``. The output is always the last
   parameter, and is also a pointer to an array. In the example above,
   ``e`` is a scalar, so only ``e[0]`` is set.

   Then given this class definition, the auto differentiated cost
   function for it can be constructed as follows.

   .. code-block:: c++

    CostFunction* cost_function
        = new AutoDiffCostFunction<MyScalarCostFunctor, 1, 2, 2>(
            new MyScalarCostFunctor(1.0));              ^  ^  ^
                                                        |  |  |
                            Dimension of residual ------+  |  |
                            Dimension of x ----------------+  |
                            Dimension of y -------------------+


   In this example, there is usually an instance for each measurement
   of ``k``.

   In the instantiation above, the template parameters following
   ``MyScalarCostFunction``, ``<1, 2, 2>`` describe the functor as
   computing a 1-dimensional output from two arguments, both
   2-dimensional.

   The framework can currently accommodate cost functions of up to 10
   independent variables, and there is no limit on the dimensionality
   of each of them.

   **WARNING 1** Since the functor will get instantiated with
   different types for ``T``, you must convert from other numeric
   types to ``T`` before mixing computations with other variables
   of type ``T``. In the example above, this is seen where instead of
   using ``k_`` directly, ``k_`` is wrapped with ``T(k_)``.

   **WARNING 2** A common beginner's error when first using
   :class:`AutoDiffCostFunction` is to get the sizing wrong. In particular,
   there is a tendency to set the template parameters to (dimension of
   residual, number of parameters) instead of passing a dimension
   parameter for *every parameter block*. In the example above, that
   would be ``<MyScalarCostFunction, 1, 2>``, which is missing the 2
   as the last template argument.


:class:`DynamicAutoDiffCostFunction`
------------------------------------

.. class:: DynamicAutoDiffCostFunction

   :class:`AutoDiffCostFunction` requires that the number of parameter
   blocks and their sizes be known at compile time, e.g., Bezier curve
   fitting, Neural Network training etc. It also has an upper limit of
   10 parameter blocks. In a number of applications, this is not
   enough.

     .. code-block:: c++

      template <typename CostFunctor, int Stride = 4>
      class DynamicAutoDiffCostFunction : public CostFunction {
      };

   In such cases :class:`DynamicAutoDiffCostFunction` can be
   used. Like :class:`AutoDiffCostFunction` the user must define a
   templated functor, but the signature of the functor differs
   slightly. The expected interface for the cost functors is:

     .. code-block:: c++

       struct MyCostFunctor {
         template<typename T>
         bool operator()(T const* const* parameters, T* residuals) const {
         }
       }

   Since the sizing of the parameters is done at runtime, you must
   also specify the sizes after creating the dynamic autodiff cost
   function. For example:

     .. code-block:: c++

       DynamicAutoDiffCostFunction<MyCostFunctor, 4> cost_function(
           new MyCostFunctor());
       cost_function.AddParameterBlock(5);
       cost_function.AddParameterBlock(10);
       cost_function.SetNumResiduals(21);

   Under the hood, the implementation evaluates the cost function
   multiple times, computing a small set of the derivatives (four by
   default, controlled by the ``Stride`` template parameter) with each
   pass. There is a performance tradeoff with the size of the passes;
   Smaller sizes are more cache efficient but result in larger number
   of passes, and larger stride lengths can destroy cache-locality
   while reducing the number of passes over the cost function. The
   optimal value depends on the number and sizes of the various
   parameter blocks.

   As a rule of thumb, try using :class:`AutoDiffCostFunction` before
   you use :class:`DynamicAutoDiffCostFunction`.

:class:`NumericDiffCostFunction`
--------------------------------

.. class:: NumericDiffCostFunction

  In some cases, its not possible to define a templated cost functor,
  for example when the evaluation of the residual involves a call to a
  library function that you do not have control over.  In such a
  situation, `numerical differentiation
  <http://en.wikipedia.org/wiki/Numerical_differentiation>`_ can be
  used.

    .. code-block:: c++

      template <typename CostFunctionNoJacobian,
                NumericDiffMethod method = CENTRAL, int M = 0,
                int N0 = 0, int N1 = 0, int N2 = 0, int N3 = 0, int N4 = 0,
                int N5 = 0, int N6 = 0, int N7 = 0, int N8 = 0, int N9 = 0>
      class NumericDiffCostFunction
        : public SizedCostFunction<M, N0, N1, N2, N3, N4, N5, N6, N7, N8, N9> {
      };

   To get a numerically differentiated :class:`CostFunction`, you must
   define a class with a ``operator()`` (a functor) that computes the
   residuals. The functor must write the computed value in the last
   argument (the only non-``const`` one) and return ``true`` to
   indicate success.  Please see :class:`CostFunction` for details on
   how the return value may be used to impose simple constraints on
   the parameter block. e.g., an object of the form

   .. code-block:: c++

     struct ScalarFunctor {
      public:
       bool operator()(const double* const x1,
                       const double* const x2,
                       double* residuals) const;
     }

   For example, consider a scalar error :math:`e = k - x'y`, where
   both :math:`x` and :math:`y` are two-dimensional column vector
   parameters, the prime sign indicates transposition, and :math:`k`
   is a constant. The form of this error, which is the difference
   between a constant and an expression, is a common pattern in least
   squares problems. For example, the value :math:`x'y` might be the
   model expectation for a series of measurements, where there is an
   instance of the cost function for each measurement :math:`k`.

   To write an numerically-differentiable class:`CostFunction` for the
   above model, first define the object

   .. code-block::  c++

     class MyScalarCostFunctor {
       MyScalarCostFunctor(double k): k_(k) {}

       bool operator()(const double* const x,
                       const double* const y,
                       double* residuals) const {
         residuals[0] = k_ - x[0] * y[0] + x[1] * y[1];
         return true;
       }

      private:
       double k_;
     };

   Note that in the declaration of ``operator()`` the input parameters
   ``x`` and ``y`` come first, and are passed as const pointers to
   arrays of ``double`` s. If there were three input parameters, then
   the third input parameter would come after ``y``. The output is
   always the last parameter, and is also a pointer to an array. In
   the example above, the residual is a scalar, so only
   ``residuals[0]`` is set.

   Then given this class definition, the numerically differentiated
   :class:`CostFunction` with central differences used for computing
   the derivative can be constructed as follows.

   .. code-block:: c++

     CostFunction* cost_function
         = new NumericDiffCostFunction<MyScalarCostFunctor, CENTRAL, 1, 2, 2>(
             new MyScalarCostFunctor(1.0));                    ^     ^  ^  ^
                                                               |     |  |  |
                                   Finite Differencing Scheme -+     |  |  |
                                   Dimension of residual ------------+  |  |
                                   Dimension of x ----------------------+  |
                                   Dimension of y -------------------------+

   In this example, there is usually an instance for each measurement
   of `k`.

   In the instantiation above, the template parameters following
   ``MyScalarCostFunctor``, ``1, 2, 2``, describe the functor as
   computing a 1-dimensional output from two arguments, both
   2-dimensional.

   The framework can currently accommodate cost functions of up to 10
   independent variables, and there is no limit on the dimensionality
   of each of them.

   The ``CENTRAL`` difference method is considerably more accurate at
   the cost of twice as many function evaluations than forward
   difference. Consider using central differences begin with, and only
   after that works, trying forward difference to improve performance.

   **WARNING** A common beginner's error when first using
   NumericDiffCostFunction is to get the sizing wrong. In particular,
   there is a tendency to set the template parameters to (dimension of
   residual, number of parameters) instead of passing a dimension
   parameter for *every parameter*. In the example above, that would
   be ``<MyScalarCostFunctor, 1, 2>``, which is missing the last ``2``
   argument. Please be careful when setting the size parameters.


   **Alternate Interface**

   For a variety of reason, including compatibility with legacy code,
   :class:`NumericDiffCostFunction` can also take
   :class:`CostFunction` objects as input. The following describes
   how.

   To get a numerically differentiated cost function, define a
   subclass of :class:`CostFunction` such that the
   :func:`CostFunction::Evaluate` function ignores the ``jacobians``
   parameter. The numeric differentiation wrapper will fill in the
   jacobian parameter if necessary by repeatedly calling the
   :func:`CostFunction::Evaluate` with small changes to the
   appropriate parameters, and computing the slope. For performance,
   the numeric differentiation wrapper class is templated on the
   concrete cost function, even though it could be implemented only in
   terms of the :class:`CostFunction` interface.

   The numerically differentiated version of a cost function for a
   cost function can be constructed as follows:

   .. code-block:: c++

     CostFunction* cost_function
         = new NumericDiffCostFunction<MyCostFunction, CENTRAL, 1, 4, 8>(
             new MyCostFunction(...), TAKE_OWNERSHIP);

   where ``MyCostFunction`` has 1 residual and 2 parameter blocks with
   sizes 4 and 8 respectively. Look at the tests for a more detailed
   example.

:class:`NumericDiffFunctor`
---------------------------

.. class:: NumericDiffFunctor

   Sometimes parts of a cost function can be differentiated
   automatically or analytically but others require numeric
   differentiation. :class:`NumericDiffFunctor` is a wrapper class
   that takes a variadic functor evaluating a function, numerically
   differentiates it and makes it available as a templated functor so
   that it can be easily used as part of Ceres' automatic
   differentiation framework.

   For example, let us assume that

   .. code-block:: c++

     struct IntrinsicProjection
       IntrinsicProjection(const double* observations);
       bool operator()(const double* calibration,
                       const double* point,
                       double* residuals);
     };

   is a functor that implements the projection of a point in its local
   coordinate system onto its image plane and subtracts it from the
   observed point projection.

   Now we would like to compose the action of this functor with the
   action of camera extrinsics, i.e., rotation and translation, which
   is given by the following templated function

   .. code-block:: c++

     template<typename T>
     void RotateAndTranslatePoint(const T* rotation,
                                  const T* translation,
                                  const T* point,
                                  T* result);

   To compose the extrinsics and intrinsics, we can construct a
   ``CameraProjection`` functor as follows.

   .. code-block:: c++

    struct CameraProjection {
       typedef NumericDiffFunctor<IntrinsicProjection, CENTRAL, 2, 5, 3>
          IntrinsicProjectionFunctor;

      CameraProjection(double* observation) {
        intrinsic_projection_.reset(
            new IntrinsicProjectionFunctor(observation)) {
      }

      template <typename T>
      bool operator()(const T* rotation,
                      const T* translation,
                      const T* intrinsics,
                      const T* point,
                      T* residuals) const {
        T transformed_point[3];
        RotateAndTranslatePoint(rotation, translation, point, transformed_point);
        return (*intrinsic_projection_)(intrinsics, transformed_point, residual);
      }

     private:
      scoped_ptr<IntrinsicProjectionFunctor> intrinsic_projection_;
    };

   Here, we made the choice of using ``CENTRAL`` differences to compute
   the jacobian of ``IntrinsicProjection``.

   Now, we are ready to construct an automatically differentiated cost
   function as

   .. code-block:: c++

    CostFunction* cost_function =
        new AutoDiffCostFunction<CameraProjection, 2, 3, 3, 5>(
           new CameraProjection(observations));

   ``cost_function`` now seamlessly integrates automatic
   differentiation of ``RotateAndTranslatePoint`` with a numerically
   differentiated version of ``IntrinsicProjection``.


:class:`CostFunctionToFunctor`
------------------------------

.. class:: CostFunctionToFunctor

   Just like :class:`NumericDiffFunctor` allows numeric
   differentiation to be mixed with automatic differentiation,
   :class:`CostFunctionToFunctor` provides an even more general
   mechanism.  :class:`CostFunctionToFunctor` is an adapter class that
   allows users to use :class:`CostFunction` objects in templated
   functors which are to be used for automatic differentiation.  This
   allows the user to seamlessly mix analytic, numeric and automatic
   differentiation.

   For example, let us assume that

   .. code-block:: c++

     class IntrinsicProjection : public SizedCostFunction<2, 5, 3> {
       public:
         IntrinsicProjection(const double* observations);
         virtual bool Evaluate(double const* const* parameters,
                               double* residuals,
                               double** jacobians) const;
     };

   is a :class:`CostFunction` that implements the projection of a
   point in its local coordinate system onto its image plane and
   subtracts it from the observed point projection. It can compute its
   residual and either via analytic or numerical differentiation can
   compute its jacobians.

   Now we would like to compose the action of this
   :class:`CostFunction` with the action of camera extrinsics, i.e.,
   rotation and translation. Say we have a templated function

   .. code-block:: c++

      template<typename T>
      void RotateAndTranslatePoint(const T* rotation,
                                   const T* translation,
                                   const T* point,
                                   T* result);


   Then we can now do the following,

   .. code-block:: c++

    struct CameraProjection {
      CameraProjection(double* observation) {
        intrinsic_projection_.reset(
            new CostFunctionToFunctor<2, 5, 3>(new IntrinsicProjection(observation_)));
      }
      template <typename T>
      bool operator()(const T* rotation,
                      const T* translation,
                      const T* intrinsics,
                      const T* point,
                      T* residual) const {
        T transformed_point[3];
        RotateAndTranslatePoint(rotation, translation, point, transformed_point);

        // Note that we call intrinsic_projection_, just like it was
        // any other templated functor.
        return (*intrinsic_projection_)(intrinsics, transformed_point, residual);
      }

     private:
      scoped_ptr<CostFunctionToFunctor<2,5,3> > intrinsic_projection_;
    };



:class:`ConditionedCostFunction`
--------------------------------

.. class:: ConditionedCostFunction

   This class allows you to apply different conditioning to the residual
   values of a wrapped cost function. An example where this is useful is
   where you have an existing cost function that produces N values, but you
   want the total cost to be something other than just the sum of these
   squared values - maybe you want to apply a different scaling to some
   values, to change their contribution to the cost.

   Usage:

   .. code-block:: c++

       //  my_cost_function produces N residuals
       CostFunction* my_cost_function = ...
       CHECK_EQ(N, my_cost_function->num_residuals());
       vector<CostFunction*> conditioners;

       //  Make N 1x1 cost functions (1 parameter, 1 residual)
       CostFunction* f_1 = ...
       conditioners.push_back(f_1);

       CostFunction* f_N = ...
       conditioners.push_back(f_N);
       ConditionedCostFunction* ccf =
         new ConditionedCostFunction(my_cost_function, conditioners);


   Now ``ccf`` 's ``residual[i]`` (i=0..N-1) will be passed though the
   :math:`i^{\text{th}}` conditioner.

   .. code-block:: c++

      ccf_residual[i] = f_i(my_cost_function_residual[i])

   and the Jacobian will be affected appropriately.


:class:`NormalPrior`
--------------------

.. class:: NormalPrior

   .. code-block:: c++

     class NormalPrior: public CostFunction {
      public:
       // Check that the number of rows in the vector b are the same as the
       // number of columns in the matrix A, crash otherwise.
       NormalPrior(const Matrix& A, const Vector& b);

       virtual bool Evaluate(double const* const* parameters,
                             double* residuals,
                             double** jacobians) const;
      };

   Implements a cost function of the form

   .. math::  cost(x) = ||A(x - b)||^2

   where, the matrix A and the vector b are fixed and x is the
   variable. In case the user is interested in implementing a cost
   function of the form

  .. math::  cost(x) = (x - \mu)^T S^{-1} (x - \mu)

  where, :math:`\mu` is a vector and :math:`S` is a covariance matrix,
  then, :math:`A = S^{-1/2}`, i.e the matrix :math:`A` is the square
  root of the inverse of the covariance, also known as the stiffness
  matrix. There are however no restrictions on the shape of
  :math:`A`. It is free to be rectangular, which would be the case if
  the covariance matrix :math:`S` is rank deficient.



:class:`LossFunction`
---------------------

.. class:: LossFunction

   For least squares problems where the minimization may encounter
   input terms that contain outliers, that is, completely bogus
   measurements, it is important to use a loss function that reduces
   their influence.

   Consider a structure from motion problem. The unknowns are 3D
   points and camera parameters, and the measurements are image
   coordinates describing the expected reprojected position for a
   point in a camera. For example, we want to model the geometry of a
   street scene with fire hydrants and cars, observed by a moving
   camera with unknown parameters, and the only 3D points we care
   about are the pointy tippy-tops of the fire hydrants. Our magic
   image processing algorithm, which is responsible for producing the
   measurements that are input to Ceres, has found and matched all
   such tippy-tops in all image frames, except that in one of the
   frame it mistook a car's headlight for a hydrant. If we didn't do
   anything special the residual for the erroneous measurement will
   result in the entire solution getting pulled away from the optimum
   to reduce the large error that would otherwise be attributed to the
   wrong measurement.

   Using a robust loss function, the cost for large residuals is
   reduced. In the example above, this leads to outlier terms getting
   down-weighted so they do not overly influence the final solution.

   .. code-block:: c++

    class LossFunction {
     public:
      virtual void Evaluate(double s, double out[3]) const = 0;
    };


   The key method is :func:`LossFunction::Evaluate`, which given a
   non-negative scalar ``s``, computes

   .. math:: out = \begin{bmatrix}\rho(s), & \rho'(s), & \rho''(s)\end{bmatrix}

   Here the convention is that the contribution of a term to the cost
   function is given by :math:`\frac{1}{2}\rho(s)`, where :math:`s
   =\|f_i\|^2`. Calling the method with a negative value of :math:`s`
   is an error and the implementations are not required to handle that
   case.

   Most sane choices of :math:`\rho` satisfy:

   .. math::

      \rho(0) &= 0\\
      \rho'(0) &= 1\\
      \rho'(s) &< 1 \text{ in the outlier region}\\
      \rho''(s) &< 0 \text{ in the outlier region}

   so that they mimic the squared cost for small residuals.

   **Scaling**

   Given one robustifier :math:`\rho(s)` one can change the length
   scale at which robustification takes place, by adding a scale
   factor :math:`a > 0` which gives us :math:`\rho(s,a) = a^2 \rho(s /
   a^2)` and the first and second derivatives as :math:`\rho'(s /
   a^2)` and :math:`(1 / a^2) \rho''(s / a^2)` respectively.


   The reason for the appearance of squaring is that :math:`a` is in
   the units of the residual vector norm whereas :math:`s` is a squared
   norm. For applications it is more convenient to specify :math:`a` than
   its square.

Instances
^^^^^^^^^

Ceres includes a number of predefined loss functions. For simplicity
we described their unscaled versions. The figure below illustrates
their shape graphically. More details can be found in
``include/ceres/loss_function.h``.

.. figure:: loss.png
   :figwidth: 500px
   :height: 400px
   :align: center

   Shape of the various common loss functions.

.. class:: TrivialLoss

      .. math:: \rho(s) = s

.. class:: HuberLoss

   .. math:: \rho(s) = \begin{cases} s & s \le 1\\ 2 \sqrt{s} - 1 & s > 1 \end{cases}

.. class:: SoftLOneLoss

   .. math:: \rho(s) = 2 (\sqrt{1+s} - 1)

.. class:: CauchyLoss

   .. math:: \rho(s) = \log(1 + s)

.. class:: ArctanLoss

   .. math:: \rho(s) = \arctan(s)

.. class:: TolerantLoss

   .. math:: \rho(s,a,b) = b \log(1 + e^{(s - a) / b}) - b \log(1 + e^{-a / b})

.. class:: ComposedLoss

   Given two loss functions ``f`` and ``g``, implements the loss
   function ``h(s) = f(g(s))``.

   .. code-block:: c++

      class ComposedLoss : public LossFunction {
       public:
        explicit ComposedLoss(const LossFunction* f,
                              Ownership ownership_f,
                              const LossFunction* g,
                              Ownership ownership_g);
      };

.. class:: ScaledLoss

   Sometimes you want to simply scale the output value of the
   robustifier. For example, you might want to weight different error
   terms differently (e.g., weight pixel reprojection errors
   differently from terrain errors).

   Given a loss function :math:`\rho(s)` and a scalar :math:`a`, :class:`ScaledLoss`
   implements the function :math:`a \rho(s)`.

   Since we treat the a ``NULL`` Loss function as the Identity loss
   function, :math:`rho` = ``NULL``: is a valid input and will result
   in the input being scaled by :math:`a`. This provides a simple way
   of implementing a scaled ResidualBlock.

.. class:: LossFunctionWrapper

   Sometimes after the optimization problem has been constructed, we
   wish to mutate the scale of the loss function. For example, when
   performing estimation from data which has substantial outliers,
   convergence can be improved by starting out with a large scale,
   optimizing the problem and then reducing the scale. This can have
   better convergence behavior than just using a loss function with a
   small scale.

   This templated class allows the user to implement a loss function
   whose scale can be mutated after an optimization problem has been
   constructed. e.g,

   .. code-block:: c++

     Problem problem;

     // Add parameter blocks

     CostFunction* cost_function =
         new AutoDiffCostFunction < UW_Camera_Mapper, 2, 9, 3>(
             new UW_Camera_Mapper(feature_x, feature_y));

     LossFunctionWrapper* loss_function(new HuberLoss(1.0), TAKE_OWNERSHIP);
     problem.AddResidualBlock(cost_function, loss_function, parameters);

     Solver::Options options;
     Solver::Summary summary;
     Solve(options, &problem, &summary);

     loss_function->Reset(new HuberLoss(1.0), TAKE_OWNERSHIP);
     Solve(options, &problem, &summary);


Theory
^^^^^^

Let us consider a problem with a single problem and a single parameter
block.

.. math::

 \min_x \frac{1}{2}\rho(f^2(x))


Then, the robustified gradient and the Gauss-Newton Hessian are

.. math::

        g(x) &= \rho'J^\top(x)f(x)\\
        H(x) &= J^\top(x)\left(\rho' + 2 \rho''f(x)f^\top(x)\right)J(x)

where the terms involving the second derivatives of :math:`f(x)` have
been ignored. Note that :math:`H(x)` is indefinite if
:math:`\rho''f(x)^\top f(x) + \frac{1}{2}\rho' < 0`. If this is not
the case, then its possible to re-weight the residual and the Jacobian
matrix such that the corresponding linear least squares problem for
the robustified Gauss-Newton step.


Let :math:`\alpha` be a root of

.. math:: \frac{1}{2}\alpha^2 - \alpha - \frac{\rho''}{\rho'}\|f(x)\|^2 = 0.


Then, define the rescaled residual and Jacobian as

.. math::

        \tilde{f}(x) &= \frac{\sqrt{\rho'}}{1 - \alpha} f(x)\\
        \tilde{J}(x) &= \sqrt{\rho'}\left(1 - \alpha
                        \frac{f(x)f^\top(x)}{\left\|f(x)\right\|^2} \right)J(x)


In the case :math:`2 \rho''\left\|f(x)\right\|^2 + \rho' \lesssim 0`,
we limit :math:`\alpha \le 1- \epsilon` for some small
:math:`\epsilon`. For more details see [Triggs]_.

With this simple rescaling, one can use any Jacobian based non-linear
least squares algorithm to robustified non-linear least squares
problems.


:class:`LocalParameterization`
------------------------------

.. class:: LocalParameterization

   .. code-block:: c++

     class LocalParameterization {
      public:
       virtual ~LocalParameterization() {}
       virtual bool Plus(const double* x,
                         const double* delta,
                         double* x_plus_delta) const = 0;
       virtual bool ComputeJacobian(const double* x, double* jacobian) const = 0;
       virtual int GlobalSize() const = 0;
       virtual int LocalSize() const = 0;
     };

   Sometimes the parameters :math:`x` can overparameterize a
   problem. In that case it is desirable to choose a parameterization
   to remove the null directions of the cost. More generally, if
   :math:`x` lies on a manifold of a smaller dimension than the
   ambient space that it is embedded in, then it is numerically and
   computationally more effective to optimize it using a
   parameterization that lives in the tangent space of that manifold
   at each point.

   For example, a sphere in three dimensions is a two dimensional
   manifold, embedded in a three dimensional space. At each point on
   the sphere, the plane tangent to it defines a two dimensional
   tangent space. For a cost function defined on this sphere, given a
   point :math:`x`, moving in the direction normal to the sphere at
   that point is not useful. Thus a better way to parameterize a point
   on a sphere is to optimize over two dimensional vector
   :math:`\Delta x` in the tangent space at the point on the sphere
   point and then "move" to the point :math:`x + \Delta x`, where the
   move operation involves projecting back onto the sphere. Doing so
   removes a redundant dimension from the optimization, making it
   numerically more robust and efficient.

   More generally we can define a function

   .. math:: x' = \boxplus(x, \Delta x),

   where :math:`x'` has the same size as :math:`x`, and :math:`\Delta
   x` is of size less than or equal to :math:`x`. The function
   :math:`\boxplus`, generalizes the definition of vector
   addition. Thus it satisfies the identity

   .. math:: \boxplus(x, 0) = x,\quad \forall x.

   Instances of :class:`LocalParameterization` implement the
   :math:`\boxplus` operation and its derivative with respect to
   :math:`\Delta x` at :math:`\Delta x = 0`.


.. function:: int LocalParameterization::GlobalSize()

   The dimension of the ambient space in which the parameter block
   :math:`x` lives.

.. function:: int LocalParamterization::LocaLocalSize()

   The size of the tangent space
   that :math:`\Delta x` lives in.

.. function:: bool LocalParameterization::Plus(const double* x, const double* delta, double* x_plus_delta) const

    :func:`LocalParameterization::Plus` implements :math:`\boxplus(x,\Delta x)`.

.. function:: bool LocalParameterization::ComputeJacobian(const double* x, double* jacobian) const

   Computes the Jacobian matrix

   .. math:: J = \left . \frac{\partial }{\partial \Delta x} \boxplus(x,\Delta x)\right|_{\Delta x = 0}

   in row major form.

Instances
^^^^^^^^^

.. class:: IdentityParameterization

   A trivial version of :math:`\boxplus` is when :math:`\Delta x` is
   of the same size as :math:`x` and

   .. math::  \boxplus(x, \Delta x) = x + \Delta x

.. class:: SubsetParameterization

   A more interesting case if :math:`x` is a two dimensional vector,
   and the user wishes to hold the first coordinate constant. Then,
   :math:`\Delta x` is a scalar and :math:`\boxplus` is defined as

   .. math::

      \boxplus(x, \Delta x) = x + \left[ \begin{array}{c} 0 \\ 1
                                  \end{array} \right] \Delta x

   :class:`SubsetParameterization` generalizes this construction to
   hold any part of a parameter block constant.

.. class:: QuaternionParameterization

   Another example that occurs commonly in Structure from Motion
   problems is when camera rotations are parameterized using a
   quaternion. There, it is useful only to make updates orthogonal to
   that 4-vector defining the quaternion. One way to do this is to let
   :math:`\Delta x` be a 3 dimensional vector and define
   :math:`\boxplus` to be

    .. math:: \boxplus(x, \Delta x) = \left[ \cos(|\Delta x|), \frac{\sin\left(|\Delta x|\right)}{|\Delta x|} \Delta x \right] * x
      :label: quaternion

   The multiplication between the two 4-vectors on the right hand side
   is the standard quaternion
   product. :class:`QuaternionParameterization` is an implementation
   of :eq:`quaternion`.



:class:`AutoDiffLocalParameterization`
--------------------------------------

.. class:: AutoDiffLocalParameterization

  :class:`AutoDiffLocalParameterization` does for
  :class:`LocalParameterization` what :class:`AutoDiffCostFunction`
  does for :class:`CostFunction`. It allows the user to define a
  templated functor that implements the
  :func:`LocalParameterization::Plus` operation and it uses automatic
  differentiation to implement the computation of the Jacobian.

  To get an auto differentiated local parameterization, you must
  define a class with a templated operator() (a functor) that computes

     .. math:: x' = \boxplus(x, \Delta x),

  For example, Quaternions have a three dimensional local
  parameterization. It's plus operation can be implemented as (taken
  from `internal/ceres/auto_diff_local_parameterization_test.cc
  <https://ceres-solver.googlesource.com/ceres-solver/+/master/include/ceres/local_parameterization.h>`_
  )

    .. code-block:: c++

      struct QuaternionPlus {
        template<typename T>
        bool operator()(const T* x, const T* delta, T* x_plus_delta) const {
          const T squared_norm_delta =
              delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];

          T q_delta[4];
          if (squared_norm_delta > T(0.0)) {
            T norm_delta = sqrt(squared_norm_delta);
            const T sin_delta_by_delta = sin(norm_delta) / norm_delta;
            q_delta[0] = cos(norm_delta);
            q_delta[1] = sin_delta_by_delta * delta[0];
            q_delta[2] = sin_delta_by_delta * delta[1];
            q_delta[3] = sin_delta_by_delta * delta[2];
          } else {
            // We do not just use q_delta = [1,0,0,0] here because that is a
            // constant and when used for automatic differentiation will
            // lead to a zero derivative. Instead we take a first order
            // approximation and evaluate it at zero.
            q_delta[0] = T(1.0);
            q_delta[1] = delta[0];
            q_delta[2] = delta[1];
            q_delta[3] = delta[2];
          }

          Quaternionproduct(q_delta, x, x_plus_delta);
          return true;
        }
      };

  Then given this struct, the auto differentiated local
  parameterization can now be constructed as

  .. code-block:: c++

     LocalParameterization* local_parameterization =
         new AutoDiffLocalParameterization<QuaternionPlus, 4, 3>;
                                                           |  |
                                Global Size ---------------+  |
                                Local Size -------------------+

  **WARNING:** Since the functor will get instantiated with different
  types for ``T``, you must to convert from other numeric types to
  ``T`` before mixing computations with other variables of type
  ``T``. In the example above, this is seen where instead of using
  ``k_`` directly, ``k_`` is wrapped with ``T(k_)``.


:class:`Problem`
----------------

.. class:: Problem

   :class:`Problem` holds the robustified non-linear least squares
   problem :eq:`ceresproblem`. To create a least squares problem, use
   the :func:`Problem::AddResidualBlock` and
   :func:`Problem::AddParameterBlock` methods.

   For example a problem containing 3 parameter blocks of sizes 3, 4
   and 5 respectively and two residual blocks of size 2 and 6:

   .. code-block:: c++

     double x1[] = { 1.0, 2.0, 3.0 };
     double x2[] = { 1.0, 2.0, 3.0, 5.0 };
     double x3[] = { 1.0, 2.0, 3.0, 6.0, 7.0 };

     Problem problem;
     problem.AddResidualBlock(new MyUnaryCostFunction(...), x1);
     problem.AddResidualBlock(new MyBinaryCostFunction(...), x2, x3);

   :func:`Problem::AddResidualBlock` as the name implies, adds a
   residual block to the problem. It adds a :class:`CostFunction`, an
   optional :class:`LossFunction` and connects the
   :class:`CostFunction` to a set of parameter block.

   The cost function carries with it information about the sizes of
   the parameter blocks it expects. The function checks that these
   match the sizes of the parameter blocks listed in
   ``parameter_blocks``. The program aborts if a mismatch is
   detected. ``loss_function`` can be ``NULL``, in which case the cost
   of the term is just the squared norm of the residuals.

   The user has the option of explicitly adding the parameter blocks
   using :func:`Problem::AddParameterBlock`. This causes additional
   correctness checking; however, :func:`Problem::AddResidualBlock`
   implicitly adds the parameter blocks if they are not present, so
   calling :func:`Problem::AddParameterBlock` explicitly is not
   required.

   :func:`Problem::AddParameterBlock` explicitly adds a parameter
   block to the :class:`Problem`. Optionally it allows the user to
   associate a :class:`LocalParameterization` object with the
   parameter block too. Repeated calls with the same arguments are
   ignored. Repeated calls with the same double pointer but a
   different size results in undefined behavior.

   You can set any parameter block to be constant using
   :func:`Problem::SetParameterBlockConstant` and undo this using
   :func:`SetParameterBlockVariable`.

   In fact you can set any number of parameter blocks to be constant,
   and Ceres is smart enough to figure out what part of the problem
   you have constructed depends on the parameter blocks that are free
   to change and only spends time solving it. So for example if you
   constructed a problem with a million parameter blocks and 2 million
   residual blocks, but then set all but one parameter blocks to be
   constant and say only 10 residual blocks depend on this one
   non-constant parameter block. Then the computational effort Ceres
   spends in solving this problem will be the same if you had defined
   a problem with one parameter block and 10 residual blocks.

   **Ownership**

   :class:`Problem` by default takes ownership of the
   ``cost_function``, ``loss_function`` and ``local_parameterization``
   pointers. These objects remain live for the life of the
   :class:`Problem`. If the user wishes to keep control over the
   destruction of these objects, then they can do this by setting the
   corresponding enums in the :class:`Problem::Options` struct.

   Note that even though the Problem takes ownership of ``cost_function``
   and ``loss_function``, it does not preclude the user from re-using
   them in another residual block. The destructor takes care to call
   delete on each ``cost_function`` or ``loss_function`` pointer only
   once, regardless of how many residual blocks refer to them.

.. function:: ResidualBlockId Problem::AddResidualBlock(CostFunction* cost_function, LossFunction* loss_function, const vector<double*> parameter_blocks)

   Add a residual block to the overall cost function. The cost
   function carries with it information about the sizes of the
   parameter blocks it expects. The function checks that these match
   the sizes of the parameter blocks listed in parameter_blocks. The
   program aborts if a mismatch is detected. loss_function can be
   NULL, in which case the cost of the term is just the squared norm
   of the residuals.

   The user has the option of explicitly adding the parameter blocks
   using AddParameterBlock. This causes additional correctness
   checking; however, AddResidualBlock implicitly adds the parameter
   blocks if they are not present, so calling AddParameterBlock
   explicitly is not required.

   The Problem object by default takes ownership of the
   cost_function and loss_function pointers. These objects remain
   live for the life of the Problem object. If the user wishes to
   keep control over the destruction of these objects, then they can
   do this by setting the corresponding enums in the Options struct.

   Note: Even though the Problem takes ownership of cost_function
   and loss_function, it does not preclude the user from re-using
   them in another residual block. The destructor takes care to call
   delete on each cost_function or loss_function pointer only once,
   regardless of how many residual blocks refer to them.

   Example usage:

   .. code-block:: c++

      double x1[] = {1.0, 2.0, 3.0};
      double x2[] = {1.0, 2.0, 5.0, 6.0};
      double x3[] = {3.0, 6.0, 2.0, 5.0, 1.0};

      Problem problem;

      problem.AddResidualBlock(new MyUnaryCostFunction(...), NULL, x1);
      problem.AddResidualBlock(new MyBinaryCostFunction(...), NULL, x2, x1);


.. function:: void Problem::AddParameterBlock(double* values, int size, LocalParameterization* local_parameterization)

   Add a parameter block with appropriate size to the problem.
   Repeated calls with the same arguments are ignored. Repeated calls
   with the same double pointer but a different size results in
   undefined behavior.

.. function:: void Problem::AddParameterBlock(double* values, int size)

   Add a parameter block with appropriate size and parameterization to
   the problem. Repeated calls with the same arguments are
   ignored. Repeated calls with the same double pointer but a
   different size results in undefined behavior.

.. function:: void Problem::RemoveResidualBlock(ResidualBlockId residual_block)

   Remove a residual block from the problem. Any parameters that the residual
   block depends on are not removed. The cost and loss functions for the
   residual block will not get deleted immediately; won't happen until the
   problem itself is deleted.

   **WARNING:** Removing a residual or parameter block will destroy
   the implicit ordering, rendering the jacobian or residuals returned
   from the solver uninterpretable. If you depend on the evaluated
   jacobian, do not use remove! This may change in a future release.
   Hold the indicated parameter block constant during optimization.

.. function:: void Problem::RemoveParameterBlock(double* values)

   Remove a parameter block from the problem. The parameterization of
   the parameter block, if it exists, will persist until the deletion
   of the problem (similar to cost/loss functions in residual block
   removal). Any residual blocks that depend on the parameter are also
   removed, as described above in RemoveResidualBlock().  If
   Problem::Options::enable_fast_parameter_block_removal is true, then
   the removal is fast (almost constant time). Otherwise, removing a
   parameter block will incur a scan of the entire Problem object.

   **WARNING:** Removing a residual or parameter block will destroy
   the implicit ordering, rendering the jacobian or residuals returned
   from the solver uninterpretable. If you depend on the evaluated
   jacobian, do not use remove! This may change in a future release.

.. function:: void Problem::SetParameterBlockConstant(double* values)

   Hold the indicated parameter block constant during optimization.

.. function:: void Problem::SetParameterBlockVariable(double* values)

   Allow the indicated parameter to vary during optimization.

.. function:: void Problem::SetParameterization(double* values, LocalParameterization* local_parameterization)

   Set the local parameterization for one of the parameter blocks.
   The local_parameterization is owned by the Problem by default. It
   is acceptable to set the same parameterization for multiple
   parameters; the destructor is careful to delete local
   parameterizations only once. The local parameterization can only be
   set once per parameter, and cannot be changed once set.

.. function:: int Problem::NumParameterBlocks() const

   Number of parameter blocks in the problem. Always equals
   parameter_blocks().size() and parameter_block_sizes().size().

.. function:: int Problem::NumParameters() const

   The size of the parameter vector obtained by summing over the sizes
   of all the parameter blocks.

.. function:: int Problem::NumResidualBlocks() const

   Number of residual blocks in the problem. Always equals
   residual_blocks().size().

.. function:: int Problem::NumResiduals() const

   The size of the residual vector obtained by summing over the sizes
   of all of the residual blocks.

.. function int Problem::ParameterBlockSize(const double* values) const;

   The size of the parameter block.

.. function int Problem::ParameterBlockLocalSize(const double* values) const;

  The size of local parameterization for the parameter block. If
  there is no local parameterization associated with this parameter
  block, then ``ParameterBlockLocalSize`` = ``ParameterBlockSize``.


.. function void Problem::GetParameterBlocks(vector<double*>* parameter_blocks) const;

  Fills the passed ``parameter_blocks`` vector with pointers to the
  parameter blocks currently in the problem. After this call,
  ``parameter_block.size() == NumParameterBlocks``.

.. function:: bool Problem::Evaluate(const Problem::EvaluateOptions& options, double* cost, vector<double>* residuals, vector<double>* gradient, CRSMatrix* jacobian)

   Evaluate a :class:`Problem`. Any of the output pointers can be
   `NULL`. Which residual blocks and parameter blocks are used is
   controlled by the :class:`Problem::EvaluateOptions` struct below.

   .. code-block:: c++

     Problem problem;
     double x = 1;
     problem.Add(new MyCostFunction, NULL, &x);

     double cost = 0.0;
     problem.Evaluate(Problem::EvaluateOptions(), &cost, NULL, NULL, NULL);

   The cost is evaluated at `x = 1`. If you wish to evaluate the
   problem at `x = 2`, then

   .. code-block:: c++

      x = 2;
      problem.Evaluate(Problem::EvaluateOptions(), &cost, NULL, NULL, NULL);

   is the way to do so.

   **NOTE** If no local parameterizations are used, then the size of
   the gradient vector is the sum of the sizes of all the parameter
   blocks. If a parameter block has a local parameterization, then
   it contributes "LocalSize" entries to the gradient vector.

.. class:: Problem::EvaluateOptions

   Options struct that is used to control :func:`Problem::Evaluate`.

.. member:: vector<double*> Problem::EvaluateOptions::parameter_blocks

   The set of parameter blocks for which evaluation should be
   performed. This vector determines the order in which parameter
   blocks occur in the gradient vector and in the columns of the
   jacobian matrix. If parameter_blocks is empty, then it is assumed
   to be equal to a vector containing ALL the parameter
   blocks. Generally speaking the ordering of the parameter blocks in
   this case depends on the order in which they were added to the
   problem and whether or not the user removed any parameter blocks.

   **NOTE** This vector should contain the same pointers as the ones
   used to add parameter blocks to the Problem. These parameter block
   should NOT point to new memory locations. Bad things will happen if
   you do.

.. member:: vector<ResidualBlockId> Problem::EvaluateOptions::residual_blocks

   The set of residual blocks for which evaluation should be
   performed. This vector determines the order in which the residuals
   occur, and how the rows of the jacobian are ordered. If
   residual_blocks is empty, then it is assumed to be equal to the
   vector containing all the parameter blocks.

``rotation.h``
--------------

Many applications of Ceres Solver involve optimization problems where
some of the variables correspond to rotations. To ease the pain of
work with the various representations of rotations (angle-axis,
quaternion and matrix) we provide a handy set of templated
functions. These functions are templated so that the user can use them
within Ceres Solver's automatic differentiation framework.

.. function:: void AngleAxisToQuaternion<T>(T const* angle_axis, T* quaternion)

   Convert a value in combined axis-angle representation to a
   quaternion.

   The value ``angle_axis`` is a triple whose norm is an angle in radians,
   and whose direction is aligned with the axis of rotation, and
   ``quaternion`` is a 4-tuple that will contain the resulting quaternion.

.. function:: void QuaternionToAngleAxis<T>(T const* quaternion, T* angle_axis)

   Convert a quaternion to the equivalent combined axis-angle
   representation.

   The value ``quaternion`` must be a unit quaternion - it is not
   normalized first, and ``angle_axis`` will be filled with a value
   whose norm is the angle of rotation in radians, and whose direction
   is the axis of rotation.

.. function:: void RotationMatrixToAngleAxis<T, row_stride, col_stride>(const MatrixAdapter<const T, row_stride, col_stride>& R, T * angle_axis)
.. function:: void AngleAxisToRotationMatrix<T, row_stride, col_stride>(T const * angle_axis, const MatrixAdapter<T, row_stride, col_stride>& R)
.. function:: void RotationMatrixToAngleAxis<T>(T const * R, T * angle_axis)
.. function:: void AngleAxisToRotationMatrix<T>(T const * angle_axis, T * R)

   Conversions between 3x3 rotation matrix with given column and row strides and
   axis-angle rotation representations. The functions that take a pointer to T instead
   of a MatrixAdapter assume a column major representation with unit row stride and a column stride of 3.

.. function:: void EulerAnglesToRotationMatrix<T, row_stride, col_stride>(const T* euler, const MatrixAdapter<T, row_stride, col_stride>& R)
.. function:: void EulerAnglesToRotationMatrix<T>(const T* euler, int row_stride, T* R)

   Conversions between 3x3 rotation matrix with given column and row strides and
   Euler angle (in degrees) rotation representations.

   The {pitch,roll,yaw} Euler angles are rotations around the {x,y,z}
   axes, respectively.  They are applied in that same order, so the
   total rotation R is Rz * Ry * Rx.

   The function that takes a pointer to T as the rotation matrix assumes a row
   major representation with unit column stride and a row stride of 3.
   The additional parameter row_stride is required to be 3.

.. function:: void QuaternionToScaledRotation<T, row_stride, col_stride>(const T q[4], const MatrixAdapter<T, row_stride, col_stride>& R)
.. function:: void QuaternionToScaledRotation<T>(const T q[4], T R[3 * 3])

   Convert a 4-vector to a 3x3 scaled rotation matrix.

   The choice of rotation is such that the quaternion
   :math:`\begin{bmatrix} 1 &0 &0 &0\end{bmatrix}` goes to an identity
   matrix and for small :math:`a, b, c` the quaternion
   :math:`\begin{bmatrix}1 &a &b &c\end{bmatrix}` goes to the matrix

   .. math::

     I + 2 \begin{bmatrix} 0 & -c & b \\ c & 0 & -a\\ -b & a & 0
           \end{bmatrix} + O(q^2)

   which corresponds to a Rodrigues approximation, the last matrix
   being the cross-product matrix of :math:`\begin{bmatrix} a& b&
   c\end{bmatrix}`. Together with the property that :math:`R(q1 * q2)
   = R(q1) * R(q2)` this uniquely defines the mapping from :math:`q` to
   :math:`R`.

   In the function that accepts a pointer to T instead of a MatrixAdapter,
   the rotation matrix ``R`` is a row-major matrix with unit column stride
   and a row stride of 3.

   No normalization of the quaternion is performed, i.e.
   :math:`R = \|q\|^2  Q`, where :math:`Q` is an orthonormal matrix
   such that :math:`\det(Q) = 1` and :math:`Q*Q' = I`.


.. function:: void QuaternionToRotation<T>(const T q[4], const MatrixAdapter<T, row_stride, col_stride>& R)
.. function:: void QuaternionToRotation<T>(const T q[4], T R[3 * 3])

   Same as above except that the rotation matrix is normalized by the
   Frobenius norm, so that :math:`R R' = I` (and :math:`\det(R) = 1`).

.. function:: void UnitQuaternionRotatePoint<T>(const T q[4], const T pt[3], T result[3])

   Rotates a point pt by a quaternion q:

   .. math:: \text{result} = R(q)  \text{pt}

   Assumes the quaternion is unit norm. If you pass in a quaternion
   with :math:`|q|^2 = 2` then you WILL NOT get back 2 times the
   result you get for a unit quaternion.


.. function:: void QuaternionRotatePoint<T>(const T q[4], const T pt[3], T result[3])

   With this function you do not need to assume that q has unit norm.
   It does assume that the norm is non-zero.

.. function:: void QuaternionProduct<T>(const T z[4], const T w[4], T zw[4])

   .. math:: zw = z * w

   where :math:`*` is the Quaternion product between 4-vectors.


.. function:: void CrossProduct<T>(const T x[3], const T y[3], T x_cross_y[3])

   .. math:: \text{x_cross_y} = x \times y

.. function:: void AngleAxisRotatePoint<T>(const T angle_axis[3], const T pt[3], T result[3])

   .. math:: y = R(\text{angle_axis}) x
