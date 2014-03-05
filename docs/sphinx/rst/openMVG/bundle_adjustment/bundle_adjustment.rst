**************************
bundle_adjustment
**************************

Bundle Adjustment (ajustement de faisceaux), is  a non linear optimization problem.
It looks to minimizing the residual error of the observed tracks (the reprojection errors of the structure :math:`X_j` to the images measures :math:`x_j^i`.
According:

* :math:`X_j` the Jnth 3D point of the structure of the scene,
* :math:`x_j^i` the observation of the projection of the 3D point :math:`X_j` in the image :math:`i`,
* :math:`P_i` the projection matrix of the image :math:`i`

For initial guess of the vector of parameters: :math:`\{X_j,P_i\}_{i,j}`: camera parameters :math:`\{P_i\}_i` and the scene structure :math:`\{X_j\}_j`.
An iterative algorithm *Levenberg-Marquardt* updates the parameters vector according a gradient descent to minimizes the residual reprojection cost:

.. math::
  \underset{ \{P_i\}_i, \{X_j\}_j}{minimize} \left\| \sum_{j=0}^{m} \sum_{i=0}^{n} x_j^i - P_i X_j \right\|_2

When subtle changes are observed on the cost function or on the norm of the parameters vector the algorithm is stopped.

* Pros :

  * convergence is observed when the initial vector is close to the optimal solution.
  
* Cons :

  * solution could be a local solution, if the cost function is not convex and the provided solution is not in the optimal convex set.

openMVG bundle_adjustment framework
=====================================

OpenMVG relies on the [Ceres-solver]_ Google library to perform the Bundle Adjustment.
In order to ease its usage openMVG provides:

* data container to setup the problem,
* functor (metric for various camera models),
* samples to show how to use them.

bundle_adjustment container
______________________________

Two containers are defined in order to refine :

* Cameras that have nothing in common:

  * A container to refine structure and cameras that do not share any common data.

    * The camera vector parameter is a 7 length vector: [Rotation(angle,axis), t, focal]
    * ``template<unsigned char NCamParam = 7> class BA_Problem_data``

* Cameras that share common data (intrinsics parameters):

  * A container to refine structure and cameras that share common data

    * To be used in the case of intrinsic parameters are shared between some cameras.
    * The camera euclidean motion is ``NExternalParam`` values long: [R(angle axis)|t]
    * The number intrinsic parameters ``NIntrinsicParam`` (i.e. 1): [focal]
    * ``template<unsigned char NExternalParam = 6,unsigned char NIntrinsicParam = 1> class BA_Problem_data_camMotionAndIntrinsic``

.. [Ceres-solver] Sameer Agarwal and Keir Mierle. Ceres Solver: Tutorial & Reference. Google Inc.
