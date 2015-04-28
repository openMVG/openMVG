*******************
geometry
*******************

Pose
===============

:class:`Pose3` defines the 3D Pose as a 3D transform: [R|C] t = -RC

.. code-block:: c++

  // Define two poses and combine them
  Pose3 pose1(RotationAroundX(0.02), Vec3(0,0,-2));
  Pose3 pose2(RotationAroundX(0.06), Vec3(-1,0,-2));
  Pose3 combinedPose = pose1 * pose2;

  // Apply the pose to a 3DPoint (World to local coordinates):
  const Vec3 pt = combinedPose(Vec3(2.6453,3.32,6.3));

Frustum & Frustum intersection
==============================

Define a camera Frustum from Pose3 and intrinsic parameters as:

- an infinite Frustum (4 Half Spaces) (a pyramid);
- a truncated Frustum (6 Half Spaces) (a truncated pyramid).

This structure is used for testing frustum intersection (see if two camera can share some visual content).

.. code-block:: c++
  
  // Build two truncated frustum and test their intersection
  Frustum frustum1(somedata, minDepth, maxDepth);
  Frustum frustum2(somedata, minDepth, maxDepth);
  bool bIntersect = frustum1.intersect(frustum2);

  // Build two infinite frustum and test their intersection
  Frustum frustum1(somedata);
  Frustum frustum2(somedata);
  bool bIntersect = frustum1.intersect(frustum2);

7DoF registration between point set
====================================

Find the rigid registration between point set using a scale, rotation and translation model.

.. code-block:: c++ 

  // Simulate two point set, apply a known transformation and estimate it back:
  const int nbPoints = 10;
  const Mat x1 = Mat::Random(3,nbPoints);
  Mat x2 = x1;

  const double scale = 2.0;
  const Mat3 rot = (Eigen::AngleAxis<double>(.2, Vec3::UnitX())
    * Eigen::AngleAxis<double>(.3, Vec3::UnitY())
    * Eigen::AngleAxis<double>(.6, Vec3::UnitZ())).toRotationMatrix();
  const Vec3 t(0.5,-0.3,.38);

  for (int i=0; i < nbPoints; ++i)
  {
    const Vec3 pt = x1.col(i);
    x2.col(i) = (scale * rot * pt + t);
  }

  // Compute the Similarity transform
  double Sc;
  Mat3 Rc;
  Vec3 tc;
  FindRTS(x1, x2, &Sc, &tc, &Rc);
  // Optional non linear refinement of the found parameters
  Refine_RTS(x1,x2,&Sc,&tc,&Rc);

  std::cout << "\n"
    << "Scale " << Sc << "\n"
    << "Rot \n" << Rc << "\n"
    << "t " << tc.transpose();

  std::cout << "\nGT\n"
    << "Scale " << scale << "\n"
    << "Rot \n" << rot << "\n"
    << "t " << t.transpose();

  

