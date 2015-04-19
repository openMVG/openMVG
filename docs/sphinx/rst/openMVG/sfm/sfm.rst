*************************
sfm
*************************

sfm is the module related to Structure from Motion.
It handles storage of SfM related data and method to solve SfM problems (camera pose estimation, structure triangulation, bundle_adjustment).

A generic SfM data container
=============================

:class:`SfM_Data` class contains all the data used to describe input of a SfM problem:

* a collection of **View**

  * the used images

* a collection of **camera extrinsic**

  * the camera poses

* a collection of **camera intrinsic**

  * the camera internal projection parameters

* a **structure**

  * the collection of landmark (3D points associated with 2d view observations)

.. code-block:: c++

  struct SfM_Data
  {
    /// Considered views
    Views views;

    /// Considered poses (indexed by view.id_pose)
    Poses poses;

    /// Considered camera intrinsics (indexed by view.id_cam)
    Intrinsics intrinsics;

    /// Structure (3D points with their 2D observations)
    Landmarks structure;

    // ...
  }

View concept
--------------

The view store information related to an image:

* image filename
* id_view (must be unique)
* id_pose
* id_intrinsic
* image size

Note that thanks to the usage of ids we can defined shared poses & shared intrinsics.

View type is **abstract** and provide a way to add new custom View type: i.e. GeoLocatedView (add GPS position, ...)

Camera Poses concept
---------------------

The camera pose store a 3D pose that define a camera rotation and position (camera rotation and center).

Camera Intrinsic concept
--------------------------

Define the parameter of a camera. It can be shared or not.
Intrinsics parameter are **abstract** and provide a way to easily add new custom camera type.

Structure/Landmarks concept
----------------------------

It defines the structure:

* 3D point with 2D view features observations.

SfM_Data cleaning
==================

Generic interface are defined to remove outlier observations:

* use a given residual pixel error to discard outlier,
* use a minimal angle along the track bearing vectors.

Triangulation
==================

Once the SfM_Data is filled with some landmark observations and poses we can compute their 3D location.

Two method are proposed:

* A blind method:

  * Triangulate tracks using all observations,

  * Inlier/Outlier classification is done with a cheirality test,

* A robust method:

  * Triangulate tracks using a RANSAC scheme,

  * Check cheirality and a pixel residual error.

Non linear refinement, Bundle Adjustment
==========================================

OpenMVG provides a generic bundle_adjustment framework to refine or keep as constant the following parameters:

* internal parameters,
* external parameters,
* 3D structure.

.. code-block:: c++

  SfM_Data sfm_data;
  // initialize the data
  // ...

  const double dResidual_before = RMSE(sfm_data);

  // Bundle adjustement over all the parameters:
  std::shared_ptr<Bundle_Adjustment> ba_object = std::make_shared<Bundle_Adjustment_Ceres>();
  ba_object->Adjust(sfm_data);

  const double dResidual_after = RMSE(sfm_data);

SfM Pipelines
==============

OpenMVG provides ready to use and customizable pipelines for:

* solving sequential SfM,
* solving global SfM,
* computing a Structure from known camera poses.

Sequential SfM
-------------------------

Global SfM
-------------------------

Structure computation from known camera poses
----------------------------------------------
