.. _chapter-users:

=====
Users
=====

* At `Google <http://www.google.com>`_, Ceres is used to:

  * Estimate the pose of `Street View`_ cars, aircrafts, and satellites.
  * Build 3D models for `PhotoTours`_.
  * Estimate satellite image sensor characteristics.
  * Stitch `panoramas`_ on Android and iOS.
  * Apply `Lens Blur`_ on Android.
  * Solve `bundle adjustment`_ and `SLAM`_ problems in `Project
    Tango`_.

* `Willow Garage`_ uses Ceres to solve `SLAM`_ problems.
* `Southwest Research Institute <http://www.swri.org/>`_ uses Ceres for
  `calibrating robot-camera systems`_.
* `Blender <http://www.blender.org>`_ uses Ceres for `planar
  tracking`_ and `bundle adjustment`_.
* `OpenMVG <http://imagine.enpc.fr/~moulonp/openMVG/>`_ an open source
  multi-view geometry library uses Ceres for `bundle adjustment`_.
* `Microsoft Research <http://research.microsoft.com/en-us/>`_ uses
  Ceres for nonlinear optimization of objectives involving subdivision
  surfaces under `skinned control meshes`_.
* `Matterport <http://www.matterport.com>`_, uses Ceres for global
  alignment of 3D point clouds and for pose graph optimization.
* `Obvious Engineering <http://obviousengine.com/>`_ uses Ceres for
  bundle adjustment for their 3D photography app `Seene
  <http://seene.co/>`_.
* The `Autonomous Systems Lab <http://www.asl.ethz.ch/>`_ at ETH
  Zurich uses Ceres for

  * Camera and Camera/IMU Calibration.
  * Large scale optimization of visual, inertial, gps and
    wheel-odometry data for long term autonomy.

* `OpenPTrack <http://openptrack.org/>`_ uses Ceres for camera
  calibration.
* The `Intelligent Autonomous System Lab <http://robotics.dei.unipd.it/>`_
  at University of Padova, Italy, uses Ceres for

  * Camera/depth sensors network calibration.
  * Depth sensor distortion map estimation.

* `Theia <http://cs.ucsb.edu/~cmsweeney/theia>`_ is an open source
  Structure from Motion library that uses Ceres for `bundle adjustment`_
  and camera pose estimation.

* The `Applied Research Laboratory <https://www.arl.psu.edu/>`_ at
  Pennsylvania State University uses in their synthetic aperture Sonar
  beamforming engine, called ASASIN , for estimating platform
  kinematics.

.. _bundle adjustment: http://en.wikipedia.org/wiki/Structure_from_motion
.. _Street View: http://youtu.be/z00ORu4bU-A
.. _PhotoTours: http://google-latlong.blogspot.com/2012/04/visit-global-landmarks-with-photo-tours.html
.. _panoramas: http://www.google.com/maps/about/contribute/photosphere/
.. _Project Tango: https://www.google.com/atap/projecttango/
.. _planar tracking: http://mango.blender.org/development/planar-tracking-preview/
.. _Willow Garage: https://www.willowgarage.com/blog/2013/08/09/enabling-robots-see-better-through-improved-camera-calibration
.. _Lens Blur: http://googleresearch.blogspot.com/2014/04/lens-blur-in-new-google-camera-app.html
.. _SLAM: http://en.wikipedia.org/wiki/Simultaneous_localization_and_mapping
.. _calibrating robot-camera systems:
   http://rosindustrial.org/news/2014/9/24/industrial-calibration-library-update-and-presentation
.. _skinned control meshes: http://research.microsoft.com/en-us/projects/handmodelingfrommonoculardepth/
