
*******************************
Global Structure from Motion
*******************************

[GlobalACSfM]_ is based on the paper "Global Fusion of Relative Motions for Robust, Accurate and Scalable Structure from Motion."  published at ICCV 2013.

Multi-view structure from motion (SfM) estimates the position and orientation of pictures in a common 3D coordinate frame. When views are treated incrementally, this external calibration can be subject to drift, contrary to global methods that distribute residual errors evenly. Here the method propose a new global calibration approach based on the fusion of relative motions between image pairs. 

.. code-block:: c++

  Require: internal camera calibration (possibly from EXIF data)
  Require: pairwise geometry consistent point correspondences
  Ensure: 3D point cloud
  Ensure: camera poses

  compute relative pairwise rotations
  detect and remove false relative pairwise rotations
    - using composition error of triplet of relative rotations
  compute the global rotation
    - using a dense least square and approximated rotations
  compute relative translations
    - using triplet of views for stability and colinear motion support
  compute the global translation
    - integration of the relative translation directions using a l-∞ method. 
  final structure and motion
    - link tracks validated per triplets and compute global structure by triangulation,
    - refine estimated parameter in a 2 step Bundle Adjustment
      - refine structure and translations
      - refine structure and camera parameters (rotations, translations).

Information and usage
========================

The chain is designed for only one intrinsic group, so all the pictures provided to the chain must have the same focal length.

The chain is used in a three step process: 
  1. Intrinsic analysis,
  2. Geometric features correspondences computation,
  3. Global Structure from Motion.

1. Intrinsic analysis
-----------------------

See intrinsic SfM entry.

2. Geometric features correspondences computation
--------------------------------------------------

As approximative or precise focal length is know we can use the essential matrix to filter putative matches.

.. code-block:: c++

  $ openMVG_main_computeMatches -g e -p 0.01 -i Dataset/images/ -o Dataset/matches/

**Tips**

  As a dense image network is required to perform global SfM it is required to detect more SIFT points per image to ensure a high probability of matching. To have even a more dense image connection graph you can use ``-g e -p 0.01 -r 0.8 -s 1 -i ...``

3. Global Structure from Motion
--------------------------------------------------

.. code-block:: c++

  $ openMVG_main_GlobalSfM -i Dataset/images/ -m Dataset/matches/ -o Dataset/outGlobalSfM




