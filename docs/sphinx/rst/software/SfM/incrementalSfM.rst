
*************************************
Incremental Structure from Motion
*************************************

The [ACSfM]_ SfM is an evolution of the implementation used for the paper "Adaptive Structure from Motion with a contrario model estimation"  published at ACCV 2012.

The incremental pipeline is a growing reconstruction process.
It starts from an initial two-view reconstruction (the seed) that is iteratively extended by adding new views and 3D points, using pose estimation and triangulation.
Due to the incremental nature of the process, successive steps of non-linear refinement, like Bundle Adjustment (BA) and Levenberg-Marquardt steps, are performed to minimize the accumulated error (drift).

Incremental Structure from Motion

.. code-block:: c++

  Require: internal camera calibration (possibly from EXIF data)
  Require: pairwise geometry consistent point correspondences
  Ensure: 3D point cloud
  Ensure: camera poses
  compute correspondence tracks t
  compute connectivity graph G (1 node per view, 1 edge when enough matches)
  pick an edge e in G with sufficient baseline
  * robustly estimate essential matrix from images of e
  triangulate validated tracks, which provides an initial reconstruction
  contract edge e
  while G contains an edge do
    pick edge e in G that maximizes union(track(e),3D points)
    * robustly estimate pose (external orientation/resection)
    triangulate new tracks
    contract edge e
    perform bundle adjustment
  end while

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

As no intrinsic is supposed to be known the fundamental matrix is used to filter putative matches.

.. code-block:: c++

  $ openMVG_main_computeMatches -g f -i Dataset/images/ -o Dataset/matches/

Refers to the SfM general help to know more about optional argument to detect more points per images pair.

**Tips**

  ``-g f -p 0.01 -r 0.8 -i ...`` helps having more matches.

3. Incremental Structure from Motion
--------------------------------------------------

- If you want refine intrinsics (focal, principal point and radial distortion) for each focal group

  .. code-block:: c++
  
    $ openMVG_main_IncrementalSfM -i Dataset/images/ -m Dataset/matches/ -o Dataset/outReconstruction/
  
- If you want only refine the focal (to use with image were the distortion have been already removed)

  .. code-block:: c++
  
    $ openMVG_main_IncrementalSfM -d 0 -i Dataset/images/ -m Dataset/matches/ -o Dataset/outReconstruction/

openMVG_main_IncrementalSfM displays to you some initial pairs that share an important number of common point.
  **Please select two image indexes that are convergent and the 3D reconstruction will start.
  The initial pair must be choosen with numerous correspondences while keeping a wide enough baseline.**

