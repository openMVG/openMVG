
*************************************
openMVG_main_IncrementalSfM
*************************************

The [ACSfM]_ SfM is an evolution of the implementation used for the paper "Adaptive Structure from Motion with a contrario model estimation"  published at ACCV 2012.

The incremental pipeline is a growing reconstruction process.
It starts from an initial two-view reconstruction (the seed) that is iteratively extended by adding new views and 3D points, using pose estimation and triangulation.
Due to the incremental nature of the process, successive steps of non-linear refinement, like Bundle Adjustment (BA) and Levenberg-Marquardt steps, are performed to minimize the accumulated error (drift).

Algorithm of the Incremental/Sequential Structure from Motion

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

The chain is designed to run on a sfm_data.json file and some pre-computed matches.

  .. code-block:: c++
  
    $ openMVG_main_IncrementalSfM -i Dataset/matches/sfm_data.json -m Dataset/matches/ -o Dataset/out_Incremental_Reconstruction/

openMVG_main_IncrementalSfM displays to you some initial pairs that share an important number of common point.
  **Please select two image index that are convergent and the 3D reconstruction will start.
  The initial pair must be choosen with numerous correspondences while keeping a wide enough baseline.**

Arguments description:

**Required parameters:**

  - **[-i|--input_file]**

    - a SfM_Data file

  - **[-m|--matchdir]**

    - path were geometric matches were stored
  
  - **[-o|--outdir]**

    - path where the output data will be stored

**Optional parameters:**

  - **[-a|--initialPairA NAME]**

    - the filename image to use (i.e. 100_7001.JPG)

  - **[-b|--initialPairB NAME]**

    - the filename image to use (i.e. 100_7002.JPG)

  - **[-c|--camera_model]**

    - The camera model type that will be used for views with unknown intrinsic:

      - 1: Pinhole
      - 2: Pinhole radial 1
      - 3: Pinhole radial 3 (default)
      - 4: Pinhole radial 3 + tangential 2
      - 5: Pinhole fisheye

  - **[-f|--refineIntrinsics]**
      User can control exactly which parameter will be considered as constant/variable and combine them by using the '|' operator.
      
    - ADJUST_ALL -> refine all existing parameters (default)
    - NONE -> intrinsic parameters are held as constant
    - ADJUST_FOCAL_LENGTH -> refine only the focal length
    - ADJUST_PRINCIPAL_POINT -> refine only the principal point position
    - ADJUST_DISTORTION -> refine only the distortion coefficient(s) (if any)
    
    - NOTE Options can be combined thanks to '|':

      - ADJUST_FOCAL_LENGTH|ADJUST_PRINCIPAL_POINT
        -> refine the focal length & the principal point position
      
      - ADJUST_FOCAL_LENGTH|ADJUST_DISTORTION
        -> refine the focal length & the distortion coefficient(s) (if any)
      
      - ADJUST_PRINCIPAL_POINT|ADJUST_DISTORTION
        -> refine the principal point position & the distortion coefficient(s) (if any)

