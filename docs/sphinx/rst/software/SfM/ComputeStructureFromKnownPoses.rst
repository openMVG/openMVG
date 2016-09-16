
********************************************
openMVG_main_ComputeStructureFromKnownPoses
********************************************

This application compute corresponding features and robustly triangulate them according the geometry of the known camera intrinsics & poses.

Algorithm of the application

.. code-block:: c++

  Require: internal + external camera calibration
  Require: image description regions (features + descriptors)
  Ensure: 3D point cloud
  compute image visibility
    list all the pair that share common visual content
      - camera frustum based
      - or structure visbility (SfM tracks) based
  list triplets of view from pairs
   for each triplets compute 3 view tracks
     if tracks triangulable add correspondences to p
  link 3 views validated matches (p) as tracks
    robustly triangulate them

Information and usage
========================

The chain is designed to run on a sfm_data.json file and some pre-computed matches.
The sfm_data file should contains:
- valid view with some defined intrinsics and camera poses,
- (optional existing structure).

  .. code-block:: c++
  
    $ openMVG_main_ComputeStructureFromKnownPoses -i Dataset/out_Reconstruction/sfm_data.json -o Dataset/out_Reconstruction/robustFitting.json

Arguments description:

**Required parameters:**

  - **[-i|--input_file]**

    - a SfM_Data file with valid intrinsics and poses and optional structure

  - **[-m|--matchdir]**

    - path were image descriptions were stored

  - **[-o|--outdir]**

    - path where the updated scene data will be stored

**Optional parameters:**

  - **[-f|--match_file]**

    - path to a matches file (pairs of the match files will be listed and used)

  - **[-p|--pair_file]**

    - path to a pairs file (only those pairs will be considered to compute the structure)
      The pair file is a list of view indexes, one pair on each line

  - **[-b|--bundle_adjustment]**

    - perform a bundle adjustment on the scene (OFF by default)

  - **[-r|--residual_threshold]**

    - maximal pixels reprojection error that will be considered for triangulations (4.0 by default)


