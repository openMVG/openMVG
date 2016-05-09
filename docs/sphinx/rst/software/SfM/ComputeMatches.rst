
***************************************************************
openMVG_main_ComputeMatches
***************************************************************

This binary compute images that have a visual overlap. Using image descriptions computed by **openMVG_main_ComputeFeatures**, we establish the corresponding putative photometric matches and filter the resulting correspondences using some robust geometric filters.

.. code-block:: c++

  $ openMVG_main_ComputeMatches -i [..\matches\sfm_data.json] -o [...\matches]

Arguments description:

**Required parameters:**

  - **[-i|--input_file]**

    - a SfM_Data file

  - **[-o|--out_dir path]**

    - path were putative and geometric matches will be stored

**Optional parameters:**
 
  - **[-f|--force: Force to recompute data]**

    - 0: (default) reload previously computed data (useful when you have kill the process and want to continue to compute)
    - 1: useful when you change have changed a command line parameter, force recomputing and re-saving.

  - **[-r|-ratio]**

    - (Nearest Neighbor distance ratio, default value is set to 0.8).
        Using 0.6 is more restrictive => provides less false positive.

  - **[-g|-geometric_model]**

    - type of model used for robust estimation from the photometric putative matches

      - f: Fundamental matrix filtering
      - e: Essential matrix filtering
      - h: Homography matrix filtering

  - **[-n|--nearest_matching_method]**

    - AUTO: auto choice from regions type,
    - For Scalar based descriptor you can use:
    
      - BRUTEFORCEL2: BruteForce L2 matching for Scalar based regions descriptor,
      - ANNL2: Approximate Nearest Neighbor L2 matching for Scalar based regions descriptor,
      - CASCADEHASHINGL2: L2 Cascade Hashing matching,
      - FASTCASCADEHASHINGL2: (default).
          L2 Cascade Hashing with precomputed hashed regions,
          (faster than CASCADEHASHINGL2 but use more memory).
    - For Binary based descriptor you must use:
    
      - BRUTEFORCEHAMMING: BruteForce Hamming matching for binary based regions descriptor,

  - **[-v|--video_mode_matching]**
  
    - (sequence matching with an overlap of X images)

      - X: with match 0 with (1->X), ...]
      - 2: will match 0 with (1,2), 1 with (2,3), ...
      - 3: will match 0 with (1,2,3), 1 with (2,3,4), ...]


  - **[-l|--pair_list]**

    - file that explicitly list the View pair that must be compared
     
Once matches have been computed you can, at your choice, you can display detected, matches as SVG files:

* **Detected keypoints**: openMVG_main_exportKeypoints
*	**Putative, Geometric matches**: openMVG_main_exportMatches
*	**Tracks**: openMVG_main_exportTracks

**Or start the 3D reconstruction:**

.. toctree::
   :maxdepth: 1

   ./IncrementalSfM.rst
   ./GlobalSfM.rst

