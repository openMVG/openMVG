
***************************************************************
Geometric consistent matches for unordered image collection
***************************************************************

Computing relative image matches is often required in computer vision.
OpenMVG module provide a binary to compute putative photometric matches and filter them according a choosen robust geometric filter.

Computation of geometry-consistent pairwise correspondences:
--------------------------------------------------------------

.. code-block:: c++

  Require: image set
  Ensure: pairwise geometrically consistent point correspondences 
  a. Compute putative photometric matches:
    detect features in each image and build their descriptor
    descriptors matching (brute force or approximate nearest neighbor)
  b. Filter geometric-consistent matches:
    robust estimation of the choosen geometric model

In order to compute pairwise geometric consistent matches use the following binary:

  .. code-block:: c++

    $ openMVG_main_computeMatches [opt. args] -i Dataset/images/ -o Dataset/matches/

  Arguments are the following:

  - **-i|-imadir** the path where image are stored.
  - **-o|-outdir** path where features, descriptors, putative and geometric matches will be exported.

  - Optional arguments
 
    - **-r|-distratio** optional argument (Nearest Neighbor distance ratio, default value is set to 0.6).

      * -r 0.6 is restrictive, -r 0.8 is less restrictive

    - **-s|-octminus1** optional argument
    
      - -s 0 Use octave 0
      - -s 1 Use the octave -1 option of SIFT Upscale the image x2.

    - **-p|-peakThreshold** optional argument (Peak threshold for SIFT detection).

        - -p [0.04 (default) to 0.01] range of possible value. 0.01 detects more SIFT keypoints.

    - **-g|-geometricModel** type of model used for robust estimation from the photometric putative matches

      - -g [f(default) or e or h] 
      - -g f Fundamental matrix filtering
      - -g e Essential matrix filtering (all the image must have the same known focal length)
      - -g h Homography matrix filtering

  Once matches have been computed you can, at your choice, you can display detected points, matches or
  start the 3D reconstruction.


Point, matching, tracks visualization:
----------------------------------------

* **Detected keypoints**: openMVG_main_exportKeypoints
*	**Putative, Geometric matches**: openMVG_main_exportMatches
*	**Tracks**: openMVG_main_exportTracks


