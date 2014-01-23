*******************************
SfM: Structure-from-Motion
*******************************

Structure from Motion computes an external camera pose per image (the motion) and a 3D point cloud (the structure) representing the pictured scene.
Inputs are images and internal camera calibration information (intrinsic parameters).
Feature points are detected in each image (e.g., SIFT) and matched between image pairs and then the SfM pipeline compute the scene and camera motion.
There are three main approaches to solve the SfM problem:
  - the incremental/sequential pipeline,
  - the hierarchical pipeline,
  - the global one.

.. figure:: imagesInput.png
   :align: center
   
   Figure : Input images, estimated camera location and structure.


openMVG proposes a customizable implementation of an Incremental Structure from Motion chain, it is an evolution of the implementation used for the paper "Adaptive Structure from Motion with a contrario model estimation" [ACSfM]_ published at ACCV 2012.

.. figure:: structureFromMotion.png
   :align: center

   Figure: Structure from Motion illustration, from pictures to 3D.

The incremental pipeline is a growing reconstruction process (i.e algorithm 2).
It starts from an initial two-view reconstruction (the seed) that is iteratively extended by adding new views and 3D points, using pose estimation and triangulation.
Due to the incremental nature of the process, successive steps of non-linear refinement, like Bundle Adjustment (BA) and Levenberg-Marquardt steps, are performed to minimize the accumulated error (drift).

The general feature correspondence and SfM processes are described in algorithms 1 and 2.
The first algorithm outputs pairwise correspondences that are consistent with the estimated fundamental matrix.
The initial pair must be choosen with numerous correspondences while keeping a wide enough baseline.
The second algorithm takes these correspondences as input and yields a 3D point cloud as well as the camera poses.
Steps marked with a star (*) are estimated within the a contrario framework.
This allows critical thresholds to be automatically adapted to the input images (and remove the choice of an empiric T threshold value).

Algorithm 1 Computation of geometry-consistent pairwise correspondences

.. code-block:: c++

	Require: image set
	Ensure: pairwise point correspondences that are geometrically consistent
	Compute putative matches:
	detect features in each image and build their descriptor
	match descriptors (brute force or approximate nearest neighbor)
	Filter geometric-consistent matches:
	* estimate fundamental matrix F

Algorithm 2 Incremental Structure from Motion

.. code-block:: c++

	Require: internal camera calibration (matrix K, possibly from EXIF data)
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

	
Once camera position and orientation have been computed, Multiple View Stereo-vision algorithms could be used 
to compute a dense scene representation .
OpenMVG exports a [PMVS]_ ready to use project (images, projection matrices and pmvs_options.txt files).

.. figure:: resultOutput.png
   :align: center

   Figure: Multiple View Stereo-vision point cloud densification on the estimated scene using [PMVS]_.

Tools to draw detected Keypoints, putative matches, geometric matches and tracks are provided
along the SfM binary directory.

Structure from Motion chain usage
=====================================

The a contrario Structure from Motion chain take as input an image collection.
Using a 3 directories based data organisation structure is suggested:

* images

  - your image sequence.

* matches

  * the image information (lists.txt), images points and matches information will be saved here

* outReconstruction

  * directory where result and log of the 3D reconstruction will be exported

1. Intrinsic analysis and export:

  The process export in outputDirectory/lists.txt file the extracted camera intrinsics in three way:
    - no information can be extracted (only the principal point position is exported)
    - the image contain EXIF Jpeg approximated focal length
      - if the camera sensor is saved in the openMVG database the focal length and principal point is exported
      - else no focal camera can be computed, only the image size is exported.
  The focal is computed as follow (with a database of knowed camera sensor size in mm for different camera model):
  double ccdw = datasheet._sensorSize; // In mm
  focal = max ( width, height ) * focalmm / ccdw;

  .. code-block:: c++
    $ cd ./software/SfM/Release/
    $ openMVG_main_CreateList [-i|--imageDirectory] [-d|--sensorWidthDatabase] [-o|--outputDirectory] [-f|--focal]

  - Usage of the automatic chain (with JPEG images)
  
  .. code-block:: c++
  
    $ openMVG_main_CreateList -i Dataset/images/ -o Dataset/matches/ -d ./software/SfM/cameraSensorWidth/cameraGenerated.txt

  - If all the camera have the same focal length and you know it exactly
  
  .. code-block:: c++
  
    $ openMVG_main_CreateList Dataset/images/ -o Dataset/matches/ -f FOCAL_LENGTH(in mm, i.e 2750)

2. Point matching:

  The first step consists in computing relative image matches (i.e algorithm 2): You have to use the openMVG_main_computeMatches software in the software/SfM openMVG module.

  .. code-block:: c++

    $ openMVG_main_computeMatches -i Dataset/images/ -o Dataset/matches/

  Arguments are the following:

  - -i|-imadir the path where image are stored.
  - -o|-outdir path where features, descriptors, putative and geometric matches will be exported.
  - -r|-distratio optional argument (Nearest Neighbor distance ratio, default value is set to 0.6).
  - -s|-octminus1 optional argument (Use the octave -1 option of SIFT or not, default value is set to false: 0).

  Once matches have been computed you can, at your choice, display detected points, matches or
  start the 3D reconstruction.

3. Point, matching visualization:

  Three softwares are available to display:

  * **Detected keypoints**: openMVG_main_exportKeypoints
  *	**Putative, Geometric matches**: openMVG_main_exportMatches
  *	**Tracks**: openMVG_main_exportTracks


4. SfM, 3D structure and camera calibration:

  The main binary in order to run the SfM process is openMVG_main_IncrementalSfM, it use previous
  computed data and is implemented as explained in algorithm 2.

  - If you want refine intrinsics (focal, principal point and radial distortion) for each focal group
  .. code-block:: c++
  
    $ openMVG_main_IncrementalSfM -i Dataset/images/ -m Dataset/matches/ -o Dataset/outReconstruction/
  
  - If you want only refine the focal (to use with image were the distortion have been already removed)
  .. code-block:: c++
  
  
    $ openMVG_main_IncrementalSfM -i Dataset/images/ -m Dataset/matches/ -o Dataset/outReconstruction/ -d 0

  openMVG_main_IncrementalSfM displays to you some initial pairs that share an important number of common point.
  Please select two image indexes that are convergent and the 3D reconstruction will start.


.. [ACSfM] Adaptive structure from motion with a contrario model estimation.
    Pierre Moulon, Pascal Monasse, and Renaud Marlet.
    In ACCV, 2012.

.. [PMVS] Accurate, dense, and robust multi-view stereopsis.
    Yasutaka Furukawa and Jean Ponce.
    IEEE Trans. on Pattern Analysis and Machine Intelligence, 32(8):1362-1376, 2010.
