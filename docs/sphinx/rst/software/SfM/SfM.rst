*******************************
SfM: Structure-from-Motion
*******************************

Structure from Motion computes an external camera pose per image (the motion) and a 3D point cloud (the structure) representing the pictured scene.
Inputs are images and internal camera calibration information (intrinsic parameters).
Feature points are detected in each image (e.g., SIFT) and matched between image pairs.
There are two main approaches to correlate detected features and solve the SfM problem: the incremental pipeline and the global method. 

.. figure:: imagesInput.png
   :align: center
   
   Figure : Input images, estimated camera location and structure.


openMVG proposes a customizable implementation of an Incremental Structure from Motion chain, it is the implementation used for the paper "Adaptive Structure from Motion with a contrario model estimation" [ACSfM]_ published at ACCV 2012.

.. figure:: structureFromMotion.png
   :align: center

   Figure: Structure from Motion illustration, from pictures to 3D.

The incremental pipeline is a growing reconstruction process (i.e algorithm 2).
It starts from an initial two-view reconstruction (the seed) that is iteratively extended by adding new views and 3D points, using pose estimation and triangulation.
Due to the incremental nature of the process, successive steps of non-linear refinement, like Bundle Adjustment (BA) and Levenberg-Marquardt steps, are performed to minimize the accumulated error (drift).

The general feature correspondence and SfM processes are described in algorithms 1 and 2.
The first algorithm outputs pairwise correspondences that are consistent with the estimated fundamental matrix.
The initial pair must be chosen with numerous correspondences while keeping a wide enough baseline.
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

The a contrario Structure from Motion chain take as input a sequence of distortion corrected
images and the intrinsic associated calibration matrix K. This calibration matrix is stored as a
"K.txt" file as a raw ascii 3 x 3 matrix in the same directory as the pictures.
Using a 3 directories based data organisation structure is suggested:

* images

  - your image sequence,
  - K.txt

* matches

  * the points and matches information will be saved here

* outReconstruction

  * directory where result and log of the 3D reconstruction will be exported

1. Point matching:

  The first step consists in computing relative image matches (i.e algorithm 2): You have to use the openMVG_main_computeMatches software in the software/SfM openMVG module.

  .. code-block:: c++

    $ openMVG_main_computeMatches -i /home/pierre/Pictures/Dataset/images -e *.JPG -o /home/pierre/Pictures/Dataset/matches

  Arguments are the following:

  - -i|-imadir the path where image are stored.
  - -e|-ext image extension i.e "*.jpg" or "*.png". Case sensitive.
  - -o|-outdir path where features, descriptors, putative and geometric matches will be exported.
  - -r|-distratio optional argument (Nearest Neighbor distance ratio, default value is set to 0.6).
  - -s|-octminus1 optional argument (Use the octave -1 option of SIFT or not, default value is set to false: 0).

  Once matches have been computed you can, at your choice, display detected points, matches or
  start the 3D reconstruction.

2. Point, matching visualization:

  Three softwares are available to display:

  * **Detected keypoints**: openMVG_main_exportKeypoints
  *	**Putative, Geometric matches**: openMVG_main_exportMatches
  *	**Tracks**: openMVG_main_exportTracks


3. SfM, 3D structure and camera calibration:

  The main binary in order to run the SfM process is openMVG_main_IncrementalSfM, it use previous
  computed data and is implemented as explained in algorithm 2.

  .. code-block:: c++

    $ openMVG_main_IncrementalSfM -i /home/pierre/Pictures/Dataset/images/ -m /home/pierre/Pictures/Dataset/matches/ -o /home/pierre/Pictures/Dataset/outReconstruction /

  openMVG_main_IncrementalSfM displays to you some initial pairs that share an important number of common point.
  Please select two image indexes that are convergent and the 3D reconstruction will start.


.. [ACSfM] Adaptive structure from motion with a contrario model estimation.
    Pierre Moulon, Pascal Monasse, and Renaud Marlet.
    In ACCV, 2012.

.. [PMVS] Accurate, dense, and robust multi-view stereopsis.
    Yasutaka Furukawa and Jean Ponce.
    IEEE Trans. on Pattern Analysis and Machine Intelligence, 32(8):1362-1376, 2010.
