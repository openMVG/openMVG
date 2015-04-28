***************************
SfM: Structure-from-Motion
***************************

SfM in a nutshell
===================

Structure from Motion computes an external camera pose per image (the motion) and a 3D point cloud (the structure) from:

- images,
- some intrinsic camera parameters,
- corresponding geometric valid features accross images.

.. figure:: imagesInput.png
   :align: center
   
   Figure : Input images, estimated camera location and structure.

openMVG tools
====================

  - 2 **Structure from Motion pipeline**:

    - an Incremental Structure from Motion chain [ACSfM]_ (ACCV 2012),
    - a Global Structure from Motion chain [GlobalACSfM]_ (ICCV 2013).

  - **tools** to visualize:

    - features,
    - photometric/geometric matches correspondences,
    - features tracks.

  - **export to existing Multiple View Stereovision pipeline**:

    - [PMVS]_, CMPMVS.

OpenMVG SfM pipelines
======================

OpenMVG SfM pipelines are used in a three step process:

1. Intrinsic image analysis & view listing:

  - describe images with internal camera calibration information (intrinsic parameters) if any.

2. features extraction & geometric correspondences computation,

  - Geometric feature matching across photo collection:

    - feature points and corresponding descriptors are detected in each image (e.g., SIFT).
    - descriptor matching between image pairs allows to build an initial corresponding photometric feature graph.
    - this graph is then geometrically filtered using robust estimation of fundamental or essential or homography matrix based.

3. SfM solving:

  - the corresponding features graph is send to the chosen SfM pipeline and it computes the scene and camera motion.


Structure from Motion chains usage
=====================================

Using a 3 directories based data organisation structure for openMVG SfM pipeline is suggested:

* **images**

  - your image sequence.

* **matches**

  * directory used to store image information, images features, descriptors and matches information.

* **outReconstruction**

  * directory used to store the SfM result and process log.

To know more about each tool visit the following link and read the doc below:

.. toctree::
   :maxdepth: 1

   ./intrinsicGroups.rst
   ./geometricMatches.rst
   ./incrementalSfM.rst
   ./globalSfM.rst
   ./SfM_OutputFormat.rst
   ./MVS.rst
