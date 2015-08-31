***************************
SfM: Structure from Motion
***************************

Structure from Motion computes an external camera pose per image (the motion) and a 3D point cloud (the structure) from:

- images,
- some intrinsic camera parameters,
- corresponding geometric valid features accross images.

.. figure:: imagesInput.png
   :align: center
   
   Figure : Input images, estimated camera location and structure.

openMVG SfM tools
====================

  - 2 **Structure from Motion (SfM) pipeline**:

    - an Incremental Structure from Motion chain [ACSfM]_ (ACCV 2012),
    - a Global Structure from Motion chain [GlobalACSfM]_ (ICCV 2013),

  - 1 **Structure from known Motion (SfM) pipeline**:

    - Structure computation from known camera poses and features.

  - **tools** to visualize:

    - features,
    - photometric/geometric matches correspondences,
    - features tracks.

  - **tools to export to existing Multiple View Stereovision (MVS) pipeline**:

    - [PMVS]_, CMPMVS.

OpenMVG SfM pipelines
======================

OpenMVG SfM pipelines run as a 4 step process:

1. Image listing
-----------------------

.. toctree::
   :maxdepth: 1

   ./SfMInit_ImageListing.rst


2. Image description computation
----------------------------------------------

.. toctree::
   :maxdepth: 1

   ./ComputeFeatures.rst

3. Corresponding images and correspondences computation
---------------------------------------------------------------------

.. toctree::
   :maxdepth: 1

   ./ComputeMatches.rst

4. SfM solving (2 methods)
---------------------------

.. toctree::
   :maxdepth: 1

   ./IncrementalSfM.rst
   ./GlobalSfM.rst

5. Optional further processing 
----------------------------------------

.. toctree::
   :maxdepth: 1

   ./ComputeSfM_DataColor.rst
   ./ComputeStructureFromKnownPoses.rst
   ./ExportUndistortedImages.rst


5. Optional further processing (3rd party)
-------------------------------------------

.. toctree::
   :maxdepth: 1

   ./MVS.rst

**You can either run by hand all the process or use pre-defined python scripts (that are using some default options).**

OpenMVG SfM pipelines demo
===========================

A complete ready to use tutorial demo is exported in your build directory. It clones an image dataset and run the SfM pipelines on it:

- openMVG_Build/software/SfM/tutorial_demo.py

In order to use easily the Sequential or the Global pipeline, ready to used script are exported in your build directory:

- openMVG_Build/software/SfM/SfM_SequentialPipeline.py
- openMVG_Build/software/SfM/SfM_GlobalPipeline.py

To use them simply run:

.. code-block:: c++

  $ cd openMVG_Build/software/SfM/
  $ python SfM_SequentialPipeline.py [full path image directory] [resulting directory]
  $ python SfM_SequentialPipeline.py ~/home/user/data/ImageDataset_SceauxCastle/images ~/home/user/data/ImageDataset_SceauxCastle/Castle_Incremental_Reconstruction

  $ python SfM_GlobalPipeline.py [full path image directory] [resulting directory]

More details about openMVG tools
=====================================

To know more about each tool visit the following link and read the doc below:

.. toctree::
   :maxdepth: 1

   ./SfMInit_ImageListing.rst
   ./ComputeFeatures.rst
   ./ComputeMatches.rst
   ./IncrementalSfM.rst
   ./GlobalSfM.rst
   ./ComputeSfM_DataColor.rst
   ./ComputeStructureFromKnownPoses.rst
   ./ExportUndistortedImages.rst

.. toctree::
   :maxdepth: 1

   ./SfM_OutputFormat.rst
   ./MVS.rst

PS: We strongly advise to use a 3 directories based data organisation structure

* **images**

  - your image sequence.

* **matches**

  * directory used to store image information, images features, descriptors and matches information.

* **outReconstruction**

  * directory used to store the SfM result and process log.

.. 
  1. Image & view listing:

    - describe images parameters:
       - image name,
       - image size,
       - internal camera calibration information (intrinsic parameters) (if any).

  2. Features & descriptors extraction:

      - Describe each view with feature points and their corresponding descriptors (e.g., SIFT).

  3. Putative matches & geometric correspondences computation,

    - Geometric feature matching across photo collection:

      - descriptor matching between image pairs allows to build an initial corresponding photometric feature graph.
      - this graph is then geometrically filtered using robust estimation of fundamental or essential or homography matrix based.

  4. SfM solving:

    - the corresponding features graph is send to the chosen SfM pipeline and it computes the scene and camera motion.
