
*************************************
Multiple View Stereovision
*************************************

Once camera position and orientation have been computed, Multiple View Stereo-vision algorithms could be used 
to compute a dense scene representation, such as:

- dense point cloud (PMVS),
- surface and texture (UMVS, CMPMVS).


Export to PMVS/CMVS
========================

OpenMVG exports [PMVS]_ ready to use project (images, projection matrices and pmvs_options.txt files).

Once a 3D calibration have been computed you can convert the SfM_Ouput files to a PMVS project.

.. code-block:: c++

  $ openMVG_main_openMVG2PMVS -i Dataset/outReconstruction/SfM_Output/ -o Dataset/outReconstruction/SfM_Output/
  $ pmvs Dataset/outReconstruction/SfM_Output/PMVS/ pmvs_options.txt

.. figure:: resultOutput.png
   :align: center

   Figure: Multiple View Stereo-vision point cloud densification on the estimated scene using [PMVS]_.

In order to use CMVS for large scene openMVG2PMVS export also the scene in the Bundler output format.

.. code-block:: c++

  $ openMVG_main_openMVG2PMVS -i Dataset/outReconstruction/SfM_Output/ -o Dataset/outReconstruction/SfM_Output/
  $ cmvs Dataset/outReconstruction/SfM_Output/PMVS/ [MaxImageCountByCluster=100]
  $ cmvs Dataset/outReconstruction/SfM_Output/PMVS/ 30
  $ genOption Dataset/outReconstruction/SfM_Output/PMVS/
  $ sh Dataset/outReconstruction/SfM_Output/PMVS/pmvs.sh


Export to MVE (Multi-View Environment)
=========================================

`MVE <http://www.gris.informatik.tu-darmstadt.de/projects/multiview-environment>`_ can import an openMVG SfM scene and use it to create dense depth map and complete dense 3D models.

.. code-block:: c++

  $ makescene SfmOutput_DIR MVE_SCENE_DIR
  $ dmrecon -s2 MVE_SCENE_DIR
  $ scene2pset -ddepth-L2 -iundist-L2 -n -s -c MVE_SCENE_DIR OUTPUT.ply
  $ fssr_octree OUTPUT.ply OUTPUT.fssr
  $ fssr_surface -t10 -c100000 OUTPUT.fssr surface.ply
  
  Call any tool without arguments to see the help.
  
You will need to compile MVE tools and `FSSR <http://www.gris.informatik.tu-darmstadt.de/projects/floating-scale-surface-reco/>`_.

Export to MVS Texturing
=======================

If you don't want to use the full MVE pipeline but only `MVS Texturing <http://www.gris.tu-darmstadt.de/projects/mvs-texturing/>`_ [Waechter2014]_ to project a set of oriented images on a mesh, one solution is to use the openMVG_main_openMVG2MVSTEXTURING binary. This binary converts your SfM_Data file into one format used by MVS Texturing. In addition, you may need to undistort your images with openMVG_main_ExportUndistortedImages as it's not handled by the openMVG_main_openMVG2MVSTEXTURING tool.

Export to CMPMVS
========================

OpenMVG exports [CMPMVS]_ ready to use project (images, projection matrices and ini configuration file).

.. code-block:: c++

  $ openMVG_main_openMVG2CMPMVS -i Dataset/outReconstruction/SfM_Output/ -o Dataset/outReconstruction/SfM_Output/
