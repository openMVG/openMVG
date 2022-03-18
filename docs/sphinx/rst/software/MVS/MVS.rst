
*************************************
MVS: Multiple View Stereovision
*************************************

Once camera position and orientation have been computed by OpenMVG, Multiple View Stereo-vision algorithms could be used to compute a dense representation of the scene by generating point clouds or mesh surfaces.

Since OpenMVG does not have itself a MVS module it proposes exporter to various existing solutions:

.. toctree::
   :maxdepth: 1

   ./AgisoftPhotoscan.rst
   ./CMPMVS.rst
   ./MVE.rst
   ./OpenMVS.rst
   ./PMVS.rst
   ./Colmap.rst
