
Multiple View Stereovision
====================================

Multiple View Stereo-vision
*****************************

Once camera position and orientation have been computed, Multiple View Stereo-vision algorithms could be used 
to compute a dense scene representation.


OpenMVG exports [PMVS]_ and [CMPMVS]_ ready to use project (images, projection matrices and pmvs_options.txt files).

.. figure:: resultOutput.png
   :align: center

   Figure: Multiple View Stereo-vision point cloud densification on the estimated scene using [PMVS]_.

.. [PMVS] Accurate, dense, and robust multi-view stereopsis.
    Yasutaka Furukawa and Jean Ponce.
    IEEE Trans. on Pattern Analysis and Machine Intelligence, 32(8):1362-1376, 2010.

.. [CMPMVS] Multi-View Reconstruction Preserving Weakly-Supported Surfaces, M. Jancosek, T. Pajdla, IEEE Conference on Computer Vision and Pattern Recognition 2011.


