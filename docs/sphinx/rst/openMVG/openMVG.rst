############################
openMVG
############################
  
The openMVG library core module
====================================

The core openMVG library provides data-structure and functionnalities to:
 - open, write, manipulate images,
 - manipulate and match local descriptions of images (features, descriptors, matching),
 - compute corresponding points between image pairs and collections,
 - estimate Multiple View Geometry relations between image pairs,
 - refine camera and structure parameters (Bundle Adjustment),
 - handle points, cameras and projections matrices.

Use within the robust_estimation libraries it allows estimation of Multiple View Geometry relations from point matches.

Here the list of libraries that:

.. toctree::
  :maxdepth: 1
  
  numeric/numeric.rst
  features/features.rst
  cameras/cameras.rst
  multiview/multiview.rst
  image/image.rst
  matching/matching.rst
  tracks/tracks.rst
  robust_estimation/robust_estimation.rst
  bundle_adjustment/bundle_adjustment.rst

