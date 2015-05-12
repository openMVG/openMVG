############################
patented
############################


The openMVG licence is MPL2 but some code used to provide interesting functionnalities have different licences or are subject to patent.

Such code takes place in the patented directory.

SIFT
=============

In order to compute corresponding points between images pairs, openMVG uses natively SIFT keypoints and their associated descriptors.
The used code is a subset of the [VLFEAT]_ library.

You can replace this keypoints and descriptors provided by any version of your choice to use openMVG in a non-research context.
Suggestions for features points:

  - CORNERS: HARRIS, FAST,
  - REGIONS: MSER,
  - BLOBS: AKAZE.

Descriptors:
  - BINARY: M-LDB (see AKAZE paper), BRIEF, Nested shape descriptor,
  - FLOATING POINT: DAISY, LIOP descriptors.


