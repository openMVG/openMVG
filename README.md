=====================================
OpenMVG (open Multiple View Geometry)
=====================================

![Logo](https://github.com/openMVG/openMVG/raw/master/logo/openMVG_Logo.png)

------------
Introduction
------------

[OpenMVG (Multiple View Geometry)](http://imagine.enpc.fr/~moulonp/openMVG/) "open Multiple View Geometry" is a library for computer-vision scientists and especially targeted to the Multiple View Geometry community. It is designed to provide an easy access to the classical problem solvers in Multiple View Geometry and solve them accurately.

The openMVG credo is: "Keep it simple, keep it maintainable". OpenMVG targets readable code that is easy to use and modify by the community.

All the features and modules are unit tested. This test driven development ensures that the code works as it should and enables more consistent repeatability. Furthermore, it makes it easier for the user to understand and learn the given features.

Specifically, this repository is aimed at integrating GPU-oriented features and matchers into the photogrammetry pipeline in order to allow for optimised reconstruction times. Currently, this repository uses CUDA >= 7.0. However, future plans may involve
using OpenCL/FPGA descriptors for faster detections and matching.

Current additional descriptors and matchers involve a CUDA implementation of LATCH and a CUDA implementation of a brute force Hamming matcher. To change the parameters of these matchers, there are files located under 
[cudaLATCH's NUM SM parameter](https://github.com/mdaiter/cudaLATCH/blob/cf05a8fdf19b83519e68cc0c184e334f83be18e5/params.hpp) and [bitMatcher's THRESHOLD parameter](https://github.com/mdaiter/openMVG/blob/custom/src/openMVG/matching_image_collection/gpu/params.hpp). These numbers allow the user to tune exactly how many CUDA SMs are used during the descriptor detection and matching, as well as the threshold for a possible match. WARNING: these numbers are multiplied by 512 to give a maximum
amount of keypoints to the CUDA FAST detector.

--------
Building
--------

See [BUILD](https://github.com/mdaiter/openMVG/raw/custom/BUILD) text file

Continuous integration:
 - linux 64 bits/GCC (Build + tests): [![Build Status](https://travis-ci.org/openMVG/openMVG.png?branch=develop)](https://travis-ci.org/openMVG/openMVG)
 - VStudio 2015 64 bits (Build): [![Build status](https://ci.appveyor.com/api/projects/status/3nv6rt41yxqx5v7i?svg=true)](https://ci.appveyor.com/project/pmoulon/openmvg)
 - Unit test coverage: [![Coverage Status](https://coveralls.io/repos/openMVG/openMVG/badge.png?branch=develop)](https://coveralls.io/r/openMVG/openMVG?branch=develop)

-------
License
-------

See [LICENSE MPL2](https://github.com/openMVG/openMVG/raw/master/license.openMVG) text file

-------
Authors
-------

See [Authors](https://github.com/openMVG/openMVG/raw/master/AUTHORS) text file

-------
Documentation
-------

See [documentation](http://openmvg.readthedocs.org/en/latest)

----------------
Acknowledgements
----------------

openMVG authors would like to thanks libmv authors for providing an inspiring 
base to design openMVG. Authors also would like to thanks Mikros Image [1] 
and LIGM-Imagine laboratory [2] for support and authorization to make this
library an opensource project.

[1] [http://www.mikrosimage.eu/](http://www.mikrosimage.eu/)
[2] [http://imagine.enpc.fr/](http://imagine.enpc.fr/)

---------
Contact
---------

openmvg-team[AT]googlegroups.com


---------
Citations
---------

If you find the library or some part of it useful, then following
publications are relevant:

[3] Moulon Pierre, Monasse Pascal and Marlet Renaud. ACCV 2012.
[Adaptive Structure from Motion with a contrario model estimation.](http://hal.archives-ouvertes.fr/index.php?halsid=1n2qdqiv2a0l5eq7qpos9us752&view_this_doc=hal-00769266&version=1)

[4] Moulon Pierre and Monasse Pascal. CVMP 2012.
[Unordered feature tracking made fast and easy.](http://hal.archives-ouvertes.fr/index.php?halsid=ggdarhl8cv1j6ohq2073eok8q3&view_this_doc=hal-00769267&version=1)

[5] Moisan Lionel, Moulon Pierre and Monasse Pascal. IPOL 2012.
[Automatic Homographic Registration of a Pair of Images, with A Contrario Elimination of Outliers.](http://dx.doi.org/10.5201/ipol.2012.mmm-oh)

[6] Moulon Pierre, Monasse Pascal and Marlet Renaud. ICCV 2013.
[Global Fusion of Relative Motions for Robust, Accurate and Scalable Structure from Motion.](http://imagine.enpc.fr/~moulonp/publis/iccv2013/index.html)

[7] Christopher Parker, Matthew Daiter, Kareem Omar, Gil Levi and Tal Hassner. ECCV Workshop 2016. [The CUDA LATCH Binary Descriptor: Because Sometimes Faster Means Better, Workshop on Local Features: State of the art, open problems and performance evaluation.](http://www.openu.ac.il/home/hassner/projects/LATCH/)

[8] Gil Levi and Tal Hassner. WACV, 2016. [LATCH: Learned Arrangements of Three Patch Codes](http://www.openu.ac.il/home/hassner/projects/LATCH/)

or cite it as:

```
  @misc{openMVG,
    author = "Pierre Moulon and Pascal Monasse and Renaud Marlet and Others",
     title = "OpenMVG. An Open Multiple View Geometry library.",
    howpublished = "\url{https://github.com/openMVG/openMVG}",
  }
```
