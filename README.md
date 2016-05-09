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

--------
Building
--------

See [BUILD](https://github.com/openMVG/openMVG/raw/master/BUILD) text file

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

or cite it as:

```
  @misc{openMVG,
    author = "Pierre Moulon and Pascal Monasse and Renaud Marlet and Others",
     title = "OpenMVG. An Open Multiple View Geometry library.",
    howpublished = "\url{https://github.com/openMVG/openMVG}",
  }
```

