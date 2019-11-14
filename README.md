OpenMVG (open Multiple View Geometry)
=====================================

![Logo](https://github.com/openMVG/openMVG/raw/master/logo/openMVG_Logo.png)

Introduction
------------

[OpenMVG (Multiple View Geometry)](http://imagine.enpc.fr/~moulonp/openMVG/) "open Multiple View Geometry" is a library for computer-vision scientists and especially targeted to the Multiple View Geometry community. It is designed to provide easy access to the classical problem solvers in Multiple View Geometry and solve them accurately.

The openMVG credo is: "Keep it simple, keep it maintainable". OpenMVG targets readable code that is easy to use and modify by the community.

All the features and modules are unit tested. This test-driven development ensures that the code works as it should and enables more consistent repeatability. Furthermore, it makes it easier for the user to understand and learn the given features.

[![GitHub license](https://img.shields.io/badge/license-MPL2-blue)](https://github.com/openMVG/openMVG/blob/master/LICENSE)

- [Build](#build)
- [Continuous integration](#integration)
- [Authors](#authors)
- [Contact](#contact)
- [Citations](#citations)
- [Acknowledgements](#acknowledgements)

## Build

Please follow this [build tutorial ](https://github.com/openMVG/openMVG/blob/master/BUILD.md) to build and use OpenMVG locally or in a Docker.

## Continuous integration

 - Build : Linux/Mac (GCC/CLANG): [![Build Status](https://travis-ci.org/openMVG/openMVG.png?branch=develop)](https://travis-ci.org/openMVG/openMVG), Windows (VStudio): [![Build status](https://ci.appveyor.com/api/projects/status/3nv6rt41yxqx5v7i?svg=true)](https://ci.appveyor.com/project/pmoulon/openmvg)
 - Code Quality [![Codacy Badge](https://api.codacy.com/project/badge/Grade/e067bc979aef48f5a96818714a5b33b9)](https://www.codacy.com/manual/pmoulon/openMVG?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=openMVG/openMVG&amp;utm_campaign=Badge_Grade) [![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/openMVG/openMVG.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/openMVG/openMVG/context:cpp)  [![CodeFactor](https://www.codefactor.io/repository/github/openmvg/openmvg/badge)](https://www.codefactor.io/repository/github/openmvg/openmvg)
<!-- - Unit test coverage: [![Coverage Status](https://coveralls.io/repos/openMVG/openMVG/badge.png?branch=develop)](https://coveralls.io/r/openMVG/openMVG?branch=develop) -->



## Authors

See [Authors](https://github.com/openMVG/openMVG/raw/master/AUTHORS) text file

## Documentation

See [documentation](http://openmvg.readthedocs.org/en/latest)

## Contact

openmvg-team[AT]googlegroups.com


## Citations

We are recommending citing `OpenMVG` if you are using the whole library or the adequate paper if you use only a submodule `AContrario Ransac, AContrario
SfM, GlobalSfM or Tracks`:
```
@inproceedings{moulon2016openmvg,
  title={Openmvg: Open multiple view geometry},
  author={Moulon, Pierre and Monasse, Pascal and Perrot, Romuald and Marlet, Renaud},
  booktitle={International Workshop on Reproducible Research in Pattern Recognition},
  pages={60--74},
  year={2016},
  organization={Springer}
}
```

[3] Moulon Pierre, Monasse Pascal and Marlet Renaud. ACCV 2012.
[Adaptive Structure from Motion with a contrario model estimation.](http://hal.archives-ouvertes.fr/index.php?halsid=1n2qdqiv2a0l5eq7qpos9us752&view_this_doc=hal-00769266&version=1)
```
@inproceedings{Moulon2012,
  doi = {10.1007/978-3-642-37447-0_20},
  year  = {2012},
  publisher = {Springer Berlin Heidelberg},
  pages = {257--270},
  author = {Pierre Moulon and Pascal Monasse and Renaud Marlet},
  title = {Adaptive Structure from Motion with a~Contrario Model Estimation},
  booktitle = {Proceedings of the Asian Computer Vision Conference (ACCV 2012)}
}
```

[4] Moulon Pierre and Monasse Pascal. CVMP 2012.
[Unordered feature tracking made fast and easy.](http://hal.archives-ouvertes.fr/index.php?halsid=ggdarhl8cv1j6ohq2073eok8q3&view_this_doc=hal-00769267&version=1)
```
@inproceedings{moulon2012unordered,
  title={Unordered feature tracking made fast and easy},
  author={Moulon, Pierre and Monasse, Pascal},
  booktitle={CVMP 2012},
  pages={1},
  year={2012}
}
```

[5] Moisan Lionel, Moulon Pierre and Monasse Pascal. IPOL 2012.
[Automatic Homographic Registration of a Pair of Images, with A Contrario Elimination of Outliers.](http://dx.doi.org/10.5201/ipol.2012.mmm-oh)
```
@article{moisan2012automatic,
  title={Automatic homographic registration of a pair of images, with a contrario elimination of outliers},
  author={Moisan, Lionel and Moulon, Pierre and Monasse, Pascal},
  journal={Image Processing On Line},
  volume={2},
  pages={56--73},
  year={2012}
}
```

[6] Moulon Pierre, Monasse Pascal, and Marlet Renaud. ICCV 2013.
[Global Fusion of Relative Motions for Robust, Accurate and Scalable Structure from Motion.](http://imagine.enpc.fr/~moulonp/publis/iccv2013/index.html)

```
@inproceedings{moulon2013global,
  title={Global fusion of relative motions for robust, accurate and scalable structure from motion},
  author={Moulon, Pierre and Monasse, Pascal and Marlet, Renaud},
  booktitle={Proceedings of the IEEE International Conference on Computer Vision},
  pages={3248--3255},
  year={2013}
}
```

## Acknowledgements

openMVG authors would like to thanks libmv authors for providing an inspiring
base to design openMVG. Authors also would like to thanks Mikros Image [1]
and LIGM-Imagine laboratory [2] for support and authorization to make this
library an opensource project.

[1] [http://www.mikrosimage.eu/](http://www.mikrosimage.eu/)
[2] [http://imagine.enpc.fr/](http://imagine.enpc.fr/)
