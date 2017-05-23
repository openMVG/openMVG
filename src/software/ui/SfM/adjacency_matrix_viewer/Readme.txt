// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

============
Introduction
============

This binary allow to navigate and view the pairwise matches from an OpenMVG
 matching file.
Press file Open, and then choose a sfm_data file and then a matches.X.Y file.
Then the matching matrix is displayed, just click on a pair and see the matches.

===========
Compilation
===========

This sample requires QT5 and OpenMVG

MAC instruction:

install qt5:
 $ sudo brew install qt5

In the cmake command line you will certainly need to tell where is located QT:
 & cmake ... -DCMAKE_PREFIX_PATH=/usr/local/Cellar/qt5/5.7.0/
