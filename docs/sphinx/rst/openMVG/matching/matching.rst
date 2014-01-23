*******************
matching
*******************

Method
=============

A generic interface to perform K-Nearest Neighbor search:

* Brute force,
* Approximate Nearest Neighbor [FLANN]_.

This module works for any size of data, it could be use to match:

* 2 or 3 vector long features (points),
* 128, 64, vector long features (like SIFT, SURF descriptors).

.. code-block:: c++


  // Setup the matcher 
  ArrayMatcherBruteForce<float> matcher;
  // The reference array
  float array[] = {0, 1, 2, 3, 4};
  // Setup the reference array of the matcher 
  matcher.Build(array, 5, 1);

  //--
  // Looking for the nearest neighbor:
  //--
  // Perform a query to look which point is closest to 1.8
  float query[] = {1.8f};
  int nIndice = -1;
  float fDistance = -1.0f;
  matcher.SearchNeighbour(query, &nIndice, &fDistance);

  //  nIndice == 2 ; // index of the found nearest neighbor
  //  fDistance == 0.4; // squared distance

  //--
  // Looking for the K=2 nearest neighbor
  //--
  vector<int> vec_nIndice;
  vector<float> vec_fDistance;
  const int K = 2;
  matcher.SearchNeighbours(query, 1, &vec_nIndice, &vec_fDistance, K);
  
  // vec_nIndice = {2,1};

Metric
=============

Used metric is customizable and enable matching under:

* L2,
* or an user customized distance (L1, ...).

Filtering
=============

Once putatives matches have been found they can be filtered thanks to filters:

* Symmetric distance (Left-Right check)
* "Nearest Neighbor Distance Ratio" distance check can be performed to remove repetitive elements.

  * when K>=2 nearest points are asked by query points.

.. [FLANN] Muja, Marius, and David G. Lowe. "Fast Approximate Nearest Neighbors with Automatic Algorithm Configuration." VISAPP (1). 2009.
