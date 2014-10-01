// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/image_features/image_features.hpp"
#include "third_party/fast/fast.h"

//
// Bibliography
//
// [1] Machine learning for high-speed corner detection
// Authors: Edward Rosten and Tom Drummond
// Conference: ECCV 2006, European Conference on Computer Vision.
//

namespace openMVG{
namespace image_features{

typedef xy* (*FastDetectorCall)(
  const unsigned char *, int, int, int, int, int *);

class FastCornerDetector
{
  int _threshold;  // Threshold called barrier in Fast paper (cf. [1]).
  int _size;  // In pixels {9,10,11,12}.

public:

  /**
	   * Creates a detector that uses the FAST detection algorithm.
	   *
	   * \param size      The size of features to detect in pixels {9,10,11,12}.
	   * \param threshold Threshold for detecting features (barrier). See the FAST
	   *                  paper for details [1].
  **/
  FastCornerDetector(int size = 9, int threshold = 30):
  _size(size), _threshold(threshold)
  {
  }

  void detect(
    const Image<unsigned char> & ima,
    std::vector<SIOPointFeature> & regions)
  {
    FastDetectorCall detector = NULL;
	  if (_size ==  9) detector =  fast9_detect_nonmax;
	  if (_size == 10) detector = fast10_detect_nonmax;
	  if (_size == 11) detector = fast11_detect_nonmax;
	  if (_size == 12) detector = fast12_detect_nonmax;
	  if (!detector) {
	    std::cout << "Invalid size for FAST detector: " << _size << std::endl;
      return;
	  }

    int num_corners = 0;
    xy* detections = detector(ima.data(),
       ima.Width(), ima.Height(), ima.Width(),
       _threshold, &num_corners);
    for (int i = 0; i < num_corners; ++i) {
     regions.push_back(SIOPointFeature(detections[i].x, detections[i].y, 3.0f));
    }
	  free( detections );
  }
};

} // image_features
} // namespace openMVG

