// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/fast/fast_detector.hpp"

#include "openMVG/image/image_container.hpp"
#include "third_party/fast/fast.h"

//
// Bibliography
//
// [1] Machine learning for high-speed corner detection
// Authors: Edward Rosten and Tom Drummond
// Conference: ECCV 2006, European Conference on Computer Vision.
//

namespace openMVG
{
namespace features
{

FastCornerDetector::FastCornerDetector
(
  int size,
  int threshold
)
:threshold_(threshold), size_(size)
{
}

void FastCornerDetector::detect
(
  const image::Image<unsigned char> & ima,
  std::vector<PointFeature> & regions
)
{
  using FastDetectorCall =
    xy* (*) (const unsigned char *, int, int, int, int, int *);

  FastDetectorCall detector = nullptr;
  if (size_ ==  9) detector =  fast9_detect_nonmax;
  if (size_ == 10) detector = fast10_detect_nonmax;
  if (size_ == 11) detector = fast11_detect_nonmax;
  if (size_ == 12) detector = fast12_detect_nonmax;
  if (!detector)
  {
    std::cout << "Invalid size for FAST detector: " << size_ << std::endl;
    return;
  }

  int num_corners = 0;
  xy* detections = detector(ima.data(),
     ima.Width(), ima.Height(), ima.Width(),
     threshold_, &num_corners);
  regions.clear();
  regions.reserve(num_corners);
  for (int i = 0; i < num_corners; ++i)
  {
    regions.emplace_back(detections[i].x, detections[i].y);
  }
  free( detections );
}

} // namespace features
} // namespace openMVG
