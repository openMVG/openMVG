// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_FAST_FAST_DETECTOR_HPP
#define OPENMVG_FEATURES_FAST_FAST_DETECTOR_HPP

#include <vector>

#include "openMVG/features/feature.hpp"
#include "openMVG/image/image_container.hpp"

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

class FastCornerDetector
{
  int threshold_;  // Threshold called barrier in the Fast paper (cf. [1]).
  int size_;  // In pixels {9,10,11,12}.

public:

  /**
   * Creates a detector that uses the FAST detection algorithm.
   *
   * \param size      The size of features to detect in pixels {9,10,11,12}.
   * \param threshold Threshold for detecting features (barrier). See the FAST
   *                  paper for details [1].
  **/
  FastCornerDetector
  (
    int size = 9,
    int threshold = 30
  );

  void detect
  (
    const image::Image<unsigned char> & ima,
    std::vector<PointFeature> & regions
  );
};

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_FAST_FAST_DETECTOR_HPP
