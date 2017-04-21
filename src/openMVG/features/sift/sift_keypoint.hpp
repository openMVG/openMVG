// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_SIFT_SIFT_KEYPOINT_HPP
#define OPENMVG_FEATURES_SIFT_SIFT_KEYPOINT_HPP

namespace openMVG{
namespace features{
namespace sift{

struct Keypoint
{
  float x; // x coordinate  (subpixel)
  float y; // y coordinate  (subpixel)
  float sigma; // level of blur (it includes the assumed image blur)
  float theta; // orientation

  // discrete coordinates
  int o; // octave level
  int s; // scale level (slice index)
  int i; // x (pixel)
  int j; // y (pixel)

  float val;      // normalized operator value (independent of the scale-space sampling)
  float edgeResp; // edge response

  // Descriptor
  openMVG::Vecf descr;
};

} // namespace sift
} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_SIFT_SIFT_KEYPOINT_HPP
