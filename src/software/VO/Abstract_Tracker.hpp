// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace openMVG { namespace Image { template<typename T> class Image; } }

namespace openMVG  {
namespace VO  {

struct Abstract_Tracker
{
  Abstract_Tracker() = default;

  // Try to track the pt_to_track features
  virtual bool track
  (
    const image::Image<unsigned char> & ima,
    const std::vector<features::PointFeature> & pt_to_track,
    std::vector<features::PointFeature> & pt_tracked,
    std::vector<bool> & status
  ) = 0;

  // suggest new feature point for tracking (count point are kept)
  virtual bool detect
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & pt_to_track,
    const size_t count
  ) const = 0;
};

} // namespace VO
} // namespace openMVG
