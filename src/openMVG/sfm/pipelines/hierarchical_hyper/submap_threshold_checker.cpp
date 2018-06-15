// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "submap_threshold_checker.hpp"

namespace openMVG{
namespace sfm{

SubmapTracksThresholdChecker::SubmapTracksThresholdChecker(int tracks_threshold)
  : SubmapThresholdChecker(), tracks_threshold_(tracks_threshold)
{}


bool SubmapTracksThresholdChecker::operator()(const HsfmSubmap & smap) const
{
  return ((!smap.is_parent)
          && smap.sfm_data.GetViews().size() > VIEWS_PER_SUBMAP_LOWER_BOUND
          && (smap.track_ids.size() > tracks_threshold_));
}

bool SubmapTracksThresholdChecker::operator()(const std::pair<IndexT, HsfmSubmap> & smap_with_id) const
{
  return operator ()(smap_with_id.second);
}

SubmapViewThresholdChecker::SubmapViewThresholdChecker(int views_threshold)
  : SubmapThresholdChecker(), views_threshold_(views_threshold)
{}

bool SubmapViewThresholdChecker::operator()(const HsfmSubmap & smap) const
{
  return ((!smap.is_parent)
          && smap.sfm_data.GetViews().size() > std::max(VIEWS_PER_SUBMAP_LOWER_BOUND, views_threshold_));
}

bool SubmapViewThresholdChecker::operator()(const std::pair<IndexT, HsfmSubmap> & smap_with_id) const
{
  return operator ()(smap_with_id.second);
}

} // namespace sfm
} // namespace openMVG
