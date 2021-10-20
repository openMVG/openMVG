// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_REGIONS_SCALE_SORT_HPP
#define OPENMVG_FEATURES_REGIONS_SCALE_SORT_HPP

#include <numeric>      // std::iota
#include <algorithm>    // std::sort

namespace openMVG {
namespace features {

/// Get the reordering indexes of the feature according the feature scale
/// reordered by increasing scale.
template <typename T>
std::vector<size_t> sort_indexes_by_scale(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing scale values of the features in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1].scale() < v[i2].scale();});

  return idx;
}

/// Sort the features and descriptor according the scale of the feature
/// A `keep_count` filter can be used to select only the feature with the largest scale
template
<
typename FeatT,
typename DescsT
>
bool SortAndSelectByRegionScale
(
  std::vector<FeatT> & feats,
  DescsT & descs,
  int keep_count = -1
)
{
  return false;
}

/// SortAndSelectByRegionScale specialization for the SIOPointFeature type.
template
<
typename FeatT,
typename DescsT
>
bool SortAndSelectByRegionScale
(
  std::vector<SIOPointFeature> & feats,
  DescsT & descs,
  int keep_count = -1
)
{
  // Reorder features and descriptors
  // a. get the re-ordering sequence
  // b. apply the re-ordering sequence
  const std::vector<size_t> re_ordering = sort_indexes_by_scale(feats);
  auto feats_reordered = feats;
  auto descs_reordered = descs;
  for (size_t i = 0; i < re_ordering.size(); ++i)
  {
    feats_reordered[i] = feats[re_ordering[i]];
    descs_reordered[i] = descs[re_ordering[i]];
  }

  // Keep only the regions with the largest scale:
  {
    const int region_to_keep =
      (keep_count == -1) ?
        feats_reordered.size() :
        std::min(keep_count, static_cast<int>(feats_reordered.size()));

    auto start_feat_it = std::prev(feats_reordered.end(), region_to_keep);
    feats.assign(start_feat_it, feats_reordered.end());

    auto start_desc_it = std::prev(descs_reordered.end(), region_to_keep);
    descs.assign(start_desc_it, descs_reordered.end());

  }
  return true;
}

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_REGIONS_SCALE_SORT_HPP
