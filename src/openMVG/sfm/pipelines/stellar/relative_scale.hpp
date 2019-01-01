// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_STELLAR_RELATIVE_SCALE_HPP
#define OPENMVG_SFM_STELLAR_RELATIVE_SCALE_HPP

#include <openMVG/types.hpp>

#include <array>
#include <vector>

namespace openMVG{
namespace sfm{

/// Structure used to store the depth ratio factor between two pair (a triplet)
struct Relative_Scale
{
  std::array<Pair,2> pairs;
  double ratio;

  Relative_Scale
  (
    const Pair & pair_a = {0,0},
    const Pair & pair_b = {0,0},
    const double ratio = 1.0
  ): ratio(ratio), pairs({{pair_a, pair_b}})
  {
  }

  static Pair_Set Get_pairs
  (
    const std::vector<Relative_Scale> & relative_scales
  )
  {
    Pair_Set pair_set;
    auto insert_pairs = [&](const Relative_Scale & rhs) {
      pair_set.insert(rhs.pairs[0]);
      pair_set.insert(rhs.pairs[1]); };
    std::for_each(relative_scales.cbegin(), relative_scales.cend(), insert_pairs);
    return std::move(pair_set);
  }
};


} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_STELLAR_RELATIVE_SCALE_HPP
