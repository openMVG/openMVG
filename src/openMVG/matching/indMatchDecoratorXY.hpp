// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IND_MATCH_DECORATOR_XY_HPP
#define OPENMVG_MATCHING_IND_MATCH_DECORATOR_XY_HPP

#include <algorithm>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <vector>

#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/stl/hash.hpp"

namespace openMVG {
namespace matching {


/// IndMatch decorator design to use the corresponding feature position for
///  deduplication purpose.
template<class T = float>
struct IndMatchPositionDecorator
{
  IndMatchPositionDecorator
  (
    T xa, T ya,
    T xb, T yb,
    const matching::IndMatch & ind
  ): x1(xa), y1(ya), x2(xb), y2(yb), match(ind)
  {
  }

  /// Comparison Operator (use only the point position for de-duplication)
  friend bool operator==
  (
    const IndMatchPositionDecorator& m1,
    const IndMatchPositionDecorator& m2
  )
  {
    // Only use the point position
    return (m1.x1 == m2.x1
         && m1.y1 == m2.y1
         && m1.x2 == m2.x2
         && m1.y2 == m2.y2);
  }

  T x1, y1, x2, y2;
  matching::IndMatch match;

  // Hash the match according the features positions.
  // Match index is ignored
  struct Hasher
  {
    std::size_t operator()(const IndMatchPositionDecorator & rhs) const
    {
      std::size_t hash_value(0);
      stl::hash_combine(hash_value, rhs.x1);
      stl::hash_combine(hash_value, rhs.x2);
      stl::hash_combine(hash_value, rhs.y1);
      stl::hash_combine(hash_value, rhs.y2);
      return hash_value;
    }
  };
};

/// Duplicates removal based on feature position.
/// Duplicate rule: if the left and the right used feature positions are repeated.
/// Duplicate: {{x1,y1}, {x2,y2}} == {{x1,y1}, {x2,y2}}
static bool DeduplicateLeftAndRight
(
  const std::vector<features::PointFeature> & left_features,
  const std::vector<features::PointFeature> & righ_features,
  matching::IndMatches & matches
)
{
  const std::size_t size_before = matches.size();
  std::unordered_set<IndMatchPositionDecorator<float>,
                     IndMatchPositionDecorator<float>::Hasher> deduplicated;
  for ( const auto & match : matches)
  {
   deduplicated.insert(
     {left_features.at(match.i_).x(),left_features.at(match.i_).y(),
      righ_features.at(match.j_).x(), righ_features.at(match.j_).y(),
      match
     });
  }
  matches.clear();
  matches.reserve(deduplicated.size());
  for (const auto & kept_match : deduplicated)
  {
    matches.push_back(kept_match.match);
  }

  return size_before != matches.size();
}

/// Duplicates removal based on feature position.
/// Duplicate rule: if the left or the right used feature position is repeated.
/// Duplicate: {{x1,y1}, {x2,y2}} == {{x1,y1}, {x2,y2}}
/// Duplicate: {{x1,y1}, {x2,y2}} == {{x1,y1}, {x3,y3}}
/// Duplicate: {{x1,y1}, {x2,y2}} == {{x3,y3}, {x2,y2}}
/// Duplicate: {{x1,y1}, {x2,y2}} == {{x1,y1}, {x2,y3}} // {x1,y1} , {x1,y1}
/// Not a duplicate: {{x1,y1}, {x2,y2}} == {{x1,y3}, {x2,y4}} // Y coords are different
static bool DeduplicateLeftOrRight
(
  const std::vector<features::PointFeature> & left_features,
  const std::vector<features::PointFeature> & right_features,
  matching::IndMatches & matches
)
{
  using FloatPair = std::pair<float, float>;
  struct PairHasher
  {
    std::size_t operator()(const FloatPair & rhs) const
    {
      std::size_t hash_value(0);
      stl::hash_combine(hash_value, rhs.first);
      stl::hash_combine(hash_value, rhs.second);
      return hash_value;
    }
  };

  const std::size_t size_before = matches.size();
  // Hash the left and right feature position in order to detect duplicate
  std::unordered_map<FloatPair, bool, PairHasher> left, right;
  for (const auto & match : matches)
  {
    const auto & left_feature = left_features.at(match.i_);
    const auto & right_feature = right_features.at(match.j_);
    if (left.count({left_feature.x(), left_feature.y()}) == 0)
      left[{left_feature.x(), left_feature.y()}] = false;    // not yet a duplicate (one occurence so far)
    else
      left.at({left_feature.x(), left_feature.y()}) = true; // mark as duplicate

    if (right.count({right_feature.x(), right_feature.y()}) == 0)
      right[{right_feature.x(), right_feature.y()}] = false; // not yet a duplicate (one occurence so far)
    else
      right.at({right_feature.x(), right_feature.y()}) = true; // mark as duplicate
  }
  // Test if some feature position are repeated
  if (left.size() != matches.size() || right.size() != matches.size())
  {
    std::vector<IndMatch> cleaned_matches;
    cleaned_matches.reserve(std::min(left.size(), right.size()));
    for (const auto & match : matches)
    {
      const auto & left_feature = left_features.at(match.i_);
      const auto & right_feature = right_features.at(match.j_);
      // If the right and left feature positions are not repeated, we keep the match
      if (!left.at({left_feature.x(), left_feature.y()}) &&
          !right.at({right_feature.x(), right_feature.y()}))
        cleaned_matches.push_back(match);
    }
    matches = std::move(cleaned_matches);
  }
  // Else there is no change to do

  return size_before != matches.size();
}

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_IND_MATCH_DECORATOR_XY_HPP
