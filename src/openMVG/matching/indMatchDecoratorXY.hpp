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
#include <set>
#include <vector>

#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"

namespace openMVG {
namespace matching {

/// IndMatch decorator.
/// Use sorting over x,y coordinates.
template<class T = float>
class IndMatchDecorator
{
  struct IndMatchDecoratorStruct
  {
    IndMatchDecoratorStruct
    (
      T xa, T ya,
      T xb, T yb,
      const matching::IndMatch & ind
    ): x1(xa), y1(ya), x2(xb), y2(yb), index(ind)
    {
    }

    /// Lexicographical ordering of matches. Used to remove duplicates.
    friend bool operator<
    (
      const IndMatchDecoratorStruct& m1,
      const IndMatchDecoratorStruct& m2
    )
    {
      if (m1.x1 < m2.x1)
        return true;

      if (m1.x1 > m2.x1)
         return false;

      return m1.y1 < m2.y1;
    }

    /// Comparison Operator (use only the point position for de-duplication)
    friend bool operator==
    (
      const IndMatchDecoratorStruct& m1,
      const IndMatchDecoratorStruct& m2
    )
    {
      // Only use the point position
      return (m1.x1==m2.x1 && m1.y1==m2.y1 &&
        m1.x2==m2.x2 && m1.y2==m2.y2);
    }

    T x1, y1, x2, y2;
    matching::IndMatch index;
  };
public:

  IndMatchDecorator
  (
    const matching::IndMatches & matches,
    const std::vector<features::SIOPointFeature> & leftFeat,
    const std::vector<features::SIOPointFeature> & rightFeat
  )
  :matches_(matches)
  {
    for ( const auto & cur_vec_match : matches_ )
    {
      const auto I = cur_vec_match.i_;
      const auto J = cur_vec_match.j_;
      const IndMatchDecoratorStruct current_match(
        leftFeat[I].x(),leftFeat[I].y(),
        rightFeat[J].x(), rightFeat[J].y(),
        cur_vec_match);
      deduplicated_matches_.insert(current_match);
    }
  }

  IndMatchDecorator
  (
    const matching::IndMatches & matches,
    const std::vector<features::PointFeature> & leftFeat,
    const std::vector<features::PointFeature> & rightFeat
  )
  :matches_(matches)
  {
    for ( const auto & cur_vec_match : matches_ )
    {
      const auto I = cur_vec_match.i_;
      const auto J = cur_vec_match.j_;
      const IndMatchDecoratorStruct current_match(
        leftFeat[I].x(),leftFeat[I].y(),
        rightFeat[J].x(), rightFeat[J].y(),
        cur_vec_match);
      deduplicated_matches_.insert(current_match);
    }
  }

  IndMatchDecorator
  (
    const matching::IndMatches & matches,
    const Mat2X & leftFeat,
    const Mat2X & rightFeat
  )
  :matches_(matches)
  {
    for ( const auto & cur_vec_match : matches_ )
    {
      const auto I = cur_vec_match.i_;
      const auto J = cur_vec_match.j_;
      const IndMatchDecoratorStruct current_match(
        leftFeat.col(I).x(),leftFeat.col(I).y(),
        rightFeat.col(J).x(), rightFeat.col(J).y(),
        cur_vec_match);
      deduplicated_matches_.insert(current_match);
    }
  }

  /// Remove duplicates (same (x,y) coords that appears multiple times)
  bool getDeduplicated(matching::IndMatches & matches) const
  {
    const auto size_before = matches_.size();

    matches.clear();
    matches.reserve(deduplicated_matches_.size());
    for (const auto & idx : deduplicated_matches_)  {
      matches.emplace_back(idx.index);
    }

    return size_before != matches.size();
  }

  friend std::ostream& operator<<(std::ostream& os, const IndMatchDecoratorStruct & m)
  {
    return os << m.x1 << ' ' << m.y1 << ' ' << m.x2 << ' ' << m.y2 << '\n';
  }

private:
  matching::IndMatches matches_;
  std::set<IndMatchDecoratorStruct> deduplicated_matches_;
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_IND_MATCH_DECORATOR_XY_HPP
