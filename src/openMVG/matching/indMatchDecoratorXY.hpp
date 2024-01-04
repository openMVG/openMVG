// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IND_MATCH_DECORATOR_XY_HPP
#define OPENMVG_MATCHING_IND_MATCH_DECORATOR_XY_HPP

#include <algorithm>
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
    explicit IndMatchDecoratorStruct
    (
      T xa, T ya,
      T xb, T yb,
      const IndMatch & ind
    )
    {
      x1 = xa; y1 = ya;
      x2 = xb; y2 = yb;
      index = ind;
    }

    /// Lexicographical ordering of matches. Used to remove duplicates.
    friend bool operator<
    (
      const IndMatchDecoratorStruct& m1,
      const IndMatchDecoratorStruct& m2
    )
    {
      if (m1 == m2) return false;

      if (m1.x1 < m2.x1)
        return m1.y1 < m2.y1;
      else
        if (m1.x1 > m2.x1)
          return m1.y1 < m2.y1;
      return m1.x1 < m2.x1;
    }

    /// Comparison Operator
    friend bool operator==
    (
      const IndMatchDecoratorStruct& m1,
      const IndMatchDecoratorStruct& m2
    )
    {
      return (m1.x1==m2.x1 && m1.y1==m2.y1 &&
        m1.x2==m2.x2 && m1.y2==m2.y2);
    }

    std::ostream & operator<<(std::ostream &os) const
    {
      return os << x1 << " " << y1 << " " << x2 << " " << y2 << "\n";
    }

    T x1,y1, x2,y2;
    IndMatch index;
  };
public:

  template<typename Feature_t>
  explicit IndMatchDecorator
  (
    const std::vector<IndMatch> & vec_matches,
    const std::vector<Feature_t> & leftFeat,
    const std::vector<Feature_t> & rightFeat
  )
  :vec_matches_(vec_matches)
  {
    for ( const auto & cur_vec_match : vec_matches )
    {
      const size_t I = cur_vec_match.i_;
      const size_t J = cur_vec_match.j_;
      vec_decoredMatches_.emplace_back(
        leftFeat[I].x(),leftFeat[I].y(),
        rightFeat[J].x(), rightFeat[J].y(), cur_vec_match);
    }
  }

  explicit IndMatchDecorator
  (
    const std::vector<IndMatch> & vec_matches,
    const Mat & leftFeat,
    const Mat & rightFeat
  )
  :vec_matches_(vec_matches)
  {
    for ( const auto & cur_vec_match : vec_matches )
    {
      const size_t I = cur_vec_match.i_;
      const size_t J = cur_vec_match.j_;
      vec_decoredMatches_.emplace_back(
        leftFeat.col(I)(0),leftFeat.col(I)(1),
        rightFeat.col(J)(0), rightFeat.col(J)(1), cur_vec_match);
    }
  }

  /// Remove duplicates (same (x1,y1) coords that appears multiple times)
  size_t getDeduplicated(std::vector<IndMatch> & vec_matches)
  {
    const size_t sizeBefore = vec_decoredMatches_.size();
    const std::set<IndMatchDecoratorStruct> set_deduplicated(
      vec_decoredMatches_.cbegin(),vec_decoredMatches_.cend());
    vec_decoredMatches_.assign(set_deduplicated.cbegin(), set_deduplicated.cend());

    vec_matches.resize(vec_decoredMatches_.size());
    for (size_t i = 0; i < vec_decoredMatches_.size(); ++i)  {
      const IndMatch & idxM = vec_decoredMatches_[i].index;
      vec_matches[i] = idxM;
    }

    return sizeBefore != vec_matches.size();
  }

  /**
  * Save the corresponding matches to file.
  * \param nameFile   The file where matches will be saved.
  * \param vec_match  The matches that we want to save.
  * \return bool True if everything was ok, otherwise false.
  */
  bool saveMatch(const char* nameFile) const
  {
    std::ofstream f(nameFile);
    if (f.is_open() ) {
      std::copy(vec_decoredMatches_.cbegin(), vec_decoredMatches_.cend(),
        std::ostream_iterator<IndMatchDecoratorStruct>(f, ""));
    }
    return f.is_open();
  }

private:
  std::vector<IndMatch> vec_matches_;
  std::vector<IndMatchDecoratorStruct> vec_decoredMatches_;
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_IND_MATCH_DECORATOR_XY_HPP
