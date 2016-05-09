
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IND_MATCH_DECORATOR_XY_H
#define OPENMVG_MATCHING_IND_MATCH_DECORATOR_XY_H

#include <iostream>
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/features/features.hpp"

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

    T x1,y1, x2,y2;
    IndMatch index;
  };
public:

  IndMatchDecorator
  (
    const std::vector<IndMatch> & vec_matches,
    const std::vector<features::SIOPointFeature> & leftFeat,
    const std::vector<features::SIOPointFeature> & rightFeat
  )
  :vec_matches_(vec_matches)
  {
    for ( const auto & cur_vec_match : vec_matches )
    {
      const size_t I = cur_vec_match.i_;
      const size_t J = cur_vec_match.j_;
      vecDecoredMatches_.push_back(
        IndMatchDecoratorStruct(leftFeat[I].x(),leftFeat[I].y(),
        rightFeat[J].x(), rightFeat[J].y(), cur_vec_match));
    }
  }

  IndMatchDecorator
  (
    const std::vector<IndMatch> & vec_matches,
    const std::vector<features::PointFeature> & leftFeat,
    const std::vector<features::PointFeature> & rightFeat
  )
  :vec_matches_(vec_matches)
  {
    for ( const auto & cur_vec_match : vec_matches )
    {
      const size_t I = cur_vec_match.i_;
      const size_t J = cur_vec_match.j_;
      vecDecoredMatches_.push_back(
        IndMatchDecoratorStruct(leftFeat[I].x(),leftFeat[I].y(),
        rightFeat[J].x(), rightFeat[J].y(), cur_vec_match));
    }
  }

  IndMatchDecorator
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
      vecDecoredMatches_.push_back(
        IndMatchDecoratorStruct(leftFeat.col(I)(0),leftFeat.col(I)(1),
        rightFeat.col(J)(0), rightFeat.col(J)(1), cur_vec_match));
    }
  }

  /// Remove duplicates (same (x1,y1) coords that appears multiple times)
  size_t getDeduplicated(std::vector<IndMatch> & vec_matches)
  {
    const size_t sizeBefore = vecDecoredMatches_.size();
    std::set<IndMatchDecoratorStruct> set_deduplicated(
      vecDecoredMatches_.begin(),vecDecoredMatches_.end());
    vecDecoredMatches_.assign(set_deduplicated.begin(), set_deduplicated.end());

    vec_matches.resize(vecDecoredMatches_.size());
    for (size_t i = 0; i < vecDecoredMatches_.size(); ++i)  {
      const IndMatch & idxM = vecDecoredMatches_[i].index;
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
      if( f.is_open() ) {
        std::copy(vecDecoredMatches_.begin(), vecDecoredMatches_.end(),
          std::ostream_iterator<IndMatchDecoratorStruct>(f, ""));
      }
      return f.is_open();
    }

    friend std::ostream& operator<<(std::ostream& os, const IndMatchDecoratorStruct & m)
    {
      return os << m.x1 << " " << m.y1 << " " << m.x2 << " " << m.y2 << "\n";
    }

private :
  std::vector<IndMatch> vec_matches_;
  std::vector<IndMatchDecoratorStruct> vecDecoredMatches_;
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_IND_MATCH_DECORATOR_XY_H
