
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
    IndMatchDecoratorStruct(
      T xa, T ya,
      T xb, T yb,
      const IndMatch & ind) {

      x1 = xa; y1 = ya;
      x2 = xb; y2 = yb;
      index = ind;
    }

    /// Lexicographical ordering of matches. Used to remove duplicates.
    friend bool operator<(const IndMatchDecoratorStruct& m1,
      const IndMatchDecoratorStruct& m2)  {

      if(m1.x1 < m2.x1) return true;
      if(m1.x1 > m2.x1) return false;

      if(m1.y1 < m2.y1) return true;
      if(m1.y1 > m2.y1) return false;

      if(m1.x2 < m2.x2) return true;
      if(m1.x2 > m2.x2) return false;

      return (m1.y2 < m2.y2);
    }

    /// Comparison Operator
    friend bool operator==(const IndMatchDecoratorStruct& m1,
      const IndMatchDecoratorStruct& m2)  {

      return (m1.x1==m2.x1 && m1.y1==m2.y1 &&
        m1.x2==m2.x2 && m1.y2==m2.y2);
    }

    T x1,y1, x2,y2;
    IndMatch index;
  };
public:
  //TODO must use x,y,x1,y1, to not depends on SIOPointFeature class
  IndMatchDecorator(const std::vector<IndMatch> & vec_matches,
    const std::vector<SIOPointFeature> & rightFeat,
    const std::vector<SIOPointFeature> & leftFeat)
    :_vec_matches(vec_matches)
  {
    for (size_t i = 0; i < vec_matches.size(); ++i) {
      const size_t I = vec_matches[i]._i;
      const size_t J = vec_matches[i]._j;
      _vecDecoredMatches.push_back(
        IndMatchDecoratorStruct(rightFeat[I].x(),rightFeat[I].y(),
        leftFeat[J].x(), leftFeat[J].y(), vec_matches[i]));
    }
  }

  void getDeduplicated(std::vector<IndMatch> & vec_matches,
    bool bVerbose = false)
  {
    std::sort(_vecDecoredMatches.begin(), _vecDecoredMatches.end());
    typename std::vector<IndMatchDecoratorStruct >::iterator end =
      std::unique(_vecDecoredMatches.begin(), _vecDecoredMatches.end());
    if(end != _vecDecoredMatches.end()) {
      if (bVerbose) {
        std::cout << "Remove " << std::distance(end, _vecDecoredMatches.end())
          << "/" << _vecDecoredMatches.size() << " duplicate matches, "
          << " keeping " << std::distance(_vecDecoredMatches.begin(), end)
          <<std::endl;
      }
      _vecDecoredMatches.erase(end, _vecDecoredMatches.end());
    }

    vec_matches.resize(_vecDecoredMatches.size());
    for (size_t i = 0; i < _vecDecoredMatches.size(); ++i)  {
      const IndMatch & idxM = _vecDecoredMatches[i].index;
      vec_matches[i] = idxM;
    }
  }

  /**
    * Save the corresponding matches to file.
    * \param nameFile   The file where matches will be saved.
    * \param vec_match  The matches that we want to save.
    * \return bool True if everything was ok, otherwise false.
    */
    bool saveMatch(const char* nameFile) const  {
      std::ofstream f(nameFile);
      if( f.is_open() ) {
        std::copy(_vecDecoredMatches.begin(), _vecDecoredMatches.end(),
          std::ostream_iterator<IndMatchDecoratorStruct>(f, ""));
      }
      return f.is_open();
    }

    friend std::ostream& operator<<(std::ostream& os, const IndMatchDecoratorStruct & m)
    {
      return os << m.x1 << " " << m.y1 << " " << m.x2 << " " << m.y2 << "\n";
    }


private :
  std::vector<IndMatch> _vec_matches;
  std::vector<IndMatchDecoratorStruct> _vecDecoredMatches;
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_IND_MATCH_DECORATOR_XY_H
