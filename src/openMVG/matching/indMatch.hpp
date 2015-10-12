
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IND_MATCH_H
#define OPENMVG_MATCHING_IND_MATCH_H

#include "openMVG/types.hpp"

#include <cereal/cereal.hpp> // Serialization

#include <iostream>
#include <set>
#include <map>
#include <vector>

namespace openMVG {
namespace matching {

/// Structure in order to save pairwise indexed references.
/// A sort operator exist in order to remove duplicates of IndMatch series.
struct IndMatch
{
  IndMatch(IndexT i = 0, IndexT j = 0)  {
    _i = i;
    _j = j;
  }

  friend bool operator==(const IndMatch& m1, const IndMatch& m2)  {
    return (m1._i == m2._i && m1._j == m2._j);
  }

  friend bool operator!=(const IndMatch& m1, const IndMatch& m2)  {
    return !(m1 == m2);
  }

  // Lexicographical ordering of matches. Used to remove duplicates.
  friend bool operator<(const IndMatch& m1, const IndMatch& m2) {
    return (m1._i < m2._i || (m1._i == m2._i && m1._j < m2._j));
  }

  /// Remove duplicates ((_i, _j) that appears multiple times)
  static bool getDeduplicated(std::vector<IndMatch> & vec_match)  {

    const size_t sizeBefore = vec_match.size();
    std::set<IndMatch> set_deduplicated( vec_match.begin(), vec_match.end());
    vec_match.assign(set_deduplicated.begin(), set_deduplicated.end());
    return sizeBefore != vec_match.size();
  }

  // Serialization
  template <class Archive>
  void serialize( Archive & ar )  {
    ar(_i, _j);
  }

  IndexT _i, _j;  // Left, right index
};

static std::ostream& operator<<(std::ostream & out, const IndMatch & obj) {
  return out << obj._i << " " << obj._j;
}

static inline std::istream& operator>>(std::istream & in, IndMatch & obj) {
  return in >> obj._i >> obj._j;
}

typedef std::vector<matching::IndMatch> IndMatches;
//--
/// Pairwise matches (indexed matches for a pair <I,J>)
/// The structure used to store corresponding point indexes per images pairs
typedef std::map< Pair, IndMatches > PairWiseMatches;

static Pair_Set getPairs(const PairWiseMatches & matches)
{
  Pair_Set pairs;
  for(PairWiseMatches::const_iterator it = matches.begin(); it != matches.end(); ++it)
    pairs.insert(it->first);
  return pairs;
}

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_IND_MATCH_H
