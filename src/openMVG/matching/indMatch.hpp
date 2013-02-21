
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IND_MATCH_H
#define OPENMVG_MATCHING_IND_MATCH_H

#include <iostream>
#include <set>
#include <vector>

namespace openMVG {
namespace matching {

/// Structure in order to save paiwise indexed references.
/// A sort operator exist in order to remove duplicates of IndMatch series.
struct IndMatch
{
  IndMatch(size_t i = 0, size_t j = 0)  {
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
    if (m1._i < m2._i)
      return m1._j < m2._j;
    else
      if (m1._i > m2._i)
        return m1._j < m2._j;
    return m1._i < m2._i;
  }

  /// Remove duplicates (same _i or _j that appears multiple times)
  static size_t getDeduplicated(std::vector<IndMatch> & vec_match){

    size_t sizeBefore = vec_match.size();
    std::set<IndMatch> set_deduplicated( vec_match.begin(), vec_match.end());
    vec_match.assign(set_deduplicated.begin(), set_deduplicated.end());
    return sizeBefore != vec_match.size();
  }

  size_t _i, _j;  // Left, right index
};

static std::ostream& operator<<(std::ostream & out, const IndMatch & obj) {
  return out << obj._i << " " << obj._j << std::endl;
}

static inline std::istream& operator>>(std::istream & in, IndMatch & obj) {
  return in >> obj._i >> obj._j;
}

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_IND_MATCH_H
