// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IND_MATCH_HPP
#define OPENMVG_MATCHING_IND_MATCH_HPP

#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "openMVG/types.hpp"

namespace openMVG {
namespace matching {

/// Structure in order to save pairwise indexed references.
/// A sort operator exist in order to remove duplicates of IndMatch series.
struct IndMatch
{
  IndMatch(IndexT i = 0, IndexT j = 0) : i_(i), j_(j)  {}

  friend bool operator==(const IndMatch& m1, const IndMatch& m2)  {
    return (m1.i_ == m2.i_ && m1.j_ == m2.j_);
  }

  friend bool operator!=(const IndMatch& m1, const IndMatch& m2)  {
    return !(m1 == m2);
  }

  // Lexicographical ordering of matches. Used to remove duplicates
  friend bool operator<(const IndMatch& m1, const IndMatch& m2)  {
    return (m1.i_ < m2.i_ || (m1.i_ == m2.i_ && m1.j_ < m2.j_));
  }

  /// Remove duplicates ((i_, j_) that appears multiple times)
  static bool getDeduplicated(std::vector<IndMatch> & vec_match)  {

    const size_t sizeBefore = vec_match.size();
    const std::set<IndMatch> set_deduplicated( vec_match.begin(), vec_match.end());
    vec_match.assign(set_deduplicated.begin(), set_deduplicated.end());
    return sizeBefore != vec_match.size();
  }

  // Serialization
  template <class Archive>
  void serialize( Archive & ar );

  IndexT i_, j_;  // Left, right index
};

inline std::ostream& operator<<(std::ostream & out, const IndMatch & obj) {
  return out << obj.i_ << " " << obj.j_;
}

inline std::istream& operator>>(std::istream & in, IndMatch & obj) {
  return in >> obj.i_ >> obj.j_;
}

using IndMatches = std::vector<matching::IndMatch>;

/// Pairwise matches (indexed matches for a pair <I,J>)
/// The interface used to store corresponding point indexes per images pairs
class PairWiseMatchesContainer
{
public:
  virtual ~PairWiseMatchesContainer() = default;
  virtual void insert(std::pair<Pair, IndMatches>&& pairWiseMatches) = 0;
};

//--
/// Pairwise matches (indexed matches for a pair <I,J>)
/// A structure used to store corresponding point indexes per images pairs
struct PairWiseMatches :
  public PairWiseMatchesContainer,
  public std::map<Pair, IndMatches>
{
  void insert(std::pair<Pair, IndMatches> && pairWiseMatches)override
  {
    std::map<Pair, IndMatches>::insert(
      std::forward<std::pair<Pair, IndMatches>>(pairWiseMatches));
  }

  // Serialization
  template <class Archive>
  void serialize( Archive & ar )  {
    ar(static_cast<std::map<Pair, IndMatches>&>(*this));
  }
};

inline Pair_Set getPairs(const PairWiseMatches & matches)
{
  Pair_Set pairs;
  for ( const auto & cur_pair : matches )
    pairs.insert(cur_pair.first);
  return pairs;
}

}  // namespace matching
}  // namespace openMVG


#endif // OPENMVG_MATCHING_IND_MATCH_HPP
