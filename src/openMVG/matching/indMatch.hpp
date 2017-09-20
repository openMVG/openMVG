// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IND_MATCH_HPP
#define OPENMVG_MATCHING_IND_MATCH_HPP

#include <algorithm>
#include <iostream>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "openMVG/stl/hash.hpp"

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

  // Used for lexicographical ordering.
  friend bool operator<(const IndMatch& m1, const IndMatch& m2)  {
    if (m1.i_ < m2.i_)
      return true;

    if (m1.i_ > m2.i_)
      return false;

    return m1.j_ < m2.j_;
  }

  struct IndMatchHasher
  {
    std::size_t operator()(const IndMatch& rhs) const
    {
      std::size_t hash_value(0);
      stl::hash_combine(hash_value, rhs.i_);
      stl::hash_combine(hash_value, rhs.j_);
      return hash_value;
    }
  };

  /// Remove duplicates ((i_, j_) that appears multiple times)
  static bool DeduplicateLeftAndRight(std::vector<IndMatch> & vec_match)  {
    const size_t size_before = vec_match.size();
    const std::unordered_set<IndMatch, IndMatchHasher> deduplicated(
      vec_match.cbegin(), vec_match.cend());
    vec_match.assign(deduplicated.cbegin(), deduplicated.cend());
    std::sort(vec_match.begin(), vec_match.end());
    return size_before != vec_match.size();
  }

  static bool DeduplicateLeftOrRight(std::vector<IndMatch> & matches)
  {
    const std::size_t size_before = matches.size();
    // Hash the left and right indices in order to detect duplicate
    std::unordered_map<IndexT, bool> left, right;
    for (const auto & match : matches)
    {
      if (left.count(match.i_) == 0)
        left[match.i_] = false;    // not yet a duplicate (one occurence so far)
      else
        left.at(match.i_) = true; // mark as duplicate

      if (right.count(match.j_) == 0)
        right[match.j_] = false; // not yet a duplicate (one occurence so far)
      else
        right.at(match.j_) = true; // mark as duplicate
    }

    // Test if some index are repeated
    if (left.size() != matches.size() || right.size() != matches.size())
    {
      std::vector<IndMatch> cleaned_matches;
      cleaned_matches.reserve(std::min(left.size(), right.size()));
      for (const auto & match : matches)
      {
        // If the right and left indexes are not repeated, we keep the match
        if (!left.at(match.i_) && !right.at(match.j_))
          cleaned_matches.push_back(match);
      }
      matches = std::move(cleaned_matches);
    }
    // Else there is no change to do

    return size_before != matches.size();
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
  public std::map< Pair, IndMatches >
{
  void insert(std::pair<Pair, IndMatches> && pairWiseMatches) override
  {
    std::map<Pair, IndMatches>::insert(
      std::forward<std::pair<Pair, IndMatches>>(pairWiseMatches));
  }

  // Serialization
  template <class Archive>
  void serialize( Archive & ar )  {
    ar(static_cast<std::map< Pair, IndMatches >&>(*this));
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
