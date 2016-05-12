// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_TYPES_H_
#define OPENMVG_TYPES_H_

#include <Eigen/Core>

#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <vector>

#ifdef __clang__

#include <utility>
#include "openMVG/stl/hash.hpp"
#include <unordered_map>

#define OPENMVG_STD_UNORDERED_MAP 1

namespace std {
  template<typename T1, typename T2>
  struct hash<std::pair<T1, T2>> {
    std::size_t operator()(std::pair<T1, T2> const &p) const {
      std::size_t seed1(0);
      stl::hash_combine(seed1, p.first);
      stl::hash_combine(seed1, p.second);

      std::size_t seed2(0);
      stl::hash_combine(seed2, p.second);
      stl::hash_combine(seed2, p.first);

      return std::min(seed1, seed2);
    }
  };
}

#endif // __clang__

/**
* @brief Main namespace of openMVG API
*/
namespace openMVG
{

/// Portable type used to store an index
typedef uint32_t IndexT;

/// Portable value used to save an undefined index value
static const IndexT UndefinedIndexT = std::numeric_limits<IndexT>::max();

/// Standard Pair of IndexT
typedef std::pair<IndexT, IndexT> Pair;

/// Set of Pair
typedef std::set<Pair> Pair_Set;

/// Vector of Pair
typedef std::vector<Pair> Pair_Vec;

#if defined OPENMVG_STD_UNORDERED_MAP

/**
* @brief Standard Hash_Map class
* @tparam K type of the keys
* @tparam V type of the values
*/
template<typename Key, typename Value>
struct Hash_Map : std::unordered_map<Key, Value> {};

#else

/**
* @brief Standard Hash_Map class
* @tparam K type of the keys
* @tparam V type of the values
*/
template<typename K, typename V>
struct Hash_Map : std::map<K, V, std::less<K>,
  Eigen::aligned_allocator<std::pair<const K, V> > > {};

#endif // OPENMVG_STD_UNORDERED_MAP

} // namespace openMVG

#endif  // OPENMVG_TYPES_H_
