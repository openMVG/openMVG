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

#if defined OPENMVG_STD_UNORDERED_MAP
  #include <unordered_map>
#endif

/**
* @brief Main namespace of openMVG API
*/
namespace openMVG
{

/// Portable type used to store an index
typedef uint32_t IndexT;

/// Portable value used to save an undefined index value
static const IndexT UndefinedIndexT = std::numeric_limits<IndexT>::max();

/// Standard Pair of int
typedef std::pair<IndexT, IndexT> Pair;

/// A set of int-pair
typedef std::set<Pair> Pair_Set;

/// Vector of Pairs
typedef std::vector<Pair> Pair_Vec;

#define OPENMVG_NO_UNORDERED_MAP 1

#if defined OPENMVG_NO_UNORDERED_MAP

/**
* @brief Standard Hash_Map class
* @tparam K type of the keys
* @tparam V type of the values
*/
template<typename K, typename V>
struct Hash_Map : std::map<K, V, std::less<K>,
  Eigen::aligned_allocator<std::pair<K, V> > > {};
#endif

#if defined OPENMVG_STD_UNORDERED_MAP

/**
* @brief Standard Hash_Map class
* @tparam K type of the keys
* @tparam V type of the values
*/
template<typename Key, typename Value>
struct Hash_Map : std::unordered_map<Key, Value> {};
#endif

} // namespace openMVG

#endif  // OPENMVG_TYPES_H_
