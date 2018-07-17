// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_TYPES_HPP
#define OPENMVG_TYPES_HPP

#ifndef OPENMVG_STD_UNORDERED_MAP

#include <Eigen/Core>

#endif

#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <vector>

#ifdef OPENMVG_STD_UNORDERED_MAP

#include <algorithm>
#include <unordered_map>
#include <utility>

#include "openMVG/stl/hash.hpp"
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

#endif // OPENMVG_STD_UNORDERED_MAP

/**
* @brief Main namespace of openMVG API
*/
namespace openMVG
{

/// Portable type used to store an index
using IndexT = uint32_t;

/// Portable value used to save an undefined index value
static const IndexT UndefinedIndexT = std::numeric_limits<IndexT>::max();

/// Standard Pair of IndexT
using Pair = std::pair<IndexT, IndexT>;

/// Set of Pair
using Pair_Set = std::set<Pair>;

/// Vector of Pair
using Pair_Vec = std::vector<Pair>;

#if defined OPENMVG_STD_UNORDERED_MAP

/**
* @brief Standard Hash_Map class
* @tparam K type of the keys
* @tparam V type of the values
*/
template<typename Key, typename Value>
using Hash_Map = std::unordered_map<Key, Value>;

#else

/**
* @brief Standard Hash_Map class
* @tparam K type of the keys
* @tparam V type of the values
*/
template<typename Key, typename Value>
using Hash_Map = std::map<Key, Value, std::less<Key>,
  Eigen::aligned_allocator<std::pair<const Key, Value>>>;

#endif // OPENMVG_STD_UNORDERED_MAP

} // namespace openMVG

#endif  // OPENMVG_TYPES_HPP
