// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_TYPES_HPP
#define OPENMVG_TYPES_HPP

#include <Eigen/Core>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <vector>

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

/**
* @brief Standard Hash_Map class
* @tparam K type of the keys
* @tparam V type of the values
*/
template<typename Key, typename Value>
using Hash_Map = std::map<Key, Value, std::less<Key>,
  Eigen::aligned_allocator<std::pair<const Key, Value>>>;

} // namespace openMVG

#endif  // OPENMVG_TYPES_HPP
