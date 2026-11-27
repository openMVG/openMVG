// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <functional>
#include <string>
#include <type_traits>
#include <utility>

// Client headers may inspect hash support before including OpenMVG.
using ClientPair = std::pair<int, int>;
constexpr bool kClientPairHashAvailable =
  std::is_default_constructible<std::hash<ClientPair>>::value;

#include "openMVG/types.hpp"
#include "third_party/CppUnitLite/TestHarness.h"

static_assert(kClientPairHashAvailable ==
  std::is_default_constructible<std::hash<ClientPair>>::value,
  "Including OpenMVG must not change std::hash support for client types");

TEST(TYPES, NonPairMapsKeepTheStandardHashType)
{
  using Map = openMVG::Hash_Map<int, std::string>;
  static_assert(std::is_same<Map, std::unordered_map<int, std::string>>::value,
    "Non-pair maps must retain their original type");
  Map values;
  values[3] = "three";
  CHECK(values.at(3) == "three");
}

TEST(TYPES, PairMapsKeepSymmetricHashValues)
{
  openMVG::Hash_Map<openMVG::Pair, int> values;
  const openMVG::Pair forward(7, 19), reverse(19, 7);
  CHECK(values.hash_function()(forward) == values.hash_function()(reverse));
}

TEST(TYPES, PairOrientationRemainsDistinct)
{
  openMVG::Hash_Map<openMVG::Pair, int> values;
  values[{7, 19}] = 2;
  values[{19, 7}] = 5;
  CHECK(values.size() == 2);
  CHECK(values.at({7, 19}) == 2);
  CHECK(values.at({19, 7}) == 5);
}

TEST(TYPES, HeterogeneousPairKeysRemainSupported)
{
  openMVG::Hash_Map<std::pair<std::string, int>, int> values;
  values[{"camera", 12}] = 8;
  CHECK(values.at({"camera", 12}) == 8);
  CHECK(values.find({"camera", 13}) == values.end());
}

TEST(TYPES, NestedPairKeysRemainSupported)
{
  using Key = std::pair<openMVG::Pair, openMVG::IndexT>;
  openMVG::Hash_Map<Key, int> values;
  values[{{7, 19}, 3}] = 8;
  CHECK(values.at({{7, 19}, 3}) == 8);
  CHECK(values.find({{7, 19}, 4}) == values.end());
}

int main()
{
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
