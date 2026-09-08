// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2023 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: jodebo_beck@gmx.de (Johannes Beck)
//         sergiu.deitsch@gmail.com (Sergiu Deitsch)

#include "ceres/internal/integer_sequence_algorithm.h"

#include <type_traits>
#include <utility>

#include "ceres/internal/jet_traits.h"

namespace ceres::internal {

// Unit tests for exclusive scan of integer sequence.
static_assert(std::is_same<ExclusiveScan<std::integer_sequence<int>>,
                           std::integer_sequence<int>>::value,
              "Unit test of calculating the exclusive scan of an integer "
              "sequence failed.");
static_assert(std::is_same<ExclusiveScan<std::integer_sequence<int, 2>>,
                           std::integer_sequence<int, 0>>::value,
              "Unit test of calculating the exclusive scan of an integer "
              "sequence failed.");
static_assert(std::is_same<ExclusiveScan<std::integer_sequence<int, 2, 1>>,
                           std::integer_sequence<int, 0, 2>>::value,
              "Unit test of calculating the exclusive scan of an integer "
              "sequence failed.");
static_assert(std::is_same<ExclusiveScan<std::integer_sequence<int, 2, 1, 10>>,
                           std::integer_sequence<int, 0, 2, 3>>::value,
              "Unit test of calculating the exclusive scan of an integer "
              "sequence failed.");

using Ranks001 = Ranks_t<Jet<double, 0>, double, Jet<double, 1>>;
using Ranks1 = Ranks_t<Jet<double, 1>>;
using Ranks110 = Ranks_t<Jet<double, 1>, Jet<double, 1>, double>;
using Ranks023 = Ranks_t<double, Jet<double, 2>, Jet<double, 3>>;
using EmptyRanks = Ranks_t<>;

// Remove zero from the ranks integer sequence
using NonZeroRanks001 = RemoveValue_t<Ranks001, 0>;
using NonZeroRanks1 = RemoveValue_t<Ranks1, 0>;
using NonZeroRanks110 = RemoveValue_t<Ranks110, 0>;
using NonZeroRanks023 = RemoveValue_t<Ranks023, 0>;

static_assert(std::is_same<RemoveValue_t<EmptyRanks, 0>,
                           std::integer_sequence<int>>::value,
              "filtered sequence does not match an empty one");
static_assert(std::is_same<RemoveValue_t<std::integer_sequence<int, 2, 2>, 2>,
                           std::integer_sequence<int>>::value,
              "filtered sequence does not match an empty one");
static_assert(
    std::is_same<RemoveValue_t<std::integer_sequence<int, 0, 0, 2>, 2>,
                 std::integer_sequence<int, 0, 0>>::value,
    "filtered sequence does not match the expected one");
static_assert(
    std::is_same<RemoveValue_t<std::make_integer_sequence<int, 6>, 7>,
                 std::make_integer_sequence<int, 6>>::value,
    "sequence not containing the element to remove must not be transformed");
static_assert(
    std::is_same<NonZeroRanks001, std::integer_sequence<int, 1>>::value,
    "sequences do not match");
static_assert(std::is_same<NonZeroRanks1, std::integer_sequence<int, 1>>::value,
              "sequences do not match");
static_assert(
    std::is_same<NonZeroRanks110, std::integer_sequence<int, 1, 1>>::value,
    "sequences do not match");
static_assert(
    std::is_same<NonZeroRanks023, std::integer_sequence<int, 2, 3>>::value,
    "sequences do not match");
static_assert(std::is_same<RemoveValue_t<std::integer_sequence<long>, -1>,
                           std::integer_sequence<long>>::value,
              "sequences do not match");
static_assert(
    std::is_same<RemoveValue_t<std::integer_sequence<short, -2, -3, -1>, -1>,
                 std::integer_sequence<short, -2, -3>>::value,
    "sequences do not match");

using J = Jet<double, 2>;
template <typename T>
using J0 = Jet<T, 0>;
using J0d = J0<double>;

// Ensure all types match
static_assert(AreAllSame_v<int, int>, "types must be the same");
static_assert(AreAllSame_v<long, long, long>, "types must be the same");
static_assert(AreAllSame_v<J0d, J0d, J0d>, "types must be the same");
static_assert(!AreAllSame_v<double, int>, "types must not be the same");
static_assert(!AreAllSame_v<int, short, char>, "types must not be the same");

// Ensure all values in the integer sequence match
static_assert(AreAllEqual_v<int, 1, 1>,
              "integer sequence must contain same values");
static_assert(AreAllEqual_v<long, 2>,
              "integer sequence must contain one value");
static_assert(!AreAllEqual_v<short, 3, 4>,
              "integer sequence must not contain the same values");
static_assert(!AreAllEqual_v<unsigned, 3, 4, 3>,
              "integer sequence must not contain the same values");
static_assert(!AreAllEqual_v<int, 4, 4, 3>,
              "integer sequence must not contain the same values");

static_assert(IsEmptyOrAreAllEqual_v<std::integer_sequence<short>>,
              "expected empty sequence is not");
static_assert(IsEmptyOrAreAllEqual_v<std::integer_sequence<unsigned, 7, 7, 7>>,
              "expected all equal sequence is not");
static_assert(IsEmptyOrAreAllEqual_v<std::integer_sequence<int, 1>>,
              "expected all equal sequence is not");
static_assert(
    IsEmptyOrAreAllEqual_v<std::integer_sequence<long, 111, 111, 111, 111>>,
    "expected all equal sequence is not");

}  // namespace ceres::internal
