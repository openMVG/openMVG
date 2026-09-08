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
// Author: sergiu.deitsch@gmail.com (Sergiu Deitsch)

#include "ceres/internal/jet_traits.h"

#include <Eigen/Core>
#include <type_traits>
#include <utility>

namespace ceres::internal {

using J = Jet<double, 2>;
// Don't care about the dual part for scalar part categorization and comparison
// tests
template <typename T>
using J0 = Jet<T, 0>;
using J0d = J0<double>;

// Extract the ranks of given types
using Ranks001 = Ranks_t<Jet<double, 0>, double, Jet<double, 1>>;
using Ranks1 = Ranks_t<Jet<double, 1>>;
using Ranks110 = Ranks_t<Jet<double, 1>, Jet<double, 1>, double>;
using Ranks023 = Ranks_t<double, Jet<double, 2>, Jet<double, 3>>;
using EmptyRanks = Ranks_t<>;

// Ensure extracted ranks match the expected integer sequence
static_assert(
    std::is_same<Ranks001, std::integer_sequence<int, 0, 0, 1>>::value,
    "ranks do not match");
static_assert(std::is_same<Ranks1, std::integer_sequence<int, 1>>::value,
              "ranks do not match");
static_assert(
    std::is_same<Ranks110, std::integer_sequence<int, 1, 1, 0>>::value,
    "ranks do not match");
static_assert(
    std::is_same<Ranks023, std::integer_sequence<int, 0, 2, 3>>::value,
    "ranks do not match");
static_assert(std::is_same<EmptyRanks, std::integer_sequence<int>>::value,
              "ranks sequence is not empty");

// Extract the underlying floating-point type
static_assert(std::is_same<UnderlyingScalar_t<double>, double>::value,
              "underlying type is not a double");
static_assert(std::is_same<UnderlyingScalar_t<J0d>, double>::value,
              "underlying type is not a double");
static_assert(std::is_same<UnderlyingScalar_t<J0<J0d>>, double>::value,
              "underlying type is not a double");
static_assert(std::is_same<UnderlyingScalar_t<J0<J0<J0d>>>, double>::value,
              "underlying type is not a double");

static_assert(CompatibleJetOperands_v<Jet<double, 1>, Jet<double, 1>>,
              "Jets must be compatible");
static_assert(CompatibleJetOperands_v<Jet<double, 1>, double>,
              "Jet and scalar must be compatible");
static_assert(CompatibleJetOperands_v<Jet<double, 2>>,
              "single Jet must be compatible");
static_assert(!CompatibleJetOperands_v<Jet<double, 1>, double, Jet<double, 2>>,
              "Jets and scalar must not be compatible");
static_assert(!CompatibleJetOperands_v<double, double>,
              "scalars must not be compatible");
static_assert(!CompatibleJetOperands_v<double>,
              "single scalar must not be compatible");
static_assert(!CompatibleJetOperands_v<>,
              "empty arguments must not be compatible");

static_assert(!PromotableJetOperands_v<double>,
              "single scalar must not be Jet promotable");
static_assert(!PromotableJetOperands_v<double, float, int>,
              "multiple scalars must not be Jet promotable");
static_assert(PromotableJetOperands_v<J0d, float, int>,
              "Jet and several scalars must be promotable");
static_assert(PromotableJetOperands_v<J0<J0d>, float, int>,
              "nested Jet and several scalars must be promotable");
static_assert(!PromotableJetOperands_v<Eigen::Array<double, 2, 3>, float, int>,
              "Eigen::Array must not be Jet promotable");
static_assert(!PromotableJetOperands_v<Eigen::Matrix<double, 3, 2>, float, int>,
              "Eigen::Matrix must not be Jet promotable");

}  // namespace ceres::internal
