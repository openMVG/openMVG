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
// Author: vitus@google.com (Michael Vitus)

#include "ceres/parallel_for.h"

#include <atomic>
#include <cmath>
#include <condition_variable>
#include <mutex>
#include <numeric>
#include <random>
#include <thread>
#include <tuple>
#include <vector>

#include "ceres/context_impl.h"
#include "ceres/internal/config.h"
#include "ceres/parallel_vector_ops.h"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace ceres::internal {

using testing::ElementsAreArray;
using testing::UnorderedElementsAreArray;

// Tests the parallel for loop computes the correct result for various number of
// threads.
TEST(ParallelFor, NumThreads) {
  ContextImpl context;
  context.EnsureMinimumThreads(/*num_threads=*/2);

  const int size = 16;
  std::vector<int> expected_results(size, 0);
  for (int i = 0; i < size; ++i) {
    expected_results[i] = std::sqrt(i);
  }

  for (int num_threads = 1; num_threads <= 8; ++num_threads) {
    std::vector<int> values(size, 0);
    ParallelFor(&context, 0, size, num_threads, [&values](int i) {
      values[i] = std::sqrt(i);
    });
    EXPECT_THAT(values, ElementsAreArray(expected_results));
  }
}

// Tests parallel for loop with ranges
TEST(ParallelForWithRange, NumThreads) {
  ContextImpl context;
  context.EnsureMinimumThreads(/*num_threads=*/2);

  const int size = 16;
  std::vector<int> expected_results(size, 0);
  for (int i = 0; i < size; ++i) {
    expected_results[i] = std::sqrt(i);
  }

  for (int num_threads = 1; num_threads <= 8; ++num_threads) {
    std::vector<int> values(size, 0);
    ParallelFor(
        &context, 0, size, num_threads, [&values](std::tuple<int, int> range) {
          auto [start, end] = range;
          for (int i = start; i < end; ++i) values[i] = std::sqrt(i);
        });
    EXPECT_THAT(values, ElementsAreArray(expected_results));
  }
}

// Tests parallel for loop with ranges and lower bound on minimal range size
TEST(ParallelForWithRange, MinimalSize) {
  ContextImpl context;
  constexpr int kNumThreads = 4;
  constexpr int kMinBlockSize = 5;
  context.EnsureMinimumThreads(kNumThreads);

  for (int size = kMinBlockSize; size <= 25; ++size) {
    std::atomic<bool> failed(false);
    ParallelFor(
        &context,
        0,
        size,
        kNumThreads,
        [&failed, kMinBlockSize](std::tuple<int, int> range) {
          auto [start, end] = range;
          if (end - start < kMinBlockSize) failed = true;
        },
        kMinBlockSize);
    EXPECT_EQ(failed, false);
  }
}

// Tests the parallel for loop with the thread ID interface computes the correct
// result for various number of threads.
TEST(ParallelForWithThreadId, NumThreads) {
  ContextImpl context;
  context.EnsureMinimumThreads(/*num_threads=*/2);

  const int size = 16;
  std::vector<int> expected_results(size, 0);
  for (int i = 0; i < size; ++i) {
    expected_results[i] = std::sqrt(i);
  }

  for (int num_threads = 1; num_threads <= 8; ++num_threads) {
    std::vector<int> values(size, 0);
    ParallelFor(
        &context, 0, size, num_threads, [&values](int thread_id, int i) {
          values[i] = std::sqrt(i);
        });
    EXPECT_THAT(values, ElementsAreArray(expected_results));
  }
}

// Tests nested for loops do not result in a deadlock.
TEST(ParallelFor, NestedParallelForDeadlock) {
  ContextImpl context;
  context.EnsureMinimumThreads(/*num_threads=*/2);

  // Increment each element in the 2D matrix.
  std::vector<std::vector<int>> x(3, {1, 2, 3});
  ParallelFor(&context, 0, 3, 2, [&x, &context](int i) {
    std::vector<int>& y = x.at(i);
    ParallelFor(&context, 0, 3, 2, [&y](int j) { ++y.at(j); });
  });

  const std::vector<int> results = {2, 3, 4};
  for (const std::vector<int>& value : x) {
    EXPECT_THAT(value, ElementsAreArray(results));
  }
}

// Tests nested for loops do not result in a deadlock for the parallel for with
// thread ID interface.
TEST(ParallelForWithThreadId, NestedParallelForDeadlock) {
  ContextImpl context;
  context.EnsureMinimumThreads(/*num_threads=*/2);

  // Increment each element in the 2D matrix.
  std::vector<std::vector<int>> x(3, {1, 2, 3});
  ParallelFor(&context, 0, 3, 2, [&x, &context](int thread_id, int i) {
    std::vector<int>& y = x.at(i);
    ParallelFor(&context, 0, 3, 2, [&y](int thread_id, int j) { ++y.at(j); });
  });

  const std::vector<int> results = {2, 3, 4};
  for (const std::vector<int>& value : x) {
    EXPECT_THAT(value, ElementsAreArray(results));
  }
}

TEST(ParallelForWithThreadId, UniqueThreadIds) {
  // Ensure the hardware supports more than 1 thread to ensure the test will
  // pass.
  const int num_hardware_threads = std::thread::hardware_concurrency();
  if (num_hardware_threads <= 1) {
    LOG(ERROR)
        << "Test not supported, the hardware does not support threading.";
    return;
  }

  ContextImpl context;
  context.EnsureMinimumThreads(/*num_threads=*/2);
  // Increment each element in the 2D matrix.
  std::vector<int> x(2, -1);
  std::mutex mutex;
  std::condition_variable condition;
  int count = 0;
  ParallelFor(&context,
              0,
              2,
              2,
              [&x, &mutex, &condition, &count](int thread_id, int i) {
                std::unique_lock<std::mutex> lock(mutex);
                x[i] = thread_id;
                ++count;
                condition.notify_all();
                condition.wait(lock, [&]() { return count == 2; });
              });

  EXPECT_THAT(x, UnorderedElementsAreArray({0, 1}));
}

// Helper function for partition tests
bool BruteForcePartition(
    int* costs, int start, int end, int max_partitions, int max_cost);
// Basic test if MaxPartitionCostIsFeasible and BruteForcePartition agree on
// simple test-cases
TEST(GuidedParallelFor, MaxPartitionCostIsFeasible) {
  std::vector<int> costs, cumulative_costs, partition;
  costs = {1, 2, 3, 5, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0};
  cumulative_costs.resize(costs.size());
  std::partial_sum(costs.begin(), costs.end(), cumulative_costs.begin());
  const auto dummy_getter = [](const int v) { return v; };

  // [1, 2, 3] [5], [0 ... 0, 7, 0, ... 0]
  EXPECT_TRUE(MaxPartitionCostIsFeasible(0,
                                         costs.size(),
                                         3,
                                         7,
                                         0,
                                         cumulative_costs.data(),
                                         dummy_getter,
                                         &partition));
  EXPECT_TRUE(BruteForcePartition(costs.data(), 0, costs.size(), 3, 7));
  // [1, 2, 3, 5, 0 ... 0, 7, 0, ... 0]
  EXPECT_TRUE(MaxPartitionCostIsFeasible(0,
                                         costs.size(),
                                         3,
                                         18,
                                         0,
                                         cumulative_costs.data(),
                                         dummy_getter,
                                         &partition));
  EXPECT_TRUE(BruteForcePartition(costs.data(), 0, costs.size(), 3, 18));
  // Impossible since there is item of cost 7
  EXPECT_FALSE(MaxPartitionCostIsFeasible(0,
                                          costs.size(),
                                          3,
                                          6,
                                          0,
                                          cumulative_costs.data(),
                                          dummy_getter,
                                          &partition));
  EXPECT_FALSE(BruteForcePartition(costs.data(), 0, costs.size(), 3, 6));
  // Impossible
  EXPECT_FALSE(MaxPartitionCostIsFeasible(0,
                                          costs.size(),
                                          2,
                                          10,
                                          0,
                                          cumulative_costs.data(),
                                          dummy_getter,
                                          &partition));
  EXPECT_FALSE(BruteForcePartition(costs.data(), 0, costs.size(), 2, 10));
}

// Randomized tests for MaxPartitionCostIsFeasible
TEST(GuidedParallelFor, MaxPartitionCostIsFeasibleRandomized) {
  std::vector<int> costs, cumulative_costs, partition;
  const auto dummy_getter = [](const int v) { return v; };

  // Random tests
  const int kNumTests = 1000;
  const int kMaxElements = 32;
  const int kMaxPartitions = 16;
  const int kMaxElCost = 8;
  std::mt19937 rng;
  std::uniform_int_distribution<int> rng_N(1, kMaxElements);
  std::uniform_int_distribution<int> rng_M(1, kMaxPartitions);
  std::uniform_int_distribution<int> rng_e(0, kMaxElCost);
  for (int t = 0; t < kNumTests; ++t) {
    const int N = rng_N(rng);
    const int M = rng_M(rng);
    int total = 0;
    costs.clear();
    for (int i = 0; i < N; ++i) {
      costs.push_back(rng_e(rng));
      total += costs.back();
    }

    cumulative_costs.resize(N);
    std::partial_sum(costs.begin(), costs.end(), cumulative_costs.begin());

    std::uniform_int_distribution<int> rng_seg(0, N - 1);
    int start = rng_seg(rng);
    int end = rng_seg(rng);
    if (start > end) std::swap(start, end);
    ++end;

    int first_admissible = 0;
    for (int threshold = 1; threshold <= total; ++threshold) {
      const bool bruteforce =
          BruteForcePartition(costs.data(), start, end, M, threshold);
      if (bruteforce && !first_admissible) {
        first_admissible = threshold;
      }
      const bool binary_search =
          MaxPartitionCostIsFeasible(start,
                                     end,
                                     M,
                                     threshold,
                                     start ? cumulative_costs[start - 1] : 0,
                                     cumulative_costs.data(),
                                     dummy_getter,
                                     &partition);
      EXPECT_EQ(bruteforce, binary_search);
      EXPECT_LE(partition.size(), M + 1);
      // check partition itself
      if (binary_search) {
        ASSERT_GT(partition.size(), 1);
        EXPECT_EQ(partition.front(), start);
        EXPECT_EQ(partition.back(), end);

        const int num_partitions = partition.size() - 1;
        EXPECT_LE(num_partitions, M);
        for (int j = 0; j < num_partitions; ++j) {
          int total = 0;
          for (int k = partition[j]; k < partition[j + 1]; ++k) {
            EXPECT_LT(k, end);
            EXPECT_GE(k, start);
            total += costs[k];
          }
          EXPECT_LE(total, threshold);
        }
      }
    }
  }
}

TEST(GuidedParallelFor, PartitionRangeForParallelFor) {
  std::vector<int> costs, cumulative_costs, partition;
  const auto dummy_getter = [](const int v) { return v; };

  // Random tests
  const int kNumTests = 1000;
  const int kMaxElements = 32;
  const int kMaxPartitions = 16;
  const int kMaxElCost = 8;
  std::mt19937 rng;
  std::uniform_int_distribution<int> rng_N(1, kMaxElements);
  std::uniform_int_distribution<int> rng_M(1, kMaxPartitions);
  std::uniform_int_distribution<int> rng_e(0, kMaxElCost);
  for (int t = 0; t < kNumTests; ++t) {
    const int N = rng_N(rng);
    const int M = rng_M(rng);
    int total = 0;
    costs.clear();
    for (int i = 0; i < N; ++i) {
      costs.push_back(rng_e(rng));
      total += costs.back();
    }

    cumulative_costs.resize(N);
    std::partial_sum(costs.begin(), costs.end(), cumulative_costs.begin());

    std::uniform_int_distribution<int> rng_seg(0, N - 1);
    int start = rng_seg(rng);
    int end = rng_seg(rng);
    if (start > end) std::swap(start, end);
    ++end;

    int first_admissible = 0;
    for (int threshold = 1; threshold <= total; ++threshold) {
      const bool bruteforce =
          BruteForcePartition(costs.data(), start, end, M, threshold);
      if (bruteforce) {
        first_admissible = threshold;
        break;
      }
    }
    EXPECT_TRUE(first_admissible != 0 || total == 0);
    partition = PartitionRangeForParallelFor(
        start, end, M, cumulative_costs.data(), dummy_getter);
    ASSERT_GT(partition.size(), 1);
    EXPECT_EQ(partition.front(), start);
    EXPECT_EQ(partition.back(), end);

    const int num_partitions = partition.size() - 1;
    EXPECT_LE(num_partitions, M);
    for (int j = 0; j < num_partitions; ++j) {
      int total = 0;
      for (int k = partition[j]; k < partition[j + 1]; ++k) {
        EXPECT_LT(k, end);
        EXPECT_GE(k, start);
        total += costs[k];
      }
      EXPECT_LE(total, first_admissible);
    }
  }
}

// Recursively try to partition range into segements of total cost
// less than max_cost
bool BruteForcePartition(
    int* costs, int start, int end, int max_partitions, int max_cost) {
  if (start == end) return true;
  if (start < end && max_partitions == 0) return false;
  int total_cost = 0;
  for (int last_curr = start + 1; last_curr <= end; ++last_curr) {
    total_cost += costs[last_curr - 1];
    if (total_cost > max_cost) break;
    if (BruteForcePartition(
            costs, last_curr, end, max_partitions - 1, max_cost))
      return true;
  }
  return false;
}

// Tests if guided parallel for loop computes the correct result for various
// number of threads.
TEST(GuidedParallelFor, NumThreads) {
  ContextImpl context;
  context.EnsureMinimumThreads(/*num_threads=*/2);

  const int size = 16;
  std::vector<int> expected_results(size, 0);
  for (int i = 0; i < size; ++i) {
    expected_results[i] = std::sqrt(i);
  }

  std::vector<int> costs, cumulative_costs;
  for (int i = 1; i <= size; ++i) {
    int cost = i * i;
    costs.push_back(cost);
    if (i == 1) {
      cumulative_costs.push_back(cost);
    } else {
      cumulative_costs.push_back(cost + cumulative_costs.back());
    }
  }

  for (int num_threads = 1; num_threads <= 8; ++num_threads) {
    std::vector<int> values(size, 0);
    ParallelFor(
        &context,
        0,
        size,
        num_threads,
        [&values](int i) { values[i] = std::sqrt(i); },
        cumulative_costs.data(),
        [](const int v) { return v; });
    EXPECT_THAT(values, ElementsAreArray(expected_results));
  }
}

TEST(ParallelAssign, D2MulX) {
  const int kVectorSize = 1024 * 1024;
  const int kMaxNumThreads = 8;
  const double kEpsilon = 1e-16;

  const Vector D_full = Vector::Random(kVectorSize * 2);
  const ConstVectorRef D(D_full.data() + kVectorSize, kVectorSize);
  const Vector x = Vector::Random(kVectorSize);
  const Vector y_expected = D.array().square() * x.array();
  ContextImpl context;
  context.EnsureMinimumThreads(kMaxNumThreads);

  for (int num_threads = 1; num_threads <= kMaxNumThreads; ++num_threads) {
    Vector y_observed(kVectorSize);
    ParallelAssign(
        &context, num_threads, y_observed, D.array().square() * x.array());

    // We might get non-bit-exact result due to different precision in scalar
    // and vector code. For example, in x86 mode mingw might emit x87
    // instructions for scalar code, thus making bit-exact check fail
    EXPECT_NEAR((y_expected - y_observed).squaredNorm(),
                0.,
                kEpsilon * y_expected.squaredNorm());
  }
}

TEST(ParallelAssign, SetZero) {
  const int kVectorSize = 1024 * 1024;
  const int kMaxNumThreads = 8;

  ContextImpl context;
  context.EnsureMinimumThreads(kMaxNumThreads);

  for (int num_threads = 1; num_threads <= kMaxNumThreads; ++num_threads) {
    Vector x = Vector::Random(kVectorSize);
    ParallelSetZero(&context, num_threads, x);

    CHECK_EQ(x.squaredNorm(), 0.);
  }
}

}  // namespace ceres::internal
