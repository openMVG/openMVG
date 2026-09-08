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

#include "benchmark/benchmark.h"
#include "ceres/context_impl.h"
#include "ceres/internal/eigen.h"
#include "ceres/parallel_for.h"
#include "glog/logging.h"

namespace ceres::internal {

// Parallel for with very small amount of work per iteration and small amount of
// iterations benchmarks performance of task scheduling
static void SchedulerBenchmark(benchmark::State& state) {
  const int vector_size = static_cast<int>(state.range(0));
  const int num_threads = static_cast<int>(state.range(1));
  ContextImpl context;
  context.EnsureMinimumThreads(num_threads);

  Vector x = Vector::Random(vector_size);
  for (auto _ : state) {
    ParallelFor(
        &context, 0, vector_size, num_threads, [&x](int id) { x[id] = 0.; });
  }
  CHECK_EQ(x.squaredNorm(), 0.);
}
BENCHMARK(SchedulerBenchmark)
    ->Args({128, 1})
    ->Args({128, 2})
    ->Args({128, 4})
    ->Args({128, 8})
    ->Args({128, 16})
    ->Args({256, 1})
    ->Args({256, 2})
    ->Args({256, 4})
    ->Args({256, 8})
    ->Args({256, 16})
    ->Args({1024, 1})
    ->Args({1024, 2})
    ->Args({1024, 4})
    ->Args({1024, 8})
    ->Args({1024, 16})
    ->Args({4096, 1})
    ->Args({4096, 2})
    ->Args({4096, 4})
    ->Args({4096, 8})
    ->Args({4096, 16});

}  // namespace ceres::internal

BENCHMARK_MAIN();
