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
// Author: alex@karatarakis.com (Alexander Karatarakis)

#include <array>

#include "benchmark/benchmark.h"
#include "ceres/jet.h"

namespace ceres {

// Cycle the Jets to avoid caching effects in the benchmark.
template <class JetType>
class JetInputData {
  using T = typename JetType::Scalar;
  static constexpr std::size_t SIZE = 20;

 public:
  JetInputData() {
    for (int i = 0; i < static_cast<int>(SIZE); i++) {
      const T ti = static_cast<T>(i + 1);

      a_[i].a = T(1.1) * ti;
      a_[i].v.setRandom();

      b_[i].a = T(2.2) * ti;
      b_[i].v.setRandom();

      c_[i].a = T(3.3) * ti;
      c_[i].v.setRandom();

      d_[i].a = T(4.4) * ti;
      d_[i].v.setRandom();

      e_[i].a = T(5.5) * ti;
      e_[i].v.setRandom();

      scalar_a_[i] = T(1.1) * ti;
      scalar_b_[i] = T(2.2) * ti;
      scalar_c_[i] = T(3.3) * ti;
      scalar_d_[i] = T(4.4) * ti;
      scalar_e_[i] = T(5.5) * ti;
    }
  }

  void advance() { index_ = (index_ + 1) % SIZE; }

  const JetType& a() const { return a_[index_]; }
  const JetType& b() const { return b_[index_]; }
  const JetType& c() const { return c_[index_]; }
  const JetType& d() const { return d_[index_]; }
  const JetType& e() const { return e_[index_]; }
  T scalar_a() const { return scalar_a_[index_]; }
  T scalar_b() const { return scalar_b_[index_]; }
  T scalar_c() const { return scalar_c_[index_]; }
  T scalar_d() const { return scalar_d_[index_]; }
  T scalar_e() const { return scalar_e_[index_]; }

 private:
  std::size_t index_{0};
  std::array<JetType, SIZE> a_{};
  std::array<JetType, SIZE> b_{};
  std::array<JetType, SIZE> c_{};
  std::array<JetType, SIZE> d_{};
  std::array<JetType, SIZE> e_{};
  std::array<T, SIZE> scalar_a_;
  std::array<T, SIZE> scalar_b_;
  std::array<T, SIZE> scalar_c_;
  std::array<T, SIZE> scalar_d_;
  std::array<T, SIZE> scalar_e_;
};

template <std::size_t JET_SIZE, class Function>
static void JetBenchmarkHelper(benchmark::State& state, const Function& func) {
  using JetType = Jet<double, JET_SIZE>;
  JetInputData<JetType> data{};
  JetType out{};
  const int iterations = static_cast<int>(state.range(0));
  for (auto _ : state) {
    for (int i = 0; i < iterations; i++) {
      func(data, out);
      data.advance();
    }
  }
  benchmark::DoNotOptimize(out);
}

template <std::size_t JET_SIZE>
static void Addition(benchmark::State& state) {
  using JetType = Jet<double, JET_SIZE>;
  JetBenchmarkHelper<JET_SIZE>(
      state, [](const JetInputData<JetType>& d, JetType& out) {
        out += +d.a() + d.b() + d.c() + d.d() + d.e();
      });
}
BENCHMARK_TEMPLATE(Addition, 3)->Arg(1000);
BENCHMARK_TEMPLATE(Addition, 10)->Arg(1000);
BENCHMARK_TEMPLATE(Addition, 15)->Arg(1000);
BENCHMARK_TEMPLATE(Addition, 25)->Arg(1000);
BENCHMARK_TEMPLATE(Addition, 32)->Arg(1000);
BENCHMARK_TEMPLATE(Addition, 200)->Arg(160);

template <std::size_t JET_SIZE>
static void AdditionScalar(benchmark::State& state) {
  using JetType = Jet<double, JET_SIZE>;
  JetBenchmarkHelper<JET_SIZE>(
      state, [](const JetInputData<JetType>& d, JetType& out) {
        out +=
            d.scalar_a() + d.scalar_b() + d.c() + d.scalar_d() + d.scalar_e();
      });
}
BENCHMARK_TEMPLATE(AdditionScalar, 3)->Arg(1000);
BENCHMARK_TEMPLATE(AdditionScalar, 10)->Arg(1000);
BENCHMARK_TEMPLATE(AdditionScalar, 15)->Arg(1000);
BENCHMARK_TEMPLATE(AdditionScalar, 25)->Arg(1000);
BENCHMARK_TEMPLATE(AdditionScalar, 32)->Arg(1000);
BENCHMARK_TEMPLATE(AdditionScalar, 200)->Arg(160);

template <std::size_t JET_SIZE>
static void Subtraction(benchmark::State& state) {
  using JetType = Jet<double, JET_SIZE>;
  JetBenchmarkHelper<JET_SIZE>(
      state, [](const JetInputData<JetType>& d, JetType& out) {
        out -= -d.a() - d.b() - d.c() - d.d() - d.e();
      });
}
BENCHMARK_TEMPLATE(Subtraction, 3)->Arg(1000);
BENCHMARK_TEMPLATE(Subtraction, 10)->Arg(1000);
BENCHMARK_TEMPLATE(Subtraction, 15)->Arg(1000);
BENCHMARK_TEMPLATE(Subtraction, 25)->Arg(1000);
BENCHMARK_TEMPLATE(Subtraction, 32)->Arg(1000);
BENCHMARK_TEMPLATE(Subtraction, 200)->Arg(160);

template <std::size_t JET_SIZE>
static void SubtractionScalar(benchmark::State& state) {
  using JetType = Jet<double, JET_SIZE>;
  JetBenchmarkHelper<JET_SIZE>(
      state, [](const JetInputData<JetType>& d, JetType& out) {
        out -=
            -d.scalar_a() - d.scalar_b() - d.c() - d.scalar_d() - d.scalar_e();
      });
}
BENCHMARK_TEMPLATE(SubtractionScalar, 3)->Arg(1000);
BENCHMARK_TEMPLATE(SubtractionScalar, 10)->Arg(1000);
BENCHMARK_TEMPLATE(SubtractionScalar, 15)->Arg(1000);
BENCHMARK_TEMPLATE(SubtractionScalar, 25)->Arg(1000);
BENCHMARK_TEMPLATE(SubtractionScalar, 32)->Arg(1000);
BENCHMARK_TEMPLATE(SubtractionScalar, 200)->Arg(160);

template <std::size_t JET_SIZE>
static void Multiplication(benchmark::State& state) {
  using JetType = Jet<double, JET_SIZE>;
  JetBenchmarkHelper<JET_SIZE>(
      state, [](const JetInputData<JetType>& d, JetType& out) {
        out *= d.a() * d.b() * d.c() * d.d() * d.e();
      });
}
BENCHMARK_TEMPLATE(Multiplication, 3)->Arg(1000);
BENCHMARK_TEMPLATE(Multiplication, 10)->Arg(1000);
BENCHMARK_TEMPLATE(Multiplication, 15)->Arg(1000);
BENCHMARK_TEMPLATE(Multiplication, 25)->Arg(1000);
BENCHMARK_TEMPLATE(Multiplication, 32)->Arg(1000);
BENCHMARK_TEMPLATE(Multiplication, 200)->Arg(160);

template <std::size_t JET_SIZE>
static void MultiplicationLeftScalar(benchmark::State& state) {
  using JetType = Jet<double, JET_SIZE>;
  JetBenchmarkHelper<JET_SIZE>(
      state, [](const JetInputData<JetType>& d, JetType& out) {
        out += d.scalar_a() *
               (d.scalar_b() * (d.scalar_c() * (d.scalar_d() * d.e())));
      });
}
BENCHMARK_TEMPLATE(MultiplicationLeftScalar, 3)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplicationLeftScalar, 10)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplicationLeftScalar, 15)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplicationLeftScalar, 25)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplicationLeftScalar, 32)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplicationLeftScalar, 200)->Arg(160);

template <std::size_t JET_SIZE>
static void MultiplicationRightScalar(benchmark::State& state) {
  using JetType = Jet<double, JET_SIZE>;
  JetBenchmarkHelper<JET_SIZE>(
      state, [](const JetInputData<JetType>& d, JetType& out) {
        out += (((d.a() * d.scalar_b()) * d.scalar_c()) * d.scalar_d()) *
               d.scalar_e();
      });
}
BENCHMARK_TEMPLATE(MultiplicationRightScalar, 3)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplicationRightScalar, 10)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplicationRightScalar, 15)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplicationRightScalar, 25)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplicationRightScalar, 32)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplicationRightScalar, 200)->Arg(160);

template <std::size_t JET_SIZE>
static void Division(benchmark::State& state) {
  using JetType = Jet<double, JET_SIZE>;
  JetBenchmarkHelper<JET_SIZE>(
      state, [](const JetInputData<JetType>& d, JetType& out) {
        out /= d.a() / d.b() / d.c() / d.d() / d.e();
      });
}
BENCHMARK_TEMPLATE(Division, 3)->Arg(1000);
BENCHMARK_TEMPLATE(Division, 10)->Arg(1000);
BENCHMARK_TEMPLATE(Division, 15)->Arg(1000);
BENCHMARK_TEMPLATE(Division, 25)->Arg(1000);
BENCHMARK_TEMPLATE(Division, 32)->Arg(1000);
BENCHMARK_TEMPLATE(Division, 200)->Arg(160);

template <std::size_t JET_SIZE>
static void DivisionLeftScalar(benchmark::State& state) {
  using JetType = Jet<double, JET_SIZE>;
  JetBenchmarkHelper<JET_SIZE>(
      state, [](const JetInputData<JetType>& d, JetType& out) {
        out += d.scalar_a() /
               (d.scalar_b() / (d.scalar_c() / (d.scalar_d() / d.e())));
      });
}
BENCHMARK_TEMPLATE(DivisionLeftScalar, 3)->Arg(1000);
BENCHMARK_TEMPLATE(DivisionLeftScalar, 10)->Arg(1000);
BENCHMARK_TEMPLATE(DivisionLeftScalar, 15)->Arg(1000);
BENCHMARK_TEMPLATE(DivisionLeftScalar, 25)->Arg(1000);
BENCHMARK_TEMPLATE(DivisionLeftScalar, 32)->Arg(1000);
BENCHMARK_TEMPLATE(DivisionLeftScalar, 200)->Arg(160);

template <std::size_t JET_SIZE>
static void DivisionRightScalar(benchmark::State& state) {
  using JetType = Jet<double, JET_SIZE>;
  JetBenchmarkHelper<JET_SIZE>(
      state, [](const JetInputData<JetType>& d, JetType& out) {
        out += (((d.a() / d.scalar_b()) / d.scalar_c()) / d.scalar_d()) /
               d.scalar_e();
      });
}
BENCHMARK_TEMPLATE(DivisionRightScalar, 3)->Arg(1000);
BENCHMARK_TEMPLATE(DivisionRightScalar, 10)->Arg(1000);
BENCHMARK_TEMPLATE(DivisionRightScalar, 15)->Arg(1000);
BENCHMARK_TEMPLATE(DivisionRightScalar, 25)->Arg(1000);
BENCHMARK_TEMPLATE(DivisionRightScalar, 32)->Arg(1000);
BENCHMARK_TEMPLATE(DivisionRightScalar, 200)->Arg(160);

template <std::size_t JET_SIZE>
static void MultiplyAndAdd(benchmark::State& state) {
  using JetType = Jet<double, JET_SIZE>;
  JetBenchmarkHelper<JET_SIZE>(
      state, [](const JetInputData<JetType>& d, JetType& out) {
        out += d.scalar_a() * d.a() + d.scalar_b() * d.b() +
               d.scalar_c() * d.c() + d.scalar_d() * d.d() +
               d.scalar_e() * d.e();
      });
}
BENCHMARK_TEMPLATE(MultiplyAndAdd, 3)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplyAndAdd, 10)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplyAndAdd, 15)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplyAndAdd, 25)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplyAndAdd, 32)->Arg(1000);
BENCHMARK_TEMPLATE(MultiplyAndAdd, 200)->Arg(160);

}  // namespace ceres

BENCHMARK_MAIN();
