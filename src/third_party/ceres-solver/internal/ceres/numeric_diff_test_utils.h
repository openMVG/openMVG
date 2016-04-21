// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
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
// Author: sameeragarwal@google.com (Sameer Agarwal)

#ifndef CERES_INTERNAL_NUMERIC_DIFF_TEST_UTILS_H_
#define CERES_INTERNAL_NUMERIC_DIFF_TEST_UTILS_H_

#include "ceres/cost_function.h"
#include "ceres/sized_cost_function.h"
#include "ceres/types.h"

namespace ceres {
namespace internal {

// Noise factor for randomized cost function.
static const double kNoiseFactor = 0.01;

// Default random seed for randomized cost function.
static const unsigned int kRandomSeed = 1234;

// y1 = x1'x2      -> dy1/dx1 = x2,               dy1/dx2 = x1
// y2 = (x1'x2)^2  -> dy2/dx1 = 2 * x2 * (x1'x2), dy2/dx2 = 2 * x1 * (x1'x2)
// y3 = x2'x2      -> dy3/dx1 = 0,                dy3/dx2 = 2 * x2
class EasyFunctor {
 public:
  bool operator()(const double* x1, const double* x2, double* residuals) const;
  void ExpectCostFunctionEvaluationIsNearlyCorrect(
      const CostFunction& cost_function,
      NumericDiffMethodType method) const;
};

class EasyCostFunction : public SizedCostFunction<3, 5, 5> {
 public:
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** /* not used */) const {
    return functor_(parameters[0], parameters[1], residuals);
  }

 private:
  EasyFunctor functor_;
};

// y1 = sin(x1'x2)
// y2 = exp(-x1'x2 / 10)
//
// dy1/dx1 =  x2 * cos(x1'x2),            dy1/dx2 =  x1 * cos(x1'x2)
// dy2/dx1 = -x2 * exp(-x1'x2 / 10) / 10, dy2/dx2 = -x2 * exp(-x1'x2 / 10) / 10
class TranscendentalFunctor {
 public:
  bool operator()(const double* x1, const double* x2, double* residuals) const;
  void ExpectCostFunctionEvaluationIsNearlyCorrect(
      const CostFunction& cost_function,
      NumericDiffMethodType method) const;
};

class TranscendentalCostFunction : public SizedCostFunction<2, 5, 5> {
 public:
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** /* not used */) const {
    return functor_(parameters[0], parameters[1], residuals);
  }
 private:
  TranscendentalFunctor functor_;
};

// y = exp(x), dy/dx = exp(x)
class ExponentialFunctor {
 public:
  bool operator()(const double* x1, double* residuals) const;
  void ExpectCostFunctionEvaluationIsNearlyCorrect(
      const CostFunction& cost_function) const;
};

class ExponentialCostFunction : public SizedCostFunction<1, 1> {
 public:
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** /* not used */) const {
    return functor_(parameters[0], residuals);
  }

 private:
  ExponentialFunctor functor_;
};

// Test adaptive numeric differentiation by synthetically adding random noise
// to a functor.
// y = x^2 + [random noise], dy/dx ~ 2x
class RandomizedFunctor {
 public:
  RandomizedFunctor(double noise_factor, unsigned int random_seed)
      : noise_factor_(noise_factor), random_seed_(random_seed) {
  }

  bool operator()(const double* x1, double* residuals) const;
  void ExpectCostFunctionEvaluationIsNearlyCorrect(
      const CostFunction& cost_function) const;

 private:
  double noise_factor_;
  unsigned int random_seed_;
};

class RandomizedCostFunction : public SizedCostFunction<1, 1> {
 public:
  RandomizedCostFunction(double noise_factor, unsigned int random_seed)
      : functor_(noise_factor, random_seed) {
  }

  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** /* not used */) const {
    return functor_(parameters[0], residuals);
  }

 private:
  RandomizedFunctor functor_;
};


}  // namespace internal
}  // namespace ceres

#endif  // CERES_INTERNAL_NUMERIC_DIFF_TEST_UTILS_H_
