// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/robust_estimation/robust_estimator_Ransac.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"

#include "CppUnitLite/TestHarness.h"

using namespace openMVG::robust;

TEST(RansacIterations, MeetsRequestedConfidence)
{
  const std::size_t count = IterationsRequired(2, 0.01, 0.5);
  CHECK_EQUAL(std::size_t{17}, count);
  CHECK(std::pow(0.75, count) <= 0.01);
  CHECK(std::pow(0.75, count - 1) > 0.01);
}

TEST(RansacIterations, SmallSuccessProbability)
{
  CHECK(IterationsRequired(5, 0.01, 0.0005) > 4096);
}

TEST(RansacIterations, SaturatesBeforeIntegerConversion)
{
  CHECK_EQUAL(std::numeric_limits<std::size_t>::max(),
    IterationsRequired(5, 0.01, 1e-10));
}

TEST(RansacIterations, UnderflowingSuccessProbability)
{
  CHECK_EQUAL(std::numeric_limits<std::size_t>::max(),
    IterationsRequired(5, 0.01, 1e-100));
}

TEST(RansacIterations, NoInliers)
{
  CHECK_EQUAL(std::numeric_limits<std::size_t>::max(),
    IterationsRequired(5, 0.01, 0.0));
}

TEST(RansacIterations, AllInliers)
{
  CHECK_EQUAL(std::size_t{1}, IterationsRequired(5, 0.01, 1.0));
}

namespace {

// Exercise the production sampler, scorer and RANSAC loop with a controlled
// minimal consensus. Every fitted model contains exactly its five samples.
struct MinimalConsensusKernel
{
  using Model = std::vector<uint32_t>;
  enum { MINIMUM_SAMPLES = 5 };

  explicit MinimalConsensusKernel(std::size_t count) : count_(count) {}
  std::size_t NumSamples() const { return count_; }
  void Fit(const std::vector<uint32_t>& samples, std::vector<Model>* models) const
  {
    ++fits;
    models->push_back(samples);
  }
  double Error(uint32_t sample, const Model& model) const
  {
    return std::find(model.begin(), model.end(), sample) == model.end() ? 1.0 : 0.0;
  }

  mutable std::size_t fits = 0;
  const std::size_t count_;
};

} // namespace

TEST(RansacIterations, LowInlierRatioRetainsIterationBudget)
{
  MinimalConsensusKernel kernel(10000);
  std::vector<uint32_t> inliers;
  std::size_t score = 0;
  RANSAC(kernel, ScorerEvaluator<MinimalConsensusKernel>(0.5), &inliers, &score);
  CHECK_EQUAL(std::size_t{4096}, kernel.fits);
  CHECK_EQUAL(std::size_t{5}, inliers.size());
  CHECK_EQUAL(std::size_t{5}, score);
}

TEST(RansacIterations, OrdinaryRatioRoundsBudgetUp)
{
  MinimalConsensusKernel kernel(10);
  RANSAC(kernel, ScorerEvaluator<MinimalConsensusKernel>(0.5));
  CHECK_EQUAL(std::size_t{146}, kernel.fits);
}

TEST(RansacIterations, CompleteConsensusStopsAfterOneFit)
{
  MinimalConsensusKernel kernel(5);
  RANSAC(kernel, ScorerEvaluator<MinimalConsensusKernel>(0.5));
  CHECK_EQUAL(std::size_t{1}, kernel.fits);
}

int main() { TestResult tr; return TestRegistry::runAllTests(tr); }
