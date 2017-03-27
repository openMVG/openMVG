// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/numeric/l1_solver_admm.hpp"
#include "openMVG/numeric/numeric.h"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include <iostream>
#include <random>

using namespace openMVG;
using namespace std;

TEST(L1Solver_ADMM, Decoding)
{
  const double dTolerance = 1e-8;
  // source length
  const unsigned int N = 256;
  // codeword length
  const unsigned int M = 4*N;
  // number of perturbations
  const unsigned int T (0.2f*M);
  // coding matrix
  const Eigen::MatrixXd G = Eigen::MatrixXd::Random(M,N);
  // source word
  const Eigen::VectorXd x = Eigen::VectorXd::Random(N);
  // code word
  const Eigen::VectorXd y = G*x;
  // channel: perturb T randomly chosen entries
  // Apply random pertubations to the initial guess
  Eigen::VectorXd observation = y;
  const Eigen::VectorXd pertubations = Eigen::VectorXd::Random(T);

  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(0, observation.size()-1);
  for (unsigned int i = 0; i < T; ++i) {
    observation(distribution(generator)) = pertubations(i);
  }

  // recover (The initial noisy guess)
  Eigen::VectorXd solution = (G.transpose()*G).inverse()*G.transpose()*observation;

  // Recover the code word
  L1Solver<Eigen::MatrixXd>::Options options;
  options.absolute_tolerance = dTolerance;
  L1Solver<Eigen::MatrixXd> l1_solver(options, G);
  l1_solver.Solve(observation, &solution);

  // Check the solution
  const Eigen::VectorXd residuals = G * solution - y;
  for (unsigned int i = 0; i < residuals.size(); ++i) {
    EXPECT_NEAR(residuals(i), 0.0, dTolerance);
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
