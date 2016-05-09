
// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "CppUnitLite/TestHarness.h"
#include "openMVG/numeric/numeric.h"
#include "testing/testing.h"

#include "openMVG/numeric/l1_solver_admm.hpp"
#include "openMVG/numeric/l1_solver_decode_pd.hpp"

#include <iostream>
#include <random>

using namespace openMVG;
using namespace std;

// L1 magic sample: l1_decode_example
//
// Let A be a M x N matrix with full rank. Given y of R^M, the problem
// (PA) minx ||y - Ax||1
// finds the vector x of R^N such that the error y - Ax has minimum l1 norm
// (i.e. we are asking that the difference between Ax and y be sparse).
//
// Suppose we have a channel code that produces a codeword c = Ax for a message x. The
// message travels over the channel, and has an unknown number of its entries corrupted.
// The decoder observes y = c + e, where e is the corruption. If e is sparse enough, then
// the decoder can use (PA) to recover x exactly. When x, A, y have real-valued entries,
// (PA) can be recast as an LP.
//
// This problem arises in the context of channel coding
// (see E. J. Candes and T. Tao. "Decoding by linear programming". in IEEE Trans.
// Inform. Theory, December 2005).
//
TEST(L1Solver_L1Magic, Decoding)
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
  // Apply random pertubations to the initial guess.
  Eigen::VectorXd observation = y;
  const Eigen::VectorXd pertubations = Eigen::VectorXd::Random(T);

  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(0, observation.size()-1);
  for (int i = 0; i < T; ++i) {
    observation(distribution(generator)) = pertubations(i);
  }

  // recover (The initial noisy guess)
  Eigen::VectorXd solution = (G.transpose()*G).inverse()*G.transpose()*observation;

  // Recover the code word
  TRobustRegressionL1PD(G, observation, solution, dTolerance, 30);

  // Check the solution
  const Eigen::VectorXd residuals = G * solution - y;
  for (unsigned int i = 0; i < residuals.size(); ++i) {
    EXPECT_NEAR(residuals(i), 0.0, dTolerance);
  }
}

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
  for (int i = 0; i < T; ++i) {
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
