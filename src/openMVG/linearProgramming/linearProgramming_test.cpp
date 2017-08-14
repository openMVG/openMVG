// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "testing/testing.h"

#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"

#include <algorithm>
#include <limits>
#include <iostream>
#include <vector>


using namespace openMVG;
using namespace openMVG::linearProgramming;

// Setup :
// max(143x + 60y)
//  s.t.
//  120x + 210y <= 15000
//  110x + 30y <= 4000
//  x + y <= 75
//  x >= 0
//  y >= 0
void BuildLinearProblem(LP_Constraints & cstraint)
{
  cstraint.nbParams_ = 2;
  cstraint.bminimize_ = false;

  //Configure objective
  cstraint.vec_cost_.push_back(143);
  cstraint.vec_cost_.push_back(60);

  cstraint.constraint_mat_ = Mat(5,2);
  cstraint.constraint_mat_ <<
    120, 210,
    110, 30,
    1, 1,
    1, 0,
    0, 1;

  cstraint.constraint_objective_ = Vec(5);
  cstraint.constraint_objective_ << 15000, 4000, 75, 0, 0;

  cstraint.vec_sign_.resize(5);
  std::fill_n(cstraint.vec_sign_.begin(), 3, LP_Constraints::LP_LESS_OR_EQUAL);
  std::fill_n(cstraint.vec_sign_.begin()+3, 2, LP_Constraints::LP_GREATER_OR_EQUAL);

  cstraint.vec_bounds_.assign(cstraint.nbParams_,
    {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()});
}

TEST(linearProgramming, osiclp_dense_sample) {

  LP_Constraints cstraint;
  BuildLinearProblem(cstraint);

  //Solve
  std::vector<double> vec_solution(2);
  OSI_CLP_SolverWrapper solver(2);
  solver.setup(cstraint);

  EXPECT_TRUE(solver.solve());
  solver.getSolution(vec_solution);

  EXPECT_NEAR( 21.875000, vec_solution[0], 1e-6);
  EXPECT_NEAR( 53.125000, vec_solution[1], 1e-6);
}

// Setup example from MOSEK API
// maximize :
// 3 x0 + 1 x1 + 5 x2 + 1 x3
// subject to
// 3 x0 + 1 x1 + 2 x2         = 30
// 2 x0 + 1 x1 + 3 x2 + 1 x3 >= 15
//        2 w1        + 3 x3 <= 25
// bounds
// 0 <= x0, x2, x3 < infinity
// 0 <= x1 <= 10
void BuildSparseLinearProblem(LP_Constraints_Sparse & cstraint)
{
  // Number of variable we are looking for
  cstraint.nbParams_ = 4; // {x0, x1, x2, x3}

  // Constraint coefficient
  sRMat & A = cstraint.constraint_mat_;
  A.resize(3,4);
  A.coeffRef(0,0) = 3;
  A.coeffRef(0,1) = 1;
  A.coeffRef(0,2) = 2;

  A.coeffRef(1,0) = 2;
  A.coeffRef(1,1) = 1;
  A.coeffRef(1,2) = 3;
  A.coeffRef(1,3) = 1;

  A.coeffRef(2,1) = 2;
  A.coeffRef(2,3) = 3;

  // Constraint objective
  Vec & C = cstraint.constraint_objective_;
  C.resize(3, 1);
  C[0] = 30;
  C[1] = 15;
  C[2] = 25;

  // Constraint sign
  std::vector<LP_Constraints::eLP_SIGN> & vec_sign = cstraint.vec_sign_;
  vec_sign.resize(3);
  vec_sign[0] = LP_Constraints::LP_EQUAL;
  vec_sign[1] = LP_Constraints::LP_GREATER_OR_EQUAL;
  vec_sign[2] = LP_Constraints::LP_LESS_OR_EQUAL;

  // Variable bounds
  cstraint.vec_bounds_.assign(4, {0.0, std::numeric_limits<double>::max()});
  cstraint.vec_bounds_[1].second = 10;

  // Objective to maximize
  cstraint.bminimize_ = false;
  cstraint.vec_cost_.resize(4);
  cstraint.vec_cost_[0] = 3;
  cstraint.vec_cost_[1] = 1;
  cstraint.vec_cost_[2] = 5;
  cstraint.vec_cost_[3] = 1;
}

TEST(linearProgramming, osiclp_sparse_sample) {

  LP_Constraints_Sparse cstraint;
  BuildSparseLinearProblem(cstraint);

  //Solve
  std::vector<double> vec_solution(4);
  OSI_CLP_SolverWrapper solver(4);
  solver.setup(cstraint);

  EXPECT_TRUE(solver.solve());
  solver.getSolution(vec_solution);

  EXPECT_NEAR( 0.00, vec_solution[0], 1e-2);
  EXPECT_NEAR( 0.00, vec_solution[1], 1e-2);
  EXPECT_NEAR( 15, vec_solution[2], 1e-2);
  EXPECT_NEAR( 8.33, vec_solution[3], 1e-2);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
