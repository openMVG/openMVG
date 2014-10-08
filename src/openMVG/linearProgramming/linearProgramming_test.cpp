// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <algorithm>
#include <iostream>
#include <vector>

#include "testing/testing.h"
#include "openMVG/linearProgramming/linearProgrammingMOSEK.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"


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
  cstraint._nbParams = 2;
  cstraint._bminimize = false;

  //Configure objective
  cstraint._vec_cost.push_back(143);
  cstraint._vec_cost.push_back(60);

  cstraint._constraintMat = Mat(5,2);
  cstraint._constraintMat <<
    120, 210,
    110, 30,
    1, 1,
    1, 0,
    0, 1;

  cstraint._Cst_objective = Vec(5);
  cstraint._Cst_objective << 15000, 4000, 75, 0, 0;

  cstraint._vec_sign.resize(5);
  std::fill_n(cstraint._vec_sign.begin(), 3, LP_Constraints::LP_LESS_OR_EQUAL);
  std::fill_n(cstraint._vec_sign.begin()+3, 2, LP_Constraints::LP_GREATER_OR_EQUAL);

  cstraint._vec_bounds = std::vector< std::pair<double,double> >(cstraint._nbParams);
  fill(cstraint._vec_bounds.begin(),cstraint._vec_bounds.end(),
      std::make_pair((double)-1e+30, (double)1e+30));
}

#ifdef OPENMVG_HAVE_MOSEK
// LP_Solve website example solving with the HighLevelFramework
TEST(linearProgramming, MOSEK_dense_sample) {

  LP_Constraints cstraint;
  BuildLinearProblem(cstraint);

  //Solve
  std::vector<double> vec_solution(2);
  MOSEK_SolveWrapper solver(2);
  solver.setup(cstraint);

  EXPECT_TRUE(solver.solve());
  solver.getSolution(vec_solution);

  EXPECT_NEAR( 21.875000, vec_solution[0], 1e-6);
  EXPECT_NEAR( 53.125000, vec_solution[1], 1e-6);

  std::cout << "Solution : " << vec_solution[0] << " " << vec_solution[1] << std::endl;
}
#endif // OPENMVG_HAVE_MOSEK

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
  cstraint._nbParams = 4; // {x0, x1, x2, x3}

  // Constraint coefficient
  sRMat & A = cstraint._constraintMat;
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
  Vec & C = cstraint._Cst_objective;
  C.resize(3, 1);
  C[0] = 30;
  C[1] = 15;
  C[2] = 25;

  // Constraint sign
  std::vector<LP_Constraints::eLP_SIGN> & vec_sign = cstraint._vec_sign;
  vec_sign.resize(3);
  vec_sign[0] = LP_Constraints::LP_EQUAL;
  vec_sign[1] = LP_Constraints::LP_GREATER_OR_EQUAL;
  vec_sign[2] = LP_Constraints::LP_LESS_OR_EQUAL;

  // Variable bounds
  cstraint._vec_bounds = std::vector< std::pair<double,double> >(4);
  fill(cstraint._vec_bounds.begin(),cstraint._vec_bounds.end(),
      std::make_pair(0.0, (double)1e+30));
  cstraint._vec_bounds[1].second = 10;

  // Objective to maximize
  cstraint._bminimize = false;
  cstraint._vec_cost.resize(4);
  cstraint._vec_cost[0] = 3;
  cstraint._vec_cost[1] = 1;
  cstraint._vec_cost[2] = 5;
  cstraint._vec_cost[3] = 1;
}

#ifdef OPENMVG_HAVE_MOSEK
// Unit test on mosek Sparse constraint
TEST(linearProgramming, mosek_sparse_sample) {

  LP_Constraints_Sparse cstraint;
  BuildSparseLinearProblem(cstraint);

  //Solve
  std::vector<double> vec_solution(4);
  MOSEK_SolveWrapper solver(4);
  solver.setup(cstraint);

  EXPECT_TRUE(solver.solve());
  solver.getSolution(vec_solution);

  EXPECT_NEAR( 0.00, vec_solution[0], 1e-2);
  EXPECT_NEAR( 0.00, vec_solution[1], 1e-2);
  EXPECT_NEAR( 15, vec_solution[2], 1e-2);
  EXPECT_NEAR( 8.33, vec_solution[3], 1e-2);
}
#endif // #ifdef OPENMVG_HAVE_MOSEK

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

#ifdef OPENMVG_HAVE_MOSEK
TEST(linearProgramming, osi_mosek_sparse_sample) {

  LP_Constraints_Sparse cstraint;
  BuildSparseLinearProblem(cstraint);

  //Solve
  std::vector<double> vec_solution(4);
  OSI_MOSEK_SolverWrapper solver(4);
  solver.setup(cstraint);

  EXPECT_TRUE(solver.solve());
  solver.getSolution(vec_solution);

  EXPECT_NEAR( 0.00, vec_solution[0], 1e-2);
  EXPECT_NEAR( 0.00, vec_solution[1], 1e-2);
  EXPECT_NEAR( 15, vec_solution[2], 1e-2);
  EXPECT_NEAR( 8.33, vec_solution[3], 1e-2);
}
#endif // #ifdef OPENMVG_HAVE_MOSEK



/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
