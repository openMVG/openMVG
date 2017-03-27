// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINEAR_PROGRAMMING_INTERFACE_HPP
#define OPENMVG_LINEAR_PROGRAMMING_INTERFACE_HPP

#include <utility>
#include <vector>

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG   {
namespace linearProgramming  {

/// Generic container for LP (Linear Programming problems).
/// Embed :
///  - objective function
///  - Constraints (coefficients, Sign, objective value),
///  - Bounds over parameter (<=, =, >=).
///  - minimize or maximize
///
struct LP_Constraints
{
  enum eLP_SIGN
  {
    LP_LESS_OR_EQUAL    = 1,  // (<=)
    LP_GREATER_OR_EQUAL = 2,  // (>=)
    LP_EQUAL            = 3   // (=)
  };

  LP_Constraints() {
    bminimize_ = false;
  }

  int nbParams_; // The number of parameter/variable in constraint.
  Mat constraint_mat_; // Constraint under Matrix form.
  Vec constraint_objective_; // Constraint objective value.
  std::vector<eLP_SIGN> vec_sign_; // Constraint sign.
  std::vector< std::pair<double, double> > vec_bounds_; // parameter/variable bounds.

  bool bminimize_; // minimize is true or maximize is false.
  std::vector<double> vec_cost_; // Objective function
};

/// Generic Sparse container for LP (Linear Programming problems).
/// Embed :
///  - Constraints (coefficients, Sign, objective value),
///  - Bounds over parameter (<=, =, >=).
/// Implementation differ from LP_Constraints, here constraints are
///  stored as a Sparse matrix.
///
struct LP_Constraints_Sparse
{
  LP_Constraints_Sparse() {
    bminimize_ = false;
  }

  // Variable part
  int nbParams_; // The number of parameter/variable in constraint.
  std::vector< std::pair<double, double> > vec_bounds_; // parameter/variable bounds.

  // Constraint part
  sRMat constraint_mat_; // Constraint under Matrix form.
  Vec constraint_objective_; // Constraint objective value.
  std::vector<LP_Constraints::eLP_SIGN> vec_sign_; // Constraint sign.

  bool bminimize_; // minimize is true or maximize is false.
  std::vector<double> vec_cost_; // Objective function
};

/// Generic LP solver (Linear Programming)
/// It's an interface to setup constraint and objective of a Linear Program.
/// Embed constraint setup, problem solving, and parameters getter.
class LP_Solver
{
public:

  LP_Solver(int nbParams):nbParams_(nbParams){}

  /// Setup constraint for the given library.
  virtual bool setup(const LP_Constraints & constraints) = 0;
  virtual bool setup(const LP_Constraints_Sparse & constraints) = 0;

  /// Setup the feasibility and found the solution that best fit the constraint.
  virtual bool solve() = 0;

  /// Get back solution. Call it after solve.
  virtual bool getSolution(std::vector<double> & estimatedParams) = 0;

protected :
  int nbParams_; // The number of parameter considered in constraint formulation.
};

} // namespace linearProgramming
} // namespace openMVG


#endif // OPENMVG_LINEAR_PROGRAMMING_INTERFACE_HPP
