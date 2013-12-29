// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINEAR_PROGRAMMING_INTERFACE_H_
#define OPENMVG_LINEAR_PROGRAMMING_INTERFACE_H_

#include <vector>
#include <utility>
#include "openMVG/numeric/numeric.h"

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
    LP_EQUAL            = 3,   // (=)
    LP_FREE             = 4 //only supported in MOSEK
  };

  LP_Constraints() {
    _bminimize = false;
  }

  int _nbParams; // The number of parameter/variable in constraint.
  Mat _constraintMat; // Constraint under Matrix form.
  Vec _Cst_objective; // Constraint objective value.
  std::vector<eLP_SIGN> _vec_sign; // Constraint sign.
  std::vector< std::pair<double, double> > _vec_bounds; // parameter/variable bounds.

  bool _bminimize; // minimize is true or maximize is false.
  std::vector<double> _vec_cost; // Objective function
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
    _bminimize = false;
  }

  // Variable part
  int _nbParams; // The number of parameter/variable in constraint.
  std::vector< std::pair<double, double> > _vec_bounds; // parameter/variable bounds.

  // Constraint part
  sRMat _constraintMat; // Constraint under Matrix form.
  Vec _Cst_objective; // Constraint objective value.
  std::vector<LP_Constraints::eLP_SIGN> _vec_sign; // Constraint sign.

  bool _bminimize; // minimize is true or maximize is false.
  std::vector<double> _vec_cost; // Objective function
};

/// Generic LP solver (Linear Programming)
/// It's an interface to setup constraint and objective of a Linear Program.
/// Embed constraint setup, problem solving, and parameters getter.
class LP_Solver
{
public:

  LP_Solver(int nbParams):_nbParams(nbParams){};

  /// Setup constraint for the given library.
  virtual bool setup(const LP_Constraints & constraints) = 0;
  virtual bool setup(const LP_Constraints_Sparse & constraints) = 0;

  /// Setup the feasibility and found the solution that best fit the constraint.
  virtual bool solve() = 0;

  /// Get back solution. Call it after solve.
  virtual bool getSolution(std::vector<double> & estimatedParams) = 0;

protected :
  int _nbParams; // The number of parameter considered in constraint formulation.
};

} // namespace linearProgramming
} // namespace openMVG


#endif // OPENMVG_LINEAR_PROGRAMMING_INTERFACE_H_
