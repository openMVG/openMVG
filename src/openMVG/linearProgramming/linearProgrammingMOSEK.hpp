// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINEAR_PROGRAMMING_LINEAR_PROGRAMMING_MOSEK_HPP
#define OPENMVG_LINEAR_PROGRAMMING_LINEAR_PROGRAMMING_MOSEK_HPP

#ifdef OPENMVG_HAVE_MOSEK

#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/numeric/numeric.h"
extern "C"{
#include "mosek.h"
}

#include <vector>

namespace openMVG   {
namespace linearProgramming  {

/// MOSEK wrapper for the LP_Solver
class MOSEK_SolveWrapper : public LP_Solver
{
public :
  MOSEK_SolveWrapper(int nbParams);

  ~MOSEK_SolveWrapper();

  //--
  // Inherited functions :
  //--

  bool setup(const LP_Constraints & constraints) override ;
  bool setup(const LP_Constraints_Sparse & constraints) override ;

  bool solve() override ;

  bool getSolution(std::vector<double> & estimatedParams) override ;

private :
  //MSKenv_t     env;
  MSKtask_t    task; // Solver object.
};


} // namespace linearProgramming
} // namespace openMVG

#endif // OPENMVG_HAVE_MOSEK

#endif // OPENMVG_LINEAR_PROGRAMMING_LINEAR_PROGRAMMING_MOSEK_HPP
