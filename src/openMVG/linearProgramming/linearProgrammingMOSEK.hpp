// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef MIMATTE_LINEAR_PROGRAMMING_INTERFACE_MOSEK_H_
#define MIMATTE_LINEAR_PROGRAMMING_INTERFACE_MOSEK_H_

#ifdef OPENMVG_HAVE_MOSEK

#include "openMVG/numeric/numeric.h"
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
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


#endif // MIMATTE_LINEAR_PROGRAMMING_INTERFACE_MOSEK_H_

#endif // OPENMVG_HAVE_MOSEK
