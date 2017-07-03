// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINEAR_PROGRAMMING_LINEAR_PROGRAMMING_OSI_X_HPP
#define OPENMVG_LINEAR_PROGRAMMING_LINEAR_PROGRAMMING_OSI_X_HPP

#include <memory>
#include <vector>

#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"

class OsiClpSolverInterface;

namespace openMVG   {
namespace linearProgramming  {

/// OSI_Clp wrapper for the LP_Solver
class OSI_X_SolverWrapper : public LP_Solver
{
public:
  explicit OSI_X_SolverWrapper(int nbParams);

  //--
  // Inherited functions:
  //--

  bool setup(const LP_Constraints & constraints) override;
  bool setup(const LP_Constraints_Sparse & constraints) override;

  bool solve() override;

  bool getSolution(std::vector<double> & estimatedParams) override;

private:
  std::shared_ptr<OsiClpSolverInterface> si;
};

using OSI_CLP_SolverWrapper = OSI_X_SolverWrapper;

} // namespace linearProgramming
} // namespace openMVG


#endif // OPENMVG_LINEAR_PROGRAMMING_LINEAR_PROGRAMMING_OSI_X_HPP
