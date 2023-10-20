// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
//:\file
//\author Ricardo Fabbri Rio de Janeiro State U. (rfabbri.github.io) 
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Pierre MOULON
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRIFOCAL_UTIL_HPP
#define OPENMVG_MULTIVIEW_TRIFOCAL_UTIL_HPP

#include <minus/chicago-default.h>
#include <minus/internal-util.h>
#include "openMVG/multiview/trifocal/trifocal_model.hpp"

namespace openMVG {
namespace trifocal {


// Converts a trifocal_model to quaternion-translation format
// Assumes tt[0] is identity
inline void
tt2qt(const trifocal_model_t &tt, double tt_qt[M::nve])
{
  typedef MiNuS::minus_util<double> util;
  util::rotm2quat(tt[1].transpose().data(), tt_qt);
  util::rotm2quat(tt[2].transpose().data(), tt_qt+4);
  for (unsigned i=0; i < 3; ++i) {
    tt_qt[8+i]   = tt[1](i,3);
    tt_qt[8+3+i] = tt[2](i,3);
  }
}

// Finds a ground truth camera among a std::vector of possible ones
bool
probe_solutions(
    const std::vector<trifocal_model_t> &solutions, 
    trifocal_model_t &gt, 
    unsigned *solution_index);

} // namespace trifocal
} // namespace OpenMVG

#endif  // OPENMVG_GEOMETRY_TRIFOCAL_UTIL_HPP
