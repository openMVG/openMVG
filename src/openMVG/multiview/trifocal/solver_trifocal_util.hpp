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

namespace openMVG {
namespace trifocal {


void
invert_intrinsics(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double px_coords[2], 
    double normalized_coords[2])
{
  const double *px = px_coords;
  double *nrm = normalized_coords;
  nrm[1] = (px[1] - K[1][2]) /K[1][1];
  nrm[0] = (px[0] - K[0][1]*nrm[1] - K[0][2])/K[0][0];
}

void
invert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double px_tgt_coords[2], 
    double normalized_tgt_coords[2])
{
  const double *tp = px_tgt_coords;
  double *t = normalized_tgt_coords;
  t[1] = tp[1]/K[1][1];
  t[0] = (tp[0] - K[0][1]*t[1])/K[0][0];
  // normalize -- works as a cache for angle computations / dot products
  // TODO: check inside minus if we are normalizing
  double n = hypot(t[0],t[1]);
  t[0] /= n; t[1] /= n;
}

// Converts a trifocal_model to quaternion-translation format
// Assumes tt[0] is identity
void
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
    unsigned *solution_index)
{
  double cameras_quat[M::nsols][M::nve];

  // translate trifocal_model from RC to QT (quaternion-translation)
  // - for each solution
  // -   translate to internal matrix form
  // -   call RC to QT

  tt2qt(gt, data::cameras_gt_quat_);
  for (unsigned s=0; s < solutions.size(); ++s)
    tt2qt(solutions[s], cameras_quat[s]);
  
  return io14::probe_all_solutions_quat(cameras_quat, data::cameras_gt_quat_, solutions.size(), solution_index);
}

} // namespace trifocal
} // namespace OpenMVG

#endif  // OPENMVG_GEOMETRY_TRIFOCAL_UTIL_HPP
