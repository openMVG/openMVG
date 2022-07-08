// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
// Copyright (c) 2019 Pierre MOULON.
// Copyright (c) 2022 Ricardo Fabbri and Gabriel Andrade
//
//\author Pierre MOULON
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Ricardo Fabbri Rio de Janeiro State U. (rfabbri.github.io)

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/trifocal/solver_trifocal_three_point.hpp"

#include <minus/minus.h>
#include <minus/chicago-default.h>

namespace openMVG {
namespace trifocal {

static unsigned constexpr max_solve_tries = 5;
using namespace MiNuS;

void Trifocal3PointPositionTangentialSolver::
Solve(const Mat &datum_0,
      const Mat &datum_1,
      const Mat &datum_2,
      std::vector<trifocal_model_t> *trifocal_tensor)
{
  double p[io::pp::nviews][io::pp::npoints][io::ncoords2d];
  double tgt[io::pp::nviews][io::pp::npoints][io::ncoords2d];

  // pack into solver's efficient representation
  for (unsigned ip=0; ip < io::pp::npoints; ++ip) {
      p[0][ip][0] = datum_0(0,ip);
      p[0][ip][1] = datum_0(1,ip);
    tgt[0][ip][0] = datum_0(2,ip);
    tgt[0][ip][1] = datum_0(3,ip);

      p[1][ip][0] = datum_1(0,ip);
      p[1][ip][1] = datum_1(1,ip);
    tgt[1][ip][0] = datum_1(2,ip);
    tgt[1][ip][1] = datum_1(3,ip);

      p[2][ip][0] = datum_2(0,ip);
      p[2][ip][1] = datum_2(1,ip);
    tgt[2][ip][0] = datum_2(2,ip);
    tgt[2][ip][1] = datum_2(3,ip);
  }

  unsigned nsols_final = 0;
  unsigned id_sols[M::nsols];
  double  cameras[M::nsols][io::pp::nviews-1][4][3];  // first camera is always [I | 0]
  for (unsigned i = 0; i < max_solve_tries; ++i)
    if (MiNuS::minus<chicago>::solve(p, tgt, cameras, id_sols, &nsols_final))
      break;

  std::vector<trifocal_model_t> &tt = *trifocal_tensor;
  tt.resize(nsols_final);
  for (unsigned s = 0; s < nsols_final; ++s) {
    tt[s][0] = Mat34::Identity();
    for (unsigned v=1; v < io::pp::nviews; ++v) {
        // eigen is col-major but minus is row-major, so memcpy cannot be used.
        for (unsigned ir = 0; ir < 3; ++ir)
          for (unsigned ic = 0; ic < 3; ++ic)
            tt[s][v](ir, ic) = cameras[id_sols[s]][v-1][ir][ic];
        for (unsigned r=0; r < 3; ++r)
          tt[s][v](r,3) = cameras[id_sols[s]][v-1][3][r];
    }
  }
  // TODO: filter the solutions by:
  // - positive depth and
  // - using tangent at 3rd point
}

} // namespace trifocal
} // namespace OpenMVG
