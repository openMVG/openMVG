
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_ROTATION_AVERAGING_COMMON_H_
#define OPENMVG_MULTIVIEW_ROTATION_AVERAGING_COMMON_H_

#include "openMVG/numeric/numeric.h"

namespace openMVG   {
namespace rotation_averaging  {

// Representation of relative rotations data
struct RelRotationData {
  size_t i, j; // view's indices
  Mat3 Rij; // view's relative rotation
  float weight;

  RelRotationData(size_t i_=0, size_t j_=0, const	Mat3 & Rij_=Mat3::Identity(), float weight_=1.0f):
  i(i_), j(j_), Rij(Rij_), weight(weight_)
  {}
};

} // namespace rotation_averaging
} // namespace openMVG

#endif //OPENMVG_MULTIVIEW_ROTATION_AVERAGING_COMMON_H_

