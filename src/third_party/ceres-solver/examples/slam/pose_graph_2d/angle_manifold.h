// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2023 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: vitus@google.com (Michael Vitus)
//         sameeragarwal@google.com (Sameer Agarwal)

#ifndef CERES_EXAMPLES_POSE_GRAPH_2D_ANGLE_MANIFOLD_H_
#define CERES_EXAMPLES_POSE_GRAPH_2D_ANGLE_MANIFOLD_H_

#include "ceres/autodiff_manifold.h"
#include "normalize_angle.h"

namespace ceres::examples {

// Defines a manifold for updating the angle to be constrained in [-pi to pi).
class AngleManifold {
 public:
  template <typename T>
  bool Plus(const T* x_radians,
            const T* delta_radians,
            T* x_plus_delta_radians) const {
    *x_plus_delta_radians = NormalizeAngle(*x_radians + *delta_radians);
    return true;
  }

  template <typename T>
  bool Minus(const T* y_radians,
             const T* x_radians,
             T* y_minus_x_radians) const {
    *y_minus_x_radians =
        NormalizeAngle(*y_radians) - NormalizeAngle(*x_radians);

    return true;
  }

  static ceres::Manifold* Create() {
    return new ceres::AutoDiffManifold<AngleManifold, 1, 1>;
  }
};

}  // namespace ceres::examples

#endif  // CERES_EXAMPLES_POSE_GRAPH_2D_ANGLE_MANIFOLD_H_
