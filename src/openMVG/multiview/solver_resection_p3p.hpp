// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_RESECTION_P3P_HPP
#define OPENMVG_MULTIVIEW_RESECTION_P3P_HPP

#include <vector>

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {
namespace euclidean_resection {

struct P3PSolver {
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 4};
  // Solve the problem of camera pose.
  static void Solve
  (
    const Mat &pt2D,
    const Mat &pt3D,
    std::vector<Mat34> *models
  );

  // Compute the residual of the projection distance(pt2D, Project(P,pt3D))
  static double Error
  (
    const Mat34 & P,
    const Vec2 & pt2D,
    const Vec3 & pt3D
  );
};

class P3P_ResectionKernel_K {
 public:
  using Model = Mat34;
  enum { MINIMUM_SAMPLES = 3 };

  P3P_ResectionKernel_K
  (
    const Mat2X &x_camera,
    const Mat3X &X,
    const Mat3 &K = Mat3::Identity()
  );

  void Fit
  (
    const std::vector<uint32_t> &samples,
    std::vector<Model> *models
  ) const;

  double Error
  (
    size_t sample,
    const Model &model
  ) const;

  size_t NumSamples() const;

 private:
  Mat2X x_image_;   // camera coordinates
  Mat3X x_camera_;  // camera coordinates (normalized)
  Mat3X X_;         // 3D points
  Mat3 K_;
};

}  // namespace euclidean_resection
}  // namespace openMVG

#endif // OPENMVG_MULTIVIEW_RESECTION_P3P_HPP
