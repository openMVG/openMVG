// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2014, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINFINITY_COMPUTER_VISION_TRIPLET_TIJS_AND_KIS_KERNEL_HPP
#define OPENMVG_LINFINITY_COMPUTER_VISION_TRIPLET_TIJS_AND_KIS_KERNEL_HPP

#include <vector>
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {
namespace trifocal {
namespace kernel {

/// A trifocal tensor seen as 3 projective cameras
struct TrifocalTensorModel {

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Mat34 P1, P2, P3;

  static double Error
  (
    const TrifocalTensorModel & t,
    const Vec2 & pt1,
    const Vec2 & pt2,
    const Vec2 & pt3
  );
};

}  // namespace kernel
}  // namespace trifocal
}  // namespace openMVG

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST(openMVG::trifocal::kernel::TrifocalTensorModel)

namespace openMVG{

/// Solve the translations and the structure of a view-triplet that have known rotations
struct translations_Triplet_Solver {
  enum { MINIMUM_SAMPLES = 4 };
  enum { MAX_MODELS = 1 };

  /// Solve the computation of the "tensor".
  static void Solve
  (
    const Mat &pt0,
    const Mat & pt1,
    const Mat & pt2,
    const std::vector<Mat3> & vec_KR,
    std::vector<trifocal::kernel::TrifocalTensorModel> *P,
    const double ThresholdUpperBound
  );

  // Compute the residual of reprojections
  static double Error
  (
    const trifocal::kernel::TrifocalTensorModel & Tensor,
    const Vec2 & pt0,
    const Vec2 & pt1,
    const Vec2 & pt2
  );
};

} // namespace openMVG

#endif // OPENMVG_LINFINITY_COMPUTER_VISION_TRIPLET_TIJS_AND_KIS_KERNEL_HPP
