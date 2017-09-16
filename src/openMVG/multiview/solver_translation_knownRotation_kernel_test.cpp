// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/solver_translation_knownRotation_kernel.hpp"
#include "openMVG/multiview/test_data_sets.hpp"

#include "testing/testing.h"

#include <vector>

using namespace openMVG;

// Estimate the translation for a pair of view for which the relative rotation is known
// Use a 2 correspondences based solver
TEST(Translation_knownRotation_Kernel, Multiview) {

  const int nViews = 10;
  const int nbPoints = 2;
  const NViewDataSet d = NRealisticCamerasRing(nViews, nbPoints,
    nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  // Solve the problem and check that fitted value are good enough
  for (int nCameraIndex = 2; nCameraIndex < nViews; ++nCameraIndex)
  {
    const Mat x0 = d._x[0];
    const Mat xCam = d._x[nCameraIndex];
    // coordinates does not need to be normalized since we have used a unit K matrix.

    // Compute GT (Ground Truth) motion
    Mat3 R_GT;
    Vec3 t_GT;
    RelativeCameraMotion(d._R[0], d._t[0], d._R[nCameraIndex], d._t[nCameraIndex], &R_GT, &t_GT);

    openMVG::translation::kernel::TranslationFromKnowRotationKernel kernel(x0, xCam, R_GT);

    const std::vector<uint32_t> samples = {0,1};
    std::vector<Vec3> vec_t;
    kernel.Fit(samples, &vec_t);

    CHECK_EQUAL(1, vec_t.size());

    // Check that the fitted model is compatible with the data
    // Here the distance to the epipolar line is used
    for (Mat::Index i = 0; i < x0.cols(); ++i) {
      EXPECT_NEAR(0.0, kernel.Error(i, vec_t[0]), 1e-8);
    }

    // Check that the GT translation and the estimated one are equal
    EXPECT_MATRIX_NEAR(t_GT.normalized(), vec_t[0].normalized(), 1e-8);
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
