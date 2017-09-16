// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/solver_resection_p3p.hpp"
#include "openMVG/multiview/test_data_sets.hpp"

#include "testing/testing.h"

using namespace openMVG;

TEST(Resection_Kernel_DLT, Multiview) {

  const int nViews = 3;
  const int nbPoints = 10;
  const NViewDataSet d =
    NRealisticCamerasRing(nViews, nbPoints,
                          nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  // Solve the problem and check that fitted value are good enough
  for (int nResectionCameraIndex = 0; nResectionCameraIndex < nViews; ++nResectionCameraIndex)
  {
    const Mat x = d._x[nResectionCameraIndex];
    const Mat X = d._X;
    openMVG::resection::kernel::PoseResectionKernel kernel(x, X);

    std::vector<Mat34> Ps;
    kernel.Fit({0,1,2,3,4,5}, &Ps);
    for (Mat::Index i = 0; i < x.cols(); ++i) {
      EXPECT_NEAR(0.0, kernel.Error(i, Ps[0]), 1e-8);
    }

    CHECK_EQUAL(1, Ps.size());

    // Check that Projection matrix is near to the GT:
    Mat34 GT_ProjectionMatrix = d.P(nResectionCameraIndex).array()
                                / d.P(nResectionCameraIndex).norm();
    Mat34 COMPUTED_ProjectionMatrix = Ps[0].array() / Ps[0].norm();
    EXPECT_MATRIX_NEAR(GT_ProjectionMatrix, COMPUTED_ProjectionMatrix, 1e-8);
  }
}

TEST(P3P_Kneip_CVPR11, Multiview) {

  const int nViews = 3;
  const int nbPoints = 3;
  const NViewDataSet d =
    NRealisticCamerasRing(nViews, nbPoints,
                          nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  // Solve the problem and check that fitted value are good enough
  for (int nResectionCameraIndex = 0; nResectionCameraIndex < nViews; ++nResectionCameraIndex)
  {
    const Mat x = d._x[nResectionCameraIndex];
    const Mat bearing_vectors = (d._K[0].inverse() * x.colwise().homogeneous()).colwise().normalized();
    const Mat X = d._X;
    openMVG::euclidean_resection::PoseResectionKernel_P3P_Kneip kernel(bearing_vectors, X);

    std::vector<Mat34> Ps;
    kernel.Fit({0,1,2}, &Ps); // 3 points sample are required, lets take the first three

    bool bFound = false;
    char index = -1;
    for (size_t i = 0; i < Ps.size(); ++i)  {
      Mat34 GT_ProjectionMatrix = d.P(nResectionCameraIndex).array()
                                / d.P(nResectionCameraIndex).norm();
      Mat34 COMPUTED_ProjectionMatrix = Ps[i].array() / Ps[i].norm();
      if ( NormLInfinity(GT_ProjectionMatrix - COMPUTED_ProjectionMatrix) < 1e-8 )
      {
        bFound = true;
        index = i;
      }
    }
    EXPECT_TRUE(bFound);

    // Check that for the found matrix the residual is small
    for (Mat::Index i = 0; i < x.cols(); ++i) {
      EXPECT_NEAR(0.0, kernel.Error(i, Ps[index]), 1e-8);
    }
  }
}

TEST(P3P_Ke_CVPR17, Multiview) {

  const int nViews = 3;
  const int nbPoints = 3;
  const NViewDataSet d =
    NRealisticCamerasRing(nViews, nbPoints,
                          nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K



  // Solve the problem and check that fitted value are good enough
  for (int nResectionCameraIndex = 0; nResectionCameraIndex < nViews; ++nResectionCameraIndex)
  {
    const Mat x = d._x[nResectionCameraIndex];
    const Mat bearing_vectors = (d._K[0].inverse() * x.colwise().homogeneous()).colwise().normalized();
    const Mat X = d._X;
    openMVG::euclidean_resection::PoseResectionKernel_P3P_Ke kernel(bearing_vectors, X);

    std::vector<Mat34> Ps;
    kernel.Fit({0,1,2}, &Ps); // 3 points sample are required, lets take the first three

    bool bFound = false;
    char index = -1;
    for (size_t i = 0; i < Ps.size(); ++i)  {
      Mat34 GT_ProjectionMatrix = d.P(nResectionCameraIndex).array()
      / d.P(nResectionCameraIndex).norm();
      Mat34 COMPUTED_ProjectionMatrix = Ps[i].array() / Ps[i].norm();
      if ( NormLInfinity(GT_ProjectionMatrix - COMPUTED_ProjectionMatrix) < 1e-8 )
      {
        bFound = true;
        index = i;
      }
    }
    EXPECT_TRUE(bFound);

    // Check that for the found matrix the residual is small
    for (Mat::Index i = 0; i < x.cols(); ++i) {
      EXPECT_NEAR(0.0, kernel.Error(i, Ps[index]), 1e-8);
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
