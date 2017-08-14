// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/linearProgramming/lInfinityCV/resection_kernel.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/robust_estimation/robust_estimator_MaxConsensus.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include <iostream>
#include <vector>

using namespace openMVG;
using namespace openMVG::robust;

TEST(Resection_L_Infinity, Robust_OutlierFree) {

  const int nViews = 3;
  const int nbPoints = 10;
  const NViewDataSet d = NRealisticCamerasRing(nViews, nbPoints,
    nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  //-- Modify a dataset (set to 0 and parse new value) (Assert good values)
  NViewDataSet d2 = d;

  const int nResectionCameraIndex = 2;
  //-- Set to 0 the future computed data to be sure of computation results :
  d2._R[nResectionCameraIndex].fill(0.0);
  d2._t[nResectionCameraIndex].fill(0.0);

  // Solve the problem and check that fitted value are good enough
  {
    using KernelType = lInfinityCV::kernel::l1PoseResectionKernel;
    const Mat & pt2D = d2._x[nResectionCameraIndex];
    const Mat & pt3D = d2._X;
    KernelType kernel(pt2D, pt3D);
    ScorerEvaluator<KernelType> scorer(2*Square(0.6));
    const Mat34 P = MaxConsensus(kernel, scorer, nullptr, 128);

    // Check that Projection matrix is near to the GT :
    const Mat34 GT_ProjectionMatrix = d.P(nResectionCameraIndex).array()
      / d.P(nResectionCameraIndex).norm();
    const Mat34 COMPUTED_ProjectionMatrix = P.array() / P.norm();

    // Extract K[R|t]
    Mat3 R,K;
    Vec3 t;
    KRt_From_P(P, &K, &R, &t);

    d2._R[nResectionCameraIndex] = R;
    d2._t[nResectionCameraIndex] = t;

    //Check matrix to GT, and residual
    EXPECT_NEAR( 0.0, FrobeniusDistance(GT_ProjectionMatrix, COMPUTED_ProjectionMatrix), 1e-2 );
    EXPECT_NEAR( 0.0, RootMeanSquareError(pt2D, pt3D.colwise().homogeneous(), COMPUTED_ProjectionMatrix), 1e-2);
  }
}

TEST(Resection_L_Infinity, Robust_OneOutlier) {

  const int nViews = 3;
  const int nbPoints = 20;
  const NViewDataSet d = NRealisticCamerasRing(nViews, nbPoints,
    nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  d.ExportToPLY("test_Before_Infinity.ply");
  //-- Modify a dataset (set to 0 and parse new value) (Assert good values)
  NViewDataSet d2 = d;

  const int nResectionCameraIndex = 2;
  //-- Set to 0 the future computed data to be sure of computation results :
  d2._R[nResectionCameraIndex].fill(0.0);
  d2._t[nResectionCameraIndex].fill(0.0);

  // Set 20% of the 3D point as outlier
  const int nbOutlier = nbPoints * 0.2;
  for (int i=0; i < nbOutlier; ++i)
  {
    d2._X.col(i)(0) += 120.0;
    d2._X.col(i)(1) += -60.0;
    d2._X.col(i)(2) += 80.0;
  }

  // Solve the problem and check that fitted value are good enough
  {
    using KernelType = lInfinityCV::kernel::l1PoseResectionKernel;
    const Mat & pt2D = d2._x[nResectionCameraIndex];
    const Mat & pt3D = d2._X;
    KernelType kernel(pt2D, pt3D);
    ScorerEvaluator<KernelType> scorer(Square(0.1)); //Highly intolerant for the test
    const Mat34 P = MaxConsensus(kernel, scorer, nullptr, 128);

    // Check that Projection matrix is near to the GT :
    const Mat34 GT_ProjectionMatrix = d.P(nResectionCameraIndex).array()
      / d.P(nResectionCameraIndex).norm();
    const Mat34 COMPUTED_ProjectionMatrix = P.array() / P.norm();

    // Extract K[R|t]
    Mat3 R,K;
    Vec3 t;
    KRt_From_P(P, &K, &R, &t);

    d2._R[nResectionCameraIndex] = R;
    d2._t[nResectionCameraIndex] = t;

    //Check matrix to GT, and residual
    EXPECT_NEAR( 0.0, FrobeniusDistance(GT_ProjectionMatrix, COMPUTED_ProjectionMatrix), 1e-3 );
    // Check that all the inlier have a very small residual reprojection error
    EXPECT_NEAR( 0.0, RootMeanSquareError(
                        pt2D.rightCols(nbPoints - nbOutlier) ,
                        pt3D.colwise().homogeneous().rightCols(nbPoints - nbOutlier),
                        COMPUTED_ProjectionMatrix),
                 1e-1);
  }
  d2.ExportToPLY("test_After_Infinity.ply");
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
