// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/solver_resection_p3p.hpp"
#include "openMVG/multiview/solver_resection_p2pt_fabbri.hpp"
#include "openMVG/multiview/solver_resection_up2p_kukelova.hpp"
#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/system/logger.hpp"

#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::cameras;

#if 0

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


TEST(P3P_Nordberg_ECCV18, Multiview) {

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
    openMVG::euclidean_resection::PoseResectionKernel_P3P_Nordberg kernel(bearing_vectors, X);

    std::vector<Mat34> Ps;
    kernel.Fit({0,1,2}, &Ps); // 3 points sample are required, lets take the first three

    bool bFound = false;
    size_t index = -1;
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

TEST(UP2PSolver_Kukelova, Multiview) {

  const int nViews = 6;
  const int nbPoints = 2;
  const NViewDataSet d =
    NRealisticCamerasRing(nViews, nbPoints,
                          nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  // Solve the problem and check that fitted value are good enough
  for (int nResectionCameraIndex = 0; nResectionCameraIndex < nViews; ++nResectionCameraIndex)
  {
    const Mat x = d._x[nResectionCameraIndex];
    const Mat bearing_vectors = (d._K[0].inverse() * x.colwise().homogeneous()).colwise().normalized();
    const Mat X = d._X;
    openMVG::euclidean_resection::PoseResectionKernel_UP2P_Kukelova kernel(bearing_vectors, X);

    std::vector<Mat34> Ps;
    kernel.Fit({0,1}, &Ps); // 2 points sample are required, lets take the first three

    bool bFound = false;
    size_t index = -1;
      for (size_t i = 0; i < Ps.size(); ++i)  {
      const Mat34 GT_ProjectionMatrix = d.P(nResectionCameraIndex).array()
        / d.P(nResectionCameraIndex).norm();
      const Mat34 COMPUTED_ProjectionMatrix = Ps[i].array() / Ps[i].norm();
      std::cout << "GT:\n " << GT_ProjectionMatrix << std::endl;
      std::cout << "Computed:\n " << COMPUTED_ProjectionMatrix << std::endl;
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
#endif

TEST(P2Pt_Fabbri_ECCV12, Multiview) 
{

  constexpr unsigned nViews = 4;
  constexpr unsigned nbPoints = 5;
  NViewOrientedDataSet d;
  nViewDatasetConfigurator config;
  NOrientedPointsCamerasSphere(nViews, nbPoints, &d, &config); // need to use config since K_ is
                                                               // ignored by getInputScene below
                                                               // (this comes from code
                                                               // initializing K not from a
                                                               // matrix)
  Mat point_tangents_2d, point_tangents_3d;
  point_tangents_2d.resize(6, nbPoints);  // x y 1 tx ty 0
  point_tangents_3d.resize(6, nbPoints);  // X Y Z TX TY TZ

  // Solve the problem and check that fitted value are good enough
  for (int nResectionCameraIndex = 0; nResectionCameraIndex < nViews; ++nResectionCameraIndex) {
    if (nResectionCameraIndex == 1)
      continue;
    OPENMVG_LOG_INFO << "View " << nResectionCameraIndex << "------------------------------------------------------";
  // unsigned nResectionCameraIndex = 0;
    const Mat &x = d._x[nResectionCameraIndex];
    Mat bearing_vectors = (d._K[0].inverse() * x.colwise().homogeneous()).colwise().hnormalized();
    const Mat &tgt = d._tgt2d[nResectionCameraIndex];

    assert (bearing_vectors.cols() == nbPoints);
    bearing_vectors.conservativeResize(3, nbPoints);
    bearing_vectors.row(2).setOnes();


    for (unsigned ip=0; ip < nbPoints; ++ip) {
      point_tangents_2d.col(ip).head(3) = bearing_vectors.col(ip);
      Pinhole_Intrinsic::invert_intrinsics_tgt((double (*)[3]) ((double *)d._K[0].data()), tgt.col(ip).data(), point_tangents_2d.col(ip).data()+3);
      point_tangents_2d.col(ip)(5) = 0;
    }

    const Mat &X = d._X;
    const Mat &T = d._Tgt3d;
    for (unsigned ip=0; ip < nbPoints; ++ip) {
      point_tangents_3d.col(ip).head(3) = X.col(ip);
      point_tangents_3d.col(ip).tail(3) = T.col(ip);
    }

    openMVG::euclidean_resection::PoseResectionKernel_P2Pt_Fabbri kernel(point_tangents_2d, point_tangents_3d);

    std::vector<Mat34> Ps;
    for (unsigned nr=0; nr < 1000; ++nr)
      kernel.Fit({2,3}, &Ps); // 2 points sample are required, lets take the first two
    OPENMVG_LOG_INFO << "Number of returned models: " << Ps.size();


    bool bFound = false;
    size_t index = -1;
    for (size_t i = 0; i < Ps.size(); ++i)  {
      Mat34 GT_ProjectionMatrix = d.Rt(nResectionCameraIndex);
      Mat34 COMPUTED_ProjectionMatrix = Ps[i].array();
//       OPENMVG_LOG_INFO << "gt:\n" << GT_ProjectionMatrix <<  "\ncomputed:\n" << COMPUTED_ProjectionMatrix;
//       OPENMVG_LOG_INFO << "NormLinf:" << NormLInfinity(GT_ProjectionMatrix - COMPUTED_ProjectionMatrix);
      if ( NormLInfinity(GT_ProjectionMatrix - COMPUTED_ProjectionMatrix) < 1e-6 ) {
        bFound = true;
        index = i;
      }
    }
    EXPECT_TRUE(bFound);

    // Check that for the found matrix the residual is small
    for (Mat::Index i = 0; i < x.cols(); ++i)
      EXPECT_NEAR(0.0, kernel.Error(i, Ps[index]), 1e-8);
   }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
