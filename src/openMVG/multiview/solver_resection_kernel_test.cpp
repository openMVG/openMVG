// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/solver_resection_p3p.hpp"
#include "openMVG/multiview/solver_resection_up2p_kukelova.hpp"
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


static void
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

static void
invert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3],
    const double px_tgt_coords[2],
    double normalized_tgt_coords[2])
{
  const double *tp = px_tgt_coords;
  double *t = normalized_tgt_coords;
  t[1] = tp[1]/K[1][1];
  t[0] = (tp[0] - K[0][1]*t[1])/K[0][0];
  double n = hypot(t[0], t[1]);
  t[0] /= n; t[1] /= n;
}

const double K_[2][3] = {
  {2584.9325098195013197, 0, 249.77137587221417903},
  {0, 2584.9325098195013197, 278.31267937919352562}
  //  0 0 1 
};

static constexpr unsigned synth_nviews_ = 1;
static constexpr unsigned synth_npts_ = 2;

// camera format: just like a 3x4 [R|T] but transposed to better fit row-major:
// | R |
// | - |
// | T'|
//
// In this case, this data is from the synthcurves multiview dataset,
// so that instead of T, C is stored:
// | R |
// | - |
// | C'|
const double
cameras_gt_[synth_nviews_][4][3] = {
// extrinsics for frame 42
{
{0.032488343069021832776, 0.14118885304658673752, 0.98944945062394118462},
{-0.28679990702465507635, 0.94965585180689460199, -0.12609352267095536027},
{-0.95743946069465868387, -0.27967744082122686367, 0.071345695037685272211},
{1071.484198049582119, 320.76850873549409471, -85.986368935484179588}
}
};

// Input points and tangents extracted from the synthcurves dataset:
//
// http://github.com/rfabbri/synthcurves-multiview-3d-dataset
//
// sub-dataset:
// spherical-ascii-100_views-perturb-radius_sigma10-normal_sigma0_01rad-minsep_15deg-no_two_cams_colinear_with_object-aspect_ratio1
//      
//
// Frame files: 42, 54, 62, 07
//
// Points and tangents:
// 
// 620
// 1009
// 3011
// 3389
// 0-based ids. +1 to get file line
// 
// Extracting this data from those files: use scripts/getlines from minus
// or just by hand
//
// This is in pixel image coordinates

const double
p_gt_[synth_nviews_][synth_npts_][2] = {
// 2D points for frame 42
{
{181.53712861181011817, 382.84265753788542952},
{361.39404494701216208, 353.17104859076766843},
{342.08123422244102585, 137.87448982117351193},
{243.09524202092040923, 298.35373638828008325},
{285.88966212823157775, 251.48973104783391364}
}
};

// 2D tangents for frame 54
// 2D tangents for frame 62
const double 
t_gt_[synth_nviews_][synth_npts_][2] = {
// 2D tangents for frame 42
{
{0.99164213923671906681, -0.12901886563609127334},
{-0.27407482149753986667, 0.96170837171207557148},
{-0.99085473006575885968, 0.1349329607853933799},
{0.21688457417910370073, -0.97619725541672619507},
{-0.88923826058271226991, 0.45744433094731007383}
};

const double pts3d_gt_[synth_npts_][3] = {
{-39.999960000000001514, 40.000016999999999712, -39.999979999999993652},
{-28.799995999999886465, 40.000010000000003174, 40.000010000000003174},
{16.241229516856364512, -45.185185185185176238, 41.368080573302677294},
{-83.024179089510298013, -7.2456979436932478222, -4.412526863075626693},
{-15.900733289698191442, -6.9926202045388530237, 12.583214033874593696}
};

const double t3d_gt_[synth_npts_][3] = {
{0, 0, 1},
{-1, 0, 0},
{-0.34011103186525448727, -0.10551104075352309153, -0.93444738015720296698},
{-0.71448475613149997621, -0.66437364477423688225, 0.21936087478196092393},
{-0.38185692861567399614, 0.23333310127926898403, -0.89428236587534382096}
};



TEST(P2Pt_Fabbri_ECCV12, Multiview) 
{

  NViewOrientedDataSet d;
  nViewDatasetConfigurator config;
  NOrientedPointsCamerasSphere(&d, &config); // need to use config since K_ is
                                              // ignored by getInputScene below
                                              // (this comes from code
                                              // initializing K not from a
                                              // matrix)
  
  // Solve the problem and check that fitted value are good enough
  for (int nResectionCameraIndex = 0; nResectionCameraIndex < synth_nviews_; ++nResectionCameraIndex) {
    const Mat &x = d._x[nResectionCameraIndex];
    const Mat bearing_vectors = (d._K[0].inverse() * x.colwise().homogeneous()).colwise().hnormalized();
    const Mat &tgt = d._tgt2d[v];

    Matrix 2d_point_tangents;
    2d_point_tangents.resize(6, synth_npts_);
    for (unsigned ip=0; ip < synth_npts_; ++ip) {
      2d_point_tangents.col(ip).head(3) = bearing_vectors.col(ip);
      invert_intrinsics_tgt(d.K[0].data(), tgt.col(ip).data()+3, 2d_point_tangents.col(ip).data()+3);
      2d_point_tangents.col(ip).back() = 0;
    }

    const Mat &X = d._X;
    const Mat &T = d._Tgt3d;
    Mat 3D_point_tangents;
    3D_point_tangents.resize(6, synth_npts_);
    for (unsigned ip=0; ip < synth_npts_; ++ip) {
      3d_point_tangents.col(ip).head(3) = X.col(ip);
      3d_point_tangents.col(ip).tail(3) = T.col(ip);
    }

    openMVG::euclidean_resection::PoseResectionKernel_P2Pt_Fabbri kernel(2d_point_tangents, 3D_point_tangents);

    std::vector<Mat34> Ps;
    kernel.Fit({0,1}, &Ps); // 2 points sample are required, lets take the first two

    bool bFound = false;
    size_t index = -1;
    for (size_t i = 0; i < Ps.size(); ++i)  {
      Mat34 GT_ProjectionMatrix = d.P(nResectionCameraIndex).array() / d.P(nResectionCameraIndex).norm();
      Mat34 COMPUTED_ProjectionMatrix = Ps[i].array() / Ps[i].norm();
      if ( NormLInfinity(GT_ProjectionMatrix - COMPUTED_ProjectionMatrix) < 1e-8 ) {
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

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
