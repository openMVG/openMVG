
// Copyright (c) 2010 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/solver_resection_p3p.hpp"
#include "openMVG/multiview/test_data_sets.hpp"

#include "testing/testing.h"

#include <vector>

using namespace openMVG;

TEST(Resection_Kernel, Multiview) {

  const int nViews = 3;
  const int nbPoints = 10;
  const NViewDataSet d = NRealisticCamerasRing(nViews, nbPoints,
    nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  const int nResectionCameraIndex = 2;

  // Solve the problem and check that fitted value are good enough
  {
    Mat x = d._x[nResectionCameraIndex];
    Mat X = d._X;
    openMVG::resection::kernel::PoseResectionKernel kernel(x, X);

    const std::vector<uint32_t> samples = {0,1,2,3,4,5};
    std::vector<Mat34> Ps;
    kernel.Fit(samples, &Ps);
    for (Mat::Index i = 0; i < x.cols(); ++i) {
      EXPECT_NEAR(0.0, kernel.Error(i, Ps[0]), 1e-8);
    }

    CHECK_EQUAL(1, Ps.size());

    // Check that Projection matrix is near to the GT :
    Mat34 GT_ProjectionMatrix = d.P(nResectionCameraIndex).array()
                                / d.P(nResectionCameraIndex).norm();
    Mat34 COMPUTED_ProjectionMatrix = Ps[0].array() / Ps[0].norm();
    EXPECT_MATRIX_NEAR(GT_ProjectionMatrix, COMPUTED_ProjectionMatrix, 1e-8);
  }
}

TEST(P3P_Kneip_CVPR11, Multiview) {

  const int nViews = 3;
  const int nbPoints = 3;
  const NViewDataSet d = NRealisticCamerasRing(nViews, nbPoints,
    nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  const int nResectionCameraIndex = 2;

  // Solve the problem and check that fitted value are good enough
  {
    const Mat x = d._x[nResectionCameraIndex];
    const Mat X = d._X;
    openMVG::euclidean_resection::P3P_ResectionKernel_K kernel(x, X, d._K[0]);

    const std::vector<uint32_t> samples = {0,1,2};
    std::vector<Mat34> Ps;
    kernel.Fit(samples, &Ps);

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

    // Check that for the found matrix residual is small
    for (Mat::Index i = 0; i < x.cols(); ++i) {
      EXPECT_NEAR(0.0, kernel.Error(i,Ps[index]), 1e-8);
    }
  }
}

// Create a new synthetic dataset for the EPnP implementation.
// It seems it do not perform well on translation like t = (0,0,x)

// Generates all necessary inputs and expected outputs for EuclideanResection.
void CreateCameraSystem(const Mat3& KK,
                        const Mat3X& x_image,
                        const Vec& X_distances,
                        const Mat3& R_input,
                        const Vec3& T_input,
                        Mat2X *x_camera,
                        Mat3X *X_world,
                        Mat3  *R_expected,
                        Vec3  *T_expected) {
  int num_points = x_image.cols();

  Mat3X x_unit_cam(3, num_points);
  x_unit_cam = KK.inverse() * x_image;

  // Create normalized camera coordinates to be used as an input to the PnP
  // function, instead of using NormalizeColumnVectors(&x_unit_cam).
  *x_camera = x_unit_cam.block(0, 0, 2, num_points);
  for (int i = 0; i < num_points; ++i){
    x_unit_cam.col(i).normalize();
  }

  // Create the 3D points in the camera system.
  Mat X_camera(3, num_points);
  for (int i = 0; i < num_points; ++i) {
    X_camera.col(i) = X_distances(i) * x_unit_cam.col(i);
  }

  // Apply the transformation to the camera 3D points
  Mat translation_matrix(3, num_points);
  translation_matrix.row(0).setConstant(T_input(0));
  translation_matrix.row(1).setConstant(T_input(1));
  translation_matrix.row(2).setConstant(T_input(2));

  *X_world = R_input * X_camera + translation_matrix;

  // Create the expected result for comparison.
  *R_expected = R_input.transpose();
  *T_expected = *R_expected * ( - T_input);
};

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
