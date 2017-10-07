
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

#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/motion_from_essential.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/solver_essential_five_point.hpp"
#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/numeric/numeric.h"

#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/multiview/triangulation.hpp"

#include "testing/testing.h"

#include <iostream>
#include <random>

using namespace openMVG;
using namespace openMVG::sfm;

struct TestData {
  //-- Dataset that encapsulate :
  // 3D points and their projection given P1 and P2
  // Link between the two camera [R|t]
  Mat3X X;
  Mat3 R;
  Vec3 t;
  Mat3 E;
  Mat34 P1, P2;
  Mat2X x1, x2;
};

TestData SomeTestData() {
  TestData d;

  // --
  // Modeling random 3D points,
  // Consider first camera as [R=I|t=0],
  // Second camera as [R=Rx*Ry*Rz|t=random],
  // Compute projection of the 3D points onto image plane.
  // --
  d.X = Mat3X::Random(3,5);

  //-- Make point in front to the cameras.
  d.X.row(0).array() -= .5;
  d.X.row(1).array() -= .5;
  d.X.row(2).array() += 1.0;

  d.R = RotationAroundZ(0.3) * RotationAroundX(0.1) * RotationAroundY(0.2);
  d.t.setRandom();

  EssentialFromRt(Mat3::Identity(), Vec3::Zero(), d.R, d.t, &d.E);

  P_From_KRt(Mat3::Identity(), Mat3::Identity(), Vec3::Zero(), &d.P1);
  P_From_KRt(Mat3::Identity(), d.R, d.t, &d.P2);

  Project(d.P1, d.X, &d.x1);
  Project(d.P2, d.X, &d.x2);

  return d;
}

TEST(FivePointsNullspaceBasis, SatisfyEpipolarConstraint) {

  const TestData d = SomeTestData();

  const Mat E_basis = FivePointsNullspaceBasis(d.x1.colwise().homogeneous(),
                                               d.x2.colwise().homogeneous());

  for (const int s : {0, 1, 2, 3}) {
    Mat3 E;
    for (const int i : {0, 1, 2}) {
      for (const int j : {0, 1, 2}) {
        E(i, j) = E_basis(3 * i + j, s);
      }
    }
    for (int i = 0; i < d.x1.cols(); ++i) {
      const Vec3 x1(d.x1(0,i), d.x1(1,i), 1);
      const Vec3 x2(d.x2(0,i), d.x2(1,i), 1);
      EXPECT_NEAR(0, x2.dot(E * x1), 1e-6);
    }
  }
}

double EvalPolynomial(Vec p, double x, double y, double z) {
  return p(coef_xxx) * x * x * x
       + p(coef_xxy) * x * x * y
       + p(coef_xxz) * x * x * z
       + p(coef_xyy) * x * y * y
       + p(coef_xyz) * x * y * z
       + p(coef_xzz) * x * z * z
       + p(coef_yyy) * y * y * y
       + p(coef_yyz) * y * y * z
       + p(coef_yzz) * y * z * z
       + p(coef_zzz) * z * z * z
       + p(coef_xx)  * x * x
       + p(coef_xy)  * x * y
       + p(coef_xz)  * x * z
       + p(coef_yy)  * y * y
       + p(coef_yz)  * y * z
       + p(coef_zz)  * z * z
       + p(coef_x)   * x
       + p(coef_y)   * y
       + p(coef_z)   * z
       + p(coef_1)   * 1;
}

TEST(o1, Evaluation) {

  Vec p1 = Vec::Zero(20), p2 = Vec::Zero(20);
  std::mt19937 gen(std::mt19937::default_seed);
  std::uniform_real_distribution<double> dis(0, 1);

  p1(coef_x) = dis(gen);
  p1(coef_y) = dis(gen);
  p1(coef_z) = dis(gen);
  p1(coef_1) = dis(gen);
  p2(coef_x) = dis(gen);
  p2(coef_y) = dis(gen);
  p2(coef_z) = dis(gen);
  p2(coef_1) = dis(gen);

  Vec p3 = o1(p1, p2);

  for (double z = -5; z < 5; ++z) {
    for (double y = -5; y < 5; ++y) {
      for (double x = -5; x < 5; ++x) {
        EXPECT_NEAR(EvalPolynomial(p3, x, y, z),
                    EvalPolynomial(p1, x, y, z) * EvalPolynomial(p2, x, y, z),
                    1e-8);
      }
    }
  }
}

TEST(o2, Evaluation) {

  Vec p1 = Vec::Zero(20), p2 = Vec::Zero(20);
  std::mt19937 gen(std::mt19937::default_seed);
  std::uniform_real_distribution<double> dis(0, 1);
  p1(coef_xx) = dis(gen);
  p1(coef_xy) = dis(gen);
  p1(coef_xz) = dis(gen);
  p1(coef_yy) = dis(gen);
  p1(coef_yz) = dis(gen);
  p1(coef_zz) = dis(gen);
  p1(coef_x)  = dis(gen);
  p1(coef_y)  = dis(gen);
  p1(coef_z)  = dis(gen);
  p1(coef_1)  = dis(gen);
  p2(coef_x)  = dis(gen);
  p2(coef_y)  = dis(gen);
  p2(coef_z)  = dis(gen);
  p2(coef_1)  = dis(gen);

  Vec p3 = o2(p1, p2);

  for (double z = -5; z < 5; ++z) {
    for (double y = -5; y < 5; ++y) {
      for (double x = -5; x < 5; ++x) {
        EXPECT_NEAR(EvalPolynomial(p3, x, y, z),
                    EvalPolynomial(p1, x, y, z) * EvalPolynomial(p2, x, y, z),
                    1e-8);
      }
    }
  }
}

/// Check that the E matrix fit the Essential Matrix properties
/// Determinant is 0
///
#define EXPECT_ESSENTIAL_MATRIX_PROPERTIES(E, expectedPrecision) { \
  EXPECT_NEAR(0, E.determinant(), expectedPrecision); \
  const Mat3 O = 2 * E * E.transpose() * E - (E * E.transpose()).trace() * E; \
  const Mat3 zero3x3 = Mat3::Zero(); \
  EXPECT_MATRIX_NEAR(zero3x3, O, expectedPrecision);\
}

TEST(FivePointsRelativePose, Random) {

  const TestData d = SomeTestData();

  std::vector<Mat3> Es;
  FivePointsRelativePose(d.x1.colwise().homogeneous(),
                         d.x2.colwise().homogeneous(),
                         &Es);

  // Recover the relative pose from E.
  std::vector<geometry::Pose3> relative_poses;
  for (size_t s = 0; s < Es.size(); ++s) {

    // Check that E holds the essential matrix constraints.
    EXPECT_ESSENTIAL_MATRIX_PROPERTIES(Es[s], 1e-8);

    std::vector<uint32_t> index(d.x1.cols());
    std::iota(index.begin(), index.end(), 0);
    geometry::Pose3 relative_pose;
    if(RelativePoseFromEssential(
      d.x1.colwise().homogeneous(),
      d.x2.colwise().homogeneous(),
      Es[s],
      index,
      &relative_pose))
    {
      relative_poses.push_back(relative_pose);
    }
  }
  CHECK(!relative_poses.empty());

  bool bsolution_found = false;
  for (size_t i = 0; i < relative_poses.size(); ++i) {

    // Check that we find the correct relative orientation.
    if (FrobeniusDistance(d.R, relative_poses[i].rotation()) < 1e-3
      && (d.t / d.t.norm()
          - relative_poses[i].translation()
            / relative_poses[i].translation().norm()).norm() < 1e-3 ) {
        bsolution_found = true;
    }
  }
  //-- Almost one solution must find the correct relative orientation
  CHECK(bsolution_found);
}


TEST(FivePointsRelativePose, test_data_sets) {

  //-- Setup a circular camera rig and assert that 5PT relative pose works.
  const int iNviews = 5;
  const NViewDataSet d = NRealisticCamerasRing(iNviews, 5,
    nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  // Compute pose [R|t] from 0 to [1;..;iNviews]
  for (int i=1; i < iNviews; ++i)
  {
    std::vector<Mat3> Es;  // Essential matrices.
    const Mat3X bearing0 = d._K[0].inverse() * d._x[0].colwise().homogeneous();
    const Mat3X bearingI = d._K[i].inverse() * d._x[i].colwise().homogeneous();
    FivePointsRelativePose(bearing0,
                           bearingI,
                           &Es);
    CHECK(!Es.empty());

    // Recover the relative pose from E.
    std::vector<geometry::Pose3> relative_poses;
    for (size_t s = 0; s < Es.size(); ++s) {

      // Check that E holds the essential matrix constraints.
      EXPECT_ESSENTIAL_MATRIX_PROPERTIES(Es[s], 1e-8);
      std::vector<uint32_t> index(d._x[0].cols());
      std::iota(index.begin(), index.end(), 0);
      geometry::Pose3 relative_pose;
      if (RelativePoseFromEssential(
        bearing0,
        bearingI,
        Es[s],
        index,
        &relative_pose))
      {
        relative_poses.push_back(relative_pose);
      }
    }

    CHECK(!relative_poses.empty());
    //-- Compute Ground Truth motion
    Mat3 R;
    Vec3 t;
    RelativeCameraMotion(d._R[0], d._t[0], d._R[i], d._t[i], &R, &t);

    // Assert that found relative motion is correct for almost one model.
    bool bsolution_found = false;
    for (size_t nModel = 0; nModel < relative_poses.size(); ++nModel) {

      // Check that we find the correct relative orientation.
      if (FrobeniusDistance(R, relative_poses[nModel].rotation()) < 1e-3
        && (t / t.norm()
            - relative_poses[nModel].translation()
              / relative_poses[nModel].translation().norm()).norm() < 1e-3 ) {
          bsolution_found = true;
      }
    }
    //-- Almost one solution must find the correct relative orientation
    CHECK(bsolution_found);
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
