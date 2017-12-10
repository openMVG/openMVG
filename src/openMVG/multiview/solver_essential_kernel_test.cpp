
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
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/solver_essential_eight_point.hpp"
#include "openMVG/multiview/test_data_sets.hpp"

#include "testing/testing.h"

#include <numeric>
using namespace openMVG;

/// Check that the E matrix fit the Essential Matrix properties
/// Determinant is 0
///
#define EXPECT_ESSENTIAL_MATRIX_PROPERTIES(E, expectedPrecision) { \
  EXPECT_NEAR(0, E.determinant(), expectedPrecision); \
  const Mat3 O = 2 * E * E.transpose() * E - (E * E.transpose()).trace() * E; \
  const Mat3 zero3x3 = Mat3::Zero(); \
  EXPECT_MATRIX_NEAR(zero3x3, O, expectedPrecision);\
}

TEST(EightPointsRelativePose, EightPointsRelativePose_Kernel_IdFocal) {

  //-- Setup a circular camera rig and assert that 8PT relative pose solver is working.
  const int iNviews = 5;
  const NViewDataSet d = NRealisticCamerasRing(iNviews, 8,
    nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  for (int i=0; i <iNviews; ++i)
  {
    const int I = i;
    const int J = (i+1)%iNviews;
    std::vector<Mat3> Es; // Essential,
    EightPointRelativePoseSolver::Solve(
      d._K[I].inverse() * d._x[I].colwise().homogeneous(),
      d._K[J].inverse() * d._x[J].colwise().homogeneous(),
      &Es);

    CHECK(!Es.empty());

    // Recover the relative pose from E.
    std::vector<geometry::Pose3> relative_poses;
    for (size_t s = 0; s < Es.size(); ++s) {

      // Check that E holds the essential matrix constraints.
      EXPECT_ESSENTIAL_MATRIX_PROPERTIES(Es[s], 1e-8);
      std::vector<uint32_t> index(d._x[I].cols());
      std::iota(index.begin(), index.end(), 0);
      geometry::Pose3 relative_pose;
      if(RelativePoseFromEssential(
        d._x[I].colwise().homogeneous(),
        d._x[J].colwise().homogeneous(),
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
    RelativeCameraMotion(d._R[I], d._t[I], d._R[J], d._t[J], &R, &t);

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

TEST(FivePointKernelTest, FivePointsRelativePose_Kernel) {

  using Kernel = essential::kernel::FivePointKernel;

  const int focal = 1000;
  const int principal_Point = 500;

  //-- Setup a circular camera rig and assert that 5PT relative pose works.
  const int iNviews = 8;
  const NViewDataSet d = NRealisticCamerasRing(iNviews, Kernel::MINIMUM_SAMPLES,
    nViewDatasetConfigurator(focal,focal,principal_Point,principal_Point,5,0)); // Suppose a camera with Unit matrix as K

  size_t found = 0;
  for (int i=1; i <iNviews; ++i)
  {
    std::vector<Mat3> Es, Rs;  // Essential, Rotation matrix.
    std::vector<Vec3> ts;      // Translation matrix.

    // Direct value do not work.
    // As we use reference, it cannot convert Mat2X& to Mat&
    const Mat x0 = d._x[0];
    const Mat x1 = d._x[i];

    Kernel kernel(x0, x1, d._K[0], d._K[i]);
    std::vector<uint32_t> samples(Kernel::MINIMUM_SAMPLES);
    std::iota(samples.begin(), samples.end(), 0);
    kernel.Fit(samples, &Es);

    // Recover the relative pose from E.
    std::vector<geometry::Pose3> relative_poses;
    for (size_t s = 0; s < Es.size(); ++s) {
      // Check that E holds the essential matrix constraints.
      EXPECT_ESSENTIAL_MATRIX_PROPERTIES(Es[s], 1e-8);
      std::vector<uint32_t> index(d._x[0].cols());
      std::iota(index.begin(), index.end(), 0);
      geometry::Pose3 relative_pose;
      if (RelativePoseFromEssential(
        d._K[0].inverse() * d._x[0].colwise().homogeneous(),
        d._K[i].inverse() * d._x[i].colwise().homogeneous(),
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
    if (bsolution_found)
      found++;
  }
  CHECK_EQUAL(iNviews-1, found);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
