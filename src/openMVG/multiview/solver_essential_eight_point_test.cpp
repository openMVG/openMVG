// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/motion_from_essential.hpp"
#include "openMVG/multiview/projection.hpp"
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

TEST(EightPointsRelativePose, RelativePose_Kernel_IdFocal) {

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

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
