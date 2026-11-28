// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp"
#include "testing/testing.h"

using namespace openMVG;

TEST(SphericalBundleResidual, LongitudeSeamValueAndDerivative)
{
  const cameras::Intrinsic_Spherical camera(2000, 1000);
  const Vec2 observed(1, 503);
  const Vec3 point = camera(Vec2(1999, 500));
  sfm::ResidualErrorFunctor_Intrinsic_Spherical functor(observed.data(), 2000, 1000);

  using Jet = ceres::Jet<double, 1>;
  Jet pose[6] = {};
  Jet position[3] = {Jet(point.x(), 0), Jet(point.y()), Jet(point.z())};
  Jet residual[2];
  EXPECT_TRUE(functor(pose, position, residual));
  EXPECT_NEAR(-2, residual[0].a, 1e-10);
  EXPECT_NEAR(-3, residual[1].a, 1e-10);
  EXPECT_NEAR(2000 / (2 * M_PI) * point.z() /
    (point.x() * point.x() + point.z() * point.z()), residual[0].v[0], 1e-10);

  const Vec2 opposite_observed(1999, 500);
  const Vec3 opposite_point = camera(Vec2(1, 500));
  const double pose_double[6] = {};
  double opposite_residual[2];
  sfm::ResidualErrorFunctor_Intrinsic_Spherical opposite(
    opposite_observed.data(), 2000, 1000);
  EXPECT_TRUE(opposite(pose_double, opposite_point.data(), opposite_residual));
  EXPECT_NEAR(2, opposite_residual[0], 1e-10);
  EXPECT_NEAR(0, opposite_residual[1], 1e-10);
}

int main() { TestResult tr; return TestRegistry::runAllTests(tr); }
