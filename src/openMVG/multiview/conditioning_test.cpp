// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/conditioning.hpp"
#include "testing/testing.h"

#include <cmath>

using namespace openMVG;

TEST(Conditioning, CentersConstantAxis) {
  Mat points(2, 3);
  points << 5, 5, 5, -2, 0, 2;
  Mat normalized;
  Mat3 transform;
  NormalizePoints(points, &normalized, &transform);
  EXPECT_TRUE(normalized.allFinite());
  EXPECT_NEAR(0.0, normalized.row(0).norm(), 1e-12);
  EXPECT_NEAR(0.0, normalized.row(1).mean(), 1e-12);
  EXPECT_NEAR(2.0, normalized.row(1).squaredNorm() / 3, 1e-12);
  EXPECT_EQ(1.0, transform(0, 0));
  EXPECT_EQ(-5.0, transform(0, 2));
}

TEST(Conditioning, CentersCoincidentPoints) {
  Mat points(2, 3);
  points << 5, 5, 5, -7, -7, -7;
  Mat normalized;
  Mat3 transform;
  NormalizePoints(points, &normalized, &transform);
  EXPECT_TRUE(normalized.allFinite());
  EXPECT_NEAR(0.0, normalized.norm(), 1e-12);
  Mat restored;
  ApplyTransformationToPoints(normalized, transform.inverse(), &restored);
  EXPECT_MATRIX_NEAR(points, restored, 1e-12);
}

TEST(Conditioning, CentersNearConstantAxis) {
  Mat points(2, 3);
  points << 0, 1, 2, 2 - 1e-5, 2, 2 + 1e-5;
  Mat normalized;
  Mat3 transform;
  NormalizePoints(points, &normalized, &transform);
  EXPECT_TRUE(normalized.allFinite());
  EXPECT_NEAR(0.0, normalized.row(1).mean(), 1e-12);
  EXPECT_EQ(1.0, transform(1, 1));
  EXPECT_NEAR(-1e-5, normalized(1, 0), 1e-12);
  EXPECT_NEAR(1e-5, normalized(1, 2), 1e-12);
}

TEST(Conditioning, NormalizesImageRanges) {
  for (const auto & dimensions : {Vec2(640, 480), Vec2(60000, 40000)}) {
    const int width = dimensions(0), height = dimensions(1);
    Mat points(2, 3);
    points << 0, width / 2.0, width, 0, height / 2.0, height;
    Mat normalized;
    Mat3 transform;
    NormalizePoints(points, &normalized, &transform, width, height);
    EXPECT_TRUE(transform.allFinite());
    EXPECT_TRUE(normalized.allFinite());
    const double scale = 1.0 / std::sqrt(double(width) * double(height));
    EXPECT_NEAR(scale, transform(0, 0), 1e-15);
    EXPECT_NEAR(scale, transform(1, 1), 1e-15);
    EXPECT_NEAR(0.0, normalized.col(1).norm(), 1e-12);
    EXPECT_NEAR(-0.5 * width * scale, normalized(0, 0), 1e-12);
    EXPECT_NEAR(0.5 * height * scale, normalized(1, 2), 1e-12);
  }
}

int main() {
  TestResult result;
  return TestRegistry::runAllTests(result);
}
