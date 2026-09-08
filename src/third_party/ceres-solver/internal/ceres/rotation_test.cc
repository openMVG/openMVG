// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2023 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: sameeragarwal@google.com (Sameer Agarwal)

#include "ceres/rotation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <utility>

#include "ceres/constants.h"
#include "ceres/internal/eigen.h"
#include "ceres/internal/euler_angles.h"
#include "ceres/internal/export.h"
#include "ceres/is_close.h"
#include "ceres/jet.h"
#include "ceres/stringprintf.h"
#include "ceres/test_util.h"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

inline constexpr double kPi = constants::pi;
const double kHalfSqrt2 = 0.707106781186547524401;

// A tolerance value for floating-point comparisons.
static double const kTolerance = std::numeric_limits<double>::epsilon() * 10;

// Looser tolerance used for numerically unstable conversions.
static double const kLooseTolerance = 1e-9;

// Use as:
// double quaternion[4];
// EXPECT_THAT(quaternion, IsNormalizedQuaternion());
MATCHER(IsNormalizedQuaternion, "") {
  double norm2 =
      arg[0] * arg[0] + arg[1] * arg[1] + arg[2] * arg[2] + arg[3] * arg[3];
  if (fabs(norm2 - 1.0) > kTolerance) {
    *result_listener << "squared norm is " << norm2;
    return false;
  }

  return true;
}

// Use as:
// double expected_quaternion[4];
// double actual_quaternion[4];
// EXPECT_THAT(actual_quaternion, IsNearQuaternion(expected_quaternion));
MATCHER_P(IsNearQuaternion, expected, "") {
  // Quaternions are equivalent upto a sign change. So we will compare
  // both signs before declaring failure.
  bool is_near = true;
  // NOTE: near (and far) can be defined as macros on the Windows platform (for
  // ancient pascal calling convention). Do not use these identifiers.
  for (int i = 0; i < 4; i++) {
    if (fabs(arg[i] - expected[i]) > kTolerance) {
      is_near = false;
      break;
    }
  }

  if (is_near) {
    return true;
  }

  is_near = true;
  for (int i = 0; i < 4; i++) {
    if (fabs(arg[i] + expected[i]) > kTolerance) {
      is_near = false;
      break;
    }
  }

  if (is_near) {
    return true;
  }

  // clang-format off
  *result_listener << "expected : "
                   << expected[0] << " "
                   << expected[1] << " "
                   << expected[2] << " "
                   << expected[3] << " "
                   << "actual : "
                   << arg[0] << " "
                   << arg[1] << " "
                   << arg[2] << " "
                   << arg[3];
  // clang-format on
  return false;
}

// Use as:
// double expected_axis_angle[3];
// double actual_axis_angle[3];
// EXPECT_THAT(actual_axis_angle, IsNearAngleAxis(expected_axis_angle));
MATCHER_P(IsNearAngleAxis, expected, "") {
  Eigen::Vector3d a(arg[0], arg[1], arg[2]);
  Eigen::Vector3d e(expected[0], expected[1], expected[2]);
  const double e_norm = e.norm();

  double delta_norm = std::numeric_limits<double>::max();
  if (e_norm > 0) {
    // Deal with the sign ambiguity near PI. Since the sign can flip,
    // we take the smaller of the two differences.
    if (fabs(e_norm - kPi) < kLooseTolerance) {
      delta_norm = std::min((a - e).norm(), (a + e).norm()) / e_norm;
    } else {
      delta_norm = (a - e).norm() / e_norm;
    }
  } else {
    delta_norm = a.norm();
  }

  if (delta_norm <= kLooseTolerance) {
    return true;
  }

  // clang-format off
  *result_listener << " arg:"
                   << " " << arg[0]
                   << " " << arg[1]
                   << " " << arg[2]
                   << " was expected to be:"
                   << " " << expected[0]
                   << " " << expected[1]
                   << " " << expected[2];
  // clang-format on
  return false;
}

// Use as:
// double matrix[9];
// EXPECT_THAT(matrix, IsOrthonormal());
MATCHER(IsOrthonormal, "") {
  for (int c1 = 0; c1 < 3; c1++) {
    for (int c2 = 0; c2 < 3; c2++) {
      double v = 0;
      for (int i = 0; i < 3; i++) {
        v += arg[i + 3 * c1] * arg[i + 3 * c2];
      }
      double expected = (c1 == c2) ? 1 : 0;
      if (fabs(expected - v) > kTolerance) {
        *result_listener << "Columns " << c1 << " and " << c2
                         << " should have dot product " << expected
                         << " but have " << v;
        return false;
      }
    }
  }

  return true;
}

// Use as:
// double matrix1[9];
// double matrix2[9];
// EXPECT_THAT(matrix1, IsNear3x3Matrix(matrix2));
MATCHER_P(IsNear3x3Matrix, expected, "") {
  for (int i = 0; i < 9; i++) {
    if (fabs(arg[i] - expected[i]) > kTolerance) {
      *result_listener << "component " << i << " should be " << expected[i];
      return false;
    }
  }

  return true;
}

// Transforms a zero axis/angle to a quaternion.
TEST(Rotation, ZeroAngleAxisToQuaternion) {
  double axis_angle[3] = {0, 0, 0};
  double quaternion[4];
  double expected[4] = {1, 0, 0, 0};
  AngleAxisToQuaternion(axis_angle, quaternion);
  EXPECT_THAT(quaternion, IsNormalizedQuaternion());
  EXPECT_THAT(quaternion, IsNearQuaternion(expected));
}

// Test that exact conversion works for small angles.
TEST(Rotation, SmallAngleAxisToQuaternion) {
  // Small, finite value to test.
  double theta = 1.0e-2;
  double axis_angle[3] = {theta, 0, 0};
  double quaternion[4];
  double expected[4] = {cos(theta / 2), sin(theta / 2.0), 0, 0};
  AngleAxisToQuaternion(axis_angle, quaternion);
  EXPECT_THAT(quaternion, IsNormalizedQuaternion());
  EXPECT_THAT(quaternion, IsNearQuaternion(expected));
}

// Test that approximate conversion works for very small angles.
TEST(Rotation, TinyAngleAxisToQuaternion) {
  // Very small value that could potentially cause underflow.
  double theta = pow(std::numeric_limits<double>::min(), 0.75);
  double axis_angle[3] = {theta, 0, 0};
  double quaternion[4];
  double expected[4] = {cos(theta / 2), sin(theta / 2.0), 0, 0};
  AngleAxisToQuaternion(axis_angle, quaternion);
  EXPECT_THAT(quaternion, IsNormalizedQuaternion());
  EXPECT_THAT(quaternion, IsNearQuaternion(expected));
}

// Transforms a rotation by pi/2 around X to a quaternion.
TEST(Rotation, XRotationToQuaternion) {
  double axis_angle[3] = {kPi / 2, 0, 0};
  double quaternion[4];
  double expected[4] = {kHalfSqrt2, kHalfSqrt2, 0, 0};
  AngleAxisToQuaternion(axis_angle, quaternion);
  EXPECT_THAT(quaternion, IsNormalizedQuaternion());
  EXPECT_THAT(quaternion, IsNearQuaternion(expected));
}

// Transforms a unit quaternion to an axis angle.
TEST(Rotation, UnitQuaternionToAngleAxis) {
  double quaternion[4] = {1, 0, 0, 0};
  double axis_angle[3];
  double expected[3] = {0, 0, 0};
  QuaternionToAngleAxis(quaternion, axis_angle);
  EXPECT_THAT(axis_angle, IsNearAngleAxis(expected));
}

// Transforms a quaternion that rotates by pi about the Y axis to an axis angle.
TEST(Rotation, YRotationQuaternionToAngleAxis) {
  double quaternion[4] = {0, 0, 1, 0};
  double axis_angle[3];
  double expected[3] = {0, kPi, 0};
  QuaternionToAngleAxis(quaternion, axis_angle);
  EXPECT_THAT(axis_angle, IsNearAngleAxis(expected));
}

// Transforms a quaternion that rotates by pi/3 about the Z axis to an axis
// angle.
TEST(Rotation, ZRotationQuaternionToAngleAxis) {
  double quaternion[4] = {sqrt(3) / 2, 0, 0, 0.5};
  double axis_angle[3];
  double expected[3] = {0, 0, kPi / 3};
  QuaternionToAngleAxis(quaternion, axis_angle);
  EXPECT_THAT(axis_angle, IsNearAngleAxis(expected));
}

// Test that exact conversion works for small angles.
TEST(Rotation, SmallQuaternionToAngleAxis) {
  // Small, finite value to test.
  double theta = 1.0e-2;
  double quaternion[4] = {cos(theta / 2), sin(theta / 2.0), 0, 0};
  double axis_angle[3];
  double expected[3] = {theta, 0, 0};
  QuaternionToAngleAxis(quaternion, axis_angle);
  EXPECT_THAT(axis_angle, IsNearAngleAxis(expected));
}

// Test that approximate conversion works for very small angles.
TEST(Rotation, TinyQuaternionToAngleAxis) {
  // Very small value that could potentially cause underflow.
  double theta = pow(std::numeric_limits<double>::min(), 0.75);
  double quaternion[4] = {cos(theta / 2), sin(theta / 2.0), 0, 0};
  double axis_angle[3];
  double expected[3] = {theta, 0, 0};
  QuaternionToAngleAxis(quaternion, axis_angle);
  EXPECT_THAT(axis_angle, IsNearAngleAxis(expected));
}

TEST(Rotation, QuaternionToAngleAxisAngleIsLessThanPi) {
  double quaternion[4];
  double angle_axis[3];

  const double half_theta = 0.75 * kPi;

  quaternion[0] = cos(half_theta);
  quaternion[1] = 1.0 * sin(half_theta);
  quaternion[2] = 0.0;
  quaternion[3] = 0.0;
  QuaternionToAngleAxis(quaternion, angle_axis);
  const double angle = std::hypot(angle_axis[0], angle_axis[1], angle_axis[2]);
  EXPECT_LE(angle, kPi);
}

static constexpr int kNumTrials = 10000;

// Takes a bunch of random axis/angle values, converts them to quaternions,
// and back again.
TEST(Rotation, AngleAxisToQuaterionAndBack) {
  std::mt19937 prng;
  std::uniform_real_distribution<double> uniform_distribution{-1.0, 1.0};
  for (int i = 0; i < kNumTrials; i++) {
    double axis_angle[3];
    // Make an axis by choosing three random numbers in [-1, 1) and
    // normalizing.
    double norm = 0;
    for (double& coeff : axis_angle) {
      coeff = uniform_distribution(prng);
      norm += coeff * coeff;
    }
    norm = sqrt(norm);

    // Angle in [-pi, pi).
    double theta = uniform_distribution(
        prng, std::uniform_real_distribution<double>::param_type{-kPi, kPi});
    for (double& coeff : axis_angle) {
      coeff = coeff * theta / norm;
    }

    double quaternion[4];
    double round_trip[3];
    // We use ASSERTs here because if there's one failure, there are
    // probably many and spewing a million failures doesn't make anyone's
    // day.
    AngleAxisToQuaternion(axis_angle, quaternion);
    ASSERT_THAT(quaternion, IsNormalizedQuaternion());
    QuaternionToAngleAxis(quaternion, round_trip);
    ASSERT_THAT(round_trip, IsNearAngleAxis(axis_angle));
  }
}

// Takes a bunch of random quaternions, converts them to axis/angle,
// and back again.
TEST(Rotation, QuaterionToAngleAxisAndBack) {
  std::mt19937 prng;
  std::uniform_real_distribution<double> uniform_distribution{-1.0, 1.0};
  for (int i = 0; i < kNumTrials; i++) {
    double quaternion[4];
    // Choose four random numbers in [-1, 1) and normalize.
    double norm = 0;
    for (double& coeff : quaternion) {
      coeff = uniform_distribution(prng);
      norm += coeff * coeff;
    }
    norm = sqrt(norm);

    for (double& coeff : quaternion) {
      coeff = coeff / norm;
    }

    double axis_angle[3];
    double round_trip[4];
    QuaternionToAngleAxis(quaternion, axis_angle);
    AngleAxisToQuaternion(axis_angle, round_trip);
    ASSERT_THAT(round_trip, IsNormalizedQuaternion());
    ASSERT_THAT(round_trip, IsNearQuaternion(quaternion));
  }
}

// Transforms a zero axis/angle to a rotation matrix.
TEST(Rotation, ZeroAngleAxisToRotationMatrix) {
  double axis_angle[3] = {0, 0, 0};
  double matrix[9];
  double expected[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  AngleAxisToRotationMatrix(axis_angle, matrix);
  EXPECT_THAT(matrix, IsOrthonormal());
  EXPECT_THAT(matrix, IsNear3x3Matrix(expected));
}

TEST(Rotation, NearZeroAngleAxisToRotationMatrix) {
  double axis_angle[3] = {1e-24, 2e-24, 3e-24};
  double matrix[9];
  double expected[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  AngleAxisToRotationMatrix(axis_angle, matrix);
  EXPECT_THAT(matrix, IsOrthonormal());
  EXPECT_THAT(matrix, IsNear3x3Matrix(expected));
}

// Transforms a rotation by pi/2 around X to a rotation matrix and back.
TEST(Rotation, XRotationToRotationMatrix) {
  double axis_angle[3] = {kPi / 2, 0, 0};
  double matrix[9];
  // The rotation matrices are stored column-major.
  double expected[9] = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  AngleAxisToRotationMatrix(axis_angle, matrix);
  EXPECT_THAT(matrix, IsOrthonormal());
  EXPECT_THAT(matrix, IsNear3x3Matrix(expected));
  double round_trip[3];
  RotationMatrixToAngleAxis(matrix, round_trip);
  EXPECT_THAT(round_trip, IsNearAngleAxis(axis_angle));
}

// Transforms an axis angle that rotates by pi about the Y axis to a
// rotation matrix and back.
TEST(Rotation, YRotationToRotationMatrix) {
  double axis_angle[3] = {0, kPi, 0};
  double matrix[9];
  double expected[9] = {-1, 0, 0, 0, 1, 0, 0, 0, -1};
  AngleAxisToRotationMatrix(axis_angle, matrix);
  EXPECT_THAT(matrix, IsOrthonormal());
  EXPECT_THAT(matrix, IsNear3x3Matrix(expected));

  double round_trip[3];
  RotationMatrixToAngleAxis(matrix, round_trip);
  EXPECT_THAT(round_trip, IsNearAngleAxis(axis_angle));
}

TEST(Rotation, NearPiAngleAxisRoundTrip) {
  double in_axis_angle[3];
  double matrix[9];
  double out_axis_angle[3];

  std::mt19937 prng;
  std::uniform_real_distribution<double> uniform_distribution{-1.0, 1.0};
  for (int i = 0; i < kNumTrials; i++) {
    // Make an axis by choosing three random numbers in [-1, 1) and
    // normalizing.
    double norm = 0;
    for (double& coeff : in_axis_angle) {
      coeff = uniform_distribution(prng);
      norm += coeff * coeff;
    }
    norm = sqrt(norm);

    // Angle in [pi - kMaxSmallAngle, pi).
    constexpr double kMaxSmallAngle = 1e-8;
    double theta =
        uniform_distribution(prng,
                             std::uniform_real_distribution<double>::param_type{
                                 kPi - kMaxSmallAngle, kPi});

    for (double& coeff : in_axis_angle) {
      coeff *= (theta / norm);
    }
    AngleAxisToRotationMatrix(in_axis_angle, matrix);
    RotationMatrixToAngleAxis(matrix, out_axis_angle);
    EXPECT_THAT(in_axis_angle, IsNearAngleAxis(out_axis_angle));
  }
}

TEST(Rotation, AtPiAngleAxisRoundTrip) {
  // A rotation of kPi about the X axis;
  // clang-format off
  static constexpr double kMatrix[3][3] = {
    {1.0,  0.0,  0.0},
    {0.0,  -1.0,  0.0},
    {0.0,  0.0,  -1.0}
  };
  // clang-format on

  double in_matrix[9];
  // Fill it from kMatrix in col-major order.
  for (int j = 0, k = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i, ++k) {
      in_matrix[k] = kMatrix[i][j];
    }
  }

  const double expected_axis_angle[3] = {kPi, 0, 0};

  double out_matrix[9];
  double axis_angle[3];
  RotationMatrixToAngleAxis(in_matrix, axis_angle);
  AngleAxisToRotationMatrix(axis_angle, out_matrix);

  LOG(INFO) << "AngleAxis = " << axis_angle[0] << " " << axis_angle[1] << " "
            << axis_angle[2];
  LOG(INFO) << "Expected AngleAxis = " << kPi << " 0 0";
  double out_rowmajor[3][3];
  for (int j = 0, k = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i, ++k) {
      out_rowmajor[i][j] = out_matrix[k];
    }
  }
  LOG(INFO) << "Rotation:";
  LOG(INFO) << "EXPECTED        |        ACTUAL";
  for (int i = 0; i < 3; ++i) {
    std::string line;
    for (int j = 0; j < 3; ++j) {
      StringAppendF(&line, "%g ", kMatrix[i][j]);
    }
    line += "         |        ";
    for (int j = 0; j < 3; ++j) {
      StringAppendF(&line, "%g ", out_rowmajor[i][j]);
    }
    LOG(INFO) << line;
  }

  EXPECT_THAT(axis_angle, IsNearAngleAxis(expected_axis_angle));
  EXPECT_THAT(out_matrix, IsNear3x3Matrix(in_matrix));
}

// Transforms an axis angle that rotates by pi/3 about the Z axis to a
// rotation matrix.
TEST(Rotation, ZRotationToRotationMatrix) {
  double axis_angle[3] = {0, 0, kPi / 3};
  double matrix[9];
  // This is laid-out row-major on the screen but is actually stored
  // column-major.
  // clang-format off
  double expected[9] = { 0.5, sqrt(3) / 2, 0,   // Column 1
                         -sqrt(3) / 2, 0.5, 0,  // Column 2
                         0, 0, 1 };             // Column 3
  // clang-format on
  AngleAxisToRotationMatrix(axis_angle, matrix);
  EXPECT_THAT(matrix, IsOrthonormal());
  EXPECT_THAT(matrix, IsNear3x3Matrix(expected));
  double round_trip[3];
  RotationMatrixToAngleAxis(matrix, round_trip);
  EXPECT_THAT(round_trip, IsNearAngleAxis(axis_angle));
}

// Takes a bunch of random axis/angle values, converts them to rotation
// matrices, and back again.
TEST(Rotation, AngleAxisToRotationMatrixAndBack) {
  std::mt19937 prng;
  std::uniform_real_distribution<double> uniform_distribution{-1.0, 1.0};
  for (int i = 0; i < kNumTrials; i++) {
    double axis_angle[3];
    // Make an axis by choosing three random numbers in [-1, 1) and
    // normalizing.
    double norm = 0;
    for (double& i : axis_angle) {
      i = uniform_distribution(prng);
      norm += i * i;
    }
    norm = sqrt(norm);

    // Angle in [-pi, pi).
    double theta = uniform_distribution(
        prng, std::uniform_real_distribution<double>::param_type{-kPi, kPi});
    for (double& i : axis_angle) {
      i = i * theta / norm;
    }

    double matrix[9];
    double round_trip[3];
    AngleAxisToRotationMatrix(axis_angle, matrix);
    ASSERT_THAT(matrix, IsOrthonormal());
    RotationMatrixToAngleAxis(matrix, round_trip);

    for (int i = 0; i < 3; ++i) {
      EXPECT_NEAR(round_trip[i], axis_angle[i], kLooseTolerance);
    }
  }
}

// Takes a bunch of random axis/angle values near zero, converts them
// to rotation matrices, and back again.
TEST(Rotation, AngleAxisToRotationMatrixAndBackNearZero) {
  std::mt19937 prng;
  std::uniform_real_distribution<double> uniform_distribution{-1.0, 1.0};
  for (int i = 0; i < kNumTrials; i++) {
    double axis_angle[3];
    // Make an axis by choosing three random numbers in [-1, 1) and
    // normalizing.
    double norm = 0;
    for (double& i : axis_angle) {
      i = uniform_distribution(prng);
      norm += i * i;
    }
    norm = sqrt(norm);

    // Tiny theta.
    constexpr double kScale = 1e-16;
    double theta =
        uniform_distribution(prng,
                             std::uniform_real_distribution<double>::param_type{
                                 -kScale * kPi, kScale * kPi});
    for (double& i : axis_angle) {
      i = i * theta / norm;
    }

    double matrix[9];
    double round_trip[3];
    AngleAxisToRotationMatrix(axis_angle, matrix);
    ASSERT_THAT(matrix, IsOrthonormal());
    RotationMatrixToAngleAxis(matrix, round_trip);

    for (int i = 0; i < 3; ++i) {
      EXPECT_NEAR(
          round_trip[i], axis_angle[i], std::numeric_limits<double>::epsilon());
    }
  }
}

// Transposes a 3x3 matrix.
static void Transpose3x3(double m[9]) {
  std::swap(m[1], m[3]);
  std::swap(m[2], m[6]);
  std::swap(m[5], m[7]);
}

// Convert Euler angles from radians to degrees.
static void ToDegrees(double euler_angles[3]) {
  for (int i = 0; i < 3; ++i) {
    euler_angles[i] *= 180.0 / kPi;
  }
}

// Compare the 3x3 rotation matrices produced by the axis-angle
// rotation 'aa' and the Euler angle rotation 'ea' (in radians).
static void CompareEulerToAngleAxis(double aa[3], double ea[3]) {
  double aa_matrix[9];
  AngleAxisToRotationMatrix(aa, aa_matrix);
  Transpose3x3(aa_matrix);  // Column to row major order.

  double ea_matrix[9];
  ToDegrees(ea);  // Radians to degrees.
  const int kRowStride = 3;
  EulerAnglesToRotationMatrix(ea, kRowStride, ea_matrix);

  EXPECT_THAT(aa_matrix, IsOrthonormal());
  EXPECT_THAT(ea_matrix, IsOrthonormal());
  EXPECT_THAT(ea_matrix, IsNear3x3Matrix(aa_matrix));
}

// Test with rotation axis along the x/y/z axes.
// Also test zero rotation.
TEST(EulerAnglesToRotationMatrix, OnAxis) {
  int n_tests = 0;
  for (double x = -1.0; x <= 1.0; x += 1.0) {
    for (double y = -1.0; y <= 1.0; y += 1.0) {
      for (double z = -1.0; z <= 1.0; z += 1.0) {
        if ((x != 0) + (y != 0) + (z != 0) > 1) continue;
        double axis_angle[3] = {x, y, z};
        double euler_angles[3] = {x, y, z};
        CompareEulerToAngleAxis(axis_angle, euler_angles);
        ++n_tests;
      }
    }
  }
  CHECK_EQ(7, n_tests);
}

// Test that a random rotation produces an orthonormal rotation
// matrix.
TEST(EulerAnglesToRotationMatrix, IsOrthonormal) {
  std::mt19937 prng;
  std::uniform_real_distribution<double> uniform_distribution{-180.0, 180.0};
  for (int trial = 0; trial < kNumTrials; ++trial) {
    double euler_angles_degrees[3];
    for (double& euler_angles_degree : euler_angles_degrees) {
      euler_angles_degree = uniform_distribution(prng);
    }
    double rotation_matrix[9];
    EulerAnglesToRotationMatrix(euler_angles_degrees, 3, rotation_matrix);
    EXPECT_THAT(rotation_matrix, IsOrthonormal());
  }
}

static double sample_euler[][3] = {{0.5235988, 1.047198, 0.7853982},
                                   {0.5235988, 1.047198, 0.5235988},
                                   {0.7853982, 0.5235988, 1.047198}};

// ZXY Intrinsic Euler Angle to rotation matrix conversion test from
// scipy/spatial/transform/test/test_rotation.py
TEST(EulerAngles, IntrinsicEulerSequence312ToRotationMatrixCanned) {
  // clang-format off
  double const expected[][9] =
      {{0.306186083320088, -0.249999816228639,  0.918558748402491,
        0.883883627842492,  0.433012359189203, -0.176776777947208,
       -0.353553128699351,  0.866025628186053,  0.353553102817459},
      { 0.533493553519713, -0.249999816228639,  0.808012821828067,
        0.808012821828067,  0.433012359189203, -0.399519181705765,
       -0.249999816228639,  0.866025628186053,  0.433012359189203},
      { 0.047366781483451, -0.612372449482883,  0.789149143778432,
        0.659739427618959,  0.612372404654096,  0.435596057905909,
       -0.750000183771249,  0.500000021132493,  0.433012359189203}};
  // clang-format on

  for (int i = 0; i < 3; ++i) {
    double results[9];
    EulerAnglesToRotation<IntrinsicZXY>(sample_euler[i], results);
    ASSERT_THAT(results, IsNear3x3Matrix(expected[i]));
  }
}

// ZXY Extrinsic Euler Angle to rotation matrix conversion test from
// scipy/spatial/transform/test/test_rotation.py
TEST(EulerAngles, ExtrinsicEulerSequence312ToRotationMatrix) {
  // clang-format off
  double const expected[][9] =
      {{0.918558725988105,  0.176776842651999,  0.353553128699352,
        0.249999816228639,  0.433012359189203, -0.866025628186053,
       -0.306186150563275,  0.883883614901527,  0.353553102817459},
      { 0.966506404215301, -0.058012606358071,  0.249999816228639,
        0.249999816228639,  0.433012359189203, -0.866025628186053,
       -0.058012606358071,  0.899519223970752,  0.433012359189203},
      { 0.659739424151467, -0.047366829779744,  0.750000183771249,
        0.612372449482883,  0.612372404654096, -0.500000021132493,
       -0.435596000136163,  0.789149175666285,  0.433012359189203}};
  // clang-format on

  for (int i = 0; i < 3; ++i) {
    double results[9];
    EulerAnglesToRotation<ExtrinsicZXY>(sample_euler[i], results);
    ASSERT_THAT(results, IsNear3x3Matrix(expected[i]));
  }
}

// ZXZ Intrinsic Euler Angle to rotation matrix conversion test from
// scipy/spatial/transform/test/test_rotation.py
TEST(EulerAngles, IntrinsicEulerSequence313ToRotationMatrix) {
  // clang-format off
  double expected[][9] =
      {{0.435595832832961, -0.789149008363071,  0.433012832394307,
        0.659739379322704, -0.047367454164077, -0.750000183771249,
        0.612372616786097,  0.612372571957297,  0.499999611324802},
      { 0.625000065470068, -0.649518902838302,  0.433012832394307,
        0.649518902838302,  0.124999676794869, -0.750000183771249,
        0.433012832394307,  0.750000183771249,  0.499999611324802},
      {-0.176777132429787, -0.918558558684756,  0.353553418477159,
        0.883883325123719, -0.306186652473014, -0.353553392595246,
        0.433012832394307,  0.249999816228639,  0.866025391583588}};
  // clang-format on
  for (int i = 0; i < 3; ++i) {
    double results[9];
    EulerAnglesToRotation<IntrinsicZXZ>(sample_euler[i], results);
    ASSERT_THAT(results, IsNear3x3Matrix(expected[i]));
  }
}

// ZXZ Extrinsic Euler Angle to rotation matrix conversion test from
// scipy/spatial/transform/test/test_rotation.py
TEST(EulerAngles, ExtrinsicEulerSequence313ToRotationMatrix) {
  // clang-format off
  double expected[][9] =
      {{0.435595832832961, -0.659739379322704,  0.612372616786097,
        0.789149008363071, -0.047367454164077, -0.612372571957297,
        0.433012832394307,  0.750000183771249,  0.499999611324802},
      { 0.625000065470068, -0.649518902838302,  0.433012832394307,
        0.649518902838302,  0.124999676794869, -0.750000183771249,
        0.433012832394307,  0.750000183771249,  0.499999611324802},
      {-0.176777132429787, -0.883883325123719,  0.433012832394307,
        0.918558558684756, -0.306186652473014, -0.249999816228639,
        0.353553418477159,  0.353553392595246,  0.866025391583588}};
  // clang-format on
  for (int i = 0; i < 3; ++i) {
    double results[9];
    EulerAnglesToRotation<ExtrinsicZXZ>(sample_euler[i], results);
    ASSERT_THAT(results, IsNear3x3Matrix(expected[i]));
  }
}

template <typename T>
struct GeneralEulerAngles : public ::testing::Test {
 public:
  static constexpr bool kIsParityOdd = T::kIsParityOdd;
  static constexpr bool kIsProperEuler = T::kIsProperEuler;
  static constexpr bool kIsIntrinsic = T::kIsIntrinsic;

  template <typename URBG>
  static void RandomEulerAngles(double* euler, URBG& prng) {
    using ParamType = std::uniform_real_distribution<double>::param_type;
    std::uniform_real_distribution<double> uniform_distribution{-kPi, kPi};
    // Euler angles should be in
    //   [-pi,pi) x [0,pi) x [-pi,pi])
    // if the outer axes are repeated and
    //   [-pi,pi) x [-pi/2,pi/2) x [-pi,pi])
    // otherwise
    euler[0] = uniform_distribution(prng);
    euler[2] = uniform_distribution(prng);
    if constexpr (kIsProperEuler) {
      euler[1] = uniform_distribution(prng, ParamType{0, kPi});
    } else {
      euler[1] = uniform_distribution(prng, ParamType{-kPi / 2, kPi / 2});
    }
  }

  static void CheckPrincipalRotationMatrixProduct(double angles[3]) {
    // Convert Shoemake's Euler angle convention into 'apparent' rotation axes
    // sequences, i.e. the alphabetic code (ZYX, ZYZ, etc.) indicates in what
    // sequence rotations about different axes are applied
    constexpr int i = T::kAxes[0];
    constexpr int j = (3 + (kIsParityOdd ? (i - 1) % 3 : (i + 1) % 3)) % 3;
    constexpr int k = kIsProperEuler ? i : 3 ^ i ^ j;
    constexpr auto kSeq =
        kIsIntrinsic ? std::array{k, j, i} : std::array{i, j, k};

    double aa_matrix[9];
    Eigen::Map<Eigen::Matrix3d, 0, Eigen::Stride<1, 3>> aa(aa_matrix);
    aa.setIdentity();
    for (int i = 0; i < 3; ++i) {
      Eigen::Vector3d angle_axis;
      if constexpr (kIsIntrinsic) {
        angle_axis = -angles[i] * Eigen::Vector3d::Unit(kSeq[i]);
      } else {
        angle_axis = angles[i] * Eigen::Vector3d::Unit(kSeq[i]);
      }
      Eigen::Matrix3d m;
      AngleAxisToRotationMatrix(angle_axis.data(), m.data());
      aa = m * aa;
    }
    if constexpr (kIsIntrinsic) {
      aa.transposeInPlace();
    }

    double ea_matrix[9];
    EulerAnglesToRotation<T>(angles, ea_matrix);

    EXPECT_THAT(aa_matrix, IsOrthonormal());
    EXPECT_THAT(ea_matrix, IsOrthonormal());
    EXPECT_THAT(ea_matrix, IsNear3x3Matrix(aa_matrix));
  }
};

using EulerSystemList = ::testing::Types<ExtrinsicXYZ,
                                         ExtrinsicXYX,
                                         ExtrinsicXZY,
                                         ExtrinsicXZX,
                                         ExtrinsicYZX,
                                         ExtrinsicYZY,
                                         ExtrinsicYXZ,
                                         ExtrinsicYXY,
                                         ExtrinsicZXY,
                                         ExtrinsicZXZ,
                                         ExtrinsicZYX,
                                         ExtrinsicZYZ,
                                         IntrinsicZYX,
                                         IntrinsicXYX,
                                         IntrinsicYZX,
                                         IntrinsicXZX,
                                         IntrinsicXZY,
                                         IntrinsicYZY,
                                         IntrinsicZXY,
                                         IntrinsicYXY,
                                         IntrinsicYXZ,
                                         IntrinsicZXZ,
                                         IntrinsicXYZ,
                                         IntrinsicZYZ>;
TYPED_TEST_SUITE(GeneralEulerAngles, EulerSystemList);

TYPED_TEST(GeneralEulerAngles, EulerAnglesToRotationMatrixAndBack) {
  std::mt19937 prng;
  std::uniform_real_distribution<double> uniform_distribution{-1.0, 1.0};
  for (int i = 0; i < kNumTrials; ++i) {
    double euler[3];
    TestFixture::RandomEulerAngles(euler, prng);

    double matrix[9];
    double round_trip[3];
    EulerAnglesToRotation<TypeParam>(euler, matrix);
    ASSERT_THAT(matrix, IsOrthonormal());
    RotationMatrixToEulerAngles<TypeParam>(matrix, round_trip);
    for (int j = 0; j < 3; ++j)
      ASSERT_NEAR(euler[j], round_trip[j], 128.0 * kLooseTolerance);
  }
}

// Check that the rotation matrix converted from euler angles is equivalent to
// product of three principal axis rotation matrices
//     R_euler = R_a2(euler_2) * R_a1(euler_1) * R_a0(euler_0)
TYPED_TEST(GeneralEulerAngles, PrincipalRotationMatrixProduct) {
  std::mt19937 prng;
  double euler[3];
  for (int i = 0; i < kNumTrials; ++i) {
    TestFixture::RandomEulerAngles(euler, prng);
    TestFixture::CheckPrincipalRotationMatrixProduct(euler);
  }
}

// Gimbal lock (euler[1] == +/-pi) handling test. If a rotation matrix
// represents a gimbal-locked configuration, then converting this rotation
// matrix to euler angles and back must produce the same rotation matrix.
//
// From scipy/spatial/transform/test/test_rotation.py, but additionally covers
// gimbal lock handling for proper euler angles, which scipy appears to fail to
// do properly.
TYPED_TEST(GeneralEulerAngles, GimbalLocked) {
  constexpr auto kBoundaryAngles = TestFixture::kIsProperEuler
                                       ? std::array{0.0, kPi}
                                       : std::array{-kPi / 2, kPi / 2};
  constexpr double gimbal_locked_configurations[4][3] = {
      {0.78539816, kBoundaryAngles[1], 0.61086524},
      {0.61086524, kBoundaryAngles[0], 0.34906585},
      {0.61086524, kBoundaryAngles[1], 0.43633231},
      {0.43633231, kBoundaryAngles[0], 0.26179939}};
  double angle_estimates[3];
  double mat_expected[9];
  double mat_estimated[9];
  for (const auto& euler_angles : gimbal_locked_configurations) {
    EulerAnglesToRotation<TypeParam>(euler_angles, mat_expected);
    RotationMatrixToEulerAngles<TypeParam>(mat_expected, angle_estimates);
    EulerAnglesToRotation<TypeParam>(angle_estimates, mat_estimated);
    ASSERT_THAT(mat_expected, IsNear3x3Matrix(mat_estimated));
  }
}

// Tests using Jets for specific behavior involving auto differentiation
// near singularity points.

using J3 = Jet<double, 3>;
using J4 = Jet<double, 4>;

namespace {

// Converts an array of N real numbers (doubles) to an array of jets
template <int N>
void ArrayToArrayOfJets(double const* const src, Jet<double, N>* dst) {
  for (int i = 0; i < N; ++i) {
    dst[i] = Jet<double, N>(src[i], i);
  }
}

// Generically initializes a Jet with type T and a N-dimensional dual part
// N is explicitly given (instead of inferred from sizeof...(Ts)) so that the
// dual part can be initialized from Eigen expressions
template <int N, typename T, typename... Ts>
Jet<T, N> MakeJet(T a, const T& v0, Ts&&... v) {
  Jet<T, N> j;
  j.a = a;                                  // Real part
  ((j.v << v0), ..., std::forward<Ts>(v));  // Fill dual part with N components
  return j;
}

J3 MakeJ3(double a, double v0, double v1, double v2) {
  J3 j;
  j.a = a;
  j.v[0] = v0;
  j.v[1] = v1;
  j.v[2] = v2;
  return j;
}

J4 MakeJ4(double a, double v0, double v1, double v2, double v3) {
  J4 j;
  j.a = a;
  j.v[0] = v0;
  j.v[1] = v1;
  j.v[2] = v2;
  j.v[3] = v3;
  return j;
}

}  // namespace

// Use EXPECT_THAT(x, testing::PointWise(JetClose(prec), y); to achieve Jet
// array comparison
MATCHER_P(JetClose, relative_precision, "") {
  using internal::IsClose;
  using LHSJetType = std::remove_reference_t<std::tuple_element_t<0, arg_type>>;
  using RHSJetType = std::remove_reference_t<std::tuple_element_t<1, arg_type>>;

  constexpr int kDualPartDimension = LHSJetType::DIMENSION;
  static_assert(
      kDualPartDimension == RHSJetType::DIMENSION,
      "Can only compare Jets with dual parts having equal dimensions");
  auto&& [x, y] = arg;
  double relative_error;
  double absolute_error;
  if (!IsClose(
          x.a, y.a, relative_precision, &relative_error, &absolute_error)) {
    *result_listener << "Real part mismatch: x.a = " << x.a
                     << " and y.a = " << y.a
                     << " where the relative error between them is "
                     << relative_error
                     << " and the absolute error between them is "
                     << absolute_error;
    return false;
  }
  for (int i = 0; i < kDualPartDimension; i++) {
    if (!IsClose(x.v[i],
                 y.v[i],
                 relative_precision,
                 &relative_error,
                 &absolute_error)) {
      *result_listener << "Dual part mismatch: x.v[" << i << "] = " << x.v[i]
                       << " and y.v[" << i << "] = " << y.v[i]
                       << " where the relative error between them is "
                       << relative_error
                       << " and the absolute error between them is "
                       << absolute_error;
      return false;
    }
  }
  return true;
}

// Log-10 of a value well below machine precision.
static const int kSmallTinyCutoff = static_cast<int>(
    2 * log(std::numeric_limits<double>::epsilon()) / log(10.0));

// Log-10 of a value just below values representable by double.
static const int kTinyZeroLimit =
    static_cast<int>(1 + log(std::numeric_limits<double>::min()) / log(10.0));

// Test that exact conversion works for small angles when jets are used.
TEST(Rotation, SmallAngleAxisToQuaternionForJets) {
  // Examine small x rotations that are still large enough
  // to be well within the range represented by doubles.
  for (int i = -2; i >= kSmallTinyCutoff; i--) {
    double theta = pow(10.0, i);
    J3 axis_angle[3] = {J3(theta, 0), J3(0, 1), J3(0, 2)};
    J3 quaternion[4];
    J3 expected[4] = {
        MakeJ3(cos(theta / 2), -sin(theta / 2) / 2, 0, 0),
        MakeJ3(sin(theta / 2), cos(theta / 2) / 2, 0, 0),
        MakeJ3(0, 0, sin(theta / 2) / theta, 0),
        MakeJ3(0, 0, 0, sin(theta / 2) / theta),
    };
    AngleAxisToQuaternion(axis_angle, quaternion);
    EXPECT_THAT(quaternion, testing::Pointwise(JetClose(kTolerance), expected));
  }
}

// Test that conversion works for very small angles when jets are used.
TEST(Rotation, TinyAngleAxisToQuaternionForJets) {
  // Examine tiny x rotations that extend all the way to where
  // underflow occurs.
  for (int i = kSmallTinyCutoff; i >= kTinyZeroLimit; i--) {
    double theta = pow(10.0, i);
    J3 axis_angle[3] = {J3(theta, 0), J3(0, 1), J3(0, 2)};
    J3 quaternion[4];
    // To avoid loss of precision in the test itself,
    // a finite expansion is used here, which will
    // be exact up to machine precision for the test values used.
    J3 expected[4] = {
        MakeJ3(1.0, 0, 0, 0),
        MakeJ3(0, 0.5, 0, 0),
        MakeJ3(0, 0, 0.5, 0),
        MakeJ3(0, 0, 0, 0.5),
    };
    AngleAxisToQuaternion(axis_angle, quaternion);
    EXPECT_THAT(quaternion, testing::Pointwise(JetClose(kTolerance), expected));
  }
}

// Test that derivatives are correct for zero rotation.
TEST(Rotation, ZeroAngleAxisToQuaternionForJets) {
  J3 axis_angle[3] = {J3(0, 0), J3(0, 1), J3(0, 2)};
  J3 quaternion[4];
  J3 expected[4] = {
      MakeJ3(1.0, 0, 0, 0),
      MakeJ3(0, 0.5, 0, 0),
      MakeJ3(0, 0, 0.5, 0),
      MakeJ3(0, 0, 0, 0.5),
  };
  AngleAxisToQuaternion(axis_angle, quaternion);
  EXPECT_THAT(quaternion, testing::Pointwise(JetClose(kTolerance), expected));
}

// Test that exact conversion works for small angles.
TEST(Rotation, SmallQuaternionToAngleAxisForJets) {
  // Examine small x rotations that are still large enough
  // to be well within the range represented by doubles.
  for (int i = -2; i >= kSmallTinyCutoff; i--) {
    double theta = pow(10.0, i);
    double s = sin(theta);
    double c = cos(theta);
    J4 quaternion[4] = {J4(c, 0), J4(s, 1), J4(0, 2), J4(0, 3)};
    J4 axis_angle[3];
    // clang-format off
    J4 expected[3] = {
        MakeJ4(2*theta, -2*s, 2*c,  0,         0),
        MakeJ4(0,        0,   0,    2*theta/s, 0),
        MakeJ4(0,        0,   0,    0,         2*theta/s),
    };
    // clang-format on
    QuaternionToAngleAxis(quaternion, axis_angle);
    EXPECT_THAT(axis_angle, testing::Pointwise(JetClose(kTolerance), expected));
  }
}

// Test that conversion works for very small angles.
TEST(Rotation, TinyQuaternionToAngleAxisForJets) {
  // Examine tiny x rotations that extend all the way to where
  // underflow occurs.
  for (int i = kSmallTinyCutoff; i >= kTinyZeroLimit; i--) {
    double theta = pow(10.0, i);
    double s = sin(theta);
    double c = cos(theta);
    J4 quaternion[4] = {J4(c, 0), J4(s, 1), J4(0, 2), J4(0, 3)};
    J4 axis_angle[3];
    // To avoid loss of precision in the test itself,
    // a finite expansion is used here, which will
    // be exact up to machine precision for the test values used.
    // clang-format off
    J4 expected[3] = {
        MakeJ4(2*theta, -2*s, 2.0, 0,   0),
        MakeJ4(0,        0,   0,   2.0, 0),
        MakeJ4(0,        0,   0,   0,   2.0),
    };
    // clang-format on
    QuaternionToAngleAxis(quaternion, axis_angle);
    EXPECT_THAT(axis_angle, testing::Pointwise(JetClose(kTolerance), expected));
  }
}

// Test that conversion works for no rotation.
TEST(Rotation, ZeroQuaternionToAngleAxisForJets) {
  J4 quaternion[4] = {J4(1, 0), J4(0, 1), J4(0, 2), J4(0, 3)};
  J4 axis_angle[3];
  J4 expected[3] = {
      MakeJ4(0, 0, 2.0, 0, 0),
      MakeJ4(0, 0, 0, 2.0, 0),
      MakeJ4(0, 0, 0, 0, 2.0),
  };
  QuaternionToAngleAxis(quaternion, axis_angle);
  EXPECT_THAT(axis_angle, testing::Pointwise(JetClose(kTolerance), expected));
}

// The following 4 test cases cover the conversion of Euler Angles to rotation
// matrices for Jets
//
// The dual parts (with dimension 3) of the resultant matrix of Jets contain the
// derivative of each matrix element w.r.t. the input Euler Angles. In other
// words, for each element in R = EulerAnglesToRotationMatrix(angles), we have
// R_ij.v = jacobian(R_ij, angles)
//
// The test data (dual parts of the Jets) is generated by analytically
// differentiating the formulas for Euler Angle to Rotation Matrix conversion

// Test ZXY/312 Intrinsic Euler Angles to rotation matrix conversion using Jets
// The two ZXY test cases specifically cover handling of Tait-Bryan angles
// i.e. last axis of rotation is different from the first
TEST(EulerAngles, Intrinsic312EulerSequenceToRotationMatrixForJets) {
  J3 euler_angles[3];
  J3 rotation_matrix[9];

  ArrayToArrayOfJets(sample_euler[0], euler_angles);
  EulerAnglesToRotation<IntrinsicZXY>(euler_angles, rotation_matrix);
  {
    // clang-format off
    const J3 expected[] = {
      MakeJ3( 0.306186083320, -0.883883627842, -0.176776571821, -0.918558748402),  // NOLINT
      MakeJ3(-0.249999816229, -0.433012359189,  0.433012832394,  0.000000000000),  // NOLINT
      MakeJ3( 0.918558748402,  0.176776777947,  0.176776558880,  0.306186083320),  // NOLINT
      MakeJ3( 0.883883627842,  0.306186083320,  0.306185986727,  0.176776777947),  // NOLINT
      MakeJ3( 0.433012359189, -0.249999816229, -0.750000183771,  0.000000000000),  // NOLINT
      MakeJ3(-0.176776777947,  0.918558748402, -0.306185964313,  0.883883627842),  // NOLINT
      MakeJ3(-0.353553128699,  0.000000000000,  0.612372616786, -0.353553102817),  // NOLINT
      MakeJ3( 0.866025628186,  0.000000000000,  0.499999611325,  0.000000000000),  // NOLINT
      MakeJ3( 0.353553102817,  0.000000000000, -0.612372571957, -0.353553128699)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(rotation_matrix,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_euler[1], euler_angles);
  EulerAnglesToRotation<IntrinsicZXY>(euler_angles, rotation_matrix);
  {
    // clang-format off
    const J3 expected[] = {
      MakeJ3( 0.533493553520, -0.808012821828, -0.124999913397, -0.808012821828),  // NOLINT
      MakeJ3(-0.249999816229, -0.433012359189,  0.433012832394,  0.000000000000),  // NOLINT
      MakeJ3( 0.808012821828,  0.399519181706,  0.216506188745,  0.533493553520),  // NOLINT
      MakeJ3( 0.808012821828,  0.533493553520,  0.216506188745,  0.399519181706),  // NOLINT
      MakeJ3( 0.433012359189, -0.249999816229, -0.750000183771,  0.000000000000),  // NOLINT
      MakeJ3(-0.399519181706,  0.808012821828, -0.374999697927,  0.808012821828),  // NOLINT
      MakeJ3(-0.249999816229,  0.000000000000,  0.433012832394, -0.433012359189),  // NOLINT
      MakeJ3( 0.866025628186,  0.000000000000,  0.499999611325,  0.000000000000),  // NOLINT
      MakeJ3( 0.433012359189,  0.000000000000, -0.750000183771, -0.249999816229)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(rotation_matrix,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_euler[2], euler_angles);
  EulerAnglesToRotation<IntrinsicZXY>(euler_angles, rotation_matrix);
  {
    // clang-format off
    const J3 expected[] = {
      MakeJ3( 0.047366781483, -0.659739427619, -0.530330235247, -0.789149143778),  // NOLINT
      MakeJ3(-0.612372449483, -0.612372404654,  0.353553418477,  0.000000000000),  // NOLINT
      MakeJ3( 0.789149143778, -0.435596057906,  0.306185986727,  0.047366781483),  // NOLINT
      MakeJ3( 0.659739427619,  0.047366781483,  0.530330196424, -0.435596057906),  // NOLINT
      MakeJ3( 0.612372404654, -0.612372449483, -0.353553392595,  0.000000000000),  // NOLINT
      MakeJ3( 0.435596057906,  0.789149143778, -0.306185964313,  0.659739427619),  // NOLINT
      MakeJ3(-0.750000183771,  0.000000000000,  0.433012832394, -0.433012359189),  // NOLINT
      MakeJ3( 0.500000021132,  0.000000000000,  0.866025391584,  0.000000000000),  // NOLINT
      MakeJ3( 0.433012359189,  0.000000000000, -0.249999816229, -0.750000183771)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(rotation_matrix,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }
}

// Test ZXY/312 Extrinsic Euler Angles to rotation matrix conversion using Jets
TEST(EulerAngles, Extrinsic312EulerSequenceToRotationMatrixForJets) {
  J3 euler_angles[3];
  J3 rotation_matrix[9];

  ArrayToArrayOfJets(sample_euler[0], euler_angles);
  EulerAnglesToRotation<ExtrinsicZXY>(euler_angles, rotation_matrix);
  {
    // clang-format off
    const J3 expected[] = {
      MakeJ3( 0.918558725988,  0.176776842652,  0.176776571821, -0.306186150563),  // NOLINT
      MakeJ3( 0.176776842652, -0.918558725988,  0.306185986727,  0.883883614902),  // NOLINT
      MakeJ3( 0.353553128699,  0.000000000000, -0.612372616786,  0.353553102817),  // NOLINT
      MakeJ3( 0.249999816229,  0.433012359189, -0.433012832394,  0.000000000000),  // NOLINT
      MakeJ3( 0.433012359189, -0.249999816229, -0.750000183771,  0.000000000000),  // NOLINT
      MakeJ3(-0.866025628186,  0.000000000000, -0.499999611325,  0.000000000000),  // NOLINT
      MakeJ3(-0.306186150563,  0.883883614902,  0.176776558880, -0.918558725988),  // NOLINT
      MakeJ3( 0.883883614902,  0.306186150563,  0.306185964313, -0.176776842652),  // NOLINT
      MakeJ3( 0.353553102817,  0.000000000000, -0.612372571957, -0.353553128699)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(rotation_matrix,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_euler[1], euler_angles);
  EulerAnglesToRotation<ExtrinsicZXY>(euler_angles, rotation_matrix);
  {
    // clang-format off
    const J3 expected[] = {
      MakeJ3( 0.966506404215, -0.058012606358,  0.124999913397, -0.058012606358),  // NOLINT
      MakeJ3(-0.058012606358, -0.966506404215,  0.216506188745,  0.899519223971),  // NOLINT
      MakeJ3( 0.249999816229,  0.000000000000, -0.433012832394,  0.433012359189),  // NOLINT
      MakeJ3( 0.249999816229,  0.433012359189, -0.433012832394,  0.000000000000),  // NOLINT
      MakeJ3( 0.433012359189, -0.249999816229, -0.750000183771,  0.000000000000),  // NOLINT
      MakeJ3(-0.866025628186,  0.000000000000, -0.499999611325,  0.000000000000),  // NOLINT
      MakeJ3(-0.058012606358,  0.899519223971,  0.216506188745, -0.966506404215),  // NOLINT
      MakeJ3( 0.899519223971,  0.058012606358,  0.374999697927,  0.058012606358),  // NOLINT
      MakeJ3( 0.433012359189,  0.000000000000, -0.750000183771, -0.249999816229)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(rotation_matrix,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_euler[2], euler_angles);
  EulerAnglesToRotation<ExtrinsicZXY>(euler_angles, rotation_matrix);
  {
    // clang-format off
    const J3 expected[] = {
      MakeJ3( 0.659739424151, -0.047366829780,  0.530330235247, -0.435596000136),  // NOLINT
      MakeJ3(-0.047366829780, -0.659739424151,  0.530330196424,  0.789149175666),  // NOLINT
      MakeJ3( 0.750000183771,  0.000000000000, -0.433012832394,  0.433012359189),  // NOLINT
      MakeJ3( 0.612372449483,  0.612372404654, -0.353553418477,  0.000000000000),  // NOLINT
      MakeJ3( 0.612372404654, -0.612372449483, -0.353553392595,  0.000000000000),  // NOLINT
      MakeJ3(-0.500000021132,  0.000000000000, -0.866025391584,  0.000000000000),  // NOLINT
      MakeJ3(-0.435596000136,  0.789149175666,  0.306185986727, -0.659739424151),  // NOLINT
      MakeJ3( 0.789149175666,  0.435596000136,  0.306185964313,  0.047366829780),  // NOLINT
      MakeJ3( 0.433012359189,  0.000000000000, -0.249999816229, -0.750000183771)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(rotation_matrix,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }
}

// Test ZXZ/313 Intrinsic Euler Angles to rotation matrix conversion using Jets
// The two ZXZ test cases specifically cover handling of proper Euler Sequences
// i.e. last axis of rotation is same as the first
TEST(EulerAngles, Intrinsic313EulerSequenceToRotationMatrixForJets) {
  J3 euler_angles[3];
  J3 rotation_matrix[9];

  ArrayToArrayOfJets(sample_euler[0], euler_angles);
  EulerAnglesToRotation<IntrinsicZXZ>(euler_angles, rotation_matrix);
  {
    // clang-format off
    const J3 expected[] = {
      MakeJ3( 0.435595832833, -0.659739379323,  0.306186321334, -0.789149008363),  // NOLINT
      MakeJ3(-0.789149008363,  0.047367454164,  0.306186298920, -0.435595832833),  // NOLINT
      MakeJ3( 0.433012832394,  0.750000183771,  0.249999816229,  0.000000000000),  // NOLINT
      MakeJ3( 0.659739379323,  0.435595832833, -0.530330235247, -0.047367454164),  // NOLINT
      MakeJ3(-0.047367454164, -0.789149008363, -0.530330196424, -0.659739379323),  // NOLINT
      MakeJ3(-0.750000183771,  0.433012832394, -0.433012359189,  0.000000000000),  // NOLINT
      MakeJ3( 0.612372616786,  0.000000000000,  0.353553128699,  0.612372571957),  // NOLINT
      MakeJ3( 0.612372571957,  0.000000000000,  0.353553102817, -0.612372616786),  // NOLINT
      MakeJ3( 0.499999611325,  0.000000000000, -0.866025628186,  0.000000000000)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(rotation_matrix,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_euler[1], euler_angles);
  EulerAnglesToRotation<IntrinsicZXZ>(euler_angles, rotation_matrix);
  {
    // clang-format off
    const J3 expected[] = {
      MakeJ3( 0.625000065470, -0.649518902838,  0.216506425348, -0.649518902838),  // NOLINT
      MakeJ3(-0.649518902838, -0.124999676795,  0.375000107735, -0.625000065470),  // NOLINT
      MakeJ3( 0.433012832394,  0.750000183771,  0.249999816229,  0.000000000000),  // NOLINT
      MakeJ3( 0.649518902838,  0.625000065470, -0.375000107735,  0.124999676795),  // NOLINT
      MakeJ3( 0.124999676795, -0.649518902838, -0.649519202838, -0.649518902838),  // NOLINT
      MakeJ3(-0.750000183771,  0.433012832394, -0.433012359189,  0.000000000000),  // NOLINT
      MakeJ3( 0.433012832394,  0.000000000000,  0.249999816229,  0.750000183771),  // NOLINT
      MakeJ3( 0.750000183771,  0.000000000000,  0.433012359189, -0.433012832394),  // NOLINT
      MakeJ3( 0.499999611325,  0.000000000000, -0.866025628186,  0.000000000000)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(rotation_matrix,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_euler[2], euler_angles);
  EulerAnglesToRotation<IntrinsicZXZ>(euler_angles, rotation_matrix);
  {
    // clang-format off
    const J3 expected[] = {
      MakeJ3(-0.176777132430, -0.883883325124,  0.306186321334, -0.918558558685),  // NOLINT
      MakeJ3(-0.918558558685,  0.306186652473,  0.176776571821,  0.176777132430),  // NOLINT
      MakeJ3( 0.353553418477,  0.353553392595,  0.612372449483,  0.000000000000),  // NOLINT
      MakeJ3( 0.883883325124, -0.176777132430, -0.306186298920, -0.306186652473),  // NOLINT
      MakeJ3(-0.306186652473, -0.918558558685, -0.176776558880, -0.883883325124),  // NOLINT
      MakeJ3(-0.353553392595,  0.353553418477, -0.612372404654,  0.000000000000),  // NOLINT
      MakeJ3( 0.433012832394,  0.000000000000,  0.750000183771,  0.249999816229),  // NOLINT
      MakeJ3( 0.249999816229,  0.000000000000,  0.433012359189, -0.433012832394),  // NOLINT
      MakeJ3( 0.866025391584,  0.000000000000, -0.500000021132,  0.000000000000)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(rotation_matrix,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }
}

// Test ZXZ/313 Extrinsic Euler Angles to rotation matrix conversion using Jets
TEST(EulerAngles, Extrinsic313EulerSequenceToRotationMatrixForJets) {
  J3 euler_angles[3];
  J3 rotation_matrix[9];

  ArrayToArrayOfJets(sample_euler[0], euler_angles);
  EulerAnglesToRotation<ExtrinsicZXZ>(euler_angles, rotation_matrix);
  {
    // clang-format off
    const J3 expected[] = {
      MakeJ3( 0.435595832833, -0.659739379323,  0.306186321334, -0.789149008363),  // NOLINT
      MakeJ3(-0.659739379323, -0.435595832833,  0.530330235247,  0.047367454164),  // NOLINT
      MakeJ3( 0.612372616786,  0.000000000000,  0.353553128699,  0.612372571957),  // NOLINT
      MakeJ3( 0.789149008363, -0.047367454164, -0.306186298920,  0.435595832833),  // NOLINT
      MakeJ3(-0.047367454164, -0.789149008363, -0.530330196424, -0.659739379323),  // NOLINT
      MakeJ3(-0.612372571957,  0.000000000000, -0.353553102817,  0.612372616786),  // NOLINT
      MakeJ3( 0.433012832394,  0.750000183771,  0.249999816229,  0.000000000000),  // NOLINT
      MakeJ3( 0.750000183771, -0.433012832394,  0.433012359189,  0.000000000000),  // NOLINT
      MakeJ3( 0.499999611325,  0.000000000000, -0.866025628186,  0.000000000000)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(rotation_matrix,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_euler[1], euler_angles);
  EulerAnglesToRotation<ExtrinsicZXZ>(euler_angles, rotation_matrix);
  {
    // clang-format off
    const J3 expected[] = {
      MakeJ3( 0.625000065470, -0.649518902838,  0.216506425348, -0.649518902838),  // NOLINT
      MakeJ3(-0.649518902838, -0.625000065470,  0.375000107735, -0.124999676795),  // NOLINT
      MakeJ3( 0.433012832394,  0.000000000000,  0.249999816229,  0.750000183771),  // NOLINT
      MakeJ3( 0.649518902838,  0.124999676795, -0.375000107735,  0.625000065470),  // NOLINT
      MakeJ3( 0.124999676795, -0.649518902838, -0.649519202838, -0.649518902838),  // NOLINT
      MakeJ3(-0.750000183771,  0.000000000000, -0.433012359189,  0.433012832394),  // NOLINT
      MakeJ3( 0.433012832394,  0.750000183771,  0.249999816229,  0.000000000000),  // NOLINT
      MakeJ3( 0.750000183771, -0.433012832394,  0.433012359189,  0.000000000000),  // NOLINT
      MakeJ3( 0.499999611325,  0.000000000000, -0.866025628186,  0.000000000000)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(rotation_matrix,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_euler[2], euler_angles);
  EulerAnglesToRotation<ExtrinsicZXZ>(euler_angles, rotation_matrix);
  {
    // clang-format off
    const J3 expected[] = {
      MakeJ3(-0.176777132430, -0.883883325124,  0.306186321334, -0.918558558685),  // NOLINT
      MakeJ3(-0.883883325124,  0.176777132430,  0.306186298920,  0.306186652473),  // NOLINT
      MakeJ3( 0.433012832394,  0.000000000000,  0.750000183771,  0.249999816229),  // NOLINT
      MakeJ3( 0.918558558685, -0.306186652473, -0.176776571821, -0.176777132430),  // NOLINT
      MakeJ3(-0.306186652473, -0.918558558685, -0.176776558880, -0.883883325124),  // NOLINT
      MakeJ3(-0.249999816229,  0.000000000000, -0.433012359189,  0.433012832394),  // NOLINT
      MakeJ3( 0.353553418477,  0.353553392595,  0.612372449483,  0.000000000000),  // NOLINT
      MakeJ3( 0.353553392595, -0.353553418477,  0.612372404654,  0.000000000000),  // NOLINT
      MakeJ3( 0.866025391584,  0.000000000000, -0.500000021132,  0.000000000000)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(rotation_matrix,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }
}

using J9 = Jet<double, 9>;

// The following 4 tests Tests the conversion of rotation matrices to Euler
// Angles for Jets.
//
// The dual parts (with dimension 9) of the resultant array of Jets contain the
// derivative of each Euler angle w.r.t. each of the 9 elements of the rotation
// matrix, or a 9-by-1 array formed from flattening the rotation matrix. In
// other words, for each element in angles = RotationMatrixToEulerAngles(R), we
// have angles.v = jacobian(angles, [R11 R12 R13 R21 ... R32 R33]);
//
// Note: the order of elements in v depend on row/column-wise flattening of
// the rotation matrix
//
// The test data (dual parts of the Jets) is generated by analytically
// differentiating the formulas for Rotation Matrix to Euler Angle conversion

// clang-format off
static double sample_matrices[][9] = {
  { 0.433012359189,  0.176776842652,  0.883883614902,  0.249999816229,  0.918558725988, -0.306186150563, -0.866025628186,  0.353553128699,  0.353553102817},  // NOLINT
  { 0.433012359189, -0.058012606358,  0.899519223971,  0.249999816229,  0.966506404215, -0.058012606358, -0.866025628186,  0.249999816229,  0.433012359189},  // NOLINT
  { 0.612372404654, -0.047366829780,  0.789149175666,  0.612372449483,  0.659739424151, -0.435596000136, -0.500000021132,  0.750000183771,  0.433012359189}   // NOLINT
};
// clang-format on

// Test rotation matrix to ZXY/312 Intrinsic Euler Angles conversion using Jets
// The two ZXY test cases specifically cover handling of Tait-Bryan angles
// i.e. last axis of rotation is different from the first
TEST(EulerAngles, RotationMatrixToIntrinsic312EulerSequenceForJets) {
  J9 euler_angles[3];
  J9 rotation_matrix[9];

  ArrayToArrayOfJets(sample_matrices[0], rotation_matrix);
  RotationMatrixToEulerAngles<IntrinsicZXY>(rotation_matrix, euler_angles);
  {
    // clang-format off
    const J9 expected[] = {
      MakeJet<9>(-0.190125743401,  0.000000000000, -1.049781178951,  0.000000000000,  0.000000000000,  0.202030634558,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000),  // NOLINT
      MakeJet<9>( 0.361366843930,  0.000000000000, -0.066815309609,  0.000000000000,  0.000000000000, -0.347182270882,  0.000000000000,  0.000000000000,  0.935414445680,  0.000000000000),  // NOLINT
      MakeJet<9>( 1.183200015636,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000, -0.404060603418,  0.000000000000, -0.989743365598)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(euler_angles,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_matrices[1], rotation_matrix);
  RotationMatrixToEulerAngles<IntrinsicZXY>(rotation_matrix, euler_angles);
  {
    // clang-format off
    const J9 expected[] = {
      MakeJet<9>( 0.059951064811,  0.000000000000, -1.030940063452,  0.000000000000,  0.000000000000, -0.061880107384,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000),  // NOLINT
      MakeJet<9>( 0.252680065344,  0.000000000000,  0.014978778808,  0.000000000000,  0.000000000000, -0.249550684831,  0.000000000000,  0.000000000000,  0.968245884001,  0.000000000000),  // NOLINT
      MakeJet<9>( 1.107149138016,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000, -0.461879804532,  0.000000000000, -0.923760579526)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(euler_angles,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_matrices[2], rotation_matrix);
  RotationMatrixToEulerAngles<IntrinsicZXY>(rotation_matrix, euler_angles);
  {
    // clang-format off
    const J9 expected[] = {
      MakeJet<9>( 0.071673287221,  0.000000000000, -1.507976776767,  0.000000000000,  0.000000000000, -0.108267107713,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000),  // NOLINT
      MakeJet<9>( 0.848062356818,  0.000000000000,  0.053708966648,  0.000000000000,  0.000000000000, -0.748074610289,  0.000000000000,  0.000000000000,  0.661437619389,  0.000000000000),  // NOLINT
      MakeJet<9>( 0.857072360427,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000, -0.989743158900,  0.000000000000, -1.142857911244)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(euler_angles,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }
}

// Test rotation matrix to ZXY/312 Extrinsic Euler Angles conversion using Jets
TEST(EulerAngles, RotationMatrixToExtrinsic312EulerSequenceForJets) {
  J9 euler_angles[3];
  J9 rotation_matrix[9];

  ArrayToArrayOfJets(sample_matrices[0], rotation_matrix);
  RotationMatrixToEulerAngles<ExtrinsicZXY>(rotation_matrix, euler_angles);
  {
    // clang-format off
    const J9 expected[] = {
      MakeJet<9>( 0.265728912717,  0.000000000000,  0.000000000000,  0.000000000000,  1.013581996386, -0.275861853641,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000),  // NOLINT
      MakeJet<9>( 0.311184173598,  0.000000000000,  0.000000000000, -0.284286741927,  0.000000000000,  0.000000000000, -0.951971659874,  0.000000000000,  0.000000000000, -0.113714586405),  // NOLINT
      MakeJet<9>( 1.190290284357,  0.000000000000,  0.000000000000,  0.390127543992,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000, -0.975319806582)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(euler_angles,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_matrices[1], rotation_matrix);
  RotationMatrixToEulerAngles<ExtrinsicZXY>(rotation_matrix, euler_angles);
  {
    // clang-format off
    const J9 expected[] = {
      MakeJet<9>( 0.253115668605,  0.000000000000,  0.000000000000,  0.000000000000,  0.969770129215, -0.250844022378,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000),  // NOLINT
      MakeJet<9>( 0.058045195612,  0.000000000000,  0.000000000000, -0.052271487648,  0.000000000000,  0.000000000000, -0.998315850572,  0.000000000000,  0.000000000000, -0.025162553041),  // NOLINT
      MakeJet<9>( 1.122153748896,  0.000000000000,  0.000000000000,  0.434474567050,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000, -0.902556744846)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(euler_angles,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_matrices[2], rotation_matrix);
  RotationMatrixToEulerAngles<ExtrinsicZXY>(rotation_matrix, euler_angles);
  {
    // clang-format off
    const J9 expected[] = {
      MakeJet<9>( 0.748180444286,  0.000000000000,  0.000000000000,  0.000000000000,  0.814235652244, -0.755776390750,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000),  // NOLINT
      MakeJet<9>( 0.450700288478,  0.000000000000,  0.000000000000, -0.381884322045,  0.000000000000,  0.000000000000, -0.900142280234,  0.000000000000,  0.000000000000, -0.209542930950),  // NOLINT
      MakeJet<9>( 1.068945699497,  0.000000000000,  0.000000000000,  0.534414175972,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000, -0.973950275281)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(euler_angles,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }
}

// Test rotation matrix to ZXZ/313 Intrinsic Euler Angles conversion using Jets
//// The two ZXZ test cases specifically cover handling of proper Euler
/// Sequences
// i.e. last axis of rotation is same as the first
TEST(EulerAngles, RotationMatrixToIntrinsic313EulerSequenceForJets) {
  J9 euler_angles[3];
  J9 rotation_matrix[9];

  ArrayToArrayOfJets(sample_matrices[0], rotation_matrix);
  RotationMatrixToEulerAngles<IntrinsicZXZ>(rotation_matrix, euler_angles);
  {
    // clang-format off
    const J9 expected[] = {
      MakeJet<9>( 1.237323270947,  0.000000000000,  0.000000000000,  0.349926947837,  0.000000000000,  0.000000000000,  1.010152467826,  0.000000000000,  0.000000000000,  0.000000000000),  // NOLINT
      MakeJet<9>( 1.209429510533,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000, -0.327326615680,  0.133630397662, -0.935414455462),  // NOLINT
      MakeJet<9>(-1.183199990019,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.404060624546,  0.989743344897,  0.000000000000)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(euler_angles,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_matrices[1], rotation_matrix);
  RotationMatrixToEulerAngles<IntrinsicZXZ>(rotation_matrix, euler_angles);
  {
    // clang-format off
    const J9 expected[] = {
      MakeJet<9>( 1.506392616830,  0.000000000000,  0.000000000000,  0.071400104821,  0.000000000000,  0.000000000000,  1.107100178948,  0.000000000000,  0.000000000000,  0.000000000000),  // NOLINT
      MakeJet<9>( 1.122964310061,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000, -0.416024849727,  0.120095910090, -0.901387983495),  // NOLINT
      MakeJet<9>(-1.289761690216,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.307691969119,  1.065877306886,  0.000000000000)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(euler_angles,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_matrices[2], rotation_matrix);
  RotationMatrixToEulerAngles<IntrinsicZXZ>(rotation_matrix, euler_angles);
  {
    // clang-format off
    const J9 expected[] = {
      MakeJet<9>( 1.066432836578,  0.000000000000,  0.000000000000,  0.536117958181,  0.000000000000,  0.000000000000,  0.971260169116,  0.000000000000,  0.000000000000,  0.000000000000),  // NOLINT
      MakeJet<9>( 1.122964310061,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000, -0.240192006893,  0.360288083393, -0.901387983495),  // NOLINT
      MakeJet<9>(-0.588002509965,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.923076812076,  0.615384416607,  0.000000000000)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(euler_angles,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }
}

// Test rotation matrix to ZXZ/313 Extrinsic Euler Angles conversion using Jets
TEST(EulerAngles, RotationMatrixToExtrinsic313EulerSequenceForJets) {
  J9 euler_angles[3];
  J9 rotation_matrix[9];

  ArrayToArrayOfJets(sample_matrices[0], rotation_matrix);
  RotationMatrixToEulerAngles<ExtrinsicZXZ>(rotation_matrix, euler_angles);
  {
    // clang-format off
    const J9 expected[] = {
      MakeJet<9>(-1.183199990019,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.404060624546,  0.989743344897,  0.000000000000),  // NOLINT
      MakeJet<9>( 1.209429510533,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000, -0.327326615680,  0.133630397662, -0.935414455462),  // NOLINT
      MakeJet<9>( 1.237323270947,  0.000000000000,  0.000000000000,  0.349926947837,  0.000000000000,  0.000000000000,  1.010152467826,  0.000000000000,  0.000000000000,  0.000000000000)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(euler_angles,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_matrices[1], rotation_matrix);
  RotationMatrixToEulerAngles<ExtrinsicZXZ>(rotation_matrix, euler_angles);
  {
    // clang-format off
    const J9 expected[] = {
      MakeJet<9>(-1.289761690216,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.307691969119,  1.065877306886,  0.000000000000),  // NOLINT
      MakeJet<9>( 1.122964310061,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000, -0.416024849727,  0.120095910090, -0.901387983495),  // NOLINT
      MakeJet<9>( 1.506392616830,  0.000000000000,  0.000000000000,  0.071400104821,  0.000000000000,  0.000000000000,  1.107100178948,  0.000000000000,  0.000000000000,  0.000000000000)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(euler_angles,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }

  ArrayToArrayOfJets(sample_matrices[2], rotation_matrix);
  RotationMatrixToEulerAngles<ExtrinsicZXZ>(rotation_matrix, euler_angles);
  {
    // clang-format off
    const J9 expected[] = {
      MakeJet<9>(-0.588002509965,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.923076812076,  0.615384416607,  0.000000000000),  // NOLINT
      MakeJet<9>( 1.122964310061,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000, -0.240192006893,  0.360288083393, -0.901387983495),  // NOLINT
      MakeJet<9>( 1.066432836578,  0.000000000000,  0.000000000000,  0.536117958181,  0.000000000000,  0.000000000000,  0.971260169116,  0.000000000000,  0.000000000000,  0.000000000000)   // NOLINT
    };
    // clang-format on
    EXPECT_THAT(euler_angles,
                testing::Pointwise(JetClose(kLooseTolerance), expected));
  }
}

TEST(Quaternion, RotatePointGivesSameAnswerAsRotationByMatrixCanned) {
  // Canned data generated in octave.
  double const q[4] = {
      +0.1956830471754074,
      -0.0150618562474847,
      +0.7634572982788086,
      -0.3019454777240753,
  };
  double const Q[3][3] = {
      // Scaled rotation matrix.
      {-0.6355194033477252, +0.0951730541682254, +0.3078870197911186},
      {-0.1411693904792992, +0.5297609702153905, -0.4551502574482019},
      {-0.2896955822708862, -0.4669396571547050, -0.4536309793389248},
  };
  double const R[3][3] = {
      // With unit rows and columns.
      {-0.8918859164053080, +0.1335655625725649, +0.4320876677394745},
      {-0.1981166751680096, +0.7434648665444399, -0.6387564287225856},
      {-0.4065578619806013, -0.6553016349046693, -0.6366242786393164},
  };

  // Compute R from q and compare to known answer.
  double Rq[3][3];
  QuaternionToScaledRotation<double>(q, Rq[0]);
  ExpectArraysClose(9, Q[0], Rq[0], kTolerance);

  // Now do the same but compute R with normalization.
  QuaternionToRotation<double>(q, Rq[0]);
  ExpectArraysClose(9, R[0], Rq[0], kTolerance);
}

TEST(Quaternion, RotatePointGivesSameAnswerAsRotationByMatrix) {
  // Rotation defined by a unit quaternion.
  double const q[4] = {
      +0.2318160216097109,
      -0.0178430356832060,
      +0.9044300776717159,
      -0.3576998641394597,
  };
  double const p[3] = {
      +0.11,
      -13.15,
      1.17,
  };

  double R[3 * 3];
  QuaternionToRotation(q, R);

  double result1[3];
  UnitQuaternionRotatePoint(q, p, result1);

  double result2[3];
  VectorRef(result2, 3) = ConstMatrixRef(R, 3, 3) * ConstVectorRef(p, 3);
  ExpectArraysClose(3, result1, result2, kTolerance);
}

// Verify that (a * b) * c == a * (b * c).
TEST(Quaternion, MultiplicationIsAssociative) {
  std::mt19937 prng;
  std::uniform_real_distribution<double> uniform_distribution{-1.0, 1.0};
  double a[4];
  double b[4];
  double c[4];
  for (int i = 0; i < 4; ++i) {
    a[i] = uniform_distribution(prng);
    b[i] = uniform_distribution(prng);
    c[i] = uniform_distribution(prng);
  }

  double ab[4];
  double ab_c[4];
  QuaternionProduct(a, b, ab);
  QuaternionProduct(ab, c, ab_c);

  double bc[4];
  double a_bc[4];
  QuaternionProduct(b, c, bc);
  QuaternionProduct(a, bc, a_bc);

  ASSERT_NEAR(ab_c[0], a_bc[0], kTolerance);
  ASSERT_NEAR(ab_c[1], a_bc[1], kTolerance);
  ASSERT_NEAR(ab_c[2], a_bc[2], kTolerance);
  ASSERT_NEAR(ab_c[3], a_bc[3], kTolerance);
}

TEST(AngleAxis, RotatePointGivesSameAnswerAsRotationMatrix) {
  std::mt19937 prng;
  std::uniform_real_distribution<double> uniform_distribution{-1.0, 1.0};
  double angle_axis[3];
  double R[9];
  double p[3];
  double angle_axis_rotated_p[3];
  double rotation_matrix_rotated_p[3];

  for (int i = 0; i < 10000; ++i) {
    double theta = (2.0 * i * 0.0011 - 1.0) * kPi;
    for (int j = 0; j < 50; ++j) {
      for (int k = 0; k < 3; ++k) {
        angle_axis[k] = uniform_distribution(prng);
        p[k] = uniform_distribution(prng);
      }

      const double inv_norm =
          theta / std::hypot(angle_axis[0], angle_axis[1], angle_axis[2]);
      for (double& angle_axi : angle_axis) {
        angle_axi *= inv_norm;
      }

      AngleAxisToRotationMatrix(angle_axis, R);
      rotation_matrix_rotated_p[0] = R[0] * p[0] + R[3] * p[1] + R[6] * p[2];
      rotation_matrix_rotated_p[1] = R[1] * p[0] + R[4] * p[1] + R[7] * p[2];
      rotation_matrix_rotated_p[2] = R[2] * p[0] + R[5] * p[1] + R[8] * p[2];

      AngleAxisRotatePoint(angle_axis, p, angle_axis_rotated_p);
      for (int k = 0; k < 3; ++k) {
        // clang-format off
        EXPECT_NEAR(rotation_matrix_rotated_p[k],
                    angle_axis_rotated_p[k],
                    kTolerance) << "p: " << p[0]
                                << " " << p[1]
                                << " " << p[2]
                                << " angle_axis: " << angle_axis[0]
                                << " " << angle_axis[1]
                                << " " << angle_axis[2];
        // clang-format on
      }
    }
  }
}

TEST(Quaternion, UnitQuaternion) {
  using Jet = ceres::Jet<double, 4>;
  std::array<Jet, 4> quaternion = {
      Jet(1.0, 0), Jet(0.0, 1), Jet(0.0, 2), Jet(0.0, 3)};
  std::array<Jet, 3> point = {Jet(0.0), Jet(0.0), Jet(0.0)};
  std::array<Jet, 3> rotated_point;
  QuaternionRotatePoint(quaternion.data(), point.data(), rotated_point.data());
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(rotated_point[i], point[i]);
    EXPECT_FALSE(rotated_point[i].v.array().isNaN().any());
  }
}

TEST(AngleAxis, NearZeroRotatePointGivesSameAnswerAsRotationMatrix) {
  std::mt19937 prng;
  std::uniform_real_distribution<double> uniform_distribution{-1.0, 1.0};
  double angle_axis[3];
  double R[9];
  double p[3];
  double angle_axis_rotated_p[3];
  double rotation_matrix_rotated_p[3];

  for (int i = 0; i < 10000; ++i) {
    double norm2 = 0.0;
    for (int k = 0; k < 3; ++k) {
      angle_axis[k] = uniform_distribution(prng);
      p[k] = uniform_distribution(prng);
      norm2 = angle_axis[k] * angle_axis[k];
    }

    double theta = (2.0 * i * 0.0001 - 1.0) * 1e-16;
    const double inv_norm = theta / sqrt(norm2);
    for (double& angle_axi : angle_axis) {
      angle_axi *= inv_norm;
    }

    AngleAxisToRotationMatrix(angle_axis, R);
    rotation_matrix_rotated_p[0] = R[0] * p[0] + R[3] * p[1] + R[6] * p[2];
    rotation_matrix_rotated_p[1] = R[1] * p[0] + R[4] * p[1] + R[7] * p[2];
    rotation_matrix_rotated_p[2] = R[2] * p[0] + R[5] * p[1] + R[8] * p[2];

    AngleAxisRotatePoint(angle_axis, p, angle_axis_rotated_p);
    for (int k = 0; k < 3; ++k) {
      // clang-format off
      EXPECT_NEAR(rotation_matrix_rotated_p[k],
                  angle_axis_rotated_p[k],
                  kTolerance) << "p: " << p[0]
                              << " " << p[1]
                              << " " << p[2]
                              << " angle_axis: " << angle_axis[0]
                              << " " << angle_axis[1]
                              << " " << angle_axis[2];
      // clang-format on
    }
  }
}

TEST(MatrixAdapter, RowMajor3x3ReturnTypeAndAccessIsCorrect) {
  double array[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  const float const_array[9] = {
      1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f};
  MatrixAdapter<double, 3, 1> A = RowMajorAdapter3x3(array);
  MatrixAdapter<const float, 3, 1> B = RowMajorAdapter3x3(const_array);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      // The values are integers from 1 to 9, so equality tests are appropriate
      // even for float and double values.
      EXPECT_EQ(A(i, j), array[3 * i + j]);
      EXPECT_EQ(B(i, j), const_array[3 * i + j]);
    }
  }
}

TEST(MatrixAdapter, ColumnMajor3x3ReturnTypeAndAccessIsCorrect) {
  double array[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  const float const_array[9] = {
      1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f};
  MatrixAdapter<double, 1, 3> A = ColumnMajorAdapter3x3(array);
  MatrixAdapter<const float, 1, 3> B = ColumnMajorAdapter3x3(const_array);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      // The values are integers from 1 to 9, so equality tests are
      // appropriate even for float and double values.
      EXPECT_EQ(A(i, j), array[3 * j + i]);
      EXPECT_EQ(B(i, j), const_array[3 * j + i]);
    }
  }
}

TEST(MatrixAdapter, RowMajor2x4IsCorrect) {
  const int expected[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  int array[8];
  MatrixAdapter<int, 4, 1> M(array);
  // clang-format off
  M(0, 0) = 1; M(0, 1) = 2; M(0, 2) = 3; M(0, 3) = 4;
  M(1, 0) = 5; M(1, 1) = 6; M(1, 2) = 7; M(1, 3) = 8;
  // clang-format on
  for (int k = 0; k < 8; ++k) {
    EXPECT_EQ(array[k], expected[k]);
  }
}

TEST(MatrixAdapter, ColumnMajor2x4IsCorrect) {
  const int expected[8] = {1, 5, 2, 6, 3, 7, 4, 8};
  int array[8];
  MatrixAdapter<int, 1, 2> M(array);
  // clang-format off
  M(0, 0) = 1; M(0, 1) = 2; M(0, 2) = 3; M(0, 3) = 4;
  M(1, 0) = 5; M(1, 1) = 6; M(1, 2) = 7; M(1, 3) = 8;
  // clang-format on
  for (int k = 0; k < 8; ++k) {
    EXPECT_EQ(array[k], expected[k]);
  }
}

TEST(RotationMatrixToAngleAxis, NearPiExampleOneFromTobiasStrauss) {
  // Example from Tobias Strauss
  // clang-format off
  const double rotation_matrix[] = {
    -0.999807135425239,    -0.0128154391194470,   -0.0148814136745799,
    -0.0128154391194470,   -0.148441438622958,     0.988838158557669,
    -0.0148814136745799,    0.988838158557669,     0.148248574048196
  };
  // clang-format on

  double angle_axis[3];
  RotationMatrixToAngleAxis(RowMajorAdapter3x3(rotation_matrix), angle_axis);
  double round_trip[9];
  AngleAxisToRotationMatrix(angle_axis, RowMajorAdapter3x3(round_trip));
  EXPECT_THAT(rotation_matrix, IsNear3x3Matrix(round_trip));
}

static void CheckRotationMatrixToAngleAxisRoundTrip(const double theta,
                                                    const double phi,
                                                    const double angle) {
  double angle_axis[3];
  angle_axis[0] = angle * sin(phi) * cos(theta);
  angle_axis[1] = angle * sin(phi) * sin(theta);
  angle_axis[2] = angle * cos(phi);

  double rotation_matrix[9];
  AngleAxisToRotationMatrix(angle_axis, rotation_matrix);

  double angle_axis_round_trip[3];
  RotationMatrixToAngleAxis(rotation_matrix, angle_axis_round_trip);
  EXPECT_THAT(angle_axis_round_trip, IsNearAngleAxis(angle_axis));
}

TEST(RotationMatrixToAngleAxis, ExhaustiveRoundTrip) {
  constexpr double kMaxSmallAngle = 1e-8;
  std::mt19937 prng;
  std::uniform_real_distribution<double> uniform_distribution1{
      kPi - kMaxSmallAngle, kPi};
  std::uniform_real_distribution<double> uniform_distribution2{
      -1.0, 2.0 * kMaxSmallAngle - 1.0};
  const int kNumSteps = 1000;
  for (int i = 0; i < kNumSteps; ++i) {
    const double theta = static_cast<double>(i) / kNumSteps * 2.0 * kPi;
    for (int j = 0; j < kNumSteps; ++j) {
      const double phi = static_cast<double>(j) / kNumSteps * kPi;
      // Rotations of angle Pi.
      CheckRotationMatrixToAngleAxisRoundTrip(theta, phi, kPi);
      // Rotation of angle approximately Pi.
      CheckRotationMatrixToAngleAxisRoundTrip(
          theta, phi, uniform_distribution1(prng));
      // Rotations of angle approximately zero.
      CheckRotationMatrixToAngleAxisRoundTrip(
          theta, phi, uniform_distribution2(prng));
    }
  }
}

}  // namespace internal
}  // namespace ceres
