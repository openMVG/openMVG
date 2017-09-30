// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
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
// Author: keir@google.com (Keir Mierle)

#include "ceres/jet.h"

#include <algorithm>
#include <cmath>

#include "glog/logging.h"
#include "gtest/gtest.h"
#include "ceres/fpclassify.h"
#include "ceres/stringprintf.h"
#include "ceres/test_util.h"

#define VL VLOG(1)

namespace ceres {
namespace internal {

const double kE = 2.71828182845904523536;

typedef Jet<double, 2> J;

// Convenient shorthand for making a jet.
J MakeJet(double a, double v0, double v1) {
  J z;
  z.a = a;
  z.v[0] = v0;
  z.v[1] = v1;
  return z;
}

// On a 32-bit optimized build, the mismatch is about 1.4e-14.
double const kTolerance = 1e-13;

void ExpectJetsClose(const J &x, const J &y) {
  ExpectClose(x.a, y.a, kTolerance);
  ExpectClose(x.v[0], y.v[0], kTolerance);
  ExpectClose(x.v[1], y.v[1], kTolerance);
}

TEST(Jet, Jet) {
  // Pick arbitrary values for x and y.
  J x = MakeJet(2.3, -2.7, 1e-3);
  J y = MakeJet(1.7,  0.5, 1e+2);

  VL << "x = " << x;
  VL << "y = " << y;

  { // Check that log(exp(x)) == x.
    J z = exp(x);
    J w = log(z);
    VL << "z = " << z;
    VL << "w = " << w;
    ExpectJetsClose(w, x);
  }

  { // Check that (x * y) / x == y.
    J z = x * y;
    J w = z / x;
    VL << "z = " << z;
    VL << "w = " << w;
    ExpectJetsClose(w, y);
  }

  { // Check that sqrt(x * x) == x.
    J z = x * x;
    J w = sqrt(z);
    VL << "z = " << z;
    VL << "w = " << w;
    ExpectJetsClose(w, x);
  }

  { // Check that sqrt(y) * sqrt(y) == y.
    J z = sqrt(y);
    J w = z * z;
    VL << "z = " << z;
    VL << "w = " << w;
    ExpectJetsClose(w, y);
  }

  { // Check that cos(2*x) = cos(x)^2 - sin(x)^2
    J z = cos(J(2.0) * x);
    J w = cos(x)*cos(x) - sin(x)*sin(x);
    VL << "z = " << z;
    VL << "w = " << w;
    ExpectJetsClose(w, z);
  }

  { // Check that sin(2*x) = 2*cos(x)*sin(x)
    J z = sin(J(2.0) * x);
    J w = J(2.0)*cos(x)*sin(x);
    VL << "z = " << z;
    VL << "w = " << w;
    ExpectJetsClose(w, z);
  }

  { // Check that cos(x)*cos(x) + sin(x)*sin(x) = 1
    J z = cos(x) * cos(x);
    J w = sin(x) * sin(x);
    VL << "z = " << z;
    VL << "w = " << w;
    ExpectJetsClose(z + w, J(1.0));
  }

  { // Check that atan2(r*sin(t), r*cos(t)) = t.
    J t = MakeJet(0.7, -0.3, +1.5);
    J r = MakeJet(2.3, 0.13, -2.4);
    VL << "t = " << t;
    VL << "r = " << r;

    J u = atan2(r * sin(t), r * cos(t));
    VL << "u = " << u;

    ExpectJetsClose(u, t);
  }

  { // Check that tan(x) = sin(x) / cos(x).
    J z = tan(x);
    J w = sin(x) / cos(x);
    VL << "z = " << z;
    VL << "w = " << w;
    ExpectJetsClose(z, w);
  }

  { // Check that tan(atan(x)) = x.
    J z = tan(atan(x));
    J w = x;
    VL << "z = " << z;
    VL << "w = " << w;
    ExpectJetsClose(z, w);
  }

  { // Check that cosh(x)*cosh(x) - sinh(x)*sinh(x) = 1
    J z = cosh(x) * cosh(x);
    J w = sinh(x) * sinh(x);
    VL << "z = " << z;
    VL << "w = " << w;
    ExpectJetsClose(z - w, J(1.0));
  }

  { // Check that tanh(x + y) = (tanh(x) + tanh(y)) / (1 + tanh(x) tanh(y))
    J z = tanh(x + y);
    J w = (tanh(x) + tanh(y)) / (J(1.0) + tanh(x) * tanh(y));
    VL << "z = " << z;
    VL << "w = " << w;
    ExpectJetsClose(z, w);
  }

  { // Check that pow(x, 1) == x.
    VL << "x = " << x;

    J u = pow(x, 1.);
    VL << "u = " << u;

    ExpectJetsClose(x, u);
  }

  { // Check that pow(x, 1) == x.
    J y = MakeJet(1, 0.0, 0.0);
    VL << "x = " << x;
    VL << "y = " << y;

    J u = pow(x, y);
    VL << "u = " << u;

    ExpectJetsClose(x, u);
  }

  { // Check that pow(e, log(x)) == x.
    J logx = log(x);

    VL << "x = " << x;
    VL << "y = " << y;

    J u = pow(kE, logx);
    VL << "u = " << u;

    ExpectJetsClose(x, u);
  }

  { // Check that pow(e, log(x)) == x.
    J logx = log(x);
    J e = MakeJet(kE, 0., 0.);
    VL << "x = " << x;
    VL << "log(x) = " << logx;

    J u = pow(e, logx);
    VL << "u = " << u;

    ExpectJetsClose(x, u);
  }

  { // Check that pow(e, log(x)) == x.
    J logx = log(x);
    J e = MakeJet(kE, 0., 0.);
    VL << "x = " << x;
    VL << "logx = " << logx;

    J u = pow(e, logx);
    VL << "u = " << u;

    ExpectJetsClose(x, u);
  }

  { // Check that pow(x,y) = exp(y*log(x)).
    J logx = log(x);
    J e = MakeJet(kE, 0., 0.);
    VL << "x = " << x;
    VL << "logx = " << logx;

    J u = pow(e, y*logx);
    J v = pow(x, y);
    VL << "u = " << u;
    VL << "v = " << v;

    ExpectJetsClose(v, u);
  }

  { // Check that pow(0, y) == 0 for y > 1, with both arguments Jets.
    // This tests special case handling inside pow().
    J a = MakeJet(0, 1, 2);
    J b = MakeJet(2, 3, 4);
    VL << "a = " << a;
    VL << "b = " << b;

    J c = pow(a, b);
    VL << "a^b = " << c;
    ExpectJetsClose(c, MakeJet(0, 0, 0));
  }

  { // Check that pow(0, y) == 0 for y == 1, with both arguments Jets.
    // This tests special case handling inside pow().
    J a = MakeJet(0, 1, 2);
    J b = MakeJet(1, 3, 4);
    VL << "a = " << a;
    VL << "b = " << b;

    J c = pow(a, b);
    VL << "a^b = " << c;
    ExpectJetsClose(c, MakeJet(0, 1, 2));
  }

  { // Check that pow(0, <1) is not finite, with both arguments Jets.
    for (int i = 1; i < 10; i++) {
      J a = MakeJet(0, 1, 2);
      J b = MakeJet(i*0.1, 3, 4);       // b = 0.1 ... 0.9
      VL << "a = " << a;
      VL << "b = " << b;

      J c = pow(a, b);
      VL << "a^b = " << c;
      EXPECT_EQ(c.a, 0.0);
      EXPECT_FALSE(IsFinite(c.v[0]));
      EXPECT_FALSE(IsFinite(c.v[1]));
    }
    for (int i = -10; i < 0; i++) {
      J a = MakeJet(0, 1, 2);
      J b = MakeJet(i*0.1, 3, 4);       // b = -1,-0.9 ... -0.1
      VL << "a = " << a;
      VL << "b = " << b;

      J c = pow(a, b);
      VL << "a^b = " << c;
      EXPECT_FALSE(IsFinite(c.a));
      EXPECT_FALSE(IsFinite(c.v[0]));
      EXPECT_FALSE(IsFinite(c.v[1]));
    }

    {
      // The special case of 0^0 = 1 defined by the C standard.
      J a = MakeJet(0, 1, 2);
      J b = MakeJet(0, 3, 4);
      VL << "a = " << a;
      VL << "b = " << b;

      J c = pow(a, b);
      VL << "a^b = " << c;
      EXPECT_EQ(c.a, 1.0);
      EXPECT_FALSE(IsFinite(c.v[0]));
      EXPECT_FALSE(IsFinite(c.v[1]));
    }
  }

  { // Check that pow(<0, b) is correct for integer b.
    // This tests special case handling inside pow().
    J a = MakeJet(-1.5, 3, 4);

    // b integer:
    for (int i = -10; i <= 10; i++) {
      J b = MakeJet(i, 0, 5);
      VL << "a = " << a;
      VL << "b = " << b;

      J c = pow(a, b);
      VL << "a^b = " << c;
      ExpectClose(c.a, pow(-1.5, i), kTolerance);
      EXPECT_TRUE(IsFinite(c.v[0]));
      EXPECT_FALSE(IsFinite(c.v[1]));
      ExpectClose(c.v[0], i * pow(-1.5, i - 1) * 3.0, kTolerance);
    }
  }

  { // Check that pow(<0, b) is correct for noninteger b.
    // This tests special case handling inside pow().
    J a = MakeJet(-1.5, 3, 4);
    J b = MakeJet(-2.5, 0, 5);
    VL << "a = " << a;
    VL << "b = " << b;

    J c = pow(a, b);
    VL << "a^b = " << c;
    EXPECT_FALSE(IsFinite(c.a));
    EXPECT_FALSE(IsFinite(c.v[0]));
    EXPECT_FALSE(IsFinite(c.v[1]));
  }

  {
    // Check that pow(0,y) == 0 for y == 2, with the second argument a
    // Jet.  This tests special case handling inside pow().
    double a = 0;
    J b = MakeJet(2, 3, 4);
    VL << "a = " << a;
    VL << "b = " << b;

    J c = pow(a, b);
    VL << "a^b = " << c;
    ExpectJetsClose(c, MakeJet(0, 0, 0));
  }

  {
    // Check that pow(<0,y) is correct for integer y. This tests special case
    // handling inside pow().
    double a = -1.5;
    for (int i = -10; i <= 10; i++) {
      J b = MakeJet(i, 3, 0);
      VL << "a = " << a;
      VL << "b = " << b;

      J c = pow(a, b);
      VL << "a^b = " << c;
      ExpectClose(c.a, pow(-1.5, i), kTolerance);
      EXPECT_FALSE(IsFinite(c.v[0]));
      EXPECT_TRUE(IsFinite(c.v[1]));
      ExpectClose(c.v[1], 0, kTolerance);
    }
  }

  {
    // Check that pow(<0,y) is correct for noninteger y. This tests special
    // case handling inside pow().
    double a = -1.5;
    J b = MakeJet(-3.14, 3, 0);
    VL << "a = " << a;
    VL << "b = " << b;

    J c = pow(a, b);
    VL << "a^b = " << c;
    EXPECT_FALSE(IsFinite(c.a));
    EXPECT_FALSE(IsFinite(c.v[0]));
    EXPECT_FALSE(IsFinite(c.v[1]));
  }

  { // Check that 1 + x == x + 1.
    J a = x + 1.0;
    J b = 1.0 + x;
    J c = x;
    c += 1.0;

    ExpectJetsClose(a, b);
    ExpectJetsClose(a, c);
  }

  { // Check that 1 - x == -(x - 1).
    J a = 1.0 - x;
    J b = -(x - 1.0);
    J c = x;
    c -= 1.0;

    ExpectJetsClose(a, b);
    ExpectJetsClose(a, -c);
  }

  { // Check that (x/s)*s == (x*s)/s.
    J a = x / 5.0;
    J b = x * 5.0;
    J c = x;
    c /= 5.0;
    J d = x;
    d *= 5.0;

    ExpectJetsClose(5.0 * a, b / 5.0);
    ExpectJetsClose(a, c);
    ExpectJetsClose(b, d);
  }

  { // Check that x / y == 1 / (y / x).
    J a = x / y;
    J b = 1.0 / (y / x);
    VL << "a = " << a;
    VL << "b = " << b;

    ExpectJetsClose(a, b);
  }

  { // Check that abs(-x * x) == sqrt(x * x).
    ExpectJetsClose(abs(-x), sqrt(x * x));
  }

  { // Check that cos(acos(x)) == x.
    J a = MakeJet(0.1, -2.7, 1e-3);
    ExpectJetsClose(cos(acos(a)), a);
    ExpectJetsClose(acos(cos(a)), a);

    J b = MakeJet(0.6,  0.5, 1e+2);
    ExpectJetsClose(cos(acos(b)), b);
    ExpectJetsClose(acos(cos(b)), b);
  }

  { // Check that sin(asin(x)) == x.
    J a = MakeJet(0.1, -2.7, 1e-3);
    ExpectJetsClose(sin(asin(a)), a);
    ExpectJetsClose(asin(sin(a)), a);

    J b = MakeJet(0.4,  0.5, 1e+2);
    ExpectJetsClose(sin(asin(b)), b);
    ExpectJetsClose(asin(sin(b)), b);
  }

  {
    J zero = J(0.0);

    // Check that J0(0) == 1.
    ExpectJetsClose(BesselJ0(zero), J(1.0));

    // Check that J1(0) == 0.
    ExpectJetsClose(BesselJ1(zero), zero);

    // Check that J2(0) == 0.
    ExpectJetsClose(BesselJn(2, zero), zero);

    // Check that J3(0) == 0.
    ExpectJetsClose(BesselJn(3, zero), zero);

    J z = MakeJet(0.1, -2.7, 1e-3);

    // Check that J0(z) == Jn(0,z).
    ExpectJetsClose(BesselJ0(z), BesselJn(0, z));

    // Check that J1(z) == Jn(1,z).
    ExpectJetsClose(BesselJ1(z), BesselJn(1, z));

    // Check that J0(z)+J2(z) == (2/z)*J1(z).
    // See formula http://dlmf.nist.gov/10.6.E1
    ExpectJetsClose(BesselJ0(z) + BesselJn(2, z), (2.0 / z) * BesselJ1(z));
  }

  { // Check that floor of a positive number works.
    J a = MakeJet(0.1, -2.7, 1e-3);
    J b = floor(a);
    J expected = MakeJet(floor(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }

  { // Check that floor of a negative number works.
    J a = MakeJet(-1.1, -2.7, 1e-3);
    J b = floor(a);
    J expected = MakeJet(floor(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }

  { // Check that floor of a positive number works.
    J a = MakeJet(10.123, -2.7, 1e-3);
    J b = floor(a);
    J expected = MakeJet(floor(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }

  { // Check that ceil of a positive number works.
    J a = MakeJet(0.1, -2.7, 1e-3);
    J b = ceil(a);
    J expected = MakeJet(ceil(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }

  { // Check that ceil of a negative number works.
    J a = MakeJet(-1.1, -2.7, 1e-3);
    J b = ceil(a);
    J expected = MakeJet(ceil(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }

  { // Check that ceil of a positive number works.
    J a = MakeJet(10.123, -2.7, 1e-3);
    J b = ceil(a);
    J expected = MakeJet(ceil(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }
}

TEST(Jet, JetsInEigenMatrices) {
  J x = MakeJet(2.3, -2.7, 1e-3);
  J y = MakeJet(1.7,  0.5, 1e+2);
  J z = MakeJet(5.3, -4.7, 1e-3);
  J w = MakeJet(9.7,  1.5, 10.1);

  Eigen::Matrix<J, 2, 2> M;
  Eigen::Matrix<J, 2, 1> v, r1, r2;

  M << x, y, z, w;
  v << x, z;

  // Check that M * v == (v^T * M^T)^T
  r1 = M * v;
  r2 = (v.transpose() * M.transpose()).transpose();

  ExpectJetsClose(r1(0), r2(0));
  ExpectJetsClose(r1(1), r2(1));
}

TEST(JetTraitsTest, ClassificationMixed) {
  Jet<double, 3> a(5.5, 0);
  a.v[0] = std::numeric_limits<double>::quiet_NaN();
  a.v[1] = std::numeric_limits<double>::infinity();
  a.v[2] = -std::numeric_limits<double>::infinity();
  EXPECT_FALSE(IsFinite(a));
  EXPECT_FALSE(IsNormal(a));
  EXPECT_TRUE(IsInfinite(a));
  EXPECT_TRUE(IsNaN(a));
}

TEST(JetTraitsTest, ClassificationNaN) {
  Jet<double, 3> a(5.5, 0);
  a.v[0] = std::numeric_limits<double>::quiet_NaN();
  a.v[1] = 0.0;
  a.v[2] = 0.0;
  EXPECT_FALSE(IsFinite(a));
  EXPECT_FALSE(IsNormal(a));
  EXPECT_FALSE(IsInfinite(a));
  EXPECT_TRUE(IsNaN(a));
}

TEST(JetTraitsTest, ClassificationInf) {
  Jet<double, 3> a(5.5, 0);
  a.v[0] = std::numeric_limits<double>::infinity();
  a.v[1] = 0.0;
  a.v[2] = 0.0;
  EXPECT_FALSE(IsFinite(a));
  EXPECT_FALSE(IsNormal(a));
  EXPECT_TRUE(IsInfinite(a));
  EXPECT_FALSE(IsNaN(a));
}

TEST(JetTraitsTest, ClassificationFinite) {
  Jet<double, 3> a(5.5, 0);
  a.v[0] = 100.0;
  a.v[1] = 1.0;
  a.v[2] = 3.14159;
  EXPECT_TRUE(IsFinite(a));
  EXPECT_TRUE(IsNormal(a));
  EXPECT_FALSE(IsInfinite(a));
  EXPECT_FALSE(IsNaN(a));
}

// ScalarBinaryOpTraits is only supported on Eigen versions >= 3.3
#if EIGEN_VERSION_AT_LEAST(3, 3, 0)
TEST(JetTraitsTest, MatrixScalarUnaryOps) {
  const J x = MakeJet(2.3, -2.7, 1e-3);
  const J y = MakeJet(1.7,  0.5, 1e+2);
  Eigen::Matrix<J, 2, 1> a;
  a << x, y;

  const J sum = a.sum();
  const J sum2 = a(0) + a(1);
  ExpectJetsClose(sum, sum2);
}

TEST(JetTraitsTest, MatrixScalarBinaryOps) {
  const J x = MakeJet(2.3, -2.7, 1e-3);
  const J y = MakeJet(1.7,  0.5, 1e+2);
  const J z = MakeJet(5.3, -4.7, 1e-3);
  const J w = MakeJet(9.7,  1.5, 10.1);

  Eigen::Matrix<J, 2, 2> M;
  Eigen::Vector2d v;

  M << x, y, z, w;
  v << 0.6, -2.1;

  // Check that M * v == M * v.cast<J>().
  const Eigen::Matrix<J, 2, 1> r1 = M * v;
  const Eigen::Matrix<J, 2, 1> r2 = M * v.cast<J>();

  ExpectJetsClose(r1(0), r2(0));
  ExpectJetsClose(r1(1), r2(1));

  // Check that M * a == M * T(a).
  const double a = 3.1;
  const Eigen::Matrix<J, 2, 2> r3 = M * a;
  const Eigen::Matrix<J, 2, 2> r4 = M * J(a);

  ExpectJetsClose(r3(0, 0), r4(0, 0));
  ExpectJetsClose(r3(1, 0), r4(1, 0));
  ExpectJetsClose(r3(0, 1), r4(0, 1));
  ExpectJetsClose(r3(1, 1), r4(1, 1));
}

TEST(JetTraitsTest, ArrayScalarUnaryOps) {
  const J x = MakeJet(2.3, -2.7, 1e-3);
  const J y = MakeJet(1.7,  0.5, 1e+2);
  Eigen::Array<J, 2, 1> a;
  a << x, y;

  const J sum = a.sum();
  const J sum2 = a(0) + a(1);
  ExpectJetsClose(sum, sum2);
}

TEST(JetTraitsTest, ArrayScalarBinaryOps) {
  const J x = MakeJet(2.3, -2.7, 1e-3);
  const J y = MakeJet(1.7,  0.5, 1e+2);

  Eigen::Array<J, 2, 1> a;
  Eigen::Array2d b;

  a << x, y;
  b << 0.6, -2.1;

  // Check that a * b == a * b.cast<T>()
  const Eigen::Array<J, 2, 1> r1 = a * b;
  const Eigen::Array<J, 2, 1> r2 = a * b.cast<J>();

  ExpectJetsClose(r1(0), r2(0));
  ExpectJetsClose(r1(1), r2(1));

  // Check that a * c == a * T(c).
  const double c = 3.1;
  const Eigen::Array<J, 2, 1> r3 = a * c;
  const Eigen::Array<J, 2, 1> r4 = a * J(c);

  ExpectJetsClose(r3(0), r3(0));
  ExpectJetsClose(r4(1), r4(1));
}
#endif   // EIGEN_VERSION_AT_LEAST(3, 3, 0)

}  // namespace internal
}  // namespace ceres
