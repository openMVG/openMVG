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
// Author: keir@google.com (Keir Mierle)

#include "ceres/jet.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cfenv>
#include <cmath>

#include "ceres/stringprintf.h"
#include "ceres/test_util.h"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

// The floating-point environment access and modification is only meaningful
// with the following pragma.
#ifdef _MSC_VER
#pragma float_control(precise, on, push)
#pragma fenv_access(on)
#elif !(defined(__ARM_ARCH) && __ARM_ARCH >= 8) && !defined(__MINGW32__)
// NOTE: FENV_ACCESS cannot be set to ON when targeting arm(v8) and MinGW
#pragma STDC FENV_ACCESS ON
#else
#define CERES_NO_FENV_ACCESS
#endif

namespace ceres::internal {

namespace {

constexpr double kE = 2.71828182845904523536;

using J = Jet<double, 2>;
// Don't care about the dual part for scalar part categorization and comparison
// tests
template <typename T>
using J0 = Jet<T, 0>;
using J0d = J0<double>;

// Convenient shorthand for making a jet.
J MakeJet(double a, double v0, double v1) {
  J z;
  z.a = a;
  z.v[0] = v0;
  z.v[1] = v1;
  return z;
}

double const kTolerance = 1e-13;

// Stores the floating-point environment containing active floating-point
// exceptions, rounding mode, etc., and restores it upon destruction.
//
// Useful for avoiding side-effects.
class Fenv {
 public:
  Fenv() { std::fegetenv(&e); }
  ~Fenv() { std::fesetenv(&e); }

  Fenv(const Fenv&) = delete;
  Fenv& operator=(const Fenv&) = delete;

 private:
  std::fenv_t e;
};

bool AreAlmostEqual(double x, double y, double max_abs_relative_difference) {
  if (std::isnan(x) && std::isnan(y)) {
    return true;
  }

  if (std::isinf(x) && std::isinf(y)) {
    return (std::signbit(x) == std::signbit(y));
  }

  Fenv env;  // Do not leak floating-point exceptions to the caller
  double absolute_difference = std::abs(x - y);
  double relative_difference =
      absolute_difference / std::max(std::abs(x), std::abs(y));

  if (std::fpclassify(x) == FP_ZERO || std::fpclassify(y) == FP_ZERO) {
    // If x or y is exactly zero, then relative difference doesn't have any
    // meaning. Take the absolute difference instead.
    relative_difference = absolute_difference;
  }
  return std::islessequal(relative_difference, max_abs_relative_difference);
}

MATCHER_P(IsAlmostEqualTo, y, "") {
  const bool result = (AreAlmostEqual(arg.a, y.a, kTolerance) &&
                       AreAlmostEqual(arg.v[0], y.v[0], kTolerance) &&
                       AreAlmostEqual(arg.v[1], y.v[1], kTolerance));
  if (!result) {
    *result_listener << "\nexpected - actual : " << y - arg;
  }
  return result;
}

const double kStep = 1e-8;
const double kNumericalTolerance = 1e-6;  // Numeric derivation is quite inexact

// Differentiate using Jet and confirm results with numerical derivation.
template <typename Function>
void NumericalTest(const char* name, const Function& f, const double x) {
  const double exact_dx = f(MakeJet(x, 1.0, 0.0)).v[0];
  const double estimated_dx =
      (f(J(x + kStep)).a - f(J(x - kStep)).a) / (2.0 * kStep);
  VLOG(1) << name << "(" << x << "), exact dx: " << exact_dx
          << ", estimated dx: " << estimated_dx;
  ExpectClose(exact_dx, estimated_dx, kNumericalTolerance);
}

// Same as NumericalTest, but given a function taking two arguments.
template <typename Function>
void NumericalTest2(const char* name,
                    const Function& f,
                    const double x,
                    const double y) {
  const J exact_delta = f(MakeJet(x, 1.0, 0.0), MakeJet(y, 0.0, 1.0));
  const double exact_dx = exact_delta.v[0];
  const double exact_dy = exact_delta.v[1];

  // Sanity check - these should be equivalent:
  EXPECT_EQ(exact_dx, f(MakeJet(x, 1.0, 0.0), MakeJet(y, 0.0, 0.0)).v[0]);
  EXPECT_EQ(exact_dx, f(MakeJet(x, 0.0, 1.0), MakeJet(y, 0.0, 0.0)).v[1]);
  EXPECT_EQ(exact_dy, f(MakeJet(x, 0.0, 0.0), MakeJet(y, 1.0, 0.0)).v[0]);
  EXPECT_EQ(exact_dy, f(MakeJet(x, 0.0, 0.0), MakeJet(y, 0.0, 1.0)).v[1]);

  const double estimated_dx =
      (f(J(x + kStep), J(y)).a - f(J(x - kStep), J(y)).a) / (2.0 * kStep);
  const double estimated_dy =
      (f(J(x), J(y + kStep)).a - f(J(x), J(y - kStep)).a) / (2.0 * kStep);
  VLOG(1) << name << "(" << x << ", " << y << "), exact dx: " << exact_dx
          << ", estimated dx: " << estimated_dx;
  ExpectClose(exact_dx, estimated_dx, kNumericalTolerance);
  VLOG(1) << name << "(" << x << ", " << y << "), exact dy: " << exact_dy
          << ", estimated dy: " << estimated_dy;
  ExpectClose(exact_dy, estimated_dy, kNumericalTolerance);
}

}  // namespace

// Pick arbitrary values for x and y.
const J x = MakeJet(2.3, -2.7, 1e-3);
const J y = MakeJet(1.7, 0.5, 1e+2);
const J z = MakeJet(1e-6, 1e-4, 1e-2);

TEST(Jet, Elementary) {
  EXPECT_THAT((x * y) / x, IsAlmostEqualTo(y));
  EXPECT_THAT(sqrt(x * x), IsAlmostEqualTo(x));
  EXPECT_THAT(sqrt(y) * sqrt(y), IsAlmostEqualTo(y));

  NumericalTest("sqrt", sqrt<double, 2>, 0.00001);
  NumericalTest("sqrt", sqrt<double, 2>, 1.0);

  EXPECT_THAT(x + 1.0, IsAlmostEqualTo(1.0 + x));
  {
    J c = x;
    c += 1.0;
    EXPECT_THAT(c, IsAlmostEqualTo(1.0 + x));
  }

  EXPECT_THAT(-(x - 1.0), IsAlmostEqualTo(1.0 - x));
  {
    J c = x;
    c -= 1.0;
    EXPECT_THAT(c, IsAlmostEqualTo(x - 1.0));
  }

  EXPECT_THAT((x * 5.0) / 5.0, IsAlmostEqualTo((x / 5.0) * 5.0));
  EXPECT_THAT((x * 5.0) / 5.0, IsAlmostEqualTo(x));
  EXPECT_THAT((x / 5.0) * 5.0, IsAlmostEqualTo(x));

  {
    J c = x;
    c /= 5.0;
    J d = x;
    d *= 5.0;
    EXPECT_THAT(c, IsAlmostEqualTo(x / 5.0));
    EXPECT_THAT(d, IsAlmostEqualTo(5.0 * x));
  }

  EXPECT_THAT(1.0 / (y / x), IsAlmostEqualTo(x / y));
}

TEST(Jet, Trigonometric) {
  EXPECT_THAT(cos(2.0 * x), IsAlmostEqualTo(cos(x) * cos(x) - sin(x) * sin(x)));
  EXPECT_THAT(sin(2.0 * x), IsAlmostEqualTo(2.0 * sin(x) * cos(x)));
  EXPECT_THAT(sin(x) * sin(x) + cos(x) * cos(x), IsAlmostEqualTo(J(1.0)));

  {
    J t = MakeJet(0.7, -0.3, +1.5);
    J r = MakeJet(2.3, 0.13, -2.4);
    EXPECT_THAT(atan2(r * sin(t), r * cos(t)), IsAlmostEqualTo(t));
  }

  EXPECT_THAT(sin(x) / cos(x), IsAlmostEqualTo(tan(x)));
  EXPECT_THAT(tan(atan(x)), IsAlmostEqualTo(x));

  {
    J a = MakeJet(0.1, -2.7, 1e-3);
    EXPECT_THAT(cos(acos(a)), IsAlmostEqualTo(a));
    EXPECT_THAT(acos(cos(a)), IsAlmostEqualTo(a));

    J b = MakeJet(0.6, 0.5, 1e+2);
    EXPECT_THAT(cos(acos(b)), IsAlmostEqualTo(b));
    EXPECT_THAT(acos(cos(b)), IsAlmostEqualTo(b));
  }

  {
    J a = MakeJet(0.1, -2.7, 1e-3);
    EXPECT_THAT(sin(asin(a)), IsAlmostEqualTo(a));
    EXPECT_THAT(asin(sin(a)), IsAlmostEqualTo(a));

    J b = MakeJet(0.4, 0.5, 1e+2);
    EXPECT_THAT(sin(asin(b)), IsAlmostEqualTo(b));
    EXPECT_THAT(asin(sin(b)), IsAlmostEqualTo(b));
  }
}

TEST(Jet, Hyperbolic) {
  // cosh(x)*cosh(x) - sinh(x)*sinh(x) = 1
  EXPECT_THAT(cosh(x) * cosh(x) - sinh(x) * sinh(x), IsAlmostEqualTo(J(1.0)));

  // tanh(x + y) = (tanh(x) + tanh(y)) / (1 + tanh(x) tanh(y))
  EXPECT_THAT(
      tanh(x + y),
      IsAlmostEqualTo((tanh(x) + tanh(y)) / (J(1.0) + tanh(x) * tanh(y))));
}

TEST(Jet, Abs) {
  EXPECT_THAT(abs(-x * x), IsAlmostEqualTo(x * x));
  EXPECT_THAT(abs(-x), IsAlmostEqualTo(sqrt(x * x)));

  {
    J a = MakeJet(-std::numeric_limits<double>::quiet_NaN(), 2.0, 4.0);
    J b = abs(a);
    EXPECT_TRUE(std::signbit(b.v[0]));
    EXPECT_TRUE(std::signbit(b.v[1]));
  }
}

TEST(Jet, Bessel) {
  J zero = J(0.0);

  EXPECT_THAT(BesselJ0(zero), IsAlmostEqualTo(J(1.0)));
  EXPECT_THAT(BesselJ1(zero), IsAlmostEqualTo(zero));
  EXPECT_THAT(BesselJn(2, zero), IsAlmostEqualTo(zero));
  EXPECT_THAT(BesselJn(3, zero), IsAlmostEqualTo(zero));

  J z = MakeJet(0.1, -2.7, 1e-3);

  EXPECT_THAT(BesselJ0(z), IsAlmostEqualTo(BesselJn(0, z)));
  EXPECT_THAT(BesselJ1(z), IsAlmostEqualTo(BesselJn(1, z)));

  //  See formula http://dlmf.nist.gov/10.6.E1
  EXPECT_THAT(BesselJ0(z) + BesselJn(2, z),
              IsAlmostEqualTo((2.0 / z) * BesselJ1(z)));
}

TEST(Jet, Floor) {
  {  // floor of a positive number works.
    J a = MakeJet(0.1, -2.7, 1e-3);
    J b = floor(a);
    J expected = MakeJet(floor(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }

  {  // floor of a negative number works.
    J a = MakeJet(-1.1, -2.7, 1e-3);
    J b = floor(a);
    J expected = MakeJet(floor(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }

  {  // floor of a positive number works.
    J a = MakeJet(10.123, -2.7, 1e-3);
    J b = floor(a);
    J expected = MakeJet(floor(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }
}

TEST(Jet, Ceil) {
  {  // ceil of a positive number works.
    J a = MakeJet(0.1, -2.7, 1e-3);
    J b = ceil(a);
    J expected = MakeJet(ceil(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }

  {  // ceil of a negative number works.
    J a = MakeJet(-1.1, -2.7, 1e-3);
    J b = ceil(a);
    J expected = MakeJet(ceil(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }

  {  // ceil of a positive number works.
    J a = MakeJet(10.123, -2.7, 1e-3);
    J b = ceil(a);
    J expected = MakeJet(ceil(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }
}

TEST(Jet, Erf) {
  {  // erf works.
    J a = MakeJet(10.123, -2.7, 1e-3);
    J b = erf(a);
    J expected = MakeJet(erf(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }
  NumericalTest("erf", erf<double, 2>, -1.0);
  NumericalTest("erf", erf<double, 2>, 1e-5);
  NumericalTest("erf", erf<double, 2>, 0.5);
  NumericalTest("erf", erf<double, 2>, 100.0);
}

TEST(Jet, Erfc) {
  {  // erfc works.
    J a = MakeJet(10.123, -2.7, 1e-3);
    J b = erfc(a);
    J expected = MakeJet(erfc(a.a), 0.0, 0.0);
    EXPECT_EQ(expected, b);
  }
  NumericalTest("erfc", erfc<double, 2>, -1.0);
  NumericalTest("erfc", erfc<double, 2>, 1e-5);
  NumericalTest("erfc", erfc<double, 2>, 0.5);
  NumericalTest("erfc", erfc<double, 2>, 100.0);
}

TEST(Jet, Cbrt) {
  EXPECT_THAT(cbrt(x * x * x), IsAlmostEqualTo(x));
  EXPECT_THAT(cbrt(y) * cbrt(y) * cbrt(y), IsAlmostEqualTo(y));
  EXPECT_THAT(cbrt(x), IsAlmostEqualTo(pow(x, 1.0 / 3.0)));

  NumericalTest("cbrt", cbrt<double, 2>, -1.0);
  NumericalTest("cbrt", cbrt<double, 2>, -1e-5);
  NumericalTest("cbrt", cbrt<double, 2>, 1e-5);
  NumericalTest("cbrt", cbrt<double, 2>, 1.0);
}

TEST(Jet, Log1p) {
  EXPECT_THAT(log1p(expm1(x)), IsAlmostEqualTo(x));
  EXPECT_THAT(log1p(x), IsAlmostEqualTo(log(J{1} + x)));

  {  // log1p(x) does not loose precision for small x
    J x = MakeJet(1e-16, 1e-8, 1e-4);
    EXPECT_THAT(log1p(x),
                IsAlmostEqualTo(MakeJet(9.9999999999999998e-17, 1e-8, 1e-4)));
    // log(1 + x) collapses to 0
    J v = log(J{1} + x);
    EXPECT_TRUE(v.a == 0);
  }
}

TEST(Jet, Expm1) {
  EXPECT_THAT(expm1(log1p(x)), IsAlmostEqualTo(x));
  EXPECT_THAT(expm1(x), IsAlmostEqualTo(exp(x) - 1.0));

  {  // expm1(x) does not loose precision for small x
    J x = MakeJet(9.9999999999999998e-17, 1e-8, 1e-4);
    EXPECT_THAT(expm1(x), IsAlmostEqualTo(MakeJet(1e-16, 1e-8, 1e-4)));
    // exp(x) - 1 collapses to 0
    J v = exp(x) - J{1};
    EXPECT_TRUE(v.a == 0);
  }
}

TEST(Jet, Exp2) {
  EXPECT_THAT(exp2(x), IsAlmostEqualTo(exp(x * log(2.0))));
  NumericalTest("exp2", exp2<double, 2>, -1.0);
  NumericalTest("exp2", exp2<double, 2>, -1e-5);
  NumericalTest("exp2", exp2<double, 2>, -1e-200);
  NumericalTest("exp2", exp2<double, 2>, 0.0);
  NumericalTest("exp2", exp2<double, 2>, 1e-200);
  NumericalTest("exp2", exp2<double, 2>, 1e-5);
  NumericalTest("exp2", exp2<double, 2>, 1.0);
}

TEST(Jet, Log) { EXPECT_THAT(log(exp(x)), IsAlmostEqualTo(x)); }

TEST(Jet, Log10) {
  EXPECT_THAT(log10(x), IsAlmostEqualTo(log(x) / log(10)));
  NumericalTest("log10", log10<double, 2>, 1e-5);
  NumericalTest("log10", log10<double, 2>, 1.0);
  NumericalTest("log10", log10<double, 2>, 98.76);
}

TEST(Jet, Log2) {
  EXPECT_THAT(log2(x), IsAlmostEqualTo(log(x) / log(2)));
  NumericalTest("log2", log2<double, 2>, 1e-5);
  NumericalTest("log2", log2<double, 2>, 1.0);
  NumericalTest("log2", log2<double, 2>, 100.0);
}

TEST(Jet, Norm) {
  EXPECT_THAT(norm(x), IsAlmostEqualTo(x * x));
  EXPECT_THAT(norm(-x), IsAlmostEqualTo(x * x));
}

TEST(Jet, Pow) {
  EXPECT_THAT(pow(x, 1.0), IsAlmostEqualTo(x));
  EXPECT_THAT(pow(x, MakeJet(1.0, 0.0, 0.0)), IsAlmostEqualTo(x));
  EXPECT_THAT(pow(kE, log(x)), IsAlmostEqualTo(x));
  EXPECT_THAT(pow(MakeJet(kE, 0., 0.), log(x)), IsAlmostEqualTo(x));
  EXPECT_THAT(pow(x, y),
              IsAlmostEqualTo(pow(MakeJet(kE, 0.0, 0.0), y * log(x))));

  // Specially cases

  // pow(0, y) == 0 for y > 1, with both arguments Jets.
  EXPECT_THAT(pow(MakeJet(0, 1, 2), MakeJet(2, 3, 4)),
              IsAlmostEqualTo(MakeJet(0, 0, 0)));

  // pow(0, y) == 0 for y == 1, with both arguments Jets.
  EXPECT_THAT(pow(MakeJet(0, 1, 2), MakeJet(1, 3, 4)),
              IsAlmostEqualTo(MakeJet(0, 1, 2)));

  // pow(0, <1) is not finite, with both arguments Jets.
  {
    for (int i = 1; i < 10; i++) {
      J a = MakeJet(0, 1, 2);
      J b = MakeJet(i * 0.1, 3, 4);  // b = 0.1 ... 0.9
      J c = pow(a, b);
      EXPECT_EQ(c.a, 0.0) << "\na: " << a << "\nb: " << b << "\na^b: " << c;
      EXPECT_FALSE(isfinite(c.v[0]))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
      EXPECT_FALSE(isfinite(c.v[1]))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
    }

    for (int i = -10; i < 0; i++) {
      J a = MakeJet(0, 1, 2);
      J b = MakeJet(i * 0.1, 3, 4);  // b = -1,-0.9 ... -0.1
      J c = pow(a, b);
      EXPECT_FALSE(isfinite(c.a))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
      EXPECT_FALSE(isfinite(c.v[0]))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
      EXPECT_FALSE(isfinite(c.v[1]))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
    }

    // The special case of 0^0 = 1 defined by the C standard.
    {
      J a = MakeJet(0, 1, 2);
      J b = MakeJet(0, 3, 4);
      J c = pow(a, b);
      EXPECT_EQ(c.a, 1.0) << "\na: " << a << "\nb: " << b << "\na^b: " << c;
      EXPECT_FALSE(isfinite(c.v[0]))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
      EXPECT_FALSE(isfinite(c.v[1]))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
    }
  }

  // pow(<0, b) is correct for integer b.
  {
    J a = MakeJet(-1.5, 3, 4);

    // b integer:
    for (int i = -10; i <= 10; i++) {
      J b = MakeJet(i, 0, 5);
      J c = pow(a, b);

      EXPECT_TRUE(AreAlmostEqual(c.a, pow(-1.5, i), kTolerance))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
      EXPECT_TRUE(isfinite(c.v[0]))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
      EXPECT_FALSE(isfinite(c.v[1]))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
      EXPECT_TRUE(
          AreAlmostEqual(c.v[0], i * pow(-1.5, i - 1) * 3.0, kTolerance))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
    }
  }

  // pow(<0, b) is correct for noninteger b.
  {
    J a = MakeJet(-1.5, 3, 4);
    J b = MakeJet(-2.5, 0, 5);
    J c = pow(a, b);
    EXPECT_FALSE(isfinite(c.a))
        << "\na: " << a << "\nb: " << b << "\na^b: " << c;
    EXPECT_FALSE(isfinite(c.v[0]))
        << "\na: " << a << "\nb: " << b << "\na^b: " << c;
    EXPECT_FALSE(isfinite(c.v[1]))
        << "\na: " << a << "\nb: " << b << "\na^b: " << c;
  }

  // pow(0,y) == 0 for y == 2, with the second argument a Jet.
  EXPECT_THAT(pow(0.0, MakeJet(2, 3, 4)), IsAlmostEqualTo(MakeJet(0, 0, 0)));

  // pow(<0,y) is correct for integer y.
  {
    double a = -1.5;
    for (int i = -10; i <= 10; i++) {
      J b = MakeJet(i, 3, 0);
      J c = pow(a, b);
      ExpectClose(c.a, pow(-1.5, i), kTolerance);
      EXPECT_FALSE(isfinite(c.v[0]))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
      EXPECT_TRUE(isfinite(c.v[1]))
          << "\na: " << a << "\nb: " << b << "\na^b: " << c;
      ExpectClose(c.v[1], 0, kTolerance);
    }
  }

  // pow(<0,y) is correct for noninteger y.
  {
    double a = -1.5;
    J b = MakeJet(-3.14, 3, 0);
    J c = pow(a, b);
    EXPECT_FALSE(isfinite(c.a))
        << "\na: " << a << "\nb: " << b << "\na^b: " << c;
    EXPECT_FALSE(isfinite(c.v[0]))
        << "\na: " << a << "\nb: " << b << "\na^b: " << c;
    EXPECT_FALSE(isfinite(c.v[1]))
        << "\na: " << a << "\nb: " << b << "\na^b: " << c;
  }
}

TEST(Jet, Hypot2) {
  // Resolve the ambiguity between two and three argument hypot overloads
  using Hypot2 = J(const J&, const J&);
  auto* const hypot2 = static_cast<Hypot2*>(&hypot<double, 2>);

  // clang-format off
  NumericalTest2("hypot2", hypot2,  0.0,   1e-5);
  NumericalTest2("hypot2", hypot2, -1e-5,  0.0);
  NumericalTest2("hypot2", hypot2,  1e-5,  1e-5);
  NumericalTest2("hypot2", hypot2,  0.0,   1.0);
  NumericalTest2("hypot2", hypot2,  1e-3,  1.0);
  NumericalTest2("hypot2", hypot2,  1e-3, -1.0);
  NumericalTest2("hypot2", hypot2, -1e-3,  1.0);
  NumericalTest2("hypot2", hypot2, -1e-3, -1.0);
  NumericalTest2("hypot2", hypot2,  1.0,   2.0);
  // clang-format on

  J zero = MakeJet(0.0, 2.0, 3.14);
  EXPECT_THAT(hypot(x, y), IsAlmostEqualTo(sqrt(x * x + y * y)));
  EXPECT_THAT(hypot(x, x), IsAlmostEqualTo(sqrt(2.0) * abs(x)));

  // The derivative is zero tangentially to the circle:
  EXPECT_THAT(hypot(MakeJet(2.0, 1.0, 1.0), MakeJet(2.0, 1.0, -1.0)),
              IsAlmostEqualTo(MakeJet(sqrt(8.0), std::sqrt(2.0), 0.0)));

  EXPECT_THAT(hypot(zero, x), IsAlmostEqualTo(x));
  EXPECT_THAT(hypot(y, zero), IsAlmostEqualTo(y));

  // hypot(x, 0, 0) == x, even when x * x underflows:
  EXPECT_EQ(
      std::numeric_limits<double>::min() * std::numeric_limits<double>::min(),
      0.0);  // Make sure it underflows
  J tiny = MakeJet(std::numeric_limits<double>::min(), 2.0, 3.14);
  EXPECT_THAT(hypot(tiny, J{0}), IsAlmostEqualTo(tiny));

  // hypot(x, 0, 0) == x, even when x * x overflows:
  EXPECT_EQ(
      std::numeric_limits<double>::max() * std::numeric_limits<double>::max(),
      std::numeric_limits<double>::infinity());
  J huge = MakeJet(std::numeric_limits<double>::max(), 2.0, 3.14);
  EXPECT_THAT(hypot(huge, J{0}), IsAlmostEqualTo(huge));
}

TEST(Jet, Hypot3) {
  J zero = MakeJet(0.0, 2.0, 3.14);

  // hypot(x, y, z) == sqrt(x^2 + y^2 + z^2)
  EXPECT_THAT(hypot(x, y, z), IsAlmostEqualTo(sqrt(x * x + y * y + z * z)));

  // hypot(x, x) == sqrt(3) * abs(x)
  EXPECT_THAT(hypot(x, x, x), IsAlmostEqualTo(sqrt(3.0) * abs(x)));

  // The derivative is zero tangentially to the circle:
  EXPECT_THAT(hypot(MakeJet(2.0, 1.0, 1.0),
                    MakeJet(2.0, 1.0, -1.0),
                    MakeJet(2.0, -1.0, 0.0)),
              IsAlmostEqualTo(MakeJet(sqrt(12.0), 1.0 / std::sqrt(3.0), 0.0)));

  EXPECT_THAT(hypot(x, zero, zero), IsAlmostEqualTo(x));
  EXPECT_THAT(hypot(zero, y, zero), IsAlmostEqualTo(y));
  EXPECT_THAT(hypot(zero, zero, z), IsAlmostEqualTo(z));
  EXPECT_THAT(hypot(x, y, z), IsAlmostEqualTo(hypot(hypot(x, y), z)));
  EXPECT_THAT(hypot(x, y, z), IsAlmostEqualTo(hypot(x, hypot(y, z))));

  // The following two tests are disabled because the three argument hypot is
  // broken in the libc++ shipped with CLANG as of January 2022.

#if !defined(_LIBCPP_VERSION)
  // hypot(x, 0, 0) == x, even when x * x underflows:
  EXPECT_EQ(
      std::numeric_limits<double>::min() * std::numeric_limits<double>::min(),
      0.0);  // Make sure it underflows
  J tiny = MakeJet(std::numeric_limits<double>::min(), 2.0, 3.14);
  EXPECT_THAT(hypot(tiny, J{0}, J{0}), IsAlmostEqualTo(tiny));

  // hypot(x, 0, 0) == x, even when x * x overflows:
  EXPECT_EQ(
      std::numeric_limits<double>::max() * std::numeric_limits<double>::max(),
      std::numeric_limits<double>::infinity());
  J huge = MakeJet(std::numeric_limits<double>::max(), 2.0, 3.14);
  EXPECT_THAT(hypot(huge, J{0}, J{0}), IsAlmostEqualTo(huge));
#endif
}

#ifdef CERES_HAS_CPP20

TEST(Jet, Lerp) {
  EXPECT_THAT(lerp(x, y, J{0}), IsAlmostEqualTo(x));
  EXPECT_THAT(lerp(x, y, J{1}), IsAlmostEqualTo(y));
  EXPECT_THAT(lerp(x, x, J{1}), IsAlmostEqualTo(x));
  EXPECT_THAT(lerp(y, y, J{0}), IsAlmostEqualTo(y));
  EXPECT_THAT(lerp(x, y, J{0.5}), IsAlmostEqualTo((x + y) / J{2.0}));
  EXPECT_THAT(lerp(x, y, J{2}), IsAlmostEqualTo(J{2.0} * y - x));
  EXPECT_THAT(lerp(x, y, J{-2}), IsAlmostEqualTo(J{3.0} * x - J{2} * y));
}

TEST(Jet, Midpoint) {
  EXPECT_THAT(midpoint(x, y), IsAlmostEqualTo((x + y) / J{2}));
  EXPECT_THAT(midpoint(x, x), IsAlmostEqualTo(x));

  {
    // midpoint(x, y) = (x + y) / 2 while avoiding overflow
    J x = MakeJet(std::numeric_limits<double>::min(), 1, 2);
    J y = MakeJet(std::numeric_limits<double>::max(), 3, 4);
    EXPECT_THAT(midpoint(x, y), IsAlmostEqualTo(x + (y - x) / J{2}));
  }

  {
    // midpoint(x, x) = x while avoiding overflow
    J x = MakeJet(std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max());
    EXPECT_THAT(midpoint(x, x), IsAlmostEqualTo(x));
  }

  {  // midpoint does not overflow for very large values
    constexpr double a = 0.75 * std::numeric_limits<double>::max();
    J x = MakeJet(a, a, -a);
    J y = MakeJet(a, a, a);
    EXPECT_THAT(midpoint(x, y), IsAlmostEqualTo(MakeJet(a, a, 0)));
  }
}

#endif  // defined(CERES_HAS_CPP20)

TEST(Jet, Fma) {
  J v = fma(x, y, z);
  J w = x * y + z;
  EXPECT_THAT(v, IsAlmostEqualTo(w));
}

TEST(Jet, FmaxJetWithJet) {
  Fenv env;
  // Clear all exceptions to ensure none are set by the following function
  // calls.
  std::feclearexcept(FE_ALL_EXCEPT);

  EXPECT_THAT(fmax(x, y), IsAlmostEqualTo(x));
  EXPECT_THAT(fmax(y, x), IsAlmostEqualTo(x));

  // Average the Jets on equality (of scalar parts).
  const J scalar_part_only_equal_to_x = J(x.a, 2 * x.v);
  const J average = (x + scalar_part_only_equal_to_x) * 0.5;
  EXPECT_THAT(fmax(x, scalar_part_only_equal_to_x), IsAlmostEqualTo(average));
  EXPECT_THAT(fmax(scalar_part_only_equal_to_x, x), IsAlmostEqualTo(average));

  // Follow convention of fmax(): treat NANs as missing values.
  const J nan_scalar_part(std::numeric_limits<double>::quiet_NaN(), 2 * x.v);
  EXPECT_THAT(fmax(x, nan_scalar_part), IsAlmostEqualTo(x));
  EXPECT_THAT(fmax(nan_scalar_part, x), IsAlmostEqualTo(x));

#ifndef CERES_NO_FENV_ACCESS
  EXPECT_EQ(std::fetestexcept(FE_ALL_EXCEPT & ~FE_INEXACT), 0);
#endif
}

TEST(Jet, FmaxJetWithScalar) {
  Fenv env;
  // Clear all exceptions to ensure none are set by the following function
  // calls.
  std::feclearexcept(FE_ALL_EXCEPT);

  EXPECT_THAT(fmax(x, y.a), IsAlmostEqualTo(x));
  EXPECT_THAT(fmax(y.a, x), IsAlmostEqualTo(x));
  EXPECT_THAT(fmax(y, x.a), IsAlmostEqualTo(J{x.a}));
  EXPECT_THAT(fmax(x.a, y), IsAlmostEqualTo(J{x.a}));

  // Average the Jet and scalar cast to a Jet on equality (of scalar parts).
  const J average = (x + J{x.a}) * 0.5;
  EXPECT_THAT(fmax(x, x.a), IsAlmostEqualTo(average));
  EXPECT_THAT(fmax(x.a, x), IsAlmostEqualTo(average));

  // Follow convention of fmax(): treat NANs as missing values.
  EXPECT_THAT(fmax(x, std::numeric_limits<double>::quiet_NaN()),
              IsAlmostEqualTo(x));
  EXPECT_THAT(fmax(std::numeric_limits<double>::quiet_NaN(), x),
              IsAlmostEqualTo(x));
  const J nan_scalar_part(std::numeric_limits<double>::quiet_NaN(), 2 * x.v);
  EXPECT_THAT(fmax(nan_scalar_part, x.a), IsAlmostEqualTo(J{x.a}));
  EXPECT_THAT(fmax(x.a, nan_scalar_part), IsAlmostEqualTo(J{x.a}));

#ifndef CERES_NO_FENV_ACCESS
  EXPECT_EQ(std::fetestexcept(FE_ALL_EXCEPT & ~FE_INEXACT), 0);
#endif
}

TEST(Jet, FminJetWithJet) {
  Fenv env;
  // Clear all exceptions to ensure none are set by the following function
  // calls.
  std::feclearexcept(FE_ALL_EXCEPT);

  EXPECT_THAT(fmin(x, y), IsAlmostEqualTo(y));
  EXPECT_THAT(fmin(y, x), IsAlmostEqualTo(y));

  // Average the Jets on equality (of scalar parts).
  const J scalar_part_only_equal_to_x = J(x.a, 2 * x.v);
  const J average = (x + scalar_part_only_equal_to_x) * 0.5;
  EXPECT_THAT(fmin(x, scalar_part_only_equal_to_x), IsAlmostEqualTo(average));
  EXPECT_THAT(fmin(scalar_part_only_equal_to_x, x), IsAlmostEqualTo(average));

  // Follow convention of fmin(): treat NANs as missing values.
  const J nan_scalar_part(std::numeric_limits<double>::quiet_NaN(), 2 * x.v);
  EXPECT_THAT(fmin(x, nan_scalar_part), IsAlmostEqualTo(x));
  EXPECT_THAT(fmin(nan_scalar_part, x), IsAlmostEqualTo(x));

#ifndef CERES_NO_FENV_ACCESS
  EXPECT_EQ(std::fetestexcept(FE_ALL_EXCEPT & ~FE_INEXACT), 0);
#endif
}

TEST(Jet, FminJetWithScalar) {
  Fenv env;
  // Clear all exceptions to ensure none are set by the following function
  // calls.
  std::feclearexcept(FE_ALL_EXCEPT);

  EXPECT_THAT(fmin(x, y.a), IsAlmostEqualTo(J{y.a}));
  EXPECT_THAT(fmin(y.a, x), IsAlmostEqualTo(J{y.a}));
  EXPECT_THAT(fmin(y, x.a), IsAlmostEqualTo(y));
  EXPECT_THAT(fmin(x.a, y), IsAlmostEqualTo(y));

  // Average the Jet and scalar cast to a Jet on equality (of scalar parts).
  const J average = (x + J{x.a}) * 0.5;
  EXPECT_THAT(fmin(x, x.a), IsAlmostEqualTo(average));
  EXPECT_THAT(fmin(x.a, x), IsAlmostEqualTo(average));

  // Follow convention of fmin(): treat NANs as missing values.
  EXPECT_THAT(fmin(x, std::numeric_limits<double>::quiet_NaN()),
              IsAlmostEqualTo(x));
  EXPECT_THAT(fmin(std::numeric_limits<double>::quiet_NaN(), x),
              IsAlmostEqualTo(x));
  const J nan_scalar_part(std::numeric_limits<double>::quiet_NaN(), 2 * x.v);
  EXPECT_THAT(fmin(nan_scalar_part, x.a), IsAlmostEqualTo(J{x.a}));
  EXPECT_THAT(fmin(x.a, nan_scalar_part), IsAlmostEqualTo(J{x.a}));

#ifndef CERES_NO_FENV_ACCESS
  EXPECT_EQ(std::fetestexcept(FE_ALL_EXCEPT & ~FE_INEXACT), 0);
#endif
}

TEST(Jet, Fdim) {
  Fenv env;
  // Clear all exceptions to ensure none are set by the following function
  // calls.
  std::feclearexcept(FE_ALL_EXCEPT);

  const J zero{};
  const J diff = x - y;
  const J diffx = x - J{y.a};
  const J diffy = J{x.a} - y;

  EXPECT_THAT(fdim(x, y), IsAlmostEqualTo(diff));
  EXPECT_THAT(fdim(y, x), IsAlmostEqualTo(zero));
  EXPECT_THAT(fdim(x, y.a), IsAlmostEqualTo(diffx));
  EXPECT_THAT(fdim(y.a, x), IsAlmostEqualTo(J{zero.a}));
  EXPECT_THAT(fdim(x.a, y), IsAlmostEqualTo(diffy));
  EXPECT_THAT(fdim(y, x.a), IsAlmostEqualTo(zero));
  EXPECT_TRUE(isnan(fdim(x, std::numeric_limits<J>::quiet_NaN())));
  EXPECT_TRUE(isnan(fdim(std::numeric_limits<J>::quiet_NaN(), x)));
  EXPECT_TRUE(isnan(fdim(x, std::numeric_limits<double>::quiet_NaN())));
  EXPECT_TRUE(isnan(fdim(std::numeric_limits<double>::quiet_NaN(), x)));

#ifndef CERES_NO_FENV_ACCESS
  EXPECT_EQ(std::fetestexcept(FE_ALL_EXCEPT & ~FE_INEXACT), 0);
#endif
}

TEST(Jet, CopySign) {
  {  // copysign(x, +1)
    J z = copysign(x, J{+1});
    EXPECT_THAT(z, IsAlmostEqualTo(x));
    EXPECT_TRUE(isfinite(z.v[0])) << z;
    EXPECT_TRUE(isfinite(z.v[1])) << z;
  }
  {  // copysign(x, -1)
    J z = copysign(x, J{-1});
    EXPECT_THAT(z, IsAlmostEqualTo(-x));
    EXPECT_TRUE(isfinite(z.v[0])) << z;
    EXPECT_TRUE(isfinite(z.v[1])) << z;
  }
  {  // copysign(-x, +1)

    J z = copysign(-x, J{+1});
    EXPECT_THAT(z, IsAlmostEqualTo(x));
    EXPECT_TRUE(isfinite(z.v[0])) << z;
    EXPECT_TRUE(isfinite(z.v[1])) << z;
  }
  {  // copysign(-x, -1)
    J z = copysign(-x, J{-1});
    EXPECT_THAT(z, IsAlmostEqualTo(-x));
    EXPECT_TRUE(isfinite(z.v[0])) << z;
    EXPECT_TRUE(isfinite(z.v[1])) << z;
  }
  {  // copysign(-0, +1)
    J z = copysign(MakeJet(-0, 1, 2), J{+1});
    EXPECT_THAT(z, IsAlmostEqualTo(MakeJet(+0, 1, 2)));
    EXPECT_FALSE(std::signbit(z.a)) << z;
    EXPECT_TRUE(isfinite(z.v[0])) << z;
    EXPECT_TRUE(isfinite(z.v[1])) << z;
  }
  {  // copysign(-0, -1)
    J z = copysign(MakeJet(-0, 1, 2), J{-1});
    EXPECT_THAT(z, IsAlmostEqualTo(MakeJet(-0, -1, -2)));
    EXPECT_TRUE(std::signbit(z.a)) << z;
    EXPECT_TRUE(isfinite(z.v[0])) << z;
    EXPECT_TRUE(isfinite(z.v[1])) << z;
  }
  {  // copysign(+0, -1)
    J z = copysign(MakeJet(+0, 1, 2), J{-1});
    EXPECT_THAT(z, IsAlmostEqualTo(MakeJet(-0, -1, -2)));
    EXPECT_TRUE(std::signbit(z.a)) << z;
    EXPECT_TRUE(isfinite(z.v[0])) << z;
    EXPECT_TRUE(isfinite(z.v[1])) << z;
  }
  {  // copysign(+0, +1)
    J z = copysign(MakeJet(+0, 1, 2), J{+1});
    EXPECT_THAT(z, IsAlmostEqualTo(MakeJet(+0, 1, 2)));
    EXPECT_FALSE(std::signbit(z.a)) << z;
    EXPECT_TRUE(isfinite(z.v[0])) << z;
    EXPECT_TRUE(isfinite(z.v[1])) << z;
  }
  {  // copysign(+0, +0)
    J z = copysign(MakeJet(+0, 1, 2), J{+0});
    EXPECT_FALSE(std::signbit(z.a)) << z;
    EXPECT_TRUE(isnan(z.v[0])) << z;
    EXPECT_TRUE(isnan(z.v[1])) << z;
  }
  {  // copysign(+0, -0)
    J z = copysign(MakeJet(+0, 1, 2), J{-0});
    EXPECT_FALSE(std::signbit(z.a)) << z;
    EXPECT_TRUE(isnan(z.v[0])) << z;
    EXPECT_TRUE(isnan(z.v[1])) << z;
  }
  {  // copysign(-0, +0)
    J z = copysign(MakeJet(-0, 1, 2), J{+0});
    EXPECT_FALSE(std::signbit(z.a)) << z;
    EXPECT_TRUE(isnan(z.v[0])) << z;
    EXPECT_TRUE(isnan(z.v[1])) << z;
  }
  {  // copysign(-0, -0)
    J z = copysign(MakeJet(-0, 1, 2), J{-0});
    EXPECT_FALSE(std::signbit(z.a)) << z;
    EXPECT_TRUE(isnan(z.v[0])) << z;
    EXPECT_TRUE(isnan(z.v[1])) << z;
  }
  {  // copysign(1, -nan)
    J z = copysign(MakeJet(1, 2, 3),
                   -J{std::numeric_limits<double>::quiet_NaN()});
    EXPECT_TRUE(std::signbit(z.a)) << z;
    EXPECT_TRUE(std::signbit(z.v[0])) << z;
    EXPECT_TRUE(std::signbit(z.v[1])) << z;
    EXPECT_FALSE(isnan(z.v[0])) << z;
    EXPECT_FALSE(isnan(z.v[1])) << z;
  }
  {  // copysign(1, +nan)
    J z = copysign(MakeJet(1, 2, 3),
                   +J{std::numeric_limits<double>::quiet_NaN()});
    EXPECT_FALSE(std::signbit(z.a)) << z;
    EXPECT_FALSE(std::signbit(z.v[0])) << z;
    EXPECT_FALSE(std::signbit(z.v[1])) << z;
    EXPECT_FALSE(isnan(z.v[0])) << z;
    EXPECT_FALSE(isnan(z.v[1])) << z;
  }
}

TEST(Jet, JetsInEigenMatrices) {
  J x = MakeJet(2.3, -2.7, 1e-3);
  J y = MakeJet(1.7, 0.5, 1e+2);
  J z = MakeJet(5.3, -4.7, 1e-3);
  J w = MakeJet(9.7, 1.5, 10.1);

  Eigen::Matrix<J, 2, 2> M;
  Eigen::Matrix<J, 2, 1> v, r1, r2;

  M << x, y, z, w;
  v << x, z;

  // M * v == (v^T * M^T)^T
  r1 = M * v;
  r2 = (v.transpose() * M.transpose()).transpose();

  EXPECT_THAT(r1(0), IsAlmostEqualTo(r2(0)));
  EXPECT_THAT(r1(1), IsAlmostEqualTo(r2(1)));
}

TEST(Jet, ScalarComparison) {
  Jet<double, 1> zero{0.0};
  zero.v << std::numeric_limits<double>::infinity();

  Jet<double, 1> one{1.0};
  one.v << std::numeric_limits<double>::quiet_NaN();

  Jet<double, 1> two{2.0};
  two.v << std::numeric_limits<double>::min() / 2;

  Jet<double, 1> three{3.0};

  auto inf = std::numeric_limits<Jet<double, 1>>::infinity();
  auto nan = std::numeric_limits<Jet<double, 1>>::quiet_NaN();
  inf.v << 1.2;
  nan.v << 3.4;

  std::feclearexcept(FE_ALL_EXCEPT);

  EXPECT_FALSE(islessgreater(zero, zero));
  EXPECT_FALSE(islessgreater(zero, zero.a));
  EXPECT_FALSE(islessgreater(zero.a, zero));

  EXPECT_TRUE(isgreaterequal(three, three));
  EXPECT_TRUE(isgreaterequal(three, three.a));
  EXPECT_TRUE(isgreaterequal(three.a, three));

  EXPECT_TRUE(isgreater(three, two));
  EXPECT_TRUE(isgreater(three, two.a));
  EXPECT_TRUE(isgreater(three.a, two));

  EXPECT_TRUE(islessequal(one, one));
  EXPECT_TRUE(islessequal(one, one.a));
  EXPECT_TRUE(islessequal(one.a, one));

  EXPECT_TRUE(isless(one, two));
  EXPECT_TRUE(isless(one, two.a));
  EXPECT_TRUE(isless(one.a, two));

  EXPECT_FALSE(isunordered(inf, one));
  EXPECT_FALSE(isunordered(inf, one.a));
  EXPECT_FALSE(isunordered(inf.a, one));

  EXPECT_TRUE(isunordered(nan, two));
  EXPECT_TRUE(isunordered(nan, two.a));
  EXPECT_TRUE(isunordered(nan.a, two));

  EXPECT_TRUE(isunordered(inf, nan));
  EXPECT_TRUE(isunordered(inf, nan.a));
  EXPECT_TRUE(isunordered(inf.a, nan.a));

  EXPECT_EQ(std::fetestexcept(FE_ALL_EXCEPT & ~FE_INEXACT), 0);
}

TEST(Jet, Nested2XScalarComparison) {
  Jet<J0d, 1> zero{J0d{0.0}};
  zero.v << std::numeric_limits<J0d>::infinity();

  Jet<J0d, 1> one{J0d{1.0}};
  one.v << std::numeric_limits<J0d>::quiet_NaN();

  Jet<J0d, 1> two{J0d{2.0}};
  two.v << std::numeric_limits<J0d>::min() / J0d{2};

  Jet<J0d, 1> three{J0d{3.0}};

  auto inf = std::numeric_limits<Jet<J0d, 1>>::infinity();
  auto nan = std::numeric_limits<Jet<J0d, 1>>::quiet_NaN();
  inf.v << J0d{1.2};
  nan.v << J0d{3.4};

  std::feclearexcept(FE_ALL_EXCEPT);

  EXPECT_FALSE(islessgreater(zero, zero));
  EXPECT_FALSE(islessgreater(zero, zero.a));
  EXPECT_FALSE(islessgreater(zero.a, zero));
  EXPECT_FALSE(islessgreater(zero, zero.a.a));
  EXPECT_FALSE(islessgreater(zero.a.a, zero));

  EXPECT_TRUE(isgreaterequal(three, three));
  EXPECT_TRUE(isgreaterequal(three, three.a));
  EXPECT_TRUE(isgreaterequal(three.a, three));
  EXPECT_TRUE(isgreaterequal(three, three.a.a));
  EXPECT_TRUE(isgreaterequal(three.a.a, three));

  EXPECT_TRUE(isgreater(three, two));
  EXPECT_TRUE(isgreater(three, two.a));
  EXPECT_TRUE(isgreater(three.a, two));
  EXPECT_TRUE(isgreater(three, two.a.a));
  EXPECT_TRUE(isgreater(three.a.a, two));

  EXPECT_TRUE(islessequal(one, one));
  EXPECT_TRUE(islessequal(one, one.a));
  EXPECT_TRUE(islessequal(one.a, one));
  EXPECT_TRUE(islessequal(one, one.a.a));
  EXPECT_TRUE(islessequal(one.a.a, one));

  EXPECT_TRUE(isless(one, two));
  EXPECT_TRUE(isless(one, two.a));
  EXPECT_TRUE(isless(one.a, two));
  EXPECT_TRUE(isless(one, two.a.a));
  EXPECT_TRUE(isless(one.a.a, two));

  EXPECT_FALSE(isunordered(inf, one));
  EXPECT_FALSE(isunordered(inf, one.a));
  EXPECT_FALSE(isunordered(inf.a, one));
  EXPECT_FALSE(isunordered(inf, one.a.a));
  EXPECT_FALSE(isunordered(inf.a.a, one));

  EXPECT_TRUE(isunordered(nan, two));
  EXPECT_TRUE(isunordered(nan, two.a));
  EXPECT_TRUE(isunordered(nan.a, two));
  EXPECT_TRUE(isunordered(nan, two.a.a));
  EXPECT_TRUE(isunordered(nan.a.a, two));

  EXPECT_TRUE(isunordered(inf, nan));
  EXPECT_TRUE(isunordered(inf, nan.a));
  EXPECT_TRUE(isunordered(inf.a, nan));
  EXPECT_TRUE(isunordered(inf, nan.a.a));
  EXPECT_TRUE(isunordered(inf.a.a, nan));

  EXPECT_EQ(std::fetestexcept(FE_ALL_EXCEPT & ~FE_INEXACT), 0);
}

TEST(JetTraitsTest, ClassificationNaN) {
  Jet<double, 1> a(std::numeric_limits<double>::quiet_NaN());
  a.v << std::numeric_limits<double>::infinity();
  EXPECT_EQ(fpclassify(a), FP_NAN);
  EXPECT_FALSE(isfinite(a));
  EXPECT_FALSE(isinf(a));
  EXPECT_FALSE(isnormal(a));
  EXPECT_FALSE(signbit(a));
  EXPECT_TRUE(isnan(a));
}

TEST(JetTraitsTest, ClassificationInf) {
  Jet<double, 1> a(-std::numeric_limits<double>::infinity());
  a.v << std::numeric_limits<double>::quiet_NaN();
  EXPECT_EQ(fpclassify(a), FP_INFINITE);
  EXPECT_FALSE(isfinite(a));
  EXPECT_FALSE(isnan(a));
  EXPECT_FALSE(isnormal(a));
  EXPECT_TRUE(signbit(a));
  EXPECT_TRUE(isinf(a));
}

TEST(JetTraitsTest, ClassificationFinite) {
  Jet<double, 1> a(-5.5);
  a.v << std::numeric_limits<double>::quiet_NaN();
  EXPECT_EQ(fpclassify(a), FP_NORMAL);
  EXPECT_FALSE(isinf(a));
  EXPECT_FALSE(isnan(a));
  EXPECT_TRUE(signbit(a));
  EXPECT_TRUE(isfinite(a));
  EXPECT_TRUE(isnormal(a));
}

TEST(JetTraitsTest, ClassificationScalar) {
  EXPECT_EQ(fpclassify(J0d{+0.0}), FP_ZERO);
  EXPECT_EQ(fpclassify(J0d{-0.0}), FP_ZERO);
  EXPECT_EQ(fpclassify(J0d{1.234}), FP_NORMAL);
  EXPECT_EQ(fpclassify(J0d{std::numeric_limits<double>::min() / 2}),
            FP_SUBNORMAL);
  EXPECT_EQ(fpclassify(J0d{std::numeric_limits<double>::quiet_NaN()}), FP_NAN);
}

TEST(JetTraitsTest, Nested2XClassificationScalar) {
  EXPECT_EQ(fpclassify(J0<J0d>{J0d{+0.0}}), FP_ZERO);
  EXPECT_EQ(fpclassify(J0<J0d>{J0d{-0.0}}), FP_ZERO);
  EXPECT_EQ(fpclassify(J0<J0d>{J0d{1.234}}), FP_NORMAL);
  EXPECT_EQ(fpclassify(J0<J0d>{J0d{std::numeric_limits<double>::min() / 2}}),
            FP_SUBNORMAL);
  EXPECT_EQ(fpclassify(J0<J0d>{J0d{std::numeric_limits<double>::quiet_NaN()}}),
            FP_NAN);
}

// The following test ensures that Jets have all the appropriate Eigen
// related traits so that they can be used as part of matrix
// decompositions.
TEST(Jet, FullRankEigenLLTSolve) {
  Eigen::Matrix<J, 3, 3> A;
  Eigen::Matrix<J, 3, 1> b, x;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      A(i, j) = MakeJet(0.0, i, j * j);
    }
    b(i) = MakeJet(i, i, i);
    x(i) = MakeJet(0.0, 0.0, 0.0);
    A(i, i) = MakeJet(1.0, i, i * i);
  }
  x = A.llt().solve(b);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(x(i).a, b(i).a);
  }
}

TEST(Jet, FullRankEigenLDLTSolve) {
  Eigen::Matrix<J, 3, 3> A;
  Eigen::Matrix<J, 3, 1> b, x;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      A(i, j) = MakeJet(0.0, i, j * j);
    }
    b(i) = MakeJet(i, i, i);
    x(i) = MakeJet(0.0, 0.0, 0.0);
    A(i, i) = MakeJet(1.0, i, i * i);
  }
  x = A.ldlt().solve(b);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(x(i).a, b(i).a);
  }
}

TEST(Jet, FullRankEigenLUSolve) {
  Eigen::Matrix<J, 3, 3> A;
  Eigen::Matrix<J, 3, 1> b, x;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      A(i, j) = MakeJet(0.0, i, j * j);
    }
    b(i) = MakeJet(i, i, i);
    x(i) = MakeJet(0.0, 0.0, 0.0);
    A(i, i) = MakeJet(1.0, i, i * i);
  }

  x = A.lu().solve(b);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(x(i).a, b(i).a);
  }
}

// ScalarBinaryOpTraits is only supported on Eigen versions >= 3.3
TEST(JetTraitsTest, MatrixScalarUnaryOps) {
  const J x = MakeJet(2.3, -2.7, 1e-3);
  const J y = MakeJet(1.7, 0.5, 1e+2);
  Eigen::Matrix<J, 2, 1> a;
  a << x, y;

  const J sum = a.sum();
  const J sum2 = a(0) + a(1);
  EXPECT_THAT(sum, IsAlmostEqualTo(sum2));
}

TEST(JetTraitsTest, MatrixScalarBinaryOps) {
  const J x = MakeJet(2.3, -2.7, 1e-3);
  const J y = MakeJet(1.7, 0.5, 1e+2);
  const J z = MakeJet(5.3, -4.7, 1e-3);
  const J w = MakeJet(9.7, 1.5, 10.1);

  Eigen::Matrix<J, 2, 2> M;
  Eigen::Vector2d v;

  M << x, y, z, w;
  v << 0.6, -2.1;

  // M * v == M * v.cast<J>().
  const Eigen::Matrix<J, 2, 1> r1 = M * v;
  const Eigen::Matrix<J, 2, 1> r2 = M * v.cast<J>();

  EXPECT_THAT(r1(0), IsAlmostEqualTo(r2(0)));
  EXPECT_THAT(r1(1), IsAlmostEqualTo(r2(1)));

  // M * a == M * T(a).
  const double a = 3.1;
  const Eigen::Matrix<J, 2, 2> r3 = M * a;
  const Eigen::Matrix<J, 2, 2> r4 = M * J(a);

  EXPECT_THAT(r3(0, 0), IsAlmostEqualTo(r4(0, 0)));
  EXPECT_THAT(r3(0, 1), IsAlmostEqualTo(r4(0, 1)));
  EXPECT_THAT(r3(1, 0), IsAlmostEqualTo(r4(1, 0)));
  EXPECT_THAT(r3(1, 1), IsAlmostEqualTo(r4(1, 1)));
}

TEST(JetTraitsTest, ArrayScalarUnaryOps) {
  const J x = MakeJet(2.3, -2.7, 1e-3);
  const J y = MakeJet(1.7, 0.5, 1e+2);
  Eigen::Array<J, 2, 1> a;
  a << x, y;

  const J sum = a.sum();
  const J sum2 = a(0) + a(1);
  EXPECT_THAT(sum, sum2);
}

TEST(JetTraitsTest, ArrayScalarBinaryOps) {
  const J x = MakeJet(2.3, -2.7, 1e-3);
  const J y = MakeJet(1.7, 0.5, 1e+2);

  Eigen::Array<J, 2, 1> a;
  Eigen::Array2d b;

  a << x, y;
  b << 0.6, -2.1;

  // a * b == a * b.cast<T>()
  const Eigen::Array<J, 2, 1> r1 = a * b;
  const Eigen::Array<J, 2, 1> r2 = a * b.cast<J>();

  EXPECT_THAT(r1(0), r2(0));
  EXPECT_THAT(r1(1), r2(1));

  // a * c == a * T(c).
  const double c = 3.1;
  const Eigen::Array<J, 2, 1> r3 = a * c;
  const Eigen::Array<J, 2, 1> r4 = a * J(c);

  EXPECT_THAT(r3(0), r3(0));
  EXPECT_THAT(r4(1), r4(1));
}

TEST(Jet, Nested3X) {
  using JJ = Jet<J, 2>;
  using JJJ = Jet<JJ, 2>;

  JJJ x;
  x.a = JJ(J(1, 0), 0);
  x.v[0] = JJ(J(1));

  JJJ y = x * x * x;

  ExpectClose(y.a.a.a, 1, kTolerance);
  ExpectClose(y.v[0].a.a, 3., kTolerance);
  ExpectClose(y.v[0].v[0].a, 6., kTolerance);
  ExpectClose(y.v[0].v[0].v[0], 6., kTolerance);

  JJJ e = exp(x);

  ExpectClose(e.a.a.a, kE, kTolerance);
  ExpectClose(e.v[0].a.a, kE, kTolerance);
  ExpectClose(e.v[0].v[0].a, kE, kTolerance);
  ExpectClose(e.v[0].v[0].v[0], kE, kTolerance);
}

#if GTEST_HAS_TYPED_TEST

using Types = testing::Types<std::int16_t,
                             std::uint16_t,
                             std::int32_t,
                             std::uint32_t,
                             std::int64_t,
                             std::uint64_t,
                             float,
                             double,
                             long double>;

template <typename T>
class JetTest : public testing::Test {};

TYPED_TEST_SUITE(JetTest, Types);

TYPED_TEST(JetTest, Comparison) {
  using Scalar = TypeParam;

  EXPECT_EQ(J0<Scalar>{0}, J0<Scalar>{0});
  EXPECT_GE(J0<Scalar>{3}, J0<Scalar>{3});
  EXPECT_GT(J0<Scalar>{3}, J0<Scalar>{2});
  EXPECT_LE(J0<Scalar>{1}, J0<Scalar>{1});
  EXPECT_LT(J0<Scalar>{1}, J0<Scalar>{2});
  EXPECT_NE(J0<Scalar>{1}, J0<Scalar>{2});
}

TYPED_TEST(JetTest, ScalarComparison) {
  using Scalar = TypeParam;

  EXPECT_EQ(J0d{0.0}, Scalar{0});
  EXPECT_GE(J0d{3.0}, Scalar{3});
  EXPECT_GT(J0d{3.0}, Scalar{2});
  EXPECT_LE(J0d{1.0}, Scalar{1});
  EXPECT_LT(J0d{1.0}, Scalar{2});
  EXPECT_NE(J0d{1.0}, Scalar{2});

  EXPECT_EQ(Scalar{0}, J0d{0.0});
  EXPECT_GE(Scalar{1}, J0d{1.0});
  EXPECT_GT(Scalar{2}, J0d{1.0});
  EXPECT_LE(Scalar{3}, J0d{3.0});
  EXPECT_LT(Scalar{2}, J0d{3.0});
  EXPECT_NE(Scalar{2}, J0d{1.0});
}

TYPED_TEST(JetTest, Nested2XComparison) {
  using Scalar = TypeParam;

  EXPECT_EQ(J0<J0d>{J0d{0.0}}, Scalar{0});
  EXPECT_GE(J0<J0d>{J0d{3.0}}, Scalar{3});
  EXPECT_GT(J0<J0d>{J0d{3.0}}, Scalar{2});
  EXPECT_LE(J0<J0d>{J0d{1.0}}, Scalar{1});
  EXPECT_LT(J0<J0d>{J0d{1.0}}, Scalar{2});
  EXPECT_NE(J0<J0d>{J0d{1.0}}, Scalar{2});

  EXPECT_EQ(Scalar{0}, J0<J0d>{J0d{0.0}});
  EXPECT_GE(Scalar{1}, J0<J0d>{J0d{1.0}});
  EXPECT_GT(Scalar{2}, J0<J0d>{J0d{1.0}});
  EXPECT_LE(Scalar{3}, J0<J0d>{J0d{3.0}});
  EXPECT_LT(Scalar{2}, J0<J0d>{J0d{3.0}});
  EXPECT_NE(Scalar{2}, J0<J0d>{J0d{1.0}});
}

TYPED_TEST(JetTest, Nested3XComparison) {
  using Scalar = TypeParam;

  EXPECT_EQ(J0<J0<J0d>>{J0<J0d>{J0d{0.0}}}, Scalar{0});
  EXPECT_GE(J0<J0<J0d>>{J0<J0d>{J0d{3.0}}}, Scalar{3});
  EXPECT_GT(J0<J0<J0d>>{J0<J0d>{J0d{3.0}}}, Scalar{2});
  EXPECT_LE(J0<J0<J0d>>{J0<J0d>{J0d{1.0}}}, Scalar{1});
  EXPECT_LT(J0<J0<J0d>>{J0<J0d>{J0d{1.0}}}, Scalar{2});
  EXPECT_NE(J0<J0<J0d>>{J0<J0d>{J0d{1.0}}}, Scalar{2});

  EXPECT_EQ(Scalar{0}, J0<J0<J0d>>{J0<J0d>{J0d{0.0}}});
  EXPECT_GE(Scalar{1}, J0<J0<J0d>>{J0<J0d>{J0d{1.0}}});
  EXPECT_GT(Scalar{2}, J0<J0<J0d>>{J0<J0d>{J0d{1.0}}});
  EXPECT_LE(Scalar{3}, J0<J0<J0d>>{J0<J0d>{J0d{3.0}}});
  EXPECT_LT(Scalar{2}, J0<J0<J0d>>{J0<J0d>{J0d{3.0}}});
  EXPECT_NE(Scalar{2}, J0<J0<J0d>>{J0<J0d>{J0d{1.0}}});
}

#endif  // GTEST_HAS_TYPED_TEST

}  // namespace ceres::internal

#ifdef _MSC_VER
#pragma float_control(pop)
#endif
