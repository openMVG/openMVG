
// Copyright (c) 2007, 2008 libmv authors.
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


#ifndef TESTING_TESTING_H_
#define TESTING_TESTING_H_

#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/numeric/numeric.h"
#include "third_party/CppUnitLite/TestHarness.h"

#define EXPECT_MATRIX_NEAR(a, b, tolerance) \
do { \
  bool dims_match = (a.rows() == b.rows()) && (a.cols() == b.cols()); \
  CHECK_EQUAL(a.rows(),b.rows()); \
  CHECK_EQUAL(a.cols(),b.cols()); \
  if (dims_match) { \
    for (int r = 0; r < a.rows(); ++r) { \
      for (int c = 0; c < a.cols(); ++c) { \
        DOUBLES_EQUAL(a(r, c), b(r, c), tolerance); \
      } \
    } \
  } \
} while(false);

#define EXPECT_MATRIX_NEAR_ZERO(a, tolerance) \
do { \
  for (int r = 0; r < a.rows(); ++r) { \
    for (int c = 0; c < a.cols(); ++c) { \
      DOUBLES_EQUAL(0.0, a(r, c), tolerance) \
    } \
  } \
} while(false);

#define EXPECT_MATRIX_EQ(a, b) \
do { \
  bool dims_match = (a.rows() == b.rows()) && (a.cols() == b.cols()); \
  CHECK_EQUAL(a.rows(),b.rows()); \
  CHECK_EQUAL(a.cols(),b.cols()); \
  if (dims_match) { \
    for (int r = 0; r < a.rows(); ++r) { \
      for (int c = 0; c < a.cols(); ++c) { \
        CHECK_EQUAL(a(r, c), b(r, c)) \
      } \
    } \
  } \
} while(false);

#define EXPECT_NEAR(expected,actual,threshold) \
  DOUBLES_EQUAL(expected,actual,threshold);

#define EXPECT_TRUE(condition) \
  CHECK(condition);

#define EXPECT_FALSE(condition) \
  CHECK(!condition);

#define EXPECT_EQ(a, b) CHECK_EQUAL(a,b);

// Check that sin(angle(a, b)) < tolerance.
#define EXPECT_MATRIX_PROP(a, b, tolerance) \
do { \
  bool dims_match = (a.rows() == b.rows()) && (a.cols() == b.cols()); \
  CHECK_EQUAL(a.rows(), b.rows()); \
  CHECK_EQUAL(a.cols(), b.cols()); \
  if (dims_match) { \
    double c = CosinusBetweenMatrices(a, b); \
    if (c * c < 1) { \
      double s = sqrt(1 - c * c); \
      DOUBLES_EQUAL(0.0, s, tolerance); \
    } \
  } \
} while(false);

#endif  // TESTING_TESTING_H_
