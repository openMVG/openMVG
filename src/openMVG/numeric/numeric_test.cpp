
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

#include "openMVG/numeric/numeric.h"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include <iostream>

using namespace openMVG;
using namespace std;

//-- Assert that stream interface is available
TEST ( TinyMatrix, print )
{
  Mat3 testMatrix = Mat3::Identity();
  std::cout << testMatrix;
}

TEST ( TinyMatrix, checkIdentity )
{
  Mat3 testMatrix, expected;

  // build expected matrix and the test matrix
  expected.fill(0);
  expected(0,0) = expected(1,1) = expected(2,2) = 1.0;

  testMatrix.setIdentity();
  std::cout << std::endl << testMatrix;
  //-- Compare expected to the testMatrix.
  EXPECT_MATRIX_NEAR( expected, testMatrix, 1e-8);
}

TEST ( TinyMatrix, product )
{
  Mat3 a, b, expected;

  // build expected matrix and the test matrix
  a(0,0) = 1.0; a(0,1) = 2.0; a(0,2) = 3.0;
  a(1,0) = 4.0; a(1,1) = 5.0; a(1,2) = 6.0;
  a(2,0) = 7.0; a(2,1) = 8.0; a(2,2) = 9.0;

  b(0,0) = 10.0; b(0,1) = 11.0; b(0,2) = 12.0;
  b(1,0) = 13.0; b(1,1) = 14.0; b(1,2) = 15.0;
  b(2,0) = 16.0; b(2,1) = 17.0; b(2,2) = 18.0;

  Mat3 resAxB = a*b;
  Mat3 expected_resAxB;
  {
    Mat3 & t = expected_resAxB;
    t(0,0) = 84.0;  t(0,1) = 90.0;    t(0,2) = 96.0;
    t(1,0) = 201.0; t(1,1) = 216.0;  t(1,2) = 231.0;
    t(2,0) = 318.0; t(2,1) = 342.0;  t(2,2) = 366.0;
  }

  Mat3 resBxA = b*a;
  Mat3 expected_resBxA;
  {
    Mat3 & t = expected_resBxA;
    t(0,0) = 138; t(0,1) = 171;  t(0,2) = 204;
    t(1,0) = 174; t(1,1) = 216;  t(1,2) = 258;
    t(2,0) = 210; t(2,1) = 261;  t(2,2) = 312;
  }

  //-- Tests
  EXPECT_MATRIX_NEAR( expected_resAxB, resAxB, 1e-8);
  EXPECT_MATRIX_NEAR( expected_resBxA, resBxA, 1e-8);
}

TEST(TinyMatrix, LookAt) {
  // Simple orthogonality check.
  Vec3 e; e[0]= 1; e[1] = 2; e[2] = 3;
  Mat3 R = LookAt(e);
  Mat3 I = Mat3::Identity();
  Mat3 RRT = R*R.transpose();
  Mat3 RTR = R.transpose()*R;

  EXPECT_MATRIX_NEAR(I, RRT, 1e-15);
  EXPECT_MATRIX_NEAR(I, RTR, 1e-15);
}

TEST(Numeric, MeanAndVarianceAlongRows) {
  int n = 4;
  Mat points(2,n);
  points << 0, 0, 1, 1,
    0, 2, 1, 3;

  Vec mean, variance;
  MeanAndVarianceAlongRows(points, &mean, &variance);

  EXPECT_NEAR(0.5, mean(0), 1e-8);
  EXPECT_NEAR(1.5, mean(1), 1e-8);
  EXPECT_NEAR(0.25, variance(0), 1e-8);
  EXPECT_NEAR(1.25, variance(1), 1e-8);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
