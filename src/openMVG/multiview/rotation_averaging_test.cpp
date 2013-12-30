// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "CppUnitLite/TestHarness.h"
#include "openMVG/multiview/rotation_averaging.hpp"
#include "openMVG/multiview/essential.hpp"
#include "testing/testing.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>

using namespace openMVG;
using namespace rotation_averaging;

TEST ( rotation_averaging, ClosestSVDRotationMatrix )
{
  Mat3 rotx = RotationAroundX(0.3);

  Mat3 Approximative_rotx = rotation_averaging::ClosestSVDRotationMatrix(rotx);

  // Check that SVD have rebuilt the matrix correctly
  EXPECT_MATRIX_NEAR( rotx, Approximative_rotx, 1e-8);
  // Check the Frobenius distance between the approximated rot matrix and the GT
  EXPECT_NEAR( 0.0, FrobeniusDistance( rotx, Approximative_rotx), 1e-8);
  // Check that the Matrix is a rotation matrix (determinant == 1)
  EXPECT_NEAR( 1.0, Approximative_rotx.determinant(), 1e-8);
}


TEST ( rotation_averaging, ClosestSVDRotationMatrixNoisy )
{
  Mat3 rotx = RotationAroundX(0.3);

  //-- Set a little of noise in the rotMatrix :
  rotx(2,2) -= 0.02;
  Mat3 Approximative_rotx = rotation_averaging::ClosestSVDRotationMatrix(rotx);

  // Check the Frobenius distance between the approximated rot matrix and the GT
  CHECK( FrobeniusDistance( rotx, Approximative_rotx) < 0.02);
  // Check that the Matrix is a rotation matrix (determinant == 1)
  EXPECT_NEAR( 1.0, Approximative_rotx.determinant(), 1e-8);
}

// Rotation averaging in a triplet:
// 0_______2
//  \     /
//   \   /
//    \ /
//     1
TEST ( rotation_averaging, RotationLeastSquare_3_Camera)
{
  using namespace std;

  //--
  // Setup 3 camera that have a relative orientation of 120°
  // Set Z axis as UP Vector for the rotation
  // They are in the same plane and looking in O={0,0,0}
  //--
  const int nCamera = 3;
  Mat3 R01 = RotationAroundZ(2.*M_PI/3.0); //120°
  Mat3 R12 = RotationAroundZ(2.*M_PI/3.0); //120°
  Mat3 R20 = RotationAroundZ(2.*M_PI/3.0); //120°
  Mat3 Id = Mat3::Identity();

  std::vector<std::pair<std::pair<size_t, size_t>, Mat3> > vec_relativeRotEstimate;
  vec_relativeRotEstimate.push_back( make_pair(make_pair(0,1), RotationAroundZ( 2.*M_PI/3.0)));
  vec_relativeRotEstimate.push_back( make_pair(make_pair(1,2), RotationAroundZ( 2.*M_PI/3.0)));
  vec_relativeRotEstimate.push_back( make_pair(make_pair(2,0), RotationAroundZ( 2.*M_PI/3.0)));

  //- Solve the global rotation estimation problem :
  std::vector<Mat3> vec_globalR;
  L2RotationAveraging(3, vec_relativeRotEstimate, vec_globalR);
  EXPECT_EQ(3, vec_globalR.size());
  // Check that the loop is closing
  EXPECT_MATRIX_NEAR(Mat3::Identity(), (vec_globalR[0]*vec_globalR[1]*vec_globalR[2]), 1e-8);

  //--
  // Check that the found relative rotation matrix give the expected rotation.
  //  -> the started relative rotation (used in the action matrix).
  //// /!\ Translation are not checked they are 0 by default.
  //--
  Mat3 R;
  Vec3 t, t0 = Vec3::Zero(), t1 = Vec3::Zero();
  RelativeCameraMotion(vec_globalR[0], t0, vec_globalR[1], t1, &R, &t);
  EXPECT_NEAR( 0, FrobeniusDistance( R01, R), 1e-2);
  std::cout << "\n" << R << "\n\n\n";

  RelativeCameraMotion(vec_globalR[1], t0, vec_globalR[2], t1, &R, &t);
  EXPECT_NEAR( 0, FrobeniusDistance( R12, R), 1e-2);
  std::cout << "\n" << R << "\n\n\n";

  RelativeCameraMotion(vec_globalR[2], t0, vec_globalR[0], t1, &R, &t);
  EXPECT_NEAR( 0, FrobeniusDistance( R20, R), 1e-2);
  std::cout << "\n" << R << "\n\n\n";
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
