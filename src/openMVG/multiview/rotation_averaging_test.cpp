// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "CppUnitLite/TestHarness.h"
#include "openMVG/multiview/rotation_averaging.hpp"
#include "openMVG/multiview/essential.hpp"
#include "testing/testing.h"

#include "openMVG/multiview/test_data_sets.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>

using namespace openMVG;
using namespace openMVG::rotation_averaging;
using namespace openMVG::rotation_averaging::l1;
using namespace openMVG::rotation_averaging::l2;

TEST ( rotation_averaging, ClosestSVDRotationMatrix )
{
  Mat3 rotx = RotationAroundX(0.3);

  Mat3 Approximative_rotx = rotation_averaging::l2::ClosestSVDRotationMatrix(rotx);

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
  Mat3 Approximative_rotx = rotation_averaging::l2::ClosestSVDRotationMatrix(rotx);

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
  Mat3 R01 = RotationAroundZ(2.*M_PI/3.0); //120°
  Mat3 R12 = RotationAroundZ(2.*M_PI/3.0); //120°
  Mat3 R20 = RotationAroundZ(2.*M_PI/3.0); //120°
  Mat3 Id = Mat3::Identity();

  std::vector<RelativeRotation > vec_relativeRotEstimate;
  vec_relativeRotEstimate.push_back( RelativeRotation(0,1, R01));
  vec_relativeRotEstimate.push_back( RelativeRotation(1,2, R12));
  vec_relativeRotEstimate.push_back( RelativeRotation(2,0, R20));

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

  RelativeCameraMotion(vec_globalR[1], t0, vec_globalR[2], t1, &R, &t);
  EXPECT_NEAR( 0, FrobeniusDistance( R12, R), 1e-2);

  RelativeCameraMotion(vec_globalR[2], t0, vec_globalR[0], t1, &R, &t);
  EXPECT_NEAR( 0, FrobeniusDistance( R20, R), 1e-2);
}

TEST ( rotation_averaging, RefineRotationsAvgL1IRLS_SimpleTriplet)
{
  using namespace std;

  //--
  // Setup 3 camera that have a relative orientation of 120°
  // Set Z axis as UP Vector for the rotation
  // They are in the same plane and looking in O={0,0,0}
  //--
  Mat3 R01 = RotationAroundZ(2.*M_PI/3.0); //120°
  Mat3 R12 = RotationAroundZ(2.*M_PI/3.0); //120°
  Mat3 R20 = RotationAroundZ(2.*M_PI/3.0); //120°
  Mat3 Id = Mat3::Identity();

  // Setup the relative motions (relative rotations)
  RelativeRotations vec_relativeRotEstimate;
  vec_relativeRotEstimate.push_back(RelativeRotation(0, 1, R01, 1));
  vec_relativeRotEstimate.push_back(RelativeRotation(1, 2, R12, 1));
  vec_relativeRotEstimate.push_back(RelativeRotation(2, 0, R20, 1));

  //- Solve the global rotation estimation problem :
  Matrix3x3Arr vec_globalR(3);
  size_t nMainViewID = 0;
  EXPECT_TRUE(GlobalRotationsRobust(vec_relativeRotEstimate, vec_globalR, nMainViewID, 0.0f, NULL));

  // Check that the loop is closing
  EXPECT_MATRIX_NEAR(Mat3::Identity(), (vec_globalR[0]*vec_globalR[1]*vec_globalR[2]), 1e-4);

  //--
  // Check that the found relative rotation matrix give the expected rotation.
  //  -> the started relative rotation (used in the action matrix).
  //// /!\ Translation are not checked they are 0 by default.
  //--
  Mat3 R;
  Vec3 t, t0 = Vec3::Zero(), t1 = Vec3::Zero();
  RelativeCameraMotion(vec_globalR[0], t0, vec_globalR[1], t1, &R, &t);
  EXPECT_NEAR( 0, FrobeniusDistance( R01, R), 1e-8);

  RelativeCameraMotion(vec_globalR[1], t0, vec_globalR[2], t1, &R, &t);
  EXPECT_NEAR( 0, FrobeniusDistance( R12, R), 1e-8);

  RelativeCameraMotion(vec_globalR[2], t0, vec_globalR[0], t1, &R, &t);
  EXPECT_NEAR( 0, FrobeniusDistance( R20, R), 1e-8);
}

// Test over a loop of cameras
TEST ( rotation_averaging, RefineRotationsAvgL1IRLS_CompleteGraph)
{
  //-- Setup a circular camera rig
  const int iNviews = 5;
  NViewDataSet d = NRealisticCamerasRing(iNviews, 5,
    nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  //Link each camera to the two next ones
  RelativeRotations vec_relativeRotEstimate;
  for (size_t i = 0; i < iNviews; ++i)
  {
    size_t index0 = i;
    size_t index1 = (i+1)%iNviews;
    size_t index2 = (i+2)%iNviews;

    //-- compute true relative rotations of the triplet
    // (index0->index1), (index1,index2), (index0->index2)
    Mat3 Rrel;
    Vec3 trel;
    RelativeCameraMotion(d._R[index0], d._t[index0], d._R[index1], d._t[index1], &Rrel, &trel);
    vec_relativeRotEstimate.push_back(RelativeRotation(index0, index1, Rrel, 1));

    RelativeCameraMotion(d._R[index1], d._t[index1], d._R[index2], d._t[index2], &Rrel, &trel);
    vec_relativeRotEstimate.push_back(RelativeRotation(index1, index2, Rrel, 1));

    RelativeCameraMotion(d._R[index0], d._t[index0], d._R[index2], d._t[index2], &Rrel, &trel);
    vec_relativeRotEstimate.push_back(RelativeRotation(index0, index2, Rrel, 1));
  }

  //- Solve the global rotation estimation problem :
  Matrix3x3Arr vec_globalR(iNviews);
  size_t nMainViewID = 0;
  bool bTest = GlobalRotationsRobust(vec_relativeRotEstimate, vec_globalR, nMainViewID, 0.0f, NULL);
  EXPECT_TRUE(bTest);

  // Check that the loop is closing
  Mat3 rotCumul = Mat3::Identity();
  for (size_t i = 0; i < iNviews; ++i)
  {
    rotCumul*= vec_globalR[i];
  }
  // Fix the sign of the rotations (put the global rotation in the same rotation axis as GT)
  if ( SIGN(vec_globalR[0](0,0)) != SIGN( d._R[0](0,0) ))
  {
    Mat3 Rrel;
    Vec3 trel;
    RelativeCameraMotion(vec_globalR[0], Vec3::Zero(), d._R[0], Vec3::Zero(), &Rrel, &trel);
    for (size_t i = 0; i < iNviews; ++i)
      vec_globalR[i] *= Rrel;
  }

  // Check that each global rotations is near the true ones
  for (size_t i = 0; i < iNviews; ++i)
  {
    EXPECT_NEAR(0.0, FrobeniusDistance(d._R[i], vec_globalR[i]), 1e-8);
  }
  EXPECT_MATRIX_NEAR(Mat3::Identity(), rotCumul, 1e-4);
}


// Test over a loop of camera that have 2 relative outliers rotations
TEST ( rotation_averaging, RefineRotationsAvgL1IRLS_CompleteGraph_outliers)
{
  //-- Setup a circular camera rig
  const int iNviews = 5;
  NViewDataSet d = NRealisticCamerasRing(iNviews, 5,
    nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  //Link each camera to the two next ones
  RelativeRotations vec_relativeRotEstimate;
  std::vector<std::pair<size_t,size_t> > vec_unique;
  for (size_t i = 0; i < iNviews; ++i)
  {
    size_t index0 = i;
    size_t index1 = (i+1)%iNviews;
    size_t index2 = (i+2)%iNviews;

    //-- compute true relative rotations of the triplet
    // (index0->index1), (index1,index2), (index0->index2)

    Mat3 Rrel;
    Vec3 trel;
    // Add the relative motion if do not exist yet
    if ( std::find(vec_unique.begin(), vec_unique.end(), std::make_pair(index0, index1)) == vec_unique.end()
      && std::find(vec_unique.begin(), vec_unique.end(), std::make_pair(index1, index0)) == vec_unique.end())
    {
      RelativeCameraMotion(d._R[index0], d._t[index0], d._R[index1], d._t[index1], &Rrel, &trel);
      vec_relativeRotEstimate.push_back(RelativeRotation(index0, index1, Rrel, 1));
      vec_unique.push_back(make_pair(index0, index1));
    }

    if ( std::find(vec_unique.begin(), vec_unique.end(), std::make_pair(index1, index2)) == vec_unique.end()
      && std::find(vec_unique.begin(), vec_unique.end(), std::make_pair(index2, index1)) == vec_unique.end())
    {
      RelativeCameraMotion(d._R[index1], d._t[index1], d._R[index2], d._t[index2], &Rrel, &trel);
      vec_relativeRotEstimate.push_back(RelativeRotation(index1, index2, Rrel, 1));
      vec_unique.push_back(make_pair(index1, index2));
    }

    if ( std::find(vec_unique.begin(), vec_unique.end(), std::make_pair(index0, index2)) == vec_unique.end()
      && std::find(vec_unique.begin(), vec_unique.end(), std::make_pair(index2, index0)) == vec_unique.end())
    {
      RelativeCameraMotion(d._R[index0], d._t[index0], d._R[index2], d._t[index2], &Rrel, &trel);
      vec_relativeRotEstimate.push_back(RelativeRotation(index0, index2, Rrel, 1));
      vec_unique.push_back(make_pair(index0, index2));
    }
  }

  // Add 2 outliers rotations between (0->1), (2->3)
  // (use a smaller weight since those rotations are less accurate)
  for (size_t i = 0; i < vec_relativeRotEstimate.size(); ++i)
  {
    if( vec_relativeRotEstimate[i].i == 0 && vec_relativeRotEstimate[i].j == 1)
      vec_relativeRotEstimate[i] = RelativeRotation(0, 1, RotationAroundX(D2R(0.1)), 0.5);
    if( vec_relativeRotEstimate[i].i == 2 && vec_relativeRotEstimate[i].j == 3)
      vec_relativeRotEstimate[i] = RelativeRotation(2, 3, RotationAroundX(D2R(0.6)), 0.5);
  }

  //- Solve the global rotation estimation problem :
  Matrix3x3Arr vec_globalR(iNviews);
  vec_globalR = d._R;
  size_t nMainViewID = 0;
  std::vector<bool> vec_inliers;
  bool bTest = GlobalRotationsRobust(vec_relativeRotEstimate, vec_globalR, nMainViewID, 0.0f, &vec_inliers);
  EXPECT_TRUE(bTest);

  std::cout << "Inliers: " << std::endl;
  std::copy(vec_inliers.begin(), vec_inliers.end(), ostream_iterator<bool>(std::cout, " "));
  std::cout << std::endl;

  // Check inlier list
  CHECK(std::accumulate(vec_inliers.begin(), vec_inliers.end(), 0) == 8);
  // Check outlier have been found
  CHECK(vec_inliers[0]== 0);
  CHECK(vec_inliers[3] == 0);

  // Remove outliers and refine
  RelativeRotations vec_relativeRotEstimateTemp;
  for (size_t i = 0; i < vec_inliers.size(); ++i)
  {
    if( vec_inliers[i] == 1)
      vec_relativeRotEstimateTemp.push_back(vec_relativeRotEstimate[i]);
  }
  vec_relativeRotEstimate.swap(vec_relativeRotEstimateTemp);
  EXPECT_TRUE( GlobalRotationsRobust(vec_relativeRotEstimate, vec_globalR, nMainViewID, 0.0f, &vec_inliers));

  // Check that the loop is closing
  Mat3 rotCumul = vec_globalR[0];
  for (size_t i = 1; i < iNviews; ++i)
  {
    rotCumul*= vec_globalR[i];
  }

  EXPECT_MATRIX_NEAR(Mat3::Identity(), rotCumul, 1e-4);
  // Fix the sign of the rotations (put the global rotation in the same rotation axis as GT)
  if ( SIGN(vec_globalR[0](0,0)) != SIGN( d._R[0](0,0) ))
  {
    Mat3 Rrel;
    Vec3 trel;
    RelativeCameraMotion(vec_globalR[0], Vec3::Zero(), d._R[0], Vec3::Zero(), &Rrel, &trel);
    for (size_t i = 0; i < iNviews; ++i)
      vec_globalR[i] *= Rrel;
  }

  // Check that each global rotations is near the true ones
  for (size_t i = 0; i < iNviews; ++i)
  {
    EXPECT_NEAR(0.0, FrobeniusDistance(d._R[i], vec_globalR[i]), 1e-8);
  }
}

/*
template<typename TYPE, int N>
inline REAL ComputePSNR(const Eigen::Matrix<REAL, N,1>& x0, const Eigen::Matrix<REAL, N,1>& x)
{
  REAL ret = std::numeric_limits<REAL>::infinity();

  REAL err = (x0 - x).squaredNorm() / N;

  TYPE max1 = 0;
  TYPE max2 = 0;
  for (unsigned i=0; i<N; ++i) {
    max1 = std::max(max1, x0(i));
    max2 = std::max(max2, x(i));
  }
  TYPE maxBoth = std::max(max1, max2);

  if ((maxBoth<1e-8) && err > 0)
    ret = 0;
  else if (err > 0)
    ret = 10.0 * std::log10(static_cast<REAL>(maxBoth*maxBoth) / err);

  return ret;
} // ComputePSNR

bool TestRunRobustRegressionL1PD()
{
  typedef Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  typedef Eigen::Matrix<REAL, Eigen::Dynamic, 1> Vector;

  // source length
  const unsigned N = 256;
  // codeword length
  const unsigned M = 4*N;
  // number of perturbations
  const unsigned T = int(0.2f*M);
  // coding matrix
  const Matrix G = Matrix::Random(M,N);
  // source word
  const Vector x = Vector::Random(N);
  // code word
  Vector y = G*x;
  // channel: perturb T randomly chosen entries
  for (unsigned i=0; i<T; ++i)
    y((int) ((rand()/RAND_MAX)* (M-1))) = Vector::Random(1)(0);
  // recover
  Vector x0 = (G.transpose()*G).inverse()*G.transpose()*y;

  Vector& xp(x0);
  RobustRegressionL1PD(G, y, x0, 1e-4, 30);

  REAL psnr = ComputePSNR<REAL,N>(*((Eigen::Matrix<REAL, N,1>*)x.data()), *((Eigen::Matrix<REAL, N,1>*)xp.data()));
  bool bPassed = (psnr > 100);
  return bPassed;
}

bool TestRobustRegressionL1PD()
{
  unsigned nIters = 12;
  unsigned nTotalCorrectSolutions = 0;
  for (unsigned i=0; i<nIters; ++i) {
    if (TestRunRobustRegressionL1PD())
      nTotalCorrectSolutions++;
  }
  const bool bPassed = (nTotalCorrectSolutions >= nIters);
  std::cout
    << "Test robust regression " << (bPassed?"passed":"FAILED")
    << " : " << nIters << " " << nTotalCorrectSolutions
    << " iterations (" << (float)nTotalCorrectSolutions/nIters << " correct) "
    << std::endl;
  return bPassed;

}

TEST ( rotation_averaging, RobustRegressionL1PD)
{
  TestRobustRegressionL1PD();
}*/

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

