// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "rigid_transformation3D_srt.hpp"

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"

#include "CppUnitLite/TestHarness.h"

#include "testing/testing.h"

#include <iostream>

using namespace openMVG;
using namespace openMVG::geometry;
using namespace std;

TEST(SRT_precision, Experiment_ScaleOnly)
{

  const std::size_t nbPoints = 10;
  Mat x1 = Mat::Random(3, nbPoints);
  Mat x2 = x1;

  const double scale = 2;
  Mat3 rot = Mat3::Identity();
  Vec3 t(0, 0, 0);

  for(std::size_t i = 0; i < nbPoints; ++i)
  {
    Vec3 pt = x1.col(i);
    x2.col(i) = (scale * rot * pt + t);
  }

  // Compute the Similarity transform
  double Sc = 0;
  Mat3 Rc;
  Vec3 tc;
  FindRTS(x1, x2, Sc, tc, Rc);
  Refine_RTS(x1, x2, Sc, tc, Rc);

  std::cout << "\n"
          << "Scale " << Sc << "\n"
          << "Rot \n" << Rc << "\n"
          << "t " << tc.transpose();
}

TEST(SRT_precision, Experiment_ScaleAndRot)
{

  const std::size_t nbPoints = 10;
  Mat x1 = Mat::Random(3, nbPoints);
  Mat x2 = x1;

  const double scale = 2;
  Mat3 rot = (Eigen::AngleAxis<double>(.2, Vec3::UnitX())
          * Eigen::AngleAxis<double>(.3, Vec3::UnitY())
          * Eigen::AngleAxis<double>(.6, Vec3::UnitZ())).toRotationMatrix();
  Vec3 t(0, 0, 0);

  for(std::size_t i = 0; i < nbPoints; ++i)
  {
    Vec3 pt = x1.col(i);
    x2.col(i) = (scale * rot * pt + t);
  }

  // Compute the Similarity transform
  double Sc = 0;
  Mat3 Rc;
  Vec3 tc;
  FindRTS(x1, x2, Sc, tc, Rc);
  Refine_RTS(x1, x2, Sc, tc, Rc);

  std::cout << "\n"
          << "Scale " << Sc << "\n"
          << "Rot \n" << Rc << "\n"
          << "t " << tc.transpose();

  std::cout << "\nGT\n"
          << "Scale " << scale << "\n"
          << "Rot \n" << rot << "\n"
          << "t " << t.transpose();
}

TEST(SRT_precision, Experiment_ScaleRotTranslation)
{

  const std::size_t nbPoints = 10;
  Mat x1 = Mat::Random(3, nbPoints);
  Mat x2 = x1;

  const double scale = 2;
  Mat3 rot = (Eigen::AngleAxis<double>(.2, Vec3::UnitX())
          * Eigen::AngleAxis<double>(.3, Vec3::UnitY())
          * Eigen::AngleAxis<double>(.6, Vec3::UnitZ())).toRotationMatrix();
  Vec3 t(0.5, -0.3, .38);

  for(std::size_t i = 0; i < nbPoints; ++i)
  {
    Vec3 pt = x1.col(i);
    x2.col(i) = (scale * rot * pt + t);
  }

  // Compute the Similarity transform
  double Sc = 0;
  Mat3 Rc;
  Vec3 tc;
  FindRTS(x1, x2, Sc, tc, Rc);
  Refine_RTS(x1, x2, Sc, tc, Rc);

  std::cout << "\n"
          << "Scale " << Sc << "\n"
          << "Rot \n" << Rc << "\n"
          << "t " << tc.transpose();

  std::cout << "\nGT\n"
          << "Scale " << scale << "\n"
          << "Rot \n" << rot << "\n"
          << "t " << t.transpose();
}

TEST(SRT_precision, ACRANSAC_noNoise)
{
  const std::size_t nbPoints = 100;
  Mat x1 = Mat::Random(3, nbPoints);
  Mat x2 = x1;

  const double scale = 2;
  Mat3 rot = (Eigen::AngleAxis<double>(.2, Vec3::UnitX())
          * Eigen::AngleAxis<double>(.3, Vec3::UnitY())
          * Eigen::AngleAxis<double>(.6, Vec3::UnitZ())).toRotationMatrix();
  Vec3 t(0.5, -0.3, .38);

  for(std::size_t i = 0; i < nbPoints; ++i)
  {
    const Vec3 &pt = x1.col(i);
    x2.col(i) = (scale * rot * pt + t);
  }

  // Compute the Similarity transform
  double Sc;
  Mat3 Rc;
  Vec3 tc;

  std::vector<std::size_t> vec_inliers;
  const bool result = ACRansac_FindRTS(x1, x2, Sc, tc, Rc, vec_inliers, true);

  EXPECT_TRUE(result);
  EXPECT_TRUE(vec_inliers.size() == nbPoints);

  std::cout << "\n"
          << "Scale " << Sc << "\n"
          << "Rot \n" << Rc << "\n"
          << "t " << tc.transpose();

  std::cout << "\nGT\n"
          << "Scale " << scale << "\n"
          << "Rot \n" << rot << "\n"
          << "t " << t.transpose();

  EXPECT_NEAR(scale, Sc, 1e-9);
  
  Mat4 RTS;
  composeRTS(Sc, tc, Rc, RTS);

  for(std::size_t i = 0; i < nbPoints; ++i)
  {
    const double error = geometry::RTSSquaredResidualError::Error(RTS, x1.col(i), x2.col(i));
    EXPECT_NEAR(error, 0.0, 1e-9);
  }
}

TEST(SRT_precision, ACRANSAC_noiseByShuffling)
{
  // it generates some points x1, it only generates the corresponding 
  // transformed points x2 for nbPoints-nbShuffles of them while the rest
  // are again taken randomly in order to generate outliers
  const std::size_t nbPoints = 100;
  const std::size_t nbNoisy = 50;

  Mat x1 = Mat::Random(3, nbPoints);
  Mat x2 = Mat::Random(3, nbPoints);

  const double scale = 2;
  Mat3 rot = (Eigen::AngleAxis<double>(.2, Vec3::UnitX())
          * Eigen::AngleAxis<double>(.3, Vec3::UnitY())
          * Eigen::AngleAxis<double>(.6, Vec3::UnitZ())).toRotationMatrix();
  Vec3 t(0.5, -0.3, .38);

  for(std::size_t i = 0; i < nbPoints - nbNoisy; ++i)
  {
    const Vec3 &pt = x1.col(i);
    x2.col(i) = (scale * rot * pt + t);
  }

  // Compute the Similarity transform
  double Sc;
  Mat3 Rc;
  Vec3 tc;

  std::vector<std::size_t> vec_inliers;
  const bool result = ACRansac_FindRTS(x1, x2, Sc, tc, Rc, vec_inliers, true);

  std::cout << "\n"
          << "Scale " << Sc << "\n"
          << "Rot \n" << Rc << "\n"
          << "t " << tc.transpose();

  std::cout << "\nGT\n"
          << "Scale " << scale << "\n"
          << "Rot \n" << rot << "\n"
          << "t " << t.transpose() << "\n";
  
  EXPECT_TRUE(result);
  // all the points must be inliers (no noise)
  const std::size_t nbInliers = vec_inliers.size();
  EXPECT_TRUE(nbInliers == nbPoints - nbNoisy);

  Mat inliers1 = Mat3X(3, nbInliers);
  Mat inliers2 = Mat3X(3, nbInliers);

  for(std::size_t i = 0; i < nbInliers; ++i)
  {
    inliers1.col(i) = x1.col(vec_inliers[i]);
    inliers2.col(i) = x2.col(vec_inliers[i]);
  }

  // check scale
  EXPECT_NEAR(scale, Sc, 1e-9);
  
  Mat4 RTS;
  composeRTS(Sc, tc, Rc, RTS);

  // check the residuals for the inliers
  for(std::size_t i = 0; i < nbInliers; ++i)
  {
    const double error = geometry::RTSSquaredResidualError::Error(RTS, inliers1.col(i), inliers2.col(i));
    EXPECT_NEAR(error, 0.0, 1e-9);
  }
}

/* ************************************************************************* */
int main()
{
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */
