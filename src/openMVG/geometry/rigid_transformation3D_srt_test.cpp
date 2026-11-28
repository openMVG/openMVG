// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/rigid_transformation3D_srt.hpp"
#include "openMVG/geometry/Similarity3_Kernel.hpp"
#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"
#include <iostream>

using namespace openMVG;
using namespace openMVG::geometry;

TEST(SRT_precision, CoincidentSourcePointsAreRejected) {
  Mat source(3, 4);
  source.colwise() = Vec3(2, -3, 5);
  Mat target(3, 4);
  target << 0, 1, 0, 1,
            0, 0, 1, 1,
            0, 0, 0, 1;

  double scale = 1;
  Vec3 translation = Vec3::Zero();
  Mat3 rotation = Mat3::Identity();
  EXPECT_FALSE(FindRTS(source, target, &scale, &translation, &rotation));

  // The robust-estimation solver must not retain a non-finite hypothesis.
  std::vector<Similarity3> models;
  kernel::Similarity3Solver::Solve(source, target, &models);
  EXPECT_TRUE(models.empty());
}

TEST(SRT_precision, MinimalNonCollinearPointsRemainValid) {
  Mat source(3, 3);
  source << 0, 1, 0,
            0, 0, 1,
            0, 0, 0;
  const Mat3 expected_rotation =
    Eigen::AngleAxis<double>(0.4, Vec3::UnitZ()).toRotationMatrix();
  const Vec3 expected_translation(2, -3, 5);
  const Mat target = (2.5 * expected_rotation * source).colwise() + expected_translation;
  double scale;
  Vec3 translation;
  Mat3 rotation;
  EXPECT_TRUE(FindRTS(source, target, &scale, &translation, &rotation));
  EXPECT_NEAR(2.5, scale, 1e-12);
  EXPECT_MATRIX_NEAR(expected_rotation, rotation, 1e-12);
  EXPECT_MATRIX_NEAR(expected_translation, translation, 1e-12);
}

TEST(SRT_precision, Experiment_ScaleOnly) {

  int nbPoints = 10;
  Mat x1 = Mat::Random(3,nbPoints);
  Mat x2 = x1;

  double scale = 2;
  Mat3 rot = Mat3::Identity();
  Vec3 t(0,0,0);

  for (int i=0; i < nbPoints; ++i)
  {
    Vec3 pt = x1.col(i);
    x2.col(i) = (scale * rot * pt + t);
  }

  // Compute the Similarity transform
  double Sc;
  Mat3 Rc;
  Vec3 tc;
  FindRTS(x1, x2, &Sc, &tc, &Rc);
  Refine_RTS(x1,x2,&Sc,&tc,&Rc);

  std::cout << "\n"
    << "Scale " << Sc << "\n"
    << "Rot \n" << Rc << "\n"
    << "t " << tc.transpose();
}

TEST(SRT_precision, Experiment_ScaleAndRot) {

  int nbPoints = 10;
  Mat x1 = Mat::Random(3,nbPoints);
  Mat x2 = x1;

  double scale = 2;
  Mat3 rot = (Eigen::AngleAxis<double>(.2, Vec3::UnitX())
      * Eigen::AngleAxis<double>(.3, Vec3::UnitY())
      * Eigen::AngleAxis<double>(.6, Vec3::UnitZ())).toRotationMatrix();
  Vec3 t(0,0,0);

  for (int i=0; i < nbPoints; ++i)
  {
    Vec3 pt = x1.col(i);
    x2.col(i) = (scale * rot * pt + t);
  }

  // Compute the Similarity transform
  double Sc;
  Mat3 Rc;
  Vec3 tc;
  FindRTS(x1, x2, &Sc, &tc, &Rc);
  Refine_RTS(x1,x2,&Sc,&tc,&Rc);

  std::cout << "\n"
    << "Scale " << Sc << "\n"
    << "Rot \n" << Rc << "\n"
    << "t " << tc.transpose();

  std::cout << "\nGT\n"
    << "Scale " << scale << "\n"
    << "Rot \n" << rot << "\n"
    << "t " << t.transpose();
}

TEST(SRT_precision, Experiment_ScaleRotTranslation) {

  int nbPoints = 10;
  Mat x1 = Mat::Random(3,nbPoints);
  Mat x2 = x1;

  double scale = 2;
  Mat3 rot = (Eigen::AngleAxis<double>(.2, Vec3::UnitX())
      * Eigen::AngleAxis<double>(.3, Vec3::UnitY())
      * Eigen::AngleAxis<double>(.6, Vec3::UnitZ())).toRotationMatrix();
  Vec3 t(0.5,-0.3,.38);

  for (int i=0; i < nbPoints; ++i)
  {
    Vec3 pt = x1.col(i);
    x2.col(i) = (scale * rot * pt + t);
  }

  // Compute the Similarity transform
  double Sc;
  Mat3 Rc;
  Vec3 tc;
  FindRTS(x1, x2, &Sc, &tc, &Rc);
  Refine_RTS(x1,x2,&Sc,&tc,&Rc);

  std::cout << "\n"
    << "Scale " << Sc << "\n"
    << "Rot \n" << Rc << "\n"
    << "t " << tc.transpose();

  std::cout << "\nGT\n"
    << "Scale " << scale << "\n"
    << "Rot \n" << rot << "\n"
    << "t " << t.transpose();
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
