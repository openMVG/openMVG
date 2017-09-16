// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/rigid_transformation3D_srt.hpp"
#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"
#include <iostream>

using namespace openMVG;
using namespace openMVG::geometry;
using namespace std;

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
