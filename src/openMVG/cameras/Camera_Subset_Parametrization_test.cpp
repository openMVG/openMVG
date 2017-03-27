// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/cameras.hpp"
using namespace openMVG;
using namespace openMVG::cameras;

#include "testing/testing.h"

//-----------------
// Test summary:
//-----------------
// - Test subset parametrization:
//   Parameter index filter to know which parameters must be held constant for
//    a given parametrization.
//-----------------
TEST(Cameras, subset_parametrization_pinhole)
{
  Pinhole_Intrinsic cam;
  // Parametrized by 3 parameters:
  // 0    -> the focal,
  // 1,2  -> the principal point
  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::NONE);
    EXPECT_EQ(3, param_filter.size());
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::ADJUST_ALL);
    EXPECT_EQ(0, param_filter.size());
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH);
    EXPECT_EQ(2, param_filter.size()); // PPX & PPY must be held constant
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 1) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 2) );
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT);
    EXPECT_EQ(1, param_filter.size()); // FOCAL must be held constant
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 0) );
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::ADJUST_DISTORTION);
    // Since the camera model does not have distortion, it is equivalent to NONE
    EXPECT_EQ(3, param_filter.size());
  }

  {
    const std::vector<int> param_filter = cam.subsetParameterization(
        Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH |
        Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT);
    // Equivalent to adjust_all
    EXPECT_EQ(0, param_filter.size()); // This camera model does not have distortion
  }

  {
    const std::vector<int> param_filter = cam.subsetParameterization(
      Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT |
      Intrinsic_Parameter_Type::ADJUST_DISTORTION);
    // Must be equivalent to Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT
    EXPECT_EQ(1, param_filter.size()); // FOCAL must be held constant
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 0) );
  }

  {
    const std::vector<int> param_filter = cam.subsetParameterization(
      Intrinsic_Parameter_Type::NONE);
    // No parameters are held constant
    EXPECT_EQ(3, param_filter.size());
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 0) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 1) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 2) );
  }
}

TEST(Cameras, subset_parametrization_pinhole_radial_k1)
{
  const Pinhole_Intrinsic_Radial_K1 cam;
  // Parametrized by 4 parameters:
  // 0    -> the focal,
  // 1,2  -> the principal point
  // 3    -> the distortion coefficient

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::NONE);
    EXPECT_EQ(4, param_filter.size());
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::ADJUST_ALL);
    EXPECT_EQ(0, param_filter.size());
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::ADJUST_DISTORTION);
    EXPECT_EQ(3, param_filter.size());
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 0) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 1) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 2) );
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(
        Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH |
        Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT);
    EXPECT_EQ(1, param_filter.size());
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 3) );
  }
}

TEST(Cameras, subset_parametrization_pinhole_radial_k3)
{
  const Pinhole_Intrinsic_Radial_K3 cam;
  // Parametrized by 6 parameters:
  // 0    -> the focal,
  // 1,2  -> the principal point
  // 3,4,5-> the distortion coefficient

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::NONE);
    EXPECT_EQ(6, param_filter.size());
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::ADJUST_ALL);
    EXPECT_EQ(0, param_filter.size());
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::ADJUST_DISTORTION);
    EXPECT_EQ(3, param_filter.size());
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 0) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 1) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 2) );
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(
        Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH |
        Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT);
    EXPECT_EQ(3, param_filter.size());
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 3) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 4) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 5) );
  }
}

TEST(Cameras, subset_parametrization_pinhole_brown)
{
  const Pinhole_Intrinsic_Brown_T2 cam;
  // Parametrized by 8 parameters:
  // 0          -> the focal,
  // 1,2        -> the principal point
  // 3,4,5,6,7  -> the distortion coefficient

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::NONE);
    EXPECT_EQ(8, param_filter.size());
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::ADJUST_ALL);
    EXPECT_EQ(0, param_filter.size());
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::ADJUST_DISTORTION);
    EXPECT_EQ(3, param_filter.size());
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 0) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 1) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 2) );
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(
        Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH |
        Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT);
    EXPECT_EQ(5, param_filter.size());
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 3) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 4) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 5) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 6) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 7) );
  }

    {
    const std::vector<int> param_filter =
      cam.subsetParameterization(
        Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH |
        Intrinsic_Parameter_Type::ADJUST_DISTORTION);
    EXPECT_EQ(2, param_filter.size());
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 1) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 2) );
  }
}

TEST(Cameras, subset_parametrization_fisheye)
{
  const Pinhole_Intrinsic_Fisheye cam;
  // Parametrized by 7 parameters:
  // 0          -> the focal,
  // 1,2        -> the principal point
  // 3,4,5,6  -> the distortion coefficient

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::NONE);
    EXPECT_EQ(7, param_filter.size());
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::ADJUST_ALL);
    EXPECT_EQ(0, param_filter.size());
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(Intrinsic_Parameter_Type::ADJUST_DISTORTION);
    EXPECT_EQ(3, param_filter.size());
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 0) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 1) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 2) );
  }

  {
    const std::vector<int> param_filter =
      cam.subsetParameterization(
        Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH |
        Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT);
    EXPECT_EQ(4, param_filter.size());
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 3) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 4) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 5) );
    EXPECT_TRUE( std::count(param_filter.begin(), param_filter.end(), 6) );
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
