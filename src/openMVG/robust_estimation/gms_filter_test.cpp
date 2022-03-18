// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/robust_estimation/gms_filter.hpp"

#include "testing/testing.h"
#include <random>

using namespace openMVG;
using namespace openMVG::robust;
using namespace std;

TEST(GMSFilter, borderCases)
{
  std::vector<Eigen::Vector2f> vec_points;
  matching::IndMatches matches;

  robust::GMSFilter gms(
    vec_points, {0, 0},
    vec_points, {0 ,0},
    matches);
  const bool with_scale_invariance = false;
  const bool with_rotation_invariance = false;
  std::vector<bool> inlier_flags;
  EXPECT_EQ(0, gms.GetInlierMask(inlier_flags,
                                 with_scale_invariance,
                                 with_rotation_invariance));
}

void GenerateCorrespondingPoints
(
  const int image_size,
  const int border_size,
  const int nb_point,
  std::vector<Eigen::Vector2f> & vec_point_left,
  std::vector<Eigen::Vector2f> & vec_point_right,
  matching::IndMatches & matches
)
{
  const std::pair<int, int> kImageSize = {image_size, image_size};

  // Setup a uniform int distribution to compute some random point position in
  // the image domain.
  std::mt19937 random_generator(std::mt19937::default_seed);
  std::uniform_int_distribution<int> dist(0 + border_size, image_size - border_size);

  vec_point_left.resize(nb_point);
  vec_point_right.resize(nb_point);
  matches.resize(nb_point);
  for (int i = 0; i < nb_point; ++i)
  {
    vec_point_left[i] << dist(random_generator), dist(random_generator);
    vec_point_right[i] = vec_point_left[i];
    matches[i] = {static_cast<IndexT>(i), static_cast<IndexT>(i)};
  }
}

TEST(GMSFilter, NoFlowField)
{
  const int kImageSize = 200;
  const int kBorder = 20;
  // Generate a lot of point (GMS is working only if dense motion statistic can be computed)
  const int kNbPoints = (kImageSize * kImageSize) * 0.1;
  std::vector<Eigen::Vector2f> vec_point_left(kNbPoints), vec_point_right(kNbPoints);
  matching::IndMatches matches;
  GenerateCorrespondingPoints(
    kImageSize, kBorder, kNbPoints, vec_point_left, vec_point_right, matches);

  robust::GMSFilter gms(
    vec_point_left,  {kImageSize, kImageSize},
    vec_point_right, {kImageSize, kImageSize},
    matches);
  const bool with_scale_invariance = false;
  const bool with_rotation_invariance = false;
  std::vector<bool> inlier_flags;
  const int num_inliers = gms.GetInlierMask(inlier_flags, with_scale_invariance, with_rotation_invariance);

  EXPECT_EQ(kNbPoints, num_inliers);
  EXPECT_EQ(kNbPoints, inlier_flags.size());
  EXPECT_EQ(kNbPoints, std::count(inlier_flags.cbegin(), inlier_flags.cend(), true));
}

TEST(GMSFilter, FixFlowField)
{
  const int kImageSize = 200;
  const int kBorder = 20;
  // Generate a lot of point (GMS is working only if dense motion statistic can be computed)
  const int kNbPoints = (kImageSize * kImageSize) * 0.1;
  std::vector<Eigen::Vector2f> vec_point_left(kNbPoints), vec_point_right(kNbPoints);
  matching::IndMatches matches;
  GenerateCorrespondingPoints(
    kImageSize, kBorder, kNbPoints, vec_point_left, vec_point_right, matches);

  // generate a fix motion for the right points
  for (auto & pt : vec_point_right)
  {
    pt += Vec2f(kBorder / 2, kBorder / 4);
  }

  robust::GMSFilter gms(
    vec_point_left,  {kImageSize,kImageSize},
    vec_point_right, {kImageSize,kImageSize},
    matches);
  const bool with_scale_invariance = false;
  const bool with_rotation_invariance = false;
  std::vector<bool> inlier_flags;
  const int num_inliers = gms.GetInlierMask(inlier_flags, with_scale_invariance, with_rotation_invariance);

  CHECK(num_inliers == kNbPoints);
  EXPECT_EQ(kNbPoints, inlier_flags.size());
}

TEST(GMSFilter, TwoFlowField)
{
  //--
  // GMS is capable of keeping various motion fields.
  // Since some point can be on the grid border some point can be misclassified,
  //  so the detection rate will not be 100%.
  //--
  const int kImageSize = 200;
  const int kBorder = 20;
  // Generate a lot of point (GMS is working only if dense motion statistic can be computed)
  const int kNbPoints = (kImageSize * kImageSize) * 0.1;
  std::vector<Eigen::Vector2f> vec_point_left(kNbPoints), vec_point_right(kNbPoints);
  matching::IndMatches matches;
  GenerateCorrespondingPoints(
    kImageSize, kBorder, kNbPoints, vec_point_left, vec_point_right, matches);

  // generate a fix motion for the right points
  for (int i = 0; i < kNbPoints/2; ++i)
  {
    vec_point_right[i] += Vec2f(kBorder / 2, kBorder / 4);
  }
  for (int i =  kNbPoints/2; i < kNbPoints; ++i)
  {
    vec_point_right[i] += Vec2f(kBorder / 4, kBorder / 2);
  }

  robust::GMSFilter gms(
    vec_point_left,  {kImageSize,kImageSize},
    vec_point_right, {kImageSize,kImageSize},
    matches);
  const bool with_scale_invariance = false;
  const bool with_rotation_invariance = false;
  std::vector<bool> inlier_flags;
  const int num_inliers = gms.GetInlierMask(inlier_flags, with_scale_invariance, with_rotation_invariance);

  CHECK(num_inliers > .98 * kNbPoints);
  EXPECT_EQ(false, inlier_flags.empty());
}

TEST(GMSFilter, OutlierDetection)
{
  const int kImageSize = 200;
  const int kBorder = 20;
  // Generate a lot of point (GMS is working only if dense motion statistic can be computed)
  const int kNbPoints = (kImageSize * kImageSize) * 0.1;
  std::vector<Eigen::Vector2f> vec_point_left(kNbPoints), vec_point_right(kNbPoints);
  matching::IndMatches matches;
  GenerateCorrespondingPoints(
    kImageSize, kBorder, kNbPoints, vec_point_left, vec_point_right, matches);

  // generate some invalid correspondences (6 outliers)
  for (const int i : {0,1,2,3,4,5})
  {
    matches[i].j_ = matches[i + 5].j_;
  }

  robust::GMSFilter gms(
    vec_point_left,  {kImageSize,kImageSize},
    vec_point_right, {kImageSize,kImageSize},
    matches);
  const bool with_scale_invariance = false;
  const bool with_rotation_invariance = false;
  std::vector<bool> inlier_flags;
  const int num_inliers = gms.GetInlierMask(inlier_flags, with_scale_invariance, with_rotation_invariance);

  // Assert that 6 outliers have been detected
  EXPECT_EQ(kNbPoints - 6, num_inliers);
  EXPECT_EQ(kNbPoints, inlier_flags.size());
  EXPECT_EQ(kNbPoints - 6, std::count(inlier_flags.cbegin(), inlier_flags.cend(), true));
  EXPECT_EQ(6, std::count(inlier_flags.cbegin(), inlier_flags.cbegin() + 6, false));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
